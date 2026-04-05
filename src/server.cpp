#ifdef COOLSOLVE_GUI

#include "coolsolve/server.h"
#include "coolsolve/runner.h"
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/solver.h"
#include "coolsolve/solution_checker.h"
#include "coolsolve/constants.h"
#include "coolsolve/fluids.h"
#include "coolsolve/variable_inference.h"
#include "coolsolve/units.h"
#include "coolsolve/latex_report.h"

#include <httplib.h>
#include "CoolProp.h"
#include <nlohmann/json.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <mutex>
#include <thread>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include <condition_variable>
#include <random>
#include <set>
#include <numeric>
#include <algorithm>

#ifdef COOLSOLVE_EMBEDDED_ASSETS
#include "embedded_assets.h"
#endif

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace coolsolve {

// ============================================================================
// C++17 helpers (no starts_with/ends_with)
// ============================================================================
static bool endsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}
static bool startsWith(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

// ============================================================================
// Base64 decoder (for plot images uploaded as data-URIs)
// ============================================================================
static std::string base64Decode(const std::string& encoded) {
    static constexpr unsigned char T[256] = {
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,62,64,64,64,63,
        52,53,54,55,56,57,58,59,60,61,64,64,64,64,64,64,
        64, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,
        15,16,17,18,19,20,21,22,23,24,25,64,64,64,64,64,
        64,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
        41,42,43,44,45,46,47,48,49,50,51,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,
        64,64,64,64,64,64,64,64,64,64,64,64,64,64,64,64
    };
    std::string out;
    out.reserve(encoded.size() * 3 / 4);
    unsigned int buf = 0;
    int bits = 0;
    for (unsigned char c : encoded) {
        if (T[c] >= 64) continue;  // skip whitespace, padding, illegal
        buf = (buf << 6) | T[c];
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<char>((buf >> bits) & 0xFF));
        }
    }
    return out;
}

// ============================================================================
// MIME type helper
// ============================================================================
static std::string getMimeType(const std::string& path) {
    if (endsWith(path, ".html")) return "text/html";
    if (endsWith(path, ".css"))  return "text/css";
    if (endsWith(path, ".js"))   return "application/javascript";
    if (endsWith(path, ".mjs"))  return "application/javascript";
    if (endsWith(path, ".json")) return "application/json";
    if (endsWith(path, ".svg"))  return "image/svg+xml";
    if (endsWith(path, ".png"))  return "image/png";
    if (endsWith(path, ".ico"))  return "image/x-icon";
    if (endsWith(path, ".woff")) return "font/woff";
    if (endsWith(path, ".woff2")) return "font/woff2";
    if (endsWith(path, ".ttf")) return "font/ttf";
    if (endsWith(path, ".map")) return "application/json";
    return "application/octet-stream";
}

// ============================================================================
// Inferred variable data (extracted from IR after solve)
// ============================================================================
struct InferredVariable {
    std::string name;
    double value = 0.0;
    std::string inferredProperty;  // "T", "P", "H", "S", "D", "Q", etc.
    std::string inferredFluid;     // "Water", "R134a", etc.
    std::string units;             // "[K]", "[Pa]", etc.
    bool isArray = false;
};

// ============================================================================
// Snapshot for Back button (one-level undo for model loading)
// ============================================================================
struct SessionSnapshot {
    std::string eescodeContent;
    std::string initialsContent;
    std::string solContent;
    std::string confContent;
    std::string openFilePath;
    std::string modelName;
    std::string debugDir;
    bool hasResult = false;
    SolveResult lastResult;
    CoolSolveRunner::PipelineTiming lastTiming;
    std::vector<InferredVariable> inferredVariables;
    std::vector<std::string> modelFluids;
    std::string latexReportContent;
    bool latexReportAvailable = false;
};

// ============================================================================
// Session state
// ============================================================================
struct Session {
    std::string eescodeContent;                 // Current .eescode source
    std::string initialsContent;                // Current .initials content
    std::string solContent;                     // Last .sol output
    std::string confContent;                    // coolsolve.conf content
    
    std::string openFilePath;                   // Path of currently open .eescode file
    std::string modelName;                      // User-facing model name (for ZIP filename and page title)
    
    fs::path tempDir;                               // Session temp directory
    std::string debugDir;                           // Last debug output directory
    
    // Snapshot for Back button
    std::unique_ptr<SessionSnapshot> previousState;
    
    // Solve state
    std::atomic<bool> solving{false};
    std::atomic<bool> solveFinished{true};
    std::atomic<bool> cancelRequested{false};
    
    // Last solve result (protected by mutex)
    std::mutex resultMutex;
    SolveResult lastResult;
    CoolSolveRunner::PipelineTiming lastTiming;
    DiagnosticCollector lastDiagnostics;
    bool hasResult = false;
    
    // Inferred variable data (for thermodynamic diagrams)
    std::vector<InferredVariable> inferredVariables;
    std::vector<std::string> modelFluids;  // unique fluids in the model
    
    // SSE progress events (protected by mutex + condition variable)
    std::mutex progressMutex;
    std::condition_variable progressCV;
    std::vector<std::string> progressEvents;
    
    // Last parametric study result (JSON string, protected by mutex)
    std::mutex parametricMutex;
    std::string lastParametricResult;
    
    // LaTeX report content (generated after solve when enableLatexReport is true)
    std::string latexReportContent;
    bool latexReportAvailable = false;  // true when a LaTeX report was generated for the current solve
    
    // Save current state as snapshot (for Back button)
    void saveSnapshot() {
        auto snap = std::make_unique<SessionSnapshot>();
        snap->eescodeContent = eescodeContent;
        snap->initialsContent = initialsContent;
        snap->solContent = solContent;
        snap->confContent = confContent;
        snap->openFilePath = openFilePath;
        snap->modelName = modelName;
        snap->debugDir = debugDir;
        snap->inferredVariables = inferredVariables;
        snap->modelFluids = modelFluids;
        snap->latexReportContent = latexReportContent;
        snap->latexReportAvailable = latexReportAvailable;
        {
            std::lock_guard<std::mutex> lock(resultMutex);
            snap->hasResult = hasResult;
            if (hasResult) {
                snap->lastResult = lastResult;
                snap->lastTiming = lastTiming;
            }
        }
        previousState = std::move(snap);
    }
    
    // Restore state from snapshot (Back button)
    bool restoreSnapshot() {
        if (!previousState) return false;
        eescodeContent = previousState->eescodeContent;
        initialsContent = previousState->initialsContent;
        solContent = previousState->solContent;
        confContent = previousState->confContent;
        openFilePath = previousState->openFilePath;
        modelName = previousState->modelName;
        debugDir = previousState->debugDir;
        inferredVariables = previousState->inferredVariables;
        modelFluids = previousState->modelFluids;
        latexReportContent = previousState->latexReportContent;
        latexReportAvailable = previousState->latexReportAvailable;
        {
            std::lock_guard<std::mutex> lock(resultMutex);
            hasResult = previousState->hasResult;
            if (previousState->hasResult) {
                lastResult = previousState->lastResult;
                lastTiming = previousState->lastTiming;
            } else {
                lastResult = SolveResult{};
                lastTiming = CoolSolveRunner::PipelineTiming{};
            }
        }
        previousState.reset();
        return true;
    }
    
    void addProgressEvent(const std::string& event) {
        {
            std::lock_guard<std::mutex> lock(progressMutex);
            progressEvents.push_back(event);
        }
        progressCV.notify_all();
    }
    
    std::vector<std::string> consumeProgressEvents() {
        std::lock_guard<std::mutex> lock(progressMutex);
        auto events = std::move(progressEvents);
        progressEvents.clear();
        return events;
    }
    
    // Wait for new events or solve completion, with timeout
    std::vector<std::string> waitForEvents(int timeoutMs = 100) {
        std::unique_lock<std::mutex> lock(progressMutex);
        if (progressEvents.empty()) {
            progressCV.wait_for(lock, std::chrono::milliseconds(timeoutMs));
        }
        auto events = std::move(progressEvents);
        progressEvents.clear();
        return events;
    }
};

// ============================================================================
// Session management for multi-user support
// ============================================================================

class SessionManager {
    std::map<std::string, std::shared_ptr<Session>> sessions_;
    std::mutex mutex_;
public:
    std::shared_ptr<Session> getOrCreate(const std::string& sessionId) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto it = sessions_.find(sessionId);
        if (it != sessions_.end()) return it->second;
        auto s = std::make_shared<Session>();
        s->tempDir = fs::temp_directory_path() / "coolsolve_sessions" / sessionId;
        fs::create_directories(s->tempDir);
        sessions_[sessionId] = s;
        return s;
    }
};

static std::string getCookieValue(const httplib::Request& req, const std::string& name) {
    auto range = req.headers.equal_range("Cookie");
    for (auto it = range.first; it != range.second; ++it) {
        size_t pos = 0;
        const auto& c = it->second;
        while (pos < c.size()) {
            while (pos < c.size() && (c[pos] == ' ' || c[pos] == ';')) pos++;
            auto eq = c.find('=', pos);
            if (eq == std::string::npos) break;
            auto semi = c.find(';', eq);
            if (semi == std::string::npos) semi = c.size();
            if (c.substr(pos, eq - pos) == name)
                return c.substr(eq + 1, semi - eq - 1);
            pos = semi + 1;
        }
    }
    return "";
}

static std::string generateSessionId() {
    static const char hex[] = "0123456789abcdef";
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 15);
    std::string id;
    id.reserve(32);
    for (int i = 0; i < 32; i++) id += hex[dis(gen)];
    return id;
}

// ============================================================================
// ZIP helpers (minimal uncompressed archive implementation)
// ============================================================================

static uint32_t zipCrc32(const uint8_t* data, size_t length) {
    uint32_t crc = 0xFFFFFFFF;
    for (size_t i = 0; i < length; i++) {
        crc ^= data[i];
        for (int j = 0; j < 8; j++)
            crc = (crc >> 1) ^ (0xEDB88320u & (~(crc & 1) + 1));
    }
    return ~crc;
}

static void zipW16(std::string& s, uint16_t v) {
    s += static_cast<char>(v & 0xFF);
    s += static_cast<char>((v >> 8) & 0xFF);
}
static void zipW32(std::string& s, uint32_t v) {
    s += static_cast<char>(v & 0xFF);
    s += static_cast<char>((v >> 8) & 0xFF);
    s += static_cast<char>((v >> 16) & 0xFF);
    s += static_cast<char>((v >> 24) & 0xFF);
}
static uint16_t zipR16(const char* p) {
    return static_cast<uint8_t>(p[0]) | (static_cast<uint8_t>(p[1]) << 8);
}
static uint32_t zipR32(const char* p) {
    return static_cast<uint8_t>(p[0]) | (static_cast<uint8_t>(p[1]) << 8)
         | (static_cast<uint8_t>(p[2]) << 16) | (static_cast<uint8_t>(p[3]) << 24);
}

static std::string createZipBundle(const std::vector<std::pair<std::string, std::string>>& files) {
    std::string out;
    struct E { std::string name; uint32_t crc, size, off; };
    std::vector<E> entries;
    for (const auto& [name, data] : files) {
        E e{name, zipCrc32(reinterpret_cast<const uint8_t*>(data.data()), data.size()),
            static_cast<uint32_t>(data.size()), static_cast<uint32_t>(out.size())};
        zipW32(out, 0x04034b50); zipW16(out, 20); zipW16(out, 0); zipW16(out, 0);
        zipW16(out, 0); zipW16(out, 0);
        zipW32(out, e.crc); zipW32(out, e.size); zipW32(out, e.size);
        zipW16(out, static_cast<uint16_t>(name.size())); zipW16(out, 0);
        out += name; out += data;
        entries.push_back(e);
    }
    uint32_t cdOff = static_cast<uint32_t>(out.size());
    for (const auto& e : entries) {
        zipW32(out, 0x02014b50); zipW16(out, 20); zipW16(out, 20);
        zipW16(out, 0); zipW16(out, 0); zipW16(out, 0); zipW16(out, 0);
        zipW32(out, e.crc); zipW32(out, e.size); zipW32(out, e.size);
        zipW16(out, static_cast<uint16_t>(e.name.size()));
        zipW16(out, 0); zipW16(out, 0); zipW16(out, 0); zipW16(out, 0);
        zipW32(out, 0); zipW32(out, e.off);
        out += e.name;
    }
    uint32_t cdSize = static_cast<uint32_t>(out.size()) - cdOff;
    zipW32(out, 0x06054b50); zipW16(out, 0); zipW16(out, 0);
    zipW16(out, static_cast<uint16_t>(entries.size()));
    zipW16(out, static_cast<uint16_t>(entries.size()));
    zipW32(out, cdSize); zipW32(out, cdOff); zipW16(out, 0);
    return out;
}

static std::map<std::string, std::string> extractZipBundle(const std::string& data) {
    std::map<std::string, std::string> files;
    size_t pos = 0;
    while (pos + 30 <= data.size()) {
        if (zipR32(data.data() + pos) != 0x04034b50) break;
        uint16_t nameLen = zipR16(data.data() + pos + 26);
        uint16_t extraLen = zipR16(data.data() + pos + 28);
        uint32_t compSize = zipR32(data.data() + pos + 18);
        std::string name(data.data() + pos + 30, nameLen);
        pos += 30 + nameLen + extraLen;
        if (pos + compSize <= data.size())
            files[name] = std::string(data.data() + pos, compSize);
        pos += compSize;
    }
    return files;
}

// ============================================================================
// Helper: Read file to string
// ============================================================================
static std::string readFileToString(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return "";
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

// ============================================================================
// Helper: Write string to file
// ============================================================================
static bool writeStringToFile(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    if (!f.is_open()) return false;
    f << content;
    return f.good();
}

// ============================================================================
// Helper: Auto-discover companion files for an .eescode file
// ============================================================================
static void discoverCompanionFiles(Session& session) {
    if (session.openFilePath.empty()) return;
    
    fs::path eesPath(session.openFilePath);
    fs::path dir = eesPath.parent_path();
    std::string stem = eesPath.stem().string();
    
    // .initials
    fs::path initialsPath = dir / (stem + ".initials");
    if (fs::exists(initialsPath)) {
        session.initialsContent = readFileToString(initialsPath.string());
    }
    
    // .sol
    fs::path solPath = dir / (stem + ".sol");
    if (fs::exists(solPath)) {
        session.solContent = readFileToString(solPath.string());
    }
    
    // coolsolve.conf
    fs::path confPath = dir / "coolsolve.conf";
    if (fs::exists(confPath)) {
        session.confContent = readFileToString(confPath.string());
    }
}

// ============================================================================
// Helper: Build JSON from solve result
// ============================================================================
static json solveResultToJSON(const SolveResult& result, const CoolSolveRunner::PipelineTiming& timing,
                              const DiagnosticCollector* runnerDiag = nullptr) {
    json j;
    j["success"] = result.success;
    j["status"] = statusToString(result.status);
    j["errorMessage"] = result.errorMessage;
    j["totalIterations"] = result.totalIterations;
    j["blocksEvaluated"] = result.blocksEvaluated;
    j["totalTimeMs"] = std::chrono::duration_cast<std::chrono::milliseconds>(result.totalTime).count();
    
    // Timing breakdown
    j["timing"] = {
        {"coolprop_warmup_ms", timing.coolprop_warmup_time_ms},
        {"parse_ms", timing.parse_time_ms},
        {"ir_ms", timing.ir_time_ms},
        {"infer_ms", timing.infer_time_ms},
        {"analysis_ms", timing.analysis_time_ms},
        {"solve_ms", timing.solve_time_ms},
        {"total_ms", timing.total_time_ms}
    };
    
    // Variables
    json vars = json::object();
    for (const auto& [name, val] : result.variables) {
        vars[name] = val;
    }
    j["variables"] = vars;
    
    // String variables
    json strVars = json::object();
    for (const auto& [name, val] : result.stringVariables) {
        strVars[name] = val;
    }
    j["stringVariables"] = strVars;
    
    // Per-block results
    json blocks = json::array();
    for (const auto& br : result.blockResults) {
        json blockObj;
        blockObj["id"] = br.id;
        blockObj["success"] = br.success;
        blockObj["status"] = statusToString(br.status);
        blockObj["iterations"] = br.iterations;
        blockObj["maxResidual"] = br.maxResidual;
        blockObj["errorMessage"] = br.errorMessage;
        blocks.push_back(blockObj);
    }
    j["blockResults"] = blocks;
    
    if (!result.detailedError.empty()) {
        j["detailedError"] = result.detailedError;
    }
    
    // Diagnostics from all phases
    if (runnerDiag && runnerDiag->size() > 0) {
        json diagArray = json::array();
        // Summarize C001 CoolProp warnings instead of listing each one
        std::map<std::string, int> c001Counts;
        int totalC001 = 0;
        for (const auto& d : runnerDiag->items()) {
            if (d.code == "C001") {
                c001Counts[d.message]++;
                totalC001++;
                continue;  // Don't add individually
            }
            json dj;
            dj["severity"] = severityToString(d.severity);
            dj["code"] = d.code;
            dj["message"] = d.message;
            dj["source"] = d.source;
            if (d.line > 0) dj["line"] = d.line;
            if (d.column > 0) dj["column"] = d.column;
            diagArray.push_back(dj);
        }
        // Add a single summary diagnostic for C001 warnings
        if (totalC001 > 0) {
            json dj;
            dj["severity"] = "warning";
            dj["code"] = "C001";
            dj["message"] = std::to_string(totalC001) + " CoolProp warning(s) during solving (" 
                          + std::to_string(c001Counts.size()) + " unique)";
            dj["source"] = "evaluator";
            diagArray.push_back(dj);
        }
        j["diagnostics"] = diagArray;
    }
    
    return j;
}

// ============================================================================
// Open browser helper (cross-platform)
// ============================================================================
static void openBrowser(const std::string& url) {
#if defined(__linux__)
    std::string cmd = "xdg-open \"" + url + "\" 2>/dev/null &";
#elif defined(__APPLE__)
    std::string cmd = "open \"" + url + "\"";
#elif defined(_WIN32)
    std::string cmd = "start \"\" \"" + url + "\"";
#else
    std::string cmd = "";
#endif
    if (!cmd.empty()) {
        (void)std::system(cmd.c_str());
    }
}

// ============================================================================
// Server implementation
// ============================================================================
int startServer(const ServerOptions& options) {
    httplib::Server svr;
    SessionManager sessionMgr;
    
    // getSession: extract session from cookie, create if new
    auto getSession = [&sessionMgr](const httplib::Request& req, httplib::Response& res) -> std::shared_ptr<Session> {
        std::string sid = getCookieValue(req, "coolsolve_session");
        if (sid.empty()) {
            sid = generateSessionId();
            res.set_header("Set-Cookie", "coolsolve_session=" + sid + "; Path=/; SameSite=Lax");
        }
        return sessionMgr.getOrCreate(sid);
    };
    
    // CoolProp warmup state
    std::atomic<bool> coolpropReady{false};
    
    // --- CORS middleware (for dev mode with Vite on different port) ---
    svr.set_pre_routing_handler([](const httplib::Request& req, httplib::Response& res) {
        res.set_header("Access-Control-Allow-Origin", "*");
        res.set_header("Access-Control-Allow-Methods", "GET, POST, PUT, DELETE, OPTIONS");
        res.set_header("Access-Control-Allow-Headers", "Content-Type, Authorization");
        
        if (req.method == "OPTIONS") {
            res.status = 204;
            return httplib::Server::HandlerResponse::Handled;
        }
        return httplib::Server::HandlerResponse::Unhandled;
    });
    
    // ================================================================
    // Health check
    // ================================================================
    svr.Get("/api/v1/health", [&](const httplib::Request&, httplib::Response& res) {
        json j = {
            {"status", "ok"},
            {"coolpropReady", coolpropReady.load()}
        };
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Warmup endpoint
    // ================================================================
    svr.Post("/api/v1/warmup", [&](const httplib::Request&, httplib::Response& res) {
        if (!coolpropReady.load()) {
            double ms = warmupCoolProp();
            coolpropReady.store(true);
            json j = {{"warmupMs", ms}};
            res.set_content(j.dump(), "application/json");
        } else {
            json j = {{"warmupMs", 0.0}, {"message", "already warmed up"}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // File Operations
    // ================================================================
    
    // GET /api/v1/files/eescode - Get current .eescode content
    svr.Get("/api/v1/files/eescode", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        json j = {
            {"content", session.eescodeContent},
            {"filePath", session.openFilePath},
            {"modelName", session.modelName}
        };
        res.set_content(j.dump(), "application/json");
    });
    
    // PUT /api/v1/files/eescode - Save .eescode content
    svr.Put("/api/v1/files/eescode", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        try {
            auto body = json::parse(req.body);
            session.eescodeContent = body.value("content", "");
            json j = {{"success", true}};
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // GET /api/v1/files/initials - Get current .initials content
    svr.Get("/api/v1/files/initials", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        json j = {{"content", session.initialsContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // PUT /api/v1/files/initials - Save .initials content  
    svr.Put("/api/v1/files/initials", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        try {
            auto body = json::parse(req.body);
            session.initialsContent = body.value("content", "");
            json j = {{"success", true}};
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // GET /api/v1/files/sol - Get .sol content
    svr.Get("/api/v1/files/sol", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        json j = {{"content", session.solContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // GET /api/v1/files/conf - Get coolsolve.conf content
    svr.Get("/api/v1/files/conf", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        json j = {{"content", session.confContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // PUT /api/v1/files/conf - Save coolsolve.conf content
    svr.Put("/api/v1/files/conf", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        try {
            auto body = json::parse(req.body);
            session.confContent = body.value("content", "");
            json j = {{"success", true}};
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // POST /api/v1/files/open - Open a file from disk (used for examples)
    svr.Post("/api/v1/files/open", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        try {
            auto body = json::parse(req.body);
            std::string filePath = body.value("path", "");
            
            if (filePath.empty() || !fs::exists(filePath)) {
                res.status = 404;
                json j = {{"error", "File not found: " + filePath}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            // Save snapshot for Back button
            session.saveSnapshot();
            
            session.openFilePath = filePath;
            session.modelName = fs::path(filePath).stem().string();
            session.eescodeContent = readFileToString(filePath);
            session.initialsContent.clear();
            session.solContent.clear();
            session.confContent.clear();
            session.debugDir.clear();
            {
                std::lock_guard<std::mutex> lock(session.resultMutex);
                session.hasResult = false;
                session.lastResult = SolveResult{};
            }
            
            // Auto-discover companion files
            discoverCompanionFiles(session);
            
            json j = {
                {"success", true},
                {"filePath", filePath},
                {"modelName", session.modelName},
                {"hasInitials", !session.initialsContent.empty()},
                {"hasSol", !session.solContent.empty()},
                {"hasConf", !session.confContent.empty()}
            };
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // New model — clear session state
    // ================================================================
    svr.Post("/api/v1/new", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        bool hadContent = !session.eescodeContent.empty();
        
        // Save snapshot for Back button (if there was content)
        if (hadContent) {
            session.saveSnapshot();
        }
        
        // Clear all model state
        session.eescodeContent.clear();
        session.initialsContent.clear();
        session.solContent.clear();
        session.confContent.clear();
        session.openFilePath.clear();
        session.modelName.clear();
        
        // Clear debug
        if (!session.debugDir.empty()) {
            std::error_code ec;
            fs::remove_all(session.debugDir, ec);
            session.debugDir.clear();
        }
        
        // Clear solve result
        {
            std::lock_guard<std::mutex> lock(session.resultMutex);
            session.hasResult = false;
            session.lastResult = SolveResult{};
            session.lastTiming = CoolSolveRunner::PipelineTiming{};
        }
        
        // Clear inference data
        session.inferredVariables.clear();
        session.modelFluids.clear();
        
        // Clear LaTeX report
        session.latexReportContent.clear();
        session.latexReportAvailable = false;
        
        json j = {{"success", true}, {"hadContent", hadContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Back — restore previous model state
    // ================================================================
    svr.Post("/api/v1/back", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        if (!session.previousState) {
            res.status = 409;
            json j = {{"error", "No previous model to restore"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        session.restoreSnapshot();
        json j = {{"success", true}, {"modelName", session.modelName}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Model name — get or set the user-facing model name
    // ================================================================
    svr.Get("/api/v1/model-name", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        json j = {{"modelName", session.modelName}};
        res.set_content(j.dump(), "application/json");
    });
    
    svr.Put("/api/v1/model-name", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        try {
            auto body = json::parse(req.body);
            session.modelName = body.value("modelName", "");
            json j = {{"success", true}};
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // Parse endpoint (for real-time editor feedback)
    // ================================================================
    svr.Post("/api/v1/parse", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        try {
            auto body = json::parse(req.body);
            std::string source = body.value("source", session.eescodeContent);
            
            EESParser parser;
            auto parseResult = parser.parse(source);
            
            json j;
            j["success"] = parseResult.success;
            j["equationCount"] = parseResult.equationCount;
            j["totalLines"] = parseResult.totalLines;
            
            json errors = json::array();
            for (const auto& err : parseResult.errors) {
                json errObj;
                errObj["line"] = err.line;
                errObj["column"] = err.column;
                errObj["message"] = err.message;
                errObj["context"] = err.context;
                errors.push_back(errObj);
            }
            j["errors"] = errors;
            
            // If parse succeeded, build IR for variable info
            if (parseResult.success) {
                auto ir = IR::fromAST(parseResult.program);
                
                // Run variable inference to detect fluids used in the model
                try {
                    inferVariables(ir);
                } catch (...) {
                    // Non-fatal: inference failure shouldn't block parse results
                }
                
                // Detect imposed variables: equations of the form `var = number`
                // These are candidates for parametric sensitivity analysis.
                std::set<std::string, CaseInsensitiveLess> imposedVars;
                std::map<std::string, double, CaseInsensitiveLess> imposedValues;
                for (const auto& eq : ir.getEquations()) {
                    if (!eq.lhs || !eq.rhs) continue;
                    // Check: LHS is a plain variable, RHS is a number (or -number)
                    if (eq.lhs->is<Variable>()) {
                        const auto& var = eq.lhs->as<Variable>();
                        if (var.indices.empty()) {  // scalar only
                            double val = 0.0;
                            bool isImposed = false;
                            if (eq.rhs->is<NumberLiteral>()) {
                                val = eq.rhs->as<NumberLiteral>().value;
                                isImposed = true;
                            } else if (eq.rhs->is<UnaryOp>()) {
                                const auto& uop = eq.rhs->as<UnaryOp>();
                                if (uop.op == "-" && uop.operand && uop.operand->is<NumberLiteral>()) {
                                    val = -uop.operand->as<NumberLiteral>().value;
                                    isImposed = true;
                                }
                            }
                            if (isImposed) {
                                imposedVars.insert(var.name);
                                imposedValues[var.name] = val;
                            }
                        }
                    }
                }
                
                // Extract model fluids from inference
                std::set<std::string> fluidSet;
                
                json variables = json::array();
                for (const auto& [name, info] : ir.getVariables()) {
                    if (!Constants::isConstant(name)) {
                        json varObj;
                        varObj["name"] = name;
                        varObj["units"] = info.units;
                        varObj["isArray"] = (name.find('[') != std::string::npos);
                        // Unit source: "code" if from EES source, "inferred" if from
                        // thermodynamic inference, empty otherwise
                        if (!info.units.empty()) {
                            varObj["unitSource"] = "code";
                        } else if (!info.inferredProperty.empty()) {
                            // Infer units from thermodynamic property type
                            std::string inferred;
                            if (info.inferredProperty == "T") inferred = "[°C]";
                            else if (info.inferredProperty == "P") inferred = "[Pa]";
                            else if (info.inferredProperty == "H") inferred = "[J/kg]";
                            else if (info.inferredProperty == "S") inferred = "[J/(kg·K)]";
                            else if (info.inferredProperty == "D") inferred = "[kg/m³]";
                            else if (info.inferredProperty == "Q") inferred = "[-]";
                            else if (info.inferredProperty == "V") inferred = "[m³/kg]";
                            else if (info.inferredProperty == "C") inferred = "[J/(kg·K)]";
                            else if (info.inferredProperty == "L") inferred = "[W/(m·K)]";
                            if (!inferred.empty()) {
                                varObj["units"] = inferred;
                                varObj["unitSource"] = "inferred";
                            } else {
                                varObj["unitSource"] = "";
                            }
                        } else {
                            varObj["unitSource"] = "";
                        }
                        // Is this variable directly imposed (var = constant)?
                        varObj["isImposed"] = imposedVars.count(name) > 0;
                        if (imposedVars.count(name) > 0) {
                            varObj["imposedValue"] = imposedValues[name];
                        }
                        variables.push_back(varObj);
                        if (!info.inferredFluid.empty()) {
                            fluidSet.insert(info.inferredFluid);
                        }
                    }
                }
                session.modelFluids.assign(fluidSet.begin(), fluidSet.end());
                
                j["variables"] = variables;
                j["variableCount"] = ir.getNonConstantVariableCount();
                j["equationCount"] = ir.getEquationCount();
                j["isSquare"] = ir.isSquare();
            }
            
            // Include parse-time diagnostics (P002, P003 warnings)
            if (parseResult.diagnostics.size() > 0) {
                json diagArray = json::array();
                for (const auto& d : parseResult.diagnostics.items()) {
                    json dj;
                    dj["severity"] = severityToString(d.severity);
                    dj["code"] = d.code;
                    dj["message"] = d.message;
                    dj["source"] = d.source;
                    if (d.line > 0) dj["line"] = d.line;
                    if (d.column > 0) dj["column"] = d.column;
                    diagArray.push_back(dj);
                }
                j["diagnostics"] = diagArray;
            }
            
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // Solve endpoint (async — returns immediately, results via SSE)
    // ================================================================
    svr.Post("/api/v1/solve", [&](const httplib::Request& req, httplib::Response& res) {
        auto sessionPtr = getSession(req, res);
        auto& session = *sessionPtr;
        if (session.solving.load()) {
            res.status = 409;
            json j = {{"error", "A solve is already in progress"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        try {
            // Parse optional body overrides
            std::string eesSource = session.eescodeContent;
            std::string initials = session.initialsContent;
            bool enableTracing = false;
            
            if (!req.body.empty()) {
                auto body = json::parse(req.body);
                if (body.contains("eescode")) eesSource = body["eescode"].get<std::string>();
                if (body.contains("initials")) initials = body["initials"].get<std::string>();
                if (body.contains("debug")) enableTracing = body["debug"].get<bool>();
            }
            
            if (eesSource.empty()) {
                res.status = 400;
                json j = {{"error", "No source code to solve"}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            // Clear previous debug output and .sol before new solve
            session.solContent.clear();
            session.latexReportContent.clear();
            session.latexReportAvailable = false;
            if (!session.debugDir.empty()) {
                std::error_code ec;
                fs::remove_all(session.debugDir, ec);
                session.debugDir.clear();
            }
            {
                std::lock_guard<std::mutex> lock(session.resultMutex);
                session.hasResult = false;
                session.lastResult = SolveResult{};
            }
            
            // Write temp files for the runner (session-isolated temp dir)
            auto tmpDir = session.tempDir;
            fs::create_directories(tmpDir);
            
            auto tmpEes = tmpDir / "model.eescode";
            writeStringToFile(tmpEes.string(), eesSource);
            
            if (!initials.empty()) {
                auto tmpInitials = tmpDir / "model.initials";
                writeStringToFile(tmpInitials.string(), initials);
            }
            
            if (!session.confContent.empty()) {
                auto tmpConf = tmpDir / "coolsolve.conf";
                writeStringToFile(tmpConf.string(), session.confContent);
            }
            
            // Clear previous progress events
            {
                std::lock_guard<std::mutex> lock(session.progressMutex);
                session.progressEvents.clear();
            }
            
            session.solving.store(true);
            session.solveFinished.store(false);
            session.cancelRequested.store(false);
            
            std::string confContent = session.confContent;
            
            // Launch solve in background thread
            std::thread([sessionPtr, tmpEes, tmpDir, enableTracing, confContent, eesSource]() {
                auto& session = *sessionPtr;
                try {
                    session.addProgressEvent("{\"type\":\"start\",\"message\":\"Solve started\"}");
                    
                    CoolSolveRunner runner(tmpEes.string());
                    
                    SolverOptions solverOpts;
                    solverOpts.verbose = false;
                    
                    // Load conf if present
                    auto confPath = tmpDir / "coolsolve.conf";
                    if (fs::exists(confPath)) {
                        loadSolverOptionsFromFile(confPath.string(), solverOpts);
                    }
                    
                    // Wire cancellation token
                    solverOpts.cancelToken = &session.cancelRequested;
                    
                    // Set up progress callback for real-time SSE events
                    solverOpts.progressCallback = [&session](
                        int blockIndex, int totalBlocks,
                        const std::string& event, int iterations, double residualNorm
                    ) {
                        json evt;
                        evt["type"] = "block";
                        evt["block"] = blockIndex;
                        evt["totalBlocks"] = totalBlocks;
                        evt["event"] = event;
                        evt["iterations"] = iterations;
                        evt["residualNorm"] = residualNorm;
                        
                        std::string msg = "Block " + std::to_string(blockIndex + 1) + "/" 
                            + std::to_string(totalBlocks);
                        if (event == "done") {
                            msg += " solved (" + std::to_string(iterations) + " iter)";
                        } else if (event == "fail") {
                            msg += " FAILED (res=" + std::to_string(residualNorm) + ")";
                        } else {
                            msg += " solving...";
                        }
                        evt["message"] = msg;
                        session.addProgressEvent(evt.dump());
                    };
                    
                    // Apply CoolProp config from solver options
                    applyCoolPropConfig(solverOpts.coolpropConfig);
                    
                    session.addProgressEvent("{\"type\":\"progress\",\"message\":\"Parsing...\"}");
                    
                    bool success = runner.run(solverOpts, enableTracing);
                    
                    // Check for parse failure and report immediately
                    if (!runner.isParseSuccess()) {
                        const auto& parseResult = runner.getParseResult();
                        std::ostringstream errMsg;
                        errMsg << "Parse failed:";
                        for (const auto& err : parseResult.errors) {
                            errMsg << "\n  Line " << err.line << ": " << err.message;
                            if (!err.context.empty()) {
                                errMsg << "\n    > " << err.context;
                            }
                        }
                        
                        // Store a minimal result
                        {
                            std::lock_guard<std::mutex> lock(session.resultMutex);
                            session.lastResult = runner.getSolveResult();
                            session.lastTiming = runner.getTiming();
                            session.lastDiagnostics = runner.getDiagnostics();
                            session.hasResult = true;
                        }
                        
                        json resultJson = solveResultToJSON(runner.getSolveResult(), runner.getTiming(),
                                                        &runner.getDiagnostics());
                        json finalEvt;
                        finalEvt["type"] = "error";
                        finalEvt["message"] = errMsg.str();
                        finalEvt["result"] = resultJson;
                        session.addProgressEvent(finalEvt.dump());
                        
                        session.solving.store(false);
                        session.solveFinished.store(true);
                        session.progressCV.notify_all();
                        return;  // Exit the solve thread early
                    }
                    
                    // Generate debug output if tracing enabled
                    if (enableTracing) {
                        auto debugPath = tmpDir / "debug_output";
                        runner.generateDebugOutput(debugPath.string(), eesSource);
                        session.debugDir = debugPath.string();
                        session.addProgressEvent("{\"type\":\"progress\",\"message\":\"Debug output generated\"}");
                    }
                    
                    // Store result
                    {
                        std::lock_guard<std::mutex> lock(session.resultMutex);
                        session.lastResult = runner.getSolveResult();
                        session.lastTiming = runner.getTiming();
                        session.lastDiagnostics = runner.getDiagnostics();
                        session.hasResult = true;
                    }
                    
                    // Extract inference data from IR for thermodynamic diagrams
                    if (runner.isIRSuccess()) {
                        const auto& ir = runner.getIR();
                        session.inferredVariables.clear();
                        std::set<std::string> fluidSet;
                        for (const auto& [name, info] : ir.getVariables()) {
                            if (Constants::isConstant(name)) continue;
                            auto it = runner.getSolveResult().variables.find(name);
                            if (it == runner.getSolveResult().variables.end()) continue;
                            InferredVariable iv;
                            iv.name = name;
                            iv.value = it->second;
                            iv.inferredProperty = info.inferredProperty;
                            iv.inferredFluid = info.inferredFluid;
                            iv.units = info.units;
                            iv.isArray = (name.find('[') != std::string::npos);
                            session.inferredVariables.push_back(iv);
                            if (!info.inferredFluid.empty()) {
                                fluidSet.insert(info.inferredFluid);
                            }
                        }
                        session.modelFluids.assign(fluidSet.begin(), fluidSet.end());
                    }
                    
                    // Build .sol content if successful
                    if (runner.isSolveSuccess()) {
                        std::ostringstream solStream;
                        solStream << std::scientific << std::setprecision(12);
                        const auto& solveResult = runner.getSolveResult();
                        const auto& ir = runner.getIR();
                        for (const auto& [name, val] : solveResult.variables) {
                            solStream << name << " = " << val;
                            const auto* info = ir.getVariable(name);
                            if (info && !info->units.empty()) {
                                solStream << " \"" << info->units << "\"";
                            }
                            solStream << "\n";
                        }
                        for (const auto& [name, val] : solveResult.stringVariables) {
                            solStream << name << " = '" << val << "'\n";
                        }
                        session.solContent = solStream.str();
                    }
                    
                    // Generate LaTeX report if enabled (lightweight string generation only)
                    if (runner.isSolveSuccess() && solverOpts.enableLatexReport) {
                        try {
                            session.latexReportContent = runner.generateLatexReportContent(session.modelName);
                            session.latexReportAvailable = !session.latexReportContent.empty();
                        } catch (const std::exception& ex) {
                            // Non-fatal: log but don't fail the solve
                            session.latexReportContent.clear();
                            session.latexReportAvailable = false;
                        }
                    }
                    
                    // Verify solution: check for CoolProp errors in final values
                    // (catches unphysical inputs that were clamped during iteration)
                    bool solutionValid = success;
                    if (runner.isSolveSuccess() && runner.isIRSuccess()) {
                        auto checkResult = checkSolution(
                            runner.getIR(),
                            runner.getSolveResult().variables,
                            runner.getSolveResult().stringVariables,
                            solverOpts.coolpropConfig);
                        if (!checkResult.allSatisfied) {
                            solutionValid = false;
                            // Merge checker diagnostics into runner diagnostics
                            // so they appear in the response
                            {
                                std::lock_guard<std::mutex> lock(session.resultMutex);
                                session.lastDiagnostics.merge(checkResult.diagnostics);
                            }
                        }
                    }
                    
                    // Build final result JSON for the done event
                    // Use session diagnostics (which may include checker diagnostics)
                    const DiagnosticCollector* diagSource = &runner.getDiagnostics();
                    DiagnosticCollector mergedDiag;
                    {
                        std::lock_guard<std::mutex> lock(session.resultMutex);
                        if (session.lastDiagnostics.size() > 0) {
                            mergedDiag = runner.getDiagnostics();
                            mergedDiag.merge(session.lastDiagnostics);
                            diagSource = &mergedDiag;
                        }
                    }
                    json resultJson = solveResultToJSON(runner.getSolveResult(), runner.getTiming(),
                                                       diagSource);
                    if (!solutionValid && runner.isSolveSuccess()) {
                        // Override success to false if solution verification failed
                        resultJson["success"] = false;
                        resultJson["status"] = "SolutionCheckFailed";
                        resultJson["errorMessage"] = "Solution verification failed: CoolProp evaluation errors in final solution";
                    }
                    if (runner.isParseSuccess()) {
                        const auto& ir = runner.getIR();
                        resultJson["equationCount"] = ir.getEquationCount();
                        resultJson["variableCount"] = ir.getNonConstantVariableCount();
                        resultJson["isSquare"] = ir.isSquare();
                    }
                    if (runner.isAnalysisSuccess()) {
                        const auto& analysis = runner.getAnalysisResult();
                        resultJson["totalBlocks"] = analysis.totalBlocks;
                        resultJson["largestBlock"] = analysis.largestBlockSize;
                    } else if (runner.isParseSuccess()) {
                        // Analysis failed (e.g. non-square system) — include analysis error
                        const auto& analysis = runner.getAnalysisResult();
                        if (!analysis.errorMessage.empty()) {
                            resultJson["errorMessage"] = analysis.errorMessage;
                        }
                    }
                    resultJson["latexReportAvailable"] = session.latexReportAvailable;
                    
                    // Send final event with embedded result
                    json finalEvt;
                    finalEvt["type"] = solutionValid ? "done" : "error";
                    finalEvt["message"] = solutionValid ? "Solve completed successfully" : "Solve failed";
                    finalEvt["result"] = resultJson;
                    session.addProgressEvent(finalEvt.dump());
                    
                } catch (const std::exception& e) {
                    json errEvt;
                    errEvt["type"] = "error";
                    errEvt["message"] = std::string("Solve error: ") + e.what();
                    session.addProgressEvent(errEvt.dump());
                }
                
                session.solving.store(false);
                session.solveFinished.store(true);
                session.progressCV.notify_all();
                
            }).detach();
            
            json j = {{"status", "started"}};
            res.set_content(j.dump(), "application/json");
            
        } catch (const std::exception& e) {
            session.solving.store(false);
            session.solveFinished.store(true);
            res.status = 500;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // SSE progress stream (real-time chunked streaming)
    // ================================================================
    svr.Get("/api/v1/solve/stream", [&](const httplib::Request& req, httplib::Response& res) {
        auto sessionPtr = getSession(req, res);
        res.set_header("Cache-Control", "no-cache");
        res.set_header("Connection", "keep-alive");
        res.set_header("X-Accel-Buffering", "no");  // Disable nginx buffering
        
        res.set_chunked_content_provider(
            "text/event-stream",
            [sessionPtr](size_t /*offset*/, httplib::DataSink& sink) {
                auto& session = *sessionPtr;
                // Wait for events (blocks up to 100ms)
                auto events = session.waitForEvents(100);
                
                for (const auto& event : events) {
                    std::string data = "data: " + event + "\n\n";
                    sink.write(data.c_str(), data.size());
                }
                
                // If solve is finished and no more events, close the stream
                if (session.solveFinished.load() && !session.solving.load()) {
                    // Drain any remaining events
                    auto remaining = session.consumeProgressEvents();
                    for (const auto& event : remaining) {
                        std::string data = "data: " + event + "\n\n";
                        sink.write(data.c_str(), data.size());
                    }
                    sink.done();
                    return false;
                }
                
                return true;
            }
        );
    });
    
    // ================================================================
    // Get last solve result
    // ================================================================
    svr.Get("/api/v1/solve/result", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        std::lock_guard<std::mutex> lock(session.resultMutex);
        if (!session.hasResult) {
            res.status = 404;
            json j = {{"error", "No solve result available"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        json j = solveResultToJSON(session.lastResult, session.lastTiming, &session.lastDiagnostics);
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Cancel solve
    // ================================================================
    svr.Post("/api/v1/solve/cancel", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        if (!session.solving.load()) {
            res.status = 409;
            json j = {{"error", "No solve is in progress"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        session.cancelRequested.store(true);
        json j = {{"success", true}, {"message", "Cancel requested"}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Parametric study endpoint (synchronous — runs all grid points)
    // ================================================================
    svr.Post("/api/v1/parametric", [&](const httplib::Request& req, httplib::Response& res) {
        auto sessionPtr = getSession(req, res);
        auto& session = *sessionPtr;
        
        if (session.solving.load()) {
            res.status = 409;
            json j = {{"error", "A solve is already in progress"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        try {
            auto body = json::parse(req.body);
            
            // Required parameters
            std::string eescode = body.value("eescode", session.eescodeContent);
            std::string initials = body.value("initials", session.initialsContent);
            std::string confContent = body.value("conf", session.confContent);
            
            // Options
            int timeoutPerPoint = body.value("timeout", 0); // seconds, 0 = no timeout
            bool updateGuesses = body.value("updateGuesses", false);

            // Swept variable(s): array of {name, values[]}
            auto sweepVars = body.at("sweepVariables");
            if (sweepVars.size() < 1 || sweepVars.size() > 2) {
                res.status = 400;
                json j = {{"error", "Provide 1 or 2 sweep variables"}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            struct SweepVar {
                std::string name;
                std::vector<double> values;
            };
            std::vector<SweepVar> sweeps;
            for (const auto& sv : sweepVars) {
                SweepVar s;
                s.name = sv.at("name").get<std::string>();
                s.values = sv.at("values").get<std::vector<double>>();
                if (s.values.empty()) {
                    res.status = 400;
                    json j = {{"error", "Sweep variable '" + s.name + "' has no values"}};
                    res.set_content(j.dump(), "application/json");
                    return;
                }
                sweeps.push_back(s);
            }
            
            // Build grid: for 1D just the values, for 2D cartesian product
            struct GridPoint {
                std::vector<std::pair<std::string, double>> overrides;
            };
            std::vector<GridPoint> grid;
            
            if (sweeps.size() == 1) {
                for (double v : sweeps[0].values) {
                    GridPoint gp;
                    gp.overrides.push_back({sweeps[0].name, v});
                    grid.push_back(gp);
                }
            } else {
                // 2D: iterate var1 as outer, var2 as inner
                for (double v1 : sweeps[0].values) {
                    for (double v2 : sweeps[1].values) {
                        GridPoint gp;
                        gp.overrides.push_back({sweeps[0].name, v1});
                        gp.overrides.push_back({sweeps[1].name, v2});
                        grid.push_back(gp);
                    }
                }
            }
            
            // For each grid point, we modify the eescode by replacing the
            // imposed equation "var = old_value" with "var = new_value".
            // This is done via regex replacement in the source text.
            
            // ---- Spiral reordering when updateGuesses is enabled ----
            // Strategy: start from the grid point closest to the original
            // imposed values, then expand outward (alternating sides for 1D,
            // BFS-by-distance for 2D).  This ensures each successive point
            // is close to an already-solved neighbour, so warm-starting from
            // the previous solution is effective.
            if (updateGuesses) {
                // Parse the original initials to find centre values
                auto parseInit = [](const std::string& text) {
                    std::map<std::string, double> m;
                    std::istringstream iss(text);
                    std::string line;
                    while (std::getline(iss, line)) {
                        auto eq = line.find('=');
                        if (eq == std::string::npos) continue;
                        std::string name = line.substr(0, eq);
                        name.erase(0, name.find_first_not_of(" \t"));
                        name.erase(name.find_last_not_of(" \t") + 1);
                        if (!name.empty() && name.back() == '$') continue;
                        std::string vs = line.substr(eq + 1);
                        auto hash = vs.find('#'); if (hash != std::string::npos) vs = vs.substr(0, hash);
                        auto sq = vs.find('\''); if (sq != std::string::npos) continue;
                        vs.erase(0, vs.find_first_not_of(" \t"));
                        vs.erase(vs.find_last_not_of(" \t") + 1);
                        try { m[name] = std::stod(vs); } catch (...) {}
                    }
                    return m;
                };
                auto initMap = parseInit(initials);

                if (sweeps.size() == 1) {
                    // 1D spiral: find index closest to initial value, then alternate L/R
                    double centre = 0.0;
                    auto it = initMap.find(sweeps[0].name);
                    if (it != initMap.end()) centre = it->second;
                    else centre = (sweeps[0].values.front() + sweeps[0].values.back()) / 2.0;

                    // grid is already sorted (linspace).  Find closest index.
                    int bestIdx = 0;
                    double bestDist = std::abs(grid[0].overrides[0].second - centre);
                    for (int k = 1; k < (int)grid.size(); ++k) {
                        double d = std::abs(grid[k].overrides[0].second - centre);
                        if (d < bestDist) { bestDist = d; bestIdx = k; }
                    }
                    std::vector<GridPoint> reordered;
                    reordered.push_back(grid[bestIdx]);
                    int lo = bestIdx - 1, hi = bestIdx + 1;
                    while (lo >= 0 || hi < (int)grid.size()) {
                        if (lo >= 0) reordered.push_back(grid[lo--]);
                        if (hi < (int)grid.size()) reordered.push_back(grid[hi++]);
                    }
                    grid = std::move(reordered);
                } else {
                    // 2D spiral: BFS ordered by Euclidean distance from centre
                    double c0 = 0.0, c1 = 0.0;
                    auto it0 = initMap.find(sweeps[0].name);
                    if (it0 != initMap.end()) c0 = it0->second;
                    else c0 = (sweeps[0].values.front() + sweeps[0].values.back()) / 2.0;
                    auto it1 = initMap.find(sweeps[1].name);
                    if (it1 != initMap.end()) c1 = it1->second;
                    else c1 = (sweeps[1].values.front() + sweeps[1].values.back()) / 2.0;

                    // Normalise axes so range spans [0,1] for distance calc
                    double range0 = sweeps[0].values.back() - sweeps[0].values.front();
                    double range1 = sweeps[1].values.back() - sweeps[1].values.front();
                    if (range0 == 0.0) range0 = 1.0;
                    if (range1 == 0.0) range1 = 1.0;

                    // Build index list sorted by distance to centre
                    std::vector<int> indices(grid.size());
                    std::iota(indices.begin(), indices.end(), 0);
                    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
                        double da0 = (grid[a].overrides[0].second - c0) / range0;
                        double da1 = (grid[a].overrides[1].second - c1) / range1;
                        double db0 = (grid[b].overrides[0].second - c0) / range0;
                        double db1 = (grid[b].overrides[1].second - c1) / range1;
                        return (da0*da0 + da1*da1) < (db0*db0 + db1*db1);
                    });
                    std::vector<GridPoint> reordered;
                    reordered.reserve(grid.size());
                    for (int idx : indices) reordered.push_back(grid[idx]);
                    grid = std::move(reordered);
                }
            }
            
            // Set up temp dir and solver options
            auto tmpDir = session.tempDir / "parametric";
            fs::create_directories(tmpDir);
            
            auto tmpEes = tmpDir / "model.eescode";
            auto tmpConf = tmpDir / "coolsolve.conf";
            
            if (!confContent.empty()) {
                writeStringToFile(tmpConf.string(), confContent);
            }
            
            // Set solving flag so UI shows progress
            session.solving.store(true);
            session.solveFinished.store(false);
            session.cancelRequested.store(false);
            {
                std::lock_guard<std::mutex> lock(session.progressMutex);
                session.progressEvents.clear();
            }
            
            // Launch parametric study in background thread
            std::thread([sessionPtr, eescode, initials, confContent, 
                         grid, sweeps, tmpDir, tmpEes, tmpConf,
                         timeoutPerPoint, updateGuesses]() {
                auto& session = *sessionPtr;
                
                session.addProgressEvent("{\"type\":\"start\",\"message\":\"Parametric study started\"}");
                
                int totalPoints = static_cast<int>(grid.size());
                json allResults = json::array();
                int successCount = 0;
                int failCount = 0;
                
                // Store the last successful solution to warm-start next point
                std::map<std::string, double, CaseInsensitiveLess> lastSolution;
                
                // Parse initials into a map for warm-starting
                auto parseInitialsStr = [](const std::string& text) {
                    std::map<std::string, double> map;
                    std::istringstream iss(text);
                    std::string line;
                    while (std::getline(iss, line)) {
                        auto eq = line.find('=');
                        if (eq == std::string::npos) continue;
                        std::string name = line.substr(0, eq);
                        // Trim
                        name.erase(0, name.find_first_not_of(" \t"));
                        name.erase(name.find_last_not_of(" \t") + 1);
                        // Skip string variables (names ending with $)
                        if (!name.empty() && name.back() == '$') continue;
                        std::string valStr = line.substr(eq + 1);
                        auto hash = valStr.find('#');
                        if (hash != std::string::npos) valStr = valStr.substr(0, hash);
                        // Remove units annotation like "[kJ/kg]" or quotes
                        auto quote = valStr.find('"');
                        if (quote != std::string::npos) valStr = valStr.substr(0, quote);
                        auto sq = valStr.find('\'');
                        if (sq != std::string::npos) continue; // skip string vars
                        valStr.erase(0, valStr.find_first_not_of(" \t"));
                        valStr.erase(valStr.find_last_not_of(" \t") + 1);
                        try {
                            map[name] = std::stod(valStr);
                        } catch (...) {}
                    }
                    return map;
                };
                
                auto baseInitials = parseInitialsStr(initials);
                
                // Helper: replace imposed value in eescode
                // Uses case-insensitive regex matching for variable name.
                // We use regex_search + manual string surgery instead of
                // regex_replace to avoid backreference ambiguity (e.g. "$31"
                // would be interpreted as group 31 instead of group 3 + "1").
                auto replaceImposedValue = [](const std::string& source, 
                                              const std::string& varName, 
                                              double newValue) -> std::string {
                    std::string escapedName;
                    for (char c : varName) {
                        if (c == '[' || c == ']' || c == '(' || c == ')' || 
                            c == '.' || c == '+' || c == '*' || c == '?') {
                            escapedName += '\\';
                        }
                        escapedName += c;
                    }
                    std::string pattern = "(^|\\n)(\\s*)" + escapedName + 
                                          "(\\s*=\\s*)-?[0-9]+\\.?[0-9]*([eE][+-]?[0-9]+)?";
                    std::regex re(pattern, std::regex_constants::icase);
                    
                    std::ostringstream valOss;
                    valOss << std::setprecision(15) << newValue;
                    std::string valStr = valOss.str();
                    
                    // Manual match-and-replace to avoid $N ambiguity
                    std::string result;
                    auto it = source.cbegin();
                    std::smatch m;
                    while (std::regex_search(it, source.cend(), m, re)) {
                        // Append text before match
                        result.append(it, it + m.position(0));
                        // Rebuild the matched region with the new value
                        result += m[1].str() + m[2].str() + varName + m[3].str() + valStr;
                        it += m.position(0) + m[0].length();
                    }
                    // Append remainder
                    result.append(it, source.cend());
                    return result;
                };
                
                // ============================================================
                // Performance optimisation: factor out invariant pipeline
                // stages so they run *once* rather than per-point.
                //
                //  1. Parse the original eescode → AST
                //  2. Build IR from AST
                //  3. inferVariables + initializeVariables  (~10 ms saved/pt)
                //  4. Structural analysis                   (~0.5 ms saved/pt)
                //
                // Per-point we only re-parse (constants change), re-build
                // the IR, load initials (which provides all solutionValues),
                // and solve using the *cached* analysis result.
                // ============================================================
                
                // Prepare solver options once (invariant across points)
                SolverOptions solverOpts;
                solverOpts.verbose = false;
                if (!confContent.empty() && fs::exists(tmpConf)) {
                    loadSolverOptionsFromFile(tmpConf.string(), solverOpts);
                }
                solverOpts.cancelToken = &session.cancelRequested;
                if (timeoutPerPoint > 0) {
                    solverOpts.timeoutSeconds = timeoutPerPoint;
                }
                applyCoolPropConfig(solverOpts.coolpropConfig);
                
                // One-time CoolProp warmup (first call pays library init cost)
                warmupCoolProp();
                
                // Run the full pipeline once on the original eescode to
                // obtain the structural analysis result.
                StructuralAnalysisResult cachedAnalysis;
                {
                    writeStringToFile(tmpEes.string(), eescode);
                    EESParser parser;
                    auto parseResult = parser.parseFile(tmpEes.string());
                    if (parseResult.success) {
                        try {
                            auto templateIR = IR::fromAST(parseResult.program);
                            inferVariables(templateIR);
                            initializeVariables(templateIR);
                            cachedAnalysis = StructuralAnalyzer::analyze(templateIR);
                        } catch (...) {
                            // Fall through — cachedAnalysis.success stays false
                        }
                    }
                    if (!cachedAnalysis.success) {
                        std::cerr << "[Parametric] WARNING: pre-analysis failed, "
                                     "falling back to full pipeline per point\n";
                    }
                }
                
                for (int i = 0; i < totalPoints; ++i) {
                    if (session.cancelRequested.load()) {
                        session.addProgressEvent("{\"type\":\"progress\",\"message\":\"Parametric study cancelled\"}");
                        break;
                    }
                    
                    const auto& gp = grid[i];
                    
                    // Progress event
                    {
                        json evt;
                        evt["type"] = "progress";
                        std::ostringstream msg;
                        msg << "Point " << (i+1) << "/" << totalPoints << ": ";
                        for (size_t j = 0; j < gp.overrides.size(); ++j) {
                            if (j > 0) msg << ", ";
                            msg << gp.overrides[j].first << " = " << gp.overrides[j].second;
                        }
                        evt["message"] = msg.str();
                        session.addProgressEvent(evt.dump());
                    }
                    
                    // Modify the eescode: replace each imposed value
                    std::string modifiedEescode = eescode;
                    for (const auto& [varName, newVal] : gp.overrides) {
                        modifiedEescode = replaceImposedValue(modifiedEescode, varName, newVal);
                    }
                    
                    // Write modified eescode
                    writeStringToFile(tmpEes.string(), modifiedEescode);
                    
                    // Build initial guesses
                    std::map<std::string, double, CaseInsensitiveLess> currentInitials(baseInitials.begin(), baseInitials.end());
                    if (updateGuesses && !lastSolution.empty()) {
                        // Warm-start from last solution (only when updateGuesses is on)
                        for (const auto& [name, val] : lastSolution) {
                            currentInitials[name] = val;
                        }
                    }
                    // Override sweep variables in initials too
                    for (const auto& [varName, newVal] : gp.overrides) {
                        currentInitials[varName] = newVal;
                    }
                    
                    auto tPointStart = std::chrono::high_resolution_clock::now();
                    bool success = false;
                    SolveResult pointSolveResult;
                    double parseMs = 0, irMs = 0, inferMs = 0, analysisMs = 0, solveMs = 0;
                    
                    if (cachedAnalysis.success) {
                        // ---- Fast path: reuse cached structural analysis ----
                        try {
                            auto t1 = std::chrono::high_resolution_clock::now();
                            EESParser parser;
                            auto parseResult = parser.parseFile(tmpEes.string());
                            auto t2 = std::chrono::high_resolution_clock::now();
                            parseMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
                            
                            if (!parseResult.success) throw std::runtime_error("parse failed");
                            
                            t1 = std::chrono::high_resolution_clock::now();
                            IR ir = IR::fromAST(parseResult.program);
                            t2 = std::chrono::high_resolution_clock::now();
                            irMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
                            
                            // Skip inferVariables + initializeVariables (the big win).
                            // The initials provide solutionValue for all variables,
                            // which the Solver uses in preference to guessValue.
                            inferMs = 0;
                            analysisMs = 0;
                            
                            // Load initials directly into the IR
                            for (const auto& [name, val] : currentInitials) {
                                if (!name.empty() && name.back() == '$') continue;
                                auto* vinfo = ir.getVariableMutable(name);
                                if (vinfo) {
                                    vinfo->solutionValue = val;
                                }
                            }
                            
                            t1 = std::chrono::high_resolution_clock::now();
                            Solver solver(ir, cachedAnalysis, solverOpts.coolpropConfig);
                            pointSolveResult = solver.solve(solverOpts, false);
                            t2 = std::chrono::high_resolution_clock::now();
                            solveMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
                            
                            success = pointSolveResult.success;
                        } catch (...) {
                            success = false;
                            pointSolveResult.errorMessage = "Exception in fast-path solve";
                        }
                    } else {
                        // ---- Fallback: full pipeline via CoolSolveRunner ----
                        // Build initials string and write file
                        std::ostringstream initOss;
                        initOss << std::scientific << std::setprecision(12);
                        for (const auto& [name, val] : currentInitials) {
                            if (!name.empty() && name.back() == '$') continue;
                            initOss << name << " = " << val << "\n";
                        }
                        auto tmpInit = tmpDir / "model.initials";
                        writeStringToFile(tmpInit.string(), initOss.str());
                        
                        CoolSolveRunner runner(tmpEes.string());
                        success = runner.run(solverOpts, false);
                        pointSolveResult = runner.getSolveResult();
                        
                        auto& t = runner.getTiming();
                        parseMs = t.parse_time_ms;
                        irMs = t.ir_time_ms;
                        inferMs = t.infer_time_ms;
                        analysisMs = t.analysis_time_ms;
                        solveMs = t.solve_time_ms;
                    }
                    
                    auto tPointEnd = std::chrono::high_resolution_clock::now();
                    
                    // Profiling: log per-phase timing to stderr
                    {
                        double totalMs = std::chrono::duration<double, std::milli>(tPointEnd - tPointStart).count();
                        std::cerr << "[Parametric] Point " << (i+1) << "/" << totalPoints
                                  << " total=" << std::fixed << std::setprecision(1) << totalMs << "ms"
                                  << " parse=" << parseMs
                                  << " ir=" << irMs
                                  << " infer=" << inferMs
                                  << " analysis=" << analysisMs
                                  << " solve=" << solveMs
                                  << (success ? " OK" : " FAIL") << "\n";
                    }
                    
                    json pointResult;
                    pointResult["index"] = i;
                    pointResult["success"] = success;
                    
                    // Store override values
                    json overrides = json::object();
                    for (const auto& [varName, newVal] : gp.overrides) {
                        overrides[varName] = newVal;
                    }
                    pointResult["overrides"] = overrides;
                    
                    if (success) {
                        successCount++;
                        json vars = json::object();
                        for (const auto& [name, val] : pointSolveResult.variables) {
                            vars[name] = val;
                        }
                        pointResult["variables"] = vars;
                        lastSolution = pointSolveResult.variables;
                    } else {
                        failCount++;
                        pointResult["errorMessage"] = pointSolveResult.errorMessage;
                        // Do NOT update lastSolution from failed runs — keep the
                        // last successful solution so the next point has good guesses.
                    }
                    
                    allResults.push_back(pointResult);
                }
                
                // Build final response
                json response;
                response["success"] = (failCount == 0);
                response["totalPoints"] = totalPoints;
                response["successCount"] = successCount;
                response["failCount"] = failCount;
                response["results"] = allResults;
                
                // Sweep variable metadata
                json sweepMeta = json::array();
                for (const auto& s : sweeps) {
                    json sm;
                    sm["name"] = s.name;
                    sm["values"] = s.values;
                    sweepMeta.push_back(sm);
                }
                response["sweepVariables"] = sweepMeta;
                
                // Send final event
                json finalEvt;
                finalEvt["type"] = failCount == 0 ? "done" : "error";
                finalEvt["message"] = "Parametric study completed: " + 
                    std::to_string(successCount) + "/" + std::to_string(totalPoints) + " points solved";
                finalEvt["result"] = response;
                session.addProgressEvent(finalEvt.dump());
                
                // Store parametric result in session
                {
                    std::lock_guard<std::mutex> lock(session.parametricMutex);
                    session.lastParametricResult = response.dump();
                }
                
                session.solving.store(false);
                session.solveFinished.store(true);
                session.progressCV.notify_all();
                
            }).detach();
            
            json j = {{"status", "started"}, {"totalPoints", static_cast<int>(grid.size())}};
            res.set_content(j.dump(), "application/json");
            
        } catch (const std::exception& e) {
            session.solving.store(false);
            session.solveFinished.store(true);
            res.status = 400;
            json j = {{"error", e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // Get last parametric study result
    // ================================================================
    svr.Get("/api/v1/parametric/result", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        std::lock_guard<std::mutex> lock(session.parametricMutex);
        if (session.lastParametricResult.empty()) {
            res.status = 404;
            json j = {{"error", "No parametric study result available"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        res.set_content(session.lastParametricResult, "application/json");
    });
    
    // ================================================================
    // Update guesses: copy .sol -> .initials
    // ================================================================
    svr.Post("/api/v1/update-guesses", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        if (session.solContent.empty()) {
            res.status = 400;
            json j = {{"error", "No solution available to use as guesses"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        session.initialsContent = session.solContent;
        json j = {{"success", true}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Variables endpoint (from last solve)
    // ================================================================
    svr.Get("/api/v1/variables", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        std::lock_guard<std::mutex> lock(session.resultMutex);
        if (!session.hasResult) {
            res.status = 404;
            json j = {{"error", "No solve result available"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        json vars = json::array();
        for (const auto& [name, val] : session.lastResult.variables) {
            json varObj;
            varObj["name"] = name;
            varObj["value"] = val;
            varObj["isArray"] = (name.find('[') != std::string::npos);
            vars.push_back(varObj);
        }
        for (const auto& [name, val] : session.lastResult.stringVariables) {
            json varObj;
            varObj["name"] = name;
            varObj["value"] = val;
            varObj["isString"] = true;
            vars.push_back(varObj);
        }
        
        json j = {{"variables", vars}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Inferred variables endpoint (for thermodynamic diagrams)
    // ================================================================
    svr.Get("/api/v1/variables/inferred", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        
        json vars = json::array();
        for (const auto& iv : session.inferredVariables) {
            json varObj;
            varObj["name"] = iv.name;
            varObj["value"] = iv.value;
            varObj["inferredProperty"] = iv.inferredProperty;
            varObj["inferredFluid"] = iv.inferredFluid;
            varObj["units"] = iv.units;
            varObj["isArray"] = iv.isArray;
            vars.push_back(varObj);
        }
        
        json j = {{"variables", vars}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // CoolProp fluids listing
    // ================================================================
    svr.Get("/api/v1/coolprop/fluids", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        
        json fluidsArr = json::array();
        auto allFluids = FluidRegistry::getAllFluids();
        // Sort by name for consistent output
        std::sort(allFluids.begin(), allFluids.end(),
            [](const auto& a, const auto& b) { return a->getName() < b->getName(); });
        
        for (const auto& fluid : allFluids) {
            std::string typeStr;
            switch (fluid->getType()) {
                case FluidType::Real: typeStr = "Real"; break;
                case FluidType::IdealGas: typeStr = "IdealGas"; break;
                case FluidType::HumidAir: typeStr = "HumidAir"; break;
                case FluidType::Incompressible: typeStr = "Incompressible"; break;
                case FluidType::Mixture: typeStr = "Mixture"; break;
                default: typeStr = "Unknown"; break;
            }
            fluidsArr.push_back({
                {"name", fluid->getName()},
                {"coolpropName", fluid->getCoolPropName()},
                {"type", typeStr},
                {"hasDome", fluid->getType() == FluidType::Real}
            });
        }
        
        json j;
        j["fluids"] = fluidsArr;
        j["modelFluids"] = session.modelFluids;
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // CoolProp saturation dome endpoint
    // ================================================================
    svr.Post("/api/v1/coolprop/saturation", [&](const httplib::Request& req, httplib::Response& res) {
        try {
            auto body = json::parse(req.body);
            if (!body.contains("fluid")) {
                res.status = 400;
                json j = {{"error", "Missing 'fluid' parameter"}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            std::string fluidName = body["fluid"].get<std::string>();
            int nPoints = body.value("nPoints", 200);
            
            // Validate fluid
            auto fluid = FluidRegistry::getFluid(fluidName);
            if (!fluid) {
                res.status = 400;
                json j = {{"error", "Unknown fluid: " + fluidName}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            if (fluid->getType() != FluidType::Real) {
                res.status = 400;
                std::string typeStr;
                switch (fluid->getType()) {
                    case FluidType::IdealGas: typeStr = "ideal gas"; break;
                    case FluidType::HumidAir: typeStr = "humid air"; break;
                    case FluidType::Incompressible: typeStr = "incompressible fluid"; break;
                    case FluidType::Mixture: typeStr = "mixture"; break;
                    default: typeStr = "unsupported fluid type"; break;
                }
                json j = {{"error", "No saturation curve for " + typeStr + ": " + fluidName}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            std::string cpFluid = fluid->getCoolPropName();
            
            auto t_start = std::chrono::high_resolution_clock::now();
            
            // Get critical and triple point temperatures
            double Tcrit, Pcrit, Ttriple;
            try {
                Tcrit = CoolProp::PropsSI("Tcrit", "", 0, "", 0, cpFluid);
                Pcrit = CoolProp::PropsSI("pcrit", "", 0, "", 0, cpFluid);
                Ttriple = CoolProp::PropsSI("T_triple", "", 0, "", 0, cpFluid);
            } catch (const std::exception& e) {
                res.status = 500;
                json j = {{"error", std::string("CoolProp error getting critical/triple point: ") + e.what()}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            // Get critical point properties
            double Hcrit = 0, Scrit = 0, Dcrit = 0;
            try {
                Hcrit = CoolProp::PropsSI("H", "T", Tcrit - 0.01, "Q", 0, cpFluid);
                Scrit = CoolProp::PropsSI("S", "T", Tcrit - 0.01, "Q", 0, cpFluid);
                Dcrit = CoolProp::PropsSI("D", "T", Tcrit - 0.01, "Q", 0, cpFluid);
            } catch (...) {
                // Non-fatal: critical point properties may fail for some fluids
            }
            
            // Generate temperature points from Ttriple+1 to Tcrit-0.1
            double Tmin = Ttriple + 1.0;
            double Tmax = Tcrit - 0.1;
            if (Tmin >= Tmax) {
                res.status = 400;
                json j = {{"error", "Temperature range too narrow for saturation curve"}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            // Use logarithmic spacing in (Tcrit - T) for better resolution near critical point
            std::vector<double> temps(nPoints);
            double logMin = std::log(Tcrit - Tmax);  // small
            double logMax = std::log(Tcrit - Tmin);   // large
            for (int i = 0; i < nPoints; i++) {
                double frac = static_cast<double>(i) / (nPoints - 1);
                double logVal = logMax - frac * (logMax - logMin);  // from far to near critical
                temps[i] = Tcrit - std::exp(logVal);
            }
            
            // Compute saturation properties using PropsSImulti for efficiency
            // CoolProp reuses the AbstractState object internally when given vectors
            json liquid, vapor;
            std::vector<double> lT, lP, lH, lS, lD;
            std::vector<double> vT, vP, vH, vS, vD;
            
            for (int i = 0; i < nPoints; i++) {
                double T = temps[i];
                try {
                    double pL = CoolProp::PropsSI("P", "T", T, "Q", 0, cpFluid);
                    double hL = CoolProp::PropsSI("H", "T", T, "Q", 0, cpFluid);
                    double sL = CoolProp::PropsSI("S", "T", T, "Q", 0, cpFluid);
                    double dL = CoolProp::PropsSI("D", "T", T, "Q", 0, cpFluid);
                    
                    double pV = CoolProp::PropsSI("P", "T", T, "Q", 1, cpFluid);
                    double hV = CoolProp::PropsSI("H", "T", T, "Q", 1, cpFluid);
                    double sV = CoolProp::PropsSI("S", "T", T, "Q", 1, cpFluid);
                    double dV = CoolProp::PropsSI("D", "T", T, "Q", 1, cpFluid);
                    
                    // Skip non-finite values
                    if (!std::isfinite(pL) || !std::isfinite(hL) || !std::isfinite(sL) || !std::isfinite(dL) ||
                        !std::isfinite(pV) || !std::isfinite(hV) || !std::isfinite(sV) || !std::isfinite(dV)) {
                        continue;
                    }
                    
                    lT.push_back(T); lP.push_back(pL); lH.push_back(hL); lS.push_back(sL); lD.push_back(dL);
                    vT.push_back(T); vP.push_back(pV); vH.push_back(hV); vS.push_back(sV); vD.push_back(dV);
                } catch (...) {
                    // Skip points where CoolProp fails (near triple/critical point)
                    continue;
                }
            }
            
            auto t_end = std::chrono::high_resolution_clock::now();
            double computeTime = std::chrono::duration<double, std::milli>(t_end - t_start).count();
            
            // Convert temperatures from CoolProp SI (K) to CoolSolve units (°C)
            UnitSystem defaultUnits;
            for (auto& t : lT) t = UnitConverter::fromSI(t, UnitType::Temperature, defaultUnits.temperature);
            for (auto& t : vT) t = UnitConverter::fromSI(t, UnitType::Temperature, defaultUnits.temperature);
            double TcritCS = UnitConverter::fromSI(Tcrit, UnitType::Temperature, defaultUnits.temperature);
            double TtripleCS = UnitConverter::fromSI(Ttriple, UnitType::Temperature, defaultUnits.temperature);
            
            liquid["T"] = lT; liquid["P"] = lP; liquid["H"] = lH; liquid["S"] = lS; liquid["D"] = lD;
            vapor["T"] = vT; vapor["P"] = vP; vapor["H"] = vH; vapor["S"] = vS; vapor["D"] = vD;
            
            json j;
            j["fluid"] = fluidName;
            j["critical"] = {{"T", TcritCS}, {"P", Pcrit}, {"H", Hcrit}, {"S", Scrit}, {"D", Dcrit}};
            j["triplePoint"] = {{"T", TtripleCS}};
            j["liquid"] = liquid;
            j["vapor"] = vapor;
            j["nPoints"] = static_cast<int>(lT.size());
            j["computeTime_ms"] = computeTime;
            
            res.set_content(j.dump(), "application/json");
            
        } catch (const json::exception& e) {
            res.status = 400;
            json j = {{"error", std::string("Invalid JSON: ") + e.what()}};
            res.set_content(j.dump(), "application/json");
        } catch (const std::exception& e) {
            res.status = 500;
            json j = {{"error", std::string("Server error: ") + e.what()}};
            res.set_content(j.dump(), "application/json");
        }
    });
    
    // ================================================================
    // Examples listing
    // ================================================================
    svr.Get("/api/v1/examples", [&](const httplib::Request&, httplib::Response& res) {
        json examples = json::array();
        
        // Look for examples in the examples/ directory relative to the binary
        // or from the source tree
        std::vector<std::string> searchPaths = {
            "examples",
            "../examples",
            fs::current_path().string() + "/examples"
        };
        
        for (const auto& searchPath : searchPaths) {
            if (fs::exists(searchPath) && fs::is_directory(searchPath)) {
                for (const auto& entry : fs::directory_iterator(searchPath)) {
                    if (entry.path().extension() == ".eescode") {
                        json exObj;
                        exObj["name"] = entry.path().stem().string();
                        exObj["path"] = fs::canonical(entry.path()).string();
                        examples.push_back(exObj);
                    }
                }
                break;  // Use first found directory
            }
        }
        
        json j = {{"examples", examples}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Bundle download (ZIP of all model files + debug if present)
    // ================================================================
    svr.Get("/api/v1/files/bundle", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        std::string stem = session.modelName.empty() ? "model" : session.modelName;
        
        std::vector<std::pair<std::string, std::string>> files;
        if (!session.eescodeContent.empty())
            files.push_back({stem + ".eescode", session.eescodeContent});
        if (!session.initialsContent.empty())
            files.push_back({stem + ".initials", session.initialsContent});
        if (!session.solContent.empty())
            files.push_back({stem + ".sol", session.solContent});
        if (!session.confContent.empty())
            files.push_back({"coolsolve.conf", session.confContent});
        
        // Include LaTeX report if generated
        if (session.latexReportAvailable && !session.latexReportContent.empty())
            files.push_back({stem + "_report.tex", session.latexReportContent});
        
        // Include debug files if they exist
        if (!session.debugDir.empty() && fs::exists(session.debugDir)) {
            for (const auto& entry : fs::directory_iterator(session.debugDir)) {
                if (entry.is_regular_file()) {
                    std::string content = readFileToString(entry.path().string());
                    files.push_back({"debug_output/" + entry.path().filename().string(), content});
                }
            }
        }
        
        // Include parametric study results if present
        {
            std::lock_guard<std::mutex> lock(session.parametricMutex);
            if (!session.lastParametricResult.empty()) {
                files.push_back({"parametric_studies.json", session.lastParametricResult});
            }
        }
        
        if (files.empty()) {
            res.status = 400;
            json j = {{"error", "No files to bundle"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        std::string zip = createZipBundle(files);
        res.set_header("Content-Disposition", "attachment; filename=\"" + stem + ".zip\"");
        res.set_content(zip, "application/zip");
    });
    
    // ================================================================
    // File upload (ZIP only)
    // ================================================================
    svr.Post("/api/v1/files/upload", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        
        // Find the uploaded ZIP file
        std::string zipContent;
        std::string zipFilename;
        for (const auto& [name, file] : req.files) {
            if (endsWith(file.filename, ".zip") && file.content.size() > 4) {
                zipContent = file.content;
                zipFilename = file.filename;
                break;
            }
        }
        
        if (zipContent.empty()) {
            res.status = 400;
            json j = {{"error", "Please upload a .zip file"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        auto extracted = extractZipBundle(zipContent);
        
        // Verify there is an .eescode file
        bool hasEescode = false;
        for (const auto& [zname, _] : extracted) {
            if (endsWith(zname, ".eescode") && !startsWith(zname, "debug_output/")) {
                hasEescode = true; break;
            }
        }
        if (!hasEescode) {
            res.status = 400;
            json j = {{"error", "ZIP must contain an .eescode file"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        // Save snapshot for Back button
        session.saveSnapshot();
        
        // Clear previous state
        session.eescodeContent.clear();
        session.initialsContent.clear();
        session.solContent.clear();
        session.confContent.clear();
        session.openFilePath.clear();
        session.modelName.clear();
        // Set model name from ZIP filename
        if (!zipFilename.empty()) {
            session.modelName = fs::path(zipFilename).stem().string();
        }
        if (!session.debugDir.empty()) {
            std::error_code ec;
            fs::remove_all(session.debugDir, ec);
            session.debugDir.clear();
        }
        {
            std::lock_guard<std::mutex> lock(session.resultMutex);
            session.hasResult = false;
            session.lastResult = SolveResult{};
        }
        
        // Extract files from ZIP
        json fileList = json::array();
        bool hasDebug = false;
        for (const auto& [zname, zcontent] : extracted) {
            if (startsWith(zname, "debug_output/") && zname.size() > 13) {
                // Write debug files to session temp dir
                auto debugDir = session.tempDir / "debug_output";
                fs::create_directories(debugDir);
                std::string debugFileName = zname.substr(13);
                writeStringToFile((debugDir / debugFileName).string(), zcontent);
                session.debugDir = debugDir.string();
                hasDebug = true;
            } else if (endsWith(zname, ".eescode")) {
                session.eescodeContent = zcontent;
                fileList.push_back(zname);
            } else if (endsWith(zname, ".initials")) {
                session.initialsContent = zcontent;
                fileList.push_back(zname);
            } else if (endsWith(zname, ".sol")) {
                session.solContent = zcontent;
                fileList.push_back(zname);
            } else if (zname == "coolsolve.conf" || endsWith(zname, ".conf")) {
                session.confContent = zcontent;
                fileList.push_back(zname);
            } else if (zname == "parametric_studies.json") {
                std::lock_guard<std::mutex> lock(session.parametricMutex);
                session.lastParametricResult = zcontent;
                fileList.push_back(zname);
            } else if (endsWith(zname, "_report.tex")) {
                session.latexReportContent = zcontent;
                session.latexReportAvailable = true;
                fileList.push_back(zname);
            }
        }
        if (hasDebug) fileList.push_back("debug_output/");
        
        json j = {{"success", true}, {"files", fileList}, {"modelName", session.modelName}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Debug output files
    // ================================================================
    svr.Get("/api/v1/debug/files", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        if (session.debugDir.empty() || !fs::exists(session.debugDir)) {
            json j = {{"files", json::array()}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        json files = json::array();
        for (const auto& entry : fs::directory_iterator(session.debugDir)) {
            if (entry.is_regular_file()) {
                json f;
                f["name"] = entry.path().filename().string();
                f["size"] = entry.file_size();
                files.push_back(f);
            }
        }
        json j = {{"files", files}};
        res.set_content(j.dump(), "application/json");
    });
    
    svr.Get("/api/v1/debug/file", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        if (session.debugDir.empty()) {
            res.status = 404;
            json j = {{"error", "No debug output available. Run a Debug Solve first."}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        auto name = req.get_param_value("name");
        if (name.empty() || name.find('/') != std::string::npos ||
            name.find('\\') != std::string::npos || name.find("..") != std::string::npos) {
            res.status = 400;
            json j = {{"error", "Invalid or missing filename"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        fs::path filePath = fs::path(session.debugDir) / name;
        if (!fs::exists(filePath)) {
            res.status = 404;
            json j = {{"error", "File not found: " + name}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        std::string content = readFileToString(filePath.string());
        json j = {{"name", name}, {"content", content}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // LaTeX report
    // ================================================================
    
    // GET /api/v1/latex/report – return the generated LaTeX source
    svr.Get("/api/v1/latex/report", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        json j;
        j["available"] = session.latexReportAvailable;
        j["content"] = session.latexReportAvailable ? session.latexReportContent : "";
        res.set_content(j.dump(), "application/json");
    });
    
    // POST /api/v1/latex/compile – compile the report to PDF and return it.
    // Body JSON: { compiler?: string, plots?: [{name: string, data: string}] }
    //   compiler  – LaTeX command (default: pdflatex, or from coolsolve.conf)
    //   plots     – array of {name, data} where data is base64-encoded PNG
    //              (the data URI prefix "data:image/png;base64," is stripped automatically)
    svr.Post("/api/v1/latex/compile", [&](const httplib::Request& req, httplib::Response& res) {
        auto& session = *getSession(req, res);
        std::cerr << "[latex/compile] Request received, latexReportAvailable=" 
                  << session.latexReportAvailable 
                  << ", contentLen=" << session.latexReportContent.size() << std::endl;
        
        if (!session.latexReportAvailable || session.latexReportContent.empty()) {
            res.status = 400;
            json j = {{"error", "No LaTeX report available. Solve with enableLatexReport = true first."}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        // Parse request body
        std::string compiler = "pdflatex";
        std::vector<std::pair<std::string, std::string>> plots; // name → decoded PNG bytes
        
        if (!req.body.empty()) {
            try {
                auto body = json::parse(req.body);
                if (body.contains("compiler") && body["compiler"].is_string()) {
                    std::string c = body["compiler"].get<std::string>();
                    if (!c.empty()) compiler = c;
                }
                if (body.contains("plots") && body["plots"].is_array()) {
                    for (const auto& p : body["plots"]) {
                        if (!p.contains("name") || !p.contains("data")) continue;
                        std::string name = p["name"].get<std::string>();
                        std::string data = p["data"].get<std::string>();
                        // Strip data URI prefix if present
                        const std::string prefix = "data:image/png;base64,";
                        if (data.size() > prefix.size() &&
                            data.compare(0, prefix.size(), prefix) == 0) {
                            data = data.substr(prefix.size());
                        }
                        std::string decoded = base64Decode(data);
                        if (!decoded.empty()) {
                            plots.push_back({name, std::move(decoded)});
                        }
                    }
                }
            } catch (...) {
                // Proceed with defaults if parsing fails
            }
        }
        
        // Also try to read latexCompiler from the session's conf
        if (!session.confContent.empty() && compiler == "pdflatex") {
            // Quick grep for latexCompiler in the conf text
            std::istringstream confStream(session.confContent);
            std::string line;
            while (std::getline(confStream, line)) {
                auto trimmed = line;
                while (!trimmed.empty() && (trimmed[0] == ' ' || trimmed[0] == '\t'))
                    trimmed.erase(trimmed.begin());
                if (trimmed.empty() || trimmed[0] == '#') continue;
                auto eq = trimmed.find('=');
                if (eq == std::string::npos) continue;
                std::string key = trimmed.substr(0, eq);
                std::string val = trimmed.substr(eq + 1);
                // Trim
                while (!key.empty() && key.back() == ' ') key.pop_back();
                while (!val.empty() && val.front() == ' ') val.erase(val.begin());
                if (key == "latexCompiler" && !val.empty()) {
                    compiler = val;
                }
            }
        }
        
        // Create temporary build directory
        auto buildDir = session.tempDir / "latex_build";
        std::error_code ec;
        fs::create_directories(buildDir, ec);
        
        // Write .tex source
        std::string texName = "report.tex";
        {
            std::ofstream f(buildDir / texName, std::ios::binary);
            if (!f.is_open()) {
                res.status = 500;
                json j = {{"error", "Failed to write .tex file"}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            f << session.latexReportContent;
        }
        
        // Write plot images
        for (const auto& [name, data] : plots) {
            // Sanitise filename: only allow alphanumeric, underscore, dot
            std::string safeName;
            for (char c : name) {
                if (std::isalnum(c) || c == '_' || c == '.') safeName += c;
            }
            if (safeName.empty()) continue;
            std::ofstream f(buildDir / safeName, std::ios::binary);
            if (f.is_open()) {
                f.write(data.data(), static_cast<std::streamsize>(data.size()));
            }
        }
        
        // Run LaTeX compiler (two passes for TOC / cross-references)
        std::string cmd = compiler + " -interaction=nonstopmode -halt-on-error"
            " -output-directory=\"" + buildDir.string() + "\""
            " \"" + (buildDir / texName).string() + "\""
            " 2>&1";
        
        std::string compilerOutput;
        for (int pass = 0; pass < 2; ++pass) {
            FILE* pipe = popen(cmd.c_str(), "r");
            if (!pipe) {
                res.status = 500;
                json j = {{"error", "Failed to run " + compiler}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            compilerOutput.clear();
            char buf[4096];
            while (fgets(buf, sizeof(buf), pipe)) {
                compilerOutput += buf;
            }
            int status = pclose(pipe);
            // Only fail on the second pass (first pass may have unresolved refs)
            if (pass == 1 && status != 0) {
                // Try to extract a meaningful error from the log
                std::string errMsg = "LaTeX compilation failed";
                auto excl = compilerOutput.find('!');
                if (excl != std::string::npos) {
                    auto eol = compilerOutput.find('\n', excl);
                    errMsg = compilerOutput.substr(excl, std::min(eol - excl, size_t(200)));
                }
                res.status = 422;
                json j = {{"error", errMsg}, {"log", compilerOutput.substr(0, 4000)}};
                res.set_content(j.dump(), "application/json");
                return;
            }
        }
        
        // Read the resulting PDF
        fs::path pdfPath = buildDir / "report.pdf";
        if (!fs::exists(pdfPath)) {
            res.status = 500;
            json j = {{"error", "PDF was not generated"}, {"log", compilerOutput.substr(0, 4000)}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        std::ifstream pdfFile(pdfPath, std::ios::binary);
        std::ostringstream pdfStream;
        pdfStream << pdfFile.rdbuf();
        std::string pdfContent = pdfStream.str();
        
        // Clean up build directory (keep it small)
        // We leave it — it's in the session temp dir, cleaned on session destruction.
        
        std::string downloadName = (session.modelName.empty() ? "model" : session.modelName) + "_report.pdf";
        std::cerr << "[latex/compile] PDF generated: " << pdfContent.size() << " bytes → " << downloadName << std::endl;
        res.set_header("Content-Disposition", "attachment; filename=\"" + downloadName + "\"");
        res.set_content(pdfContent, "application/pdf");
    });
    
    // ================================================================
    // Static file serving (SPA)
    // ================================================================
    
#ifdef COOLSOLVE_EMBEDDED_ASSETS
    // Embedded mode: serve assets from compiled-in data
    // Serve known static assets from embedded data
    svr.Get("/.*", [](const httplib::Request& req, httplib::Response& res) {
        std::string path = req.path;
        if (startsWith(path, "/api/")) return; // Let API handlers deal with it
        
        // Try exact path match
        const EmbeddedAsset* asset = getEmbeddedAsset(path);
        
        // Serve asset if found
        if (asset) {
            res.set_content(reinterpret_cast<const char*>(asset->data), asset->size, asset->mimeType);
            // Cache hashed assets for 1 year
            if (path.find("/assets/") != std::string::npos) {
                res.set_header("Cache-Control", "public, max-age=31536000, immutable");
            }
            return;
        }
        
        // SPA fallback: serve index.html for unknown routes
        const EmbeddedAsset* index = getEmbeddedAsset("/index.html");
        if (index) {
            res.set_content(reinterpret_cast<const char*>(index->data), index->size, "text/html");
            return;
        }
    });
    std::cout << "Serving GUI from embedded assets\n";
#else
    // Filesystem mode: serve from dist directory or dev server
    std::string staticDir;
    
#ifdef COOLSOLVE_GUI_DEV_DIR
    // Dev mode: serve from filesystem
    staticDir = COOLSOLVE_GUI_DEV_DIR;
    if (!options.guiDevDir.empty()) {
        staticDir = options.guiDevDir;
    }
#endif

    if (!staticDir.empty() && fs::exists(staticDir)) {
        svr.set_mount_point("/", staticDir);
        std::cout << "Serving GUI from: " << staticDir << "\n";
    }
    
    // SPA fallback: serve index.html for any non-API, non-file route
    svr.set_error_handler([&](const httplib::Request& req, httplib::Response& res) {
        if (res.status == 404 && !startsWith(req.path, "/api/")) {
            // Try to serve index.html for SPA routing
            if (!staticDir.empty()) {
                auto indexPath = fs::path(staticDir) / "index.html";
                if (fs::exists(indexPath)) {
                    std::string content = readFileToString(indexPath.string());
                    res.set_content(content, "text/html");
                    res.status = 200;
                    return;
                }
            }
            
            // If no static files, serve a minimal redirect page
            res.set_content(R"html(<!DOCTYPE html>
<html>
<head><title>CoolSolve GUI</title></head>
<body>
<h1>CoolSolve GUI</h1>
<p>The frontend is not built yet. Run <code>cd gui && npm install && npm run build</code> first.</p>
<p>Or in dev mode: <code>cd gui && npm run dev</code> (on port 5173) while the server runs on this port.</p>
<h2>API Status</h2>
<p>The API server is running. Try <a href="/api/v1/health">/api/v1/health</a></p>
</body>
</html>)html", "text/html");
            res.status = 200;
        }
    });
#endif

    // ================================================================
    // Start server
    // ================================================================
    int port = options.port;
    
    // Try to find an available port
    for (int attempt = 0; attempt < 10; ++attempt) {
        if (svr.bind_to_port("0.0.0.0", port)) {
            break;
        }
        std::cerr << "Port " << port << " is in use, trying " << (port + 1) << "\n";
        port++;
    }
    
    std::string url = "http://localhost:" + std::to_string(port);
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════╗\n";
    std::cout << "║         CoolSolve GUI Server             ║\n";
    std::cout << "╠══════════════════════════════════════════╣\n";
    std::cout << "║  URL: " << std::left << std::setw(35) << url << "║\n";
    std::cout << "║  Press Ctrl+C to stop                    ║\n";
    std::cout << "╚══════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    // Warmup CoolProp in background thread
    std::thread warmupThread([&]() {
        double ms = warmupCoolProp();
        coolpropReady.store(true);
        std::cout << "CoolProp warmed up in " << std::fixed << std::setprecision(0) << ms << " ms\n";
    });
    warmupThread.detach();
    
    // Open browser
    if (options.openBrowser) {
        // Small delay to let server start
        std::thread([url]() {
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
            openBrowser(url);
        }).detach();
    }
    
    // Listen (blocking)
    svr.listen_after_bind();
    
    return 0;
}

} // namespace coolsolve

#endif // COOLSOLVE_GUI
