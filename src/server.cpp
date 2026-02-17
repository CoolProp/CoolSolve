#ifdef COOLSOLVE_GUI

#include "coolsolve/server.h"
#include "coolsolve/runner.h"
#include "coolsolve/parser.h"
#include "coolsolve/ir.h"
#include "coolsolve/structural_analysis.h"
#include "coolsolve/evaluator.h"
#include "coolsolve/solver.h"
#include "coolsolve/constants.h"

#include <httplib.h>
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
// Snapshot for Back button (one-level undo for model loading)
// ============================================================================
struct SessionSnapshot {
    std::string eescodeContent;
    std::string initialsContent;
    std::string solContent;
    std::string confContent;
    std::string openFilePath;
    std::string debugDir;
    bool hasResult = false;
    SolveResult lastResult;
    CoolSolveRunner::PipelineTiming lastTiming;
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
    bool hasResult = false;
    
    // SSE progress events (protected by mutex + condition variable)
    std::mutex progressMutex;
    std::condition_variable progressCV;
    std::vector<std::string> progressEvents;
    
    // Save current state as snapshot (for Back button)
    void saveSnapshot() {
        auto snap = std::make_unique<SessionSnapshot>();
        snap->eescodeContent = eescodeContent;
        snap->initialsContent = initialsContent;
        snap->solContent = solContent;
        snap->confContent = confContent;
        snap->openFilePath = openFilePath;
        snap->debugDir = debugDir;
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
        debugDir = previousState->debugDir;
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
static json solveResultToJSON(const SolveResult& result, const CoolSolveRunner::PipelineTiming& timing) {
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
            {"filePath", session.openFilePath}
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
        json j = {{"success", true}};
        res.set_content(j.dump(), "application/json");
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
                
                json variables = json::array();
                for (const auto& [name, info] : ir.getVariables()) {
                    if (!Constants::isConstant(name)) {
                        json varObj;
                        varObj["name"] = name;
                        varObj["units"] = info.units;
                        varObj["isArray"] = (name.find('[') != std::string::npos);
                        variables.push_back(varObj);
                    }
                }
                j["variables"] = variables;
                j["variableCount"] = ir.getNonConstantVariableCount();
                j["equationCount"] = ir.getEquationCount();
                j["isSquare"] = ir.isSquare();
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
                    
                    // Apply CoolProp config
                    CoolPropConfig cpConfig;
                    applyCoolPropConfig(cpConfig);
                    
                    session.addProgressEvent("{\"type\":\"progress\",\"message\":\"Parsing...\"}");
                    
                    bool success = runner.run(solverOpts, enableTracing);
                    
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
                        session.hasResult = true;
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
                    
                    // Build final result JSON for the done event
                    json resultJson = solveResultToJSON(runner.getSolveResult(), runner.getTiming());
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
                    }
                    
                    // Send final event with embedded result
                    json finalEvt;
                    finalEvt["type"] = success ? "done" : "error";
                    finalEvt["message"] = success ? "Solve completed successfully" : "Solve failed";
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
        json j = solveResultToJSON(session.lastResult, session.lastTiming);
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
        std::string stem = "model";
        if (!session.openFilePath.empty())
            stem = fs::path(session.openFilePath).stem().string();
        
        std::vector<std::pair<std::string, std::string>> files;
        if (!session.eescodeContent.empty())
            files.push_back({stem + ".eescode", session.eescodeContent});
        if (!session.initialsContent.empty())
            files.push_back({stem + ".initials", session.initialsContent});
        if (!session.solContent.empty())
            files.push_back({stem + ".sol", session.solContent});
        if (!session.confContent.empty())
            files.push_back({"coolsolve.conf", session.confContent});
        
        // Include debug files if they exist
        if (!session.debugDir.empty() && fs::exists(session.debugDir)) {
            for (const auto& entry : fs::directory_iterator(session.debugDir)) {
                if (entry.is_regular_file()) {
                    std::string content = readFileToString(entry.path().string());
                    files.push_back({"debug_output/" + entry.path().filename().string(), content});
                }
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
        for (const auto& [name, file] : req.files) {
            if (endsWith(file.filename, ".zip") && file.content.size() > 4) {
                zipContent = file.content;
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
            }
        }
        if (hasDebug) fileList.push_back("debug_output/");
        
        json j = {{"success", true}, {"files", fileList}};
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
