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
// Session state (single session for local mode)
// ============================================================================
struct Session {
    std::string eescodeContent;                 // Current .eescode source
    std::string initialsContent;                // Current .initials content
    std::string solContent;                     // Last .sol output
    std::string confContent;                    // coolsolve.conf content
    
    std::string openFilePath;                   // Path of currently open .eescode file
    
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
    Session session;
    
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
    svr.Get("/api/v1/files/eescode", [&](const httplib::Request&, httplib::Response& res) {
        json j = {
            {"content", session.eescodeContent},
            {"filePath", session.openFilePath}
        };
        res.set_content(j.dump(), "application/json");
    });
    
    // PUT /api/v1/files/eescode - Save .eescode content
    svr.Put("/api/v1/files/eescode", [&](const httplib::Request& req, httplib::Response& res) {
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
    svr.Get("/api/v1/files/initials", [&](const httplib::Request&, httplib::Response& res) {
        json j = {{"content", session.initialsContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // PUT /api/v1/files/initials - Save .initials content  
    svr.Put("/api/v1/files/initials", [&](const httplib::Request& req, httplib::Response& res) {
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
    svr.Get("/api/v1/files/sol", [&](const httplib::Request&, httplib::Response& res) {
        json j = {{"content", session.solContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // GET /api/v1/files/conf - Get coolsolve.conf content
    svr.Get("/api/v1/files/conf", [&](const httplib::Request&, httplib::Response& res) {
        json j = {{"content", session.confContent}};
        res.set_content(j.dump(), "application/json");
    });
    
    // PUT /api/v1/files/conf - Save coolsolve.conf content
    svr.Put("/api/v1/files/conf", [&](const httplib::Request& req, httplib::Response& res) {
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
    
    // POST /api/v1/files/open - Open a file from disk (local mode)
    svr.Post("/api/v1/files/open", [&](const httplib::Request& req, httplib::Response& res) {
        try {
            auto body = json::parse(req.body);
            std::string filePath = body.value("path", "");
            
            if (filePath.empty() || !fs::exists(filePath)) {
                res.status = 404;
                json j = {{"error", "File not found: " + filePath}};
                res.set_content(j.dump(), "application/json");
                return;
            }
            
            session.openFilePath = filePath;
            session.eescodeContent = readFileToString(filePath);
            session.initialsContent.clear();
            session.solContent.clear();
            session.confContent.clear();
            
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
    
    // POST /api/v1/files/save - Save to current file path (local mode)
    svr.Post("/api/v1/files/save", [&](const httplib::Request&, httplib::Response& res) {
        if (session.openFilePath.empty()) {
            res.status = 400;
            json j = {{"error", "No file is currently open"}};
            res.set_content(j.dump(), "application/json");
            return;
        }
        
        bool ok = writeStringToFile(session.openFilePath, session.eescodeContent);
        
        // Also save .initials if content exists
        if (!session.initialsContent.empty()) {
            fs::path dir = fs::path(session.openFilePath).parent_path();
            std::string stem = fs::path(session.openFilePath).stem().string();
            writeStringToFile((dir / (stem + ".initials")).string(), session.initialsContent);
        }
        
        // Save .conf if content exists
        if (!session.confContent.empty()) {
            fs::path dir = fs::path(session.openFilePath).parent_path();
            writeStringToFile((dir / "coolsolve.conf").string(), session.confContent);
        }
        
        json j = {{"success", ok}};
        res.set_content(j.dump(), "application/json");
    });
    
    // ================================================================
    // Parse endpoint (for real-time editor feedback)
    // ================================================================
    svr.Post("/api/v1/parse", [&](const httplib::Request& req, httplib::Response& res) {
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
            
            // Write temp files for the runner
            auto tmpDir = fs::temp_directory_path() / "coolsolve_gui";
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
            std::thread([&session, tmpEes, tmpDir, enableTracing, confContent]() {
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
    svr.Get("/api/v1/solve/stream", [&](const httplib::Request&, httplib::Response& res) {
        res.set_header("Cache-Control", "no-cache");
        res.set_header("Connection", "keep-alive");
        res.set_header("X-Accel-Buffering", "no");  // Disable nginx buffering
        
        res.set_chunked_content_provider(
            "text/event-stream",
            [&session](size_t /*offset*/, httplib::DataSink& sink) {
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
    svr.Get("/api/v1/solve/result", [&](const httplib::Request&, httplib::Response& res) {
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
    // Update guesses: copy .sol -> .initials
    // ================================================================
    svr.Post("/api/v1/update-guesses", [&](const httplib::Request&, httplib::Response& res) {
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
    svr.Get("/api/v1/variables", [&](const httplib::Request&, httplib::Response& res) {
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
