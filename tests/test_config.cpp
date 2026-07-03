/**
 * Tests for coolsolve.conf loading (loadSolverOptionsFromFile).
 * Run from build directory; examples are expected at ../examples/.
 */

#include <catch2/catch_test_macros.hpp>
#include "coolsolve/solver.h"
#include <cstdlib>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

static fs::path getExamplesDir() {
    const char* env = std::getenv("COOLSOLVE_EXAMPLES_DIR");
    if (env && env[0] != '\0') return fs::path(env);
    return fs::path("../examples");
}

TEST_CASE("Load examples/coolsolve.conf", "[config][solver]") {
    fs::path examplesDir = getExamplesDir();
    fs::path configPath = examplesDir / "coolsolve.conf";
    if (!fs::exists(configPath)) {
        SKIP("Examples coolsolve.conf not found: " << configPath.string());
    }
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    REQUIRE(loaded);
    // File is all comments, so defaults are unchanged (e.g. maxIterations)
    REQUIRE(options.maxIterations == 100);
}

TEST_CASE("Load non-existent config returns false", "[config][solver]") {
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile("/nonexistent/coolsolve.conf", options);
    REQUIRE_FALSE(loaded);
}

TEST_CASE("Config file options are applied", "[config][solver]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_config.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "# test\nmaxIterations = 99\ntolerance = 1e-6\nverbose = true\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.maxIterations == 99);
    REQUIRE(options.tolerance == 1e-6);
    REQUIRE(options.verbose == true);
}

TEST_CASE("Tearing options are loaded from config", "[config][solver][tearing]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_tearing.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "enableTearing = true\n";
    f << "tearingMaxIterations = 200\n";
    f << "tearingMinBlockSize = 2\n";
    f << "tearingInnerIterations = 8\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.enableTearing == true);
    REQUIRE(options.tearingMaxIterations == 200);
    REQUIRE(options.tearingMinBlockSize == 2);
    REQUIRE(options.tearingInnerIterations == 8);
}

TEST_CASE("CoolProp integration options are loaded from config", "[config][solver][coolprop]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_coolprop.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "coolpropBackend = TTSE&HEOS\n";
    f << "coolpropUseAbstractState = true\n";
    f << "coolpropEnableAnalyticalDerivatives = false\n";
    f << "coolpropCacheEnabled = true\n";
    f << "coolpropEnableSuperancillaries = false\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.coolpropConfig.backend == "TTSE&HEOS");
    REQUIRE(options.coolpropConfig.useAbstractState == true);
    REQUIRE(options.coolpropConfig.enableAnalyticalDerivatives == false);
    REQUIRE(options.coolpropConfig.cacheEnabled == true);
    REQUIRE(options.coolpropConfig.enableSuperancillaries == false);
}

TEST_CASE("CoolProp config defaults from file with no CoolProp keys", "[config][solver][coolprop]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_cpdefault.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "maxIterations = 50\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    // CoolProp defaults should be unchanged
    REQUIRE(options.coolpropConfig.backend == "HEOS");
    REQUIRE(options.coolpropConfig.useAbstractState == true);
    REQUIRE(options.coolpropConfig.enableAnalyticalDerivatives == true);
    REQUIRE(options.coolpropConfig.cacheEnabled == true);
    REQUIRE(options.coolpropConfig.enableSuperancillaries == true);
}

TEST_CASE("TrustRegion hybrd Broyden options are loaded from config", "[config][solver][trustregion][hybrd]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_tr_hybrd.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "trBroydenRecomputeInterval = 5\n";
    f << "trBroydenRestartRejects = 3\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.trBroydenRecomputeInterval == 5);
    REQUIRE(options.trBroydenRestartRejects == 3);
}

TEST_CASE("TrustRegion hybrd Broyden options default when not set", "[config][solver][trustregion][hybrd]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_tr_hybrd_default.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "maxIterations = 50\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.trBroydenRecomputeInterval == 0);
    REQUIRE(options.trBroydenRestartRejects == 2);
}

TEST_CASE("Multi-start options are loaded from config", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_multistart.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "multiStartEnabled = false\n";
    f << "multiStartMaxRestarts = 6\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.multiStartEnabled == false);
    REQUIRE(options.multiStartMaxRestarts == 6);
}

TEST_CASE("Multi-start options default when not set", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_multistart_default.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "maxIterations = 50\n";  // unrelated key
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.multiStartEnabled == true);
    REQUIRE(options.multiStartMaxRestarts == 4);
}

TEST_CASE("Multi-start negative max restarts falls back to default", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_multistart_neg.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "multiStartMaxRestarts = -3\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    // Negative values are rejected and reset to the default (4).
    REQUIRE(options.multiStartMaxRestarts == 4);
}

TEST_CASE("Multi-start numCores is loaded from config", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_ms_cores.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "multiStartNumCores = 0\n";  // 0 = auto
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.multiStartNumCores == 0);
}

TEST_CASE("Multi-start numCores negative falls back to sequential default", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_ms_cores_neg.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "multiStartNumCores = -2\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.multiStartNumCores == 1);
}

