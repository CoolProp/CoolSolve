/**
 * Tests for coolsolve.conf loading (loadSolverOptionsFromFile).
 * Run from build directory; examples are expected at ../examples/.
 */

#include <catch2/catch_test_macros.hpp>
#include "coolsolve/solver.h"
#include "coolsolve/integral/integrator.h"
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

TEST_CASE("Integral options are loaded from config", "[config][solver][integral]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_integral.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "integralMethod = RK45\n";
    f << "integralFixedStep = 0.0\n";
    f << "integralMaxSteps = 5000\n";
    f << "integralRelTol = 1e-8\n";
    f << "integralAbsTol = 1e-10\n";
    f << "integralMinStep = 1e-6\n";
    f << "integralMaxStep = 0.5\n";
    f << "integralRichardson = true\n";
    f << "integralOutputInterval = 0.25\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.integralMethod == "rk45");
    REQUIRE(options.integralFixedStep == 0.0);
    REQUIRE(options.integralMaxSteps == 5000);
    REQUIRE(options.integralRelTol == 1e-8);
    REQUIRE(options.integralAbsTol == 1e-10);
    REQUIRE(options.integralMinStep == 1e-6);
    REQUIRE(options.integralMaxStep == 0.5);
    REQUIRE(options.integralRichardson == true);
    REQUIRE(options.integralOutputInterval == 0.25);

    // The method string round-trips through parseIntegralMethod.
    coolsolve::IntegratorOptions::Method m;
    REQUIRE(coolsolve::parseIntegralMethod(options.integralMethod, m));
    REQUIRE(m == coolsolve::IntegratorOptions::RK45);
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
    // Legacy `multiStartEnabled = false` must map to Never.
    REQUIRE(options.multiStartMode == coolsolve::MultiStartMode::Never);
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

TEST_CASE("Multi-start numCores negative falls back to default", "[config][solver][multistart]") {
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
    // Negative values are rejected and reset to the default (4).
    REQUIRE(options.multiStartNumCores == 4);
}

TEST_CASE("Multi-start numCores default is 4", "[config][solver][multistart]") {
    // The built-in default must be 4 (multi-start is failure-only, so the
    // extra threads are amortised across expensive candidate re-solves).
    coolsolve::SolverOptions options;
    REQUIRE(options.multiStartNumCores == 4);
}

TEST_CASE("Multi-start mode is loaded from config", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_ms_mode.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "multiStartMode = always\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.multiStartMode == coolsolve::MultiStartMode::Always);
    // The legacy bool is kept in sync.
    REQUIRE(options.multiStartEnabled == true);
}

TEST_CASE("Multi-start mode accepts deepsearch/never (and aliases)", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();

    auto loadMode = [&](const std::string& val) {
        fs::path p = tmpDir / "coolsolve_test_ms_mode_variants.conf";
        std::ofstream f(p);
        f << "multiStartMode = " << val << "\n";
        f.close();
        coolsolve::SolverOptions options;
        coolsolve::loadSolverOptionsFromFile(p.string(), options);
        fs::remove(p);
        return options.multiStartMode;
    };

    REQUIRE(loadMode("deepsearch")  == coolsolve::MultiStartMode::InDeepSearch);
    REQUIRE(loadMode("InDeepSearch") == coolsolve::MultiStartMode::InDeepSearch);
    REQUIRE(loadMode("in deep search") == coolsolve::MultiStartMode::InDeepSearch);
    REQUIRE(loadMode("never")       == coolsolve::MultiStartMode::Never);
    REQUIRE(loadMode("off")         == coolsolve::MultiStartMode::Never);
    // Invalid value falls back to InDeepSearch (the default).
    REQUIRE(loadMode("banana")      == coolsolve::MultiStartMode::InDeepSearch);
}

TEST_CASE("Multi-start legacy bool maps to mode", "[config][solver][multistart]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_ms_legacy.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    // Legacy boolean: `false` should map to Never for backward compatibility.
    f << "multiStartEnabled = false\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.multiStartMode == coolsolve::MultiStartMode::Never);
    REQUIRE(options.multiStartEnabled == false);
}

TEST_CASE("Multi-start default mode is InDeepSearch", "[config][solver][multistart]") {
    coolsolve::SolverOptions options;
    REQUIRE(options.multiStartMode == coolsolve::MultiStartMode::InDeepSearch);
    // isMultiStartActive() must be false in a normal solve (deepSearch off).
    REQUIRE_FALSE(options.isMultiStartActive());
    options.deepSearch = true;
    REQUIRE(options.isMultiStartActive());
}

TEST_CASE("Deep search pipeline is loaded from config", "[config][solver][deepsearch]") {
    fs::path tmpDir = fs::temp_directory_path();
    fs::path configPath = tmpDir / "coolsolve_test_deepsearch.conf";
    std::ofstream f(configPath);
    REQUIRE(f.is_open());
    f << "deepSearchPipeline = Newton, LevenbergMarquardt, BisectionND\n";
    f << "deepSearchPipelineMode = parallel\n";
    f.close();
    coolsolve::SolverOptions options;
    bool loaded = coolsolve::loadSolverOptionsFromFile(configPath.string(), options);
    fs::remove(configPath);
    REQUIRE(loaded);
    REQUIRE(options.deepSearchPipeline.size() == 3);
    REQUIRE(options.deepSearchPipeline[0] == coolsolve::SolverStrategy::Newton);
    REQUIRE(options.deepSearchPipeline[1] == coolsolve::SolverStrategy::LevenbergMarquardt);
    REQUIRE(options.deepSearchPipeline[2] == coolsolve::SolverStrategy::BisectionND);
    REQUIRE(options.deepSearchPipelineMode == coolsolve::SolverPipelineMode::Parallel);
}

TEST_CASE("Deep search pipeline default is the full sequential chain", "[config][solver][deepsearch]") {
    coolsolve::SolverOptions options;
    REQUIRE(options.deepSearchPipeline.size() == 7);
    REQUIRE(options.deepSearchPipelineMode == coolsolve::SolverPipelineMode::Sequential);
    REQUIRE(options.deepSearch == false);  // never persisted from conf
}

