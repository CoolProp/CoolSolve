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
