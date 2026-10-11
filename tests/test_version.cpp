#include <catch2/catch_test_macros.hpp>
#include "coolsolve/version.h"
#include <cctype>
#include <string>

using namespace coolsolve;

TEST_CASE("Version info reports CoolSolve and the linked CoolProp", "[version]") {
    const VersionInfo info = getVersionInfo();

    REQUIRE(info.coolsolve == COOLSOLVE_VERSION);
    REQUIRE_FALSE(info.coolprop.empty());

    // CoolProp commit hash (written by CoolProp's header generation at build time: the full
    // SHA-1 of the FetchContent checkout; a system CoolProp may report an abbreviation).
    REQUIRE(info.coolpropRevision.size() >= 7);
    for (char c : info.coolpropRevision)
        REQUIRE(std::isxdigit(static_cast<unsigned char>(c)));
}

TEST_CASE("versionText has the two-line format parsed by the library tools", "[version]") {
    VersionInfo info{"1.2.3", "8.1.0dev", "75af78166fe6467e9e4ce53ffe083cdc96863363"};
    REQUIRE(versionText(info) ==
            "CoolSolve 1.2.3\n"
            "CoolProp 8.1.0dev (git revision 75af78166fe6467e9e4ce53ffe083cdc96863363)\n");

    // A CoolProp without git information must not break the format.
    info.coolpropRevision.clear();
    REQUIRE(versionText(info) ==
            "CoolSolve 1.2.3\n"
            "CoolProp 8.1.0dev (git revision unknown)\n");
}
