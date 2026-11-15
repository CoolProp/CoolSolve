/**
 * Tests for the ZIP bundle (zip_bundle.h): the archives of the zip tools users have (the fixtures of
 * tools/zip_fixtures.py in tests/fixtures/zip/: Info-ZIP, 7-Zip, bsdtar, Python, the layouts of macOS and
 * Windows), the round trip of the bundle *Save* writes, and the messages of the archives that cannot be
 * read (encrypted, another compression method, not a ZIP, damaged, unsafe names, too large, not UTF-8).
 * The upload route itself is tested by tests/test_bundle_api.py.
 *
 * Run with: ./coolsolve_tests "[zip]"
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "coolsolve/zip_bundle.h"

#include <filesystem>
#include <fstream>
#include <sstream>

using namespace coolsolve;
using Catch::Matchers::ContainsSubstring;
namespace fs = std::filesystem;

namespace {

fs::path fixtureDir() {
    const fs::path here = fs::path(__FILE__).parent_path() / "fixtures" / "zip";
    if (fs::exists(here)) return here;
    return fs::path("..") / "tests" / "fixtures" / "zip";
}

std::string readFixture(const std::string& name) {
    std::ifstream file(fixtureDir() / name, std::ios::binary);
    REQUIRE(file.is_open());
    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

/// The files of every fixture (MODEL and INITIALS of tools/zip_fixtures.py).
std::string fixtureModel() {
    std::string s = "\"The model of the ZIP fixtures (tests/test_zip_bundle.cpp): T in \xC2\xB0" "C\"\n";
    for (int i = 1; i <= 20; ++i) s += "y_" + std::to_string(i) + " = x^" + std::to_string(i) + " + " + std::to_string(i) + "\n";
    return s + "x = 2\n";
}

std::string fixtureInitials() {
    std::string s;
    for (int i = 1; i <= 20; ++i) s += "y_" + std::to_string(i) + "=1\n";
    return s;
}

}  // namespace

TEST_CASE("ZIP bundle: archives of the usual zip tools", "[zip]") {
    // Deflate with extra fields (the zip command, as the archive of the bug report), data descriptors
    // (streamed), 7-Zip, libarchive and ZIP64 records: the same two files
    for (const char* name : {"infozip_deflate.zip", "infozip_stream.zip", "7zip_deflate.zip", "bsdtar.zip",
                             "python_zip64.zip"}) {
        INFO(name);
        const auto files = extractZipBundle(readFixture(name));
        REQUIRE(files.size() == 2);
        CHECK(files.at("model.eescode") == fixtureModel());
        CHECK(files.at("model.initials") == fixtureInitials());
    }
}

TEST_CASE("ZIP bundle: the ZIP of a folder", "[zip]") {
    SECTION("zip -r: the folder entry and the folder of the names are left out") {
        const auto files = extractZipBundle(readFixture("infozip_folder.zip"));
        REQUIRE(files.size() == 2);
        CHECK(files.at("model.eescode") == fixtureModel());
        CHECK(files.count("model.initials") == 1);
    }
    SECTION("macOS' Compress: __MACOSX/, ._* and .DS_Store are not files of the model") {
        const auto files = extractZipBundle(readFixture("macos_folder.zip"));
        REQUIRE(files.size() == 2);
        CHECK(files.at("model.eescode") == fixtureModel());
        CHECK(files.count("model.initials") == 1);
    }
    SECTION("files in several folders keep their names") {
        const auto files = extractZipBundle(createZipBundle({{"a/model.eescode", "x=1"}, {"b/model.initials", "x=1"}}));
        CHECK(files.count("a/model.eescode") == 1);
        CHECK(files.count("b/model.initials") == 1);
    }
    SECTION("the debug folder of a bundle is not a folder of the ZIP") {
        const auto files = extractZipBundle(createZipBundle({{"debug_output/report.md", "# Report"}}));
        CHECK(files.count("debug_output/report.md") == 1);
    }
}

TEST_CASE("ZIP bundle: names", "[zip]") {
    SECTION("Windows' Compressed folder: code page 437 names become UTF-8") {
        const auto files = extractZipBundle(readFixture("windows_cp437.zip"));
        CHECK(files.count("\xC3\x89" "changeur.eescode") == 1);  // Échangeur
        CHECK(files.count("\xC3\x89" "changeur.initials") == 1);
    }
    SECTION("UTF-8 names are kept") {
        const std::string name = "\xC3\xA9" "change.eescode";  // échange
        CHECK(extractZipBundle(createZipBundle({{name, "x=1"}})).count(name) == 1);
    }
    SECTION("backslashes and ./ are separators of folders") {
        const auto files = extractZipBundle(createZipBundle({{"./model.eescode", "x=1"}, {"debug_output\\report.md", "#"}}));
        CHECK(files.count("model.eescode") == 1);
        CHECK(files.count("debug_output/report.md") == 1);
    }
}

TEST_CASE("ZIP bundle: round trip of the bundle Save writes", "[zip]") {
    const std::vector<std::pair<std::string, std::string>> bundle = {
        {"model.eescode", fixtureModel()},
        {"model.initials", fixtureInitials()},
        {"empty.sol", ""},
        {"debug_output/report.md", "# Report\n"},
    };
    const auto files = extractZipBundle(createZipBundle(bundle));
    REQUIRE(files.size() == bundle.size());
    for (const auto& [name, content] : bundle) CHECK(files.at(name) == content);
    CHECK(extractZipBundle(createZipBundle({})).empty());
}

TEST_CASE("ZIP bundle: archives that cannot be read", "[zip]") {
    SECTION("password-protected") {
        CHECK_THROWS_WITH(extractZipBundle(readFixture("infozip_encrypted.zip")),
                          ContainsSubstring("'model.eescode' is password-protected"));
    }
    SECTION("another compression method than Deflate, named in the message") {
        CHECK_THROWS_WITH(extractZipBundle(readFixture("7zip_bzip2.zip")),
                          ContainsSubstring("compressed with BZIP2") && ContainsSubstring("the Deflate method"));
        CHECK_THROWS_WITH(extractZipBundle(readFixture("python_lzma.zip")), ContainsSubstring("compressed with LZMA"));
    }
    SECTION("not a ZIP archive") {
        CHECK_THROWS_WITH(extractZipBundle("x = 2\n"), ContainsSubstring("not a ZIP archive, or it is damaged"));
        CHECK_THROWS_WITH(extractZipBundle(""), ContainsSubstring("not a ZIP archive"));
        CHECK_THROWS_WITH(extractZipBundle(std::string("7z\xBC\xAF\x27\x1C", 6) + std::string(32, '\0')),
                          ContainsSubstring("is a 7z archive, not a ZIP archive"));
        CHECK_THROWS_WITH(extractZipBundle(std::string("Rar!\x1A\x07\x01\x00", 8) + std::string(32, '\0')),
                          ContainsSubstring("is a RAR archive"));
    }
    SECTION("damaged: truncated, or a byte of the compressed data changed") {
        const std::string zip = readFixture("infozip_deflate.zip");
        CHECK_THROWS_AS(extractZipBundle(zip.substr(0, zip.size() / 2)), ZipBundleError);
        // A byte in the Deflate stream of model.eescode: after the local header, its name and extra field
        auto u16 = [&](std::size_t at) { return static_cast<unsigned char>(zip[at]) | (static_cast<unsigned char>(zip[at + 1]) << 8); };
        const std::size_t data = 30 + u16(26) + u16(28);
        std::string altered = zip;
        altered[data + 10] = static_cast<char>(altered[data + 10] ^ 0x55);
        CHECK_THROWS_WITH(extractZipBundle(altered), ContainsSubstring("'model.eescode' cannot be extracted"));
    }
    SECTION("a name outside the folder of the archive") {
        CHECK_THROWS_WITH(extractZipBundle(createZipBundle({{"../evil.eescode", "x=1"}})),
                          ContainsSubstring("names a file outside its folder ('../evil.eescode')"));
        CHECK_THROWS_WITH(extractZipBundle(createZipBundle({{"debug_output/../../evil.md", "#"}})),
                          ContainsSubstring("outside its folder"));
        CHECK_THROWS_WITH(extractZipBundle(createZipBundle({{"/etc/evil.eescode", "x=1"}})),
                          ContainsSubstring("outside its folder"));
        CHECK_THROWS_WITH(extractZipBundle(createZipBundle({{"C:\\evil.eescode", "x=1"}})),
                          ContainsSubstring("outside its folder"));
        CHECK(extractZipBundle(createZipBundle({{"a..b.eescode", "x=1"}})).count("a..b.eescode") == 1);
    }
    SECTION("larger than the limit once uncompressed") {
        const std::string zip = readFixture("infozip_deflate.zip");  // 520 bytes uncompressed
        CHECK_THROWS_WITH(extractZipBundle(zip, 500), ContainsSubstring("exceed"));
        CHECK(extractZipBundle(zip, 520).size() == 2);
    }
}

TEST_CASE("ZIP bundle: text files must be UTF-8", "[zip]") {
    CHECK_NOTHROW(checkUtf8Text("model.eescode", fixtureModel()));
    CHECK_NOTHROW(checkUtf8Text("model.eescode", ""));
    CHECK_NOTHROW(checkUtf8Text("model.eescode", "\xF0\x9F\x94\xA5 = 1"));  // U+1F525, 4 bytes
    // A model saved by EES (Windows ANSI): ° is the byte 0xB0
    CHECK_THROWS_WITH(checkUtf8Text("model.eescode", "x = 2\nT = 20 \"\xB0" "C\"\n"),
                      ContainsSubstring("'model.eescode' is not UTF-8 text (byte 0xB0 on line 2)"));
    CHECK_THROWS_AS(checkUtf8Text("a", "\xC3"), ZipBundleError);              // truncated character
    CHECK_THROWS_AS(checkUtf8Text("a", "\xC0\xAF"), ZipBundleError);          // overlong '/'
    CHECK_THROWS_AS(checkUtf8Text("a", "\xED\xA0\x80"), ZipBundleError);      // surrogate U+D800
    CHECK_THROWS_AS(checkUtf8Text("a", "\xF4\x90\x80\x80"), ZipBundleError);  // above U+10FFFF
}
