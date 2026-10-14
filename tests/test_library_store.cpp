/**
 * Tests for LibraryStore: the index of the CoolSolve Library snapshot (menu of verified models,
 * loading of a model folder with its companion files). Uses an in-memory library, so it does
 * not depend on the embedded snapshot.
 */

#include <catch2/catch_test_macros.hpp>
#include "coolsolve/library_store.h"

using coolsolve::LibraryFileSource;
using coolsolve::LibraryModel;
using coolsolve::LibraryStore;

namespace {

const char* kIndex = R"({
  "snapshot": {"content_hash": "sha256:abc", "models": 4},
  "models": [
    {"id": "CSL-0001", "name": "orc_simple", "title": "Simple ORC", "category": "cycles/organic_rankine",
     "status": "verified", "level": 3, "kind": "steady", "path": "models/cycles/organic_rankine/orc_simple",
     "main_file": "orc_simple.eescode"},
    {"id": "CSL-0002", "name": "blocked_model", "title": "A blocked model", "category": "cycles/organic_rankine",
     "status": "blocked", "level": 2, "kind": "steady", "path": "models/cycles/organic_rankine/blocked_model",
     "main_file": "blocked_model.eescode"},
    {"id": "CSL-0003", "name": "comfort", "title": "Comfort", "category": "buildings",
     "status": "verified", "level": 1, "kind": "steady", "path": "models/buildings/comfort",
     "main_file": "comfort.eescode"},
    {"id": "CSL-0004", "name": "alpha", "title": "alpha first", "category": "cycles/organic_rankine",
     "status": "verified", "level": 1, "kind": "steady", "path": "models/cycles/organic_rankine/alpha",
     "main_file": "alpha.eescode"}
  ]
})";

const char* kTaxonomy = R"({
  "categories": [
    {"key": "fundamentals", "title": "Fundamentals", "children": [{"key": "properties", "title": "Properties"}]},
    {"key": "cycles", "title": "Cycles", "children": [
        {"key": "steam_power", "title": "Steam power"},
        {"key": "organic_rankine", "title": "Organic Rankine"}]},
    {"key": "buildings", "title": "Buildings"}
  ]
})";

LibraryStore makeStore() {
    return LibraryStore(LibraryFileSource::fromMap({
        {"library.json", kIndex},
        {"taxonomy.json", kTaxonomy},
        {"models/cycles/organic_rankine/orc_simple/orc_simple.eescode", "T = 1"},
        {"models/cycles/organic_rankine/orc_simple/orc_simple_coolsolve.eescode", "T = 2"},
        {"models/cycles/organic_rankine/orc_simple/orc_simple.initials", "T=5"},
        {"models/cycles/organic_rankine/orc_simple/orc_simple.sol", "sol"},
        {"models/cycles/organic_rankine/orc_simple/coolsolve.conf", "maxIterations = 50"},
        {"models/cycles/organic_rankine/orc_simple/orc_simple-map.csv", "a,b\n1,2"},
        {"models/cycles/organic_rankine/orc_simple/orc_simple-map-2.csv", "c\n3"},
        {"models/cycles/organic_rankine/orc_simple/other-table.csv", "not a table of this model"},
        {"models/cycles/organic_rankine/orc_simple/README.md", "# readme"},
        {"models/cycles/organic_rankine/orc_simple/figures/f.png", "png"},
        {"models/cycles/organic_rankine/blocked_model/blocked_model.eescode", "x = y"},
        {"models/buildings/comfort/comfort.eescode", "c = 1"},
        {"models/buildings/comfort/comfort-clo.csv", "model table"},
        {"models/buildings/comfort/alt_comfort.eescode", "c = 2"},
        {"models/buildings/comfort/alt_comfort-clo.csv", "variant table"},
    }));
}

}  // namespace

TEST_CASE("Library menu lists only verified models, by taxonomy", "[library]") {
    auto store = makeStore();
    REQUIRE(store.available());
    auto menu = store.menu();

    REQUIRE(menu["total"] == 3);
    // Taxonomy order, empty categories (fundamentals) and empty sub-categories (steam_power) left out
    REQUIRE(menu["categories"].size() == 2);
    REQUIRE(menu["categories"][0]["key"] == "cycles");
    REQUIRE(menu["categories"][0]["count"] == 2);
    REQUIRE(menu["categories"][0]["children"].size() == 1);
    auto orc = menu["categories"][0]["children"][0];
    REQUIRE(orc["key"] == "organic_rankine");
    // The blocked model is not in the menu; the others are sorted by title, case-insensitively
    REQUIRE(orc["models"].size() == 2);
    REQUIRE(orc["models"][0]["id"] == "CSL-0004");
    REQUIRE(orc["models"][1]["id"] == "CSL-0001");
    // A category without sub-categories holds its models directly
    REQUIRE(menu["categories"][1]["key"] == "buildings");
    REQUIRE(menu["categories"][1]["models"][0]["id"] == "CSL-0003");
}

TEST_CASE("Library snapshot information comes from library.json", "[library]") {
    REQUIRE(makeStore().snapshot()["content_hash"] == "sha256:abc");
}

TEST_CASE("Loading a library model brings its whole folder", "[library]") {
    auto store = makeStore();
    LibraryModel m;
    std::string err;
    REQUIRE(store.loadModel("CSL-0001", "", m, err));
    REQUIRE(m.name == "orc_simple");
    REQUIRE(m.eescode == "T = 1");
    REQUIRE(m.initials == "T=5");
    REQUIRE(m.sol == "sol");
    REQUIRE(m.conf == "maxIterations = 50");
    // <name>-<table>.csv: the table name is what follows the first hyphen; other models' tables are ignored
    REQUIRE(m.tables.size() == 2);
    REQUIRE(m.tables.at("map") == "a,b\n1,2");
    REQUIRE(m.tables.at("map-2") == "c\n3");
    REQUIRE(m.variants == std::vector<std::string>{"orc_simple_coolsolve.eescode"});
}

TEST_CASE("A library model loads by name and a variant by file name", "[library]") {
    auto store = makeStore();
    LibraryModel m;
    std::string err;
    REQUIRE(store.loadModel("orc_simple", "orc_simple_coolsolve.eescode", m, err));
    REQUIRE(m.eescode == "T = 2");
    // The variant has no .initials of its own: the model's are used, and so are its tables
    REQUIRE(m.initials == "T=5");
    REQUIRE(m.tables.size() == 2);
    REQUIRE(m.variants == std::vector<std::string>{"orc_simple.eescode"});
}

TEST_CASE("A variant's own table wins over the model's table of the same name", "[library]") {
    auto store = makeStore();
    LibraryModel m;
    std::string err;
    // "alt_comfort-clo.csv" sorts before "comfort-clo.csv": the priority must not depend on the order
    REQUIRE(store.loadModel("CSL-0003", "alt_comfort.eescode", m, err));
    REQUIRE(m.tables.size() == 1);
    REQUIRE(m.tables.at("clo") == "variant table");

    REQUIRE(store.loadModel("CSL-0003", "", m, err));
    REQUIRE(m.tables.at("clo") == "model table");
    REQUIRE(m.variants == std::vector<std::string>{"alt_comfort.eescode"});
}

TEST_CASE("A blocked model can still be opened", "[library]") {
    auto store = makeStore();
    LibraryModel m;
    std::string err;
    REQUIRE(store.loadModel("CSL-0002", "", m, err));
    REQUIRE(m.status == "blocked");
    REQUIRE(m.eescode == "x = y");
    REQUIRE(m.initials.empty());
}

TEST_CASE("Loading a library model fails with a clear error", "[library]") {
    auto store = makeStore();
    LibraryModel m;
    std::string err;

    REQUIRE_FALSE(store.loadModel("CSL-9999", "", m, err));
    REQUIRE(err.find("Unknown library model") != std::string::npos);

    // Only the .eescode files listed in the model folder can be opened: no path from the client
    REQUIRE_FALSE(store.loadModel("CSL-0001", "../comfort/comfort.eescode", m, err));
    REQUIRE(err.find("has no file") != std::string::npos);
    REQUIRE_FALSE(store.loadModel("CSL-0001", "README.md", m, err));
}

TEST_CASE("A library without library.json is reported as unavailable", "[library]") {
    LibraryStore store(LibraryFileSource::fromMap({}));
    REQUIRE_FALSE(store.available());
    REQUIRE(store.menu()["available"] == false);
    REQUIRE(store.menu()["total"] == 0);

    LibraryModel m;
    std::string err;
    REQUIRE_FALSE(store.loadModel("CSL-0001", "", m, err));
    REQUIRE(err.find("not available") != std::string::npos);

    LibraryStore noSource{LibraryFileSource{}};
    REQUIRE_FALSE(noSource.available());
}

TEST_CASE("A corrupt library.json is reported, not thrown", "[library]") {
    LibraryStore store(LibraryFileSource::fromMap({{"library.json", "{not json"}}));
    REQUIRE_FALSE(store.available());
    REQUIRE(store.unavailableReason().find("invalid") != std::string::npos);

    LibraryStore noModels(LibraryFileSource::fromMap({{"library.json", "{}"}}));
    REQUIRE_FALSE(noModels.available());
}
