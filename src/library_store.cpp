#include "coolsolve/library_store.h"

#include <algorithm>
#include <cctype>
#include <set>

#ifdef COOLSOLVE_EMBEDDED_LIBRARY
#include "embedded_library.h"
#endif

namespace coolsolve {

using json = nlohmann::json;

namespace {

bool endsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool startsWith(const std::string& s, const std::string& prefix) {
    return s.compare(0, prefix.size(), prefix) == 0;
}

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

/// "a/b/x.eescode" -> "x"
std::string stemOf(const std::string& file) {
    auto slash = file.find_last_of('/');
    std::string base = slash == std::string::npos ? file : file.substr(slash + 1);
    auto dot = base.find_last_of('.');
    return dot == std::string::npos ? base : base.substr(0, dot);
}

/// A model of the menu; `functions` is the number of its FUNCTIONs and PROCEDUREs (what an
/// `$INCLUDE` of the model imports).
json menuEntry(const json& m, int functions) {
    return {{"id", m.value("id", "")},
            {"name", m.value("name", "")},
            {"title", m.value("title", "")},
            {"level", m.value("level", 0)},
            {"kind", m.value("kind", "")},
            {"functions", functions}};
}

/// Levenshtein distance, case-insensitive.
size_t editDistance(const std::string& a, const std::string& b) {
    std::vector<size_t> row(b.size() + 1);
    for (size_t j = 0; j <= b.size(); ++j) row[j] = j;
    for (size_t i = 1; i <= a.size(); ++i) {
        size_t diagonal = row[0];
        row[0] = i;
        for (size_t j = 1; j <= b.size(); ++j) {
            const size_t above = row[j];
            const bool same = std::tolower(static_cast<unsigned char>(a[i - 1])) ==
                              std::tolower(static_cast<unsigned char>(b[j - 1]));
            row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + (same ? 0 : 1)});
            diagonal = above;
        }
    }
    return row[b.size()];
}

/// "hx_ees_function" -> {"hx", "ees", "function"}
std::set<std::string> nameWords(const std::string& name) {
    std::set<std::string> words;
    std::string word;
    for (char c : lower(name) + "_") {
        if (std::isalnum(static_cast<unsigned char>(c))) {
            word += c;
        } else if (!word.empty()) {
            words.insert(word);
            word.clear();
        }
    }
    return words;
}

void sortByTitle(json& models) {
    std::sort(models.begin(), models.end(), [](const json& a, const json& b) {
        return lower(a.value("title", "")) < lower(b.value("title", ""));
    });
}

}  // namespace

LibraryFileSource LibraryFileSource::fromMap(std::map<std::string, std::string> files) {
    auto shared = std::make_shared<std::map<std::string, std::string>>(std::move(files));
    LibraryFileSource src;
    src.read = [shared](const std::string& path, std::string& content) {
        auto it = shared->find(path);
        if (it == shared->end()) return false;
        content = it->second;
        return true;
    };
    src.list = [shared](const std::string& prefix) {
        std::vector<std::string> out;
        for (auto it = shared->lower_bound(prefix); it != shared->end() && startsWith(it->first, prefix); ++it)
            out.push_back(it->first);
        return out;
    };
    return src;
}

LibraryFileSource LibraryFileSource::embedded() {
    LibraryFileSource src;
#ifdef COOLSOLVE_EMBEDDED_LIBRARY
    src.read = [](const std::string& path, std::string& content) {
        const EmbeddedLibraryFile* f = getEmbeddedLibraryFile(path);
        if (!f) return false;
        content.assign(reinterpret_cast<const char*>(f->data), f->size);
        return true;
    };
    src.list = [](const std::string& prefix) { return listEmbeddedLibraryFiles(prefix); };
#endif
    return src;
}

const LibraryStore& embeddedLibrary() {
    static const LibraryStore store(LibraryFileSource::embedded());
    return store;
}

LibraryStore::LibraryStore(LibraryFileSource source) : source_(std::move(source)) {
    if (!source_.read || !source_.list) {
        reason_ = "this build does not embed the CoolSolve Library";
        return;
    }
    std::string text;
    if (!source_.read("library.json", text)) {
        reason_ = "this build does not embed the CoolSolve Library";
        return;
    }
    try {
        index_ = json::parse(text);
        if (!index_.contains("models") || !index_["models"].is_array())
            throw std::runtime_error("no \"models\" array");
    } catch (const std::exception& e) {
        reason_ = std::string("library.json is invalid: ") + e.what();
        return;
    }
    // The taxonomy only gives titles and order: a missing one degrades the menu, not the library.
    if (source_.read("taxonomy.json", text)) {
        try {
            taxonomy_ = json::parse(text);
        } catch (const std::exception&) {
            taxonomy_ = json::object();
        }
    }
    available_ = true;
}

json LibraryStore::snapshot() const {
    if (!available_) return json::object();
    return index_.value("snapshot", json::object());
}

json LibraryStore::menu() const {
    json out = {{"available", available_}, {"categories", json::array()}, {"total", 0}};
    if (!available_) {
        out["reason"] = reason_;
        return out;
    }

    // FUNCTIONs and PROCEDUREs per model ID
    std::map<std::string, int> functionCount;
    if (index_.contains("functions") && index_["functions"].is_array())
        for (const auto& f : index_["functions"]) ++functionCount[f.value("model_id", "")];

    // Verified models by category path ("cycles/organic_rankine" or "buildings")
    std::map<std::string, json> byCategory;
    int total = 0;
    for (const auto& m : index_["models"]) {
        if (m.value("status", "") != "verified") continue;
        const auto count = functionCount.find(m.value("id", ""));
        byCategory[m.value("category", "")].push_back(menuEntry(m, count == functionCount.end() ? 0 : count->second));
        ++total;
    }
    for (auto& kv : byCategory) sortByTitle(kv.second);

    auto take = [&](const std::string& path) {
        json models = json::array();
        auto it = byCategory.find(path);
        if (it != byCategory.end()) {
            models = it->second;
            byCategory.erase(it);
        }
        return models;
    };

    // Declared categories first, in taxonomy order
    if (taxonomy_.contains("categories")) {
        for (const auto& top : taxonomy_["categories"]) {
            const std::string topKey = top.value("key", "");
            json cat = {{"key", topKey}, {"title", top.value("title", topKey)},
                        {"models", take(topKey)}, {"children", json::array()}};
            int count = static_cast<int>(cat["models"].size());
            if (top.contains("children")) {
                for (const auto& sub : top["children"]) {
                    const std::string subKey = sub.value("key", "");
                    json models = take(topKey + "/" + subKey);
                    if (models.empty()) continue;
                    count += static_cast<int>(models.size());
                    cat["children"].push_back({{"key", subKey}, {"title", sub.value("title", subKey)},
                                               {"models", models}});
                }
            }
            if (count == 0) continue;
            cat["count"] = count;
            out["categories"].push_back(cat);
        }
    }
    // Models of a category the taxonomy does not declare (should not happen: the library
    // validates it) are kept visible rather than dropped.
    if (!byCategory.empty()) {
        json other = {{"key", "other"}, {"title", "Other"}, {"models", json::array()}, {"children", json::array()}, {"count", 0}};
        for (auto& kv : byCategory) {
            for (auto& m : kv.second) other["models"].push_back(m);
        }
        sortByTitle(other["models"]);
        other["count"] = static_cast<int>(other["models"].size());
        out["categories"].push_back(other);
    }
    out["total"] = total;
    return out;
}

bool LibraryStore::loadModel(const std::string& idOrName, const std::string& file, LibraryModel& out,
                             std::string& error) const {
    if (!available_) {
        error = "Model library not available: " + reason_;
        return false;
    }
    const json* found = findModel(idOrName);
    if (!found) {
        error = "Unknown library model: " + idOrName;
        return false;
    }

    LibraryModel model;
    model.id = found->value("id", "");
    model.name = found->value("name", "");
    model.title = found->value("title", "");
    model.status = found->value("status", "");
    model.folder = found->value("path", "");
    if (model.folder.empty() || model.folder.find("..") != std::string::npos) {
        error = "Invalid folder for library model " + idOrName;
        return false;
    }
    const std::string dir = model.folder + "/";

    // The .eescode files of the folder: only those can be opened (no path is taken from the client)
    const std::string mainFile = found->value("main_file", "");
    std::vector<std::string> eescodeFiles;
    for (const auto& path : source_.list(dir)) {
        std::string rel = path.substr(dir.size());
        if (rel.find('/') == std::string::npos && endsWith(rel, ".eescode")) eescodeFiles.push_back(rel);
    }
    model.file = file.empty() ? mainFile : file;
    if (std::find(eescodeFiles.begin(), eescodeFiles.end(), model.file) == eescodeFiles.end()) {
        error = "Library model " + model.id + " has no file '" + model.file + "'";
        return false;
    }
    for (const auto& f : eescodeFiles)
        if (f != model.file) model.variants.push_back(f);

    if (!source_.read(dir + model.file, model.eescode)) {
        error = "Cannot read " + dir + model.file;
        return false;
    }

    // Companion files. A variant ("x_coolsolve.eescode") may have its own .initials, .sol and
    // tables or share those of the model: the file's stem comes first, then the model name.
    const std::string stem = stemOf(model.file);
    std::vector<std::string> stems = {stem};
    if (model.name != stem) stems.push_back(model.name);
    for (const auto& s : stems) {
        if (model.initials.empty()) source_.read(dir + s + ".initials", model.initials);
        if (model.sol.empty()) source_.read(dir + s + ".sol", model.sol);

        // Lookup tables "<stem>-<table>.csv": the table name is what follows "<stem>-"
        const std::string prefix = dir + s + "-";
        for (const auto& path : source_.list(prefix)) {
            if (!endsWith(path, ".csv") || path.find('/', dir.size()) != std::string::npos) continue;
            std::string table = path.substr(prefix.size(), path.size() - prefix.size() - 4);
            std::string csv;
            if (!table.empty() && !model.tables.count(table) && source_.read(path, csv))
                model.tables.emplace(table, std::move(csv));
        }
    }
    source_.read(dir + "coolsolve.conf", model.conf);

    out = std::move(model);
    return true;
}

const json* LibraryStore::findModel(const std::string& idOrName) const {
    if (!available_ || idOrName.empty()) return nullptr;
    const std::string wanted = lower(idOrName);   // EES names ignore the case
    for (const auto& m : index_["models"])
        if (lower(m.value("id", "")) == wanted || lower(m.value("name", "")) == wanted) return &m;
    return nullptr;
}

bool LibraryStore::readFile(const std::string& path, std::string& content) const {
    return available_ && path.find("..") == std::string::npos && source_.read(path, content);
}

std::vector<std::string> LibraryStore::closestNames(const std::string& name, size_t count) const {
    if (!available_) return {};
    const std::set<std::string> words = nameWords(name);
    std::vector<std::pair<size_t, std::string>> scored;   // (score, name), lower is closer
    for (const auto& m : index_["models"]) {
        const std::string candidate = m.value("name", "");
        const size_t distance = editDistance(name, candidate);
        size_t shared = 0;
        for (const auto& w : nameWords(candidate)) shared += words.count(w);
        // A few typos, or a word in common: a plausible meant name
        if (distance > std::max<size_t>(3, name.size() / 3) && shared == 0) continue;
        scored.emplace_back(distance - std::min(distance, 4 * shared), candidate);
    }
    std::sort(scored.begin(), scored.end());
    std::vector<std::string> out;
    for (size_t i = 0; i < scored.size() && i < count; ++i) out.push_back(scored[i].second);
    return out;
}

}  // namespace coolsolve
