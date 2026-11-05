#include "coolsolve/include_resolver.h"

#include "coolsolve/library_store.h"
#include "coolsolve/lookup_table.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

namespace coolsolve {

namespace {

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

bool iequals(const std::string& a, const std::string& b) { return lower(a) == lower(b); }

std::string trim(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

/// The text without its comments: `{...}`, `"..."` and `// ...` (a `'...'` string is kept).
std::string withoutComments(const std::string& text) {
    std::string out;
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (c == '\'') {
            const auto end = text.find('\'', i + 1);
            const size_t stop = end == std::string::npos ? text.size() : end + 1;
            out += text.substr(i, stop - i);
            i = stop - 1;
        } else if (c == '{' || c == '"') {
            const auto end = text.find(c == '{' ? '}' : '"', i + 1);
            if (end == std::string::npos) break;
            i = end;
        } else if (c == '/' && i + 1 < text.size() && text[i + 1] == '/') {
            break;
        } else {
            out += c;
        }
    }
    return out;
}

bool endsWithIgnoreCase(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() && iequals(s.substr(s.size() - suffix.size()), suffix);
}

bool readTextFile(const fs::path& path, std::string& content) {
    std::ifstream file(path, std::ios::binary);
    if (!file) return false;
    std::ostringstream buffer;
    buffer << file.rdbuf();
    content = buffer.str();
    return true;
}

/// A file of a folder on disk with its tables `<stem>-<table>.csv`.
bool readFolderFile(const fs::path& path, IncludedSource& out) {
    if (!readTextFile(path, out.content)) return false;
    out.tables = readCompanionTableFiles(path.parent_path().string(), path.stem().string());
    return true;
}

}  // namespace

std::string IncludedSource::key() const {
    return lower(fromLibrary ? "library:" + libraryFolder + "/" + fileName : "project:" + fileName);
}

std::string IncludedSource::describe() const {
    if (!fromLibrary) return fileName;
    if (fileName == libraryName + ".eescode") return "library model " + libraryName + " (" + libraryId + ")";
    return fileName + " of library model " + libraryName + " (" + libraryId + ")";
}

bool parseIncludeName(const std::string& content, std::string& name, std::string& error) {
    std::string text = trim(withoutComments(content));
    if (text.size() >= 2 && text.front() == '\'' && text.back() == '\'') text = trim(text.substr(1, text.size() - 2));
    if (text.empty()) {
        error = "$INCLUDE needs the name of a library model or of a file of the model folder "
                "(e.g. $INCLUDE hx_ees_function)";
        return false;
    }
    if (text.find('/') != std::string::npos || text.find('\\') != std::string::npos ||
        text.find("..") != std::string::npos || text.find(':') != std::string::npos) {
        error = "$INCLUDE " + text + ": only a file of the model folder or a library model can be included, "
                "by its name without path (e.g. $INCLUDE my_functions for my_functions.eescode)";
        return false;
    }
    if (endsWithIgnoreCase(text, ".eescode")) text = trim(text.substr(0, text.size() - 8));
    for (const char* ext : {".lib", ".txt", ".ees", ".fdl", ".dll", ".lkt", ".tpl"}) {
        if (endsWithIgnoreCase(text, ext)) {
            error = "$INCLUDE " + text + ": CoolSolve includes .eescode files and library models only; "
                    "convert the EES file into " + text.substr(0, text.size() - std::string(ext).size()) +
                    ".eescode in the model folder (docs/language_reference.md, $INCLUDE)";
            return false;
        }
    }
    if (text.empty() || text.find('\'') != std::string::npos) {
        error = "$INCLUDE: invalid name '" + content + "'";
        return false;
    }
    name = text;
    return true;
}

std::shared_ptr<const IncludeResolver> IncludeResolver::forDirectory(const std::string& dir, const std::string& mainFile,
                                                                     const LibraryStore* library) {
    auto resolver = std::shared_ptr<IncludeResolver>(new IncludeResolver());
    const fs::path folder = dir.empty() ? fs::path(".") : fs::path(dir);
    resolver->project_ = [folder](const std::string& fileName, IncludedSource& out) {
        std::error_code ec;
        if (fs::is_regular_file(folder / fileName, ec)) return readFolderFile(folder / fileName, out);
        // The case of a name may differ from the file's (models written on Windows)
        for (const auto& entry : fs::directory_iterator(folder, ec)) {
            if (entry.is_regular_file(ec) && iequals(entry.path().filename().string(), fileName))
                return readFolderFile(entry.path(), out);
        }
        return false;
    };
    resolver->mainFile_ = mainFile;
    resolver->library_ = library;
    return resolver;
}

std::shared_ptr<const IncludeResolver> IncludeResolver::forFiles(std::map<std::string, std::string> files,
                                                                 std::map<std::string, std::string> tables,
                                                                 const std::string& mainFile,
                                                                 const LibraryStore* library) {
    auto resolver = std::shared_ptr<IncludeResolver>(new IncludeResolver());
    auto sharedFiles = std::make_shared<const std::map<std::string, std::string>>(std::move(files));
    auto sharedTables = std::make_shared<const std::map<std::string, std::string>>(std::move(tables));
    resolver->project_ = [sharedFiles, sharedTables](const std::string& fileName, IncludedSource& out) {
        for (const auto& [file, text] : *sharedFiles) {
            if (!iequals(file, fileName)) continue;
            out.content = text;
            const std::string prefix = lower(fs::path(file).stem().string()) + "-";
            for (const auto& [name, csv] : *sharedTables)
                if (name.size() > prefix.size() && lower(name.substr(0, prefix.size())) == prefix)
                    out.tables[name.substr(prefix.size())] = csv;
            return true;
        }
        return false;
    };
    resolver->mainFile_ = mainFile;
    resolver->library_ = library;
    return resolver;
}

std::shared_ptr<const IncludeResolver> IncludeResolver::forLibraryModel(const std::string& libraryId,
                                                                        const std::string& mainFile,
                                                                        const LibraryStore* library) {
    auto resolver = std::shared_ptr<IncludeResolver>(new IncludeResolver());
    resolver->rootLibraryId_ = libraryId;
    resolver->mainFile_ = mainFile;
    resolver->library_ = library;
    return resolver;
}

std::shared_ptr<const IncludeResolver> IncludeResolver::libraryOnly(const LibraryStore* library) {
    auto resolver = std::shared_ptr<IncludeResolver>(new IncludeResolver());
    resolver->library_ = library;
    return resolver;
}

bool IncludeResolver::readLibraryModelFile(const std::string& libraryId, const std::string& fileName,
                                           IncludedSource& out) const {
    if (!library_) return false;
    LibraryModel model;
    std::string error;
    if (!library_->loadModel(libraryId, fileName, model, error)) return false;
    IncludedSource source;
    source.fromLibrary = true;
    source.libraryId = model.id;
    source.libraryName = model.name;
    source.libraryFolder = model.folder;
    source.fileName = model.file;
    source.content = std::move(model.eescode);
    source.tables = std::move(model.tables);
    out = std::move(source);
    return true;
}

bool IncludeResolver::resolve(const std::string& name, const IncludedSource* from, IncludedSource& out,
                              std::string& error) const {
    const std::string fileName = name + ".eescode";
    // The library model whose folder holds the includer: the includer itself, or the model opened from the library
    const std::string folderModel = from ? (from->fromLibrary ? from->libraryId : "") : rootLibraryId_;
    const bool ownFile = !from && iequals(fileName, mainFile_);
    if (ownFile) {
        error = "Cannot $INCLUDE '" + name + "': " + fileName + " is the model itself";
        return false;
    }

    // 1. The folder of the includer: a file of the project, or of the library model
    if (!folderModel.empty()) {
        if (readLibraryModelFile(folderModel, fileName, out)) {
            out.name = name;
            return true;
        }
    } else if (project_) {
        IncludedSource source;
        if (project_(fileName, source)) {
            source.name = name;
            source.fileName = fileName;
            out = std::move(source);
            return true;
        }
    }

    // 2. A library model, by name or ID
    const bool libraryAvailable = library_ && library_->available();
    if (libraryAvailable && readLibraryModelFile(name, "", out)) {
        out.name = name;
        return true;
    }

    std::string where;
    if (!folderModel.empty()) {
        where = "no file " + fileName + " in the folder of library model " + folderModel;
    } else if (project_) {
        where = "no file " + fileName + " in the folder of the model";
    } else {
        where = "the model has no folder (it was not read from a file)";
    }
    error = "Cannot $INCLUDE '" + name + "': " + where;
    if (!libraryAvailable) {
        error += ", and this build does not embed the CoolSolve Library";
        return false;
    }
    error += ", and no library model is named '" + name + "'";
    const auto close = library_->closestNames(name);
    if (!close.empty()) {
        error += " (did you mean ";
        for (size_t i = 0; i < close.size(); ++i) error += (i ? ", " : "") + close[i];
        error += "?)";
    }
    return false;
}

}  // namespace coolsolve
