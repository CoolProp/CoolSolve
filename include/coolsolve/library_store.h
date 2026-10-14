#pragma once

#include <functional>
#include <memory>
#include <map>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace coolsolve {

/**
 * Read access to the files of a CoolSolve Library snapshot.
 *
 * Paths are relative to the library root ("library.json", "models/cycles/x/README.md").
 * The server feeds the store with the files embedded in the binary; tests feed it an
 * in-memory map, so the store itself knows nothing about embedding or HTTP.
 */
struct LibraryFileSource {
    /// Content of a file; false if it does not exist.
    std::function<bool(const std::string& path, std::string& content)> read;
    /// Paths (sorted) of the files whose path starts with `prefix`.
    std::function<std::vector<std::string>(const std::string& prefix)> list;

    /// Source over an in-memory map path -> content.
    static LibraryFileSource fromMap(std::map<std::string, std::string> files);
};

/** The files of a library model folder, ready to be loaded into a GUI session. */
struct LibraryModel {
    std::string id;                          ///< permanent ID, "CSL-0012"
    std::string name;                        ///< unique name, also the file stem of the main file
    std::string title;
    std::string status;                      ///< verified, runs, blocked, ...
    std::string folder;                      ///< "models/cycles/organic_rankine/orc_simple_r245fa"
    std::string file;                        ///< the .eescode file loaded (main file or a variant)
    std::string eescode;
    std::string initials;                    ///< empty if the model has none
    std::string sol;                         ///< CoolSolve solution (baseline), empty if none
    std::string conf;                        ///< coolsolve.conf, empty if none
    std::map<std::string, std::string> tables;  ///< lookup table name -> CSV content
    std::vector<std::string> variants;       ///< the other .eescode files of the folder
};

/**
 * Index of the CoolSolve Library snapshot (library.json + taxonomy.json) and loader of
 * model folders. A model is a whole folder: besides the `.eescode` file it carries
 * `.initials`, `coolsolve.conf`, `<name>-<table>.csv` lookup tables, the `.sol`
 * baseline, and a `README.md` shown by the documentation viewer.
 */
class LibraryStore {
public:
    explicit LibraryStore(LibraryFileSource source);

    /// False when the snapshot is missing or its library.json cannot be parsed.
    bool available() const { return available_; }
    /// Why the library is not available (empty when it is).
    const std::string& unavailableReason() const { return reason_; }

    /// `snapshot` block of library.json (content hash, counts, taxonomy version).
    nlohmann::json snapshot() const;

    /**
     * The model menu of the GUI: the taxonomy (categories, sub-categories) holding the
     * models whose status is `verified`, empty categories left out. Blocked, failing or
     * documented-only models stay out of the menu; the dashboard (/library/) lists them all.
     *
     *     {"categories":[{"key","title","count","models":[...],"children":[{"key","title","models":[...]}]}],
     *      "total":N}
     * with models as {"id","name","title","level","kind"}.
     */
    nlohmann::json menu() const;

    /**
     * Load a model folder by ID ("CSL-0012") or by name. `file` selects one of the folder's
     * variants (e.g. "x_coolsolve.eescode"); empty = the main file. Any status is allowed:
     * a blocked model opens too (its native file may then fail to solve).
     * On failure returns false and sets `error`.
     */
    bool loadModel(const std::string& idOrName, const std::string& file, LibraryModel& out,
                   std::string& error) const;

private:
    LibraryFileSource source_;
    bool available_ = false;
    std::string reason_;
    nlohmann::json index_;     // library.json
    nlohmann::json taxonomy_;  // taxonomy.json
    bool readOptional(const std::string& path, std::string& content) const;
};

}  // namespace coolsolve
