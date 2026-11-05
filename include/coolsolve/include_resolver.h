#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>

namespace coolsolve {

class LibraryStore;

/** The file an `$INCLUDE` directive resolved to (CS-FEAT-IMPORT). */
struct IncludedSource {
    std::string name;            ///< The name written in the directive ("hx_ees_function", "CSL-0331").
    std::string fileName;        ///< The file read: "hx_ees_function.eescode".
    bool fromLibrary = false;    ///< A file of a model folder of the CoolSolve Library (else: of the project folder).
    std::string libraryId;       ///< "CSL-0331" when `fromLibrary`.
    std::string libraryName;     ///< The library model's name when `fromLibrary`.
    std::string libraryFolder;   ///< "models/heat_exchangers/methods/hx_ees_function" when `fromLibrary`.
    std::string content;         ///< The source text.
    /// Its lookup tables, table name -> CSV text: the `<stem>-<table>.csv` files of its folder (for a
    /// library model, those `LibraryStore::loadModel` loads). Its routines read them as the model's own
    /// tables, which take precedence (addIncludedLookupTables, runner.h).
    std::map<std::string, std::string> tables;

    /// Identity of the file, the same for two includes of it: "library:<folder>/<file>" or "project:<file>".
    std::string key() const;
    /// For messages and reports: "my_props.eescode" or "library model hx_ees_function (CSL-0331)".
    std::string describe() const;
};

/**
 * Where the `$INCLUDE <name>` directives of a model find their files (CS-FEAT-IMPORT).
 *
 * `<name>` is, in this order:
 *  1. the file `<name>.eescode` of the folder of the model (the *project folder*), so that a file of
 *     the project takes precedence over a library model of the same name;
 *  2. the model of the CoolSolve Library embedded in the binary whose name or ID is `<name>`
 *     (`hx_ees_function`, `CSL-0331`): its main `.eescode` file.
 *
 * A library model resolves its own `$INCLUDE`s in its library folder and then in the library, never
 * in the project folder: a library model does not depend on the files of the model that includes it.
 * The resolution is the same in the CLI (the folder of the model file, whatever the working directory)
 * and in the GUI (the folder of the file given on the command line, the files of an uploaded ZIP, or
 * the folder of a library model), so a model includes the same files wherever it is run.
 */
class IncludeResolver {
public:
    /// Reads `fileName` ("x.eescode") and its tables from the project folder into `out.content` and
    /// `out.tables`; false when the folder has no such file.
    using ProjectReader = std::function<bool(const std::string& fileName, IncludedSource& out)>;

    /// Files of a folder on disk (the CLI: the folder of the model file). `mainFile` is the model's own file
    /// name, which cannot be included.
    static std::shared_ptr<const IncludeResolver> forDirectory(const std::string& dir, const std::string& mainFile,
                                                               const LibraryStore* library);
    /// Files held in memory, file name -> content (the GUI session of an uploaded ZIP or of a new model).
    /// The tables of file `<stem>.eescode` are the entries `<stem>-<table>` of `tables` (table name -> CSV),
    /// the names a ZIP upload gives to the files `<stem>-<table>.csv` of another model.
    static std::shared_ptr<const IncludeResolver> forFiles(std::map<std::string, std::string> files,
                                                           std::map<std::string, std::string> tables,
                                                           const std::string& mainFile, const LibraryStore* library);
    /// The folder of a library model (the GUI session of a model opened from the library).
    static std::shared_ptr<const IncludeResolver> forLibraryModel(const std::string& libraryId,
                                                                  const std::string& mainFile,
                                                                  const LibraryStore* library);
    /// No project folder: library models only (a source parsed without file).
    static std::shared_ptr<const IncludeResolver> libraryOnly(const LibraryStore* library);

    /**
     * Resolve `name` for a directive of `from` (null: of the model itself).
     * @return false with `error` set (a sentence naming what was looked for, with the closest library
     *         names) when nothing matches.
     */
    bool resolve(const std::string& name, const IncludedSource* from, IncludedSource& out, std::string& error) const;

private:
    IncludeResolver() = default;

    /// A file of the folder of library model `libraryId` (empty `fileName`: its main file) with its tables.
    bool readLibraryModelFile(const std::string& libraryId, const std::string& fileName, IncludedSource& out) const;

    ProjectReader project_;            ///< Empty: no project folder.
    std::string mainFile_;             ///< The model's own file name (refused as a target).
    std::string rootLibraryId_;        ///< Set when the model itself is a library model.
    const LibraryStore* library_ = nullptr;
};

/**
 * The name of a `$INCLUDE` directive from its content (`hx_ees_function`, `'my_props.eescode'`, with an
 * optional comment): quotes, comments and the `.eescode` extension removed.
 * @return false with `error` for an empty name, a path (`/`, `\`, `..`) or another extension (the
 *         `.lib`, `.txt`, `.ees` files of EES).
 */
bool parseIncludeName(const std::string& content, std::string& name, std::string& error);

}  // namespace coolsolve
