# embed_library.cmake
# Generates C++ source files that embed a snapshot of the CoolSolve Library
# (https://github.com/IntSusEnergySystems/CoolSolve_Library) in the binary.
#
# Usage (run by the `coolsolve_library_snapshot` target at EVERY build, and once at
# configure time -- see CMakeLists.txt):
#   cmake -DLIBRARY_DIR=<library checkout> -DOUTPUT_CPP=<path> -DOUTPUT_H=<path>
#         -DSTAMP_FILE=<path> -P embed_library.cmake
#
# The script hashes the files it would embed and compares the hash and the library
# commit with STAMP_FILE: when nothing changed it exits at once and leaves OUTPUT_CPP
# untouched, so the compiler is not re-run. A library that changed (git pull, local
# edit) is therefore always picked up by `cmake --build`.
#
# What is embedded (allow-list, a mirror of EMBED_FILES in the library's
# tools/build_index.py; never sources/, tools/, templates/ or work/):
#   library.html, library.json, taxonomy.json, docs/taxonomy.md and, in each model
#   folder (a folder holding model.json): README.md, model.json, coolsolve.conf,
#   *.eescode, *.initials, *.sol, <name>-<table>.csv and figures/*.{png,svg,jpg,jpeg,webp}
# Files are served under the URL prefix /library/ with their path in the library.

if(NOT LIBRARY_DIR OR NOT OUTPUT_CPP OR NOT OUTPUT_H OR NOT STAMP_FILE)
    message(FATAL_ERROR "embed_library.cmake requires LIBRARY_DIR, OUTPUT_CPP, OUTPUT_H, STAMP_FILE")
endif()
if(NOT EXISTS "${LIBRARY_DIR}/library.json")
    message(FATAL_ERROR "embed_library.cmake: '${LIBRARY_DIR}' is not a CoolSolve Library (no library.json)")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/embed_common.cmake")

# --- Collect the files to embed (paths relative to the library root) ---
set(LIB_FILES library.html library.json taxonomy.json docs/taxonomy.md)
set(LIB_MODEL_FILES "")   # the model files only: counted against library.json
file(GLOB_RECURSE MODEL_JSONS RELATIVE "${LIBRARY_DIR}" "${LIBRARY_DIR}/models/*/model.json")
list(SORT MODEL_JSONS)
foreach(MJ ${MODEL_JSONS})
    get_filename_component(MDIR "${MJ}" DIRECTORY)
    file(GLOB CANDIDATES RELATIVE "${LIBRARY_DIR}" "${LIBRARY_DIR}/${MDIR}/*" "${LIBRARY_DIR}/${MDIR}/figures/*")
    list(SORT CANDIDATES)
    foreach(C ${CANDIDATES})
        if(IS_DIRECTORY "${LIBRARY_DIR}/${C}")
            continue()
        endif()
        string(LENGTH "${MDIR}/" PREFIX_LEN)
        string(SUBSTRING "${C}" ${PREFIX_LEN} -1 REL)      # path inside the model folder
        string(TOLOWER "${REL}" REL_LC)
        if(REL STREQUAL "README.md" OR REL STREQUAL "model.json" OR REL STREQUAL "coolsolve.conf"
           OR REL_LC MATCHES "^[^/]+[.](eescode|initials|sol)$"
           OR REL_LC MATCHES "^[^/]+-[^/]+[.]csv$"
           OR REL_LC MATCHES "^figures/[^/]+[.](png|svg|jpg|jpeg|webp)$")
            list(APPEND LIB_MODEL_FILES "${C}")
        endif()
    endforeach()
endforeach()
list(APPEND LIB_FILES ${LIB_MODEL_FILES})
list(SORT LIB_FILES)
list(LENGTH LIB_FILES LIB_FILE_COUNT)

# --- Identity: hash of the paths and contents of everything embedded ---
set(HASH_INPUT "")
foreach(F ${LIB_FILES})
    file(SHA256 "${LIBRARY_DIR}/${F}" FILE_SHA)
    string(APPEND HASH_INPUT "${F}:${FILE_SHA}\n")
endforeach()
string(SHA256 SNAPSHOT_HASH "${HASH_INPUT}")
string(SUBSTRING "${SNAPSHOT_HASH}" 0 16 SNAPSHOT_HASH)

# --- Library commit (informative), "+local-changes" when embedded files are modified ---
set(LIBRARY_COMMIT "unknown")
find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${LIBRARY_DIR}" rev-parse --short HEAD
        OUTPUT_VARIABLE COMMIT_OUT RESULT_VARIABLE COMMIT_RES ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(COMMIT_RES EQUAL 0 AND COMMIT_OUT)
        set(LIBRARY_COMMIT "${COMMIT_OUT}")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${LIBRARY_DIR}" status --porcelain
                -- library.html library.json taxonomy.json docs/taxonomy.md models
            OUTPUT_VARIABLE DIRTY_OUT ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(DIRTY_OUT)
            set(LIBRARY_COMMIT "${LIBRARY_COMMIT}+local-changes")
        endif()
    endif()
endif()

# The commit is part of the stamp: committing the files embedded so far changes the
# identity printed by `coolsolve --version`, not only a change of content.
set(STAMP "${SNAPSHOT_HASH} ${LIBRARY_COMMIT}")
if(EXISTS "${STAMP_FILE}" AND EXISTS "${OUTPUT_CPP}" AND EXISTS "${OUTPUT_H}")
    file(READ "${STAMP_FILE}" OLD_STAMP)
    if("${OLD_STAMP}" STREQUAL "${STAMP}")
        return()
    endif()
endif()

# --- Consistency with the library's own index (warn only) ---
file(READ "${LIBRARY_DIR}/library.json" LIBRARY_JSON_HEAD LIMIT 600)
string(REGEX MATCH "\"files\": *([0-9]+)" _ "${LIBRARY_JSON_HEAD}")
set(INDEX_FILES "${CMAKE_MATCH_1}")
list(LENGTH LIB_MODEL_FILES MODEL_FILE_COUNT)
math(EXPR MODEL_FILE_COUNT_WITH_TAXONOMY "${MODEL_FILE_COUNT} + 2")   # + taxonomy.json, docs/taxonomy.md
if(NOT "${INDEX_FILES}" STREQUAL "${MODEL_FILE_COUNT_WITH_TAXONOMY}")
    message(WARNING "CoolSolve Library: ${MODEL_FILE_COUNT_WITH_TAXONOMY} files match the embedding allow-list but "
            "library.json lists ${INDEX_FILES}. Either library.json is stale (run `python3 tools/build_index.py` "
            "in the library) or the allow-list of cmake/embed_library.cmake differs from EMBED_FILES in "
            "tools/build_index.py.")
endif()

# --- Generate ---
set(HEADER_CONTENT [=[
// Auto-generated by embed_library.cmake -- do not edit
#pragma once
#include <string>
#include <vector>

/// One file of the embedded CoolSolve Library snapshot.
struct EmbeddedLibraryFile {
    const unsigned char* data;   // null-terminated (text files can be used as C strings)
    unsigned int size;           // without the terminator
    const char* mimeType;
};

/// File by path relative to the library root ("models/cycles/.../README.md"); nullptr if absent.
const EmbeddedLibraryFile* getEmbeddedLibraryFile(const std::string& path);

/// Paths (sorted) of the embedded files that start with `prefix` (e.g. "models/cycles/x/").
std::vector<std::string> listEmbeddedLibraryFiles(const std::string& prefix);

/// Git revision of the library checkout the snapshot was taken from ("unknown" without git).
const char* embeddedLibraryCommit();

/// Content hash of the embedded files (same files, same hash).
const char* embeddedLibraryHash();

]=])

set(SOURCE_CONTENT [=[
// Auto-generated by embed_library.cmake -- do not edit
#include "embedded_library.h"
#include <map>

]=])

set(MAP_ENTRIES "")
set(INDEX 0)
foreach(F ${LIB_FILES})
    coolsolve_file_as_hex_array("${LIBRARY_DIR}/${F}" FILE_HEX FILE_SIZE)
    guess_mime_type("${F}" MIME)
    if(F MATCHES "[.](md|eescode|initials|sol|conf|csv)$")
        set(MIME "text/plain; charset=utf-8")
    endif()
    string(APPEND SOURCE_CONTENT "static const unsigned char libfile_${INDEX}[] = {\n    ${FILE_HEX}\n};\n")
    string(APPEND MAP_ENTRIES "    {\"${F}\", {libfile_${INDEX}, ${FILE_SIZE}, \"${MIME}\"}},\n")
    math(EXPR INDEX "${INDEX} + 1")
endforeach()

string(APPEND SOURCE_CONTENT "\nstatic const std::map<std::string, EmbeddedLibraryFile> libraryFiles = {\n${MAP_ENTRIES}};\n")
string(APPEND SOURCE_CONTENT "
const EmbeddedLibraryFile* getEmbeddedLibraryFile(const std::string& path) {
    auto it = libraryFiles.find(path);
    return it == libraryFiles.end() ? nullptr : &it->second;
}

std::vector<std::string> listEmbeddedLibraryFiles(const std::string& prefix) {
    std::vector<std::string> out;
    for (auto it = libraryFiles.lower_bound(prefix); it != libraryFiles.end(); ++it) {
        if (it->first.compare(0, prefix.size(), prefix) != 0) break;
        out.push_back(it->first);
    }
    return out;
}

const char* embeddedLibraryCommit() { return \"${LIBRARY_COMMIT}\"; }
const char* embeddedLibraryHash() { return \"${SNAPSHOT_HASH}\"; }
")

file(WRITE "${OUTPUT_CPP}" "${SOURCE_CONTENT}")
# The header does not depend on the library: rewriting it would recompile its includers
# (server.cpp, version.cpp) at every library change.
set(OLD_HEADER "")
if(EXISTS "${OUTPUT_H}")
    file(READ "${OUTPUT_H}" OLD_HEADER)
endif()
if(NOT "${OLD_HEADER}" STREQUAL "${HEADER_CONTENT}")
    file(WRITE "${OUTPUT_H}" "${HEADER_CONTENT}")
endif()
file(WRITE "${STAMP_FILE}" "${STAMP}")

message(STATUS "Embedded ${LIB_FILE_COUNT} CoolSolve Library files (commit ${LIBRARY_COMMIT}, hash ${SNAPSHOT_HASH}) into ${OUTPUT_CPP}")
