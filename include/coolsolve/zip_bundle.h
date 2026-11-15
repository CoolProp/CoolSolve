#pragma once

#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace coolsolve {

/**
 * @file zip_bundle.h
 * @brief The ZIP bundle of the GUI: Save writes one, Open reads any ZIP a zip tool makes.
 *
 * The archives are read with miniz (built into CoolProp): stored and Deflate entries (the
 * method of every zip tool by default), data descriptors, ZIP64. Names are UTF-8, or IBM
 * code page 437 when they are not valid UTF-8 (the encoding of the ZIP format when the
 * UTF-8 flag is not set, Windows' *Compressed folder*). See docs/gui.md §7.2.
 */

/// A ZIP archive (or a file of it) that cannot be read; `what()` is the message for the user.
class ZipBundleError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// The largest uncompressed size of an uploaded bundle (all its files): 256 MB. A model bundle
/// weighs kilobytes; the limit stops an archive that expands to gigabytes.
constexpr std::size_t kZipBundleMaxBytes = std::size_t(256) << 20;

/**
 * @brief A ZIP archive of `files` (name in the archive, content), stored without compression.
 *
 * The bundle downloaded by *Save* (`GET /api/v1/files/bundle`).
 */
std::string createZipBundle(const std::vector<std::pair<std::string, std::string>>& files);

/**
 * @brief The files of the ZIP archive `data`: name in the archive -> content.
 *
 * Folders, the `__MACOSX/` folder and the `._*` files of macOS are left out. When every file
 * is in one top folder (the ZIP of a folder), that folder is removed from the names, except
 * `debug_output/` (the debug files of a bundle). `\` in names becomes `/`.
 *
 * @param maxBytes the largest uncompressed size of all files together.
 * @throws ZipBundleError when `data` is not a ZIP archive or is damaged, when a file is
 *         encrypted, compressed with a method other than Deflate, has an absolute name or a
 *         `..` in its name, or when the files exceed `maxBytes`.
 */
std::map<std::string, std::string> extractZipBundle(const std::string& data,
                                                    std::size_t maxBytes = kZipBundleMaxBytes);

/**
 * @brief Check that the file `name` of a bundle is UTF-8 text (the GUI exchanges text as JSON).
 *
 * @throws ZipBundleError naming the first invalid byte and its line, e.g. a model saved by EES
 *         (Windows ANSI) with a `°` in it.
 */
void checkUtf8Text(const std::string& name, const std::string& content);

}  // namespace coolsolve
