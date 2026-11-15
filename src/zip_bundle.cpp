#include "coolsolve/zip_bundle.h"

#include <cstdint>
#include <cstdio>
#include <cstring>

// miniz is compiled into CoolProp (externals/miniz-*); its zlib names (crc32, inflate, ...) are not wanted here.
#define MINIZ_NO_ZLIB_COMPATIBLE_NAMES
#include "miniz.h"

namespace coolsolve {

namespace {

void zipW16(std::string& s, uint16_t v) {
    s += static_cast<char>(v & 0xFF);
    s += static_cast<char>((v >> 8) & 0xFF);
}

void zipW32(std::string& s, uint32_t v) {
    s += static_cast<char>(v & 0xFF);
    s += static_cast<char>((v >> 8) & 0xFF);
    s += static_cast<char>((v >> 16) & 0xFF);
    s += static_cast<char>((v >> 24) & 0xFF);
}

/// Offset of the first byte of `s` that does not start a valid UTF-8 character (RFC 3629: no
/// overlong form, no surrogate, nothing above U+10FFFF); npos when `s` is all valid.
std::size_t invalidUtf8Offset(const std::string& s) {
    std::size_t i = 0;
    while (i < s.size()) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c < 0x80) { ++i; continue; }
        std::size_t len = 0;
        unsigned char lo = 0x80, hi = 0xBF;  // range of the second byte
        if (c >= 0xC2 && c <= 0xDF) len = 2;
        else if (c >= 0xE0 && c <= 0xEF) { len = 3; if (c == 0xE0) lo = 0xA0; if (c == 0xED) hi = 0x9F; }
        else if (c >= 0xF0 && c <= 0xF4) { len = 4; if (c == 0xF0) lo = 0x90; if (c == 0xF4) hi = 0x8F; }
        if (len == 0 || i + len > s.size()) return i;
        for (std::size_t k = 1; k < len; ++k) {
            const unsigned char cc = static_cast<unsigned char>(s[i + k]);
            if (k == 1 ? (cc < lo || cc > hi) : (cc & 0xC0) != 0x80) return i;
        }
        i += len;
    }
    return std::string::npos;
}

void appendUtf8(std::string& out, uint32_t cp) {
    if (cp < 0x80) {
        out += static_cast<char>(cp);
    } else if (cp < 0x800) {
        out += static_cast<char>(0xC0 | (cp >> 6));
        out += static_cast<char>(0x80 | (cp & 0x3F));
    } else {
        out += static_cast<char>(0xE0 | (cp >> 12));
        out += static_cast<char>(0x80 | ((cp >> 6) & 0x3F));
        out += static_cast<char>(0x80 | (cp & 0x3F));
    }
}

/// A file name of IBM code page 437 in UTF-8: the encoding of ZIP names when the UTF-8 flag is not
/// set (APPNOTE.TXT, appendix D), which Windows' Compressed folder writes for accented names.
std::string cp437ToUtf8(const std::string& s) {
    static const uint16_t upper[128] = {
        0x00C7, 0x00FC, 0x00E9, 0x00E2, 0x00E4, 0x00E0, 0x00E5, 0x00E7,
        0x00EA, 0x00EB, 0x00E8, 0x00EF, 0x00EE, 0x00EC, 0x00C4, 0x00C5,
        0x00C9, 0x00E6, 0x00C6, 0x00F4, 0x00F6, 0x00F2, 0x00FB, 0x00F9,
        0x00FF, 0x00D6, 0x00DC, 0x00A2, 0x00A3, 0x00A5, 0x20A7, 0x0192,
        0x00E1, 0x00ED, 0x00F3, 0x00FA, 0x00F1, 0x00D1, 0x00AA, 0x00BA,
        0x00BF, 0x2310, 0x00AC, 0x00BD, 0x00BC, 0x00A1, 0x00AB, 0x00BB,
        0x2591, 0x2592, 0x2593, 0x2502, 0x2524, 0x2561, 0x2562, 0x2556,
        0x2555, 0x2563, 0x2551, 0x2557, 0x255D, 0x255C, 0x255B, 0x2510,
        0x2514, 0x2534, 0x252C, 0x251C, 0x2500, 0x253C, 0x255E, 0x255F,
        0x255A, 0x2554, 0x2569, 0x2566, 0x2560, 0x2550, 0x256C, 0x2567,
        0x2568, 0x2564, 0x2565, 0x2559, 0x2558, 0x2552, 0x2553, 0x256B,
        0x256A, 0x2518, 0x250C, 0x2588, 0x2584, 0x258C, 0x2590, 0x2580,
        0x03B1, 0x00DF, 0x0393, 0x03C0, 0x03A3, 0x03C3, 0x00B5, 0x03C4,
        0x03A6, 0x0398, 0x03A9, 0x03B4, 0x221E, 0x03C6, 0x03B5, 0x2229,
        0x2261, 0x00B1, 0x2265, 0x2264, 0x2320, 0x2321, 0x00F7, 0x2248,
        0x00B0, 0x2219, 0x00B7, 0x221A, 0x207F, 0x00B2, 0x25A0, 0x00A0,
    };
    std::string out;
    for (unsigned char c : s) appendUtf8(out, c < 0x80 ? c : upper[c - 0x80]);
    return out;
}

/// The name of the compression method `method` of a ZIP entry, for the message that refuses it.
std::string methodName(unsigned method) {
    switch (method) {
        case 9: return "Deflate64";
        case 12: return "BZIP2";
        case 14: return "LZMA";
        case 93: return "Zstandard";
        case 95: return "XZ";
        case 98: return "PPMd";
        default: return "method " + std::to_string(method);
    }
}

/// What `data` is when it is not a ZIP archive but another archive format, else empty.
std::string otherArchiveFormat(const std::string& data) {
    struct Magic { const char* bytes; std::size_t size; const char* name; };
    static const Magic magics[] = {
        {"Rar!\x1A\x07", 6, "a RAR archive"},
        {"7z\xBC\xAF\x27\x1C", 6, "a 7z archive"},
        {"\x1F\x8B", 2, "a gzip file"},
        {"BZh", 3, "a bzip2 file"},
        {"\xFD" "7zXZ", 5, "an XZ file"},
    };
    for (const auto& m : magics)
        if (data.size() >= m.size && std::memcmp(data.data(), m.bytes, m.size) == 0) return m.name;
    return "";
}

/// Ends the miniz reader of an archive however extractZipBundle leaves.
struct ZipReader {
    mz_zip_archive zip;
    bool open = false;
    ZipReader() { mz_zip_zero_struct(&zip); }
    ~ZipReader() { if (open) mz_zip_reader_end(&zip); }
    std::string lastError() { return mz_zip_get_error_string(mz_zip_get_last_error(&zip)); }
};

/// macOS metadata, not files of the model: the `__MACOSX/` folder, `._*` files, `.DS_Store`.
bool isMacMetadata(const std::string& name) {
    if (name.compare(0, 9, "__MACOSX/") == 0) return true;
    const auto slash = name.rfind('/');
    const std::string base = slash == std::string::npos ? name : name.substr(slash + 1);
    return base.compare(0, 2, "._") == 0 || base == ".DS_Store";
}

/// Whether `name` would be written outside the folder it is extracted to: absolute, or with `..`.
bool isUnsafePath(const std::string& name) {
    if (!name.empty() && name[0] == '/') return true;
    if (name.size() >= 2 && name[1] == ':') return true;  // C:...
    std::size_t start = 0;
    while (start <= name.size()) {
        auto end = name.find('/', start);
        if (end == std::string::npos) end = name.size();
        if (name.compare(start, end - start, "..") == 0 && end - start == 2) return true;
        start = end + 1;
    }
    return false;
}

}  // namespace

std::string createZipBundle(const std::vector<std::pair<std::string, std::string>>& files) {
    std::string out;
    struct E { std::string name; uint32_t crc, size, off; };
    std::vector<E> entries;
    for (const auto& [name, data] : files) {
        E e{name,
            static_cast<uint32_t>(mz_crc32(MZ_CRC32_INIT, reinterpret_cast<const unsigned char*>(data.data()), data.size())),
            static_cast<uint32_t>(data.size()), static_cast<uint32_t>(out.size())};
        zipW32(out, 0x04034b50); zipW16(out, 20); zipW16(out, 0); zipW16(out, 0);
        zipW16(out, 0); zipW16(out, 0);
        zipW32(out, e.crc); zipW32(out, e.size); zipW32(out, e.size);
        zipW16(out, static_cast<uint16_t>(name.size())); zipW16(out, 0);
        out += name; out += data;
        entries.push_back(e);
    }
    uint32_t cdOff = static_cast<uint32_t>(out.size());
    for (const auto& e : entries) {
        zipW32(out, 0x02014b50); zipW16(out, 20); zipW16(out, 20);
        zipW16(out, 0); zipW16(out, 0); zipW16(out, 0); zipW16(out, 0);
        zipW32(out, e.crc); zipW32(out, e.size); zipW32(out, e.size);
        zipW16(out, static_cast<uint16_t>(e.name.size()));
        zipW16(out, 0); zipW16(out, 0); zipW16(out, 0); zipW16(out, 0);
        zipW32(out, 0); zipW32(out, e.off);
        out += e.name;
    }
    uint32_t cdSize = static_cast<uint32_t>(out.size()) - cdOff;
    zipW32(out, 0x06054b50); zipW16(out, 0); zipW16(out, 0);
    zipW16(out, static_cast<uint16_t>(entries.size()));
    zipW16(out, static_cast<uint16_t>(entries.size()));
    zipW32(out, cdSize); zipW32(out, cdOff); zipW16(out, 0);
    return out;
}

std::map<std::string, std::string> extractZipBundle(const std::string& data, std::size_t maxBytes) {
    ZipReader reader;
    if (!mz_zip_reader_init_mem(&reader.zip, data.data(), data.size(), 0)) {
        const std::string other = otherArchiveFormat(data);
        if (!other.empty())
            throw ZipBundleError("the file is " + other + ", not a ZIP archive: save the model folder as a ZIP file");
        throw ZipBundleError("the file is not a ZIP archive, or it is damaged (" + reader.lastError() + ")");
    }
    reader.open = true;

    struct Entry { mz_uint index; std::string name; std::size_t size; };
    std::vector<Entry> entries;
    std::size_t total = 0;
    const mz_uint count = mz_zip_reader_get_num_files(&reader.zip);
    for (mz_uint i = 0; i < count; ++i) {
        mz_zip_archive_file_stat st;
        if (!mz_zip_reader_file_stat(&reader.zip, i, &st))
            throw ZipBundleError("the ZIP archive is damaged (" + reader.lastError() + ")");
        std::string name = st.m_filename;
        if (invalidUtf8Offset(name) != std::string::npos) name = cp437ToUtf8(name);
        for (char& c : name) if (c == '\\') c = '/';
        while (name.compare(0, 2, "./") == 0) name.erase(0, 2);
        if (st.m_is_directory || name.empty() || name.back() == '/' || isMacMetadata(name)) continue;

        if (isUnsafePath(name))
            throw ZipBundleError("the ZIP archive names a file outside its folder ('" + name + "')");
        if (st.m_is_encrypted)
            throw ZipBundleError("'" + name + "' is password-protected: CoolSolve reads ZIP files without a password");
        if (!st.m_is_supported)
            throw ZipBundleError("'" + name + "' is compressed with " + methodName(st.m_method) +
                                 ", which CoolSolve cannot read: zip it again with the Deflate method (the default of "
                                 "zip tools, of Windows' Compressed folder and of macOS' Compress)");
        if (st.m_uncomp_size > maxBytes || total + st.m_uncomp_size > maxBytes)
            throw ZipBundleError("the files of the ZIP archive exceed " + std::to_string(maxBytes >> 20) +
                                 " MB once uncompressed");
        total += static_cast<std::size_t>(st.m_uncomp_size);
        entries.push_back({i, name, static_cast<std::size_t>(st.m_uncomp_size)});
    }

    // The ZIP of a folder: every file in one top folder, which is not part of the model's names
    std::string top;
    for (const auto& e : entries) {
        const auto slash = e.name.find('/');
        const std::string first = slash == std::string::npos ? "" : e.name.substr(0, slash + 1);
        if (first.empty() || (!top.empty() && first != top)) { top.clear(); break; }
        top = first;
    }
    if (top == "debug_output/") top.clear();

    std::map<std::string, std::string> files;
    for (const auto& e : entries) {
        std::string content(e.size, '\0');
        if (!mz_zip_reader_extract_to_mem(&reader.zip, e.index, &content[0], content.size(), 0))
            throw ZipBundleError("'" + e.name + "' cannot be extracted, the ZIP archive is damaged (" +
                                 reader.lastError() + ")");
        files[e.name.substr(top.size())] = std::move(content);
    }
    return files;
}

void checkUtf8Text(const std::string& name, const std::string& content) {
    const std::size_t at = invalidUtf8Offset(content);
    if (at == std::string::npos) return;
    std::size_t line = 1;
    for (std::size_t i = 0; i < at; ++i) if (content[i] == '\n') ++line;
    char byte[8];
    std::snprintf(byte, sizeof byte, "0x%02X", static_cast<unsigned char>(content[at]));
    throw ZipBundleError("'" + name + "' is not UTF-8 text (byte " + byte + " on line " + std::to_string(line) +
                         "): save it with the UTF-8 encoding (EES writes Windows ANSI) and zip it again");
}

}  // namespace coolsolve
