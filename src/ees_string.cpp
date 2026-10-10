#include "coolsolve/ees_string.h"

namespace coolsolve {

namespace {

/// Lower-case form used for the comparison: ASCII letters, and the capitals of the
/// Latin-1 supplement (UTF-8: C3 80..C3 9E except the multiplication sign C3 97,
/// lower case = second byte + 0x20). Locale independent.
std::string foldCase(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        if (c >= 'A' && c <= 'Z') {
            out += static_cast<char>(c + ('a' - 'A'));
        } else if (c == 0xC3 && i + 1 < s.size()) {
            unsigned char d = static_cast<unsigned char>(s[++i]);
            if (d >= 0x80 && d <= 0x9E && d != 0x97) d = static_cast<unsigned char>(d + 0x20);
            out += static_cast<char>(c);
            out += static_cast<char>(d);
        } else {
            out += static_cast<char>(c);
        }
    }
    return out;
}

}  // namespace

bool eesStringEqual(const std::string& a, const std::string& b) {
    return a.size() == b.size() && foldCase(a) == foldCase(b);
}

}  // namespace coolsolve
