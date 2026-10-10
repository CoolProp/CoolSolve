#include "coolsolve/ees_string.h"

namespace coolsolve {

std::string eesFoldCase(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) out += static_cast<char>(eesFoldByteAt(s, i));
    return out;
}

bool eesStringEqual(const std::string& a, const std::string& b) {
    return a.size() == b.size() && eesCompareFolded(a, b) == 0;
}

}  // namespace coolsolve
