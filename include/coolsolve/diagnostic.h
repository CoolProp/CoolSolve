#pragma once

#include <string>
#include <vector>
#include <mutex>

namespace coolsolve {

// ============================================================================
// Diagnostic Severity & Structure
// ============================================================================

enum class DiagnosticSeverity { Error, Warning, Info };

/**
 * @brief A single diagnostic message (error, warning, or info).
 *
 * Diagnostics are produced by the parser, evaluator, solver, and solution
 * checker.  They carry a machine-readable code (e.g. "P003") and a
 * human-readable message, plus optional source location.
 */
struct Diagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::Info;
    std::string code;       // Machine-readable ID, e.g. "P003", "C001"
    std::string message;    // Human-readable text
    int line   = 0;         // Source line (0 = not applicable)
    int column = 0;         // Source column (0 = not applicable)
    std::string source;     // Subsystem: "parser", "evaluator", "solver",
                            //            "coolprop", "checker"
};

// ============================================================================
// Diagnostic Collector
// ============================================================================

/**
 * @brief Lightweight, thread-safe collector for diagnostics.
 *
 * Default-constructed with no heap allocation; the internal vector is only
 * allocated on the first push.  Movable but not copyable.
 */
class DiagnosticCollector {
public:
    DiagnosticCollector() = default;
    ~DiagnosticCollector() = default;

    DiagnosticCollector(DiagnosticCollector&& other) noexcept {
        std::lock_guard<std::mutex> lock(other.mu_);
        items_ = std::move(other.items_);
    }
    DiagnosticCollector& operator=(DiagnosticCollector&& other) noexcept {
        if (this != &other) {
            std::lock_guard<std::mutex> lockOther(other.mu_);
            std::lock_guard<std::mutex> lockThis(mu_);
            items_ = std::move(other.items_);
        }
        return *this;
    }

    DiagnosticCollector(const DiagnosticCollector& other) {
        std::lock_guard<std::mutex> lock(other.mu_);
        items_ = other.items_;
    }
    DiagnosticCollector& operator=(const DiagnosticCollector& other) {
        if (this != &other) {
            std::lock_guard<std::mutex> lockOther(other.mu_);
            std::lock_guard<std::mutex> lockThis(mu_);
            items_ = other.items_;
        }
        return *this;
    }

    void push(Diagnostic d) {
        std::lock_guard<std::mutex> lock(mu_);
        items_.push_back(std::move(d));
    }

    void push(DiagnosticSeverity sev, const std::string& code,
              const std::string& msg, const std::string& src,
              int line = 0, int col = 0) {
        Diagnostic d;
        d.severity = sev;
        d.code     = code;
        d.message  = msg;
        d.source   = src;
        d.line     = line;
        d.column   = col;
        push(std::move(d));
    }

    bool hasErrors() const {
        std::lock_guard<std::mutex> lock(mu_);
        for (const auto& d : items_)
            if (d.severity == DiagnosticSeverity::Error) return true;
        return false;
    }

    bool hasWarnings() const {
        std::lock_guard<std::mutex> lock(mu_);
        for (const auto& d : items_)
            if (d.severity == DiagnosticSeverity::Warning) return true;
        return false;
    }

    bool empty() const {
        std::lock_guard<std::mutex> lock(mu_);
        return items_.empty();
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(mu_);
        return items_.size();
    }

    const std::vector<Diagnostic>& items() const { return items_; }

    /// Merge another collector (move semantics).
    void merge(DiagnosticCollector&& other) {
        std::lock_guard<std::mutex> lockOther(other.mu_);
        std::lock_guard<std::mutex> lockThis(mu_);
        items_.insert(items_.end(),
                      std::make_move_iterator(other.items_.begin()),
                      std::make_move_iterator(other.items_.end()));
        other.items_.clear();
    }

    /// Merge another collector (copy semantics).
    void merge(const DiagnosticCollector& other) {
        std::lock_guard<std::mutex> lockOther(other.mu_);
        std::lock_guard<std::mutex> lockThis(mu_);
        items_.insert(items_.end(), other.items_.begin(), other.items_.end());
    }

    /// Clear all diagnostics.
    void clear() {
        std::lock_guard<std::mutex> lock(mu_);
        items_.clear();
    }

private:
    mutable std::mutex mu_;
    std::vector<Diagnostic> items_;
};

// ============================================================================
// Severity helpers
// ============================================================================

inline const char* severityToString(DiagnosticSeverity sev) {
    switch (sev) {
        case DiagnosticSeverity::Error:   return "error";
        case DiagnosticSeverity::Warning: return "warning";
        case DiagnosticSeverity::Info:    return "info";
    }
    return "info";
}

}  // namespace coolsolve
