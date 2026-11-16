#include "coolsolve/chem_equil.h"
#include "coolsolve/nasa_glenn.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace coolsolve::chem_equil {

namespace {

// The method is documented in chem_equil.h. The notation: species j, elements c, a_cj the atoms of element c in species j,
// b_c the abundance of element c (b_O = 1, the ratios otherwise), g_j = mu0_j/RT + ln(P/P0), the element potentials pi_c and
// nu = ln N; the moles of the species n_j = exp(nu - g_j + sum_c a_cj pi_c).

/// The elements; the index of an element + 1 is that of its ratio in `Inputs` (O has no input, its abundance is 1).
enum Element { kO = 0, kAr = 1, kC = 2, kH = 3, kN = 4 };
constexpr int kElements = 5;
constexpr int kMaxSize = kElements + 1;   ///< the potentials and nu, for the derivatives

constexpr double kInnerTolerance = 1e-13;   ///< on the element balances, relative to the abundance of each element
constexpr double kOuterTolerance = 1e-13;   ///< on ln(sum x): not tighter than the balances, which carry sum n
constexpr double kMaxPotentialStep = 5.0;   ///< the largest change of a potential in one Newton step
constexpr double kMaxNuStep = 2.0;          ///< the largest change of ln N before the root is bracketed
constexpr double kRidge = 1e-13;            ///< relative regularisation of the scaled Hessian (its null direction, below)
constexpr int kMaxInnerSteps = 500;
constexpr int kMaxOuterSteps = 100;
constexpr int kMaxBacktracks = 40;

std::string numberText(double v) {
    std::ostringstream oss;
    oss << std::setprecision(10) << v;
    return oss.str();
}

/// The NASA-Glenn record of each species and the atoms of each element in it, found once.
struct SpeciesTable {
    const nasa_glenn::Species* data[kSpeciesCount];
    int atoms[kSpeciesCount][kElements];
};

const SpeciesTable& speciesTable() {
    static const SpeciesTable table = [] {
        SpeciesTable t{};
        for (int j = 0; j < kSpeciesCount; ++j) {
            t.data[j] = nasa_glenn::find(kSpeciesNames[j]);
            if (!t.data[j]) throw std::logic_error(std::string("CHEM_EQUIL: no NASA-Glenn data for ") + kSpeciesNames[j]);
            const nasa_glenn::Composition& c = t.data[j]->composition;
            t.atoms[j][kO] = c.O;
            t.atoms[j][kAr] = c.Ar;
            t.atoms[j][kC] = c.C;
            t.atoms[j][kH] = c.H;
            t.atoms[j][kN] = c.N;
        }
        return t;
    }();
    return table;
}

/// The abundance of each element, from the inputs (O: 1).
void abundances(const Inputs& in, double b[kElements]) {
    b[kO] = 1.0;
    for (int e = 1; e < kElements; ++e) b[e] = in[e + 1];
}

/// The present species: those whose elements all have an abundance above 0.
void presentSpecies(const double b[kElements], bool present[kSpeciesCount]) {
    const SpeciesTable& t = speciesTable();
    for (int j = 0; j < kSpeciesCount; ++j) {
        present[j] = true;
        for (int e = 0; e < kElements; ++e) present[j] = present[j] && (t.atoms[j][e] == 0 || b[e] > 0.0);
    }
}

TemperatureRange rangeOf(const bool present[kSpeciesCount]) {
    const SpeciesTable& t = speciesTable();
    TemperatureRange r{0.0, INFINITY};
    for (int j = 0; j < kSpeciesCount; ++j) {
        if (!present[j]) continue;
        r.tMin = std::max(r.tMin, t.data[j]->tMin());
        r.tMax = std::min(r.tMax, t.data[j]->tMax());
    }
    return r;
}

/// The problem at one state: the elements present (numbered 0 ... ne-1 for the potentials) and the species present.
struct Problem {
    double P = 0.0, T = 0.0;
    int ne = 0;                    ///< the number of elements present
    int element[kElements] = {};   ///< the Element of each potential
    double b[kElements] = {};      ///< the abundance of each element present
    int ns = 0;                    ///< the number of species present
    int species[kSpeciesCount] = {};   ///< their index in the outputs
    double a[kElements][kSpeciesCount] = {};   ///< atoms of element c in present species k
    double g[kSpeciesCount] = {};  ///< mu0/RT + ln(P/P0) of present species k
    double u[kSpeciesCount] = {};  ///< h/(R T^2) = -dg/dT of present species k [1/K]
};

void setup(const Inputs& in, Problem& p) {
    if (!(in[0] > 0.0) || !std::isfinite(in[0]))
        throw std::domain_error("the pressure P = " + numberText(in[0]) + " Pa is not positive");
    if (!std::isfinite(in[1])) throw std::domain_error("the temperature T = " + numberText(in[1]) + " K is not finite");
    for (int i = 2; i < kInputCount; ++i)
        if (!(in[i] >= 0.0) || !std::isfinite(in[i]))
            throw std::domain_error(std::string("the ratio ") + kInputNames[i] + " = " + numberText(in[i]) +
                                    " is negative or not finite (the atoms of the element per oxygen atom)");
    double b[kElements];
    abundances(in, b);
    // the carbon must fit in CO, CO2 and CH4: at most one per oxygen atom and one per four hydrogen atoms
    if (b[kC] >= b[kO] + 0.25 * b[kH])
        throw std::domain_error("C\\O = " + numberText(b[kC]) + " is not below 1 + H\\O/4 = " + numberText(1.0 + 0.25 * b[kH]) +
                                ": the carbon cannot be held by CO, CO2 and CH4, the mixture would form graphite, which the "
                                "14 gaseous species of the equilibrium do not include");
    bool present[kSpeciesCount];
    presentSpecies(b, present);
    const TemperatureRange range = rangeOf(present);
    if (!(in[1] >= range.tMin && in[1] <= range.tMax))
        throw std::out_of_range("the temperature T = " + numberText(in[1]) + " K is outside the range of the NASA-Glenn data of "
                                "the species of the mixture, " + numberText(range.tMin) + " K to " + numberText(range.tMax) + " K");
    p.P = in[0];
    p.T = in[1];
    p.ne = 0;
    for (int e = 0; e < kElements; ++e)
        if (b[e] > 0.0) {
            p.element[p.ne] = e;
            p.b[p.ne] = b[e];
            ++p.ne;
        }
    const SpeciesTable& t = speciesTable();
    const double R = nasa_glenn::kGasConstant;
    const double lnP = std::log(p.P / nasa_glenn::kStandardPressure);
    p.ns = 0;
    for (int j = 0; j < kSpeciesCount; ++j) {
        if (!present[j]) continue;
        const int k = p.ns++;
        p.species[k] = j;
        for (int c = 0; c < p.ne; ++c) p.a[c][k] = t.atoms[j][p.element[c]];
        const nasa_glenn::Point q = nasa_glenn::evaluate(*t.data[j], p.T);
        p.g[k] = q.h / (R * p.T) - q.s / R + lnP;
        p.u[k] = q.h / (R * p.T * p.T);
    }
}

/// The moles of the species present at (pi, nu) [lnn: their logarithms], their sum, the residuals of the element balances
/// r_c = b_c - sum_k a_ck n_k and the matrix M = sum_k n_k a_k a_k^T. False when a value is not finite.
struct State {
    double n[kSpeciesCount], lnn[kSpeciesCount];
    double sum;
    double r[kElements];
    double M[kElements][kElements];
    double residual;   ///< max_c |r_c| / b_c
};

bool evaluate(const Problem& p, const double* pi, double nu, State& s) {
    s.sum = 0.0;
    for (int k = 0; k < p.ns; ++k) {
        double e = nu - p.g[k];
        for (int c = 0; c < p.ne; ++c) e += p.a[c][k] * pi[c];
        s.lnn[k] = e;
        s.n[k] = std::exp(e);
        s.sum += s.n[k];
    }
    if (!std::isfinite(s.sum)) return false;
    s.residual = 0.0;
    for (int c = 0; c < p.ne; ++c) {
        double sum = 0.0;
        for (int k = 0; k < p.ns; ++k) sum += p.a[c][k] * s.n[k];
        s.r[c] = p.b[c] - sum;
        s.residual = std::max(s.residual, std::abs(s.r[c]) / p.b[c]);
        for (int d = 0; d <= c; ++d) {
            double m = 0.0;
            for (int k = 0; k < p.ns; ++k) m += p.a[c][k] * p.a[d][k] * s.n[k];
            s.M[c][d] = s.M[d][c] = m;
        }
    }
    return true;
}

/// Solves (M + ridge) x = rhs for the symmetric positive semi-definite M of size n: Cholesky factorisation of the Jacobi-scaled
/// matrix D^-1/2 M D^-1/2 with kRidge added to its unit diagonal. The ridge bounds the step in the direction M has (almost) no
/// weight in: the stoichiometric mixture at a low temperature, where the species that tell an excess of oxygen from a lack of
/// it are 1E-40 of the mixture and the balances of O, C and H are, to the rounding, one condition. False when M has a zero or
/// non-finite diagonal.
bool solveScaled(const double M[kElements][kElements], int n, double* x) {
    double d[kElements], L[kElements][kElements];
    for (int i = 0; i < n; ++i) {
        if (!(M[i][i] > 0.0) || !std::isfinite(M[i][i])) return false;
        d[i] = 1.0 / std::sqrt(M[i][i]);
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j <= i; ++j) {
            double v = M[i][j] * d[i] * d[j] + (i == j ? kRidge : 0.0);
            for (int k = 0; k < j; ++k) v -= L[i][k] * L[j][k];
            if (i == j) {
                if (!(v > 0.0)) return false;
                L[i][i] = std::sqrt(v);
            } else {
                L[i][j] = v / L[j][j];
            }
        }
    }
    double y[kElements];
    for (int i = 0; i < n; ++i) {
        double v = x[i] * d[i];
        for (int k = 0; k < i; ++k) v -= L[i][k] * y[k];
        y[i] = v / L[i][i];
    }
    for (int i = n - 1; i >= 0; --i) {
        double v = y[i];
        for (int k = i + 1; k < n; ++k) v -= L[k][i] * y[k];
        y[i] = v / L[i][i];
    }
    for (int i = 0; i < n; ++i) x[i] = y[i] * d[i];
    return true;
}

/// The starting point: the complete combustion (lean: C to CO2, H to H2O, the rest of the oxygen to O2; rich: the oxygen to CO,
/// then to H2O and CO2, the rest of the hydrogen to H2 and of the carbon to CH4), N to N2, Ar to Ar, every present species
/// having at least a small share. The potentials make the most abundant independent species (one per element) exact.
void startingPoint(const Problem& p, double* pi, double& nu) {
    double b[kElements] = {};
    for (int c = 0; c < p.ne; ++c) b[p.element[c]] = p.b[c];
    double amount[kSpeciesCount] = {};   // by output index
    auto set = [&](const char* name, double v) {
        for (int j = 0; j < kSpeciesCount; ++j)
            if (std::string(kSpeciesNames[j]) == name) amount[j] = std::max(v, 0.0);
    };
    set("Ar", b[kAr]);
    set("N2", 0.5 * b[kN]);
    const double O = b[kO], C = b[kC], H = b[kH];
    if (O >= 2.0 * C + 0.5 * H) {   // lean or stoichiometric
        set("CO2", C);
        set("H2O", 0.5 * H);
        set("O2", 0.5 * (O - 2.0 * C - 0.5 * H));
    } else if (O >= C) {   // rich, the carbon still all in CO and CO2
        const double h2o = std::min(0.5 * H, O - C);
        const double co2 = std::min(C, O - C - h2o);
        set("CO", C - co2);
        set("CO2", co2);
        set("H2O", h2o);
        set("H2", 0.5 * H - h2o);
    } else {   // very rich: the oxygen in CO, the rest of the carbon in CH4 (4 C - 4 O <= H, see setup)
        set("CO", O);
        set("CH4", C - O);
        set("H2", 0.5 * H - 2.0 * (C - O));
    }
    double total = 0.0;
    for (double v : amount) total += v;
    for (int k = 0; k < p.ns; ++k) amount[p.species[k]] = std::max(amount[p.species[k]], 1e-6 * total);
    total = 0.0;
    for (int k = 0; k < p.ns; ++k) total += amount[p.species[k]];
    nu = std::log(total);

    // the basis: the most abundant species whose compositions are independent (elimination on the atoms)
    int order[kSpeciesCount];
    for (int k = 0; k < p.ns; ++k) order[k] = k;
    std::stable_sort(order, order + p.ns, [&](int x, int y) { return amount[p.species[x]] > amount[p.species[y]]; });
    double rows[kElements][kElements + 1];   // a_k^T pi = ln n_k - nu + g_k, reduced
    int pivotColumn[kElements];
    int rank = 0;
    for (int idx = 0; idx < p.ns && rank < p.ne; ++idx) {
        const int k = order[idx];
        double row[kElements + 1];
        for (int c = 0; c < p.ne; ++c) row[c] = p.a[c][k];
        row[p.ne] = std::log(amount[p.species[k]]) - nu + p.g[k];
        for (int r = 0; r < rank; ++r) {
            const double f = row[pivotColumn[r]] / rows[r][pivotColumn[r]];
            if (f != 0.0)
                for (int c = 0; c <= p.ne; ++c) row[c] -= f * rows[r][c];
        }
        int best = -1;
        for (int c = 0; c < p.ne; ++c)
            if (std::abs(row[c]) > 1e-9 && (best < 0 || std::abs(row[c]) > std::abs(row[best]))) best = c;
        if (best < 0) continue;   // dependent on the species kept
        for (int c = 0; c <= p.ne; ++c) rows[rank][c] = row[c];
        pivotColumn[rank] = best;
        ++rank;
    }
    if (rank < p.ne) throw std::logic_error("CHEM_EQUIL: the species do not span the elements");
    for (int r = rank - 1; r >= 0; --r) {   // back-substitution
        double v = rows[r][p.ne];
        for (int s = r + 1; s < rank; ++s) v -= rows[r][pivotColumn[s]] * pi[pivotColumn[s]];
        pi[pivotColumn[r]] = v / rows[r][pivotColumn[r]];
    }
}

/// The inner problem at nu: the potentials that maximise the concave function D(pi) = sum_c b_c pi_c - sum_k n_k(pi, nu), whose
/// gradient is the residual of the element balances and whose Hessian is -M. Its maximum exists for every nu (the elements fit
/// in the species: setup), so Newton's method with a line search on D converges from any start. `pi` and `s` are updated.
bool innerSolve(const Problem& p, double nu, double* pi, State& s, int& steps) {
    if (!evaluate(p, pi, nu, s)) return false;
    auto objective = [&](const double* v, const State& st) {
        double f = -st.sum;
        for (int c = 0; c < p.ne; ++c) f += p.b[c] * v[c];
        return f;
    };
    for (int iteration = 0; iteration < kMaxInnerSteps; ++iteration) {
        if (s.residual < kInnerTolerance) return true;
        double step[kElements];
        for (int c = 0; c < p.ne; ++c) step[c] = s.r[c];
        if (!solveScaled(s.M, p.ne, step)) return false;
        double largest = 0.0, slope = 0.0;
        for (int c = 0; c < p.ne; ++c) {
            largest = std::max(largest, std::abs(step[c]));
            slope += s.r[c] * step[c];
        }
        double lambda = largest > kMaxPotentialStep ? kMaxPotentialStep / largest : 1.0;
        const double f = objective(pi, s);
        State trial;
        double piTrial[kElements];
        bool accepted = false;
        for (int k = 0; k < kMaxBacktracks; ++k, lambda *= 0.5) {
            for (int c = 0; c < p.ne; ++c) piTrial[c] = pi[c] + lambda * step[c];
            if (!evaluate(p, piTrial, nu, trial)) continue;
            // the ascent of D, or, where the rounding of D hides it (near the maximum), a smaller residual
            if (objective(piTrial, trial) >= f + 1e-4 * lambda * slope || trial.residual < 0.5 * s.residual) {
                accepted = true;
                break;
            }
        }
        ++steps;
        if (!accepted) return false;
        for (int c = 0; c < p.ne; ++c) pi[c] = piTrial[c];
        s = trial;
    }
    return s.residual < kInnerTolerance;
}

[[noreturn]] void notConverged(const Problem& p, const std::string& what) {
    throw std::runtime_error("the equilibrium at P = " + numberText(p.P) + " Pa, T = " + numberText(p.T) + " K did not converge (" +
                             what + ")");
}

/// LU factorisation with partial pivoting of the bordered matrix [[M, b], [b^T, 0]] of the derivatives, solved in place.
struct Lu {
    double m[kMaxSize][kMaxSize];
    int perm[kMaxSize];
    int n = 0;
    bool factor() {
        for (int i = 0; i < n; ++i) perm[i] = i;
        for (int k = 0; k < n; ++k) {
            int pivot = k;
            for (int i = k + 1; i < n; ++i)
                if (std::abs(m[i][k]) > std::abs(m[pivot][k])) pivot = i;
            if (!(std::abs(m[pivot][k]) > 0.0)) return false;
            if (pivot != k) {
                std::swap(perm[k], perm[pivot]);
                for (int j = 0; j < n; ++j) std::swap(m[k][j], m[pivot][j]);
            }
            for (int i = k + 1; i < n; ++i) {
                m[i][k] /= m[k][k];
                for (int j = k + 1; j < n; ++j) m[i][j] -= m[i][k] * m[k][j];
            }
        }
        return true;
    }
    void solve(double* x) const {
        double y[kMaxSize];
        for (int i = 0; i < n; ++i) {
            y[i] = x[perm[i]];
            for (int k = 0; k < i; ++k) y[i] -= m[i][k] * y[k];
        }
        for (int i = n - 1; i >= 0; --i) {
            for (int k = i + 1; k < n; ++k) y[i] -= m[i][k] * y[k];
            y[i] /= m[i][i];
        }
        for (int i = 0; i < n; ++i) x[i] = y[i];
    }
};

}  // namespace

TemperatureRange temperatureRange(const Inputs& in) {
    double b[kElements];
    abundances(in, b);
    bool present[kSpeciesCount];
    presentSpecies(b, present);
    return rangeOf(present);
}

Result solve(const Inputs& in, bool derivatives) {
    Problem p;
    setup(in, p);
    const int ne = p.ne;

    double pi[kElements] = {};
    double nu = 0.0;
    startingPoint(p, pi, nu);

    // The outer problem: nu such that the moles of the species sum to N = exp(nu), psi(nu) = ln(sum n) - nu = 0. psi decreases,
    // with the slope -b^T M^-1 b / sum n, between -1 and 0 (the tangent of the potentials is pi' = -M^-1 b): Newton's method
    // inside the bracket found, bisection where a step leaves it; the potentials follow the tangent from one nu to the next.
    State s;
    int steps = 0;
    double lo = -INFINITY, hi = INFINITY;   // psi > 0 at lo, psi < 0 at hi
    for (int outer = 0;; ++outer) {
        if (outer == kMaxOuterSteps) notConverged(p, "the moles of the mixture");
        if (!innerSolve(p, nu, pi, s, steps)) notConverged(p, "the element potentials");
        const double psi = std::log(s.sum) - nu;
        if (std::abs(psi) < kOuterTolerance) break;
        (psi > 0.0 ? lo : hi) = nu;
        // the bracket closed on nu to the rounding (the stoichiometric mixture at a low temperature, where the balances hold
        // sum n to kInnerTolerance only): the root is found
        if (hi - lo <= 4e-16 * std::max(1.0, std::abs(nu)) && std::abs(psi) < 1e3 * kOuterTolerance) break;
        double tangent[kElements];
        for (int c = 0; c < ne; ++c) tangent[c] = p.b[c];
        if (!solveScaled(s.M, ne, tangent)) notConverged(p, "the moles of the mixture");
        double bMb = 0.0;
        for (int c = 0; c < ne; ++c) bMb += p.b[c] * tangent[c];
        const double slope = -bMb / s.sum;
        double next = slope < 0.0 ? nu - psi / slope : nu + (psi > 0.0 ? kMaxNuStep : -kMaxNuStep);
        if (!(std::isfinite(lo) && std::isfinite(hi))) next = std::max(nu - kMaxNuStep, std::min(nu + kMaxNuStep, next));
        else if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        for (int c = 0; c < ne; ++c) pi[c] -= tangent[c] * (next - nu);
        nu = next;
    }

    Result result{};
    result.iterations = steps;
    for (int k = 0; k < p.ns; ++k) result.x[p.species[k]] = s.n[k] / s.sum;
    if (!derivatives) return result;

    // The implicit function theorem at the solution of G_c = sum_k a_ck n_k - b_c = 0, G_N = sum_k n_k - exp(nu) = 0:
    // [[M, b], [b^T, 0]] d(pi, nu)/d theta = -dG/d theta, then d ln x_k = -dg_k + sum_c a_ck dpi_c (x_k = exp(-g_k + a_k^T pi)).
    // In the null direction of M (above) only the trace species move: the ridge keeps their derivatives bounded.
    Lu lu;
    lu.n = ne + 1;
    double scale[kElements];
    for (int c = 0; c < ne; ++c) scale[c] = std::sqrt(s.M[c][c]);
    for (int c = 0; c < ne; ++c) {
        for (int d = 0; d < ne; ++d) lu.m[c][d] = s.M[c][d] / (scale[c] * scale[d]) + (c == d ? kRidge : 0.0);
        lu.m[c][ne] = lu.m[ne][c] = p.b[c] / scale[c];
    }
    lu.m[ne][ne] = 0.0;
    if (!lu.factor()) notConverged(p, "the derivatives: singular system");
    for (int input = 0; input < kInputCount; ++input) {
        double rhs[kMaxSize] = {};
        double dg[kSpeciesCount] = {};   // dg_k / d input
        if (input == 0) {   // P: dn_k/dP = -n_k/P
            for (int c = 0; c < ne; ++c) rhs[c] = p.b[c] / p.P;
            rhs[ne] = s.sum / p.P;
            for (int k = 0; k < p.ns; ++k) dg[k] = 1.0 / p.P;
        } else if (input == 1) {   // T: dn_k/dT = n_k u_k
            for (int c = 0; c < ne; ++c) {
                double v = 0.0;
                for (int k = 0; k < p.ns; ++k) v += p.a[c][k] * s.n[k] * p.u[k];
                rhs[c] = -v;
            }
            double v = 0.0;
            for (int k = 0; k < p.ns; ++k) v += s.n[k] * p.u[k];
            rhs[ne] = -v;
            for (int k = 0; k < p.ns; ++k) dg[k] = -p.u[k];
        } else {   // the ratio of an element: dG_c/db_c = -1
            int c = -1;
            for (int e = 0; e < ne; ++e)
                if (p.element[e] == input - 1) c = e;
            if (c < 0) continue;   // an absent element: its species are 0, their derivatives too
            rhs[c] = 1.0;
        }
        for (int c = 0; c < ne; ++c) rhs[c] /= scale[c];
        lu.solve(rhs);
        for (int c = 0; c < ne; ++c) rhs[c] /= scale[c];   // dpi_c
        for (int k = 0; k < p.ns; ++k) {
            double dlnx = -dg[k];
            for (int c = 0; c < ne; ++c) dlnx += p.a[c][k] * rhs[c];
            result.dx[p.species[k]][input] = result.x[p.species[k]] * dlnx;
        }
    }
    return result;
}

}  // namespace coolsolve::chem_equil
