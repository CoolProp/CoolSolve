#include "coolsolve/nh3h2o.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

// Phase equilibrium and the three-property flash of the ammonia - water mixture on the IAPWS G4-01
// Helmholtz formulation (nh3h2o.cpp). Structure:
//   1. phase quantities as first-order jets in (T, rho, x), from the second-order jet of Phi;
//   2. numerical tools: small dense solver, damped Newton, Brent;
//   3. the saturation solver: 4 equations (equal p, mu_1, mu_2) in the 6 unknowns (T, p, xL, xV, rhoL, rhoV), two
//      of which are given (bubble and dew point at T or p, or the pair at T and p), with guesses from the liquid
//      fugacities against an ideal vapour and a warm-start cache;
//   4. the 15 flashes, each an outer root search over p or T around the saturation solver, then a full Newton
//      polish and the implicit function theorem on the converged system for the exact derivatives.

namespace coolsolve::nh3h2o {

double massFractionOf(double x) { return x * kM2 / molarMass(x); }
double dMassFractionDx(double x) {
    const double m = molarMass(x);
    return kM1 * kM2 / (m * m);
}
double moleFractionOf(double xi) { return xi * kM1 / ((1.0 - xi) * kM2 + xi * kM1); }

namespace {

/// The mole fraction of water 1 - x for the mass fraction of ammonia `xi`, from the mass fraction of water 1 - xi (exact
/// for xi >= 0.5), not by subtracting the rounded mole fraction of ammonia from 1.
double moleFractionComplementOf(double xi) { return (1.0 - xi) * kM2 / ((1.0 - xi) * kM2 + xi * kM1); }

using taylor::Jet;
using J1 = Jet<1>;

constexpr double kXMin = 1e-10;   ///< composition floor of the equilibrium solvers (mole fraction)
constexpr double kPRef = 1.0e5;   ///< reference pressure of the ideal-gas constants of the guesses [Pa]

double clampX(double x) { return std::min(std::max(x, kXMin), 1.0 - kXMin); }

/// Clamps a mole fraction of ammonia x and its complement xc (the mole fraction of water) together: both at least kXMin.
void clampX(double& x, double& xc) {
    if (x < kXMin) {
        x = kXMin;
        xc = 1.0 - kXMin;
    } else if (xc < kXMin) {
        xc = kXMin;
        x = 1.0 - kXMin;
    }
}

/// x1 < x2 for two mole fractions of ammonia given with their complements c1, c2, decided on the representation that
/// resolves them: the complements when both are ammonia-rich (a double near 1 has a spacing of 1.1e-16).
bool lessX(double x1, double c1, double x2, double c2) { return x1 >= 0.5 && x2 >= 0.5 ? c2 < c1 : x1 < x2; }

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error("nh3h2o: " + message); }

std::string fmt(double v) {
    std::ostringstream os;
    os.precision(10);
    os << v;
    return os.str();
}

// ---------------------------------------------------------------------------------------------------------------
// 1. Phase quantities
// ---------------------------------------------------------------------------------------------------------------

/// The quantities of one homogeneous phase at (T, rho, x) as first-order jets in (T, rho, x) (gradient indices 0, 1,
/// 2): pressure [Pa], chemical potentials of water (1) and ammonia (2) [J/mol], and the MASS-based h, s, u [J/kg,
/// J/(kg K)] and v [m3/kg]. x is the mole fraction of ammonia.
struct Phase {
    J1 p, mu1, mu2, h, s, u, v;
};

/// The phase at (T, rho, x) with the mole fraction of water xc = 1 - x given separately, evaluated without the rounding
/// of x (mu_1 contains ln(xc); see `reducedHelmholtz`). The x derivatives are those of the mixture of x and 1 - x.
Phase evaluatePhase(double T, double rho, double x, double xc) {
    clampX(x, xc);
    const Jet<2> jT = Jet<2>::variable(T, 0), jR = Jet<2>::variable(rho, 1), jX = Jet<2>::variable(x, 2);
    const Jet<2> f = kR * jT * reducedHelmholtz(jT, jR, jX, complementJet(jX, xc));
    const J1 T1 = J1::variable(T, 0), r1 = J1::variable(rho, 1), x1 = J1::variable(x, 2), xc1 = complementJet(x1, xc);
    const J1 f1(f), fT = derivative(f, 0), fR = derivative(f, 1), fX = derivative(f, 2);
    const J1 p = r1 * r1 * fR;
    const J1 u = f1 - T1 * fT;
    const J1 g = f1 + p / r1;
    const J1 invM = 1.0 / (xc1 * kM1 + x1 * kM2);
    Phase ph;
    ph.p = p;
    ph.mu1 = g - x1 * fX;
    ph.mu2 = g + xc1 * fX;
    ph.u = u * invM;
    ph.h = (u + p / r1) * invM;
    ph.s = (-fT) * invM;
    ph.v = invM / r1;
    return ph;
}

/// The phase at (T, rho, x) for a mole fraction that is itself the given value or the variable (a single phase, a
/// start): 1 - x is then its exact complement (the subtraction is exact for x >= 0.5).
Phase evaluatePhase(double T, double rho, double x) { return evaluatePhase(T, rho, x, 1.0 - x); }

struct PressureSlope {
    double p, dpdrho;
};

/// p and dp/drho at fixed T and x (the density solvers).
PressureSlope pressureSlope(double T, double rho, double x) {
    const Jet<2> jT = Jet<2>::constant(T), jR = Jet<2>::variable(rho, 1), jX = Jet<2>::constant(clampX(x));
    const Jet<2> f = kR * jT * reducedHelmholtz(jT, jR, jX);
    const J1 fR = derivative(f, 1);
    return {rho * rho * fR.value(), 2.0 * rho * fR.value() + rho * rho * fR.gradient(1)};
}

// ---------------------------------------------------------------------------------------------------------------
// 2. Numerical tools
// ---------------------------------------------------------------------------------------------------------------

/// Solves A X = B for n <= 8 (Gauss-Jordan with partial pivoting); B is n x m. False if A is singular.
bool solveLinear(int n, double a[8][8], int m, double b[8][8]) {
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int r = col + 1; r < n; ++r)
            if (std::abs(a[r][col]) > std::abs(a[pivot][col])) pivot = r;
        if (!(std::abs(a[pivot][col]) > 1e-300)) return false;
        if (pivot != col) {
            for (int j = 0; j < n; ++j) std::swap(a[pivot][j], a[col][j]);
            for (int j = 0; j < m; ++j) std::swap(b[pivot][j], b[col][j]);
        }
        const double inv = 1.0 / a[col][col];
        for (int j = 0; j < n; ++j) a[col][j] *= inv;
        for (int j = 0; j < m; ++j) b[col][j] *= inv;
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            const double factor = a[r][col];
            if (factor == 0.0) continue;
            for (int j = 0; j < n; ++j) a[r][j] -= factor * a[col][j];
            for (int j = 0; j < m; ++j) b[r][j] -= factor * b[col][j];
        }
    }
    return true;
}

struct NewtonSetup {
    int n = 0;
    double rowScale[8];
    double colScale[8];
    bool logVar[8] = {false, false, false, false, false, false, false, false};
    double maxStep[8];   ///< largest step: fraction of colScale, or |d ln v| for a logarithmic variable
    double lower[8], upper[8];
    /// The change of each Newton variable that moves the evaluated quantity by one representable step (default:
    /// one ulp of the variable). For a composition iterated as logit(x) it is set to that of x itself,
    /// ulp(x) / (x (1 - x)): the quantisation of a mole fraction held as a double, which `frozenPolish` was written for.
    /// Since the solvers evaluate x and 1 - x from the logit (`satNewton`), it overestimates the resolution of a
    /// composition near pure ammonia (up to 1e5 times ulp(z)); it is kept with `frozenPolish` (docs/nh3h2o_review.md §4).
    std::function<void(const double*, double*)> resolution;
    /// Variables that are compositions: their quantisation is in `resolution`; the evaluation jitter is measured with
    /// the other variables only.
    bool composition[8] = {false, false, false, false, false, false, false, false};
};

using Evaluator = std::function<bool(const double*, double*, double (*)[8])>;

/// The smallest residual that can be told from noise, per row and scaled by the row scale, the sum of
///  - the resolution of the representation of the variables: the change of F produced by moving any variable by one
///    representable step, sum |J_ij| u_j (u from `NewtonSetup::resolution`); and
///  - the round-off jitter of the evaluation of F, measured at the point: the largest deviation, over a few
///    perturbations of the non-composition variables by 1 to 3 ulp, of F from its linearisation (the largest sample
///    once, not their sum: each sample is one draw of the same jitter).
bool noiseFloor(const NewtonSetup& ns, const double* v, const double* F, double (*J)[8], const Evaluator& eval,
                double* floorOut) {
    const int n = ns.n;
    double u[8];
    for (int j = 0; j < n; ++j) u[j] = std::abs(v[j]) * 2.220446049250313e-16;
    if (ns.resolution) ns.resolution(v, u);
    for (int i = 0; i < n; ++i) {
        double f = 0.0;
        for (int j = 0; j < n; ++j) f += std::abs(J[i][j]) * u[j];
        floorOut[i] = f / ns.rowScale[i];
    }
    const double kUlp = 2.220446049250313e-16;
    const double pattern[6][8] = {{1, -1, 2, -2, 3, -3, 1, -1}, {-1, 2, -3, 1, -2, 3, -1, 2}, {2, 3, 1, -1, -3, -2, 3, 1},
                                  {-3, 1, -2, 2, 1, -1, -2, 3}, {1, 1, -1, -1, 2, 2, -3, -3}, {3, -2, 1, -3, 2, -1, 1, -2}};
    double jitter[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (const auto& c : pattern) {
        double vp[8], Fp[8], Jp[8][8];
        for (int j = 0; j < n; ++j) vp[j] = ns.composition[j] ? v[j] : v[j] * (1.0 + c[j] * kUlp);
        if (!eval(vp, Fp, Jp)) return false;
        for (int i = 0; i < n; ++i) {
            double lin = F[i];
            for (int j = 0; j < n; ++j) lin += J[i][j] * (vp[j] - v[j]);
            jitter[i] = std::max(jitter[i], std::abs(Fp[i] - lin) / ns.rowScale[i]);
        }
    }
    for (int i = 0; i < n; ++i) floorOut[i] += jitter[i];
    return true;
}

/**
 * Damped Newton iteration on F(v) = 0 with the Jacobian supplied by `eval` (false from `eval` = step refused).
 *
 * Convergence is a scaled residual below `tol`. The other exits are for an iteration that has stopped lowering the
 * residual (a step that does not halve it): the point is accepted if every row is within `tol` plus two units of the
 * noise floor measured there (`noiseFloor`: the resolution of the representation of the variables and the round-off
 * jitter of the evaluation), i.e. when no representable point can do better; otherwise `frozenPolish` is tried.
 */
bool newtonSolve(const NewtonSetup& ns, double* v, const Evaluator& eval, double tol, int maxIter) {
    const int n = ns.n;
    double F[8], J[8][8];
    if (!eval(v, F, J)) return false;
    auto merit = [&](const double* r) {
        double m = 0.0;
        for (int i = 0; i < n; ++i) m = std::max(m, std::abs(r[i]) / ns.rowScale[i]);
        return m;
    };
    auto atResolutionFloor = [&]() {
        double floorOut[8];
        if (!noiseFloor(ns, v, F, J, eval, floorOut)) return false;
        for (int i = 0; i < n; ++i)
            if (std::abs(F[i]) / ns.rowScale[i] > tol + 2.0 * floorOut[i]) return false;
        return true;
    };
    // Written for a variable quantised coarsely (a mole fraction near 1 held as a double, before the solvers evaluated 1 - x
    // from the logit): the full Newton step asks it for a change that a double cannot represent and the other variables
    // overshoot in response. Freeze that variable and the one row it dominates, solve the other rows strictly for the other
    // variables, and move the frozen variable by Newton on the dropped row, whose sensitivity is then the Schur complement
    // S = J[rq][jq] - J[rq][cols] M^-1 J[rows][jq]. The dropped row is accepted within two units of its quantisation floor
    // |S| u (`NewtonSetup::resolution`, plus the measured noise); all the other rows satisfy `tol`. It is tried at every
    // stall, also far from the solution (a second iteration there), and its quantisation floor is now an overestimate: kept
    // until the maintainer decides on its removal (docs/nh3h2o_review.md §4).
    auto frozenPolish = [&]() {
        double u[8];
        for (int j = 0; j < n; ++j) u[j] = std::abs(v[j]) * 2.220446049250313e-16;
        if (ns.resolution) ns.resolution(v, u);
        int jq = -1, rq = -1;
        double best = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                const double c = std::abs(J[i][j]) * u[j] / ns.rowScale[i];
                if (c > best) {
                    best = c;
                    jq = j;
                    rq = i;
                }
            }
        if (jq < 0 || best < 10.0 * tol || n < 2) return false;   // not a quantisation problem
        double w[8], Fw[8], Jw[8][8];
        for (int j = 0; j < n; ++j) w[j] = v[j];
        for (int i = 0; i < n; ++i) {
            Fw[i] = F[i];
            for (int j = 0; j < n; ++j) Jw[i][j] = J[i][j];
        }
        int rows[8], cols[8], m = 0;
        for (int i = 0; i < n; ++i)
            if (i != rq) rows[m++] = i;
        m = 0;
        for (int j = 0; j < n; ++j)
            if (j != jq) cols[m++] = j;
        auto reducedMerit = [&]() {
            double r = 0.0;
            for (int i = 0; i < n; ++i)
                if (i != rq) r = std::max(r, std::abs(Fw[i]) / ns.rowScale[i]);
            return r;
        };
        auto step = [&](int j, double d) {
            w[j] = ns.logVar[j] ? w[j] * std::exp(d / w[j]) : w[j] + d;
            w[j] = std::min(std::max(w[j], ns.lower[j]), ns.upper[j]);
        };
        double quantisation = 0.0;
        for (int outer = 0; outer < 6; ++outer) {
            // the other variables for the other rows (Newton on the reduced system)
            for (int it = 0; it < 6 && reducedMerit() >= tol; ++it) {
                double A[8][8], B[8][8];
                for (int a = 0; a < n - 1; ++a) {
                    for (int b = 0; b < n - 1; ++b) A[a][b] = Jw[rows[a]][cols[b]] / ns.rowScale[rows[a]];
                    B[a][0] = -Fw[rows[a]] / ns.rowScale[rows[a]];
                }
                if (!solveLinear(n - 1, A, 1, B)) return false;
                for (int b = 0; b < n - 1; ++b) step(cols[b], B[b][0]);
                if (!eval(w, Fw, Jw)) return false;
            }
            if (!(reducedMerit() < tol)) return false;
            // the dropped row: its Schur complement and quantisation floor
            double Mx[8][8], c[8][8];
            for (int a = 0; a < n - 1; ++a) {
                for (int b = 0; b < n - 1; ++b) Mx[a][b] = Jw[rows[a]][cols[b]];
                c[a][0] = Jw[rows[a]][jq];
            }
            if (!solveLinear(n - 1, Mx, 1, c)) return false;
            double S = Jw[rq][jq];
            for (int b = 0; b < n - 1; ++b) S -= Jw[rq][cols[b]] * c[b][0];
            double uq[8];
            for (int j = 0; j < n; ++j) uq[j] = std::abs(w[j]) * 2.220446049250313e-16;
            if (ns.resolution) ns.resolution(w, uq);
            quantisation = std::abs(S) * uq[jq] / ns.rowScale[rq];
            double floorOut[8];
            if (!noiseFloor(ns, w, Fw, Jw, eval, floorOut)) return false;
            const double allowed = tol + 2.0 * std::max(floorOut[rq], quantisation);
            if (std::abs(Fw[rq]) / ns.rowScale[rq] <= allowed) {
                for (int j = 0; j < n; ++j) v[j] = w[j];
                for (int i = 0; i < n; ++i) F[i] = Fw[i];
                return true;
            }
            if (!(std::abs(S) > 0.0)) return false;
            step(jq, -Fw[rq] / S);   // Newton on the frozen variable, the others eliminated
            if (!eval(w, Fw, Jw)) return false;
        }
        return false;
    };
    double m0 = merit(F);
    int stalled = 0;
    for (int it = 0; it < maxIter; ++it) {
        if (!(m0 >= tol)) {
            return std::isfinite(m0);
        }
        if (stalled >= 1 && (atResolutionFloor() || frozenPolish())) return true;
        double A[8][8], B[8][8];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) A[i][j] = J[i][j] / ns.rowScale[i];
            B[i][0] = -F[i] / ns.rowScale[i];
        }
        if (!solveLinear(n, A, 1, B)) return false;
        double d[8], lam = 1.0;
        for (int j = 0; j < n; ++j) {
            d[j] = B[j][0];
            const double rel = ns.logVar[j] ? std::abs(d[j]) / std::max(std::abs(v[j]), 1e-300) : std::abs(d[j]) / ns.colScale[j];
            if (rel > ns.maxStep[j]) lam = std::min(lam, ns.maxStep[j] / rel);
        }
        double vn[8], Fn[8], Jn[8][8], mn = 0.0;
        bool accepted = false;
        for (int ls = 0; ls < 12 && !accepted; ++ls, lam *= 0.5) {
            for (int j = 0; j < n; ++j) {
                vn[j] = ns.logVar[j] ? v[j] * std::exp(lam * d[j] / v[j]) : v[j] + lam * d[j];
                vn[j] = std::min(std::max(vn[j], ns.lower[j]), ns.upper[j]);
            }
            if (!eval(vn, Fn, Jn)) continue;
            mn = merit(Fn);
            if (std::isfinite(mn) && (mn < m0 || ls == 11)) accepted = true;
        }
        if (!accepted) return false;
        stalled = mn < 0.5 * m0 ? 0 : stalled + 1;   // the iteration no longer lowers the residual
        if (mn >= m0) {   // keep the better point; it is the answer if nothing representable does better
            return atResolutionFloor() || frozenPolish();
        }
        for (int j = 0; j < n; ++j) v[j] = vn[j];
        for (int i = 0; i < n; ++i) {
            F[i] = Fn[i];
            for (int j = 0; j < n; ++j) J[i][j] = Jn[i][j];
        }
        m0 = mn;
    }
    return m0 < tol;
}

/// Brent's root finder on [a, b] with f(a) and f(b) of opposite signs (or zero).
bool brent(const std::function<double(double)>& f, double a, double b, double fa, double fb, double tol, double& root,
           int maxIter = 200) {
    if (fa == 0.0) return root = a, true;
    if (fb == 0.0) return root = b, true;
    if (!(fa * fb < 0.0)) return false;
    double c = a, fc = fa, d = b - a, e = d;
    for (int it = 0; it < maxIter; ++it) {
        if (fb * fc > 0.0) {
            c = a;
            fc = fa;
            d = e = b - a;
        }
        if (std::abs(fc) < std::abs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        const double tolAct = 2.0 * std::numeric_limits<double>::epsilon() * std::abs(b) + 0.5 * tol;
        const double xm = 0.5 * (c - b);
        if (std::abs(xm) <= tolAct || fb == 0.0) return root = b, true;
        if (std::abs(e) >= tolAct && std::abs(fa) > std::abs(fb)) {
            double p, q;
            const double s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                const double qq = fa / fc, r = fb / fc;
                p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::abs(p);
            if (2.0 * p < std::min(3.0 * xm * q - std::abs(tolAct * q), std::abs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        b += std::abs(d) > tolAct ? d : (xm > 0 ? tolAct : -tolAct);
        fb = f(b);
    }
    return root = b, true;
}

// ---------------------------------------------------------------------------------------------------------------
// Single-phase density roots
// ---------------------------------------------------------------------------------------------------------------

/// A starting density of the liquid branch [mol/m3]: the saturated-liquid ancillary of water scaled to the
/// reducing density of the mixture.
double liquidDensityGuess(double T, double x) {
    const double theta = std::min(std::max(1.0 - T / reducingTemperature(x), 0.0), 1.0);
    const double delta = 1.0 + 1.99274064 * std::cbrt(theta) + 1.09965342 * std::pow(theta, 2.0 / 3.0) -
                         0.510839303 * std::pow(theta, 5.0 / 3.0);
    return delta * reducingDensity(x);
}

/// The density [mol/m3] of the liquid-like (`liquid`) or vapour-like root of p(T, rho, x) = p; false if the branch
/// has no root at this pressure.
bool densityRoot(double T, double p, double x, bool liquid, double& rho) {
    double r = liquid ? liquidDensityGuess(T, x) : p / (kR * T);
    for (int it = 0; it < 100; ++it) {
        PressureSlope ps;
        try {
            ps = pressureSlope(T, r, x);
        } catch (const std::runtime_error&) {
            return false;
        }
        if (!std::isfinite(ps.p)) return false;
        // the pressure of a liquid is a sum of terms of the size rho R T: its noise floor is relative to that
        if (std::abs(ps.p - p) <= 1e-12 * r * kR * T + 1e-13 * p && ps.dpdrho > 0.0) {
            rho = r;
            return true;
        }
        if (!(ps.dpdrho > 0.0)) {   // mechanically unstable: move away from the spinodal, towards the branch
            r *= liquid ? 1.08 : 0.92;
            continue;
        }
        double step = (p - ps.p) / ps.dpdrho;
        const double up = liquid ? 0.5 * r : 0.25 * r, down = liquid ? 0.2 * r : 0.5 * r;
        step = std::min(std::max(step, -down), up);
        r += step;
        if (!(r > 0.0)) return false;
    }
    return false;
}

// ---------------------------------------------------------------------------------------------------------------
// 3. Saturation: (T, p, xL, xV, rhoL, rhoV) with p_L = p_V = p and mu_1, mu_2 equal in both phases
// ---------------------------------------------------------------------------------------------------------------

struct Sat {
    double v[6] = {0, 0, 0, 0, 0, 0};   ///< T, p, xL, xV, rhoL, rhoV
    /// The mole fractions of water 1 - xL, 1 - xV, each evaluated to its own round-off (from the logit of the
    /// composition or from the specified value), never by subtracting the rounded xL, xV from 1.
    double w[2] = {1.0, 1.0};
    double T() const { return v[0]; }
    double p() const { return v[1]; }
    double xL() const { return v[2]; }
    double xV() const { return v[3]; }
    double rhoL() const { return v[4]; }
    double rhoV() const { return v[5]; }
    double wL() const { return w[0]; }
    double wV() const { return w[1]; }
    /// Sets the composition v[k] (k = 2 liquid, 3 vapour) and its complement, clamped together.
    void setX(int k, double x, double xc) {
        clampX(x, xc);
        v[k] = x;
        w[k - 2] = xc;
    }
};

enum SatKind { kBubbleT = 0, kDewT, kBubbleP, kDewP, kTP };
constexpr int kFixed[5][2] = {{0, 2}, {0, 3}, {1, 2}, {1, 3}, {0, 1}};
const char* const kKindName[5] = {"bubble point at T and x", "dew point at T and x", "bubble point at p and x",
                                  "dew point at p and x", "vapour-liquid equilibrium at T and p"};

/// The 4 equilibrium residuals and their 4 x 6 Jacobian at v = (T, p, xL, xV, rhoL, rhoV), with w = (1 - xL, 1 - xV)
/// the complements of the compositions.
void satResiduals(const double* v, const double* w, double F[4], double J[4][6], Phase* outL = nullptr,
                  Phase* outV = nullptr) {
    const Phase L = evaluatePhase(v[0], v[4], v[2], w[0]);
    const Phase V = evaluatePhase(v[0], v[5], v[3], w[1]);
    F[0] = L.p.value() - v[1];
    F[1] = V.p.value() - v[1];
    F[2] = L.mu1.value() - V.mu1.value();
    F[3] = L.mu2.value() - V.mu2.value();
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 6; ++j) J[i][j] = 0.0;
    J[0][0] = L.p.gradient(0);  J[0][1] = -1.0;  J[0][2] = L.p.gradient(2);  J[0][4] = L.p.gradient(1);
    J[1][0] = V.p.gradient(0);  J[1][1] = -1.0;  J[1][3] = V.p.gradient(2);  J[1][5] = V.p.gradient(1);
    J[2][0] = L.mu1.gradient(0) - V.mu1.gradient(0);
    J[2][2] = L.mu1.gradient(2);  J[2][3] = -V.mu1.gradient(2);
    J[2][4] = L.mu1.gradient(1);  J[2][5] = -V.mu1.gradient(1);
    J[3][0] = L.mu2.gradient(0) - V.mu2.gradient(0);
    J[3][2] = L.mu2.gradient(2);  J[3][3] = -V.mu2.gradient(2);
    J[3][4] = L.mu2.gradient(1);  J[3][5] = -V.mu2.gradient(1);
    if (outL) *outL = L;
    if (outV) *outV = V;
}

/// The logit ln(x / xc) of a mole fraction x with its complement xc (clamped together).
double logitOf(double x, double xc) {
    clampX(x, xc);
    return std::log(x / xc);
}
/// The mole fraction of the logit z, 1 / (1 + exp(-z)); its complement is sigmoid(-z), to its own round-off.
double sigmoid(double z) { return 1.0 / (1.0 + std::exp(-z)); }

/// The change of logit(x) that moves the double x by one representable step: the resolution of a composition that is
/// iterated as its logit.
double logitResolution(double z) {
    const double x = sigmoid(z);
    return (std::nextafter(x, 2.0) - x) / (x * (1.0 - x));
}

/// True if both phases are mechanically stable, dp/drho > 0 at fixed T and x. The equilibrium equations also have
/// roots with a phase between its spinodals (dp/drho < 0); such a root is not a physical state.
bool mechanicallyStable(const Phase& L, const Phase& V) { return L.p.gradient(1) > 0.0 && V.p.gradient(1) > 0.0; }

/// Newton on the saturation system for the free variables of `kind` (the other two keep their values in `s`). The
/// compositions are iterated as z = logit(x), which keeps the iteration well scaled up to the pure-component ends, and
/// evaluated as x = sigmoid(z) with the complement xc = sigmoid(-z): both are known to their own round-off, so ln(xc)
/// in mu_1 is resolved near pure ammonia as well. False if the iteration does not converge or converges to the
/// trivial solution, to swapped phases or to a mechanically unstable phase.
bool satNewton(Sat& s, SatKind kind, int maxIter = 60) {
    const int fixedA = kFixed[kind][0], fixedB = kFixed[kind][1];
    int free[4], nf = 0;
    for (int k = 0; k < 6; ++k)
        if (k != fixedA && k != fixedB) free[nf++] = k;
    NewtonSetup ns;
    ns.n = 4;
    // the pressure of a liquid is a sum of terms of the size rho R T: its noise is relative to that, not to p
    const double rowScale[4] = {std::max(s.v[1], s.v[4] * kR * s.v[0]), s.v[1], kR * 300.0, kR * 300.0};
    const double colScale[6] = {100.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    const double maxStep[6] = {0.05, 0.7, 1.5, 1.5, 0.25, 0.4};
    const double lower[6] = {120.0, 1e-3, -23.0, -23.0, 1.0, 1e-6};
    const double upper[6] = {1000.0, 2e9, 23.0, 23.0, 2e5, 2e5};
    auto isComposition = [](int k) { return k == 2 || k == 3; };
    for (int i = 0; i < 4; ++i) {
        ns.rowScale[i] = rowScale[i];
        ns.colScale[i] = colScale[free[i]];
        ns.logVar[i] = free[i] == 1 || free[i] == 4 || free[i] == 5;
        ns.composition[i] = isComposition(free[i]);
        ns.maxStep[i] = maxStep[free[i]];
        ns.lower[i] = lower[free[i]];
        ns.upper[i] = upper[free[i]];
    }
    ns.resolution = [&](const double* uv, double* out) {
        for (int i = 0; i < 4; ++i)
            out[i] = isComposition(free[i]) ? logitResolution(uv[i]) : std::abs(uv[i]) * 2.220446049250313e-16;
    };
    double u[4];
    for (int i = 0; i < 4; ++i)
        u[i] = isComposition(free[i]) ? logitOf(s.v[free[i]], s.w[free[i] - 2]) : s.v[free[i]];
    double full[6], fullW[2] = {s.w[0], s.w[1]};
    for (int k = 0; k < 6; ++k) full[k] = s.v[k];
    // the phases of the last evaluation, for the stability check of the root (most often the root itself)
    double lastU[4] = {0, 0, 0, 0};
    Phase lastL, lastV;
    bool haveLast = false;
    const Evaluator eval = [&](const double* uv, double* F, double (*J)[8]) {
        double dxdu[4];
        for (int i = 0; i < 4; ++i) {
            if (isComposition(free[i])) {
                full[free[i]] = sigmoid(uv[i]);
                fullW[free[i] - 2] = sigmoid(-uv[i]);
                dxdu[i] = full[free[i]] * fullW[free[i] - 2];
            } else {
                full[free[i]] = uv[i];
                dxdu[i] = 1.0;
            }
        }
        double Jf[4][6], Ff[4];
        haveLast = false;
        try {
            satResiduals(full, fullW, Ff, Jf, &lastL, &lastV);
        } catch (const std::runtime_error&) {
            return false;
        }
        for (int i = 0; i < 4; ++i) lastU[i] = uv[i];
        haveLast = true;
        for (int i = 0; i < 4; ++i) {
            F[i] = Ff[i];
            for (int j = 0; j < 4; ++j) J[i][j] = Jf[i][free[j]] * dxdu[j];
        }
        return true;
    };
    if (!newtonSolve(ns, u, eval, 1e-12, maxIter)) return false;
    for (int i = 0; i < 4; ++i) {
        if (isComposition(free[i])) {
            s.v[free[i]] = sigmoid(u[i]);
            s.w[free[i] - 2] = sigmoid(-u[i]);
        } else {
            s.v[free[i]] = u[i];
        }
    }
    // the trivial solution (one phase) satisfies every equation. Near a critical point the residuals are flat along the
    // trivial manifold, and Newton can stop on two copies of one phase (seen: ln(rhoL/rhoV) = 2.8e-5 and |xL - xV| =
    // 3.7e-8 for a dew point at 395 K and X = 0.999, 4 % above the true dew pressure; 4.2e-5 and 4.4e-6 for a bubble
    // point at 435 K and X = 0.9). A real pair is that close only within about 1e-9 of the critical point (0.2 and 2e-3
    // at 1e-4 of it, 410 K).
    if (std::abs(s.v[2] - s.v[3]) < 1e-4 && std::abs(std::log(s.v[4] / s.v[5])) < 1e-3) return false;
    if (!(s.v[4] > s.v[5])) return false;
    if (haveLast && std::equal(u, u + 4, lastU)) return mechanicallyStable(lastL, lastV);
    double F[4], J[4][6];
    Phase L, V;
    try {
        satResiduals(s.v, s.w, F, J, &L, &V);
    } catch (const std::runtime_error&) {
        return false;
    }
    return mechanicallyStable(L, V);
}

/// R T times the constants of the ideal-gas chemical potentials at the reference pressure: for the pure ideal gas i,
/// mu_i(T, p) = C_i(T) + R T ln p (water, ammonia).
void idealGasConstants(double T, double& C1, double& C2) {
    const double RT = kR * T, tau0 = kT0 / T, delta0 = kPRef / (RT * kRho0);
    C1 = RT * (idealHelmholtz(tau0, delta0, 0.0) + 1.0) - RT * std::log(kPRef);
    C2 = RT * (idealHelmholtz(tau0, delta0, 1.0) + 1.0) - RT * std::log(kPRef);
}

/// The partial pressures f_i that the liquid of composition xL at T would have in an ideal vapour, and its density.
/// The liquid is evaluated at its own bubble pressure against an ideal vapour (iterated from pStart), or at pStart
/// itself (`atStart`, the liquid of an equilibrium at a given pressure).
struct Fug {
    bool ok = false;
    double f1 = 0, f2 = 0, rhoL = 0;
};

Fug liquidFugacities(double T, double xL, double pStart = 1.0e5, bool atStart = false) {
    double C1, C2;
    idealGasConstants(T, C1, C2);
    const double RT = kR * T;
    Fug fug;
    double p = pStart;
    for (int it = 0, escalations = 0; it < (atStart ? 1 : 6); ++it) {
        double rho;
        // a liquid root can be missing at a low pressure near the critical temperature (beyond its spinodal)
        while (!densityRoot(T, p, xL, true, rho)) {
            if (atStart || ++escalations > 8) return fug;
            p *= 4.0;
        }
        if (rho < reducingDensity(xL)) return fug;   // above the critical temperature there is no liquid
        Phase ph;
        try {
            ph = evaluatePhase(T, rho, xL);
        } catch (const std::runtime_error&) {
            return fug;
        }
        fug.f1 = std::exp((ph.mu1.value() - C1) / RT);
        fug.f2 = std::exp((ph.mu2.value() - C2) / RT);
        fug.rhoL = rho;
        fug.ok = std::isfinite(fug.f1 + fug.f2) && fug.f1 + fug.f2 > 0.0;
        if (!fug.ok) return fug;
        const double pNew = fug.f1 + fug.f2;
        if (std::abs(pNew - p) < 1e-5 * p) break;
        p = pNew;
    }
    return fug;
}

double vaporDensityGuess(double T, double p, double y) {
    double rho;
    if (densityRoot(T, p, y, false, rho)) return rho;
    return p / (kR * T);
}

/// Fugacity coefficients of the vapour (T, p, y) against the ideal gas of the same composition,
/// phi_i = f_i / (y_i p) with f_i the partial pressure the formulation gives (the same quantity as `Fug` for a liquid).
/// They correct the ideal-vapour assumption of the guesses (about 2 % of the pressure of ammonia at 3 bar).
struct VaporPhi {
    bool ok = false;
    double phi1 = 1.0, phi2 = 1.0;
};

VaporPhi vaporPhi(double T, double p, double y) {
    VaporPhi out;
    double rho;
    if (!densityRoot(T, p, y, false, rho)) return out;
    Phase ph;
    try {
        ph = evaluatePhase(T, rho, y);
    } catch (const std::runtime_error&) {
        return out;
    }
    double C1, C2;
    idealGasConstants(T, C1, C2);
    const double RT = kR * T, yc = clampX(y);
    out.phi1 = std::exp((ph.mu1.value() - C1) / RT) / ((1.0 - yc) * p);
    out.phi2 = std::exp((ph.mu2.value() - C2) / RT) / (yc * p);
    out.ok = std::isfinite(out.phi1) && std::isfinite(out.phi2) && out.phi1 > 0.0 && out.phi2 > 0.0;
    return out;
}

/// Root of g(z) in the logit of the liquid composition, z in [-16, 16] (g is NaN where there is no liquid, e.g. pure
/// ammonia above its critical temperature). The ends are tried first; if one of them has no liquid the range is
/// scanned for the first sign change between two valid points.
bool solveLogitRoot(const std::function<double(double)>& g, double& z) {
    const double lo = -16.0, hi = 16.0, glo = g(lo), ghi = g(hi);
    if (std::isfinite(glo) && std::isfinite(ghi)) return brent(g, lo, hi, glo, ghi, 5e-3, z, 40);
    double zPrev = 0.0, gPrev = 0.0;
    bool havePrev = false;
    for (double zz = lo; zz <= hi + 1e-9; zz += 2.0) {
        const double gz = g(zz);
        if (!std::isfinite(gz)) {
            havePrev = false;
            continue;
        }
        if (havePrev && (gPrev == 0.0 || gPrev * gz < 0.0)) return brent(g, zPrev, zz, gPrev, gz, 5e-3, z, 40);
        zPrev = zz;
        gPrev = gz;
        havePrev = true;
    }
    return false;
}

bool guessBubbleT(double T, double xL, Sat& s, bool realVapour = true) {
    const Fug f = liquidFugacities(T, xL);
    if (!f.ok) return false;
    double p = f.f1 + f.f2, y = clampX(f.f2 / p), yc = f.f1 / p;
    for (int pass = 0; pass < (realVapour ? 3 : 0); ++pass) {   // the vapour is not an ideal gas
        const VaporPhi phi = vaporPhi(T, p, y);
        if (!phi.ok) break;
        p = f.f1 / phi.phi1 + f.f2 / phi.phi2;
        y = clampX((f.f2 / phi.phi2) / p);
        yc = (f.f1 / phi.phi1) / p;
    }
    s.v[0] = T; s.v[1] = p; s.v[4] = f.rhoL; s.v[5] = vaporDensityGuess(T, p, y);
    s.setX(2, xL, 1.0 - xL);
    s.setX(3, y, yc);
    return true;
}

bool guessDewT(double T, double yV, Sat& s, bool realVapour = true) {
    const double y = clampX(yV), target = std::log(y / (1.0 - y));
    double phi1 = 1.0, phi2 = 1.0, p = 0.0, xL = 0.0, xLc = 1.0;
    Fug f;
    for (int pass = 0; pass < (realVapour ? 3 : 1); ++pass) {   // the vapour is not an ideal gas: its fugacity coefficients follow p
        auto g = [&](double z) {
            const Fug fl = liquidFugacities(T, sigmoid(z));
            return fl.ok ? std::log((fl.f2 / phi2) / (fl.f1 / phi1)) - target : std::numeric_limits<double>::quiet_NaN();
        };
        double z;
        if (!solveLogitRoot(g, z)) return false;
        xL = sigmoid(z);
        xLc = sigmoid(-z);
        f = liquidFugacities(T, xL);
        if (!f.ok) return false;
        p = (f.f1 / phi1) / (1.0 - y);
        const VaporPhi phi = vaporPhi(T, p, y);
        if (!phi.ok) break;
        phi1 = phi.phi1;
        phi2 = phi.phi2;
    }
    s.v[0] = T; s.v[1] = p; s.v[4] = f.rhoL; s.v[5] = vaporDensityGuess(T, p, y);
    s.setX(2, xL, xLc);
    s.setX(3, y, 1.0 - y);
    return true;
}

bool guessTP(double T, double p, Sat& s, bool realVapour = true) {
    double phi1 = 1.0, phi2 = 1.0, xL = 0.0, xLc = 1.0, y = 0.0, yc = 1.0;
    bool atP = false;   // the liquid at its own ideal-vapour bubble pressure, or at p
    Fug f;
    for (int pass = 0; pass < (realVapour ? 3 : 1); ++pass) {   // the vapour is not an ideal gas: its fugacity coefficients follow the composition
        auto g = [&](double z) {
            const Fug fl = liquidFugacities(T, sigmoid(z), p, atP);
            return fl.ok ? std::log((fl.f1 / phi1 + fl.f2 / phi2) / p) : std::numeric_limits<double>::quiet_NaN();
        };
        double z;
        if (!solveLogitRoot(g, z)) {
            // within a few percent of the pressure of the pure ammonia the ideal vapour has no root: start from the
            // fugacity coefficients of an ammonia-rich vapour, as pure as the margin to the pure-ammonia pressure needs,
            // and from the liquid at p (its compression from its ideal-vapour bubble pressure to p, a fraction of a
            // percent of its fugacities, is of the size of that margin)
            if (pass > 0 || !realVapour) return false;
            bool found = false;
            atP = true;
            for (const double yRich : {0.999999, 0.9999, 0.99}) {
                const VaporPhi rich = vaporPhi(T, p, yRich);
                if (!rich.ok) continue;
                phi1 = rich.phi1;
                phi2 = rich.phi2;
                if ((found = solveLogitRoot(g, z))) break;
            }
            if (!found) return false;
        }
        xL = sigmoid(z);
        xLc = sigmoid(-z);
        f = liquidFugacities(T, xL, p, atP);
        if (!f.ok) return false;
        y = clampX((f.f2 / phi2) / (f.f1 / phi1 + f.f2 / phi2));
        yc = (f.f1 / phi1) / (f.f1 / phi1 + f.f2 / phi2);
        const VaporPhi phi = vaporPhi(T, p, y);
        if (!phi.ok) break;
        phi1 = phi.phi1;
        phi2 = phi.phi2;
    }
    s.v[0] = T; s.v[1] = p; s.v[4] = f.rhoL; s.v[5] = vaporDensityGuess(T, p, y);
    s.setX(2, xL, xLc);
    s.setX(3, y, yc);
    return true;
}

/// Lee-Kesler vapour pressure of a pure component [Pa] (a starting point only).
double leeKeslerPressure(double T, double Tc, double pc, double omega) {
    const double Tr = T / Tc, Tr6 = std::pow(Tr, 6);
    const double f0 = 5.92714 - 6.09648 / Tr - 1.28862 * std::log(Tr) + 0.169347 * Tr6;
    const double f1 = 15.2518 - 15.6875 / Tr - 13.4721 * std::log(Tr) + 0.43577 * Tr6;
    return pc * std::exp(f0 + omega * f1);
}

/// The temperature at which an ideal solution (Raoult's law) of mole fraction x has the bubble or dew pressure p.
double raoultTemperature(bool bubble, double p, double x) {
    auto pressure = [&](double T) {
        const double pw = leeKeslerPressure(T, kTc1, 22.064e6, 0.3443), pa = leeKeslerPressure(T, kTc2, 11.333e6, 0.2526);
        return bubble ? (1.0 - x) * pw + x * pa : 1.0 / ((1.0 - x) / pw + x / pa);
    };
    double lo = 120.0, hi = 900.0;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        (pressure(mid) < p ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

/// A guess for a saturation solved at the pressure p: the temperature is searched so that the guess at T has this
/// pressure (kind = bubble or dew at fixed composition), starting from Raoult's law.
bool guessAtPressure(SatKind kind, double p, double x, Sat& s, bool realVapour = true) {
    const double target = std::log(p);
    constexpr double kFailed = 10.0;
    auto g = [&](double T) {
        Sat trial;
        const bool ok = kind == kBubbleP ? guessBubbleT(T, x, trial, realVapour) : guessDewT(T, x, trial, realVapour);
        return ok ? std::log(trial.v[1]) - target : kFailed;   // beyond the critical temperature: "too hot"
    };
    const double T0 = raoultTemperature(kind == kBubbleP, p, clampX(x));
    double lo = 0.9 * T0, hi = 1.15 * T0, glo = g(lo), ghi = g(hi);
    for (int it = 0; it < 12 && !(glo < 0.0); ++it) {   // widen downwards until the guess pressure is below the target
        hi = lo; ghi = glo;
        lo *= 0.9;
        if (lo < 130.0) return false;
        glo = g(lo);
    }
    for (int it = 0; it < 12 && !(ghi > 0.0); ++it) {   // and upwards until it is above (or the guess ceases to exist)
        lo = hi; glo = ghi;
        hi = std::min(hi * 1.15, 1100.0);
        ghi = g(hi);
    }
    if (!(glo < 0.0) || !(ghi > 0.0)) return false;
    double T;
    if (!brent(g, lo, hi, glo, ghi, 0.5, T, 60)) return false;
    return kind == kBubbleP ? guessBubbleT(T, x, s, realVapour) : guessDewT(T, x, s, realVapour);
}

struct SatCache {
    bool valid = false;
    Sat s;
};
thread_local SatCache g_satCache[5];

/// The outcome of a saturation solve. `NoGuess`: every guess function failed (no liquid at this T, no bracket at this
/// p), the evidence that the state is outside the two-phase region or above the critical locus. `NotConverged`: a guess
/// existed but Newton failed from every start or its root was rejected; it says nothing about the phase of the state.
enum class SatStatus { Found, NoGuess, NotConverged };

/// The saturation state of `kind` for the specified values a, b (T and x, or p and x, or T and p) in `s`, with its
/// status (`s` is meaningful for `Found` only). For a composition b, `bc` is its complement 1 - b evaluated without the
/// rounding of b (unused for kTP). Never throws for a missing or failed saturation. The starts, in this order: the
/// warm-start cache, `hint` (a start supplied by the caller), the guesses. The stable root does not depend on the
/// start; the order only affects the cost.
SatStatus satAttempt(SatKind kind, double a, double b, double bc, Sat& s, const Sat* hint = nullptr) {
    SatCache& cache = g_satCache[kind];
    const bool pressureSpec = kind == kBubbleP || kind == kDewP;
    const bool pairTP = kind == kTP;
    auto specify = [&](Sat& t) {   // the two specified values
        t.v[kFixed[kind][0]] = a;
        if (pairTP)
            t.v[kFixed[kind][1]] = b;
        else
            t.setX(kFixed[kind][1], b, bc);
    };
    auto close = [&](const Sat& c) {
        const double ca = pressureSpec ? c.p() : c.T();
        if (std::abs(a / ca - 1.0) > 0.25) return false;
        if (pairTP) return std::abs(b / c.p() - 1.0) < 0.25;
        const double cx = kind == kBubbleT || kind == kBubbleP ? c.xL() : c.xV();
        return std::abs(b - cx) < 0.25;
    };
    if (cache.valid && close(cache.s)) {
        s = cache.s;
        specify(s);
        if (satNewton(s, kind)) {
            cache.s = s;
            return SatStatus::Found;
        }
    }
    bool anyGuess = false;
    if (hint) {
        anyGuess = true;
        s = *hint;
        specify(s);
        if (satNewton(s, kind)) {
            cache.valid = true;
            cache.s = s;
            return SatStatus::Found;
        }
    }
    // the guess with the vapour as a real gas first; the ideal-vapour guess (cruder, but it never over-corrects at high
    // pressure) when Newton does not converge from it
    double guessT = 0.0, guessP = 0.0;   // T and p of the first guess, with the specified one of them
    for (const bool realVapour : {true, false}) {
        bool ok = false;
        switch (kind) {
            case kBubbleT: ok = guessBubbleT(a, b, s, realVapour); break;
            case kDewT: ok = guessDewT(a, b, s, realVapour); break;
            case kBubbleP: ok = guessAtPressure(kBubbleP, a, b, s, realVapour); break;
            case kDewP: ok = guessAtPressure(kDewP, a, b, s, realVapour); break;
            case kTP: ok = guessTP(a, b, s, realVapour); break;
        }
        if (!ok) continue;
        anyGuess = true;
        specify(s);
        if (guessT == 0.0) {
            guessT = s.T();
            guessP = s.p();
        }
        if (satNewton(s, kind)) {
            cache.valid = true;
            cache.s = s;
            return SatStatus::Found;
        }
    }
    // near the critical locus Newton can slide from a guess to the trivial solution (rejected in satNewton): the last
    // start is the equilibrium pair at the T and p of the guess, a real pair near the one sought
    if (!pairTP && guessT > 0.0) {
        Sat pair;
        if (satAttempt(kTP, guessT, guessP, 0.0, pair) == SatStatus::Found) {
            s = pair;
            specify(s);
            if (satNewton(s, kind)) {
                cache.valid = true;
                cache.s = s;
                return SatStatus::Found;
            }
        }
    }
    return anyGuess ? SatStatus::NotConverged : SatStatus::NoGuess;
}

/// The error of a saturation solve that did not give a state. Only `NoGuess` outside a known two-phase region may say
/// that the state is outside it; `insideDome` (the state is between a bubble and a dew point that were found) says
/// that any failure is one of the solver.
[[noreturn]] void satFailure(SatKind kind, double a, double b, SatStatus status, bool insideDome = false) {
    const bool pressureSpec = kind == kBubbleP || kind == kDewP;
    const std::string first = pressureSpec ? "p = " + fmt(a) + " Pa" : "T = " + fmt(a) + " K";
    const std::string second = kind == kTP ? "p = " + fmt(b) + " Pa" : "x(mole) = " + fmt(b);
    const std::string where = " (" + first + ", " + second + ")";
    if (insideDome)
        fail(std::string("the ") + kKindName[kind] + where + " did not converge, although the state is between its bubble "
             "and dew points");
    if (status == SatStatus::NoGuess)
        fail(std::string("no ") + kKindName[kind] + " found" + where +
             ": the state is probably outside the two-phase region or above the critical locus");
    fail(std::string("the ") + kKindName[kind] + where + " did not converge to a pair of stable phases");
}

/// The saturation state of `kind` (an error naming the specification if there is none); `bc` as for `satAttempt`.
Sat satSolve(SatKind kind, double a, double b, double bc) {
    Sat s;
    const SatStatus status = satAttempt(kind, a, b, bc, s);
    if (status != SatStatus::Found) satFailure(kind, a, b, status);
    return s;
}

/// A start for the equilibrium pair inside a bracket of two saturated states `lo` and `hi` (the dew and bubble points of
/// a composition, or two points of the saturation curve): the six variables interpolated linearly at the position of
/// `value` between lo.v[k] and hi.v[k] (k = 0 for T, 1 for p).
Sat between(const Sat& lo, const Sat& hi, int k, double value) {
    const double span = hi.v[k] - lo.v[k];
    const double t = span != 0.0 ? std::min(std::max((value - lo.v[k]) / span, 0.0), 1.0) : 0.5;
    Sat s;
    for (int j = 0; j < 6; ++j) s.v[j] = lo.v[j] + t * (hi.v[j] - lo.v[j]);
    for (int j = 0; j < 2; ++j) s.w[j] = lo.w[j] + t * (hi.w[j] - lo.w[j]);
    return s;
}

/// The equilibrium pair at (T, p) of a state known to be inside the two-phase region, with a start `hint` from the
/// bracket: any failure is an error of the solver, never a claim that the state is outside the region.
Sat satInside(double T, double p, const Sat& hint) {
    Sat s;
    const SatStatus status = satAttempt(kTP, T, p, 0.0, s, &hint);
    if (status != SatStatus::Found) satFailure(kTP, T, p, status, true);
    return s;
}

/// True if the saturation of `kind` was found (the scans, where a failed node is skipped); `bc` as for `satAttempt`.
bool trySat(SatKind kind, double a, double b, double bc, Sat& out, const Sat* hint = nullptr) {
    return satAttempt(kind, a, b, bc, out, hint) == SatStatus::Found;
}

// ---------------------------------------------------------------------------------------------------------------
// 4. The flashes
// ---------------------------------------------------------------------------------------------------------------

/// Relative margin at the phase boundaries: an input that is saturated to round-off is single-phase.
constexpr double kBoundaryTol = 1e-7;

struct Spec {
    Prop prop;
    double value;
};

double specValue(const Spec spec[3], Prop prop, bool* found = nullptr) {
    for (int i = 0; i < 3; ++i)
        if (spec[i].prop == prop) {
            if (found) *found = true;
            return spec[i].value;
        }
    if (found) *found = false;
    return 0.0;
}

/// The converged state: two-phase w = (T, p, xL, xV, rhoL, rhoV, xi, q) with the complements wc = (1 - xL, 1 - xV) of
/// its compositions (as `Sat::w`), or single-phase (T, rho, x) with the phase flag.
struct Solution {
    bool twoPhase = false;
    double w[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    double wc[2] = {1.0, 1.0};
    double s[3] = {0, 0, 0};
    bool vaporLike = false;
};

double mixProp(Prop prop, const Phase& L, const Phase& V, double q) {
    auto pick = [&](const J1& a, const J1& b) { return (1.0 - q) * a.value() + q * b.value(); };
    switch (prop) {
        case Prop::H: return pick(L.h, V.h);
        case Prop::S: return pick(L.s, V.s);
        case Prop::U: return pick(L.u, V.u);
        case Prop::V: return pick(L.v, V.v);
        default: fail("internal error: not a mixture property");
    }
}

/// Quality (vapour mass fraction) of the overall mass fraction xi at the saturation pair.
double qualityOf(double xi, const Sat& s) {
    const double xiL = massFractionOf(s.xL()), xiV = massFractionOf(s.xV());
    return (xi - xiL) / (xiV - xiL);
}

Solution twoPhaseSolution(const Sat& s, double q) {
    Solution sol;
    sol.twoPhase = true;
    for (int k = 0; k < 6; ++k) sol.w[k] = s.v[k];
    for (int k = 0; k < 2; ++k) sol.wc[k] = s.w[k];
    sol.w[6] = (1.0 - q) * massFractionOf(s.xL()) + q * massFractionOf(s.xV());
    sol.w[7] = q;
    return sol;
}

Solution singleSolution(double T, double rho, double x, bool vaporLike) {
    Solution sol;
    sol.s[0] = T;
    sol.s[1] = rho;
    sol.s[2] = x;
    sol.vaporLike = vaporLike;
    return sol;
}

// ----- the value of a property row on a state and its gradient --------------------------------------------------

/// Value of `prop` on a two-phase state and its gradient with respect to w = (T, p, xL, xV, rhoL, rhoV, xi, q).
double twoPhaseRow(Prop prop, const double* w, const Phase& L, const Phase& V, double* g) {
    for (int k = 0; k < 8; ++k) g[k] = 0.0;
    switch (prop) {
        case Prop::T: g[0] = 1.0; return w[0];
        case Prop::P: g[1] = 1.0; return w[1];
        case Prop::X: g[6] = 1.0; return w[6];
        case Prop::Q: g[7] = 1.0; return w[7];
        default: break;
    }
    const double q = w[7];
    const J1 &a = prop == Prop::H ? L.h : prop == Prop::S ? L.s : prop == Prop::U ? L.u : L.v;
    const J1 &b = prop == Prop::H ? V.h : prop == Prop::S ? V.s : prop == Prop::U ? V.u : V.v;
    g[0] = (1.0 - q) * a.gradient(0) + q * b.gradient(0);
    g[2] = (1.0 - q) * a.gradient(2);
    g[3] = q * b.gradient(2);
    g[4] = (1.0 - q) * a.gradient(1);
    g[5] = q * b.gradient(1);
    g[7] = b.value() - a.value();
    return (1.0 - q) * a.value() + q * b.value();
}

/// Value of `prop` on a single-phase state (T, rho, x) and its gradient with respect to it.
double singleRow(Prop prop, const double* s, const Phase& ph, double* g) {
    for (int k = 0; k < 3; ++k) g[k] = 0.0;
    auto copy = [&](const J1& a) {
        for (int k = 0; k < 3; ++k) g[k] = a.gradient(k);
        return a.value();
    };
    switch (prop) {
        case Prop::T: g[0] = 1.0; return s[0];
        case Prop::X: g[2] = dMassFractionDx(s[2]); return massFractionOf(s[2]);
        case Prop::P: return copy(ph.p);
        case Prop::H: return copy(ph.h);
        case Prop::S: return copy(ph.s);
        case Prop::U: return copy(ph.u);
        case Prop::V: return copy(ph.v);
        default: fail("internal error: the quality is not a single-phase input");
    }
}

double scaleOf(Prop prop, double value) {
    switch (prop) {
        case Prop::T: return value;
        case Prop::P: return value;
        case Prop::V: return std::abs(value);
        case Prop::H:
        case Prop::U: return std::max(std::abs(value), 1.0e4);
        case Prop::S: return std::max(std::abs(value), 1.0e2);
        default: return 1.0;
    }
}

FlashResult finishTwoPhase(Solution& sol, const Spec spec[3], const std::string& label) {
    double* w = sol.w;
    NewtonSetup ns;
    ns.n = 8;
    for (int i = 0; i < 8; ++i) {
        ns.rowScale[i] = 1.0;
        ns.colScale[i] = 1.0;
        ns.maxStep[i] = 0.1;
        ns.lower[i] = -1e300;
        ns.upper[i] = 1e300;
    }
    const double RT = kR * w[0];
    ns.rowScale[0] = std::max(w[1], w[4] * RT);
    ns.rowScale[1] = w[1];
    ns.rowScale[2] = ns.rowScale[3] = kR * 300.0;
    for (int j = 0; j < 3; ++j) ns.rowScale[5 + j] = scaleOf(spec[j].prop, spec[j].value);
    ns.colScale[0] = 100.0;
    ns.logVar[1] = ns.logVar[4] = ns.logVar[5] = true;
    ns.maxStep[1] = 0.3; ns.maxStep[4] = ns.maxStep[5] = 0.2;
    ns.maxStep[7] = 0.3; ns.maxStep[6] = 0.1; ns.maxStep[2] = ns.maxStep[3] = 1.0; ns.maxStep[0] = 0.02;
    ns.lower[0] = 120.0; ns.upper[0] = 1000.0; ns.lower[1] = 1e-3; ns.upper[1] = 2e9;
    ns.lower[2] = ns.lower[3] = -23.0; ns.upper[2] = ns.upper[3] = 23.0;
    ns.lower[4] = 1.0; ns.lower[5] = 1e-6; ns.upper[4] = ns.upper[5] = 2e5;   // the bounds of satNewton
    ns.lower[6] = 0.0; ns.upper[6] = 1.0; ns.lower[7] = -0.5; ns.upper[7] = 1.5;

    ns.composition[2] = ns.composition[3] = true;
    ns.resolution = [](const double* u, double* out) {
        for (int k = 0; k < 8; ++k) out[k] = std::abs(u[k]) * 2.220446049250313e-16;
        out[2] = logitResolution(u[2]);
        out[3] = logitResolution(u[3]);
    };
    // w holds the compositions as mole fractions (with their complements in sol.wc); the iteration vector u holds their
    // logits z, evaluated as x = sigmoid(z) and 1 - x = sigmoid(-z)
    auto toW = [](const double* u, double* v, double* c) {
        for (int k = 0; k < 8; ++k) v[k] = u[k];
        v[2] = sigmoid(u[2]);
        v[3] = sigmoid(u[3]);
        c[0] = sigmoid(-u[2]);
        c[1] = sigmoid(-u[3]);
    };
    auto assemble = [&](const double* v, const double* c, double* F, double J[8][8], Phase* L, Phase* V) {
        double Jf[4][6], Ff[4];
        satResiduals(v, c, Ff, Jf, L, V);
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j) J[i][j] = 0.0;
        for (int i = 0; i < 4; ++i) {
            F[i] = Ff[i];
            for (int j = 0; j < 6; ++j) J[i][j] = Jf[i][j];
        }
        const double xiL = massFractionOf(v[2]), xiV = massFractionOf(v[3]), q = v[7];
        F[4] = v[6] - (1.0 - q) * xiL - q * xiV;
        J[4][6] = 1.0;
        J[4][7] = -(xiV - xiL);
        J[4][2] = -(1.0 - q) * dMassFractionDx(v[2]);
        J[4][3] = -q * dMassFractionDx(v[3]);
        for (int j = 0; j < 3; ++j) {
            double g[8];
            F[5 + j] = twoPhaseRow(spec[j].prop, v, *L, *V, g) - spec[j].value;
            for (int k = 0; k < 8; ++k) J[5 + j][k] = g[k];
        }
    };
    const Evaluator eval = [&](const double* u, double* F, double (*J)[8]) {
        double v[8], c[2];
        toW(u, v, c);
        Phase L, V;
        try {
            assemble(v, c, F, J, &L, &V);
        } catch (const std::runtime_error&) {
            return false;
        }
        for (int i = 0; i < 8; ++i) {   // chain rule to the logit variables
            J[i][2] *= v[2] * c[0];
            J[i][3] *= v[3] * c[1];
        }
        return true;
    };
    double u[8];
    for (int k = 0; k < 8; ++k) u[k] = w[k];
    u[2] = logitOf(w[2], sol.wc[0]);
    u[3] = logitOf(w[3], sol.wc[1]);
    if (!newtonSolve(ns, u, eval, 1e-12, 30)) fail("the two-phase state of " + label + " did not converge");
    toW(u, w, sol.wc);
    if (!(w[7] >= -1e-9 && w[7] <= 1.0 + 1e-9)) fail("the state of " + label + " is not two-phase (q = " + fmt(w[7]) + ")");

    FlashResult res;
    res.twoPhase = true;
    double F[8], J[8][8];
    Phase L, V;
    assemble(w, sol.wc, F, J, &L, &V);
    if (!mechanicallyStable(L, V))
        fail("the two-phase state of " + label + " has a mechanically unstable phase (dp/drho <= 0; liquid " +
             fmt(L.p.gradient(1)) + ", vapour " + fmt(V.p.gradient(1)) + " Pa m3/mol)");
    double B[8][8];
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j) B[i][j] = (i == 5 + j && j < 3) ? 1.0 : 0.0;
    if (!solveLinear(8, J, 3, B)) fail("singular Jacobian at the two-phase state of " + label);
    const Prop props[8] = {Prop::T, Prop::P, Prop::X, Prop::H, Prop::S, Prop::U, Prop::V, Prop::Q};
    for (int k = 0; k < 8; ++k) {
        double g[8];
        res.value[k] = twoPhaseRow(props[k], w, L, V, g);
        for (int j = 0; j < 3; ++j) {
            double d = 0.0;
            for (int m = 0; m < 8; ++m) d += g[m] * B[m][j];
            res.derivative[k][j] = d;
        }
    }
    return res;
}

/// Single-phase polish, outputs and derivatives; also returns the 3 x 3 matrix d(T, rho, x)/d(inputs).
FlashResult finishSingle(Solution& sol, const Spec spec[3], const std::string& label, double D[3][3] = nullptr) {
    double* s = sol.s;
    NewtonSetup ns;
    ns.n = 3;
    for (int j = 0; j < 3; ++j) {
        ns.rowScale[j] = scaleOf(spec[j].prop, spec[j].value);
        // the pressure of a liquid is a sum of terms of the size rho R T: its round-off is relative to that, not to p
        if (spec[j].prop == Prop::P) ns.rowScale[j] = std::max(spec[j].value, s[1] * kR * s[0]);
        ns.maxStep[j] = 0.1;
        ns.colScale[j] = 1.0;
    }
    ns.colScale[0] = 100.0;
    ns.logVar[1] = true;
    ns.maxStep[1] = 0.2; ns.maxStep[2] = 0.05;
    ns.lower[0] = 120.0; ns.upper[0] = 2000.0; ns.lower[1] = 1e-6; ns.upper[1] = 2e5; ns.lower[2] = 0.0; ns.upper[2] = 1.0;
    const Evaluator eval = [&](const double* v, double* F, double (*J)[8]) {
        Phase ph;
        try {
            ph = evaluatePhase(v[0], v[1], v[2]);
        } catch (const std::runtime_error&) {
            return false;
        }
        for (int j = 0; j < 3; ++j) {
            double g[3];
            F[j] = singleRow(spec[j].prop, v, ph, g) - spec[j].value;
            for (int k = 0; k < 3; ++k) J[j][k] = g[k];
        }
        return true;
    };
    if (!newtonSolve(ns, s, eval, 1e-12, 30)) fail("the single-phase state of " + label + " did not converge");
    FlashResult res;
    res.twoPhase = false;
    const Phase ph = evaluatePhase(s[0], s[1], s[2]);
    double A[8][8], B[8][8];
    for (int j = 0; j < 3; ++j) {
        double g[3];
        singleRow(spec[j].prop, s, ph, g);
        for (int k = 0; k < 3; ++k) A[j][k] = g[k];
    }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) B[i][j] = i == j ? 1.0 : 0.0;
    if (!solveLinear(3, A, 3, B)) fail("singular Jacobian at the single-phase state of " + label);
    if (D)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) D[i][j] = B[i][j];
    const Prop props[7] = {Prop::T, Prop::P, Prop::X, Prop::H, Prop::S, Prop::U, Prop::V};
    for (int k = 0; k < 7; ++k) {
        double g[3];
        res.value[k] = singleRow(props[k], s, ph, g);
        for (int j = 0; j < 3; ++j) res.derivative[k][j] = g[0] * B[0][j] + g[1] * B[1][j] + g[2] * B[2][j];
    }
    res.value[7] = sol.vaporLike ? 1.001 : -0.001;
    for (int j = 0; j < 3; ++j) res.derivative[7][j] = 0.0;
    return res;
}

// ----- the value solvers -----------------------------------------------------------------------------------------

/// Single-phase state at (T, p, x): the liquid (`liquid`) or vapour root.
Solution singleAtTP(double T, double p, double x, bool liquid, const std::string& label) {
    double rho;
    if (!densityRoot(T, p, x, liquid, rho))
        fail("no " + std::string(liquid ? "liquid" : "vapour") + " density found for " + label);
    return singleSolution(T, rho, x, !liquid);
}

/// Phase of the single-phase state at (T, p, x) when there is no two-phase region to compare with (above the
/// critical locus): the root of lowest Gibbs energy.
Solution singleSupercritical(double T, double p, double x, const std::string& label) {
    double rl = 0, rv = 0;
    const bool hasL = densityRoot(T, p, x, true, rl), hasV = densityRoot(T, p, x, false, rv);
    if (!hasL && !hasV) fail("no density found for " + label);
    double rho = hasL ? rl : rv;
    if (hasL && hasV && std::abs(rl / rv - 1.0) > 1e-6) {
        const Phase L = evaluatePhase(T, rl, x), V = evaluatePhase(T, rv, x);
        const double gL = (1.0 - x) * L.mu1.value() + x * L.mu2.value(), gV = (1.0 - x) * V.mu1.value() + x * V.mu2.value();
        rho = gL <= gV ? rl : rv;
    }
    // without a two-phase region to compare with, the label (q = -0.001 or 1.001) follows the reducing density
    return singleSolution(T, rho, x, rho < reducingDensity(x));
}

// ---------------------------------------------------------------------------------------------------------------
// flash
// ---------------------------------------------------------------------------------------------------------------

std::string propName(Prop p) {
    static const char* const names[8] = {"T", "P", "x", "h", "s", "u", "v", "q"};
    return names[static_cast<int>(p) - 1];
}

std::string describe(const Spec spec[3]) {
    std::string out = "NH3H2O(";
    for (int i = 0; i < 3; ++i) out += (i ? ", " : "") + propName(spec[i].prop) + " = " + fmt(spec[i].value);
    return out + ")";
}

/// The flashes with a given quality and an unknown overall composition: (T, h|s|u|v, q) with the pressure unknown and
/// (P, h|s|u|v, q) with the temperature unknown. The saturated states are scanned along the liquid composition (a
/// logit grid from 1e-4 to 1 - 1e-4) for a sign change of `prop` - `target` at the quality `q`, then the root is
/// refined with the equilibrium pair at (T, p). The property of a mixture is not monotone in the composition (heat
/// of mixing), so there can be two solutions; the one of the lowest ammonia content is returned.
Solution solveUnknownComposition(bool temperatureGiven, double fixedValue, Prop prop, double target, double q,
                                 const std::string& label) {
    auto at = [&](const Sat& s) {
        const Phase L = evaluatePhase(s.T(), s.rhoL(), s.xL(), s.wL()), V = evaluatePhase(s.T(), s.rhoV(), s.xV(), s.wV());
        return mixProp(prop, L, V, q) - target;
    };
    const SatKind kind = temperatureGiven ? kBubbleT : kBubbleP;
    struct Node {
        double z = 0.0, f = 0.0;
        Sat sat;
        bool ok = false;
    };
    auto evalNode = [&](double z) {
        Node n;
        n.z = z;
        n.ok = trySat(kind, fixedValue, sigmoid(z), sigmoid(-z), n.sat);
        if (n.ok) n.f = at(n.sat);
        return n;
    };
    // refine a sign change between two saturated states with the equilibrium pair at (T, p)
    auto refine = [&](const Node& A, const Node& B, Solution& out) {
        // the pressure rises with x at fixed T, the temperature falls at fixed p
        const double lo = temperatureGiven ? A.sat.p() : B.sat.T(), hi = temperatureGiven ? B.sat.p() : A.sat.T();
        const double flo = temperatureGiven ? A.f : B.f, fhi = temperatureGiven ? B.f : A.f;
        const int k = temperatureGiven ? 1 : 0;   // the coordinate that varies between the two nodes
        auto f = [&](double v) {
            Sat pair;
            const Sat hint = between(A.sat, B.sat, k, v);
            if (!trySat(kTP, temperatureGiven ? fixedValue : v, temperatureGiven ? v : fixedValue, 0.0, pair, &hint))
                return std::numeric_limits<double>::quiet_NaN();
            return at(pair);
        };
        double root;
        if (!brent(f, lo, hi, flo, fhi, 1e-13 * hi, root)) return false;
        const Sat pair =
            satInside(temperatureGiven ? fixedValue : root, temperatureGiven ? root : fixedValue, between(A.sat, B.sat, k, root));
        out = twoPhaseSolution(pair, q);
        return true;
    };
    constexpr int kPoints = 25;
    std::vector<Node> nodes;
    for (int k = 0; k < kPoints; ++k) nodes.push_back(evalNode(-9.0 + 18.0 * k / (kPoints - 1)));
    Solution out;
    for (int k = 0; k + 1 < kPoints; ++k)
        if (nodes[k].ok && nodes[k + 1].ok && (nodes[k].f == 0.0 || nodes[k].f * nodes[k + 1].f <= 0.0) &&
            refine(nodes[k], nodes[k + 1], out))
            return out;
    // no sign change at the nodes: the property may dip through the target between two of them (a double root of
    // the heat of mixing): golden-section search for the extremum around the node closest to the target
    int best = -1;
    for (int k = 0; k < kPoints; ++k)
        if (nodes[k].ok && (best < 0 || std::abs(nodes[k].f) < std::abs(nodes[best].f))) best = k;
    if (best >= 0 && best > 0 && best + 1 < kPoints && nodes[best - 1].ok && nodes[best + 1].ok) {
        const double sign = nodes[best].f > 0.0 ? 1.0 : -1.0;   // minimise sign * f
        double a = nodes[best - 1].z, b = nodes[best + 1].z;
        const double gr = 0.3819660112501051;
        double c1 = a + gr * (b - a), c2 = b - gr * (b - a);
        Node n1 = evalNode(c1), n2 = evalNode(c2);
        for (int it = 0; it < 40 && n1.ok && n2.ok; ++it) {
            for (const Node* n : {&n1, &n2})
                if (sign * n->f < 0.0) {
                    // roots on both sides of this node: take the one of lower ammonia content
                    if (refine(nodes[best - 1], *n, out)) return out;
                    if (refine(*n, nodes[best + 1], out)) return out;
                }
            if (sign * n1.f < sign * n2.f) {
                b = c2; c2 = c1; n2 = n1; c1 = a + gr * (b - a); n1 = evalNode(c1);
            } else {
                a = c1; c1 = c2; n1 = n2; c2 = b - gr * (b - a); n2 = evalNode(c2);
            }
        }
    }
    fail("no state for " + label + ": the value is outside the range at this " +
         (temperatureGiven ? "temperature" : "pressure") + " and quality");
}

/// The bubble and the dew point of the mole fraction x (complement xc = 1 - x, see `satAttempt`) at the temperature
/// (`temperatureGiven`) or the pressure `a`, with their status. A pure component (x or xc at the composition floor
/// kXMin) has no two-phase region of non-zero width, and the
/// mixture cannot represent its saturation pair from both sides (one phase would have to be purer than the floor):
/// nothing is solved and both statuses are `NoGuess`, so that the flashes use the single-phase rules of a state without
/// a two-phase region (the density root of lowest Gibbs energy for T, p, x).
struct Boundaries {
    Sat bubble, dew;
    SatKind bubbleKind = kBubbleT, dewKind = kDewT;
    SatStatus bubbleStatus = SatStatus::NoGuess, dewStatus = SatStatus::NoGuess;
    bool notConverged() const { return bubbleStatus == SatStatus::NotConverged || dewStatus == SatStatus::NotConverged; }
};

Boundaries boundariesOf(bool temperatureGiven, double a, double x, double xc) {
    Boundaries bd;
    bd.bubbleKind = temperatureGiven ? kBubbleT : kBubbleP;
    bd.dewKind = temperatureGiven ? kDewT : kDewP;
    if (x <= kXMin || xc <= kXMin) return bd;   // a pure component
    bd.bubbleStatus = satAttempt(bd.bubbleKind, a, x, xc, bd.bubble);
    bd.dewStatus = satAttempt(bd.dewKind, a, x, xc, bd.dew);
    return bd;
}

/// The error of a state whose phase is unknown because a saturation it needs did not converge.
[[noreturn]] void failUndecided(const Boundaries& bd, double a, double x) {
    if (bd.bubbleStatus == SatStatus::NotConverged) satFailure(bd.bubbleKind, a, x, bd.bubbleStatus);
    satFailure(bd.dewKind, a, x, bd.dewStatus);
}

Solution solveFlash(const Spec spec[3], const std::string& label) {
    const int code = 100 * static_cast<int>(spec[0].prop) + 10 * static_cast<int>(spec[1].prop) + static_cast<int>(spec[2].prop);
    const double a = spec[0].value, b = spec[1].value, c = spec[2].value;
    auto inRange = [&](double q) { return q >= 0.0 && q <= 1.0; };

    switch (code) {
        case 123: {   // T, p, x
            const double T = a, p = b, x = moleFractionOf(c), xc = moleFractionComplementOf(c);
            const Boundaries bd = boundariesOf(true, T, x, xc);
            const Sat &sb = bd.bubble, &sd = bd.dew;
            const bool hasB = bd.bubbleStatus == SatStatus::Found, hasD = bd.dewStatus == SatStatus::Found;
            if (hasB != hasD) {
                // one saturation point only: above the critical temperature of x (two dew points, the retrograde
                // region), near the critical locus where a guess finds no liquid or Newton does not converge, or near a
                // pure component. The equilibrium pair at (T, p), started from the point that was found, decides:
                // two-phase if its compositions bracket x (not within the boundary margin of that point)
                const Sat& known = hasB ? sb : sd;
                Sat s;
                if (std::abs(p / known.p() - 1.0) > kBoundaryTol && satAttempt(kTP, T, p, 0.0, s, &known) == SatStatus::Found &&
                    lessX(s.xL(), s.wL(), x, xc) && lessX(x, xc, s.xV(), s.wV()))
                    return twoPhaseSolution(s, qualityOf(c, s));
            }
            if (hasB && p >= sb.p() * (1.0 - kBoundaryTol)) return singleAtTP(T, p, x, true, label);
            if (hasD && p <= sd.p() * (1.0 + kBoundaryTol)) return singleAtTP(T, p, x, false, label);
            if (hasB && hasD) {
                const Sat s = satInside(T, p, between(sd, sb, 1, p));
                return twoPhaseSolution(s, qualityOf(c, s));
            }
            // a single phase only on the evidence of the saturations: found, or no guess at all (above the critical
            // locus); a saturation that did not converge leaves the phase unknown
            if (bd.notConverged()) failUndecided(bd, T, x);
            return singleSupercritical(T, p, x, label);
        }
        case 128: {   // T, p, q
            if (!inRange(c)) fail("the quality must be in [0, 1] (" + label + ")");
            return twoPhaseSolution(satSolve(kTP, a, b, 0.0), c);
        }
        case 138: {   // T, x, q
            if (!inRange(c)) fail("the quality must be in [0, 1] (" + label + ")");
            const double x = moleFractionOf(b), xc = moleFractionComplementOf(b);
            if (c == 0.0) return twoPhaseSolution(satSolve(kBubbleT, a, x, xc), 0.0);
            if (c == 1.0) return twoPhaseSolution(satSolve(kDewT, a, x, xc), 1.0);
            const Sat sb = satSolve(kBubbleT, a, x, xc), sd = satSolve(kDewT, a, x, xc);
            double p;
            auto f = [&](double pp) { return qualityOf(b, satInside(a, pp, between(sd, sb, 1, pp))) - c; };
            if (!brent(f, sd.p(), sb.p(), 1.0 - c, -c, 1e-12 * sb.p(), p)) fail("no state for " + label);
            return twoPhaseSolution(satInside(a, p, between(sd, sb, 1, p)), c);
        }
        case 137: {   // T, x, v
            const double T = a, x = moleFractionOf(b), xc = moleFractionComplementOf(b), rho = 1.0 / (c * molarMass(x));
            const Boundaries bd = boundariesOf(true, T, x, xc);
            const Sat &sb = bd.bubble, &sd = bd.dew;
            const bool hasB = bd.bubbleStatus == SatStatus::Found, hasD = bd.dewStatus == SatStatus::Found;
            if (!hasB || !hasD) {
                if (!bd.notConverged())   // above the critical locus
                    return singleSolution(T, rho, x, rho < reducingDensity(x));
                // a saturation did not converge: a single phase only on the evidence of the one that was found
                if (hasB && rho >= sb.rhoL() * (1.0 - kBoundaryTol)) return singleSolution(T, rho, x, false);
                if (hasD && rho <= sd.rhoV() * (1.0 + kBoundaryTol)) return singleSolution(T, rho, x, true);
                failUndecided(bd, T, x);
            }
            if (rho >= sb.rhoL() * (1.0 - kBoundaryTol)) return singleSolution(T, rho, x, false);
            if (rho <= sd.rhoV() * (1.0 + kBoundaryTol)) return singleSolution(T, rho, x, true);
            auto vmix = [&](double pp) {
                const Sat s = satInside(T, pp, between(sd, sb, 1, pp));
                const double q = qualityOf(b, s);
                return (1.0 - q) / (s.rhoL() * molarMass(s.xL())) + q / (s.rhoV() * molarMass(s.xV())) - c;
            };
            const double flo = vmix(sd.p()), fhi = vmix(sb.p());
            double p;
            if (!brent(vmix, sd.p(), sb.p(), flo, fhi, 1e-12 * sb.p(), p)) fail("no two-phase state for " + label);
            const Sat s = satInside(T, p, between(sd, sb, 1, p));
            return twoPhaseSolution(s, qualityOf(b, s));
        }
        case 148: case 158: case 168: case 178:   // T, h|s|u|v, q: the pressure is unknown
            if (!inRange(c)) fail("the quality must be in [0, 1] (" + label + ")");
            return solveUnknownComposition(true, a, spec[1].prop, b, c, label);
        case 234: case 235: {   // P, x, h|s
            const Prop prop = spec[2].prop;
            const double p = a, x = moleFractionOf(b), xc = moleFractionComplementOf(b);
            const Boundaries bd = boundariesOf(false, p, x, xc);
            const Sat &sb = bd.bubble, &sd = bd.dew;
            const bool hasB = bd.bubbleStatus == SatStatus::Found, hasD = bd.dewStatus == SatStatus::Found;
            auto subcooled = [&]() {   // the liquid below the bubble temperature
                auto f = [&](double T) {
                    double rho;
                    if (!densityRoot(T, p, x, true, rho)) return std::numeric_limits<double>::quiet_NaN();
                    const Phase ph = evaluatePhase(T, rho, x);
                    return mixProp(prop, ph, ph, 0.0) - c;
                };
                double T, Tlo = 180.0, flo = std::numeric_limits<double>::quiet_NaN();
                for (; Tlo < sb.T() && !(flo < 0.0); Tlo += 10.0) flo = f(Tlo);   // the liquid root needs a meaningful state
                if (!brent(f, Tlo - 10.0, sb.T(), flo, f(sb.T()), 1e-10, T)) fail("no liquid state for " + label + ": the value is below the range");
                return singleAtTP(T, p, x, true, label);
            };
            auto superheated = [&]() {   // the vapour above the dew temperature
                auto f = [&](double T) {
                    double rho;
                    if (!densityRoot(T, p, x, false, rho)) return std::numeric_limits<double>::quiet_NaN();
                    const Phase ph = evaluatePhase(T, rho, x);
                    return mixProp(prop, ph, ph, 0.0) - c;
                };
                double T;
                const double Thi = 1500.0;
                if (!brent(f, sd.T(), Thi, f(sd.T()), f(Thi), 1e-10, T)) fail("no vapour state for " + label + ": the value is above the range");
                return singleAtTP(T, p, x, false, label);
            };
            const double nan = std::numeric_limits<double>::quiet_NaN();
            double hb = nan, hd = nan;
            if (hasB) {
                const Phase Lb = evaluatePhase(sb.T(), sb.rhoL(), sb.xL(), sb.wL());
                hb = mixProp(prop, Lb, Lb, 0.0);
            }
            if (hasD) {
                const Phase Vd = evaluatePhase(sd.T(), sd.rhoV(), sd.xV(), sd.wV());
                hd = mixProp(prop, Vd, Vd, 0.0);
            }
            if (!hasB || !hasD) {
                // one saturation point at most: a single phase on its side of the one that was found
                if (hasB && c <= hb + kBoundaryTol * std::abs(hb))
                    return c >= hb - kBoundaryTol * std::abs(hb) ? singleAtTP(sb.T(), p, x, true, label) : subcooled();
                if (hasD && c >= hd - kBoundaryTol * std::abs(hd))
                    return c <= hd + kBoundaryTol * std::abs(hd) ? singleAtTP(sd.T(), p, x, false, label) : superheated();
                // otherwise a single phase only if no saturation failed: no guess at all is the evidence of a state
                // above the critical locus (single phase at any T; search T on the root of lowest Gibbs energy)
                if (bd.notConverged()) failUndecided(bd, p, x);
                auto f = [&](double T) {
                    Solution s = singleSupercritical(T, p, x, label);
                    return mixProp(prop, evaluatePhase(T, s.s[1], x), evaluatePhase(T, s.s[1], x), 0.0) - c;
                };
                double T;
                if (!brent(f, 200.0, 1500.0, f(200.0), f(1500.0), 1e-10, T)) fail("no state for " + label);
                return singleSupercritical(T, p, x, label);
            }
            const double hTol = kBoundaryTol * std::max(std::abs(hb), std::abs(hd));
            if (std::abs(c - hb) <= hTol) return singleAtTP(sb.T(), p, x, true, label);   // saturated liquid
            if (std::abs(c - hd) <= hTol) return singleAtTP(sd.T(), p, x, false, label);   // saturated vapour
            if (c <= hb) return subcooled();
            if (c >= hd) return superheated();
            auto f = [&](double T) {
                const Sat s = satInside(T, p, between(sb, sd, 0, T));
                const Phase L = evaluatePhase(T, s.rhoL(), s.xL(), s.wL()), V = evaluatePhase(T, s.rhoV(), s.xV(), s.wV());
                return mixProp(prop, L, V, qualityOf(b, s)) - c;
            };
            double T;
            if (!brent(f, sb.T(), sd.T(), hb - c, hd - c, 1e-12 * sd.T(), T)) fail("no two-phase state for " + label);
            const Sat s = satInside(T, p, between(sb, sd, 0, T));
            return twoPhaseSolution(s, qualityOf(b, s));
        }
        case 238: {   // P, x, q
            if (!inRange(c)) fail("the quality must be in [0, 1] (" + label + ")");
            const double p = a, x = moleFractionOf(b), xc = moleFractionComplementOf(b);
            if (c == 0.0) return twoPhaseSolution(satSolve(kBubbleP, p, x, xc), 0.0);
            if (c == 1.0) return twoPhaseSolution(satSolve(kDewP, p, x, xc), 1.0);
            const Sat sb = satSolve(kBubbleP, p, x, xc), sd = satSolve(kDewP, p, x, xc);
            auto f = [&](double T) { return qualityOf(b, satInside(T, p, between(sb, sd, 0, T))) - c; };
            double T;
            if (!brent(f, sb.T(), sd.T(), -c, 1.0 - c, 1e-12 * sd.T(), T)) fail("no state for " + label);
            return twoPhaseSolution(satInside(T, p, between(sb, sd, 0, T)), c);
        }
        case 248: case 258: case 268: case 278:   // P, h|s|u|v, q: the temperature is unknown
            if (!inRange(c)) fail("the quality must be in [0, 1] (" + label + ")");
            return solveUnknownComposition(false, a, spec[1].prop, b, c, label);
        default: fail("unsupported combination of properties: " + label);
    }
}

void checkInputs(const Spec spec[3], const std::string& label) {
    for (int i = 0; i < 3; ++i) {
        const double v = spec[i].value;
        if (!std::isfinite(v)) fail("an input is not finite (" + label + ")");
        switch (spec[i].prop) {
            case Prop::T: if (!(v > 0.0)) fail("the temperature must be positive (" + label + ")"); break;
            case Prop::P: if (!(v > 0.0)) fail("the pressure must be positive (" + label + ")"); break;
            case Prop::V: if (!(v > 0.0)) fail("the specific volume must be positive (" + label + ")"); break;
            case Prop::X: if (!(v >= 0.0 && v <= 1.0)) fail("the mass fraction must be in [0, 1] (" + label + ")"); break;
            case Prop::Q: if (!(v >= 0.0 && v <= 1.0)) fail("the quality must be in [0, 1] (" + label + ")"); break;
            default: break;
        }
    }
}

}  // namespace

Saturation saturation(SatSpec spec, double a, double b) {
    // a given mole fraction b is the composition itself: 1 - b is its exact complement (exact for b >= 0.5)
    const Sat s = satSolve(static_cast<SatKind>(spec), a, b, spec == SatSpec::TP ? 0.0 : 1.0 - b);
    return {s.T(), s.p(), s.xL(), s.xV(), s.rhoL(), s.rhoV(), s.wL(), s.wV()};
}

FlashResult flash(Prop first, Prop second, Prop third, double a, double b, double c) {
    const Spec spec[3] = {{first, a}, {second, b}, {third, c}};
    const std::string label = describe(spec);
    if (!(static_cast<int>(first) < static_cast<int>(second) && static_cast<int>(second) < static_cast<int>(third)))
        fail("the three properties must be distinct and in the order T, P, x, h, s, u, v, q (" + label + ")");
    checkInputs(spec, label);
    try {
        Solution sol = solveFlash(spec, label);
        return sol.twoPhase ? finishTwoPhase(sol, spec, label) : finishSingle(sol, spec, label);
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        if (message.find(label) != std::string::npos) throw;
        throw std::runtime_error(message + " [state: " + label + "]");
    }
}

SinglePhaseExtras singlePhaseExtras(Prop first, Prop second, Prop third, double a, double b, double c) {
    const Spec spec[3] = {{first, a}, {second, b}, {third, c}};
    const std::string label = describe(spec);
    checkInputs(spec, label);
    Solution sol;
    try {
        sol = solveFlash(spec, label);
    } catch (const std::runtime_error& e) {
        const std::string message = e.what();
        if (message.find(label) != std::string::npos) throw;
        throw std::runtime_error(message + " [state: " + label + "]");
    }
    if (sol.twoPhase) fail("the specific heat and the expansion coefficient are defined for single-phase states only (" + label + ")");
    double D[3][3];
    finishSingle(sol, spec, label, D);
    const double T = sol.s[0], rho = sol.s[1], x = clampX(sol.s[2]);
    const Jet<3> jT = Jet<3>::variable(T, 0), jR = Jet<3>::variable(rho, 1), jX = Jet<3>::variable(x, 2);
    const Jet<3> f3 = kR * jT * reducedHelmholtz(jT, jR, jX);
    const Jet<2> T2(jT), R2(jR), f2(f3), fT = derivative(f3, 0), fR = derivative(f3, 1);
    const Jet<2> p2 = R2 * R2 * fR;
    const Jet<2> h2 = f2 - T2 * fT + p2 / R2;
    const J1 hT = derivative(h2, 0), hR = derivative(h2, 1), pT = derivative(p2, 0), pR = derivative(p2, 1);
    const J1 x1 = J1::variable(x, 2), r1 = J1::variable(rho, 1);
    const J1 cp = (hT - hR * pT / pR) / ((1.0 - x1) * kM1 + x1 * kM2);
    const J1 beta = pT / (pR * r1);
    SinglePhaseExtras out;
    out.cp = cp.value();
    out.expansion = beta.value();
    for (int j = 0; j < 3; ++j) {
        out.dcp[j] = out.dexpansion[j] = 0.0;
        for (int k = 0; k < 3; ++k) {
            out.dcp[j] += cp.gradient(k) * D[k][j];
            out.dexpansion[j] += beta.gradient(k) * D[k][j];
        }
    }
    return out;
}

}  // namespace coolsolve::nh3h2o
