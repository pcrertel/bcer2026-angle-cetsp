// ==========================================================================
// max_angle_chain.hpp
//
// Closed-form block update of the BCD (single-disk subproblem of the paper):
// chain H -> P_i -> b -> P_k -> L with b free in the disk D_j = (c_j, r_j);
// maximize f(b) = theta_i + theta_j + theta_k, equivalently minimize the
// total turning 3*pi - f. Closed form: when the segment P_iP_k meets D_j the
// straight-passage witness is returned; otherwise the optimum is attained at
// one of the at most eight candidates of the theorem "Closed-form solution of
// the single-block subproblem" of the paper (all from roots of quadratic
// equations), and the best of them is returned.
//
// Implementation notes:
//   * Header-only, no heap allocation, no exceptions, POD types only,
//     fixed-size stack buffers (<= 20 raw candidates).
//   * solve_chain() is stateless and reentrant, so independent solves can
//     run concurrently.
//   * Threshold tests use squared distances throughout (no square roots
//     taken for comparisons). Build with -O3 -march=native -fno-math-errno.
//
// Usage:
//   #include "max_angle_chain.hpp"
//   maxangle::ChainResult res;
//   bool ok = maxangle::solve_chain(H, Pi, Pk, L, cj, rj, res);
//
// ==========================================================================
#pragma once

#include <cmath>
#include <cstdint>

// Optional fast atan2: about 88% of the time of one solve is spent in atan2, and
// the glibc atan2 is slow on Linux/Intel. Build with -DMAXANGLE_USE_SLEEF and
// link -lsleef to route the atan2 calls through SLEEF's 1-ulp scalar kernel
// Sleef_atan2d1_u10 (see build.sh and sleef_compat.hpp). The results agree to
// rounding, but they are not bit-identical to std::atan2, so a run of the
// ALNS-BCD may follow a different trajectory. The runs of the paper used the
// default, std::atan2.
#if defined(MAXANGLE_USE_SLEEF)
#include <sleef.h>
#define MAXANGLE_ATAN2(y, x) Sleef_atan2d1_u10((y), (x))
#else
#define MAXANGLE_ATAN2(y, x) std::atan2((y), (x))
#endif

namespace maxangle {

inline constexpr double kPi    = 3.14159265358979323846;
inline constexpr double kTwoPi = 6.28318530717958647692;

// ------------------------------------------------------------------ Vec2 ---

struct Vec2 { double x, y; };

inline constexpr Vec2   operator+(Vec2 a, Vec2 b) noexcept { return {a.x + b.x, a.y + b.y}; }
inline constexpr Vec2   operator-(Vec2 a, Vec2 b) noexcept { return {a.x - b.x, a.y - b.y}; }
inline constexpr Vec2   operator*(double s, Vec2 a) noexcept { return {s * a.x, s * a.y}; }
inline constexpr double dot  (Vec2 a, Vec2 b) noexcept { return a.x * b.x + a.y * b.y; }
inline constexpr double cross(Vec2 a, Vec2 b) noexcept { return a.x * b.y - a.y * b.x; }
inline constexpr Vec2   perp (Vec2 a)         noexcept { return {-a.y, a.x}; }
// sqrt(dot) instead of std::hypot: ~4x faster, safe at geometric scales.
inline double norm(Vec2 a) noexcept { return std::sqrt(dot(a, a)); }
inline double dist(Vec2 a, Vec2 b) noexcept { return norm(a - b); }
// squared distance: all threshold tests compare against squared bounds,
// which is the same predicate without the square root
inline constexpr double dist2(Vec2 a, Vec2 b) noexcept {
    return dot(a - b, a - b);
}
inline double heading(Vec2 v) noexcept { return MAXANGLE_ATAN2(v.y, v.x); }

// wrap to [-pi, pi), matching the Python reference `_wrap`
inline double wrap_pi(double a) noexcept {
    double w = std::fmod(a + kPi, kTwoPi);
    if (w < 0.0) w += kTwoPi;
    return w - kPi;
}

// angle A-P-B at the vertex P, in [0, pi]; atan2 form, well conditioned
// near 0 and pi (no arccos cancellation)
inline double angle_at(Vec2 P, Vec2 A, Vec2 B) noexcept {
    const Vec2 u = A - P, v = B - P;
    return MAXANGLE_ATAN2(std::fabs(cross(u, v)), dot(u, v));
}

// f(b) = theta_i + theta_j + theta_k  (= 3*pi - sum |tau_v|)
inline double theta_sum(Vec2 b, Vec2 H, Vec2 Pi, Vec2 Pk, Vec2 L) noexcept {
    return angle_at(Pi, H, b) + angle_at(b, Pi, Pk) + angle_at(Pk, b, L);
}

// ---------------------------------------------------------------- result ---

enum class ChainStatus : std::uint8_t {
    Point,       // bound not attained; maximizer is a closed-form candidate
    Plateau,     // bound 3*pi - |Delta| attained on a region; point = witness
    PointDisk,   // r_j = 0: the disk is the single point c_j
    Degenerate,  // P_i or P_k in D_j; point is a feasible lower bound only
    Straight     // the open segment (P_i,P_k) meets D_j: psi is globally
                 // minimised on that segment, so the whole chord D_j n line
                 // (P_i,P_k) is optimal and no enumeration is needed. This
                 // takes precedence over Degenerate, which it answers exactly.
                 // Numeric value 4 is appended, so codes 0-3 are unchanged.
};

inline const char* to_string(ChainStatus s) noexcept {
    switch (s) {
        case ChainStatus::Point:      return "point";
        case ChainStatus::Plateau:    return "plateau";
        case ChainStatus::PointDisk:  return "point_disk";
        case ChainStatus::Degenerate: return "degenerate";
        case ChainStatus::Straight:   return "straight";
    }
    return "?";
}

// `family` is the class of the closed-form theorem of the paper that produced the point: 'a' chord-angle
// tangency, 'b' tangent point from an anchor, 'c' line-circle intersection,
// 'e' the center of a disk of radius zero, 'f' the segment witness. Purely
// informational: it does not affect the returned optimum.
// `slot` distinguishes the individual points inside a family:
//   family 'a': 0 (external tangency; there may be two such points);
//   family 'b': 0,1 the two tangent points from p, 2,3 those from q;
//   family 'c': 0 the entry point E_p, 1 the entry point E_q.
// Purely informational: it does not affect the returned optimum. Note that
// duplicates are removed after construction, so a point that coincides with an
// earlier one keeps the earlier slot.
struct ChainCandidate { Vec2 b; double f; char family; std::uint8_t slot; };

struct ChainResult {
    static constexpr int kMaxCandidates = 20;
    ChainStatus status;
    Vec2   point;          // maximizer: the best candidate (on dD_j unless Straight)
    double theta_sum;      // f(point), radians
    double turning_sum;    // 3*pi - theta_sum
    double bound;          // 3*pi - |Delta|
    double delta;          // Delta in [-pi, pi)
    int    n_candidates;   // evaluated candidates (closed-form theorem of the paper)
    int    best_index;     // argmax f over candidates
    int    min_index;      // argmin f (the side of maximal turning)
    ChainCandidate candidates[kMaxCandidates];
};

// ------------------------------------------------------ geometric pieces ---

namespace detail {

// real roots of a2 x^2 + a1 x + a0 = 0 with the reference thresholds
inline int quad_roots(double a2, double a1, double a0, double r[2]) noexcept {
    double nrm = std::fabs(a2);
    if (std::fabs(a1) > nrm) nrm = std::fabs(a1);
    if (std::fabs(a0) > nrm) nrm = std::fabs(a0);
    if (nrm < 1e-300) nrm = 1e-300;
    if (std::fabs(a2) > 1e-13 * nrm) {
        const double disc = a1 * a1 - 4.0 * a2 * a0;
        if (disc < 0.0) return 0;
        const double s = std::sqrt(disc);
        r[0] = (-a1 - s) / (2.0 * a2);
        r[1] = (-a1 + s) / (2.0 * a2);
        return 2;
    }
    if (std::fabs(a1) > 1e-13 * nrm) { r[0] = -a0 / a1; return 1; }
    return 0;
}

// intersections of the line P + t u with the circle (c, r); 0 or 2 points
inline int line_circle(Vec2 P, Vec2 u, Vec2 c, double r, Vec2 out[2]) noexcept {
    const double nu = norm(u);
    if (nu == 0.0) return 0;
    const Vec2 uh = (1.0 / nu) * u;
    const Vec2 w  = P - c;
    const double bq   = dot(uh, w);
    const double cq   = dot(w, w) - r * r;
    const double disc = bq * bq - cq;
    if (disc < 0.0) return 0;
    const double s = std::sqrt(disc);
    out[0] = P + (-bq - s) * uh;
    out[1] = P + (-bq + s) * uh;
    return 2;
}


// tangent points from external b to the circle (c, r); 0 or 2
inline int tangent_points(Vec2 b, Vec2 c, double r, Vec2 out[2]) noexcept {
    const Vec2 e = b - c;
    const double d2 = dot(e, e);
    if (d2 <= r * r) return 0;
    const double fpar = r * r / d2;
    const double h    = r * std::sqrt(d2 - r * r) / d2;
    const Vec2 ep = perp(e);
    out[0] = c + fpar * e + h * ep;
    out[1] = c + fpar * e - h * ep;
    return 2;
}

// tangency points of circles through P_i, P_k with dD_j: the Problem R
// quadratic (6), tangency point (9); computed in EVERY branch
inline int tangency_circle_points(Vec2 Pi, Vec2 Pk, Vec2 cj, double rj,
                                  Vec2 out[2], std::uint8_t kind[2]) noexcept {
    const Vec2 u = Pk - Pi;
    const double c = norm(u);
    if (c == 0.0) return 0;
    const Vec2 M  = 0.5 * (Pi + Pk);
    const Vec2 th = (1.0 / c) * u;
    Vec2 nh = perp(th);
    if (dot(cj - M, nh) < 0.0) nh = -1.0 * nh;
    const double q = dot(cj - M, nh);
    const double S = dot(cj - M, cj - M) - 0.25 * c * c - rj * rj;
    double ks[2];
    const int nk = quad_roots(4.0 * (q * q - rj * rj), -4.0 * q * S,
                              S * S - rj * rj * c * c, ks);
    int n = 0;
    for (int i = 0; i < nk; ++i) {
        const Vec2  O = M + ks[i] * nh;
        const double R = std::sqrt(ks[i] * ks[i] + 0.25 * c * c);
        const double d = dist(cj, O);
        if (d > 1e-12) {
            // sign(S - 2kq) separates the external tangency from the internal
            // one: + external, - internal.
            kind[n] = (S - 2.0 * ks[i] * q >= 0.0) ? 0 : 1;
            out[n++] = O + (R / d) * (cj - O);
        }
    }
    return n;
}

// Straight-passage witness (class C_f): a point of the piece of the open
// segment (P_i,P_k) that lies inside D_j, from the segment quadratic (5).
// The WHOLE piece is optimal, so the point returned here is a tie-break. We
// return the point of maximal clearance to the two anchors,
//     argmax_z min{|z-P_i|, |z-P_k|}  over the piece,
// which keeps the iterates of a coordinate-descent smoother away from
// near-coincident waypoints, where the turning angles are ill conditioned.
// On the segment min{|z-P_i|,|z-P_k|} = |P_k-P_i| * min(t, 1-t) is unimodal
// with its peak at t = 1/2, so the maximiser over the piece [lo,hi] is simply
// t* = clamp(1/2, lo, hi): a closed form. Measured against the midpoint rule
// on 200,000 instances, this doubles the clearance at every quantile and
// halves the frequency of near-coincidences, at identical objective value.
// Fast path: t* = 1/2 whenever the midpoint of P_iP_k lies in D_j, since
// clamp(1/2,lo,hi) = 1/2 exactly in that case. The test costs one squared
// distance, so the segment quadratic below is solved only when the midpoint
// falls outside the disk.
inline bool pi_witness(Vec2 Pi, Vec2 Pk, Vec2 cj, double rj, double tol,
                       Vec2& out) noexcept {
    const Vec2 M{0.5 * (Pi.x + Pk.x), 0.5 * (Pi.y + Pk.y)};
    if (dist2(M, cj) <= rj * rj) { out = M; return true; }
    const Vec2 u = Pk - Pi;
    const double c  = norm(u);
    double scale = c;
    { const double t1 = dist(cj, Pi); if (t1 > scale) scale = t1; }
    if (rj  > scale) scale = rj;
    if (1.0 > scale) scale = 1.0;
    const Vec2 w = cj - Pi;
    const double bq = dot(u, w);
    const double aq = c * c;
    const double cq = dot(w, w) - rj * rj;
    double disc = bq * bq - aq * cq;
    double thr = bq * bq;
    if (std::fabs(aq * cq) > thr) thr = std::fabs(aq * cq);
    if (thr < 1e-300) thr = 1e-300;
    if (disc < -1e-12 * thr || aq == 0.0) return false;
    if (disc < 0.0) disc = 0.0;
    const double s = std::sqrt(disc);
    double lo = (bq - s) / aq; if (lo < 0.0) lo = 0.0;
    double hi = (bq + s) / aq; if (hi > 1.0) hi = 1.0;
    if (lo > hi) return false;
    const double tm = (0.5 < lo) ? lo : ((0.5 > hi) ? hi : 0.5);  // clamp(1/2,lo,hi)
    const Vec2 P = Pi + tm * u;
    const bool inside = dist(P, cj) <= rj + tol * scale;
    const double margin = (tm < 1.0 - tm) ? tm : 1.0 - tm;
    if (inside && margin > 1e-9) { out = P; return true; }
    return false;
}

} // namespace detail

// ----------------------------------------------------------------- solver ---

// Maximize theta_i + theta_j + theta_k over b in ball(c_j, r_j).
// Returns false on invalid input (coincident anchors, r_j < 0, or an empty
// candidate set in a pathological degenerate configuration); on success,
// fills `res`. Stateless and reentrant: safe to call from many threads.
// `short_circuit` enables the straight-passage early return. It is valid for
// the MINIMISATION of the total turning (equivalently, the maximisation of
// theta_sum) only: it returns without building the candidate set, so the
// opposite direction (maximal turning, res.min_index) is not
// available in that branch. Pass false to always build the full set.
inline bool solve_chain(Vec2 H, Vec2 Pi, Vec2 Pk, Vec2 L, Vec2 cj, double rj,
                        ChainResult& res, double tol = 1e-9,
                        bool short_circuit = true) noexcept {
    using namespace detail;

    const double dic = dist(Pi, cj), dkc = dist(Pk, cj);
    double scale = 1.0;
    if (dic > scale) scale = dic;
    if (dkc > scale) scale = dkc;
    if (rj > scale) scale = rj;
    const double ts = tol * scale, ts2 = ts * ts;
    if (dist2(H, Pi) <= ts2) return false;   // incoming heading undefined
    if (dist2(L, Pk) <= ts2) return false;   // outgoing heading undefined
    if (dist2(Pi, Pk) <= ts2) return false;  // P_i and P_k coincide
    if (rj < 0.0) return false;

    const double delta = wrap_pi(heading(L - Pk) - heading(Pi - H));
    const double bound = 3.0 * kPi - std::fabs(delta);
    res.delta = delta;
    res.bound = bound;

    const bool degenerate = (dic <= rj + ts) || (dkc <= rj + ts);

    if (rj <= tol * scale * 1e-3 && !degenerate) {          // point disk
        const double f = theta_sum(cj, H, Pi, Pk, L);
        res.status = ChainStatus::PointDisk;
        res.point = cj; res.theta_sum = f; res.turning_sum = 3.0 * kPi - f;
        res.n_candidates = 1; res.best_index = 0; res.min_index = 0;
        res.candidates[0] = {cj, f, 'e', 0};
        return true;
    }

    // ---- straight-passage branch (Corollary: the whole chord is optimal) ---
    // psi is globally minimised on the open segment (P_i,P_k); when that
    // segment meets D_j the disk constraint is inactive and one formula
    // answers the block update. The proof of this result does not use
    // P_i, P_k outside D_j -- it only needs (z-P_i) + (P_k-z) = P_k-P_i and
    // the triangle inequality on the circle of directions -- so this branch
    // is taken even when the `degenerate` flag is set. Indeed P_i or P_k
    // inside D_j FORCES the segment to meet D_j, so the configuration that
    // the flag reports is exactly one this branch answers exactly.
    if (short_circuit) {
        Vec2 wpt;
        if (pi_witness(Pi, Pk, cj, rj, tol, wpt)) {
            const double f = theta_sum(wpt, H, Pi, Pk, L);
            res.status = ChainStatus::Straight;
            res.point = wpt; res.theta_sum = f; res.turning_sum = 3.0 * kPi - f;
            res.n_candidates = 1; res.best_index = 0; res.min_index = 0;
            res.candidates[0] = {wpt, f, 'f', 0};
            return true;
        }
    }

    // ---- candidate set of the theorem: classes (a), (b), (c) -------------
    // Eight points at most, ALL on the circle dD_j (Theorem "Closed-form
    // solution of the single-block subproblem" of the paper):
    //   (a) the points where a circle through P_i and P_k touches dD_j from
    //       OUTSIDE (S - 2k*eta > 0): one point when the line P_iP_k misses
    //       D_j, two when it meets it. The internal tangency is a local
    //       maximum of psi along dD_j and is never needed;
    //   (b) the four points where a line through P_i or through P_k touches
    //       dD_j (never empty, since P_i, P_k lie off D_j);
    //   (c) the points E_p and E_q where the extensions of the edges
    //       H->P_i beyond P_i and L->P_k beyond P_k enter D_j, each taken only
    //       when that extension meets the disk. On the other half of l_p
    //       (resp. l_q) the turn at P_i (resp. P_k) equals pi, and along the
    //       extension psi is nondecreasing away from the anchor, so the exit
    //       point is never better than the entry point.
    // Class (f), the straight-passage witness, is returned by the short
    // circuit above and is empty whenever that branch does not apply.
    Vec2 raw[ChainResult::kMaxCandidates];
    char fam[ChainResult::kMaxCandidates];
    std::uint8_t slt[ChainResult::kMaxCandidates];
    int  nraw = 0;
    { Vec2 t[2]; std::uint8_t kd[2];
      const int n = tangency_circle_points(Pi, Pk, cj, rj, t, kd);
      for (int i = 0; i < n; ++i)
        if (kd[i] == 0)                                   // external tangency only
          { raw[nraw] = t[i]; slt[nraw] = 0; fam[nraw++] = 'a'; } }
    { Vec2 wpt; if (pi_witness(Pi, Pk, cj, rj, tol, wpt))
        { raw[nraw] = wpt; slt[nraw] = 0; fam[nraw++] = 'f'; } }
    { Vec2 t[2]; int n = tangent_points(Pi, cj, rj, t);
      for (int i = 0; i < n; ++i)
        { raw[nraw] = t[i]; slt[nraw] = (std::uint8_t)i; fam[nraw++] = 'b'; }
      n = tangent_points(Pk, cj, rj, t);
      for (int i = 0; i < n; ++i)
        { raw[nraw] = t[i]; slt[nraw] = (std::uint8_t)(2 + i); fam[nraw++] = 'b'; } }
    // (c): the entry points E_p and E_q, eq. (entry-point) of the paper. The
    // extension of H->P_i beyond P_i is P_i + s*u with u = P_i - H, s > 0;
    // since P_i lies off D_j, both crossings of its line with dD_j have the
    // sign of mu = <u, c_j - P_i>, so the extension meets the disk only if
    // mu > 0, and then it enters at the smaller root, which line_circle()
    // returns first. Same for E_q with P_k and u = P_k - L.
    { const Vec2 anchor[2] = {Pi, Pk};
      const Vec2 dir[2]    = {Pi - H, Pk - L};
      for (int m = 0; m < 2; ++m) {                          // (c) <= 2
        if (dot(dir[m], cj - anchor[m]) <= 0.0) continue;    // disk not ahead
        Vec2 t[2];
        if (line_circle(anchor[m], dir[m], cj, rj, t) == 2)
          { raw[nraw] = t[0]; slt[nraw] = (std::uint8_t)m; fam[nraw++] = 'c'; }
      } }

    // ---- filter (feasibility, anchor proximity) and dedupe ---------------
    // all threshold tests in squared form: same predicates, no square roots
    const double rfeas = rj * (1.0 + 1e-9) + 1e-12;
    const double rfeas2 = rfeas * rfeas;
    const double prox  = 1e-7 * scale, prox2 = prox * prox;
    int nc = 0, best = -1, worst = -1;
    double fbest = -1.0, fworst = 1e300;
    for (int i = 0; i < nraw; ++i) {
        const Vec2 X = raw[i];
        if (dist2(X, cj) > rfeas2) continue;
        if (dist2(X, Pi) <= prox2 || dist2(X, Pk) <= prox2) continue;
        bool dup = false;
        for (int j = 0; j < nc; ++j)
            if (std::fabs(X.x - res.candidates[j].b.x) <= 1e-9 &&
                std::fabs(X.y - res.candidates[j].b.y) <= 1e-9) { dup = true; break; }
        if (dup) continue;
        const double f = theta_sum(X, H, Pi, Pk, L);
        res.candidates[nc] = {X, f, fam[i], slt[i]};
        if (f > fbest)  { fbest = f;  best  = nc; }
        if (f < fworst) { fworst = f; worst = nc; }
        ++nc;
    }
    if (nc == 0) {
        // Cannot happen under the hypotheses of the theorem: outside the
        // straight-passage branch P_i and P_k lie off D_j, so the four tangent
        // points of class (b) always exist and lie on dD_j. It can only be
        // reached with the short circuit disabled on a degenerate instance;
        // the center is then used as a feasible fallback, not as a candidate
        // class of the theorem.
        const double f = theta_sum(cj, H, Pi, Pk, L);
        res.candidates[0] = {cj, f, 'e', 0};
        nc = 1; best = 0; worst = 0; fbest = f; fworst = f;
    }
    res.n_candidates = nc; res.best_index = best; res.min_index = worst;

    const Vec2  bstar = res.candidates[best].b;
    const double fstar = fbest;

    if (degenerate) {
        res.status = ChainStatus::Degenerate;
        res.point = bstar; res.theta_sum = fstar;
        res.turning_sum = 3.0 * kPi - fstar;
        return true;
    }

    // The returned point is always the best candidate of the theorem, which
    // lies on dD_j. When its value attains the bound 3*pi - |Delta| (the
    // smallest possible total turning, attained on the plateau region), the
    // status reports it, but the point is not replaced: every point of the
    // plateau inside D_j has the same value.
    res.status = (fstar >= bound - 1e-7) ? ChainStatus::Plateau : ChainStatus::Point;
    res.point = bstar; res.theta_sum = fstar;
    res.turning_sum = 3.0 * kPi - fstar;
    return true;
}

} // namespace maxangle
