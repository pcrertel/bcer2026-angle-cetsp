// ==========================================================================
// chain_c_api.cpp -- C interface to the closed-form block update
// (max_angle_chain.hpp), called from Fortran through iso_c_binding.
//
// Build (default; add -DMAXANGLE_USE_SLEEF -include sleef_compat.hpp and link -lsleef
// for the faster atan2 of SLEEF, see max_angle_chain.hpp):
//     g++ -O3 -march=native -std=c++17 -fno-math-errno -c chain_c_api.cpp
//     gfortran yourcode.f90 chain_c_api.o -lstdc++ -o yourprog
//
// One call optimizes ONE window: chain i- -> i -> j -> k -> k+ with the
// point of j free inside the disk (cj, rj); the four neighbor POINTS are
// fixed. It returns the optimal location of j and the three deviation
// angles (pi - interior angle, radians) at i, j, k evaluated there.
// Minimizing dev(1)+dev(2)+dev(3) is the objective. The call is stateless
// and thread-safe.
// ==========================================================================
#include "max_angle_chain.hpp"

extern "C"
// returns  0 = regular optimum ("point")      1 = plateau (bound attained)
//          2 = point disk (rj = 0)            3 = degenerate (i or k in disk j)
//          4 = straight (segment i-k meets disk j; the whole chord is optimal)
//         -1 = invalid input (coincident anchors or rj < 0)
// inputs : h(2) = point of i-,  p_i(2) = point of i,  p_k(2) = point of k,
//          l(2) = point of k+,  cj(2) = center of disk j,  rj = its radius
// outputs: b(2) = optimal point of j (feasible: |b - cj| <= rj)
//          dev(3) = deviations at i, j, k evaluated at b (radians)
int solve_chain_c(const double* h, const double* p_i, const double* p_k,
                  const double* l, const double* cj, double rj,
                  double* b, double* dev)
{
    using namespace maxangle;
    const Vec2 H{h[0], h[1]}, Pi{p_i[0], p_i[1]}, Pk{p_k[0], p_k[1]},
               L{l[0], l[1]}, C{cj[0], cj[1]};
    ChainResult r;
    if (!solve_chain(H, Pi, Pk, L, C, rj, r)) return -1;
    // snap onto the disk so the absolute feasibility certificate holds
    Vec2 q = r.point;
    const double d = dist(q, C);
    if (d > rj && d > 0.0) q = C + (rj / d) * (q - C);
    b[0] = q.x; b[1] = q.y;
    dev[0] = kPi - angle_at(Pi, H, q);
    dev[1] = kPi - angle_at(q, Pi, Pk);
    dev[2] = kPi - angle_at(Pk, q, L);
    switch (r.status) {
        case ChainStatus::Point:      return 0;
        case ChainStatus::Plateau:    return 1;
        case ChainStatus::PointDisk:  return 2;
        case ChainStatus::Degenerate: return 3;
        case ChainStatus::Straight:   return 4;
    }
    return -1;
}
