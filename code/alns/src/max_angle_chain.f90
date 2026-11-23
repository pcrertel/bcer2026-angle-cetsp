! ==========================================================================
! max_angle_chain.f90 -- closed-form solution of the single-disk subproblem
! of the BCD (Section 3 of the paper).
!
! Given the fixed points h, p_i, p_k, l and the disk D_j (center cj,
! radius rj), it returns a point b of D_j that minimizes the sum of the
! turning angles at p_i, b and p_k:
!   * if the open segment (p_i,p_k) meets the disk ("straight" case), the
!     point of the segment inside the disk farthest from the anchors;
!   * otherwise, the best of the candidate points of the theorem (at most
!     eight, all on the boundary of D_j): (a) points of tangency of circles
!     through p_i and p_k, (b) points of tangency of the lines through p_i
!     and through p_k, (c) points where the extensions of the edges enter
!     the disk.
! All routines are pure: no state, reentrant.
!
!   st = solve_chain_f(h, p_i, p_k, l, cj, rj, b, dev)
!   st:   0 point, 1 plateau, 2 point disk, 3 degenerate, 4 straight,
!        -1 invalid input
!   b(2): optimal point of the block (projected onto the disk);
!   dev(3): turning angles at i, j and k
! ==========================================================================
module max_angle_chain_f

    implicit none
    private
    public :: solve_chain_f

    real(kind=8), parameter :: kPi    = 3.14159265358979323846d0
    real(kind=8), parameter :: kTwoPi = 6.28318530717958647692d0
    integer, parameter :: kMaxCand = 20

contains

    pure real(kind=8) function dot2(a, b)
        real(kind=8), intent(in) :: a(2), b(2)
        dot2 = a(1) * b(1) + a(2) * b(2)
    end function dot2

    pure real(kind=8) function cross2(a, b)
        real(kind=8), intent(in) :: a(2), b(2)
        cross2 = a(1) * b(2) - a(2) * b(1)
    end function cross2

    pure function perp2(a) result(p)
        real(kind=8), intent(in) :: a(2)
        real(kind=8) :: p(2)
        p(1) = -a(2); p(2) = a(1)
    end function perp2

    pure real(kind=8) function normv(a)
        real(kind=8), intent(in) :: a(2)
        normv = sqrt(dot2(a, a))
    end function normv

    pure real(kind=8) function distv(a, b)
        real(kind=8), intent(in) :: a(2), b(2)
        distv = normv(a - b)
    end function distv

    pure real(kind=8) function dist2v(a, b)
        real(kind=8), intent(in) :: a(2), b(2)
        dist2v = dot2(a - b, a - b)
    end function dist2v

    pure real(kind=8) function heading(v)
        real(kind=8), intent(in) :: v(2)
        heading = atan2(v(2), v(1))
    end function heading

    ! reduces to [-pi, pi)
    pure real(kind=8) function wrap_pi(a)
        real(kind=8), intent(in) :: a
        real(kind=8) :: w
        w = mod(a + kPi, kTwoPi)
        if ( w < 0.0d0 ) w = w + kTwoPi
        wrap_pi = w - kPi
    end function wrap_pi

    ! angle A-P-B at the vertex P, in [0, pi], atan2 form
    pure real(kind=8) function angle_at(P, A, B)
        real(kind=8), intent(in) :: P(2), A(2), B(2)
        real(kind=8) :: u(2), v(2)
        u = A - P; v = B - P
        angle_at = atan2(abs(cross2(u, v)), dot2(u, v))
    end function angle_at

    ! f(b) = theta_i + theta_j + theta_k (= 3 pi - sum of the turning angles)
    pure real(kind=8) function theta_sum(b, H, Pi, Pk, L)
        real(kind=8), intent(in) :: b(2), H(2), Pi(2), Pk(2), L(2)
        theta_sum = angle_at(Pi, H, b) + angle_at(b, Pi, Pk) + angle_at(Pk, b, L)
    end function theta_sum

    ! real roots of a2 x^2 + a1 x + a0 = 0
    pure subroutine quad_roots(a2, a1, a0, r, n)
        real(kind=8), intent(in) :: a2, a1, a0
        real(kind=8), intent(out) :: r(2)
        integer, intent(out) :: n
        real(kind=8) :: nrm, disc, s
        r = 0.0d0
        nrm = abs(a2)
        if ( abs(a1) > nrm ) nrm = abs(a1)
        if ( abs(a0) > nrm ) nrm = abs(a0)
        if ( nrm < 1.0d-300 ) nrm = 1.0d-300
        if ( abs(a2) > 1.0d-13 * nrm ) then
            disc = a1 * a1 - 4.0d0 * a2 * a0
            if ( disc < 0.0d0 ) then
                n = 0; return
            end if
            s = sqrt(disc)
            r(1) = (-a1 - s) / (2.0d0 * a2)
            r(2) = (-a1 + s) / (2.0d0 * a2)
            n = 2; return
        end if
        if ( abs(a1) > 1.0d-13 * nrm ) then
            r(1) = -a0 / a1; n = 1; return
        end if
        n = 0
    end subroutine quad_roots

    ! intersections of the line P + t u with the circle (c, r): 0 or 2 points
    pure subroutine line_circle(P, u, c, r, out, n)
        real(kind=8), intent(in) :: P(2), u(2), c(2), r
        real(kind=8), intent(out) :: out(2, 2)
        integer, intent(out) :: n
        real(kind=8) :: nu, uh(2), w(2), bq, cq, disc, s
        out = 0.0d0
        nu = normv(u)
        if ( nu == 0.0d0 ) then
            n = 0; return
        end if
        uh = (1.0d0 / nu) * u
        w = P - c
        bq = dot2(uh, w)
        cq = dot2(w, w) - r * r
        disc = bq * bq - cq
        if ( disc < 0.0d0 ) then
            n = 0; return
        end if
        s = sqrt(disc)
        out(:, 1) = P + (-bq - s) * uh
        out(:, 2) = P + (-bq + s) * uh
        n = 2
    end subroutine line_circle

    ! points of tangency from a point b outside the circle (c, r): 0 or 2
    pure subroutine tangent_points(b, c, r, out, n)
        real(kind=8), intent(in) :: b(2), c(2), r
        real(kind=8), intent(out) :: out(2, 2)
        integer, intent(out) :: n
        real(kind=8) :: e(2), d2, fpar, h, ep(2)
        out = 0.0d0
        e = b - c
        d2 = dot2(e, e)
        if ( d2 <= r * r ) then
            n = 0; return
        end if
        fpar = r * r / d2
        h = r * sqrt(d2 - r * r) / d2
        ep = perp2(e)
        out(:, 1) = c + fpar * e + h * ep
        out(:, 2) = c + fpar * e - h * ep
        n = 2
    end subroutine tangent_points

    ! points of tangency of circles through P_i and P_k with dD_j; kind = 0 external, 1 internal
    pure subroutine tangency_circle_points(Pi, Pk, cj, rj, out, kind, n)
        real(kind=8), intent(in) :: Pi(2), Pk(2), cj(2), rj
        real(kind=8), intent(out) :: out(2, 2)
        integer, intent(out) :: kind(2), n
        real(kind=8) :: u(2), c, M(2), th(2), nh(2), q, S, ks(2), O(2), R, d
        integer :: nk, i
        out = 0.0d0; kind = 0; n = 0
        u = Pk - Pi
        c = normv(u)
        if ( c == 0.0d0 ) return
        M = 0.5d0 * (Pi + Pk)
        th = (1.0d0 / c) * u
        nh = perp2(th)
        if ( dot2(cj - M, nh) < 0.0d0 ) nh = -1.0d0 * nh
        q = dot2(cj - M, nh)
        S = dot2(cj - M, cj - M) - 0.25d0 * c * c - rj * rj
        call quad_roots(4.0d0 * (q * q - rj * rj), -4.0d0 * q * S, S * S - rj * rj * c * c, ks, nk)
        do i = 1, nk
            O = M + ks(i) * nh
            R = sqrt(ks(i) * ks(i) + 0.25d0 * c * c)
            d = distv(cj, O)
            if ( d > 1.0d-12 ) then
                n = n + 1
                if ( S - 2.0d0 * ks(i) * q >= 0.0d0 ) then
                    kind(n) = 0
                else
                    kind(n) = 1
                end if
                out(:, n) = O + (R / d) * (cj - O)
            end if
        end do
    end subroutine tangency_circle_points

    ! point of the open segment (P_i,P_k) inside D_j farthest from the anchors
    pure subroutine pi_witness(Pi, Pk, cj, rj, tol, out, achou)
        real(kind=8), intent(in) :: Pi(2), Pk(2), cj(2), rj, tol
        real(kind=8), intent(out) :: out(2)
        logical, intent(out) :: achou
        real(kind=8) :: Mp(2), u(2), c, scale, t1, w(2), bq, aq, cq, disc, thr, s, lo, hi, tm, P(2), margin
        logical :: inside
        out = 0.0d0
        achou = .false.
        Mp = 0.5d0 * (Pi + Pk)
        if ( dist2v(Mp, cj) <= rj * rj ) then
            out = Mp; achou = .true.; return
        end if
        u = Pk - Pi
        c = normv(u)
        scale = c
        t1 = distv(cj, Pi); if ( t1 > scale ) scale = t1
        if ( rj > scale ) scale = rj
        if ( 1.0d0 > scale ) scale = 1.0d0
        w = cj - Pi
        bq = dot2(u, w)
        aq = c * c
        cq = dot2(w, w) - rj * rj
        disc = bq * bq - aq * cq
        thr = bq * bq
        if ( abs(aq * cq) > thr ) thr = abs(aq * cq)
        if ( thr < 1.0d-300 ) thr = 1.0d-300
        if ( disc < -1.0d-12 * thr .or. aq == 0.0d0 ) return
        if ( disc < 0.0d0 ) disc = 0.0d0
        s = sqrt(disc)
        lo = (bq - s) / aq; if ( lo < 0.0d0 ) lo = 0.0d0
        hi = (bq + s) / aq; if ( hi > 1.0d0 ) hi = 1.0d0
        if ( lo > hi ) return
        if ( 0.5d0 < lo ) then
            tm = lo
        else if ( 0.5d0 > hi ) then
            tm = hi
        else
            tm = 0.5d0
        end if
        P = Pi + tm * u
        inside = distv(P, cj) <= rj + tol * scale
        if ( tm < 1.0d0 - tm ) then
            margin = tm
        else
            margin = 1.0d0 - tm
        end if
        if ( inside .and. margin > 1.0d-9 ) then
            out = P; achou = .true.
        end if
    end subroutine pi_witness

    ! core of solve_chain (tol = 1e-9, short_circuit = .true.): status and point
    pure subroutine solve_core(H, Pi, Pk, L, cj, rj, ok, status, point)
        real(kind=8), intent(in) :: H(2), Pi(2), Pk(2), L(2), cj(2), rj
        logical, intent(out) :: ok
        integer, intent(out) :: status
        real(kind=8), intent(out) :: point(2)
        real(kind=8), parameter :: tol = 1.0d-9
        real(kind=8) :: dic, dkc, scale, ts, ts2, delta, bound, wpt(2), t(2, 2), X(2), f, fbest
        real(kind=8) :: raw(2, kMaxCand), anchor(2, 2), dir(2, 2), rfeas, rfeas2, prox, prox2
        real(kind=8) :: cand(2, kMaxCand)
        integer :: kd(2), nraw, n, i, j, m, nc, best
        logical :: degenerate, dup, achou

        ok = .false.; status = -1; point = 0.0d0
        dic = distv(Pi, cj); dkc = distv(Pk, cj)
        scale = 1.0d0
        if ( dic > scale ) scale = dic
        if ( dkc > scale ) scale = dkc
        if ( rj > scale ) scale = rj
        ts = tol * scale; ts2 = ts * ts
        if ( dist2v(H, Pi) <= ts2 ) return
        if ( dist2v(L, Pk) <= ts2 ) return
        if ( dist2v(Pi, Pk) <= ts2 ) return
        if ( rj < 0.0d0 ) return

        delta = wrap_pi(heading(L - Pk) - heading(Pi - H))
        bound = 3.0d0 * kPi - abs(delta)
        degenerate = (dic <= rj + ts) .or. (dkc <= rj + ts)

        if ( rj <= tol * scale * 1.0d-3 .and. .not. degenerate ) then   ! disk of radius zero
            ok = .true.; status = 2; point = cj; return
        end if

        ! straight case: the open segment meets the disk
        call pi_witness(Pi, Pk, cj, rj, tol, wpt, achou)
        if ( achou ) then
            ok = .true.; status = 4; point = wpt; return
        end if

        ! candidate points of the theorem: (a), (f), (b), (c)
        nraw = 0
        call tangency_circle_points(Pi, Pk, cj, rj, t, kd, n)
        do i = 1, n
            if ( kd(i) == 0 ) then
                nraw = nraw + 1; raw(:, nraw) = t(:, i)
            end if
        end do
        call pi_witness(Pi, Pk, cj, rj, tol, wpt, achou)
        if ( achou ) then
            nraw = nraw + 1; raw(:, nraw) = wpt
        end if
        call tangent_points(Pi, cj, rj, t, n)
        do i = 1, n
            nraw = nraw + 1; raw(:, nraw) = t(:, i)
        end do
        call tangent_points(Pk, cj, rj, t, n)
        do i = 1, n
            nraw = nraw + 1; raw(:, nraw) = t(:, i)
        end do
        anchor(:, 1) = Pi; anchor(:, 2) = Pk
        dir(:, 1) = Pi - H; dir(:, 2) = Pk - L
        do m = 1, 2
            if ( dot2(dir(:, m), cj - anchor(:, m)) <= 0.0d0 ) cycle
            call line_circle(anchor(:, m), dir(:, m), cj, rj, t, n)
            if ( n == 2 ) then
                nraw = nraw + 1; raw(:, nraw) = t(:, 1)
            end if
        end do

        ! filter (feasibility, distance to the anchors) and duplicates
        rfeas = rj * (1.0d0 + 1.0d-9) + 1.0d-12
        rfeas2 = rfeas * rfeas
        prox = 1.0d-7 * scale; prox2 = prox * prox
        nc = 0; best = 0; fbest = -1.0d0
        do i = 1, nraw
            X = raw(:, i)
            if ( dist2v(X, cj) > rfeas2 ) cycle
            if ( dist2v(X, Pi) <= prox2 .or. dist2v(X, Pk) <= prox2 ) cycle
            dup = .false.
            do j = 1, nc
                if ( abs(X(1) - cand(1, j)) <= 1.0d-9 .and. abs(X(2) - cand(2, j)) <= 1.0d-9 ) then
                    dup = .true.; exit
                end if
            end do
            if ( dup ) cycle
            f = theta_sum(X, H, Pi, Pk, L)
            nc = nc + 1
            cand(:, nc) = X
            if ( f > fbest ) then
                fbest = f; best = nc
            end if
        end do
        if ( nc == 0 ) then                      ! safeguard: the center of the disk
            nc = 1; best = 1; cand(:, 1) = cj
            fbest = theta_sum(cj, H, Pi, Pk, L)
        end if

        ok = .true.
        point = cand(:, best)
        if ( degenerate ) then
            status = 3
        else if ( fbest >= bound - 1.0d-7 ) then
            status = 1
        else
            status = 0
        end if
    end subroutine solve_core

    ! entry point
    integer function solve_chain_f(h, p_i, p_k, l, cj, rj, b, dev)
        real(kind=8), intent(in) :: h(2), p_i(2), p_k(2), l(2), cj(2), rj
        real(kind=8), intent(out) :: b(2), dev(3)
        logical :: ok
        integer :: status
        real(kind=8) :: q(2), d
        call solve_core(h, p_i, p_k, l, cj, rj, ok, status, q)
        if ( .not. ok ) then
            b = 0.0d0; dev = 0.0d0
            solve_chain_f = -1; return
        end if
        ! projects onto the disk, so that b is feasible
        d = distv(q, cj)
        if ( d > rj .and. d > 0.0d0 ) q = cj + (rj / d) * (q - cj)
        b = q
        dev(1) = kPi - angle_at(p_i, h, q)
        dev(2) = kPi - angle_at(q, p_i, p_k)
        dev(3) = kPi - angle_at(p_k, q, l)
        solve_chain_f = status
    end function solve_chain_f

end module max_angle_chain_f
