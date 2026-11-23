! ==========================================================================
! subproblem_timing.f90 -- time of one closed-form block update (Section 5.1 of
! the paper): the call solve_chain_f of max_angle_chain.f90, the same call made
! by the ALNS-BCD in discmin_cf. Random configurations:
!   * disk at the origin, radius r ~ U[0.5, 10]; scale BOX ~ U[1, 10];
!   * draws with sqrt(8) BOX <= r (empty annulus) are discarded;
!   * p and q uniform in area in the annulus between r and sqrt(8) BOX;
!   * p- and q+ uniform in the full disk of radius sqrt(8) BOX (they may lie
!     in the disk);
!   * configurations in which the segment (p,q) meets the disk are discarded;
!   * N accepted configurations (default 1,000,000).
! Prints the fractions of discarded draws and the mean time per call (best of
! three passes).
! Usage: ./subproblem_timing [N] [seed]
! ==========================================================================
program subproblem_timing
    use max_angle_chain_f, only: solve_chain_f
    implicit none
    real(kind=8), parameter :: pi2 = 6.28318530717958647692d0
    integer :: n, sd, nseed, k, i, rep, st, n_invalid, straight
    integer(kind=8) :: draws, n_empty, n_touch, n_inside, c0, c1, crate
    integer, allocatable :: seed(:)
    real(kind=8), allocatable :: X(:, :)
    real(kind=8) :: r, box, rout, u(2), H(2), P(2), Q(2), L(2), b(2), dv(3), best, s, cj(2)
    character(len=32) :: arg

    n = 1000000; sd = 20260901
    if ( command_argument_count() >= 1 ) then
        call get_command_argument(1, arg); read(arg, *) n
    end if
    if ( command_argument_count() >= 2 ) then
        call get_command_argument(2, arg); read(arg, *) sd
    end if
    call random_seed(size=nseed)
    allocate(seed(nseed)); seed = sd + 37 * [(k, k = 1, nseed)]
    call random_seed(put=seed)

    allocate(X(9, n))
    draws = 0; n_empty = 0; n_touch = 0; n_inside = 0
    i = 0
    do while ( i < n )
        draws = draws + 1
        call random_number(u)
        r = 0.5d0 + 9.5d0 * u(1)
        box = 1.0d0 + 9.0d0 * u(2)
        rout = sqrt(8.0d0) * box
        if ( rout <= r ) then
            n_empty = n_empty + 1; cycle
        end if
        call annulus_point(P); call annulus_point(Q); call disk_point(H); call disk_point(L)
        if ( dist_origin_segment(P, Q) <= r ) then
            n_touch = n_touch + 1; cycle
        end if
        if ( norm2(H) < r .or. norm2(L) < r ) n_inside = n_inside + 1
        i = i + 1
        X(1:2, i) = H; X(3:4, i) = P; X(5:6, i) = Q; X(7:8, i) = L; X(9, i) = r
    end do

    cj = 0.0d0
    n_invalid = 0; straight = 0
    do i = 1, n
        st = solve_chain_f(X(1:2, i), X(3:4, i), X(5:6, i), X(7:8, i), cj, X(9, i), b, dv)
        if ( st < 0 ) n_invalid = n_invalid + 1
        if ( st == 4 ) straight = straight + 1
    end do

    best = huge(1.0d0); s = 0.0d0
    do rep = 1, 3
        call system_clock(c0, crate)
        do i = 1, n
            st = solve_chain_f(X(1:2, i), X(3:4, i), X(5:6, i), X(7:8, i), cj, X(9, i), b, dv)
            s = s + dv(2)
        end do
        call system_clock(c1)
        best = min(best, real(c1 - c0, 8) / real(crate, 8))
    end do

    write(*, '(a, i0, a, i0)') '# subproblem_timing.f90  N=', n, '  seed=', sd
    write(*, '(a, i0)') 'draws: ', draws
    write(*, '(a, i0, a, f6.2, a)') '  discarded, empty annulus (sqrt8*BOX <= r): ', n_empty, '  (', &
        100.0d0 * n_empty / draws, '% of the draws)'
    write(*, '(a, i0, a, f6.2, a)') '  discarded, segment meets the disk:         ', n_touch, '  (', &
        100.0d0 * n_touch / draws, '% of the draws)'
    write(*, '(a, i0, a, f6.2, a)') 'accepted configurations: ', n, ' ; with p- or q+ inside D: ', &
        100.0d0 * n_inside / n, '%'
    write(*, '(a, i0, a, i0)') 'invalid calls: ', n_invalid, ' ; straight status (should not occur): ', straight
    write(*, '(a, es10.3, a)') 'mean time per call (best of 3 passes): ', best / n, ' s'
    if ( s == -1.0d0 ) write(*, *) s

contains

    subroutine annulus_point(z)
        real(kind=8), intent(out) :: z(2)
        real(kind=8) :: w(2), rad
        call random_number(w)
        rad = sqrt(r * r + w(1) * (rout * rout - r * r))
        z(1) = rad * cos(pi2 * w(2)); z(2) = rad * sin(pi2 * w(2))
    end subroutine annulus_point

    subroutine disk_point(z)
        real(kind=8), intent(out) :: z(2)
        real(kind=8) :: w(2), rad
        call random_number(w)
        rad = rout * sqrt(w(1))
        z(1) = rad * cos(pi2 * w(2)); z(2) = rad * sin(pi2 * w(2))
    end subroutine disk_point

    ! distance from the origin to the segment [a, b]
    real(kind=8) function dist_origin_segment(a, b)
        real(kind=8), intent(in) :: a(2), b(2)
        real(kind=8) :: dx, dy, L2, t
        dx = b(1) - a(1); dy = b(2) - a(2); L2 = dx * dx + dy * dy
        t = 0.0d0
        if ( L2 > 0.0d0 ) t = -(a(1) * dx + a(2) * dy) / L2
        t = max(0.0d0, min(1.0d0, t))
        dist_origin_segment = sqrt((a(1) + t * dx)**2 + (a(2) + t * dy)**2)
    end function dist_origin_segment

end program subproblem_timing
