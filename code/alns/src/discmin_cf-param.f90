    ! =====================================================================
    ! discmin_cf  -  closed-form block update of the BCD (single-disk subproblem).
    !
    ! The C solver (max_angle_chain.hpp / chain_c_api.cpp) returns an optimal
    ! visit point b of the window. Acceptance and the angle bookkeeping are done
    ! with psi0 (turningangle, atan2), so that the O(1) running cost stays
    ! consistent with evalfangle. The point b is accepted unless it raises psi by
    ! more than tol_acc, so that ties are resolved in favour of b: in the segment
    ! case, b is the point of maximal clearance, which keeps consecutive visit
    ! points apart. The tolerance only absorbs rounding, since b is an exact
    ! minimizer of the window. On output, zmin and the three window angles are
    ! updated and diff = psi_new - psi_old.
    ! =====================================================================
    subroutine discmin_cf(inst,param,count,pol,zpp,zpr,zne,znn,zmin,diff,&
    anglepr,anglemin,anglene)

        use iso_c_binding
        implicit none

        ! TYPE ARGUMENTS
        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(in) :: param
        type(counter_type),intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer,intent(in) :: pol
        real(kind=8),intent(inout) :: anglepr,anglemin,anglene
        real(kind=8),intent(out) :: diff

        ! ARRAY ARGUMENTS
        real(kind=8),intent(in) :: zpp(2),zpr(2),zne(2),znn(2)
        real(kind=8),intent(inout) :: zmin(2)

        ! PARAMETERS
        real(kind=8), parameter :: tol_acc = 1.0d-8

        ! C INTERFACE
        interface
           integer(c_int) function solve_chain_c(h,p_i,p_k,l,cj,rj,b,dev) &
                 bind(C, name="solve_chain_c")
             use iso_c_binding
             real(c_double), intent(in)  :: h(2),p_i(2),p_k(2),l(2),cj(2)
             real(c_double), value       :: rj
             real(c_double), intent(out) :: b(2),dev(3)
           end function solve_chain_c
        end interface

        ! LOCAL SCALARS
        real(kind=8) :: initialpsi,newpsi,r,tmppr,tmpmin,tmpne
        real(c_double) :: h(2),p_i(2),p_k(2),l(2),cj(2),rj,b(2),dev(3)
        integer(c_int) :: st

        ! LOCAL ARRAYS
        real(kind=8) :: c(2)

        c(1:2) = inst%centercoor(1:2,pol)
        r      = inst%radii(pol)

        ! Window cost at the incumbent (fills the incumbent angles)
        call psi0(zpp,zpr,zne,znn,zmin,initialpsi,anglepr,anglemin,anglene)

        ! Optimal point b of the window, by the closed-form solver
        h=zpp; p_i=zpr; p_k=zne; l=znn; cj=c; rj=r
        st = solve_chain_c(h,p_i,p_k,l,cj,rj,b,dev)

        if ( st .lt. 0 ) then                 ! invalid input -> keep incumbent
            diff = 0.0d0
            return
        end if

        ! Evaluates b with psi0 (keeps the cost bookkeeping consistent)
        call psi0(zpp,zpr,zne,znn,b,newpsi,tmppr,tmpmin,tmpne)

        if ( newpsi .lt. initialpsi + tol_acc ) then
            zmin(1:2) = b(1:2)
            anglepr = tmppr; anglemin = tmpmin; anglene = tmpne
            diff = newpsi - initialpsi
        else
            diff = 0.0d0
        end if

    end subroutine discmin_cf
