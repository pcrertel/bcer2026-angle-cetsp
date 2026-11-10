module blockcoordinatedescent

    use instance, only: instance_type,parameter_type,counter_type,turningangle

    implicit none

    ! PARAMETERS
    real(kind=8), parameter :: pi = acos(-1.0d0), macheps = epsilon( 1.0d0 ),&
    macheps12 = sqrt( macheps ), tolcoinc = 1.0d-12

    contains

    integer function nextind(i,imin,imax)

        implicit none

        integer, intent(in) :: i,imin,imax

        if (i .lt. imax) then
            nextind = i + 1
        else
            nextind = imin
        end if

        return
    end function nextind

    ! *****************************************************************

    integer function prevind(i,imin,imax)

        implicit none

        integer, intent(in) :: i,imin,imax

        if (i .gt. imin) then
            prevind = i - 1
        else
            prevind = imax
        end if

        return
    end function prevind

    ! *****************************************************************

    ! =================================================================
    ! angleBCD  -  BCD with the queue block selection rule (Section 5.1 of the paper).
    ! A FIFO queue holds the positions whose window may still change; each position
    ! enters the queue at most once at a time (flags inq). The cost is updated in O(1)
    ! through diff. After the update of position p, if its visit point moved
    ! (displacement > macheps12*max(1,|z|)), even without decreasing the cost (ties of
    ! the segment case), the 5 positions p-2..p+2 are put back in the queue. The method
    ! stops when the queue is empty, when the cost reaches 2*pi, after npts consecutive
    ! updates without a relative decrease of macheps12 with respect to the value fref at
    ! the last progress, or after bcdmaxcic*npts block updates. At the end, the cost and
    ! the angles are recomputed with evalfangle.
    ! =================================================================
    subroutine angleBCD(inst,param,count,npts,tour,x,cost,angles)

        implicit none

        ! TYPE ARGUMENTS
        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(in) :: param
        type(counter_type), intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer, intent(in) :: npts
        real(kind=8), intent(inout) :: cost

        ! ARRAY ARGUMENTS
        integer, intent(in) :: tour(inst%nballs)
        real(kind=8), intent(inout) :: x(2,inst%nballs),angles(inst%nballs)

        call angleBCDqueue(inst,param,count,npts,tour,x,cost,angles)
        call evalfangle(npts,tour,x,cost,angles)

    end subroutine angleBCD

    ! *****************************************************************

    subroutine angleBCDqueue(inst,param,count,npts,tour,x,cost,angles)

        implicit none

        ! TYPE ARGUMENTS
        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(in) :: param
        type(counter_type), intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer, intent(in) :: npts
        real(kind=8), intent(inout) :: cost

        ! ARRAY ARGUMENTS
        integer, intent(in) :: tour(inst%nballs)
        real(kind=8), intent(inout) :: x(2,inst%nballs),angles(inst%nballs)

        ! LOCAL SCALARS
        integer :: p,v,d,kprev1,kprev2,knext1,knext2,head,tail,nq,nupd,noprogress
        real(kind=8) :: diff,curcost,fref

        ! LOCAL ARRAYS
        real(kind=8) :: zpp(2),zpr(2),zne(2),znn(2),zold(2)
        integer :: queue(npts)
        logical :: inq(npts)

        count%bcdcal = count%bcdcal + 1

        ! a route whose cost is 2*pi is already optimal (convex polygon)
        if ( abs( cost - 2.0d0 * pi) .le. macheps ) return

        curcost = cost
        fref = cost          ! value at the last progress
        noprogress = 0       ! consecutive updates without progress
        nupd = 0             ! block updates performed

        ! the queue starts with all positions, in the order of the tour
        do p = 1, npts
            queue(p) = p
            inq(p)   = .true.
        end do
        head = 1 ; tail = npts ; nq = npts

        do while ( nq > 0 )
            p = queue(head) ; head = mod(head,npts) + 1 ; nq = nq - 1
            inq(p) = .false.

            kprev1 = prevind(p     ,1,npts)
            kprev2 = prevind(kprev1,1,npts)
            knext1 = nextind(p     ,1,npts)
            knext2 = nextind(knext1,1,npts)

            zpp(1:2) = x(1:2,tour(kprev2))
            zpr(1:2) = x(1:2,tour(kprev1))
            zne(1:2) = x(1:2,tour(knext1))
            znn(1:2) = x(1:2,tour(knext2))
            zold(1:2) = x(1:2,tour(p))

            ! closed-form block update
            call discmin(inst,param,count,tour(p),zpp,zpr,zne,znn,&
            x(1:2,tour(p)),diff,angles(tour(kprev1)),&
            angles(tour(p)),angles(tour(knext1)))

            curcost = curcost + diff
            nupd = nupd + 1

            ! if the visit point moved, improving or not, the 5 windows p-2..p+2 go back
            ! to the queue
            if ( norm2( x(1:2,tour(p)) - zold(1:2) ) .gt. &
                 macheps12 * max( 1.0d0, norm2( zold(1:2) ) ) ) then
                do d = -2, 2
                    v = modulo(p-1+d, npts) + 1
                    if ( .not. inq(v) ) then
                        inq(v) = .true.
                        tail = mod(tail,npts) + 1
                        queue(tail) = v
                        nq = nq + 1
                    end if
                end do
            end if

            ! stopping criteria
            if ( curcost .lt. fref * ( 1.0d0 - macheps12 ) ) then
                fref = curcost
                noprogress = 0
            else
                noprogress = noprogress + 1
            end if

            if ( abs( curcost - 2.0d0 * pi ) .le. macheps .or. noprogress .ge. npts .or. &
                 nupd .ge. param%bcdmaxcic * npts ) then
                cost = curcost
                count%bcdupd = count%bcdupd + nupd
                return
            end if
        end do

        ! empty queue: no visit point moves
        cost = curcost
        count%bcdupd = count%bcdupd + nupd

    end subroutine angleBCDqueue

    ! *****************************************************************

    ! Cyclic BCD, used only in the comparison of the block selection rules of Section 5.1
    ! (program bcd_compare). The cost and the angles are updated in O(1) per block with
    ! the local rule of psi0, and recomputed with evalfangle at the end. The method stops
    ! when a cycle fails to decrease the cost by a relative amount of macheps12, when the
    ! cost reaches 2*pi, or after bcdmaxcic cycles.
    subroutine angleBCDcyclic(inst,param,count,npts,tour,x,cost,angles)

        implicit none

        ! TYPE ARGUMENTS
        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(in) :: param
        type(counter_type), intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer, intent(in) :: npts
        real(kind=8), intent(inout) :: cost

        ! ARRAY ARGUMENTS
        integer, intent(in) :: tour(inst%nballs)
        real(kind=8), intent(inout) :: x(2,inst%nballs),angles(inst%nballs)

        ! LOCAL SCALARS
        integer :: k,kprev1,knext1,kprev2,knext2,cic
        real(kind=8) :: prevcost,curcost,diff

        ! LOCAL ARRAYS
        real(kind=8) :: zpp(2),zpr(2),zne(2),znn(2),xprev(2,inst%nballs),&
        prevangles(inst%nballs)

        count%bcdcal = count%bcdcal + 1

        if ( abs( cost - 2.0d0 * pi) .gt. macheps ) then

            cic = 0
            curcost = cost
            prevcost = cost
            xprev = x
            prevangles = angles

            do
                cic = cic + 1

                ! one cycle: all visit points, in the order of the tour
                do k = 1, npts
                    kprev1 = prevind(k     ,1,npts)
                    kprev2 = prevind(kprev1,1,npts)
                    knext1 = nextind(k     ,1,npts)
                    knext2 = nextind(knext1,1,npts)

                    zpp(1:2) = x(1:2,tour(kprev2))
                    zpr(1:2) = x(1:2,tour(kprev1))
                    zne(1:2) = x(1:2,tour(knext1))
                    znn(1:2) = x(1:2,tour(knext2))

                    call discmin(inst,param,count,tour(k),zpp,zpr,zne,znn,&
                    x(1:2,tour(k)),diff,angles(tour(kprev1)),&
                    angles(tour(k)),angles(tour(knext1)))

                    curcost = curcost + diff
                end do
                count%bcdupd = count%bcdupd + npts

                ! a cycle that increases the cost is undone
                if ( curcost .gt. prevcost * ( 1.0d0 + macheps12 ) ) then
                    x = xprev
                    angles = prevangles
                    cost = prevcost
                    exit
                end if

                if ( abs( curcost - 2.0d0 * pi) .le. macheps .or. &
                     curcost .ge. prevcost * ( 1.0d0 - macheps12 ) .or. &
                     cic .ge. param%bcdmaxcic ) then
                    cost = curcost
                    exit
                end if

                prevcost = curcost
                xprev = x
                prevangles = angles
            end do
        end if

        call evalfangle(npts,tour,x,cost,angles)

    end subroutine angleBCDcyclic

    ! *****************************************************************

    ! closed-form solution of the single-disk subproblem (see discmin_cf-param.f90)
    subroutine discmin(inst,param,count,pol,zpp,zpr,zne,znn,zmin,diff,&
    anglepr,anglemin,anglene)
        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(in) :: param
        type(counter_type),intent(inout) :: count
        integer,intent(in) :: pol
        real(kind=8),intent(inout) :: anglepr,anglemin,anglene
        real(kind=8),intent(out) :: diff
        real(kind=8),intent(in) :: zpp(2),zpr(2),zne(2),znn(2)
        real(kind=8),intent(inout) :: zmin(2)

        call discmin_cf(inst,param,count,pol,zpp,zpr,zne,znn,zmin,diff,&
        anglepr,anglemin,anglene)
    end subroutine discmin

    ! *****************************************************************

    subroutine evalfangle(npts,tour,x,fangle,angles)

        implicit none

        ! SCALAR ARGUMENTS
        integer, intent(in) :: npts
        real(kind=8), intent(out) :: fangle

        ! ARRAY ARGUMENTS
        integer, intent(in) :: tour(npts)
        real(kind=8), intent(in) :: x(2,*)
        real(kind=8), intent(out) :: angles(*)

        ! This routine must compute the objective function of the
        ! angle-CETSP problem, i.e., the sum of direction changes
        ! angles of a given tour x.
        !
        ! Coincident visit points. A point within
        ! tolcoinc*max(1,|a|,|b|) of the first point of the current group is the
        ! same location: the corresponding disks are visited at the same point.
        ! The cost is the total turning angle of the closed polygon whose
        ! vertices are the first points of the groups; the turning angle of a
        ! group is assigned to its first point, and the others receive zero.
        ! Since the same vertices are used on both sides of every turn, the
        ! polygon closes and the cost is never below 2*pi. The BCD and the ALNS
        ! therefore recompute the cost and the angles with this routine at the
        ! end, instead of relying on the incremental bookkeeping. If all points
        ! coincide, the cost is 2*pi, the limit of closed routes with distinct
        ! points.

        ! LOCAL SCALARS
        integer :: k,nv
        real(kind=8) :: alpha

        ! LOCAL ARRAYS
        integer :: v(npts)

        fangle = 0.0d0
        do k = 1, npts
            angles(tour(k)) = 0.0d0
        end do

        ! v(1:nv): the disks whose points are the vertices of the polygon
        nv = 1
        v(1) = tour(1)
        do k = 2, npts
            if ( .not. same_point( x(1:2,tour(k)), x(1:2,v(nv)) ) ) then
                nv = nv + 1
                v(nv) = tour(k)
            end if
        end do
        do while ( nv .gt. 1 )
            if ( .not. same_point( x(1:2,v(1)), x(1:2,v(nv)) ) ) exit
            nv = nv - 1
        end do

        if ( nv .eq. 1 ) then
            fangle = 2.0d0 * pi
            angles(tour(1)) = fangle
            return
        end if

        do k = 1, nv
            call turningangle(x(1:2,v(prevind(k,1,nv))),x(1:2,v(k)),x(1:2,v(nextind(k,1,nv))),alpha)
            fangle = fangle + alpha
            angles(v(k)) = alpha
        end do

    end subroutine evalfangle

    ! *****************************************************************

    logical function same_point(a,b)

        implicit none

        real(kind=8), intent(in) :: a(2),b(2)

        same_point = norm2( a - b ) .le. tolcoinc * max( 1.0d0, norm2( a ), norm2( b ) )

    end function same_point

    ! *****************************************************************

    subroutine psi0(pprev,p,q,qnext,z,psival,anglep,anglez,angleq)

        implicit none

        ! SCALAR ARGUMENTS
        real(kind=8), intent(out) :: psival,anglep,anglez,angleq

        ! ARRAY ARGUMENTS
        real(kind=8), intent(in) :: pprev(2),p(2),q(2),qnext(2),z(2)

        ! This routine must compute the total turning angle of the polyline 
        ! pprev(2) → p(2) → z(2) → q(2) → qnext(2).

        if ( same_point( p(1:2), z(1:2) ) .or. same_point( q(1:2), z(1:2) ) ) then
            ! z = p or z = q
            call turningangle(pprev,p,q,anglep)
            anglez = 0.0d0
            call turningangle(p,q,qnext,angleq)
        else
            call turningangle(pprev,p,z,anglep)
            call turningangle(p,z,q,anglez)
            call turningangle(z,q,qnext,angleq)
        end if

        psival = anglep + anglez + angleq
    end subroutine psi0

    ! *****************************************************************

    include 'discmin_cf-param.f90'

end module blockcoordinatedescent
