! ==========================================================================
! bcd_compare.f90 -- block selection rules of the BCD (Section 5.1 of the paper)
! on a fixed visiting sequence. From the same starting points (the centers of
! the disks) and with the same closed-form block update, it runs
!
!   cyclic : angleBCDcyclic, which sweeps the visit points in the order of the
!            route and stops when a cycle does not decrease the cost by a
!            relative amount of sqrt(machine epsilon);
!   queue  : angleBCD, the queue rule used in the ALNS-BCD.
!
! The visiting sequence is the order of the disks in the instance file, and
! no cap on the number of block updates is used. Each rule is run NREP times
! from the same start, and the reported time is the smallest of the NREP runs
! (the solutions are identical).
!
! Usage:   ./bcd_compare <instance file>
! Output:  one line
!   m  psi_cyclic  psi_queue  updates_cyclic  updates_queue
!      cycles_cyclic  cycles_queue  time_cyclic_s  time_queue_s
! (updates = block updates; cycles = updates/m, rounded up for the queue)
! ==========================================================================
program bcd_compare

    use instance, only: instance_type, parameter_type, counter_type
    use blockcoordinatedescent, only: angleBCD, angleBCDcyclic, evalfangle

    implicit none

    integer, parameter :: NREP = 5
    type(instance_type) :: inst
    type(parameter_type) :: param
    type(counter_type) :: count

    integer :: m, i, allocerr, rule, rep
    integer :: upd(2), cyc(2)
    real(kind=8) :: psi(2), tbest(2), t0, t1, c
    character(len=256) :: fname

    integer, allocatable :: tour(:)
    real(kind=8), allocatable :: x(:,:), xc(:,:), angles(:)

    call get_command_argument(1,fname)
    open(unit=20,file=trim(fname),status='old')
    read(20,*) m
    allocate(inst%centercoor(2,m),inst%radii(m),tour(m),x(2,m),xc(2,m), &
             angles(m),stat=allocerr)
    do i = 1, m
        read(20,*) inst%centercoor(1,i),inst%centercoor(2,i),inst%radii(i)
    end do
    close(20)
    inst%nballs = m

    param%bcdmaxcic = 10000                           ! no cap in practice

    do i = 1, m
        tour(i) = i                                   ! sequence = order of the file
        xc(1:2,i) = inst%centercoor(1:2,i)            ! start at the centers
    end do

    do rule = 1, 2
        tbest(rule) = huge(1.0d0)
        do rep = 1, NREP
            count%bcdcal = 0
            count%bcdupd = 0
            x = xc
            call evalfangle(m,tour,x,c,angles)
            call cpu_time(t0)
            if ( rule .eq. 1 ) then
                call angleBCDcyclic(inst,param,count,m,tour,x,c,angles)
            else
                call angleBCD(inst,param,count,m,tour,x,c,angles)
            end if
            call cpu_time(t1)
            tbest(rule) = min( tbest(rule), t1 - t0 )
        end do
        psi(rule) = c
        upd(rule) = count%bcdupd
        cyc(rule) = ( count%bcdupd + m - 1 ) / m
    end do

    write(*,'(I6,2(1X,F18.10),2(1X,I9),2(1X,I7),2(1X,ES12.5))') &
        m, psi(1), psi(2), upd(1), upd(2), cyc(1), cyc(2), tbest(1), tbest(2)

    deallocate(inst%centercoor,inst%radii,tour,x,xc,angles)

end program bcd_compare
