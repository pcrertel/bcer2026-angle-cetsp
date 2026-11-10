module instance

    implicit none

    type :: instance_type
        integer :: nballs,lenfilename
        real(kind=8) :: xymax                 ! scale factor of the normalization (1 if not normalized)
        real(kind=8) :: xyshift(2)            ! translation of the normalization (0 if not normalized)
        integer, allocatable :: nn_list(:,:)  ! disks sorted by distance from each disk
        real(kind=8), allocatable :: centercoor(:,:),radii(:),distmatrix(:,:),&
        min_angle_matrix(:,:,:)               ! lower bounds (angle-bound) of the paper
        character(len=256) :: filename        ! path of the instance file, as read
        character(len=128) :: name            ! name of the instance: file name without folders and extension
    end type instance_type

    type :: parameter_type
        integer :: bcdmaxcic,a,b,nworst,alnsmaxiter,segment_size,nab_cap,tri_c
        real(kind=8) :: seed,initial_seed,w,t_end,alpha,sigma1,sigma2,&
        sigma3,rfactor,p_split,p_blink,p_worst
    end type parameter_type

    type :: counter_type
        integer :: bcdcal,bcdupd,alnsiter,nimprov,last_improv_iter,max_no_improv,&
        n_visited,hashsize
        real(kind=8) :: costhc,last_improv_time
        integer :: remov_calls(11),remov_success(11),insert_calls(4),insert_success(4)
        integer(kind=8), allocatable :: hash_table(:)
    end type counter_type

    contains

    real(kind=8) function drand(ix)
        implicit none
        ! This is the random number generator of Schrage:
        !
        ! L. Schrage, A more portable Fortran random number generator, ACM
        ! Transactions on Mathematical Software 5 (1979), 132-138.
        ! SCALAR ARGUMENT
        real(kind=8), intent(inout) :: ix

        ! LOCAL ARRAYS
        real(kind=8) :: a,p,b15,b16,xhi,xalo,leftlo,fhi,k

        data a/16807.d0/,b15/32768.d0/,b16/65536.d0/,p/2147483647.d0/

        xhi= ix/b16
        xhi= xhi - dmod(xhi,1.d0)
        xalo= (ix-xhi*b16)*a
        leftlo= xalo/b16
        leftlo= leftlo - dmod(leftlo,1.d0)
        fhi= xhi*a + leftlo
        k= fhi/b15
        k= k - dmod(k,1.d0)
        ix= (((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
        if (ix.lt.0) ix= ix + p
        drand= ix*4.656612875d-10

        return
    end function drand

    ! *****************************************************************

    subroutine turningangle(p1,p2,p3,angle)

        implicit none

        ! SCALAR ARGUMENTS
        real(kind=8),intent(out) :: angle

        ! ARRAY ARGUMENTS
        real(kind=8),intent(in) :: p1(2),p2(2),p3(2)

        ! This routine computes the turning angle between segments (p1,p2) and (p2,p3).

        ! The angle between u = p2-p1 and v = p3-p2 is computed as
        ! atan2(|u x v|, <u,v>), which is accurate over the whole range [0,pi].
        ! The former acos(<u,v>/(|u||v|)) lost about eight digits near 0 and pi.

        ! LOCAL SCALARS
        real(kind=8) :: u(2),v(2)

        u = p2 - p1
        v = p3 - p2
        if ( norm2( u ) .le. 0.0d0 .or. norm2( v ) .le. 0.0d0 ) then
            angle = 0.0d0
            return
        end if

        angle = atan2( abs( u(1)*v(2) - u(2)*v(1) ), u(1)*v(1) + u(2)*v(2) )

    end subroutine turningangle

    ! *****************************************************************

    subroutine inputdata(inst,param,count,normalize)

        implicit none

        ! TYPE ARGUMENTS
        type(instance_type), intent(inout) :: inst
        type(parameter_type), intent(inout) :: param
        type(counter_type), intent(inout) :: count

        ! SCALAR ARGUMENTS
        logical, intent(in) :: normalize

        ! This routine reads the input parameters and the instance, and builds the
        ! geometric structures (distances, neighbor lists and lower bounds).

        ! LOCAL SCALARS
        integer :: i,j,k,allocerr,key_idx
        real(kind=8) :: key_dist,dist_ij,dist_jk,alpha,beta,theta_c,theta_min

        ! =====================================================================
        ! PARAMETERS (standard input, one per line; notation of Section 5.3 of the paper)
        ! =====================================================================
        read(*,*) inst%filename
        read(*,*) param%seed
        read(*,*) param%alnsmaxiter         ! N, number of iterations
        read(*,*) param%p_worst             ! beta, exponent of the ranking rule
        read(*,*) param%w                   ! chi, initial worsening accepted with probability 1/2
        read(*,*) param%t_end               ! T_end, final temperature
        read(*,*) param%rfactor             ! lambda, reaction factor
        read(*,*) param%segment_size        ! iterations per weight segment
        read(*,*) param%sigma1              ! score of a new best solution
        read(*,*) param%sigma2              ! score of an accepted solution that does not increase the cost
        read(*,*) param%sigma3              ! score of an accepted solution that increases the cost
        read(*,*) param%nab_cap             ! largest number of removed disks
        read(*,*) param%p_split             ! splitting probability of the SISR random split
        read(*,*) param%p_blink             ! skipping probability of the SISR greedy with blinks

        inst%filename = adjustl(inst%filename)
        inst%lenfilename = len_trim(inst%filename)
        ! name of the instance = file name without the folders and the extension, whatever the path
        block
            integer :: ib,ie
            ib = index(inst%filename(1:inst%lenfilename),'/',back=.true.) + 1
            ie = index(inst%filename(ib:inst%lenfilename),'.',back=.true.)
            if ( ie .gt. 1 ) then
                ie = ib + ie - 2
            else
                ie = inst%lenfilename
            end if
            inst%name = inst%filename(ib:ie)
        end block
        param%initial_seed = param%seed
        param%bcdmaxcic = 1
        count%hashsize = 500009_8

        open(unit=10,file=trim(adjustl(inst%filename)),status='old')
        read(10,*) inst%nballs
        write(*,*) 'nballs = ',inst%nballs

        allocate(inst%centercoor(2,inst%nballs),inst%radii(inst%nballs),&
        inst%distmatrix(inst%nballs,inst%nballs),inst%nn_list(inst%nballs,inst%nballs),&
        count%hash_table(count%hashsize),stat=allocerr)

        if (allocerr.ne.0) then
            write(*,*) 'Allocation error - inputdata.'
            stop
        end if

        do i=1,inst%nballs
            read(10,*) inst%centercoor(1:2,i),inst%radii(i)
        end do
        close(10)

        ! the number of removed disks is drawn from a triangular distribution on [a,b]
        ! with mode tri_c
        param%a = 1
        param%tri_c = 1
        if ( param%nab_cap .gt. 0 ) then
            param%b = min( inst%nballs, param%nab_cap )
        else
            param%b = inst%nballs
        end if
        if ( param%b .lt. 1 ) param%b = 1

        param%nworst=param%b

        ! Normalization: translates the lower left corner of the rectangle that contains
        ! the centers to the origin and divides the coordinates AND THE RADII by the
        ! largest side of this rectangle, so that the centers lie in [0,1]^2. The
        ! objective (sum of the turning angles) is invariant under translation and
        ! scaling, and the absolute tolerances of the BCD hold on the same scale for
        ! all instances.
        if (normalize) then
            inst%xyshift(1)=minval(inst%centercoor(1,1:inst%nballs))
            inst%xyshift(2)=minval(inst%centercoor(2,1:inst%nballs))
            inst%xymax=max(maxval(inst%centercoor(1,1:inst%nballs))-inst%xyshift(1),&
                           maxval(inst%centercoor(2,1:inst%nballs))-inst%xyshift(2))
            if (inst%xymax.le.0.0d0) inst%xymax=1.0d0
            do j=1,inst%nballs
                inst%centercoor(1:2,j)=(inst%centercoor(1:2,j)-inst%xyshift(1:2))/inst%xymax
                inst%radii(j)=inst%radii(j)/inst%xymax
            end do
        else
            inst%xymax=1.0d0
            inst%xyshift=0.0d0
        end if

        do i=1,inst%nballs
            inst%distmatrix(i,i)=0.d0
            do j=i+1,inst%nballs
                inst%distmatrix(i,j)=norm2(inst%centercoor(1:2,i)-inst%centercoor(1:2,j))
                inst%distmatrix(j,i)=inst%distmatrix(i,j)
            end do
        end do 

        do i=1,inst%nballs
            do j=1,inst%nballs
                inst%nn_list(i,j)=j
            end do
            do j=2,inst%nballs
                key_idx=inst%nn_list(i,j)
                key_dist=inst%distmatrix(i,key_idx)
                k=j-1
                do while (k.ge.1)
                    if (inst%distmatrix(i,inst%nn_list(i,k)).le.key_dist) exit
                    inst%nn_list(i,k+1)=inst%nn_list(i,k)
                    k=k-1
                end do
                inst%nn_list(i,k+1)=key_idx
            end do
        end do

        allocate(inst%min_angle_matrix(inst%nballs,inst%nballs,inst%nballs),stat=allocerr)
        if (allocerr/=0) then
            write(*,*) 'Allocation error - min_angle_matrix'
            stop
        end if
        
        inst%min_angle_matrix=0.0d0 

        do i=1,inst%nballs
            do j=1,inst%nballs
                if (i.eq.j) cycle
                do k=i+1,inst%nballs
                    if (k.eq.j) cycle

                    dist_ij=inst%distmatrix(i,j)
                    dist_jk=inst%distmatrix(j,k)

                    if (inst%radii(i)+inst%radii(j).ge.dist_ij) then
                        alpha=acos(-1.0d0) 
                    else
                        alpha=asin((inst%radii(i)+inst%radii(j))/dist_ij)
                    end if

                    if (inst%radii(j)+inst%radii(k).ge.dist_jk) then
                        beta=acos(-1.0d0) 
                    else
                        beta=asin((inst%radii(j)+inst%radii(k))/dist_jk)
                    end if

                    call turningangle(inst%centercoor(1:2,i),inst%centercoor(1:2,j),inst%centercoor(1:2,k),theta_c)
                    theta_min=max(0.0d0,theta_c-alpha-beta)

                    inst%min_angle_matrix(i,j,k)=theta_min
                    inst%min_angle_matrix(k,j,i)=theta_min
                end do
            end do
        end do

        count%bcdcal=0
        count%bcdupd=0
        count%last_improv_time=0.0d0
        count%alnsiter=0
        count%nimprov=0
        count%remov_calls=0
        count%insert_calls=0
        count%remov_success=0
        count%insert_success=0
        count%last_improv_iter=0
        count%max_no_improv=0
        count%n_visited=0
        count%hash_table=-1_8
    end subroutine inputdata

end module instance
