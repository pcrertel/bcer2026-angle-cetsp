module metaheuristic

    use instance, only: instance_type,parameter_type,counter_type,drand,turningangle
    use blockcoordinatedescent, only: angleBCD,discmin,nextind,prevind,evalfangle

    implicit none

    ! Removal operators: 1 Worst angles, 2 Random, 3 Cluster, 4 Random string, 5 Bad string,
    ! 6 Good string, 7 Worst-sum string, 8 Savings, 9 SISR random split, 10 SISR spatial,
    ! 11 Lens. Insertion operators: 1 Greedy, 2 Regret-2, 3 Regret-3, 4 SISR greedy with
    ! blinks. The strings removed by operators 4 and 6 are reinserted as single elements.
    integer, parameter :: num_remov_ops = 11, num_insert_ops = 4

    ! strings reinserted in both directions / times in which the reversed direction was better
    ! (see subroutine insert), for the report at the end of the ALNS
    integer, save :: rev_tries = 0, rev_chosen = 0

    ! orders used to reinsert the disks in the SISR insertion (see sisr_reorder)
    integer, save :: sisr_orders(3) = 0

    contains

    ! =================================================================
    ! initial_solution -- initial solution of the ALNS, by convex hull + cheapest
    ! insertion (the Cheapest Insertion heuristic of Fischer et al. 2014 for the
    ! Angle-TSP), with the cost of each insertion computed as in the insertion
    ! operators of the ALNS:
    !   1. the route starts with the vertices of the convex hull of the centers,
    !      visited at the centers;
    !   2. for each disk w outside the route and each edge u -> v of the route, the
    !      visit point of w is the minimizer of the single-disk subproblem (discmin,
    !      closed form), with p = x(u), q = x(v), p- = x(prev(u)) and q+ = x(next(v));
    !      the cost of the insertion is the increase of the turning angles at u, w and
    !      v. The pair (w, edge) of smallest cost is inserted;
    !   3. after w is inserted, the visit points of u and v are updated by one block
    !      update each (as in the insertion operators);
    !   4. at the end, the BCD runs until convergence (bcdmaxcic = 200).
    ! With radius zero the visit point is the center itself, and the construction is
    ! the classical cheapest insertion. Each disk outside the route keeps its best
    ! edge. An insertion between u and v moves the points of u, w and v, which changes
    ! the cost of the edges that start at prev(prev(u)), prev(u), u, w, v and next(v).
    ! The disks whose best edge was one of these scan the route again, and the others
    ! compare their best edge with these edges.
    ! =================================================================
    subroutine initial_solution(inst,param,count,tour,x,cost,angles)

        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        type(counter_type),intent(inout) :: count
        real(kind=8),intent(out) :: cost
        integer,intent(out) :: tour(inst%nballs)
        real(kind=8),intent(out) :: x(2,inst%nballs),angles(inst%nballs)

        integer :: m,nh,k,i,j,w,wbest,u,v,maxcic0,nin,a,affected(6)
        integer :: hull(inst%nballs),nxt(inst%nballs),prv(inst%nballs),bestu(inst%nballs)
        logical :: inroute(inst%nballs),rescan
        real(kind=8) :: bestd(inst%nballs),d,dmin,zw(2)
        real(kind=8) :: c(2,inst%nballs),xx(2,inst%nballs)

        m = inst%nballs
        c(1:2,1:m) = inst%centercoor(1:2,1:m)
        xx = c
        maxcic0 = param%bcdmaxcic

        call convex_hull(m,c,hull,nh)

        inroute = .false.
        do k = 1, nh
            i = hull(k)
            j = hull(modulo(k,nh)+1)
            nxt(i) = j
            prv(j) = i
            inroute(i) = .true.
        end do
        nin = nh

        do w = 1, m
            if ( .not. inroute(w) ) call scan(w)
        end do

        do while ( nin .lt. m )
            wbest = 0
            dmin = huge(1.0d0)
            do w = 1, m
                if ( inroute(w) ) cycle
                if ( bestd(w) .lt. dmin ) then
                    dmin = bestd(w)
                    wbest = w
                end if
            end do
            w = wbest
            u = bestu(w)
            v = nxt(u)
            d = insertion_cost(w,u,zw)          ! visit point of w on this edge
            nxt(u) = w; prv(w) = u
            nxt(w) = v; prv(v) = w
            xx(1:2,w) = zw
            inroute(w) = .true.
            nin = nin + 1
            ! one block update at u and another at v
            call update_point(u)
            call update_point(v)
            ! edges whose cost changed
            affected(1) = prv(prv(u)); affected(2) = prv(u); affected(3) = u
            affected(4) = w;           affected(5) = v;      affected(6) = nxt(v)
            do i = 1, m
                if ( inroute(i) ) cycle
                rescan = any( affected .eq. bestu(i) )
                if ( rescan ) then
                    call scan(i)
                else
                    do k = 1, 6
                        a = affected(k)
                        d = insertion_cost(i,a,zw)
                        if ( d .lt. bestd(i) ) then
                            bestd(i) = d
                            bestu(i) = a
                        end if
                    end do
                end if
            end do
        end do

        tour(1) = hull(1)
        do k = 2, m
            tour(k) = nxt(tour(k-1))
        end do

        x(1:2,1:m) = xx(1:2,1:m)
        call evalfangle(m,tour,x,cost,angles)
        param%bcdmaxcic = 200
        call angleBCD(inst,param,count,m,tour,x,cost,angles)
        param%bcdmaxcic = maxcic0

    contains

        ! cost of inserting w on the edge a -> nxt(a), with the visit point of w given by
        ! the single-disk subproblem; also returns this point in z
        real(kind=8) function insertion_cost(w,a,z)
            integer, intent(in) :: w,a
            real(kind=8), intent(out) :: z(2)
            integer :: b,pa,qb
            real(kind=8) :: dif,apr,amin,ane,t4,t5
            b = nxt(a); pa = prv(a); qb = nxt(b)
            z = c(1:2,w)
            apr = 0.0d0; amin = 0.0d0; ane = 0.0d0
            call discmin(inst,param,count,w,xx(1:2,pa),xx(1:2,a),xx(1:2,b),xx(1:2,qb),&
            z,dif,apr,amin,ane)
            call turningangle(xx(1:2,pa),xx(1:2,a),xx(1:2,b),t4)
            call turningangle(xx(1:2,a),xx(1:2,b),xx(1:2,qb),t5)
            insertion_cost = apr + amin + ane - t4 - t5
        end function insertion_cost

        ! one block update of the visit point of disk a of the route
        subroutine update_point(a)
            integer, intent(in) :: a
            real(kind=8) :: dif,apr,amin,ane
            apr = 0.0d0; amin = 0.0d0; ane = 0.0d0
            call discmin(inst,param,count,a,xx(1:2,prv(prv(a))),xx(1:2,prv(a)),&
            xx(1:2,nxt(a)),xx(1:2,nxt(nxt(a))),xx(1:2,a),dif,apr,amin,ane)
        end subroutine update_point

        ! full scan: best edge of the current route for disk w
        subroutine scan(w)
            integer, intent(in) :: w
            integer :: a0,a1
            real(kind=8) :: dd,zz(2)
            bestd(w) = huge(1.0d0)
            bestu(w) = 0
            a0 = hull(1)
            a1 = a0
            do
                dd = insertion_cost(w,a1,zz)
                if ( dd .lt. bestd(w) ) then
                    bestd(w) = dd
                    bestu(w) = a1
                end if
                a1 = nxt(a1)
                if ( a1 .eq. a0 ) exit
            end do
        end subroutine scan

    end subroutine initial_solution

    ! =================================================================
    ! convex_hull -- vertices of the convex hull of the points c(:,1:m), in
    ! counterclockwise order (Andrew's monotone chain), without collinear or repeated
    ! points. If all points are collinear, returns the two extreme points.
    ! =================================================================
    subroutine convex_hull(m,c,hull,nh)

        implicit none

        integer, intent(in) :: m
        real(kind=8), intent(in) :: c(2,m)
        integer, intent(out) :: hull(m),nh

        integer :: ord(m),k,i,j,tmp,lo
        integer :: h(2*m+1)

        ! sorts the indices by (x,y) (insertion sort; m is at most of the order of 10^3)
        do i = 1, m
            ord(i) = i
        end do
        do i = 2, m
            tmp = ord(i)
            j = i - 1
            do while ( j .ge. 1 )
                if ( lex_less(tmp, ord(j)) ) then
                    ord(j+1) = ord(j)
                    j = j - 1
                else
                    exit
                end if
            end do
            ord(j+1) = tmp
        end do

        k = 0
        do i = 1, m                                   ! lower chain
            do while ( k .ge. 2 )
                if ( cross(h(k-1), h(k), ord(i)) .gt. 0.0d0 ) exit
                k = k - 1
            end do
            k = k + 1
            h(k) = ord(i)
        end do
        lo = k + 1
        do i = m - 1, 1, -1                           ! upper chain
            do while ( k .ge. lo )
                if ( cross(h(k-1), h(k), ord(i)) .gt. 0.0d0 ) exit
                k = k - 1
            end do
            k = k + 1
            h(k) = ord(i)
        end do
        nh = k - 1                                    ! the last one repeats the first
        if ( nh .lt. 2 ) then
            nh = min( m, 2 )
            hull(1) = ord(1)
            if ( m .ge. 2 ) hull(2) = ord(m)
        else
            hull(1:nh) = h(1:nh)
        end if

    contains

        logical function lex_less(a,b)
            integer, intent(in) :: a,b
            lex_less = c(1,a) .lt. c(1,b) .or. ( c(1,a) .eq. c(1,b) .and. c(2,a) .lt. c(2,b) )
        end function lex_less

        real(kind=8) function cross(o,a,b)
            integer, intent(in) :: o,a,b
            cross = (c(1,a)-c(1,o))*(c(2,b)-c(2,o)) - (c(2,a)-c(2,o))*(c(1,b)-c(1,o))
        end function cross

    end subroutine convex_hull

    subroutine alns(inst,param,count,tour,x,cost,angles,start,nworst,alnsmaxiter)

        implicit none

        ! TYPE ARGUMENTS
        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        type(counter_type),intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer,intent(in) :: nworst,alnsmaxiter
        real(kind=8),intent(in) :: start
        real(kind=8),intent(inout) :: cost

        ! ARRAY ARGUMENTS
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(inout) :: x(2,inst%nballs),angles(inst%nballs)

        ! LOCAL SCALARS
        integer :: nab,npts,op_remov,op_insert,i,cycle_length
        real(kind=8) :: currentcost,bestcost,newcost,r,T,To,To_curr,awarded_score,sum_p
        logical :: is_string

        ! LOCAL ARRAYS
        integer :: bestour(inst%nballs),currenttour(inst%nballs),newtour(inst%nballs),&
        removed_cities(nworst),use_remov(num_remov_ops),use_insert(num_insert_ops)
        real(kind=8) :: bestsol(2,inst%nballs),currentsol(2,inst%nballs),newsol(2,inst%nballs),&
        bestangles(inst%nballs),currentangles(inst%nballs),newangles(inst%nballs),&
        w_remov(num_remov_ops),w_insert(num_insert_ops),score_remov(num_remov_ops),&
        score_insert(num_insert_ops)

        ! Initialization of the incumbent solution
        bestour = tour
        bestsol = x
        bestcost = cost
        bestangles = angles

        currenttour = tour
        currentsol = x
        currentcost = cost
        currentangles = angles

        count%alnsiter = 0

        ! Initial temperature: a worsening of chi times the initial cost is accepted with
        ! probability 1/2
        To = - ( param%w * count%costhc ) / log( 0.5d0 )
        T = To
        To_curr = To

        ! Three temperature cycles of ceil(N/3) iterations; in each cycle the temperature is
        ! multiplied by alpha at every iteration so that it reaches t_end at the end of the cycle
        cycle_length = ( alnsmaxiter + 2 ) / 3
        param%alpha = exp( log( param%t_end / To ) / dble(cycle_length) )

        ! All operators start with weight 1
        w_remov = 1.0d0
        w_insert = 1.0d0
        score_remov = 0.0d0
        score_insert = 0.0d0
        use_remov = 0
        use_insert = 0

        write(*,*) '============ SIMULATED ANNEALING ==========='
        write(*,*) 'T_start     : ',T
        write(*,*) 'T_end       : ',param%t_end
        write(*,*) 'Alpha       : ',param%alpha
        write(*,*) 'Cycle length: ',cycle_length
        write(*,*) '============================================'

        do while ( count%alnsiter .lt. alnsmaxiter )

            count%alnsiter = count%alnsiter + 1

            ! The new solution always starts from the current solution
            newtour = currenttour
            newsol = currentsol
            newcost = currentcost
            newangles = currentangles

            ! Removal operator selection (roulette wheel)
            r = drand(param%seed)
            sum_p = sum( w_remov )
            r = r * sum_p
            sum_p = 0.0d0
            op_remov = num_remov_ops
            do i = 1, num_remov_ops
                sum_p = sum_p + w_remov(i)
                if ( r .le. sum_p ) then
                    op_remov = i
                    exit
                end if
            end do

            ! Insertion operator selection (roulette wheel)
            r = drand(param%seed)
            sum_p = sum( w_insert )
            r = r * sum_p
            sum_p = 0.0d0
            op_insert = num_insert_ops
            do i = 1, num_insert_ops
                sum_p = sum_p + w_insert(i)
                if ( r .le. sum_p ) then
                    op_insert = i
                    exit
                end if
            end do

            use_remov(op_remov) = use_remov(op_remov) + 1
            use_insert(op_insert) = use_insert(op_insert) + 1

            ! Removal phase
            call removal(inst,param,count,npts,nab,nworst,is_string,&
            newtour,newsol,newangles,removed_cities,op_remov)

            ! Reinsertion phase
            call insert(inst,param,count,npts,nab,nworst,newcost,is_string,&
            newtour,newsol,newangles,removed_cities,op_insert)

            ! Acceptance phase
            call acceptsol(inst,param,count,start,newcost,newtour,newsol,newangles,&
            currentcost,currenttour,currentsol,currentangles,bestcost,bestour,&
            bestsol,bestangles,op_remov,op_insert,T,awarded_score)

            ! Accumulates the scores of the operators in the current segment
            score_remov(op_remov) = score_remov(op_remov) + awarded_score
            score_insert(op_insert) = score_insert(op_insert) + awarded_score

            ! Adaptive weight update at the end of a segment
            if ( mod( count%alnsiter, param%segment_size ) .eq. 0 ) then
                do i = 1, num_remov_ops
                    if ( use_remov(i) .gt. 0 ) then
                        w_remov(i) = w_remov(i) * ( 1.0d0 - param%rfactor ) + &
                        param%rfactor * ( score_remov(i) / dble( use_remov(i) ) )
                    end if
                    w_remov(i) = max(w_remov(i), 0.001d0)
                end do
                do i = 1, num_insert_ops
                    if ( use_insert(i) .gt. 0 ) then
                        w_insert(i) = w_insert(i) * ( 1.0d0 - param%rfactor ) + &
                        param%rfactor * ( score_insert(i) / dble( use_insert(i) ) )
                    end if
                    w_insert(i) = max(w_insert(i), 0.001d0)
                end do
                score_remov = 0.0d0
                score_insert = 0.0d0
                use_remov = 0
                use_insert = 0
            end if

            ! Restart of the temperature at the end of the first and second cycles: the peak
            ! is halved, and the cooling factor is recomputed to reach t_end at the end of the
            ! next cycle
            if ( mod( count%alnsiter, cycle_length ) .eq. 0 .and. count%alnsiter .lt. alnsmaxiter ) then
                To_curr = To_curr * 0.5d0
                T = To_curr
                param%alpha = exp( log( param%t_end / To_curr ) / dble(cycle_length) )
                write(*,*) '>>> REHEAT AT ITERATION', count%alnsiter, 'NEW TEMPERATURE:', T
                call flush(6)
            end if

            ! Cooling
            T = max( T * param%alpha, param%t_end )
        end do

        write(*,*) ''
        write(*,*) '==================================================='
        write(*,*) '            ALNS OPERATOR STATISTICS               '
        write(*,*) '==================================================='
        write(*,*) 'REMOVAL OPERATORS (Calls / New best solutions):'
        write(*,*) '1 (Worst angles)       : ', count%remov_calls(1), ' / ', count%remov_success(1)
        write(*,*) '2 (Random)             : ', count%remov_calls(2), ' / ', count%remov_success(2)
        write(*,*) '3 (Cluster)            : ', count%remov_calls(3), ' / ', count%remov_success(3)
        write(*,*) '4 (Random string)      : ', count%remov_calls(4), ' / ', count%remov_success(4)
        write(*,*) '5 (Bad string)         : ', count%remov_calls(5), ' / ', count%remov_success(5)
        write(*,*) '6 (Good string)        : ', count%remov_calls(6), ' / ', count%remov_success(6)
        write(*,*) '7 (Worst-sum string)   : ', count%remov_calls(7), ' / ', count%remov_success(7)
        write(*,*) '8 (Savings)            : ', count%remov_calls(8), ' / ', count%remov_success(8)
        write(*,*) '9 (SISR random split)  : ', count%remov_calls(9), ' / ', count%remov_success(9)
        write(*,*) '10 (SISR spatial)      : ', count%remov_calls(10), ' / ', count%remov_success(10)
        write(*,*) '11 (Lens)              : ', count%remov_calls(11), ' / ', count%remov_success(11)
        write(*,*) '---------------------------------------------------'
        write(*,*) 'INSERTION OPERATORS (Calls / New best solutions):'
        write(*,*) '1 (Greedy)             : ', count%insert_calls(1), ' / ', count%insert_success(1)
        write(*,*) '2 (Regret-2)           : ', count%insert_calls(2), ' / ', count%insert_success(2)
        write(*,*) '3 (Regret-3)           : ', count%insert_calls(3), ' / ', count%insert_success(3)
        write(*,*) '4 (SISR greedy blinks) : ', count%insert_calls(4), ' / ', count%insert_success(4)
        write(*,*) '==================================================='
        write(*,*) 'Max Iterations without Improvement: ', count%max_no_improv
        write(*,*) 'Total Improvements found: ', count%nimprov
        write(*,*) 'Strings reinserted both ways / reversed chosen: ', rev_tries, ' / ', rev_chosen
        write(*,*) 'SISR insertion orders random / far / close: ', sisr_orders(1), sisr_orders(2), sisr_orders(3)
        write(*,*) '==================================================='
        write(*,*) ''

        tour = bestour
        x = bestsol
        cost = bestcost
        angles = bestangles

        ! Final BCD on the best solution
        param%bcdmaxcic = 150
        call angleBCD(inst,param,count,inst%nballs,tour,x,cost,angles)

    end subroutine alns

    ! *****************************************************************
    ! *****************************************************************

    subroutine removal(inst,param,count,npts,nab,nworst,is_string,&
        tour,x,tangle,removed_cities,op)
        
        implicit none

        ! TYPE ARGUMENTS
        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(inout) :: param
        type(counter_type), intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer, intent(in) :: nworst,op
        integer, intent(out) :: nab,npts
        logical, intent(out) :: is_string
        
        ! ARRAY ARGUMENTS
        integer, intent(inout) :: tour(inst%nballs)
        integer, intent(out) :: removed_cities(nworst)  
        real(kind=8), intent(inout) :: x(2,inst%nballs),tangle(inst%nballs)
        
        nab = sample_triangular(param, param%a, param%b, param%tri_c)
        if ( nab .gt. nworst ) nab = nworst
        if ( nab .ge. inst%nballs ) nab = inst%nballs - 1
        if ( nab .lt. 1 ) nab = 1

        npts = inst%nballs
        is_string = .false. 

        count%remov_calls(op) = count%remov_calls(op) + 1

        ! Select Removal Operator based on 'op'
        if ( op .eq. 1 ) then
            call worst_angles_removal(inst,param,npts,nab,nworst,tour,tangle,removed_cities)
        else if ( op .eq. 2 ) then
            call random_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        else if ( op .eq. 3 ) then
            call cluster_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        else if ( op .eq. 4 ) then
            is_string = .true.
            call random_string_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        else if ( op .eq. 5 ) then
            call bad_string_removal(inst,param,npts,nab,nworst,tour,tangle,removed_cities)
        else if ( op .eq. 6 ) then
            is_string = .true.
            call good_string_removal(inst,param,npts,nab,nworst,tour,tangle,removed_cities)
        else if ( op .eq. 7 ) then
            call worst_sum_string_removal(inst,param,npts,nab,nworst,tour,tangle,removed_cities)
        else if ( op .eq. 8 ) then
            call savings_removal(inst,param,npts,nab,nworst,tour,x,tangle,removed_cities)
        else if ( op .eq. 9 ) then
            call sisr_rand_split_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        else if ( op .eq. 10 ) then
            call sisr_spatial_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        else if ( op .eq. 11 ) then
            ! the two ends of the chosen edge stay in the route
            if ( nab .gt. npts - 2 ) nab = npts - 2
            call lens_removal(inst,param,npts,nab,nworst,tour,x,removed_cities)
        end if

    end subroutine removal

    ! *****************************************************************
    ! *****************************************************************

    subroutine insert(inst,param,count,npts,nab,nworst,cost,is_string,&
    tour,x,tangle,removed_cities,op)
        
        implicit none

        ! TYPE ARGUMENTS
        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(inout) :: param
        type(counter_type), intent(inout) :: count

        ! SCALAR ARGUMENTS
        integer, intent(inout) :: npts
        integer, intent(in) :: nab,nworst,op
        real(kind=8), intent(inout) :: cost
        logical, intent(in) :: is_string
        
        ! ARRAY ARGUMENTS
        integer, intent(inout) :: tour(inst%nballs)
        real(kind=8), intent(inout) :: x(2,inst%nballs),tangle(inst%nballs)
        integer, intent(in) :: removed_cities(nworst)

        ! LOCAL (second attempt, with the reversed string)
        integer :: tourb(inst%nballs),nptsb,revc(nworst)
        real(kind=8) :: xb(2,inst%nballs),tangleb(inst%nballs),costb

        count%insert_calls(op) = count%insert_calls(op) + 1

        ! When the removal produced a string (is_string: operators 4 and 6), it is
        ! reinserted as a single element, in both directions:
        ! the insertion operator is applied twice, with the string in its original and in
        ! its reversed direction, and the better of the two repaired routes is kept. The
        ! string is stored in removed_cities(1:nab) in the order of the route, so that
        ! reversing this segment traverses it in the opposite direction. Reinserted at the
        ! same place, the reversed string is a 2-opt move; elsewhere, an Or-opt move with
        ! reversal.
        if ( is_string .and. nab .ge. 2 ) then
            tourb = tour; xb = x; tangleb = tangle; nptsb = npts; costb = cost

            call apply_string(npts,cost,tour,x,tangle,removed_cities)
            call evalfangle(npts,tour,x,cost,tangle)

            revc = removed_cities
            revc(1:nab) = removed_cities(nab:1:-1)
            call apply_string(nptsb,costb,tourb,xb,tangleb,revc)
            call evalfangle(nptsb,tourb,xb,costb,tangleb)

            rev_tries = rev_tries + 1
            if ( costb .lt. cost ) then
                rev_chosen = rev_chosen + 1
                tour = tourb; x = xb; tangle = tangleb; npts = nptsb; cost = costb
            end if
        else
            call apply_string(npts,cost,tour,x,tangle,removed_cities)

            ! The cost and the angles are recomputed here for the whole route, in O(npts),
            ! because the incremental bookkeeping of the operators may drift, and with
            ! coincident visit points the angles stored in tangle may become inconsistent
            ! (see evalfangle).
            call evalfangle(npts,tour,x,cost,tangle)
        end if

    contains

        subroutine apply_string(np,c,t,xx,tg,rc)
            integer, intent(inout) :: np
            real(kind=8), intent(inout) :: c
            integer, intent(inout) :: t(inst%nballs)
            real(kind=8), intent(inout) :: xx(2,inst%nballs),tg(inst%nballs)
            integer, intent(in) :: rc(nworst)

            if ( op .eq. 1 ) then
                call greedy_insertion(inst,param,count,np,nab,nworst,c,&
                is_string,t,xx,tg,rc)
            else if ( op .eq. 2 ) then
                call regretk_insertion(inst,param,count,2,np,nab,nworst,c,&
                is_string,t,xx,tg,rc)
            else if ( op .eq. 3 ) then
                call regretk_insertion(inst,param,count,3,np,nab,nworst,c,&
                is_string,t,xx,tg,rc)
            else if ( op .eq. 4 ) then
                call sisr_greedy_blink_insertion(inst,param,count,np,nab,nworst,c,&
                is_string,t,xx,tg,rc)
            end if
        end subroutine apply_string

    end subroutine insert

    ! *****************************************************************
    ! *****************************************************************
    integer function sample_triangular(param, a, b, c)
        implicit none
        type(parameter_type), intent(inout) :: param
        integer, intent(in) :: a, b, c
        real(kind=8) :: u, fc, x, da, db, dc

        if ( b .le. a ) then
            sample_triangular = a
            return
        end if

        da = dble(a); db = dble(b); dc = dble(c)
        u = drand(param%seed)
        fc = (dc - da) / (db - da)

        if ( u .lt. fc ) then
            x = da + sqrt( u * (db - da) * (dc - da) )
        else
            x = db - sqrt( (1.0d0 - u) * (db - da) * (db - dc) )
        end if

        sample_triangular = nint(x)
        if ( sample_triangular .lt. a ) sample_triangular = a
        if ( sample_triangular .gt. b ) sample_triangular = b
    end function sample_triangular

    ! *****************************************************************
    integer function select_by_power(param, list_size)
        implicit none
        type(parameter_type), intent(inout) :: param
        integer, intent(in) :: list_size
        real(kind=8) :: y

        if ( list_size .le. 1 ) then
            select_by_power = 1
            return
        end if

        y = drand(param%seed)
        select_by_power = 1 + int( (y ** param%p_worst) * dble(list_size) )
        if ( select_by_power .gt. list_size ) select_by_power = list_size
        if ( select_by_power .lt. 1 ) select_by_power = 1
    end function select_by_power

    ! *****************************************************************
    ! *****************************************************************
    subroutine worst_angles_removal(inst,param,npts,nab,nworst,tour,angles,removed_cities)
        
        implicit none

        ! TYPE ARGUMENTS
        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param

        ! SCALAR ARGUMENTS
        integer,intent(inout) :: npts 
        integer,intent(in) :: nab,nworst

        ! ARRAY ARGUMENTS
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(in) :: angles(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)

        ! LOCAL SCALARS
        integer :: j,k,sel_idx,temp_city,keep_count
        real(kind=8) :: temp_angle
        
        ! LOCAL ARRAYS
        integer :: ranked_cities(npts)
        real(kind=8) :: ranked_angles(npts)
        logical :: is_removed(inst%nballs)

        removed_cities = 0
        is_removed = .false.

        do k=1,npts
            ranked_cities(k) = tour(k)
            ranked_angles(k) = angles(tour(k))
        end do
        do k=2,npts
            temp_angle = ranked_angles(k)
            temp_city  = ranked_cities(k)
            j = k - 1
            do while ( j .ge. 1 )
                if ( ranked_angles(j) .ge. temp_angle ) exit
                ranked_angles(j+1) = ranked_angles(j)
                ranked_cities(j+1) = ranked_cities(j)
                j = j - 1
            end do
            ranked_angles(j+1) = temp_angle
            ranked_cities(j+1) = temp_city
        end do

        do k=1,nab
            sel_idx = (k-1) + select_by_power(param, npts-k+1)
            if ( sel_idx .gt. npts ) sel_idx = npts

            removed_cities(k) = ranked_cities(sel_idx)
            is_removed(ranked_cities(sel_idx)) = .true.

            temp_city = ranked_cities(k)
            ranked_cities(k) = ranked_cities(sel_idx)
            ranked_cities(sel_idx) = temp_city
        end do

        keep_count = 0
        do k=1,npts
            if (.not. is_removed(tour(k))) then
                keep_count = keep_count + 1
                tour(keep_count) = tour(k)
            end if
        end do
        npts = keep_count

    end subroutine worst_angles_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine random_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: k,rand_idx,pool_size,selected_pos,keep_count
        real(kind=8) :: r
        integer :: pool(inst%nballs)
        logical :: is_removed(inst%nballs)

        removed_cities = 0
        is_removed = .false.
        pool_size = npts
        
        pool(1:pool_size) = (/ (k, k=1,pool_size) /)

        do k=1,nab
            r = drand(param%seed)
            rand_idx = 1 + int( r * dble( pool_size ) )
            rand_idx = min( rand_idx, pool_size ) 

            selected_pos = pool(rand_idx)
            
            removed_cities(k) = tour(selected_pos)
            is_removed(removed_cities(k)) = .true.

            pool(rand_idx) = pool(pool_size)
            pool_size = pool_size - 1
        end do

        keep_count = 0
        do k=1,npts
            if (.not. is_removed(tour(k))) then
                keep_count = keep_count + 1
                tour(keep_count) = tour(k)
            end if
        end do

        if (keep_count .lt. inst%nballs) then
            tour(keep_count+1:inst%nballs) = 0
        end if

        npts = keep_count

    end subroutine random_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine cluster_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: k,seed_idx,seed_city,curr_city,keep_count
        real(kind=8) :: r
        logical :: is_removed(inst%nballs)

        is_removed = .false.
        removed_cities = 0

        r = drand(param%seed)
        seed_idx = 1 + int( r * dble( npts ) )
        seed_idx = min( seed_idx, npts ) 
        seed_city = tour(seed_idx)

        do k=1,nab
            curr_city = inst%nn_list(seed_city,k) 
            removed_cities(k) = curr_city
            is_removed(curr_city) = .true.
        end do

        keep_count = 0
        do k=1,npts
            if (.not. is_removed(tour(k))) then
                keep_count = keep_count + 1
                tour(keep_count) = tour(k)
            end if
        end do

        if (keep_count .lt. inst%nballs) then
            tour(keep_count+1:inst%nballs) = 0
        end if

        npts = keep_count

    end subroutine cluster_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine random_string_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: seed_idx,start_idx,end_idx,wrap_end
        real(kind=8) :: r

        removed_cities = 0

        r = drand(param%seed)
        seed_idx = 1 + int(r * dble(npts))
        seed_idx = min(seed_idx,npts) 
        
        start_idx = seed_idx - ( nab / 2 )
        if ( start_idx .lt. 1 ) then
            start_idx = start_idx + npts
        end if

        end_idx = start_idx + nab - 1

        if ( end_idx .le. npts ) then
            removed_cities(1:nab) = tour(start_idx:end_idx)
            if ( end_idx .lt. npts ) then
                tour(start_idx:npts-nab) = tour(end_idx+1:npts)
            end if
        else
            wrap_end = end_idx - npts
            removed_cities(1:npts-start_idx+1) = tour(start_idx:npts)
            removed_cities(npts-start_idx+2:nab) = tour(1:wrap_end)
            tour(1:start_idx-wrap_end-1) = tour(wrap_end+1:start_idx-1)
        end if

        npts = npts - nab

        if ( npts .lt. inst%nballs ) then
            tour(npts+1:inst%nballs) = 0
        end if

    end subroutine random_string_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine bad_string_removal(inst,param,npts,nab,nworst,tour,angles,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(inout) :: angles(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: j,k,seed_idx,start_idx,end_idx,wrap_end,temp_idx
        real(kind=8) :: temp_angle
        integer :: ranked_pos(npts)
        real(kind=8) :: ranked_angles(npts)

        removed_cities = 0

        do k=1,npts
            ranked_pos(k)    = k          
            ranked_angles(k) = angles(tour(k))
        end do
        do k=2,npts
            temp_angle = ranked_angles(k)
            temp_idx   = ranked_pos(k)
            j = k - 1
            do while ( j .ge. 1 )
                if ( ranked_angles(j) .ge. temp_angle ) exit
                ranked_angles(j+1) = ranked_angles(j)
                ranked_pos(j+1)    = ranked_pos(j)
                j = j - 1
            end do
            ranked_angles(j+1) = temp_angle
            ranked_pos(j+1)    = temp_idx
        end do

        seed_idx = ranked_pos( select_by_power(param, npts) )
        
        start_idx = seed_idx - (nab / 2)
        if (start_idx .lt. 1) then
            start_idx = start_idx + npts
        end if
        end_idx = start_idx + nab - 1

        if (end_idx .le. npts) then
            removed_cities(1:nab) = tour(start_idx:end_idx)
            if (end_idx .lt. npts) then
                tour(start_idx : npts - nab) = tour(end_idx + 1 : npts)
            end if
        else
            wrap_end = end_idx - npts
            removed_cities(1 : npts - start_idx + 1) = tour(start_idx:npts)
            removed_cities(npts - start_idx + 2 : nab) = tour(1:wrap_end)
            tour(1 : start_idx - wrap_end - 1) = tour(wrap_end + 1 : start_idx - 1)
        end if

        npts = npts - nab
        if (npts .lt. inst%nballs) then
            tour(npts+1:inst%nballs) = 0
        end if

    end subroutine bad_string_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine good_string_removal(inst,param,npts,nab,nworst,tour,angles,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(inout) :: angles(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: j,k,seed_idx,start_idx,end_idx,wrap_end,temp_idx
        real(kind=8) :: temp_angle
        integer :: ranked_pos(npts)
        real(kind=8) :: ranked_angles(npts)

        removed_cities = 0

        do k=1,npts
            ranked_pos(k)    = k
            ranked_angles(k) = angles(tour(k))
        end do
        do k=2,npts
            temp_angle = ranked_angles(k)
            temp_idx   = ranked_pos(k)
            j = k - 1
            do while ( j .ge. 1 )
                if ( ranked_angles(j) .le. temp_angle ) exit
                ranked_angles(j+1) = ranked_angles(j)
                ranked_pos(j+1)    = ranked_pos(j)
                j = j - 1
            end do
            ranked_angles(j+1) = temp_angle
            ranked_pos(j+1)    = temp_idx
        end do

        seed_idx = ranked_pos( select_by_power(param, npts) )

        start_idx = seed_idx - (nab / 2)
        if (start_idx .lt. 1) then
            start_idx = start_idx + npts
        end if
        
        end_idx = start_idx + nab - 1

        if (end_idx .le. npts) then
            removed_cities(1:nab) = tour(start_idx:end_idx)
            if (end_idx .lt. npts) then
                tour(start_idx : npts - nab) = tour(end_idx + 1 : npts)
            end if
        else
            wrap_end = end_idx - npts
            removed_cities(1 : npts - start_idx + 1) = tour(start_idx:npts)
            removed_cities(npts - start_idx + 2 : nab) = tour(1:wrap_end)
            tour(1 : start_idx - wrap_end - 1) = tour(wrap_end + 1 : start_idx - 1)
        end if

        npts = npts - nab
        if (npts .lt. inst%nballs) then
            tour(npts+1:inst%nballs) = 0
        end if

    end subroutine good_string_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine worst_sum_string_removal(inst,param,npts,nab,nworst,tour,angles,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(in) :: angles(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: j,k,start_idx,end_idx,wrap_end,leaving_idx,entering_idx,temp_idx
        real(kind=8) :: current_sum,temp_sum
        integer :: ranked_starts(npts)
        real(kind=8) :: window_sums(npts)

        removed_cities = 0

        current_sum = 0.0d0
        do k=1,nab
            current_sum = current_sum + angles(tour(k))
        end do
        window_sums(1)   = current_sum
        ranked_starts(1) = 1

        do k=2,npts
            leaving_idx = k - 1
            entering_idx = k + nab - 1
            if (entering_idx .gt. npts) entering_idx = entering_idx - npts
            current_sum = current_sum - angles(tour(leaving_idx)) + angles(tour(entering_idx))
            window_sums(k)   = current_sum
            ranked_starts(k) = k
        end do

        do k=2,npts
            temp_sum = window_sums(k)
            temp_idx = ranked_starts(k)
            j = k - 1
            do while ( j .ge. 1 )
                if ( window_sums(j) .ge. temp_sum ) exit
                window_sums(j+1)   = window_sums(j)
                ranked_starts(j+1) = ranked_starts(j)
                j = j - 1
            end do
            window_sums(j+1)   = temp_sum
            ranked_starts(j+1) = temp_idx
        end do

        start_idx = ranked_starts( select_by_power(param, npts) )

        end_idx = start_idx + nab - 1

        if (end_idx .le. npts) then
            removed_cities(1:nab) = tour(start_idx:end_idx)

            if (end_idx .lt. npts) then
                tour(start_idx : npts - nab) = tour(end_idx + 1 : npts)
            end if
        else
            wrap_end = end_idx - npts
            removed_cities(1 : npts - start_idx + 1) = tour(start_idx:npts)
            removed_cities(npts - start_idx + 2 : nab) = tour(1:wrap_end)
            tour(1 : start_idx - wrap_end - 1) = tour(wrap_end + 1 : start_idx - 1)
        end if

        npts = npts - nab
        if (npts .lt. inst%nballs) then
            tour(npts+1:inst%nballs) = 0
        end if

    end subroutine worst_sum_string_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine savings_removal(inst,param,npts,nab,nworst,tour,x,angles,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(in) :: x(2,inst%nballs)
        real(kind=8),intent(in) :: angles(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: j,k,prev,next,prev2,next2,sel_idx,temp_city,keep_count
        real(kind=8) :: new_angle_prev,new_angle_next,current_savings,temp_sav
        integer :: ranked_cities(npts)
        real(kind=8) :: ranked_savings(npts)
        logical :: is_removed(inst%nballs)

        removed_cities = 0
        is_removed = .false.

        if (npts .le. 3) return

        do k=1,npts
            prev = k - 1
            if (prev .lt. 1) prev = prev + npts
            
            next = k + 1
            if (next .gt. npts) next = next - npts
            
            prev2 = prev - 1
            if (prev2 .lt. 1) prev2 = prev2 + npts
            
            next2 = next + 1
            if (next2 .gt. npts) next2 = next2 - npts

            call turningangle(x(1:2,tour(prev2)), x(1:2,tour(prev)), x(1:2,tour(next)), new_angle_prev)
            call turningangle(x(1:2,tour(prev)), x(1:2,tour(next)), x(1:2,tour(next2)), new_angle_next)

            current_savings = angles(tour(prev)) + angles(tour(k)) + &
            angles(tour(next)) - new_angle_prev - new_angle_next

            ranked_cities(k)  = tour(k)
            ranked_savings(k) = current_savings
        end do

        do k=2,npts
            temp_sav  = ranked_savings(k)
            temp_city = ranked_cities(k)
            j = k - 1
            do while ( j .ge. 1 )
                if ( ranked_savings(j) .ge. temp_sav ) exit
                ranked_savings(j+1) = ranked_savings(j)
                ranked_cities(j+1)  = ranked_cities(j)
                j = j - 1
            end do
            ranked_savings(j+1) = temp_sav
            ranked_cities(j+1)  = temp_city
        end do

        do k=1,nab
            sel_idx = (k-1) + select_by_power(param, npts-k+1)
            if ( sel_idx .gt. npts ) sel_idx = npts

            removed_cities(k) = ranked_cities(sel_idx)
            is_removed(ranked_cities(sel_idx)) = .true.

            temp_city = ranked_cities(k)
            ranked_cities(k) = ranked_cities(sel_idx)
            ranked_cities(sel_idx) = temp_city
        end do

        keep_count = 0
        do k=1,npts
            if (.not. is_removed(tour(k))) then
                keep_count = keep_count + 1
                tour(keep_count) = tour(k)
            end if
        end do

        if (keep_count .lt. inst%nballs) then
            tour(keep_count+1:inst%nballs) = 0
        end if
        
        npts = keep_count

    end subroutine savings_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine sisr_rand_split_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)
        
        integer :: start_idx, l1, l2, m, m_max, k, real_idx, keep_count
        real(kind=8) :: r, p_split
        logical :: is_removed(inst%nballs)

        removed_cities = 0
        is_removed = .false.

        p_split = param%p_split

        r = drand(param%seed)
        start_idx = 1 + int( r * dble( npts ) )
        start_idx = min( start_idx, npts ) 

        r = drand(param%seed)

        if ( r .lt. p_split .and. nab .gt. 1 ) then
            r = drand(param%seed)
            l1 = 1 + int( r * dble( nab - 1 ) )
            l1 = min(l1, nab - 1) 
            l2 = nab - l1

            m_max = 4
            r = drand(param%seed)
            m = 1 + int( r * dble( m_max ) )
            m = min( m, m_max ) 

            if (nab + m .ge. npts) m = max( 1, npts - nab - 1 )

            do k = 0, l1 - 1
                real_idx = mod( start_idx - 1 + k, npts ) + 1
                removed_cities(k + 1) = tour(real_idx)
                is_removed(tour(real_idx)) = .true.
            end do

            do k = 0, l2 - 1
                real_idx = mod( start_idx - 1 + l1 + m + k, npts ) + 1
                removed_cities(l1 + k + 1) = tour(real_idx)
                is_removed(tour(real_idx)) = .true.
            end do

        else
            do k = 0, nab - 1
                real_idx = mod( start_idx - 1 + k, npts ) + 1
                removed_cities(k + 1) = tour(real_idx)
                is_removed(tour(real_idx)) = .true.
            end do
        end if

        keep_count = 0
        do k = 1, npts
            if ( .not. is_removed(tour(k)) ) then
                keep_count = keep_count + 1
                tour(keep_count) = tour(k)
            end if
        end do

        npts = keep_count
        if (npts .lt. inst%nballs) then
            tour(npts+1:inst%nballs) = 0
        end if

    end subroutine sisr_rand_split_removal

    ! *****************************************************************
    ! *****************************************************************
    subroutine sisr_spatial_removal(inst,param,npts,nab,nworst,tour,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        integer,intent(out) :: removed_cities(nworst)

        integer :: idx_A,idx_B,k,real_idx,keep_count,rem_count,city_A,city_k,&
        min_top_dist,top_dist,l1,j,idx_k
        real(kind=8) :: r
        logical :: is_removed(inst%nballs)
        integer :: pos_in_tour(inst%nballs)

        removed_cities = 0
        is_removed = .false.
        rem_count = 0
        pos_in_tour = 0 

        do k = 1, npts
            pos_in_tour(tour(k)) = k
        end do

        r = drand(param%seed)
        idx_A = 1 + int(r * dble(npts))
        idx_A = min(idx_A, npts) 
        city_A = tour(idx_A)

        min_top_dist = max( int( npts * 0.1d0 ), 2 )

        idx_B = -1
        do j = 2, inst%nballs
            city_k = inst%nn_list(city_A, j)
            idx_k = pos_in_tour(city_k)
            
            if (idx_k .eq. 0) cycle
            
            top_dist = min(abs(idx_k - idx_A), npts - abs(idx_k - idx_A))
            
            if (top_dist .ge. min_top_dist) then
                idx_B = idx_k
                exit 
            end if
        end do

        if (idx_B .eq. -1) then
            idx_B = mod(idx_A + npts/2 - 1, npts) + 1
        end if

        l1 = max(1, nab / 2)

        k = 0
        do while (rem_count .lt. l1 .and. k .lt. npts)
            real_idx = mod(idx_A - 1 + k, npts) + 1
            city_k = tour(real_idx)
            
            if (.not. is_removed(city_k)) then
                rem_count = rem_count + 1
                removed_cities(rem_count) = city_k
                is_removed(city_k) = .true.
            end if
            k = k + 1
        end do

        k = 0
        do while (rem_count .lt. nab .and. k .lt. npts)
            real_idx = mod(idx_B - 1 + k, npts) + 1
            city_k = tour(real_idx)
            
            if (.not. is_removed(city_k)) then
                rem_count = rem_count + 1
                removed_cities(rem_count) = city_k
                is_removed(city_k) = .true.
            end if
            k = k + 1
        end do

        keep_count = 0
        do k = 1, npts
            if (.not. is_removed(tour(k))) then
                keep_count = keep_count + 1
                tour(keep_count) = tour(k)
            end if
        end do

        npts = keep_count
        if (npts .lt. inst%nballs) then
            tour(npts+1:inst%nballs) = 0
        end if

    end subroutine sisr_spatial_removal

    ! *****************************************************************
    ! *****************************************************************
    ! lens_removal -- lens removal (inspired by the LENS procedure of Stanek et al.
    ! 2019, also used by Cavagnini et al. 2024). Points close to a straight leg of the
    ! route fit in it with a small turn: inserted between u and v, a point z turns by
    ! turn(x_u, z, x_v), and the points with a turn of at most theta form the "lens" of
    ! thickness theta of the edge u -> v. The operator
    !   1. sorts the edges of the current route by decreasing length and draws one by
    !      the power rule (select_by_power), favoring the long legs;
    !   2. evaluates each disk w outside the edge by the lower bound on the turn it
    !      would make between u and v, as in eq. (angle-bound): turn(x_u, c_w, x_v) minus
    !      arcsin(r_w/|c_w - x_u|) and arcsin(r_w/|c_w - x_v|) (pi when the argument
    !      exceeds 1, and at least zero);
    !   3. removes the nab disks with the smallest turn, those that best fit the edge,
    !      that is, those of the thinnest lens that contains nab disks. There is no
    !      thickness parameter.
    ! Cost O(m log m) per call.
    subroutine lens_removal(inst,param,npts,nab,nworst,tour,x,removed_cities)

        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(in) :: x(2,inst%nballs)
        integer,intent(out) :: removed_cities(nworst)

        integer :: k,e,u,v,w,t,kmin,nr,ord(npts)
        real(kind=8) :: comp(npts),score(npts),p(2),q(2),a,dp,dq,bp,bq,smin
        logical :: leaves(npts)
        real(kind=8), parameter :: pi = acos(-1.0d0)

        removed_cities = 0

        ! 1. edges k -> k+1 (positions), sorted by decreasing length
        do k = 1, npts
            comp(k) = norm2( x(1:2,tour(nextind(k,1,npts))) - x(1:2,tour(k)) )
            ord(k) = k
        end do
        call sort_decreasing(npts,comp,ord)
        e = ord( select_by_power(param, npts) )
        u = tour(e)
        v = tour(nextind(e,1,npts))
        p = x(1:2,u)
        q = x(1:2,v)

        ! 2. lower bound on the turn of each disk inserted between u and v
        do k = 1, npts
            w = tour(k)
            if ( w .eq. u .or. w .eq. v ) then
                score(k) = huge(1.0d0)
                cycle
            end if
            call turningangle(p,inst%centercoor(1:2,w),q,a)
            dp = norm2( inst%centercoor(1:2,w) - p )
            dq = norm2( inst%centercoor(1:2,w) - q )
            bp = pi
            bq = pi
            if ( dp .gt. inst%radii(w) ) bp = asin( inst%radii(w) / dp )
            if ( dq .gt. inst%radii(w) ) bq = asin( inst%radii(w) / dq )
            score(k) = max( 0.0d0, a - bp - bq )
        end do

        ! 3. removes the nab disks with the smallest turn (partial selection, O(nab m))
        leaves = .false.
        do t = 1, nab
            smin = huge(1.0d0)
            kmin = 0
            do k = 1, npts
                if ( leaves(k) ) cycle
                if ( score(k) .lt. smin ) then
                    smin = score(k)
                    kmin = k
                end if
            end do
            leaves(kmin) = .true.
            removed_cities(t) = tour(kmin)
        end do

        nr = 0
        do k = 1, npts
            if ( leaves(k) ) cycle
            nr = nr + 1
            tour(nr) = tour(k)
        end do
        npts = nr
        if ( npts .lt. inst%nballs ) tour(npts+1:inst%nballs) = 0

    end subroutine lens_removal

    ! *****************************************************************
    ! sort_decreasing -- heapsort of the indices ord(1:n) by decreasing key(ord)
    subroutine sort_decreasing(n,chave,ord)

        implicit none

        integer,intent(in) :: n
        real(kind=8),intent(in) :: chave(n)
        integer,intent(inout) :: ord(n)

        integer :: i,iend,tmp

        do i = n / 2, 1, -1
            call sift_down(i,n)
        end do
        do iend = n, 2, -1
            tmp = ord(1); ord(1) = ord(iend); ord(iend) = tmp
            call sift_down(1,iend-1)
        end do

    contains

        ! min-heap on the key, so that the result is in decreasing order
        subroutine sift_down(root_node,last)
            integer,intent(in) :: root_node,last
            integer :: r,f,t
            r = root_node
            do
                f = 2 * r
                if ( f .gt. last ) exit
                if ( f + 1 .le. last ) then
                    if ( chave(ord(f+1)) .lt. chave(ord(f)) ) f = f + 1
                end if
                if ( chave(ord(r)) .le. chave(ord(f)) ) exit
                t = ord(r); ord(r) = ord(f); ord(f) = t
                r = f
            end do
        end subroutine sift_down

    end subroutine sort_decreasing

    ! *****************************************************************
    ! *****************************************************************
    subroutine greedy_insertion(inst,param,count,npts,nab,nworst,cost,is_string,tour,&
    x,tangle,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        type(counter_type),intent(inout) :: count
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        real(kind=8),intent(inout) :: cost
        logical,intent(in) :: is_string
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(inout) :: x(2,inst%nballs),tangle(inst%nballs)
        integer,intent(in) :: removed_cities(nworst)

        logical :: pruned
        integer :: i,j,k,ell,kprev,knext,jbest,pospr,pospp,posne,posnn,numins,&
        bsize,ncbcd,natoc
        real(kind=8) :: diff,costbest,costtrial,strincost,cost_lower_bound
        real(kind=8) :: costcmp,costbestcmp
        integer :: tourtrial(inst%nballs),pos(5),pbest(4),idabcd(8),cityabcdbest(8)
        real(kind=8) :: xtmp(2,inst%nballs),xbest(2,4),zpp(2),zpr(2),&
        zne(2),znn(2),tangletmp(8),tanglebest(8)

        if ( is_string ) then
            numins = 1      
            bsize = nab     
            ncbcd = 4       
            natoc = 8       
        else
            numins = nab
            bsize = 1
            ncbcd = 3
            natoc = 5
        end if

        cost = 0.0d0
        do k = 1, npts
            kprev = prevind(k,1,npts)
            knext = nextind(k,1,npts)
            call turningangle(x(1:2,tour(kprev)),x(1:2,tour(k)),&
            x(1:2,tour(knext)),tangle(tour(k)))
            cost = cost + tangle(tour(k))
        end do

        strincost = 0.0d0
        if ( is_string .and. bsize .gt. 2 ) then
            do k = 2, bsize - 1
                call turningangle(x(1:2,removed_cities(k-1)),x(1:2,removed_cities(k)),&
                x(1:2,removed_cities(k+1)),tangle(removed_cities(k)))
                strincost = strincost + tangle(removed_cities(k))
            end do
        end if


        do i = 1, numins
            costbest = huge(1.0d0)
            costbestcmp = huge(1.0d0)
            jbest = 0
            do j = 0, npts - 1

                pos(1) = npts - j
                pos(2) = npts + 1 - j
                if ( is_string ) then
                    pos(3) = npts + bsize - j
                    pos(4) = nextind(pos(3),1,npts+bsize)
                    pos(5) = nextind(pos(4),1,npts+bsize)
                else
                    pos(3) = nextind(pos(2),1,npts+bsize)   
                    pos(4) = nextind(pos(3),1,npts+bsize)   
                end if

                if ( j .eq. 0 ) then
                    tourtrial(1:npts) = tour(1:npts)
                    if ( is_string ) then
                        tourtrial(npts+1:npts+bsize) = removed_cities(1:bsize)
                    else
                        tourtrial(npts+bsize) = removed_cities(i)
                    end if
                    costtrial = cost - tangle(tour(npts)) - tangle(tour(1))
                else
                    tourtrial(npts+bsize-j+1) = tourtrial(npts+1-j)
                    if ( is_string ) then
                        tourtrial(npts+1-j:npts+bsize-j) = removed_cities(1:bsize)
                        costtrial = cost - tangle(tourtrial(pos(1))) - tangle(tourtrial(pos(4)))
                    else
                        tourtrial(npts+1-j) = removed_cities(i)
                        costtrial = cost - tangle(tourtrial(pos(1))) - tangle(tourtrial(pos(3)))
                    end if
                end if

                if ( is_string ) then
                    costtrial = costtrial + strincost
                    idabcd(1) = prevind(pos(1),1,npts+bsize)
                    idabcd(2) = pos(1)
                    idabcd(3) = pos(2)
                    idabcd(4) = nextind(pos(2),1,npts+bsize)    
                    idabcd(5) = prevind(pos(3),1,npts+bsize)    
                    idabcd(6) = pos(3)
                    idabcd(7) = pos(4)
                    idabcd(8) = nextind(pos(4),1,npts+bsize)
                    
                    costtrial = costtrial - tangle(tourtrial(idabcd(1))) - tangle(tourtrial(idabcd(8)))
                    if ( bsize .gt. 2 ) then
                        costtrial = costtrial - tangle(tourtrial(idabcd(4)))
                        if ( bsize .gt. 3 ) then
                            costtrial = costtrial - tangle(tourtrial(idabcd(5)))
                        end if
                    end if
                else
                    idabcd(1) = prevind(pos(1),1,npts+bsize)
                    idabcd(2) = pos(1)
                    idabcd(3) = pos(2)
                    idabcd(4) = pos(3)
                    idabcd(5) = nextind(pos(3),1,npts+bsize)
                    costtrial = costtrial - tangle(tourtrial(idabcd(1))) - tangle(tourtrial(idabcd(5)))
                end if

                cost_lower_bound = costtrial
                
                if ( is_string ) then
                    cost_lower_bound = cost_lower_bound &
                        + inst%min_angle_matrix(tourtrial(prevind(idabcd(1),1,npts+bsize)),&
                        tourtrial(idabcd(1)),tourtrial(idabcd(2))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(1)),tourtrial(idabcd(2)),&
                        tourtrial(idabcd(3))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(2)),tourtrial(idabcd(3)),&
                        tourtrial(nextind(idabcd(3),1,npts+bsize)))
                    
                    cost_lower_bound = cost_lower_bound &
                        + inst%min_angle_matrix(tourtrial(prevind(idabcd(6),1,npts+bsize)),&
                        tourtrial(idabcd(6)), tourtrial(idabcd(7))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(6)),tourtrial(idabcd(7)),&
                        tourtrial(idabcd(8))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(7)),tourtrial(idabcd(8)),&
                        tourtrial(nextind(idabcd(8),1,npts+bsize)))

                    if ( bsize .gt. 2 ) then
                        cost_lower_bound = cost_lower_bound &
                            + inst%min_angle_matrix(tourtrial(prevind(idabcd(4),1,npts+bsize)),&
                            tourtrial(idabcd(4)),tourtrial(nextind(idabcd(4),1,npts+bsize)))
                        if ( bsize .gt. 3 ) then
                            cost_lower_bound = cost_lower_bound &
                                + inst%min_angle_matrix(tourtrial(prevind(idabcd(5),1,npts+bsize)),&
                                tourtrial(idabcd(5)), tourtrial(nextind(idabcd(5),1,npts+bsize)))
                        end if
                    end if
                else
                    cost_lower_bound = cost_lower_bound &
                        + inst%min_angle_matrix(tourtrial(prevind(idabcd(1),1,npts+bsize)),&
                        tourtrial(idabcd(1)), tourtrial(idabcd(2))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(1)),tourtrial(idabcd(2)),&
                        tourtrial(idabcd(3))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(2)),tourtrial(idabcd(3)),&
                        tourtrial(idabcd(4))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(3)),tourtrial(idabcd(4)),&
                        tourtrial(idabcd(5))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(4)),tourtrial(idabcd(5)),&
                        tourtrial(nextind(idabcd(5),1,npts+bsize)))
                end if

                if ( cost_lower_bound .ge. costbestcmp ) then
                    pruned = .true.
                else
                    pruned = .false.
                end if

                if ( .not. pruned ) then

                    xtmp(1:2,tourtrial(pos(1:ncbcd))) = x(1:2,tourtrial(pos(1:ncbcd)))
                    do k = 1, natoc
                        tangletmp(k) = tangle(tourtrial(idabcd(k)))
                    end do

                    do ell = 1, ncbcd
                        pospr = prevind(pos(ell),1,npts+bsize)
                        pospp = prevind(pospr   ,1,npts+bsize)
                        posne = nextind(pos(ell),1,npts+bsize)
                        posnn = nextind(posne   ,1,npts+bsize)
                    
                        zpp(1:2) = x(1:2,tourtrial(pospp))
                        zpr(1:2) = x(1:2,tourtrial(pospr))
                        zne(1:2) = x(1:2,tourtrial(posne))
                        znn(1:2) = x(1:2,tourtrial(posnn))
                
                        call discmin(inst,param,count,tourtrial(pos(ell)),zpp,zpr,zne,znn,&
                        x(1:2,tourtrial(pos(ell))),diff,tangle(tourtrial(pospr)),&
                        tangle(tourtrial(pos(ell))),tangle(tourtrial(posne))) 
                    end do

                    if ( is_string ) then
                        costtrial = costtrial + tangle(tourtrial(idabcd(1))) + tangle(tourtrial(idabcd(2))) &
                                    + tangle(tourtrial(idabcd(3))) + tangle(tourtrial(idabcd(6))) &
                                    + tangle(tourtrial(idabcd(7))) + tangle(tourtrial(idabcd(8)))
                        
                        if ( bsize .gt. 2 ) then
                            costtrial = costtrial + tangle(tourtrial(idabcd(4)))
                            if ( bsize .gt. 3 ) then
                                costtrial = costtrial + tangle(tourtrial(idabcd(5)))
                            end if
                        end if
                    else
                        costtrial = costtrial + tangle(tourtrial(idabcd(1))) + tangle(tourtrial(idabcd(2))) &
                                    + tangle(tourtrial(idabcd(3))) + tangle(tourtrial(idabcd(4))) &
                                    + tangle(tourtrial(idabcd(5)))
                    end if

                    ! the positions are compared by their cost
                    costcmp = costtrial

                    if ( j .eq. 0 .or. costcmp .lt. costbestcmp ) then
                        xbest(1:2,1:ncbcd) = x(1:2,tourtrial(pos(1:ncbcd)))
                        pbest(1:ncbcd) = pos(1:ncbcd)
                        do k = 1, natoc
                            cityabcdbest(k) = tourtrial(idabcd(k))
                            tanglebest(k) = tangle(tourtrial(idabcd(k)))
                        end do
                        costbest = costtrial
                        costbestcmp = costcmp
                        jbest = j
                    end if

                    x(1:2,tourtrial(pos(1:ncbcd))) = xtmp(1:2,tourtrial(pos(1:ncbcd)))
                    do k = 1, natoc
                        tangle(tourtrial(idabcd(k))) = tangletmp(k)
                    end do
                end if
            end do
            
            if ( jbest .gt. 0 ) then
                tour(npts+1-jbest+bsize:npts+bsize) = tour(npts+1-jbest:npts)
            end if
    
            if ( is_string ) then
                tour(npts+1-jbest:npts+bsize-jbest) = removed_cities(1:bsize)
            else
                tour(npts+1-jbest) = removed_cities(i)
            end if

            x(1:2,tour(pbest(1:ncbcd))) = xbest(1:2,1:ncbcd)

            do k = 1, natoc
                tangle(cityabcdbest(k)) = tanglebest(k)
            end do

            cost = costbest
            npts = npts + bsize
        end do

    end subroutine greedy_insertion

    subroutine regretk_insertion(inst,param,count,regk,npts,nab,nworst,cost,is_string,&
    tour,x,tangle,removed_cities)
        
        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        type(counter_type),intent(inout) :: count
        integer,intent(in) :: regk,nab,nworst
        integer,intent(inout) :: npts
        real(kind=8),intent(inout) :: cost
        logical,intent(in) :: is_string
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(inout) :: x(2,inst%nballs),tangle(inst%nballs)
        integer,intent(in) :: removed_cities(nworst)

        integer :: step,i,j,k,ell,kprev,knext,pospr,pospp,posne,posnn,&
        best_city_idx,local_jbest,global_jbest,city_in
        real(kind=8) :: diff,costtrial,cost1,cost2,cost3,regret,max_regret,&
        cost_lower_bound,threshold,ang_pr,ang_in,ang_ne

        logical :: inserted(nworst)
        integer :: tourtrial(inst%nballs),pos(3),idabcd(5)
        real(kind=8) :: local_tanglebest(3),global_tanglebest(3)
        real(kind=8) :: zpp(2),zpr(2),zne(2),znn(2)

        ! Indivisible strings do not compete through the regret; they go to the greedy insertion
        if ( is_string ) then
            call greedy_insertion(inst,param,count,npts,nab,nworst,cost,is_string,&
            tour,x,tangle,removed_cities)
            return
        end if

        inserted = .false.


        cost = 0.0d0
        do k = 1, npts
            kprev = prevind(k,1,npts)
            knext = nextind(k,1,npts)
            call turningangle(x(1:2,tour(kprev)),x(1:2,tour(k)),&
            x(1:2,tour(knext)),tangle(tour(k)))
            cost = cost + tangle(tour(k))
        end do

        ! Main loop: inserts exactly one disk per step
        do step = 1, nab
            
            max_regret = -huge(1.0d0)
            best_city_idx = 0

            ! Evaluates all disks not yet inserted
            do i = 1, nab
                if ( inserted(i) ) cycle
                
                city_in = removed_cities(i)
                cost1 = huge(1.0d0) / 10.0d0
                cost2 = huge(1.0d0) / 10.0d0
                cost3 = huge(1.0d0) / 10.0d0

                do j = 0, npts - 1
                    pos(1) = npts - j
                    pos(2) = npts + 1 - j
                    pos(3) = nextind(pos(2),1,npts+1)

                    if ( j .eq. 0 ) then
                        tourtrial(1:npts) = tour(1:npts)
                        tourtrial(npts+1) = city_in
                        costtrial = cost - tangle(tour(npts)) - tangle(tour(1))
                    else
                        tourtrial(npts+2-j) = tourtrial(npts+1-j)
                        tourtrial(npts+1-j) = city_in
                        costtrial = cost - tangle(tourtrial(pos(1))) - tangle(tourtrial(pos(3)))
                    end if

                    idabcd(1) = prevind(pos(1),1,npts+1)
                    idabcd(2) = pos(1)
                    idabcd(3) = pos(2)
                    idabcd(4) = pos(3)
                    idabcd(5) = nextind(pos(3),1,npts+1)

                    if ( regk .eq. 2 ) then
                        threshold = cost2
                    else
                        threshold = cost3
                    end if

                    ! O(1) pruning: evaluates the lower bound only for the 3 angles formed with the disk
                    cost_lower_bound = costtrial &
                        + inst%min_angle_matrix(tourtrial(idabcd(1)), tourtrial(idabcd(2)), tourtrial(idabcd(3))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(2)), tourtrial(idabcd(3)), tourtrial(idabcd(4))) &
                        + inst%min_angle_matrix(tourtrial(idabcd(3)), tourtrial(idabcd(4)), tourtrial(idabcd(5)))

                    if ( cost_lower_bound .ge. threshold ) then
                        cycle
                    end if

                    ! FAST EVALUATION: uses the center of the disk only. No BCD routine is called.
                    call turningangle(x(1:2, tourtrial(idabcd(1))), x(1:2, tourtrial(idabcd(2))), &
                    inst%centercoor(1:2, city_in), ang_pr)
                    call turningangle(x(1:2, tourtrial(idabcd(2))), inst%centercoor(1:2, city_in), &
                    x(1:2, tourtrial(idabcd(4))), ang_in)
                    call turningangle(inst%centercoor(1:2, city_in), x(1:2, tourtrial(idabcd(4))), &
                    x(1:2, tourtrial(idabcd(5))), ang_ne)

                    costtrial = costtrial + ang_pr + ang_in + ang_ne

                    if ( costtrial .lt. cost1 ) then
                        cost3 = cost2
                        cost2 = cost1
                        cost1 = costtrial
                        local_jbest = j
                        
                        ! Stores the 3 evaluated angles
                        local_tanglebest(1) = ang_pr
                        local_tanglebest(2) = ang_in
                        local_tanglebest(3) = ang_ne
                        
                    else if ( costtrial .lt. cost2 ) then
                        cost3 = cost2
                        cost2 = costtrial
                    else if ( costtrial .lt. cost3 ) then
                        if ( regk .eq. 3 ) cost3 = costtrial
                    end if
                end do

                if ( regk .eq. 2 ) then
                    regret = cost2 - cost1
                else
                    regret = ( cost2 - cost1 ) + ( cost3 - cost1 )
                end if
                
                if ( regret .gt. max_regret ) then
                    max_regret = regret
                    best_city_idx = i
                    global_jbest = local_jbest
                    global_tanglebest = local_tanglebest
                end if
            end do

            ! =========================================================
            ! INSERTION OF THE CHOSEN DISK
            ! =========================================================
            inserted(best_city_idx) = .true.
            city_in = removed_cities(best_city_idx)
            
            if ( global_jbest .gt. 0 ) then
                tour(npts+2-global_jbest:npts+1) = tour(npts+1-global_jbest:npts)
            end if
            tour(npts+1-global_jbest) = city_in

            x(1:2, city_in) = inst%centercoor(1:2, city_in)

            pos(1) = npts + 1 - global_jbest - 1
            if (pos(1) .lt. 1) pos(1) = pos(1) + npts + 1
            pos(2) = npts + 1 - global_jbest
            pos(3) = nextind(pos(2), 1, npts+1)

            idabcd(1) = prevind(pos(1), 1, npts+1)
            idabcd(2) = pos(1)
            idabcd(3) = pos(2)
            idabcd(4) = pos(3)
            idabcd(5) = nextind(pos(3), 1, npts+1)

            ! Breaks the edge and subtracts the two old angles
            cost = cost - tangle(tour(idabcd(2))) - tangle(tour(idabcd(4)))
            
            ! Adds the angles stored in the evaluation step
            tangle(tour(idabcd(2))) = global_tanglebest(1)
            tangle(tour(idabcd(3))) = global_tanglebest(2)
            tangle(tour(idabcd(4))) = global_tanglebest(3)

            cost = cost + global_tanglebest(1) + global_tanglebest(2) + global_tanglebest(3)

            ! Finally, block updates of the inserted point and its neighbors
            do ell = 1, 3
                pospr = prevind(pos(ell), 1, npts+1)
                pospp = prevind(pospr   , 1, npts+1)
                posne = nextind(pos(ell), 1, npts+1)
                posnn = nextind(posne   , 1, npts+1)
            
                zpp(1:2) = x(1:2, tour(pospp))
                zpr(1:2) = x(1:2, tour(pospr))
                zne(1:2) = x(1:2, tour(posne))
                znn(1:2) = x(1:2, tour(posnn))
        
                call discmin(inst,param,count,tour(pos(ell)),zpp,zpr,zne,znn,&
                x(1:2,tour(pos(ell))),diff,tangle(tour(pospr)),&
                tangle(tour(pos(ell))),tangle(tour(posne)))

                cost = cost + diff
            end do

            npts = npts + 1
        end do

    end subroutine regretk_insertion
    
    ! *****************************************************************
    ! *****************************************************************
    subroutine sisr_greedy_blink_insertion(inst,param,count,npts,nab,nworst,cost,is_string,&
    tour,x,tangle,removed_in)

        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        type(counter_type),intent(inout) :: count
        integer,intent(inout) :: npts
        integer,intent(in) :: nab,nworst
        real(kind=8),intent(inout) :: cost
        logical,intent(in) :: is_string
        integer,intent(inout) :: tour(inst%nballs)
        real(kind=8),intent(inout) :: x(2,inst%nballs),tangle(inst%nballs)
        integer,intent(in) :: removed_in(nworst)

        integer :: i,j,k,ell,kprev,knext,pospr,pospp,posne,posnn,jbest,numins,bsize,ncbcd,natoc
        real(kind=8) :: diff,costbest,costtrial,strincost,p_blink,r,cost_lower_bound
        logical :: blinked,evaluated_at_least_one,pruned
        integer :: tourtrial(inst%nballs),pos(4),idabcd(8),pbest(4),cityabcdbest(8)
        real(kind=8) :: xtmp(2,inst%nballs),xbest(2,4),zpp(2),zpr(2),&
        zne(2),znn(2),tangletmp(8),tanglebest(8)
        integer :: removed_cities(nworst)

        p_blink = param%p_blink

        ! Reinsertion order (Christiaens and Vanden Berghe 2020): the removed disks are
        ! sorted before they are reinserted one at a time. Without demands, the orders are
        ! random (weight 4), farthest first (weight 2) and closest first (weight 1) from
        ! the centroid of the centers, which plays the role of the depot. Strings are
        ! reinserted as single elements and are not reordered.
        removed_cities = removed_in
        if ( .not. is_string .and. nab .ge. 2 ) then
            call sisr_reorder(inst,param,nab,removed_cities(1:nab))
        end if

        if ( is_string ) then
            numins = 1      
            bsize = nab     
            ncbcd = 4       
            natoc = 8       
        else
            numins = nab
            bsize = 1
            ncbcd = 3
            natoc = 5
        end if

        cost = 0.0d0
        do k = 1, npts
            kprev = prevind(k,1,npts)
            knext = nextind(k,1,npts)
            call turningangle(x(1:2,tour(kprev)),x(1:2,tour(k)),&
            x(1:2,tour(knext)),tangle(tour(k)))
            cost = cost + tangle(tour(k))
        end do

        strincost = 0.0d0
        if ( is_string .and. bsize .gt. 2 ) then
            do k = 2, bsize - 1
                call turningangle(x(1:2,removed_cities(k-1)),x(1:2,removed_cities(k)),&
                x(1:2,removed_cities(k+1)),tangle(removed_cities(k)))
                strincost = strincost + tangle(removed_cities(k))
            end do
        end if

        do i = 1, numins
            
            costbest = huge(1.0d0)
            jbest = 0
            evaluated_at_least_one = .false.

            do j = 0, npts - 1
                
                r = drand(param%seed)
                
                if ( r .lt. p_blink .and. (j .lt. npts - 1 .or. evaluated_at_least_one) ) then
                    blinked = .true.
                else
                    blinked = .false.
                    evaluated_at_least_one = .true.
                end if

                pos(1) = npts - j           
                pos(2) = npts + 1 - j       
                if ( is_string ) then
                    pos(3) = npts + bsize - j               
                    pos(4) = nextind(pos(3),1,npts+bsize)   
                else
                    pos(3) = nextind(pos(2),1,npts+bsize)   
                    pos(4) = nextind(pos(3),1,npts+bsize)   
                end if

                if ( j .eq. 0 ) then
                    tourtrial(1:npts) = tour(1:npts)
                    if ( is_string ) then
                        tourtrial(npts+1:npts+bsize) = removed_cities(1:bsize)
                    else
                        tourtrial(npts+bsize) = removed_cities(i)
                    end if
                    costtrial = cost - tangle(tour(npts)) - tangle(tour(1))
                else
                    tourtrial(npts+bsize-j+1) = tourtrial(npts+1-j)
                    if ( is_string ) then
                        tourtrial(npts+1-j:npts+bsize-j) = removed_cities(1:bsize)
                        costtrial = cost - tangle(tourtrial(pos(1))) - tangle(tourtrial(pos(4)))
                    else
                        tourtrial(npts+1-j) = removed_cities(i)
                        costtrial = cost - tangle(tourtrial(pos(1))) - tangle(tourtrial(pos(3)))
                    end if
                end if

                if ( .not. blinked ) then
                    if ( is_string ) then
                        costtrial = costtrial + strincost
                        idabcd(1) = prevind(pos(1),1,npts+bsize)
                        idabcd(2) = pos(1)
                        idabcd(3) = pos(2)
                        idabcd(4) = nextind(pos(2),1,npts+bsize)  
                        idabcd(5) = prevind(pos(3),1,npts+bsize)  
                        idabcd(6) = pos(3)
                        idabcd(7) = pos(4)
                        idabcd(8) = nextind(pos(4),1,npts+bsize)
                        
                        costtrial = costtrial - tangle(tourtrial(idabcd(1))) - tangle(tourtrial(idabcd(8)))
                        if ( bsize .gt. 2 ) then
                            costtrial = costtrial - tangle(tourtrial(idabcd(4)))
                            if ( bsize .gt. 3 ) costtrial = costtrial - tangle(tourtrial(idabcd(5)))
                        end if
                    else
                        idabcd(1) = prevind(pos(1),1,npts+bsize)
                        idabcd(2) = pos(1)
                        idabcd(3) = pos(2)
                        idabcd(4) = pos(3)
                        idabcd(5) = nextind(pos(3),1,npts+bsize)
                        costtrial = costtrial - tangle(tourtrial(idabcd(1))) - tangle(tourtrial(idabcd(5)))
                    end if

                    cost_lower_bound = costtrial

                    if ( is_string ) then
                        cost_lower_bound = cost_lower_bound &
                            + inst%min_angle_matrix(tourtrial(prevind(idabcd(1),1,npts+bsize)),&
                            tourtrial(idabcd(1)), tourtrial(idabcd(2))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(1)), tourtrial(idabcd(2)),&
                            tourtrial(idabcd(3))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(2)), tourtrial(idabcd(3)),&
                            tourtrial(nextind(idabcd(3),1,npts+bsize)))
                        
                        cost_lower_bound = cost_lower_bound &
                            + inst%min_angle_matrix(tourtrial(prevind(idabcd(6),1,npts+bsize)),&
                            tourtrial(idabcd(6)), tourtrial(idabcd(7))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(6)), tourtrial(idabcd(7)),&
                            tourtrial(idabcd(8))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(7)), tourtrial(idabcd(8)),&
                            tourtrial(nextind(idabcd(8),1,npts+bsize)))

                        if ( bsize .gt. 2 ) then
                            cost_lower_bound = cost_lower_bound &
                                + inst%min_angle_matrix(tourtrial(prevind(idabcd(4),1,npts+bsize)),&
                                tourtrial(idabcd(4)), tourtrial(nextind(idabcd(4),1,npts+bsize)))
                            if ( bsize .gt. 3 ) then
                                cost_lower_bound = cost_lower_bound &
                                    + inst%min_angle_matrix(tourtrial(prevind(idabcd(5),1,npts+bsize)),&
                                    tourtrial(idabcd(5)), tourtrial(nextind(idabcd(5),1,npts+bsize)))
                            end if
                        end if
                    else
                        cost_lower_bound = cost_lower_bound &
                            + inst%min_angle_matrix(tourtrial(prevind(idabcd(1),1,npts+bsize)),&
                            tourtrial(idabcd(1)), tourtrial(idabcd(2))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(1)),tourtrial(idabcd(2)),&
                            tourtrial(idabcd(3))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(2)),tourtrial(idabcd(3)),&
                            tourtrial(idabcd(4))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(3)),tourtrial(idabcd(4)),&
                            tourtrial(idabcd(5))) &
                            + inst%min_angle_matrix(tourtrial(idabcd(4)),tourtrial(idabcd(5)),&
                            tourtrial(nextind(idabcd(5),1,npts+bsize)))
                    end if

                    if ( cost_lower_bound .ge. costbest ) then
                        pruned = .true.
                    else
                        pruned = .false.
                    end if

                    if ( .not. pruned ) then

                        xtmp(1:2,tourtrial(pos(1:ncbcd))) = x(1:2,tourtrial(pos(1:ncbcd)))
                        do k = 1, natoc
                            tangletmp(k) = tangle(tourtrial(idabcd(k)))
                        end do

                        do ell = 1, ncbcd
                            pospr = prevind(pos(ell),1,npts+bsize)
                            pospp = prevind(pospr   ,1,npts+bsize)
                            posne = nextind(pos(ell),1,npts+bsize)
                            posnn = nextind(posne   ,1,npts+bsize)

                            zpp(1:2) = x(1:2,tourtrial(pospp))
                            zpr(1:2) = x(1:2,tourtrial(pospr))
                            zne(1:2) = x(1:2,tourtrial(posne))
                            znn(1:2) = x(1:2,tourtrial(posnn))

                            call discmin(inst,param,count,tourtrial(pos(ell)),zpp,zpr,zne,znn,&
                            x(1:2,tourtrial(pos(ell))),diff,tangle(tourtrial(pospr)),&
                            tangle(tourtrial(pos(ell))),tangle(tourtrial(posne))) 
                        end do

                        if ( is_string ) then
                            costtrial = costtrial + tangle(tourtrial(idabcd(1))) + tangle(tourtrial(idabcd(2))) &
                                        + tangle(tourtrial(idabcd(3))) + tangle(tourtrial(idabcd(6))) &
                                        + tangle(tourtrial(idabcd(7))) + tangle(tourtrial(idabcd(8)))
                            if ( bsize .gt. 2 ) then
                                costtrial = costtrial + tangle(tourtrial(idabcd(4)))
                                if ( bsize .gt. 3 ) costtrial = costtrial + tangle(tourtrial(idabcd(5)))
                            end if
                        else
                            costtrial = costtrial + tangle(tourtrial(idabcd(1))) + tangle(tourtrial(idabcd(2))) &
                                        + tangle(tourtrial(idabcd(3))) + tangle(tourtrial(idabcd(4))) &
                                        + tangle(tourtrial(idabcd(5)))
                        end if

                        if ( costtrial .lt. costbest ) then
                            xbest(1:2,1:ncbcd) = x(1:2,tourtrial(pos(1:ncbcd)))
                            pbest(1:ncbcd) = pos(1:ncbcd)
                            do k = 1, natoc 
                                cityabcdbest(k) = tourtrial(idabcd(k))
                                tanglebest(k) = tangle(tourtrial(idabcd(k)))
                            end do
                            costbest = costtrial
                            jbest = j
                        end if

                        x(1:2,tourtrial(pos(1:ncbcd))) = xtmp(1:2,tourtrial(pos(1:ncbcd)))
                        do k = 1, natoc
                            tangle(tourtrial(idabcd(k))) = tangletmp(k)
                        end do
                    end if 
                end if 
            end do
            
            if ( jbest .gt. 0 ) then
                tour(npts+1-jbest+bsize:npts+bsize) = tour(npts+1-jbest:npts)
            end if
            
            if ( is_string ) then
                tour(npts+1-jbest:npts+bsize-jbest) = removed_cities(1:bsize)
            else
                tour(npts+1-jbest) = removed_cities(i)
            end if

            x(1:2,tour(pbest(1:ncbcd))) = xbest(1:2,1:ncbcd)
            do k = 1, natoc
                tangle(cityabcdbest(k)) = tanglebest(k)
            end do

            cost = costbest
            npts = npts + bsize
        end do

    end subroutine sisr_greedy_blink_insertion

    ! *****************************************************************
    ! sisr_reorder -- reinsertion order of the disks in the SISR insertion: random (with
    ! probability 4/7), from the farthest to the closest to the centroid of the centers
    ! (2/7), or from the closest to the farthest (1/7), the weights 4:2:1 of the random,
    ! "far" and "close" orders of Christiaens and Vanden Berghe (2020).
    subroutine sisr_reorder(inst,param,n,v)

        implicit none

        type(instance_type),intent(in) :: inst
        type(parameter_type),intent(inout) :: param
        integer,intent(in) :: n
        integer,intent(inout) :: v(n)

        integer :: i,j,tmp,modo
        real(kind=8) :: r,dk,cen(2),d(n)
        logical :: decrescente

        r = 7.0d0 * drand(param%seed)
        if ( r .lt. 4.0d0 ) then
            modo = 1
        else if ( r .lt. 6.0d0 ) then
            modo = 2
        else
            modo = 3
        end if
        sisr_orders(modo) = sisr_orders(modo) + 1

        if ( modo .eq. 1 ) then
            ! Fisher-Yates shuffle
            do i = n, 2, -1
                j = 1 + int( drand(param%seed) * dble(i) )
                if ( j .gt. i ) j = i
                tmp = v(i); v(i) = v(j); v(j) = tmp
            end do
            return
        end if

        cen(1) = sum( inst%centercoor(1,1:inst%nballs) ) / dble(inst%nballs)
        cen(2) = sum( inst%centercoor(2,1:inst%nballs) ) / dble(inst%nballs)
        do i = 1, n
            d(i) = norm2( inst%centercoor(1:2,v(i)) - cen )
        end do
        decrescente = ( modo .eq. 2 )
        do i = 2, n                              ! insertion sort (n is small)
            tmp = v(i); dk = d(i)
            j = i - 1
            do while ( j .ge. 1 )
                if ( decrescente ) then
                    if ( d(j) .ge. dk ) exit
                else
                    if ( d(j) .le. dk ) exit
                end if
                v(j+1) = v(j); d(j+1) = d(j)
                j = j - 1
            end do
            v(j+1) = tmp; d(j+1) = dk
        end do

    end subroutine sisr_reorder

    ! *****************************************************************
    ! *****************************************************************
    function hash_tour(n, tour) result(h)
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: tour(n)
        integer(kind=8) :: h
        
        integer :: i, u, v, kprev
        integer(kind=8) :: edge_hash
        integer(kind=8), parameter :: p1 = 2654435761_8
        integer(kind=8), parameter :: p2 = 2246822519_8

        h = 0_8
        do i = 1, n
            if (i .eq. 1) then
                kprev = n
            else
                kprev = i - 1
            end if
            
            u = tour(kprev)
            v = tour(i)
            
            if (u .gt. v) then
                u = tour(i)
                v = tour(kprev)
            end if
            
            edge_hash = ieor( int(u, 8) * p1, int(v, 8) * p2 )
            
            h = h + edge_hash
        end do
    end function hash_tour

    ! *****************************************************************
    ! *****************************************************************
    subroutine acceptsol(inst,param,count,start,newcost,newtour,newsol,newangles,&
    currentcost,currenttour,currentsol,currentangles,bestcost,bestour,bestsol,&
    bestangles,op_remov,op_insert,T,awarded_score)
        
        implicit none

        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(inout) :: param
        type(counter_type), intent(inout) :: count
        real(kind=8), intent(in) :: start,T
        real(kind=8), intent(inout) :: bestcost,currentcost,newcost
        integer, intent(in) :: op_remov,op_insert
        real(kind=8), intent(out) :: awarded_score

        integer, intent(in) :: newtour(inst%nballs)
        integer, intent(inout) :: bestour(inst%nballs),currenttour(inst%nballs)
        real(kind=8), intent(inout) :: currentsol(2,inst%nballs),newsol(2,inst%nballs),&
        currentangles(inst%nballs),newangles(inst%nballs),bestsol(2,inst%nballs),&
        bestangles(inst%nballs)
        
        integer :: iter_diff
        integer(kind=8) :: hash_key, h_idx, initial_idx
        real(kind=8) :: finish,time,delta,p,r
        logical :: accept,already_visited

        awarded_score = 0.0d0

        delta = newcost - currentcost   
        accept = .false.

        if ( delta .le. 0.0d0 ) then
            accept = .true.
        else
            if ( (-delta / T) .lt. -700.0d0 ) then
                p = 0.0d0
            else
                p = exp( - delta / T )
            end if
            r = drand(param%seed)
            if ( r .lt. p ) accept = .true.
        end if

        if ( accept ) then
            
            hash_key = hash_tour(inst%nballs, newtour)
            
            initial_idx = int( mod( abs( hash_key ), int( count%hashsize, 8 ) ) + 1_8 )
            h_idx = initial_idx
            already_visited = .false.

            do
                if ( count%hash_table(h_idx) .eq. -1_8 ) then
                    exit
                else if ( count%hash_table(h_idx) .eq. hash_key ) then
                    already_visited = .true.
                    exit
                end if
                h_idx = h_idx + 1
                if ( h_idx .gt. count%hashsize ) h_idx = 1
                if ( h_idx .eq. initial_idx ) exit
            end do

            if ( newcost .lt. bestcost ) then
                awarded_score = param%sigma1  
            else if ( .not. already_visited ) then
                if ( delta .le. 0.0d0 ) then
                    awarded_score = param%sigma2  
                else
                    awarded_score = param%sigma3  
                end if
            else
                awarded_score = 0.0d0 
            end if

            if ( .not. already_visited .and. count%hash_table(h_idx) .eq. -1_8 ) then
                count%hash_table(h_idx) = hash_key
                count%n_visited = count%n_visited + 1
            end if

            currentcost = newcost
            currenttour = newtour
            currentsol = newsol
            currentangles = newangles
            
            if ( newcost .lt. bestcost ) then
                count%remov_success(op_remov) = count%remov_success(op_remov) + 1
                count%insert_success(op_insert) = count%insert_success(op_insert) + 1

                iter_diff = count%alnsiter - count%last_improv_iter
                if ( iter_diff .gt. count%max_no_improv ) count%max_no_improv = iter_diff
                count%last_improv_iter = count%alnsiter
                count%nimprov = count%nimprov + 1

                bestcost = newcost
                bestour = newtour
                bestsol = newsol
                bestangles = newangles

                ! BCD on the new best solution, with at most bcdmaxcic*m = m block updates
                call angleBCD(inst,param,count,inst%nballs,bestour,bestsol,bestcost,bestangles)

                currentcost = bestcost
                currentsol = bestsol
                currentangles = bestangles

                call cpu_time(finish)
                time = finish - start
                count%last_improv_time = time 

                call save_best(inst,param,count,bestcost,time,'improvement',bestour,bestsol)
                write(*,*) ' iter = ',count%alnsiter,' Best sol = ',bestcost,' Time = ',time
                write(*,*) '      >>> IMPROVED BY: Removal OP = ',op_remov,', Insert OP = ',op_insert
                call flush(6)
            end if
        end if
    end subroutine acceptsol

    ! *****************************************************************
    ! *****************************************************************
    ! Writes the solution (cost, visiting order and visit points) to the file
    ! solutions/best-<instance>-seed-<seed>.txt. It is called with the initial
    ! solution, at each new best solution of the ALNS (after the BCD) and with the
    ! final solution, so that, if the run is interrupted (e.g., by the time limit),
    ! the file holds the last improvement. The file is first written to a temporary
    ! file, which is then renamed, so that an interruption during the writing never
    ! leaves a truncated file. The points are written in the original coordinates of
    ! the instance (the normalization is undone).
    !
    ! Format:  line 1: instance seed state cost iteration time(s)
    !          line 2: m
    !          m lines: disk (1..m, in the order of the instance file) x y
    subroutine save_best(inst,param,count,cost,time,state,tour,x)

        implicit none

        type(instance_type), intent(in) :: inst
        type(parameter_type), intent(in) :: param
        type(counter_type), intent(in) :: count
        real(kind=8), intent(in) :: cost,time
        character(len=*), intent(in) :: state
        integer, intent(in) :: tour(inst%nballs)
        real(kind=8), intent(in) :: x(2,inst%nballs)

        integer :: k,un,ios
        character(len=512) :: fn,fntmp

        write(fn,'(A,A,A,I0,A)') 'solutions/best-',trim(inst%name),'-seed-',&
        int(param%initial_seed),'.txt'
        fntmp = trim(fn) // '.tmp'

        open(newunit=un,file=trim(fntmp),status='replace',action='write',iostat=ios)
        if ( ios .ne. 0 ) then
            call execute_command_line('mkdir -p solutions')
            open(newunit=un,file=trim(fntmp),status='replace',action='write',iostat=ios)
            if ( ios .ne. 0 ) return
        end if

        write(un,'(A,1X,I0,1X,A,1X,ES24.16,1X,I0,1X,F14.3)') trim(inst%name),&
        int(param%initial_seed),state,cost,count%alnsiter,time
        write(un,'(I0)') inst%nballs
        do k = 1, inst%nballs
            write(un,'(I0,2(1X,ES24.16))') tour(k),&
            x(1:2,tour(k)) * inst%xymax + inst%xyshift(1:2)
        end do
        close(un)

        call rename(trim(fntmp),trim(fn))

    end subroutine save_best

end module metaheuristic
