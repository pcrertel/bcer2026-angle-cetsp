program alns_cetsp

    use instance, only: instance_type,parameter_type,counter_type,inputdata
    use metaheuristic, only: alns,initial_solution,save_best

    implicit none

    ! LOCAL SCALARS
    logical :: normalize
    integer :: allocerr
    real(kind=8) :: ccost,cost0,start,finish,time
    type(instance_type) :: inst
    type(parameter_type) :: param
    type(counter_type) :: count

    ! LOCAL ARRAYS
    integer, allocatable :: ctour(:),tour0(:)
    real(kind=8), allocatable :: cangles(:),angles0(:),csol(:,:),xsol0(:,:)

    normalize = .true.    ! normalized coordinates and radii (see inputdata)

    call inputdata(inst,param,count,normalize)

    allocate(ctour(inst%nballs),csol(2,inst%nballs),&
    cangles(inst%nballs),angles0(inst%nballs),&
    tour0(inst%nballs),xsol0(2,inst%nballs),stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error - Main.'
        stop
    end if

    ! output folder of the best solutions, if it does not exist yet
    call execute_command_line('mkdir -p solutions')

    ! initial solution: convex hull + cheapest insertion, with the visit point of each
    ! inserted disk given by the single-disk subproblem, followed by the BCD
    write(*,*) 'Initial solution (convex hull + cheapest insertion + BCD)'
    call initial_solution(inst,param,count,tour0,xsol0,cost0,angles0)
    count%costhc = cost0
    write(*,*) 'Initial cost = ',cost0
    call save_best(inst,param,count,cost0,0.0d0,'initial',tour0,xsol0)

    ! ALNS-BCD
    ctour = tour0
    csol = xsol0
    ccost = cost0
    cangles = angles0

    call cpu_time(start)
    call alns(inst,param,count,ctour,csol,ccost,cangles,start,&
    param%nworst,param%alnsmaxiter)
    call cpu_time(finish)

    time = finish - start
    call save_best(inst,param,count,ccost,time,'final',ctour,csol)

    write(*,*) '=================================='
    write(*,*) 'Finished ALNS.'
    write(*,*) 'Cost = ',ccost,' Time = ',time
    write(*,*) '=================================='
    call flush(6)

    deallocate(inst%centercoor,inst%radii,inst%distmatrix,inst%nn_list,&
    inst%min_angle_matrix,count%hash_table,ctour,csol,tour0,xsol0,cangles,angles0,&
    stat=allocerr)

    if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error - Main.'
        stop
    end if

end program alns_cetsp
