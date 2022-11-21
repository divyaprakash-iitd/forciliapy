subroutine ibmc(MM,PP,b,boundary,co,kval,aelem,FN)
    use mod_pressure,       only: generate_laplacian_sparse, calculate_pressure_sparse
    ! use mod_amgx,           only: calculate_pressure_amgx
    use mod_mesh
    use mod_time
    ! use mod_time_stepping
    ! use mod_boundary
    use mod_io
    use mod_ibm
    use mod_vec
    use fem2d
    implicit none

    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: MM(:,:)
    real(dp), intent(inout) :: PP(:,:)
    real(dp), intent(in) :: b(:,:,:)
    logical, intent(in) :: boundary(:,:)
    real(dp), intent(in) :: co, kval
    real(dp), intent(in) :: aelem(:)
    real(dp), intent(inout) :: FN(:,:)

    ! Parameters
    integer, parameter :: PI = 3.141592653589793
    ! Computational Domain
    real(dp)    :: Lx
    real(dp)    :: Ly
    ! Mesh Paramaters
    integer  :: Nx
    integer  :: Ny
    ! Simulation time Paramaters
    real(dp)    :: tsim
    real(dp)    :: dt  
    real(dp)    :: t
    integer  :: it_save
    ! Physical Constants
    real(dp)    :: nu 
    real(dp)    :: rho   
    ! Boundary values
    real(dp)    :: utop, vtop, ubottom, vbottom, &
                       uleft, vleft, uright, vright, ulid
    ! Matrices to store fields
    real(dp), allocatable :: u(:,:), v(:,:), us(:,:), vs(:,:), R(:,:), &
                                 P(:,:), A(:,:,:), Fx(:,:), Fy(:,:), umid(:,:), vmid(:,:)
    ! Temporary/Miscellaneous variable
    integer  :: it, NN, il, ip, err
    real(dp)    :: tp  ! Time period
    ! Mesh
    type(mesh)      :: M

    ! Namelists for input
    namelist /time/ dt, it_save, tsim
    namelist /grid/ Nx, Ny, Lx, Ly
    namelist /flow/ nu, rho, utop, TP
    
    type(festruct) :: cilium

    !---------------------- Begin Calculations ------------------------------------!

    ! Read input data from file
    open(1004,file="input_params.dat",form='formatted')
    READ(unit=1004,nml=grid,iostat=err)
    close(1004)

    ! Construct and write Mesh data
    M = mesh('M',Lx,Ly,Nx,Ny)
    call write_mesh(M,'u')
    call write_mesh(M,'v')
    call write_mesh(M,'p')

    ! Allocate matrices ofr u,v,us,vs,rhs
    allocate(u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(umid(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(vmid(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    allocate(R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub))
    allocate(Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub))    
    allocate(Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub))
    NN = 5
    allocate(A(1:Nx,1:Ny,NN))

    ! Initialize
    u   = 0.0d0
    v   = 0.0d0
    us  = 0.0d0
    vs  = 0.0d0
    R   = 0.0d0
    A   = 0.0d0
    Fx  = 0.0d0
    Fy  = 0.0d0

    ! Define boundary conditions for velocity
    utop    = 0.2d0
    vtop    = 0.0d0
    ubottom = 0.0d0
    vbottom = 0.0d0
    uleft   = 0.0d0
    vleft   = 0.0d0
    uright  = 0.0d0
    vright  = 0.0d0
    
    ! Generate Laplacian matrix
    call generate_laplacian_sparse(A,M%dx,M%dy)

    ! Read input data from file
    open(1004,file="input_params.dat",form='formatted')
    READ(unit=1004,nml=time,iostat=err)
    READ(unit=1004,nml=flow,iostat=err)
    close(1004)

    write(*,'(2(A,I8),2(A,1p1e15.6))') "Nx = ",Nx, " Ny = ",Ny, &
    " Lx = ",Lx, " Ly = ",Ly

    ! Define an festruct object before
    ! Initialize the festruct

    cilium = festruct(MM,PP,b,boundary,aelem,FN,co,kval)

    call apply_parabolic_initialization(M,u,uleft)

    ! Start time loop
    t = 0.0d0
    it = 0

    ! init_status = .False. ! AmgX initialization status

    ! Lid velocity
    ulid = utop

    do while (t.lt.tsim)
        ! ulid is the amplitude
        utop = 2*3.1415/tp*ulid*cos(2*3.1415*t/tp)
        t = t + dt
        it = it + 1
        

        ! RK2: Step 1
        ! Apply velocity boundary conditions
        call apply_boundary_channel(M,u,v,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
    
        call cilium%calculate_forces()
        Fx = 0.0d0 ! Initialize the forces at every time-step
        Fy = 0.0d0
        call spread_force(M,cilium,Fx,Fy)

        call cdu(M,u,v,us,nu,Fx)
        us = u + us*dt/2
        call cdv(M,u,v,vs,nu,Fy)
        vs = v + vs*dt/2

        ! Apply velocity boundary conditions to us and vs
        call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

        ! Form the RHS of the pressure poisson equation
        call calculate_rhs(M,us,vs,R,rho,0.5d0*dt)

        ! Solve for pressure
        ! call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
        call calculate_pressure_sparse(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:))

        ! Perform the corrector step to obtain the velocity
        call corrector(M,umid,vmid,us,vs,P,rho,0.5d0*dt)

        ! Interpolate velocity
        call interpolate_velocity(M,cilium,umid,vmid)

        ! Update structure
        call cilium%update_position(dt/2)

        ! RK2: Step 2
       
        ! Calculate forces
        call cilium%calculate_forces()
        call apply_boundary_channel(M,umid,vmid,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
        
        Fx = 0.0d0 ! Initialize the forces at every time-step
        Fy = 0.0d0
        call spread_force(M,cilium,Fx,Fy)

        call cdu(M,umid,vmid,us,nu,Fx)
        us = u + us*dt
        call cdv(M,umid,vmid,vs,nu,Fy)
        vs = v + vs*dt
        
        ! Apply velocity boundary conditions to us and vs
        call apply_boundary_channel(M,us,vs,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)

        ! Form the RHS of the pressure poisson equation
        call calculate_rhs(M,us,vs,R,rho,dt)

        ! Solve for pressure
        ! call calculate_pressure_amgx(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:),init_status)
        call calculate_pressure_sparse(A(2:M%Nx-1,:,:),P(2:M%Nx-1,:),R(2:M%Nx-1,:))

        ! Perform the corrector step to obtain the velocity
        call corrector(M,u,v,us,vs,p,rho,dt)
    
        ! Interpolate velocity
        call interpolate_velocity(M,cilium,u,v)

        ! Update structure
        call cilium%update_position(dt)

        write(*,'(A,F10.5)') 'time = ', t
        
        ! Write files every Nth timestep
        if (mod(it,it_save).eq.0) then 
            call write_field(u,'u',it) 
            call write_field(v,'v',it) 
        end if

    end do

end subroutine ibmc
