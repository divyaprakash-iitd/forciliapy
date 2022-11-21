module mod_time_stepping
    use iso_fortran_env, only: int32, real64
    use mod_pressure,       only: generate_laplacian_sparse, calculate_pressure_sparse 
    use mod_mesh
    use mod_time
    use mod_boundary
    use mod_io
    use mod_ibm
    use mod_vec
    implicit none
    
contains

    subroutine time_loop(FP,BC,M,u,v,us,vs,Fx,Fy,A,P,R,tsim,dt,it_save,tp)
        real(real64), intent(in)          :: FP(:)
        real(real64), intent(in)          :: BC(:)
        class(mesh), intent(inout)        :: M
        real(real64), intent(inout)       :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: us(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vs(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(inout)       :: Fx(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub)
        real(real64), intent(inout)       :: Fy(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)
        real(real64), intent(in)          :: A(:,:,:)
        real(real64), intent(inout)       :: P(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(inout)       :: R(M%xp%lb:M%xp%ub,M%yp%lb:M%yp%ub)
        real(real64), intent(in)          :: tsim
        real(real64), intent(in)          :: dt
        integer(int32), intent(in)        :: it_save
        real(real64), intent(in)          :: tp

        ! Intermediate velocities
        real(real64)       :: umid(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), vmid(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Iteration 
        integer(int32) :: it

        ! Time
        real(real64) :: t

        ! Fluid properties
        real(real64) :: nu, rho        

        ! Boundary conditions
        real(real64)    :: utop, vtop, ubottom, vbottom, &
                               uleft, vleft, uright, vright, ulid
        
        ! AmgX
        logical         :: init_status

        ! Assign boundary values
        utop    = BC(1)
        vtop    = BC(2)
        ubottom = BC(3)
        vbottom = BC(4)
        uleft   = BC(5)
        vleft   = BC(6)
        uright  = BC(7)
        vright  = BC(8)


        ! Get fluid properties
        nu  = FP(1)
        rho = FP(2)

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
        
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0

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

            ! RK2: Step 2
            
            call apply_boundary_channel(M,umid,vmid,utop,ubottom,uleft,uright,vtop,vbottom,vleft,vright)
            
            Fx = 0.0d0 ! Initialize the forces at every time-step
            Fy = 0.0d0

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
        
            write(*,'(A,F10.5)') 'time = ', t
            
            ! Write files every Nth timestep
            if (mod(it,it_save).eq.0) then 
                call write_field(u,'u',it) 
                call write_field(v,'v',it) 
            end if

        end do
    end subroutine time_loop

end module mod_time_stepping
