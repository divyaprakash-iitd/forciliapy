module mod_ibm
    use iso_fortran_env, only: int32, real64, real64
    use mod_mesh
    use fem2d
    implicit none
   
    real(real64), parameter :: PI = 3.141592653589793

    private
    public :: spread_force, interpolate_velocity

contains
   
    subroutine spread_force(M,B,Fx,Fy)
        class(mesh), intent(in) :: M
        ! class(ib), intent(in) :: B
        class(festruct), intent(in) :: B
        real(real64), intent(in out) :: Fx(:,:), Fy(:,:)

        ! Define variables to store locations
        real(real64) :: Lx, Ly      ! Lagrangian locations
        real(real64) :: Ex, Ey      ! Eulerian locations
        real(real64) :: Flx, Fly    ! Forces at Lagrangian locations
        integer(int32) :: np ! Number of Lagrangian particles

        ! Indices
        integer(int32) :: i, j, inp, LeftIndex, RightIndex, BottomIndex, TopIndex, StencilSize

        np = size(B%P,1)
        StencilSize = 4; ! Only on one side

        ! Calculate the x-direction force on the u-velocity cells
        ! Iterating over all the grid points including the boundary values
        do inp = 1,np
            Lx = B%P(inp,1)
            Ly = B%P(inp,2)
            Flx = B%FN(inp,1)
            ! Make sure that the following indices lie inside the domain            
            LeftIndex   = 0 + M%xu%lb + floor((Lx-StencilSize*M%dx)/M%dx) 
            BottomIndex = 1 + M%yu%lb + floor((Ly-StencilSize*M%dy)/M%dy)   ! Increment by 1 since the first value is negative and not zero

                    ! print *, floor((Ly-StencilSize*M%dy)/M%dy)
            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%u_mesh(i,j)%x
                    Ey = M%u_mesh(i,j)%y
                    
                    Fx(i,j) = Fx(i,j) + Flx * dirac( [(Ex-Lx), (Ey-Ly)], M%dx)
                end do
            end do
        end do

        ! Calculate the y-direction force on the v-velocity cells
        do inp = 1,np
            Lx = B%P(inp,1)
            Ly = B%P(inp,2)
            Fly = B%FN(inp,2)

            LeftIndex   = 1 + M%xv%lb + floor((Lx-StencilSize*M%dx)/M%dx) ! Increment by 1 since the first value is negative and not zero 
            BottomIndex = 0 + M%yv%lb + floor((Ly-StencilSize*M%dy)/M%dy)

            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%v_mesh(i,j)%x
                    Ey = M%v_mesh(i,j)%y
                    
                    Fy(i,j) = Fy(i,j) + Fly * dirac( [(Ex-Lx), (Ey-Ly)], M%dy)
                end do
            end do
        end do

        contains

        pure function dirac(x,h)
            ! Defined for a uniform grid
            real(real64), intent(in) :: x(2)
            real(real64), intent(in) :: h
            real(real64) :: dirac

            integer(int32) :: ii
            real(real64) :: phi, r

            dirac = 1.0d0
            do ii = 1,2
                r = x(ii)/h 
                if (abs(r).le.2) then
                    phi = 0.25d0 * (1 + cos(PI*r/2))
                else 
                    phi = 0.0d0
                end if
            
                dirac = dirac*phi
            end do        

            dirac = (1/h**2) * dirac

        end function dirac

    end subroutine spread_force

    subroutine interpolate_velocity(M,B,u,v)
        class(mesh), intent(in) :: M
        ! class(ib), intent(in out) :: B
        class(festruct), intent(inout) :: B
        real(real64), intent(in) :: u(M%xu%lb:M%xu%ub,M%yu%lb:M%yu%ub), v(M%xv%lb:M%xv%ub,M%yv%lb:M%yv%ub)

        ! Define variables to store locations
        real(real64) :: Lx, Ly      ! Lagrangian locations
        real(real64) :: Ex, Ey      ! Eulerian locations
        real(real64) :: UL, VL      ! Velocity at Lagrangian locations
        real(real64) :: UE, VE      ! Velocity at Eulerian locations
        integer(int32) :: np        ! Number of Lagrangian particles
        
        ! Indices
        integer(int32) :: i, j, inp, LeftIndex, RightIndex, BottomIndex, TopIndex, StencilSize

        np = size(B%P,1)
        StencilSize = 4

        ! Calculate the u velocity of Lagrangian points
        ! Iterating over all the grid points including the boundary values for u-velocity cells
        do inp = 1,np
            Lx = B%P(inp,1)
            Ly = B%P(inp,2)
        
            ! Make sure that the following indices lie inside the domain            
            LeftIndex   = 0 + M%xu%lb + floor((Lx-StencilSize*M%dx)/M%dx) 
            BottomIndex = 1 + M%yu%lb + floor((Ly-StencilSize*M%dy)/M%dy)   ! Increment by 1 since the first value is negative and not zero

            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%u_mesh(i,j)%x ! u-cell x location
                    Ey = M%u_mesh(i,j)%y ! u-cell y location
                    UE = u(i,j)
                    B%U(inp,1) = B%U(inp,1) + UE * dirac( [(Ex-Lx), (Ey-Ly)], M%dx) * M%dx**2
                end do
            end do
        end do

        ! Calculate the v velocity of Lagrangian points
        ! Iterating over all the grid points including the boundary values for v-velocity cells
        
        do inp = 1,np
            Lx = B%P(inp,1)
            Ly = B%P(inp,2)
       
            LeftIndex   = 1 + M%xv%lb + floor((Lx-StencilSize*M%dx)/M%dx) ! Increment by 1 since the first value is negative and not zero 
            BottomIndex = 0 + M%yv%lb + floor((Ly-StencilSize*M%dy)/M%dy)

            RightIndex  = LeftIndex     +   2*StencilSize
            TopIndex    = BottomIndex   +   2*StencilSize

            do j = BottomIndex,TopIndex
                do i = LeftIndex,RightIndex
                    Ex = M%v_mesh(i,j)%x ! v-cell x location
                    Ey = M%v_mesh(i,j)%y ! v-cell y location
                    VE = v(i,j)
                    B%U(inp,2) = B%P(inp,2) + VE * dirac( [(Ex-Lx), (Ey-Ly)], M%dy) * M%dy**2
                end do
            end do
        end do

        contains 

        function dirac(x,h)
            ! Defined for a uniform grid
            real(real64), intent(in) :: x(2)
            real(real64), intent(in) :: h
            real(real64) :: dirac

            integer(int32) :: ii
            real(real64) :: phi, r

            dirac = 1.0d0
            do ii = 1,2
                r = x(ii)/h 
                if (abs(r).le.2) then
                    phi = 0.25d0 * (1 + cos(PI*r/2))
                else 
                    phi = 0.0d0
                end if
            
                dirac = dirac*phi
            end do        

            dirac = (1/h**2) * dirac

        end function dirac

    end subroutine interpolate_velocity

    end module mod_ibm

