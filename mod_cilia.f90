module mod_cilia
    use iso_fortran_env, only: int32, real64
    use calculate_forces
    implicit none
    
    type :: cilia
        integer(int32), allocatable :: M(:,:) ! Indices providing connectivity information
        real(real64), allocatable :: P(:,:) ! Coordinates of points
        real(real64), allocatable :: b(:,:,:) ! Shape function coefficients
        logical, allocatable :: boundary(:,:) ! Top, bottom, left and right boundaries
        real(real64), allocatable :: aelem(:)
        real(real64), allocatable :: FN(:,:)
        real(real64) :: co, kval

        contains
            procedure :: calculate_forces, update_position

    end type cilia

    contains

    subroutine update_position(self, dt, mu)
        class(cilia), intent(inout) :: self
        real(real64), intent(in) :: dt, mu

        integer(int32) :: i, npoints

        npoints = size(self%P,1)

        do i = 1,size(self%P,1)
            self%P(i,1) = self%P(i,1) + dt/mu * self%FN(i,1)
            self%P(i,2) = self%P(i,2) + dt/mu * self%FN(i,2)
        end do
    
        
    end subroutine update_position

    subroutine calculate_forces(self)
        class(cilia), intent(inout) :: self
        real(real64) :: D(size(self%M,1),2,2)
        real(real64) :: didx(size(self%M,1),2,3)
        real(real64) :: djdx(size(self%M,1),2,3)
        real(real64) :: jac(size(self%M,1))
        real(real64) :: F(size(self%M,1),2,3)
        integer :: j, k, nelem, nidx

        nelem = size(self%M,1)

        ! Calculate deformation gradient
        D = deformationgradient(self%M,self%P,self%b)

        ! Calculate the invariant derivative term
        didx = invariantderivative(D,self%b)

        ! Calculate the jacobian derivative term
        djdx = jacobianderivative(D,self%b)

        ! Calculate the jacobian
        jac = jacobian(D)

        ! Calculate forces at each node per element 
        ! (3 nodes, 2 components of force per element)
        F = 0.0d0
        do k = 1, nelem
            F(k,:,:) = -0.5d0 * self%co * self%aelem(k) * &
                        ( didx(k,:,:) - 2.0d0/jac(k) * djdx(k,:,:) * &
                        (1.0d0 - self%kval/self%co * log(jac(k)))) 
        end do

        ! Adding forces for all nodes
        do k = 1, nelem
            do j = 1, 3
                nidx = self%M(k,j)
               self%FN(nidx,1) = self%FN(nidx,1) + F(k,1,j)            
               self%FN(nidx,2) = self%FN(nidx,2) + F(k,2,j)            
            end do
        end do 

    end subroutine calculate_forces


end module mod_cilia