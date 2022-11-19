module fem2d
    ! TO-DO: 1. Name the module fem2d
    !        2. Move all the required subroutines here
    !           (Not necessarily type bound methods)
    
    use iso_fortran_env, only: int32, real64
    ! use calculate_forces
    implicit none
    
    type :: festruct
        integer(int32), allocatable :: M(:,:) ! Indices providing connectivity information
        real(real64), allocatable :: P(:,:) ! Coordinates of points
        real(real64), allocatable :: b(:,:,:) ! Shape function coefficients
        logical, allocatable :: boundary(:,:) ! Top, bottom, left and right boundaries
        real(real64), allocatable :: aelem(:)
        real(real64), allocatable :: FN(:,:)
        real(real64) :: co, kval

        contains
            procedure :: calculate_forces
            procedure :: update_position

    end type festruct

    contains
    
    subroutine calculate_forces(self)
        class(festruct), intent(inout) :: self
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
    
    subroutine update_position(self, dt, mu)
        class(festruct), intent(inout) :: self
        real(real64), intent(in) :: dt, mu

        integer(int32) :: i, npoints

        npoints = size(self%P,1)

        do i = 1,size(self%P,1)
            self%P(i,1) = self%P(i,1) + dt/mu * self%FN(i,1)
            self%P(i,2) = self%P(i,2) + dt/mu * self%FN(i,2)
        end do
            
    end subroutine update_position

    ! Helper functions
    ! Jacobian derivative term
    function jacobianderivative(D,b) result(djdx)
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: D(:,:,:)
        real(dp), intent(in) :: b(:,:,:)
        real(dp) :: djdx(size(b,1),2,3)

        integer :: k, nelem
        real(dp) :: MD(2,2)

        nelem = size(b,1)

        do concurrent (k = 1:nelem)
            ! Rearrange the elements of the deformation matrix
            MD(1,1) = D(k,2,2)
            MD(1,2) = -D(k,2,1)
            MD(2,1) = -D(k,1,2)
            MD(2,2) = D(k,1,1)

            djdx(k,:,:) = matmul(MD, b(k,:,:))
        end do

    end function jacobianderivative 

    ! Invariant derivative term
    function invariantderivative(D,b) result(didx)
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: D(:,:,:)
        real(dp), intent(in) :: b(:,:,:)
        real(dp) :: didx(size(b,1),2,3)

        integer :: k, nelem

        nelem = size(b,1)
        
        do concurrent (k = 1:nelem)
            didx(k,:,:) = 2 * matmul(D(k,:,:), b(k,:,:)) 
        end do

    end function invariantderivative

    function deformationgradient(M,P,b) result(D)
        integer, parameter :: dp = kind(0.d0)
        integer, intent(in) :: M(:,:)  ! Triangulation
        real(dp), intent(in) :: P(:,:) ! Points
        real(dp), intent(in) :: b(:,:,:) ! Shape function coefficient
        real(dp) :: D(size(M,1),2,2) ! Deformation gradient

        integer :: i, j, nelem, idxelem, nvertex

        real(dp) :: X(2,3) ! Coordinates matrix of nodes of a triangle

        nelem = size(M,1) ! No of elements 
        nvertex = 3 ! No of nodes in an element\triangle 

        do i = 1, nelem
            do concurrent (j = 1:nvertex)
                idxelem = M(i,j)
                X(:,j) = P(idxelem,:)
            end do
            D(i,:,:) = matmul(X, transpose(b(i,:,:)))
        end do
    end function deformationgradient

    function jacobian(F) result(J)
    ! subroutine jacobian(F,J)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: F(:,:,:)
        real(dp)             :: J(size(F,1))
        ! real(dp), intent(out)             :: J(size(F,1))

        integer             :: k, nelem

        nelem = size(F,1)
        do k = 1, nelem
            J(k) = determinant(F(k,:,:))
        end do

    end function jacobian
    ! end subroutine jacobian

    function determinant(M) result(D)
        implicit none
        integer, parameter :: dp = kind(0.d0)
        real(dp), intent(in) :: M(:,:)
        real(dp) :: D

        D = M(1,1) * M(2,2) - M(1,2) * M(2,1)

        ! D = M(1,1) * (M(2,2) * M(3,3) - M(3,2) * M(2,3)) - &
        !     M(1,2) * (M(2,1) * M(3,3) - M(3,1) * M(2,3)) + &
        !     M(1,3) * (M(2,1) * M(3,2) - M(2,2) * M(3,1))

    end function determinant

end module fem2d