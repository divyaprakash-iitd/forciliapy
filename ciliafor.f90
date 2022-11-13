module calculate_forces
    implicit none
    
    private

    public :: deformationgradient, jacobian, invariantderivative, jacobianderivative
    
contains

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

! subroutine deformationgradient(M,P,b,D)
function deformationgradient(M,P,b) result(D)
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: M(:,:)
    real(dp), intent(in) :: P(:,:)
    real(dp), intent(in) :: b(:,:,:)
    ! real(dp), intent(out) :: D(size(M,1),2,2)
    real(dp) :: D(size(M,1),2,2)

    integer :: i, j, nelem, idxelem, nvertex

    real(dp) :: X(2,3)

    nelem = size(M,1)
    nvertex = 3

    do i = 1, nelem
        do concurrent (j = 1:nvertex)
            idxelem = M(i,j)
            X(:,j) = P(idxelem,:)
        end do
        D(i,:,:) = matmul(X, transpose(b(i,:,:)))
    end do
    
! end subroutine deformationgradient
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
    
end module calculate_forces