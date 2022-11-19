module solvefem
    use calculate_forces
    use mod_cilia
    implicit none

contains

subroutine computeforces(M,P,b,co,Kval,aelem,FN)
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: M(:,:)
    real(dp), intent(in) :: P(:,:)
    real(dp), intent(in) :: b(:,:,:)
    real(dp), intent(in) :: co, kval
    real(dp), intent(in) :: aelem(size(M,1))
    real(dp), intent(inout) :: FN(:,:)

    real(dp) :: D(size(M,1),2,2)
    real(dp) :: didx(size(M,1),2,3)
    real(dp) :: djdx(size(M,1),2,3)
    real(dp) :: jac(size(M,1))
    real(dp) :: F(size(M,1),2,3)
    integer :: j, k, nelem, nidx

    nelem = size(M,1)

    ! Calculate deformation gradient
    D = deformationgradient(M,P,b)

    ! Calculate the invariant derivative term
    didx = invariantderivative(D,b)

    ! Calculate the jacobian derivative term
    djdx = jacobianderivative(D,b)

    ! Calculate the jacobian
    jac = jacobian(D)

    ! Calculate forces at each node per element 
    ! (3 nodes, 2 components of force per element)
    F = 0.0d0
    do k = 1, nelem
        F(k,:,:) = -0.5d0 * co * aelem(k) * &
                    ( didx(k,:,:) - 2.0d0/jac(k) * djdx(k,:,:) * &
                    (1.0d0 - kval/co * log(jac(k)))) 
    end do

    ! Adding forces for all nodes
    do k = 1, nelem
        do j = 1, 3
            nidx = M(k,j)
            FN(nidx,1) = FN(nidx,1) + F(k,1,j)            
            FN(nidx,2) = FN(nidx,2) + F(k,2,j)            
        end do
    end do 

end subroutine computeforces
   
subroutine solve(M,P,b,co,Kval,aelem,boundary,FN,mu,dt)
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: M(:,:)
    real(dp), intent(inout) :: P(:,:)
    real(dp), intent(in) :: b(:,:,:)
    logical, intent(in) :: boundary(:,:)
    real(dp), intent(in) :: co, kval, mu, dt
    real(dp), intent(in) :: aelem(:)
    real(dp), intent(inout) :: FN(:,:)

    integer(int32) :: iter, niter, ipoints, npoints
    real(real64) :: ftip
    type(cilia) :: cilium

    niter = 10000
    npoints = size(P,1)
    ftip = 100.0d0
    cilium = cilia(M,P,b,boundary,aelem,FN,co,kval)

    print *, cilium%co

    do iter = 1, niter
        cilium%FN = 0.0d0
        ! forall (ipoints = 1:npoints, boundary(ipoints,1).eqv..True.) &
            ! FN(ipoints,1) = ftip
        where (boundary(:,1).eqv..True.)
            cilium%FN(:,1) = ftip/3.0d0
        end where
       

        call cilium%calculate_forces()

        ! Fix the bottom nodes (Make force zero)
        forall (ipoints = 1:npoints, boundary(ipoints,2).eqv..True.) &
            cilium%FN(ipoints,2) = 0.0d0

        call cilium%update_position(dt,mu)

        ! Fix the bottom nodes (Make displacement zero)
        forall (ipoints = 1:npoints, boundary(ipoints,2).eqv..True.) &
            cilium%P(ipoints,2) = 0.0d0

    end do

    FN = cilium%FN
    P = cilium%P

end subroutine solve

end module solvefem