module solvefem
    use iso_fortran_env, only: int32, real64
    use fem2d
    implicit none

contains
   
subroutine solve(M,P,b,co,Kval,aelem,boundary,FN,niter)
    ! TO-DO: 1. Create N festruct given nfestruct as input and store them in an array
    !        2. Define elemental functions/subroutines so that all the calculation is
    !           done simultaneously across all cilium in the festruct array
    integer, parameter :: dp = kind(0.d0)
    integer, intent(in) :: M(:,:)
    real(dp), intent(inout) :: P(:,:)
    real(dp), intent(in) :: b(:,:,:)
    logical, intent(in) :: boundary(:,:)
    real(dp), intent(in) :: co, kval
    real(dp), intent(in) :: aelem(:)
    real(dp), intent(inout) :: FN(:,:)
    integer, intent(in) :: niter
   
    real(dp) :: U(size(FN,1),size(FN,2))
    real(dp) :: mu, dt
    integer :: iter, ipoints, npoints
    real(dp) :: ftip
    type(festruct) :: cilium

    U = 0.0d0
    mu = 1000
    dt = 0.001
    npoints = size(P,1)
    ftip = 100.0d0
    cilium = festruct(M,P,b,boundary,aelem,FN,U,co,kval)

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