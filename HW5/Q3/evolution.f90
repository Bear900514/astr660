subroutine evolution()
    use Simulation_data
    use IO, only : output
    implicit none

    integer :: n
    integer :: interval
    real    :: dt, time


    n        = 0
    time     = 0.0

    dt = abs(dx/c)*cfl

    do while(time .le. tend)

        ! reset boundary condition

        !call boundary(u)
        if (mod(n,io_interval) .eq. 0) then
            print *, "n =", n ," Time =", time
            call output(n,time)
        endif
        call update(time, dt)
        
        n    = n + 1
        time = time + dt
    enddo

end subroutine evolution


!!
!! 
!!
subroutine update(time, dt)
    use Simulation_data
    implicit none
    real, intent(in) :: time ,dt
    integer :: i, j
    real    :: FL, FR

    ! 1st order in time

    call boundary_x(u)
    uold  = u
    do j = istart, iend
        do i = istart, iend
            call flux_x(i,j,dt,FL,FR)
            u(i,j) = uold(i,j) - dt/dx*(FR-FL)
        enddo
    enddo
    call boundary_y(u)
    uold  = u
    do j = istart, iend
        do i = istart, iend
            call flux_y(i,j,dt,FL,FR)
            u(i,j) = uold(i,j) - dt/dx*(FR-FL)
        enddo
    enddo
end subroutine update

!
! Routine to evalue flux the cell edge
!
subroutine flux_x(i,j,dt,FL, FR)
    use Simulation_data
    implicit none
    integer, intent(in) :: i,j
    real, intent(in)    :: dt
    real, intent(out)   :: FL, FR

    real :: sig, a, b, qL, qR


    ! Arithmetic average method
    !FL = 0.5*c*(uold(i-1)+uold(i))
    !FR = 0.5*c*(uold(i+1)+uold(i))

    ! The Lax-Friedrichs Method
    !FL = 0.5*c*(uold(i-1)+uold(i)) -0.5*dx/dt*(uold(i)-uold(i-1))
    !FR = 0.5*c*(uold(i+1)+uold(i)) -0.5*dx/dt*(uold(i+1)-uold(i))

    ! The upwind method
    !FL = ! TODO
    !FR = ! TODO

    ! The Lax-Wendroff Method
    !FL = ! TODO
    !FR = ! TODO

    !! Use piecewise linear and slope limiter

    !! left state
    call get_slope(dx,uold(i-2,j),uold(i-1,j),uold(i,j),sig) ! compute sig(i-1)
    FL = uold(i-1,j)*c+0.5*c*(dx-dt*c)*sig! TODO

    !! right state
    call get_slope(dx,uold(i-1,j),uold(i,j),uold(i+1,j),sig) ! compute sig(i)
    FR = uold(i,j)*c+0.5*c*(dx-dt*c)*sig! TODO

    return

end subroutine flux_x

subroutine flux_y(i,j,dt,FL, FR)
    use Simulation_data
    implicit none
    integer, intent(in) :: i,j
    real, intent(in)    :: dt
    real, intent(out)   :: FL, FR

    real :: sig, a, b, qL, qR

    !! Use piecewise linear and slope limiter
    !! left state
    call get_slope(dx,uold(i,j-2),uold(i,j-1),uold(i,j),sig) ! compute sig(i-1)
    FL = uold(i,j-1)*c+0.5*c*(dx-dt*c)*sig! TODO

    !! right state
    call get_slope(dx,uold(i,j-1),uold(i,j),uold(i,j+1),sig) ! compute sig(i)
    FR = uold(i,j)*c+0.5*c*(dx-dt*c)*sig! TODO

    return

end subroutine flux_y

subroutine get_slope(dx,l,m,r,sig)
    implicit none
    real, intent(in)  :: dx
    real, intent(in)  :: l   ! left 
    real, intent(in)  :: m   ! middle
    real, intent(in)  :: r   ! right
    real, intent(out) :: sig ! the slope
    real :: a, b

    ! compute a and b as the left/right slopes 
    a = (m-l)/dx
    b = (r-m)/dx

    ! TODO: implement the minmod limiter
    if (abs(a).lt.abs(b) .and. (a*b).gt.0) then
        sig = a
    else if (abs(a).gt.abs(b) .and. (a*b).gt.0) then
        sig = b
    else if ((a*b) .le. 0) then
        sig = 0.0
    endif
    return
end subroutine get_slope
