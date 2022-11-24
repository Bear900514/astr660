subroutine evolution()
    use Simulation_data
    use IO, only : output
    implicit none

    integer :: n
    real    :: dt, time

    n        = 0
    time     = 0.0

    dt = abs(dx/c)*cfl

    do while(time .le. tend)

        ! reset boundary condition
        call boundary(u)

        ! dump output times with frequency set by io_interval
        if (mod(n,io_interval) .eq. 0) then
            print *, "n =", n ," Time =", time
            call output(n,time)
        endif

        ! update the solution
        call update(time, dt)
        
        n    = n + 1
        time = time + dt
    enddo

end subroutine evolution

subroutine update(time, dt)
    use Simulation_data
    implicit none
    real, intent(in) :: time ,dt
    integer :: i
    real :: FL,FR

    uold = u

    do i = istart, iend
        ! upwind
        !u(i)=uold(i)-c*dt*(uold(i)-uold(i-1))/(dx)
        ! TODO: implement the upwind or FTCS methods here
        ! FTCS
        !u(i)=uold(i)-c*dt*(uold(i+1)-uold(i-1))/(dx*2)
        call flux(i,dt,FL,FR)
        u(i)=uold(i)-dt/dx*(FR-FL)
    enddo

end subroutine update

subroutine flux(i,dt,FL,FR)
    use Simulation_data
    implicit none
    real, intent(in) :: dt
    integer, intent(in) :: i
    real, intent(out) :: FL,FR
    !Arithmetic
    !FL=0.5*c*(uold(i-1)+uold(i))
    !FR=0.5*c*(uold(i)+uold(i+1))
    !The Lax-Friedrichs
    FL=0.5*c*(uold(i-1)+uold(i))-dx*(uold(i)-uold(i-1))/(2*dt)
    FR=0.5*c*(uold(i)+uold(i+1))-dx*(uold(i+1)-uold(i))/(2*dt)

end subroutine flux
