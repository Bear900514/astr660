!---------------------------------------------------
!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.21
!
! Problem:
!
!        Simulating binary orbits
!
program binary

    use constants
    use Simulation_data
    use physics, only : initial, update
    use output, only : record

    implicit none

    integer :: i, interval, noutput, step
    real :: dt, time, tmax

    ! initial setup
    step     = 0             ! start from 0th step
    dt       = 0.00001 * yr  ! sec, time step
    time     = 0.0           ! sec, start from t=0
    tmax     = 10.0 * yr     ! sec, max simulation time

    ! output intervals to record trajectory
    noutput  = 500
    interval = (tmax/dt)/noutput

    ! initialize the stars
    call initial()

    ! the main time loop
    do while (time .le. tmax)

        ! update stars
        call update(dt)

        ! update time
        time = time + dt

        ! output data
        if (mod(step,interval) .eq. 0) then
            call record(time)
        endif

        ! update step
        step = step + 1
    end do
    
    print *, "DONE!"

end program binary

