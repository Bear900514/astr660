!---------------------------------------------------
! The physics module
!
module physics
    use Simulation_data
    implicit none
    contains
        
        subroutine initial()
            !
            ! setup initial conditions of each stars
            ! in this example we only have two stars
            !
            use constants, only : au, msun, pi, G
            implicit none
            integer :: i
            real :: m1, m2, force

            m1 = 1.0 * msun
            m2 = 2.0 * msun

            !
            ! Use Kepler's law to evaluate the orbital period
            ! and use Newton's law to evaluate the force between
            ! two stars
            !

            separation = 3.0*au !TODO
            period     = ((4.0*pi**2*separation**3)/(G*(m1+m2)))**(1.0/2.0) !TODO
            force      = G*m1*m2/(separation**2.0) !TODO

            !
            ! setup initial conditions of star M1
            stars(1)%mass = m1
            stars(1)%x    = -2.0*au !TODO
            stars(1)%y    = 0.0 !TODO
            stars(1)%vx   = 0.0 !TODO
            stars(1)%vy   = 2.0*pi*2*au/period
            !stars(1)%vy   = 2000000.0 !TODO
            stars(1)%ax   = force/stars(1)%mass !TODO
            stars(1)%ay   = 0.0 !TODO

            !
            ! setup initial conditions of star M2
            stars(2)%mass = m2
            stars(2)%x    = 1.0*au !TODO
            stars(2)%y    = 0.0 !TODO
            stars(2)%vx   = 0.0 !TODO
            stars(2)%vy   = -2.0*pi*1*au/period
            !stars(2)%vy   = -1000000.0 !TODO
            stars(2)%ax   = -force/stars(2)%mass !TODO
            stars(2)%ay   = 0.0 !TODO


        end subroutine initial

        subroutine update(dt)
            use constants
            implicit none
            real, intent(in)  :: dt
            integer :: i,j
            real    :: x, y, rsq, fx, fy
            real    :: radius, force, angle ,distance
            do i=1,2
                !angle = 
                !
                ! In this example we use a first order scheme (Euler method)
                ! we approximate dx/dt = v  --> x^(n+1) - x^n = v * dt
                ! therefore, x at step n+1 = x^(n+1) = x^n + v * dt
                !
                ! the same approximation can be applied to dv/dt = a
                !
                x=stars(i)%x
                y=stars(i)%y
                ! update position to t = t + dt
                !TODO
                stars(i)%x = stars(i)%x + stars(i)%vx * dt!TODO
                stars(i)%y = stars(i)%y + stars(i)%vy * dt!TODO
                ! update velocity to t = t + dt
                !TODO
                stars(i)%vx = stars(i)%vx + stars(i)%ax * dt
                stars(i)%vy = stars(i)%vy + stars(i)%ay * dt
                ! update accelerations to t = t + dt
                !TODO
               
                !fx = force*(x/(2.0*au/i))
                !fx = force*(y/(2.0*au/i))!(y/(2.0*au/i)) 
            enddo
            distance=((stars(2)%x-stars(1)%x)**2.0+(stars(2)%y-stars(1)%y)**2.0)**(1.0/2.0)
            !print *,distance
            force = G*stars(1)%mass*stars(2)%mass/distance**2 
            fx = force*(stars(2)%x-stars(1)%x)/distance
            fy = force*(stars(2)%y-stars(1)%y)/distance
            stars(1)%ax = fx/stars(1)%mass
            stars(1)%ay = fy/stars(1)%mass
            stars(2)%ax = -fx/stars(2)%mass
            stars(2)%ay = -fy/stars(2)%mass
            return
        end subroutine update

end module physics

