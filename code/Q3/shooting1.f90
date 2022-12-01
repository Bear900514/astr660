!---------------------------------------------------
!
! National Tsing Hua University
!
! ASTR 660 Computational Astrophysics
!
! Created:  Kuo-Chuan Pan 2020
! Modified: Karen Yang 2022.10.25
!
! Problem:
!
!        Solving boundary value problems

program shooting1

    use Solver , only : rk4
    implicit none
    external :: my_func
    real, parameter :: pi = 4.0*atan(1.0)
    real :: h,t,tend
    integer :: i,n
    real :: try_y2,theta,vx,x,ans,fa,fx,err
    real, save :: a = 0.0001   ! bracking interval [a,b]
    real, save :: b = pi/6.0
    real, dimension(2) :: y, ynext
    err=10.0
    do while(err .ge. 0.0000001)
        ans= (a+b)/2.0
        fa=(2.0*30**2*cos(a)*sin(a)/9.8)-50
        fx=(2.0*30**2*cos(ans)*sin(ans)/9.8)-50
        err= abs(fx)
        if(sign(1.,fx)/=sign(1.,fa)) then
            b=ans
        else
            a=ans
        endif
    enddo
    print *,ans, fx ,fa, err
    ! a trial value
    theta=ans
    try_y2 = 30.0*sin(theta)
    vx = 30.0*cos(theta)
    h    = 0.01  ! step size
    t    = 0.0   ! initial t
    x    = 0.0
    n    = 2   ! number of ODEs

    ! initial conditions
    y(1) =  0.0      ! y(1) = u
    y(2) =  try_y2   ! y(2) = u' 

    do while(y(2) .Gt. -try_y2)
        if (y(2) .Le. -try_y2) then
            h = 0
        endif
        call rk4(n, y, ynext, t, h, my_func)
        y = ynext
        t = t + h
        x = x + h*vx
        print *, t, x, ynext(1),ynext(2)
    enddo
    print *, "x = ", x , "theta = " , theta
    print *, "The desired value is 50."

end program shooting1

!-----------------------------------------------
!
!  Solving Boundary Value problem
!
!  u'' = 6t      0 < t < 1
!
!  with BC
!
!         u(t=0) = 1 and
!         u(t=1) = 1
!
!-----------------------------------------------

subroutine my_func(n, t, yin, k)
    implicit none
    integer, intent(in) :: n  ! number of ODEs
    real, intent(in)    :: t
    real, dimension(n), intent(in)  :: yin
    real, dimension(n), intent(out) :: k    ! dydt

    ! in this example  n = 2
    k(1) =  yin(2) ! TODO
    k(2) =  -9.8  ! TODO
    return
end subroutine my_func







