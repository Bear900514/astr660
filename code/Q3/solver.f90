module Solver 

    implicit none
    contains
    subroutine bisection(func, xs, err)
        implicit none
        real, external    :: func    ! the function to solve
        real, intent(out) :: xs      ! solution
        real, intent(out) :: err     ! error
        real, save :: a = -5.0        ! bracking interval [a,b]
        real, save :: b = 0.0        ! bracking interval [a,b]
        real  :: fa, fx              ! f(a) and f(x)

        xs= (a+b)/2
        fa=func(a)
        fx=func(xs)
        err= abs(func(xs))
        if(sign(1.,fx)/=sign(1.,fa)) then
            b=xs
        else
            a=xs
        endif
        return
        
        end subroutine bisection

    subroutine rk2(n, yin, ynext, t, h, func)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: t, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: func
        integer            :: i
        real, dimension(n) :: k1, k2
        real,dimension(n)  :: y2

        ! compute k1 = func(t, yin)
        call func(n,t,yin,k1)
        ! compute y2 = yin + h*k1
        do i=1,n
            y2(i)=yin(i) +h*k1(i)
        enddo
        ! compute k2 = func(t+h, y2)
        call func(n,t+h,y2,k2)
        ! compute ynext 
         do i=1,n
            ynext(i)=yin(i)+h*(k1(i)+k2(i))/2.0
         enddo
    end subroutine rk2
    subroutine rk4(n, yin, ynext, t, h, func)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: t, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: func
        integer :: i
        real              :: h2
        real,dimension(n) :: k1,k2,k3,k4
        real,dimension(n) :: y2,y3,y4

        ! compute k1 = func(t, yin)
        call func(n,t,yin,k1)
        ! compute y2 = yin + h*k1/2
        do i=1,n
            y2(i)=yin(i) +h*k1(i)/2.0
        enddo
        ! compute k2 = func(t+h/2, y2)
        call func(n,t+h/2.0,y2,k2)
         ! compute y3 = yin + h*k2/2
        do i=1,n
            y3(i)=yin(i) +h*k2(i)/2.0
        enddo
        ! compute k3 = func(t+h/2, y3)
        call func(n,t+h/2.0,y3,k3)
         ! compute y4 = yin + h*k3
        do i=1,n
            y4(i)=yin(i) +h*k3(i)
        enddo
        ! compute k4 = func(t+h, y4)
        call func(n,t+h,y4,k4)
        ! compute ynext 
         do i=1,n
            ynext(i)=yin(i)+h*(k1(i)+2*k2(i)+2*k3(i)+k4(i))/6.0
         enddo
    end subroutine rk4

end module Solver
