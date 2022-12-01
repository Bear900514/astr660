module Solver 

    implicit none
    contains

    subroutine euler(n, yin, ynext, t, h, func)
        implicit none
        integer, intent(in) :: n      ! number of ODEs
        real, intent(in)    :: t, h
        real, dimension(n), intent(in)  :: yin
        real, dimension(n), intent(out)  :: ynext
        external      :: func
        integer            :: i
        real, dimension(n) :: k1

        ! call func to obtain the values of dydt
        call func(n,t,yin,k1)
        
        ynext(1) = yin(1) + k1(1)
        ynext(2) = yin(2) + k1(2)
        ynext(3) = yin(3) + k1(3)
        ynext(4) = yin(4) + k1(4)

        ! compute ynext using the Euler's method

    end subroutine euler

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
            ynext(i)=yin(i)+h*(k1(i)+2.0*k2(i)+2.0*k3(i)+k4(i))/6.0
         enddo
    end subroutine rk4
end module Solver
