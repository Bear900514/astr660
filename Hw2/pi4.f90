program pi

    implicit none
    real*8 :: dA ,A
    real*8 :: error
    integer :: i

    real*8,parameter:: p=4.0*atan(1.0)
    integer,parameter :: Nmax=7
    integer,dimension(Nmax) :: N
    N=(/10,100,1000,10000,&
        100000,1000000,10000000/)
    
    
    open(unit=8,file="pi_error4-41.dat")
    write(8, *)
    !8 Format('N,error,A')
    
    do i =1,Nmax 
      call calculate_function_integral(N(i), A)
      error=(p-2*A)/(p)
      write(8, *) N(i),error,2*A
      !8 Format(I16,F12.9,F10.7)
      !print *, N(i),error
    enddo
      
    !-- print out the result
    !print *, "PI = ", 2.*area
    
    CLOSE(unit=8)

end program pi

subroutine calculate_function_integral(N, A)
    implicit none
    integer, intent(in)  :: N

    real*8, intent(out) :: A ! area of the function

    integer :: i
    real    :: my_func
    real*8    :: x, h, dx, dA,area

    dx = 2./N
    x=2.
    area  = 0.
    do i = 1, N-1
        x = -1.+ dx*(i-1)
        !h = my_func(x+(dx/2.))
        !h = (my_func(x)+my_func(x+dx))/2.
        h = (my_func(x)+4*my_func(x+(dx/2.))+my_func(x+dx))/6.
        dA = h * dx
        area  = area + dA
    enddo
    A=area
    return
end subroutine calculate_function_integral

real function my_func(x)
    ! the function return the y values of a half circle with radius = 1
    real*8 :: x
    my_func = sqrt(1.0 - x**2.)
    return
end function
