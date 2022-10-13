module solver

    ! ---------------------------------------------------------
    ! ---------------------------------------------------------

    implicit none
    contains 

        !--------------------------------------------------
        !
        ! Implment the bisection method to solve the func
        !
        !
        ! Inputs:   func
        !
        ! Outputs:
        !           xs   : solution of func
        !           error: relative error
        !
        !--------------------------------------------------


        subroutine bisection(func, xs, err)
        implicit none
        real, external    :: func    ! the function to solve
        real, intent(out) :: xs      ! solution
        real, intent(out) :: err     ! error
        real, save :: a = 1.0        ! bracking interval [a,b]
        real, save :: b = 3.0        ! bracking interval [a,b]
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


        !--------------------------------------------------
        !
        ! Implment the Newton's method to solve my_func
        !
        !
        ! Inputs: func, dfunc
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine newton(func, dfunc, xs, err)
        implicit none
        real, external    :: func  ! the function to solve
        real, external    :: dfunc ! the first derivative of the function to solve
        real, intent(out) :: xs    ! solution
        real, intent(out) :: err   ! error
        real, save :: x = 2.0      ! trial value
        real :: fx
        real :: dfdx

        dfdx = dfunc(x)
        fx= func(x)   
        xs = x-fx/dfdx
        err=abs(func(xs))
        x=xs
        

        end subroutine newton


        !--------------------------------------------------
        !
        ! Implment the Secant method to solve my_func
        !
        !
        ! Inputs: None
        !
        ! Outputs:
        !           xs   : solution of my_func
        !           error: relative error
        !
        !--------------------------------------------------
        subroutine secant(func, xs, err)
        implicit none
        real, external    :: func  ! the function to solve
        real, intent(out) :: xs    ! solution
        real, intent(out) :: err   ! error

        real,save :: x0 = 1.0  ! initial guess
        real,save :: x1 = 3.0  ! initial guess

        real :: fx0, fx1       ! f(x0) and f(x1)

        fx0= func(x0)
        fx1= func(x1)
        xs = x0-fx0/((fx0-fx1)/(x0-x1))
        err=abs(func(xs))
        x1=x0
        x0=xs
        end subroutine secant


end module solver
