module IO
    implicit none
    contains
        subroutine output(n,time)
            use Simulation_data
            implicit none
            integer, intent(in) :: n
            real, intent(in)    :: time
            real, parameter :: small = 1.0e-99

            integer      :: i
            character(7) :: cycstr
            character(58)      :: filename

            real  ,dimension(istart-ibuf:iend+ibuf)      :: ua ! analytical solution
            real         :: pa, pi, pe
            
            write(cycstr,10) n,'.d'
10  format(i5,a2)

            ! fill remaing zeros
            do i=1,5
                if(cycstr(i:i).eq.' ') cycstr(i:i)='0'
            enddo

            filename = 'advection_LW'//cycstr
            open(100,file=filename,status='unknown')

            ! write header
            write(100,29) "# i", "x", "U(x)", "UA(x)"

            pa = 0.15 + c * time
            pi = pa-0.05
            pe = pa+0.05

            do i=1,iend
                ua(i)=max(exp(-1000.0*(x(i)-0.1)**2.0),small)
                !if ((x(i) .ge. pi) .and. (x(i) .le. pe)) then
                    !ua = 1.0
                !else
                    !ua = 0.01
                !endif
                write(100,30) i, x(i), u(i), ua(i)
            enddo

29 format(a6,3a24)
30 format(i6,3e24.12)

            close(100)
            return
        end subroutine output

end module IO
