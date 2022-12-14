subroutine boundary_x(v)
    !
    ! apply BC on array v
    !
    use Simulation_data
    implicit none
    real, dimension(istart-ibuf:iend+ibuf,istart-ibuf:iend+ibuf), intent(inout) :: v
    integer :: i,j


    ! apply boundary condition
    do j = 1 ,imax
    ! left boundary (period)
        do i = 1, ibuf
            v(istart-i,j) = v(iend-i+1,j)
        enddo
        ! right boundary (period)
        do i = 1, ibuf
            v(iend+i,j) = v(istart+i-1,j)
        enddo
    enddo

end subroutine boundary_x

subroutine boundary_y(v)
    !
    ! apply BC on array v
    !
    use Simulation_data
    implicit none
    real, dimension(istart-ibuf:iend+ibuf,istart-ibuf:iend+ibuf), intent(inout) :: v
    integer :: i,j


    ! apply boundary condition
    do j = 1 ,ibuf
    ! left boundary (period)
        do i = 1,imax
            v(i,istart-j) = v(i,iend-j+1)
        enddo
        ! right boundary (period)
        do i = 1,imax
            v(i,iend+j) = v(i,istart+j-1)
        enddo
    enddo

end subroutine boundary_y
