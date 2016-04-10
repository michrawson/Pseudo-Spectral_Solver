
program main
    use variable_coeff_wave_eq_pseudo_rk
    use zero_finder

    implicit none

    integer, parameter           :: n = 128
    integer, parameter           :: tmax = 8
    integer                      :: j, k
    real ( kind = 8 ), parameter :: tplot = 0.15

    real ( kind = 8 ), dimension(int(tmax/tplot)+1)      :: tdata
!    real ( kind = 8 ), dimension(int(tmax/tplot)+1,n,n)  :: result

    real ( kind = 8 )                 :: result=0
    real ( kind = 8 )                 :: p_x(n)=0
    real ( kind = 8 )                 :: p_y(n,n)=0
    real ( kind = 8 )                 :: pi

    pi = 4.*atan(1.)

    p_x = [ (j, j = 1,n) ]
    p_x = p_x*2.0/n*pi
    do k = 1, n
        do j = 1, n
            p_y(k,j) = cos(p_x(k)) + cos(p_x(j))
        end do
    end do

    result = triginterp([pi-sqrt(2.0), pi-sqrt(2.0)],p_x,p_y)
    print *, "result", result

!    call variable_coeff_wave_eq_pseudo_rk_run(x,y,tdata,result)
!    PRINT *,"x",x
!    PRINT *,"tdata",tdata
!    PRINT *,"result",result(1,:,1)

end program main
