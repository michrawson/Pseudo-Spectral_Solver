
program main
    use variable_coeff_wave_eq_pseudo_rk
    use zero_finder

    implicit none

    integer, parameter           :: n = 128
    integer, parameter           :: tmax = 8
    integer                      :: j
    real ( kind = 8 ), parameter :: tplot = 0.15

    real ( kind = 8 ), dimension(N)                      :: x, y
    real ( kind = 8 ), dimension(int(tmax/tplot)+1)      :: tdata
    real ( kind = 8 ), dimension(int(tmax/tplot)+1,n,n)    :: result

    real ( kind = 8 )                 :: p_xi=0
    real ( kind = 8 )                 :: p_x(n)=0
    real ( kind = 8 )                 :: p_y(n)=0
    real ( kind = 8 )                   :: pi

    pi = 4.*atan(1.)

    p_xi = pi

    p_x = [ (j, j = 1,n) ]
    p_x = p_x/n*2*pi
    p_y = cos(p_x)
    p_xi = triginterp(p_xi,p_x,p_y)
    print *, "p_xi", p_xi


!    call variable_coeff_wave_eq_pseudo_rk_run(x,y,tdata,result)
!    PRINT *,"x",x
!    PRINT *,"tdata",tdata
!    PRINT *,"result",result(1,:,1)

end program main
