
program main
    use variable_coeff_wave_eq_pseudo_rk

    implicit none
    integer, parameter           :: n = 42
    integer, parameter           :: nplots = 2

    real ( kind = 8 ), dimension(nplots+1)      :: tdata
    real ( kind = 8 ), dimension(nplots+1,n,n)  :: result

    real ( kind = 8 )                 :: x(n), y(n)

    integer :: plotgap = 200
    real ( kind = 8 ) :: sigma1 = .1
    real ( kind = 8 ) :: sigma2 = 1.0
    real ( kind = 8 ) :: dt, pi, delta2

    pi = 4.*atan(1.)
    dt = 2.0*pi/n/4./25.
    delta2 = 4.0

    call variable_coeff_wave_eq_pseudo_rk_run(nplots, plotgap, &
            sigma1, sigma2, dt, delta2,n,x,y,tdata,result)

    PRINT *,"x",x(1:10)
    PRINT *,"tdata",tdata(1:10)
    PRINT *,"result",result(1,1:10,1:10)

end program main
