
program main
    use variable_coeff_wave_eq_pseudo_rk

    implicit none
    integer, parameter           :: n = 16
    integer, parameter           :: nplots = 1

    real ( kind = 8 ), dimension(nplots+1)      :: tdata
    real ( kind = 8 ), dimension(1:nplots+1,1:n,1:n)  :: result

    real ( kind = 8 )                 :: x(n), y(n)

    integer           :: plotgap = 100
    real ( kind = 8 ) :: sigma1 = 2.52
    real ( kind = 8 ) :: sigma2 = 13.4
    real ( kind = 8 ) :: dt, pi, delta2, h

    pi = 4.*atan(1.)
    dt = 1./2.**6.
    delta2 = 0.8
    h = 44.0/n

    call variable_coeff_wave_eq_pseudo_rk_run(nplots, plotgap, &
        sigma1, sigma2, dt, delta2,h,n,x,y,tdata,result)

    PRINT *,"x",x(1)
    PRINT *,"tdata",tdata(1)
    PRINT *,"result",result(1,1,1)

end program main
