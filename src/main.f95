
program main
    use finite_differences
    use variable_coeff_wave_eq
    implicit none

    double precision                                        :: error(9+1) = 0

    integer, parameter                          :: N = 128, tmax = 8
    double precision, parameter                         :: tplot = 0.15
    double precision, dimension(N)                      :: x = 0
    double precision, dimension(N,N)    :: result
    double precision, dimension(int(tmax/tplot))      :: tdata

    call finite_differences_run(error)
    PRINT *,"error",error(1:10)

    call variable_coeff_wave_eq_run(x,tdata,result)
    PRINT *,"x",x(1:10)
    PRINT *,"tdata",tdata(1:10)
    PRINT *,"result",result(1:10,1)

end program main
