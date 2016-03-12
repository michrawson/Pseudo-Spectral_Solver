
program main
    use variable_coeff_wave_eq_pseudo_rk
    implicit none

    integer, parameter          :: n = 128
    integer, parameter          :: tmax = 8
    real ( kind = 8 ), parameter :: tplot = 0.15

    real ( kind = 8 ), dimension(N)                      :: x
    real ( kind = 8 ), dimension(int(tmax/tplot)+1)      :: tdata
    real ( kind = 8 ), dimension(int(tmax/tplot)+1,N)    :: result

    call variable_coeff_wave_eq_pseudo_rk_run(x,tdata,result)
    PRINT *,"x",x
    PRINT *,"tdata",tdata
    PRINT *,"result",result(1,:)

end program main
