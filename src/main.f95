
program main
    use finite_differences
!    use variable_coeff_wave_eq
    use variable_coeff_wave_eq_pseudo
    implicit none

    double precision                                        :: error(9+1) = 0

    integer, parameter                          :: N = 128, tmax = 8
    double precision, parameter                         :: tplot = 0.15
    double precision, dimension(N)                      :: x = 0
    double precision, dimension(int(tmax/tplot)+1,N)    :: result
    double precision, dimension(int(tmax/tplot))      :: tdata

!    call finite_differences_run(error)
!    PRINT *,"error",error(1:10)

    call variable_coeff_wave_eq_pseudo_run(x,tdata,result)
    PRINT *,"x",x
    PRINT *,"tdata",tdata
    PRINT *,"result",result(1,:)
    PRINT *,"result",result(2,:)
!    PRINT *,"result",result(:,10)
!    PRINT *,"result",result(:,100)
!    PRINT *,"result",result(:,N)

end program main
