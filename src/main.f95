
program main
    use finite_differences
    use variable_coeff_wave_eq
    implicit none

    real                    :: error(9+1) = 0

    integer, parameter      :: N = 128
    real, dimension(N)                  :: x = 0
    real, dimension(:,:), allocatable   :: result
    real, dimension(:), allocatable     :: tdata

    call finite_differences_run(error)
    PRINT *,"error",error(1:10)

    call variable_coeff_wave_eq_run(x,tdata,result)
    PRINT *,"x",x(1:10)
    PRINT *,"tdata",tdata(1:10)
    PRINT *,"result",result(1:10,1)

end program main
