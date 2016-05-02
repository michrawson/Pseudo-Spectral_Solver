module solve
    use variable_coeff_wave_eq_pseudo_rk
    implicit none

contains

    subroutine solve_run( nplots, plotgap, sigma1, sigma2, dt, delta2, &
        h, n, x, y, tdata, result,int_res, int2_res, max_res, min_res)
        implicit none
        integer, intent(in)          ::  nplots, plotgap
        real ( kind = 8 ), intent(in) :: sigma1, sigma2, dt, delta2, h
        integer, intent(in)          :: n
        real ( kind = 8 ), dimension(n), intent(out)                   :: x, y
        real ( kind = 8 ), dimension(nplots+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(nplots+1,n,n), intent(out)        :: result
        real ( kind = 8 ), dimension(1:nplots), intent(out) :: int_res, int2_res, max_res, min_res

        call variable_coeff_wave_eq_pseudo_rk_run( nplots, plotgap, sigma1, sigma2, dt, &
            delta2, h, n, x, y, tdata, result,int_res, int2_res, max_res, min_res)

    end subroutine solve_run

end module solve
