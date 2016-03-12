module variable_coeff_wave_eq_pseudo_rk
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run(x,tdata,result)
    use fft_prime
    implicit none
        integer, parameter          :: n = 128
        integer, parameter          :: tmax = 8
        real ( kind = 8 ), parameter :: tplot = 0.15

        real ( kind = 8 ), dimension(n), intent(out)                   :: x
        real ( kind = 8 ), dimension(int(tmax/tplot)+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(int(tmax/tplot)+1,n), intent(out) :: result

        real ( kind = 8 ), dimension(n)              :: c, v, k1,k2,k3,k4,prime

        integer               :: i, j, plotgap, nplots
        real ( kind = 8 )      :: pi, h, t, dt

        pi = 4.*atan(1.)
        h = 2.*pi/N
        x = h*(/ (j,j=1,N) /)
        t = 0
        dt = h/4.
        c = .2 + sin(x-1.)**2.
        v = exp(-100*((x-1)**2.))

        plotgap = nint(tplot/dt)
        dt = tplot/plotgap
        nplots = nint(tmax/tplot)

        tdata = -1
        tdata(1) = t

        result = 0
        result(1,1:N) = v

        do i=1,nplots
            do j = 1,plotgap
                t = t+dt

                call fft_prime_run(v, prime)
                k1 = -1.0*prime

                call fft_prime_run(v + dt*k1/2, prime)
                k2 = -1.0*prime

                call fft_prime_run(v + dt*k2/2, prime)
                k3 = -1.0*prime

                call fft_prime_run(v + dt*k3, prime)
                k4 = -1.0*prime

                v = v + dt/6.0*(k1 + 2.*k2 + 2.*k3 + k4)

            end do
            result(i+1,:) = v
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_rk_run

end module variable_coeff_wave_eq_pseudo_rk
