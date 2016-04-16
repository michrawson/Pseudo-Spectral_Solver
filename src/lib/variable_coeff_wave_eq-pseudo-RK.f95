module variable_coeff_wave_eq_pseudo_rk
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run(x, y, tdata, result)
    use fftw_prime
    use max_finder
    implicit none
        integer, parameter          :: n = 128
        integer, parameter          :: tmax = 4
        integer, parameter          :: nplots = 10
        real ( kind = 8 ), parameter :: tplot = 0.2

        real ( kind = 8 ), dimension(n), intent(out)                   :: x, y
        real ( kind = 8 ), dimension(int(tmax/tplot)+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(10+1,n,n), intent(out) :: result

        real ( kind = 8 ), dimension(n,n) :: v, k1,k2,k3,k4
        real ( kind = 8 ), dimension(n,n) :: vx, vy
        real ( kind = 8 ), dimension(n,n) :: prime_x,prime_y

        integer               :: i, j, plotgap
!        , nplots
        real ( kind = 8 )      :: pi, h, t, dt, curr_max, prev_max, delta2

        pi = 4.*atan(1.)
        h = 2.*pi/N
        x = h*(/ (j-n/2,j=1,N) /)
        y = h*(/ (j-n/2,j=1,N) /)
        t = 0
        dt = h/4./20.
        delta2 = 4.0

        do i=1,n
            do j = 1,n
                v(i,j) = exp(-(((x(i)-pi/6)/0.1)**2.)-((y(j)/1.0)**2.))
!                v(i,j) = exp(-(((x(i)-pi/6)/0.1)**2.)-((y(j)/0.1)**2.)) + exp(-(((x(i)+pi/6)/0.1)**2.)-((y(j)/0.1)**2.))
            end do
        end do

!        nplots = 20
        plotgap = nint(tplot/dt)

        tdata = -1
        tdata(1) = t

        result = 0
        result(1,1:n,1:n) = v

!        prev_max = max_finder_run(maxloc(v), v)

        do i=1, nplots
            print *,"nplot #",i
            do j = 1,plotgap
                print *,"plotgap #",j
                t = t+dt

                call fft_prime_2d_partial_x_run(v, prime_x)
                call fft_prime_2d_partial_y_run(v, prime_y)
                call poisson2df(n,n,h,h,v,vx,vy,1)

                k1 = delta2*(-1.0 * (vy * prime_x) - (vx * prime_y))

                call fft_prime_2d_partial_x_run(v + dt*k1/2, prime_x)
                call fft_prime_2d_partial_y_run(v + dt*k1/2, prime_y)
                call poisson2df(n,n,h,h,v,vx,vy,1)

                k2 = delta2*(-1.0 * (vy * prime_x) - (vx * prime_y))

                call fft_prime_2d_partial_x_run(v + dt*k2/2, prime_x)
                call fft_prime_2d_partial_y_run(v + dt*k2/2, prime_y)
                call poisson2df(n,n,h,h,v,vx,vy,1)

                k3 = delta2*(-1.0 * (vy * prime_x) - (vx * prime_y))

                call fft_prime_2d_partial_x_run(v + dt*k3, prime_x)
                call fft_prime_2d_partial_y_run(v + dt*k3, prime_y)
                call poisson2df(n,n,h,h,v,vx,vy,1)

                k4 = delta2*(-1.0 * (vy * prime_x) - (vx * prime_y))

                v = v + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)

                if (maxval(-result(i,1:n,1:n)) < 0) then
                    PRINT *,"error: minval: ", maxval(-result(i,1:n,1:n))
                end if

!                PRINT *,"error: volume: diff", (abs( SUM(MATMUL(v,(/ (1,j=1,N) /))) &
!                                             - SUM(MATMUL(result(1,1:n,1:n),(/ (1,j=1,N) /))) ))
!
!                PRINT *,"error: maxval: diff", ( maxval(result(i,1:n,1:n)) - maxval(v) )

            end do

!            curr_max = max_finder_run(maxloc(v), v)
!            PRINT *,"error: better maxval: diff", ( prev_max - curr_max )
!
!            prev_max = curr_max

            result(i+1,1:n,1:n) = v
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_rk_run

end module variable_coeff_wave_eq_pseudo_rk
