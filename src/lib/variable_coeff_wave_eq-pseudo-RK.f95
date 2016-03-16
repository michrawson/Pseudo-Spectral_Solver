module variable_coeff_wave_eq_pseudo_rk
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run(x, y, tdata, result)
    use fft_prime
    implicit none
        integer, parameter          :: n = 128
        integer, parameter          :: tmax = 8
        real ( kind = 8 ), parameter :: tplot = 0.15

        real ( kind = 8 ), dimension(n), intent(out)                   :: x, y
        real ( kind = 8 ), dimension(int(tmax/tplot)+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(int(tmax/tplot)+1,n,n), intent(out) :: result

        real ( kind = 8 ), dimension(n,n)              :: v, k1,k2,k3,k4,prime_x,prime_y,prime_phi_x,prime_phi_y

        integer               :: i, j, k, plotgap, nplots
        real ( kind = 8 )      :: pi, h, t, dt

        pi = 4.*atan(1.)
        h = 2.*pi/N
        x = h*(/ (j-n/2,j=1,N) /)
        y = h*(/ (j-n/2,j=1,N) /)
        t = 0
        dt = h/4.

        do i=1,n
            do j = 1,n
                v(i,j) = exp(-(((x(i)-pi/6)/0.1)**2.)-((y(j)/0.1)**2.)) + exp(-(((x(i)+pi/6)/0.1)**2.)-((y(j)/0.1)**2.))
            end do
        end do

        do i=1,n
!            prime_phi_x(:,i) = 1
            prime_phi_x(:,i) = 2. * x
!            prime_phi_y(i,:) = 1
            prime_phi_y(i,:) = 2. * y
        end do

        plotgap = nint(tplot/dt)
        dt = tplot/plotgap
        nplots = nint(tmax/tplot)

        tdata = -1
        tdata(1) = t

        result = 0
        result(1,1:n,1:n) = v

        do i=1,nplots
            do j = 1,plotgap
                t = t+dt

                call fft_prime_2d_partial_x_run(v, prime_x)
                call fft_prime_2d_partial_y_run(v, prime_y)

                k1 = prime_phi_y*prime_x - prime_phi_x*prime_y

!                k1 = 2. * spread(y, dim=1, ncopies=n) * prime_x - 2. * spread(x, dim=1, ncopies=n) * prime_y

                call fft_prime_2d_partial_x_run(v + dt*k1/2, prime_x)
                call fft_prime_2d_partial_y_run(v + dt*k1/2, prime_y)
                k2 = prime_phi_y*prime_x - prime_phi_x*prime_y

!                k2 = 2. * spread(y, dim=1, ncopies=n) * prime_x - 2. * spread(x, dim=1, ncopies=n) * prime_y

                call fft_prime_2d_partial_x_run(v + dt*k2/2, prime_x)
                call fft_prime_2d_partial_y_run(v + dt*k2/2, prime_y)
                k3 = prime_phi_y*prime_x - prime_phi_x*prime_y

!                k3 = 2. * spread(y, dim=1, ncopies=n) * prime_x - 2. * spread(x, dim=1, ncopies=n) * prime_y

                call fft_prime_2d_partial_x_run(v + dt*k3, prime_x)
                call fft_prime_2d_partial_y_run(v + dt*k3, prime_y)
                k4 = prime_phi_y*prime_x - prime_phi_x*prime_y

!                k4 = 2. * spread(y, dim=1, ncopies=n) * prime_x - 2. * spread(x, dim=1, ncopies=n) * prime_y

                v = v + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)

            end do
            result(i+1,1:n,1:n) = v
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_rk_run

end module variable_coeff_wave_eq_pseudo_rk
