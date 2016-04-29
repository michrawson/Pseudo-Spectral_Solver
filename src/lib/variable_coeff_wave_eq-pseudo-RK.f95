module variable_coeff_wave_eq_pseudo_rk
use, intrinsic :: iso_c_binding
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run( nplots, plotgap, sigma1, sigma2, dt, delta2, &
                                                     h, n, x, y, tdata, result)
    use fftw_prime
    use max_finder
    implicit none
        integer, intent(in)          ::  nplots, plotgap
        real ( kind = 8 ), intent(in) :: sigma1, sigma2, dt, delta2, h
        integer, intent(in)          :: n
        real ( kind = 8 ), dimension(n), intent(out)                   :: x, y
        real ( kind = 8 ), dimension(nplots+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(nplots+1,n,n), intent(out)        :: result


        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: temp
        real ( kind = 8 ), dimension(n,n) :: v, k1,k2,k3,k4,k5,k6
        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: vx, vy
        real ( kind = 8 ), dimension(n,n) :: prime_x,prime_y

        integer               :: i, j
        real ( kind = 8 )      :: pi, t, curr_max, prev_max, start, finish

        pi = 4.*atan(1.)
        x = h*(/ (j-n/2,j=1,N) /)
        y = h*(/ (j-n/2,j=1,N) /)
        t = 0

!        v = 0
        do i = 1,n
            do j = 1,n
                v(i,j) = exp(-((x(i)**2.)/sigma1)-((y(j)**2.)/sigma2))
!                if (x(i) <= 15 .AND. x(i) >= -15) then
!                    v(i,j) = exp(-(y(j)**2.)*(-4.0*log(10.0**(-3)))/(1 + 0.025*cos(1.68*x(j)))**2.0)
!                end if
            end do
        end do

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

                temp = v
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_y_run(n, temp, prime_y)
                vx=0
                vy=0
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k1 = delta2*(-(vy * prime_x) + (vx * prime_y))

                temp = v + dt*k1/3.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_y_run(n, temp, prime_y)
!                temp = v + dt*k1/2.0
                vx=0
                vy=0
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k2 = delta2*(-(vy * prime_x) + (vx * prime_y))

                temp = v + dt*(4.0*k1+6.0*k2)/25.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_y_run(n, temp, prime_y)

!                temp = v + dt*k2/2.0
                vx=0
                vy=0
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k3 = delta2*(-(vy * prime_x) + (vx * prime_y))

                temp = v + dt*(k1-12.0*k2+15.0*k3)/4.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_y_run(n, temp, prime_y)

!                temp = v + dt*k3
                vx=0
                vy=0
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k4 = delta2*(-(vy * prime_x) + (vx * prime_y))

                temp = v + dt*(6.0*k1+90.0*k2-50.0*k3+8.0*k4)/81.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_y_run(n, temp, prime_y)
                vx=0
                vy=0
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k5 = delta2*(-(vy * prime_x) + (vx * prime_y))

                temp = v + dt*(6.0*k1+36.0*k2+10.0*k3+8.0*k4)/75.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_y_run(n, temp, prime_y)
                vx=0
                vy=0
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k6 = delta2*(-(vy * prime_x) + (vx * prime_y))

                v = v + dt/192.0*(23.0*k1 + 125.0*k3 - 81.0*k5 + 125.0*k6)
!                v = v + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
!                v = v + dt/6.0*(k1 + 4.0*k2 + k3)

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
