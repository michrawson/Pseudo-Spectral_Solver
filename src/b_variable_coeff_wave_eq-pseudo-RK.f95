module variable_coeff_wave_eq_pseudo_rk
use, intrinsic :: iso_c_binding
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run( nplots, plotgap, sigma1, sigma2, dt, delta2, &
                                                     h, n, x, y, tdata, result)
    use fftw_prime
    use max_finder
    use int_finder
    implicit none
        integer, intent(in)          ::  nplots, plotgap
        real ( kind = 8 ), intent(in) :: sigma1, sigma2, dt, delta2, h
        integer, intent(in)          :: n
        real ( kind = 8 ), dimension(n), intent(out)                   :: x, y
        real ( kind = 8 ), dimension(1:nplots+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(1:nplots+1,1:n,1:n), intent(out)        :: result

        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: temp
        real ( kind = 8 ), dimension(n,n) :: v, k1,k2,k3,k4,k5,k6
        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: vx, vy
        real ( kind = 8 ), dimension(n,n) :: prime_x,prime_y

        integer                :: i, j, k
        real ( kind = 8 )      :: pi, t, curr_max, prev_max, start, finish
        real ( kind = 8 )      :: prev_int_v,prev_int_v2,curr_int_v,curr_int_v2

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

        call int_finder_run(n, h, v, prev_int_v)

        call int_finder_run(n, h, v*v, prev_int_v2)

        tdata(1) = t

        result(1,:,:) = v(:,:)

!        prev_max = max_finder_run(n, maxloc(v), v)

        do i=1, nplots
!            print *,"nplot #",i
            do j = 1,plotgap
!                print *,"plotgap #",j

                t = t+dt

                temp = v
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_x_run(n, transpose(temp), prime_y)
                prime_y = transpose(prime_y)
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k1 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*k1/3.0
                temp = v + dt*k1/2.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_x_run(n, transpose(temp), prime_y)
                prime_y = transpose(prime_y)
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k2 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*(4.0*k1+6.0*k2)/25.0
                temp = v + dt*k2/2.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_x_run(n, transpose(temp), prime_y)
                prime_y = transpose(prime_y)

                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k3 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*(k1-12.0*k2+15.0*k3)/4.0
                temp = v + dt*k3
                call fft_prime_2d_partial_x_run(n, temp, prime_x)
                call fft_prime_2d_partial_x_run(n, transpose(temp), prime_y)
                prime_y = transpose(prime_y)
                call poisson2df(n,n,h,h,temp,vx,vy,1)

                k4 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*(6.0*k1+90.0*k2-50.0*k3+8.0*k4)/81.0
!                call fft_prime_2d_partial_x_run(n, temp, prime_x)
!                call fft_prime_2d_partial_x_run(n, transpose(temp), prime_y)
!                prime_y = transpose(prime_y)
!                call poisson2df(n,n,h,h,temp,vx,vy,1)
!
!                k5 = delta2*(-(vy * prime_x) + (vx * prime_y))
!
!                temp = v + dt*(6.0*k1+36.0*k2+10.0*k3+8.0*k4)/75.0
!                call fft_prime_2d_partial_x_run(n, temp, prime_x)
!                call fft_prime_2d_partial_x_run(n, transpose(temp), prime_y)
!                prime_y = transpose(prime_y)
!                call poisson2df(n,n,h,h,temp,vx,vy,1)
!
!                k6 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                v = v + dt/192.0*(23.0*k1 + 125.0*k3 - 81.0*k5 + 125.0*k6)
                v = v + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
!                v = v + dt/6.0*(k1 + 4.0*k2 + k3)

                if (maxval(-result(i,1:n,1:n)) < 0) then
                    PRINT *,"error: minval: ", maxval(-result(i,1:n,1:n))
                end if

                call int_finder_run(n, h, v, curr_int_v)
                PRINT *,"error: integral n^1: diff", prev_int_v-curr_int_v
                prev_int_v = curr_int_v

                call int_finder_run(n, h, v*v, curr_int_v2)
                PRINT *,"error: integral n^2: diff", prev_int_v2-curr_int_v2
                prev_int_v2 = curr_int_v2

!                PRINT *,"error: volume: diff", (abs( SUM(MATMUL(v,(/ (1,j=1,N) /))) &
!                                             - SUM(MATMUL(result(1,1:n,1:n),(/ (1,j=1,N) /))) ))
!
!                PRINT *,"error: maxval: diff", ( maxval(result(i,1:n,1:n)) - maxval(v) )

            end do

!            curr_max = max_finder_run(n, maxloc(v), v)
!            PRINT *,"error: maxval: diff", ( prev_max - curr_max )
!            prev_max = curr_max

            result(i+1,:,:) = v(:,:)
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_rk_run

end module variable_coeff_wave_eq_pseudo_rk
