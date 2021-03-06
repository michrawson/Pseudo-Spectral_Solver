module variable_coeff_wave_eq_pseudo_rk
    use, intrinsic :: iso_c_binding
    use fftw_prime
    use max_finder
    use int_finder
    implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run( nplots, plotgap, sigma1, sigma2, dt, delta2, &
        h, n, x, y, tdata, result, int_res, int2_res, max_res, min_res)
        implicit none
        integer, intent(in)          ::  nplots, plotgap
        real ( kind = 8 ), intent(in) :: sigma1, sigma2, dt, delta2, h
        integer, intent(in)          :: n
        real ( kind = 8 ), dimension(n), intent(out)                   :: x, y
        real ( kind = 8 ), dimension(1:nplots+1), intent(out)   :: tdata
        real ( kind = 8 ), dimension(1:nplots+1,1:n,1:n), intent(out)   :: result
        real ( kind = 8 ), dimension(1:nplots+1), intent(out) :: int_res, int2_res, max_res, min_res

        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: temp, temp_t
        real ( kind = 8 ), dimension(n,n) :: v, v_temp, k1,k2,k3,k4,k5,k6
        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: vx, vy
        complex (C_DOUBLE_COMPLEX), dimension(3*n,3*n) :: f2, fhat2, fhatx, fhaty
        real ( kind = 8 ), dimension(n,n) :: prime_x,prime_y

        integer                :: i, j, v_maxloc(2)
        integer ( kind = 8 )   :: poisson2df_fftw_plan_forward,poisson2df_fftw_plan_back
        type(C_PTR)            :: fft_prime_plan_forward, fft_prime_plan_back
        real ( kind = 8 )      :: pi, t, step_size

        integer       :: is_converged, iter_count

        poisson2df_fftw_plan_forward=0
        poisson2df_fftw_plan_back=0
        fft_prime_plan_forward=c_null_ptr
        fft_prime_plan_back=c_null_ptr

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

        tdata(1) = t

        result(:,:,:) = 0

        result(1,:,:) = v(:,:)

        int_res=0
        int2_res=0
        max_res=0
        min_res=0

        max_res(1) = -min_finder_approx(n, -v, is_converged)
        min_res(1) = minval(v)
        int2_res(1) = int_finder_run(n, h, v*v)
        int_res(1) = int_finder_run(n, h, v)

        do i=1, nplots
            print *,"nplot #",i
            do j = 1,plotgap
!                print *,"plotgap #",j

                t = t+dt

                temp = v
                call fft_prime_2d_partial_x_run(n, temp, prime_x, fft_prime_plan_forward, fft_prime_plan_back)
                temp_t = transpose(temp)
                call fft_prime_2d_partial_x_run(n, temp_t, prime_y, fft_prime_plan_forward, fft_prime_plan_back)
                prime_y = transpose(prime_y)

                call poisson2df(n,n,h,h,temp,vx,vy,1,poisson2df_fftw_plan_forward, &
                    poisson2df_fftw_plan_back,f2, fhat2, fhatx, fhaty)

                k1 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*k1/3.0
                temp = v + dt*k1/2.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x, fft_prime_plan_forward, fft_prime_plan_back)
                temp_t = transpose(temp)
                call fft_prime_2d_partial_x_run(n, temp_t, prime_y, fft_prime_plan_forward, fft_prime_plan_back)
                prime_y = transpose(prime_y)

                call poisson2df(n,n,h,h,temp,vx,vy,1,poisson2df_fftw_plan_forward, &
                                poisson2df_fftw_plan_back,f2, fhat2, fhatx, fhaty)

                k2 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*(4.0*k1+6.0*k2)/25.0
                temp = v + dt*k2/2.0
                call fft_prime_2d_partial_x_run(n, temp, prime_x, fft_prime_plan_forward, fft_prime_plan_back)
                temp_t = transpose(temp)
                call fft_prime_2d_partial_x_run(n, temp_t, prime_y, fft_prime_plan_forward, fft_prime_plan_back)
                prime_y = transpose(prime_y)

                call poisson2df(n,n,h,h,temp,vx,vy,1,poisson2df_fftw_plan_forward, &
                                poisson2df_fftw_plan_back,f2, fhat2, fhatx, fhaty)

                k3 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*(k1-12.0*k2+15.0*k3)/4.0
                temp = v + dt*k3
                call fft_prime_2d_partial_x_run(n, temp, prime_x, fft_prime_plan_forward, fft_prime_plan_back)
                temp_t = transpose(temp)
                call fft_prime_2d_partial_x_run(n, temp_t, prime_y, fft_prime_plan_forward, fft_prime_plan_back)
                prime_y = transpose(prime_y)

                call poisson2df(n,n,h,h,temp,vx,vy,1,poisson2df_fftw_plan_forward, &
                                poisson2df_fftw_plan_back,f2, fhat2, fhatx, fhaty)

                k4 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                temp = v + dt*(6.0*k1+90.0*k2-50.0*k3+8.0*k4)/81.0
!                call fft_prime_2d_partial_x_run(n, temp, prime_x, fft_prime_plan_forward, fft_prime_plan_back)
!                temp_t = transpose(temp)
!                call fft_prime_2d_partial_x_run(n, temp_t, prime_y, fft_prime_plan_forward, fft_prime_plan_back)
!                prime_y = transpose(prime_y)
!                call poisson2df(n,n,h,h,temp,vx,vy,1,poisson2df_fftw_plan_forward, &
!                                poisson2df_fftw_plan_back,f2, fhat2, fhatx, fhaty)
!
!                k5 = delta2*(-(vy * prime_x) + (vx * prime_y))
!
!                temp = v + dt*(6.0*k1+36.0*k2+10.0*k3+8.0*k4)/75.0
!                call fft_prime_2d_partial_x_run(n, temp, prime_x, fft_prime_plan_forward, fft_prime_plan_back)
!                temp_t = transpose(temp)
!                call fft_prime_2d_partial_x_run(n, temp_t, prime_y, fft_prime_plan_forward, fft_prime_plan_back)
!                prime_y = transpose(prime_y)
!                call poisson2df(n,n,h,h,temp,vx,vy,1,poisson2df_fftw_plan_forward, &
!                                poisson2df_fftw_plan_back,f2, fhat2, fhatx, fhaty)
!
!                k6 = delta2*(-(vy * prime_x) + (vx * prime_y))

!                v = v + dt/192.0*(23.0*k1 + 125.0*k3 - 81.0*k5 + 125.0*k6)
                v = v + dt/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
!                v = v + dt/6.0*(k1 + 4.0*k2 + k3)

            end do


            print *,"maxval ",maxval(v)

            max_res(i+1) = -min_finder_approx(n, -v, is_converged)
            min_res(i+1) = minval(v)
            int2_res(i+1) = int_finder_run(n, h, v*v)
            int_res(i+1) = int_finder_run(n, h, v)

            result(i+1,:,:) = v(:,:)
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_rk_run

end module variable_coeff_wave_eq_pseudo_rk
