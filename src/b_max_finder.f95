module max_finder
    use dft_mod
    implicit none


contains

    real ( kind = 8 ) function max_finder_run_caller(p_x1, p_x2)
        implicit none

        integer, parameter           :: n = 128
        integer, intent(in)                :: p_x1, p_x2

        real ( kind = 8 )                 :: p_x(n)
        real ( kind = 8 )                 :: p_v(n,n)
        real ( kind = 8 )                   :: pi
        integer                             :: j,k, is_converged

        pi = 4.*atan(1.)

        p_x = [ (j, j = 1,n) ]
        p_x = p_x*2.0/n*pi

        do k = 1, n
            do j = 1, n
                p_v(k,j) = -cos(p_x(k)) - cos(p_x(j))
            end do
        end do

        max_finder_run_caller = min_finder_approx(n, -p_v, is_converged)

    end function max_finder_run_caller

    real ( kind = 8 ) function min_finder_approx(n, p_v, is_converged)
        implicit none
        integer, intent(in)                 :: n
        real ( kind = 8 ), intent(in)       :: p_v(n,n)
        integer, intent(inout)              :: is_converged

        real ( kind = 8 )                   :: step_size, start, finish, minimum
        complex ( kind = 8 ), dimension(n, n)  ::  fft_2d_prime_x_of_v, fft_2d_prime_y_of_v, &
            fft_2d_prime_x_x_of_v, fft_2d_prime_y_y_of_v, fft_2d_prime_x_y_of_v, fft_2d_prime_y_x_of_v
        integer                             :: step_size_iter_count, x0_index(2)

        call cpu_time(start)

        x0_index = minloc(p_v)

        fft_2d_prime_x_of_v = p_v
        call fft_prime_2d_partial_x(n, fft_2d_prime_x_of_v)

        fft_2d_prime_y_of_v = p_v
        call fft_prime_2d_partial_y(n, fft_2d_prime_y_of_v)


        fft_2d_prime_x_x_of_v = p_v
        call fft_prime_prime_2d_partial_x_x(n, fft_2d_prime_x_x_of_v)

        fft_2d_prime_y_y_of_v = p_v
        call fft_prime_prime_2d_partial_y_y(n, fft_2d_prime_y_y_of_v)

        fft_2d_prime_x_y_of_v = p_v
        call fft_prime_prime_2d_partial_x_y(n, fft_2d_prime_x_y_of_v)

        fft_2d_prime_y_x_of_v = fft_2d_prime_x_y_of_v
!        fft_2d_prime_y_x_of_v = p_v
!        call fft_prime_prime_2d_partial_y_x(n, fft_2d_prime_y_x_of_v)

        step_size = 1.0
        do step_size_iter_count = 1, 5

            minimum = min_finder_Newton(n, x0_index, p_v, step_size, is_converged, &
                fft_2d_prime_x_of_v, fft_2d_prime_y_of_v, fft_2d_prime_x_x_of_v, &
                fft_2d_prime_y_y_of_v, fft_2d_prime_x_y_of_v, fft_2d_prime_y_x_of_v)

            print *,"min_finder_run",minimum
            if (minimum > minval(p_v) .and. is_converged == 1) then
                exit
            else
                print *,"FAIL"
                step_size = step_size/2.0
            end if
        end do

        min_finder_approx = min(minimum, minval(p_v))

        call cpu_time(finish)
        print *,"min_finder_run Seconds: ",finish-start

    end function min_finder_approx

    real ( kind = 8 ) function min_finder_Newton(n, x0_index, p_v, step_size, is_converged, &
            fft_2d_prime_x_of_v, fft_2d_prime_y_of_v, fft_2d_prime_x_x_of_v, &
            fft_2d_prime_y_y_of_v, fft_2d_prime_x_y_of_v, fft_2d_prime_y_x_of_v)
        implicit none
        integer, intent(in)                 :: n, x0_index(2)
        real ( kind = 8 ), intent(in)       :: p_v(n,n), step_size
        integer, intent(inout)              :: is_converged
        complex ( kind = 8 ), dimension(n, n), intent(in)  ::  fft_2d_prime_x_of_v, fft_2d_prime_y_of_v, &
            fft_2d_prime_x_x_of_v, fft_2d_prime_y_y_of_v, fft_2d_prime_x_y_of_v, fft_2d_prime_y_x_of_v

        real ( kind = 8 )                   :: x(2), grad_y(2), x0(2), &
                                               hess_y(2,2), hess_inv_y(2,2), pi
        integer                             :: iter_count

        is_converged=0

        pi = 4.*atan(1.)
        x0 = x0_index*2.0*pi/n
        x = x0

        do iter_count = 1, 30
            grad_y(1) = triginterp_fft(n, x, fft_2d_prime_x_of_v)
            grad_y(2) = triginterp_fft(n, x, fft_2d_prime_y_of_v)

            hess_y(1,1) = triginterp_fft(n, x,fft_2d_prime_x_x_of_v)
            hess_y(1,2) = triginterp_fft(n, x,fft_2d_prime_x_y_of_v)
            hess_y(2,1) = triginterp_fft(n, x,fft_2d_prime_y_x_of_v)
            hess_y(2,2) = triginterp_fft(n, x,fft_2d_prime_y_y_of_v)

            if (abs((hess_y(1,1)*hess_y(2,2))-(hess_y(1,2)*hess_y(2,1)))<0.000001) then
                print *,"hess Singular!"
            end if

            hess_inv_y(1,1) = hess_y(2,2)
            hess_inv_y(1,2) = -1.*hess_y(1,2)
            hess_inv_y(2,1) = -1.*hess_y(2,1)
            hess_inv_y(2,2) = hess_y(1,1)
            hess_inv_y = 1.0/((hess_y(1,1)*hess_y(2,2))-(hess_y(1,2)*hess_y(2,1)))*hess_inv_y

!            print *,"Newton grad_y",MATMUL(hess_inv_y, grad_y)

!            print *,"Newton grad_y angle",ACOS(DOT_PRODUCT(-MATMUL(hess_inv_y, grad_y),-grad_y) &
!                                /(norm2(MATMUL(hess_inv_y, grad_y))*norm2(grad_y)))

!            print *,"hess_y",hess_y
!            print *,"hess_inv_y",hess_inv_y

            x = x - step_size*MATMUL(hess_inv_y, grad_y)

!            print *,"step_size",step_size
!            print *,"x",x
!            print *,"triginterp",triginterp(n, x, p_v)

!            print *,"desc count",iter_count

            if (maxval(abs(grad_y)) < 10d-13) then
                is_converged=1
                exit
            end if
        end do

        min_finder_Newton = triginterp(n, x, p_v)

    end function min_finder_Newton

    real ( kind = 8 ) function triginterp_caller(p_x1, p_x2)
        implicit none

        integer, parameter           :: n = 128
        real ( kind = 8 ), intent(in)                :: p_x1, p_x2

        real ( kind = 8 )                 :: p_x(n)
        real ( kind = 8 )                 :: p_v(n,n)
        real ( kind = 8 )                   :: pi
        integer                             :: j,k

        pi = 4.*atan(1.)

        p_x = [ (j, j = 1,n) ]
        p_x = p_x*2.0/n*pi

        do k = 1, n
            do j = 1, n
                p_v(k,j) = cos(p_x(k)) + cos(p_x(j))
            end do
        end do

        triginterp_caller = triginterp(n, [p_x1, p_x2], p_v)

    end function triginterp_caller

    real ( kind = 8 ) function triginterp(n, p_xi,p_v)
        integer, intent(in)          :: n

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        real ( kind = 8 ), intent(in)                :: p_v(n,n)

        real ( kind = 8 )                   :: pi
        complex ( kind = 8 )                ::          y(n,n)

        pi = 4.*atan(1.)

        y = p_v
        call dft_2d(n,y)
        if (minval(real(y)) > -0.01 .and. maxval(real(y)) < 0.01) then
            print *,"dft_2d of y = zero ",y
        end if
        triginterp = triginterp_fft(n,p_xi,y)

    end function triginterp


    real ( kind = 8 ) function triginterp_fft(n, p_xi,p_v_hat)
        integer, intent(in)          :: n

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        complex ( kind = 8 ), intent(in)             ::   p_v_hat(n,n)

        real ( kind = 8 )                   :: pi
        complex ( kind = 8 )                :: result
        integer                             :: k, j

        pi = 4.*atan(1.)

        result = 0
        do k = 1, n
            do j = 1, n
                result = result + exp((0., 1.) * (p_xi(1) * (-n/2+k) + p_xi(2) * (-n/2+j))) &
                            * p_v_hat(k,j)
            end do
        end do

        triginterp_fft = real(result/2.0/pi/n)
    end function triginterp_fft

end module max_finder
