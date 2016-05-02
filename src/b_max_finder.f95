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
        integer                             :: j,k

        pi = 4.*atan(1.)

        p_x = [ (j, j = 1,n) ]
        p_x = p_x*2.0/n*pi

        do k = 1, n
            do j = 1, n
                p_v(k,j) = -cos(p_x(k)) - cos(p_x(j))
            end do
        end do

        max_finder_run_caller = max_finder_run(n, [p_x1, p_x2], p_v)

    end function max_finder_run_caller

    real ( kind = 8 ) function max_finder_run(n, x0_index, p_v)
        implicit none
        integer, intent(in)                  :: n
        integer, intent(in)       :: x0_index(2)
        real ( kind = 8 ), intent(in)       :: p_v(n,n)
        real ( kind = 8 )                   :: pi, x(2), grad_y(2), x0(2), &
            hess_y(2,2), hess_inv_y(2,2)
        complex ( kind = 8 ), dimension(n, n)  ::  fft_2d_prime_x_of_y, fft_2d_prime_y_of_y, &
            fft_2d_prime_x_x_of_y, fft_2d_prime_y_y_of_y, fft_2d_prime_x_y_of_y, fft_2d_prime_y_x_of_y
        integer                             :: iter_count

        pi = 4.*atan(1.)

        x0 = x0_index*2.0*pi/n

        x = x0

        fft_2d_prime_x_of_y = p_v
        call fft_prime_2d_partial_x(n, fft_2d_prime_x_of_y)

        fft_2d_prime_y_of_y = p_v
        call fft_prime_2d_partial_y(n, fft_2d_prime_y_of_y)


        fft_2d_prime_x_x_of_y = p_v
        call fft_prime_prime_2d_partial_x_x(n, fft_2d_prime_x_x_of_y)

        fft_2d_prime_y_y_of_y = p_v
        call fft_prime_prime_2d_partial_y_y(n, fft_2d_prime_y_y_of_y)

        fft_2d_prime_x_y_of_y = p_v
        call fft_prime_prime_2d_partial_x_y(n, fft_2d_prime_x_y_of_y)

        fft_2d_prime_y_x_of_y = fft_2d_prime_x_y_of_y
!        fft_2d_prime_y_x_of_y = p_v
!        call fft_prime_prime_2d_partial_y_x(n, fft_2d_prime_y_x_of_y)

        grad_y = [  triginterp_fft(n, x, fft_2d_prime_x_of_y), &
            triginterp_fft(n, x, fft_2d_prime_y_of_y) ]

        iter_count = 0
        do while (maxval(abs(grad_y)) > 10d-12)
            grad_y(1) = triginterp_fft(n, x, fft_2d_prime_x_of_y)
            grad_y(2) = triginterp_fft(n, x, fft_2d_prime_y_of_y)

            hess_y(1,1) = triginterp_fft(n, x,fft_2d_prime_x_x_of_y)
            hess_y(1,2) = triginterp_fft(n, x,fft_2d_prime_x_y_of_y)
            hess_y(2,1) = triginterp_fft(n, x,fft_2d_prime_y_x_of_y)
            hess_y(2,2) = triginterp_fft(n, x,fft_2d_prime_y_y_of_y)
            hess_inv_y(1,1) = hess_y(2,2)
            hess_inv_y(1,2) = -1*hess_y(1,2)
            hess_inv_y(2,1) = -1*hess_y(2,1)
            hess_inv_y(2,2) = hess_y(1,1)
            hess_inv_y = 1.0/((hess_y(1,1)*hess_y(2,2))-(hess_y(1,2)*hess_y(2,1)))*hess_inv_y

            x = x - MATMUL(hess_inv_y, grad_y)

            iter_count = iter_count + 1
        !            print *,"desc count",iter_count
        end do

        max_finder_run = triginterp(n, x, p_v)

    end function max_finder_run

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

        triginterp_caller = triginterp(n, [p_x1, p_x2],p_v)

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
                result = result + exp((0., 1.) * (p_xi(1) * (-n/2+k) + p_xi(2) * (-n/2+j))) * p_v_hat(k,j)
            end do
        end do

        triginterp_fft = real(result/2.0/pi/n)
    end function triginterp_fft

end module max_finder
