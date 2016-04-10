module max_finder
use dft_mod
implicit none



contains

    subroutine max_finder_run(x0_index, p_y, result)
        implicit none
        integer, parameter                  :: n = 128
        integer, intent(in)       :: x0_index(2)
        real ( kind = 8 ), intent(in)       :: p_y(n,n)
        real ( kind = 8 ), intent(inout)    :: result
        real ( kind = 8 )                   :: pi, step_size, x(2), grad_y(2), x0(2), f_x
        complex ( kind = 8 )                ::    fft_2d_prime_x_of_y(n,n), fft_2d_prime_y_of_y(n,n)
        integer                             :: iter_count=0

        pi = 4.*atan(1.)

        x0 = x0_index*2.0*pi/n

        x = x0

        print *,"x",x
        print *,"max",triginterp(x,p_y)

        step_size = 1.0

        fft_2d_prime_x_of_y = p_y
        call fft_prime_2d_partial_x(fft_2d_prime_x_of_y)

        fft_2d_prime_y_of_y = p_y
        call fft_prime_2d_partial_y(fft_2d_prime_y_of_y)

        grad_y = [ triginterp_fft(x, fft_2d_prime_x_of_y), triginterp_fft(x, fft_2d_prime_y_of_y) ]
        do while (maxval(abs(grad_y))>0.000001)
            step_size = step_size*(1.5**2)
            grad_y = [ triginterp_fft(x, fft_2d_prime_x_of_y), triginterp_fft(x, fft_2d_prime_y_of_y) ]
            print *,"grad",grad_y
            f_x = triginterp(x, p_y)
            do while (f_x > triginterp(x + step_size * grad_y, p_y))
                step_size = step_size/1.5
            end do
            print *,"step_size",step_size
            print *,"triginterp(x)",f_x
            x = x + step_size * grad_y
            iter_count = iter_count + 1
        end do

        print *,"x",x
        print *,"max",triginterp(x,p_y)
        print *,"iter_count",iter_count

        result = triginterp(x,p_y)
    end subroutine max_finder_run

    real ( kind = 8 ) function triginterp_caller(p_x1, p_x2)
        implicit none

        integer, parameter           :: n = 128
        real ( kind = 8 ), intent(in)                :: p_x1, p_x2

        real ( kind = 8 )                 :: p_x(n)=0
        real ( kind = 8 )                 :: p_y(n,n)=0
        real ( kind = 8 )                   :: pi
        integer                             :: j,k

        pi = 4.*atan(1.)

        p_x = [ (j, j = 1,n) ]
        p_x = p_x*2.0/n*pi

        do k = 1, n
            do j = 1, n
                p_y(k,j) = cos(p_x(k)) + cos(p_x(j))
            end do
        end do

        triginterp_caller = triginterp([p_x1, p_x2],p_y)

    end function triginterp_caller

    real ( kind = 8 ) function triginterp(p_xi,p_y)
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        real ( kind = 8 ), intent(in)                :: p_y(n,n)

        real ( kind = 8 )                   :: pi
        complex ( kind = 8 )                ::          y(n,n)

        pi = 4.*atan(1.)

        y = p_y
        call dft_2d(y)
        triginterp = triginterp_fft(p_xi,y)

    end function triginterp


    real ( kind = 8 ) function triginterp_fft(p_xi,p_y_hat)
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        complex ( kind = 8 ), intent(in)             ::   p_y_hat(n,n)

        real ( kind = 8 )                   :: pi
        complex ( kind = 8 )                ::          temp, result
        integer                             :: k, j

        pi = 4.*atan(1.)

        result = 0
        do k = 1, n
            do j = 1, n
                temp = exp((0., 1.) * (p_xi(1) * (-n/2+k) + p_xi(2) * (-n/2+j))) * p_y_hat(k,j)
                result = result + temp
            end do
        end do
        result = result/2.0/pi/n

        triginterp_fft = real(result)
    end function triginterp_fft

end module max_finder
