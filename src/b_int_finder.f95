module int_finder
    use dft_mod
    use, intrinsic :: iso_c_binding
    implicit none

contains

    real ( kind = 8 ) function int_finder_run_caller(n)
        implicit none

        integer, intent(in)                 :: n

        real ( kind = 8 )                   :: x(n)
        real ( kind = 8 )                   :: y(n,n)
        real ( kind = 8 )                   :: pi,h,hx,xa,xb,ret
        integer                             :: j,k

        pi = 4.*atan(1.)
        h = 45.0/n
        xa = ceiling((1 - n/2)*h)
        xb = floor((n/2)*h)
        hx=0.5

        x = [ (j,j=1,n) ]
        x = x-32.
        print *,"x",x

        do k = 1, n
            do j = 1, n
                y(k,j) = sin(x(k)) + sin(x(j))
            end do
        end do

        int_finder_run_caller = int_finder_run(n, h, y)

    end function int_finder_run_caller

    real (kind=8) function int_finder_run(n, h, p_v)
        implicit none
        integer, intent(in)                 :: n
        real ( kind = 8 ), intent(in)       :: h, p_v(n,n)
        complex (C_DOUBLE_COMPLEX), dimension(n,n) :: v_hat

        v_hat = p_v
        call fft2_F95(n, v_hat)

        int_finder_run = real(v_hat(1,1))*h*h

    end function int_finder_run

    real ( kind = 8 ) function triginterp_caller(p_x1, p_x2, n)
        implicit none

        real ( kind = 8 ), intent(in)                :: p_x1, p_x2
        integer, intent(in)           :: n

        real ( kind = 8 )                 :: p_x(n)
        real ( kind = 8 )                 :: p_y(n,n)
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

        triginterp_caller = triginterp([p_x1, p_x2],n,p_y)

    end function triginterp_caller

    real ( kind = 8 ) function triginterp(p_xi,n,p_y)

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        integer, intent(in)           :: n
        real ( kind = 8 ), intent(in)                :: p_y(n,n)

        real ( kind = 8 )                   :: pi
        complex (C_DOUBLE_COMPLEX)                ::          y(n,n)

        pi = 4.*atan(1.)

        y = p_y
        call dft_2d(n,y)
        triginterp = triginterp_fft(p_xi,n,y)

    end function triginterp


    real ( kind = 8 ) function triginterp_fft(p_xi,n,p_y_hat)

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        integer, intent(in)           :: n
        complex (C_DOUBLE_COMPLEX), intent(in)             ::   p_y_hat(n,n)

        real ( kind = 8 )                   :: pi
        complex (C_DOUBLE_COMPLEX)                ::          temp, result
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

end module int_finder
