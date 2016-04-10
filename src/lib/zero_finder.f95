module zero_finder
implicit none

contains

    real ( kind = 8 ) function triginterp_caller(p_xi)
        implicit none

        integer, parameter           :: n = 128
        real ( kind = 8 ), intent(in)                :: p_xi

        real ( kind = 8 )                 :: p_x(n)=0
        real ( kind = 8 )                 :: p_y(n)=0
        real ( kind = 8 )                   :: pi
        integer                             :: j

        pi = 4.*atan(1.)

        p_x = [ (j, j = 1,n) ]
        p_x = p_x*2.0/n*pi
        p_y = cos(p_x)
        triginterp_caller = triginterp(p_xi,p_x,p_y)

    end function triginterp_caller

    subroutine fft(x, y)
    implicit none
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(inout)                :: y(n)
        real ( kind = 8 ), intent(in)                   :: x(n)

        complex*16                          :: y_hat(n)
        integer                             :: k, j

        if (abs(x(2)-x(1)) > 0) then
            y_hat = 0;
            do k = 1, n
                do j = 1, n
                    y_hat(k) = y_hat(k) + exp((0., 1.) * (-1) * (k-n/2) * x(j))*y(j);
                end do
            end do
            y_hat = abs(x(2)-x(1))*y_hat;
            y = y_hat
        end if

    end subroutine fft


    real ( kind = 8 ) function triginterp(p_xi,p_x,p_y)
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(in)                :: p_xi
        real ( kind = 8 ), intent(in)                :: p_x(n)
        real ( kind = 8 ), intent(in)                :: p_y(n)

        real ( kind = 8 )                   :: y_hat(n+1), pi, y(n)
        complex*16                          :: temp, result
        integer                             :: k

        pi = 4.*atan(1.)

        y = p_y
        call fft(p_x, y)
        y_hat(1:n) = y(1:n)
        y_hat(1:n+1) = [y_hat(n), y_hat(1:n)]

        result = 0
        do k = 1, n+1
            temp = exp((0., 1.) * p_xi * (-n/2-1+k)) * y_hat(k)
            if (k==1 .OR. k==n+1) then
                result = result + temp/2.0
            else
                result = result + temp
            end if
        end do
        result = result/2.0/pi

        triginterp = result
    end function triginterp

end module zero_finder
