module zero_finder
implicit none

contains

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

        triginterp_caller = triginterp([p_x1, p_x2],p_x,p_y)

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

    subroutine fft_2d(x, y)
    implicit none
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(inout)                :: y(n,n)
        real ( kind = 8 ), intent(in)                   :: x(n)

        complex*16                          :: y_hat(n,n)
        integer                             :: k, j, p, u

        if (abs(x(2)-x(1)) > 0) then
            y_hat = 0;
            do p = 1, n
                do u = 1, n
                    do k = 1, n
                        do j = 1, n
                            y_hat(p,u) = y_hat(p,u) + y(k,j) &
         * exp((0., 1.) * (-1) * ((p-n/2) * x(k) + (u-n/2) * x(j)))
                        end do
                    end do
                end do
            end do
            y_hat = abs(x(2)-x(1))*y_hat
            y = y_hat
        end if
    end subroutine fft_2d

    real ( kind = 8 ) function triginterp(p_xi,p_x,p_y)
        integer, parameter          :: n = 128

        real ( kind = 8 ), intent(in)                :: p_xi(2)
        real ( kind = 8 ), intent(in)                :: p_x(n)
        real ( kind = 8 ), intent(in)                :: p_y(n,n)

        real ( kind = 8 )                   :: y_hat(n,n), pi, y(n,n)
        complex*16                          :: temp, result
        integer                             :: k, j

        pi = 4.*atan(1.)

        y = p_y
        call fft_2d(p_x, y)
        y_hat(1:n,1:n) = y(1:n,1:n)

        result = 0
        do k = 1, n
            do j = 1, n
                temp = exp((0., 1.) * (p_xi(1) * (-n/2+k) + p_xi(2) * (-n/2+j))) * y_hat(k,j)
                result = result + temp
            end do
        end do
        result = result/2.0/pi/n

        triginterp = result
    end function triginterp

end module zero_finder
