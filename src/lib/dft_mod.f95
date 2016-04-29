module dft_mod
implicit none

contains

    subroutine fft_prime_prime_2d_partial_y_y(v)
    implicit none
        integer, parameter                                 :: n = 128
        complex (kind = 8), dimension(n, n), intent(inout)      :: v
        integer               :: i, j

        call dft_2d(v)
        do i=1,n
            v(i,:) = (-1)*(/ (j,j=-N/2+1,N/2-1), 0 /) * (/ (j,j=-N/2+1,N/2-1), 0 /) * v(i,:)
        end do
    end subroutine fft_prime_prime_2d_partial_y_y

    subroutine fft_prime_prime_2d_partial_x_x(v)
    implicit none
        integer, parameter                                 :: n = 128
        complex (kind = 8), dimension(n, n), intent(inout)      :: v
        integer               :: i, j

        call dft_2d(v)
        do i=1,n
            v(:,i) = (-1)*(/ (j,j=-N/2+1,N/2-1), 0 /) * (/ (j,j=-N/2+1,N/2-1), 0 /) * v(:,i)
        end do
    end subroutine fft_prime_prime_2d_partial_x_x

    subroutine fft_prime_prime_2d_partial_x_y(v)
    implicit none
        integer, parameter                                 :: n = 128
        complex (kind = 8), dimension(n, n), intent(inout)      :: v
        integer               :: i, j

        call dft_2d(v)
        do i=1,n
            v(:,i) = (0., 1.) * (/ (j,j=-N/2+1,N/2-1), 0 /) * v(:,i)
        end do
        do i=1,n
            v(i,:) = (0., 1.) * (/ (j,j=-N/2+1,N/2-1), 0 /) * v(i,:)
        end do

    end subroutine fft_prime_prime_2d_partial_x_y

    subroutine fft_prime_prime_2d_partial_y_x(v)
    implicit none
        integer, parameter                                 :: n = 128
        complex (kind = 8), dimension(n, n), intent(inout)      :: v
        integer               :: i, j

        call dft_2d(v)
        do i=1,n
            v(i,:) = (0., 1.) * (/ (j,j=-N/2+1,N/2-1), 0 /) * v(i,:)
        end do
        do i=1,n
            v(:,i) = (0., 1.) * (/ (j,j=-N/2+1,N/2-1), 0 /) * v(:,i)
        end do

    end subroutine fft_prime_prime_2d_partial_y_x

!------------------------------------------------------------------------------------------
    subroutine fft_prime_2d_partial_x(v)
    implicit none
        integer, parameter                                 :: n = 128
        complex (kind = 8), dimension(n, n), intent(inout)      :: v
        integer               :: i, j

        call dft_2d(v)
        do i=1,n
            v(:,i) = (0., 1.) * (/ (j,j=-N/2+1,-1), (j,j=0,N/2-1), 0 /) * v(:,i)
        end do
    end subroutine fft_prime_2d_partial_x

    subroutine fft_prime_2d_partial_y(v)
    implicit none
        integer, parameter                                 :: n = 128
        complex (kind = 8), dimension(n, n), intent(inout)      :: v
        integer              :: i, j

        call dft_2d(v)
        do i=1,n
            v(i,:) = (0., 1.) * (/ (j,j=-N/2+1,-1), (j,j=0,N/2-1), 0 /) * v(i,:)
        end do
    end subroutine fft_prime_2d_partial_y
!------------------------------------------------------------------------------------------
    subroutine dft(y)
    implicit none
        integer, parameter                  :: n = 128

        complex (kind = 8), intent(inout)    :: y(n)
        real ( kind = 8 )                   :: x(n), pi

        complex (kind = 8)                 :: y_hat(n)
        integer                             :: k, j

        pi = 4.*atan(1.)

        x = [ (j, j = 1,n) ]
        x = x*2.0/n*pi

        y_hat = 0;
        do k = 1, n
            do j = 1, n
                y_hat(k) = y_hat(k) + exp((0., 1.) * (-1.0) * (k-n/2) * x(j))*y(j);
            end do
        end do
        y_hat = abs(x(2)-x(1))*y_hat;
        y = y_hat
    end subroutine dft

    subroutine dft_2d(y)
    implicit none
        integer, parameter          :: n = 128

        complex (kind = 8), intent(inout)                :: y(n,n)
        real ( kind = 8 )                   :: x(n), pi

        complex (kind = 8)                          :: y_hat(n,n)
        integer                             :: k, j, p, u


        pi = 4.*atan(1.)

        x = [ (j, j = 1,n) ]
        x = x*2.0/n*pi

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
    end subroutine dft_2d

end module dft_mod
