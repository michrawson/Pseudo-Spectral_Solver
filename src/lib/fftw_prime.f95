module fftw_prime
implicit none

contains

    subroutine fft_prime_2d_partial_x_run(n, v, v_prime)
    implicit none
        integer, intent(in)                                :: n
        real ( kind = 8 ), dimension(n, n), intent(in)     :: v
        real ( kind = 8 ), dimension(n, n), intent(inout)  :: v_prime
        integer               :: i

        do i=1,n
            call fft_prime_run( n, v(:,i), v_prime(:,i))
        end do

    end subroutine fft_prime_2d_partial_x_run

    subroutine fft_prime_2d_partial_y_run(n, v, v_prime)
    implicit none
        integer, intent(in)                                :: n
        real ( kind = 8 ), dimension(n, n), intent(in)     :: v
        real ( kind = 8 ), dimension(n, n), intent(inout)  :: v_prime
        integer               :: i

        do i=1,n
            call fft_prime_run( n, v(i, :), v_prime(i, :))
        end do

    end subroutine fft_prime_2d_partial_y_run

    subroutine fft_prime_2d_run(n, v, v_prime)
    implicit none
        integer, intent(in)                                :: n
        real ( kind = 8 ), dimension(n, n), intent(in)     :: v
        real ( kind = 8 ), dimension(n, n), intent(inout)  :: v_prime
        real ( kind = 8 ), dimension(n, n)                 :: v_temp_x, v_temp_y
        integer               :: i, j

        do i=1,n
            call fft_prime_run( n, v(:,i), v_temp_x(:,i))
            call fft_prime_run( n, v(i,:), v_temp_y(i,:))
        end do

        do i=1,n
            do j = 1,n
                v_prime(i,j) = v_temp_x(i,j) + v_temp_y(i,j)
            end do
        end do

    end subroutine fft_prime_2d_run

    subroutine fft_prime_run(n, v, v_prime)
    implicit none
        integer, intent(in)                                :: n
        real ( kind = 8 ), dimension(n), intent(in)     :: v
        real ( kind = 8 ), dimension(n), intent(inout)  :: v_prime
        integer ( kind = 4 ), parameter :: fftw_forward = -1
        integer ( kind = 4 ), parameter :: fftw_backward = +1
        integer ( kind = 4 ), parameter :: fftw_estimate = 64
        integer ( kind = 4 ), parameter :: FFTW_MEASURE = 0
        integer                         :: j, plan_forward, plan_backward
        complex (kind = 8) :: w_hat(n), ifft(n), v2(n), v_hat(n)

        call dfftw_plan_dft_1d_ ( plan_forward, N, v2, v_hat, FFTW_FORWARD, fftw_estimate )
        v2(1:n) = v(1:n)
        call dfftw_execute_ ( plan_forward )
!        call dfftw_destroy_plan_(plan_forward)

        call dfftw_plan_dft_1d_ ( plan_backward, N, w_hat, ifft, FFTW_BACKWARD, fftw_estimate )

        w_hat(1:n) = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat(1:n)

        call dfftw_execute_ ( plan_backward )
!        call dfftw_destroy_plan_(plan_backward)
        ifft(1:n) = ifft(1:n)/n

        v_prime(1:n) = real(ifft(1:n))

    end subroutine fft_prime_run

end module fftw_prime
