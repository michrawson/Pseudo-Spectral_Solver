module fftw_prime
use, intrinsic :: iso_c_binding
implicit none

include 'fftw3.f03'

contains

    subroutine fft_prime_2d_partial_x_run(n, v, v_prime)
    implicit none
        integer, intent(in)                                 :: n
        complex (C_DOUBLE_COMPLEX), dimension(n, n), intent(in)   :: v
        real ( kind = 8 ), dimension(n, n), intent(inout)   :: v_prime
        integer                                             :: i

        do i=1,n
            call fft_prime_run( n, v(:, i), v_prime(:, i))
        end do

    end subroutine fft_prime_2d_partial_x_run

    subroutine fft_prime_2d_partial_y_run(n, v, v_prime)
    implicit none
        integer, intent(in)                              :: n
        complex (C_DOUBLE_COMPLEX), dimension(n, n), intent(in)  :: v
        real ( kind = 8 ), dimension(n, n), intent(inout)   :: v_prime
        integer                                             :: i
        complex (C_DOUBLE_COMPLEX)                          :: v2(n,n)

        v2 = transpose(v)
        do i=1,n
            call fft_prime_run( n, v2(:, i), v_prime(:, i))
        end do
        v_prime = transpose(v_prime)

    end subroutine fft_prime_2d_partial_y_run

    subroutine fft_prime_2d_run(n, v, v_prime)
    implicit none
        integer, intent(in)                                :: n
        complex (C_DOUBLE_COMPLEX), dimension(n, n), intent(inout)     :: v
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
        complex (C_DOUBLE_COMPLEX), dimension(n), intent(in)     :: v
        real ( kind = 8 ), dimension(n), intent(inout)  :: v_prime
        integer                         :: j
        type(C_PTR)                     :: plan_forward, plan_backward
        complex(C_DOUBLE_COMPLEX)       :: w_hat(n), ifft(n), v2(n), v_hat(n)

        plan_forward = fftw_plan_dft_1d ( N, v2, v_hat, FFTW_FORWARD, FFTW_MEASURE )
        v2(1:n) = v(1:n)
        call fftw_execute_dft ( plan_forward, v2, v_hat )
        call fftw_destroy_plan(plan_forward)

        plan_backward = fftw_plan_dft_1d ( N, w_hat, ifft, FFTW_BACKWARD, FFTW_MEASURE )
        w_hat(1:n) = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat(1:n)

        call fftw_execute_dft ( plan_backward, w_hat, ifft )
        call fftw_destroy_plan(plan_backward)

        ifft = ifft/n

        v_prime = real(ifft)

    end subroutine fft_prime_run

!    subroutine fft_prime_run(n, v, v_prime)
!    implicit none
!        integer, intent(in)                                :: n
!        complex *16, dimension(n), intent(inout)     :: v
!        real ( kind = 8 ), dimension(n), intent(inout)        :: v_prime
!        integer ( kind = 4 ), parameter :: fftw_forward = -1
!        integer ( kind = 4 ), parameter :: fftw_backward = +1
!        integer ( kind = 4 ), parameter :: FFTW_ESTIMATE = 64
!        integer ( kind = 4 ), parameter :: FFTW_MEASURE = 0
!        integer ( kind = 4 ), parameter :: FFTW_PATIENT = 32
!        integer                         :: j, plan_forward, plan_backward
!        complex *16, allocatable        :: w_hat(:), ifft(:), v2(:), v_hat(:)
!
!        allocate(w_hat(n))
!        allocate(ifft(n))
!        allocate(v2(n))
!        allocate(v_hat(n))
!
!        v2 = v
!        call dfftw_plan_dft_1d ( plan_forward, N, v2, v_hat, FFTW_FORWARD, FFTW_ESTIMATE )
!        call dfftw_execute_dft ( plan_forward, v2, v_hat )
!        call dfftw_destroy_plan(plan_forward)
!
!        w_hat = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat(1:n)
!        call dfftw_plan_dft_1d ( plan_backward, N, w_hat, ifft, FFTW_BACKWARD, FFTW_ESTIMATE )
!        call dfftw_execute_dft ( plan_backward, w_hat, ifft )
!        call dfftw_destroy_plan(plan_backward)
!
!        ifft = ifft/n
!
!        v_prime = real(ifft)
!
!        deallocate(w_hat)
!        deallocate(ifft)
!        deallocate(v2)
!        deallocate(v_hat)
!
!    end subroutine fft_prime_run

end module fftw_prime
