module fft_prime
implicit none

contains

    subroutine fft_prime_run(v, v_prime)
    implicit none
        integer, parameter                              :: n = 128
        real ( kind = 8 ), dimension(n), intent(in)     :: v
        real ( kind = 8 ), dimension(n), intent(inout)  :: v_prime
        integer ( kind = 4 ), parameter :: fftw_forward = -1
        integer ( kind = 4 ), parameter :: fftw_backward = +1
        integer ( kind = 4 ), parameter :: fftw_estimate = 64
        integer                         :: j, plan_forward, plan_backward
        complex*16 w_hat(n), ifft(n), v2(n), v_hat(n)

        w_hat=0
        ifft=0
        v2=0
        v_hat=0

        v2(1:n) = v(1:n)
        call dfftw_plan_dft_1d_ ( plan_forward, N, v2, v_hat, FFTW_FORWARD, FFTW_ESTIMATE )
        call dfftw_execute_ ( plan_forward )

        w_hat(1:n) = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat(1:n)

        call dfftw_plan_dft_1d_ ( plan_backward, N, w_hat, ifft, FFTW_BACKWARD, FFTW_ESTIMATE )
        call dfftw_execute_ ( plan_backward )
        ifft(1:n) = ifft(1:n)/n

        v_prime(1:n) = real(ifft(1:n))

    end subroutine fft_prime_run

end module fft_prime