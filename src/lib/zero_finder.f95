module zero_finder
implicit none

contains

    subroutine zero_finder_run(x0, f, result)
    use fft_prime
    implicit none
        integer, parameter          :: n = 128

        integer, intent(in)                   :: x0(2)
        real ( kind = 8 ), intent(in)                   :: f(n,n)
        real ( kind = 8 ), intent(inout) :: result

        real ( kind = 8 ), dimension(n)                 :: f_prime_x, f_prime_prime_x, f_prime_prime_prime_x
        real ( kind = 8 ), dimension(n)                 :: f_prime_y, f_prime_prime_y, f_prime_prime_prime_y
        real ( kind = 8 ), dimension(2)                 :: f_prime, f_prime_prime, f_prime_prime_prime    ! at x0
        real ( kind = 8 ), dimension(2)                 :: h

        call fft_prime_run(f(x0(1),:), f_prime_x)
        call fft_prime_run(f(:,x0(2)), f_prime_y)

        f_prime = [f_prime_x(x0(1)), f_prime_y(x0(2))]

        call fft_prime_run(f_prime_x, f_prime_prime_x)
        call fft_prime_run(f_prime_y, f_prime_prime_y)

        f_prime_prime = [f_prime_prime_x(x0(1)), f_prime_prime_y(x0(2))]

        call fft_prime_run(f_prime_prime_x, f_prime_prime_prime_x)
        call fft_prime_run(f_prime_prime_y, f_prime_prime_prime_y)

        f_prime_prime_prime = [f_prime_prime_prime_x(x0(1)), f_prime_prime_prime_y(x0(2))]

        h=0
        if (f_prime_prime_prime(1) /= 0) then
        h(1) = ((- f_prime_prime(1) + sqrt((f_prime_prime(1)**2.0) - 2.0*f_prime_prime_prime(1)*f_prime(1))) &
                /f_prime_prime_prime(1))
        end if

        if (f_prime_prime_prime(2) /= 0) then
        h(2) = ((- f_prime_prime(2) + sqrt((f_prime_prime(2)**2.0) - 2.0*f_prime_prime_prime(2)*f_prime(2))) &
                /f_prime_prime_prime(2))
        end if

        PRINT *,"f_prime:", f_prime
        PRINT *,"f_prime_prime:", f_prime_prime
        PRINT *,"f_prime_prime_prime:", f_prime_prime_prime

        PRINT *,"h:", h

        result = f(x0(1),x0(2)) + h(1)*f_prime(1) + h(1)**2.0 * f_prime_prime(1) / 2.0 &
                                + h(2)*f_prime(2) + h(2)**2.0 * f_prime_prime(2) / 2.0

    end subroutine zero_finder_run

end module zero_finder
