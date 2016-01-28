subroutine run(error)
    integer, parameter      :: dp = kind(0.d0)
    integer, parameter      :: MAX_N = 2**(9+2)
    real(dp), intent(inout) :: error(9)
    real, dimension(MAX_N)  :: x, u, u_prime, e
    real, dimension(9)      :: Nvec
    real, dimension(MAX_N,MAX_N) :: D, tmp1, tmp2
    integer                 :: i, j, N
    real                    :: pi, h

    pi = 4.*atan(1.)
    x = 0
    u = 0
    error = 0
    u_prime=0
    Nvec=0
    e=1
    D=0
    tmp1=0
    tmp2=0

    Nvec = (/ (2**(i+2),i=1,9) /)

    do i = 1, 9
        N = Nvec(i)
        h = 2*pi/N
        x(1:N) = -pi + ((/ (j,j=1,N) /)*h)
        u(1:N) = exp(sin(x(1:N)))
        u_prime(1:N) = (/ (cos(x(j)),j=1,N) /)

        do j = 1, N-1
            tmp1(j,j+1) = 2*e(j)/3
        end do
        tmp1(N,1) = 2*e(N)/3

        do j = 1, N-2
            tmp2(j,j+2) = e(j)/12
        end do
        tmp2(N-1,1) = e(N-1)/12
        tmp2(N,2) = e(N)/12

        D(1:N, 1:N) = tmp1(1:N, 1:N) - tmp2(1:N, 1:N)

        D(1:N, 1:N) = (D(1:N, 1:N) - transpose(D(1:N, 1:N))) / h

        error(i) = maxval(abs( MATMUL(D(1:N, 1:N), u(1:N)) - u_prime(1:N)))

    end do

end subroutine

program finite_differences
    implicit none
    interface
        subroutine run(error)
            integer, parameter      :: dp = kind(0.d0)
            integer, parameter      :: MAX_N = 2**(9+2)
            real(dp), intent(inout) :: error(9)
            real, dimension(MAX_N)  :: x, u, u_prime, e
            real, dimension(9)      :: Nvec
            real, dimension(MAX_N,MAX_N) :: D, tmp1, tmp2
            integer                 :: i, j, N
            real                    :: pi, h
        end subroutine run
    end interface

    integer, parameter      :: dp = kind(0.d0)
    integer, parameter      :: MAX_N = 2**(9+2)
    real(dp)                :: error(9)

    error = 0

    call run(error)

    PRINT *,error

end program finite_differences
