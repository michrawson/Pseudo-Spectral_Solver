module finite_differences
        implicit none

        contains

            double precision function inner_loop(N)
                implicit none
                integer, intent(in)     :: N
                double precision                :: error
                double precision, dimension(:), allocatable      :: x, u, u_prime, e
                double precision, dimension(:,:), allocatable    :: D, tmp1, tmp2
                integer                 :: j
                double precision                    :: pi, h

                allocate(x(N))
                allocate(u(N))
                allocate(u_prime(N))
                allocate(e(N))

                allocate(D(N,N))
                allocate(tmp1(N,N))
                allocate(tmp2(N,N))

                pi = 4.*atan(1.)
                error = 0
                e=1

                x = 0
                u = 0
                u_prime=0
                D=0
                tmp1=0
                tmp2=0

        !        PRINT *,"N",N

                h = 2*pi/N
        !        PRINT *,"h",h
                x = -pi + ((/ (j,j=1,N) /)*h)
        !        PRINT *,"x",x(1:N)
                u = exp(sin(x))
        !        PRINT *,"u",u(1:N)
                u_prime = cos(x) * u
        !        PRINT *,"u_prime",u_prime

                do j = 1, N-1
                    tmp1(j,j+1) = 2*e(j)/3
                end do
                tmp1(N,1) = 2*e(N)/3

        !        PRINT *,"tmp1",tmp1(1:N,1:N)

                do j = 1, N-2
                    tmp2(j,j+2) = e(j)/12
                end do
                tmp2(N-1,1) = e(N-1)/12
                tmp2(N,2) = e(N)/12

        !        PRINT *,"tmp2",tmp2

                D = tmp1 - tmp2

        !        PRINT *,"D",D

                D = (D - transpose(D)) / h

        !        PRINT *,"D",D

                error = maxval(abs( MATMUL(D, u) - u_prime))

                inner_loop = error

            end function inner_loop

            subroutine finite_differences_run(error)
                implicit none
                double precision, dimension(9+1), intent(inout) :: error
                integer, dimension(9+1)    :: Nvec
                integer                 :: i, N

                error = 0
                Nvec=0

                Nvec = (/ ((2**(i+2)),i=1,9+1) /)

                do i = 1, 9+1
                    N = Nvec(i)
                    error(i) = inner_loop(N)
                end do
            end subroutine finite_differences_run

end module finite_differences
