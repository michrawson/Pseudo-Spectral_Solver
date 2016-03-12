module variable_coeff_wave_eq_pseudo_rk
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_rk_run(x,tdata,result)
    implicit none
        integer, parameter          :: n = 128
        integer, parameter          :: tmax = 8
        real ( kind = 8 ), parameter :: tplot = 0.15

        real ( kind = 8 ), dimension(n), intent(out)                   :: x
        real ( kind = 8 ), dimension(int(tmax/tplot)+1), intent(out)     :: tdata
        real ( kind = 8 ), dimension(int(tmax/tplot)+1,n), intent(out) :: result

        real ( kind = 8 ), dimension(n)              :: c, w, v, k1,k2,k3,k4
        complex ( kind = 8 ), dimension(n)          :: w_hat, ifft, v2, v_hat
        real ( kind = 8 ), dimension(n,n)            :: D

        integer               :: i, j, plotgap, nplots, plan_backward, plan_forward
        real ( kind = 8 )      :: pi, h, t, dt

        integer ( kind = 4 ), parameter :: fftw_forward = -1
        integer ( kind = 4 ), parameter :: fftw_backward = +1
        integer ( kind = 4 ), parameter :: fftw_estimate = 64

        pi = 4.*atan(1.)
        h = 2.*pi/N
        x = h*(/ (j,j=1,N) /)
        t = 0
        dt = h/4.
        c = 1.0
        !.2 + sin(x-1.)**2.
        v = exp(-100*((x-1)**2.))

        plotgap = nint(tplot/dt)
        dt = tplot/plotgap
        nplots = nint(tmax/tplot)

        tdata = -1
        tdata(1) = t

        result = 0
        result(1,1:N) = v

        do i=0,n-1
            do j=1,n
                if (i==0) then
                    D(j,j) = 0
                else
                    if (j+i <= n) then
                        D(j+i,j) = ((-1.)**(i))/2./tan(i*h/2.)
                    end if
                    if (j-i >= 1) then
                        D(j-i,j) = ((-1.)**(i+1))/2./tan(i*h/2.)
                    end if
                end if
            end do
        end do


        do i=1,nplots
            do j = 1,plotgap
                t = t+dt

                k1 = -1.0*MATMUL(D, v)
                k2 = -1.0*MATMUL(D, v + dt*k1/2)
                k3 = -1.0*MATMUL(D, v + dt*k2/2)
                k4 = -1.0*MATMUL(D, v + dt*k3)
                v = v + dt/6.0*(k1 + 2.*k2 + 2.*k3 + k4)

            end do
            result(i+1,:) = v
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_rk_run

end module variable_coeff_wave_eq_pseudo_rk
