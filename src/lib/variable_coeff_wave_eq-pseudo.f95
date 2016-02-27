module variable_coeff_wave_eq_pseudo
implicit none

contains

    subroutine variable_coeff_wave_eq_pseudo_run(x,tdata,result)
    implicit none
        integer, parameter          :: n = 128
        integer, parameter          :: tmax = 8
        double precision, parameter :: tplot = 0.15

        double precision, dimension(n), intent(out)                   :: x
        double precision, dimension(int(tmax/tplot)+1), intent(out)     :: tdata
        double precision, dimension(int(tmax/tplot)+1,n), intent(out) :: result

        double precision, dimension(n)              :: c, w, v, vold, vnew
        complex ( kind = 8 ), dimension(N)          :: w_hat, ifft, v2, v_hat
        double precision, dimension(n,n)            :: D


        integer               :: i, j, plotgap, nplots, k, plan_backward, plan_forward
        double precision      :: pi, h, t, dt

        integer ( kind = 4 ), parameter :: fftw_forward = -1
        integer ( kind = 4 ), parameter :: fftw_backward = +1
        integer ( kind = 4 ), parameter :: fftw_estimate = 64

        pi = 4.*atan(1.)
        h = 2.*pi/N
        x = h*(/ (j,j=1,N) /)
        t = 0
        dt = h/4.
        c = .2 + sin(x-1.)**2.
        v = exp(-100*((x-1)**2.))
        vold = exp(-100*((x-.2*dt-1)**2.))

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

                w = MATMUL(D, v)

                if (1==0) then
                    v_hat=0
                    v2 = v
                    call dfftw_plan_dft_1d_ ( plan_forward, N, v2, v_hat, FFTW_FORWARD, FFTW_ESTIMATE )
                    call dfftw_execute_ ( plan_forward )

                    w_hat = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat

                    call dfftw_plan_dft_1d_ ( plan_backward, N, w_hat, ifft, FFTW_BACKWARD, FFTW_ESTIMATE )
                    call dfftw_execute_ ( plan_backward )
                    ifft = ifft/n

                    w = real(ifft)
                end if

                vnew = vold - 2*dt*c*w
                vold = v
                v = vnew

            end do
            result(i+1,:) = v
            tdata(i+1) = t
        end do

    end subroutine variable_coeff_wave_eq_pseudo_run

end module variable_coeff_wave_eq_pseudo
