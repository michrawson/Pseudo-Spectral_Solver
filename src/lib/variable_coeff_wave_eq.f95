module variable_coeff_wave_eq
        implicit none
        integer, parameter:: dp=kind(0.d0)

        contains
            subroutine variable_coeff_wave_eq_run(x,tdata,result)
                implicit none
                
                integer, parameter          :: N = 128
                integer, parameter          :: tmax = 8
                double precision, parameter :: tplot = 0.15
                
                double precision, dimension(N), intent(out)                     :: x
                double precision, dimension(int(tmax/tplot)), intent(out)       :: tdata
                double precision, dimension(int(tmax/tplot)+1,N), intent(out)   :: result

                double precision, dimension(N)              :: c, v, vold, v_hat, vnew, w
                double complex, dimension(N)                :: w_hat, ifft
                double precision, dimension(:), allocatable :: wsave, work

                integer                         :: i, j, plotgap, nplots, ier
                double precision                :: pi, h, t, dt

                pi = 4.*atan(1.)
                h = 2*pi/N
                x = h*(/ (j,j=1,N) /)
                t = 0
                dt = h/4
                c = .2 + sin(x-1)**2.
                v = exp(-100*((x-1)**2.))
                vold = exp(-100*((x-.2*dt-1)**2.))

                plotgap = nint(tplot/dt)
                dt = tplot/plotgap
                nplots = nint(tmax/tplot)

                tdata = -1

                result = 0
                result(1,1:N) = v

                w_hat = 0

                allocate(wsave(3*N))
                wsave = 0
                allocate(work(3*N))
                work = 0
                
                do i=1,nplots
                    do j = 1,plotgap
                        t = t+dt

                        call cfft1i ( N, wsave, 3*N, ier )

                        v_hat = v
                        call cfft1f( N, 1, v_hat, N, wsave, 3*N, work, 3*N, ier)

                        w_hat = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat
                        
                        ifft = w_hat
                        call cfft1b( N, 1, ifft, N, wsave, 3*N, work, 3*N, ier)

                        w = real(ifft)
                        vnew = vold - 2*dt*c*w
                        vold = v
                        v = vnew
                    end do
                    result(i+1,:) = v
                    tdata(i) = t
                end do

            end subroutine variable_coeff_wave_eq_run

end module variable_coeff_wave_eq
