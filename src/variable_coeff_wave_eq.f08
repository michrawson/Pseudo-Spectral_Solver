module variable_coeff_wave_eq
        implicit none

        contains

            subroutine variable_coeff_wave_eq_run(x,tdata,result)
                implicit none
                
                integer, parameter      :: N = 128, tmax = 8
                
                real, dimension(N), intent(out)                 :: x
                real, dimension(N)      :: c, v, vold, v_hat, vnew, w
                complex, dimension(N)   :: w_hat, ifft
                real, dimension(:,:), allocatable, intent(out)  :: result
                real, dimension(:), allocatable, intent(out)    :: tdata
                real, dimension(:), allocatable                 :: wsave, work

                integer                 :: i, j, plotgap, nplots, ier
                real                    :: pi, h, t, dt, tplot

                pi = 4.*atan(1.)
                h = 2*pi/N
                x = h*(/ (j,j=1,N) /)
                t = 0
                dt = h/4
                c = .2 + sin(x-1)**2.
                v = exp(-100*((x-1)**2.))
                vold = exp(-100*((x-.2*dt-1)**2.))

                tplot = .15
                plotgap = nint(tplot/dt)
                dt = tplot/plotgap
                nplots = nint(tmax/tplot)

                allocate(tdata(nplots))
                tdata = -1

                allocate(result(nplots+1,N))
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