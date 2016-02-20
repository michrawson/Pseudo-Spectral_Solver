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
           
                complex ( kind = 8 ) in2(n)
                complex ( kind = 8 ) in3(n)
                complex ( kind = 8 ) out2(n)
                integer ( kind = 8 ) plan_backward
                integer ( kind = 8 ) plan_forward

  integer ( kind = 4 ), parameter :: fftw_r2hc = 0
  integer ( kind = 4 ), parameter :: fftw_hc2r = 1
  integer ( kind = 4 ), parameter :: fftw_dht = 2
  integer ( kind = 4 ), parameter :: fftw_redft00 = 3
  integer ( kind = 4 ), parameter :: fftw_redft01 = 4
  integer ( kind = 4 ), parameter :: fftw_redft10 = 5
  integer ( kind = 4 ), parameter :: fftw_redft11 = 6
  integer ( kind = 4 ), parameter :: fftw_rodft00 = 7
  integer ( kind = 4 ), parameter :: fftw_rodft01 = 8
  integer ( kind = 4 ), parameter :: fftw_rodft10 = 9
  integer ( kind = 4 ), parameter :: fftw_rodft11 = 10
  integer ( kind = 4 ), parameter :: fftw_forward = -1
  integer ( kind = 4 ), parameter :: fftw_backward = +1
  integer ( kind = 4 ), parameter :: fftw_measure = 0
  integer ( kind = 4 ), parameter :: fftw_destroy_input = 1
  integer ( kind = 4 ), parameter :: fftw_unaligned = 2
  integer ( kind = 4 ), parameter :: fftw_conserve_memory = 4
  integer ( kind = 4 ), parameter :: fftw_exhaustive = 8
  integer ( kind = 4 ), parameter :: fftw_preserve_input = 16
  integer ( kind = 4 ), parameter :: fftw_patient = 32
  integer ( kind = 4 ), parameter :: fftw_estimate = 64
  integer ( kind = 4 ), parameter :: fftw_estimate_patient = 128
  integer ( kind = 4 ), parameter :: fftw_believe_pcost = 256
  integer ( kind = 4 ), parameter :: fftw_dft_r2hc_icky = 512
  integer ( kind = 4 ), parameter :: fftw_nonthreaded_icky = 1024
  integer ( kind = 4 ), parameter :: fftw_no_buffering = 2048
  integer ( kind = 4 ), parameter :: fftw_no_indirect_op = 4096
  integer ( kind = 4 ), parameter :: fftw_allow_large_generic = 8192
  integer ( kind = 4 ), parameter :: fftw_no_rank_splits = 16384
  integer ( kind = 4 ), parameter :: fftw_no_vrank_splits = 32768
  integer ( kind = 4 ), parameter :: fftw_no_vrecurse = 65536
  integer ( kind = 4 ), parameter :: fftw_no_simd = 131072

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

                        in2 = 1
                        out2 = 0
                        in3 = 0
                        call dfftw_plan_dft_1d_ ( plan_forward, N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE )

                        call dfftw_execute_ ( plan_forward )

                        call dfftw_plan_dft_1d_ ( plan_backward, N, out2, in3, FFTW_BACKWARD, FFTW_ESTIMATE )

                        call dfftw_execute_ ( plan_backward )

                        PRINT *,"diff",abs(abs(in2) - abs(in3))

                        PRINT *,"rel diff",abs(in2/in3)

                        return

                        call cfft1i ( N, wsave, 3*N, ier )
                        PRINT *,"ier",ier

                        v_hat = v
                        call cfft1f( N, 1, v_hat, N, wsave, 3*N, work, 3*N, ier)
                        PRINT *,"ier",ier

                        PRINT *,"diff",abs(v_hat - v)

                        call cfft1b( N, 1, v_hat, N, wsave, 3*N, work, 3*N, ier)
                        PRINT *,"ier",ier

                        PRINT *,"diff",abs(v_hat - v)
                        return

                        w_hat = (0., 1.) * (/ (j,j=0,N-1) /) * v_hat
                        
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
