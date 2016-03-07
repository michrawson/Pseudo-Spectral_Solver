
program main
    use finite_differences
!    use variable_coeff_wave_eq
    use variable_coeff_wave_eq_pseudo
    implicit none

    real ( kind = 8 )                                        :: error(9+1) = 0

    integer, parameter                          :: N = 30, tmax = 8
    real ( kind = 8 ), parameter                         :: tplot = 0.15
!    real ( kind = 8 ), dimension(N)                      :: x = 0
    real ( kind = 8 ), dimension(int(tmax/tplot)+1,N)    :: result
    real ( kind = 8 ), dimension(int(tmax/tplot))      :: tdata

!    call finite_differences_run(error)
!    PRINT *,"error",error(1:10)

!    call variable_coeff_wave_eq_pseudo_run(x,tdata,result)
!    PRINT *,"x",x
!    PRINT *,"tdata",tdata
!    PRINT *,"result",result(1,:)
!    PRINT *,"result",result(2,:)
!    PRINT *,"result",result(:,10)
!    PRINT *,"result",result(:,100)
!    PRINT *,"result",result(:,N)

    result = test()

        contains
            real ( kind = 8 ) function test()
                implicit none

                real ( kind = 8 ), dimension(n)     :: x, y, diff, vxp, vxp2, vxp3, vyp, vyp2, vyp3, temp
                complex ( kind = 8 ), dimension(n)  :: w_hat, ifft, v2, v_hat
                real ( kind = 8 ), dimension(n,n)   :: v
                real ( kind = 8 ), dimension(n,n)   :: D
                integer                             :: j,i,plan_backward, plan_forward
                real ( kind = 8 )                   :: pi, h

                integer ( kind = 4 ), parameter :: fftw_forward = -1
                integer ( kind = 4 ), parameter :: fftw_backward = +1
                integer ( kind = 4 ), parameter :: fftw_estimate = 64

                pi = 4*atan(1.)
                h=2*pi/n;
                x=h*(/ (j,j=1,N) /);
                y=x;

                d=0
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

                v = 0;

                diff=0;
                do i=1,n
                    v(:,i) = cos(x)*exp(sin(y(i)));

                    vxp= -sin(x)*exp(sin(y(i)));
                    vxp2 = MATMUL(D,v(:,i));

                    v_hat=0
                    v2 = v(:,i)
                    call dfftw_plan_dft_1d_ ( plan_forward, N, v2, v_hat, FFTW_FORWARD, FFTW_ESTIMATE )
                    call dfftw_execute_ ( plan_forward )

                    w_hat = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat

                    call dfftw_plan_dft_1d_ ( plan_backward, N, w_hat, ifft, FFTW_BACKWARD, FFTW_ESTIMATE )
                    call dfftw_execute_ ( plan_backward )
                    ifft = ifft/n

                    vxp3 = real(ifft)

                    diff = MAX(MAXVAL(diff), MAXVAL(abs(vxp - vxp2)), MAXVAL(abs(vxp3 - vxp)))
                end do

                do i=1,n

                    vyp= cos(x(i))*exp(sin(y))*cos(y);

                    vyp2 = MATMUL(D,v(i,:));

                    v_hat=0
                    v2 = v(i,:)
                    call dfftw_plan_dft_1d_ ( plan_forward, N, v2, v_hat, FFTW_FORWARD, FFTW_ESTIMATE )
                    call dfftw_execute_ ( plan_forward )

                    w_hat = (0., 1.) * (/ (j,j=0,N/2-1), 0, (j,j=-N/2+1,-1) /) * v_hat

                    call dfftw_plan_dft_1d_ ( plan_backward, N, w_hat, ifft, FFTW_BACKWARD, FFTW_ESTIMATE )
                    call dfftw_execute_ ( plan_backward )
                    ifft = ifft/n

                    vyp3 = real(ifft)

                    diff = MAX(MAXVAL(diff), MAXVAL(abs(vyp - vyp2)), MAXVAL(abs(vyp3 - vyp)));
                end do

                PRINT *,"linf diff", MAXVAL(abs(diff))

                test = MAXVAL(abs(diff))

            end function test

end program main
