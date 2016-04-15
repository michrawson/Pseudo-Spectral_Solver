
program main
    use variable_coeff_wave_eq_pseudo_rk

    implicit none

    integer, parameter           :: n = 128
    integer, parameter           :: tmax = 8
    real ( kind = 8 ), parameter :: tplot = 0.15

    real ( kind = 8 ), dimension(int(tmax/tplot)+1)      :: tdata
    real ( kind = 8 ), dimension(int(tmax/tplot)+1,n,n)  :: result

    real ( kind = 8 )                 :: x(n), y(n)

    call variable_coeff_wave_eq_pseudo_rk_run(x,y,tdata,result)
    PRINT *,"x",x(1:10)
    PRINT *,"tdata",tdata(1:10)
    PRINT *,"result",result(1,1:10,1:10)

!    nx=64
!    ny=64
!    hx=1
!    hy=1
!
!    allocate( ux(nx,ny) )
!    allocate( uy(nx,ny) )
!    allocate( uxt(nx,ny) )
!    allocate( uyt(nx,ny) )
!    allocate( f0(nx,ny) )
!
!    f0=0
!    ux=0
!    uy=0
!
!    call poisson2df(nx,ny,hx,hy,f0,ux,uy,1)

end program main
