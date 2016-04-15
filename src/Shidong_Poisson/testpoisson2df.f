 

        implicit real *8 (a-h,o-z)
        real *8, allocatable :: x(:),y(:)
        complex *16, allocatable :: ux(:,:),f0(:,:),uxt(:,:)
        complex *16, allocatable :: uy(:,:),uyt(:,:)
        complex *16 ima,d,dn,df
        data ima/(0.0d0,1.0d0)/
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
c
        print *, 'Enter nx'
        read *, nx
        call prinf('nx=*',nx,1)

c        print *, 'Enter iprec'
c        read *, iprec
c        call prinf('iprec=*',iprec,1)

c        nx=2**6
c        ny=2**6
        ny=nx

c       iprec=0 - 6 digits; iprec=1 - 12 digits
        iprec=1

        D=6
        a1=1.0d0
        a1=1.3d0
        
        call testfunini(a1)
        call truesolini(a1)

        hx=2*D*sqrt(a1)/nx
        hy=hx

        allocate( x(nx) )
        allocate( y(ny) )

        x=(/ (-nx/2+i, i=0,nx-1)/)
        y=(/ (-ny/2+i, i=0,ny-1)/)

        x=x*hx
        y=y*hy

        allocate( ux(nx,ny) )
        allocate( uy(nx,ny) )
        allocate( uxt(nx,ny) )
        allocate( uyt(nx,ny) )
        allocate( f0(nx,ny) )
c
c       compute the right hand side function on the uniform grid
c
        call testfun(x,y,nx,ny,f0)

c
c       compute the numerical solution
c     
        call poisson2df(nx,ny,hx,hy,f0,ux,uy,iprec)
        call prin2('ux=*',dble(ux),60)
        call prin2('uy=*',dble(uy),60)

c
c       compute the analytical solution
c
        call truesol(x,y,nx,ny,uxt,uyt)
        call prin2('uxt=*',dble(uxt),60)
        call prin2('uyt=*',dble(uyt),60)

c
c       compute the relative l2 error
c
        d=0
        dn=0
        do j=1,ny
           do i=1,nx
              df=uxt(i,j)-ux(i,j)
              d=d+df*conjg(df)
              dn=dn+uxt(i,j)*conjg(uxt(i,j))
           enddo
        enddo

        d=sqrt(dble(d))
        dn=sqrt(dble(dn))

        call prin2('relative l2 error=*',d/dn,1)

        stop
        end
c
c
c
c
        subroutine truesol(x,y,nx,ny,ux,uy)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny)
        real *8, allocatable :: x2(:),y2(:)
        complex *16 ux(nx,ny),uy(nx,ny),u
        save

        allocate( x2(nx) )
        allocate( y2(ny) )

        x2=x*x
        y2=y*y

        do j=1,ny
           do i=1,nx
!              r=sqrt(x2(i)+y2(j))
!              call uspec(x(i),y(j),a,u,ux(i,j),uy(i,j))
              ux(i,j) = -2.0*x(i)*exp(-(x2(i)+y2(j)))
              uy(i,j) = -2.0*y(j)*exp(-(x2(i)+y2(j)))
           enddo
        enddo

	deallocate(x2)
	deallocate(y2)

        return
        entry truesolini(a7)

        a=a7

        return
        end
c
c
c
c
        subroutine uspec(x,y,a,u,ux,uy)
c
c       returns the solution to -\Delta u = exp(-r^2/a) in 2D.
c
        implicit real *8 (a-h,o-z)
        data gamma/0.5772156649015328606065120900824024310421593d0/
        complex *16 u,ux,uy,du
c
        rho=sqrt(x**2+y**2)
        r=rho**2/a

        if (r .gt. 1d-3) then
           u=-a/4*(eone(r)+2*log(rho))
           du=-(1-exp(-r))/r/2
           ux=du*x
           uy=du*y
        else
           r2=r*r
           r3=r2*r
           r4=r3*r
           r5=r4*r
           r6=r5*r

           u=-a/4*(-gamma+r-r2/4+r3/18-r4/96+r5/600-r6/4320+log(a))
           du=-(1-r/2+r2/6-r3/24+r4/120-r5/720)/2
           ux=du*x
           uy=du*y
        endif

        return
        end
c
c
c
c
        subroutine testfun(x,y,nx,ny,f)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny)
        real *8, allocatable :: x2(:),y2(:)
        complex *16 f(nx,ny)
        save
 
        allocate( x2(nx) )
        allocate( y2(ny) )

        x2=x*x
        y2=y*y

        do j=1,ny
           do i=1,nx
         f(i,j)= 2.0*exp(-(x2(i)+y2(j))) - 4.0*x2(i)*exp(-(x2(i)+y2(j)))
     &         + 2.0*exp(-(x2(i)+y2(j))) - 4.0*y2(j)*exp(-(x2(i)+y2(j)))
!                f(i,j) = cos(x(i)) + cos(y(j))
           enddo
        enddo

	deallocate(x2)
	deallocate(y2)

        return
        entry testfunini(a7)

        a=a7

        return
        end
c
c
c
c
