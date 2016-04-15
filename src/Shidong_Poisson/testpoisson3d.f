 

        implicit real *8 (a-h,o-z)
        real *8, allocatable :: x(:),y(:),z(:)
        complex *16, allocatable :: uc(:,:,:),f0(:,:,:),ut(:,:,:)
        complex *16 ima,d,dn,df
        integer *8 count
        integer *4 t(2),s,pid,iseed(12)
        dimension dran(10)
        data ima/(0.0d0,1.0d0)/
c
c       SET ALL PARAMETERS
c        
        call prini(6,13)
c
        print *, 'Enter nx'
        read *, nx
        call prinf('nx=*',nx,1)

        ny=nx
        nz=nx
c
c       set random unit normal vectors
c
        call system_clock(count)
c        print *, 'count=', count
        t = transfer(count, t)
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 
        s = ieor(s, pid)
        iseed(1) = t(1) + 36269
        iseed(2) = t(2) + 72551
        iseed(3) = pid
        iseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        call random_seed(put=iseed)
        call random_number(dran)

        rn1x=dran(1)-0.5
        rn1y=dran(2)-0.5
        rn1z=dran(3)-0.5

        rn1=sqrt(rn1x**2+rn1y**2+rn1z**2)

        rn1x=rn1x/rn1
        rn1y=rn1y/rn1
        rn1z=rn1z/rn1

        rn2x=dran(4)-0.5
        rn2y=dran(5)-0.5
        rn2z=dran(6)-0.5

        rn2=sqrt(rn2x**2+rn2y**2+rn2z**2)

        rn2x=rn2x/rn2
        rn2y=rn2y/rn2
        rn2z=rn2z/rn2

	call prin2('rn1x=*',rn1x,1)
	call prin2('rn1y=*',rn1y,1)
	call prin2('rn1z=*',rn1z,1)

	call prin2('rn2x=*',rn2x,1)
	call prin2('rn2y=*',rn2y,1)
	call prin2('rn2z=*',rn2z,1)

c       iprec=0 - 6 digits; iprec=1 - 12 digits
        iprec=1

        ifder=1

        D=6
        a1=1.0d0
        a2=1.0d0
        a3=1.0d0

        call truesolini(a1,a2,a3,ifder)
        call testfunini(a1,a2,a3)

        hx=2*D*sqrt(a1)/nx
        hy=2*D*sqrt(a2)/ny
        hz=2*D*sqrt(a3)/nz

        allocate( x(nx) )
        allocate( y(ny) )
        allocate( z(nz) )

        x=(/ (-nx/2+i, i=0,nx-1)/)
        y=(/ (-ny/2+i, i=0,ny-1)/)
        z=(/ (-nz/2+i, i=0,nz-1)/)

        x=x*hx
        y=y*hy
        z=z*hz

        allocate( uc(nx,ny,nz) )
        allocate( f0(nx,ny,nz) )
        allocate( ut(nx,ny,nz) )
c
c       compute the right hand side function on the uniform grid
c
        call testfun(x,y,z,nx,ny,nz,f0)
c
c       compute the numerical solution
c     
        call poisson3d(nx,ny,nz,hx,hy,hz,rn1x,rn1y,rn1z,
     &    rn2x,rn2y,rn2z,f0,uc,iprec,ifder)

        call truesol(x,y,z,nx,ny,nz,rn1x,rn1y,rn1z,rn2x,rn2y,rn2z,ut)
c
c       compute the relative l2 error
c

        d=0
        dn=0
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 df=ut(i,j,k)-uc(i,j,k)
                 d=d+df*conjg(df)
                 dn=dn+ut(i,j,k)*conjg(ut(i,j,k))
              enddo
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
        subroutine truesol(x,y,z,nx,ny,nz,rn1x,rn1y,rn1z,
     &    rn2x,rn2y,rn2z,unn)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny),z(nz)
        real *8, allocatable :: x2(:),y2(:),z2(:)
        complex *16 unn(nx,ny,nz)
        save

        allocate( x2(nx) )
        allocate( y2(ny) )
        allocate( z2(nz) )
        
        x2=x*x
        y2=y*y
        z2=z*z

        if (ifder .eq. 0) then

        do k=1,nz
           do j=1,ny
              do i=1,nx
                 unn(i,j,k)=-exp(-(x2(i)/ax+y2(j)/ay+z2(k)/az))
              enddo
           enddo
        enddo

        elseif (ifder .eq. 1) then
        
        dn1n2=rn1x*rn2x/ax+rn1y*rn2y/ay+rn1z*rn2z/az
        do k=1,nz
           do j=1,ny
              do i=1,nx
                 rn1=x(i)*rn1x/ax+y(j)*rn1y/ay+z(k)*rn1z/az
                 rn2=x(i)*rn2x/ax+y(j)*rn2y/ay+z(k)*rn2z/az

                 unn(i,j,k)=-(4*rn1*rn2-2*dn1n2)
     &               *exp(-(x2(i)/ax+y2(j)/ay+z2(k)/az))
              enddo
           enddo
        enddo

        endif

	deallocate(x2)
	deallocate(y2)
	deallocate(z2)

        return
        entry truesolini(ax7,ay7,az7,ifder7)

        ax=ax7
        ay=ay7
        az=az7

        ifder=ifder7

        return
        end
c
c
c
c       
        subroutine testfun(x,y,z,nx,ny,nz,f)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny),z(nz)
        real *8, allocatable :: x2(:),y2(:),z2(:)
        complex *16 f(nx,ny,nz)
        save

        allocate( x2(nx) )
        allocate( y2(ny) )
        allocate( z2(nz) )

        x2=x*x
        y2=y*y
        z2=z*z

        ax2=ax*ax
        ay2=ay*ay
        az2=az*az

        rtmp=-2.0d0/ax-2.0d0/ay-2.0d0/az

        do k=1,nz
           do j=1,ny
              do i=1,nx
                 f(i,j,k)=(rtmp+4*(x2(i)/ax2+y2(j)/ay2+z2(k)/az2))
     &               *exp(-(x2(i)/ax+y2(j)/ay+z2(k)/az))
              enddo
           enddo
        enddo

	deallocate(x2)
	deallocate(y2)
	deallocate(z2)

        return
        entry testfunini(ax7,ay7,az7)

        ax=ax7
        ay=ay7
        az=az7

        return
        end
c
c
c
