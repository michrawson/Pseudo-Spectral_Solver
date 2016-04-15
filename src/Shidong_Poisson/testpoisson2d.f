 

        implicit real *8 (a-h,o-z)
        real *8, allocatable :: x(:),y(:)
        complex *16, allocatable :: uc(:,:),f0(:,:),ut(:,:)
        complex *16 ima,d,dn,df
        integer *8 count
        integer *4 t(2),s,pid,iseed(12)
        dimension dran(4)
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

c        rn1x=hkrand(0)-0.5
c        rn1y=hkrand(0)-0.5

        rn1x=dran(1)-0.5
        rn1y=dran(2)-0.5

        rn1=sqrt(rn1x**2+rn1y**2)

        rn1x=rn1x/rn1
        rn1y=rn1y/rn1

c        rn2x=hkrand(0)-0.5
c        rn2y=hkrand(0)-0.5

        rn2x=dran(3)-0.5
        rn2y=dran(4)-0.5

        rn2=sqrt(rn2x**2+rn2y**2)

        rn2x=rn2x/rn2
        rn2y=rn2y/rn2


c        rn2x=rn1x
c        rn2y=rn1y

	call prin2('rn1x=*',rn1x,1)
	call prin2('rn1y=*',rn1y,1)

	call prin2('rn2x=*',rn2x,1)
	call prin2('rn2y=*',rn2y,1)

c       iprec=0 - 6 digits; iprec=1 - 12 digits
        iprec=1

        ifder=1
        iffrac=0


        D=6
        a1=1.0d0

        a2=1.5d0

        call testfunini(a1,iffrac)
        call truesolini(a1,ifder,iffrac)

        call testfun0ini(a1,a2)
        call truesol0ini(a1,a2)

        hx=2*D*sqrt(a1)/nx

        hy=hx
        if (iffrac .eq. 0) hy=2*D*sqrt(a2)/ny

        allocate( x(nx) )
        allocate( y(ny) )

        x=(/ (-nx/2+i, i=0,nx-1)/)
        y=(/ (-ny/2+i, i=0,ny-1)/)

        x=x*hx
        y=y*hy


        allocate( uc(nx,ny) )
        allocate( ut(nx,ny) )
        allocate( f0(nx,ny) )
c
c       compute the right hand side function on the uniform grid
c
        call testfun(x,y,nx,ny,f0)

c
c       compute the numerical solution
c     
        call poisson2d(nx,ny,hx,hy,rn1x,rn1y,rn2x,rn2y,
     &      f0,uc,iprec,ifder,iffrac)
        call prin2('uc=*',dble(uc),60)

c
c       compute the analytical solution
c
        call truesol(x,y,nx,ny,rn1x,rn1y,rn2x,rn2y,ut)
        call prin2('ut=*',dble(ut),60)
c
c       compute the relative l2 error
c
        d=0
        dn=0
        do j=1,ny
           do i=1,nx
              df=ut(i,j)-uc(i,j)
              d=d+df*conjg(df)
              dn=dn+ut(i,j)*conjg(ut(i,j))
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
        subroutine truesol(x,y,nx,ny,rn1x,rn1y,rn2x,rn2y,unn)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny)
        real *8, allocatable :: x2(:),y2(:)
        complex *16 unn(nx,ny)
        data pi/3.141592653589793238462643383279502884197169399d0/
        save


        if ((ifder .eq. 1) .and. (iffrac .eq. 0)) then
           call truesol0(x,y,nx,ny,rn1x,rn1y,rn2x,rn2y,unn)
           return
        endif
           
        allocate( x2(nx) )
        allocate( y2(ny) )
        
        x2=x*x
        y2=y*y

        if ((ifder .eq. 1) .and. (iffrac .eq. 1)) then
        do j=1,ny
           do i=1,nx
              r2=x2(i)+y2(j)
              rho=r2/a/2
              di0=besei0(rho)
              di1=besei1(rho)

              dn1n2=rn1x*rn2x+rn1y*rn2y
              rn1=x(i)*rn1x+y(j)*rn1y
              rn2=x(i)*rn2x+y(j)*rn2y
              

              unn(i,j)=c*((di1-di0)*dn1n2+rn1*rn2/a
     &            *(2*di0-2*di1-di1/rho))
           enddo
        enddo

        unn(nx/2+1,ny/2+1)=-c*dn1n2
        
        elseif ((ifder .eq. 0) .and. (iffrac .eq. 1)) then
        do j=1,ny
           do i=1,nx
              r2=x2(i)+y2(j)
              rho=r2/a/2
              di0=besei0(rho)

              unn(i,j)=c1*di0

           enddo
        enddo

        endif

        do j=1,ny
           do i=1,nx
              unn(i,j)=exp(-(x2(i) + y2(j)))
           enddo
        enddo


	deallocate(x2)
	deallocate(y2)

        return
        entry truesolini(a7,ifder7,iffrac7)

        a=a7

        ifder=ifder7

        iffrac=iffrac7

        c=1.0d0/4/sqrt(pi*a)

        c1=sqrt(a/pi)/4

        return
        end
c
c
c
c       
        subroutine truesol0(x,y,nx,ny,rn1x,rn1y,rn2x,rn2y,unn)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny)
        real *8, allocatable :: x2(:),y2(:)
        complex *16 unn(nx,ny)
        save

        allocate( x2(nx) )
        allocate( y2(ny) )
        
        x2=x*x
        y2=y*y

        dn1n2=rn1x*rn2x/ax+rn1y*rn2y/ay
        do j=1,ny
           do i=1,nx
              rn1=x(i)*rn1x/ax+y(j)*rn1y/ay
              rn2=x(i)*rn2x/ax+y(j)*rn2y/ay
              unn(i,j)=-(4*rn1*rn2-2*dn1n2)*exp(-(x2(i)/ax+y2(j)/ay))
           enddo
        enddo

        do j=1,ny
           do i=1,nx
              unn(i,j)=exp(-(x2(i)+y2(j)))
           enddo
        enddo

	deallocate(x2)
	deallocate(y2)

        return
        entry truesol0ini(ax7,ay7)

        ax=ax7
        ay=ay7

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
        
        if (iffrac .eq. 0) then
           call testfun0(x,y,nx,ny,f)
           return
        endif

        allocate( x2(nx) )
        allocate( y2(ny) )

        x2=x*x
        y2=y*y

        do j=1,ny
           do i=1,nx
!              f(i,j)=exp(-(x2(i)+y2(j))/a)
              f(i,j)=exp(-(x2(i)+y2(j))) + exp(-(x2(i)+y2(j)))
!              f(i,j)= cos(x(i)) + cos(y(j))
           enddo
        enddo

	deallocate(x2)
	deallocate(y2)

        return
        entry testfunini(a7,iffrac7)

        iffrac=iffrac7

        a=a7

        return
        end
c
c
c
        subroutine testfun0(x,y,nx,ny,f)
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny)
        real *8, allocatable :: x2(:),y2(:)
        complex *16 f(nx,ny)
        save

        allocate( x2(nx) )
        allocate( y2(ny) )

        x2=x*x
        y2=y*y

        ax2=ax*ax
        ay2=ay*ay

        rtmp=-2.0d0/ax-2.0d0/ay

        do j=1,ny
           do i=1,nx
!              f(i,j)=(rtmp+4*(x2(i)/ax2+y2(j)/ay2))
!     &            *exp(-(x2(i)/ax+y2(j)/ay))
              f(i,j)=exp(-(x2(i)+y2(j))) + exp(-(x2(i)+y2(j)))
           enddo
        enddo

	deallocate(x2)
	deallocate(y2)

        return
        entry testfun0ini(ax7,ay7)

        ax=ax7
        ay=ay7

        return
        end
c
c
c
