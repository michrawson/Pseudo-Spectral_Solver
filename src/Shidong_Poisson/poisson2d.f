c
c
c
c
        subroutine poisson2d(nx,ny,hx,hy,rn1x,rn1y,rn2x,rn2y,
     &    f0,uc,iprec,ifder,iffrac)
c
c       input parameters
c       
c       nx, ny: number of points in x, y directions
c       hx, hy: spatial step size in x, y directions
c       rn1x, rn1y, rn2x, rn2y: two unit normal vectors specifying two dipole orientations
c       f0: the right hand side over the uniform grid in R^2
c       iprec: precision control, 0 gives 6 digit  accuracy, 1 gives 12 digit accuracy
c       ifder: 0 - computes the solution and rn1x,rn1y,rn2x,rn2y are NOT used
c              1 - computes the dipolar interactions uc_{n_1 n_2}
c
c       iffrac: 0 - computes the solution to the Poisson equation -\Delta uc = f0
c               1 - computes the solution to the fractional Poisson equation \sqrt(-\Dealta)uc = f0
c
c       output parameters
c
c       uc : the computed solution on the same spatial grid as f0
c
c       WARNING: only works even nx and ny now, ifder=0, iffrac=0 case not implemented!
c
        implicit real *8 (a-h,o-z)
        complex *16 uc(nx,ny),f0(nx,ny)
        complex *16, allocatable :: u2(:,:)

        dk=min(1.0d0/hx/nx,1.0d0/hy/ny)

c        R0=0.1d0
c        R1=0.6d0

        R0=0.8d0*dk
        R1=10*dk
        call prin2('R1=*',R1,1)

        if (iprec .eq. 1) then
           kfac=3
        else
           kfac=2
        endif

        call ufftpart(hx,hy,nx,ny,rn1x,rn1y,rn2x,rn2y,kfac,
     &      R0,R1,f0,uc,ifder,iffrac)

        if (iprec .eq. 1) then
           nr=60
           nphi=60
        else
           nr=40
           nphi=40
        endif

c	nr=50
c	nphi=50

        allocate( u2(nx,ny) )
        call uradpart(nr,nphi,hx,hy,nx,ny,rn1x,rn1y,rn2x,rn2y,
     &      R0,R1,f0,u2,ifder,iffrac)

        uc=uc+u2

        if (ifder .eq. 1) then
           uc=-uc
        endif

        if (ifder .eq. 0) then
           d=1/(2*pi)**2
           uc =uc*d
        endif

        return
        end
c
c
c
c
        subroutine fft2(nx,ny,fval,fhat)
        implicit real *8 (a-h,o-z)
        include "fftw3.f"
        integer *8 plan
        complex *16 fval(nx,ny),fhat(nx,ny)

        t1=second()
        call dfftw_plan_dft_2d(plan, nx,ny, fval,fhat,
     &      FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, fval, fhat)
        call dfftw_destroy_plan(plan)
        t2=second()

        call prin2('after forward fft, time (sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine ifft2(nx,ny,fhat,fval)
        implicit real *8 (a-h,o-z)
        include "fftw3.f"
        integer *8 plan
        complex *16 fval(nx,ny),fhat(nx,ny)

        t1=second()
        call dfftw_plan_dft_2d(plan, nx,ny, fhat,fval,
     &      FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, fhat, fval)
        call dfftw_destroy_plan(plan)

        nt=nx*ny
        fval=fval/nt

        t2=second()

        call prin2('after backward fft, time (sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine ufftpart(hx,hy,nx,ny,rn1x,rn1y,rn2x,rn2y,kfac,
     &    R0,R1,f0,uc,ifder,iffrac)
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: rk1(:),rk2(:)
        parameter (pi=3.141592653589793d0)
        complex *16 uc(nx,ny),f0(nx,ny)
        complex *16, allocatable :: f1(:,:),fhat(:,:)
	complex *16, allocatable :: f2(:,:),fhat2(:,:)


        n1=nx*kfac
        n2=ny*kfac

        allocate( f1(n1,n2) )
        allocate( f2(n1,n2) )
        allocate( fhat(n1,n2) )
        allocate( fhat2(n1,n2) )

        call zeropad(nx,ny,kfac,f0,f1)
        call ifftshift2(n1,n2,f1,f2)
        call fft2(n1,n2,f2,fhat2)
        call ifftshift2(n1,n2,fhat2,fhat)

        allocate( rk1(n1) )
        allocate( rk2(n2) )

        rk1=(/ (k, k=-n1/2,n1/2-1) /)
        rk2=(/ (k, k=-n2/2,n2/2-1) /)

        rk1=rk1/hx/n1
        rk2=rk2/hy/n2

        if ((ifder .eq. 1) .and. (iffrac .eq. 1)) then
        do j=1,n2
           do i=1,n1
              r2=rk1(i)**2+rk2(j)**2
              rkn=(rk1(i)*rn1x+rk2(j)*rn1y)*(rk1(i)*rn2x+rk2(j)*rn2y)
              r=sqrt(r2)
              call erfcbase(r,R0,R1,v2)
              fhat(i,j)=fhat(i,j)*rkn/r*(1-v2)
           enddo
        enddo
        
        elseif ((ifder .eq. 1) .and. (iffrac .eq. 0)) then

        do j=1,n2
           do i=1,n1
              r2=rk1(i)**2+rk2(j)**2
              rkn=(rk1(i)*rn1x+rk2(j)*rn1y)*(rk1(i)*rn2x+rk2(j)*rn2y)
              r=sqrt(r2)
              call erfcbase(r,R0,R1,v2)
              fhat(i,j)=fhat(i,j)*rkn/r2*(1-v2)
           enddo
        enddo

        elseif ((ifder .eq. 0) .and. (iffrac .eq. 1)) then
        
        do j=1,n2
           do i=1,n1
              r2=rk1(i)**2+rk2(j)**2
              r=sqrt(r2)
              call erfcbase(r,R0,R1,v2)
              fhat(i,j)=fhat(i,j)/r*(1-v2)
           enddo
        enddo

        endif

        fhat(n1/2+1,n2/2+1)=0

        call ifftshift2(n1,n2,fhat,fhat2)
        call ifft2(n1,n2,fhat2,fhat)
        call ifftshift2(n1,n2,fhat,fhat2)

        ncx=nx*(kfac-1)/2
        ncy=ny*(kfac-1)/2

        uc=fhat2((/(ncx+i, i=1,nx)/),(/(ncy+i, i=1,ny)/))
        
        return
        end
c
c
c
c
        subroutine ifftshift2(n1,n2,fin,fout)
        implicit real *8 (a-h,o-z)
        complex *16 fin(n1,n2), fout(n1,n2)
        
        m1=n1/2
        m2=n2/2

        fout=fin((/ (/(i,i=m1+1,n1)/), (/(i,i=1,m1)/)/),
     &      (/ (/(i,i=m2+1,n2)/), (/(i,i=1,m2)/)/))

        return
        end
c
c
c
c
        subroutine uradpart(nr,nphi,hx,hy,nx,ny,rn1x,rn1y,rn2x,rn2y,
     &    R0,R1,f0,u2,ifder,iffrac)
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: r(:),phi(:),wr(:),wphi(:)
        real *8, allocatable :: rk1(:,:),rk2(:,:)
	real *8, allocatable :: x1(:),y1(:)
        parameter (pi=3.141592653589793d0)
        complex *16 u2(nx,ny),f0(nx,ny)
        complex *16, allocatable :: fhat(:,:),f1(:)

        allocate( r(nr) )
        allocate( wr(nr) )

        allocate( phi(nphi) )
        allocate( wphi(nphi) )

	nt=nphi*nr
	allocate(x1(nt))
	allocate(y1(nt))

	allocate(f1(nt))


        allocate( rk1(nphi,nr) )
        allocate( rk2(nphi,nr) )
        allocate( fhat(nphi,nr) )

        call rnodes2(nr,nphi,R1,r,wr,phi,wphi)
        call fhateval(r,phi,nr,nphi,nx,ny,f0,hx,hy,fhat)

        cx=2*pi*hx
        cy=2*pi*hy

        if ((ifder .eq. 1) .and. (iffrac .eq. 1)) then

        do j=1,nr
           call erfcbase(r(j),R0,R1,v2)
           do i=1,nphi
              rk1(i,j)=r(j)*cos(phi(i))
              rk2(i,j)=r(j)*sin(phi(i))
              v3=(cos(phi(i))*rn1x+sin(phi(i))*rn1y)*
     &            (cos(phi(i))*rn2x+sin(phi(i))*rn2y)*r(j)**2
              fhat(i,j)=fhat(i,j)*wphi(i)*wr(j)*v2*v3
           enddo
        enddo

        elseif ((ifder .eq. 1) .and. (iffrac .eq. 0)) then

        do j=1,nr
           call erfcbase(r(j),R0,R1,v2)
           do i=1,nphi
              rk1(i,j)=r(j)*cos(phi(i))
              rk2(i,j)=r(j)*sin(phi(i))
              v3=(cos(phi(i))*rn1x+sin(phi(i))*rn1y)*
     &            (cos(phi(i))*rn2x+sin(phi(i))*rn2y)*r(j)
              fhat(i,j)=fhat(i,j)*wphi(i)*wr(j)*v2*v3
           enddo
        enddo

        elseif ((ifder .eq. 0) .and. (iffrac .eq. 1)) then
        do j=1,nr
           call erfcbase(r(j),R0,R1,v2)
           do i=1,nphi
              rk1(i,j)=r(j)*cos(phi(i))
              rk2(i,j)=r(j)*sin(phi(i))
              fhat(i,j)=fhat(i,j)*wphi(i)*wr(j)*v2
           enddo
        enddo

        endif

        rk1=cx*rk1
        rk2=cy*rk2


	x1=reshape(rk1,(/nt/))
	y1=reshape(rk2,(/nt/))

	f1=reshape(fhat,(/nt/))

        iflag=1
        eps=1e-12
	t1=second()
        call nufft2d1f90(nt,rk1,rk2,fhat,iflag,eps,nx,ny,u2,ier)
	t2=second()

        u2=u2*nt

	call prin2('after nufft2d1, time(sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine rnodes2(nr,nphi,R1,r,wr,phi,wphi)
        implicit real *8 (a-h,o-z)
        dimension r(nr),phi(nphi),wr(nr),wphi(nphi)
        parameter (pi=3.141592653589793d0)

        itype=1
        call legeexps(itype,nr,r,u,v,wr)

        r2=R1/2
        r=r2+r2*r
        wr=r2*wr

        h=2*pi/nphi
        phi=(/ (i-1, i=1,nphi) /)
        phi=h*phi

        wphi=(/ (h, i=1,nphi) /)

        return
        end
c
c
c
c
        subroutine fhateval(r,phi,nr,nphi,nx,ny,fx,hx,hy,fhat)
        implicit real *8 (a-h,o-z)
        dimension r(nr),phi(nphi)
        real *8, allocatable :: rk1(:,:),rk2(:,:)
	real *8, allocatable :: x1(:),y1(:)
        parameter (pi=3.141592653589793d0)
        complex *16 fhat(nphi,nr),fx(nx,ny)
	complex *16, allocatable :: f1(:)

        cx=2*pi*hx
        cy=2*pi*hy

        nt=nr*nphi
	allocate(x1(nt))
	allocate(y1(nt))
	allocate(f1(nt))

        allocate( rk1(nphi,nr) )
        allocate( rk2(nphi,nr) )

        do j=1,nr
           do i=1,nphi
              rk1(i,j)=cx*r(j)*cos(phi(i))
              rk2(i,j)=cy*r(j)*sin(phi(i))
           enddo
        enddo

	x1=reshape(rk1,(/nt/))
	y1=reshape(rk2,(/nt/))

        iflag=-1
        eps=1e-12
	t1=second()
        call nufft2d2f90(nt,x1,y1,f1,iflag,eps,nx,ny,fx,ier)
        t2=second()

	fhat=reshape(f1,(/nphi,nr/))
        fhat=fhat*hx*hy

	call prin2('after nufft2d2, time (sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine zeropad(nx,ny,nfac,fin,fout)
        implicit real *8 (a-h,o-z)
        complex *16 fin(nx,ny),fout(nx*nfac,ny*nfac)

        fout=0

        ncx=nx*(nfac-1)/2;
        ncy=ny*(nfac-1)/2;

        fout((/(ncx+i,i=1,nx)/),(/(ncy+i,i=1,ny)/))=fin

        return
        end
c
c
c
c
        subroutine erfcbase(r,t0,t1,f)
        implicit real *8 (a-h,o-z)
        
        x=(r-(t0+t1)/2)/(t1-t0)*12
        call errfunc(x,f)
        f=f/2

        return
        end
c
c
c
c
