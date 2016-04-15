c
c
c
c
        subroutine poisson3d(nx,ny,nz,hx,hy,hz,rn1x,rn1y,rn1z,
     &    rn2x,rn2y,rn2z,f0,uc,iprec,ifder)
c
c       input parameters
c       
c       nx, ny, nz: number of points in x, y, and z directions
c       hx, hy, hz: spatial step size in x, y, and z directions
c       rn1x, rn1y, rn1z, rn2x, rn2y, rn2z: two unit normal vectors 
c                        specifying two dipole orientations
c       f0: the right hand side over the equispace grid in R^2
c       iprec: precision control, 0 gives 6 digit  accuracy, 1 gives 12 digit accuracy
c       ifder: 0 - computes the solution and rn1x,rn1y,rn1z,rn2x,rn2y,rn2z are NOT used
c              1 - computes the dipolar interactions
c
c       output parameters
c
c       uc : the computed solution on the same spatial grid as f0
c
c       WARNING: only works even nx, ny, and nz now!
c

        implicit real *8 (a-h,o-z)
        complex *16 uc(nx,ny,nz),f0(nx,ny,nz)
        complex *16, allocatable :: u2(:,:,:)
        data pi/3.141592653589793238462643383279502884197169399d0/

c        R1=1d-2
c        dkmax=2.0d0/hx
        dk=min(1.0d0/hx/nx,1.0d0/hy/ny,1.0d0/hz/nz)
        call prin2('dk=*',dk,1)

        R1=2*dk**2
        call prin2('R1=*',R1,1)
        kfac=2
c        R1=1d-3
        call ufftpart(hx,hy,hz,nx,ny,nz,rn1x,rn1y,rn1z,
     &      rn2x,rn2y,rn2z,kfac,R1,f0,uc,ifder) 

        if (iprec .eq. 1) then
           nr=40
           nzt=40
           nphi=40
        else
           nr=24
           nzt=24
           nphi=24
        endif
c	nr=30
c	nzt=30
c	nphi=30

        allocate( u2(nx,ny,nz) )
        call uradpart(nr,nzt,nphi,hx,hy,hz,x,y,z,nx,ny,nz,
     &    rn1x,rn1y,rn1z,rn2x,rn2y,rn2z,R1,f0,u2,ifder) 

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
        subroutine fft3(nx,ny,nz,fval,fhat)
        implicit real *8 (a-h,o-z)
        include "fftw3.f"
        integer *8 plan
        complex *16 fval(nx,ny,nz),fhat(nx,ny,nz)

        t1=second()
        call dfftw_plan_dft_3d(plan, nx,ny,nz, fval,fhat,
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
        subroutine ifft3(nx,ny,nz,fhat,fval)
        implicit real *8 (a-h,o-z)
        include "fftw3.f"
        integer *8 plan
        complex *16 fval(nx,ny,nz),fhat(nx,ny,nz)

        t1=second()
        call dfftw_plan_dft_3d(plan, nx,ny,nz, fhat,fval,
     &      FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, fhat, fval)
        call dfftw_destroy_plan(plan)

        nt=nx*ny*nz
        fval=fval/nt

        t2=second()

        call prin2('after backward fft, time (sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine ufftpart(hx,hy,hz,nx,ny,nz,rn1x,rn1y,rn1z,
     &      rn2x,rn2y,rn2z,kfac,R1,f0,uc,ifder) 
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: x(:),y(:),z(:)
        real *8, allocatable :: rk1(:),rk2(:),rk3(:)
        parameter (pi=3.141592653589793d0)
        complex *16 uc(nx,ny,nz),f0(nx,ny,nz)
        complex *16, allocatable :: f1(:,:,:),fhat(:,:,:)
	complex *16, allocatable :: f2(:,:,:),fhat2(:,:,:)


        n1=nx*kfac
        n2=ny*kfac
        n3=nz*kfac

        allocate( f1(n1,n2,n3) )
        allocate( f2(n1,n2,n3) )
        allocate( fhat(n1,n2,n3) )
        allocate( fhat2(n1,n2,n3) )

        call zeropad(nx,ny,nz,kfac,f0,f1)
        call ifftshift3(n1,n2,n3,f1,f2)
        call fft3(n1,n2,n3,f2,fhat2)
        call ifftshift3(n1,n2,n3,fhat2,fhat)


        allocate( rk1(n1) )
        allocate( rk2(n2) )
        allocate( rk3(n3) )

        rk1=(/ (k, k=-n1/2,n1/2-1) /)
        rk2=(/ (k, k=-n2/2,n2/2-1) /)
        rk3=(/ (k, k=-n3/2,n3/2-1) /)

        rk1=rk1/hx/n1
        rk2=rk2/hy/n2
        rk3=rk3/hz/n3

        if (ifder .eq. 0) then

        do k=1,n3
           do j=1,n2
              do i=1,n1
                 r2=rk1(i)**2+rk2(j)**2+rk3(k)**2
                 fhat(i,j,k)=fhat(i,j,k)/r2*(1-exp(-r2/R1))
              enddo
           enddo
        enddo

        elseif (ifder .eq. 1) then

        do k=1,n3
           do j=1,n2
              do i=1,n1
                 r2=rk1(i)**2+rk2(j)**2+rk3(k)**2
                 rn1=rk1(i)*rn1x+rk2(j)*rn1y+rk3(k)*rn1z
                 rn2=rk1(i)*rn2x+rk2(j)*rn2y+rk3(k)*rn2z

                 fhat(i,j,k)=fhat(i,j,k)*rn1*rn2/r2*(1-exp(-r2/R1))
              enddo
           enddo
        enddo

        endif           

        fhat(n1/2+1,n2/2+1,n3/2+1)=0
        
        call ifftshift3(n1,n2,n3,fhat,fhat2)
        call ifft3(n1,n2,n3,fhat2,fhat)
        call ifftshift3(n1,n2,n3,fhat,fhat2)

        ncx=nx*(kfac-1)/2
        ncy=ny*(kfac-1)/2
        ncz=nz*(kfac-1)/2

        uc=fhat2((/(ncx+i, i=1,nx)/),(/(ncy+i, i=1,ny)/),
     &      (/(ncz+i, i=1,nz)/))
        
        return
        end
c
c
c
c
        subroutine ifftshift3(n1,n2,n3,fin,fout)
        implicit real *8 (a-h,o-z)
        complex *16 fin(n1,n2,n3), fout(n1,n2,n3)
        
        m1=n1/2
        m2=n2/2
        m3=n3/2

        fout=fin((/ (/(i,i=m1+1,n1)/), (/(i,i=1,m1)/)/),
     &      (/ (/(i,i=m2+1,n2)/), (/(i,i=1,m2)/)/),
     &      (/ (/(i,i=m3+1,n3)/), (/(i,i=1,m3)/)/))

        return
        end
c
c
c
c
        subroutine uradpart(nr,nzt,nphi,hx,hy,hz,x,y,z,nx,ny,nz,
     &    rn1x,rn1y,rn1z,rn2x,rn2y,rn2z,R1,f0,u2,ifder) 
        implicit real *8 (a-h,o-z)
        dimension x(nx),y(ny),z(nz)
        real *8, allocatable :: r(:),zt(:),phi(:),wr(:),wz(:),wphi(:)
        real *8, allocatable :: rk1(:,:,:),rk2(:,:,:),rk3(:,:,:)
	real *8, allocatable :: x1(:),y1(:),z1(:)
        parameter (pi=3.141592653589793d0)
        complex *16 u2(nx,ny,nz),f0(nx,ny,nz)
        complex *16, allocatable :: fhat(:,:,:),f1(:)

        allocate( r(nr) )
        allocate( wr(nr) )

        allocate( zt(nzt) )
        allocate( wz(nzt) )

        allocate( phi(nphi) )
        allocate( wphi(nphi) )

	nt=nphi*nzt*nr
	allocate(x1(nt))
	allocate(y1(nt))
	allocate(z1(nt))
	allocate(f1(nt))


        allocate( rk1(nphi,nzt,nr) )
        allocate( rk2(nphi,nzt,nr) )
        allocate( rk3(nphi,nzt,nr) )
        allocate( fhat(nphi,nzt,nr) )

        call rnodes3(nr,nzt,nphi,6*sqrt(R1),r,wr,zt,wz,phi,wphi)
        call fhateval(r,zt,phi,nr,nzt,nphi,nx,ny,nz,f0,hx,hy,hz,fhat)

        cx=2*pi*hx
        cy=2*pi*hy
        cz=2*pi*hz

        if (ifder .eq. 0) then

        do k=1,nr
           v2=exp(-r(k)**2/R1)
           do j=1,nzt
              do i=1,nphi
                 rk1(i,j,k)=r(k)*sqrt(1-zt(j)**2)*cos(phi(i))
                 rk2(i,j,k)=r(k)*sqrt(1-zt(j)**2)*sin(phi(i))
                 rk3(i,j,k)=r(k)*zt(j)

                 fhat(i,j,k)=fhat(i,j,k)*wphi(i)*wz(j)*wr(k)*v2
              enddo
           enddo
        enddo

        elseif (ifder .eq. 1) then

        do k=1,nr
           v2=exp(-r(k)**2/R1)
           do j=1,nzt
              do i=1,nphi
                 rk1(i,j,k)=r(k)*sqrt(1-zt(j)**2)*cos(phi(i))
                 rk2(i,j,k)=r(k)*sqrt(1-zt(j)**2)*sin(phi(i))
                 rk3(i,j,k)=r(k)*zt(j)

                 rn1=rk1(i,j,k)*rn1x+rk2(i,j,k)*rn1y+rk3(i,j,k)*rn1z
                 rn2=rk1(i,j,k)*rn2x+rk2(i,j,k)*rn2y+rk3(i,j,k)*rn2z
                 fhat(i,j,k)=fhat(i,j,k)*wphi(i)*wz(j)*wr(k)*v2*rn1*rn2
              enddo
           enddo
        enddo

        endif

        rk1=cx*rk1
        rk2=cy*rk2
        rk3=cz*rk3


	x1=reshape(rk1,(/nt/))
	y1=reshape(rk2,(/nt/))
	z1=reshape(rk3,(/nt/))
	f1=reshape(fhat,(/nt/))

        iflag=1
        eps=1e-12
	t1=second()
        call nufft3d1f90(nt,x1,y1,z1,f1,iflag,eps,nx,ny,nz,u2,ier)
	t2=second()

        u2=u2*nt

	call prin2('after nufft3d1, time(sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine rnodes3(nr,nz,nphi,R1,r,wr,z,wz,phi,wphi)
        implicit real *8 (a-h,o-z)
        dimension r(nr),z(nz),phi(nphi),wr(nr),wz(nz),wphi(nphi)
        parameter (pi=3.141592653589793d0)

        itype=1
        call legeexps(itype,nz,z,u,v,wz)

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
        subroutine fhateval(r,z,phi,nr,nzt,nphi,nx,ny,nz,fx,hx,hy,hz,
     &    fhat)
        implicit real *8 (a-h,o-z)
        dimension r(nr),z(nzt),phi(nphi)
        real *8, allocatable :: rk1(:,:,:),rk2(:,:,:),rk3(:,:,:)
	real *8, allocatable :: x1(:),y1(:),z1(:)
        parameter (pi=3.141592653589793d0)
        complex *16 fhat(nphi,nzt,nr),fx(nx,ny,nz)
	complex *16, allocatable :: f1(:)

        cx=2*pi*hx
        cy=2*pi*hy
        cz=2*pi*hz

        nt=nr*nzt*nphi
	allocate(x1(nt))
	allocate(y1(nt))
	allocate(z1(nt))
	allocate(f1(nt))

        allocate( rk1(nphi,nzt,nr) )
        allocate( rk2(nphi,nzt,nr) )
        allocate( rk3(nphi,nzt,nr) )

        do k=1,nr
           do j=1,nzt
              do i=1,nphi
                 rk1(i,j,k)=cx*r(k)*sqrt(1-z(j)**2)*cos(phi(i))
                 rk2(i,j,k)=cy*r(k)*sqrt(1-z(j)**2)*sin(phi(i))
                 rk3(i,j,k)=cz*r(k)*z(j)
              enddo
           enddo
        enddo

	x1=reshape(rk1,(/nt/))
	y1=reshape(rk2,(/nt/))
	z1=reshape(rk3,(/nt/))

        iflag=-1
        eps=1e-12
	t1=second()
        call nufft3d2f90(nt,x1,y1,z1,f1,iflag,eps,nx,ny,nz,fx,ier)
        t2=second()

	fhat=reshape(f1,(/nphi,nzt,nr/))
        fhat=fhat*hx*hy*hz

	call prin2('after nufft3d2, time (sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine zeropad(nx,ny,nz,nfac,fin,fout)
        implicit real *8 (a-h,o-z)
        complex *16 fin(nx,ny,nz),fout(nx*nfac,ny*nfac,nz*nfac)

        fout=0

        ncx=nx*(nfac-1)/2;
        ncy=ny*(nfac-1)/2;
        ncz=nz*(nfac-1)/2;

        fout((/(ncx+i,i=1,nx)/),(/(ncy+i,i=1,ny)/),
     &	   (/(ncz+i,i=1,nz)/))=fin

        return
        end
