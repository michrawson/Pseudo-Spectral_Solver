cc Copyright (C) 2013-2014: Shidong Jiang
cc Contact: shidong.jiang@njit.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
c
c
c
c
c
        subroutine poisson2df(nx,ny,hx,hy,f0,ux,uy,iprec,
     &      aplan, bplan, f2, fhat2, fhatx, fhaty)
c
c       This subroutine solves the 2D Poisson equation -\Delta u = f  
c       in free space and returns the user u_x and u_y on a uniform grid.
c
c       input parameters
c       
c       nx, ny: number of points in x, y directions
c       hx, hy: spatial step size in x, y directions
c       f0: the right hand side over a unifrom grid in R^2
c       iprec: precision control, 0 gives 6 digit  accuracy, 1 gives 12 digit accuracy
c
c       output parameters
c
c       ux : du/dx on the same spatial grid as f0
c       uy : du/dy on the same spatial grid as f0
c
c       WARNING: only works even nx and ny now!
c
        implicit real *8 (a-h,o-z)
        complex *16 ux(nx,ny),uy(nx,ny),f0(nx,ny)
        complex *16, allocatable :: ux2(:,:),uy2(:,:)
        integer *8  :: aplan, bplan
        complex *16 :: fhatx(3*nx,3*ny),fhaty(3*nx,3*ny)
        complex *16 :: f2(3*nx,3*ny),fhat2(3*nx,3*ny)
        data pi/3.141592653589793238462643383279502884197169399d0/

!        fftw_plan_forward=0
!        fftw_plan_back=0

        dk=min(1.0d0/hx/nx,1.0d0/hy/ny)

c        R0=0.1d0
c        R1=0.6d0

        R0=0.8d0*dk
        R1=10*dk
c        call prin2('R1=*',R1,1)

        if (iprec .eq. 1) then
           kfac=3
        else
           kfac=2
        endif

        call ufftpart(hx,hy,nx,ny,kfac,nx*kfac,ny*kfac,R0,R1,f0,ux,uy,
     &     aplan, bplan, f2, fhat2, fhatx, fhaty)

        if (iprec .eq. 1) then
           nr=60
           nphi=60
        else
           nr=40
           nphi=40
        endif

c	nr=50
c	nphi=50

        allocate( ux2(nx,ny) )
        allocate( uy2(nx,ny) )
        call uradpart(nr,nphi,hx,hy,nx,ny,R0,R1,f0,ux2,uy2)

        ux=ux+ux2
        uy=uy+uy2

        d=1/(2*pi)

        ux = ux*d
        uy = uy*d

        return
        end
c
c
c
c
        subroutine fft2(nx,ny,fval,fhat,plan)
        implicit real *8 (a-h,o-z)
        include "fftw3.f"
        integer *8 plan
        complex *16 fval(nx,ny),temp(nx,ny),fhat(nx,ny)

        t1=second()
!        plan=0
        if (plan .eq. 0) then
            temp = fval
!            print *,"comp plan"
            call dfftw_plan_dft_2d(plan, nx,ny, fval,fhat,
     &            FFTW_FORWARD, FFTW_MEASURE)
            fval = temp
        end if
        call dfftw_execute_dft(plan, fval, fhat)
!        call dfftw_destroy_plan(plan)
        t2=second()

        call prin2('after forward fft, time (sec)=*',t2-t1,1)

        return
        end
c
c
c
c
        subroutine ifft2(nx,ny,fhat,fval,plan)
        implicit real *8 (a-h,o-z)
        include "fftw3.f"
        integer *8 plan
        complex *16 fval(nx,ny),temp(nx,ny),fhat(nx,ny)

        t1=second()
!        plan=0
        if (plan .eq. 0) then
            temp = fhat
!            print *,"comp plan"
            call dfftw_plan_dft_2d(plan, nx,ny, fhat,fval,
     &            FFTW_BACKWARD, FFTW_MEASURE)
            fhat = temp
        end if
        call dfftw_execute_dft(plan, fhat, fval)
!        call dfftw_destroy_plan(plan)

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
        subroutine ufftpart(hx,hy,nx,ny,kfac,n1,n2,R0,R1,f0,ux,uy,
     &      aplan, bplan, f2, fhat2, fhatx, fhaty)
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: rk1(:),rk2(:)
        parameter (pi=3.141592653589793d0)
        complex *16 ux(nx,ny),uy(nx,ny),f0(nx,ny),ima
        complex *16, allocatable :: f1(:,:),fhat(:,:)
        complex *16 :: fhatx(n1,n2),fhaty(n1,n2)
	    complex *16 :: f2(n1,n2),fhat2(n1,n2)
        integer *8  :: aplan, bplan
        data ima/(0.0d0,1.0d0)/

!        n1=nx*kfac
!        n2=ny*kfac

        allocate( f1(n1,n2) )
!        allocate( f2(n1,n2) )
        allocate( fhat(n1,n2) )
!        allocate( fhatx(n1,n2) )
!        allocate( fhaty(n1,n2) )
!        allocate( fhat2(n1,n2) )

        call zeropad(nx,ny,kfac,f0,f1)
        call ifftshift2(n1,n2,f1,f2)
!        fftw_plan_forward=0
        call fft2(n1,n2,f2,fhat2,aplan)
        call ifftshift2(n1,n2,fhat2,fhat)

        allocate( rk1(n1) )
        allocate( rk2(n2) )

        rk1=(/ (k, k=-n1/2,n1/2-1) /)
        rk2=(/ (k, k=-n2/2,n2/2-1) /)

        rk1=rk1/hx/n1
        rk2=rk2/hy/n2

        do j=1,n2
           do i=1,n1
              r2=rk1(i)**2+rk2(j)**2
              r=sqrt(r2)
              call erfcbase(r,R0,R1,v2)
              fhatx(i,j)=ima*fhat(i,j)*rk1(i)/r2*(1-v2)
              fhaty(i,j)=ima*fhat(i,j)*rk2(j)/r2*(1-v2)
           enddo
        enddo

        ncx=nx*(kfac-1)/2
        ncy=ny*(kfac-1)/2

        fhatx(n1/2+1,n2/2+1)=0
        call ifftshift2(n1,n2,fhatx,fhat2)
!        fftw_plan_back=0
        call ifft2(n1,n2,fhat2,fhatx,bplan)
        call ifftshift2(n1,n2,fhatx,fhat2)

        ux=fhat2((/(ncx+i, i=1,nx)/),(/(ncy+i, i=1,ny)/))
        

        fhaty(n1/2+1,n2/2+1)=0
        call ifftshift2(n1,n2,fhaty,fhat2)
!        fftw_plan_back=0
        call ifft2(n1,n2,fhat2,fhaty,bplan)
        call ifftshift2(n1,n2,fhaty,fhat2)

        uy=fhat2((/(ncx+i, i=1,nx)/),(/(ncy+i, i=1,ny)/))
        
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
        subroutine uradpart(nr,nphi,hx,hy,nx,ny,R0,R1,f0,ux2,uy2)
        implicit real *8 (a-h,o-z)
        real *8, allocatable :: r(:),phi(:),wr(:),wphi(:)
        real *8, allocatable :: rk1(:,:),rk2(:,:)
	real *8, allocatable :: x1(:),y1(:)
        parameter (pi=3.141592653589793d0)
        complex *16 ux2(nx,ny),uy2(nx,ny),f0(nx,ny),ima
        complex *16, allocatable :: fhat(:,:),f1(:)
        complex *16, allocatable :: fhatx(:,:),fhaty(:,:)
        data ima/(0.0d0,1.0d0)/

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
        allocate( fhatx(nphi,nr) )
        allocate( fhaty(nphi,nr) )

        call rnodes2(nr,nphi,R1,r,wr,phi,wphi)
        call fhateval(r,phi,nr,nphi,nx,ny,f0,hx,hy,fhat)

        cx=2*pi*hx
        cy=2*pi*hy

        do j=1,nr
           call erfcbase(r(j),R0,R1,v2)
           do i=1,nphi
              rk1(i,j)=r(j)*cos(phi(i))
              rk2(i,j)=r(j)*sin(phi(i))
              fhatx(i,j)=fhat(i,j)*wphi(i)*wr(j)*v2*ima*rk1(i,j)/r(j)
              fhaty(i,j)=fhat(i,j)*wphi(i)*wr(j)*v2*ima*rk2(i,j)/r(j)
           enddo
        enddo

        rk1=cx*rk1
        rk2=cy*rk2


	x1=reshape(rk1,(/nt/))
	y1=reshape(rk2,(/nt/))

	f1=reshape(fhat,(/nt/))

        iflag=1
        eps=2e-13
	t1=second()
        call nufft2d1f90(nt,rk1,rk2,fhatx,iflag,eps,nx,ny,ux2,ier)
        call nufft2d1f90(nt,rk1,rk2,fhaty,iflag,eps,nx,ny,uy2,ier)
	t2=second()

        ux2=ux2*nt
        uy2=uy2*nt

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
        eps=2e-13
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

        ncx=nx*(nfac-1)/2
        ncy=ny*(nfac-1)/2

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
