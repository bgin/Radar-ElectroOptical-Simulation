c
c     file tcud3cr.f
c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                  copyright (c) 2008 by UCAR                   *
c     *                                                               *
c     *       University Corporation for Atmospheric Research         *
c     *                                                               *
c     *                      all rights reserved                      *
c     *                                                               *
c     *                     MUDPACK  version 5.0.1                    *
c     *                                                               *
c     *                 A Fortran Package of Multigrid                *
c     *                                                               *
c     *                Subroutines and Example Programs               *
c     *                                                               *
c     *      for Solving Elliptic Partial Differential Equations      *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *                         John Adams                            *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the National Center for Atmospheric Research          *
c     *                                                               *
c     *                Boulder, Colorado  (80307)  U.S.A.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the National Science Foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c ... purpose
c
c     test program for the MUDPACK solver cud3cr
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cud3cr.f, cudcom.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for cud3cr (see cud3cr.d)
c
c **********************************************************
c **********************************************************
c
c     a sample program/test driver for cud3cr is below. it can be
c     executed as an initial test.  The output from executing
c     the code in this file is listed after the problem description.
c
c     Problem Description:
c
c     test cud3cr by solving the complex nonseparable 3-d elliptic pde
c     with cross derivative terms:
c
c          cxx*pxx + cyy*pyy + czz*pzz + cx*px + cy*py + cz*pz + ce*pe +
c
c          cxy*pxy + cxz*pxz + cyz*pyz = r(x,y,z).
c
c     on the region 0 < x < 1, 0 < y < (pi+pi), 0 < z < 1.
c     let s = sin(y), c = cos(y).  the complex coefficients are
c     given by:
c
c          cxx = cmplx(1.0+0.5*s*z,1.0+0.5*c*z)
c
c          cyy = (1.,1.)+cmplx(x,z)
c
c          czz = cmplx(1.0+0.5*s*x,1.0+0.5*c*x)
c
c          cx = cy = cz = (0.0,0.0)
c
c          ce = -cmplx(z,x)
c
c          cxy = cmplx(sqrt(real(cxx)*real(cyy),sqrt(aimag(cxx)*aimag(cyy))
c
c          cxz = (0.0,0.0)
c
c          cyz = cmplx(sqrt(real(cyy)*real(czz),sqrt(aimag(cyy)*aimag(czz))
c
c     assume the solution is periodic in y and is specified at x=0 and z=0.
c     further assume mixed oblique derivative conditions of the form
c
c          px + cmplx(z,1-z)*py + cmplx(s,c)*pz - p(1,y,z) = g(y,z)
c
c     at x=1 and mixed normal derivative conditions of the form
c
c          pz + cmplx(x,s*s)*p(x,y,1)  = h(x,y)
c
c     at z=1.  for testing purposes, use the exact solution
c
c          p(x,y,z) = exp(x*z)*cmplx(c,s)
c
c     to set boundary conditions, the right hand side, and compute
c     exact error.  results from approximating this problem on a
c     25 by 65 by 25 x-y-z grid using cud3cr are given below.  point
c     relaxation and an error tolerance of .001 are used.
c
c ******************************************************
c     output (32 bit floating point arithmetic)
c *******************************************************
c
c      cud3cr test
c
c      input arguments
c      intl =  0
c      nxa =  1 nxb =  2
c      nyc =  0 nyd =  0
c      nze =  1 nzf =  2
c      ixp =  3 jyq =  2 kzr =  3
c      iex =  3 jey =  6 kez =  3
c      nx =  13 ny =  65 nz =  13
c      iguess =  0 maxcy =  2
c      method =  0 work space input =  173394
c      xa =  0.00 xb =  1.00
c      yc =  0.00 yd =  6.28
c      ze =  0.00 zf =  1.00
c      tolmax =  0.000E+00
c
c      multigrid options
c      mgopt(1)  =  0
c
c      new cud3cr arguments
c      icrs
c         1    0    1
c      tol = 0.0010
c      maxit =  10
c
c      initial call
c      intl =  0
c      iguess =  0
c      ierror =  0
c      minimum required work space length =   173394
c
c      approximation call
c      intl =  1
c      iguess =  0
c      ierror =   0
c      number of outer iterations executed =   9
c      relative difference profile:
c      0.6567  0.1462  0.0638  0.0203  0.0099  0.0051  0.0016  0.0011  0.0004
c      exact least squares error =  0.127E-02
c
c ************************************************************
c     end of output
c ************************************************************
c
      program tcd3cr
      implicit none
c
c     set grid sizes and predetermined minimal required equired work
c     with parameter statements
c
      integer iixp,jjyq,kkzr,iiex,jjey,kkez,nnx,nny,nnz,llengt,mmaxit
      parameter(iixp=3,jjyq=2,kkzr=3)
      parameter (iiex=3,jjey=6,kkez=3)
      parameter(nnx = iixp*2**(iiex-1)+1)
      parameter(nny = jjyq*2**(jjey-1)+1)
      parameter(nnz = kkzr*2**(kkez-1)+1)
      parameter(llengt=173394)
      parameter(mmaxit = 10)
      complex rhs(nnx,nny,nnz),phi(nnx,nny,nnz),work(llengt)
      integer iparm(23),mgopt(4),icrs(3)
      integer i,j,k,ierror,maxit,iouter
      real fparm(8),rmax(mmaxit)

c
c     use labelled common to identify integer and floating point arguments
c
      integer intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nx,ny,nz,iguess,maxcy,method,meth2,nwork,lwkmin,itero
      common /iprm/intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,
     +kez,nx,ny,nz,iguess,maxcy,method,meth2,nwork,lwkmin,itero
      real xa,xb,yc,yd,ze,zf,tolmax,relmax
      common/fprm/xa,xb,yc,yd,ze,zf,tolmax,relmax
      real pi,dx,dy,dz,c,s,x,y,z,tol,err2
      complex cxx,cyy,czz,cx,cy,cz,ce,cxy,cyz
      complex pxx,pyy,pzz,px,py,pz,pxy,pxz,pyz,pe
      equivalence(iparm,intl)
      equivalence(fparm,xa)
      external cof,bd3cr,cxyf,cxzf,cyzf
      pi = 4.*atan(1.)
c
c     set interval end points
c
      xa = 0.0
      xb = 1.0
      yc = 0.0
      yd = pi+pi
      ze = 0.0
      zf = 1.0
c
c     set required no error control within multigrid cycling
c
      tolmax = 0.0
c
c     set integer input arguments
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 1
      nxb = 2
      nyc = 0
      nyd = 0
      nze = 1
      nzf = 2
c
c     set grid size arguments from parameter statements
c
      ixp = iixp
      jyq = jjyq
      kzr = kkzr
      iex = iiex
      jey = jjey
      kez = kkez
      nx = nnx
      ny = nny
      nz = nnz
c
c     flag nonzero xy and yz and zero xz cross derivatives terms
c
      icrs(1) = 1
      icrs(2) = 0
      icrs(3) = 1
c
c     set two multigrid cycles per outer iteration
c
      maxcy = 2
c
c     set point relaxation
c
      method = 0
      meth2 = 0
c
c     set work space length input
c
      nwork = llengt
c
c     set uniform grid interval lengths in each dimension
c
      dx = (xb-xa)/(nx-1)
      dy = (yd-yc)/(ny-1)
      dz = (zf-ze)/(nz-1)
c
c     set right hand side and preset solution to zero
c     this also sets specified (Dirchlet) b.c. in phi
c
      do j=1,ny
	y = (j-1)*dy
	s = sin(y)
	do k=1,nz
	  z = (k-1)*dz
	  do i=1,nx
	    x = (i-1)*dx
	    call exact(x,y,z,pe,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz)
	    call cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
	    call cxyf(x,y,z,cxy)
	    call cyzf(x,y,z,cyz)
	    rhs(i,j,k) = cxx*pxx+cyy*pyy+czz*pzz+cxy*pxy+cyz*pyz+
     +                   cx*px+cy*py+cz*pz+ce*pe
	    phi(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     set default multigrid options
c
      mgopt(1) = 0
c
c     set error control tolerance and outer iteration limit of 10
c
      tol = .001
      maxit = mmaxit
c
c     print input arguments shared with cud3
c
      write(*,50)
   50 format(//' cud3cr test ' )
      write(*,100)intl,nxa,nxb,nyc,nyd,nze,nzf,ixp,jyq,kzr,iex,jey,kez,
     +            nx,ny,nz,iguess,maxcy,method,nwork,xa,xb,yc,yd,ze,zf,
     +            tolmax,mgopt(1)
  100 format(/' input arguments ',
     +/' intl = ',i2,
     +/' nxa = ',i2,' nxb = ',i2,
     +/' nyc = ',i2,' nyd = ',i2,
     +/' nze = ',i2,' nzf = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' kzr = ',i2,
     +/' iex = ',i2, ' jey = ',i2, ' kez = ',i2,
     +/' nx = ',i3,' ny = ',i3,' nz = ',i3,
     +/' iguess = ',i2,' maxcy = ',i2,
     +/' method = ',i2, ' work space input = ',i7,
     +/' xa = ',f5.2,' xb = ',f5.2,
     +/' yc = ',f5.2,' yd = ',f5.2,
     +/' ze = ',f5.2,' zf = ',f5.2,
     +/' tolmax = ' ,e10.3
     +//' multigrid options '
     +/' mgopt(1)  = ',i2)
c
c     print new cud3cr arguments
c
      write(*,101) (icrs(i),i=1,3),tol,maxit
  101 format(/' new cud3cr arguments ', / ' icrs '/ 3i5,
     +/' tol = ',f6.4 /' maxit = ',i3)
c
c *** initialization call
c
      intl = 0
      iguess = 0
      call cud3cr(iparm,fparm,work,cof,bd3cr,rhs,phi,mgopt,icrs,
     +cxyf,cxzf,cyzf,tol,maxit,iouter,rmax,ierror)
      write (*,200) intl,iguess,ierror,iparm(22)
  200 format(/' initial call', /' intl = ',i2, /' iguess = ',i2,
     +/' ierror = ',i2, /' minimum required work space length = ',i8)
      if (ierror.gt.0) call exit(0)
c
c *** noninitial call
c
      intl = 1
      iguess = 0
      call cud3cr(iparm,fparm,work,cof,bd3cr,rhs,phi,mgopt,icrs,
     +cxyf,cxzf,cyzf,tol,maxit,iouter,rmax,ierror)
      write (*,201) intl,iguess,ierror,iouter,(rmax(i),i=1,iouter)
  201 format(/' approximation call ', /' intl = ',i2, /' iguess = ',i2,
     +/' ierror = ',i3,
     +/' number of outer iterations executed = ',i3,
     +/' relative difference profile:', /(10(f6.4,2x)))
c
c     compute and print exact least squares error after iouter iterations
c
      err2 = 0.
      do j=1,ny
	y = (j-1)*dy
	c = cos(y)
	do k=1,nz
	  z = (k-1)*dz
	  do i=1,nx
	    x = (i-1)*dx
	    pe = cmplx((x*c*z)**3,0.0)
	    err2 = err2 + cabs(pe - phi(i,j,k))**2
	  end do
	end do
      end do
      err2 = sqrt(err2/(nx*ny*nz))
      write(*,202) err2
  202 format(' exact least squares error = ',e10.3)
      end

      subroutine exact(x,y,z,pe,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz)
c
c     set exact solution and partial derivatives: p(x,y,z) = (x*cos(y)*z)**3
c
      implicit none
      real x,y,z
      complex pe,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz
      real s,c,exz
      real xc,xz,cz,xcz
      c = cos(y)
      s = sin(y)
      xc = x*c
      xz = x*z
      cz = c*z
      xcz = xc*z
      pe = cmplx(xcz**3,0.0)
      px = cmplx(3.*x*x*cz**3,0.0)
      pxx = cmplx(6.*x*cz**3,0.0)
      py =cmplx( -3.*c*c*s*xz**3,0.0)
      pyy = cmplx(-3.*xz**3*(c**3-2.*c*s**2),0.0)
      pz = cmplx(3.*z*z*xc**3,0.0)
      pzz = cmplx(6.*z*xc**3,0.0)
      pxy = cmplx(-9.*x**2*c*c*s*z**3,0.0)
      pxz = cmplx(9.*xz**2*c**3,0.0)
      pyz =cmplx( -9.*c*c*s*x**3*z**2,0.0)
      return
      exz = exp(x*z)
      pe = exz*cmplx(c,z)
      px = z*pe
      py = -exz*s
      pz = x*pe +exz*cmplx(c,1.0)
      pxx = z*px
      pyy = -exz*c
      pzz =  x*(pz+exz*cmplx(c,1.0))
      return
      end

      subroutine cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
c
c     set noncross derivative pde coefficients
c
      implicit none
      real x,y,z,s,c
      complex cxx,cyy,czz,cx,cy,cz,ce
      s = sin(y)
      c = cos(y)
      cxx = cmplx(1.+0.5*s*z, 1.+0.5*c*z)
      cyy = cmplx(1.+x,1.+z)
      czz = cmplx(1+0.5*s*x,1.+0.5*c*x)
      cx = (0.,0.)
      cy = cx
      cz = cz
      ce = -cmplx(x+z,0.0)
      return

      cxx = cmplx(1.+0.5*s*z,0.0)
      cyy = cmplx(1.+x*z,0.0)
      czz = cmplx(1.+0.5*x*s,0.0)
      cx = (0.0,0.0)
      cy = cmplx(-x*z,0.0)
      cz = (0.0,0.0)
      ce = cmplx(-(x+z),0.0)
      return

      cxx = (1.0,1.0)+0.5*cmplx(s*z,c*z)
      cyy = (1.0,1.0)+cmplx(x,z)
      czz = (1.0,1.0)+0.5*cmplx(x*s,x*c)
      cx = (0.0,0.0)
      cy = cmplx(-x*z,x*z)
      cz = (0.0,0.0)
      ce = -cmplx((x+z),0.0)
      return
      end

      subroutine cxyf(x,y,z,cxy)
c
c     set x-y cross term coefficient at (x,y,z)
c
      implicit none
      real x,y,z,cxxr,cxxi,cyyr,cyyi,cxyr,cxyi
      complex cxx,cyy,czz,cx,cy,cz,ce,cxy
      call cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
      cxxr = cxx
      cxxi = aimag(cxx)
      cyyr = cyy
      cyyi = aimag(cyy)
      cxyr = sqrt(cxxr*cyyr)
      cxyi = sqrt (cxxi*cyyi)
      cxy = cmplx(cxyr,cxyi)
      return
      end

      subroutine cxzf(x,y,z,cxz)
c
c     this is a dummy subroutine since cxz=0.0 for all (x,y,z)
c
      return
      end

      subroutine cyzf(x,y,z,cyz)
c
c     set y-z cross term coefficient at (x,y,z)
c
      real x,y,z,czzr,czzi,cyyr,cyyi,cyzr,cyzi
      complex cxx,cyy,czz,cx,cy,cz,ce,cyz
      call cof(x,y,z,cxx,cyy,czz,cx,cy,cz,ce)
      cyyr = cyy
      cyyi = aimag(cyy)
      czzr = czz
      czzi = aimag(czz)
      cyzr = sqrt(cyyr*cyyi)
      cyzi = sqrt (cyyi*czzi)
      cyz = cmplx(cyzr,cyzi)
      return
      end

      subroutine bd3cr(kbdy,xory,yorz,a,b,c,g)
c
c     pass mixed derivative boundary conditions to cud3cr
c
      implicit none
      integer kbdy
      real x,y,z,xory,yorz
      complex a,b,c,g
      complex pe,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz
      if (kbdy.eq.2) then
c
c     upper x boundary (mixed oblique)
c
	x = 1.0
	y = xory
	z = yorz
	a = cmplx(z,1.-z)
	b = cmplx(sin(y),cos(y))
	c = - (1.0,1.0)
	call exact(x,y,z,pe,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz)
	g = px + a*py + b*pz + c*pe
	return
      else if (kbdy.eq.6) then
c
c     upper z boundary (mixed normal)
c
	z = 1.0
	x = xory
	y = yorz
	c = cmplx(x,sin(y)**2)
	call exact(x,y,z,pe,px,py,pz,pxx,pyy,pzz,pxy,pxz,pyz)
	g = pz + a*px + b*py + c*pe
	return
      end if
      end

