c
c     file tcud24sp.f
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
c     test program for the complex MUDPACK solver cud24sp
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cud24sp.f, cud2sp.f, cudcom.f
c
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for cud24sp
c
c **********************************************************
c **********************************************************
c
c
c     a sample program/test driver for cud24sp is listed below.  it
c     can be executed as an initial test.  the output is listed for
c     the test case described.
c
c     test the driver below by solving the complex separable elliptic pde
c
c          cmplx(1.+x*x,1-x*x)*pxx + cmplx(exp(1.-y),exp(y))*(pyy-py) -
c
c          (cmplx(x,x) + cmplx(y,y))*pe(x,y) = r(x,y)
c
c     on the (x,y) region
c
c         [1/4,3/4] X [1/3,2/3]
c
c     with specified boundaries at upper x boundry xb = 3/4 and
c     the lower y boundary yc = 1/3 and mixed boundary conditions
c     at the lower x boundry xa = 1/4
c
c          dp/dx - p = g(y)
c
c     and upper y boundary yd = 2/3
c
c          dp/dy + p = h(x)
c
c     Use point relaxation, the default multigrid options and
c     the exact solution
c
c          pe(x,y) = cmplx(y**5,x**4+1.0)
c
c     to set the righthand side, boundary conditions, and compute
c     error.  the approximation is generated on a 129 by 97 grid.
c     First cud2sp is called to render a second-order approximation.
c     then cud24sp is called for a fourth-order estimate
c
c
c *************************************************************
c     output (64 bit floating point arithmetic)
c *************************************************************
c
c     cud2sp test
c
c     integer input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
c     ixp =  2 jyq =  3 iex =  7 jey =  6
c     nx = 129 ny =  97 iguess =  0 maxcy =  3
c     method =  0 work space estimate =   48976
c
c     multigrid option arguments
c     kcycle =  0
c     iprer =  0
c     ipost =  0
c     intpol =  0
c
c     floating point input parameters
c     xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
c     tolerance (error control) =    0.000E+00
c
c     discretization call to cud2sp intl =  0
c     ierror =  0 minimum work space =   48976
c
c     approximation call to cud2sp
c     intl =  1 method =  0 iguess =  0
c     ierror =  0
c     maximum error  =   0.317E-03
c
c     cud24sp test  ierror =  0
c     maximum error  =   0.621E-06
c
c ***************************************************************
c      end of output
c ***************************************************************
c
      program tcud24sp
      implicit none
      integer iixp,jjyq,iiex,jjey,nnx,nny,llw
c
c     set grid sizes with parameter statements
c
      parameter (iixp = 2 , jjyq = 3 , iiex =7, jjey = 6)
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     set minimum required work space (see tmud2sp.f)
c
      parameter (llw = 48976)
c
c     dimension solution,right hand side, and work arrays
c
      complex p(nnx,nny),r(nnx,nny),w(llw)
c     put integer and floating point parameter names in contiguous
c     storeage for labelling purposes
c
      integer iprm(16),mgopt(4)
      real fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itcud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftcud2sp/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      complex cxx,cyy,cx,cy,ce,pxx,pyy,px,py,pe
      complex cex,cey
      real dlx,dly,x,y,errmax
c
c     declare coefficient and boundary condition input subroutines external
c
      external cfx,cfy,bndc
c
c
c     set input integer arguments
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 2
      nxb = 1
      nyc = 1
      nyd = 2
c
c     set grid sizes from parameter statements
c
      ixp = iixp
      jyq = jjyq
      iex = iiex
      jey = jjey
      nx = nnx
      ny = nny
c
c     set three multigrid cycles
c
      maxcy = 3
c
c     set work space length approximation from parameter statement
c
      nwork = llw
c
c     set point relaxation
c
      method = 0
c
c     flag no initial guess (this sets full multigrid cycling)
c
      iguess = 0
c
c     set end points of solution rectangle in (x,y) space
c
      xa = 0.0
      xb = 1.0
      yc = 0.0
      yd = 1.0
c
c     set mesh increments
c
      dlx = (xb-xa)/float(nx-1)
      dly = (yd-yc)/float(ny-1)
c
c     set for no error control
c
      tolmax = 0.0
c
c     set right hand side in r
c     initialize p to zero
c
      do i=1,nx
	x = xa+float(i-1)*dlx
	call cfx(x,cxx,cx,cex)
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call cfy(y,cyy,cy,cey)
	  ce = cex+cey
	  call exact(x,y,pxx,pyy,px,py,pe)
	  r(i,j) = cxx*pxx+cyy*pyy+cx*px+cy*py+ce*pe
	  p(i,j) = (0.0,0.0)
	end do
      end do
c
c     set specified boundaries in p
c
      x = xb
      do j=1,ny
	y = yc+float(j-1)*dly
	call exact(x,y,pxx,pyy,px,py,pe)
	p(nx,j) = pe
      end do
      y = yc
      do i=1,nx
	x = xa+float(i-1)*dlx
	call exact(x,y,pxx,pyy,px,py,pe)
	p(i,1) = pe
      end do
c
c     set default multigrid opitons
c
      mgopt(1) = 0
c
c     print input parameters (except multigrid options which are default)
c
      write(6,100)
  100 format(//' cud2sp test ')
      write (*,101) (iprm(i),i=1,15)
  101 format(/' integer input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2
     +/' nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2,
     +/' method = ',i2, ' work space estimate = ',i7)
      write (*,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option arguments ',
     +/' kcycle = ',i2,
     +/' iprer = ',i2,
     +/' ipost = ',i2
     +/' intpol = ',i2)
      write(*,103) xa,xb,yc,yd,tolmax
  103 format(/' floating point input parameters ',
     +/' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3,
     +/' tolerance (error control) =   ',e10.3)
c
c     intiialization call
c
      write(*,104) intl
  104 format(/' discretization call to cud2sp', ' intl = ', i2)
      call cud2sp(iprm,fprm,w,cfx,cfy,bndc,r,p,mgopt,ierror)
c
c     print error parameter and minimum work space requirement
c
      write (*,105) ierror,iprm(16)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
      write(*,106) intl,method,iguess
  106 format(/' approximation call to cud2sp ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      call cud2sp(iprm,fprm,w,cfx,cfy,bndc,r,p,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
c
c     compute and print exact maximum error
c
      errmax = 0.0
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call exact(x,y,pxx,pyy,px,py,pe)
	  errmax = amax1(errmax,cabs((p(i,j)-pe)))
	end do
      end do
      write(*,108) errmax
  108 format(' maximum error  =  ',e10.3)
c
c     fourth-order estimate
c
      call cud24sp(w,p,ierror)
      write(*,109) ierror
  109 format(/' cud24sp test ', ' ierror = ',i2)
      errmax = 0.0
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call exact(x,y,pxx,pyy,px,py,pe)
	  errmax = amax1(errmax,cabs((p(i,j)-pe)))
	end do
      end do
      write(*,108) errmax
      end

      subroutine cfx(x,cxx,cx,cex)
c
c     input x dependent complex coefficients
c
      implicit none
      real x
      complex cxx,cx,cex
      cxx = cmplx(1.+x*x,1.-x*x)
      cx = (0.,0.)
      cex = -cmplx(x,x)
      return
      end

      subroutine cfy(y,cyy,cy,cey)
c
c     input y dependent complex coefficients
c
      implicit none
      real y
      complex cyy,cy,cey
      cyy = cmplx(exp(1.0-y),exp(y))
      cy = -cyy
      cey = -cmplx(y,y)
      return
      end

      subroutine bndc(kbdy,xory,cons,gbdy)
c
c     input mixed complex derivative b.c. to cud2sp
c
      implicit none
      integer kbdy
      real xory,x,y
      complex cons,gbdy,pe,px,py,pxx,pyy
      real xa,xb,yc,yd,tolmax,relmax
      common/ftcud2sp/xa,xb,yc,yd,tolmax,relmax
      if (kbdy.eq.1) then
c
c     x=xa boundary
c
	y = xory
	x = xa
	cons =(-1.0,0.0)
	call exact(x,y,pxx,pyy,px,py,pe)
	gbdy = px + cons*pe
	return
      end if
      if (kbdy.eq.4) then
c
c     y=yd boundary
c
	y = yd
	x = xory
	cons = (1.0,0.0)
	call exact(x,y,pxx,pyy,px,py,pe)
	gbdy = py + cons*pe
	return
      end if
      end

      subroutine exact(x,y,pxx,pyy,px,py,pe)
c
c     set exact solution used for testing cud2sp
c
      implicit none
      real x,y
      complex pe,px,py,pxx,pyy
      pe = cmplx(y**5.,x**4+1.)
      px = cmplx(0.,4.*x**3)
      pxx = cmplx(0.,12.*x**2)
      py = cmplx(5.*y**4,0.)
      pyy = cmplx(20.*y**3,0.)
      return
      end
