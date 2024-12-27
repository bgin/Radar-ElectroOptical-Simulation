c
c     file tcud2.f
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
c     test program for the complex MUDPACK solver cud2
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cud2.f, cudcom.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for cud2
c
c **********************************************************
c **********************************************************
c
c     a sample program/test driver for cud2 is below.  it can be
c     executed as an initial test.  the output is listed for the
c     test case described.
c
c     test the driver below by solving the complex linear elliptic pde
c
c        cmplx(1.+x*x,1.+y*y)*pxx +
c
c        cmplx(exp(-x),exp(-y))*(pyy - py) +
c
c        cmplx(y,x)*p(x,y) = r(x,y)
c
c     on the unit square with specified boundary conditions at
c     xb = 1.0, yc = 1.0 and mixed boundary conditions
c
c           dp/dx - cmplx(y,y)*p(xa,y) = g(y) at x = xa
c     and
c
c           dp/dy + cmplx(x,x)*p(x,yd) = h(x) at y = yd.
c
c     the exact solution
c
c          p(x,y) = cmplx(x**5,y**5) + 1.0
c
c     is used for testing.  one full multigrid cycle (no initial guess)
c     with red/black gauss-seidel point relaxation and the default multigrid
c     options is sufficient to reach discretization level error on a
c     49 x 65 grid.
c
c
c ******************************************************
c     output (32 bit floating point arithmetic)
c *******************************************************
c
c     cud2 test
c
c     integer input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
c     ixp =  3 jyq =  2 iex =  5 jey =  6
c     nx =  49 ny =  65 iguess =  0 maxcy =  1
c     method =  0 work space estimate =   42776
c
c     multigrid option arguments
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     floating point input parameters
c     xa =  0.000 xb =  1.000 yc =  0.000 yd =  1.000
c     tolerance (error control) =    0.000E+00
c
c     discretization call to cud2 intl =  0
c     ierror =  0 minimum work space =   34192
c
c     approximation call to cud2
c     intl =  1 method =  0 iguess =  0
c     ierror =  0
c     maximum error  =   0.615E-03
c
c ************************************************************
c     end of output
c ************************************************************
c
      program tcud2
      implicit none
      integer iixp,jjyq,iiex,jjey,nnx,nny,llw
c
c     set grid sizes with parameter statements
c
      parameter (iixp = 3 , jjyq = 2 , iiex =5, jjey = 6)
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     estimate work length approximation for method=0 (see cud2.d)
c
      parameter (llw=(40*nnx*nny+8*(nnx+nny+2))/3)
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
      common/itcud2/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftcud2/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      complex pe,px,py,pxx,pyy,cxx,cyy,cx,cy,ce
      real dlx,dly,x,y,errmax
c
c     declare coefficient and boundary condition input subroutines external
c
      external cof,bndc
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
c     set for one multigrid cycle
c
      maxcy = 1
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
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call cof(x,y,cxx,cyy,cx,cy,ce)
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
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     print input parameters (except multigrid options which are default)
c
      write(6,100)
  100 format(//' cud2 test ')
      write (*,101) (iprm(i),i=1,15)
  101 format(/' integer input arguments ',
     +/'intl = ',i2,' nxa = ',i2,' nxb = ',i2,' nyc = ',i2,' nyd = ',i2,
     +/' ixp = ',i2,' jyq = ',i2,' iex = ',i2,' jey = ',i2,
     +/' nx = ',i3,' ny = ',i3,' iguess = ',i2,' maxcy = ',i2,
     +/' method = ',i2, ' work space estimate = ',i7)
      write (*,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option arguments ',
     +/' kcycle = ',i2,
     +/' iprer = ',i2,
     +/' ipost = ',i2,
     +/' intpol = ',i2)
      write(*,103) xa,xb,yc,yd,tolmax
  103 format(/' floating point input parameters ',
     +/' xa = ',f6.3,' xb = ',f6.3,' yc = ',f6.3,' yd = ',f6.3,
     +/' tolerance (error control) =   ',e10.3)
c
c     intiialization call
c
      write(*,104) intl
  104 format(/' discretization call to cud2', ' intl = ', i2)
      call cud2(iprm,fprm,w,cof,bndc,r,p,mgopt,ierror)
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
  106 format(/' approximation call to cud2',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      call cud2(iprm,fprm,w,cof,bndc,r,p,mgopt,ierror)
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
      end

      subroutine cof(x,y,cxx,cyy,cx,cy,ce)
c
c     input pde coefficients at any grid point (x,y) in the solution region
c     (xa.le.x.le.xb,yc.le.y.le.yd) to cud2
c
      implicit none
      real x,y
      complex cxx,cyy,cx,cy,ce
      cxx = cmplx(1.+x*x,1.+y*y)
      cyy = cmplx(exp(-x),exp(-y))
      cx = (0.,0.)
      cy = -cyy
      ce = -cmplx(y,x)
      return
      end

      subroutine bndc(kbdy,xory,alfa,gbdy)
c
c     input mixed derivative b.c. to cud2
c
      implicit none
      integer kbdy
      real xory,x,y
      complex alfa,gbdy
      real xa,xb,yc,yd,tolmax,relmax
      common/ftcud2/xa,xb,yc,yd,tolmax,relmax
      complex pe,px,py,pxx,pyy
      if (kbdy.eq.1) then
c
c     x=xa boundary (nxa must equal 2)
c     b.c. has the form px + alfxa(y)*pe = gbdxa(y)
c     alfa and gbdy corresponding to alfxa(y),gbdxa(y)
c     must be output
c
	y = xory
	x = xa
	call exact(x,y,pxx,pyy,px,py,pe)
	alfa = -cmplx(y,y)
	gbdy = px + alfa*pe
	return
      end if
      if (kbdy.eq.4) then
c
c     y = yd boundary (nyd must equal 2)
c     b.c. has the form py + alfyd(x)*pe = gbdyd(x)
c     alfa and gbdy corresponding to alfyd(x),gbdyd(x)
c     must be output
c
	y = yd
	x = xory
	call exact(x,y,pxx,pyy,px,py,pe)
	alfa = cmplx(x,x)
	gbdy = py + alfa*pe
      return
      end if
      end

      subroutine exact(x,y,pxx,pyy,px,py,pe)
c
c     this subroutine is used to set an exact solution for testing cud2
c
      implicit none
      real x,y
      complex pxx,pyy,px,py,pe
      pe = cmplx(x**5,y**5)+1.0
      px = cmplx(5*x**4,0.)
      py = cmplx(0.,5*y**4)
      pxx = cmplx(20.*x**3,0.)
      pyy = cmplx(0.,20.*y**3)
      return
      end
