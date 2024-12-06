c
c     file tcuh24cr.f
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
c     test program for the complex MUDPACK solver cuh24cr
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     cuh24cr.f, cuh2cr.f, cudcom.f
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for cuh24cr
c
c **********************************************************
c **********************************************************
c
c     a sample program/test driver for cuh24cr is below.  it can
c     be executed as an initial test.  the output is listed for
c     the test case described.
c
c     test the driver below by solving the elliptic pde
c
c      cmplx(1.+y*y,1.-y*y)*d(dp/dx)/dx + cmplx(x*y,-x*y)*d(dp/dx)/dy
c
c      cmplx(1.+x*x,1.-x*x)*d(dp/dy)/dy + cmplx(y,-y)*dp/dx +
c
c      cmplx(x,-x)*dp/dy + cmplx(x+y,-x-y)*p(x,y) = r(x,y)
c
c     on the region
c
c       0.25 < x < 1.0, 0.0 < y < 0.50
c
c     with specified boundary conditions at xa=0.25, xb=1.0, yc=0.0
c     and mixed boundary conditions at yd=0.5 of the form:
c
c     cmplx(-x,x)*dp/dx+cmplx(1.+x,1.-x)*dp/dy+cmplx(-x,x)*p(x,yd)=gbdyd(x)
c
c     choose a 45 x 73 grid and use line relaxation in the y direction.
c     forr testing purposes use the exact solution
c
c       p(x,y) = cmplx((x*y)**3,-(x*y)**3) + (1.,1.)
c
c     Choosing grid arguments
c
c       ixp = 11, iex = 3, jyq = 9, jey = 4
c
c     fits the required 45 X 73 grid exactly with the coarsening
c
c       45 X 73 > 23 X 37 > 12 X 19 > 12 X 10
c
c     The 12 X 10 coarsest grid has too many points for effective
c     error reduction with relaxation alone.  cuh2cr maintains
c     multigrid convergence efficiency by utilitzing a direct method
c     (Gaussian elimination) whenever the coarsest grid is encountered
c     within multigrid cycling.  Two cycles with cuh2cr are executed
c     and then cuh24cr is called for a fourth-order estimate
c
c
c ******************************************************
c     output (64 bit floating point arithmetic)
c *******************************************************
c
c     cuh2cr test
c
c     integer input arguments
c     intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
c     ixp = 11 jyq =  9 iex =  3 jey =  4
c     nx =  45 ny =  73 iguess =  0 maxcy =  2
c     method =  2 work space estimate =   69661
c
c     multigrid option arguments
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     floating point input parameters
c     xa =  0.250 xb =  1.000 yc =  0.000 yd =  0.500
c     tolerance (error control) =    0.000E+00
c
c     discretization call to cuh2cr intl =  0
c     ierror =  0 minimum work space =   69661
c
c     approximation call to cuh2cr
c     intl =  1 method =  2 iguess =  0 maxcy =  2
c     ierror =  0
c     maximum error  =   0.114E-04
c
c     cuh24cr test  ierror =  0
c     maximum error  =   0.739E-07
c
c **********************************************************
c     end of output
c **********************************************************
c
      program tcuh24cr
      implicit none
      integer iixp,jjyq,iiex,jjey,nnx,nny,llw,mmx,mmy
c
c     set grid sizes with parameter statements
c
      parameter (iixp =11 , jjyq = 9 , iiex =3, jjey = 4)
      parameter (mmx = iixp+1,mmy = jjyq+1)
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     set exact minimal required work space (see tcuh2cr.f)
c
      parameter (llw = 69661)
c
c     dimension solution,right hand side, and work arrays
c
      complex p(nnx,nny),r(nnx,nny),w(llw)
      integer iw(mmx,mmy)
c
c     put integer and floating point parameter names in contiguous
c     storeage for labelling purposes
c
      integer iprm(16),mgopt(4)
      real fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itcuh2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftcur2/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      complex cxx,cxy,cyy,cx,cy,ce,pxx,pxy,pyy,px,py,pe
      real dlx,dly,x,y,errmax
c
c     declare coefficient and boundary condition input subroutines external
c
      external cofcr,bndcr
c
c     set input integer arguments
c
      intl = 0
c
c     set boundary condition flags
c
      nxa = 1
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
c     set two multigrid cycles
c
      maxcy = 2
c
c     set work space length approximation from parameter statement
c
      nwork = llw
c
c     set point relaxation
c
      method = 2
c
c     flag no initial guess (this sets full multigrid cycling)
c
      iguess = 0
c
c     set end points of solution rectangle in (x,y) space
c
      xa = 0.25
      xb = 1.0
      yc = 0.0
      yd = 0.50
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
	  call cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
	  call exact(x,y,pxx,pxy,pyy,px,py,pe)
	  r(i,j) = cxx*pxx+cxy*pxy+cyy*pyy+cx*px+cy*py+ce*pe
	  p(i,j) = (0.0,0.0)
	end do
      end do
c
c     set specified boundaries in p
c
      x = xa
      do j=1,ny
	y = yc+float(j-1)*dly
	call exact(x,y,pxx,pxy,pyy,px,py,pe)
	p(1,j) = pe
      end do
      x = xb
      do j=1,ny
	y = yc+float(j-1)*dly
	call exact(x,y,pxx,pxy,pyy,px,py,pe)
	p(nx,j) = pe
      end do
      y = yc
      do i=1,nx
	x = xa+float(i-1)*dlx
	call exact(x,y,pxx,pxy,pyy,px,py,pe)
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
c     print input arguments
c
      write(6,100)
  100 format(//' cuh2cr test ')
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
  104 format(/' discretization call to cuh2cr', ' intl = ', i2)
      call cuh2cr(iprm,fprm,w,iw,cofcr,bndcr,r,p,mgopt,ierror)
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
      write(*,106) intl,method,iguess,maxcy
  106 format(/' approximation call to cuh2cr',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      call cuh2cr(iprm,fprm,w,iw,cofcr,bndcr,r,p,mgopt,ierror)
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
	  call exact(x,y,pxx,pxy,pyy,px,py,pe)
	  errmax = amax1(errmax,cabs((p(i,j)-pe)))
	end do
      end do
      write(*,108) errmax
  108 format(' maximum error  =  ',e10.3)
c
c     attempt fourth-order estimate with difference corrections
c
      call cuh24cr(w,iw,cofcr,bndcr,p,ierror)
      write (*,109) ierror
  109 format(/ ' cuh24cr test ', ' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
      errmax = 0.0
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call exact(x,y,pxx,pxy,pyy,px,py,pe)
	  errmax = amax1(errmax,cabs((p(i,j)-pe)))
	end do
      end do
      write(*,108) errmax
      end

      subroutine cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
c
c     input pde coefficients at any point (x,y) in the solution region
c     (xa.le.x.le.xb,yc.le.y.le.yd) to cuh2cr
c
      implicit none
      complex cxx,cxy,cyy,cx,cy,ce
      real x,y,xa,xb,yc,yd,tolmax,relmax
      common/ftcur2/xa,xb,yc,yd,tolmax,relmax
      cxx = cmplx(1.+y**2,1.-y**2)
      cxy = cmplx(x*y,-x*y)
      cyy = cmplx(1.+x**2,1.-x*x)
      cx = cmplx(y,-y)
      cy = cmplx(x,-x)
      ce = -cmplx((x+y),-x-y)
      return
      end

      subroutine bndcr(kbdy,xory,alfa,beta,gama,gbdy)
c
c     input mixed "oblique" derivative b.c. to cuh2cr
c
      implicit none
      integer kbdy
      real xory,x,y,xa,xb,yc,yd,tolmax,relmax
      complex alfa,beta,gama,gbdy,pe,px,py,pxx,pxy,pyy
      common/ftcur2/xa,xb,yc,yd,tolmax,relmax
      if (kbdy.eq.4) then
c
c     y=yd boundary (nyd must equal 2 if this code is to be executed).
c     b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
c     where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
c     betyd(x),gamyd(x),gbdyd(y) must be output.
c
	y = yd
	x = xory
	alfa = -cmplx(x,-x)
	beta = cmplx(1.+x,1.-x)
	gama = -cmplx(x,-x)
	call exact(x,y,pxx,pxy,pyy,px,py,pe)
	gbdy = alfa*px + beta*py + gama*pe
	return
      end if
      end

      subroutine exact(x,y,pxx,pxy,pyy,px,py,pe)
c
c     this subroutine is used for setting an exact solution in order
c     to test subroutine cuh2cr.
c
      implicit none
      real x,y,xy,xy2,xy3
      complex pxx,pxy,pyy,px,py,pe
      xy = x*y
      xy2 = xy*xy
      xy3 = xy2*xy
      pe = cmplx(xy3,-xy3)
      px = 3.*y*cmplx(xy2,-xy2)
      py = 3.*x*cmplx(xy2,-xy2)
      pxx = 6.*y*y*cmplx(xy,-xy)
      pxy = 9.*xy*cmplx(xy,-xy)
      pyy = 6.*x*x*cmplx(xy,-xy)
      return
      end

