c
c     file tmuh2cr.f
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
c     test program for the mudpack solver muh2cr
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     muh2cr.f, mudcom.f
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for muh2cr
c
c **********************************************************
c **********************************************************
c
c
c     a sample program/test driver for muh2cr is listed below.  it
c     can be executed as an initial test.  the output is listed
c     for the test case described.
c
c     test muh2cr below by solving the nonseparable elliptic pde
c     with cross derivative term
c
c          (1.+y**2)*pxx + (1.+x**2)*pyy + 2.*x*y*pxy +
c
c          y*px + x*py - (x*y)*pe = r(x,y)
c
c     on a grid as close to 60 by 70 as muh2cr size constraints
c     allow.  the solution region is the unit square.  assume a
c     mixed derivative boundary condition at y=1 of the form
c
c          -x * dp/dx + (1+x) * dp/dy - x * pe = gbdyd(x).
c
c     and specified (Dirchlet) boundary conditions elsewhere.  the
c     exact solution
c
c          p(x,y) = (x*y)**5
c
c     is used to set the right hand side, boundary conditions, and
c     compute the error.
c
c     red/black gauss-seidel point relaxation is used along with the
c     the default multigrid options.  one full multigrid cycle reaches
c     discretization level error for this problem.  Choosing grid
c     size parameters
c
c       ixp = 15, jyq=9, iex=3,jey=4
c
c     yields a 61 by 73 grid.  The grid coarsening is
c
c       61 X 73 > 31 X 37 > 16 X 19 > 16 X 10
c
c     The coarsest 16 X 10 grid has too much resolution for effective
c     error reduction with relaxation only.  muh2cr uses a direct
c     method whenever the 16 X 10 grid is encountered which maintains
c     multigrid convergence efficiency
c
c ******************************************************
c     output (32 bit floating point arithmetic)
c *******************************************************
c
c     muh2cr test
c
c     integer input arguments
c     intl =  0 nxa =  1 nxb =  1 nyc =  1 nyd =  2
c     ixp = 15 jyq =  9 iex =  3 jey =  5
c     nx =  61 ny = 145 iguess =  0 maxcy =  1
c     method =  0 work space estimate =  151335
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
c     discretization call to muh2cr intl =  0
c     ierror =  0 minimum work space =  149055
c
c     approximation call to muh2cr
c     intl =  1 method =  0 iguess =  0
c     ierror =  0
c     maximum error  =   0.412E-03
c
c **********************************************************
c      end of output
c **********************************************************
c
      program tmuh2cr
      implicit none
c
c     set grid size params
c
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
      integer iixp1,jjyq1
      parameter (iixp = 15 , jjyq = 9, iiex = 3, jjey = 5 )
      parameter(iixp1 = iixp+1, jjyq1 = jjyq+1)
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     estimate work space for point relaxation (see muh2cr.d)
c
      parameter (llwork=(7*(nnx+2)*(nny+2)+44*nnx*nny)/3 )
      real phi(nnx,nny),rhs(nnx,nny),work(llwork)
      integer iwork(iixp1,jjyq1)
c
c     put integer and floating point argument names in contiguous
c     storeage for labelling in vectors iprm,fprm
c
      integer iprm(16),mgopt(4)
      real fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmuh2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmuh2cr/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      real dlx,dly,x,y,cxx,cxy,cyy,cx,cy,ce,pxx,pxy,pyy,px,py,pe,errmax
c
c     declare coefficient and boundary condition input subroutines external
c
      external cofcr,bndcr
c
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
c     set multigrid arguments (w(2,1) cycling with fully weighted
c     residual restriction and cubic prolongation)
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     set for one cycle
c
      maxcy = 1
c
c     set no initial guess forcing full multigrid cycling
c
      iguess = 0
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set point relaxation
c
      method = 0
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
c     set for no error control flag
c
      tolmax = 0.0
c
c     set right hand side in rhs
c     initialize phi to zero
c
      do i=1,nx
	x = xa+float(i-1)*dlx
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
	  call exacr(x,y,pxx,pxy,pyy,px,py,pe)
	  rhs(i,j) = cxx*pxx+cxy*pxy+cyy*pyy+cx*px+cy*py+ce*pe
	  phi(i,j) = 0.0
	end do
      end do
c
c     set specified boundaries in phi at x=xa,xb and y=yc
c
      do j=1,ny
	y = yc+float(j-1)*dly
	call exacr(xa,y,pxx,pxy,pyy,px,py,pe)
	phi(1,j) = pe
	call exacr(xb,y,pxx,pxy,pyy,px,py,pe)
	phi(nx,j) = pe
      end do
      do i=1,nx
	x = xa+float(i-1)*dlx
	call exacr(x,yc,pxx,pxy,pyy,px,py,pe)
	phi(i,1) = pe
      end do
      write(*,100)
  100 format(//' muh2cr test ')
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
  104 format(/' discretization call to muh2cr', ' intl = ', i2)
      call muh2cr(iprm,fprm,work,iwork,cofcr,bndcr,rhs,phi,mgopt,ierror)
      write (*,200) ierror,iprm(16)
  200 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
      write(*,106) intl,method,iguess
  106 format(/' approximation call to muh2cr',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      call muh2cr(iprm,fprm,work,iwork,cofcr,bndcr,rhs,phi,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
      if (ierror .le. 0) then
c
c     compute and print maximum norm of error
c
      errmax = 0.0
      do j=1,ny
	y = yc+(j-1)*dly
	do i=1,nx
	  x = xa+(i-1)*dlx
	  call exacr(x,y,pxx,pxy,pyy,px,py,pe)
	  errmax = amax1(errmax,abs((phi(i,j)-pe)))
	end do
      end do
      write(*,201) errmax
  201 format(' maximum error  =  ',e10.3)
      end if
      end

      subroutine cofcr(x,y,cxx,cxy,cyy,cx,cy,ce)
c
c     input pde coefficients at any grid point (x,y) in the solution region
c     (xa.le.x.le.xb,yc.le.y.le.yd) to muh2cr
c
      implicit none
      real x,y,cxx,cxy,cyy,cx,cy,ce
      cxx = 1.+y**2
      cxy = 2.*x*y
      cyy = 1.+x**2
      cx = y
      cy = x
      ce = -(x*y)
      return
      end

      subroutine bndcr(kbdy,xory,alfa,beta,gama,gbdy)
c
c     input mixed "oblique" derivative b.c. to muh2cr
c     at upper y boundary
c
      implicit none
      integer kbdy
      real xory,alfa,beta,gama,gbdy
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmuh2cr/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmuh2cr/xa,xb,yc,yd,tolmax,relmax
      real x,y,pxx,pxy,pyy,px,py,pe
      if (kbdy.eq.4) then
c
c     y=yd boundary (nyd must equal 2 if this code is to be executed).
c     b.c. has the form alfyd(x)*px+betyd(x)*py+gamyd(x)*pe = gbdyd(x)
c     where x = yorx.   alfa,beta,gama,gbdy corresponding to alfyd(x),
c     betyd(x),gamyd(x),gbdyd(y) must be output.
c
      y = yd
      x = xory
      alfa = -x
      beta = 1.+x
      gama = -x
      call exacr(x,y,pxx,pxy,pyy,px,py,pe)
      gbdy = alfa*px + beta*py + gama*pe
      return
      end if
      end

      subroutine exacr(x,y,pxx,pxy,pyy,px,py,pe)
c
c     this subroutine is used for setting an exact solution
c     to test subroutine muh2cr.
c
      implicit none
      real x,y,pxx,pxy,pyy,px,py,pe
      pe = (x*y)**5
      px = 5.*(x*y)**4*y
      py = 5.*(x*y)**4*x
      pxx = 20.*(x*y)**3*y*y
      pxy = 25.*(x*y)**4
      pyy = 20.*(x*y)**3*x*x
      return
      end
