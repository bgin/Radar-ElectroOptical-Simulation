c
c     file tmud2sp.f
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
c     test program for the mudpack solver mud2sp
c
c ... see documentation and test files provided in this distribution
c
c ... required mudpack files
c
c     mud2sp.f, mudcom.f
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud2sp
c
c **********************************************************
c **********************************************************
c
c
c     a sample program/test driver for mud2sp is listed below.  it
c     can be executed as an initial test.  the output is listed for
c     listed for the test case described.
c
c     test the driver below by solving the separable elliptic pde
c
c     (1.+x**2)*pxx + exp(1.-y)*(pyy-py) - (x+y)*pe = r(x,y)
c
c     on the unit square with specified boundary conditions at
c     xb = 1.0, yc = 0.0 and mixed boundary conditions
c
c          dp/dx - pe(0.0,y) = ga(y)  (at x = 0.0)
c
c          dp/dy + pe(x,1.0) = gd(x)  (at y = 1.0)
c
c     use point relaxation and choose a grid as close to 60 x 50
c     as the grid size constraints allow.  use the exact solution
c
c          pe(x,y) = (x**3+y**3+1.0)/3
c
c     for testing.
c
c *************************************************************
c     output (32 bit floating point arithmetic)
c *************************************************************
c
c      mud2sp test
c
c     integer input arguments
c     intl =  0 nxa =  2 nxb =  1 nyc =  1 nyd =  2
c     ixp =  2 jyq =  3 iex =  6 jey =  5
c     nx =  65 ny =  49 iguess =  0 maxcy =  1
c     method =  0 work space estimate =   17065
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
c     discretization call to mud2sp intl =  0
c     ierror =  0 minimum work space =   13264
c
c     approximation call to mud2sp
c     intl =  1 method =  0 iguess =  0
c     ierror =  0
c     maximum error  =   0.520E-04
c
c ***************************************************************
c      end of output
c ***************************************************************
c
      program tmud2sp
      implicit none
      integer iixp,jjyq,iiex,jjey,nnx,nny,llwork
c
c     set grid sizes with parameter statements
c
      parameter (iixp = 2 , jjyq = 3 , iiex = 6, jjey = 5)
      parameter (nnx=iixp*2**(iiex-1)+1, nny=jjyq*2**(jjey-1)+1)
c
c     set work space length approximation for point relaxation
c     see (mud2sp.d)
c
      parameter (llwork=5*(nnx*nny+2*(nnx+nny)))
      real phi(nnx,nny),rhs(nnx,nny),work(llwork)
c
c     put integer and floating point argument names in contiguous
c     storeage for labelling in vectors iprm,fprm
c
      integer iprm(16),mgopt(4)
      real fprm(6)
      integer intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      common/itmud2sp/intl,nxa,nxb,nyc,nyd,ixp,jyq,iex,jey,nx,ny,
     +              iguess,maxcy,method,nwork,lwrkqd,itero
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
      equivalence(intl,iprm)
      equivalence(xa,fprm)
      integer i,j,ierror
      real dlx,dly,x,y,cxx,cyy,cx,cy,ce,pxx,pyy,px,py,pe,errmax
      real cex,cey
c
c     declare coefficient and boundary condition input subroutines external
c
      external cofx,cofy,bndsp
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
	call cofx(x,cxx,cx,cex)
	do j=1,ny
	  y = yc+float(j-1)*dly
	  call cofy(y,cyy,cy,cey)
	  ce = cex+cey
	  call exact(x,y,pxx,pyy,px,py,pe)
	  rhs(i,j) = cxx*pxx+cyy*pyy+cx*px+cy*py+ce*pe
	  phi(i,j) = 0.0
	end do
      end do
c
c     set specified boundaries in phi
c
      x = xb
      do j=1,ny
	y = yc+float(j-1)*dly
	call exact(x,y,pxx,pyy,px,py,pe)
	phi(nx,j) = pe
      end do
      y = yc
      do i=1,nx
	x = xa+float(i-1)*dlx
	call exact(x,y,pxx,pyy,px,py,pe)
	phi(i,1) = pe
      end do
      write(*,100)
  100 format(//' mud2sp test ')
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
  104 format(/' discretization call to mud2sp', ' intl = ', i2)
      call mud2sp(iprm,fprm,work,cofx,cofy,bndsp,rhs,phi,mgopt,ierror)
c
c     print error parameter and minimum work space requirement
c
      write (*,200) ierror,iprm(16)
  200 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
      write(*,106) intl,method,iguess
  106 format(/' approximation call to mud2sp',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2)
      call mud2sp(iprm,fprm,work,cofx,cofy,bndsp,rhs,phi,mgopt,ierror)
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
	  call exact(x,y,pxx,pyy,px,py,pe)
	  errmax = amax1(errmax,abs((phi(i,j)-pe)))
	end do
      end do
      write(*,201) errmax
  201 format(' maximum error  =  ',e10.3)
      end if
      end

      subroutine cofx(x,cxx,cx,cex)
c
c     input x dependent coefficients
c
      implicit none
      real x,cxx,cx,cex
      cxx = 1.0+x*x
      cx = 0.0
      cex = -x
      return
      end

      subroutine cofy(y,cyy,cy,cey)
c
c     input y dependent coefficients
c
      implicit none
      real y,cyy,cy,cey
      cyy = exp(1.0-y)
      cy = -cyy
      cey = -y
      return
      end

      subroutine bndsp(kbdy,xory,alfa,gbdy)
c
c     input mixed derivative b.c. to mud2sp
c
      implicit none
      integer kbdy
      real xory,alfa,gbdy,x,y,pe,px,py,pxx,pyy
      real xa,xb,yc,yd,tolmax,relmax
      common/ftmud2sp/xa,xb,yc,yd,tolmax,relmax
      if (kbdy.eq.1) then  ! x=xa boundary
	y = xory
	x = xa
	call exact(x,y,pxx,pyy,px,py,pe)
	alfa = -1.0
	gbdy = px + alfa*pe
	return
      end if
      if (kbdy.eq.4) then  ! y=yd boundary
	y = yd
	x = xory
	call exact(x,y,pxx,pyy,px,py,pe)
	alfa = 1.0
	gbdy = py + alfa*pe
	return
      end if
      end

      subroutine exact(x,y,pxx,pyy,px,py,pe)
c
c     set an exact solution for testing mud2sp
c
      implicit none
      real x,y,pxx,pyy,px,py,pe
      pe = (x**3+y**3+1.0)/3.0
      px = x*x
      py = y*y
      pxx = x+x
      pyy = y+y
      return
      end

