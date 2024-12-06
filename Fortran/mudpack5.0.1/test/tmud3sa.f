c
c     file tmud3sa.f
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
c     test program for the MUDPACK solver mud3sa
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mud3sa.f, mudcom.f, mud3ln.f, mud3pn.f
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud3sa
c
c **********************************************************
c **********************************************************
c
c
c
c     a sample program/test driver for mud3sa is below. it can be
c     executed as an initial test.  the output is listed for the
c     the test case described.
c
c     solve the following helmholtz equation (divergence form) in
c     spherical coordinates.
c
c          div(grad(u(r,t,p)) - u(r,t,p) = f(r,t,p)  (original pde)
c
c     after multiplying the left and right hand sides by r*r*sin(t) this
c     becomes
c
c          r*r*sin(t)*div(grad(u(r,t,p)) - r*r*sin(t)*u(r,t,p)
c
c          = r*r*sin(t)*f(r,t,p)
c
c     or (expanded form)
c
c
c          d(r*r*sin(t)*du/dr)/dr + d(sin(t)*du/dt)/dt + d(1/sin(t)*du/dp)/dp
c
c          - r*r*sin(t)*u(r,t,p)  =  r*r*sin(t)*f(r,t,p).
c
c
c     this form is suitable for using mud3sa.  the solution region is:
c
c          (pi = 4.*atan(1.0) radians)
c
c          0.5  <  r  <  1.0       r is the "radius" coordinate
c
c          pi/4  <  t  <  3*pi/4   t is the colatitude coordinate
c
c          0  <  p  <  pi+pi       p is the longitude coordinate
c
c     asssume the solution is specified at the r,t boundaries and is periodic
c     in p.  use line relaxation in the r direction and choose a solution
c     grid as close to 50 x 30 x 60 as the mudpack size constraints allow.
c     the exact solution
c
c          u(r,t,p) = (r*sin(t)*sin(p)*r*sin(t)*cos(p)*r*cos(t))**2
c
c                   = (x*y*z)**2   (in cartesian coordinates)
c
c     is used for testing purposes. one full multigrid cycle (no initial
c     guess) with the default multigrid options is executed.  error control
c     is not used.
c
c ******************************************************
c     output (32 bit floating point arithmetic
c *******************************************************
c
c     mud3sa test
c
c     input arguments
c     intl =  0 nra =  1 nrb =  1 ntc =  1 ntd =  1
c     npe =  0 npf =  0
c     irp =  3 jtq =  2 kpr =  2
c     ier =  5 jet =  5 kep =  6
c     nr =  49 nt =  33 np =  65 iguess =  0 maxcy =  1
c     method =  1 work space length input = 1598164
c     ra =  0.50 rb =  1.00
c     tc =  0.79 td =  2.36
c     pe =  0.00 pf =  6.28
c     tolmax =  0.000E+00
c
c     multigrid options
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     discretization call to mud3sa intl =  0
c     ierror =  0 minimum work space = 1598164
c
c     approximation call to mud3sa
c     intl =  1 method =  1 iguess =  0 maxcy =  1
c     ierror =  0
c     maximum error  =   0.369E-03
c
c ************************************************************
c     end of output
c ************************************************************
c
      program tmud3sa
      implicit none
c
c     set grid sizes with parameter statements
c
      integer iirp,jjtq,kkpr,iier,jjet,kkep,nnr,nnt,nnp
      parameter(iirp=3, jjtq=2, kkpr=2, iier=5, jjet=5, kkep=6)
      parameter (nnr = iirp*2**(iier-1)+1)
      parameter (nnt = jjtq*2**(jjet-1)+1)
      parameter (nnp = kkpr*2**(kkep-1)+1)
c
c     set predetermined minimum work space required
c
      integer llwork
      parameter (llwork = 1598164)
c
c     dimension solution,right hand side, and work arrays
c
      real uso(nnr,nnt,nnp),rhs(nnr,nnt,nnp),work(llwork)
c
c     put integer and floating point parameter names in contiguous
c     storeage for labelling in equivalenced vectors iprm,fprm
c
      real fprm(8)
      integer iprm(23),mgopt(5)
      integer intl,nra,nrb,ntc,ntd,npe,npf,irp,jtq,kpr,ier,jet,
     +              kep,nr,nt,np,iguess,maxcy,method,method2,nwork,
     +              lwrkqd,itero
      common/itmud3/intl,nra,nrb,ntc,ntd,npe,npf,irp,jtq,kpr,ier,jet,
     +              kep,nr,nt,np,iguess,maxcy,method,method2,nwork,
     +              lwrkqd,itero
      real ra,rb,tc,td,pe,pf,tolmax,relmax
      common/ftmud3/ra,rb,tc,td,pe,pf,tolmax,relmax
      real sinat,cosat,sinap,cosap,dlr,dlt,dlp
      integer i,j,k,ierror
      real urr,utt,upp,ur,ut,up,ue,r,r6,t,p,pi
      real st,ct,sp,cp,ep,et,sint,cost,errm
c
c     vectors for saving sin,cos on lat-long points
c
      common/sincos/sinat(33),cosat(33),sinap(65),cosap(65),dlr,dlt,dlp
      equivalence(intl,iprm)
      equivalence(ra,fprm)
c
c     declare coefficient input functions external
c
      real lam
      external sigr,sigt,sigp,bndc,lam
c
c
c     set input integer arguments
c
      intl = 0
c
c     set boundary condition flags
c
      nra = 1
      nrb = 1
      ntc = 1
      ntd = 1
      npe = 0
      npf = 0
c
c     set grid sizes from arguments statements
c
      irp = iirp
      jtq = jjtq
      kpr = kkpr
      ier = iier
      jet = jjet
      kep = kkep
      nr = nnr
      nt = nnt
      np = nnp
c
c     set one cycle limit
c
      maxcy = 1
c
c     set work space length approximation from parameter statement
c
      nwork = llwork
c
c     set line relaxation in the r direction
c
      method = 1
      method2 = 0
c
c     flag no initial guess which forces full multigrid cycling
c
      iguess = 0
c
c     set ends of solution "cube" in (r,t,p) space  (actually the sphere)
c
      pi = 4.*atan(1.)
      ra = 0.5
      rb = 1.0
      tc = 0.25*pi
      td = 0.75*pi
      pe = 0.0
      pf = pi+pi
c
c     set mesh increments
c
      dlr = (rb-ra)/float(nr-1)
      dlt = (td-tc)/float(nt-1)
      dlp = (pf-pe)/float(np-1)
c
c     preset sin,cos on lat-long grid to save computation
c
      do k=1,np
	p = pe +(k-1)*dlp
	sinap(k) = sin(p)
	cosap(k) = cos(p)
      end do
      do j=1,nt
	t = tc+(j-1)*dlt
	sinat(j) = sin(t)
	cosat(j) = cos(t)
      end do
c
c     set for no error control
c
      tolmax = 0.0
c
c     set right hand side in rhs and initialize solution to zero
c
      do j=1,nt
	t = tc+(j-1)*dlt
	sint = sinat(j)
	cost = cosat(j)
	do k=1,np
	  p = pe+(k-1)*dlp
	  do i=1,nr
	    r = ra+float(i-1)*dlr
	    call exact(r,t,p,urr,utt,upp,ur,ut,up,ue)
	    rhs(i,j,k) = r*sint*(r*urr+2.*ur)+sint*utt+cost*ut+upp/sint
     +                 - lam(r,t,p)*ue
	    uso(i,j,k) = 0.0
	  end do
	end do
      end do
c
c     set specified values in uso at r and t boundaries
c
c     r = ra
      r6 = ra**6
      do k=1,np
	sp = sinap(k)
	cp = cosap(k)
	ep = (sp*cp)**2
	do j=1,nt
	  st = sinat(j)
	  ct = cosat(j)
	  et = st*st*(st*ct)**2
	  uso(1,j,k) = r6*et*ep
	end do
      end do
c
c     r = rb
      r6 = rb**6
      do k=1,np
	sp = sinap(k)
	cp = cosap(k)
	ep = (sp*cp)**2
	do j=1,nt
	  st = sinat(j)
	  ct = cosat(j)
	  et = st*st*(st*ct)**2
	  uso(nr,j,k) = r6*et*ep
	end do
      end do

c     t = tc
      st = sinat(1)
      ct = cosat(1)
      et = st*st*(st*ct)**2
      do k=1,np
	sp = sinap(k)
	cp = cosap(k)
	ep = (cp*sp)**2
	do i=1,nr
	  r6 = (ra+float(i-1)*dlr)**6
	  uso(i,1,k) = r6*et*ep
	end do
      end do

c     t = td
      st = sinat(nt)
      ct = cosat(nt)
      et = st*st*(st*ct)**2
      do k=1,np
	sp = sinap(k)
	cp = cosap(k)
	ep = (cp*sp)**2
	do i=1,nr
	  r6 = (ra+float(i-1)*dlr)**6
	  uso(i,nt,k) = r6*et*ep
	end do
      end do
c
c     set default multigrid options
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     write input arguments
c
      write(*,50)
   50 format(//' mud3sa test')
      write(*,100)intl,nra,nrb,ntc,ntd,npe,npf,irp,jtq,kpr,ier,jet,kep,
     +            nr,nt,np,iguess,maxcy,method,nwork,ra,rb,tc,td,pe,pf,
     +            tolmax,mgopt(1),mgopt(2),mgopt(3),mgopt(4)
  100 format(/' input arguments',
     +/'intl = ',i2,' nra = ',i2,' nrb = ',i2,' ntc = ',i2,' ntd = ',i2,
     +/' npe = ',i2, ' npf = ',i2,
     +/' irp = ',i2,' jtq = ',i2,' kpr = ',i2,
     +/' ier = ',i2, ' jet = ',i2, ' kep = ',i2,
     +/' nr = ',i3,' nt = ',i3,' np = ',i3, ' iguess = ',i2,
     +/' maxcy = ',i2,
     +/' method = ',i2, ' work space length input = ',i7,
     +/' ra = ',f5.2,' rb = ',f5.2,
     +/' tc = ',f5.2,' td = ',f5.2,
     +/' pe = ',f5.2,' pf = ',f5.2,
     +/' tolmax = ' ,e10.3
     +//' multigrid options '
     +/' kcycle = ',i2
     +/' iprer = ',i2
     +/' ipost = ',i2
     +/' intpol = ',i2 )
      write(*,104) intl
  104 format(/' discretization call to mud3sa', ' intl = ', i2)
      call mud3sa(iprm,fprm,work,sigr,sigt,sigp,lam,bndc,rhs,uso,
     +            mgopt,ierror)
      write (*,105) ierror,iprm(22)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     attempt solution
c
      intl = 1
      write(*,106) intl,method,iguess,maxcy
  106 format(/' approximation call to mud3sa ',
     +/' intl = ',i2, ' method = ',i2,' iguess = ',i2, ' maxcy = ',i2)
      call mud3sa(iprm,fprm,work,sigr,sigt,sigp,lam,bndc,rhs,uso,
     +            mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2)
      if (ierror.gt.0) call exit(0)
c
c     compute maximum error
c
      errm = 0.0
      do k=1,np
	sp = sinap(k)
	cp = cosap(k)
	ep = (cp*sp)**2
	do j=1,nt
	  st = sinat(j)
	  ct = cosat(j)
	  et = st*st*(st*ct)**2
	  do i=1,nr
	    r6 = (ra+float(i-1)*dlr)**6
c     exact continuous solution
	    ue = r6*et*ep
	    errm = amax1(errm,abs(uso(i,j,k)-ue))
	  end do
	end do
      end do
      write(*,108) errm
  108 format(' maximum error  =  ',e10.3)
      end

      function sigr(r,t,p)
c
c     coefficient for r derivative (mud3sa will call sigr off grid)
c
      implicit none
      real sigr,r,t,p
      sigr = r*r*sin(t)
      return
      end

      function sigt(r,t,p)
c
c     coefficient for theta derivative (mud3sa will call sigt off grid)
c
      implicit none
      real sigt,r,t,p
      sigt = sin(t)
      return
      end

      function sigp(r,t,p)
c
c     coefficient for phi derivative (mud3sa will call sigp off grid)
c
      implicit none
      real sigp,r,t,p
      sigp = 1.0/sin(t)
      return
      end

      function lam(r,t,p)
c
c     input zero order coefficient in self adjoint pde at (r,t,p) to mud3sa
c     (only grid values needed)
c
      implicit none
      real lam,r,t,p
      integer j
      integer intl,nra,nrb,ntc,ntd,npe,npf,irp,jtq,kpr,ier,jet,
     +              kep,nr,nt,np,iguess,maxcy,method,method2,nwork,
     +              lwrkqd,itero
      common/itmud3/intl,nra,nrb,ntc,ntd,npe,npf,irp,jtq,kpr,ier,jet,
     +              kep,nr,nt,np,iguess,maxcy,method,method2,nwork,
     +              lwrkqd,itero
      real ra,rb,tc,td,pe,pf,tolmax,relmax
      real sinat,cosat,sinap,cosap,dlr,dlt,dlp
      common/ftmud3/ra,rb,tc,td,pe,pf,tolmax,relmax
      common/sincos/sinat(33),cosat(33),sinap(65),cosap(65),dlr,dlt,dlp
      j = int((t-tc)/dlt+0.5)+1
      lam = r*r*sinat(j)
      return
      end
c
c
      subroutine bndc(kbdy,rort,torp,alfa,gbdy)
c
c     a dummy routine since there is no mixed derivative b.c. for this problem
c
      return
      end

      subroutine exact(r,t,p,urr,utt,upp,ur,ut,up,ue)
c
c     set exact solution taken from u(x,y,z) = (x*y*z)**2 transformed
c     with spherical coordinates
c
      common/itmud3/intl,nra,nrb,ntc,ntd,npe,npf,irp,jtq,kpr,ier,jet,
     +              kep,nr,nt,np,iguess,maxcy,method,method2,nwork,
     +              lwrkqd,itero
      common/ftmud3/ra,rb,tc,td,pe,pf,tolmax,relmax
      common/sincos/sinat(33),cosat(33),sinap(65),cosap(65),dlr,dlt,dlp
c
c     set subscripts for current lat-lon grid point
c
      j = int((t-tc)/dlt+0.5)+1
      k = int((p-pe)/dlp+0.5)+1
c
c     set sin, cos from pre-computed vectors
c
      st = sinat(j)
      ct = cosat(j)
      sp = sinap(k)
      cp = cosap(k)
c
c     set intermediate quantities
c
      r6 = r**6
      ep = (cp*sp)**2
      dep = 2.*(cp*sp)*(cp*cp-sp*sp)
      ddep = 2.*((cp**2-sp**2)**2-4.*(sp*cp)**2)
      et = st*st*(st*ct)**2
      det = 2.*(2.*(st*ct)**3 -ct*st**5)
      ddet = 12.*(ct*st)**2*(ct**2-st**2)+2.*st**4*(st**2-4.*ct**2)
c
c     set exact values of continuous solution and its partial derivatives
c
      ue = r6*et*ep
      ur = 6.*r**5*et*ep
      urr = 30.*r**4*et*ep
      ut = r6*det*ep
      utt = r6*ddet*ep
      up = r6*et*dep
      upp = r6*et*ddep
      return
      end

