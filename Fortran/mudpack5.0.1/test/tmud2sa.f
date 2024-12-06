c
c     file tmud2sa.f
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
c     test program for the MUDPACK solver mud2
c
c ... see documentation and test files provided in this distribution
c
c ... required MUDPACK files
c
c     mud2sa.f, mudcom.f
c
c
c
c *********************************************************
c *********************************************************
c
c     sample program/test driver for mud2sa
c
c **********************************************************
c **********************************************************
c
c
c     a sample program/test driver for mud2sa is below. it can be
c     executed as an initial test.  output is listed for the test
c     case described.
c
c     test the nonseparable self-adjoint solver mud2sa by approximating
c     the solution to the following elliptic pde in divergence form on
c     the full surface of a sphere of radius one.
c
c          div(sigma(t,p)*grad(u(t,p))) - lambda(t,p)*u(t,p) = f(t,p)
c
c     t and p are colatitude and longitude in radians.  choose a grid
c     as close to 2.5 degrees radians (73 x 145) as the mudpack size
c     constraints allow. multiplying thru by sin(t) puts the pde in the
c     following self-adjoint form suitable for mud2sa:
c
c          d(sin(t)*sigma*du/dt)/dt + d(sigma/sin(t)*du/dp)/dp -
c
c          sin(t)*lambda*u(t,p) = sin(t)*f(t,p).
c
c
c     for testing use the coefficients and exact solution:
c
c          sigma(t,p) = 1.5 + (sin(t)*cos(p))**2
c
c          lambda(t,p) = - sigma(t,p)
c
c          u(t,p) =  [sin(t)*(sin(t)*cos(t)*sin(p)*cos(p))]**2
c
c     (the exact solution is the restriction of the solution u(x,y,z) =
c     (x*y*z)**2 in cartesian coordinates to the surface of the sphere
c     in spherical coordinates).  assume the solution u(t,p) is specified
c     at the poles and is periodic in longitude.
c
c ************************************************************
c     output (32 bit floating point arithmetic)
c ***********************************************************
c
c     mud2sa test
c
c     integer input parameters
c     intl =  0 nta =  1 ntb =  1 npc =  0 npd =  0
c     itp =  5 jpq =  5 iet =  5 jep =  6
c     nt =  81 np = 161 iguess =  0 maxcy =  1
c      method =  2 work space estimate =  260820

c
c     multigrid option parameters
c     kcycle =  2
c     iprer =  2
c     ipost =  1
c     intpol =  3
c
c     floating point input parameters
c     ta =  0.000 tb =  3.142 pc =  0.000 pd =  6.283
c     tolerance (error control) =   0.000E+00
c
c     discretization call to mud2sa  intl =  0
c     ierror =  0 minimum work space =  225200
c
c     approximation call to mud2sa
c     intl =  1 method =  2 iguess =  0 maxcy =  1
c     tolmax =  0.00
c     ierror =  0
c     maximum error  =  0.655E-04
c
c ************************************************************
c ************************************************************
c
      program tmud2sa
      implicit none
      integer iitp,jjpq,iiet,jjep,nnt,nnp,llwork
c
c     set grid size with parameter statements
c
      parameter(iitp = 5, jjpq = 5, iiet = 5, jjep = 6)
      parameter ( nnt = iitp*2**(iiet-1) + 1)
      parameter ( nnp = jjpq*2**(jjep-1) + 1)
c
c     set work space for longitude relaxation only (see mud2sa.d)
c
      parameter (llwork = 20*nnt*nnp)
c
c     dimension solution,right hand side, and work arrays
c
      real u(nnt,nnp),r(nnt,nnp),w(llwork)
c
c     dimension input argument vectors and set up continguous storage
c     for labelling entries
c
      integer iparm(17),mgopt(4)
      real fparm(6)
      integer intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np,
     +                iguess,maxcy,method,nwork,lwork,iter
      common / iprm / intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np,
     +                iguess,maxcy,method,nwork,lwork,iter
      real ta,tb,pc,pd,tolmax,sinat,cosat,sinap,cosap,dlt,dlp
      common / fprm / ta,tb,pc,pd,tolmax,sinat(81),cosat(81),
     +                sinap(161),cosap(161),dlt,dlp
      integer i,j,ierror
      real pi,sint,cost,sinp,cosp,sigma,dsigdt,dsigdp
      real ctt,cpp,ct,cp,ce,tmp,dt,dp,dtt,dpp,ue,ut,up,utt,upp
      real sigt,sigp,xlmb,errm,p,t
c
c     equivlance iparm,fparm with labelled commons iprm,fprm
c
      equivalence (intl,iparm)
      equivalence (ta,fparm)
c
c     declare coefficient and boundary condition input subroutines external
c
      external sigt,sigp,xlmb,bndc
c
c     set input integer parameters
c
c     initialization
c
      intl = 0
c
c     set boundary condition flags: poles specified, longitude periodic
c
      nta = 1
      ntb = 1
      npc = 0
      npd = 0
c
c     set grid sizes from parameter statements
c
      itp = iitp
      jpq = jjpq
      iet = iiet
      jep = jjep
      nt = nnt
      np = nnp
c
c     full multigrid cycling (no initial guess at finest grid)
c
      iguess = 0
c
c     set one multigrid cycle
c
      maxcy = 1
c
c     set line relaxation in the longitude direction
c
      method = 2
c
c     set work space estimate
c
      nwork = llwork
c
c     set multigrid parameters (w(2,1) cycling with fully weighted
c     residual restriction and cubic prolongation).  this can also
c     be set by inputting mgopt(1) = 0 to mud2sa.
c
      mgopt(1) = 2
      mgopt(2) = 2
      mgopt(3) = 1
      mgopt(4) = 3
c
c     set floating point input parameters
c
      pi = 4.0*atan(1.0)
c
c     interval end points (in radians)
c
      ta = 0.0
      tb = pi
      pc = 0.
      pd = pi+pi
c
c     no error control
c
      tolmax = 0.0
c
c     set mesh increments
c
      dlt = (tb-ta)/float(nt-1)
      dlp = (pd-pc)/float(np-1)
c
c     preset sin,cos vectors to save computation on grid points
c
      do i=1,nt
	t = ta+(i-1)*dlt
	sinat(i) = sin(t)
	cosat(i) = cos(t)
      end do
      do j=1,np
	p = pc+(j-1)*dlp
	sinap(j) = sin(p)
	cosap(j) = cos(p)
      end do
c
c     initialize right hand side and solution array except at poles
c
      do  i=2,nt-1
	sint = sinat(i)
	cost = cosat(i)
	do j=1,np
	  sinp = sinap(j)
	  cosp = cosap(j)
	  sigma = 1.5 + (sint*cosp)**2
	  dsigdt = 2.0*cosp*cosp*sint*cost
	  dsigdp = -2.0*sint*sint*cosp*sinp
c
c         compute expanded pde coefficients
c
	  ctt = sint*sigma
	  cpp = sigma/sint
	  ct = sint*dsigdt + cost*sigma
	  cp = dsigdp/sint
	  ce = -sint*sigma
c
c         set intermediate variables for exact solution
c
	  tmp = (sint*cosp*sint*sinp*cost)
	  dt = (2.*sint*cost*cost-sint**3)*(cosp*sinp)
	  dp = (cosp**2-sinp**2)*(sint**2*cost)
	  dtt = (2.*cost**3-4.*cost*sint**2-3.*sint**2*cost)*(cosp*sinp)
	  dpp = (-4.*cosp*sinp)*(sint**2*cost)
c
c         set continuous solution and partial derivatives
c
	  ue = tmp*tmp
	  ut = 2.*tmp*dt
	  up = 2.*tmp*dp
	  utt = 2.*(dt*dt+tmp*dtt)
	  upp = 2.*(dp*dp+tmp*dpp)
c
c         set right hand side of continuous pde on grid
c
	  r(i,j) = ctt*utt+cpp*upp+ct*ut+cp*up+ce*ue
c
c         initialize solution array to zero
c
	  u(i,j) = 0.0
	end do
      end do
c
c     set u, r(unused) at poles
c
      do j=1,np
	u(1,j) = 0.0
	r(1,j) = 0.0
	u(nt,j) = 0.0
	r(nt,j) = 0.0
      end do
c
c     print input parameters
c
      write(*,100)
  100 format(//' mud2sa test' )

      write (*,101) (iparm(i),i=1,15)
  101 format(/' integer input parameters ',
     +/'intl = ',i2,' nta = ',i2,' ntb = ',i2,' npc = ',i2,' npd = ',i2,
     +/' itp = ',i2,' jpq = ',i2,' iet = ',i2,' jep = ',i2
     +/' nt = ',i3,' np = ',i3,' iguess = ',i2,' maxcy = ',i2,
     +/' method = ',i2, ' work space estimate = ',i7)

      write (*,102) (mgopt(i),i=1,4)
  102 format(/' multigrid option parameters ',
     +/' kcycle = ',i2,
     +/' iprer = ',i2,
     +/' ipost = ',i2
     +/' intpol = ',i2)

      write(*,103) (fparm(i),i=1,5)
  103 format(/' floating point input parameters ',
     +/' ta = ',f6.3,' tb = ',f6.3,' pc = ',f6.3,' pd = ',f6.3,
     +/' tolerance (error control) =  ' ,e10.3)
c
c     discretization call to mud2sa
c
      write(*,104) intl
  104 format(/' discretization call to mud2sa ', ' intl = ',i2)
      call mud2sa(iparm,fparm,w,sigt,sigp,xlmb,bndc,r,u,mgopt,ierror)
      write (*,105) ierror,iparm(16)
  105 format(' ierror = ',i2, ' minimum work space = ',i7)
      if (ierror.gt.0) call exit(0)
c
c     aprroximation call to mud2sa
c
      intl = 1
      write(*,106) intl, method, iguess, maxcy, tolmax
  106 format(/' approximation call to mud2sa ',
     +/' intl = ',i2, ' method = ',i2 , ' iguess = ',i2, ' maxcy = ',i2
     +/' tolmax = ',f5.2)
      call mud2sa(iparm,fparm,w,sigt,sigp,xlmb,bndc,r,u,mgopt,ierror)
      write (*,107) ierror
  107 format(' ierror = ',i2 )
      if (ierror.gt.0) call exit(0)
      if (ierror .le. 0) then
c
c     compute and print exact maximum error
c
      errm = 0.0
      do j=1,np
	sinp = sinap(j)
	cosp = cosap(j)
	do i=1,nt
	  sint = sinat(i)
	  cost = cosat(i)
	  ue = (sint*cosp*sint*sinp*cost)**2
	  errm = amax1(errm,abs((u(i,j)-ue)))
	end do
      end do
      write(*,108) errm
  108 format(' maximum error  = ',e10.3 )
      end if
      end

      function sigt(t,p)
c
c     colatitude coefficient (this will be called off grid by mud2sa)
c
      implicit none
      real t,p,sigt,sint,cosp
      sint = sin(t)
      cosp = cos(p)
      sigt = sint*(1.5+(sint*cosp)**2)
      return
      end

      function sigp(t,p)
c
c     longitude coefficient (this will be called off grid by mud2sa)
c
      implicit none
      real t,p,sigp,sint,cosp
      sint = sin(t)
      cosp = cos(p)
c     avoid polar singularities and division by zero
      if (abs(sint) .gt. 1.e-7) then
      sigp = (1.5+(sint*cosp)**2)/sint
      else
      sigp = 1.0
      end if
      return
      end

      function xlmb(t,p)
c
c     zero order coefficient evaluated only on grid points
c
      implicit none
      real t,p,xlmb
      integer i,j
      integer intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np,
     +                iguess,maxcy,method,nwork,lwork,iter
      real ta,tb,pc,pd,tolmax,sinat,cosat,sinap,cosap,dlt,dlp
      common / iprm / intl,nta,ntb,npc,npd,itp,jpq,iet,jep,nt,np,
     +                iguess,maxcy,method,nwork,lwork,iter
      common / fprm / ta,tb,pc,pd,tolmax,sinat(81),cosat(81),
     +                sinap(161),cosap(161),dlt,dlp
      i = int((t-ta)/dlt+0.5)+1
      j = int((p-pc)/dlp+0.5)+1
      xlmb = sinat(i)*(1.5 + (sinat(i)*cosap(j))**2)
      return
      end

      subroutine bndc(kbdy,torp,alfa,gbdy)
c
c     this subroutine must be provided as a dummy argument even though
c     there are no mixed derivative b.c.
c
      implicit none
      integer kbdy
      real torp,alfa,gbdy
      return
      end


