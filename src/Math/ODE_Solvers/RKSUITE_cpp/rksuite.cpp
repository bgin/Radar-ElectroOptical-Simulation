#include "rksuite.h"
#include <stdlib.h>


inline double sign(double a, double b)
{
  return (b >= 0.0 ? fabs(a) : -fabs(a));
}


inline double max(double a, double b)
{
  return (a >= b ? a : b);
}

inline double min(double a, double b)
{
  return (a <= b ? a : b);
}

inline int sign(int a, int b)
{
  return (b >= 0 ? abs(a) : -abs(a));
}

inline int max(int a, int b)
{
  return (a >= b ? a : b);
}

inline int min(int a, int b)
{
  return (a <= b ? a : b);
}


void RKSUITE::setup(int neq, double tstart, double ystart[], double tend,
            double tol, double thres[], int method, char* task,
            bool errass, double hstart, double work[], int lenwrk,
            bool mesage)
{
//      SUBROUTINE SETUP(NEQ,TSTART,YSTART,TEND,TOL,THRES,METHOD,TASK,
//     &                 ERRASS,HSTART,WORK,LENWRK,MESAGE)
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code SETUP and how it is used in
//C  conjunction with UT or CT to solve initial value problems, you should study
//C  the document file rksuite.doc carefully before attempting to use the code.
//C  The following "Brief Reminder" is intended only to remind you of the
//C  meaning, type, and size requirements of the arguments.
//C
//C  The environmental parameters OUTCH, MCHEPS, and DWARF are used in the
//C  following description.  To find out their values
//C
//C       CALL ENVIRN(OUTCH,MCHEPS,DWARF)
//C
//C  INPUT VARIABLES
//C
//C     NEQ       - INTEGER
//C                 The number of differential equations in the system.
//C                 Constraint: NEQ >= 1
//C     TSTART    - double PRECISION
//C                 The initial value of the independent variable.
//C     YSTART(*) - double PRECISION array of length NEQ
//C                 The vector of initial values of the solution components.
//C     TEND      - double PRECISION
//C                 The integration proceeds from TSTART in the direction of
//C                 TEND. You cannot go past TEND.
//C                 Constraint: TEND must be clearly distinguishable from TSTART
//C                 in the precision available.
//C     TOL       - double PRECISION
//C                 The relative error tolerance.
//C                 Constraint: 0.01D0 >= TOL >= 10*MCHEPS
//C     THRES(*)  - double PRECISION array of length NEQ
//C                 THRES(L) is the threshold for the Ith solution component.
//C                 Constraint: THRES(L) >= SQRT(DWARF)
//C     METHOD    - INTEGER
//C                 Specifies which Runge-Kutta pair is to be used.
//C                  = 1 - use the (2,3) pair
//C                  = 2 - use the (4,5) pair
//C                  = 3 - use the (7,8) pair
//C     TASK      - CHARACTER*(*)
//C                 Only the first character of TASK is significant.
//C                 TASK(1:1) = 'U' or 'u' - UT is to be used
//C                           = 'C' or 'c' - CT is to be used
//C                 Constraint: TASK(1:1) = 'U'or 'u' or'C' or 'c'
//C     ERRASS    - LOGICAL
//C                 = .FALSE. - do not attempt to assess the true error.
//C                 = .TRUE.  - assess the true error. Costs roughly twice
//C                             as much as the integration with METHODs 2 and
//C                             3, and three times with METHOD = 1.
//C     HSTART    - double PRECISION
//C                 0.0D0     - select automatically the first step size.
//C                 non-zero  - try HSTART for the first step.
//C
//C  WORKSPACE
//C
//C     WORK(*) - double PRECISION array of length LENWRK
//C               Do not alter the contents of this array after calling SETUP.
//C
//C  INPUT VARIABLES
//C
//C     LENWRK  - INTEGER
//C               Length of WORK(*): How big LENWRK must be depends
//C               on the task and how it is to be solved.
//C
//C               LENWRK = 32*NEQ is sufficient for all cases.
//C
//C               If storage is a problem, the least storage possible
//C               in the various cases is:
//C
//C                 If TASK = 'U' or 'u', then
//C                   if ERRASS = .FALSE. and
//C                     METHOD = 1, LENWRK must be at least 10*NEQ
//C                            = 2                          20*NEQ
//C                            = 3                          16*NEQ
//C                   if ERRASS = .TRUE. and
//C                     METHOD = 1, LENWRK must be at least 15*NEQ
//C                            = 2                          32*NEQ
//C                            = 3                          21*NEQ
//C
//C                 If TASK = 'C' or 'c', then
//C                   if ERRASS = .FALSE. and
//C                     METHOD = 1, LENWRK must be at least 10*NEQ
//C                            = 2                          14*NEQ
//C                            = 3                          16*NEQ
//C                   if ERRASS = .TRUE. and
//C                     METHOD = 1, LENWRK must be at least 15*NEQ
//C                            = 2                          26*NEQ
//C                            = 3                          21*NEQ
//C
//C                 Warning:  To exploit the interpolation capability
//C                 of METHODs 1 and 2, you have to call INTRP.  This
//C                 subroutine requires working storage in addition to
//C                 that specified here.
//C
//C     MESAGE    - LOGICAL
//C                 Specifies whether you want informative messages written to
//C                 the standard output channel OUTCH.
//C                 = .TRUE.   - provide messages
//C                 = .FALSE.  - do not provide messages
//C
//C  In the event of a "catastrophic" failure to call SETUP correctly, the
//C  nature of the catastrophe is reported on the standard output channel,
//C  regardless of the value of MESAGE.  Unless special provision was made
//C  in advance (see rksuite.doc), the computation then comes to a STOP.
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block for General Workspace Pointers ..
//      INTEGER           PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      COMMON /RKCOM3/   PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      SAVE   /RKCOM3/
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Global Error Assessment ..
//      double PRECISION  MAXERR, LOCMAX
//      INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
//     &                  PRZYNU
//      LOGICAL           ERASON, ERASFL
//      COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
//     &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
//      SAVE   /RKCOM6/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Common Block for Integrator Options ..
//      LOGICAL           MSG, UTASK
//      COMMON /RKCOM8/   MSG, UTASK
//      SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      CHARACTER*6       SRNAME
//      PARAMETER         (SRNAME='SETUP')
//      INTEGER           MINUS1
//      LOGICAL           TELL
//      PARAMETER         (MINUS1=-1,TELL=.FALSE.)
//      double PRECISION  ONE, ZERO, PT01
//      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,PT01=0.01D+0)
   const char srname[] = "SETUP";
   const int minus1 = -1;
   const bool tell = false;
   const double one = 1.0;
   const double zero = 0.0;
   const double pt01 = 0.01;
//C     .. Local Scalars ..
   double hmin;
   int flag, freepr, ier, l, lintpl, lreq, nrec, vecstg;
   bool legalt, reqstg;
   char task1;
//C     .. External Subroutines ..
//      EXTERNAL          CONST, MCONST, RKMSG, RKSIT
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX, SIGN
//C     .. Executable Statements ..
//C
//C  Clear previous flag values of subprograms in the suite.
//C
   ier = minus1;
   rksit(tell,srname,ier);
//C
   ier = 1;
   nrec = 0;
//C
//C  Fetch output channel and machine constants; initialise common
//C  block /RKCOM7/
//C
   mconst(method);
//C
//C  Check for valid input of trivial arguments
   task1 = task[0];
   legalt = task1 == 'U' || task1 == 'u' ||
           task1 == 'C' || task1 == 'c';
   if (!legalt) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** You have set the first character of");
      sprintf(&rkcom9.rec[1][0]," ** TASK to be '%c'. It must be one of",task1);
      sprintf(&rkcom9.rec[2][0]," ** 'U','u','C' or 'c'.");
   }
   else if (neq < 1) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have set NEQ = %6d which is less than 1.",neq);
   }
   else if (method < 1 || method > 3) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have set METHOD = %6d which is not 1, 2, or 3.",method);
   }
   else if (tstart == tend) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have set TSTART = TEND = %13.5lf.",tstart);
   }
   else if ((tol > pt01) || (tol < rkcom7.rndoff)) {
      ier = 911;
      nrec = 2;
      sprintf(&rkcom9.rec[0][0]," ** You have set TOL = %13.5lf which is not permitted. The",tol);
      sprintf(&rkcom9.rec[1][0]," ** range of permitted values is (%13.5lf,0.01D0).",rkcom7.rndoff);
   }
   else {
      l = 0;
      label20:
      if (thres[l] < rkcom7.tiny) {
         ier = 911;
         nrec = 2;
         sprintf(&rkcom9.rec[0][0]," ** You have set THRES(%6d) to be %13.5lf which is",l,thres[l]);
         sprintf(&rkcom9.rec[1][0]," ** less than the permitted minimum,%13.5lf.",rkcom7.tiny);
      }
      l++;
      if (ier != 911 && l < neq) goto label20;
   }
//C
//C  Return if error detected
//C
   if (ier != 1) goto label80;
//C
//C  Set formula definitions and characteristics by means of arguments
//C  in the call list and COMMON blocks /RKCOM4/ and /RKCOM5/
//C
   rkconst(method,vecstg,reqstg,lintpl);
//C
//C  Set options in /RKCOM8/
   rkcom8.utask = task1 == 'U' || task1 == 'u';
   rkcom8.msg = mesage;
//C
//C  Initialise problem status in /RKCOM1/ and /RKCOM2/
   rkcom1.neqn = neq;
   rkcom1.tstrt = tstart;
   rkcom1.tnd = tend;
   rkcom2.t = tstart;
   rkcom2.told = tstart;
   rkcom1.dir = sign(one,tend-tstart);
//C
//C  In CT the first step taken will have magnitude H.  If HSTRT = ABS(HSTART)
//C  is not equal to zero, H = HSTRT.  If HSTRT is equal to zero, the code is
//C  to find an on-scale initial step size H.  To start this process, H is set
//C  here to an upper bound on the first step size that reflects the scale of
//C  the independent variable.  UT has some additional information, namely the
//C  first output point, that is used to refine this bound in UT when UTASK
//C  is .TRUE..  If HSTRT is not zero, but it is either too big or too small,
//C  the input HSTART is ignored and HSTRT is set to zero to activate the
//C  automatic determination of an on-scale initial step size.
//C
   rkcom1.hstrt = fabs(hstart);
   hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(tstart),fabs(tend)));
   if (rkcom1.hstrt > fabs(tend-tstart) || rkcom1.hstrt < hmin) rkcom1.hstrt = zero;
   if (rkcom1.hstrt == zero) {
     rkcom2.h = max(fabs(tend-tstart)/rkcom5.rs3,hmin);
   }
   else {
     rkcom2.h = rkcom1.hstrt;
   }
   rkcom2.hold = zero;
   rkcom1.tolr = tol;
   rkcom2.nfcn = 0;
   rkcom2.svnfcn = 0;
   rkcom2.okstp = 0;
   rkcom2.flstp = 0;
   rkcom2.first = true;
   rkcom2.last = false;
//C
//C  WORK(*) is partioned into a number of arrays using pointers. These
//C  pointers are set in /RKCOM3/.
   rkcom3.prthrs = 0;
//C                           the threshold values
   rkcom3.prerst = rkcom3.prthrs + neq;
//C                           the error estimates
   rkcom3.prwt = rkcom3.prerst + neq;
//C                           the weights used in the local error test
   rkcom3.pryold = rkcom3.prwt + neq;
//C                           the previous value of the solution
   rkcom3.prscr = rkcom3.pryold + neq;
//C                           scratch array used for the higher order
//C                           approximate solution and for the previous
//C                           value of the derivative of the solution
   rkcom3.pry = rkcom3.prscr + neq;
//C                           the dependent variables
   rkcom3.pryp = rkcom3.pry + neq;
//C                           the derivatives
   rkcom3.prstgs = rkcom3.pryp + neq;
//C                           intermediate stages held in an internal
//C                           array STAGES(NEQ,VECSTG)
//C
   freepr = rkcom3.prstgs + vecstg*neq;
//C
//C  Allocate storage for interpolation if the TASK = 'U' or 'u' was
//C  specified. INTP and LINTPL returned by CONST indicate whether there
//C  is an interpolation scheme associated with the pair and how much
//C  storage is required.
//C
   rkcom3.printp = 1;
   rkcom3.lnintp = 1;
   if (rkcom8.utask) {
      if (rkcom4.intp) {
         rkcom3.lnintp = lintpl*neq;
         if (reqstg) {
            rkcom3.printp = freepr;
            freepr = rkcom3.printp + rkcom3.lnintp;
         }
         else {
            rkcom3.printp = rkcom3.prstgs;
            freepr = max(rkcom3.printp+vecstg*neq,rkcom3.printp+rkcom3.lnintp);
         }
      }
   }
//C
//C  Initialise state and allocate storage for global error assessment
//C  using /RKCOM6/
   rkcom6.gnfcn = 0;
   rkcom6.maxerr = zero;
   rkcom6.locmax = tstart;
   rkcom6.erason = errass;
   rkcom6.erasfl = false;
   if (errass) {
//C
//C  Storage is required for the stages of a secondary integration. The
//C  stages of the primary intergration can only be overwritten in the
//C  cases where there is no interpolant or the interpolant does not
//C  require information about the stages (e.g. METHOD 3 and METHOD 1,
//C  respectively).
     if (!reqstg) {
          rkcom6.przstg = rkcom3.prstgs;
      }
      else {
         rkcom6.przstg = freepr;
         freepr = rkcom6.przstg + vecstg*neq;
      }
      rkcom6.przy = freepr;
      rkcom6.przyp = rkcom6.przy + neq;
      rkcom6.przers = rkcom6.przyp + neq;
      rkcom6.przerr = rkcom6.przers + neq;
      rkcom6.przynu = rkcom6.przerr + neq;
      freepr = rkcom6.przynu + neq;
   }
   else {
      rkcom6.przstg = 1;
      rkcom6.przy = 1;
      rkcom6.przyp = 1;
      rkcom6.przers = 1;
      rkcom6.przerr = 1;
      rkcom6.przynu = 1;
   }
//C
   lreq = freepr - 1;
//C
//C  Check for enough workspace and suitable range of integration
//C
   if (lenwrk < lreq) {
      ier = 911;
      nrec = 2;
      sprintf(&rkcom9.rec[0][0]," ** You have not supplied enough workspace. You gave LENWRK");
      sprintf(&rkcom9.rec[1][0]," ** as%6d, but it must be at least %6d.",lenwrk,lreq);
   }
   else {
      hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(tstart),fabs(tend)));
      if (fabs(tend-tstart) < hmin) {
         ier = 911;
         nrec = 4;
         sprintf(&rkcom9.rec[0][0]," ** You have set values for TEND and TSTART that are not");
         sprintf(&rkcom9.rec[1][0]," ** clearly distinguishable for the method and the precision");
         sprintf(&rkcom9.rec[2][0]," ** of the computer being used. ABS(TEND-TSTART) is %13.5lf",fabs(tend-tstart));
         sprintf(&rkcom9.rec[3][0]," ** but should be at least %13.5lf.",hmin);
      }
   }
//C
//C  Return if error detected
//C
   if (ier != 1) goto label80;
//C
//C  Initialize elements of the workspace
   for (l = 0; l < neq; l++) {
      work[rkcom3.prthrs+l] = thres[l];
      work[rkcom3.pry+l] = ystart[l];
   }
//C
//C  Initialize the global error to zero when ERRASS = .TRUE.
   if (errass) {
      for (l = 0; l < neq; l++) {
         work[rkcom6.przerr+l] = zero;
      }
   }
//C
   label80:
//C
   rkmsg(ier,srname,nrec,flag);
//C
   return;
}

void RKSUITE::ut(void (*f)(double, double*, double*), double twant, double& tgot,
     double ygot[], double ypgot[], double ymax[], double work[], int& uflag)
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code UT and how it is used in
//C  conjunction with SETUP to solve initial value problems, you should study
//C  the document file rksuite.doc carefully before proceeding further.  The
//C  following "Brief Reminder" is intended only to remind you of the meaning,
//C  type, and size requirements of the arguments.
//C
//C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM:
//C
//C     F         - name of the subroutine for evaluating the differential
//C                 equations.
//C
//C  The subroutine F must have the form
//C
//C  SUBROUTINE F(T,Y,YP)
//C  double PRECISION T,Y(*),YP(*)
//C     Given input values of the independent variable T and the solution
//C     components Y(*), for each L = 1,2,...,NEQ evaluate the differential
//C     equation for the derivative of the Ith solution component and place the
//C     value in YP(L).  Do not alter the input values of T and Y(*).
//C  RETURN
//C  END
//C
//C  INPUT VARIABLE
//C
//C     TWANT     - double PRECISION
//C                 The next value of the independent variable where a
//C                 solution is desired.
//C
//C                 Constraints: TWANT must lie between the previous value
//C                 of TGOT (TSTART on the first call) and TEND. TWANT can be
//C                 equal to TEND, but it must be clearly distinguishable from
//C                 the previous value of TGOT (TSTART on the first call) in
//C                 the precision available.
//C
//C  OUTPUT VARIABLES
//C
//C     TGOT      - double PRECISION
//C                 A solution has been computed at this value of the
//C                 independent variable.
//C     YGOT(*)   - double PRECISION array of length NEQ
//C                 Approximation to the true solution at TGOT. Do not alter
//C                 the contents of this array
//C     YPGOT(*)  - double PRECISION array of length NEQ
//C                 Approximation to the first derivative of the true
//C                 solution at TGOT.
//C     YMAX(*)   - double PRECISION array of length NEQ
//C                 YMAX(L) is the largest magnitude of YGOT(L) computed at any
//C                 time in the integration from TSTART to TGOT. Do not alter
//C                 the contents of this array.
//C
//C  WORKSPACE
//C
//C     WORK(*)   - double PRECISION array as used in SETUP
//C                 Do not alter the contents of this array.
//C
//C  OUTPUT VARIABLE
//C
//C     UFLAG     - INTEGER
//C
//C                       SUCCESS.  TGOT = TWANT.
//C                 = 1 - Complete success.
//C
//C                       "SOFT" FAILURES
//C                 = 2 - Warning:  You are using METHOD = 3 inefficiently
//C                       by computing answers at many values of TWANT.  If
//C                       you really need answers at so many specific points,
//C                       it would be more efficient to compute them with
//C                       METHOD = 2.  To do this you would need to restart
//C                       from TGOT, YGOT(*) by a call to SETUP.  If you wish
//C                       to continue as you are, you may.
//C                 = 3 - Warning:  A considerable amount of work has been
//C                       expended.  If you wish to continue on to TWANT, just
//C                       call UT again.
//C                 = 4 - Warning:  It appears that this problem is "stiff".
//C                       You really should change to another code that is
//C                       intended for such problems, but if you insist, you can
//C                       continue with UT by calling it again.
//C
//C                       "HARD" FAILURES
//C                 = 5 - You are asking for too much accuracy. You cannot
//C                       continue integrating this problem.
//C                 = 6 - The global error assessment may not be reliable beyond
//C                       the current point in the integration.  You cannot
//C                       continue integrating this problem.
//C
//C                       "CATASTROPHIC" FAILURES
//C                 = 911 - The nature of the catastrophe is reported on
//C                         the standard output channel. Unless special
//C                         provision was made in advance (see rksuite.doc),
//C                         the computation then comes to a STOP.
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//      double PRECISION  TGOT, TWANT
//      INTEGER           UFLAG
//C     .. Array Arguments ..
//      double PRECISION  WORK(*), YGOT(*), YMAX(*), YPGOT(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block for General Workspace Pointers ..
//      INTEGER           PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      COMMON /RKCOM3/   PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      SAVE   /RKCOM3/
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Common Block for Integrator Options ..
//      LOGICAL           MSG, UTASK
//      COMMON /RKCOM8/   MSG, UTASK
//      SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      CHARACTER*6       SRNAME
//      PARAMETER         (SRNAME='UT')
   const char srname[] = "UT";
//      LOGICAL           ASK, TELL
//      PARAMETER         (ASK=.TRUE.,TELL=.FALSE.)
   const bool ask = true, tell = false;
//      INTEGER           MINUS1, MINUS2
//      PARAMETER         (MINUS1=-1,MINUS2=-2)
   const int minus1 = -1, minus2 = -2;
//      double PRECISION  ZERO
//      PARAMETER         (ZERO=0.0D+0)
   const double zero = 0.0;
//C     .. Local Scalars ..
   double hmin, tnow;
   int cflag, ier, l, nrec, state;
   bool baderr, goback;
//C     .. External Subroutines ..
//      EXTERNAL          CHKFL, CT, INTRP, RESET, RKMSG, RKSIT
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX, MIN
//C     .. Save statement ..
//      SAVE              UTEND, TLAST
   static double utend, tlast;
//C     .. Executable Statements ..
   ier = 1;
   nrec = 0;
   goback = false;
   baderr = false;
//C
//C  Is it permissible to call UT?
//C
   rksit(ask,"SETUP",state);
   if (state == 911) {
      ier = 912;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** A catastrophic error has already been detected elsewhere.");
      goto label100;
   }
   if (state == minus1) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have not called SETUP, so you cannot use UT.");
      goto label100;
   }
   if (!rkcom8.utask) {
      ier = 911;
      nrec = 2;
      sprintf(&rkcom9.rec[0][0]," ** You have called UT after you specified in SETUP that");
      sprintf(&rkcom9.rec[1][0]," ** you were going to use CT. This is not permitted.");
      goto label100;
   }
   rksit(ask,srname,state);
   if (state == 5 || state == 6) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** This routine has already returned with a hard failure.");
      sprintf(&rkcom9.rec[1][0]," ** You must call SETUP to start another problem.");
      goto label100;
   }
   state = minus2;
   rksit(tell,srname,state);
//C
   if (rkcom2.first) {
//C
//C  First call.
//C
//C  A value of TND is specified in SETUP. When INTP = .FALSE., as with
//C  METHD = 3, output is obtained at the specified TWANT by resetting TND
//C  to TWANT.  At this point, before the integration gets started, this can
//C  be done with a simple assignment.  Later it is done with a call to RESET.
//C  The original TND is SAVEd as a local variable UTEND.
//C
      utend = rkcom1.tnd;
      if (!rkcom4.intp) rkcom1.tnd = twant;
//C
//C  The last TGOT returned is SAVEd in the variable TLAST.  T (a variable
//C  passed through the common block RKCOM2) records how far the integration
//C  has advanced towards the specified TND.  When output is obtained by
//C  interpolation, the integration goes past the TGOT returned (T is closer
//C  to the specified TND than TGOT).  Initialize these variables and YMAX(*).
      tlast = rkcom1.tstrt;
      tgot = rkcom1.tstrt;
      for (l = 0; l < rkcom1.neqn; l++) {
         ymax[l] = fabs(work[rkcom3.pry+l]);
      }
//C
//C  If the code is to find an on-scale initial step size H, a bound was placed
//C  on H in SETUP.  Here the first output point is used to refine this bound.
      if (rkcom1.hstrt == zero) {
         rkcom2.h = min(fabs(rkcom2.h),fabs(twant-rkcom1.tstrt));
         hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(rkcom1.tstrt),fabs(rkcom1.tnd)));
         rkcom2.h = max(rkcom2.h,hmin);
      }
//C
   }
   else {
//C
//C  Subsequent call.
//C
      if (tlast == utend) {
         ier = 911;
         nrec = 3;
         sprintf(&rkcom9.rec[0][0]," ** You have called UT after reaching TEND. (Your last");
         sprintf(&rkcom9.rec[1][0]," ** call to UT resulted in TGOT = TEND.)  To start a new");
         sprintf(&rkcom9.rec[2][0]," ** problem, you will need to call SETUP.");
         goto label100;
      }
//C
   }
//C
//C  Check for valid TWANT.
//C
   if (rkcom1.dir*(twant-tlast) <= zero) {
      ier = 911;
      nrec = 4;
      sprintf(&rkcom9.rec[0][0]," ** You have made a call to UT with a TWANT that does");
      sprintf(&rkcom9.rec[1][0]," ** not lie between the previous value of TGOT (TSTART");
      sprintf(&rkcom9.rec[2][0]," ** on the first call) and TEND. This is not permitted.");
      sprintf(&rkcom9.rec[3][0]," ** Check your program carefully.");
      goto label100;
   }
   if (rkcom1.dir*(twant-utend) > zero) {
      hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(twant),fabs(utend)));
      if (fabs(twant-utend) < hmin) {
         ier = 911;
         nrec = 5;
         sprintf(&rkcom9.rec[0][0]," ** You have made a call to UT with a TWANT that does");
         sprintf(&rkcom9.rec[1][0]," ** not lie between the previous value of TGOT (TSTART on");
         sprintf(&rkcom9.rec[2][0]," ** the first call) and TEND. This is not permitted. TWANT");
         sprintf(&rkcom9.rec[3][0]," ** is very close to TEND, so you may have meant to set");
         sprintf(&rkcom9.rec[4][0]," ** it to be TEND exactly.  Check your program carefully.");
      }
      else {
         ier = 911;
         nrec = 4;
         sprintf(&rkcom9.rec[0][0]," ** You have made a call to UT with a TWANT that does");
         sprintf(&rkcom9.rec[1][0]," ** not lie between the previous value of TGOT (TSTART");
         sprintf(&rkcom9.rec[2][0]," ** on the first call) and TEND. This is not permitted.");
         sprintf(&rkcom9.rec[3][0]," ** Check your program carefully.");
      }
      goto label100;
   }
   if (!rkcom4.intp) {
      hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(tlast),fabs(twant)));
      if (fabs(twant-tlast) < hmin) {
         ier = 911;
         nrec = 4;
         sprintf(&rkcom9.rec[0][0]," ** You have made a call to UT with a TWANT that is not");
         sprintf(&rkcom9.rec[1][0]," ** sufficiently different from the last value of TGOT");
         sprintf(&rkcom9.rec[2][0]," ** (TSTART on the first call).  When using METHOD = 3,");
         sprintf(&rkcom9.rec[3][0]," ** it must differ by at least %13.5lf.",hmin);
         goto label100;
      }
//C
//C  We have a valid TWANT. There is no interpolation with this METHD and
//C  therefore we step to TWANT exactly by resetting TND with a call to RESET.
//C  On the first step this matter is handled differently as explained above.
//C
      if (!rkcom2.first) {
         reset(twant);
         chkfl(ask,baderr);
         if (baderr) goto label100;
      }
   }
//C
//C  Process output, decide whether to take another step.
//C
   label40:
//C
   if (rkcom4.intp) {
//C
//C  Interpolation is possible with this METHD.  The integration has
//C  already reached T. If this is past TWANT, GOBACK is set .TRUE. and
//C  the answers are obtained by interpolation.
//C
      goback = rkcom1.dir*(rkcom2.t-twant) >= zero;
      if (goback) {
         intrp(twant,"Both solution and derivative",rkcom1.neqn,ygot,
            ypgot,f,work,&work[rkcom3.printp],rkcom3.lnintp);
         chkfl(ask,baderr);
         if (baderr) goto label100;
         tgot = twant;
      }
   }
   else {
//C
//C  Interpolation is not possible with this METHD, so output is obtained
//C  by integrating to TWANT = TND.  Both YGOT(*) and YPGOT(*) are then
//C  already loaded with the solution at TWANT by CT.
//C
      goback = rkcom2.t == twant;
      if (goback) tgot = twant;
   }
//C
//C  Updating of YMAX(*) is done here to account for the fact that when
//C  interpolation is done, the integration goes past TGOT.  Note that YGOT(*)
//C  is not defined until CT is called.  YMAX(*) was initialized at TSTRT
//C  from values stored in WORK(*), so only needs to be updated for T
//C  different from TSTRT.
   if (rkcom2.t != rkcom1.tstrt) {
      for (l = 0; l < rkcom1.neqn; l++) {
         ymax[l] = max(ymax[l],fabs(ygot[l]));
      }
   }
//C
//C  If done, go to the exit point.
   if (goback) goto label100;
//C
//C  Take a step with CT in the direction of TND.  On exit, the solution is
//C  advanced to TNOW.  The way CT is written, the approximate solution at
//C  TNOW is available in both YGOT(*) and in WORK(*).  If output is obtained by
//C  stepping to the end (TNOW = TWANT = TND), YGOT(*) can be returned directly.
//C  If output is obtained by interpolation, the subroutine INTRP that does this
//C  uses the values in WORK(*) for its computations and places the approximate
//C  solution at TWANT in the array YGOT(*) for return to the calling program.
//C  The approximate derivative is handled in the same way. TNOW is output from
//C  CT and is actually a copy of T declared above in a common block.
//C
   ct(f,tnow,ygot,ypgot,work,cflag);
   ier = cflag;
//C
//C  A successful step by CT is indicated by CFLAG = 1 or = 2.
   if (cflag == 1) {
      goto label40;
   }
   else if (cflag == 2) {
//C
//C  Supplement the warning message written in CT.
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** The last message was produced on a call to CT from UT.");
      sprintf(&rkcom9.rec[1][0]," ** In UT the appropriate action is to change to METHOD = 2,");
      sprintf(&rkcom9.rec[2][0]," ** or, if insufficient memory is available, to METHOD = 1.");
   }
   else if (cflag <= 6) {
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** The last message was produced on a call to CT from UT.");
   }
   else {
      baderr = true;
   }
   tgot = rkcom2.t;
//C
//C  Update YMAX(*) before the return.
   for (l = 0; l < rkcom1.neqn; l++) {
      ymax[l] = max(ymax[l],fabs(ygot[l]));
   }
//C
//C  Exit point for UT.
//C
   label100:
//C
   if (baderr) {
      ier = 911;
      nrec = 4;
      sprintf(&rkcom9.rec[0][0]," ** An internal call by UT to a subroutine resulted in an");
      sprintf(&rkcom9.rec[1][0]," ** error that should not happen.  Check your program");
      sprintf(&rkcom9.rec[2][0]," ** carefully for array sizes, correct number of arguments,");
      sprintf(&rkcom9.rec[3][0]," ** type mismatches ... .");
   }
//C
   tlast = tgot;
//C
//C  All exits are done here after a call to RKMSG to report
//C  what happened and set UFLAG.
//C
   rkmsg(ier,srname,nrec,uflag);
//C
   return;
}

void RKSUITE::stat(int& totfcn, int& stpcst, double& waste, int& stpsok, double& hnext)
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code STAT and how it is used in
//C  conjunction with the integrators CT and UT, you should study the
//C  document file rksuite.doc carefully before attempting to use the code.
//C  The following "Brief Reminder" is intended only to remind you of the
//C  meaning, type, and size requirements of the arguments.
//C
//C  STAT is called to obtain some details about the integration.
//C
//C  OUTPUT VARIABLES
//C
//C     TOTFCN    - INTEGER
//C                 Total number of calls to F in the integration so far --
//C                 a measure of the cost of the integration.
//C     STPCST    - INTEGER
//C                 Cost of a typical step with this METHOD measured in
//C                 calls to F.
//C     WASTE     - double PRECISION
//C                 The number of attempted steps that failed to meet the
//C                 local error requirement divided by the total number of
//C                 steps attempted so far in the integration.
//C     STPSOK    - INTEGER
//C                 The number of accepted steps.
//C     HNEXT     - double PRECISION
//C                 The step size the integrator plans to use for the next step.
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//      double PRECISION  HNEXT, WASTE
//      INTEGER           STPCST, STPSOK, TOTFCN
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Global Error Assessment ..
//      double PRECISION  MAXERR, LOCMAX
//      INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
//     &                  PRZYNU
//      LOGICAL           ERASON, ERASFL
//      COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
//     &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
//      SAVE   /RKCOM6/
//C     .. Common Block for Integrator Options ..
//      LOGICAL           MSG, UTASK
//      COMMON /RKCOM8/   MSG, UTASK
//      SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      CHARACTER*6       SRNAME
//      PARAMETER         (SRNAME='STAT')
   const char srname[] = "STAT";
//      LOGICAL           ASK
//      INTEGER           MINUS1, MINUS2
//      PARAMETER         (ASK=.TRUE.,MINUS1=-1,MINUS2=-2)
   const bool ask = true;
   const int minus1 = -1, minus2 = -2;
//      double PRECISION  ZERO
//      PARAMETER         (ZERO=0.0D+0)
   const double zero = 0.0;
//C     .. Local Scalars ..
   int flag, ier, nrec, state;
//C     .. External Subroutines ..
//      EXTERNAL          RKMSG, RKSIT
//C     .. Intrinsic Functions ..
//      INTRINSIC         DBLE
//C     .. Executable Statements ..
//C
   ier = 1;
   nrec = 0;
//C
//C  Is it permissible to call STAT?
//C
   rksit(ask,srname,state);
   if (state == 911) {
      ier = 912;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** A catastrophic error has already been detected elsewhere.");
      goto label20;
   }
   if (state == minus2) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** You have already made a call to STAT after a hard");
      sprintf(&rkcom9.rec[1][0]," ** failure was reported from the integrator. You cannot");
      sprintf(&rkcom9.rec[2][0]," ** call STAT again.");
      goto label20;
   }
   rksit(ask,"CT",state);
   if (state == minus1) {
      ier = 911;
      nrec = 1;
      if (rkcom8.utask) {
         sprintf(&rkcom9.rec[0][0]," ** You have not called UT, so you cannot use STAT.");
      }
      else {
         sprintf(&rkcom9.rec[0][0]," ** You have not called CT, so you cannot use STAT.");
      }
      goto label20;
   }
//C
//C  Set flag so that the routine can only be called once after a hard 
//C  failure from the integrator.
   if (state == 5 || state == 6) ier = minus2;
//C
   totfcn = rkcom2.svnfcn + rkcom2.nfcn;
   if (rkcom6.erason) totfcn = totfcn + rkcom6.gnfcn;
   stpcst = int(rkcom5.cost);
   stpsok = rkcom2.okstp;
   if (rkcom2.okstp <= 1) {
      waste = zero;
   }
   else {
      waste = double(rkcom2.flstp)/double(rkcom2.flstp+rkcom2.okstp);
   }
   hnext = rkcom2.h;
//C
   label20:
//C
   rkmsg(ier,srname,nrec,flag);
//C
   return;
}

void RKSUITE::glberr(double rmserr[], double& errmax, double& terrmx, double work[])
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code GLBERR and how it is used in
//C  conjunction with UT and CT to solve initial value problems, you should
//C  study the document file rksuite.doc carefully before attempting to use 
//C  the code.  The following "Brief Reminder" is intended only to remind you 
//C  of the meaning, type, and size requirements of the arguments.
//C
//C  If ERRASS was set .TRUE. in the call to SETUP, then after any call to UT
//C  or CT to advance the integration to TNOW or TWANT, the subroutine GLBERR
//C  may be called to obtain an assessment of the true error of the integration.
//C  At each step and for each solution component Y(L), a more accurate "true"
//C  solution YT(L), an average magnitude "size(L)" of its size, and its error
//C                abs(Y(L) - YT(L))/max("size(L)",THRES(L))
//C  are computed.  The assessment returned in RMSERR(L) is the RMS (root-mean-
//C  square) average of the error in the Lth solution component over all steps
//C  of the integration from TSTART through TNOW.
//C
//C  OUTPUT VARIABLES
//C
//C     RMSERR(*) - double PRECISION array of length NEQ
//C                 RMSERR(L) approximates the RMS average of the true error 
//C                 of the numerical solution for the Ith solution component,
//C                 L = 1,2,...,NEQ.  The average is taken over all steps from
//C                 TSTART to TNOW.
//C     ERRMAX    - double PRECISION
//C                 The maximum (approximate) true error taken over all
//C                 solution components and all steps from TSTART to TNOW.
//C     TERRMX    - double PRECISION
//C                 First value of the independent variable where the
//C                 (approximate) true error attains the maximum value ERRMAX.
//C
//C  WORKSPACE
//C
//C     WORK(*)   - double PRECISION array as used in SETUP and UT or CT
//C                 Do not alter the contents of this array.
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//    double PRECISION  ERRMAX, TERRMX
//C     .. Array Arguments ..
//    double PRECISION  RMSERR(*), WORK(*)
//C     .. Common Block for Problem Definition ..
//    double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//    INTEGER           NEQN
//    COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//    SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//    double PRECISION  T, H, TOLD, HOLD
//    INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//    LOGICAL           FIRST, LAST
//    COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//   &                  FIRST, LAST
//    SAVE   /RKCOM2/
//C     .. Common Block for Global Error Assessment ..
//    double PRECISION  MAXERR, LOCMAX
//    INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
//   &                  PRZYNU
//    LOGICAL           ERASON, ERASFL
//    COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
//   &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
//    SAVE   /RKCOM6/
//C     .. Common Block for Integrator Options ..
//    LOGICAL           MSG, UTASK
//    COMMON /RKCOM8/   MSG, UTASK
//    SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//    CHARACTER*80      REC(10)
//    COMMON /RKCOM9/   REC
//    SAVE   /RKCOM9/
//C     .. Parameters ..
//    CHARACTER*6       SRNAME
//    PARAMETER         (SRNAME='GLBERR')
   const char srname[] = "GLBERR";
//    LOGICAL           ASK
//    PARAMETER         (ASK=.TRUE.)
   const bool ask = true;
//    INTEGER           MINUS1, MINUS2
//    PARAMETER         (MINUS1=-1,MINUS2=-2)
   const int minus1 = -1, minus2 = -2;
//C     .. Local Scalars ..
   int flag, ier, l, nrec, state;
//C     .. External Subroutines ..
//      EXTERNAL          RKMSG, RKSIT
//C     .. Intrinsic Functions ..
//    INTRINSIC         SQRT
//C     .. Executable Statements ..
//C
   ier = 1;
   nrec = 0;
//C
//C  Is it permissible to call GLBERR?
//C
   rksit(ask,srname,state);
   if (state == 911) {
      ier = 912;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** A catastrophic error has already been detected elsewhere.");
      goto label40;
   }
   if (state == minus2) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** You have already made a call to GLBERR after a hard");
      sprintf(&rkcom9.rec[1][0]," ** failure was reported from the integrator. You cannot");
      sprintf(&rkcom9.rec[2][0]," ** call GLBERR again.");
      goto label40;
   }
   rksit(ask,"CT",state);
   if (state == minus1) {
      ier = 911;
      nrec = 1;
      if (rkcom8.utask) {
         sprintf(&rkcom9.rec[0][0]," ** You have not yet called UT, so you cannot call GLBERR.");
      }
      else {
         sprintf(&rkcom9.rec[0][0]," ** You have not yet called CT, so you cannot call GLBERR.");
      }
      goto label40;
   }
//C
//C  Set flag so that the routine can only be called once after a hard 
//C  failure from the integrator.
   if (state == 5 || state == 6) ier = minus2;
//C
//C  Check that ERRASS was set properly for error assessment in SETUP.
//C
   if (!rkcom6.erason) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** No error assessment is available since you did not");
      sprintf(&rkcom9.rec[1][0]," ** ask for it in your call to the routine SETUP.");
      sprintf(&rkcom9.rec[2][0]," ** Check your program carefully.");
      goto label40;
   }
//C
//C  Check to see if the integrator has not actually taken a step.
//C
   if (rkcom2.okstp == 0) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** The integrator has not actually taken any successful");
      sprintf(&rkcom9.rec[1][0]," ** steps.  You cannot call GLBERR in this circumstance.");
      sprintf(&rkcom9.rec[2][0]," ** Check your program carefully.");
      goto label40;
   }
//C
//C  Compute RMS error and set output variables.
//C
   errmax = rkcom6.maxerr;
   terrmx = rkcom6.locmax;
   for (l = 0; l < rkcom1.neqn; l++) {
      rmserr[l] = sqrt(work[rkcom6.przerr+l]/rkcom2.okstp);
   }
//C
   label40:
//C
   rkmsg(ier,srname,nrec,flag);
//C
   return;
}

void RKSUITE::ct(void (*f)(double, double*, double*), double& tnow, double ynow[],
   double ypnow[], double work[],int& cflag)
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code CT and how it is used in
//C  conjunction with SETUP to solve initial value problems, you should study
//C  the document file rksuite.doc carefully before attempting to use the code.
//C  The following "Brief Reminder" is intended only to remind you of the
//C  meaning, type, and size requirements of the arguments.
//C
//C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE CALLING PROGRAM:
//C
//C     F         - name of the subroutine for evaluating the differential
//C                 equations.
//C
//C  The subroutine F must have the form
//C
//C  SUBROUTINE F(T,Y,YP)
//C  double PRECISION T,Y(*),YP(*)
//C     Using the input values of the independent variable T and the solution
//C     components Y(*), for each L = 1,2,...,NEQ evaluate the differential
//C     equation for the derivative of the Lth solution component and place the
//C     value in YP(L).  Do not alter the input values of T and Y(*).
//C  RETURN
//C  END
//C
//C  OUTPUT VARIABLES
//C
//C     TNOW      - double PRECISION
//C                 Current value of the independent variable.
//C     YNOW(*)   - double PRECISION array of length NEQ
//C                 Approximation to the true solution at TNOW.
//C     YPNOW(*)  - double PRECISION array of length NEQ
//C                 Approximation to the first derivative of the
//C                 true solution at TNOW.
//C
//C  WORKSPACE
//C
//C     WORK(*)   - double PRECISION array as used in SETUP
//C                 Do not alter the contents of this array.
//C
//C  OUTPUT VARIABLE
//C
//C     CFLAG     - INTEGER
//C
//C                       SUCCESS.  A STEP WAS TAKEN TO TNOW.
//C                 = 1 - Complete success.
//C
//C                       "SOFT" FAILURES
//C                 = 2 - Warning:  You have obtained an answer by integrating
//C                       to TEND (TNOW = TEND).  You have done this at least
//C                       100 times, and monitoring of the computation reveals
//C                       that this way of getting output has degraded the
//C                       efficiency of the code. If you really need answers at
//C                       so many specific points, it would be more efficient to
//C                       get them with INTRP.  (If METHOD = 3, you would need
//C                       to change METHOD and restart from TNOW, YNOW(*) by a
//C                       call to SETUP.)  If you wish to continue as you are,
//C                       you may.
//C                 = 3 - Warning:  A considerable amount of work has been
//C                       expended. To continue the integration, just call
//C                       CT again.
//C                 = 4 - Warning:  It appears that this problem is "stiff".
//C                       You really should change to another code that is
//C                       intended for such problems, but if you insist, you 
//C                       can continue with CT by calling it again.
//C
//C                       "HARD" FAILURES
//C                 = 5 - You are asking for too much accuracy. You cannot
//C                       continue integrating this problem.
//C                 = 6 - The global error assessment may not be reliable beyond
//C                       the current point in the integration.  You cannot
//C                       continue integrating this problem.
//C
//C                       "CATASTROPHIC" FAILURES
//C                 = 911 - The nature of the catastrophe is reported on
//C                         the standard output channel. Unless special
//C                         provision was made in advance (see rksuite.doc),
//C                         the computation then comes to a STOP.
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//    double PRECISION  TNOW
//    INTEGER           CFLAG
//C     .. Array Arguments ..
//    double PRECISION  WORK(*), YNOW(*), YPNOW(*)
//C     .. Subroutine Arguments ..
//    EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//    double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//    INTEGER           NEQN
//    COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//    SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//    double PRECISION  T, H, TOLD, HOLD
//    INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//    LOGICAL           FIRST, LAST
//    COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//   &                  FIRST, LAST
//    SAVE   /RKCOM2/
//C     .. Common Block for General Workspace Pointers ..
//    INTEGER           PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//   &                  PRSTGS, PRINTP, LNINTP
//    COMMON /RKCOM3/   PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//   &                  PRSTGS, PRINTP, LNINTP
//    SAVE   /RKCOM3/
//C     .. Common Block to hold Formula Characterisitcs ..
//    double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//   &                  RS, RS1, RS2, RS3, RS4
//    INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//    LOGICAL           FSAL
//    COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//   &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//   &                  NSEC, FSAL
//    SAVE   /RKCOM5/
//C     .. Common Block for Global Error Assessment ..
//    double PRECISION  MAXERR, LOCMAX
//    INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
//   &                  PRZYNU
//    LOGICAL           ERASON, ERASFL
//    COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
//   &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
//    SAVE   /RKCOM6/
//C     .. Common Block for Environment Parameters ..
//    double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//    INTEGER           OUTCH
//    COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//   &                  OUTCH
//    SAVE   /RKCOM7/
//C     .. Common Block for Integrator Options ..
//    LOGICAL           MSG, UTASK
//    COMMON /RKCOM8/   MSG, UTASK
//    SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//    CHARACTER*80      REC(10)
//    COMMON /RKCOM9/   REC
//    SAVE   /RKCOM9/
//C     .. Parameters ..
//    CHARACTER*6       SRNAME
//    PARAMETER         (SRNAME='CT')
   const char srname[] = "CT";
//    LOGICAL           ASK, TELL
//    PARAMETER         (ASK=.TRUE.,TELL=.FALSE.)
   const bool ask = true, tell = false;
//    INTEGER           MINUS1, MINUS2
//    PARAMETER         (MINUS1=-1,MINUS2=-2)
   const int minus1 = -1, minus2 = -2;
//    INTEGER           MAXFCN
//    PARAMETER         (MAXFCN=5000)
   const int maxfcn = 5000;
//    double PRECISION  ZERO, PT1, PT9, ONE, TWO, HUNDRD
//    PARAMETER         (ZERO=0.0D+0,PT1=0.1D+0,PT9=0.9D+0,ONE=1.0D+0,
//   &                  TWO=2.0D+0,HUNDRD=100.0D+0)
   const double zero = 0.0, pt1 = 0.1, pt9 = 0.9, one = 1.0, two = 2.0,
      hundrd = 100.0;
//C     .. Local Scalars ..
   double alpha, beta, err, hmin, htry, tau,
      temp1, temp2, ypnorm;
   int ier, l, nrec, point, state;
   bool failed, main, phase1, phase3, toomch;
//C     .. External Subroutines ..
//      EXTERNAL          RKMSG, RKSIT, STEP, STIFF, TRUERR
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX, MIN, SIGN
//C     .. Save statement ..
   static int jflstp, ntend, ynew, ypold;
   static double errold, havg;
   static bool phase2, chkeff;
//      SAVE              JFLSTP, NTEND, ERROLD, HAVG, PHASE2, YNEW,
//     &                  YPOLD, CHKEFF
//C     .. Executable Statements ..
//C
   ier = 1;
   nrec = 0;
//C
//C  Is it permissible to call CT?
//C
   rksit(ask,"SETUP",state);
   if (state == 911) {
      ier = 912;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** A catastrophic error has already been detected elsewhere.");
      goto label100;
   }
   if (state == minus1) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have not called SETUP, so you cannot use CT.");
      goto label100;
   }
   if (rkcom8.utask) {
      rksit(ask,"UT",state);
      if (state != minus2) {
         ier = 911;
         nrec = 2;
         sprintf(&rkcom9.rec[0][0]," ** You have called CT after you specified in SETUP that");
         sprintf(&rkcom9.rec[1][0]," ** you were going to use UT. This is not permitted.");
         rkcom8.utask = false;
         goto label180;
      }
   }
   rksit(ask,srname,state);
   if (state == 5 || state == 6) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** CT has already returned with a flag value of 5 or 6.");
      sprintf(&rkcom9.rec[1][0]," ** You cannot continue integrating this problem. You must");
      sprintf(&rkcom9.rec[2][0]," ** call SETUP to start another problem.");
      goto label180;
   }
//C
   if (rkcom2.first) {
//C
//C  First call in an integration -- initialize everything.
//C
      chkeff = false;
      ntend = 0;
      jflstp = 0;
//C
//C  A scratch area of WORK(*) starting at PRSCR is used to hold two
//C  arrays in this subroutine: the higher order approximate solution at
//C  the end of a step and the approximate derivative of the solution at
//C  the end of the last step. To make this clear, local pointers YNEW and
//C  YPOLD are used.
      ynew = rkcom3.prscr;
      ypold = rkcom3.prscr;
//C
//C  For this first step T was initialized to TSTRT in SETUP and the
//C  starting values YSTART(*) were loaded into the area of WORK(*) reserved
//C  for the current solution approximation starting at location PRY. The
//C  derivative is now computed and stored in WORK(*) starting at PRYP.
//C  Subsequently these arrays are copied to the output vectors YNOW(*)
//C  and YPNOW(*).
      f(rkcom2.t,&work[rkcom3.pry],&work[rkcom3.pryp]);
      rkcom2.nfcn++;
      for (l = 0; l < rkcom1.neqn; l++) {
         ynow[l] = work[rkcom3.pry+l];
         ypnow[l] = work[rkcom3.pryp+l];
      }
//C
//C  Set dependent variables for error assessment.
      if (rkcom6.erason) {
         for (l = 0; l < rkcom1.neqn; l++) {
            work[rkcom6.przy+l] = ynow[l];
            work[rkcom6.przyp+l] = ypnow[l];
         }
      }
//C
//C  The weights for the control of the error depend on the size of the
//C  solution at the beginning and at the end of the step. On the first
//C  step we do not have all this information. Whilst determining the
//C  initial step size we initialize the weight vector to the larger of
//C  abs(Y(L)) and the threshold for this component.
      for (l = 0; l < rkcom1.neqn; l++) {
         work[rkcom3.prwt+l] = max(fabs(ynow[l]),work[rkcom3.prthrs+l]);
      }
//C
//C  If HSTRT is equal to zero, the code is to find an on-scale initial step
//C  size H.  CT has an elaborate scheme of three phases for finding such an H,
//C  and some preparations are made earlier.  In SETUP an upper bound is placed
//C  on H that reflects the scale of the independent variable. When UTASK is
//C  .TRUE., UT refines this bound using the first output point.  Here in CT
//C  PHASE1 applies a rule of thumb based on the error control, the order of the
//C  the formula, and the size of the initial slope to get a crude approximation
//C  to an on-scale H.  PHASE2 may reduce H in the course of taking the first
//C  step.  PHASE3 repeatedly adjusts H and retakes the first step until H is
//C  on scale.
//C
//C  A guess for the magnitude of the first step size H can be provided to SETUP
//C  as HSTART.  If it is too big or too small, it is ignored and the automatic
//C  determination of an on-scale initial step size is activated.  If it is
//C  acceptable, H is set to HSTART in SETUP.  Even when H is supplied to CT,
//C  PHASE3 of the scheme for finding an on-scale initial step size is made
//C  active so that the code can deal with a bad guess.
//C
      phase1 = rkcom1.hstrt == zero;
      phase2 = phase1;
      phase3 = true;
      if (phase1) {
         rkcom2.h = fabs(rkcom2.h);
         ypnorm = zero;
         for (l = 0; l < rkcom1.neqn; l++) {
            if (fabs(ynow[l]) != zero) {
               ypnorm = max(ypnorm,fabs(ypnow[l])/work[rkcom3.prwt+l]);
            }
         }
         tau = pow(rkcom1.tolr, rkcom5.expon);
         if (rkcom2.h*ypnorm > tau) rkcom2.h = tau/ypnorm;
         hmin = max(rkcom7.tiny,
            rkcom5.toosml*max(fabs(rkcom1.tstrt),fabs(rkcom1.tnd)));
         rkcom2.h = rkcom1.dir*max(rkcom2.h,hmin);
         phase1 = false;
      }
   }
//C
   else {
//C
//C Continuation call
//C
      if (rkcom2.last) {
         ier = 911;
         nrec = 3;
         sprintf(&rkcom9.rec[0][0]," ** You have already reached TEND ( = %13.5lf).",rkcom1.tnd);
         sprintf(&rkcom9.rec[1][0]," ** To integrate further with the same problem you must");
         sprintf(&rkcom9.rec[2][0]," ** call the routine RESET with a new value of TEND.");
         goto label180;
      }
   }
//C
//C  Begin computation of a step here.
//C
   failed = false;
//C
   label100:
   rkcom2.h = sign(fabs(rkcom2.h),rkcom1.dir);
//C
//C  Reduce the step size if necessary so that the code will not step
//C  past TND.  "Look ahead" to prevent unnecessarily small step sizes.
   rkcom2.last = rkcom1.dir*((rkcom2.t+rkcom2.h)-rkcom1.tnd) >= zero;
   if (rkcom2.last) {
      rkcom2.h = rkcom1.tnd - rkcom2.t;
   }
   else if (rkcom1.dir*((rkcom2.t+two*rkcom2.h)-rkcom1.tnd) >= zero) {
      rkcom2.h = (rkcom1.tnd-rkcom2.t)/two;
   }
//C
//C  When the integrator is at T and attempts a step of H, the function
//C  defining the differential equations will be evaluated at a number of
//C  arguments between T and T+H.  If H is too small, these arguments cannot
//C  be clearly distinguished in the precision available.
//C
   hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(rkcom2.t),fabs(rkcom2.t+rkcom2.h)));
   if (fabs(rkcom2.h) < hmin) {
      ier = 5;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** In order to satisfy your error requirements CT would");
      sprintf(&rkcom9.rec[1][0]," ** have to use a step size of %13.5lf at TNOW = %13.5lf",rkcom2.h,rkcom2.t);
      sprintf(&rkcom9.rec[2][0]," ** This is too small for the machine precision.");
      goto label180;
   }
//C
//C  Monitor the impact of output on the efficiency of the integration.
//C
   if (chkeff) {
      ntend++;
      if (ntend >= 100 && ntend >= rkcom2.okstp/3) {
         ier = 2;
         nrec = 6;
         sprintf(&rkcom9.rec[0][0]," ** More than 100 output points have been obtained by");
         sprintf(&rkcom9.rec[1][0]," ** integrating to TEND.  They have been sufficiently close");
         sprintf(&rkcom9.rec[2][0]," ** to one another that the efficiency of the integration has");
         sprintf(&rkcom9.rec[3][0]," ** been degraded. It would probably be (much) more efficient");
         sprintf(&rkcom9.rec[4][0]," ** to obtain output by interpolating with INTRP (after");
         sprintf(&rkcom9.rec[5][0]," ** changing to METHOD=2 if you are using METHOD = 3).");
         ntend = 0;
         goto label180;
      }
   }
//C
//C  Check for stiffness and for too much work.  Stiffness can be
//C  checked only after a successful step.
//C
   if (!failed) {
//C
//C  Check for too much work.
      toomch = rkcom2.nfcn > maxfcn;
      if (toomch) {
         ier = 3;
         nrec = 3;
         sprintf(&rkcom9.rec[0][0]," ** Approximately %6d function evaluations have been",maxfcn);
         sprintf(&rkcom9.rec[1][0]," ** used to compute the solution since the integration");
         sprintf(&rkcom9.rec[2][0]," ** started or since this message was last printed.");
//C
//C  After this warning message, NFCN is reset to permit the integration
//C  to continue.  The total number of function evaluations in the primary
//C  integration is SVNFCN + NFCN.
         rkcom2.svnfcn = rkcom2.svnfcn + rkcom2.nfcn;
         rkcom2.nfcn = 0;
      }
//C
//C  Check for stiffness.  NREC is passed on to STIFF because when
//C  TOOMCH = .TRUE. and stiffness is diagnosed, the message about too
//C  much work is augmented inside STIFF to explain that it is due to
//C  stiffness.
      stiff(f,havg,jflstp,toomch,maxfcn,work,ier,nrec);
//C
      if (ier != 1) goto label180;
   }
//C
//C  Take a step.  Whilst finding an on-scale H (PHASE2 = .TRUE.), the input
//C  value of H might be reduced (repeatedly), but it will not be reduced
//C  below HMIN.  The local error is estimated, a weight vector is formed,
//C  and a weighted maximum norm, ERR, of the local error is returned.
//C  The variable MAIN is input as .TRUE. to tell STEP that this is the
//C  primary, or "main", integration.
//C
//C  H resides in the common block /RKCOM2/ which is used by both CT and STEP;
//C  since it may be changed inside STEP, a local copy is made to ensure
//C  portability of the code.
//C
   main = true;
   htry = rkcom2.h;
   step(f,rkcom1.neqn,rkcom2.t,&work[rkcom3.pry],&work[rkcom3.pryp],&work[rkcom3.prstgs],rkcom1.tolr,htry,
      &work[rkcom3.prwt],&work[ynew],&work[rkcom3.prerst],err,main,hmin,
      &work[rkcom3.prthrs],phase2);
   rkcom2.h = htry;
//C
//C  Compare the norm of the local error to the tolerance.
//C
   if (err > rkcom1.tolr) {
//C
//C  Failed step.  Reduce the step size and try again.
//C
//C  First step:  Terminate PHASE3 of the search for an on-scale step size.
//C               The step size is not on scale, so ERR may not be accurate;
//C               reduce H by a fixed factor.  Failed attempts to take the
//C               first step are not counted.
//C  Later step:  Use ERR to compute an "optimal" reduction of H.  More than
//C               one failure indicates a difficulty with the problem and an
//C               ERR that may not be accurate, so reduce H by a fixed factor.
//C
      if (rkcom2.first) {
         phase3 = false;
         alpha = rkcom5.rs1;
      }
      else {
         rkcom2.flstp++;
         jflstp++;
         if (failed) {
            alpha = rkcom5.rs1;
         }
         else {
            alpha = pow(rkcom5.safety*(rkcom1.tolr/err),rkcom5.expon);
            alpha = max(alpha,rkcom5.rs1);
         }
      }
      rkcom2.h = alpha*rkcom2.h;
      failed = true;
      goto label100;
   }
//C
//C  Successful step.
//C
//C  Predict a step size appropriate for the next step.  After the first
//C  step the prediction can be refined using an idea of H.A. Watts that
//C  takes account of how well the prediction worked on the previous step.
   beta = pow((err/rkcom1.tolr),rkcom5.expon);
   if (!rkcom2.first) {
      temp1 = (pow(err,rkcom5.expon))/rkcom2.h;
      temp2 = (pow(errold,rkcom5.expon))/rkcom2.hold;
      if (temp1 < temp2*hundrd && temp2 < temp1*hundrd) {
         beta = beta*(temp1/temp2);
      }
   }
   alpha = rkcom5.rs3;
   if (rkcom5.safety < beta*alpha) alpha = rkcom5.safety/beta;
//C
//C  On the first step a search is made for an on-scale step size.  PHASE2
//C  of the scheme comes to an end here because a step size has been found
//C  that is both successful and has a credible local error estimate. Except
//C  in the special case that the first step is also the last, the step is
//C  repeated in PHASE3 as long as an increase greater than RS2 appears
//C  possible.  An increase as big as RS3 is permitted.  A step failure
//C  terminates PHASE3.
//C
   if (rkcom2.first) {
      phase2 = false;
      phase3 = phase3 && !rkcom2.last && (alpha > rkcom5.rs2);
      if (phase3) {
         rkcom2.h = alpha*rkcom2.h;
         goto label100;
      }
   }
//C
//C  After getting on scale, step size changes are more restricted.
   alpha = min(alpha,rkcom5.rs);
   if (failed) alpha = min(alpha,one);
   alpha = max(alpha,rkcom5.rs1);
   rkcom2.hold = rkcom2.h;
   rkcom2.h = alpha*rkcom2.h;
//C
//C  For the diagnosis of stiffness, an average accepted step size, HAVG,
//C  must be computed and SAVEd.
   if (rkcom2.first) {
      havg = rkcom2.hold;
   }
   else {
      havg = pt9*havg + pt1*rkcom2.hold;
   }
//C
   rkcom2.first = false;
   errold = err;
   rkcom2.told = rkcom2.t;
//C
//C  Take care that T is set to precisely TND when the end of the
//C  integration is reached.
   if (rkcom2.last) {
      rkcom2.t = rkcom1.tnd;
   }
   else {
      rkcom2.t += rkcom2.hold;
   }
//C
//C  Increment counter on accepted steps.  Note that successful steps
//C  that are repeated whilst getting on scale are not counted.
   rkcom2.okstp++;
//C
//C  Advance the current solution and its derivative.  (Stored in WORK(*)
//C  with the first location being PRY and PRYP, respectively.)  Update the
//C  previous solution and its derivative.  (Stored in WORK(*) with the first
//C  location being PRYOLD and YPOLD, respectively.)  Note that the previous
//C  derivative will overwrite YNEW(*).
//C
   for (l = 0; l < rkcom1.neqn; l++) {
      work[rkcom3.pryold+l] = work[rkcom3.pry+l];
      work[rkcom3.pry+l] = work[ynew+l];
      work[ypold+l] = work[rkcom3.pryp+l];
   }
//C
   if (rkcom5.fsal) {
//C
//C  When FSAL = .TRUE., YP(*) is the last stage of the step.
      point = rkcom3.prstgs + (rkcom5.lststg-1)*rkcom1.neqn;
      for (l=0; l < rkcom1.neqn; l++) {
         work[rkcom3.pryp+l] = work[point+l];
      }
   }
   else {
//C
//C  Call F to evaluate YP(*).
      f(rkcom2.t,&work[rkcom3.pry],&work[rkcom3.pryp]);
      rkcom2.nfcn++;
   }
//C
//C  If global error assessment is desired, advance the secondary
//C  integration from TOLD to T.
//C
   if (rkcom6.erason) {
      truerr(f,rkcom1.neqn,&work[rkcom3.pry],rkcom1.tolr,&work[rkcom3.prwt],&work[rkcom6.przy],
         &work[rkcom6.przyp],&work[rkcom6.przerr],&work[rkcom6.przynu],&work[rkcom6.przers],
         &work[rkcom6.przstg],ier);
      if (ier == 6) {
//C
//C  The global error estimating procedure has broken down. Treat it as a
//C  failed step. The solution and derivative are reset to their values at
//C  the beginning of the step since the last valid error assessment refers
//C  to them.
         rkcom2.okstp--;
         rkcom6.erasfl = true;
         rkcom2.last = false;
         rkcom2.t = rkcom2.told;
         rkcom2.h = rkcom2.hold;
         for (l = 0; l < rkcom1.neqn; l++) {
            work[rkcom3.pry+l] = work[rkcom3.pryold+l];
            work[rkcom3.pryp+l] = work[ypold+l];
         }
         if (rkcom2.okstp > 1) {
            nrec = 2;
            sprintf(&rkcom9.rec[0][0]," ** The global error assessment may not be reliable for T past");
            sprintf(&rkcom9.rec[1][0]," ** TNOW = %13.5lf.  The integration is being terminated.",rkcom2.t);
         }
         else {
            nrec = 2;
            sprintf(&rkcom9.rec[0][0]," ** The global error assessment algorithm failed at the start");
            sprintf(&rkcom9.rec[1][0]," ** the integration.  The integration is being terminated.");
         }
         goto label180;
      }
   }
//C
//C
//C  Exit point for CT
//C
  label180:
//C
//C  Set the output variables and flag that interpolation is permitted
//C
   if (ier < 911) {
      tnow = rkcom2.t;
      rkcom2.last = tnow == rkcom1.tnd;
      chkeff = rkcom2.last;
      for (l = 0; l < rkcom1.neqn; l++) {
         ynow[l] = work[rkcom3.pry+l];
         ypnow[l] = work[rkcom3.pryp+l];
      }
      if (ier == 1) {
         state = minus2;
         rksit(tell,"INTRP",state);
      }
   }
//C
//C  Call RKMSG to report what happened and set CFLAG.
//C
   rkmsg(ier,srname,nrec,cflag);
//C
   return;
}

void RKSUITE::intrp(double twant, char reqest[], int nwant, double ywant[], double ypwant[],
   void (*f)(double, double*, double*), double work[], double wrkint[],
   int lenint)
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code INTRP and how it is used in
//C  conjunction with CT to solve initial value problems, you should study the
//C  document file rksuite.doc carefully before attempting to use the code. The
//C  following "Brief Reminder" is intended only to remind you of the meaning,
//C  type, and size requirements of the arguments.
//C
//C  When integrating with METHOD = 1 or 2, answers may be obtained inexpensively
//C  between steps by interpolation. INTRP is called after a step by CT from a
//C  previous value of T, called TOLD below, to the current value of T to get
//C  an answer at TWANT. You can specify any value of TWANT you wish, but
//C  specifying a value outside the interval [TOLD,T] is likely to yield
//C  answers with unsatisfactory accuracy.
//C
//C  INPUT VARIABLE
//C
//C     TWANT     - double PRECISION
//C                 The value of the independent variable where a solution
//C                 is desired.
//C
//C  The interpolant is to be evaluated at TWANT to approximate the solution
//C  and/or its first derivative there.  There are three cases:
//C
//C  INPUT VARIABLE
//C
//C     REQEST    - CHARACTER*(*)
//C                 Only the first character of REQEST is significant.
//C                 REQEST(1:1) = 'S' or 's'- compute approximate 'S'olution
//C                                           only.
//C                             = 'D' or 'd'- compute approximate first
//C                                           'D'erivative of the solution only.
//C                             = 'B' or 'b'- compute 'B'oth approximate solution
//C                                           and first derivative.
//C                 Constraint: REQEST(1:1) must be 'S','s','D','d','B' or 'b'.
//C
//C  If you intend to interpolate at many points, you should arrange for the
//C  the interesting components to be the first ones because the code
//C  approximates only the first NWANT components.
//C
//C  INPUT VARIABLE
//C
//C     NWANT     - INTEGER
//C                 Only the first NWANT components of the answer are to be
//C                 computed.
//C                 Constraint:  NEQ >= NWANT >= 1
//C
//C  OUTPUT VARIABLES
//C
//C     YWANT(*)  - double PRECISION array of length NWANT
//C                 Approximation to the first NWANT components of the true
//C                 solution at TWANT when REQESTed.
//C     YPWANT(*) - double PRECISION array of length NWANT
//C                 Approximation to the first NWANT components of the first
//C                 derivative of the true solution at TWANT when REQESTed.
//C
//C  NAME DECLARED IN AN EXTERNAL STATEMENT IN THE PROGRAM CALLING INTRP:
//C
//C     F         - name of the subroutine for evaluating the differential
//C                 equations as provided to CT.
//C
//C  WORKSPACE
//C
//C     WORK(*)   - double PRECISION array as used in SETUP and CT
//C                 Do not alter the contents of this array.
//C
//C     WRKINT(*) - double PRECISION array of length LENINT
//C                 Do not alter the contents of this array.
//C
//C     LENINT    - INTEGER
//C                 Length of WRKINT. If
//C                 METHOD = 1, LENINT must be at least 1
//C                        = 2, LENINT must be at least NEQ+MAX(NEQ,5*NWANT)
//C                        = 3--CANNOT BE USED WITH THIS SUBROUTINE
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//      double PRECISION  TWANT
//      INTEGER           LENINT, NWANT
//      CHARACTER*(*)     REQEST
//C     .. Array Arguments ..
//      double PRECISION  WORK(*), WRKINT(*), YPWANT(*), YWANT(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block for General Workspace Pointers ..
//      INTEGER           PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      COMMON /RKCOM3/   PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      SAVE   /RKCOM3/
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Common Block for Integrator Options ..
//      LOGICAL           MSG, UTASK
//      COMMON /RKCOM8/   MSG, UTASK
//      SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      CHARACTER*6       SRNAME
//      PARAMETER         (SRNAME='INTRP')
   const char srname[] = "INTRP";
//      LOGICAL           ASK
//      INTEGER           PLUS1, MINUS1, MINUS2
//      PARAMETER         (ASK=.TRUE.,PLUS1=1,MINUS1=-1,MINUS2=-2)
   const bool ask = true;
   const int plus1 = 1, minus1 = -1, minus2 = -2;
//C     .. Local Scalars ..
   int flag, ichk, ier, nrec, state, state1;
   bool legalr;
   char reqst1;
//C     .. External Subroutines ..
//      EXTERNAL          EVALI, FORMI, RKMSG, RKSIT
//C     .. Intrinsic Functions ..
//      INTRINSIC         MAX
//C     .. Save statement ..
//      SAVE              NWNTSV, ININTP, STARTP
//C     .. Data statements ..
//      DATA              NWNTSV/MINUS1/
   static int nwntsv = minus1, startp;
   static bool inintp;
//C     .. Executable Statements ..
//C
   ier = 1;
   nrec = 0;
//C
//C  Is it permissible to call INTRP?
//C
   rksit(ask,"CT",state);
   if (state == 911) {
      ier = 912;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** A catastrophic error has already been detected elsewhere.");
      goto label20;
   }
   if (rkcom8.utask) {
      rksit(ask,"UT",state1);
      if (state1 != minus2) {
         ier = 911;
         nrec = 2;
         sprintf(&rkcom9.rec[0][0]," ** You have called INTRP after you specified to SETUP");
         sprintf(&rkcom9.rec[1][0]," ** that you were going to use UT. This is not permitted.");
         goto label20;
      }
   }
   if (state == minus1) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have not called CT, so you cannot use INTRP.");
      goto label20;
   }
   if (state > plus1) {
      ier = 911;
      nrec = 2;
      sprintf(&rkcom9.rec[0][0]," ** CT has returned with a flag value greater than 1.");
      sprintf(&rkcom9.rec[1][0]," ** You cannot call INTRP in this circumstance.");
      goto label20;
   }
//C
//C  Check input
//C
   reqst1 = reqest[0];
   legalr = reqst1 == 'S' || reqst1 == 's' ||
            reqst1 == 'D' || reqst1 == 'd' ||
            reqst1 == 'B' || reqst1 == 'b';
   if (!legalr) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** You have set the first character of");
      sprintf(&rkcom9.rec[1][0]," ** REQEST to be '%c'. It must be one of",reqst1);
      sprintf(&rkcom9.rec[2][0]," ** 'S','s','D','d','B' or 'b'.");
      goto label20;
   }
//C
   if (nwant > rkcom1.neqn) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** You have specified the value of NWANT to be %6d. This",nwant);
      sprintf(&rkcom9.rec[1][0]," ** is greater than %6d, which is the number of equations",rkcom1.neqn);
      sprintf(&rkcom9.rec[2][0]," ** in the system being integrated.");
      goto label20;
   }
   else if (nwant < 1) {
      ier = 911;
      nrec = 3;
      sprintf(&rkcom9.rec[0][0]," ** You have specified the value of NWANT to be %6d, but",nwant);
      sprintf(&rkcom9.rec[1][0]," ** this is less than 1. You cannot interpolate a zero or");
      sprintf(&rkcom9.rec[2][0]," ** negative number of components.");
      goto label20;
   }
//C
   if (rkcom4.methd == 1) {
      if (lenint < 1) {
         ier = 911;
         nrec = 2;
         sprintf(&rkcom9.rec[0][0]," ** You have specified LENINT to be %6d.",lenint);
         sprintf(&rkcom9.rec[1][0]," ** This is too small. LENINT must be at least 1.");
         goto label20;
      }
      startp = /*1*/0;
   }
   else if (rkcom4.methd == 2) {
      ichk = rkcom1.neqn + max(rkcom1.neqn,5*nwant);
      if (lenint < ichk) {
         ier = 911;
         nrec = 3;
         sprintf(&rkcom9.rec[0][0]," ** You have specified LENINT to be %6d. This is too",lenint);
         sprintf(&rkcom9.rec[1][0]," ** small. NINT must be at least NEQ + MAX(NEQ, 5*NWANT)");
         sprintf(&rkcom9.rec[2][0]," ** which is %6d.",ichk);
         goto label20;
      }
      startp = rkcom1.neqn + /*1*/0;
   }
   else if (rkcom4.methd == 3) {
      ier = 911;
      nrec = 5;
      sprintf(&rkcom9.rec[0][0]," ** You have been using CT with METHOD = 3 to integrate your");
      sprintf(&rkcom9.rec[1][0]," ** equations. You have just called INTRP, but interpolation");
      sprintf(&rkcom9.rec[2][0]," ** is not available for this METHOD. Either use METHOD = 2,");
      sprintf(&rkcom9.rec[3][0]," ** for which interpolation is available, or use RESET to make");
      sprintf(&rkcom9.rec[4][0]," ** CT step exactly to the points where you want output.");
      goto label20;
   }
//C
//C  Has the interpolant been initialised for this step?
//C
   rksit(ask,srname,state);
   inintp = (state != minus2);
//C
//C  Some initialization must be done before interpolation is possible.
//C  To reduce the overhead, the interpolating polynomial is formed for
//C  the first NWANT components.  In the unusual circumstance that NWANT
//C  is changed while still interpolating within the span of the current
//C  step, the scheme must be reinitialized to accomodate the additional
//C  components.
//C
   if (!inintp || nwant != nwntsv) {
   	bool notinintp = !inintp;
//C
//C  At present the derivative of the solution at the previous step, YPOLD(*),
//C  is stored in the scratch array area starting at PRSCR. In the case of
//C  METHD = 1 we can overwrite the stages.
//C
      if (rkcom4.methd == 1) {
         formi(f,rkcom1.neqn,nwant,&work[rkcom3.pry],&work[rkcom3.pryp],&work[rkcom3.pryold],
            &work[rkcom3.prscr],&work[rkcom3.prstgs],notinintp,
            &work[rkcom3.prstgs],&work[rkcom3.prstgs]);
      }
      else {
         formi(f,rkcom1.neqn,nwant,&work[rkcom3.pry],&work[rkcom3.pryp],&work[rkcom3.pryold],
            &work[rkcom3.prscr],&work[rkcom3.prstgs],notinintp,wrkint,
            &wrkint[startp]);
      }
//C
//C  Set markers to show that interpolation has been initialized for
//C  NWANT components.
      nwntsv = nwant;
      inintp = true;
   }
//C
//C  The actual evaluation of the interpolating polynomial and/or its first
//C  derivative is done in EVALI.
//C
   if (rkcom4.methd == 1) {
      evali(&work[rkcom3.pry],&work[rkcom3.pryp],&work[rkcom3.prstgs],twant,reqest,
         nwant,ywant,ypwant);
   }
   else {
      evali(&work[rkcom3.pry],&work[rkcom3.pryp],&wrkint[startp],twant,reqest,
         nwant,ywant,ypwant);
   }
//C
   label20:
//C
   rkmsg(ier,srname,nrec,flag);
//C
   return;
}

void RKSUITE::reset(double tendnu)
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  If you are not familiar with the code RESET and how it is used in
//C  conjunction with CT to solve initial value problems, you should study the
//C  document file rksuite.doc carefully before attempting to use the code. The
//C  following "Brief Reminder" is intended only to remind you of the meaning,
//C  type, and size requirements of the arguments.
//C
//C  The integration using CT proceeds from TSTART in the direction of TEND, and
//C  is now at TNOW.  To reset TEND to a new value TENDNU, just call RESET with
//C  TENDNU as the argument.  You must continue integrating in the same
//C  direction, so the sign of (TENDNU - TNOW) must be the same as the sign of
//C  (TEND - TSTART). To change direction you must restart by a call to SETUP.
//C
//C  INPUT VARIABLE
//C
//C     TENDNU    - double PRECISION
//C                 The new value of TEND.
//C                 Constraint: TENDNU and TNOW must be clearly distinguishable
//C                 in the precision used.  The sign of (TENDNU - TNOW) must be
//C                 the same as the sign of (TEND - TSTART).
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C     .. Scalar Arguments ..
//      double PRECISION  TENDNU
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Common Block for Integrator Options ..
//      LOGICAL           MSG, UTASK
//      COMMON /RKCOM8/   MSG, UTASK
//      SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      CHARACTER*6       SRNAME
//      PARAMETER         (SRNAME='RESET')
   const char srname[] = "RESET";
//      LOGICAL           ASK
//      INTEGER           MINUS1, MINUS2
//      PARAMETER         (ASK=.TRUE.,MINUS1=-1,MINUS2=-2)
   const bool ask = true;
   const int minus1 = -1, minus2 = -2;
//      double PRECISION  ZERO
//      PARAMETER         (ZERO=0.0D+0)
   const double zero = 0.0;
//C     .. Local Scalars ..
   double hmin, tdiff;
   int flag, ier, nrec, state, state1;
//C     .. External Subroutines ..
//      EXTERNAL          RKMSG, RKSIT
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX
//C     .. Executable Statements ..
   ier = 1;
   nrec = 0;
//C
//C  Is it permissible to call RESET?
//C
   rksit(ask,"CT",state);
   if (state == 911) {
      ier = 912;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** A catastrophic error has already been detected elsewhere.");
      goto label20;
   }
   if (rkcom8.utask) {
      rksit(ask,"UT",state1);
      if (state1 != minus2) {
         ier = 911;
         nrec = 2;
         sprintf(&rkcom9.rec[0][0]," ** You have called RESET after you specified to SETUP that");
         sprintf(&rkcom9.rec[1][0]," ** you were going to use UT. This is not permitted.");
         goto label20;
      }
   }
   if (state == minus1) {
      ier = 911;
      nrec = 1;
      sprintf(&rkcom9.rec[0][0]," ** You have not called CT, so you cannot use RESET.");
      goto label20;
   }
   if (state == 5 || state == 6) {
      ier = 911;
      nrec = 2;
      sprintf(&rkcom9.rec[0][0]," ** CT has returned with CFLAG =  %6d.",state);
      sprintf(&rkcom9.rec[1][0]," ** You cannot call RESET in this circumstance.");
      goto label20;
   }
//C
//C  Check value of TENDNU
//C
   if (rkcom1.dir > zero && tendnu <= rkcom2.t) {
      ier = 911;
      nrec = 4;
      sprintf(&rkcom9.rec[0][0]," ** Integration is proceeding in the positive direction. The");
      sprintf(&rkcom9.rec[1][0]," ** current value for the independent variable is %13.5lf",rkcom2.t);
      sprintf(&rkcom9.rec[2][0]," ** and you have set TENDNU = %13.5lf.  TENDNU must be",tendnu);
      sprintf(&rkcom9.rec[3][0]," ** greater than T.");
   }
   else if (rkcom1.dir < zero && tendnu >= rkcom2.t) {
      ier = 911;
      nrec = 4;
      sprintf(&rkcom9.rec[0][0]," ** Integration is proceeding in the negative direction. The");
      sprintf(&rkcom9.rec[1][0]," ** current value for the independent variable is %13.5lf",rkcom2.t);
      sprintf(&rkcom9.rec[2][0]," ** and you have set TENDNU = %13.5lf.  TENDNU must be",tendnu);
      sprintf(&rkcom9.rec[3][0]," ** less than T.");
   }
   else {
      hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(rkcom2.t),fabs(tendnu)));
      tdiff = fabs(tendnu-rkcom2.t);
      if (tdiff < hmin) {
         ier = 911;
         nrec = 4;
         sprintf(&rkcom9.rec[0][0]," ** The current value of the independent variable T is %13.5lf.",rkcom2.t);
         sprintf(&rkcom9.rec[1][0]," ** The TENDNU you supplied has ABS(TENDNU-T) = %13.5lf.",tdiff);
         sprintf(&rkcom9.rec[2][0]," ** For the METHOD and the precision of the computer being");
         sprintf(&rkcom9.rec[3][0]," ** used, this difference must be at least %13.5lf.",hmin);
      }
   }
   if (ier == 911) goto label20;
//C
   rkcom1.tnd = tendnu;
   rkcom2.last = false;
//C
   label20:
//C
   rkmsg(ier,srname,nrec,flag);
//C
   return;
}

void RKSUITE::mconst(int method)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:   Sets machine-dependent global quantities
//C
//C  Common:    Initializes:    /RKCOM7/ OUTCH, MCHEPS, DWARF, RNDOFF,
//C                                      SQRRMC, CUBRMC, TINY
//C             Reads:          none
//C             Alters:         none
//C
//C  Comments:
//C  =========
//C  OUTCH, MCHEPS, DWARF are pure environmental parameters with values
//C  obtained from a call to ENVIRN. The other quantities set depend on
//C  the environmental parameters, the implementation, and, possibly,
//C  METHOD. At present the METHODs implemented in the RK suite do not
//C  influence the values of these quantities.
//C  OUTCH  - Standard output channel
//C  MCHEPS - Largest positive number such that 1.0D0 + MCHEPS = 1.0D0
//C  DWARF  - Smallest positive number
//C  RNDOFF - 10 times MCHEPS
//C  SQRRMC - Square root of MCHEPS
//C  CUBRMC - Cube root of MCHEPS
//C  TINY   - Square root of DWARF
//C
//C     .. Scalar Arguments ..
//      INTEGER           METHOD
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Parameters ..
//      double PRECISION  TEN, THIRD
//      PARAMETER         (TEN=10.0D+0,THIRD=1.0D+0/3.0D+0)
   const double ten = 10.0, third = 1.0/3.0;
//C     .. External Subroutines ..
//      EXTERNAL          ENVIRN
//C     .. Intrinsic Functions ..
//      INTRINSIC         SQRT
//C     .. Executable Statements ..
//C
   envirn(rkcom7.outch,rkcom7.mcheps,rkcom7.dwarf);
//C
   rkcom7.rndoff = ten*rkcom7.mcheps;
   rkcom7.sqrrmc = sqrt(rkcom7.mcheps);
   rkcom7.cubrmc = pow(rkcom7.mcheps,third);
   rkcom7.tiny = sqrt(rkcom7.dwarf);
//C
   return;
}

void RKSUITE::envirn(int& outch, double& mcheps, double& dwarf)
{
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C  The RK suite requires some environmental parameters that are provided by
//C  this subroutine.  The values provided with the distribution codes are those
//C  appropriate to the IEEE standard.  They must be altered, if necessary, to
//C  those appropriate to the computing system you are using before calling the 
//C  codes of the suite.  
//C
//C       ================================================================
//C       ================================================================
//C        TO MAKE SURE THAT THESE MACHINE AND INSTALLATION DEPENDENT 
//C        QUANTITIES ARE SPECIFIED PROPERLY, THE DISTRIBUTION VERSION 
//C        WRITES A MESSAGE ABOUT THE MATTER TO THE STANDARD OUTPUTCHANNEL
//C        AND TERMINATES THE RUN.  THE VALUES PROVIDED IN THEDISTRIBUTION
//C        VERSION SHOULD BE ALTERED, IF NECESSARY, AND THE "WRITE" AND 
//C        "STOP" STATEMENTS COMMENTED OUT.
//C       ================================================================
//C       ================================================================
//C
//C  OUTPUT VARIABLES 
//C
//C     OUTCH     - INTEGER
//C                 Standard output channel
//C     MCHEPS    - double PRECISION
//C                 MCHEPS is the largest positive number such that
//C                 1.0D0 + MCHEPS = 1.0D0. 
//C     DWARF     - double PRECISION
//C                 DWARF is the smallest positive number.
//C
//C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//C
//C     .. Scalar Arguments ..  
//      INTEGER           OUTCH
//      double PRECISION  DWARF, MCHEPS
//C     .. Executable Statements ..      
//C
//C  The following six statements are to be Commented out after verification that
//C  the machine and installation dependent quantities are specified correctly.
//C  If you pass copies of RKSUITE on to others, please give them the whole
//C  distribution version of RKSUITE, and in particular, give them a version 
//C  of ENVIRN that does not have the following six statements Commented out.
//      WRITE(*,*) ' Before using RKSUITE, you must verify that the  '
//      WRITE(*,*) ' machine- and installation-dependent quantities  '
//      WRITE(*,*) ' specified in the subroutine ENVIRN are correct, '
//      WRITE(*,*) ' and then Comment these WRITE statements and the '
//      WRITE(*,*) ' STOP statement out of ENVIRN.                   '
//      STOP
//C
//C  The following values are appropriate to IEEE arithmetic with the typical
//C  standard output channel.
//C
      outch = 6;
//    mcheps = 1.11e-16;
//    dwarf = 2.23e-308;
      mcheps = DBL_EPSILON;
      dwarf = DBL_MIN;
//      sprintf(&rkcom9.rec[0][0],"mcheps = %13.5lf dwarf = %13.5lf",mcheps,dwarf);
//C      
//C------------------------------------------------------------------------------
//C  If you have the routines D1MACH and I1MACH on your system, you could
//C  replace the preceding statements by the following ones to obtain the 
//C  appropriate machine dependent numbers. The routines D1MACH and I1MACH 
//C  are public domain software.  They are available from NETLIB.
//C      .. Scalar Arguments ..  
//C      INTEGER           OUTCH
//C      double PRECISION  DWARF, MCHEPS
//C      .. External Functions ..
//C      INTEGER           I1MACH
//C      double PRECISION  D1MACH
//C      .. Executable Statements ..
//C
//C      OUTCH = I1MACH(2)
//C      MCHEPS = D1MACH(3)
//C      DWARF = D1MACH(1)
//C
//C  If you have the NAG Fortran Library available on your system, you could 
//C  replace the preceding statements by the following ones to obtain the 
//C  appropriate machine dependent numbers.
//C
//C      .. Scalar Arguments ..  
//C      INTEGER           OUTCH
//C      double PRECISION  DWARF, MCHEPS
//C      .. External Functions ..
//C      double PRECISION  X02AJF, X02AMF
//C      .. Executable Statements ..
//C
//C      CALL X04AAF(0,OUTCH)
//C      MCHEPS = X02AJF()
//C      DWARF = X02AMF()
//C
//C  If you have the IMSL MATH/LIBRARY available on your system, you could
//C  replace the preceding statements by the following ones to obtain the
//C  appropriate machine dependent numbers.
//C
//C      .. Scalar Arguments ..  
//C      INTEGER           OUTCH
//C      double PRECISION  DWARF, MCHEPS
//C      .. External Functions ..
//C      double PRECISION  DMACH
//C      .. Executable Statements ..
//C
//C      CALL UMACH(2,OUTCH)
//C      MCHEPS = DMACH(4)
//C      DWARF = DMACH(1)
//C------------------------------------------------------------------------------
//C
   return;
}

void RKSUITE::rkconst(int method, int& vecstg, bool& reqstg, int& lintpl)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//*************************************************
//C
//C  Purpose:   Set formula definitions and formula characteristics for
//C             selected method. Return storage requirements for the
//C             selected method.
//C
//C  Input:     METHOD
//C  Output:    VECSTG, REQSTG, LINTPL
//C
//C  Common:    Initializes:    /RKCOM4/ A(*,*), B(*), C(*), BHAT(*), R(*),
//C                                      E(*), PTR(*), NSTAGE, METHD, INTP, MINTP
//C                             /RKCOM5/ TOOSML, COST, SAFETY, EXPON, STBRAD,
//C                                      TANANG, RS, RS1, RS2, RS3, RS4, ORDER,
//C                                      LSTSTG, MAXTRY, NSEC, FSAL
//C             Reads:          /RKCOM7/ RNDOFF
//C             Alters:         none
//C
//C  Comments:
//C  =========
//C  Runge-Kutta formula pairs are described by a set of coefficients
//C  and by the setting of a number of parameters that describe the
//C  characteristics of the pair.  The input variable METHD indicates
//C  which of the three pairs available is to be set. In the case of
//C  METHD = 2 additional coefficients are defined that make interpolation
//C  of the results possible and provide an additional error estimator.
//C  VECSTG is the number of columns of workspace required to compute the
//C  stages of a METHD. For interpolation purposes the routine returns via
//C  COMMON the logical variable INTP indicating whether interpolation is
//C  possible and via the call list:
//C  REQSTG - whether the stages are required to form the
//C           interpolant (set .FALSE. if INTP=.FALSE.)
//C  LINTPL - the number of extra columns of storage required for use
//C           with UT (set 0 if INTP=.FALSE.)
//C
//C  Quantities set in common blocks:
//C  METHD - copy of METHOD
//C  A, B, C, BHAT - coefficients of the selected method
//C  R      - extra coefficents for interpolation with METHD = 2
//C  E      - extra coefficients for additional local error estimate
//C           with METHD = 2
//C  PTR    - vector of pointers indicating how individual stages are to
//C           be stored.  With it zero coefficients of the formulas can
//C           be exploited to reduce the storage required
//C  NSTAGE - number of stages for the specified METHD
//C  INTP   - indicates whether there is an associated interpolant
//C           (depending on the method, the user may have to supply
//C           extra workspace)
//C  MINTP  - the degree of the interpolating polynomial, if one exists
//C  FSAL   - indicates whether the last stage of a step can be used as
//C           the first stage of the following step
//C  LSTSTG - pointer to location of last stage for use with FSAL=.TRUE.
//C  ORDER  - the lower order of the pair of Runge-Kutta formulas that
//C           constitute a METHD
//C  TANANG, 
//C  STBRAD - the stability region of the formula used to advance
//C           the integration is approximated by a sector in the left half
//C           complex plane.  TANANG is the tangent of the interior angle
//C           of the sector and STBRAD is the radius of the sector.
//C  COST   - cost of a successful step in function evaluations
//C  MAXTRY - limit on the number of iterations in the stiffness check. As
//C           set, no more than 24 function evaluations are made in the check.
//C  NSEC   - each step of size H in the primary integration corresponds to
//C           NSEC steps of size H/NSEC in the secondary integration when
//C           global error assessment is done.
//C  EXPON  - used to adjust the step size; this code implements an error
//C           per step control for which EXPON = 1/(ORDER + 1).
//C  SAFETY - quantity used in selecting the step size
//C  TOOSML - quantity used to determine when a step size is too small for
//C           the precision available
//C  RS, RS1,
//C  RS2, RS3,
//C  RS4    - quantities used in determining the maximum and minimum change
//C           change in step size (set independently of METHD)
//C
//C  Further comments on SAFETY:
//C ============================
//C  The code estimates the largest step size that will yield the specified
//C  accuracy.  About half the time this step size would result in a local
//C  error that is a little too large, and the step would be rejected.  For
//C  this reason a SAFETY factor is used to reduce the "optimal" value to one
//C  more likely to succeed.  Unfortunately, there is some art in choosing this
//C  value. The more expensive a failed step, the smaller one is inclined to 
//C  take this factor. However, a small factor means that if the prediction were
//C  good, more accuracy than desired would be obtained and the behavior of the
//C  error would then be irregular.  The more stringent the tolerance, the better
//C  the prediction, except near limiting precision. Thus the general range of 
//C  tolerances expected influences the choice of SAFETY.
//C
//C     .. Scalar Arguments ..
//      INTEGER           LINTPL, METHOD, VECSTG
//      LOGICAL           REQSTG
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Parameters ..
//      double PRECISION  ONE, ZERO, TWO, FIFTY, FIVEPC
//      PARAMETER         (ONE=1.0D+0,ZERO=0.0D+0,TWO=2.0D+0,FIFTY=50.D+0,
//     &                  FIVEPC=0.05D+0)
   const double one = 1.0, zero = 0.0, two = 2.0, fifty = 50.0, fivepc = 0.05;
//C     .. Local Scalars ..
   double cdiff, diff;
   int i, j;
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, DBLE, INT, MAX, MIN
//C     .. Executable Statements ..
//C
   rkcom4.methd = method;
//C
//    GO TO (20,40,100) rkcom4.methd
   switch (rkcom4.methd) {
   case 1:
//C
//C  METHD = 1.
//C    This pair is from "A 3(2) Pair of Runge-Kutta Formulas" by P. Bogacki
//C    and L.F. Shampine, Appl. Math. Lett., 2, pp. 321-325, 1989.  The authors
//C    are grateful to P. Bogacki for his assistance in implementing the pair.
//C
//   20 CONTINUE
      rkcom4.nstage = 4;
      rkcom5.fsal = true;
      rkcom5.order = 2;
      rkcom5.tanang = 8.9e0;
      rkcom5.stbrad = 2.3e0;
      rkcom5.safety = 0.8e0;
      rkcom4.intp = true;
      rkcom4.mintp = 3;
      reqstg = false;
      lintpl = 2;
      rkcom5.nsec = 3;
//C
      rkcom4.ptr[1] = 0;
      rkcom4.ptr[2] = 1;
      rkcom4.ptr[3] = 2;
      rkcom4.ptr[4] = 3;
//C
      rkcom4.a[2][1] = 1.0e0/2.0e0;
      rkcom4.a[3][1] = 0.0e0;
      rkcom4.a[3][2] = 3.0e0/4.0e0;
      rkcom4.a[4][1] = 2.0e0/9.0e0;
      rkcom4.a[4][2] = 1.0e0/3.0e0;
      rkcom4.a[4][3] = 4.0e0/9.0e0;
//C
//C  The coefficients BHAT(*) refer to the formula used to advance the
//C  integration, here the one of order 3.  The coefficients B(*) refer
//C  to the other formula, here the one of order 2. For this pair, BHAT(*)
//C  is not needed since FSAL = .TRUE.
//C
      rkcom4.b[1] = 7.0e0/24.0e0;
      rkcom4.b[2] = 1.0e0/4.0e0;
      rkcom4.b[3] = 1.0e0/3.0e0;
      rkcom4.b[4] = 1.0e0/8.0e0;
//C
      rkcom4.c[1] = 0.0e0;
      rkcom4.c[2] = 1.0e0/2.0e0;
      rkcom4.c[3] = 3.0e0/4.0e0;
      rkcom4.c[4] = 1.0e0;
//C
//      GO TO 120
      break;
   case 2:
//C
//C  METHD = 2
//C    This pair is from "An Efficient Runge-Kutta (4,5) Pair" by P. Bogacki
//C    and L.F. Shampine, Rept. 89-20, Math. Dept., Southern Methodist
//C    University, Dallas, Texas, USA, 1989.  The authors are grateful to
//C    P. Bogacki for his assistance in implementing the pair.  Shampine and
//C    Bogacki subsequently modified the formula to enhance the reliability of
//C    the pair.  The original fourth order formula is used in an estimate of
//C    the local error.  If the step fails, the computation is broken off.  If
//C    the step is acceptable, the first evaluation of the next step is done,
//C    i.e., the pair is implemented as FSAL and the local error of the step
//C    is again estimated with a fourth order formula using the additional data.
//C    The step must succeed with both estimators to be accepted.  When the
//C    second estimate is formed, it is used for the subsequent adjustment of
//C    the step size because it is of higher quality.  The two fourth order
//C    formulas are well matched to leading order, and only exceptionally do
//C    the estimators disagree -- problems with discontinuous coefficients are
//C    handled more reliably by using two estimators as is global error
//C    estimation.
//C
//   40 CONTINUE
      rkcom4.nstage = 8;
      rkcom5.fsal = true;
      rkcom5.order = 4;
      rkcom5.tanang = 5.2e0;
      rkcom5.stbrad = 3.9e0;
      rkcom5.safety = 0.8e0;
      rkcom4.intp = true;
      reqstg = true;
      rkcom4.mintp = 6;
      lintpl = 6;
      rkcom5.nsec = 2;
//C
      rkcom4.ptr[1] = 0;
      rkcom4.ptr[2] = 1;
      rkcom4.ptr[3] = 2;
      rkcom4.ptr[4] = 3;
      rkcom4.ptr[5] = 4;
      rkcom4.ptr[6] = 5;
      rkcom4.ptr[7] = 6;
      rkcom4.ptr[8] = 7;
//C
      rkcom4.a[2][1] = 1.0e0/6.0e0;
      rkcom4.a[3][1] = 2.e0/27.e0;
      rkcom4.a[3][2] = 4.e0/27.e0;
      rkcom4.a[4][1] = 183.e0/1372.e0;
      rkcom4.a[4][2] = -162.e0/343.e0;
      rkcom4.a[4][3] = 1053.e0/1372.e0;
      rkcom4.a[5][1] = 68.e0/297.e0;
      rkcom4.a[5][2] = -4.e0/11.e0;
      rkcom4.a[5][3] = 42.e0/143.e0;
      rkcom4.a[5][4] = 1960.e0/3861.e0;
      rkcom4.a[6][1] = 597.e0/22528.e0;
      rkcom4.a[6][2] = 81.e0/352.e0;
      rkcom4.a[6][3] = 63099.e0/585728.e0;
      rkcom4.a[6][4] = 58653.e0/366080.e0;
      rkcom4.a[6][5] = 4617.e0/20480.e0;
      rkcom4.a[7][1] = 174197.e0/959244.e0;
      rkcom4.a[7][2] = -30942.e0/79937.e0;
      rkcom4.a[7][3] = 8152137.e0/19744439.e0;
      rkcom4.a[7][4] = 666106.e0/1039181.e0;
      rkcom4.a[7][5] = -29421.e0/29068.e0;
      rkcom4.a[7][6] = 482048.e0/414219.e0;
      rkcom4.a[8][1] = 587.e0/8064.e0;
      rkcom4.a[8][2] = 0.e0;
      rkcom4.a[8][3] = 4440339.e0/15491840.e0;
      rkcom4.a[8][4] = 24353.e0/124800.e0;
      rkcom4.a[8][5] = 387.e0/44800.e0;
      rkcom4.a[8][6] = 2152.e0/5985.e0;
      rkcom4.a[8][7] = 7267.e0/94080.e0;
//C
//C  The coefficients B(*) refer to the formula of order 4.
//C
      rkcom4.b[1] = 2479.e0/34992.e0;
      rkcom4.b[2] = 0.e0;
      rkcom4.b[3] = 123.e0/416.e0;
      rkcom4.b[4] = 612941.e0/3411720.e0;
      rkcom4.b[5] = 43.e0/1440.e0;
      rkcom4.b[6] = 2272.e0/6561.e0;
      rkcom4.b[7] = 79937.e0/1113912.e0;
      rkcom4.b[8] = 3293.e0/556956.e0;
//C
//C  The coefficients E(*) refer to an estimate of the local error based on
//C  the first formula of order 4.  It is the difference of the fifth order
//C  result, here located in A(8,*), and the fourth order result.  By
//C  construction both E(2) and E(7) are zero.
//C
      rkcom4.e[1] = -3.e0/1280.e0;
      rkcom4.e[2] = 0.e0;
      rkcom4.e[3] = 6561.e0/632320.e0;
      rkcom4.e[4] = -343.e0/20800.e0;
      rkcom4.e[5] = 243.e0/12800.e0;
      rkcom4.e[6] = -1.e0/95.e0;
      rkcom4.e[7] = 0.e0;
//C
      rkcom4.c[1] = 0.e0;
      rkcom4.c[2] = 1.e0/6.e0;
      rkcom4.c[3] = 2.e0/9.e0;
      rkcom4.c[4] = 3.e0/7.e0;
      rkcom4.c[5] = 2.e0/3.e0;
      rkcom4.c[6] = 3.e0/4.e0;
      rkcom4.c[7] = 1.e0;
      rkcom4.c[8] = 1.e0;
//C
//C  To do interpolation with this pair, some extra stages have to be computed.
//C  The following additional A(*,*) and C(*) coefficients are for this purpose.
//C  In addition there is an array R(*,*) that plays a role for interpolation
//C  analogous to that of BHAT(*) for the basic step.
//C
      rkcom4.c[9] = 1.e0/2.e0;
      rkcom4.c[10] = 5.e0/6.e0;
      rkcom4.c[11] = 1.e0/9.e0;
//C
      rkcom4.a[9][1] = 455.e0/6144.e0;
      rkcom4.a[10][1] = -837888343715.e0/13176988637184.e0;
      rkcom4.a[11][1] = 98719073263.e0/1551965184000.e0;
      rkcom4.a[9][2] = 0.e0;
      rkcom4.a[10][2] = 30409415.e0/52955362.e0;
      rkcom4.a[11][2] = 1307.e0/123552.e0;
      rkcom4.a[9][3] = 10256301.e0/35409920.e0;
      rkcom4.a[10][3] = -48321525963.e0/759168069632.e0;
      rkcom4.a[11][3] = 4632066559387.e0/70181753241600.e0;
      rkcom4.a[9][4] = 2307361.e0/17971200.e0;
      rkcom4.a[10][4] = 8530738453321.e0/197654829557760.e0;
      rkcom4.a[11][4] = 7828594302389.e0/382182512025600.e0;
      rkcom4.a[9][5] = -387.e0/102400.e0;
      rkcom4.a[10][5] = 1361640523001.e0/1626788720640.e0;
      rkcom4.a[11][5] = 40763687.e0/11070259200.e0;
      rkcom4.a[9][6] = 73.e0/5130.e0;
      rkcom4.a[10][6] = -13143060689.e0/38604458898.e0;
      rkcom4.a[11][6] = 34872732407.e0/224610586200.e0;
      rkcom4.a[9][7] = -7267.e0/215040.e0;
      rkcom4.a[10][7] = 18700221969.e0/379584034816.e0;
      rkcom4.a[11][7] = -2561897.e0/30105600.e0;
      rkcom4.a[9][8] = 1.e0/32.e0;
      rkcom4.a[10][8] = -5831595.e0/847285792.e0;
      rkcom4.a[11][8] = 1.e0/10.e0;
      rkcom4.a[10][9] = -5183640.e0/26477681.e0;
      rkcom4.a[11][9] = -1.e0/10.e0;
      rkcom4.a[11][10] = -1403317093.e0/11371610250.e0;
//C
      for (i = 1; i <= 11; i++) {
         rkcom4.r[i][1] = 0.e0;
      }
      for (i = 1; i <= 6; i++) {
         rkcom4.r[2][i] = 0.e0;
      }
      rkcom4.r[1][6] = -12134338393.e0/1050809760.e0;
      rkcom4.r[1][5] = -1620741229.e0/50038560.e0;
      rkcom4.r[1][4] = -2048058893.e0/59875200.e0;
      rkcom4.r[1][3] = -87098480009.e0/5254048800.e0;
      rkcom4.r[1][2] = -11513270273.e0/3502699200.e0;
//C
      rkcom4.r[3][6] = -33197340367.e0/1218433216.e0;
      rkcom4.r[3][5] = -539868024987.e0/6092166080.e0;
      rkcom4.r[3][4] = -39991188681.e0/374902528.e0;
      rkcom4.r[3][3] = -69509738227.e0/1218433216.e0;
      rkcom4.r[3][2] = -29327744613.e0/2436866432.e0;
//C
      rkcom4.r[4][6] = -284800997201.e0/19905339168.e0;
      rkcom4.r[4][5] = -7896875450471.e0/165877826400.e0;
      rkcom4.r[4][4] = -333945812879.e0/5671036800.e0;
      rkcom4.r[4][3] = -16209923456237.e0/497633479200.e0;
      rkcom4.r[4][2] = -2382590741699.e0/331755652800.e0;
//C
      rkcom4.r[5][6] = -540919.e0/741312.e0;
      rkcom4.r[5][5] = -103626067.e0/43243200.e0;
      rkcom4.r[5][4] = -633779.e0/211200.e0;
      rkcom4.r[5][3] = -32406787.e0/18532800.e0;
      rkcom4.r[5][2] = -36591193.e0/86486400.e0;
//C
      rkcom4.r[6][6] = 7157998304.e0/374350977.e0;
      rkcom4.r[6][5] = 30405842464.e0/623918295.e0;
      rkcom4.r[6][4] = 183022264.e0/5332635.e0;
      rkcom4.r[6][3] = -3357024032.e0/1871754885.e0;
      rkcom4.r[6][2] = -611586736.e0/89131185.e0;
//C
      rkcom4.r[7][6] = -138073.e0/9408.e0;
      rkcom4.r[7][5] = -719433.e0/15680.e0;
      rkcom4.r[7][4] = -1620541.e0/31360.e0;
      rkcom4.r[7][3] = -385151.e0/15680.e0;
      rkcom4.r[7][2] = -65403.e0/15680.e0;
//C
      rkcom4.r[8][6] = 1245.e0/64.e0;
      rkcom4.r[8][5] = 3991.e0/64.e0;
      rkcom4.r[8][4] = 4715.e0/64.e0;
      rkcom4.r[8][3] = 2501.e0/64.e0;
      rkcom4.r[8][2] = 149.e0/16.e0;
      rkcom4.r[8][1] = 1.e0;
//C
      rkcom4.r[9][6] = 55.e0/3.e0;
      rkcom4.r[9][5] = 71.e0;
      rkcom4.r[9][4] = 103.e0;
      rkcom4.r[9][3] = 199.e0/3.e0;
      rkcom4.r[9][2] = 16.0e0;
//C
      rkcom4.r[10][6] = -1774004627.e0/75810735.e0;
      rkcom4.r[10][5] = -1774004627.e0/25270245.e0;
      rkcom4.r[10][4] = -26477681.e0/359975.e0;
      rkcom4.r[10][3] = -11411880511.e0/379053675.e0;
      rkcom4.r[10][2] = -423642896.e0/126351225.e0;
//C
      rkcom4.r[11][6] = 35.e0;
      rkcom4.r[11][5] = 105.e0;
      rkcom4.r[11][4] = 117.e0;
      rkcom4.r[11][3] = 59.e0;
      rkcom4.r[11][2] = 12.e0;
//C
//      GO TO 120
      break;
//C
//C  METHD = 3
//C    This pair is from "High Order Embedded Runge-Kutta Formulae" by P.J.
//C    Prince and J.R. Dormand, J. Comp. Appl. Math.,7, pp. 67-75, 1981.  The
//C    authors are grateful to P. Prince and J. Dormand for their assistance in
//C    implementing the pair.
//C
//  100 CONTINUE
   case 3:
      rkcom4.nstage = 13;
      rkcom5.fsal = false;
      rkcom5.order = 7;
      rkcom5.tanang = 11.0e0;
      rkcom5.stbrad = 5.2e0;
      rkcom5.safety = 0.8e0;
      rkcom4.intp = false;
      reqstg = false;
      rkcom4.mintp = 0;
      lintpl = 0;
      rkcom5.nsec = 2;
//C
      rkcom4.ptr[1] = 0;
      rkcom4.ptr[2] = 1;
      rkcom4.ptr[3] = 2;
      rkcom4.ptr[4] = 1;
      rkcom4.ptr[5] = 3;
      rkcom4.ptr[6] = 2;
      rkcom4.ptr[7] = 4;
      rkcom4.ptr[8] = 5;
      rkcom4.ptr[9] = 6;
      rkcom4.ptr[10] = 7;
      rkcom4.ptr[11] = 8;
      rkcom4.ptr[12] = 9;
      rkcom4.ptr[13] = 1;
//C
      rkcom4.a[2][1] = 5.55555555555555555555555555556e-2;
      rkcom4.a[3][1] = 2.08333333333333333333333333333e-2;
      rkcom4.a[3][2] = 6.25e-2;
      rkcom4.a[4][1] = 3.125e-2;
      rkcom4.a[4][2] = 0.e0;
      rkcom4.a[4][3] = 9.375e-2;
      rkcom4.a[5][1] = 3.125e-1;
      rkcom4.a[5][2] = 0.e0;
      rkcom4.a[5][3] = -1.171875e0;
      rkcom4.a[5][4] = 1.171875e0;
      rkcom4.a[6][1] = 3.75e-2;
      rkcom4.a[6][2] = 0.e0;
      rkcom4.a[6][3] = 0.e0;
      rkcom4.a[6][4] = 1.875e-1;
      rkcom4.a[6][5] = 1.5e-1;
      rkcom4.a[7][1] = 4.79101371111111111111111111111e-2;
      rkcom4.a[7][2] = 0.e0;
      rkcom4.a[7][3] = 0.0e0;
      rkcom4.a[7][4] = 1.12248712777777777777777777778e-1;
      rkcom4.a[7][5] = -2.55056737777777777777777777778e-2;
      rkcom4.a[7][6] = 1.28468238888888888888888888889e-2;
      rkcom4.a[8][1] = 1.6917989787292281181431107136e-2;
      rkcom4.a[8][2] = 0.e0;
      rkcom4.a[8][3] = 0.e0;
      rkcom4.a[8][4] = 3.87848278486043169526545744159e-1;
      rkcom4.a[8][5] = 3.59773698515003278967008896348e-2;
      rkcom4.a[8][6] = 1.96970214215666060156715256072e-1;
      rkcom4.a[8][7] = -1.72713852340501838761392997002e-1;
      rkcom4.a[9][1] = 6.90957533591923006485645489846e-2;
      rkcom4.a[9][2] = 0.e0;
      rkcom4.a[9][3] = 0.e0;
      rkcom4.a[9][4] = -6.34247976728854151882807874972e-1;
      rkcom4.a[9][5] = -1.61197575224604080366876923982e-1;
      rkcom4.a[9][6] = 1.38650309458825255419866950133e-1;
      rkcom4.a[9][7] = 9.4092861403575626972423968413e-1;
      rkcom4.a[9][8] = 2.11636326481943981855372117132e-1;
      rkcom4.a[10][1] = 1.83556996839045385489806023537e-1;
      rkcom4.a[10][2] = 0.e0;
      rkcom4.a[10][3] = 0.e0;
      rkcom4.a[10][4] = -2.46876808431559245274431575997e0;
      rkcom4.a[10][5] = -2.91286887816300456388002572804e-1;
      rkcom4.a[10][6] = -2.6473020233117375688439799466e-2;
      rkcom4.a[10][7] = 2.84783876419280044916451825422e0;
      rkcom4.a[10][8] = 2.81387331469849792539403641827e-1;
      rkcom4.a[10][9] = 1.23744899863314657627030212664e-1;
      rkcom4.a[11][1] = -1.21542481739588805916051052503e0;
      rkcom4.a[11][2] = 0.e0;
      rkcom4.a[11][3] = 0.e0;
      rkcom4.a[11][4] = 1.66726086659457724322804132886e1;
      rkcom4.a[11][5] = 9.15741828416817960595718650451e-1;
      rkcom4.a[11][6] = -6.05660580435747094755450554309e0;
      rkcom4.a[11][7] = -1.60035735941561781118417064101e1;
      rkcom4.a[11][8] = 1.4849303086297662557545391898e1;
      rkcom4.a[11][9] = -1.33715757352898493182930413962e1;
      rkcom4.a[11][10] = 5.13418264817963793317325361166e0;
      rkcom4.a[12][1] = 2.58860916438264283815730932232e-1;
      rkcom4.a[12][2] = 0.e0;
      rkcom4.a[12][3] = 0.e0;
      rkcom4.a[12][4] = -4.77448578548920511231011750971e0;
      rkcom4.a[12][5] = -4.3509301377703250944070041181e-1;
      rkcom4.a[12][6] = -3.04948333207224150956051286631e0;
      rkcom4.a[12][7] = 5.57792003993609911742367663447e0;
      rkcom4.a[12][8] = 6.15583158986104009733868912669e0;
      rkcom4.a[12][9] = -5.06210458673693837007740643391e0;
      rkcom4.a[12][10] = 2.19392617318067906127491429047e0;
      rkcom4.a[12][11] = 1.34627998659334941535726237887e-1;
      rkcom4.a[13][1] = 8.22427599626507477963168204773e-1;
      rkcom4.a[13][2] = 0.e0;
      rkcom4.a[13][3] = 0.e0;
      rkcom4.a[13][4] = -1.16586732572776642839765530355e1;
      rkcom4.a[13][5] = -7.57622116690936195881116154088e-1;
      rkcom4.a[13][6] = 7.13973588159581527978269282765e-1;
      rkcom4.a[13][7] = 1.20757749868900567395661704486e1;
      rkcom4.a[13][8] = -2.12765911392040265639082085897e0;
      rkcom4.a[13][9] = 1.99016620704895541832807169835e0;
      rkcom4.a[13][10] = -2.34286471544040292660294691857e-1;
      rkcom4.a[13][11] = 1.7589857770794226507310510589e-1;
      rkcom4.a[13][12] = 0.e0;
//C
//C  The coefficients BHAT(*) refer to the formula used to advance the
//C  integration, here the one of order 8.  The coefficients B(*) refer
//C  to the other formula, here the one of order 7.
//C
      rkcom4.bhat[1] = 4.17474911415302462220859284685e-2;
      rkcom4.bhat[2] = 0.e0;
      rkcom4.bhat[3] = 0.e0;
      rkcom4.bhat[4] = 0.e0;
      rkcom4.bhat[5] = 0.e0;
      rkcom4.bhat[6] = -5.54523286112393089615218946547e-2;
      rkcom4.bhat[7] = 2.39312807201180097046747354249e-1;
      rkcom4.bhat[8] = 7.0351066940344302305804641089e-1;
      rkcom4.bhat[9] = -7.59759613814460929884487677085e-1;
      rkcom4.bhat[10] = 6.60563030922286341461378594838e-1;
      rkcom4.bhat[11] = 1.58187482510123335529614838601e-1;
      rkcom4.bhat[12] = -2.38109538752862804471863555306e-1;
      rkcom4.bhat[13] = 2.5e-1;
//C
      rkcom4.b[1] = 2.9553213676353496981964883112e-2;
      rkcom4.b[2] = 0.e0;
      rkcom4.b[3] = 0.e0;
      rkcom4.b[4] = 0.e0;
      rkcom4.b[5] = 0.e0;
      rkcom4.b[6] = -8.28606276487797039766805612689e-1;
      rkcom4.b[7] = 3.11240900051118327929913751627e-1;
      rkcom4.b[8] = 2.46734519059988698196468570407e0;
      rkcom4.b[9] = -2.54694165184190873912738007542e0;
      rkcom4.b[10] = 1.44354858367677524030187495069e0;
      rkcom4.b[11] = 7.94155958811272872713019541622e-2;
      rkcom4.b[12] = 4.44444444444444444444444444445e-2;
      rkcom4.b[13] = 0.e0;
//C
      rkcom4.c[1] = 0.e0;
      rkcom4.c[2] = 5.55555555555555555555555555556e-2;
      rkcom4.c[3] = 8.33333333333333333333333333334e-2;
      rkcom4.c[4] = 1.25e-1;
      rkcom4.c[5] = 3.125e-1;
      rkcom4.c[6] = 3.75e-1;
      rkcom4.c[7] = 1.475e-1;
      rkcom4.c[8] = 4.65e-1;
      rkcom4.c[9] = 5.64865451382259575398358501426e-1;
      rkcom4.c[10] = 6.5e-1;
      rkcom4.c[11] = 9.24656277640504446745013574318e-1;
      rkcom4.c[12] = 1.e0;
      rkcom4.c[13] = rkcom4.c[12];
//C
//      GO TO 120
      break;
//C
//C  The definitions of all pairs come here for the calculation of
//C  LSTSTG, RS1, RS2, RS3, RS4, COST, MAXTRY, EXPON, TOOSML, and VECSTG.
//C
//  120 CONTINUE
   }
   rkcom5.lststg = rkcom4.ptr[rkcom4.nstage];
   if (rkcom5.fsal) {
      rkcom5.cost = double(rkcom4.nstage-1);
   }
   else {
      rkcom5.cost = double(rkcom4.nstage);
   }
//C
//C  MAXTRY - limit on the number of iterations of a computation made in
//C  diagnosing stiffness.  There are at most Q = 3 function calls per
//C  iteration. MAXTRY is determined so that  Q*MAXTRY <= 5% of the cost of
//C  50 steps and 1 <= MAXTRY <= 8. This limits the number of calls to FCN
//C  in each diagnosis of stiffness to 24 calls.
//C
   rkcom5.maxtry = min(8,max(1,int(fivepc*rkcom5.cost*fifty)));
//C
   rkcom5.expon = one/(rkcom5.order+one);
//C
//C     In calculating CDIFF it is assumed that there will be a non-zero
//C     difference |C(I) - C(J)| less than one. If C(I) = C(J) for any I not
//C     equal to J, they should be made precisely equal by assignment.
//C
   cdiff = one;
   for (i = 1; i < rkcom4.nstage; i++) {
      for (j = i + 1; j <= rkcom4.nstage; j++) {
         diff = fabs(rkcom4.c[i]-rkcom4.c[j]);
         if (diff != zero) cdiff = min(cdiff,diff);
      }
   }
   rkcom5.toosml = rkcom7.rndoff/cdiff;
//C
//C  Determine the number of columns needed in STAGES(1:NEQ,*) (this will be
//C  at most NSTAGE-1 since the first stage is held in a separate array).
//C  The PTR array contains the column positions of the stages.
//C
   vecstg = 0;
   for (i = 2; i <= rkcom4.nstage; i++) {
      vecstg = max(rkcom4.ptr[i],vecstg);
   }
//C
   rkcom5.rs = two;
   rkcom5.rs1 = one/rkcom5.rs;
   rkcom5.rs2 = rkcom5.rs*rkcom5.rs;
   rkcom5.rs3 = rkcom5.rs*rkcom5.rs*rkcom5.rs;
   rkcom5.rs4 = one/rkcom5.rs3;
//C
   return;
}

void RKSUITE::formi(void (*f)(double, double*, double*), int neq, int nwant, double y[], double yp[],
   double yold[], double ypold[], double stages/*[neq][]*/[], bool calstg,
   double xstage[], double p[])
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:    Forms an interpolating polynomial for use with
//C              METHDs 1 or 2.
//C
//C  Input:      NEQ, NWANT, T, Y(*), YP(*), HOLD, YOLD(*), YPOLD(*),
//C              STAGES(NEQ,*), CALSTG
//C  Output:     P(*), XSTAGE(NEQ)
//C  External:   F
//C
//C  Common:     Initializes:    none
//C              Reads:          /RKCOM4/ A(*,*), C(*), R(*), METHD, MINTP
//C                              /RKCOM2/ T, TOLD, HOLD
//C              Alters:         /RKCOM2/ NFCN
//C
//C  Comments:
//C  =========
//C  The integration has reached T with a step HOLD from TOLD = T-HOLD.
//C  Y(*),YP(*) and YOLD(*),YPOLD(*) approximate the solution and its
//C  derivative at T and TOLD respectively.  STAGES(NEQ,*) holds the stages
//C  computed in taking this step. In the case of METHD = 2 it is necessary 
//C  to compute some more stages in this subroutine. CALSTG indicates whether
//C  or not the extra stages need to be computed. A(*,*) and C(*) are used in 
//C  computing these stages. The extra stages are stored in STAGES(NEQ,*) and 
//C  XSTAGE(*).  The coefficients of the interpolating polynomials for the first
//C  NWANT components of the solution are returned in the array P(*). The 
//C  polynomial is of degree MINTP = 3 for METHD = 1 and of degree MINTP = 6 
//C  for METHD = 2. The vector R(*) is used for workspace when METHD = 2.
//C
//C     .. Scalar Arguments ..
//      INTEGER           NEQ, NWANT
//      LOGICAL           CALSTG
//C     .. Array Arguments ..
//      double PRECISION  P(*), STAGES(NEQ,*), XSTAGE(*), Y(*), YOLD(*),
//     &                  YP(*), YPOLD(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Local Scalars ..
   double d1, d2, d3, d4, hyp, hypold;
   int i, j, k, l;
//C     .. Executable Statements ..
//C
   if (rkcom4.methd == 1) {
//C
//C  METHD = 1.  Use the cubic Hermite interpolant that is is fully
//C  specified by the values and slopes at the two ends of the step.
//C
      for (l = 0; l < nwant; l++) {
         d1 = y[l] - yold[l];
         hyp = rkcom2.hold*yp[l];
         hypold = rkcom2.hold*ypold[l];
         d2 = hyp - d1;
         d3 = d1 - hypold;
         d4 = d2 - d3;
         p[l] = d2 + d4;
         p[nwant+l] = d4;
      }
//C
   }
   else {
//C
//C  METHD = 2.
//C
      if (calstg) {
//C
//C  Compute the extra stages needed for interpolation using the facts that
//C       1. Stage 1 is YPOLD(*).
//C       2. Stage i (i>1) is stored in STAGES(1:NEQ,i).
//C       3. This pair is FSAL, i.e. STAGES(1:NEQ,7)=YP(1:NEQ), which frees
//C          up STAGES(1:NEQ,7) for use by stage 9.
//C       4. XSTAGE(1:NEQ) is used for stage 10.
//C       5. The coefficient of stage 2 in the interpolant is always 0, so
//C          STAGES(1:NEQ,1) is used for stage 11.
//C  The vector P(1:NEQ) is used as workspace for computing the stages.
//C
         for (i = 9; i <= 11; i++) {
            for (l = 0; l < neq; l++) {
               p[l] = rkcom4.a[i][1]*ypold[l];
            }
            for (j = 2; j < i; j++) {
               if (j <= 7) {
                  for (l = 0; l < neq; l++) {
                     p[l] = p[l] + rkcom4.a[i][j]*
                        stages/*[l][j-1]*/[l+neq*(j-1-1)];
                  }
               }
               else if (j == 8) {
                  for (l = 0; l < neq; l++) {
                     p[l] = p[l] + rkcom4.a[i][j]*yp[l];
                  }
               }
               else if (j == 9) {
                  for (l = 0; l < neq; l++) {
                     p[l] = p[l] + rkcom4.a[i][j]*
                        stages/*[l][7]*/[l+neq*(7-1)];
                  }
               }
               else if (j == 10) {
                  for (l = 0; l < neq; l++) {
                     p[l] = p[l] + rkcom4.a[i][j]*xstage[l];
                  }
               }
            }
            for (l = 0; l < neq; l++) {
               p[l] = yold[l] + rkcom2.hold*p[l];
            }
            if (i == 9) {
               f(rkcom2.told+rkcom4.c[i]*rkcom2.hold,p,
                  &stages/*[1][7]*/[neq*(7-1)]);
               rkcom2.nfcn++;
            }
            else if (i == 10) {
               f(rkcom2.told+rkcom4.c[i]*rkcom2.hold,p,xstage);
               rkcom2.nfcn++;
            }
            else {
               f(rkcom2.told+rkcom4.c[i]*rkcom2.hold,p,
                  &stages/*[1][1]*/[neq*(1-1)]);
               rkcom2.nfcn++;
            }
         }
      }
//C
//C  Form the coefficients of the interpolating polynomial in its shifted
//C  and scaled form.  The transformation from the form in which the
//C  polynomial is derived can be somewhat ill-conditioned.  The terms 
//C  are grouped so as to minimize the errors of the transformation.
//C
//C  Coefficient of SIGMA**6
      k = 4*nwant;
      for (l = 0; l < nwant; l++) {
         p[k+l] = rkcom4.r[5][6]*stages/*[l][4]*/[l+neq*(4-1)] +
                  ((rkcom4.r[10][6]*xstage[l]+rkcom4.r[8][6]*yp[l])+
                  (rkcom4.r[7][6]*stages/*[l][6]*/[l+neq*5]+rkcom4.r[6][6]*stages/*[l][5]*/[l+neq*4])) +
                  ((rkcom4.r[4][6]*stages/*[l][3]*/[l+neq*2]+rkcom4.r[9][6]*stages/*[l][7]*/[l+neq*6])+
                  (rkcom4.r[3][6]*stages/*[l][2]*/[l+neq]+rkcom4.r[11][6]*stages/*[l][1]*/[l])+
                  rkcom4.r[1][6]*ypold[l]);
      }
//C
//C  Coefficient of SIGMA**5
      k = 3*nwant;
      for (l = 0; l < nwant; l++) {
         p[k+l] = (rkcom4.r[10][5]*xstage[l]+rkcom4.r[9][5]*stages/*[l][7]*/[l+neq*6]) +
                  ((rkcom4.r[7][5]*stages/*[l][6]*/[l+neq*5]+rkcom4.r[6][5]*stages/*[l][5]*/[l+neq*4])+
                  rkcom4.r[5][5]*stages/*[l][4]*/[l+neq*3]) + ((rkcom4.r[4][5]*stages/*[l][3]*/[l+neq*2]+
                  rkcom4.r[8][5]*yp[l])+(rkcom4.r[3][5]*stages/*[l][2]*/[l+neq]+rkcom4.r[11][5]*
                  stages/*[l][1]*/[l])+rkcom4.r[1][5]*ypold[l]);
      }
//C
//C  Coefficient of SIGMA**4
      k = 2*nwant;
      for (l = 0; l < nwant; l++) {
         p[k+l] = ((rkcom4.r[4][4]*stages/*[l][3]*/[l+neq*2]+rkcom4.r[8][4]*yp[l])+
                  (rkcom4.r[7][4]*stages/*[l][6]*/[l+neq*5]+rkcom4.r[6][4]*stages/*[l][5]*/[l+neq*4])+
                  rkcom4.r[5][4]*stages/*[l][4]*/[l+neq*3]) + ((rkcom4.r[10][4]*xstage[l]+
                  rkcom4.r[9][4]*stages/*[l][7]*/[l+neq*6])+(rkcom4.r[3][4]*stages/*[l][2]*/[l+neq]+
                  rkcom4.r[11][4]*stages/*[l][1]*/[l])+rkcom4.r[1][4]*ypold[l]);
      }
//C
//C  Coefficient of SIGMA**3
      k = nwant;
      for (l = 0; l < nwant; l++) {
         p[k+l] = rkcom4.r[5][3]*stages/*[l][4]*/[l+neq*3] + rkcom4.r[6][3]*stages/*[l][5]*/[l+neq*4] +
                  ((rkcom4.r[3][3]*stages/*[l][2]*/[l+neq]+rkcom4.r[9][3]*stages/*[l][7]*/[l+neq*6])+
                  (rkcom4.r[10][3]*xstage[l]+rkcom4.r[8][3]*yp[l])+rkcom4.r[1][3]*
                  ypold[l])+((rkcom4.r[4][3]*stages/*[l][3]*/[l+neq*2]+rkcom4.r[11][3]*
                  stages/*[l][1]*/[l])+rkcom4.r[7][3]*stages/*[l][6]*/[l+neq*5]);
      }
//C
//C  Coefficient of SIGMA**2
//C
      for (l = 0; l < nwant; l++) {
         p[l] = rkcom4.r[5][2]*stages/*[l][4]*/[l+neq*3] + ((rkcom4.r[6][2]*stages/*[l][5]*/[l+neq*4]+
                rkcom4.r[8][2]*yp[l])+rkcom4.r[1][2]*ypold[l]) +
                ((rkcom4.r[3][2]*stages/*[l][2]*/[l+neq]+rkcom4.r[9][2]*stages/*[l][7]*/[l+neq*6])+
                rkcom4.r[10][2]*xstage[l]) + ((rkcom4.r[4][2]*stages/*[l][3]*/[l+neq*2]+
                rkcom4.r[11][2]*stages/*[l][1]*/[l])+rkcom4.r[7][2]*stages/*[l][6]*/[l+neq*5]);
      }
//C
//C  Scale all the coefficients by the step size.
//C
      for (l = 0; l < nwant*(rkcom4.mintp-1); l++) {
         p[l] = rkcom2.hold*p[l];
      }
//C
   }
//C
   return;
}

void RKSUITE::evali(double y[], double yp[], double p /* [nwant][] */ [], double twant,
   char reqest[], int nwant, double ywant[], double ypwant[])
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//*************************************************
//C
//C  Purpose:    Evaluation of an interpolating polynomial and/or its
//C              first derivative.
//C
//C  Input:      Y(*), YP(*), P(NWANT,*), TWANT, REQEST, NWANT
//C  Output:     YWANT(*), YPWANT(*)
//C
//C  Common:     Initializes:    none
//C              Reads:          /RKCOM2/ HOLD, T
//C                              /RKCOM4/ MINTP
//C              Alters:         none
//C
//C  Comments:
//C  =========
//C  The interpolant is evaluated at TWANT to approximate the solution,
//C  YWANT, and/or its first derivative there, YPWANT. Only the first
//C  NWANT components of the answer are computed. There are three cases
//C  that are indicated by the first character of REQEST:
//C    REQEST(1:1) = 'S' or 's'- compute approximate 'S'olution only.
//C                = 'D' or 'd'- compute approximate first 'D'erivative
//C                              of the solution only.
//C                = 'B' or 'b'- compute 'B'oth approximate solution and
//C                              first derivative.
//C  The coefficents of the polynomial are contained in Y(*), YP(*) and
//C  P(NWANT,*).
//C
//C     .. Scalar Arguments ..
//      double PRECISION  TWANT
//      INTEGER           NWANT
//      CHARACTER*(*)     REQEST
//C     .. Array Arguments ..
//      double PRECISION  P(NWANT,*), Y(*), YP(*), YPWANT(*), YWANT(*)
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Local Scalars ..
   double sigma;
   int k, l;
   char reqst1;
//C     .. Executable Statements ..
//C
//C  Evaluate the interpolating polynomial of degree MINTP in terms of the
//C  shifted and scaled independent variable SIGMA.
//C
   sigma = (twant-rkcom2.t)/rkcom2.hold;
//C
   reqst1 = reqest[0];
   if (reqst1 == 'S' || reqst1 == 's' ||
      reqst1 == 'B' || reqst1 == 'b') {
//C
      for (l = 0; l < nwant; l++) {
         ywant[l] = p /* [l][rkcom4.mintp-1] */ [l+nwant*(rkcom4.mintp-2)]*sigma;
      }
      for (k = rkcom4.mintp - 2; k >= 1; k--) {
         for (l = 0; l < nwant; l++) {
            ywant[l] = (ywant[l]+p/*[l][k]*/[l+nwant*(k-1)])*sigma;
         }
      }
      for (l = 0; l < nwant; l++) {
         ywant[l] = (ywant[l]+rkcom2.hold*yp[l])*sigma + y[l];
      }
   }
//C
//C  Evaluate the derivative of the interpolating polynomial.
//C
   if (reqst1 == 'D' || reqst1 == 'd' ||
      reqst1 == 'B' || reqst1 == 'b') {
//C
//C  The derivative of the interpolating polynomial with respect to TWANT 
//C  is the derivative with respect to S divided by HOLD.
//C
      for (l = 0; l < nwant; l++) {
         ypwant[l] = rkcom4.mintp*
         p/*[l][rkcom4.mintp-1]*/[l+nwant*(rkcom4.mintp-2)]*sigma;
      }
      for (k = rkcom4.mintp -1; k >= 2; k--) {
         for (l = 0; l < nwant; l++) {
            ypwant[l] = (ypwant[l]+k*p/*[l][k-1]*/[l+nwant*(k-2)])*sigma;
         }
      }
      for (l = 0; l < nwant; l++) {
         ypwant[l] = (ypwant[l]+rkcom2.hold*yp[l])/rkcom2.hold;
      }
   }
//C
   return;
}

void RKSUITE::rkmsg(int ier, const char* srname, int nrec, int& flag)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      To process error messages and terminate the program
//C                in the event of a "catastrophic" failure.
//C
//C  Input:        IER, SRNAME, NREC
//C  Output:       FLAG
//C
//C  Common:       Initializes:    none
//C                Reads:          /RKCOM7/ OUTCH
//C                                /RKCOM8/ MSG, UTASK
//C                                /RKCOM9/ REC
//C                Alters:         none
//C
//C  Comments:
//C  =========
//C  The output variable FLAG is assigned the value of the input variable IER.
//C
//C  IER = -2  reports a successful call of the subroutine SRNAME and
//C            indicates that special action is to be taken elsewhere
//C            in the suite.  FLAG is set and a return is effected.
//C
//C  IER = 1   reports a successful call of the subroutine SRNAME.  FLAG
//C            is set and a return is effected.
//C
//C  1 < IER < 911 and MSG = .TRUE.: a message of NREC records contained in
//C            the array REC(*) is written to the standard output channel, 
//C            OUTCH.  FLAG is set and a return is effected.
//C
//C  IER = 911 reports a "catastrophic" error was detected in SRNAME.  A
//C            message is written to OUTCH regardless of the value of MSG and
//C            normally the execution of the program is terminated.  The
//C            execution is not terminated if the error is the result of an
//C            indirect call to CT, RESET, or INTRP through UT (UTASK = .TRUE.).
//C            Termination can be prevented by using the subroutine SOFTFL.
//C
//C  IER = 912 reports that a "catastrophic" error was detected earlier and
//C            termination was prevented, but the user has failed to take
//C            appropriate remedial action.  Execution is terminated.
//C
//C     .. Scalar Arguments ..
//      INTEGER           FLAG, IER, NREC
//      CHARACTER*(*)     SRNAME
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Common Block for Integrator Options ..
//      LOGICAL           MSG, UTASK
//      COMMON /RKCOM8/   MSG, UTASK
//      SAVE   /RKCOM8/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      INTEGER           PLUS1
//      LOGICAL           ASK, TELL
//      PARAMETER         (PLUS1=1,ASK=.TRUE.,TELL=.FALSE.)
   const int plus1 = 1;
   const bool ask = true, tell = false;
//C     .. Local Scalars ..
   int i;
   bool baderr, ok, on, utcall;
//C     .. External Subroutines ..
//      EXTERNAL          CHKFL, RKSIT, SOFTFL
//C     .. Executable Statements ..
//C
//C  Check where the call came from - if it is an indirect call from UT,
//C  the run is not STOPped.
   on = false;
   utcall = (!strcmp(srname,"RESET") || !strcmp(srname,"CT") ||
      !strcmp(srname,"INTRP")) && rkcom8.utask;
//C
//C  Check if can continue with integrator.
   ok = (!strcmp(srname,"CT") || !strcmp(srname,"UT")) &&
      (ier == 2 || ier == 3 || ier == 4);
//C
//C  Check if program termination has been overridden.
//   softel(ask,on);
//C
   if ((rkcom8.msg && ier > plus1) || ier >= 911) {
         printf("\n **\n");
         for (i=1; i<=nrec; i++) printf("%s\n",rkcom9.rec[i-1]);
         if (ier >= 911) {
            printf(" **\n");
            printf(" ** Catastrophic error detected in %s.\n",srname);
            printf(" **\n");
            if ((!(utcall || on) && ier == 911) ||
               ier == 912) {
               printf(" **\n");
               printf(" ** Execution of your program is being terminated.\n");
               printf(" **\n");
               exit(-1);
            }
         }
         else if (ok) {
            printf(" **\n");
            printf(" ** Warning from routine %s with flag set %6d.\n",srname,ier);
            printf(" ** You can continue integrating this problem.\n");
            printf(" **\n");
         }
         else {
            printf(" **\n");
            printf(" ** Warning from routine %s with flag set %6d.\n",srname,ier);
            printf(" ** You cannot continue integrating this problem.\n");
            printf(" **\n");
         }
      }
      for (i = nrec; i < 10; i++) {
         rkcom9.rec[i][0] = '\0';
   }
      flag = ier;
//C
//C  TELL RKSIT the status of the routine associated with SRNAME
      rksit(tell,srname,flag);
//C
//C  Indicate that a catastrophic error has been detected
      baderr = flag >= 911;
      chkfl(tell,baderr);
//C
      return;
//C
}

void RKSUITE::rksit(bool ask, const char* srname, int& state)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      To save or enquire about the status of each
//C                subprogram in the suite.
//C
//C  Input:        ASK, SRNAME
//C  Input/output: STATE
//C
//C
//C  Comments:
//C  =========
//C  SRNAME indicates which routine is of interest in the call to RKSIT.
//C
//C  If ASK=.FALSE., then the value of STATE (which as used is the error
//C  flag) for the routine SRNAME is saved internally. This value of STATE
//C  is usually positive.  There are two exceptions:
//C     1. SRNAME='SETUP' and STATE=-1 indicates a completely new problem,
//C        so all SAVEd states are cleared.
//C     2. STATE=-2 is used by some routines in the suite to indicate
//C        circumstances which require special action.
//C
//C  If ASK=.TRUE., then RKSIT first checks to see if there were any
//C  catastrophic errors, that is, a SAVEd state has value 911. This should
//C  happen only when the user has overridden program termination in the event
//C  of catastrophic failures from routines in the package but has failed to
//C  take appropriate action. If this is the case, then RKSIT returns a value
//C  of STATE = 911 which forces a termination of execution inside RKMSG. If 
//C  no catastrophic errors are flagged, then STATE returns the saved state 
//C  value for the routine specified by SRNAME.
//C
//C     .. Scalar Arguments ..
//      INTEGER           STATE
//      LOGICAL           ASK
//      CHARACTER*(*)     SRNAME
//C     .. Parameters ..
//      INTEGER           STATES, MINUS1
//      PARAMETER         (STATES=7,MINUS1=-1)
   const int states = 7, minus1 = -1;
//C     .. Local Scalars ..
   int i, name;
//C     .. Local Arrays ..
   static int svsta[states] = {minus1,minus1,minus1,minus1,minus1,minus1,minus1};
//C     .. Save statement ..
//      SAVE              SVSTA
//C     .. Data statements ..
//      DATA              SVSTA/STATES*MINUS1/
//C     .. Executable Statements ..
//C
   if (!strcmp(srname,"SETUP")) {
      name = 1;
   }
   else if (!strcmp(srname,"UT")) {
      name = 2;
   }
   else if (!strcmp(srname,"STAT")) {
      name = 3;
   }
   else if (!strcmp(srname,"GLBERR")) {
      name = 4;
   }
   else if (!strcmp(srname,"CT")) {
      name = 5;
   }
   else if (!strcmp(srname,"INTRP")) {
      name = 6;
   }
   else if (!strcmp(srname,"RESET")) {
      name = 7;
   }
   else {
      name = 0;
   }
//C
//C  (Re)initialize if SETUP is telling RKSIT to do so.
   if (!ask && name == 1 && state == minus1) {
      for (i = 0; i < states; i++) {
         svsta[i] = minus1;
      }
      goto label60;
   }
//C
//C  Check for 911 on exit from a previous call.
   if (ask) {
      for (i = 0; i < states; i++) {
         if (svsta[i] == 911) {
            state = 911;
            goto label60;
         }
      }
   }
//C
   if (ask) {
      state = svsta[name-1];
   }
   else {
      svsta[name-1] = state;
   }
//C
   label60:
//C
   return;
}

void RKSUITE::truerr(void (*f)(double, double*, double*), int neq, double y[],
   double tol, double weight[], double zy[], double zyp[], double zerror[],
   double zynew[], double zerres[], double zstage /* [neq][] */ [], int& ier)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      Compute a running RMS measure of the true (global) error
//C                for a general Runge-Kutta pair.
//C
//C
//C  Input:        NEQ, Y(*), TOL, WEIGHT(*),
//C  Input/output: ZY(*), ZYP(*), ZERROR(*)
//C  Workspace:    ZYNEW(*), ZERRES(*), ZSTAGE(NEQ,*)
//C  Output:       IER
//C  External:     F
//C
//C  Common:       Initializes:    none
//C                Reads:          /RKCOM2/ T, HOLD
//C                                /RKCOM5/ TOOSML, ORDER, NSEC
//C                                /RKCOM7/ TINY
//C                Alters:         /RKCOM6/ MAXERR, LOCMAX, GNFCN
//C
//C  Comments:
//C  =========
//C  A secondary integration is performed using a fraction of the step size 
//C  of the primary integration. ZY(*) and ZYP(*) are the approximate solution
//C  and first derivative of this secondary integration. ZERRES(*) contains the 
//C  error estimates for the secondary integration. ZYNEW(*) and ZSTAGE(*,*) are
//C  workspace for taking a step. The error assessment is computed using the
//C  difference of the primary and secondary solutions at the primary
//C  integration points as an estimate of the true error there.  The weights 
//C  used are those of the error test of the primary integration. This error 
//C  assessment is maintained in the vector ZERROR(*).  MAXERR and LOCMAX 
//C  contain the maximum contribution to the assessment and its location,
//C  respectively.  The number of calls to F is counted by GNFCN.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  TOL
//      INTEGER           IER, NEQ
//C     .. Array Arguments ..
//      double PRECISION  WEIGHT(*), Y(*), ZERRES(*), ZERROR(*),
//     &                  ZSTAGE(NEQ,*), ZY(*), ZYNEW(*), ZYP(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Global Error Assessment ..
//      double PRECISION  MAXERR, LOCMAX
//      INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
//     &                  PRZYNU
//      LOGICAL           ERASON, ERASFL
//      COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
//     &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
//      SAVE   /RKCOM6/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Parameters ..
//      double PRECISION  PT1, TEN, DUMMY
//      PARAMETER         (PT1=0.1D0,TEN=10.0D0,DUMMY=1.0D0)
   const double pt1 = 0.1, ten = 10.0, dummy = 1.0;
//C     .. Local Scalars ..
   double diff, errmax, hmin, hsec, mxerlc, tsec, zlerr,
      ztest1, ztest2;
   int istep, l, level;
   bool ldummy, main;
//C     .. Local Arrays ..
   double dumarr[1];
//C     .. External Subroutines ..
//      EXTERNAL          STEP
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, DBLE, MAX
//C     .. Executable Statements ..
   tsec = rkcom2.t - rkcom2.hold;
   hsec = rkcom2.hold/double(rkcom5.nsec);
   hmin = max(rkcom7.tiny,rkcom5.toosml*max(fabs(tsec),fabs(rkcom2.t)));
   if (fabs(hsec) < hmin) {
      ier = 6;
      goto label120;
   }
   ztest1 = tol/double(rkcom5.nsec);
   ztest2 = tol/ten;
   level = 0;
//C
//C  The subroutine STEP is used to take a step.  In its use in the primary
//C  integration provision is made for getting on scale in the first step.
//C  In this situation only the subroutine might reduce the step size.  By
//C  setting MAIN = .FALSE., the subroutine will take a step of the size input.
//C  In this use of the subroutine, all items of the call list appearing after
//C  MAIN are dummy variables.
//C
//C  Perform secondary integration.
   main = false;
   ldummy = false;
   for (istep = 1; istep <= rkcom5.nsec; istep++) {
//C
//C  Take a step.
      step(f,neq,tsec,zy,zyp,zstage,ztest1,hsec,weight,zynew,
         zerres,zlerr,main,dummy,dumarr,ldummy);
//C
//C  The primary integration is using a step size of HUSED and the secondary
//C  integration is using the smaller step size HSEC = HUSED/NSEC.  If steps
//C  of this size were taken from the same starting point and the asymptotic
//C  behavior were evident, the smaller step size would result in a local error
//C  that is considerably smaller, namely by a factor of 1/(NSEC**(ORDER+1)).
//C  If the two approximate solutions are close and TOLR is neither too large nor
//C  too small, this should be approximately true.  The step size is chosen in
//C  the primary integration so that the local error ERR is no larger than TOLR.
//C  The local error, ZLERR, of the secondary integration is compared to TOLR in
//C  an attempt to diagnose a secondary integration that is not rather more
//C  accurate than the primary integration.
//C
      if (zlerr >= ztest1) {
         level = 2;
      }
      else if (zlerr > ztest2) {
         level++;
      }
      if (level >= 2) {
         ier = 6;
         goto label120;
      }
//C
//C  Advance TSEC and the dependent variables ZY(*) and ZYP(*).
      tsec = rkcom2.t - double(rkcom5.nsec-istep)*hsec;
      for (l = 0; l < neq; l++) {
         zy[l] = zynew[l];
      }
//C
      if (rkcom5.fsal) {
//C
//C  When FSAL = .TRUE., the derivative ZYP(*) is the last stage of the step.
         for (l = 0; l < neq; l++) {
            zyp[l] = zstage /* [l][rkcom5.lststg] */ [l+neq*(rkcom5.lststg-1)];
         }
      }
      else {
//C
//C  Call F to evaluate ZYP(*).
         f(tsec,zy,zyp);
         rkcom6.gnfcn++;
      }
//C
   }
//C
//C  Update the maximum error seen, MAXERR, and its location, LOCMAX.
//C  Use local variables ERRMAX and MXERLC.
//C
   errmax = rkcom6.maxerr;
   mxerlc = rkcom6.locmax;
   for (l = 0; l < neq; l++) {
      diff = fabs(zy[l]-y[l])/weight[l];
      if (diff > errmax) {
         errmax = diff;
         mxerlc = rkcom2.t;
      }
   }
//C
//C  If the global error is greater than 0.1D0, the solutions have diverged so
//C  far that comparing them may not provide a reliable estimate of the global
//C  error. The test is made before ZERROR(*) and MAXERR, LCMXER are updated so
//C  that on a failure, they refer to the last reliable results.
//C
   if (errmax > pt1) {
      ier = 6;
      goto label120;
   }
   else {
      rkcom6.maxerr = errmax;
      rkcom6.locmax = mxerlc;
      for (l = 0; l < neq; l++) {
         diff = fabs(zy[l]-y[l])/weight[l];
         zerror[l] = zerror[l] + diff * diff;
      }
      ier = 1;
   }
//C
//C  Exit point for TRUERR
   label120:
//C
   return;
}

void RKSUITE::step(void (*f)(double, double*, double*), int neq, double tnow,
   double* y, double* yp, double stages /* [neq][] */ [], double tol, double& htry,
   double* weight, double* ynew, double* errest, double& err, bool main,
   double hmin, double* thres, bool& phase2)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      To compute a step of an explicit Runge-Kutta
//C                method and estimate the local error of the step.
//C
//C  Input:        NEQ, TNOW, Y(*), YP(*), TOL, MAIN, HMIN, THRES(*)
//C  Input/output: HTRY, PHASE2, LAST, WEIGHT(*)
//C  Output:       STAGES(NEQ,*), YNEW(*), ERREST(*), ERR
//C
//C  Common:       Initializes:    none
//C                Reads:          /RKCOM1/ TND
//C                                /RKCOM2/ LAST
//C                                /RKCOM4/ A, B, C, BHAT, PTR, NSTAGE, METHD
//C                                /RKCOM5/ FSAL
//C                Alters:         /RKCOM2/ NFCN, LAST
//C                                /RKCOM6/ GNFCN
//C
//C  Comments:
//C  =========
//C  From an approximate solution Y(*) at TNOW and first derivative there,
//C  YP(*) = F(TNOW,Y,YP), a step is taken to get an approximation YNEW(*)
//C  at TNOW + HTRY. The Runge-Kutta method and how it is used are defined
//C  by A, B, C, BHAT, PTR, NSTAGE, METHD and FSAL. Intermediate stages
//C  of the method are stored in the array STAGES(NEQ,*). The error in
//C  each solution component is estimated and returned in ERREST(*). A
//C  weighted maximum norm of the local error, ERR, is formed. For some
//C  methods an intermediate error estimate can be computed before completion
//C  of the step (see routine STEPB); if the estimate is greater than the
//C  specified tolerance TOL, the computation of the step is terminated.
//C
//C  When global error estimation is desired, two integrations are done.
//C  The usual integration is referred to as the "primary", or "main",
//C  integration (MAIN=.TRUE.).  For global error estimation another,
//C  "secondary" integration (MAIN=.FALSE.) is carried out with a smaller
//C  step size.  The weight vector WEIGHT(*) used in computing ERR is
//C  determined by the main integration.  Thus this argument is output when
//C  MAIN = .TRUE. and input when MAIN = .FALSE..
//C
//C  When taking the first step in an integration, the logical variable
//C  PHASE2 may be input as .TRUE. and if the first step is the whole of
//C  the range of integration, then LAST will be .TRUE.. When PHASE2=.TRUE.,
//C  the first three stages are monitored to help assure that the step
//C  size H is small enough for the integration to be stable and for the
//C  estimate of the error of the step to be credible. Calls are made to
//C  the subroutine STEPA for this purpose. If necessary, H will be
//C  reduced in STEPA (and LAST altered accordingly) and the step retried
//C  in STEP until an acceptable value is found.
//C
//C  In the primary integration the number of calls to F is counted by
//C  NFCN, and in the secondary integration, by GNFCN.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  ERR, HMIN, HTRY, TNOW, TOL
//      INTEGER           NEQ
//      LOGICAL           MAIN, PHASE2
//C     .. Array Arguments ..
//      double PRECISION  ERREST(*), STAGES(NEQ,*), THRES(*), WEIGHT(*),
//     &                  Y(*), YNEW(*), YP(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Global Error Assessment ..
//      double PRECISION  MAXERR, LOCMAX
//      INTEGER           GNFCN, PRZSTG, PRZY, PRZYP, PRZERS, PRZERR,
//     &                  PRZYNU
//      LOGICAL           ERASON, ERASFL
//      COMMON /RKCOM6/   MAXERR, LOCMAX, GNFCN, PRZSTG, PRZY, PRZYP,
//     &                  PRZERS, PRZERR, PRZYNU, ERASON, ERASFL
//      SAVE   /RKCOM6/
//C     .. Parameters ..
//      double PRECISION  ZERO, HALF, ONE
//      PARAMETER         (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
   const double zero = 0.0, half = 0.5, one = 1.0;
//C     .. Local Scalars ..
   double avgy, tstg;
   int i, j, l;
   bool cutbak;
//C     .. External Subroutines ..
//      EXTERNAL          STEPA, STEPB
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX, SIGN
//C     .. Executable Statements ..
//C
//C  Many of the following loops over L = 1, NEQ have constant array values
//C  inside. The code is written with clarity in mind.  Any optimizing
//C  compiler will identify these occurrences and take appropriate action.
//C  A check for zero multipliers has been included so as to prevent
//C  needless computation resulting from the storing of zero coefficients
//C  in the arrays for the sake of clarity.  The array ERREST(*) is used
//C  for working storage in this computation.
//C
   label20:
   if (main) {
      if (phase2) {
//C
//C  Initialize weights for measuring the local error.
         for (l = 0; l < neq; l++) {
            weight[l] = max(thres[l],fabs(y[l]));
         }
      }
   }
//C
   for (i = 2; i <= rkcom4.nstage; i++) {
      for (j = 1; j < i; j++) {
         if (j == 1) {
            for (l = 0; l < neq; l++) {
               errest[l] = rkcom4.a[i][1]*yp[l];
            }
         }
         else {
            if (rkcom4.a[i][j] != zero) {
               for (l = 0; l < neq; l++) {
                  errest[l] = errest[l] + rkcom4.a[i][j]* /* stages[l][rkcom4.ptr[j]] */
                  stages[l+neq*(rkcom4.ptr[j]-1)];
               }
            }
         }
      }
      for (l = 0; l < neq; l++) {
         ynew[l] = y[l] + htry*errest[l];
      }
//C
//C  METHD = 2 is special in that an estimate of the local error can be
//C  formed before the step is completed.  If the step is a failure,
//C  return immediately.  Otherwise, complete the step and compute a more
//C  accurate error estimate.
      if (rkcom4.methd == 2 && i == 7) {
         stepb(neq,y,yp,htry,ynew,stages,thres,err,main,weight);
         if (err > tol) return;
      }
//C
      tstg = tnow + rkcom4.c[i]*htry;
      if (main && rkcom2.last && rkcom4.c[i] == one) tstg = rkcom1.tnd;
      f(tstg,ynew,/* stages[1][rkcom4.ptr[i]] */ &stages[neq*(rkcom4.ptr[i]-1)]);
//C
//C  Increment the counter for the number of function evaluations
//C  depending on whether the primary or secondary integration is taking
//C  place.
      if (main) {
         rkcom2.nfcn++;
      }
      else {
         rkcom6.gnfcn++;
      }
//C
//C----------------------------------------------------------------------
//C  When PHASE2 is .TRUE. we are in the second phase of the automatic
//C  selection of the initial step size.  The results of the first three
//C  stages are monitored in the subroutine STEPA for evidence that H is
//C  too large -- instability and/or an unreliable estimate of the error
//C  of the step is then possible.  When the subroutine believes H to be
//C  too large, it returns CUTBAK = .TRUE. and a suitably reduced H for
//C  another try.
//C
      if (main) {
         if (phase2) {
            if (i <= 3 && fabs(htry) > hmin) {
               stepa(tnow,y,yp,tstg,ynew,/* stages[1][rkcom4.ptr[i]] */
                  &stages[neq*(rkcom4.ptr[i]-1)],
                  htry,weight,cutbak);
               if (cutbak) {
                  rkcom2.last = false;
//C
//C  Make sure that STEPA does not reduce the step size below the
//C  minimum. If it does, reset H to HMIN and deactivate PHASE2.
                  if (fabs(htry) <= hmin) {
                     htry = sign(hmin,htry);
                     phase2 = false;
                  }
                  goto label20;
               }
            }
         }
      }
//C----------------------------------------------------------------------
//C
   }
//C
//C  Some formulas are constructed so that the last stage represents
//C  the result of the step (FSAL=.TRUE.), hence if the step is acceptable,
//C  it will be the first stage for the next step. When FSAL=.FALSE., we
//C  have to complete the computation of the step.
//C
   if (!rkcom5.fsal) {
      for (i = 1; i <= rkcom4.nstage; i++) {
         if (i == 1) {
            for (l = 0; l < neq; l++) {
               errest[l] = rkcom4.bhat[1]*yp[l];
            }
         }
         else {
            if (rkcom4.bhat[i] != zero) {
               for (l = 0; l < neq; l++) {
                  errest[l] = errest[l] + rkcom4.bhat[i]* /* stages[l][rkcom4.ptr[i]] */
                  stages[l+neq*(rkcom4.ptr[i]-1)];
               }
            }
         }
      }
      for (l = 0; l < neq; l++) {
         ynew[l] = y[l] + htry*errest[l];
      }
   }
//C
//C  Form an estimate of the error in the lower order formula by comparing
//C  it to the higher order formula of the pair. ERREST(*) has been used
//C  as working storage above.  The higher order approximation has been
//C  formed as YNEW(*) = Y(*) + HTRY*ERREST(*) where ERREST(*) is a linear
//C  combination of the stages of the formula. The lower order result also
//C  has the form Y(*) plus HTRY times a different linear combination of
//C  the stages. Hence, this different linear combination of stages for
//C  the lower order formula can just be subtracted from the combination
//C  stored in ERREST(*) to produce the errors. The result is then
//C  multiplied by HTRY to obtain the error estimate.
//C
   for (i = 1; i <= rkcom4.nstage; i++) {
      if (i == 1 && rkcom4.b[1] != zero) {
         for (l = 0; l < neq; l++) {
            errest[l] = errest[l] - rkcom4.b[1]*yp[l];
         }
      }
      else {
         if (rkcom4.b[i] != zero) {
            for (l = 0; l < neq; l++) {
               errest[l] = errest[l] - rkcom4.b[i]* /* stages[l][rkcom4.ptr[i]] */
               stages[l+neq*(rkcom4.ptr[i]-1)];
            }
         }
      }
   }
   for (l = 0; l < neq; l++) {
      errest[l] = htry*errest[l];
   }
//C
//C  The error in a solution component is measured relative to a weight
//C  that is the larger of a threshold and the size of the solution over
//C  the step.  Using the magnitude of a solution component at both ends
//C  of the step in the definition of "size" increases the robustness of
//C  the test. When global error estimation is specified, the weight
//C  vector WEIGHT(*) is defined by the primary integration and is then
//C  used in the secondary integration.
//C
   if (main) {
      for (l = 0; l < neq; l++) {
         avgy = half*(fabs(y[l])+fabs(ynew[l]));
         weight[l] = max(avgy,thres[l]);
      }
   }
//C
   err = zero;
   for (l = 0; l < neq; l++) {
      err = max(err,fabs(errest[l]/weight[l]));
   }
//C
   return;
}

void RKSUITE::stepa(double tnow, double y[], double yp[], double tstg, double ystg[],
   double ypstg[], double& htry, double weight[], bool& cutbak)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      To calculate an "on-scale" step size for phase 2 of
//C                the initial step size computation.
//C
//C  Input:        TNOW, Y(*), YP(*), TSTG, YSTG(*), YPSTG(*)
//C  Input/output: HTRY, WEIGHT
//C  Output:       CUTBAK
//C
//C  Common:       Initializes:    none
//C                Reads:          /RKCOM1/ TND, NEQ
//C                                /RKCOM5/ STBRAD, RS1, RS4
//C                                /RKCOM7/ RNDOFF
//C                Alters:         none
//C
//C  Comments:
//C  =========
//C  This subroutine is used during the first three stages of the first step.
//C  A Lipschitz constant L for the differential equation in autonomous form
//C  is approximated, and the product abs(HTRY)*L is compared to an approximate
//C  radius, STBRAD, of the stability region of the method. The step size is 
//C  reduced as necessary, within a range specified by the step size control 
//C  parameters RS1 and RS4, to assure stability and give some confidence in 
//C  the error estimator.  If HTRY is reduced, CUTBAK is set .TRUE..
//C
//C  Y(*) and YP(*) contain the solution and its derivative at TNOW and
//C  similarly YSTG(*) and YPSTG(*) contain approximations at TSTG.
//C
//C  Normally the weights used in the control of the error depend on the
//C  size of the solution at the beginning and at the end of the step, but
//C  at this time we do not have a solution at the end of the step.  Each
//C  stage YSTG(*) of the Runge - Kutta process represents a low order
//C  approximation to the solution at TSTG.  Because the initial value of
//C  WEIGHT(*) provided in the first phase of the scheme is based only on
//C  the solution at T and THRES(*), it is continually updated in STEPA to
//C  account for the size of the solution throughout the step as revealed
//C  by the intermediate stages YSTG(*). Inside this subroutine only, the
//C  differential equation is converted to autonomous form. After the
//C  conversion, the end of the interval of integration, TND, is used
//C  to define a suitable weight for the independent variable.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  HTRY, TNOW, TSTG
//      LOGICAL           CUTBAK
//C     .. Array Arguments ..
//      double PRECISION  WEIGHT(*), Y(*), YP(*), YPSTG(*), YSTG(*)
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Parameters ..
//      double PRECISION  ZERO
//      PARAMETER         (ZERO=0.0D0)
   const double zero = 0.0;
//C     .. Local Scalars ..
   double argdif, fdiff, scl, tdiff, twt, wt, ynrm, ystgnm;
   int l;
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX, MIN
//C     .. Executable Statements ..
//C
//C  Update the weights to account for the current intermediate solution
//C  approximation YSTG(*).  Compute the sizes of Y(*) and YSTG(*) in the
//C  new norm.  The size of the Lipschitz constant is assessed by a difference
//C  in the arguments Y(*), YSTG(*) and a difference in the function evaluated
//C  at these arguments.
//C
   ynrm = zero;
   ystgnm = zero;
   argdif = zero;
   fdiff = zero;
   for (l = 0; l < rkcom1.neqn; l++) {
      wt = max(weight[l],fabs(ystg[l]));
      weight[l] = wt;
      ynrm = max(ynrm,fabs(y[l])/wt);
      ystgnm = max(ystgnm,fabs(ystg[l])/wt);
      argdif = max(argdif,fabs(ystg[l]-y[l])/wt);
      fdiff = max(fdiff,fabs(ypstg[l]-yp[l])/wt);
   }
//C
//C  The transformation of the equation to autonomous form is done
//C  implicitly.  The difference of the arguments must take into account
//C  the difference between the values of the independent variable T and
//C  TSTG. The difference of the corresponding component of the function
//C  is zero because of the way the standard transformation is done.
//C
   tdiff = tstg - tnow;
   twt = fabs(rkcom1.tnd-tnow);
   ynrm = max(ynrm,fabs(tnow)/twt);
   ystgnm = max(ystgnm,fabs(tstg)/twt);
   argdif = max(argdif,fabs(tdiff)/twt);
//C
//C  The ratio FDIFF/ARGDIF is a lower bound for, and an approximation to, a
//C  Lipschitz constant L for the differential equation written in autonomous
//C  form.  First we must ask if the difference ARGDIF is significant in the 
//C  precision available.  If it appears to be, we insist that abs(HTRY)*L be 
//C  less than an approximate radius, STBRAD, of the stability region of the
//C  method.  This is more stringent than necessary for stability, possibly a
//C  lot more stringent, but the aim is to get an HTRY small enough that the
//C  error estimate for the step is credible.  The reduction is required to be
//C  at least as much as the step control parameter RS1. It is necessary to 
//C  limit the reduction of HTRY at any one time because we may be misled in 
//C  the size of the reduction that is appropriate due to nonlinearity of the 
//C  differential equation and to inaccurate weights caused by HTRY much too 
//C  large.  The reduction is not permitted to be more than the step control 
//C  parameter RS4.
//C
   cutbak = false;
   if (argdif > rkcom7.rndoff*max(ynrm,ystgnm)) {
      if ((fabs(htry)*fdiff) > (rkcom5.stbrad*argdif)) {
         scl = (rkcom5.stbrad*argdif)/(fabs(htry)*fdiff);
         scl = min(scl,rkcom5.rs1);
         scl = max(scl,rkcom5.rs4);
         htry = scl*htry;
         cutbak = true;
      }
   }
//C
   return;
}

void RKSUITE::stepb(int neq, double y[], double yp[], double h, double ynew[],
   double stages/*[neq][]*/[], double thres[], double& err, bool main,
   double weight[])
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      To compute an error estimate for METHD = 2 prior
//C                to completing the step.
//C
//C  Input:        NEQ, Y(*), YP(*), H, STAGES(NEQ,*), THRES(*), MAIN,
//C                WEIGHT(*)
//C  Output:       ERR
//C
//C  Common:       Initializes:    none
//C                Reads:          /RKCOM4/ E, PTR
//C                Alters:         none
//C
//C  Comments:
//C  =========
//C  If global error assessment is taking place, then MAIN = .FALSE. and
//C  the weight vector generated by the primary integration is used.  The
//C  error estimate is a linear combination (with coefficients in E(*))
//C  of the stages stored in STAGES(*,*) (located by PTR(*)).
//C
//C     .. Scalar Arguments ..
//      double PRECISION  ERR, H
//      INTEGER           NEQ
//      LOGICAL           MAIN
//C     .. Array Arguments ..
//      double PRECISION  STAGES(NEQ,*), THRES(*), WEIGHT(*), Y(*),
//     &                  YNEW(*), YP(*)
//C     .. Common Block to hold Formula Definitions ..
//      double PRECISION  A(13,13), B(13), C(13), BHAT(13), R(11,6),
//     &                  E(7)
//      INTEGER           PTR(13), NSTAGE, METHD, MINTP
//      LOGICAL           INTP
//      COMMON /RKCOM4/   A, B, C, BHAT, R, E, PTR, NSTAGE, METHD,
//     &                  MINTP, INTP
//      SAVE   /RKCOM4/
//C     .. Parameters ..
//      double PRECISION  ZERO, HALF
//      PARAMETER         (ZERO=0.0D0,HALF=0.5D0)
   const double zero = 0.0, half = 0.5;
//C     .. Local Scalars ..
   double avgy, sum, wt;
   int index, l;
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MAX
//C     .. Executable Statements ..
//C
   err = zero;
   for (l = 0; l < neq; l++) {
//C
//C  Estimate the local error of component L. The coding makes use of
//C  E(2) = 0.0D0 and E(7) = 0.0D0.
//C
      sum = rkcom4.e[1]*yp[l];
      for (index = 3; index <= 6; index++) {
         sum = sum + rkcom4.e[index]* /* stages[l][rkcom4.ptr[index]] */
         stages[l+neq*(rkcom4.ptr[index]-1)];
      }
//C
//C  The local error is H*SUM.  A weighted maximum norm of SUM is formed 
//C  and then the factor of H is taken into account.
//C  
      if (main) {
         avgy = half*(fabs(y[l])+fabs(ynew[l]));
         wt = max(avgy,thres[l]);
      }
      else {
         wt = weight[l];
      }
//C
      err = max(err,fabs(sum/wt));
   }
   err = fabs(h)*err;
//C
   return;
}

void RKSUITE::stiff(void (*f)(double, double*, double*), double havg, int& jflstp,
   bool toomch, int maxfcn, double work[], int& ier, int& nrec)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:      Diagnose stiffness.  This depends on two things: whether
//C                the step size is being restricted on grounds of stability
//C                and whether the integration to TND can be completed in no
//C                more than MAXFCN function evaluations.
//C
//C  Input:        HAVG, TOOMCH, MAXFCN, WORK(*)
//C  Input/output: JFLSTP
//C  Output:       IER, NREC
//C  Workspace:    WORK(*)
//C  External:     F
//C
//C  Common:       Initializes:    /RKCOM9/ REC
//C                Reads:          /RKCOM1/ TND, NEQN
//C                                /RKCOM2/ T, H, NFCN, SVNFCN, OKSTP
//C                                /RKCOM3/ PRY, PRYP, PRTHRS, PRWT, PRSCR,
//C                                         PRSTGS, PRYOLD
//C                                /RKCOM5/ COST
//C                Alters:         /RKCOM2/ NFCN
//C                                /RKCOM9/ REC
//C
//C     .. Scalar Arguments ..
//      double PRECISION  HAVG
//      INTEGER           IER, JFLSTP, MAXFCN, NREC
//      LOGICAL           TOOMCH
//C     .. Array Arguments ..
//      double PRECISION  WORK(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Common Block for General Workspace Pointers ..
//      INTEGER           PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      COMMON /RKCOM3/   PRTHRS, PRERST, PRWT, PRYOLD, PRSCR, PRY, PRYP,
//     &                  PRSTGS, PRINTP, LNINTP
//      SAVE   /RKCOM3/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Error Message ..
//      CHARACTER*80      REC(10)
//      COMMON /RKCOM9/   REC
//      SAVE   /RKCOM9/
//C     .. Parameters ..
//      double PRECISION  HALF
//      PARAMETER         (HALF=0.5D0)
   const double half = 0.5;
//C     .. Local Scalars ..
   double avgy, xtrawk;
   int l;
   bool lotsfl, stif, unsure;
//C     .. External Subroutines ..
//    EXTERNAL          STIFFA
//C     .. Intrinsic Functions ..
//    INTRINSIC         ABS, DBLE, MAX, MOD
//C     .. Executable Statements ..
//C
   if ((rkcom2.okstp-10) % 40 == 0) {
      lotsfl = (jflstp >> 10) != 0;
      jflstp = 0;
   }
   else {
      lotsfl = false;
   }
//C
//C  If either too much work has been done or there are lots of failed steps,
//C  test for stiffness.
//C
   if (toomch || lotsfl) {
//C
//C  Regenerate weight vector
      for (l = 0; l < rkcom1.neqn; l++) {
         avgy = half*(fabs(work[rkcom3.pry+l])+fabs(work[rkcom3.pryold+l]));
         work[rkcom3.prwt+l] = max(avgy,work[rkcom3.prthrs+l]);
      }
//C
//C  STIFFA determines whether the problem is STIFF. In some circumstances it
//C  is UNSURE.  The decision depends on two things: whether the step size is
//C  being restricted on grounds of stability and whether the integration to
//C  TND can be completed in no more than MAXFCN function evaluations.  The
//C  last four arguments of STIFFA are vectors of length NEQN used for working
//C  storage.  Some storage in WORK(*) reserved for the stages (there are a
//C  minimum of three such vectors reserved for the METHDs implemented) and
//C  the scratch vector starting at PRSCR are used for this purpose.
//C
      stiffa(f,rkcom2.t,&work[rkcom3.pry],rkcom2.h,havg,rkcom1.tnd,maxfcn,&work[rkcom3.prwt],
         &work[rkcom3.pryp],&work[rkcom3.prerst],unsure,stif,&work[rkcom3.prstgs],
         &work[rkcom3.prstgs+rkcom1.neqn],&work[rkcom3.prstgs+2*rkcom1.neqn],
         &work[rkcom3.prscr]);
      if (!unsure) {
         if (stif) {
//C
//C  Predict how much eXTRA WorK will be needed to reach TND.
            xtrawk = (rkcom5.cost*fabs((rkcom1.tnd-rkcom2.t)/havg))/double(rkcom2.svnfcn+rkcom2.nfcn);
            ier = 4;
            sprintf(&rkcom9.rec[nrec][0]," ** Your problem has been diagnosed as stiff.  If the");
            sprintf(&rkcom9.rec[nrec+1][0]," ** situation persists, it will cost roughly %13.5lf", xtrawk);
            sprintf(&rkcom9.rec[nrec+2][0]," ** times as much to reach TEND as it has cost to reach TNOW.");
            sprintf(&rkcom9.rec[nrec+3][0]," ** You should probably change to a code intended for");
            sprintf(&rkcom9.rec[nrec+4][0]," ** stiff problems.");
            nrec += 5;
         }
      }
   }
//C
   return;
}

void RKSUITE::stiffa(void (*f)(double, double*, double*),double x, double y[],
   double hnow, double havg, double xend, int maxfcn, double wt[],
   double fxy[], double v0[], bool& unsure, bool& stif, double v1[],
   double v2[], double v3[], double vtemp[])
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  External:     F
//C  Input:        X, Y(*), HNOW, HAVG, XEND, MAXFCN, WT(*), FXY(*)
//C  Input/Output  V0(*)
//C  Output:       UNSURE, STIF
//C  Workspace:    V1(*), V2(*), V3(*), VTEMP(*)
//C
//C  Common:       Initializes:    none
//C                Reads:          /RKCOM1/ TND, NEQN
//C                                /RKCOM5/ COST, STBRAD, TANANG
//C                                /RKCOM7/ SQRRMC, CUBRMC
//C                Alters:         none
//C
//C  STIFFA diagnoses stiffness for an explicit Runge-Kutta code.  When it
//C  is called, either many step failures have been observed, or a lot of
//C  work has been done.
//C
//C  The NEQ equations of the problem are defined by the subroutine F(X,Y,YP).
//C  When STIFFA is called, the integration has reached X where the approximate
//C  solution is Y(*).  The vector FXY(*) is defined by a call of F(X,Y,FXY).
//C  It is an input argument because it is usually available from the integrator.
//C
//C  The last successful step was of size HNOW, and an average step size is
//C  HAVG.  A weighted norm is used to measure the local error with the error
//C  in solution component L divided by the positive weight WT(L) provided in 
//C  the vector WT(*).
//C
//C  Explicit Runge - Kutta codes estimate the local error of Y(*) by
//C  forming the difference of two approximate solutions.  This difference
//C  must be provided in the vector V0(*).  When this difference is too
//C  small to be significant, STIFFA will replace it with a "random" vector.
//C
//C  STIF is set .TRUE. when the average step size appears to be restricted
//C  on grounds of stability.  In certain cases the variable UNSURE is set 
//C  .TRUE.; the value of STIF is then not defined.
//C
//C  The stability region of the explicit Runge-Kutta formula is described
//C  by quantities TANANG and STBRAD that are communicated by the setup routine
//C  via COMMON.  Stability regions often change sharply near the imaginary
//C  axis so that it is difficult to classify the stiffness of a problem with
//C  eigenvalues of a local Jacobian that are "near" the imaginary axis.  For
//C  this reason,  we consider only points Z in the upper left half complex
//C  plane for which TAN( IMAG(Z)/( - RE(Z))) <= TANANG. Eigenvalues outside
//C  this region are one reason for the code being UNSURE.  The stability
//C  region is approximated by the intersection of a disk with this sector.
//C  The radius of this disk is called STBRAD.
//C
//C  Working storage must be provided via the four vectors V1(*),V2(*),
//C  V3(*),VTEMP(*).  These vectors must be of length at least NEQ.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  HAVG, HNOW, X, XEND
//      INTEGER           MAXFCN
//      LOGICAL           STIF, UNSURE
//C     .. Array Arguments ..
//      double PRECISION  FXY(*), V0(*), V1(*), V2(*), V3(*), VTEMP(*),
//     &                  WT(*), Y(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Formula Characterisitcs ..
//      double PRECISION  TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4
//      INTEGER           ORDER, LSTSTG, MAXTRY, NSEC
//      LOGICAL           FSAL
//      COMMON /RKCOM5/   TOOSML, COST, SAFETY, EXPON, STBRAD, TANANG,
//     &                  RS, RS1, RS2, RS3, RS4, ORDER, LSTSTG, MAXTRY,
//     &                  NSEC, FSAL
//      SAVE   /RKCOM5/
//C     .. Common Block for Environment Parameters ..
//      double PRECISION  MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY
//      INTEGER           OUTCH
//      COMMON /RKCOM7/   MCHEPS, DWARF, RNDOFF, SQRRMC, CUBRMC, TINY,
//     &                  OUTCH
//      SAVE   /RKCOM7/
//C     .. Parameters ..
//      double PRECISION  LARGE
//      PARAMETER         (LARGE=1.0D+10)
//      double PRECISION  ZERO, P001, P9, ONE, TWO, FIVE, FIFTH
//      PARAMETER         (ZERO=0.0D+0,P001=0.001D+0,P9=0.9D+0,ONE=1.0D+0,
//     &                  TWO=2.0D+0,FIVE=5.0D+0,FIFTH=0.2D+0)
   const double large = 1.0e10;
   const double zero = 0.0, p001 = 0.001, p9 = 0.9, one = 1.0;
   const double two = 2.0, five = 5.0, fifth = 0.2;
//C     .. Local Scalars ..
   double alpha1, alpha2, beta1, beta2, d1, d2, det1,
      det2, dist, res2, rho, rho2, rold, scale, v0nrm,
      v0v0, v0v1, v0v2, v1v1, v1v2, v1v3, v2v2, v2v3,
      v3nrm, v3v3, xtrfcn, ynrm;
   int l, ntry;
   bool rootre;
//C     .. Local Arrays ..
   double r1[2], r2[2], root1[2], root2[2];
//C     .. External Functions ..
//      double PRECISION  DOTPRD
//      EXTERNAL          DOTPRD
//C     .. External Subroutines ..
//      EXTERNAL          STIFFB, STIFFC, STIFFD
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, MIN, SQRT
//C     .. Executable Statements ..
//C
//C  If the current step size differs substantially from the average,
//C  the problem is not stiff.
//C
   if (fabs(hnow/havg) > five || fabs(hnow/havg) < fifth) {
      stif = false;
      unsure = false;
      return;
   }
   else {
      unsure = true;
   }
//C
//C  The average step size is used to predict the cost in function evaluations
//C  of finishing the integration to XEND.  If this cost is no more than MAXFCN,
//C  the problem is declared not stiff: If the step size is being restricted on
//C  grounds of stability, it will stay close to HAVG.  The prediction will
//C  then be good, but the cost is too low to consider the problem stiff.  If
//C  the step size is not close to HAVG, the problem is not stiff.  Either way
//C  there is no point to testing for a step size restriction due to stability.
//C
   xtrfcn = rkcom5.cost*fabs((xend-x)/havg);
   if (xtrfcn <= maxfcn) {
      stif = false;
      unsure = false;
      return;
   }
   else {
      unsure = true;
   }
//C
//C  There have been many step failures or a lot of work has been done.  Now 
//C  we must determine if this is due to the stability characteristics of the
//C  formula.  This is done by calculating the dominant eigenvalues of the
//C  local Jacobian and then testing whether HAVG corresponds to being on the
//C  boundary of the stability region.
//C
//C  The size of Y(*) provides scale information needed to approximate
//C  the Jacobian by differences.
//C
   ynrm = sqrt(dotprd(y,y,wt,rkcom1.neqn));
   scale = ynrm*rkcom7.sqrrmc;
   if (scale == zero) {
//C
//C  Degenerate case.  Y(*) is (almost) the zero vector so the scale is not 
//C  defined.  The input vector V0(*) is the difference between Y(*) and a 
//C  lower order approximation to the solution that is within the error 
//C  tolerance.  When Y(*) vanishes, V0(*) is itself an acceptable approximate
//C  solution, so we take SCALE from it, if this is possible.
//C
      ynrm = sqrt(dotprd(v0,v0,wt,rkcom1.neqn));
      scale = ynrm*rkcom7.sqrrmc;
      if (scale == zero) {
         unsure = true;
         return;
      }
   }
//C
   v0v0 = dotprd(v0,v0,wt,rkcom1.neqn);
   if (v0v0 == zero) {
//C
//C  Degenerate case.  V0(*) is (almost) the zero vector so cannot
//C  be used to define a direction for an increment to Y(*).  Try a
//C  "random" direction.
//C
      for (l = 0; l < rkcom1.neqn; l++) {
         v0[l] = one;
      }
      v0v0 = dotprd(v0,v0,wt,rkcom1.neqn);
   }
   v0nrm = sqrt(v0v0);
   for (l = 0; l < rkcom1.neqn; l++) {
      v0[l] = v0[l]/v0nrm;
   }
   v0v0 = one;
//C
//C  Use a nonlinear power method to estimate the two dominant eigenvalues.
//C  V0(*) is often very rich in the two associated eigenvectors.  For this 
//C  reason the computation is organized with the expectation that a minimal 
//C  number of iterations will suffice.  Indeed, it is necessary to recognize 
//C  a kind of degeneracy when there is a dominant real eigenvalue.  The
//C  subroutine STIFFB does this.  In the first try, NTRY = 1, a Rayleigh 
//C  quotient for such an eigenvalue is initialized as ROLD.  After each 
//C  iteration, REROOT computes a new Rayleigh quotient and tests whether the
//C  two approximations agree to one tenth of one per cent and the eigenvalue,
//C  eigenvector pair satisfy a stringent test on the residual.  ROOTRE = .TRUE.
//C  signals that a single dominant real root has been found.
//C
   ntry = 1;
   label60:
//C
   stiffd(v0,havg,x,y,f,fxy,wt,scale,v0v0,v1,v1v1,vtemp);
//C
//C  The quantity SQRT(V1V1/V0V0) is a lower bound for the product of HAVG
//C  and a Lipschitz constant.  If it should be LARGE, stiffness is not
//C  restricting the step size to the stability region.  The principle is
//C  clear enough, but the real reason for this test is to recognize an
//C  extremely inaccurate computation of V1V1 due to finite precision
//C  arithmetic in certain degenerate circumstances.
//C
   if (sqrt(v1v1) > large*sqrt(v0v0)) {
      unsure = true;
      return;
   }
//C
   v0v1 = dotprd(v0,v1,wt,rkcom1.neqn);
   if (ntry == 1) {
      rold = v0v1/v0v0;
//C
//C  This is the first Rayleigh quotient approximating the product of HAVG
//C  and a dominant real eigenvalue.  If it should be very small, the
//C  problem is not stiff.  It is important to test for this possibility so
//C  as to prevent underflow and degeneracies in the subsequent iteration.
//C
      if (fabs(rold) < rkcom7.cubrmc) {
         unsure = false;
         stif = false;
         return;
      }
   }
   else {
      stiffb(v1v1,v0v1,v0v0,rold,rho,root1,root2,rootre);
      if (rootre) goto label100;
   }
   stiffd(v1,havg,x,y,f,fxy,wt,scale,v1v1,v2,v2v2,vtemp);
   v0v2 = dotprd(v0,v2,wt,rkcom1.neqn);
   v1v2 = dotprd(v1,v2,wt,rkcom1.neqn);
   stiffb(v2v2,v1v2,v1v1,rold,rho,root1,root2,rootre);
   if (rootre) goto label100;
//C
//C  Fit a quadratic in the eigenvalue to the three successive iterates
//C  V0(*),V1(*),V2(*) of the power method to get a first approximation to
//C  a pair of eigenvalues.  A test made earlier in STIFFB implies that
//C  the quantity DET1 here will not be too small.
//C
   det1 = v0v0*v1v1 - v0v1*v0v1;
   alpha1 = (-v0v0*v1v2+v0v1*v0v2)/det1;
   beta1 = (v0v1*v1v2-v1v1*v0v2)/det1;
//C
//C  Iterate again to get V3, test again for degeneracy, and then fit a
//C  quadratic to V1(*),V2(*),V3(*) to get a second approximation to a pair
//C  of eigenvalues.
//C
   stiffd(v2,havg,x,y,f,fxy,wt,scale,v2v2,v3,v3v3,vtemp);
   v1v3 = dotprd(v1,v3,wt,rkcom1.neqn);
   v2v3 = dotprd(v2,v3,wt,rkcom1.neqn);
   stiffb(v3v3,v2v3,v2v2,rold,rho,root1,root2,rootre);
   if (rootre) goto label100;
   det2 = v1v1*v2v2 - v1v2*v1v2;
   alpha2 = (-v1v1*v2v3+v1v2*v1v3)/det2;
   beta2 = (v1v2*v2v3-v2v2*v1v3)/det2;
//C
//C  First test the residual of the quadratic fit to see if we might
//C  have determined a pair of eigenvalues.
//C
   res2 = fabs(v3v3+v2v2*alpha2*alpha2+v1v1*beta2*beta2+two*v2v3*alpha2+
      two*v1v3*beta2+two*v1v2*alpha2*beta2);
   if (res2 <= v3v3*p001*p001) {
//C
//C  Calculate the two approximate pairs of eigenvalues.
//C
      stiffc(alpha1,beta1,r1,r2);
      stiffc(alpha2,beta2,root1,root2);
//C
//C  The test for convergence is done on the larger root of the second
//C  approximation.  It is complicated by the fact that one pair of roots 
//C  might be real and the other complex.  First calculate the spectral 
//C  radius RHO of HAVG*J as the magnitude of ROOT1.  Then see if one of 
//C  the roots R1,R2 is within one per cent of ROOT1.  A subdominant root 
//C  may be very poorly approximated if its magnitude is much smaller than 
//C  RHO -- this does not matter in our use of these eigenvalues.
//C
      rho = sqrt(root1[0]*root1[0]+root1[1]*root1[1]);
      d1 = (root1[0]-r1[0])*(root1[0]-r1[0]) +
         (root1[1]-r1[1])*(root1[1]-r1[1]);
      d2 = (root1[0]-r2[0])*(root1[0]-r2[0]) +
         (root1[1]-r2[1])*(root1[1]-r2[1]);
      dist = sqrt(min(d1,d2));
      if (dist <= p001*rho) goto label100;
   }
//C
//C  Do not have convergence yet.  Because the iterations are cheap, and
//C  because the convergence criterion is stringent, we are willing to try
//C  a few iterations.
//C
   if (ntry < rkcom5.maxtry) {
      ntry++;
      v3nrm = sqrt(v3v3);
      for (l = 0; l < rkcom1.neqn; l++) {
         v0[l] = v3[l]/v3nrm;
      }
      v0v0 = one;
      goto label60;
   }
   else {
      unsure = true;
      return;
   }
//C
//C                        **************
//C
//C  We now have the dominant eigenvalues.  Decide if the average step
//C  size is being restricted on grounds of stability.  Check the real
//C  parts of the eigenvalues.  First see if the dominant eigenvalue is
//C  in the left half plane -- there won't be a stability restriction
//C  unless it is. If there is another eigenvalue of comparable magnitude
//C  with a positive real part, the problem is not stiff. If the dominant
//C  eigenvalue is too close to the imaginary axis, we cannot diagnose
//C  stiffness.
//C
   label100:
   if (root1[0] > zero) {
      stif = false;
      unsure = false;
      return;
   }
   rho2 = sqrt(root2[0]*root2[0]+root2[1]*root2[1]);
   if (rho2 >= p9*rho && root2[0] > zero) {
      stif = false;
      unsure = false;
      return;
   }
   if (fabs(root1[1]) > fabs(root1[0])*rkcom5.tanang) {
      unsure = true;
      return;
   }
//C
//C  If the average step size corresponds to being well within the
//C  stability region, the step size is not being restricted because
//C  of stability.
//C
   stif = rho >= p9*rkcom5.stbrad;
   unsure = false;
   return;
}

void RKSUITE::stiffb(double v1v1, double v0v1, double v0v0, double& rold,
   double& rho, double root1[], double root2[], bool& rootre)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Input:        V1V1, V0V1, V0V0
//C  Input/output: ROLD
//C  Output:       RHO, ROOT1(*),ROOT2(*),ROOTRE
//C
//C  Decide if the iteration has degenerated because of a strongly
//C  dominant real eigenvalue.  Have just computed the latest iterate.
//C  V1V1 is its dot product with itself, V0V1 is the dot product
//C  of the previous iterate with the current one, and V0V0 is the
//C  dot product of the previous iterate with itself.  ROLD is a
//C  previous Rayleigh quotient approximating a dominant real
//C  eigenvalue.  It must be computed directly the first time the
//C  subroutine is called.  It is updated each call to STIFFB, hence
//C  is available for subsequent calls.
//C
//C  If there is a strongly dominant real eigenvalue, ROOTRE is set
//C  .TRUE., ROOT1(*) returns the eigenvalue, RHO returns the magnitude
//C  of the eigenvalue, and ROOT2(*) is set to zero.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  RHO, ROLD, V0V0, V0V1, V1V1
//      LOGICAL           ROOTRE
//C     .. Array Arguments ..
//      double PRECISION  ROOT1(2), ROOT2(2)
//C     .. Parameters ..
//      double PRECISION  ZERO, P001
//      PARAMETER         (ZERO=0.0D+0,P001=0.001D+0)
   const double zero = 0.0, p001 = 0.001;
//C     .. Local Scalars ..
   double det, r, res;
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS
//C     .. Executable Statements ..
//C
   r = v0v1/v0v0;
   rho = fabs(r);
   det = v0v0*v1v1 - v0v1*v0v1;
   res = fabs(det/v0v0);
   rootre = det == zero || (res <= v1v1*p001*p001 &&
      fabs(r-rold) <= p001*rho);
   if (rootre) {
      root1[0] = r;
      root1[1] = zero;
      root2[0] = zero;
      root2[1] = zero;
   }
   rold = r;
//C
   return;
}

void RKSUITE::stiffc(double alpha, double beta, double r1[], double r2[])
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Input:  ALPHA, BETA
//C  Output: R1(*), R2(*)
//C
//C  This subroutine computes the two complex roots R1 and R2 of
//C  the quadratic equation X**2 + ALPHA*X + BETA = 0.  The magnitude
//C  of R1 is greater than or equal to the magnitude of R2. R1 and R2 are
//C  returned as vectors of two components with the first being the real
//C  part of the complex number and the second being the imaginary part.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  ALPHA, BETA
//C     .. Array Arguments ..
//      double PRECISION  R1(2), R2(2)
//C     .. Parameters ..
//      double PRECISION  ZERO, TWO
//      PARAMETER         (ZERO=0.0D+0,TWO=2.0D+0)
   const double zero = 0.0, two = 2.0;
//C     .. Local Scalars ..
   double disc, sqdisc, temp;
//C     .. Intrinsic Functions ..
//      INTRINSIC         ABS, SQRT
//C     .. Executable Statements ..
   temp = alpha/two;
   disc = temp*temp - beta;
   if (disc == zero) {
//C
//C  double root.
//C
      r1[0] = -temp;
      r1[1] = zero;
      r2[0] = r1[0];
      r2[1] = r1[1];
      return;
   }
//C
   sqdisc = sqrt(fabs(disc));
   if (disc < zero) {
//C
//C  Complex conjugate roots.
//C
      r1[0] = -temp;
      r1[1] = sqdisc;
      r2[0] = r1[0];
      r2[1] = -r1[1];
   }
   else {
//C
//C  Real pair of roots.  Calculate the bigger one in R1(1).
//C
      if (temp > zero) {
         r1[0] = -temp - sqdisc;
      }
      else {
         r1[0] = -temp + sqdisc;
      }
      r1[1] = zero;
      r2[0] = beta/r1[0];
      r2[1] = zero;
   }
//C
   return;
}

void RKSUITE::stiffd(double v[], double havg, double x, double y[],
   void (*f)(double, double[], double[]), double fxy[], double wt[],
   double scale, double vdotv, double z[], double& zdotz, double vtemp[])
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  External:     F
//C  Input:        V(*), HAVG, X, Y(*), FXY(*), WT(*), SCALE, VDOTV,
//C  Output:       Z(*), ZDOTZ
//C  Workspace:    VTEMP(*)
//C
//C  For an input vector V(*) of length NEQ, this subroutine computes a vector
//C  Z(*) that approximates the product HAVG*J*V where HAVG is an input scalar 
//C  and J is the Jacobian matrix of a function F evaluated at the input 
//C  arguments (X,Y(*)).  This function is defined by a subroutine of the form
//C  F(T,U,F) that when given T and U(*), returns the value of the function in 
//C  F(*).  The input vector FXY(*) is defined by F(X,Y,FXY).  Scaling is a 
//C  delicate matter.  A weighted Euclidean norm is used with the (positive) 
//C  weights provided in WT(*).  The input scalar SCALE is the square root of 
//C  the unit roundoff times the norm of Y(*).  The square of the norm of the
//C  input vector V(*) is input as VDOTV.  The routine outputs the square of
//C  the norm of the output vector Z(*) as ZDOTZ.  The subroutine calls the
//C  double PRECISION FUNCTION DOTPRD(U,V,WT,NEQ) to compute the dot (inner)
//C  product.  The vector VTEMP(*) is used for working storage.
//C
//C     .. Scalar Arguments ..
//      double PRECISION  HAVG, SCALE, VDOTV, X, ZDOTZ
//C     .. Array Arguments ..
//      double PRECISION  FXY(*), V(*), VTEMP(*), WT(*), Y(*), Z(*)
//C     .. Subroutine Arguments ..
//      EXTERNAL          F
//C     .. Common Block for Problem Definition ..
//      double PRECISION  TSTRT, TND, DIR, HSTRT, TOLR
//      INTEGER           NEQN
//      COMMON /RKCOM1/   TSTRT, TND, DIR, HSTRT, TOLR, NEQN
//      SAVE   /RKCOM1/
//C     .. Common Block to hold Problem Status ..
//      double PRECISION  T, H, TOLD, HOLD
//      INTEGER           NFCN, SVNFCN, OKSTP, FLSTP
//      LOGICAL           FIRST, LAST
//      COMMON /RKCOM2/   T, H, TOLD, HOLD, NFCN, SVNFCN, OKSTP, FLSTP,
//     &                  FIRST, LAST
//      SAVE   /RKCOM2/
//C     .. Local Scalars ..
   double temp1, temp2;
   int l;
//C     .. External Functions ..
//      double PRECISION  DOTPRD
//      EXTERNAL          DOTPRD
//C     .. Intrinsic Functions ..
//      INTRINSIC         SQRT
//C     .. Executable Statements ..
//C
//C  Scale V(*) so that it can be used as an increment to Y(*)
//C  for an accurate difference approximation to the Jacobian.
//C
   temp1 = scale/sqrt(vdotv);
   for (l = 0; l < rkcom1.neqn; l++) {
      vtemp[l] = y[l] + temp1*v[l];
   }
//C
   f(x,vtemp,z);
   rkcom2.nfcn++;
//C
//C  Form the difference approximation.  At the same time undo
//C  the scaling of V(*) and introduce the factor of HAVG.
//C
   temp2 = havg/temp1;
   for (l = 0; l < rkcom1.neqn; l++) {
      z[l] = temp2*(z[l]-fxy[l]);
   }
//C
   zdotz = dotprd(z,z,wt,rkcom1.neqn);
//C
   return;
}

double RKSUITE::dotprd(double u[], double v[], double wt[], int neq)
{
//C************************************************
//C**** NOT A DESIGNATED USER-CALLABLE ROUTINE ****
//C************************************************
//C
//C  Purpose:   To compute a weighted Euclidean dot (inner) product of
//C             two vectors.
//C
//C  Input:     U(*), V(*), WT(*), NEQ
//C  Output:    the result DOTPRD is returned via the subprogram name
//C
//C  Comments:
//C  =========
//C  The vectors U(*), V(*), and WT(*) are of length NEQ. The components
//C  of WT(*) are weights that must be non-zero.
//C
//C     .. Scalar Arguments ..
//      INTEGER           NEQ
//C     .. Array Arguments ..
//      double PRECISION  U(*), V(*), WT(*)
//C     .. Parameters ..
//      double PRECISION  ZERO
//      PARAMETER         (ZERO=0.0D0)
   const double zero = 0.0;
//C     .. Local Scalars ..
   double sum;
   int l;
//C     .. Executable Statements ..
//C
   sum = zero;
   for (l = 0; l < neq; l++) {
      sum += (u[l]/wt[l])*(v[l]/wt[l]);
   }
//C
//   DOTPRD = SUM
//C
   return sum;
}

void RKSUITE::softfl(bool ask, bool& on)
{
//C
//C  Purpose:      To prevent a program STOP after a "catastrophic"
//C                failure when using a routine from RKSUITE.
//C
//C  Input:        ASK
//C  Input/output: ON
//C
//C  Comments:
//C  =========
//C  When a "catastrophic" failure is detected, the default action of
//C  RKSUITE is to write an explanation to the standard output channel,
//C  OUTCH, and STOP.  This subroutine can be used to prevent the STOP and
//C  so allow the main program to continue.  To do this, you call SOFTFL with
//C  ASK = .FALSE. and ON = .TRUE.  You must then call the subroutine CHKFL
//C  after every call to a user-callable routine in RKSUITE to check whether
//C  a catastrophic error occurred and take appropriate action if it did.  Of
//C  course, you may call SETUP at any time to start a new problem, but calling
//C  any other user-callable routine in RKSUITE after a catastrophic error will
//C  lead to a STOP (even when "soft failure" has been set "on").
//C
//C  When ON is set by a call to SOFTFL with ASK = .FALSE., the value of ON
//C  is SAVEd.  The subroutine RKMSG in RKSUITE calls SOFTFL with ASK = .TRUE.
//C  to find out the SAVEd value of ON.
//C
//C     .. Scalar Arguments ..
//      LOGICAL           ASK, ON
//C     .. Local Scalars ..
    static bool soft = false;
//C     .. Save statement ..
//      SAVE              SOFT
//C     .. Data statements ..
//      DATA              SOFT/.FALSE./
//C     .. Executable Statements ..
//C
   if (ask) {
      on = soft;
   }
   else {
      soft = on;
   }
//C
   return;
}

void RKSUITE::chkfl(bool ask, bool& error)
{
//C
//C  Purpose:      Enquiry routine used in conjunction with SOFTFL.
//C                Reports whether a "catastrophic" error was detected.
//C
//C  Input:        ASK
//C  Input/output: ERROR
//C
//C  Comments:
//C  =========
//C  When a "catastrophic" failure is detected, the default action of
//C  RKSUITE is to write an explanation to the standard output channel,
//C  OUTCH, and STOP.  SOFTFL can be used to prevent the STOP and so
//C  allow the main program to continue.  It is then necessary to call
//C  CHKFL with ASK = .TRUE. after every call to a user-callable routine 
//C  in RKSUITE to check whether a catastrophic error occurred and take 
//C  appropriate action if it did.  If there was a catastrophic error, 
//C  ERROR is returned .TRUE.  Of course, you may call SETUP at any time 
//C  to start a new problem, but calling any other user-callable routine 
//C  in RKSUITE after a catastrophic error will lead to a STOP (even when
//C  "soft failure" has been set "on").
//C
//C  When a catastrophic failure (IER = 911) is detected in one of
//C  the routines in RKSUITE, it calls CHKFL with ASK = .FALSE. and
//C  ERROR = .TRUE.  This value of ERROR is SAVEd.
//C
//C     .. Scalar Arguments ..
//      LOGICAL           ASK, ERROR
//C     .. Local Scalars ..
   static bool saverr = false;
//C     .. Save statement ..
//      SAVE              SAVERR
//C     .. Data statements ..
//      DATA              SAVERR/.FALSE./
//C     .. Executable Statements ..
//C
   if (ask) {
      error = saverr;
   }
   else {
      saverr = error;
   }
//C
   return;
}

