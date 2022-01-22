/*Translated by FOR_C, v3.4.2 (-), on 07/09/115 at 08:31:15 */
/*FOR_C Options SET: ftn=u io=c no=p op=aimnv s=dbov str=l x=f - prototypes */
#include <math.h>
#if
#include <omp.h>
#include "fcrt.h"
#include "diva.h"
#include <float.h>
#include <stdlib.h>
#include <string.h>
		/* PARAMETER translations */
#define	CM1	(-1.e0)
#define	KDIM	20
#define	LTXTAB	7
#define	LTXTAC	71
#define	LTXTAD	183
#define	LTXTAE	226
#define	LTXTAF	290
#define	LTXTAG	400
#define	LTXTAH	427
#define	LTXTAI	455
#define	LTXTAJ	498
#define	LTXTAK	535
#define	LTXTAL	573
#define	LTXTAM	620
#define	LTXTAN	655
#define	MAXORD	2
#define	MAXSTF	1
#define	MECONT	50
#define	MEEMES	52
#define	MEIDAT	24
#define	MEIVEC	57
#define	MEMDA1	27
#define	MENTXT	23
#define	MERET	51
#define	METEXT	53
		/* end of PARAMETER translations */
 
	/* COMMON translations */
struct t_divaev {
	double eeps2, eept75, eovep2, ovtm75, ovd10, eeps10, eeps16, erov10;
	}divaev;
#pragma omp threadprivate(divaev)

struct t_divasc {
	double tn, xi[KDIM];
	long int iopst, kordi, kqmaxd, kqmaxi, ldt, maxdif, maxint, nkdko,
	 nte, nyny, ndtf, numdt;
	}divasc;
#pragma omp threadprivate(divasc)

struct t_divamc {
	double tg[2], tgstop[2], tmark, tmarkx, tout, tolg, hc, hdec,
	 hinc, hincc, hmax, hmaxp9, hmin, alpha[KDIM], beta[KDIM + 1],
	 d[MAXORD][MAXSTF + MAXORD], g[MAXORD][KDIM], v[KDIM + MAXORD],
	 ds[MAXORD][MAXSTF + MAXORD], gs[KDIM], sigma[KDIM], rbq[KDIM],
	 dnoise, eave, eimax, eimin, emax, erep, robnd, snoise, fdat[11];
	long int icf, ics, igflg, igtype[2], igstop[2], ilgrep, ings,
	 iop3, iop4, iop5, iop6, iop7, iop8, iop9, iop10, iop11, iop12,
	 iop13, iop14, iop15, iop16, iop17, iop18, iop19, iop20, iop21,
	 iop22, iop21s, itolep, iy, kemax, kis, kmark, kord1i, kord2i,
	 kpred, kqdcon, kqicon, kqmaxs, kqmxds, kqmxil, kqmxip, kqmxis,
	 ksc, ksout, ksstrt, kstep, lex, linc, lincd, lincq, lsc, maxkqd,
	 maxkqi, method, ne, neptol, ng, ngtot, noiseq, noutko, ntolf,
	 ny, idat[6];
	}divamc;
#pragma omp threadprivate(divamc)

	/* end of COMMON translations */
void /*FUNCTION*/ diva(
double tspecs[],
double y[],
double f[],
long kord[],
long neq,
void (*divaf)(double[],double[],double[],long[]),
void (*divao)(double[],double[],double[],long[]),
long idimt,
long idimy,
long idimf,
long idimk,
long iopt[])
{
	long int ihi, ilow, j, jl, jlim, k, kferr, kgo, kqq, _i, _r;
	static long int iopiva[2];
	static char mtxtaa[3][234]={"DIVA$BThe interval [1, 10**6], bounds the allowed values for NTE=$I.$EFor option $I, the interval [$I, $I], bounds the$ allowed values for the integration order which is set to $I.$EOption 16 must be used for error control.$EF($I) = ",
	 "$F, but it must be -1.0 when skipping the error check.$EFor option $I, the interval [$I, $I] bounds the allowed values for KORD($I)=$I, which is used to specify an $Boutput type for printing.$Eoutput group for printing.$Eequation gro",
	 "up for variational equations.$Eorder for a$ differential equation.$Eequation group for diagnostic print.$Eequation group for integration order control.$Eequation group for error control.$EOption 5 argument must be .le. 0 or .gt. 4.$E"};
	static char mtxtab[1][56]={"KORD values for this option starting at KORD($M) are:$E"};
	static long mloc[12]={LTXTAI,LTXTAK,LTXTAL,LTXTAM,LTXTAJ,LTXTAG,
	 LTXTAH,LTXTAB,LTXTAC,LTXTAD,LTXTAE,LTXTAN};
	static long mact[16]={MEEMES,38,24,0,MENTXT,0,METEXT,MECONT,MEMDA1,
	 0,METEXT,MEIDAT,0,MEIVEC,0,MECONT};
	static long mact1[4]={METEXT,MEIVEC,0,MERET};
	static char text1[1][11]={"IOPT()= $B"};
	static int _aini = 1;
	/* EQUIVALENCE translations */
	long   _e0[32];
	long int * __restrict const intchk = (long*)_e0;
	long int * __restrict const nxtchk = (long*)((long*)_e0 + 1);
	/* end of EQUIVALENCE translations */
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const F = &f[0] - 1;
	double * __restrict const Fdat = &divamc.fdat[0] - 1;
	long * __restrict const Idat = &divamc.idat[0] - 1;
	long * __restrict const Iopiva = &iopiva[0] - 1;
	long * __restrict const Iopt = &iopt[0] - 1;
	long * __restrict const Kord = &kord[0] - 1;
	long * __restrict const Mact = &mact[0] - 1;
	long * __restrict const Mact1 = &mact1[0] - 1;
	long * __restrict const Mloc = &mloc[0] - 1;
	double * __restrict const Tspecs = &tspecs[0] - 1;
	double * __restrict const Y = &y[0] - 1;
		/* end of OFFSET VECTORS */
	if( _aini ){ /* Do 1 TIME INITIALIZATIONS! */
		Iopiva[1] = 1111;
		_aini = 0;
	}
 
	/* Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
	 * ALL RIGHTS RESERVED.
	 * Based on Government Sponsored Research NAS7-03001.
	 *>> 2015-03-15 DIVA  Krogh  Removed extra call divabu after noise test
	 *>> 2015-03-15 DIVA  Krogh  Forced restart needs more reduction in h.
	 *>> 2010-02-20 DIVA  Krogh  Fixed calling DIVAOP with array other than F
	 *>> 2009-11-03 DIVA  Krogh  Added option 11, more variables initialized.
	 *>> 2009-10-30 DIVA  Krogh  Gave KSSTRT and ROBND initial values.
	 *>> 2009-10-30 DIVA  Krogh  Fixed reference to undefined location in F.
	 *>> 2009-10-21 DIVA  Krogh  Got rid of NaN in diag. print when LSC=3.
	 *>> 2009-10-15 DIVA  Krogh  A few changes on how noise is handled.
	 *>> 2002-11-12 DIVA  Krogh  Fixed problem integrating to final output pt
	 *>> 2002-08-29 DIVA  Krogh  Added test for invalid HMIN/HMAX.
	 *>> 2002-07-26 DIVA  Krogh  Added KOUTKO to fully support Option 10.
	 *>> 2002-05-14 DIVA  Krogh  Fix starting prob. for Option 18.
	 *>> 2002-05-13 DIVA  Krogh  Put exponent letter in  numbers missing them
	 *>> 2002-05-12 DIVA  Krogh  Added error message for bad option 5 usage.
	 *>> 2001-09-07 DIVA  Krogh  Changes to allow user tol on G-Stops.
	 *>> 2001-05-25 DIVA  Krogh  Minor change for making .f90 version.
	 *>> 2001-05-18 DIVA  Krogh  Less computing with no error test
	 *>> 2001-05-17 DIVA  Krogh  Fixed so with no error test can't start dump
	 *>> 2001-04-24 DIVA  Krogh  Inserted comments from ivacom.
	 *>> 2000-12-01 DIVA  Krogh  Removed (some of) unused C1, MAXSTF, METEXT.
	 *>> 1999-12-28 DIVA  Krogh  Saved S in DIVACR for output consistency.
	 *>> 1999-08-19 DIVA  Krogh  Removed superfluous test above label 3520.
	 *>> 1997-04-22 DIVA  Krogh  Got rid of assigned go to's. F=0 if diag.
	 *>> 1996-08-26 DIVA  Krogh  Initialize F to 0 if dumping solution.
	 *>> 1996-08-23 DIVA  Krogh  Print TN not TSPECS(1) in error messages.
	 *>> 1996-05-30 DIVA  Krogh  Changed DERIVS/OUTPUT to  DIVAF/DIVAO.
	 *>> 1996-04-27 DIVA  Krogh  Changes to use .C. and C%%.
	 *>> 1996-03-30 DIVA  Krogh  Added external statement.
	 *>> 1996-03-25 DIVA  Krogh  Introduced TEXT1 to comply with F77.
	 *>> 1996-02-27 DIVA  Krogh  Fixed so DUMP not affected by ignored eqs.
	 *>> 1995-12-18 DIVA  Krogh  Fixed so no solution dump on 0 length integ.
	 *>> 1995-11-09 DIVA  Krogh  Fixed so char. data at col. 72 is not ' '.
	 *>> 1995-06-19 DIVA  Krogh  Fixed prob. with discon. just after restart.
	 *>> 1995-05-09 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction
	 *>> 1995-04-26 DIVA  Krogh  Use KQMAXS instead of KQMAXI when LDIS>1000.
	 *>> 1995-04-26 DIVA  Krogh  Keep current KQL on discontinutiy.
	 *>> 1994-12-16 DIVA  Krogh  Fixed option 12 with K12 < 0.
	 *>> 1994-11-11 DIVA  Krogh  Declared all vars.
	 *>> 1994-11-02 DIVA  Krogh  Changes to use M77CON
	 *>> 1994-09-08 DIVA  Krogh  Added CHGTYP code.
	 *>> 1994-07-11 DIVA  Krogh  Fix to get same state with/without var. eqs.
	 *>> 1994-03-07 DIVA  Krogh  Allow larger order in single precision.
	 *>> 1994-01-14 DIVA  Krogh  Minor change to allow changing TFINAL.
	 *>> 1993-04-27 DIVA  Krogh  Additions for Conversion to C.
	 *>> 1993-04-12 DIVA  Krogh  Converted to use slightly altered MESS.
	 *>> 1993-04-12 DIVA  Krogh  Fixed LSC so sol. saved when HMAX is small.
	 *>> 1992-10-13 DIVA  Krogh  Fixed G-Stop/discontinuity code interaction.
	 *>> 1992-09-21 DIVA  Krogh  Fixed bug in discontinuity code.
	 *>> 1992-09-09 DIVA  Krogh  Fixed bug - Var. Eqs. with discontinuities.
	 *>> 1992-08-07 DIVA  Krogh  Storage map printed only if option 10 .ne. 0
	 *>> 1992-07-16 DIVA  Krogh  Restored correct discontinuity code.
	 *>> 1992-06-16 DIVA  Krogh  Eliminate reuse of storage for option 12.
	 *>> 1992-04-08 DIVA  Krogh  Removed unused labels, 1020, 2120.
	 *>> 1992-03-30 DIVA  Krogh  Fixed bug in DIVAOP error message.
	 *>> 1992-03-12 DIVA  Krogh  Simplified DIVABU, more digits in B's.
	 *>> 1992-01-16 DIVA  Krogh  Fixed minor bug in error messages.
	 *>> 1991-12-03 DIVA  Krogh  Major change for improved error checks.
	 *>> 1991-06-17 DIVA  Krogh  Fixed bug in checking storage allocation.
	 *>> 1991-04-11 DIVA  Krogh  Fixed minor bug re. option 12 in DIVAOP.
	 *>> 1991-03-28 DIVA  Krogh  Removed check at label 650 for KORD2I<0.
	 *>> 1991-02-08 DIVA  Krogh  Changed some floats to generics
	 *>> 1990-11-08 DIVA  Krogh  Fixed bug on TSPECS on discon.
	 *>> 1990-09-14 DIVA  Krogh  Fixed bug when discon. and sol. save.
	 *>> 1990-09-13 DIVA  Krogh  Increased dimension of BETA by 1.
	 *>> 1990-09-13 DIVA  Krogh  Added one more poss. on rel. error test.
	 *>> 1990-09-11 DIVA  Krogh  Recent change messed up getting dump output.
	 *>> 1990-06-05 DIVA  Krogh  Fixed bug in noise test, comments in IVACOM.
	 *>> 1990-05-08 DIVA  Krogh  Fixed new bug when TMARK hit in DIVAG.
	 *>> 1990-04-17 DIVA  Krogh  Fixed minor problem in DIVAIN error msg.
	 *>> 1990-04-10 DIVA  Krogh  Fixed interaction between discon. & dump.
	 *>> 1990-03-23 DIVA  Krogh  Fixed bug on option "-2", see 1989-12-07.
	 *>> 1990-03-20 DIVA  Krogh  Fixed rarely occuring loop.
	 *>> 1990-01-29 DIVA  Krogh  Removed unneeded labels.
	 *>> 1989-12-14 DIVA  Krogh  Saved common block DIVAEV.
	 *>> 1989-12-07 DIVA  Krogh  Added option "2" to DIVAOP.
	 *>> 1989-11-09 DIVA  Krogh  Made GG a save var. in DIVAHC
	 *>> 1989-08-21 DIVA  Krogh  Fix out of bounds ref. to V in DIVABU
	 *>> 1989-07-26 DIVA  Krogh  Fix bug in initial dim. check
	 *>> 1989-07-21 DIVA  Krogh  Code for integrating discontinuities
	 *>> 1987-12-07 DIVA  Krogh  Initial code.
	 *
	 *--D replaces "?": ?IVA,?IVAA,?IVABU,?IVACO,?IVACR,?IVAEV,?IVAF,?IVAHC,
	 *-- & ?IVAG,?IVAIN,?IVAMC,?IVAO,?IVAOP,?IVAPR,?IVASC,?IVACE,?IVAIE,
	 *-- & ?IVAPE,?MESS
	 *
	 *
	 * Note a "*" at the start of a name is used to indicate "D" for the
	 * double precision version and "S" for the single precision version.
	 *
	 * When converting between precisions, don't forget to change the value
	 * of KDIM set in parameter statements in a variety of routines, and to
	 * adjust comments for the data statements associated with EIBND in
	 * *IVACR, and B in *IVAHC.
	 *
	 * Entries
	 *  *IVA    Main entry for starting the package.
	 *  *IVAA   Main program inside the package, calls the other routines, */
	/*          and does checks for output, and noise.  Called by the user
	 *          if reverse communication is used.
	 *  *IVABU  Back ups the solution to the current base time, if a step
	 *          that has been started must be taken over for some reason.
	 *  *IVACO  Called by user to get certain information from the common
	 *          blocks.
	 *  *IVACR  Corrects the solution, estimates errors, and selects order.
	 *  *IVADB  Subroutine to assist in debugging codes.  Called by user to
	 *          get a formatted list of all the variables used in the
	 *          integration.  Not required in usual case.
	 *  *IVADE  Needed only for delay differential equations.  This is called
	 *          by the user from the derivative subprogram.
	 *  *IVAG   Required only if the user has G-Stops, i.e. places to call
	 *          his output subroutine where certain functions have zeroes.
	 *  *IVAHC  Compute coefficients that depend on the step size history.
	 *  *IVAIN  Used to interpolate to arbitrary points.
	 *  *IVAOP  Used to process user option requests.
	 *  *IVAPR  Used to update the differences and to predict the solution
	 *          at the end of the current step.
	 *
	 * External Routines
	 *  *1MACH  Not used in the Fortran 95 version.  ("*" is "D" for double
	 *          and "R" for single precision.) This returns constants that
	 *          depend on the floating point arithmetic.  Input arguments of
	 *          1 to 4 give respectively:  underflow limit, overflow limit,
	 *          smallest relative difference between two floating point
	 *          numbers, and the largest relative difference between two
	 *          floating point numbers.
	 * DERIVS (formal) Name of subroutine to be called for computing */
	/*  OPTCHK  Used in checking storage allocation.
	 *  *MESS   Used to output error messages and diaganostic messages.
	 *          (Just MESS if no floating point is output.)
	 *  *ZERO   Called only if *IVAG is used.  Iterates to find zeros of
	 *          arbitrary (continuous) functions.
	 *
	 * Common blocks -- As a left over from the distant past, some variables
	 *   are in common so that they would be saved.
	 *  *IVAEV  Holds variables that depend on the environment.
	 *  *IVAMC  The main common block for the package.
	 *  *IVASC  The secondary common block for the package.  This contains
	 *          variables that are required for doing interpolation and is
	 *          separate to simplify saving the variables that are required
	 *          when the solution is being dumped (saved).
	 *
	 * Common variables and local variables
	 * ALPHA  (*IVAMC) Array with I-th entry = (current step size) / XI(I).
	 *   Used in computing integration coefficients.
	 * B      (*IVAHC) Array used to get started on computing integration
	 *   coefficients.  B(K) = 1. / (K*(K+1))
	 * BAKMIN (*IVADE) The largest delay at the initial point.
	 * BETA   (*IVAMC) Array with I-th entry = product (K=1,I-1) of
	 *   (current (XI(K)) / XI(K) from previous step),  BETA(1)=1.  Used in
	 *    updating the difference tables.
	 * C      (*IVAIN) Array used to hold integration/interpolation coeffs.
	 * C0     Parameter = 0. (in *IVAA,DE,CR,A,G,HC,IN,OP,PR)
	 * C1     Parameter = 1. (in *IVA,A,CR,DA,HC,IN,OP)
	 * C10    Parameter = 10. (in *IVAA,CR,OP)
	 * C1000  Parameter = 1000. (in *IVACR)
	 * C16    Parameter = 16. (in *IVAA,OP)
	 * C1M3   Parameter = .001 (in *IVAA)
	 * C1M5   Parameter = .00001 (in *IVAA)
	 * C1P125 Parameter = 1.125 (in *IVAA,HC,OP)
	 * C1P3   Parameter = 1.3 (in *IVAA)
	 * C1P4   Parameter = 1.4 (in *IVACR)
	 * C2     Parameter = 2. (in *IVAA,DE,BU,CR,IN,OP)
	 * C20    Parameter = 20. (in *IVACR)
	 * C2P5M3 Parameter = .0025 (in *IVAA)
	 * C4     Parameter = 4. (in *IVACR,OP)
	 * C40    Parameter = 40. (in *IVACR)
	 * C4096  Parameter = 4096. (in *IVAA)
	 * C6     Parameter = 6. (in *IVAA)
	 * C8M3   Parameter = .008 (in *IVAA)
	 * CM2    Parameter = -2. (in *IVACR)
	 * CM8    Parameter = -8. (in *IVACR)
	 * CMP5   Parameter = -.5 (in *IVACR)
	 * CMP75  Parameter = -.75 (in *IVAOP)
	 * CP0625 Parameter = .0625 (in *IVAA)
	 * CP1    Parameter = .1 (in *IVAA,CR,DA,HC)
	 * CP125  Parameter = .125 (in *IVACR)
	 * CP25   Parameter = .25 (in *IVAA,CR,DE,OP)
	 * CP3    Parameter = .3 (in *IVAA,OP)
	 * CP4    Parameter = .4 (in *IVAA)
	 * CP5    Parameter = .5 (in *IVAA,CR,DA,DE,HC,OP)
	 * CP5625 Parameter = .5625 (in *IVAHC)
	 * CP625  Parameter = .625 (in *IVAOP)
	 * CP75   Parameter = .75 (in *IVACR,OP)
	 * CP8    Parameter = .8 (in *IVACR)
	 * CP875  Parameter = .875 (in *IVAA, OP)
	 * CP9    Parameter = .9 (in *IVAOP)
	 * CP9375 Parameter = .9375 (in *IVACR)
	 * CQ3125 Parameter = .03125 (in *IVACR)
	 * CRBQI  Parameter = .421875 (in *IVAHC)  Initial val for computing RBQ.
	 * CSUM   (*IVAIN) Array used to contain partial sums of the integration
	 *   coefficients.  This is used to corrrect for a difference table that
	 *   has not yet been updated.
	 * D      (*IVAMC) Array to be used later to store coefficients for
	 *   integrating stiff equations.
	 *   derivatives.  Not used if option 13 is set.
	 * DISADJ (*IVAA) Value of stepsize when discontinuity is indicated. */
	/* DNOISE (*IVAMC) Used in determining if noise is limiting the
	 *   precision.  It is usually |highest difference used in correcting|
	 *   of the equation with the largest error estimate.
	 * DS     (*IVAMC) Array to be used later to store coefficients for
	 *   estimating errors when integrating stiff equations.
	 * DVC2   (*IVADB) Array used for output of variables HC to TOUT in
	 *   common block *IVAMC.
	 * E      (*IVACR) (Estimated error) / (Requested accuracy)
	 * EAVE   (*IVAMC) This is a weighted average of past values of EIMAX.
	 *   It is adjusted to account for expected changes due to step changes.
	 * EEPS10 (*IVAEV) = 10. * (machine epsilon).
	 * EEPS16 (*IVAEV) = 16. * (machine epsilon).
	 * EEPS2  (*IVAEV) =  2. * (machine epsilon).
	 * EEPT75 (*IVAEV) = (machine epsilon) ** (.75)
	 * EI     (*IVACR) Estimate for what E would be if step size increased.
	 * EIBND  (*IVACR) Array containing limits on the estimated error with
	 *   the stepsize increased.  This array tends to make the code a little
	 *   more conservative on step size increases at low order.
	 * EIMAX  (*IVAMC) Estimate of (error estimate / error requested) if the
	 *   step size should be increased.
	 * EIMIN  (*IVAMC) An error estimate is small enough to allow a step
	 *   increase if the estimate of ((error with the step size increased) /
	 *   (error requested)) is less than EIMIN.
	 * EIMINO (*IVAA) Set to C8M3 and never changed.  When step size is being
	 *   reduced if EIMIN .le. EIMINO then the reduction factor is set to
	 *   CP875.  This variable could be a parameter.
	 * EMAX   (*IVAMC) Largest value computed for (error estimate) / (error
	 *   requested).
	 * EOVEP2 (*IVAEV) = EEPS2 * (largest floating point number).
	 * EPS    (*IVACR) Current absolute error tolerance.  Also used for
	 *   temporary storage when computing the desired value of EPS.
	 * ERCOEF (*IVACR) (Error coefficient from formula) / EPS
	 * EREP   (*IVAMC) If EMAX > EREP, a step is repeated.  Ordinarily
	 *   this has the value .3.  This is set < 0 if the error tolerance is
	 *   specified improperly, and is set to a large value if the user
	 *   requests complete control over the step size.  EREP is also set
	 *   < 0 after a user specified discontinuity.
	 * EROV10 (*IVAEV) = 10. / (largest floating point number).
	 * ETA    (*IVAIN) Array used in computing integration/interp. coeffs.
	 * EVC    (*IVADB) Array used for output of variables EEPS2 to EROV10 in
	 *   common block *IVAEV.
	 * EXR    (*IVAA) Set to CP1 and never changed.  If it is estimated the
	 *   the (error estimate) / (error requested) on the next step will be
	 *   .ge. EXR then the step size is reduced.  Could be a parameter.
	 * F      (formal) Array used to store derivative values, the difference
	 *   tables, error tolerance requests, and values used by some other
	 *   options. (in *IVA,A,BU,CR,DA,DB,G,IN,PR)
	 * FDAT  (*IVAMC) Used to store data for error messages.  (Local array in
	 *   *IVAIN.)
	 * FOPT  (formal) in *IVAOP.  Passed as place to save floating point data
	 *   for options.  This package passes F in for FOPT when calling *IVAOP.
	 * G      (*IVAMC) Integration coefficients used for predicting solution.
	 *   G(I, J) gives the I-th coefficient for integrating a J-th order
	 *   differential equation.  G(1, 1) is equal to the step size.
	 * GAMMA  (*IVAIN) Array used in computing integration/interp. coeffs.
	 * GG     (*IVAHC) Array of length = max. differential equation order
	 *   allowed by code - 1.  GG(K) = (HH**(K+1)) / K!
	 * GNEW   (formal) in *IVAG.  Current value for vector function g, whose
	 *   zeroes are to be found.
	 * GOINT  (*IVACR) Used for assigned go to used in computing integration
	 *   coefficients.
	 * GOLD   (*IVAG) Previous value for element of G whose zero search is
	 *   active.
	 * GS     (*IVAMC) Integration coefficients used in estimating errors.
	 * GT     (formal) in *IVAG.  Previous value of GNEW.
	 * HC     (*IVAMC) Ratio of (new step size) / (old step size)
	 * HDEC   (*IVAMC) Default value to use for HC when reducing the step
	 *   size.  (Values closer to 1 may be used some of the time.)
	 * HH     Equivalenced to G(1,1) = current step size in *IVAA,CR,DA,G,HC.
	 * HI     (*IVAIN) Step length from the base value of the independent
	 *   variable for the interpolation.
	 * HINC   (*IVAMC) Default value to use for HC when increasing the step
	 *   size.  (Values closer to 1 may be used some of the time.)
	 * HINCC  (*IVAMC) Actual value used for default value of HC when
	 *   increasing the step size.  Set to HINC after start is considered
	 *   complete.  During the start HINCC is set to 1.125.
	 * HMAX   (*IVAMC) Largest value allowed for abs(step size).  Default
	 *   value is a very large number.
	 * HMAXP9 (*IVAMC) .9 * HMAX.
	 * HMIN   (*IVAMC) Smallest value allowed for abs(step size).  Default
	 *   value is 0.
	 * HNEW   (*IVADE) Value of step size when iterating at initial point
	 *   for delay differential equations.
	 * I      Used for temporary storage. (*IVAA,BU,CR,DA,DE,G,IN,OP,PR)
	 * IA     (*IVAOP) absolute value of first integer stored for an option.
	 * ICF    (*IVAMC) Final index for current loop in *IVACR.  Required by
	 *   option 18.
	 * ICI    (*IVAIN) Temporary index, = 0 for interpolation, 1 or 0 for
	 *   differentiation, and d-1, d-2, ... 0 for integration, where d is the
	 *   order of the differential equation.  Index of first location
	 *   in C() used is ICI + an offset.
	 * ICS    (*IVAMC) Starting index for current loop in *IVACR.
	 * ID     (formal) Array use to contain integer data from common.  Values
	 *   are returned in locations 1 to 5 as follows.
	 *   1    KEMAX = Index of equation with largest error estimate
	 *   2    KSTEP = Current step number
	 *   3    NUMDT = Number of differences used for each equation
	 *   4            Reserved for future use
	 *   5            Reserved for future use
	 * IDAT   (*IVAMC) Used to store integer for error messages.  (Also used */
	/*   in *IVAA for temporary storage of KORD(2).  (Local array in *IVAIN.)
	 * IDE    (*IVADE - formal) Array used to contain past information so
	 *   that delays can stretch back indefinitely.  If the first location is
	 *   0, then any interpolations requested must be in the range of the
	 *   current difference tables.  At present, only the value 0 is allowed
	 *   in IDE(1).  This array is intended for the support of saving long
	 *   past histories.  IDE(2) must contain the declared dimension of WDE.
	 * IDEF   (*IVADE -- formal) Flag giving indicaion of what is going on.
	 *   = 0  User should compute derivatives and return to the main
	 *        integrator.
	 *   = 1  Code is computing additional values in order to get past data
	 *        necessary for starting.  User should compute derivatives and
	 *        call *IVADE.
	 *   < 0  Indicates an error condition.  If *IVADE is called without
	 *        changing the value of IDEF, the integration is stopped and an
	 *        error message printed.  Possible error flags are:
	 *    -1  Difference tables do not span back far enough to compute the
	 *        past values of Y needed.
	 *    -2  There is not enough space in WDE to get the required starting
	 *        values.
	 * IDIMF  (formal) Declared dimension of F().
	 * IDIMK  (formal) Declared dimension of KORD().
	 * IDIMT  (formal) Declared dimension of TSPECS().
	 * IDIMY  (formal) Declared dimension of Y().
	 * IDT    (*IVAIN) Used as a base index into the difference table.
	 * IFLAG  (formal in *IVAG) Used for communication with user.
	 *   = 1  Continue as if *IVAG was not called.
	 *   = 2  Check KORD(1) as one would do at start of OUTPUT if no G-Stops
	 *        were present. (Exit if in DERIVS.)
	 *   = 3  Return to the integrator.
	 *   = 4  Compute G and return to *IVAG.
	 *   = 5  A G-Stop has been found, and NSTOP gives its index.  (If NSTOP
	 *        < 0, the stop was an extrapolating stop.)
	 *   = 6  Same as 5, but requested accuracy was not met.
	 *   = 7  Same as 5, but there is a probable error in computing G.
	 *   = 8  Fatal error of some type.  (An error message has been printed.)
	 * IG     (*IVAG)  IG = KORD(2) on the initial entry (0 for extrapolating
	 *   G-Stops, and 1 for interpolating).
	 * IGFLG  (*IVAMC) Used primarily in *ivag, but also used in *iva to keep
	 *   track of the state of GSTOP calculations.
	 *   = -2 Extrapolatory G's initialized, but not the interpolatory.
	 *   = -1 Interpolatory G's initialized, but not the extrapolatory.
	 *   =  0 Set when integration is started or restarted, or option setting
	 *        GSTOP is set.
	 *   =  1 Iterating to find a GSTOP.
	 *   =  2 User told that a GSTOP was found.
	 *   =  3 Checking G's at point where a GSTOP was located.
	 *   =  4 Checking G's at a T output point.
	 *   =  5 Usual case, no sign change detected.
	 * IGSTOP (*IVAMC) IGSTOP(k) is set in *ivag to the index of the last G
	 *   with a 0, where k is one for an interpolatory G-Stop, and k is two
	 *   for an extrapolatory G-Stop.
	 * IGTYPE (*IVAMC) Array with two elements as for IGSTOP, but this saves
	 *   a flag giving the nature of convergence to the stop.
	 *   = 0  All known G-stops completely processed.
	 *   = 4  Need to compute next value while iterating.
	 *   = 5  Got good convergence.
	 *   = 6  Got convergence, but not to desired accuracy.
	 *   = 7  Problem in getting convergence.
	 *   = 8  A fatal error of some type.
	 * IHI    (*IVA) Last location used by the current option.
	 * ILGREP (*IVAMC) Used when correction to keep track of equations that
	 *   are to use a certain error tolerance.
	 * ILGROR (*IVACR) Index of last equation in the current group of
	 *   equations grouped for selecting integration order.
	 * ILOW   (*IVA) First location used by the current option.
	 * INCOM  (*IVADE) Array equivalenced to LDT in the common block *IVASC.
	 *   Used to simplify saving information in the common block.
	 * INCOP  (*IVAOP) Array containing data giving the amount of space in
	 *   IOPT used for each of the options.
	 * INGS   Current index for G-stop being examined in DIVAG.
	 * INICAS (*IVADE) Used to track the initialization for a delay equation.
	 *   = 1  Very beginning.
	 *   = 2  Getting derivative at the very beginning.
	 *   = 3  Getting derivatives at points prior to the initial point.
	 *   = 4  Getting derivative at initial point after iteration is started.
	 * INTCHK (*IVA) Array passed to OPTCHK containing information on storage
	 *   allocation.  See comments in OPTCHK for details.
	 * INTEG  (*IVAIN) Number of integrations being done. (<0 for
	 *   differentiations and =0 for interpolation.)  Also used as counter
	 *   when computing integration coefficients.
	 *        (*IVAPR) Number of integrations being done.
	 * INTEGS (*IVAPR) = -1 for equations that are not stiff, 0 for those
	 *   that are stiff.
	 * INTEGZ (*IVAIN) min(INTEG, 0)
	 * INTERP (*IVAIN) added to the usual integration order to get the order
	 *   to be used when interpolating: 3-KQMAXI, if HI=0; 1, if
	 *   |HI| > |XI(1)| and HI * XI(1) < 0; 0, otherwise -- the usual case.
	 * IOP10  (*IVAMC) Number of times diagnostic output is to be given when
	 *   leaving *ivacr (the corrector).
	 * IOP11  (*IVAMC) Gives current step number of the method.  Tells how
	 *   many of certain coefficients must be computed. (Has nothing to do
	 *   with options.) = min(max integ order + 1, KDIM).  Also set when
	 *   starting to flag that certain memory locations must be set to 0.
	 * IOP12  (*IVAMC) Points to location in F() where user supplied values
	 *   of HINC, HDEC, HMIN, and HMAX.  (0 if option 12 not used.)
	 * IOP13  (*IVAMC) If not zero, reverse communication will be used for
	 *   getting the values of derivatives.  Associated with option 13.
	 * IOP14  (*IVAMC) If not zero, reverse communication will be used in
	 *   place of calls to the output routine.  Associated with option 14. */
	/* IOP15  (*IVAMC) If not zero, a return will be made to the user after
	 *   the initialization.  Associated with option 15.  This might be used
	 *   to overlay *iva, some of the user's code, and perhaps *ivaop.
	 * IOP16  (*IVAMC) Points to location in KORD() where information for
	 *   specifying the error tolerance is specified.  See option 16.
	 * IOP17  (*IVAMC) Used in initialization for option 17, afterwards this
	 *   cell is used by KEXIT which is equivalenced to IOP17.
	 * IOP18  (*IVAMC) Points to location in KORD() where information for
	 *   specifying a grouping of equations for derivative evaluation is
	 *   stored.  See option 18.
	 * IOP19  (*IVAMC) Points to location in KORD() where information for
	 *   specifying a grouping of equations for integration order control
	 *   is stored.  See option 19.
	 * IOP20  (*IVAMC) Used for option 20, gives first location in F where
	 *   estimated errors are to be stored.  Expected to be useful in a
	 *   program for solving boundary value problems using multiple shooting.
	 * IOP21  (*IVAMC) Was used for stiff equations option (never completely
	 *   coded).  The optional code still uses this (don't activate it!).
	 *   Now used to flag the location if F where the user has stored the
	 *    tolerance to use in finding G-Stops.
	 * IOP21S (*IVAMC) Was used for stiff equations see above.
	 * IOP22  (*IVAMC) Set aside for possible option for stiff equations.
	 * IOP3   (*IVAMC) Value set by option 3.
	 *   =  0 Interpolate to final point. (The default)
	 *   =  1 Integrate to final point.
	 *   = -1 Extrapolate to final point.
	 * IOP4   (*IVAMC) Value set by option 4.  The output routine is called
	 *   with KORD(1) = 4, every IOP4 steps.  (Default value for IOP4 is a
	 *   very large number.
	 * IOP5   (*IVAMC) Value provided by option 5, used to specify extra
	 *   output points.
	 * IOP6   (*IVAMC) Value provided by option 6.  If nonzero, the output
	 *   routine is called at the end of every step.  If > 0, there are
	 *   IOP6 interpolating G-Stops.
	 * IOP7   (*IVAMC) Value provided by option 7.  If > 0, there are K7
	 *   extrapolating G-Stops.
	 * IOP8   (*IVAMC) Value provided by option 8.  If nonzero, the output
	 *   routine is called with KORD(1)=8 whenever the step size is changed.
	 * IOP9   (*IVAMC) Value provided by option 9.  Used to specify that the
	 *   user wishes to save the solution.
	 * IOPIVA (*IVA) Used to save length of IOPT vector for error messages.
	 * IOPST  (*IVASC) Intended for possible use in stiff equations.
	 * IOPT   (formal *IVA and IVAOP) Used to specify options.
	 * IOPTC  (*IVAOP) In *IVAOP equivalenced so that IOPTC(3) is equivalent
	 *   to IOP3.
	 * IOPTS  (*IVAOP) Array containing the current default values to be
	 *   stored into IOPTC.
	 * IORD   (*IVACR) Index of first equation in the current group of
	 *   equations grouped for selecting integration order.
	 * IOUTKO (*IVADC) Used in *IVADI to point to KORD to keep track of
	 *   equation grouping for diagnostic output.
	 * ISVCOM (*IVADE) Used to save info. in the common block *IVASC.
	 * ITERS  (*IVADE) Counts iterations in starting delay differential
	 *   equations.  Max. value for this is arbitrarily 100.
	 * ITOLEP (*IVAMC) Used for temporary storage, and for the index of a
	 *   tolerance relative to the start of tolerances.
	 * IVC1   (*IVADB) Array used for output of variables IOPST to NUMDT in
	 *   common block *IVASC.
	 * IVC2   (*IVADB) Array used for output of variables ICF to NY in
	 *   common block *IVAMC.
	 * IWB    (*IVADE) Current base index for saving F values in WDE when
	 *   starting delay differential equations.
	 * IY     (*IVAMC) Used for the current index to the Y() array.  (Local
	 *   variable in *IVAIN used in computing IYI.)  Equivalenced to
	 *   IZFLAG in *IVAG.
	 * IYI    (*IVAIN) Y(IYI) is currently being computed.
	 * IYN    (*IVAIN) Y(IYN) is base Y() corresponding to Y(IYI).
	 * IYNI   (*IVAIN) Used as base index for computing IYN as IY is for INI.
	 * IYO    (*IVADE) Points to first base value of Y for current
	 *   interpolation when getting values for a delay differential equation.
	 * IZFLAG (*IVAG)  Equivalenced to IY.  Set to 0 initially, and later
	 *   set to the value returned by *ZERO.
	 *    = 0  Value set on entry at start of search.
	 *    = 1  Compute next g again.
	 *    = 2  Normal terminiation.
	 *    = 3  Normal termination -- error criterion not satisfied.
	 *    = 4  Apparent discontinuity -- no zero found.
	 *    = 5  Couldn't find a sign change.
	 *    = 6  *ZERO was called with a bad value in IZFLAG.
	 * J      For temporary storage. (In *IVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
	 * J1     (*IVAA & DA) Used for temporary storage.
	 * J2     (*IVAA) Used for temporary storage.
	 * JL     (*IVA) Used for checking storage.
	 * JLGREP (*IVACR) Contents of first location of KORD (called LGROUP in
	 *   *IVACR) for the current error tolerance rule.
	 * JLGROR (*IVACR) Contents of first location of KORD for the current
	 *   integration order control.
	 * JLIM   (*IVA) Used for checking second item in KORD list for options
	 *   16 and 19.
	 * K      For temporary storage.  (In *IVA,A,BU,CR,DA,DB,DE,HC,IN,OP,PR)
	 * KDIM   Parameter giving the largest number of differences supported.
	 *        Used in all the routines.
	 * KEMAX  (*IVAMC) Index associated with equation giving the largest
	 *   value for (estimated error) / (requested error).
	 * KEXIT  (*IVAMC) Equivalenced to IOP17 which is not used after
	 *   initialization.  Defines actions when KORD2I = -7.  (Referenced in
	 *   (*IVAA,DA,G).)
	 *   =  1  Take the step over with reduced H.
	 *   =  2  Take the step over.
	 *   =  3  Do the end of step call to OUTPUT. */
	/*   =  4  Reset TMARK, then do same as for KEXIT = 2.
	 *   =  5  Reset TMARK, then do same as for KEXIT = 3.
	 *   =  6  Give the fatal error diagnostic.
	 * KFERR  (*IVA)  Temporary storage in checking for option 16.
	 * KGO    (*IVA)  Used to tell from whence a check is being done or an
	 *   error message is being written.
	 *   = 1 Checking an equation group for variational equations.
	 *   = 2 Checking an equation group for diagnostic print.
	 *   = 3 Checking an equation group for integration order control.
	 *   = 4 Checking an equation group for error control.
	 *   = 5 Checking an equation group for specifying ODE orders.
	 *   = 6 Found a problem with output type for printing.
	 *   = 7 Found a problem with an output group for printing.
	 *   = 8 Found a problem with input NEQ.
	 *   = 9 Order specified for the ODE's in the system is out of range.
	 *   =10 Option 16 was not used (an error).
	 *   =11 Error tolerance of 0 specified without proper flags.
	 * KIS    (*IVAMC) Used to check if it is time to dump the solution.
	 *   The check involves incrementing KIS at the end of the step, and
	 *   dumping the solution if KIS is 0.
	 *   = -1  Set in *ivacr when it is time to dump solution
	 *   =  0  When starting
	 *   =  2  After being dumped.
	 *   This is set to 1000 just after a user specified discontinuity, and
	 *   counted up from that point.
	 * KMARK  (*IVAMC) Identifies the type of output associated with the next
	 *   output point specified by TSPECS.
	 * KONV   (*IVADE) Counts iterations.  Test for convergence if KONV > 1.
	 * KORD   (formal in *IVA,A,BU,CR,DA,DB,DE,G,IN,PR) KORD(1) is used to
	 *   return flags for the user to test, and KORD(2) tells what routine
	 *   the flag is associated with.  See KORD1I and KORD2I below and the
	 *   write up for the program.  KORD(3) is used for communicating extra
	 *   information to the user in some cases.  KORD(4) to KORD(NTE+3) are
	 *   used for integration order for the equations, and the rest of KORD()
	 *   is available for user options.
	 * KORD1I (*IVAMC) Helps in defining the state of the integrator.
	 *   Frequently has the same value as KORD(1).  Meaning depends on the
	 *   value of KORD(2), or the value about to be assigned to KORD(2).
	 *   <  0  Happens when preparing to give output with extrapolation.
	 *   =  0  Happens when checking F at points for noise test.
	 *   =  1  (KORD(2)=-1)  End of integration has been reached.
	 *   =  1  (KORD(2)= 0)  Computing first predicted derivative.
	 *   =  1  (KORD(2)= 1)  Output for initial point.
	 *   =  2  (KORD(2)=-1)  Giving diagnostic for noise limiting precision.
	 *   =  2  (KORD(2)= 0)  Computing corrected derivative.
	 *   =  2  (KORD(2)= 1)  Output for TSPECS(3).
	 *   =  3  (KORD(2)=-1)  Diagnostic for step size reduction too fast.
	 *   =  3  (KORD(2)= 0)  Computing variational derivative.
	 *   =  3  (KORD(2)= 1)  Output for TSPECS(4).
	 *   =  4  (KORD(2)=-1)  Error, discontinuity.
	 *   =  4  (KORD(2)= 1)  Output for certain number of steps.
	 *   =  5  (KORD(2)= 0)  Get initial derivatives for stiff equations.
	 *   =  5  (KORD(2)= 1)  Extra output from TSPECS.
	 *   =  6  (KORD(2)= 1)  End of step output.
	 *   =  7  (KORD(2)= 0)  Evaluate G before extrapolated output point.
	 *   =  7  (KORD(2)= 1)  Evaluate G before extrapolated output point.
	 *                       (Also used when checking for other G's after
	 *                        finding one.)
	 *   =  8  (KORD(2)= 1)  Tell user step size has changed.
	 *   =  9  (KORD(2)= 1)  Request for user to save solution.
	 *   = 11  (KORD(2)=-1)  Error, step size too small at end of start.
	 *   = 12  (KORD(2)=-1)  Error, step size is too small.
	 *   = 13  (KORD(2)=-1)  Error, output points specified badly.
	 *   = 21  (KORD(2)=-1)  H too small to give reasonable change when added
	 *                       to T.
	 *   = 22  (KORD(2)=-1)  Error, bad tolerance.
	 *   = 23  (KORD(2)=-1)  Set after message for a fatal error.
	 *   = 24  Set on error message in *iva, along with KORD2I = -4.
	 *   Also used as an index into MLOC in *IVAA when an error is being
	 *   processsed, see MLOC below.
	 * KORD2I (*IVAMC) Helps in defining the state of the integrator.
	 *   Frequently has the same value as KORD(2).
	 *   = -3  Set in *ivag, to get a derivative evaluation.
	 *   = -2  Set in *ivag, to get another entry to OUTPUT.
	 *   = -1  Return to calling program, done, interrupt, or got an error.
	 *   =  1  Calling OUTPUT or returning to user for OUTPUT type action.
	 *   =  0  Calling DERIVS or returning to user for DERIVS type action.
	 *   = -4  Error message in *iva and in *ivaop, along with KORD1I = 24.
	 *   = -5  Starting
	 *   = -6  Starting, getting the initial derivative value or derivatives
	 *         for the noise test.
	 *   = -7  Done some extrapolation, KEXIT defines the action to take.
	 *         Set in *ivag to activate KEXIT action in *iva.
	 *   = -8  Set when user has requested adjustment of the difference
	 *         tables for a discontinutiy.
	 * KORDI  (*IVASC) Order of differential equation being integrated.  If
	 *   all orders are the same, this set once at the beginning.
	 * KOUNT   (*IVADE) Count of number of points back from the initial point
	 *   when solving a delay differential equation.
	 * KOUNTM  (*IVADE) Largest value currrently allowed for KOUNT.
	 * KOUNTX  (*IVADE) Largest value allowed for KOUNTM.
	 * KOUTKO  Used in DIVACR to track where output is wanted.
	 * KPRED  (*IVAMC) Value assigned to KORD1I when getting a predicted
	 *   derivative.  (1 used now, 5 planned for use with stiff equations.)
	 * KQD    (*IVACR) = max(2, integration order)
	 * KQDCON (*IVAMC) Number of coefficients computed with constant step
	 *   size for stiff equations.
	 * KQICON (*IVAMC) Number of coefficients computed with constant step
	 *   size for nonstiff equations.
	 * KQL    (*IVACR) Integration order at start of (*IVACR) */
	/* KQLORD (*IVACR) Saved value of KQL when equations are grouped for
	 *   controlling the integration order.
	 * KQMAXD (*IVASC) Maximum integration order used for stiff equations.
	 * KQMAXI (*IVASC) Maximum integration order used for nonstiff equations.
	 * KQMAXS (*IVAMC) Maximum integration order for equations that have
	 *   some limit on the error that can be committed.
	 * KQMXDS (*IVAMC) Used to save KQMAXD in case step is repeated and the
	 *   solution must be dumped.
	 * KQMXI  (*IVAIN) Maximum integration order used for integration or
	 *   interpolation, = KQMAXI+INTERP-1.
	 * KQMXS  (*IVAIN) Maximum step number, = max(KQMXI, KQMAXD).
	 * KQMXIL (*IVAMC) Value of KQMAXI the last time integration coefficients
	 *   were computed.
	 * KQMXIP (*IVAMC) = KQMAXI + MAXINT, for computing integration coeffs.
	 * KQMXIS (*IVAMC) Used to save KQMAXI in case step is repeated and the
	 *   solution must be dumped.
	 * KQN    (*IVACR) Value of integration order at end of *IVACR.
	 * KQQ    Used for the integration order for current equation.  (Values
	 *   < 0 are intended for stiff equations.)  (In *IVA,BU,DA,IN,PR)
	 * KSC    (*IVAMC) Number of steps that have been taken with a constant
	 *   step size.
	 * KSOUT  (*IVAMC) When KSTEP reaches this value, the output routine is
	 *   called with KORD(1) = 4.  The default value is a very large number.
	 * KSSTRT (*IVAMC) Set when ending one derivative per step to KSTEP + 2.
	 *   Checked later in *IVAHC to decide whether to set the step changing
	 *   factors to their nominal values.
	 * KSTEP  (*IVAMC) Number of steps taken since the start of integration.
	 * L      Used for temporary storage.  In *IVAIN, L is the initial value
	 *   of LDT, except L=1 if LDT=-1, and MAXINT .ge. 0.  (Used in *IVAA,BU
	 *   CR,DA,DB,IN,PR.)
	 * LAHAG  (*IVADB) Used to get proper offset into an diagnostic message.
	 * LAIAG  (*IVADB) Used to get proper offset into an diagnostic message.
	 * LDIS   (*IVAA) Count of steps since user flagged a discontinuity.
	 * LDT    (*IVASC) Used to keep track of state of difference table.
	 *   = -5  Used only on first step to indicate that an extra iteration
	 *         is desired to get a firm estimate on the error.
	 *   = -4  Set on initialization before there is any difference table.
	 *   = -3  Set just after predicting, interpolation is not allowed when
	 *         this value is set.
	 *   = -2  Set when difference table is to be updated to the end of the
	 *         current step, but no interpolation is to be done.  (For
	 *         dumping the solution.)
	 *   =  0  Calculations for current step are complete, it is not o.k. to
	 *         update the difference table.
	 *   =  1  Difference table has been updated to the end of the current
	 *         step, e.g. by doing an interpolation.
	 *   =  2  Set when doing a special interpolation during computation of
	 *         derivatives.  (For delay equations.)
	 * LEX    (*IVAMC) Indicates how to get values at next output point:
	 *   = -1  Extrapolate
	 *   =  0  Interpolate (The usual case.)
	 *   =  1  Integrate to the output point, integration is not continued.
	 * LGO    (*IVAIN) Used as an an assigned go to.  Result is to add in
	 *   extra correction term when LDT has been set to 2.
	 * LGROUP (formal) This is a part of KORD passed into *IVACR.  The first
	 *   location is the start of the information on the grouping of
	 *   equations for error control.
	 * LINC   (*IVAMC) Used to indicate state of step size selection.
	 *   = -10 After computed derivatives at base time, after computing other
	 *         extra derivatives for the noise test.
	 *   = -9  After computed second extra derivative for noise test.
	 *   = -8  After computed first extra derivative for noise test.
	 *   = -7  Dumping the solution and then doing a user initiated restart,
	 *         or getting ready to compute extra derivatives for the noise
	 *         test.
	 *   = -6  Dumping the solution before a restart.
	 *   = -5  Set on the first step, and also set when dumping the solution
	 *         after a discontinuity.
	 *   = -4  Repeat step with no change in the step size.
	 *   = -3  Set when the error tolerance is set improperly.
	 *   = -2  User has complete control of selecting the step size.
	 *   = -1  Step is being repeated.
	 *   =  0  Step size is not to be increased on this step.
	 *   = k>0 Step size can be increased by HINCC**k.
	 * LINCD  (*IVAMC) Value of smallest k for which HINCC**k .ge. 2.
	 *   (=-2 if user is specifying all step size changes.)
	 * LINCQ  (*IVAMC) Value of smallest k for which HINCC**k .ge. 4.
	 * LIOPT  (*IVAOP) Value of the last index in IOPT on the last call.
	 *   Used so *IVA can print IOPT in error messages.
	 * LL     (*IVACR) Temporary variable used when equations are grouped
	 *   for integration order control.
	 * LNOTM1 (*IVAIN) Logical variable = L .ne. -1.  If LNOTM1 is true,
	 *   storage in Y() is different in some way lost to antiquity.  Such
	 *   a case can only arise in the case of stiff equations.
	 * LOCF1  (*IVADB) Gives packed data needed for output of tables by the
	 *   message processor MESS.  See comments there under METABL for defs.
	 * LOCF2  (*IVADB) As for LOCF1 above.
	 * LOCM   (*IVAA) Parameter = 32*256, used to unpack integers stored
	 *   in MLOC for use in error message processing.
	 * LPRINT (formal, *IVADB) Defines how much printing is to be done in
	 *   *IVADB.  Let |LPRINT| = 10*N1 + N2     (N1,N2 digits)
	 *    N1=1   Do not print any variables external to the integrator.
	 *    N1=2   Print  tspecs, current y, past y, current f, all pertinent
	 *           contents of KORD, and TOL.
	 *    N1=3   Above + difference tables up to highest difference used.
	 *    N1=4   Same as N1=1 + all in storage allocated for differences.
	 *    N2=1   Do not print any variables internal to the integrator.
	 *    N2=2   Print all scalar variables in interpolation common block.
	 *    N2=3   Above + all scalar variables in main integ. common block.
	 *    N2=4   Same as N1=3 + all used in arrays XI,BETA,ALPHA, first */
	/*           column of G, GS,RBQ,SIGMA
	 *    N2=5   Same as N1=4 + all used in arrays G,D,DS,V
	 * LSC    (*IVAMC) Indicates if starting or if noise may be present.
	 *   =k<0 -k steps have been taken for which noise appears to be limiting
	 *        the precision.
	 *   = 0  Usual case
	 *   = 1  Doing 1 derivative per step after initial part of start.
	 *   = 2  Used as flag that it is time to set LSC=0.
	 *   = 3  Third step, hold the order constant.
	 *   = 4  Second step, increase orders from 2 to 3.
	 *   = 5  First step, third time through the first step (if required).
	 *   = 6  First step, second time through.
	 *   = 7  First step, first time through.
	 *   = 8  Set on initialization.
	 * LTXT?? Names of this form are used in setting up data statements for
	 *   error messages.  These names are generated automatically by PMESS,
	 *   the program that makes up these messages.
	 * LX     (*IVAA) Used for temporary storage in computing TMARKA().
	 *        ( formal *IVADE)  An integer array containing extra
	 *   information, as follows.
	 *  LX(1) Points to a location in Y beyond those already in use.  Values
	 *        of Y requested are computed at TSPECS(1) - Y(LX(1)) and stored
	 *        starting at Y(LX(1)+1).  If this index is 0, no more extra Y
	 *        values are to be computed.
	 *  LX(2) Index of the first equation for which the Y's above are to be
	 *        computed.  Y(LX(1)+1) will correspond to this first equation
	 *        index.
	 *  LX(3) Index of the last equation for which the Y's above are to be
	 *        computed.  Thus the Y's stored starting at Y(LX(1)+1) will
	 *        require no more space than half the space ordinarily required
	 *        for the array Y(), and may require significantly less.
	 *  LX(4) Maximum number of times to integrate F to get Y.  This should
	 *        be > 0, and less than or equal to the order of the highest
	 *        order differential equation.  (= 0 is allowed, but probably
	 *        not what you want.  It would give a value only for F.)  Space
	 *        must be set aside for all integrals of F, even if not all are
	 *        requested.  For a first order system, all Y's are just the
	 *        first integrals of the corresponding F's.  For higher order
	 *        equations, the first Y associated with a given F is the d-th
	 *        integral of the corresponding F, where d is the order of the
	 *        equation, and the last Y corresponding to the F is the first
	 *        integral of that F.
	 *  LX(5) As for LX(4), but gives the index for the fewest number of
	 *        times to integrate F.  Ordinarily this should be > 0.  If 0 is
	 *        requested, an estimate for the value of F at the delay point is
	 *        computed.  This should not be 0 more than once, for equations
	 *        covering the same index, since later such requests would write
	 *        over the earlier results.
	 *  LX(5i+k) , k = 1, 2, ... 5.  Treated as for the cases above.  If i
	 *        different cases of delayed Y's are to be computed, then
	 *        LX(5i+1) must be 0.
	 * LX2    (*IVADE) Value of LX(5i+2), when working on the i-th delay.
	 * MACT   Used in the programs which call the error message program.
	 *   This array difines the actions to be taken by that program.  (In
	 *   (*IVA,A,DA,DE,G,IN,OP)
	 * MACT0  (*IVADB) Used to call the message program, see MACT.
	 * MACT?  As for MACT, in (*IVA,CR,DB)
	 * MACTFV (*IVADB) As for MACT0.
	 * MAXDIF (*IVASC) Maximum differentiations required for stiff equations.
	 * MAXINT (*IVASC) Maximum integrations required.  (= max. order of
	 *   differential equations if equations are not stiff.)
	 * MAXKQ  (*IVA, BU)e
	 * MAXKQD (*IVAMC) Largest integration order allowed for stiff equations.
	 * MAXKQI (*IVAMC) Largest integ. order allowed for nonstiff equations.
	 * ME???? Parameters defining constants used for interaction with the
	 *   error message program MESS.  See comments there for definitions.
	 *   (In *IVA,A,DA,DE,G,IN,OP)
	 * METHOD (*IVAMC) Defines kind of methods being used.
	 *   = -1  Only stiff equations are being integrated.
	 *   =  0  Only nonstiff equations are being integrated.
	 *   =  1  Both kinds of methods are required.
	 * MLOC   (*IVA,A,DE) Contains locations in MTEXT for error messages.  In
	 *   *IVAA this data is packed using MLOC??, see below.
	 * MLOC?? (*IVAA) Parameters constructed to aid in making up packed data
	 *   for processing error messages.  Low two digits give the value of
	 *   KORD1I to use for the error index and later processing, the next two
	 *   give the error severity level, and the rest point to text used for
	 *   the message.
	 * MODF2  (*IVADB) Used in constructing the same kind of packed data as
	 *   described for LOCF1 above.
	 * MULTJ  Local to DIVAOP for calls not using F.
	 * MTEXT  (*IVA,A,CR,IN,OP) Text for error messages.
	 * MTXT?? (*IVA,A,CR,DA,DB,DE,G,IN,OP) Equivalenced into MTEXT.
	 * N      Used for temporary storage.  (In *IVAHC,IN,PR)
	 * NDTF   (*IVASC) Location in F() where difference table starts.
	 * NE     (*IVAMC) Number of equations in the first group.  (=NTE if
	 *   option 18 is not used.)
	 * NEDDIG (*IVADB) Parameter = -MEDDIG.
	 * NEPTOL (*IVAMC) Used for temporary storage and to save the value of
	 *   ITOLEP for error messages.
	 * NEQ    (formal) Total number of equations being integrated.
	 * NG     (*IVAMC) Used in *ivag for the number of g's in the current
	 *   context.
	 * NGSTOP (*IVAG) Dimension 2 array equivalenced to IOP6, and IOP7.  To
	 *   get the number of interpolating and extrapolating G-Stops.
	 * NGTOT  (*IVAMC) NGTOT(1) gives the number of interpolating G-Stops,
	 *   and NGTOT(2) gives the number of extrapolating G-Stops.
	 * NKDKO  (*IVASC) If this is nonzero (option 17), it gives the location
	 *   in KORD() where a vector defining the order of each equation is
	 *   specified. */
	/* NLX    (*IVADE) Temporary index used to keep track of interpolations
	 *   being done to get Y() values for a delay differential equation.
	 * NOISEQ (*IVAMC) max(2, order of equation for which (error estimate)/
	 *   (error requested) is a maximum).
	 * NOUTKO (*IVAMC) If nonzero, gives the index in KORD where information
	 *   on what equations are to be included in the diagnostic output is
	 *   given.   See option 10.
	 * NSTOP  (formal) In *IVAG.  Index of the G-stop, see IFLAG.
	 * NTE    (*IVASC) Total number of equations being integrated = NEQ.
	 * NTEXT  (formao *IVADB) Character variable containing heading text.
	 * NTOLF  (*IVAMC) First location in F() where tolerance specifying
	 *   accuracy desired is stored.
	 * NUMDT  (*IVASC) Maximum allowed number of differences available for
	 *   doing an integration.
	 * NXTCHK (*IVA) Equivalenced to INTCHK(1), which gives the next
	 *   available location in INTCHK for storing data on storage allocation.
	 * NY     (*IVAMC) Total order of the system.
	 * NYNY   (*IVASC) Location in Y() where the base value for Y() is saved.
	 * NYNYO  (*IVADE) Equivalenced to the saved value from common of NYNY.
	 * OUTPUT (formal) Name of subroutine to be called for the output of
	 *   data or for computing G-Stops.  Not used if option 14 is set.
	 * OVD10  (*IVAEV) (largest floating point number) / 10.
	 * OVTM75 (*IVAEV) (largest floating point number) ** (-.75)
	 * RBQ    (*IVAMC) Array containing data for the preliminary noise test.
	 * RD     (formal *IVACO) Array use to contain floating point data from
	 *   common.  Values are returned in locations 1 to 3 as follows.
	 *   1    EMAX =  Max. ratio of estimated error to requested error
	 *   2            Reserved for future use
	 *   3            Reserved for future use
	 * REF    (*IVACR) Array of length 3 used for translating error tolerance
	 *   type into the factor used for exponential averaging for that type.
	 * RND    (*IVACR) Usually the current estimated error.  Used in deciding
	 *   if noise is limiting precision.
	 * RNOISE (*IVACR) Value used in comparison with RBQ() for preliminary
	 *   noise test.
	 * ROBND  (*IVAMC) Used to influence the selection of integration order.
	 *   The larger ROBND, the harder it is to increase the order and the
	 *   easier it is to decrease it.
	 * RVC2   (*IVADB) Array used for output of variables DNOISE to SNOISE in
	 *   common block *IVAMC.  These are variables that don't require a great
	 *   deal of precision.
	 * S      (*IVACR) Estimate of (step size) * eigenvalue of Jacobian.
	 * SIGMA  (*IVAMC) The k-th entry of this array contains a factor that
	 *   gives the amount the k-th difference is expected to increase if the
	 *   step size in increased.  These numbers get bigger it there is a past
	 *   history of increasing the step size.
	 * SIGMAS (*IVAA) Saved value of SIGMA(k) from the last step, where k =
	 *   integration order for equation with index KEMAX.
	 * SNOISE (*IVAMC) Value used in comparison with RBQ() on equation with
	 *   largest value for (error estimate) / (error requested).
	 * T      (formal) in *IVAIN. T(1) contains the point to be interpolated
	 *   to, and T(2) is used in a check that |HI| .le. |T(2)|.  When used by
	 *   other routines in this package, TSPECS is passed in for T.
	 * TB      (*IVADE) Base time for current interpolation.
	 * TC      (*IVADE) Original value of TN when getting past Y's for a
	 *   delay differential equation.
	 * TEMP   Used for temporary storage, in *IVAHC,PR
	 * TEMPA  (*IVACR) Array equivalenced to (TPS1,TPS2,TPS3,TPS4).
	 * TEMPAO (*IVACR) Array used to accumulate values in TEMPA.
	 * TG     (*IVAMC) TG(1) gives the last value of TSPECS(1) for which an
	 *   interpolatory G-Stop has been computed.  TG(2) is defined similarly
	 *   for extrapolatory G-Stops.
	 * TGSTOP (*IVAMC) TGSTOP(1) gives the value of TSPECS(1) where the last
	 *   0 for an interpolatory G-Stop was found.  TGSTOP(2) is defined
	 *   similarly for extrapolatory G-Stops.
	 * TMARK  (*IVAMC) Location of the next output point.
	 * TMARKA (*IVAA)  Array of length 2 equivalenced to TMARK (and TMARKX).
	 * TMARKX (*IVAMC) Location of the next output point to be found using
	 *   integration or extrapolation.  This variable must follow immediately
	 *   after TMARK in the common block.
	 * TN     (*IVASC) The value of TSPECS(1) at the conclusion of the last
	 *   step.
	 * TNEQ   (*IVADB) Array of dimension 1 equivalenced to TN so that an
	 *   array can be passed to *MESS.
	 * TOL    (formal) This is a part of F passed into *IVACR.  The first
	 *   location is the start of the information on the tolerances for error
	 *   control.
	 * TOLD   (*IVAG) Value of TSPECS(1) on one side of a zero.
	 * TOLG   (*IVAMC) Tolerance to pass to dzero when locating G-Stops.
	 * TOUT   (*IVAMC) Location of next output point defined by value of
	 *   TSPECS(3).  Such output is given with KORD(1) = 2.
	 * TP     (*IVA,A,DA,DE,HC) Used for temporary storage.
	 * TP1    (*IVAA,DA,HC,IN,PR) Used for temporary storage.
	 * TP2    (*IVAA,DA,HC,PR) Used for temporary storage.
	 * TP3    (*IVAA) Used for temporary storage.
	 * TPD    (*IVABU) Used for temporary storage.
	 * TPP    (*IVACR) Used for temporary storage.  Usually same as TPS3.
	 * TPS1   (*IVAA,CR) Used for temporary storage.  (In *IVACR is the
	 *   difference of order KQQ-2)
	 * TPS2   (*IVAA,CR) Used for temporary storage.  (In *IVACR is the
	 *   difference of order KQQ-1)
	 * TPS3   (*IVACR) Contains the difference of order KQQ.  This is the
	 *   last difference used in the corrector.
	 * TPS4   (*IVACR) Contains the difference of order KQQ+1.
	 * TPS5   (*IVACR) Temporary storage.
	 * TPS6   (*IVACR) Temporary storage.
	 * TPS7   (*IVACR) Temporary storage.
	 * TSAVE  (*IVAG) Value of TSPECS(1) before starting the search for a 0.
	 * TSPECS (formal *IVA,A,DB,DE,G)
	 *   TSPECS(1) is the current value of the independent variable. */
	/*   TSPECS(2) is the current value of the step size.
	 *   TSPECS(3) is the increment to use between output points that give
	 *             output with KORD(1) = 2.
	 *   TSPECS(4) is the "final" output point.
	 * V      (*IVAMC) Array used in computing integration coefficients.
	 * XI     (*IVASC) XI(K) = TSPECS(1) - value of TSPECS(1) K steps
	 *   previous.
	 * W      (*IVAHC) Array used in computing integration coefficients.
	 * WDE    (formal, *IVADE)  Array used for working storage.  This storage
	 *   is used to save derivative values when iterating to get started.  To
	 *   be safe one should allow as much space as is allowed for differences
	 *   in F.  In most cases the start will not require this much space
	 *   however.  This array is also intended for the support of saving long
	 *   past histories.
	 * Y      (formal, *IVA,A,CR,DA,DB,DE,G,IN,PR) Array containing the
	 *   independent variable and all derivatives up to order one less than
	 *   the order of the differential equation.  Also use to save these
	 *   values at the beginning of the current step, the base values.
	 * YN     (formal, in *IVAPR)  Base values of y, these follow the
	 *   current values of the dependent variable, y, in Y().
	 *
	 *
	 *++S Default KDIM = 16
	 *++  Default KDIM = 20
	 *++  Default MAXORD = 2, MAXSTF = 1
	 *++  Default INTEGO, VAREQ, OUTPUT, DUMP, GSTOP, EXTRAP
	 *++  Default STIFF=.F., ARGM=.F., ERRSTO=.F.
	 * */
	/*--D Next line special: P=>D, X=>Q */
 
	/* *********************** Internal Variables ***************************
	 *
	 * Comments for variables used in this package can be found in the file
	 *   IVACOM.
	 *
	 * *********************** Type Declarations ****************************
	 * */
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
	/*.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS. */
 
 
 
 
	/*                      Declarations for error message processing.
	 * */
 
 
	/* ********* Error message text ***************
	 *[Last 2 letters of Param. name]  [Text generating message.]
	 *AA DIVA$B
	 *AB The interval [1, 10**6], bounds the allowed values for NTE=$I.$E
	 *AC For option $I, the interval [$I, $I], bounds the allowed $C
	 *   values for the integration order which is set to $I.$E
	 *AD Option 16 must be used for error control.$E
	 *AE F($I) = $F, but it must be -1.0 when skipping the error check.$E
	 *AF For option $I, the interval [$I, $I] bounds the allowed $C
	 *   values for KORD($I)=$I, which is used to specify an $B
	 *AG output type for printing.$E
	 *AH output group for printing.$E
	 *AI equation group for variational equations.$E
	 *AJ order for a differential equation.$E
	 *AK equation group for diagnostic print.$E
	 *AL equation group for integration order control.$E
	 *AM equation group for error control.$E
	 *AN Option 5 argument must be .le. 0 or .gt. 4.$E
	 *   $
	 *AO KORD values for this option starting at KORD($M) are:$E */
	/* End of automatically generated error message code.
	 *
	 *        for KGO =     1      2      3      4      5      6      7 */
	/*           KGO        8      9     10     11      12
	 *
	 *                      1  2  3 4       5  6       7       8       9 10 */
	/*           11      12 13     14 15     16 */
 
	/* ************** START OF EXECUTABLE CODE ******************
	 *
	 *     **** TEST IF CONTINUING AN INTEGRATION */
	if (Kord[1] != 0)
		goto L_330;
	/*     **** INITIALIZE VARIOUS SCALARS */
	divamc.kstep = 0;
	divamc.kqmxis = 0;
	divamc.kord2i = -5;
	Kord[2] = -1;
	divasc.nte = neq;
	divamc.ne = divasc.nte;
	divamc.tolg = 0.e0;
	/*     **** SET UP OPTIONS */
	if (Iopt[1] != 0)
		divaop( iopt, f );
	divaop( iopiva, f );
	if (Iopt[1] == 0)
		Iopiva[2] = 1;
 
	if ((__builtin_expect(divamc.ne <= 0,0))) || (__builtin_expect(divamc.ne > 1000000,0)))
	{
		Idat[1] = divamc.ne;
		kgo = 8;
		goto L_650;
	}
	/*                         Set up diagnostic print on storage allocation. */
	intchk[0] = 245;
	if (divamc.iop10 != 0)
		intchk[0] = 247;
 
	/*     **** CHECK TSPECS STORAGE ALLOCATION */
	intchk[2] = idimt;
	intchk[3] = 4;
	*nxtchk = 4;
	if (divamc.iop5 != 0)
	{
		intchk[4] = 5;
		intchk[5] = 5;
		if (divamc.iop5 > 0)
		{
			intchk[6] = divamc.iop5 - 4;
			if (divamc.iop5 < 5)
			{
				kgo = 12;
				goto L_600;
			}
		}
		else
		{
			ihi = -divamc.iop5;
			jl = 4;
			for (ihi = ihi; ihi <= (idimk - 3); ihi += 3)
			{
				j = labs( Kord[ihi] );
				if (j == 0)
					goto L_20;
				if (labs( Kord[ihi + 2] ) > 1)
				{
					Idat[2] = -1;
					Idat[3] = 1;
					kgo = 6;
					goto L_600;
				}
				if ((j <= jl) || (j > Kord[ihi + 1]))
				{
					kgo = 7;
					Idat[2] = jl + 1;
					Idat[3] = Kord[ihi + 1];
					goto L_610;
				}
				jl = Kord[ihi + 1];
			}
			if (Kord[ihi] != 0)
				ihi += 3;
L_20:
			intchk[6] = jl - 4;
		}
		*nxtchk = 7;
	}
L_25:
	optchk( intchk, iopt, "DIVA / TSPECS$E" );
	if (*nxtchk < 0)
		divamc.kord2i = -4;
 
	/*     **** CHECK KORD STORAGE ALLOCATION */
	intchk[2] = idimk;
	intchk[3] = divamc.ne + 3;
	*nxtchk = 4;
	if (divamc.iop5 < 0)
	{
		intchk[4] = 5;
		intchk[5] = -divamc.iop5;
		intchk[6] = ihi + divamc.iop5;
		*nxtchk = 7;
	}
 
	/*++  Code for VAREQ is active */
	if (divamc.iop18 != 0)
	{
		divamc.ne = labs( Kord[divamc.iop18] );
		intchk[*nxtchk] = 18;
		ilow = divamc.iop18;
		kgo = 1;
		/*.       **** CHECK OPTION FOR VALID INPUT */
		goto L_430;
	}
	/*++  End */
L_30:
	;
	if (divasc.nkdko != 0)
	{
		/*                        **** STORAGE ALLOCATED FOR ODE ORDERS */
		intchk[*nxtchk] = 17;
		intchk[*nxtchk + 1] = divasc.nkdko;
		intchk[*nxtchk + 2] = divasc.nte;
		*nxtchk += 3;
	}
	/*++  Code for STIFF is inactive
	 *      IF (IOPST .ne. 0) then
	 *         INTCHK(NXTCHK) = 17
	 *         INTCHK(NXTCHK+1) = IOPST
	 *         INTCHK(NXTCHK+2) = NTE
	 *         NXTCHK = NXTCHK + 3
	 *      end if
	 *++  End
	 *
	 * **** SET INITIAL INTEGRATION ORDERS, TEST ODE ORDERS ****
	 * */
	divasc.maxint = 0;
	divasc.maxdif = 0;
	divamc.ny = 0;
	for (k = 1; k <= divasc.nte; k++)
	{
		if (divasc.nkdko != 0)
			divasc.kordi = Kord[divasc.nkdko + k - 1];
		divamc.ny += labs( divasc.kordi );
		/*++  Code for STIFF is inactive
		 *      IF (IOPST .EQ. 0) GO TO 60
		 *c.    **** CHECK FOR POTENTIAL STIFF EQUATION
		 *      JS = abs(KORD(IOPST+K-1)) - 1
		 *      IF ( JS ) 52,60,54
		 *c.    **** EQUATION IS NOT ACTIVE
		 *   52 KQQ = 0
		 *      GO TO 56
		 *c.    **** EQUATION USES IMPLICIT METHOD
		 *   54 KQQ = -1
		 *      IF (JS .GT. abs(KORDI)) then
		 *        Set up an error message.
		 *      end if
		 *      MAXINT = max(MAXINT, abs(KORDI) - JS)
		 *   56 IF (KORDI .GE. 0) GO TO 70
		 *      KORDI = -1 - KORDI
		 *      JS = JS - 1
		 *      MAXDIF = max(MAXDIF, JS, 1)
		 *      GO TO 70
		 *++  End
		 *     **** EQUATION IS TO USE AN EXPLICIT METHOD */
L_60:
		kqq = 1;
		divasc.maxint = max( divasc.maxint, divasc.kordi );
L_70:
		if ((__builtin_expect(divasc.kordi > MAXORD,0)) || (__builtin_expect(divasc.kordi <= 0,0)))
		{
			/*                    Set up error message.  KORDI is out of range. */
			Idat[1] = 17;
			Idat[2] = 1;
			Idat[3] = MAXORD;
			if (divasc.nkdko != 0)
			{
				kgo = 5;
				ilow = divasc.nkdko;
				ihi = divasc.nkdko + k - 1;
				goto L_640;
			}
			else
			{
				kgo = 9;
				Idat[4] = divasc.kordi;
				goto L_650;
			}
		}
		Kord[k + 3] = kqq;
	}
	/*     **** SET FLAGS WHICH DEPEND ON METHOD USED
	 *++  Code for STIFF is inactive
	 *      METHOD = 1
	 *      IF (MAXINT .GT. 0) IF (MAXDIF) 85,90,85
	 *      METHOD = -1
	 *   85 CONTINUE
	 *      KPRED = 5
	 *      GO TO 100
	 *++  End */
L_90:
	divamc.method = 0;
	divamc.kpred = 1;
L_100:
	;
 
	/* ******* CHECK KORD FOR DIAGNOSTIC OUTPUT CONTROL *********
	 *
	 *++  Code for OUTPUT is active */
	if (divamc.iop10 > 0)
	{
		if (divamc.noutko != 0)
		{
			intchk[*nxtchk] = 10;
			ilow = divamc.noutko;
			/*.    **** Check option for valid input */
			kgo = 2;
			goto L_430;
		}
	}
	/*++  End */
L_110:
	;
 
	/* ********** CHECK KORD FOR INTEGRATION ORDER CONTROL ******
	 *
	 *++  Code for INTEGO is active */
	if (divamc.iop19 != 0)
	{
		/*.           **** Check option for valid input */
		intchk[*nxtchk] = 19;
		ilow = divamc.iop19;
		jlim = -30;
		kgo = 3;
		goto L_430;
	}
	/*++  End */
L_120:
	;
 
	/* ********** CHECK SET UP FOR ERROR TOLERANCES *************
	 * */
	intchk[*nxtchk] = 16;
	ilow = divamc.iop16;
	jlim = -5;
	kgo = 4;
	if (divamc.iop16 != 0)
		goto L_430;
	/*.                      **** IN CURRENT CODE, IOP16=0 IS AN ERROR */
	kgo = 10;
	goto L_650;
L_150:
	;
	/*     **** CHECK KORD STORAGE ALLOCATION */
	optchk( intchk, iopt, "DIVA / KORD$E" );
	if (*nxtchk < 0)
		divamc.kord2i = -4;
 
	/*     ******** DONE CHECKING KORD STORAGE ALLOCATION *******
	 *
	 *     **** CHECK  Y  STORAGE ALLOCATION */
	intchk[2] = idimy;
	intchk[3] = divamc.ny + divamc.ny;
	*nxtchk = 4;
	divasc.nyny = divamc.ny + 1;
	optchk( intchk, iopt, "DIVA / Y$E" );
	if (*nxtchk < 0)
		divamc.kord2i = -4;
 
	/*     **** CHECK  F  STORAGE ALLOCATION */
	intchk[2] = idimf;
	intchk[3] = divasc.nte;
	*nxtchk = 4;
	if (divamc.iop16 != 0)
	{
		/*                                Error tolerance info. */
		intchk[4] = 16;
		intchk[5] = divamc.ntolf;
		intchk[6] = ihi - divamc.iop16 + 1;
		*nxtchk = 7;
	}
	if (divamc.iop12 > 0)
	{
		intchk[*nxtchk] = 12;
		intchk[*nxtchk + 1] = divamc.iop12;
		intchk[*nxtchk + 2] = 4;
		*nxtchk += 3;
	}
	if (divamc.iop21 > 0)
	{
		intchk[*nxtchk] = 21;
		intchk[*nxtchk + 1] = divamc.iop21;
		intchk[*nxtchk + 2] = 1;
		*nxtchk += 3;
	}
 
	/*++  Code for ERRSTO is inactive
	 *      IF (IOP20 .ne. 0) then
	 *c.                                Space for saving error estimates
	 *         INTCHK(NXTCHK) = 20
	 *         INTCHK(NXTCHK) = IOP20
	 *         INTCHK(NXTCHK) = NTE
	 *         NXTCHK = NXTCHK + 3
	 *      end if
	 *++  Code for STIFF is inactive
	 *      if (IOP21 .gt. 0) then
	 *c.                               Info. for stiff equations
	 *         INTCHK(NXTCHK) = 21
	 *         INTCHK(NXTCHK+1) = IOP21
	 *         INTCHK(NXTCHK+2) = IOP21S
	 *         NXTCHK = NXTCHK + 3
	 *      end if
	 *      MAXKQD = min(MAXKQI, 6)
	 *++  End
	 *                          Set aside space for the difference tables. */
	intchk[*nxtchk] = 0;
	intchk[*nxtchk + 1] = -KDIM*divasc.nte;
	intchk[*nxtchk + 2] = 0;
	*nxtchk += 3;
	intchk[*nxtchk] = -5*divasc.nte;
	optchk( intchk, iopt, "DIVA / F$E" );
	if (*nxtchk < 0)
	{
		divamc.kord2i = -4;
	}
	else if (divamc.kord2i != -4)
	{
		for (k = *nxtchk + 1; k <= intchk[*nxtchk]; k++)
		{
			if (intchk[intchk[k]] == 0)
			{
				divasc.ndtf = intchk[intchk[k] + 1];
				divasc.numdt = min( KDIM, (intchk[intchk[k] + 2] -
				 divasc.ndtf + 1)/divasc.nte );
				divamc.maxkqi = divasc.numdt - 1;
			}
			else
			{
				/*         Take a quick return if needed space was not specified by user. */
				divamc.kord2i = -4;
			}
		}
	}
	if (divamc.iop9 + labs( divamc.iop10 ) + divamc.iop11 != 0)
	{
		/* Insure user doesn't get in trouble with F not iniitalized. */
		for (k = divasc.ndtf; k <= (divasc.ndtf + divasc.nte*divasc.numdt -
		 1); k++)
		{
			F[k] = 0.e0;
		}
	}
L_320:
	;
	if ((divamc.kord2i == -4) || (divamc.iop10 != 0))
	{
		Mact1[3] = Iopiva[2];
		mess( mact1, (char*)text1,11, iopt );
		divamc.kord1i = 24;
		Kord[1] = 24;
	}
	divamc.tmark = Tspecs[1];
	divamc.tmarkx = Tspecs[4] + Tspecs[2];
 
	/*     **** DONE WITH INITIALIZATION AND CHECKING INPUTS */
	if (divamc.iop13 + divamc.iop14 + divamc.iop15 != 0)
		return;
L_330:
	divaa( tspecs, y, f, kord, divaf, divao );
	return;
 
	/* ************ LOOP TO CHECK OPTION SPECIFICATIONS *********
	 * */
L_430:
	jl = 0;
	for (ihi = ilow; ihi <= idimk; ihi++)
	{
		j = Kord[ihi];
		switch (kgo)
		{
			case 1: goto L_460;
			case 2: goto L_480;
			case 3: goto L_490;
			case 4: goto L_490;
		}
		/*     **** CHECK ON VARIATIONAL EQUATIONS */
L_460:
		;
		/*++  Code for VAREQ is active */
		switch (IARITHIF(j - divasc.nte))
		{
			case -1: goto L_470;
			case  0: goto L_565;
			case  1: goto L_620;
		}
L_470:
		if (j == 0)
			goto L_560;
		if (j <= jl)
			goto L_620;
		/*++  End
		 *     **** Check on diagnostic output option */
L_480:
		;
		/*++  Code for OUTPUT is active
		 *.    **** CHECK IF DONE */
		if (j >= divasc.nte)
			goto L_565;
		if (j <= jl)
			goto L_620;
		goto L_550;
		/*++  End */
L_490:
		;
		/*     **** Check integration order control (KGO=3) and
		 *     **** error tolerance equation grouping (KGO=4). */
		switch (IARITHIF(j - divasc.nte))
		{
			case -1: goto L_500;
			case  0: goto L_565;
			case  1: goto L_620;
		}
L_500:
		switch (IARITHIF(j))
		{
			case -1: goto L_510;
			case  0: goto L_530;
			case  1: goto L_540;
		}
L_510:
		if ((jl <= 0) && (ihi != ilow))
			goto L_620;
		if (j < jlim)
		{
			/*                         Output an error message. */
			Idat[2] = jlim;
			Idat[3] = 0;
			goto L_630;
		}
L_520:
		jl = -jl;
		goto L_560;
L_530:
		if (kgo == 3)
			goto L_520;
		kferr = divamc.ntolf + ihi - ilow;
		if (F[kferr] == CM1)
			goto L_510;
		/*                         Set up error message, TOL must be -1. */
		Idat[1] = kferr;
		kgo = 11;
		goto L_650;
L_540:
		if (labs( jl ) >= labs( j ))
			goto L_620;
L_550:
		jl = j;
L_560:
		;
	}
L_565:
	*nxtchk += 3;
	intchk[*nxtchk - 2] = ilow;
	intchk[*nxtchk - 1] = ihi - ilow + 1;
	switch (kgo)
	{
		case 1: goto L_30;
		case 2: goto L_110;
		case 3: goto L_120;
		case 4: goto L_150;
	}
 
	/*     **** AN ERROR HAS BEEN MADE
	 *                  Error in setting up TSPECS for extra output */
L_600:
	ihi += 2;
L_610:
	ilow = -divamc.iop5;
	goto L_630;
	/*                  Error in KORD indices */
L_620:
	Idat[2] = labs( jl ) + 1;
	Idat[3] = divasc.nte;
	/*                  Set up for print of message about KORD */
L_630:
	Idat[1] = intchk[*nxtchk];
L_640:
	Idat[4] = ihi;
	Idat[5] = Kord[ihi];
 
	/* ***************** Process Errors *************************************
	 * */
L_650:
	divamc.kord2i = -4;
	Mact[4] = LTXTAF;
	if (kgo >= 8)
		Mact[4] = -1;
	Mact[6] = Mloc[kgo];
	/*--D Next line special: P=>S, X=>D */
	dmess( mact, (char*)mtxtaa,234, divamc.idat, divamc.fdat );
	if (kgo < 8)
	{
		Mact[10] = ilow;
		Mact[13] = ilow;
		Mact[15] = -min( ihi + 2, idimk );
		mess( &Mact[9], (char*)mtxtab,56, kord );
		if (kgo <= 4)
			goto L_565;
	}
	/*              5   6   7    8    9   10   11  12 */
	switch (kgo - 4)
	{
		case 1: goto L_100;
		case 2: goto L_25;
		case 3: goto L_25;
		case 4: goto L_320;
		case 5: goto L_100;
		case 6: goto L_150;
		case 7: goto L_660;
		case 8: goto L_25;
	}
L_660:
	kgo = 4;
	goto L_565;
} /* end of function */
/*   End of DIVA */
 
		/* PARAMETER translations */
#define	C0	0.e0
#define	C1	1.e0
#define	C10	10.e0
#define	C16	16.e0
#define	C1M3	1.e-3
#define	C1M5	1.e-5
#define	C1P125	1.125e0
#define	C1P3	1.3e0
#define	C2	2.e0
#define	C2P5M3	2.5e-3
#define	C4096	4096.e0
#define	C6	6.e0
#define	C8M3	8.e-3
#define	CMP75	(-.75e0)
#define	CP0625	.0625e0
#define	CP1	.1e0
#define	CP25	.25e0
#define	CP3	.3e0
#define	CP4	.4e0
#define	CP5	.5e0
#define	CP875	.875e0
#define	LOCM	(32*256)
#undef	LTXTAC
#define	LTXTAC	40
#undef	LTXTAD
#define	LTXTAD	80
#undef	LTXTAE
#define	LTXTAE	129
#undef	LTXTAF
#define	LTXTAF	189
#undef	LTXTAG
#define	LTXTAG	237
#undef	LTXTAH
#define	LTXTAH	272
#undef	LTXTAI
#define	LTXTAI	302
#undef	LTXTAJ
#define	LTXTAJ	342
#undef	LTXTAK
#define	LTXTAK	390
#undef	LTXTAL
#define	LTXTAL	455
#undef	LTXTAM
#define	LTXTAM	484
#undef	LTXTAN
#define	LTXTAN	531
#define	MLOCAC	(23 + 32*(99 + 256*LTXTAC))
#define	MLOCAD	(13 + 32*(38 + 256*LTXTAD))
#define	MLOCAE	(22 + 32*(38 + 256*LTXTAE))
#define	MLOCAF	(3 + 32*(14 + 256*LTXTAF))
#define	MLOCAG	(21 + 32*(38 + 256*LTXTAG))
#define	MLOCAH	(2 + 32*(25 + 256*LTXTAH))
#define	MLOCAI	(11 + 32*(38 + 256*LTXTAI))
#define	MLOCAJ	(12 + 32*(38 + 256*LTXTAJ))
		/* end of PARAMETER translations */
 
void /*FUNCTION*/ divaa(
double tspecs[],
double y[],
double f[],
long kord[],
void (*divaf)(double[],double[],double[],long[]),
void (*divao)(double[],double[],double[],long[]))
{
	long int i, j, j1, j2, k, l, lx;
	double xp, xp1;
	static double disadj, sigmas, tp, tp1, tp2, tp3, tps1, tps2;
	static char mtxtaa[3][187]={"DIVAA$BAt: TN=$F, KSTEP=$I, with H=$F$EA previously reported error was fatal.$EPrint points not properly ordered: TSPEC($I)=$F$EAn error tolerance of 0 requires setting special flags. $ ",
	 "$BStep size reduced too fast, doing a restart.  $BH is so small that TN + H = TN.  $BError tolerance too small.  $BStep size at$ end of start < HMIN=$F, $BError estimates require a steps",
	 "ize < HMIN=$F, $B(Estimated Error) / (Requested Error) for equation $I$ is $F.  $BTolerance $I is F($I) = $F.$BTolerance $I is F($I) * F($I) = $F * $F = $F.$B  Replacing F($I) with $F.$E"};
	static long mloc[8]={MLOCAC,MLOCAD,MLOCAE,MLOCAF,MLOCAG,MLOCAH,
	 MLOCAI,MLOCAJ};
	static long mact[17]={MEEMES,0,0,0,MENTXT,0,METEXT,MENTXT,0,METEXT,
	 MENTXT,0,METEXT,MENTXT,LTXTAN,METEXT,MERET};
	static double exr = CP1;
	static double eimino = C8M3;
	static long ldis = 0;
	double * __restrict const hh = (double*)divamc.g;
	long int * __restrict const kexit = (long*)&divamc.iop17;
	double * __restrictconst tmarka = (double*)&divamc.tmark;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const Beta = &divamc.beta[0] - 1;
	double * __restrict const F = &f[0] - 1;
	double * __restrict const Fdat = &divamc.fdat[0] - 1;
	long * __restrict const Idat = &divamc.idat[0] - 1;
	long * __restrict const Kord = &kord[0] - 1;
	long * __restrict const Mact = &mact[0] - 1;
	long * __restrict const Mloc = &mloc[0] - 1;
	double * __restrict const Rbq = &divamc.rbq[0] - 1;
	double * __restrict const Sigma = &divamc.sigma[0] - 1;
	double * __restrict const Tmarka = &tmarka[0] - 1;
	double * __restrict const Tspecs = &tspecs[0] - 1;
	double * __restrict const Xi = &divasc.xi[0] - 1;
	double * __restrict const Y = &y[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1989-02-24 DIVAA  Krogh   Big error with BETA(2)=1+epsilon -- looped
	 *>> 1988-07-27 DIVAA  Krogh   Fixed to allow restart on a restart.
	 *>> 1988-03-07 DIVAA  Krogh   Initial code.
	 *
	 *  MAIN SUBROUTINE FOR VARIABLE ORDER INTEGRATION OF ORDINARY
	 *  DIFFERENTIAL EQUATIONS
	 * */
	/*--D Next line special: P=>D, X=>Q */
 
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
	/*.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS. */
 
 
 
	/*--D Next line special: P=>D, X=>Q */
 
	/*                      Declarations for error message processing.
	 * */
 
	/* ********* Error message text ***************
	 *[Last 2 letters of Param. name]  [Text generating message.]
	 *AA DIVAA$B
	 *AB At: TN=$F, KSTEP=$I, with H=$F$E
	 *AC A previously reported error was fatal.$E
	 *AD Print points not properly ordered: TSPEC($I)=$F$E
	 *AE An error tolerance of 0 requires setting special flags.  $B
	 *AF Step size reduced too fast, doing a restart.  $B
	 *AG H is so small that TN + H = TN.  $B
	 *AH Error tolerance too small.  $B
	 *AI Step size at end of start < HMIN=$F, $B
	 *AJ Error estimates require a stepsize < HMIN=$F, $B
	 *AK (Estimated Error) / (Requested Error) for equation $I is $F.  $B
	 *AL Tolerance $I is F($I) = $F.$B
	 *AM Tolerance $I is F($I) * F($I) = $F * $F = $F.$B
	 *AN   Replacing F($I) with $F.$E */
 
	/*                     KORD1I   Severity   Loc. message */
 
 
	/*                      1 2 3 4       5 6       7       8 9      10 */
	/*           11 12      13      14      15      16     17
	 * */
 
	/* ************** START OF EXECUTABLE CODE ******************************
	 * */
L_660:
	switch (IARITHIF(divamc.kord2i))
	{
		case -1: goto L_670;
		case  0: goto L_1380;
		case  1: goto L_1840;
	}
L_670:
	if (divamc.kord2i == -1)
		goto L_2140;
	if (divamc.kord2i == -5)
		goto L_720;
	switch (IARITHIF(divamc.kord2i + 8))
	{
		case -1: goto L_1840;
		case  0: goto L_710;
		case  1: goto L_1840;
	}
	/*     **** SPECIAL OUTPUT CASE (EXTRAPOLATION OR GSTOP) */
L_680:
	switch (*kexit)
	{
		case 1: goto L_1220;
		case 2: goto L_1190;
		case 3: goto L_830;
		case 4: goto L_690;
		case 5: goto L_690;
		case 6: goto L_2190;
	}
	/*     **** RESET TMARK BEFORE GOING WHERE DIRECTED BY KEXIT */
L_690:
	*kexit -= 2;
	goto L_1890;
	/*     **** USER HAS REQUESTED A RESTART */
L_700:
	divamc.kord2i = -5;
	divamc.igflg = 0;
	if ((Kord[2] >= 2) && (divamc.lsc < 3))
		divamc.kord2i = -8;
	if ((divamc.kord1i <= 3) || (divamc.kord1i == 5))
		goto L_1890;
	if (divamc.kord2i == -5)
		goto L_720;
	/*                Set up for a discontinuity */
L_710:
	xp = Tspecs[1];
	if (Kord[2] == 3)
	{
		/*                Adjust for discontinuity in Y */
		j = 0;
		for (i = 1; i <= divasc.nte; i++)
		{
			if (divasc.nkdko != 0)
				divasc.kordi = Kord[divasc.nkdko + i - 1];
			k = 1;
			j += divasc.kordi;
			xp1 = Y[j];
L_711:
			Y[divasc.nyny + j - k] += xp1;
			if (k < divasc.kordi)
			{
				xp1 = Y[j - k] + (divasc.tn - xp)*xp1/(double)( k );
				k += 1;
				goto L_711;
			}
		}
	}
	else
	{
		divamc.igflg = -3;
	}
	disadj = *hh;
	xp1 = (xp - divasc.tn)/Xi[1];
	if (xp1 < CP25)
	{
		k = 1;
		if (xp1 < CMP75)
			k = 2;
		Tspecs[1] = divasc.tn - Xi[k];
		Tspecs[2] = 2.e0*(divasc.tn - Tspecs[1]);
		divain( &Tspecs[1], y, f, kord );
		for (j = 1; j <= divamc.ny; j++)
		{
			Y[divasc.nyny + j - 1] = Y[j];
		}
		/*          Move difference tables back one step */
		divasc.tn = Tspecs[1];
L_714:
		divabu( f, kord );
		if (k == 2)
		{
			divamc.ksc = max( divamc.ksc - 1, 1 );
			for (k = max( 1, divamc.ksc ); k <= (divamc.iop11 - 1); k++)
			{
				Beta[k + 1] = Beta[k]*(Xi[k]/(Xi[k + 1] - Xi[1]));
			}
			k = 0;
			goto L_714;
		}
	}
	/*      Take step to discontinuity. */
	*hh = xp - divasc.tn;
	divamc.erep = -fabs( divamc.erep );
	divamc.kis = 1000;
	ldis = 1;
	divamc.lsc = 0;
	divamc.linc = 0;
	divamc.hincc = C1P125;
	divamc.ksstrt = divamc.kstep + 2;
	goto L_1115;
	/* ********
	 * INITIALIZE FOR STARTING AN INTEGRATION
	 * ******** */
L_720:
	*hh = Tspecs[2];
	divamc.lsc = 8;
	divamc.linc = -5;
	divamc.lincd = 64;
	divasc.ldt = -4;
	divasc.kqmaxi = 0;
	divamc.eimin = CP3;
	divamc.eave = C10;
	Xi[1] = C0;
	divamc.igflg = 0;
	divamc.kis = 0;
	divamc.ksstrt = 0;
	divamc.robnd = 0.e0;
 
 
	/* GO COMPUTE INITIAL DERIVATIVES */
L_730:
	divamc.kord2i = -6;
	/*++  Code for VAREQ is active */
	divamc.ics = 0;
	divamc.icf = divamc.ne;
	/*++  End */
L_740:
	divamc.kord1i = divamc.kpred;
	goto L_1360;
	/*   RETURN AFTER COMPUTING INITIAL (OR NOISE TEST) DERIVATIVES */
L_750:
	;
	/*++  Code for VAREQ is active */
	if (divamc.iop18 == 0)
		goto L_790;
	/*.    **** SPECIAL LOGIC TO COMPUTE VARIATIONAL DERIVATIVES */
L_760:
	if (divamc.kord1i == 3)
		goto L_770;
	divamc.kord1i = 3;
	Kord[3] = 0;
L_770:
	if (divamc.icf == divasc.nte)
		goto L_790;
	if (divamc.kord2i == -6)
	{
		if (divamc.ics == 1)
			goto L_790;
	}
	divamc.ics = divamc.icf + 1;
L_780:
	Kord[3] += 1;
	divamc.icf = divamc.iop18 + Kord[3];
	divamc.icf = labs( Kord[divamc.icf] );
	switch (IARITHIF(divamc.icf))
	{
		case -1: goto L_1360;
		case  0: goto L_780;
		case  1: goto L_1360;
	}
	/*++  End */
L_790:
	divamc.ics = 1;
	divamc.icf = divamc.ne;
	/*++  Code for VAREQ is active */
	if (divamc.kord2i == 0)
		switch (ARITHIF(divamc.erep))
		{
			case -1: goto L_2220;
			case  0: goto L_2220;
			case  1: goto L_1430;
		}
	/*++  End */
	switch (IARITHIF(divamc.linc + 5))
	{
		case -1: goto L_1490;
		case  0: goto L_810;
		case  1: goto L_1490;
	}
	/* END OF SPECIAL CODE FOR INITIALIZATION AT THE START
	 * ********
	 * UPDATE VARIABLES TO PREPARE FOR NEXT STEP
	 * ******** */
L_800:
	divasc.ldt = 0;
	divamc.eimin = C2P5M3 + divamc.eimin*(C6*divamc.eave + divamc.eimax)/
	 ((divamc.eimin + C1)*(C6*divamc.eimax + divamc.eave));
L_810:
	for (j = 1; j <= divamc.ny; j++)
	{
		Y[divasc.nyny + j - 1] = Y[j];
	}
	divasc.tn = Tspecs[1];
	/* ********
	 * TEST FOR VARIOUS TYPES OF OUTPUT
	 * ********
	 *     TEST FOR END OF STEP OUTPUT (OR IF DUMP OUTPUT TO BE TESTED FOR) */
	if (divamc.iop6 == 0)
		goto L_840;
	/*   SET UP FOR END OF STEP OUTPUT */
	divamc.kord1i = 6;
	goto L_1810;
	/*     SET UP AFTER OTHER COMPUTATIONS HAVE BEEN MADE (TSPECS(1).NE.TN) */
L_830:
	divamc.kord1i = 6;
	goto L_860;
	/*     TEST FOR AN OUTPUT POINT */
L_840:
	switch (ARITHIF(*hh*(divamc.tmark - divasc.tn)))
	{
		case -1: goto L_1760;
		case  0: goto L_1760;
		case  1: goto L_850;
	}
	/*     TEST FOR TOO MANY STEPS OUTPUT */
L_850:
	;
	if (divamc.ksout > divamc.kstep)
		goto L_890;
	divamc.kord1i = 4;
L_860:
	if (Tspecs[1] == divasc.tn)
		goto L_1810;
	/*     GO INTERPOLATE VALUES AT END OF LAST STEP */
L_870:
	Tspecs[1] = divasc.tn;
	goto L_1780;
	/*     CONTINUE AFTER TOO MANY STEPS OUTPUT */
L_880:
	;
	divamc.ksout = divamc.kstep + divamc.iop4;
L_890:
	;
	/*++  Code for DUMP is active */
	if (divamc.iop9 == 0)
		goto L_920;
	/*++  Code for DUMP & STIFF is inactive
	 *      KQMXDS=KQMAXD
	 *++  Code for DUMP is active */
	divamc.kqmxis = divasc.kqmaxi;
	divamc.kis += 1;
	if (divamc.kis != 0)
		goto L_920;
	/*.   TIME TO DUMP THE SOLUTION */
L_900:
	divamc.kord1i = 9;
	/*.    SET TO UPDATE THE DIFFERENCE TABLE */
	if (divasc.ldt == 1)
		goto L_1810;
	/* Note that doing this update can lead to very small differences in the
	 * results because of round off differences. */
	divasc.ldt = -2;
	goto L_1780;
	/*++  End
	 * RETURN AFTER DUMPING THE SOLUTION */
L_910:
	;
	/*++  Code for DUMP is active */
	divamc.kis = 2;
	/*.    TEST IF SOLUTION DUMP DUE TO RESTART, END, OR
	 *.    DROP IN INTEG. ORDER */
	if (divamc.linc < 0)
		switch (IARITHIF(divamc.linc + 6))
		{
			case -1: goto L_1860;
			case  0: goto L_1750;
			case  1: goto L_1180;
		}
	/*++  End
	 * END OF TESTING FOR VARIOUS TYPES OF OUTPUT
	 *   TEST IF STEPSIZE MAY BE INCREASED OR IF TESTS SHOULD BE MADE
	 *   FOR DECREASING THE STEPSIZE (OR IF STARTING OR USER SELECTING H) */
L_920:
	divamc.kstep += 1;
	switch (IARITHIF(divamc.linc))
	{
		case -1: goto L_930;
		case  0: goto L_980;
		case  1: goto L_940;
	}
L_930:
	;
	/*     **** ONLY POSSIBLE VALUES AT THIS POINT ARE
	 *          LINC = -2 OR -5 */
	switch (IARITHIF(divamc.linc + 5))
	{
		case -1: goto L_1120;
		case  0: goto L_1110;
		case  1: goto L_1120;
	}
	/* ********
	 * ERROR ESTIMATES INDICATE STEPSIZE CAN BE INCREASED
	 * ******** */
L_940:
	divamc.hc = divamc.hincc;
	if (divamc.linc > 1)
		divamc.hc = powi(divamc.hc,divamc.linc);
	*hh *= divamc.hc;
	/*     TEST IF NEW STEPSIZE IS TOO BIG */
	if (fabs( *hh ) > divamc.hmax)
		switch (ARITHIF(divamc.hmax))
		{
			case -1: goto L_970;
			case  0: goto L_970;
			case  1: goto L_960;
		}
L_950:
	divamc.eave = divamc.eimax;
	divamc.robnd = CP3 + divamc.hincc;
	goto L_1110;
	/*     NEW STEPSIZE IS TOO BIG */
L_960:
	if (fabs( Xi[1] ) >= divamc.hmaxp9)
		goto L_970;
	*hh = sign( divamc.hmax, *hh );
	goto L_950;
	/*     RESTORE THE OLD STEPSIZE */
L_970:
	*hh = Xi[1];
	divamc.linc = 0;
	goto L_1150;
	/* END OF CODE FOR CASE WHEN ERROR ESTIMATES INDICATE STEPSIZE INCREASE
	 * ********
	 * TEST IF ESTIMATED ERRORS INDICATE STEPSIZE SHOULD BE DECREASED
	 * ******** */
L_980:
	divamc.robnd = C1P3;
	if (divamc.eimax <= divamc.eave)
		goto L_990;
	divamc.eave += CP4*(divamc.eimax - divamc.eave);
	switch (ARITHIF((divamc.eimax*divamc.emax) - C1M3))
	{
		case -1: goto L_1000;
		case  0: goto L_1010;
		case  1: goto L_1010;
	}
L_990:
	divamc.eave = divamc.eimax;
	if ((divamc.eimax*divamc.emax) >= divamc.eimin)
		goto L_1010;
L_1000:
	divamc.robnd = CP3 + (Sigma[divamc.kqmaxs]/Sigma[divamc.kqmaxs - 1]);
	goto L_1180;
	/*     TEST IF STEPSIZE SHOULD BE REDUCED */
L_1010:
	if (divamc.emax*divamc.eimax < exr*divamc.eave)
		goto L_1180;
	/* ********
	 * ERROR ESTIMATES INDICATE STEPSIZE SHOULD BE REDUCED
	 * ******** */
	divamc.hc = divamc.hdec;
	if (divamc.eimin <= eimino)
		goto L_1030;
	divamc.eimin = eimino;
	divamc.hc = CP875;
L_1030:
	*hh = divamc.hc*Xi[1];
	switch (IARITHIF(divamc.lsc - 1))
	{
		case -1: goto L_1040;
		case  0: goto L_1080;
		case  1: goto L_1090;
	}
L_1040:
	if (fabs( *hh ) >= divamc.hmin)
		goto L_1090;
	if (fabs( CP875*Xi[1] ) <= divamc.hmin)
		switch (IARITHIF(divamc.linc))
		{
			case -1: goto L_1050;
			case  0: goto L_970;
			case  1: goto L_970;
		}
	*hh = sign( divamc.hmin, *hh );
	goto L_1090;
	/*     STEPSIZE IS TOO SMALL TO BE REDUCED
	 *     SET UP ERROR INDICATORS AND PREPARE FOR RETURN TO USER */
L_1050:
	divamc.kord1i = 8;
	goto L_2240;
	/*     PROCEED WITH CURRENT STEPSIZE DESPITE ERROR BEING TOO BIG */
L_1070:
	*hh = Xi[1];
	Tspecs[2] = *hh;
	divamc.emax = C0;
	divamc.linc = 0;
	switch (IARITHIF(divamc.kord1i - 2))
	{
		case -1: goto L_1420;
		case  0: goto L_1430;
		case  1: goto L_1430;
	}
	/*     SET LSC TO END STARTING PHASE */
L_1080:
	divamc.lsc = 2;
	/*     CHECK IF REPEATING A STEP */
L_1090:
	if (divamc.linc != -1)
		goto L_1110;
	/*   WHEN REPEATING A STEP, BACK UP THE DIFFERENCES AND STEPSIZE INFO. */
L_1100:
	divabu( f, kord );
	/*   TEST IF NOISE TEST (LINC = -7) OR IF H IS NOT
	 *     BEING CHANGED (LINC = -4) */
	switch (IARITHIF(divamc.linc + 4))
	{
		case -1: goto L_1780;
		case  0: goto L_1180;
		case  1: goto L_1110;
	}
	/* ********
	 * STEPSIZE IS TO BE CHANGED
	 * ******** */
L_1110:
	;
	/* MODIFY STEPSIZE TO REDUCE ROUNDOFF ERROR IN ACCUMULATING INDEP. VAR. */
	tp = C2*fabs( divasc.tn ) + C4096*fabs( *hh );
	tp = (tp + fabs( *hh )) - tp;
	if (tp != C0)
		*hh = sign( tp, *hh );
	/*     TEST IF NEW STEPSIZE SELECTED ACTUALLY GIVES A CHANGE */
	if (*hh == Tspecs[2])
		goto L_1140;
L_1115:
	Tspecs[2] = *hh;
	if (divamc.iop8 == 0)
		goto L_1140;
	/*     SETUP TO TELL USER ABOUT STEPSIZE CHANGE (OR TO CHANGE STEPSIZE) */
L_1120:
	divamc.kord1i = 8;
	goto L_860;
	/*     RETURN AFTER TELLING USER ABOUT STEPSIZE CHANGE */
L_1130:
	*hh = Tspecs[2];
L_1140:
	if (*hh != Xi[1])
		divamc.kqicon = -1;
L_1150:
	divamc.hc = fmin( divaev.eovep2, fabs( *hh ) )/divaev.eeps2;
	/* ********
	 * PREPARE FOR BEGINNING A NEW STEP
	 * ******** */
	if (divamc.linc > 0)
	{
		divamc.linc = min( divamc.linc, divamc.lincq ) + divamc.lincq;
		goto L_1190;
	}
L_1180:
	divamc.linc = divamc.lincd;
L_1190:
	if (divamc.hc > fabs( divasc.tn ))
		goto L_1200;
	/*     **** GIVE SINGULARITY DIAGNOSTIC */
	divamc.kord1i = 5;
	goto L_2240;
L_1200:
	Tspecs[1] = divasc.tn + *hh;
	if (divamc.lex == 0)
		goto L_1250;
	if (*hh*(Tspecs[1] - divamc.tmarkx) < C0)
		goto L_1250;
	Tspecs[1] = divamc.tmarkx;
	*hh = divamc.tmarkx - divasc.tn;
	divamc.linc = 64;
	if (divamc.lex > 0)
		goto L_1240;
	if ((divamc.lsc < 4) && (*hh/Xi[1] < CP3))
		goto L_1230;
L_1220:
	*hh *= CP875;
	goto L_1110;
	/*     **** GIVE OUTPUT AT CURRENT TMARK (WITH EXTRAPOLATION) */
L_1230:
	divamc.kord1i = -divamc.kmark;
	goto L_1770;
	/*     **** INTEGRATE TO TMARKX */
L_1240:
	divamc.kqicon = -1;
	/*   TEST IF SUBROUTINE FOR COMPUTING INTEGRATION COEFF. SHOULD BE CALLED */
L_1250:
	;
	/*++  Code for STIFF is inactive
	 *      IF ((KQMAXI .LT.KQICON) .OR. (KQMAXD.LT.KQDCON)) GO TO 1320
	 *++  Code for ~STIFF is active */
	if (divasc.kqmaxi < divamc.kqicon)
		goto L_1320;
	/*++  End
	 *   GO COMPUTE COEFFICIENTS REQUIRED FOR THE INTEGRATION
	 *     TEST IF STARTING */
	if (divamc.lsc < 7)
		goto L_1310;
L_1260:
	divasc.kqmaxi = 2;
	/*++  Code for STIFF is inactive
	 *      IF (METHOD) 1262,1270,1264
	 * 1262 KQMAXI=0
	 * 1264 KQMAXD=max(MAXDIF,2)
	 *      CALL DIVAHC
	 *c.  SET UP TO GO DO INITIALIZATION FOR CASE OF STIFF EQUATIONS
	 *      KORD1I=5
	 *      GO TO 1350
	 *++  End
	 *   INITIALIZE FOR EQUATIONS WHICH ARE NOT STIFF */
L_1270:
	divasc.kqmaxd = 0;
	divahc();
	j = divasc.ndtf;
	for (i = 1; i <= divasc.nte; i++)
	{
		/*++  Code for STIFF is inactive
		 *         if (KORD(I + 3) .le. 0) go to 1290
		 *++  End */
		Kord[i + 3] = 1;
		/*     INITIALIZE THE DIFFERENCE TABLE */
		if (divasc.ldt == -4)
			F[j] = F[i];
		F[j + 1] = C0;
		F[j + 2] = C0;
L_1290:
		;
		j += divasc.numdt;
	}
	if (divamc.lsc == 5)
		goto L_1340;
	divamc.lsc = 7;
	divasc.ldt = 1;
	goto L_1330;
	/*   INTEGRATION IS NOT BEING STARTED */
L_1310:
	k = Kord[divamc.kemax + 3];
	sigmas = Sigma[k];
	divahc();
	/*     **** ADJUST EAVE */
	tps1 = Beta[k];
	if (tps1 > C1)
		tps1 = CP5*tps1 + CP5;
	divamc.eave *= tps1*(Sigma[k]/sigmas);
	/*     TEST BELOW USED TO GET SAME RESULTS WITH/WITHOUT EXTRA EQUATIONS */
	if (k > divamc.kqicon)
		divamc.lsc = max( divamc.lsc, -3 );
	/* END OF SPECIAL LOGIC FOR CASE WHEN INTEG. COEFF. ROUTINE IS CALLED */
L_1320:
	;
	/* ********
	 * PREDICT Y
	 * ******** */
L_1330:
	;
	/*++  Code for ~ARGM is active */
	divapr( y, &Y[divasc.nyny], f, kord );
	/*++  Code for ARGM is inactive
	 *      CALL DIVAPE
	 *++  End
	 *     GO GET PREDICTED DERIVATIVES */
L_1340:
	divamc.kord1i = divamc.kpred;
	/* ********
	 * CALL DIVAF  (OR RETURN)
	 * ******** */
L_1350:
	divamc.kord2i = 0;
L_1360:
	Kord[1] = divamc.kord1i;
	Kord[2] = 0;
	if (divamc.iop13 != 0)
		return;
	(*divaf)( &Tspecs[1], y, f, &Kord[1] );
	/*     TEST FOR SPECIAL USER RETURN */
L_1380:
	if (Kord[1] < 0)
		goto L_2130;
	/*     TEST FOR SPECIAL CASE */
	if (divamc.kord2i != 0)
		goto L_660;
	/* ********
	 * TRANSFER CONTROL TO PROPER PLACE AFTER COMPUTING DERIVATIVES
	 * ******** */
	switch (IARITHIF(divamc.kord1i - 2))
	{
		case -1: goto L_1400;
		case  0: goto L_800;
		case  1: goto L_1390;
	}
L_1390:
	;
	/*++  Code for VAREQ is active */
	switch (IARITHIF(divamc.ics - divamc.icf))
	{
		case -1: goto L_1410;
		case  0: goto L_1410;
		case  1: goto L_760;
	}
	/*++  End
	 * ********
	 * PREPARE FOR CORRECTING, AND CORRECT Y
	 * ******** */
L_1400:
	divamc.itolep = 0;
	divamc.ilgrep = 0;
	divamc.iy = 1;
	divamc.eimax = C1M5;
	divamc.emax = C0;
	divasc.kqmaxi = 2;
	divamc.kqmaxs = 2;
	/*++  Code for STIFF is inactive
	 *      IF (METHOD) 1404,1410,1406
	 * 1404 KQMAXI=0
	 * 1406 KQMAXD=2
	 *++  End */
L_1410:
	;
	/*++  Code for ~ARGM is active */
	divacr( y, f, kord, &F[divamc.ntolf], &Kord[divamc.iop16] );
	/*++  Code for ARGM is inactive
	 *      CALL DIVACE
	 *++  End
	 *     TEST IF ESTIMATED ERROR IS TOO BIG (OR IF DIAGNOSTIC CALLED FOR) */
	if (divamc.emax > divamc.erep)
		switch (ARITHIF(divamc.erep))
		{
			case -1: goto L_2210;
			case  0: goto L_2210;
			case  1: goto L_1670;
		}
L_1420:
	;
	/*++  Code for VAREQ is active */
	if (divamc.iop18 != 0)
		goto L_760;
	/*++  End */
L_1430:
	divamc.kord1i = 2;
	/*     TEST IF NOISE APPEARS TO LIMIT PRECISION */
	if (divamc.emax < C0)
		goto L_1470;
	/*++  Code for ~STIFF is active */
	switch (IARITHIF(divamc.lsc))
	{
		case -1: goto L_1450;
		case  0: goto L_1360;
		case  1: goto L_1610;
	}
	/*++  Code for STIFF is inactive
	 *      IF (LSC) 1450,1460,1610
	 *++  End
	 *     SET LSC=0 IF NOISE NO LONGER APPEARS TO LIMIT PRECISION
	 *     OR IF THE END OF THE STARTING PHASE HAS BEEN REACHED */
L_1450:
	divamc.lsc = 0;
L_1460:
	switch (IARITHIF(divamc.method))
	{
		case -1: goto L_800;
		case  0: goto L_1350;
		case  1: goto L_1350;
	}
	/* ********
	 * NOISE APPEARS TO BE LIMITING PRECISION
	 * ******** */
L_1470:
	;
	if (divamc.lsc <= 0)
		divamc.lsc = max( divamc.lsc - 1, -divamc.kqmaxs );
	if (divamc.lsc == -1)
		goto L_1460;
	if (fabs( divamc.emax ) < exr)
		goto L_1590;
	divamc.linc = -7;
	tps2 = powi(C1 + Beta[divamc.noiseq - 1],divamc.noiseq);
	if (divamc.snoise < divaev.eeps10*tps2)
		goto L_1550;
	tp = sign( divaev.eept75*fabs( divasc.tn ) + divaev.ovtm75, *hh );
	if (fabs( tp ) > fabs( *hh ))
		goto L_1550;
	Tspecs[1] = divasc.tn + tp;
	divamc.kord1i = 0;
	/*     **** GO TO BACK UP THE DIFFERENCES AND GET F(TSPECS(1)) */
	goto L_1100;
	/*     **** SOLUTION HAS BEEN INTERPOLATED AND F COMPUTED */
L_1490:
	;
	divamc.kord1i = 0;
	divamc.linc -= 1;
	switch (IARITHIF(divamc.linc + 9))
	{
		case -1: goto L_1510;
		case  0: goto L_1520;
		case  1: goto L_1500;
	}
L_1500:
	Tspecs[1] = divasc.tn + (tp + tp);
	tp1 = F[divamc.kemax];
	tp2 = F[divasc.ndtf + divasc.numdt*divamc.kemax - divasc.numdt];
	goto L_1780;
	/*     **** COMPUTE 2-ND DIFFERENCE AT CLOSE SPACED T VALUES */
L_1510:
	tp2 = tp3;
L_1520:
	tp3 = F[divamc.kemax];
	tps1 = fabs( (tp3 - tp1) - (tp1 - tp2) );
	if ((C16*tps1*tps2) >= divamc.dnoise)
		switch (IARITHIF(divamc.linc + 9))
		{
			case -1: goto L_1550;
			case  0: goto L_870;
			case  1: goto L_1550;
		}
L_1530:
	;
	tps2 = CP25*divamc.snoise/Rbq[divamc.noiseq];
	for (k = 2; k <= divasc.numdt; k++)
	{
		tps1 += tps1;
		Rbq[k] = fmax( tps1, tps2*Rbq[k] );
	}
	divamc.linc = 0;
 
	/*FTK Next two lines added 2009-10-15 */
	if (fabs( divamc.emax ) < divamc.erep)
		goto L_1460;
	/*FTK  LINC = -1  And then on 2015-03-14 commented out this line */
 
	*hh *= CP875;
	goto L_1040;
	/*     **** SET UP TO GIVE NOISE DIAGNOSTIC */
L_1550:
	divamc.kord1i = 6;
	goto L_2240;
	/*     **** AFTER GIVING NOISE DIAGNOSTIC */
L_1560:
	divamc.kord1i = 2;
	if (Kord[2] >= 0)
	{
		tps1 = divaev.eeps10;
		/*FTK Next line added 2009-10-15 */
		if (tps1 < .49e0*Rbq[2])
			goto L_1530;
	}
	/*     **** SET NEW VALUE FOR OFFENDING TOL */
	F[divamc.ntolf + divamc.itolep - 1] = Fdat[7];
	switch (IARITHIF(divamc.linc + 7))
	{
		case -1: goto L_1180;
		case  0: goto L_1570;
		case  1: goto L_1180;
	}
L_1570:
	divamc.linc = 0;
L_1580:
	switch (IARITHIF(divamc.lsc))
	{
		case -1: goto L_1460;
		case  0: goto L_1460;
		case  1: goto L_1610;
	}
	/*     **** CHANGE HINCC AND ADJUST SIGMA( ) */
L_1590:
	if (divamc.lsc != -4)
		goto L_1580;
	if (divamc.hincc == C1P125)
		goto L_1580;
	tps1 = C1P125/divamc.hincc;
	tps2 = 1.0e0;
	for (k = 2; k <= divamc.iop11; k++)
	{
		tps2 *= tps1;
		Sigma[k] *= tps2;
	}
	divamc.eave *= powi(tps1,1 - Kord[divamc.kemax + 3]);
	divamc.lincd = 6;
	divamc.lincq = 12;
	divamc.hincc = C1P125;
	goto L_1460;
	/*   END OF CODE FOR CASE WHEN NOISE APPEARS TO LIMIT PRECISION
	 * ********
	 * SPECIAL LOGIC FOR STARTING THE INTEGRATION
	 * ******** */
L_1610:
	if (divamc.lsc == 1)
		goto L_800;
	divamc.lsc -= 1;
	switch (IARITHIF(divamc.lsc - 2))
	{
		case -1: goto L_1620;
		case  0: goto L_1640;
		case  1: goto L_1650;
	}
L_1620:
	if (divamc.eimax <= (CP0625*divamc.eave*(Sigma[divamc.kqmaxs]/
	 sigmas)*powi(Beta[divamc.kqmaxs + 1],2)))
		goto L_800;
L_1630:
	divamc.ksstrt = divamc.kstep + 2;
	/*   TEST IF STEPSIZE IS TOO SMALL BEFORE ENDING STARTING PHASE */
	if (fabs( *hh ) >= divamc.hmin)
		goto L_1450;
	/*     GIVE DIAGNOSTIC FOR STEPSIZE TOO SMALL AT END OF START */
	divamc.kord1i = 7;
	goto L_2240;
	/*   SET LSC TO DO ONE DERIVATIVE EVAL. PER STEP */
L_1640:
	divamc.lsc = 1;
	goto L_800;
	/*     TEST IF FIRST TIME THROUGH THE FIRST STEP */
L_1650:
	if (divamc.lsc == 6)
		goto L_1340;
	/*     END STARTING PHASE IF CONVERGENCE OF CORRECTOR ITERATES TOO SLOW */
	if (divasc.ldt == -5)
		goto L_1660;
	divamc.lsc = min( divamc.lsc, 4 );
	goto L_800;
L_1660:
	divasc.ldt = 0;
	switch (IARITHIF(divamc.lsc - 4))
	{
		case -1: goto L_1260;
		case  0: goto L_1630;
		case  1: goto L_1260;
	}
	/* END OF SPECIAL LOGIC FOR STARTING THE INTEGRATION
	 * ********
	 * ESTIMATED ERROR IS TOO BIG
	 * ******** */
L_1670:
	switch (ARITHIF(Beta[2] - C1))
	{
		case -1: goto L_1690;
		case  0: goto L_1730;
		case  1: goto L_1680;
	}
L_1680:
	divamc.hc = C1/Beta[2];
	if (Beta[2] >= C1P125)
		goto L_1740;
L_1690:
	if (Beta[2] > CP1)
		goto L_1730;
	/*   REQUIRED STEPSIZE REDUCTION IS TOO RAPID -- GIVE A DIAGNOSTIC */
	divamc.kord1i = 4;
	goto L_2240;
 
	/*     TEST KORD(2) AFTER ABOVE DIAGNOSTIC OR A DISCONTINUITY DIAGNOSTIC */
L_1700:
	;
	if (Kord[2] == 0)
		goto L_1730;
	/*  TEST IF SOLUTION MUST BE DUMPED BEFORE A RESTART */
L_1710:
	divamc.linc = -1;
	/*++  Code for DUMP is active */
	if (divamc.iop9 == 0)
		goto L_1750;
	if (divamc.kis == 2)
		goto L_1750;
	divamc.linc = -6;
	/*.    GO DUMP SOLUTION BEFORE REPEATING THE STEP */
L_1720:
	divasc.kqmaxi = divamc.kqmxis;
	/*++  Code for DUMP & STIFF is inactive
	 *      KQMAXD=KQMXDS
	 *++  Code for DUMP is active */
	divabu( f, kord );
	goto L_900;
	/*++  End
	 *   SET UP TO REPEAT THE STEP */
L_1730:
	divamc.hc = CP5;
L_1740:
	divamc.linc = -1;
	if (divamc.lsc <= 3)
		goto L_1030;
	/*   RESTART THE INTEGRATION IF ERROR IS TOO BIG ON FIRST OR SECOND STEP
	 * LOOP TO SELECT A NEW INITIAL STEPSIZE */
L_1750:
	divamc.lsc = 7;
L_1755:
	*hh *= CP5;
	divamc.emax *= CP25;
	if (divamc.emax >= CP3)
		goto L_1755;
	goto L_1090;
	/*   END OF SELECTING A NEW INITIAL STEPSIZE
	 * END OF LOGIC FOR CASE WHEN ESTIMATED ERROR IS TOO BIG
	 * ********
	 * INTEGRATION HAS REACHED AN OUTPUT POINT
	 * ******** */
L_1760:
	if (divamc.kmark == 0)
		goto L_1920;
	divamc.kord1i = min( divamc.kmark, 5 );
	Kord[3] = divamc.kmark;
	if (Tspecs[1] == divamc.tmark)
		goto L_1790;
L_1770:
	Tspecs[1] = divamc.tmark;
L_1780:
	divain( &Tspecs[1], y, f, kord );
L_1790:
	;
	switch (IARITHIF(divamc.kord1i))
	{
		case -1: goto L_1800;
		case  0: goto L_730;
		case  1: goto L_1810;
	}
	/*   OUTPUT POINT IS OBTAINED BY EXTRAPOLATION */
L_1800:
	;
	/*++  Code for EXTRAP is active */
	divamc.kord1i = -divamc.kord1i;
	divamc.kord2i = -7;
	*kexit = 4;
	/*.  TEST IF GSTOP-S ARE PRESENT
	 *++  Code for EXTRAP &  GSTOP is active */
	if (divamc.ngtot == 0)
		goto L_1820;
	divamc.igflg = 4;
	*kexit = 2;
	divamc.kord1i = 7;
	switch (IARITHIF(divamc.iop7))
	{
		case -1: goto L_740;
		case  0: goto L_1820;
		case  1: goto L_740;
	}
	/*++  End
	 * ********
	 * CALL DIVAO  (OR RETURN)
	 * ******** */
L_1810:
	divamc.kord2i = 1;
L_1820:
	Kord[1] = divamc.kord1i;
	Kord[2] = 1;
	if (divamc.iop14 != 0)
		return;
	(*divao)( &Tspecs[1], y, f, &Kord[1] );
	/*     TEST FOR SPECIAL USER RETURN OR FOR A RESTART
	 *++  Code for ~DUMP is inactive
	 * 1840 IF (KORD(1)) 2130,700,1880
	 *++  Code for DUMP is active */
L_1840:
	if (Kord[1] > 0)
		goto L_1880;
L_1850:
	if (divamc.iop9 == 0)
		goto L_1870;
	/*.    **** GO DUMP THE SOLUTION */
	divamc.linc = -7;
	divamc.itolep = Kord[1];
	Idat[1] = Kord[2];
	divamc.neptol = divamc.kord1i;
	if (divamc.lsc != 8)
		goto L_900;
L_1860:
	divamc.linc = min( 0, divamc.lincd );
	divamc.kord1i = divamc.neptol;
	Kord[1] = divamc.itolep;
	Kord[2] = Idat[1];
L_1870:
	switch (IARITHIF(Kord[1]))
	{
		case -1: goto L_2130;
		case  0: goto L_700;
		case  1: goto L_2100;
	}
	/*++  End */
L_1880:
	if (divamc.kord2i < 0)
		switch (-divamc.kord2i)
		{
			case 1: goto L_2140;
			case 2: goto L_1810;
			case 3: goto L_1350;
			case 4: goto L_2110;
			case 5: goto L_720;
			case 6: goto L_750;
			case 7: goto L_680;
			case 8: goto L_710;
		}
	if (divamc.kord2i == 0)
		goto L_1380;
	/* ********
	 * TRANSFER CONTROL TO PROPER PLACE AFTER OUTPUT
	 * ******** */
L_1890:
	switch (IARITHIF(divamc.kord1i - 5))
	{
		case -1: goto L_1910;
		case  0: goto L_1930;
		case  1: goto L_1900;
	}
L_1900:
	switch (IARITHIF(divamc.kord1i - 8))
	{
		case -1: goto L_840;
		case  0: goto L_1130;
		case  1: goto L_910;
	}
L_1910:
	switch (IARITHIF(divamc.kord1i - 3))
	{
		case -1: goto L_1920;
		case  0: goto L_1930;
		case  1: goto L_880;
	}
	/*   GET NEW TOUT */
L_1920:
	divamc.tout = Tspecs[1] + Tspecs[3];
	/* GET NEW TMARK (NEXT INDEP. VAR. OUTPUT POINT) */
L_1930:
	xp = divamc.tmark;
	k = divamc.kmark;
	divamc.tmark = divamc.tout;
	divamc.kmark = 2;
	divamc.lex = 0;
	switch (IARITHIF(divamc.iop5))
	{
		case -1: goto L_1940;
		case  0: goto L_1980;
		case  1: goto L_1970;
	}
L_1940:
	i = -divamc.iop5;
L_1950:
	i += 3;
	j1 = Kord[i - 3];
	switch (IARITHIF(j1))
	{
		case -1: goto L_1950;
		case  0: goto L_1980;
		case  1: goto L_1960;
	}
L_1960:
	j2 = Kord[i - 2];
	l = Kord[i - 1];
	goto L_1990;
L_1970:
	j1 = 5;
	j2 = divamc.iop5;
	l = 0;
	if (j2 >= j1)
		goto L_1990;
L_1980:
	j1 = 4;
	j2 = 4;
	l = divamc.iop3;
 
	/*     **** LOOP TO SET NEW TMARK (AND TMARKX) */
L_1990:
	for (j = j1; j <= j2; j++)
	{
		/*        **** TEST IF EXTRAPOLATION NOT POSSIBLE */
		if (l == 0)
			goto L_2010;
		lx = 2;
		switch (IARITHIF(divamc.lex))
		{
			case -1: goto L_2020;
			case  0: goto L_2030;
			case  1: goto L_2020;
		}
L_2000:
		divamc.lex = l;
L_2010:
		lx = 1;
L_2020:
		switch (ARITHIF(*hh*(Tspecs[j] - Tmarka[lx])))
		{
			case -1: goto L_2030;
			case  0: goto L_2060;
			case  1: goto L_2060;
		}
L_2030:
		if (j == 4)
			goto L_2050;
		switch (ARITHIF(*hh*(Tspecs[j] - xp)))
		{
			case -1: goto L_2060;
			case  0: goto L_2040;
			case  1: goto L_2050;
		}
L_2040:
		if ((k >= j) || (k == 3))
			goto L_2060;
L_2050:
		Tmarka[lx] = Tspecs[j];
		if (lx == 2)
			goto L_2000;
		divamc.kmark = j;
L_2060:
		;
	}
	if (divamc.iop5 < 0)
		goto L_1950;
	if (j1 != 4)
		goto L_1980;
	if (divamc.kmark == 4)
		divamc.kmark = 3;
	/*     **** TEST IF NEW TMARK IS ACCEPTABLE */
	switch (ARITHIF(*hh*(xp - divamc.tmark)))
	{
		case -1: goto L_2070;
		case  0: goto L_2080;
		case  1: goto L_2090;
	}
L_2070:
	switch (IARITHIF(divamc.kord2i - 1))
	{
		case -1: goto L_670;
		case  0: goto L_840;
		case  1: goto L_670;
	}
L_2080:
	if (k != divamc.kmark)
		goto L_2070;
	/*++  Code for DUMP is active */
	if (divamc.kord1i == 3)
		goto L_1850;
	/*++  Code for ~DUMP is inactive
	 *      IF (KORD1I .EQ. 3) GO TO 2100
	 *++  End */
L_2090:
	if (divamc.kord1i == 13)
		goto L_2190;
	/* SETUP TO INDICATE ERROR IN SPECIFICATION OF OUTPUT POINTS */
	divamc.kord1i = 2;
	Idat[2] = divamc.kmark;
	if (divamc.kmark <= 3)
		Idat[2] = divamc.kmark + 1;
	Fdat[3] = Tspecs[Idat[2]];
	goto L_2240;
	/*     SET KORD1I=1 TO INDICATE THAT END OF INTEGRATION HAS BEEN REACHED */
L_2100:
	divamc.kord1i = 1;
	/* ********
	 * RETURN TO USER
	 * ******** */
L_2110:
	divamc.kord2i = -1;
	Kord[1] = divamc.kord1i;
L_2130:
	Kord[2] = -1;
	return;
	/* ********
	 * TRANSFER CONTROL TO PROPER PLACE AFTER RETURN TO USER
	 * ******** */
L_2140:
	switch (IARITHIF(divamc.kord1i - 2))
	{
		case -1: goto L_2150;
		case  0: goto L_1560;
		case  1: goto L_2160;
	}
L_2150:
	divamc.kord2i = 1;
	goto L_1930;
L_2160:
	switch (IARITHIF(divamc.kord1i - 4))
	{
		case -1: goto L_1700;
		case  0: goto L_2200;
		case  1: goto L_2170;
	}
L_2170:
	switch (IARITHIF(divamc.kord1i - 13))
	{
		case -1: goto L_2180;
		case  0: goto L_1930;
		case  1: goto L_2190;
	}
L_2180:
	if (fabs( *hh ) >= divamc.hmin)
		switch (IARITHIF(divamc.kord1i - 11))
		{
			case -1: goto L_1030;
			case  0: goto L_1450;
			case  1: goto L_1030;
		}
	if (Kord[2] == 0)
		switch (IARITHIF(divamc.kord1i - 11))
		{
			case -1: goto L_1070;
			case  0: goto L_800;
			case  1: goto L_1070;
		}
	/*   ERROR MESSAGES HAVE BEEN IGNORED -- COMPUTATION CAN NOT CONTINUE */
L_2190:
	divamc.kord1i = 1;
	goto L_2240;
 
	/*        AFTER A DISCONTINUITY RETURN */
L_2200:
	divamc.linc = -4;
	switch (IARITHIF(Kord[2]))
	{
		case -1: goto L_1710;
		case  0: goto L_1730;
		case  1: goto L_1100;
	}
	/* ********
	 * PROBLEM ENCOUNTERED WHEN CORRECTING
	 * ******** */
L_2210:
	if (ldis == 0)
		goto L_2230;
	/*           Extra checks when had a user specified discontinuity.
	 *++  Code for VAREQ is active */
	if (divamc.iop18 != 0)
		goto L_760;
	/*++  End */
L_2220:
	divamc.kord1i = 2;
	ldis += 1;
	tp = disadj/ *hh;
	if (divamc.kis >= 1000)
	{
		if (ldis == 2)
		{
			if (divamc.kqmaxs <= 3)
			{
				ldis = 0;
				divamc.erep = fabs( divamc.erep );
				Tspecs[2] = *hh*fmin( fmin( tp, SQ(tp) ), pow(CP25*
				 exr/divamc.emax,.333333333e0) );
				goto L_720;
			}
			divamc.linc = -5;
			if (divamc.iop9 == 0)
				divamc.kis = 1001;
			goto L_800;
		}
		if (divamc.iop9 == 0)
			divamc.kis += 1;
		if (divamc.kqmaxs <= ldis + 2)
			divamc.kis = ldis + 1;
		divamc.linc = min( divamc.linc, ldis - 2 );
	}
	if (ldis > 2*divamc.kqmaxs)
	{
		divamc.erep = fabs( divamc.erep );
		ldis = 0;
		if (divamc.emax > divamc.erep)
			goto L_1670;
		goto L_1430;
	}
	if (tp >= powi(divamc.hincc,divamc.linc + 2))
	{
		if ((ldis != 3) && (tp > (double)( divamc.kqmaxs )))
			divamc.lsc = 1;
		divamc.eimin = CP5;
		divamc.eave *= powi(tp,8);
	}
	if (divamc.lsc == 2)
		goto L_1630;
	if (divamc.emax > exr)
		goto L_1730;
	goto L_1430;
 
L_2230:
	divamc.erep = fabs( divamc.erep );
	/*++  Code for DUMP is active */
	if (divamc.linc < -3)
		goto L_1720;
	/*++  End
	 *     BAD TOL */
	divamc.kord1i = 3;
	/* ********
	 * ERROR PROCESSING
	 * ******** */
L_2240:
	Fdat[1] = divasc.tn;
	Fdat[2] = *hh;
	Idat[1] = divamc.kstep;
	divamc.itolep = max( divamc.neptol, -divamc.neptol - 1 );
	j = 3;
	if (divamc.kord1i >= 7)
	{
		j = 4;
		Fdat[3] = divamc.hmin;
	}
	if (divamc.kord1i <= 3)
	{
		if (divamc.kord1i < 3)
		{
			k = 8;
		}
		else
		{
			Mact[9] = LTXTAL;
			Fdat[3] = C0;
			Idat[2] = divamc.itolep;
			Idat[3] = divamc.itolep + divamc.ntolf - 1;
			k = 11;
		}
	}
	else
	{
		Mact[9] = LTXTAK;
		Fdat[j] = divamc.emax;
		Idat[2] = divamc.kemax;
		Idat[3] = divamc.itolep;
		Idat[4] = divamc.itolep + divamc.ntolf - 1;
		Fdat[j + 1] = F[Idat[4]];
		k = 14;
		if (divamc.kord1i == 6)
		{
			k = 17;
			Idat[5] = Idat[4];
			Fdat[7] = 32.e0*fabs( divamc.emax )*Fdat[j + 1];
			Fdat[j + 2] = Fdat[7];
		}
		Mact[12] = LTXTAL;
		if (divamc.neptol < 0)
		{
			Mact[12] = LTXTAM;
			Idat[6] = Idat[4];
			Idat[5] = Idat[4] + 1;
			Fdat[j + 2] = F[Idat[5]];
			Fdat[j + 3] = Fdat[j + 1]*Fdat[j + 2];
		}
	}
	/* Set the location for the first part of the message that varies, set
	 * the error severity, and the index number, print the error and
	 * return or stop. */
	l = Mloc[divamc.kord1i];
	Mact[6] = l/LOCM;
	Mact[2] = (l - Mact[6]*LOCM)/32;
	divamc.kord1i = l%32;
	Mact[3] = divamc.kord1i;
	Mact[k] = MERET;
	/*--D Next line special: P=>S, X=>D */
	dmess( mact, (char*)mtxtaa,187, divamc.idat, divamc.fdat );
	Mact[k] = MENTXT;
	goto L_2110;
 
} /* end of function */
/*   End of DIVAA */
 
void /*FUNCTION*/ divabu(
double f[],
long kord[])
{
	long int i, j, k, kqq, l;
	double tpd;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const Beta = &divamc.beta[0] - 1;
	double * __restrict const F = &f[0] - 1;
	long * __restrict const Kord = &kord[0] - 1;
	double * __restrict const Xi = &divasc.xi[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1987-12-07 DIVABU Krogh   Initial code.
	 *
	 * THIS SUBROUTINE RESTORES THE DIFFERENCE TABLE TO ITS STATE
	 * AT THE BEGINNING OF THE CURRENT STEP.  IF THE INTEGRATION ORDER
	 * WAS INCREASED, IT IS REDUCED. THE COMMON ARRAY XI IS ALSO
	 * RESTORED TO ITS STATE AT THE BEGINNING OF THE STEP. IF THE
	 * STEPSIZE IS NOT BEING CHANGED, THE ARRAY V USED TO COMPUTE
	 * INTEGRATION COEFFICIENTS IS RESTORED.
	 * */
 
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
 
 
	/* ********* START OF EXECUTABLE CODE **********
	 *
	 * ********
	 * BEGIN LOOP TO BACK UP DIFFERENCE TABLES
	 * ******** */
	l = divasc.ndtf - 1;
	for (i = 1; i <= divasc.nte; i++)
	{
		kqq = Kord[i + 3];
		/*++  Code for STIFF is inactive
		 *         IF (KQQ) 2302,2400,2310
		 *c.           EQUATION IS STIFF
		 * 2302    IF (LINC.GE.0) GO TO 2310
		 *         IF (F(L+1+I)) 2306,2308,2304
		 *c.     ORDER WAS INCREASED, AND THUS MUST BE DECREASED (KQQ.LT.0)
		 * 2304    KQQ=KQQ+1
		 *         KORD(I+3) = KQQ
		 *         GO TO 2308
		 *c.     ORDER WAS DECREASED
		 * 2306    KQQ=KQQ-1
		 * 2308    KQQ=max(2,-KQQ)
		 *         GO TO 2350
		 *++  End
		 *     EQUATION IS NOT STIFF */
L_2310:
		if (kqq > 2)
		{
			if (F[l + kqq] == C0)
			{
				/*                 ORDER WAS INCREASED, AND THUS MUST BE DECREASED */
				kqq -= 1;
				Kord[i + 3] = kqq;
			}
		}
		j = min( kqq, divamc.ksc );
		divasc.kqmaxi = max( divasc.kqmaxi, kqq );
		if (kqq != 1)
			F[l + kqq + 1] = 0.e0;
		/*           BACK UP FOR BACKWARD DIFFERENCES */
		for (k = 1; k <= j; k++)
		{
			F[l + k] -= F[l + k + 1];
		}
		if (kqq > divamc.ksc)
		{
			/*           BACK UP FOR MODIFIED DIVIDED DIFFERENCES */
			for (k = j + 1; k <= kqq; k++)
			{
				F[l + k] = (F[l + k] - F[l + k + 1])/Beta[k];
			}
		}
L_2400:
		F[l + kqq + 1] /= Beta[kqq + 1];
		l += divasc.numdt;
	}
	/* END OF LOOP TO BACK UP DIFFERENCE TABLES
	 * ********
	 * BACK UP XI TO BEGINNING OF THE STEP
	 * ******** */
	i = divamc.ksc + 1;
	switch (IARITHIF(i - divamc.iop11 - 1))
	{
		case -1: goto L_2420;
		case  0: goto L_2440;
		case  1: goto L_2450;
	}
L_2420:
	tpd = Xi[1];
	/*                Check below needed when starting? */
	if (tpd == Xi[2])
		goto L_2450;
	for (k = i; k <= divamc.iop11; k++)
	{
		Xi[k - 1] = Xi[k] - tpd;
	}
L_2440:
	Xi[divamc.iop11] = C2*Xi[divamc.iop11 - 1];
	if (divamc.iop11 != 2)
		Xi[divamc.iop11] -= Xi[divamc.iop11 - 2];
L_2450:
	divamc.kqicon = -1;
	divamc.icf = divamc.ne;
	divamc.ics = 1;
	divasc.ldt = 1;
	return;
} /* end of function */
/*   End of DIVABU */
 
void /*FUNCTION*/ divaco(
long id[],
double rd[])
{
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	long * __restrict const Id = &id[0] - 1;
	double * __restrict const Rd = &rd[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1987-12-07 DIVACO Krogh   Initial code.
	 *
	 * THIS SUBROUTINE RETURNS THE FOLLOWING DATA FROM COMMON
	 * ID(1) = KEMAX  =  INDEX OF EQUATION WITH LARGEST ERROR ESTIMATE
	 * ID(2) = KSTEP  =  CURRENT STEP NUMBER
	 * ID(3) = NUMDT  =  NUMBER OF DIFFERENCES USED FOR EACH EQUATION
	 * ID(4) =           RESERVED FOR FUTURE USE
	 * ID(5) =           RESERVED FOR FUTURE USE
	 * RD(1) = EMAX   =  MAX. RATIO OF ESTIMATED ERROR TO REQUESTED ERROR
	 * RD(2) =           RESERVED FOR FUTURE USE
	 * RD(3) =           RESERVED FOR FUTURE USE
	 * */
 
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
 
 
	Id[1] = divamc.kemax;
	Id[2] = divamc.kstep;
	Id[3] = divasc.numdt;
	Rd[1] = divamc.emax;
	return;
} /* end of function */
/*   End of DIVACO */
 
		/* PARAMETER translations */
#define	C1000	1000.e0
#define	C1P4	1.4e0
#define	C20	20.e0
#define	C4	4.e0
#define	C40	40.e0
#define	CM2	(-2.e0)
#define	CM8	(-8.e0)
#define	CMP5	(-.5e0)
#define	CP125	.125e0
#define	CP75	.75e0
#define	CP8	.8e0
#define	CP9375	.9375e0
#define	CQ3125	.03125e0
#define	METABL	55
		/* end of PARAMETER translations */
 
void /*FUNCTION*/ divacr(
double y[],
double f[],
long kord[],
double tol[],
long lgroup[])
{
	long int i, ilgror, iord, itolor, j, jlgrep, jlgror, k, kqd, kql,
	 kqlord, kqn, l, ll, _i, _r;
	static long int koutko, lkqmax;
	double e, ei, eps, ercoef, rnd, rnoise, s, tempao[4], tp2, tpp,
	 tps5, tps6, tps7;
	static double ref[4];
	static char mtxtaa[1][105]={"KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) EIMIN=$(E8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G RQ=$(E11.5)$G$E"};
	static char mtxtab[1][89]={"I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$HHIGH ORDER PREDICTED$ DIFFERENCES$HRNOISE$HSTIFF$HBETA$E"};
	static long mact1[2]={METEXT,MERET};
	static long mact2[12]={METABL,1,0,14,400201,300202,801503,1507501,
	 1103504,901501,1002501,1205501};
	static double eibnd[KDIM - 1]={.1e0,.1e0,.14e0,.19e0,.26e0,.36e0,
	 .50e0,.69e0,.94e0,C1,C1,C1,C1,C1,C1,C1,C1,C1,C1};
	static int _aini = 1;
	/* EQUIVALENCE translations */
	double _e0[4];
	double * __restrict const hh = (double*)divamc.g;
	double * __restrict const tempa = (double*)_e0;
	double * __restrict const tps1 = (double*)_e0;
	double * __restrict const tps2 = (double*)((double*)_e0 + 1);
	double * __restrict const tps3 = (double*)((double*)_e0 + 2);
	double * __restrict const tps4 = (double*)((double*)_e0 + 3);
	/* end of EQUIVALENCE translations */
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const Beta = &divamc.beta[0] - 1;
	double * __restrict const Eibnd = &eibnd[0] - 1;
	double * __restrict const F = &f[0] - 1;
	double * __restrict const Fdat = &divamc.fdat[0] - 1;
	double * __restrict const Gs = &divamc.gs[0] - 1;
	long * __restrict const Idat = &divamc.idat[0] - 1;
	long * __restrict const Kord = &kord[0] - 1;
	long * __restrict const Lgroup = &lgroup[0] - 1;
	long * __restrict const Mact1 = &mact1[0] - 1;
	long * __restrict const Mact2 = &mact2[0] - 1;
	double * __restrictconst Rbq = &divamc.rbq[0] - 1;
	double * __restrict const Ref = &ref[0] - 1;
	double * __restrict const Sigma = &divamc.sigma[0] - 1;
	double * __restrict const Tempa = &tempa[0] - 1;
	double * __restrict const Tempao = &tempao[0] - 1;
	double * __restrict const Tol = &tol[0] - 1;
	double * __restrict const Y = &y[0] - 1;
		/* end of OFFSET VECTORS */
	if( _aini ){ /* Do 1 TIME INITIALIZATIONS! */
		Ref[1] = C1;
		Ref[2] = CP9375;
		Ref[3] = CP75;
		Ref[4] = CP5;
		_aini = 0;
	}
 
	/*>> 1988-08-25 DIVACR Krogh   Fix bug in relative error test.
	 *>> 1988-01-15 DIVACR Krogh   Initial code.
	 *
	 * THIS SUBROUTINE
	 *   1. CORRECTS Y FOR EQUATIONS WHICH ARE NOT STIFF
	 *   2. ESTIMATES ERRORS
	 *   3. SELECTS INTEGRATION ORDERS
	 *   4. TESTS IF NOISE LIMITS THE PRECISION
	 *
	 *     Y = VECTOR OF PREDICTED VALUES ON ENTRY, AND OF CORRECTED
	 *         VALUES WHEN THE RETURN IS MADE.
	 * LGROUP= VECTOR INDICATING HOW ERROR TOLERANCES ARE TO BE GROUPED
	 *         (AND POSSIBLY HOW INTEGRATION ORDERS ARE TO BE GROUPED).
	 *   TOL = VECTOR CONTAINING ERROR TOLERANCES (AND POSSIBLY RELATIVE
	 *         ERROR FACTORS).
	 *     F = VECTOR GIVING PREDICTED DERIVATIVE VALUES AND DIFF. TABLES.
	 *    KD = VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS
	 *         (IF EQUATIONS HAVE DIFFERENT ORDERS).
	 *    KQ = VECTOR OF INTEGRATION ORDERS.
	 * */
	/*--D Next line special: P=>D, X=>Q */
 
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
 
	/*.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS. */
 
	/*++  Code for INTEGO is active */
	/*++  End */
	/*             Parameters for Interface to MESS and DMESS */
	/* ********* Error message text ***************
	 *[Last 2 letters of Param. name]  [Text generating message.]
	 *AA KSTEP=$(I6) T=$(E15.8) H=$(E12.5) LSC=$(I3) $C
	 *   EIMIN=$(E8.2) EAVE=$G KSC=$(I2) SIGMA($J)=$G $C
	 *   RQ=$(E11.5)$G$E
	 *   $
	 *AB I$HKQ$HLI$HE$HEI$HEPS$HF$H$H$H$H
	 *   HIGH ORDER PREDICTED DIFFERENCES$HRNOISE$HSTIFF$HBETA$E */
 
	/* (rr=repeat, t=3/5 for I/E format)  wwddtrr  wwddtrr  wwddtrr */
	/*         wwddtrr  wwddtrr  wwddtrr  wwddtrr  wwddtrr
	 *          End of stuff for interface to message processor
	 * */
	/*++ Save data by elements if ~.C.
	 *++ Of next 20 lines, only the first KDIM-1 are active */
	/*     data EIBND(20) / C1 /
	 *
	 *++  Code for ARGM is inactive
	 *      RETURN
	 *      ENTRY DIVACE
	 *++  End
	 * ********
	 * START OF CODE
	 * ******** */
	l = divasc.ndtf - 1;
	if (divamc.ics != 1)
		l += (divamc.ics - 1)*divasc.numdt;
	for (i = divamc.ics; i <= divamc.icf; i++)
	{
		if (divasc.nkdko != 0)
			divasc.kordi = Kord[divasc.nkdko + i - 1];
		divamc.iy += labs( divasc.kordi );
		kql = Kord[i + 3];
		kqn = labs( kql );
		kqd = max( 2, kqn );
		/* ********
		 * OBTAIN ERROR TOLERANCE SPECIFIED BY THE USER
		 * ******** */
		if (i <= divamc.ilgrep)
			switch (IARITHIF(kql))
			{
				case -1: goto L_2600;
				case  0: goto L_3310;
				case  1: goto L_2610;
			}
		divamc.itolep = labs( divamc.itolep ) + 1;
		eps = Tol[divamc.itolep];
		divamc.ilgrep = Lgroup[divamc.itolep];
		/*   TEST IF SIMPLE ABSOLUTE ERROR TEST IS BEING USED */
		if (divamc.ilgrep > 0)
			goto L_2580;
		jlgrep = divamc.ilgrep;
		/*     GET OLD RELATIVE ERROR FACTOR */
		tps6 = Tol[divamc.itolep + 1];
		divamc.ilgrep = Lgroup[divamc.itolep + 1];
		divamc.itolep = -divamc.itolep - 1;
 
		switch (IARITHIF(jlgrep + 1))
		{
			case -1: goto L_2540;
			case  0: goto L_2570;
			case  1: goto L_2510;
		}
		/*   NO CHECK ON THE ERROR ESTIMATE IS TO BE MADE */
L_2510:
		switch (ARITHIF(eps + C1))
		{
			case -1: goto L_2520;
			case  0: goto L_2590;
			case  1: goto L_2520;
		}
		/*   ERROR TOLERANCE IS SPECIFIED IMPROPERLY */
L_2520:
		divamc.kemax = i;
		divamc.neptol = divamc.itolep;
		divamc.linc = -3;
		divamc.erep = -fabs( divamc.erep );
		return;
		/*   COMPUTE NEW RELATIVE ERROR FACTOR */
L_2540:
		;
		*tps1 = C0;
		for (j = i; j <= divamc.ilgrep; j++)
		{
			*tps1 += fabs( F[j] );
		}
		*tps1 = fabs( *hh )**tps1/(double)( divamc.ilgrep - i + 1 );
		if (divamc.lsc <= 2)
			goto L_2560;
		/*     ON FIRST 3 STEPS INCREASE TPS6 WHEN COMPUTING REL. ERROR FACTOR */
		tps6 = fmax( C4**tps1, tps6 );
		/*     ON 1-ST TIME THROUGH THE FIRST STEP, REL. ERR. FAC. IS NOT STORED */
		if (divamc.lsc == 7)
			goto L_2570;
L_2560:
		;
		tps6 = fmax( *tps1, tps6 );
		/*   STORE NEW RELATIVE ERROR FACTOR */
		Tol[-divamc.itolep] = tps6*Ref[-jlgrep - 1];
		/*   COMPUTE ABSOLUTE ERROR TOLERANCE */
L_2570:
		eps *= tps6;
L_2580:
		if (eps <= C0)
			goto L_2520;
L_2590:
		switch (IARITHIF(kql))
		{
			case -1: goto L_2600;
			case  0: goto L_3330;
			case  1: goto L_2610;
		}
		/* END OF OBTAINING ERROR TOLERANCE
		 * ********
		 * OBTAIN INFORMATION USED FOR ERROR ESTIMATION, ORDER SELECTION, ETC.
		 * ********
		 * EQUATION IS STIFF */
L_2600:
		;
		/*++  Code for STIFF is inactive
		 *      JS=abs(KORD(NJSKO+I-1))-1
		 *      JSI=JS
		 *      TPP=C0
		 *      TPS4=F(L+KQD+2)
		 *      TPS3=F(L+KQD+1)
		 *      TPS2=F(L+KQD)
		 *      TPS1=F(L+KQD-1)
		 *      IF (KQD.EQ.2) TPS1=Y(IY-1)
		 *      E=ABS(TPS3)+ABS(TPS4)
		 *      EI=E+ABS(TPS2)
		 *      RND=EI
		 *      IF (KORDI.GE.0) GO TO 2604
		 *c.    EQUATION IS IMPLICIT
		 *      JSI=JSI-1
		 *      IF (JSI.NE.0) GO TO 2604
		 *      IF (KORDI.EQ.-1) GO TO 2602
		 *      ERCOEF=GS(KQN+1)
		 *      GO TO 2606
		 * 2602 ERCOEF=.5D0*DS(KQD,1)
		 *      JSI=1
		 *      GO TO 2606
		 *c.    END OF SPECIAL CODE FOR IMPLICIT EQUATIONS
		 * 2604 ERCOEF = DS(KQD,JSI)
		 * 2606 ERCOEF = ABS(ERCOEF) / EPS
		 *      IF (LSC.LE.2)  GO TO 2710
		 *      IF (LSC-5) 2650,2710,2710
		 *c.  END OF CODE FOR STIFF EQUATIONS
		 *++  End
		 *
		 * EQUATION IS NOT STIFF */
L_2610:
		tpp = F[i] - F[l + 1];
		*tps3 = tpp;
		*tps4 = tpp - F[l + kqd + 1];
		*tps2 = tpp + F[l + kqd];
		*tps1 = tpp + F[l + kqd - 1];
		e = fabs( *tps3 ) + fabs( *tps4 );
		rnd = e;
		ei = e + fabs( *tps2 );
		ercoef = fabs( Gs[kqn + 1] )/eps;
		if (kql >= 4)
			goto L_2710;
		/*   TEST IF STARTING OR IF INTEGRATION ORDER IS ONE */
		if (divamc.lsc <= 2)
			switch (IARITHIF(kql - 2))
			{
				case -1: goto L_2660;
				case  0: goto L_2710;
				case  1: goto L_2710;
			}
		/* ********
		 * LOGIC ASSOCIATED WITH STARTING THE INTEGRATION
		 * ******** */
		*tps4 = C0;
		switch (IARITHIF(divamc.lsc - 4))
		{
			case -1: goto L_2650;
			case  0: goto L_2640;
			case  1: goto L_2620;
		}
		/* FIRST STEP */
L_2620:
		e *= CQ3125;
		*tps3 = C0;
		F[l + 4] = C0;
		s = C0;
		/*   TEST IF FIRST TIME THROUGH THE FIRST STEP */
		if (divamc.lsc == 7)
			goto L_2690;
		/*   COMPUTE S=ESTIMATE OF H * EIGENVALUE OF JACOBIAN = 2*(F(A)-F(B))/
		 *   (F(B)-F(C)) WHERE F(A)=CURRENT F(I), AND F(B) AND F(C) PRECEDING
		 *   VALUES OR ESTIMATES OF F(I) */
		tpp = F[i] - F[l + 5];
		*tps4 = tpp;
		e = C2*fabs( *tps4 );
		if (s != C0)
			s = (*tps4 + *tps4)/s;
		switch (ARITHIF(s + CP125))
		{
			case -1: goto L_2630;
			case  0: goto L_2700;
			case  1: goto L_2700;
		}
		/*     SET LDT=-5  TO INDICATE POSSIBLE PROBLEMS DUE TO INSTABILITY */
L_2630:
		divasc.ldt = -5;
		goto L_2690;
		/*   ADJUST CORRECTION MADE ON SECOND STEP */
L_2640:
		tpp *= CP8;
		/*   ADJUST ESTIMATED ERRORS ON SECOND AND THIRD STEPS */
L_2650:
		e = fabs( *tps3 );
		rnd = C4*e;
		goto L_2710;
		/* END OF SPECIAL LOGIC FOR STARTING
		 * ********
		 * INTEGRATION ORDER =1 IS TREATED AS A SPECIAL CASE
		 * ******** */
L_2660:
		tpp += F[l + 2];
		if (Beta[2] >= C1P4)
			ei *= C1000;
		/*   ESTIMATE NEW VALUE FOR S */
		s = F[l + 4];
		if (s == C0)
			goto L_2680;
		s = fmax( CM8, C2*Beta[2]*(*tps1 - *tps2 - F[l + 5])/s );
		if (s >= CMP5)
			goto L_2670;
		/*   MODIFY TPP (TO GET BETTER STABILITY CHARACTERISTICS) */
		tpp *= fmax( CP25, (CM2 - C2*s)/(s*s) );
L_2670:
		*tps4 *= fabs( s );
L_2680:
		e = CP25*(e + fabs( *tps4 ));
		ei += fabs( *tps4*s );
		/*     STORE INFORMATION REQUIRED TO ESTIMATE S ON NEXT STEP */
L_2690:
		F[l + 4] = tpp;
L_2700:
		F[l + 5] = F[i];
		/* END OF SPECIAL CODE FOR INTEGRATION ORDER =1
		 * ********
		 * CODE FOR NOISE TEST AND GETTING ERROR ESTIMATE
		 * ******** */
L_2710:
		e *= ercoef;
		rnoise = C0;
		if (eps < C0)
			goto L_2810;
		tps5 = fabs( F[l + 2] ) + fabs( F[i] );
		if (tps5 == C0)
			goto L_2760;
L_2720:
		rnoise = rnd/tps5;
		if (rnoise > Rbq[kqd])
			switch (ARITHIF(rnoise - C1))
			{
				case -1: goto L_2760;
				case  0: goto L_2750;
				case  1: goto L_2750;
			}
		/*   NOISE IS APPARENTLY SLOWING CONVERGENCE OF THE DIFFERENCES
		 *     REDUCE EI */
		ei = rnd;
		tps5 = fabs( divaev.eeps2*Y[divamc.iy - 1] )/eps;
		if (tps5 < fabs( e ))
			switch (IARITHIF(divamc.lsc))
			{
				case -1: goto L_2730;
				case  0: goto L_2730;
				case  1: goto L_2760;
			}
		e = tps5;
		rnoise = C0;
L_2730:
		e = -fabs( e );
		if (divamc.eimin > CP1)
			ei *= C10*divamc.eimin;
		/*     COMPUTE REDUCTION TO BE MADE IN EI */
		if (rnoise > (C20*Rbq[kqd]))
			goto L_2760;
		k = -6 - divamc.lsc;
L_2740:
		if (k <= 0)
			goto L_2760;
		/*     REDUCE EI WHEN NOISE APPARENTLY LIMITS PRECISION */
		k -= 1;
		ei *= CP5;
		if (ei > divamc.eimin)
			goto L_2740;
		goto L_2760;
L_2750:
		*tps4 = 1.1e0*rnd;
		*tps3 = rnd;
L_2760:
		;
		/*   TEST FOR STIFFNESS GOES HERE WHEN IMPLEMENTED
		 * *       INGREDIENTS OF TEST MAY INCLUDE --
		 * *       RNOISE, WHETHER (ABS(TPS4).GT.ABS(TPS3)),
		 * *       WHETHER EMAX IS INCREASING, RESULT OF TEST ON
		 * *       PREVIOUS STEPS, ETC.
		 *
		 * ********
		 * COMPUTE ERROR ESTIMATES AND INFORMATION FOR SELECTING THE STEPSIZE
		 * ******** */
		if (e >= fabs( divamc.emax ))
			goto L_2770;
		if (-e <= fabs( divamc.emax ))
			goto L_2780;
		divamc.snoise = rnoise;
		divamc.dnoise = rnd;
		divamc.noiseq = kqd;
		/*   STORE PARAMETERS ASSOCIATED WITH LARGEST VALUE OF E */
L_2770:
		divamc.emax = e;
		divamc.kemax = i;
		divamc.neptol = divamc.itolep;
		/*   DETERMINE HOW MUCH STEPSIZE CAN BE INCREASED */
L_2780:
		ei *= ercoef*Sigma[kqd];
		divamc.eimax = fmax( divamc.eimax, ei );
		if (divamc.linc <= 0)
			goto L_2810;
		k = 0;
L_2790:
		if (ei >= fmin( divamc.eimin, Eibnd[kqn] ))
			goto L_2800;
		k += 1;
		if (k == divamc.linc)
			goto L_2810;
		ei *= Sigma[kqd];
		goto L_2790;
L_2800:
		divamc.linc = k;
		/* END OF COMPUTING ERROR ESTIMATES */
L_2810:
		;
		/*++  Code for ERRSTO is inactive
		 *      IF (IOP20 .EQ. 0) GO TO 780
		 *c.********
		 *c.STORE ERROR ESTIMATE (OPTIONAL)
		 *c.********
		 *      F(IOP20+I-1)=TPS3*GS(KQN+1)
		 *c.END OF STORING ERROR ESTIMATE
		 *++  Code for INTEGO | ERRSTO is active */
		if (divamc.iop19 == 0)
			goto L_3090;
		/*.********
		 *.EQUATIONS ARE GROUPED TO USE SAME INTEGRATION METHOD (OPTIONAL)
		 *.********
		 *++  Code for INTEGO is active */
		if (i > 1)
			switch (IARITHIF(i - ilgror))
			{
				case -1: goto L_2900;
				case  0: goto L_2900;
				case  1: goto L_2830;
			}
		itolor = divamc.iop19;
L_2830:
		jlgror = Kord[itolor];
		itolor += 1;
		if (jlgror > 0)
			goto L_2870;
		ilgror = Kord[itolor];
		itolor += 1;
		switch (IARITHIF(jlgror + 1))
		{
			case -1: goto L_2840;
			case  0: goto L_2850;
			case  1: goto L_2890;
		}
L_2840:
		if (jlgror < -2)
			switch (IARITHIF(kqd + jlgror))
			{
				case -1: goto L_2850;
				case  0: goto L_2880;
				case  1: goto L_2880;
			}
		/*.INITIALIZE FOR ACCUMULATING VARIABLES USED IN ORDER SELECTION */
L_2850:
		iord = i;
		kqlord = kql;
		for (k = 1; k <= 4; k++)
		{
			Tempao[k] = fabs( Tempa[k] );
		}
		goto L_2930;
		/*.ORDERS IN CURRENT GROUP CAN BE DIFFERENT */
L_2870:
		ilgror = jlgror;
		goto L_3090;
		/*.ORDER IS NOT GOING TO BE CHANGED */
L_2880:
		jlgror = 0;
L_2890:
		switch (IARITHIF(kql))
		{
			case -1: goto L_3240;
			case  0: goto L_3270;
			case  1: goto L_3270;
		}
		/*.TAKE ACTION FOR EQUATION WHICH IS NOT THE FIRST IN THE GROUP */
L_2900:
		switch (IARITHIF(jlgror))
		{
			case -1: goto L_2910;
			case  0: goto L_2890;
			case  1: goto L_3090;
		}
		/*.ACCUMULATE VARIABLES USED IN ORDER SELECTION */
L_2910:
		for (k = 1; k <= 4; k++)
		{
			Tempao[k] += fabs( Tempa[k] );
		}
		/*.    TEST IF THIS IS LAST EQUATION IN THE GROUP */
L_2930:
		if (i != ilgror)
			switch (IARITHIF(kql))
			{
				case -1: goto L_3310;
				case  0: goto L_3290;
				case  1: goto L_3290;
			}
		/*.SET UP TO GO SELECT INTEGRATION ORDER */
		kql = 0;
		for (k = 1; k <= 4; k++)
		{
			Tempa[k] = Tempao[k];
		}
		goto L_3090;
		/*.INTEGRATION ORDER HAS BEEN SELECTED
		 *++  Code for INTEGO | STIFF is active */
L_2950:
		;
		/*++  Code for INTEGO is active */
		kql = kqlord;
		switch (IARITHIF(kqn - labs( kql )))
		{
			case -1: goto L_2960;
			case  0: goto L_2980;
			case  1: goto L_3020;
		}
		/*.  TEST IF ORDER CAN BE DECREASED */
L_2960:
		if (jlgror >= -2)
			switch (IARITHIF(kql))
			{
				case -1: goto L_3010;
				case  0: goto L_3040;
				case  1: goto L_3040;
			}
		/*.    INTEGRATION ORDER WAS SELECTED OUTSIDE PERMITTED RANGE */
L_2970:
		kqn = labs( kql );
		/*.    INTEGRATION ORDER IS NOT GOING TO BE CHANGED */
L_2980:
		if ((kql != 1) || (divamc.lsc > 0))
			switch (IARITHIF(kql))
			{
				case -1: goto L_3030;
				case  0: goto L_3040;
				case  1: goto L_3040;
			}
		/*.    SET  4-TH ENTRY IN DIFFERENCE TABLES SO THAT STANDARD ADAMS
		 *.    METHOD IS USED WHEN KQL=1 */
L_2990:
		for (k = iord; k <= i; k++)
		{
			F[divasc.ndtf + k*divasc.numdt - divasc.numdt + 3] = C0;
		}
		goto L_3270;
		/*.  ORDER FOR STIFF EQUATION WAS REDUCED */
L_3010:
		;
		/*++  Code for INTEGO & STIFF is inactive
		 *      IF (KQN.LT.JSI) GO TO 990
		 *      TPP=-C1
		 *      GO TO 1090
		 *c.  TEST IF ORDER CAN BE INCREASED
		 *++  Code for INTEGO is active */
L_3020:
		if (jlgror == -2)
			goto L_2970;
		/*++  Code for INTEGO & STIFF is inactive
		 *      IF (KQL.GE.0) GO TO 1140
		 *      IF ((JSI.NE.0).AND.(KQN.GT.(MAXKQD+JSI))) GO TO 990
		 *      TPP=C1
		 *c.  STORE RESULTS FOR STIFF EQUATIONS
		 *++  Code for INTEGO is active */
L_3030:
		;
		/*++  Code for INTEGO & STIFF is inactive
		 *      DO 3035 K=IORD,I
		 *      KORD(K+3) = -KQN
		 * 3035 F(NDTF+K*NUMDT-NUMDT)=TPP
		 *      GO TO 3245
		 *c.  STORE RESULTS FOR EQUATIONS WHICH ARE NOT STIFF
		 *++  Code for INTEGO is active */
L_3040:
		ll = divasc.ndtf + divasc.numdt*iord - divasc.numdt;
		for (j = iord; j <= i; j++)
		{
			Kord[j + 3] = kqn;
			switch (IARITHIF(kqn - kql))
			{
				case -1: goto L_3050;
				case  0: goto L_3070;
				case  1: goto L_3060;
			}
L_3050:
			F[ll + kqd - 1] += F[j] - F[ll];
			goto L_3080;
L_3060:
			F[ll + kqn] = F[ll + kqd];
L_3070:
			F[ll + kqd] = C0;
L_3080:
			ll += divasc.numdt;
		}
		switch (IARITHIF(kqn - 1))
		{
			case -1: goto L_3270;
			case  0: goto L_2990;
			case  1: goto L_3270;
		}
		/*++  End
		 *.********
		 *.SELECT INTEGRATION ORDER
		 *.******** */
L_3090:
		if (divamc.lsc <= 0)
			goto L_3120;
		/*. SPECIAL ORDER SELECTION WHEN STARTING */
		switch (IARITHIF(divamc.lsc - 3))
		{
			case -1: goto L_3110;
			case  0: goto L_3210;
			case  1: goto L_3100;
		}
L_3100:
		if (divamc.lsc == 5)
			switch (ARITHIF(s + .125e0))
			{
				case -1: goto L_3160;
				case  0: goto L_3130;
				case  1: goto L_3130;
			}
		switch (IARITHIF(divamc.lsc - 6))
		{
			case -1: goto L_3130;
			case  0: goto L_3130;
			case  1: goto L_3210;
		}
L_3110:
		if (C40*fmin( fabs( *tps4 ), fabs( *tps3 ) ) > fabs( *tps2 ))
		{
			if (eps != -C1)
				divamc.lsc = 2;
		}
		if (fabs( *tps4 ) < fabs( *tps3 ))
			switch (ARITHIF(C4*fabs( *tps4 ) - fabs( *tps2 )))
			{
				case -1: goto L_3130;
				case  0: goto L_3130;
				case  1: goto L_3210;
			}
		/*.  CHECK IF ORDER CAN BE INCREASED OR SHOULD BE DECREASED */
L_3120:
		tps5 = divamc.robnd*fabs( *tps4 );
		tps6 = divamc.robnd*(tps5 + fabs( *tps3 ));
		tps7 = fabs( *tps1 ) + fabs( *tps2 );
		if (tps5 >= fabs( *tps3 ))
			goto L_3140;
		if (tps6 >= tps7)
			goto L_3210;
L_3130:
		if (kqn >= divamc.maxkqi)
			goto L_3210;
		/*.    INCREASE THE INTEGRATION ORDER */
		kqn += 1;
		/*++  Code for INTEGO | STIFF is active */
		switch (IARITHIF(kql))
		{
			case -1: goto L_3230;
			case  0: goto L_2950;
			case  1: goto L_3250;
		}
		/*++  Code for ~(INTEGO | STIFF) is inactive
		 *      GO TO 3250
		 *++  End
		 *.  CHECK IF ORDER SHOULD BE DECREASED */
L_3140:
		if (tps6 < tps7)
			goto L_3210;
		if (tps5 < fabs( *tps3 - *tps4 ))
			goto L_3210;
		if ((*tps3 == *tps4) && (divamc.lsc <= 0))
			goto L_3210;
		switch (IARITHIF(kqn - 2))
		{
			case -1: goto L_3210;
			case  0: goto L_3160;
			case  1: goto L_3180;
		}
L_3160:
		kqn = 1;
		/*++  Code for INTEGO | STIFF is active */
		switch (IARITHIF(kql))
		{
			case -1: goto L_3220;
			case  0: goto L_2950;
			case  1: goto L_3170;
		}
		/*++  End
		 *.    WHEN ORDER IS REDUCED TO 1 WITH ADAMS METHOD SET F(L+4)=0 */
L_3170:
		F[l + 4] = C0;
		goto L_3260;
		/*.    DECREASE THE INTEGRATION ORDER */
L_3180:
		kqn -= 1;
		/*++  Code for INTEGO | STIFF is active */
		switch (IARITHIF(kql))
		{
			case -1: goto L_3220;
			case  0: goto L_2950;
			case  1: goto L_3200;
		}
		/*++  End */
L_3200:
		F[l + kqd] += tpp;
		goto L_3260;
		/*   NO CHANGE IN INTEGRATION ORDER IS BEING MADE */
L_3210:
		;
		/*++  Code for INTEGO is active */
		switch (IARITHIF(kql))
		{
			case -1: goto L_3240;
			case  0: goto L_2950;
			case  1: goto L_3270;
		}
		/*++  Code for ~INTEGO is inactive
		 *         TPS1 = EEPS10
		 *      GO TO 1530
		 *++  End
		 * END OF SELECTING INTEGRATION ORDER
		 * ********
		 * COMPUTE MAXIMUM INTEGRATION ORDERS AND SET NEW ONES (IF ANY)
		 * ********
		 * EQUATION IS STIFF
		 *     ORDER WAS DECREASED
		 *++  Code for INTEGO | STIFF is active */
L_3220:
		;
		/*++  Code for STIFF is inactive
		 *      IF (KQN.LT.JSI) GO TO 3236
		 *      F(L+1)=-C1
		 *      GO TO 3233
		 *c.    ORDER WAS INCREASED
		 *++  Code for INTEGO |  STIFF  is active */
L_3230:
		;
		/*++  Code for STIFF is inactive
		 *      IF ((JSI.NE.0).AND.(KQN.GT.(MAXKQD+JSI))) GO TO 3236
		 *      F(L+1)=C1
		 * 3233 KORD(I+3) = -KQN
		 *      GO TO 3245
		 *      ORDER WAS SET TO AN UNACCEPTABLE VALUE
		 * 3236 KQN=abs(KQL)
		 *      ORDER IS NOT BEING CHANGED
		 *++  Code for STIFF |  INTEGO is active */
L_3240:
		;
		/*++  Code for STIFF is inactive
		 *      F(L+1)=C0
		 * 3245 IF (JSI.NE.0) KQMAXD=max(KQN,KQMAXD)
		 *      IF (JS.LT.abs(KORDI)) KQMAXI=max(KQN,KQMAXI)
		 *      GO TO 3290
		 *++  End
		 * EQUATION IS NOT STIFF
		 *     ORDER INCREASED */
L_3250:
		F[l + kqn + 1] = -F[l + kqd + 1];
		if (divamc.lsc > 0)
			F[l + kqn + 1] = F[l + 1] - F[i];
		/*     ORDER CHANGED */
L_3260:
		Kord[i + 3] = kqn;
L_3270:
		divasc.kqmaxi = max( kqn, divasc.kqmaxi );
		if (eps > C0)
			divamc.kqmaxs = max( kqn, divamc.kqmaxs );
		F[l + kqd + 1] = C0;
L_3290:
		;
		if (kqn > divamc.kis)
			goto L_3310;
		/*.********
		 *.DETERMINE IF TIME TO STORE SOLUTION (OPTIONAL)
		 *.******** */
		if (divamc.kis >= 1000)
		{
			tp2 = fmax( 1.5e0, (double)( kqn )*powi(C2,1001 - divamc.kis) )*
			 fabs( *tps4 );
L_3295:
			if (tp2 > fabs( F[l + kqn] ))
			{
				if (kqn <= kql)
				{
					kqn -= 1;
					if (kqn > 1)
						goto L_3295;
					kqn = 1;
				}
			}
			Kord[i + 3] = kqn;
			if (i == 1)
				lkqmax = 0;
			lkqmax = max( kqn, lkqmax );
			divasc.kqmaxi = lkqmax;
			if (divamc.kis == 1000)
			{
				if (i == divamc.kemax)
					divamc.emax = (double)( 8 + SQ(kqn) )*fabs( divamc.emax );
				goto L_3325;
			}
			/*++  Code for DUMP is active */
		}
		else if ((e != C0) && (eps > C0))
		{
			if (divamc.iop9 > 0)
				switch (ARITHIF((fabs( e )*powi((double)( divamc.kis -
				 kqn + 2 ),kqn + 1)) - 1.e-2))
				{
					case -1: goto L_3310;
					case  0: goto L_3310;
					case  1: goto L_3300;
				}
L_3300:
			divamc.kis = -1;
			/*++  End */
		}
L_3310:
		;
		/* ********
		 * CORRECT
		 * ******** */
		for (k = 1; k <= divasc.kordi; k++)
		{
			/*++  Code for ~{p,x} is active */
			Y[divamc.iy - k] += divamc.g[k - 1][kql]*tpp;
			/*++  Code for {p,x} is inactive
			 *c--D Next line special: P=>D, X=>Q
			 *            Y(IY - K) = Y(IY - K) + dble(G(KQL + 1, K)) * dble(TPP)
			 *++  END */
		}
		/* END OF CORRECTING */
L_3325:
		;
		/*++  Code for OUTPUT is active */
		if (divamc.iop10 > 0)
		{
			if (i == 1)
			{
				Idat[1] = divamc.kstep;
				Idat[2] = divamc.lsc;
				Idat[3] = divamc.ksc;
				Idat[4] = divamc.iop11;
				Fdat[1] = divasc.tn;
				Fdat[2] = *hh;
				Fdat[3] = divamc.eimin;
				Fdat[4] = divamc.eave;
				Fdat[5] = Sigma[divamc.iop11];
				Fdat[6] = divamc.robnd;
				Mact2[3] = divasc.nte;
				/*--D Next line special: P=>S, X=>D */
				dmess( mact1, (char*)mtxtaa,105, divamc.idat, divamc.fdat );
				koutko = divamc.noutko;
			}
			if (koutko != 0)
			{
				if (Kord[koutko] > 0)
				{
					if (i < Kord[koutko])
						goto L_3328;
					koutko += 1;
				}
				else
				{
					if (i >= labs( Kord[koutko] ))
						koutko += 1;
				}
			}
			Idat[1] = i;
			Idat[2] = kql;
			Idat[3] = divamc.linc;
			Fdat[1] = e;
			Fdat[2] = ei;
			Fdat[3] = eps;
			Fdat[4] = F[i];
			Fdat[5] = *tps1;
			Fdat[6] = *tps2;
			Fdat[7] = *tps3;
			Fdat[8] = *tps4;
			Fdat[9] = rnoise;
			Fdat[10] = 0.e0;
			if (kql == 1)
				Fdat[10] = s;
			Fdat[11] = Beta[kqd];
			/*--D Next line special: P=>S, X=>D */
			dmess( mact2, (char*)mtxtab,89, divamc.idat, divamc.fdat );
L_3328:
			if (i == divasc.nte)
				divamc.iop10 -= 1;
		}
		/*++  End */
L_3330:
		l += divasc.numdt;
	}
	return;
} /* end of function */
/*   End of DIVACR */
 
		/* PARAMETER translations */
#define	CP5625	.5625e0
#define	CRBQI	.421875e0
		/* end of PARAMETER translations */
 
void /*FUNCTION*/ divahc()
{
	long int j, k, n;
	double temp, tp, tp1, tp2;
	 static double gg[MAXORD - 1 + 1/MAXORD], w[KDIM + MAXORD];
	 static double b[KDIM + MAXORD]={5.000000000000000000000000000000000000000e-1,
	 1.666666666666666666666666666666666666667e-1,8.333333333333333333333333333333333333333e-2,
	 5.000000000000000000000000000000000000000e-2,3.333333333333333333333333333333333333333e-2,
	 2.380952380952380952380952380952380952381e-2,1.785714285714285714285714285714285714286e-2,
	 1.388888888888888888888888888888888888889e-2,1.111111111111111111111111111111111111111e-2,
	 9.090909090909090909090909090909090909091e-3,7.575757575757575757575757575757575757576e-3,
	 6.410256410256410256410256410256410256410e-3,5.494505494505494505494505494505494505495e-3,
	 4.761904761904761904761904761904761904762e-3,4.166666666666666666666666666666666666667e-3,
	 3.676470588235294117647058823529411764706e-3,3.267973856209150326797385620915032679739e-3,
	 2.923976608187134502923976608187134502924e-3,2.631578947368421052631578947368421052632e-3,
	 2.380952380952380952380952380952380952381e-3,2.164502164502164502164502164502164502165e-3,
	 1.976284584980237154150197628458498023715e-3};
	double *__restrict const hh = (double*)divamc.g;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const Alpha = &divamc.alpha[0] - 1;
	double * __restrict const B = &b[0] - 1;
	double * __restrict const Beta = &divamc.beta[0] - 1;
	double * __restrict const Gg = &gg[0] - 1;
	double * __restrict const Gs = &divamc.gs[0] - 1;
	double * __restrict const Rbq = &divamc.rbq[0] - 1;
	double * __restrict const Sigma = &divamc.sigma[0] - 1;
	double * __restrict const V = &divamc.v[0] - 1;
	double * __restrict const W = &w[0] - 1;
	double * __restrict const Xi = &divasc.xi[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1988-05-20 DIVAHC Krogh   Initial code.
	 *
	 * SUBROUTINE TO COMPUTE COEFFICIENTS REQUIRED FOR INTEGRATING
	 * ORDINARY DIFFERENTIAL EQUATIONS
	 * */
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
 
	/*.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS. */
	/*                 K - 1 + 1 / K  is equivalent to max(1, K-1) */
	/*++  Code for STIFF is inactive
	 *      INTEGER          GODIF
	 *++  End */
 
 
	/*  B(K)= 1/(K*(K+1))
	 *++ Save data by elements if ~.C.
	 *++ Of next 23 lines, only the first KDIM+MAXORD are active */
	/*     data B(23) / 1.811594202898550724637681159420289855072D-3 /
	 *
	 * ********
	 * START OF CODE
	 * ********
	 *     SET STEP NUMBER OF METHOD
	 *++  Code for STIFF is inactive
	 *      IOP11 = MIN(max(KQMAXI,KQMAXD) + 1), KDIM)
	 *++  Code for ~STIFF is active */
	divamc.iop11 = min( divasc.kqmaxi + 1, KDIM );
	/*++  End
	 *     TEST IF STEPSIZE WAS CHANGED */
	if (divamc.kqicon >= 0)
		goto L_3510;
	/* ********
	 * STEPSIZE JUST CHANGED
	 * ********
	 *     SET CONSTANTS DEPENDING ON NEW STEPSIZE */
	divamc.kqmxil = divasc.kqmaxi;
	tp1 = *hh;
	Gg[1] = tp1*tp1;
	divamc.g[1][0] = Gg[1]*CP5;
	if (divasc.maxint <= 2)
		goto L_3450;
	for (k = 3; k <= divasc.maxint; k++)
	{
		Gg[k - 1] = divamc.g[k - 2][0]*tp1;
		divamc.g[k - 1][0] = Gg[k - 1]/(double)( k );
	}
	/*     SET CONSTANTS INDICATING STEP CHANGE */
L_3450:
	divamc.kqicon = 0;
	/*++  Code for STIFF is inactive
	 *      KQDCON=0
	 *++  End */
	divamc.kqmxip = 1;
	divamc.ksc = 1;
	if (divamc.lsc < 7)
		goto L_3490;
	/*     SPECIAL SET-UP OF CONSTANTS ON THE VERY FIRST STEP */
	divamc.hincc = C1P125;
	divamc.lincd = 6;
	divamc.lincq = 12;
	if (divamc.hinc > C0)
		goto L_3460;
	divamc.lincd = -2;
	divamc.linc = -2;
	divamc.robnd = C1;
L_3460:
	Sigma[1] = 1.0e0;
	Beta[1] = C1;
	for (n = 1; n <= divamc.iop11; n++)
	{
		/*++  Code for STIFF is inactive
		 *      D(1,N)=C0
		 *++  End */
		Xi[n] = tp1;
		Alpha[n] = C1;
		Beta[n + 1] = C1;
		Sigma[n + 1] = (double)( n + 1 )*Sigma[n]*divamc.hincc;
	}
	temp = divaev.eeps16;
	Rbq[1] = C1;
	Rbq[2] = CP1;
	tp = CRBQI;
	/*     **** IN THE LOOP BELOW RBQ(K) IS COMPUTED TO BE
	 *          APPROXIMATELY (3/4 ** ((K-1) ** 2 - 1) / 10
	 *          .5625 = (3/4) ** 2    TP = (3/4) ** (2*K -3) */
	for (k = 3; k <= KDIM; k++)
	{
		temp += temp;
		Rbq[k] = fmax( temp, Rbq[k - 1]*tp );
		tp *= CP5625;
	}
	goto L_3560;
	/*     SET-UP AFTER THE FIRST STEP */
L_3490:
	tp2 = Xi[1];
	Xi[1] = tp1;
	Beta[2] = tp1/tp2;
	k = 2;
	if (divamc.hincc == divamc.hinc)
		goto L_3540;
	if ((divamc.lsc != 0) || ((divamc.kstep - divamc.ksstrt - divamc.kqmaxs) <
	 10))
		goto L_3540;
	divamc.hincc = C1;
	divamc.lincd = 0;
L_3500:
	divamc.lincd += 1;
	divamc.hincc *= divamc.hinc;
	if (divamc.hincc < 2.e0)
		goto L_3500;
	divamc.linc = (divamc.linc*(divamc.lincd + divamc.lincd))/divamc.lincq;
	divamc.lincq = divamc.lincd + divamc.lincd;
	divamc.hincc = divamc.hinc;
	goto L_3540;
	/* END OF LOGIC FOR CASE WHEN STEPSIZE JUST CHANGED
	 *     TEST IF MAXIMUM INTEGRATION ORDER DID NOT INCREASE */
L_3510:
	if (divasc.kqmaxi > divamc.kqmxil)
	{
		/* ********
		 * INTEGRATION ORDER WAS INCREASED -- GET NEW V'S
		 * ******** */
		divamc.kqmxil = divasc.kqmaxi;
		divamc.kqmxip = divamc.kqmxil + divasc.maxint;
		k = divamc.kqmxip;
		V[k] = B[k];
		if (divamc.kqicon == 1)
			goto L_3530;
		/*     if (KQICON .eq. K) KQICON = KQICON - 1 --- Removed 1999-08-19 */
		for (n = 2; n <= divamc.kqicon; n++)
		{
			k -= 1;
			V[k] += -Alpha[n]*V[k + 1];
		}
		/* END OF GETTING NEW V'S */
	}
	else
	{
		divamc.iop11 = max( divamc.iop11, divamc.kqmxil + 1 );
	}
L_3530:
	if (divamc.iop11 <= divamc.ksc)
		goto L_3560;
	/* ********
	 * COMPUTE PARAMETERS WHICH ARE STILL CHANGING AS A RESULT OF
	 * A CHANGE IN THE STEPSIZE
	 * ******** */
	tp2 = Xi[divamc.ksc];
	/*     UPDATE CONSTANT STEP COUNTER */
	divamc.ksc += 1;
	k = divamc.ksc;
	Beta[k] = C1;
L_3540:
	;
	temp = divamc.hincc;
 
	/*   LOOP TO COMPUTE NEW VALUES OF PARAMETERS */
	for (n = k; n <= divamc.iop11; n++)
	{
		tp1 = tp2 + *hh;
		tp2 = Xi[n];
		Xi[n] = tp1;
		Alpha[n] = *hh/tp1;
		Beta[n + 1] = Beta[n]*(tp1/tp2);
		temp = fmax( temp, (double)( n )*(Alpha[n]*divamc.hincc) );
		Sigma[n] = Sigma[n - 1]*temp;
	}
	if (divamc.iop11 != KDIM)
		Xi[divamc.iop11 + 1] = tp2 + *hh;
	/* END OF CODE FOR COMPUTING PARAMETERS WHICH ARE STILL CHANGING
	 * */
L_3560:
	if (divamc.kqicon >= divamc.kqmxip)
		goto L_3690;
	/* ********
	 * COMPUTE INTEGRATION COEFFICIENTS WHICH ARE STILL CHANGING
	 * ******** */
	divamc.kqmxil = max( divasc.kqmaxi, divamc.kqmxil );
	divamc.kqmxip = divamc.kqmxil + divasc.maxint;
	j = divamc.kqmxip - divamc.kqicon;
	n = divamc.kqicon + 1;
	divamc.kqicon = n;
	if (n != 1)
		goto L_3580;
	/* INITIALIZE V AND W */
	for (k = 1; k <= j; k++)
	{
		V[k] = B[k];
		W[k] = V[k];
	}
	goto L_3600;
	/* UPDATE V AND INITIALIZE W */
L_3580:
	if (n == KDIM)
		goto L_3690;
	for (k = 1; k <= j; k++)
	{
		V[k] += -Alpha[n]*V[k + 1];
		W[k] = V[k];
	}
	/* SET TRANSFER FOR LOOP BELOW DEPENDING ON VALUE OF MAXINT */
L_3600:
	;
	goto L_3660;
 
L_3640:
	j -= 1;
	/* INNER LOOP FOR COMPUTING INTEGRATION COEFFICIENTS */
	for (k = 1; k <= j; k++)
	{
		W[k] += -Alpha[n]*W[k + 1];
	}
	/*     STORE INTEGRATION COEFFICIENTS */
L_3660:
	divamc.g[0][n] = *hh*W[1];
	Gs[n + 1] = divamc.g[0][n] - divamc.g[0][n - 1];
	/*++  Code for MAXORD >= 2 is active */
	if (divasc.maxint >= 2)
	{
		divamc.g[1][n] = Gg[1]*W[2];
		/*++  Code for MAXORD >= 3 is inactive
		 *        if (MAXINT .gt. 2) then
		 *           DO 3665 K=3,MAXINT
		 *3665          G(N+1,K)=GG(K-1)*W(K)
		 *        end if
		 *++  Code for MAXORD >= 2 is active */
	}
	/*++  End */
	n += 1;
	if (n <= divamc.kqmxil)
		goto L_3640;
	/* END OF COMPUTING INTEGRATION COEFFICIENTS
	 * */
L_3690:
	;
	/*++  Code for STIFF is inactive
	 *      IF (KQDCON.GT.KQMAXD) GO TO 4662
	 *c.********
	 *c.COMPUTE DIFFERENTIATION COEFFICIENTS WHICH ARE STILL CHANGING
	 *c.********
	 *c.SET TRANSFER FOR LOOP BELOW, DEPENDING ON VALUE OF MAXDIF
	 *++  Code for STIFF & MAXORD >= 2 is inactive
	 *      IF (MAXDIF-2) 3693,3692,3691
	 * 3691 ASSIGN 3696 TO GODIF
	 *      GO TO 3694
	 * 3692 ASSIGN 3698 TO GODIF
	 *      GO TO 3694
	 * 3693 ASSIGN 3699 TO GODIF
	 *++  Code for STIFF is inactive
	 * 3694 KQDCON=KQDCON+1
	 *c.LOOP FOR COMPUTING DIFFERENTIATION COEFFICIENTS
	 *      DO 3699 N=KQDCON,KQMAXD
	 *      DS(N+1,2)=C1/XI(N)
	 *      D(N+1,1)=DS(N+1,2)+D(N,1)
	 *      DS(N+1,1)=DS(N+1,2)/D(N+1,1)
	 *++  Code for STIFF & MAXORD >= 2 is inactive
	 *      GO TO GODIF, (3696,3698,3699)
	 * 3696 CONTINUE
	 *++  Code for STIFF & MAXORD >= 3 is inactive
	 *      DO 3697 K=3,MAXDIF
	 *      DS(N+1,K)=D(N,K-2) * (K-1)/XI(N)
	 * 3697 D(N+1,K-1)=DS(N+1,K) + D(N,K-1)
	 *++  Code for STIFF is inactive
	 * 3698 CONTINUE
	 *++  Code for STIFF & MAXORD >= 2 is inactive
	 *      D(N+1,MAXDIF)=D(N,MAXDIF) + D(N,MAXDIF-1) * (MAXDIF)/XI(N)
	 *++  Code for STIFF is inactive
	 * 3699 CONTINUE
	 *++  End
	 *
	 * END OF COMPUTING DIFFERENTIATION COEFFICIENTS */
	return;
} /* end of function */
/*   End of DIVAHC */
 
		/* PARAMETER translations */
#undef	LTXTAC
#define	LTXTAC	41
#undef	LTXTAD
#define	LTXTAD	94
		/* end of PARAMETER translations */
 
void /*FUNCTION*/ divain(
double t[],
double y[],
double f[],
long kord[])
{
	LOGICAL32 lnotm1;
	long int i, ici, idat[1], idt, integ, integz, interp, iy, iyi,
	 iyn, iyni, j, k, kqmxi, kqmxs, kqq, l, n;
	double c[KDIM + MAXORD - 1], csum[KDIM + MAXORD - 1], eta[KDIM],
	 fdat[6], gamma[KDIM], hi, tp1, xp1;
	static char mtxtaa[1][155]={"DIVAIN$BInterpolating at T(1)=$F with $BTN=$F, T(2)=$F and H=$F.  T(1) must be in [$F, $F].$Einternal variable LDT = $I.  Interpolation not allowed now.$E"};
	static long mact[8]={MEEMES,0,0,0,MENTXT,0,METEXT,MERET};
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const C = &c[0] - 1;
	double * __restrict const Csum = &csum[0] - 1;
	double * __restrict const Eta = &eta[0] - 1;
	double * __restrict const F = &f[0] - 1;
	double * __restrict const Fdat = &fdat[0] - 1;
	double * __restrict const Gamma = &gamma[0] - 1;
	long * __restrict const Idat = &idat[0] - 1;
	long * __restrict const Kord = &kord[0] - 1;
	long * __restrict const Mact = &mact[0] - 1;
	double * __restrict const T = &t[0] - 1;
	double * __restrict const Xi = &divasc.xi[0] - 1;
	double * __restrict const Y = &y[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1988-01-14 DIVAIN Krogh   Initial code.
	 *
	 *  SUBROUTINE TO DO INTERPOLATION FOR VARIABLE ORDER INTEG. ROUTINE
	 * */
	/*--D Next line special: P=>D, X=>Q */
	/*++ Substitute for KDIM, MAXORD below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
	/*              Stuff for processing error messages */
	/* ********* Error message text ***************
	 *[Last 2 letters of Param. name]  [Text generating message.]
	 *AA DIVAIN$B
	 *AB Interpolating at T(1)=$F with $B
	 *AC TN=$F, T(2)=$F and H=$F.  T(1) must be in [$F, $F].$E
	 *AD internal variable LDT = $I.  Interpolation not allowed now.$E */
 
	/*                      1 2 3 4       5 6       7      8 */
 
	/*++  Code for ARGM is inactive
	 *      ENTRY DIVAIE
	 *++  End
	 * ********
	 * START OF CODE -- CHECK ON STATE OF DIFFERENCE TABLE
	 * ******** */
	l = divasc.ldt;
	switch (IARITHIF(l))
	{
		case -1: goto L_3710;
		case  0: goto L_3730;
		case  1: goto L_3780;
	}
L_3710:
	switch (IARITHIF(l + 2))
	{
		case -1: goto L_4170;
		case  0: goto L_3730;
		case  1: goto L_3720;
	}
L_3720:
	if (divasc.maxint >= 0)
		l = 1;
	goto L_3840;
	/* ********
	 * UPDATE DIFFERENCE TABLE TO START OF NEXT STEP
	 * ******** */
L_3730:
	k = divasc.ndtf;
	for (i = 1; i <= divasc.nte; i++)
	{
		kqq = Kord[i + 3];
		if (kqq <= 0)
			goto L_3760;
		/* EQUATION IS NOT STIFF */
		tp1 = F[i] - F[k];
		/* LOOP TO DO UPDATING */
		n = k + max( labs( kqq ), 2 );
		for (j = k; j <= n; j++)
		{
			F[j] += tp1;
		}
L_3760:
		;
		k += divasc.numdt;
	}
	divasc.ldt = 1;
	if (l != 0)
		return;
	/* END OF UPDATING DIFFERENCE TABLE
	 * ********
	 * INITIALIZE FOR COMPUTATION OF COEFFICIENTS
	 * ******** */
L_3780:
	interp = 0;
	hi = T[1] - divasc.tn;
	Gamma[1] = hi/Xi[1];
	switch (ARITHIF(Gamma[1]))
	{
		case -1: goto L_3790;
		case  0: goto L_3800;
		case  1: goto L_3810;
	}
L_3790:
	if (Gamma[1] >= -C1)
		goto L_3820;
	interp = 1;
	switch (ARITHIF(fabs( hi ) - fabs( T[2] )))
	{
		case -1: goto L_3820;
		case  0: goto L_3820;
		case  1: goto L_4180;
	}
L_3800:
	interp = 2 - divasc.kqmaxi;
	goto L_3820;
L_3810:
	if (Gamma[1] > C2)
		switch (IARITHIF(divasc.ldt - 2))
		{
			case -1: goto L_4180;
			case  0: goto L_3820;
			case  1: goto L_3820;
		}
L_3820:
	kqmxi = divasc.kqmaxi + interp - 1;
	/*++  Code for STIFF is inactive
	 *      KQMXS=max(KQMXI,KQMAXD)
	 *++  Code for ~STIFF is active */
	kqmxs = kqmxi;
	/*++  End */
	for (n = 2; n <= kqmxs; n++)
	{
		Gamma[n] = (hi + Xi[n - 1])/Xi[n];
	}
L_3840:
	lnotm1 = l != -1;
	integ = divasc.maxint;
	if (integ <= 0)
		switch (IARITHIF(integ + divasc.maxdif))
		{
			case -1: goto L_4160;
			case  0: goto L_3950;
			case  1: goto L_3950;
		}
	/* ********
	 * COMPUTE INTEGRATION COEFFICIENTS
	 * ********
	 *     INITIAL SET-UP
	 *         COMPUTE INITIAL C VALUES */
	for (n = 1; n <= integ; n++)
	{
		C[n] = hi/(double)( n );
	}
	i = integ + 1;
	integ += kqmxi;
	for (n = i; n <= integ; n++)
	{
		C[n] = C[n - 1]*((double)( n - divasc.maxint )/(double)( n ));
	}
	/*         COMPUTE ETA'S */
	for (n = 1; n <= kqmxi; n++)
	{
		Eta[n] = hi/Xi[n];
	}
	/*         COMPUTE C(K)'S TO CORRESPOND TO G(K-MAXINT+1,MAXINT),
	 *         K=MAXINT, MAXINT+1,..., MAXINT+KQMXI-1 */
	i = integ;
L_3880:
	j = integ;
	integ = j - 1;
	if (integ <= divasc.maxint)
		goto L_3900;
	for (n = j; n <= i; n++)
	{
		C[n] = Eta[n - integ]*C[n] + C[n - 1];
	}
	goto L_3880;
L_3900:
	for (n = j; n <= i; n++)
	{
		C[n] *= Eta[n - integ];
	}
	/*         END OF COMPUTING  G(---,MAXINT) */
	integz = 0;
	goto L_3940;
	/*         COMPUTE C(K)-S TO CORRESPOND TO G(K-INTEG+1,INTEG),
	 *         K=INTEG+1,INTEG+2,..., INTEG+KQMXI */
L_3920:
	for (n = 1; n <= kqmxi; n++)
	{
		C[integ + n] = Gamma[n]*C[integ + n - 1] - Eta[n]*C[integ + n];
	}
L_3940:
	ici = integ - 1;
	goto L_4020;
	/* END OF COMPUTING INTEGRATION COEFFICIENTS
	 * ********
	 * COMPUTE COEFFICIENTS FOR INTERPOLATION
	 * ******** */
L_3950:
	C[1] = C1;
	ici = 0;
	for (n = 1; n <= kqmxs; n++)
	{
		C[n + 1] = Gamma[n]*C[n];
	}
	switch (IARITHIF(integ + 1))
	{
		case -1: goto L_3970;
		case  0: goto L_3990;
		case  1: goto L_4010;
	}
	/* END OF COMPUTING INTERPOLATION COEFFICIENTS
	 *
	 *     SET-UP TO COMPUTE DIFFERENTIATION COEFFICIENTS REQUIRED
	 *     IN ORDER TO GET COEFFICIENTS ACTUALLY USED */
L_3970:
	integ = 0;
	ici = 1;
L_3980:
	integ -= 1;
	if (integ == divasc.maxint)
		ici = 0;
	/* ********
	 * COMPUTE DIFFERENTIATION COEFFICIENTS
	 * ******** */
L_3990:
	interp = max( interp, 0 );
	tp1 = (double)( -integ );
	C[1] = tp1*C[1]/Xi[-integ];
	j = divasc.kqmaxd + integ;
	for (n = 1; n <= j; n++)
	{
		C[n + 1] = (tp1*C[n])/Xi[n - integ] + Gamma[n - integ]*C[n];
	}
	/*     C(N) NOW CORRESPONDS TO THE DIFFERENTIAL COEFFICIENT
	 *          D(N-INTEG,-INTEG) */
L_4010:
	integz = integ;
	if (ici != 0)
		goto L_3980;
	/* END OF COMPUTING DIFFERENTIATION COEFFICIENTS
	 * ********
	 * BEGINNING OF LOOP TO DO
	 *         INTEGRATION       (INTEG.GT.0)
	 *         INTERPOLATION     (INTEG.EQ.0)
	 *         DIFFERENTIATION   (INTEG.LT.0)
	 * TO THE POINT INDICATED BY T.
	 * ********
	 *     SET UP INITIAL INDICES */
L_4020:
	if (divasc.nyny < 0)
	{
		iy = -divasc.nyny;
		iyni = divasc.nyny + ici + 1;
		if (divasc.ldt == 2)
		{
			Csum[ici + 1] = C[ici + 1];
			for (j = ici + 2; j <= (integ + kqmxi); j++)
			{
				Csum[j] = Csum[j - 1] + C[j];
			}
		}
	}
	else
	{
		iy = 1;
		iyni = divasc.nyny + ici - 1;
	}
	idt = divasc.ndtf - integz;
	for (i = 1; i <= divasc.nte; i++)
	{
		if (divasc.nkdko != 0)
			divasc.kordi = Kord[divasc.nkdko + i - 1];
		iy += labs( divasc.kordi );
		kqq = Kord[i + 3];
		/*         GET INDEX OF HIGHEST ORDER DIFFERENCE TO BE USED */
		k = max( labs( kqq ) + interp, 2 );
		iyi = -integ;
		switch (IARITHIF(kqq))
		{
			case -1: goto L_4030;
			case  0: goto L_4130;
			case  1: goto L_4040;
		}
		/* EQUATION IS STIFF */
L_4030:
		;
		/*++  Code for STIFF is inactive
		 *      JS=abs(KORD(NJSKO+I-1))-1
		 *      IYI=IYI-JS
		 *      IF(LNOTM1) IF (IYI) 4034,4032,4130
		 *      IF (KORDI.LT.0) IYI=IYI+1
		 *      IYI=IYI+MAXINT-abs(KORDI)
		 *      IF (IYI) 4034,4130,4130
		 *c.      IF EQUATION IS IMPLICIT DO NOT COMPUTE AN F
		 * 4032 IF (KORDI.LT.0) GO TO 4130
		 *c.      TEST IF INTEG TOO BIG FOR THIS EQUATION
		 * 4034 IF (abs(KORDI).LT.-IYI) GO TO 4130
		 *      IYI=IYI+IY
		 *      IYN=IYI+IYNI
		 *c. COMPUTE INNER PRODUCT FOR STIFF EQUATIONS
		 *      IF (INTEGZ.EQ.0) GO TO ???
		 *c.    DIFFERENTIATING
		 *      TP1 = C0
		 *      DO 4036 J = K+INTEGZ, 1, -1
		 *         TP1 = TP1 + C(J) * F(IDT+J-1)
		 * 4036 CONTINUE
		 *c.    TEST WHETHER TO STORE RESULT IN Y OR F
		 *      IF (IYI-IY) 4080, 4090, 4080
		 *c.    INTEGRATING OR INTERPOLATING
		 *      TP1 = C0
		 *      DO 4037 J = ICI + K, ICI + 2, -1
		 *         TP1 = TP1 + C(J) * F(IDT+J-ICI-1)
		 * 4037 CONTINUE
		 *      IF (INTEG.EQ.0) GO TO 4120
		 *      TP1=TP1 + C(ICI+1)*Y(IYN+1)
		 *++  End */
		goto L_4100;
		/* END OF SPECIAL CODE FOR STIFF EQUATIONS
		 *
		 * EQUATION IS NOT STIFF */
L_4040:
		if (lnotm1)
			switch (IARITHIF(iyi))
			{
				case -1: goto L_4050;
				case  0: goto L_4060;
				case  1: goto L_4130;
			}
		iyi += divasc.maxint - divasc.kordi;
		if (iyi >= 0)
			goto L_4130;
		/*       TEST IF INTEG TOO BIG FOR THIS EQUATION */
L_4050:
		if (divasc.kordi < -iyi)
			goto L_4130;
L_4060:
		iyi += iy;
		iyn = iyi + iyni;
		/*  COMPUTE INNER PRODUCT FOR EQUATION WHICH IS NOT STIFF */
		xp1 = C0;
		if (divasc.ldt == 2)
		{
			if (kqq != divasc.kqmaxi)
				xp1 = Csum[k + integz + ici]*F[idt + integz + divasc.numdt - 1];
		}
		for (j = k + integz + ici; j >= (ici + 1); j--)
		{
			xp1 += C[j]*F[idt - ici - 1 + j];
		}
		switch (IARITHIF(integ))
		{
			case -1: goto L_4080;
			case  0: goto L_4090;
			case  1: goto L_4100;
		}
		/* STORE FINAL RESULT IN Y WHEN DIFFERENTIATING */
L_4080:
		;
		Y[iyi] = xp1;
		goto L_4130;
		/* STORE INTERPOLATED VALUE IN F (OR STIFF DIFFERENTIATION) */
L_4090:
		F[i] = xp1;
		goto L_4130;
		/* PICK UP EXTRA STUFF TO ADD TO INNER PRODUCT WHEN INTEGRATING */
L_4100:
		k = ici;
		if (k == 0)
			goto L_4120;
L_4110:
		;
		xp1 = C[k]*(xp1 + Y[iyn]);
		iyn -= 1;
		k -= 1;
		if (k != 0)
			goto L_4110;
		/* STORE FINAL RESULT IN Y WHEN INTEGRATING (OR STIFF INTERPOLATION) */
L_4120:
		Y[iyi] = xp1 + Y[iyn];
L_4130:
		;
		idt += divasc.numdt;
	}
 
	integ -= 1;
	if (integ >= -divasc.maxdif)
		switch (IARITHIF(integ))
		{
			case -1: goto L_3990;
			case  0: goto L_3950;
			case  1: goto L_3920;
		}
L_4160:
	return;
	/* ********
	 * ERROR PROCESSING
	 * ******** */
L_4170:
	Mact[2] = 68;
	Mact[3] = 11;
	Mact[6] = LTXTAD;
	Idat[1] = divasc.ldt;
	goto L_4190;
L_4180:
	Mact[2] = 28;
	Mact[3] = 1;
	Mact[6] = LTXTAC;
	Fdat[2] = divasc.tn;
	Fdat[3] = T[2];
	Fdat[4] = Xi[1];
	Fdat[5] = divasc.tn - T[2];
	Fdat[6] = divasc.tn + C2*Xi[1];
	if (Xi[1] < 0)
	{
		Fdat[5] = Fdat[6];
		Fdat[6] = divasc.tn - T[2];
	}
L_4190:
	Fdat[1] = T[1];
	/*--D Next line special: P=>S, X=>D */
	dmess( mact, (char*)mtxtaa,155, idat, fdat );
	if (Mact[2] < 50)
		goto L_3820;
	return;
} /* end of function */
/*   End of DIVAIN */
 
		/* PARAMETER translations */
#define	CP625	.625e0
#define	CP9	.9e0
#undef	LTXTAC
#define	LTXTAC	49
		/* end of PARAMETER translations */
 
void /*FUNCTION*/ divaop(
long iopt[],
double fopt[])
{
	long int _l0, i, ia, j, k, multj;
	static long int liopt;
	static char mtxtaa[1][77]={"DIVAOP$BError in IOPT() specifications: IOPT =$EHMIN$ = $F is > HMAX = $F.$E"};
	static long mact[7]={MEEMES,88,24,0,MEIVEC,0,MERET};
	static long mact1[5]={MEEMES,28,24,LTXTAC,MERET};
	static long iopts[23]={0,0,0,500000,0,0,0,0,0,0,0,0,0,0,0,0,1,
	 0,0,0,0,0,0};
	static long incop[22]={1,3,2,2,2,2,2,1,2,3,1,2,1,1,1,3,2,2,2,2,
	 2,2};
	long int *const ioptc = (long*)((long*)&divamc.iop3 + -2);
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const Fopt = &fopt[0] - 1;
	long * __restrict const Incop = &incop[0] - 1;
	long * __restrict const Iopt = &iopt[0] - 1;
	long * __restrict const Ioptc = &ioptc[0] - 1;
	long * __restrict const Iopts = &iopts[0] - 1;
	long * __restrict const Mact = &mact[0] - 1;
	long * __restrict const Mact1 = &mact1[0] - 1;
	double * __restrict const Xi = &divasc.xi[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1987-12-07 DIVAOP Krogh   Initial code.
	 *
	 *  SUBROUTINE TO SET UP OPTIONS FOR DIFFERENTIAL EQUATION  PACKAGE -IVA */
 
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
 
	/*.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS. */
 
 
	/*                      Declarations for error message processing.
	 * */
 
	/* ********* Error message text ***************
	 *[Last 2 letters of Param. name]  [Text generating message.]
	 *AA DIVAOP$B
	 *AB Error in IOPT() specifications: IOPT =$E
	 *AC HMIN = $F is > HMAX = $F.$E */
	/* **** End of text generated by pmess
	 *                      1   2  3   4       5  6      7 */
 
	/*                      IOP4       IOP17 */
 
	/*                  1  2    3  8  9 10 11 12  13  16   17 21 22 */
 
	/* ********* START OF EXECUTABLE CODE ***********************
	 * */
	multj = 1;
	k = 1;
L_4200:
	i = Iopt[k];
	ia = labs( i );
	/* 1 and 6 lines below need 21 changed if more options are added. */
	if (ia <= 21)
		switch (IARITHIF(i))
		{
			case -1: goto L_4220;
			case  0: goto L_4520;
			case  1: goto L_4280;
		}
	if (ia != 1111)
		goto L_4490;
	if (i < 0)
	{
		multj = -1;
		k += 1;
		goto L_4200;
	}
	Iopt[2] = liopt;
 
	/*     ****  INITIALIZE FOR STARTING A NEW INTEGRATION */
	for (j = 3; j <= 23; j++)
	{
		Ioptc[j] = Iopts[j];
	}
	divamc.ksout = Iopts[4];
	divamc.kmark = 1 - Iopts[1];
	divasc.kordi = Iopts[17];
	divasc.nkdko = max( -divasc.kordi, 0 );
	divasc.iopst = Iopts[22];
	goto L_4260;
 
	/*     **** SET A NOMINAL VALUE */
L_4220:
	Iopts[ia] = 0;
	if (ia == 12)
		goto L_4420;
	switch (IARITHIF(ia - 2))
	{
		case -1: goto L_4400;
		case  0: goto L_4240;
		case  1: goto L_4230;
	}
L_4230:
	if (ia == 4)
		Iopts[4] = 500000;
	if (ia == 21)
		divamc.tolg = 0.e0;
	goto L_4390;
 
	/*     **** SET ALL OPTIONS TO THEIR NOMINAL VALUES */
L_4240:
	ia = 1;
	Iopts[1] = 0;
	for (j = 3; j <= 22; j++)
	{
		Iopts[j] = 0;
		Ioptc[j] = 0;
	}
	Iopts[4] = 500000;
	Ioptc[4] = Iopts[4];
	Iopts[17] = 1;
	divamc.tolg = 0.e0;
L_4260:
	divamc.ngtot = Iopts[7] + max( Iopts[6], 0 );
	if (Iopts[12] == 0)
		goto L_4420;
L_4270:
	return;
 
	/*     **** SET SPECIFIED OPTION */
L_4280:
	j = Iopt[k + 1];
	switch (IARITHIF(Incop[ia] - 2))
	{
		case -1: goto L_4290;
		case  0: goto L_4330;
		case  1: goto L_4300;
	}
	/*     **** OPTION INVOLVES NO EXTRA PARAMETERS */
L_4290:
	Iopts[ia] = 1;
	switch (IARITHIF(ia - 2))
	{
		case -1: goto L_4400;
		case  0: goto L_4400;
		case  1: goto L_4390;
	}
	/*     **** TAKE CARE OF SECOND EXTRA PARAMETER */
L_4300:
	if (ia != 10)
		goto L_4310;
	divamc.noutko = Iopt[k + 2];
	switch (IARITHIF(divamc.noutko))
	{
		case -1: goto L_4500;
		case  0: goto L_4350;
		case  1: goto L_4350;
	}
L_4310:
	if (ia != 16)
		goto L_4320;
	divamc.ntolf = Iopt[k + 2];
	switch (IARITHIF(divamc.ntolf))
	{
		case -1: goto L_4500;
		case  0: goto L_4500;
		case  1: goto L_4350;
	}
L_4320:
	if (j == 3)
	{
		if (divamc.kmark != 3)
		{
			if (Xi[1]*(Fopt[Iopt[k + 2]] - divamc.tmark) >= C0)
				goto L_4400;
		}
	}
	divamc.tmark = Fopt[Iopt[k + 2]];
	divamc.kmark = j;
	goto L_4400;
	/*     **** TAKE CARE OF FIRST EXTRA PARAMETER */
L_4330:
	;
	if (ia == 12)
		goto L_4410;
	if (ia == 4)
		divamc.ksout = j;
	if (ia == 21)
		divamc.tolg = Fopt[j];
L_4350:
	Iopts[ia] = j*multj;
	if (labs( ia - 7 ) > 1)
		goto L_4360;
	/*     **** SET SPECIAL PARAMETERS FOR GSTOP-S */
	divamc.igflg = 0;
	divamc.ngtot = Iopts[7] + max( Iopts[6], 0 );
	/*     **** TEST FOR ERROR */
	if (j > 500)
		goto L_4500;
L_4360:
	if (j > 0)
		goto L_4390;
	if ((ia == 5) || (ia == 17))
		goto L_4390;
	switch (IARITHIF(j + 1))
	{
		case -1: goto L_4500;
		case  0: goto L_4370;
		case  1: goto L_4380;
	}
L_4370:
	if (ia == 7)
		goto L_4500;
L_4380:
	if (((ia == 4) || (ia == 11)) || (ia >= 16))
		goto L_4500;
	/*     **** STORE SAVED VALUE IN COMMON */
L_4390:
	Ioptc[ia] = Iopts[ia];
 
	/*     **** INCREMENT K TO GET NEXT OPTION */
L_4400:
	k += Incop[ia];
	goto L_4200;
 
	/* ******* SET UP INFORMATION FOR CHANGING STEPSIZE *********
	 *
	 *     **** TEST IF VALUES ARE ALREADY SET */
L_4410:
	if (Iopts[12] != 0)
		goto L_4430;
	/*     **** SET NOMINAL VALUES FOR VARIABLES ONLY SET ONCE */
L_4420:
	divamc.erep = CP3;
	/*     **** SET NOMINAL VALUES FOR STEPSIZE CONTROL AND ENV. CONSTANTS */
	divaev.eeps2 = DBL_EPSILON;
	divaev.eeps16 = C16*divaev.eeps2;
	divaev.eeps10 = CP625*divaev.eeps16;
	divaev.eept75 = pow(divaev.eeps2,CP75);
	divaev.eeps2 += divaev.eeps2;
	divaev.ovd10 = DBL_MAX;
	divaev.erov10 = C10/divaev.ovd10;
	divaev.eovep2 = divaev.ovd10*divaev.eeps2;
	divaev.ovtm75 = pow(divaev.ovd10,CMP75);
	divaev.ovd10 /= C10;
	divamc.hinc = C2;
	divamc.hdec = CP5;
	divamc.hmin = divaev.erov10;
	divamc.hmax = divaev.ovd10;
	if (i != 12)
		goto L_4470;
L_4430:
	Iopts[12] = j*multj;
	switch (IARITHIF(j))
	{
		case -1: goto L_4450;
		case  0: goto L_4470;
		case  1: goto L_4460;
	}
	/*     **** SET UP TO GIVE USER COMPLETE STEPSIZE CONTROL */
L_4450:
	divamc.erep = C1/divaev.erov10;
	divamc.hinc = -C2;
	divamc.iop8 = 1;
	/*## Recent code 12/16/94 */
	divamc.lincd = -2;
	divamc.linc = -2;
	divamc.robnd = C1;
	/*## Recent code 9/6/2001 */
	divamc.tolg = 0.e0;
	/*## End of recent code */
	goto L_4480;
	/*     **** SET USER VALUES FOR STEPSIZE CONTROL */
L_4460:
	if (Fopt[j] != C0)
		divamc.hinc = fmax( C1P125, fmin( Fopt[j], C4 ) );
	if (Fopt[j + 1] != C0)
		divamc.hdec = fmin( CP875, fmax( Fopt[j + 1], CP25 ) );
	if (Fopt[j + 2] != C0)
		divamc.hmin = Fopt[j + 2];
	if (Fopt[j + 3] != C0)
		divamc.hmax = Fopt[j + 3];
	if ((divamc.hmin > divamc.hmax) || (divamc.hmax <= 0.e0))
	{
		dmess( mact1, (char*)mtxtaa,77, iopt, &Fopt[j + 2] );
		divamc.kord2i = -4;
	}
L_4470:
	divamc.kqicon = -1;
L_4480:
	divamc.hmaxp9 = divamc.hmax*CP9;
	switch (IARITHIF(i - 1111))
	{
		case -1: goto L_4400;
		case  0: goto L_4270;
		case  1: goto L_4400;
	}
 
	/* ***************** ERROR  IN  IOPT ************************
	 * */
L_4490:
	ia = 1;
L_4500:
	Mact[6] = k + Incop[ia] - 1;
	mess( mact, (char*)mtxtaa,77, iopt );
	divamc.kord1i = 24;
	divamc.kord2i = -4;
	/* Usual return with no error is here. */
L_4520:
	liopt = k;
	return;
} /* end of function */
/*   End of DIVAOP */
 
void /*FUNCTION*/ divapr(
double y[],
double yn[],
double f[],
long kord[])
{
	long int i, integ, j, k, kqq, l, n;
	double temp[KDIM], tp1, xp;
	static long integs = -1;
		/* OFFSET Vectors w/subscript range: 1 to dimension */
	double * __restrict const Beta = &divamc.beta[0] - 1;
	double * __restrict const F = &f[0] - 1;
	long *  __restrict const Kord = &kord[0] - 1;
	double * __restrict const Temp = &temp[0] - 1;
	double * __restrict const Y = &y[0] - 1;
	double * __restrict const Yn = &yn[0] - 1;
		/* end of OFFSET VECTORS */
 
	/*>> 1988-01-13 DIVAPR Krogh   Initial code.
	 *
	 * THIS SUBROUTINE
	 *   1. UPDATES THE DIFFERENCE TABLE FROM THE PREVIOUS STEP (IF NOT
	 *      DONE ALREADY).
	 *   2. PREDICTS WHAT THE VALUES OF THE DEPENDENT VARIABLES, Y, AND
	 *      THE DIFFERENCE TABLE, DT, WILL BE AT THE END OF THE CURRENT STEP.
	 *
	 *   Y = VECTOR OF PREDICTED VALUES COMPUTED BY THIS SUBROUTINE.
	 *   YN= VECTOR OF VALUES OF Y COMPUTED ON THE LAST STEP.
	 *   F = VECTOR OF DERIVATIVE VALUES.
	 *   DT= ARRAY CONTAINING DIFFERENCE TABLES.
	 *   KD= VECTOR GIVING ORDERS OF THE DIFFERENTIAL EQUATIONS (IF
	 *       EQUATIONS HAVE DIFFERENT ORDERS).
	 *   KQ= VECTOR OF INTEGRATION ORDERS.
	 * */
	/*--D Next line special: P=>D, X=>Q */
 
	/*++ Substitute for KDIM, MAXORD, MAXSTF below */
	/*--D Next line special: P=>D, X=>Q */
 
	/*--D Next line special: P=>D, X=>Q */
 
 
	/*.    SPECIFICATION OF ENVIRONMENTAL CONSTANTS. */
 
 
	/*--D Next line special: P=>D, X=>Q */
 
	/*++  Code for ARGM is inactive
	 *      RETURN
	 *      ENTRY DIVAPE
	 *++  End
	 * ********
	 * START OF CODE
	 * ******** */
	divamc.iy = 0;
	l = divasc.ndtf - 1;
	for (i = 1; i <= divasc.nte; i++)
	{
		integ = divasc.kordi;
		if (divasc.nkdko != 0)
			integ = Kord[divasc.nkdko + i - 1];
		kqq = Kord[i + 3];
		k = max( labs( kqq ), 2 );
		switch (IARITHIF(kqq))
		{
			case -1: goto L_4530;
			case  0: goto L_4520;
			case  1: goto L_4540;
		}
L_4520:
		divamc.iy += labs( integ );
		goto L_4670;
		/* ********
		 * EQUATION IS STIFF, OR IMPLICIT
		 * ******** */
L_4530:
		;
		/*++  Code for STIFF is inactive
		 *      KQQ=-KQQ
		 *      N=KQQ-1
		 *      JS=abs(KORD(NJSKO+I-1))-1
		 *      IMPLIC=INTEG
		 *      INTEG=abs(IMPLIC)-JS
		 *c.    SET INTEGS FOR STIFF EQUATIONS
		 *      INTEGS=0
		 *      IF (K-KSC) 160,160,140
		 *c.END OF SET-UP FOR STIFF EQUATIONS
		 *++  End
		 * ********
		 * EQUATION IS NOT STIFF
		 * ******** */
L_4540:
		n = kqq;
		if (divasc.ldt != 0)
			switch (IARITHIF(k - divamc.ksc))
			{
				case -1: goto L_4570;
				case  0: goto L_4570;
				case  1: goto L_4550;
			}
		/*     DIFFERENCE TABLE HAS NOT BEEN UPDATED */
		tp1 = F[i] - F[l + 1];
		switch (IARITHIF(k - divamc.ksc))
		{
			case -1: goto L_4610;
			case  0: goto L_4610;
			case  1: goto L_4590;
		}
		/* END OF SET-UP FOR EQUATIONS WHICH ARE NOT STIFF
		 * ********
		 * GET PREDICTED DIFFERENCES FROM UPDATED DIFFERENCE TABLE
		 * ******** */
L_4550:
		F[l + k + 1] *= Beta[k + 1];
		Temp[k] = F[l + k]*Beta[k];
		F[l + k] = Temp[k];
		/* LOOP FOR MODIFIED DIVIDED DIFFERENCES */
L_4560:
		k -= 1;
		if (k <= divamc.ksc)
			goto L_4580;
		Temp[k] = F[l + k]*Beta[k];
		F[l + k] = Temp[k] + F[l + k + 1];
		goto L_4560;
		/* CODE FOR BACKWARD DIFFERENCES */
L_4570:
		F[l + k + 1] = F[l + k + 1];
		Temp[k] = F[l + k];
		k -= 1;
 
L_4580:
		Temp[k] = F[l + k];
		F[l + k] = Temp[k] + F[l + k + 1];
		k -= 1;
		if (k != 0)
			goto L_4580;
		goto L_4630;
		/* ********
		 * UPDATE DIFFERENCE TABLE AND GET PREDICTED DIFFERENCES
		 * ********
		 * CODE FOR MODIFIED DIVIDED DIFFERENCES */
L_4590:
		F[l + k + 1] = (F[l + k + 1] + tp1)*Beta[k + 1];
		Temp[k] = (F[l + k] + tp1)*Beta[k];
		F[l + k] = Temp[k];
L_4600:
		k -= 1;
		if (k <= divamc.ksc)
			goto L_4620;
		Temp[k] = (F[l + k] + tp1)*Beta[k];
		F[l + k] = Temp[k] + F[l + k + 1];
		goto L_4600;
		/* CODE FOR BACKWARD DIFFERENCES */
L_4610:
		F[l + k + 1] = F[l + k + 1] + tp1;
		Temp[k] = F[l + k] + tp1;
		F[l + k] = Temp[k];
		k -= 1;
 
L_4620:
		Temp[k] = F[l + k] + tp1;
		F[l + k] = Temp[k] + F[l + k + 1];
		k -= 1;
		if (k != 0)
			goto L_4620;
		/* ********
		 * COMPUTE Y-S OBTAINED USING INTEGRATION
		 * ********
		 *     TEST IF NEXT Y TO BE OBTAINED BY INTERPOLATION */
L_4630:
		;
		/*++  Code for STIFF is inactive
		 *      IF (INTEG.EQ.0) GO TO 4662
		 *++  End */
		divamc.iy += 1;
		/*     FORM INNER PRODUCT */
		xp = C0;
		for (j = integs + n + 1; j >= (integs + 2); j--)
		{
			/*++  Code for ~{p,x} is active */
			xp += divamc.g[integ - 1][j - 1]*Temp[j];
			/*++  Code for {p,x} is inactive
			 *c--D Next line special: P=>D, X=>Q
			 *            XP = XP + dble(G(J, INTEG)) * dble(TEMP(J))
			 *++  END */
		}
		k = integ + integs;
		for (j = k; j >= 1; j--)
		{
			/*++  Code for ~{p,x} is active */
			xp += divamc.g[j - 1][0]*Yn[divamc.iy + j];
			/*++  Code for {p,x} is inactive
			 *c--D Next line special: P=>D, X=>Q
			 *            XP = XP + dble(G(1, J)) * dble(YN(IY + J))
			 *++  END */
		}
		Y[divamc.iy] = Yn[divamc.iy] + xp;
		integ -= 1;
		switch (IARITHIF(k))
		{
			case -1: goto L_4670;
			case  0: goto L_4670;
			case  1: goto L_4630;
		}
		/* END OF COMPUTING Y-S OBTAINED BY INTEGRATION
		 * ********
		 * COMPUTE Y-S OBTAINED USING INTERPOLATION AND DIFFERENTIATION
		 * ********
		 *++  Code for STIFF is inactive
		 *c.    RESTORE INTEGS FOR EQUATIONS WHICH ARE NOT STIFF
		 * 4662 INTEGS=-1
		 *      IY=IY+1
		 *c.    COMPUTE Y USING INTERPOLATION
		 *      Y(IY)=YN(IY) + F(L+2)
		 *      IF (KQQ.EQ.1) Y(IY)=YN(IY)
		 * 4663 INTEG=INTEG+1
		 *      IF (INTEG.EQ.JS) IF (IMPLIC) 4680,4680,4664
		 *c.    COMPUTE INTEG-TH DERIVATIVE
		 *      XP = C0
		 * 4664 DO 4666 J = KQQ+1, INTEG+1, -1
		 *         XP = XP + D(J, INTEG) * TEMP(J)
		 * 4666 CONTINUE
		 *      IF (INTEG.EQ.JS) GO TO 4667
		 *      IY=IY+1
		 *      Y(IY)=XP
		 *      GO TO 4663
		 *c.STORE PREDICTED VALUE FOR F
		 * 4667 CONTINUE
		 *      F(L+NUMDT)=XP
		 *++  End */
L_4670:
		l += divasc.numdt;
	}
	divasc.ldt = -3;
	return;
} /* end of function */
 
