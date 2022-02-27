/*
 *  cuLsoda.h
 *
 *	File Notes:	This file is a conversion of the REAL precision Livermore Solver for
 *	Ordinary Differential Equations with automatic switching for stiff and non-stiff
 *	problems (DLSODA)
 */

/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."
 
 - From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */
#ifndef GMS_CULSODA_CU_H
#define GMS_CULSODA_CU_H

#define REAL double


/* Common Block Declarations */
struct cuLsodaCommonBlock
{
	REAL /*rowns[209],*/ CM_conit, CM_crate, CM_ccmax, CM_el0, CM_h__, CM_hmin, CM_hmxi, CM_hu, CM_rc, CM_tn, CM_uround, CM_pdest, CM_pdlast, CM_ratio, CM_hold, CM_rmax;
	REAL  CM_el[13], CM_elco[156]	/* was [13][12] */, CM_tesco[36]	/* was [3][12] */;
	REAL CM_rls[218];
	REAL CM_tsw, /*rowns2[20],*/ CM_pdnorm;
	REAL /*rownd2,*/ CM_cm1[12], CM_cm2[5];
	REAL CM_rlsa[22];
	REAL CM_sm1[12];
	int CM_init, CM_mxstep, CM_mxhnil, CM_nhnil, CM_nslast, CM_nyh, /*iowns[6],*/ CM_icf, 
	CM_ierpj, CM_iersl, CM_jcur, CM_jstart, CM_kflag, CM_l, CM_lyh, CM_lewt, CM_lacor, CM_lsavf,
	CM_lwm, CM_liwm, CM_meth, CM_miter, CM_maxord, CM_maxcor, CM_msbp, CM_mxncf, CM_n, CM_nq, 
	CM_nst, CM_nfe, CM_nje, CM_nqu;
	int /*iownd[6],*/ CM_ialth, CM_ipup, CM_lmax, /*meo,*/ CM_nqnyh, CM_nslp;
	int CM_ils[37];
	int CM_insufr, CM_insufi, CM_ixpr, /*iowns2[2],*/ CM_jtyp, CM_mused, CM_mxordn, CM_mxords; 
	int /*iownd2[3],*/ CM_icount, CM_irflag;
	int CM_ilsa[9];
};

/* End Common Block */ 


#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


//#define Fex_and_Jex_definition
//struct myFex
//{
//	__device__ void operator()(int *neq, REAL *t, REAL *y, REAL *ydot/*, void *otherData*/)
//	{
//		ydot[0] = (REAL)1.0E4 * y[1] * y[2] - (REAL).04E0 * y[0];
//		ydot[2] = (REAL)3.0E7 * y[1] * y[1];
//		ydot[1] = (REAL)-1.0 * (ydot[0] + ydot[2]);
//	}
//};

//struct myJex
//{
//	__device__ void operator()(int *neq, REAL *t, REAL *y, int ml, int mu, REAL *pd, int nrowpd/*, void *otherData*/)
//	{
//		return;
//	}
//};




/* dlsoda.f -- translated by f2c (version 20090411).
 You must link the resulting object file with libf2c:
 on Microsoft Windows system, link with libf2c.lib;
 on Linux or Unix systems, link with .../path/to/libf2c.a -lm
 or, if you install libf2c.a in a standard place, with -lf2c -lm
 -- in that order, at the end of the command line, as in
 cc *.o -lf2c -lm
 Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
 
 http://www.netlib.org/f2c/libf2c.zip
 */




/*template<typename Fex, typename Jex>
__device__ int dlsoda_(Fex, int *, REAL *, REAL *, REAL *, int *, REAL *, REAL *, int *, int *, int *, REAL *, int *, int *, int *, Jex, int *, struct cuLsodaCommonBlock *);

template<typename Fex, typename Jex> 
__device__ int dstoda_(int *neq, REAL *y, REAL *yh, int *NOT_nyh, REAL *yh1, REAL *ewt, REAL *savf, REAL *acor, REAL *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

template<typename Fex, typename Jex> 
__device__ int dprja_(int *neq, REAL *y, REAL *yh, int *NOT_nyh, REAL *ewt, REAL *ftem, REAL *savf, REAL *wm, int *iwm, Fex f, Jex jac, struct cuLsodaCommonBlock *common);

__device__ int dsolsy_(REAL *wm, int *iwm, REAL *x, REAL *tem, struct cuLsodaCommonBlock *common);
__device__ int dintdy_(REAL *t, int k, REAL *yh, int *NOT_nyh, REAL *dky, int *iflag, struct cuLsodaCommonBlock *common);
__device__ int dcfode_(int meth, REAL *DCFODE_elco, REAL *DCFODE_tesco, struct cuLsodaCommonBlock *common);
__device__ int dsolsy_(REAL *wm, int *iwm, REAL *x, REAL *tem, struct cuLsodaCommonBlock *common);
__device__ int dewset_(int *PARAM_n, int *itol, REAL *rtol, REAL *atol, REAL *ycur, REAL *ewt, struct cuLsodaCommonBlock *common);
__device__ REAL dmnorm_(int *PARAM_n, REAL *v, REAL *w, struct cuLsodaCommonBlock *common);
__device__ REAL dfnorm_(int *PARAM_n, REAL *a, REAL *w, struct cuLsodaCommonBlock *common);
__device__ REAL dbnorm_(int *PARAM_n, REAL *a, int *nra, int *ml, int *mu, REAL *w, struct cuLsodaCommonBlock *common);
__device__ int dsrcma_(REAL *rsav, int *isav, int *job, struct cuLsodaCommonBlock *common);
__device__ int dgefa_(REAL *a, int *lda, int *PARAM_n, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
__device__ int dgesl_(REAL *a, int *lda, int *PARAM_n, int *ipvt, REAL *b, int job, struct cuLsodaCommonBlock *common);
__device__ int dgbfa_(REAL *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, int *info, struct cuLsodaCommonBlock *common);
__device__ int dgbsl_(REAL *abd, int *lda, int *PARAM_n, int *ml, int *mu, int *ipvt, REAL *b, int job, struct cuLsodaCommonBlock *common);
__device__ REAL dumach_( struct cuLsodaCommonBlock *common);
//__device__ int xsetf_(int *mflag, struct CommonBlock *common);
//__device__ int xsetun_(int *lun, struct CommonBlock *common);
__device__ int ixsav_(int ipar, int *ivalue, int iset, struct cuLsodaCommonBlock *common);
__device__ int idamax_(int *PARAM_n, REAL *dx, int incx, struct cuLsodaCommonBlock *common);
__device__ int daxpy_(int *PARAM_n, REAL *da, REAL *dx, int incx, REAL *dy, int incy, struct cuLsodaCommonBlock *common);
__device__ int dumsum_(REAL a, REAL b, REAL *c__, struct cuLsodaCommonBlock *common);
__device__ int dscal_(int *PARAM_n, REAL *da, REAL *dx, int incx, struct cuLsodaCommonBlock *common);
__device__ REAL ddot_(int *PARAM_n, REAL *dx, int incx, REAL *dy, int incy, struct cuLsodaCommonBlock *common);
__device__ REAL d_sign(REAL *a, REAL *b);
__host__ __device__ void cuLsodaCommonBlockInit(struct cuLsodaCommonBlock *common);

*/



#endif



