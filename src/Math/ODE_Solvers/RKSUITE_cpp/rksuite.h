#ifndef RKSUITE_H
#define RKSUITE_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <malloc.h>
/****************************************************************************
 *  RKSUITE/C++ - Conversion of RKSUITE to C++
 *
 *  Conversion from FORTRAN to C by anonymous
 *  Conversion from C to C++ by
 *      D. F. Linton
 *      Millennium Engineering and Integration Company (MEI)
 *      600 Jackson Court
 *      Satellite Beach, FL 32937 USA
 *
 *  Changes are Open Source under the Lesser GPL.
 *
 *  Not warrantied by MEI or the author.  Use only after testing to confirm
 *  suitability for your application.
 *
 *  Please read http://www.netlib.org/ode/rksuite/rksuite.doc before using.
 *
 *  If you are not familiar with the routines of RKSUITE and how they are used
 *  to solve initial value problems, you should study the document rksuite.doc 
 *  carefully before proceeding further.  The following function documentation
 *  are simply "Brief Reminders" intended only to remind you of the meaning,
 *  type, and size requirements of the arguments.
 *
 ****************************************************************************
 */

///////////////////////////////////////////////////////////////////////////////
/// @class RKSUITE 
///
/// Used to create an integrator object.
///
/// All intermediate data is stored witin the object.
///
///////////////////////////////////////////////////////////////////////////////

class RKSUITE {
public:

    /** 
    Create an instance of an RKSUITE integrator.  The setup() method must be
    called to initialize the integrator before use.
    */
    RKSUITE() {
        work = NULL;
        lenint = 0;
        wrkint = NULL;
    }

    /**
    Destroy this instance of an HDSAPI application.
    */
    ~RKSUITE() {
        if (work != NULL)
            delete[] work;
        if (wrkint != NULL)
            delete[] wrkint;
    }

    /**
    Configures and initializes this integrator object.
    @param neq - The number of differential equations in the system.
        Constraint: neq >= 1
    @param tstart - The initial value of the independent variable.
    @param ystart - The vector of initial values of the solution components.
    @param tend - The integration proceeds from tstart in the direction of
        tend. You cannot go past tend without calling reset() or setup() again.
    @param tol - The relative error tolerance.
        Constraint: 0.01D0 >= tol >= 10*MCHEPS
    @param thres - threshold for the solution components.
        Constraint: thres(L) >= SQRT(DWARF)
    @param method - Specifies which Runge-Kutta pair is to be used.
        1 - use the (2,3) pair, 2 - use the (4,5) pair, 3 - use the (7,8) pair
    @param task - 'U' to use ut() or 'C' to use ct()
    @param errass - bool, true to attempt to assess the true error.
    @param hstart - if non-zero, use hstart as the first step size
    @param mesage - bool, true if you want informative messages written to 
        standard output
    */
    void setup(int neq, double tstart, double ystart[], double tend,
                        double tol, double thres[], int method, char* task,
                        bool errass, double hstart, bool mesage) {
        this->neq = neq;
        if (work != NULL)
            delete[] work; 
        work = new double[32*neq];
        if (wrkint != NULL)
            delete[] wrkint;
        lenint = 0;
        setup( neq, tstart, ystart, tend, tol, thres, method, task, errass, hstart,
            work, 32*neq, mesage );
    };

    /**
    UT - the "usual" task.  If all you want is to step to a time call this.
    @param f - IN derivatives function (double tin, double yin[neq], double ypout[neq])
    @param twant - IN next value of the independent variable
    @param tgot - OUT acheived value of the independent variable
    @param ygot - OUT approximate solution at tgot
    @param ypgot - OUT approximate first derivatives at tgot
    @param ymax - IN/OUT largest magnitude of ygot during integration.  Do not alter.
    @param uflag - OUT status, (1 - success, 2 - inefficient usage, 3 - work intensive,
        4 - problem probably stiff, 5 - too much accuracy requested, 6 - global error 
        assessment unreliable beyond this point, 911 - catastrophic failure reported on
        stdout.
    */
    void ut(void (*f)(double, double*, double*), double twant, double& tgot,
        double ygot[], double ypgot[], double ymax[], int& uflag) {
        ut( f, twant, tgot, ygot, ypgot, ymax, work, uflag );
    }

    /**
    Called to obtain some details about the integration.
    @param totfnc - OUT total number of calls to the derivatives function
    @param stpcst - OUT cost of a step in derivatives function calls
    @param waste - OUT number of failed steps
    @param stpsok - OUT number of accepted steps
    @param hnext - step size the integrator plans to use for the next step
    */
    void stat(int& totfcn, int& stpcst, double& waste, int& stpsok, double& hnext);

    /**
    If ERRASS was set .TRUE. in the call to SETUP, then after any call to UT
    or CT to advance the integration to TNOW or TWANT, the subroutine GLBERR
    may be called to obtain an assessment of the true error of the integration.
    At each step and for each solution component Y(L), a more accurate "true"
    solution YT(L), an average magnitude "size(L)" of its size, and its error
                abs(Y(L) - YT(L))/max("size(L)",THRES(L))
    are computed.  The assessment returned in RMSERR(L) is the RMS (root-mean-
    square) average of the error in the Lth solution component over all steps
    of the integration from TSTART through TNOW.    
    @param rmserr - OUT approximate RMS error of solution elements
    @param errmax - OUT maximum (approximate) true error taken over all
        solution components and all steps from TSTART to TNOW.
    @param terrmx - OUT First value of the independent variable where the
        (approximate) true error attains the maximum value ERRMAX.
    */
    void glberr(double rmserr[], double& errmax, double& terrmx) {
        glberr( rmserr, errmax, terrmx, work );
    }

    /**
    CT - the "complex" task.  If need to check for events along the way call this.
    @param f - IN derivatives function (double tin, double yin[neq], double ypout[neq])
    @param tnow - OUT next value of the independent variable
    @param ynow - OUT next value of the solution vector
    @param ypnow - OUT first derivatives at tnow
    @param cflag - OUT status, (1 - success, 2 - inefficient usage, 3 - work intensive,
        4 - problem probably stiff, 5 - too much accuracy requested, 6 - global error 
        assessment unreliable beyond this point, 911 - catastrophic failure reported on
        stdout.
    */
    void ct(void (*f)(double, double*, double*), double& tnow, double ynow[],
        double ypnow[], int& cflag) {
            ct( f, tnow, ynow, ypnow, work, cflag );
    }

    /**
    Obtain state and derivative values at intermediate points during the current step.
    @param twant - IN next value of the independent variable
    @param reqest - IN 'S' compute solution only, 'D' compute derivatives only, 'B' both
    @param nwant - IN number of solution vector components to interpolate (generally neq)
    @param ywant - OUT requested solution vector component elements
    @param ypwant - OUT requested solution derivative vector component elements
    @param f - IN derivatives function (double tin, double yin[neq], double ypout[neq])
    */
    void intrp(double twant, char reqest[], int nwant, double ywant[], double ypwant[],
        void (*f)(double, double*, double*)) {
            int length = neq + ((neq > 5*nwant) ? neq : 5*nwant);
            if (length != lenint) {
                if (wrkint != NULL)
                    delete[] wrkint;
                lenint = length;
                wrkint = new double[lenint];
            }
            intrp( twant, reqest, nwant, ywant, ypwant, f, work, wrkint, lenint );
    }

    /**
    Change the ending value of the independent variable
    @param tendnu - new final value
    */
    void reset(double tendnu);

    /**
    Prevent a program STOP after a "catastrophic" failure when using a routine from RKSUITE.
    When a "catastrophic" failure is detected, the default action of
    RKSUITE is to write an explanation to the standard output, and STOP.  
    This method can be used to prevent the STOP and so allow the program
    to continue.  To do this, you call SOFTFL with ON = true.  You must then 
    call the method chkfl() after every call to an RKSUITE method to check whether
    a catastrophic error occurred and take appropriate action if it did.  Of
    course, you may call setup() at any time to start a new problem, but calling
    any other user-callable routine in RKSUITE after a catastrophic error will
    lead to a STOP (even when "soft failure" has been set "on").

    @param on - true to prevent stops
    */
    void softfl(bool& on) {
        softfl( false, on );
    }

    /**
    Check if a catastrophic error has occured.
    @param error - OUT bool, true if fatal error occured.
    */
    void chkfl(bool& error) {
        chkfl( true, error );
    }

protected:

    void setup(int neq, double tstart, double ystart[], double tend,
                        double tol, double thres[], int method, char* task,
                bool errass, double hstart, double work[], int lenwrk,
                bool mesage);
    void ut(void (*f)(double, double*, double*), double twant, double& tgot,
        double ygot[], double ypgot[], double ymax[], double work[], int& uflag);
    void ct(void (*f)(double, double*, double*), double& tnow, double ynow[],
	    double ypnow[], double work[],int& cflag);
    void intrp(double twant, char reqest[], int nwant, double ywant[], double ypwant[],
	    void (*f)(double, double*, double*), double work[], double wrkint[],
        int lenint);
    void glberr(double rmserr[], double& errmax, double& terrmx, double work[]);
    void envirn(int& outch, double& mcheps, double& dwarf);
    void softfl(bool ask, bool& on);
    void chkfl(bool ask, bool& error);

    void mconst(int method);
    void rkconst(int method, int& vecstg, bool& reqstg, int& lintpl);
    void rkmsg(int ier, const char* srname, int nrec, int& flag);
    void rksit(bool ask, const char* srname, int& state);
    void step(void (*f)(double, double*, double*), int neq, double tnow,
	    double* y, double* yp, double stages /* [neq][] */ [], double tol, double& htry,
	    double* weight, double* ynew, double* errest, double& err, bool main,
	    double hmin, double* thres, bool& phase2); 
    void stepa(double tnow, double y[], double yp[], double tstg, double ystg[],
	    double ypstg[], double& htry, double weight[], bool& cutbak);
    void stepb(int neq, double y[], double yp[], double h, double ynew[],
	    double stages/*[neq][]*/[], double thres[], double& err, bool main,
	    double weight[]);
    void stiff(void (*f)(double, double*, double*), double havg, int& jflstp,
	    bool toomch, int maxfcn, double work[], int& ier, int& nrec);
    void stiffa(void (*f)(double, double*, double*),double x, double y[],
	    double hnow, double havg, double xend, int maxfcn, double wt[],
	    double fxy[], double v0[], bool& unsure, bool& stif, double v1[],
	    double v2[], double v3[], double vtemp[]);
    void stiffb(double v1v1, double v0v1, double v0v0, double& rold,
	    double& rho, double root1[], double root2[], bool& rootre);
    void stiffc(double alpha, double beta, double r1[], double r2[]);
    void stiffd(double v[], double havg, double x, double y[],
	    void (*f)(double, double[], double[]), double fxy[], double wt[],
	    double scale, double vdotv, double z[], double& zdotz, double vtemp[]);
    double dotprd(double u[], double v[], double wt[], int neq);
    void evali(double y[], double yp[], double p /* [nwant][] */ [], double twant,
	    char reqest[], int nwant, double ywant[], double ypwant[]);
    void formi(void (*f)(double, double*, double*), int neq, int nwant, double y[], double yp[],
	    double yold[], double ypold[], double stages/*[neq][]*/[], bool calstg,
    double xstage[], double p[]);
    void truerr(void (*f)(double, double*, double*), int neq, double y[],
	    double tol, double weight[], double zy[], double zyp[], double zerror[],
	    double zynew[], double zerres[], double zstage /* [neq][] */ [], int& ier);

    double*    work;
    int        neq;
    int        lenint;
    double*    wrkint;

    struct
    {
    double tstrt, tnd, dir, hstrt, tolr;
    int neqn;
    } rkcom1;

    struct
    {
    double t, h, told, hold;
    int nfcn, svnfcn, okstp, flstp;
    bool first, last;
    } rkcom2;

    struct
    {
    int prthrs, prerst, prwt, pryold, prscr,
            pry, pryp, prstgs, printp, lnintp;
    } rkcom3;

    struct
    {
    double a[14][14], r[12][7];
    double b[14], c[14], bhat[14], e[8];
    int ptr[14];
    int nstage, methd, mintp;
    bool intp;
    } rkcom4;

    struct
    {
    double toosml, cost, safety, expon, stbrad, tanang,
            rs, rs1, rs2, rs3, rs4;
    int order, lststg, maxtry, nsec;
    bool fsal;
    } rkcom5;

    struct
    {
    double maxerr, locmax;
    int gnfcn, przstg, przy, przyp, przers, przerr, przynu;
    bool erason, erasfl;
    } rkcom6;

    struct
    {
    double mcheps, dwarf, rndoff, sqrrmc, cubrmc, tiny;
    int outch;
    } rkcom7;

    struct
    {
    bool msg, utask;
    } rkcom8;

    struct
    {
        char rec[10][80];
    } rkcom9;

};
#endif
