#include <math.h>
extern "C" {
#include "rksuite.h"
}

#define N_HALTS   1 // number of halt tests
#define N_ELEMENTS  6  // number of state vector elements

// helper functions for exponentiation to integer powers
#define SQR(x) ((x)*(x))
#define P2(x) ((x)*(x))
#define P3(x) (P2(x)*(x))
#define P4(x) (P2(x)*P2(x))
#define P5(x) (P3(x)*P2(x))
#define P7(x) (P5(x)*P2(x))

/*---------------------------------------------------------------------
 *
 *  getGravity : Compute state vector derivatives
 *  
 *	getGravity( gravity[] , efg[] )
 *  
 *  efg[] - IN double ecef position[3] [m]
 *  gravity[] - OUT double ecef gravity-vector[3] [m/s/s]
 */
static
void getGravity( double gravity[3], double efg[3] ) {
	double r = (double) sqrt( SQR(efg[0]) + SQR(efg[1]) + SQR(efg[2]) );
	const double  a = 6378137;
	const double  b = 6.356752314245164e+006;
	const double  J2 = +1.082629989000000e-003;
	const double  mu = 3.986004418000000e+014;
	
	gravity[0] = (double) (efg[0]/P3(r));
	gravity[1] = (double) (efg[1]/P3(r));
	gravity[2] = (double) (efg[2]/P3(r));
	gravity[0] = gravity[0] - (double) (3.0/2.0*P2(a)*J2*(5.0*efg[0]*P2(efg[2])/P7(r)-efg[0]/P5(r)));
	gravity[1] = gravity[1] - (double) (3.0/2.0*P2(a)*J2*(5.0*efg[1]*P2(efg[2])/P7(r)-efg[1]/P5(r)));
	gravity[2] = gravity[2] - (double) (3.0/2.0*P2(a)*J2*(5.0*P3(efg[2])/P7(r)-3.0*efg[2]/P5(r)));
	gravity[0] *= (double) -mu;
	gravity[1] *= (double) -mu;
	gravity[2] *= (double) -mu;
}


/*---------------------------------------------------------------------
 *
 *  computeDerivatives : Compute state vector derivatives
 *  
 *	computeDerivatives( t , y[] , f[] )
 *  
 *  t - IN double  time [s]
 *  y[] - IN double ecef position[3], velocity[3] [m, m/s]
 *  f[] - OUT double ecef dposition/dt[3], dvelocity/dt[3] [m/s,m/s/s]
 */
static
void computeDerivatives(
				  double t,
				  double y[N_ELEMENTS],
				  double f[N_ELEMENTS] ) {

	double gravity[3];
	getGravity( gravity, y );

	f[0] = y[3];
	f[1] = y[4];
	f[2] = y[5];
	f[3] = gravity[0];
	f[4] = gravity[1];
	f[5] = gravity[2];
}

/*---------------------------------------------------------------------
 *
 *  computeHalts : Compute integration halting conditions
 *  
 *	computeHalts( t , y[] , f[] )
 *  
 *  t - IN double  time [s]
 *  y[] - IN double ecef position[3], velocity[3] [m, m/s]
 *  f[] - OUT double ecef halting function values (stops when zero)
 */
static                  
void computeHalts( double    t,
                   double    y[N_ELEMENTS],
                   double    G[N_HALTS] ) {
	G[0] = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]) - 6378140.0;
}


// This example test driver integrates a simple ballistic fall in vacuum.
// First do one second steps for the first 5 seconds and then integrate
// to ground impact.
// EXPECTED OUTPUT:
//    1.0  6379135.1     1000.0        0.0       -9.8     1000.0        0.0
//    2.0  6379120.4     2000.0        0.0      -19.6     1000.0        0.0
//    3.0  6379095.8     3000.0        0.0      -29.4     1000.0        0.0
//    4.0  6379061.5     4000.0        0.0      -39.2     1000.0        0.0
//    5.0  6379017.4     5000.0        0.0      -49.1     1000.0        0.0
//   14.4  6378123.8    14391.8        0.0     -141.2      999.8        0.0

int main( int argc, char* argv[] ) {
    RKSUITE    rksuite;  // create the integrator object instance
	double tnow = 0.0;     
    // current time
	double y[N_ELEMENTS] = { 6378140.0+1000.0, 0.0, 0.0, 0.0, 1000.0, 0.0 };
    // state vector
	double yp[N_ELEMENTS];  // storage for prior state vector 

	double TOL = 1e-4;  // relative tolerance
    double thres[N_ELEMENTS] = { 1e-2, 1e-2, 1e-2, 1e-4, 1e-4, 1e-4 };
    // absolute tolerances

	int    method = 2;  // RK(4,5)
	
    double stopTime = 1.0;
	int   cflag;  // status return 

    // first setup up the integrator
	rksuite.setup( 6, tnow, y, stopTime, TOL, thres, 
		method, "CT", false, 0.0f, false );

    // next take 1.0 second steps for five seconds

    while (tnow < 5.0) {
        // the ct() function does not guarantee to advance all the 
        // way to the stop time.  Keep stepping until it does.
        do {
	        computeDerivatives( tnow, y, yp );
	        rksuite.ct( computeDerivatives, tnow, y, yp, cflag );
	        if (cflag >= 5) {
		        printf("RKSUITE error %d\n", cflag );
		        return false;
	        }
        } while (tnow < stopTime);
        printf("%6.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
            tnow, y[0], y[1], y[2], y[3], y[4], y[5] );
        stopTime += 1.0;
        rksuite.reset( stopTime );
    };

    // next integrate until the object impacts the spherical earth.
    // to do this repeatedly call ct() checking the computeHalts()
    // output for a sign change.  when this occurs interpolate to the
    // exact impact time.
    stopTime = 10800.0;  // a gilligan
    rksuite.reset( stopTime );  

	double g[N_HALTS];  // current halt function values
    double gp[N_HALTS]; // prior halt function values
    double tp;    //prior cycle time

	computeHalts( tnow, y, gp );  // initialize halt function values
	tp = tnow;

    bool eventFound = false;
    // ct() requires the current derivatives as input initially
	computeDerivatives( tnow, y, yp );

	while(!eventFound) {

		rksuite.ct( computeDerivatives, tnow, y, yp, cflag );
		if (cflag >= 5) {  // ignore warnings
			printf("RKSUITE fatal error %d\n", cflag );
			return false;
		}
		if (tnow >= stopTime) {
			eventFound = true;  // ct() may have exited due to reaching the stop time
			break;
		}

		computeHalts( tnow, y, g );  // get the current halt vector values
		for (int k = 0; k < N_HALTS; k++) { // for each halt vector element
			double twant, ywant[6], ypwant[6];

			if (g[k]*gp[k] < 0.0) {  // if halt vector element changed sign
				double gt[N_HALTS];
                // ct() advanced beyond the sign change.  Reduce the interval
                // by a binary search until the sign change is located.
                // initial interval [tp .. tnow]
				bool   search = true;
				do {
					double dt = tnow - tp; // time in the interval
                    twant = 0.5*(tnow + tp); // mid point of interval
                    // interpolate the solution vector; we don't need the derivatives
                    // estimated to check for impact
					rksuite.intrp( twant, "S", 6, ywant, ypwant, computeDerivatives );
					computeHalts( twant, ywant, gt );  // halt vector at mid-point
					if (fabs(dt) < 1.0e-3 || fabs(gt[k]) < 1e-1) {
						search = false;  // interval too small or impact determined
                    } else if (gt[k]*g[k] >= 0) { 
                        // if intermediate point is same sign at current
                        // new interval [tp .. twant]
						tnow = twant;
                    } else {  
                        // else intermediate point is same sign as prior
                        // new interval [twant .. tnow]
						tp = twant;						
					}
				} while (search);
                // copy time and solution vector to output variables
				tnow = twant;
                for (int i = 0; i < 6; i++) {
				    y[i] = ywant[i];
                }
                computeDerivatives( tnow, y, yp ); //update the deriviates at impact
				eventFound = true;
				break;
			}
			gp[k] = g[k];
		}  // for k
		tp = tnow; // time at end of this cycle become the prior time for next cycle
	}  // while ! eventFound

    printf("%6.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf %10.1lf\n",
        tnow, y[0], y[1], y[2], y[3], y[4], y[5] );
    return 0;
}
