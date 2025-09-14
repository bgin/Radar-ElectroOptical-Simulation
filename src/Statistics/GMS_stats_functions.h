
#ifndef __GMS_STATS_FUNCTIONS_H__
#define __GMS_STATS_FUNCTIONS_H__

 // Helper statistical functions.

 //Globals
static const double eps = 1.0e-12;

static bool is_population_modified = false;

static const char writing_file[] = "Started writing results to file:";
static const char mean_accumul[] = "accumulated_mean:";
static const char population_value[] = "population_value:";
static const char index_i[] = "index i:";
static const char index_j[] = "index j:";
static const char ar_mean[] = "Arithmetic Mean:";
static const char geometric_mean[] = "Geometric Mean:";
static const char geometric_median[] = "Geometric Median:";
static const char standard_deviation[] = "Standard Deviation:";
static const char skeweness[] = "Skeweness:";
static const char writing_file_end[] = "Finished writing results to file:";
// Functions declarations.

// bool arith_mean(double **, const int, const int, double *)
/*
 * Computes arithmetic mean of sampled population.
 * Arguments: double ** sampled population, const int i, const int j, double * mean
 * Returns: double * mean by its argument and bool true on success otherwise false on failure.
 */
bool    arith_mean(double **, const int, const int, double *);

// bool geo_mean(double**, const int, const int, double *)
/*
* Computes geometric mean of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool   geo_mean(double **, const int, const int, double *);

// bool std_deviation(double **, const int, const int, double *)
/*
* Computes standard deviation of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool   std_deviation(double **, const int, const int, double *);

// bool geo_median(double **, const int, const int, double *)
/*
* Computes geometric median of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool   geomet_median(double **, const int, const int, double *);
// bool skew(double **, const int, const int, double *)
/*
* Computes skewennes of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool  skew(double **, const int, const int, double *);

//  File writing versions of functions

// bool arith_mean_wf(double **, const int, const int, double *, const char *)
/*
* Computes arithmetic mean of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean, const char * fname (file name)
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool  arith_mean_wf(double **, const int, const int, double *);

//  bool geo_mean_wf(double **, const int, const int, double *, const char *)
/*
* Computes geometric mean of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean, const char * fname (file name to be written)
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool  geo_mean_wf(double **, const int, const int, double *);

// bool std_deviation_wf(double **, const int, const int, double *, const char *)
/*
* Computes standard deviation of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean, const char * fname (file name to be writte)
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool  std_deviation_wf(double **, const int, const int, double *);
// bool geo_median_wf(double **, const int, const int, double *, const char *)
/*
* Computes geometric median of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean, const char *fname (file name to be written)
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool  geo_median_wf(double **, const int, const int, double *);

// bool skew_wf(double **, const int, const int, double *, const char *)
/*
* Computes skewennes of sampled population.
* Arguments: double ** sampled population, const int i, const int j, double * mean, const char * fname (file name to be written)
* Returns: double * mean by its argument and bool true on success otherwise false on failure.
*/
bool  skew_wf(double [][100], const int, const int, double *);

// Test functions - drivers.
void Test_skew_wf();
#endif /*__GMS_STATS_FUNCTIONS_H__*/
