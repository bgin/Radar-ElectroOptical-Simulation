void naca4_cambered ( double m, double p, double t, double c, int n, 
  double xc[], double xu[], double yu[], double xl[], double yl[] );
double *naca4_symmetric ( double t, double c, int n, double x[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
void r8mat_write ( string output_filename, int m, int n, double table[] );
double *r8vec_linspace_new ( int n, double a, double b );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );


