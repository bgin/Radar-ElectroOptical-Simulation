# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "naca.hpp"

//****************************************************************************80

void naca4_cambered ( double m, double p, double t, double c, int n, 
  double xc[], double xu[], double yu[], double xl[], double yl[] )

//****************************************************************************80
//
//  Purpose:
//
//    NACA4_CAMBERED: (xu,yu), (xl,yl) for a NACA cambered 4-digit airfoil.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
//    "The characteristics of 78 related airfoil sections from tests in
//    the variable-density wind tunnel",
//    NACA Report 460, 1933.
//
//  Parameters:
//
//    Input, double M, the maximum camber.
//    0.0 < M.
//
//    Input, double P, the location of maximum camber.
//    0.0 < P < 1.0
//
//    Input, double T, the maximum relative thickness.
//    0.0 < T <= 1.0
//
//    Input, double C, the chord length.
//    0.0 < C.
//
//    Input, int N, the number of sample points.
//
//    Input, double XC[N], points along the chord length.  
//    0.0 <= XC(*) <= C.
//
//    Output, double XU[N], YU[N], XL[N], YL[N], for each value of 
//    XC, measured along the camber line, the corresponding values (XU,YU) 
//    on the upper airfoil surface and (XL,YL) on the lower airfoil surface.
//
{
  double divisor;
  double dycdx;
  int i;
  double theta;
  double yc;
  double yt;

  for ( i = 0; i < n; i++ )
  {
    if ( 0.0 <= xc[i] / c && xc[i] / c <= p )
    {
      divisor = p * p;
    }
    else if ( p <= xc[i] / c && xc[i] / c <= 1.0 )
    {
      divisor = pow ( 1.0 - p, 2 );
    }
    else
    {
      divisor = 1.0;
    }

    dycdx = 2.0 * m * ( p - xc[i] / c ) / divisor;

    theta = atan ( dycdx );
   
    yt = 5.0 * t * c * ( 
       0.2969 * sqrt ( xc[i] / c ) 
       + (((( 
         - 0.1015 ) * ( xc[i] / c ) 
         + 0.2843 ) * ( xc[i] / c ) 
         - 0.3516 ) * ( xc[i] / c ) 
         - 0.1260 ) * ( xc[i] / c ) );

    if ( 0.0 <= xc[i] / c && xc[i] / c <= p )
    {
      yc = m * xc[i] * ( 2.0 * p - xc[i] / c ) / p / p;
    }
    else if ( p <= xc[i] / c && xc[i] / c <= 1.0 )
    {
      yc = m * ( xc[i] - c ) * ( 2.0 * p - xc[i] / c - 1.0 )
        / ( 1.0 - p ) / ( 1.0 - p );
    }
    else
    {
      yc = 0.0;
    }

    xu[i] = xc[i] - yt * sin ( theta );
    yu[i] = yc + yt * cos ( theta );
    xl[i] = xc[i] + yt * sin ( theta );
    yl[i] = yc - yt * cos ( theta );
  }
  return;
}
//****************************************************************************80

double *naca4_symmetric ( double t, double c, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    NACA4_SYMMETRIC evaluates y(x) for a NACA symmetric 4-digit airfoil.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eastman Jacobs, Kenneth Ward, Robert Pinkerton,
//    "The characteristics of 78 related airfoil sections from tests in
//    the variable-density wind tunnel",
//    NACA Report 460, 1933.
//
//  Parameters:
//
//    Input, double T, the maximum relative thickness.
//
//    Input, double C, the chord length.
//
//    Input, int N, the number of sample points.
//
//    Input, double X[N], points along the chord length.  
//    0.0 <= X(*) <= C.
//
//    Output, double NACA4_SYMMETRIC[N], for each value of X, the corresponding
//    value of Y so that (X,Y) is on the upper wing surface, and (X,-Y) is on the
//    lower wing surface.
//
{
  int i;
  double *y;

  y = new double[n];

  for ( i = 0; i < n; i++ )
  {
    y[i] = 5.0 * t * c * ( 
      0.2969 * sqrt ( x[i] / c ) 
      + (((( 
      - 0.1015 ) * ( x[i] / c ) 
      + 0.2843 ) * ( x[i] / c ) 
      - 0.3516 ) * ( x[i] / c ) 
      - 0.1260 ) * ( x[i] / c ) );
  }

  return y;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
