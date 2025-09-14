# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

#include "GMS_triangulation.h"

//****************************************************************************80

void alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *alpha_min, double *alpha_ave,
  double *alpha_area )

//****************************************************************************80
//
//  Purpose:
//
//    ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
//
//  Discusion:
//
//    The ALPHA measure evaluates the uniformity of the shapes of the triangles
//    defined by a triangulated pointset.
//
//    We compute the minimum angle among all the triangles in the triangulated
//    dataset and divide by the maximum possible value (which, in degrees,
//    is 60).  The best possible value is 1, and the worst 0.  A good
//    triangulation should have an ALPHA score close to 1.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, real ( kind = 8 ) Z(2,N), the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
//    the triangulation.
//
//    Output, double *ALPHA_MIN, the minimum value of ALPHA over all
//    triangles.
//
//    Output, double *ALPHA_AVE, the value of ALPHA averaged over
//    all triangles.
//
//    Output, double *ALPHA_AREA, the value of ALPHA averaged over
//    all triangles and weighted by area.
//
{
  double a_angle;
  int a_index;
  double a_x;
  double a_y;
  double ab_len;
  double alpha;
  double area;
  double area_total;
  double b_angle;
  int b_index;
  double b_x;
  double b_y;
  double bc_len;
  double c_angle;
  int c_index;
  double c_x;
  double c_y;
  double ca_len;
  double pi = 3.141592653589793;
  int triangle;
  double value;

  *alpha_min = r8_huge ( );
  *alpha_ave = 0.0;
  *alpha_area = 0.0;
  area_total = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    a_index = triangle_node[0+triangle*triangle_order];
    b_index = triangle_node[1+triangle*triangle_order];
    c_index = triangle_node[2+triangle*triangle_order];

    a_x = z[0+(a_index-1)*2];
    a_y = z[1+(a_index-1)*2];
    b_x = z[0+(b_index-1)*2];
    b_y = z[1+(b_index-1)*2];
    c_x = z[0+(c_index-1)*2];
    c_y = z[1+(c_index-1)*2];

    area = 0.5 * r8_abs ( a_x * ( b_y - c_y )
                        + b_x * ( c_y - a_y )
                        + c_x * ( a_y - b_y ) );

    ab_len = sqrt ( pow ( a_x - b_x, 2 ) + pow ( a_y - b_y, 2 ) );
    bc_len = sqrt ( pow ( b_x - c_x, 2 ) + pow ( b_y - c_y, 2 ) );
    ca_len = sqrt ( pow ( c_x - a_x, 2 ) + pow ( c_y - a_y, 2 ) );
//
//  Take care of a ridiculous special case.
//
    if ( ab_len == 0.0 && bc_len == 0.0 && ca_len == 0.0 )
    {
      a_angle = 2.0 * pi / 3.0;
      b_angle = 2.0 * pi / 3.0;
      c_angle = 2.0 * pi / 3.0;
    }
    else
    {
      if ( ca_len == 0.0 || ab_len == 0.0 )
      {
        a_angle = pi;
      }
      else
      {
        a_angle = arc_cosine (
          ( ca_len * ca_len + ab_len * ab_len - bc_len * bc_len )
          / ( 2.0 * ca_len * ab_len ) );
      }

      if ( ab_len == 0.0 || bc_len == 0.0 )
      {
        b_angle = pi;
      }
      else
      {
        b_angle = arc_cosine (
          ( ab_len * ab_len + bc_len * bc_len - ca_len * ca_len )
          / ( 2.0 * ab_len * bc_len ) );
      }

      if ( bc_len == 0.0 || ca_len == 0.0 )
      {
        c_angle = pi;
      }
      else
      {
        c_angle = arc_cosine (
          ( bc_len * bc_len + ca_len * ca_len - ab_len * ab_len )
          / ( 2.0 * bc_len * ca_len ) );
      }
    }
    *alpha_min = r8_min ( *alpha_min, a_angle );
    *alpha_min = r8_min ( *alpha_min, b_angle );
    *alpha_min = r8_min ( *alpha_min, c_angle );

    *alpha_ave = *alpha_ave + *alpha_min;

    *alpha_area = *alpha_area + area * *alpha_min;

    area_total = area_total + area;
  }
  *alpha_ave = *alpha_ave / ( double ) ( triangle_num );
  *alpha_area = *alpha_area / area_total;
//
//  Normalize angles from [0,pi/3] radians into qualities in [0,1].
//
  *alpha_min = *alpha_min * 3.0 / pi;
  *alpha_ave = *alpha_ave * 3.0 / pi;
  *alpha_area = *alpha_area * 3.0 / pi;

  return;
}
//****************************************************************************80

double angle_rad_2d ( double p1[2], double p2[2], double p3[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
//
//  Discussion:
//
//      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
//
//        P1
//        /
//       /
//      /
//     /
//    P2--------->P3
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double P1[2], P2[2], P3[2], define the rays
//    P1 - P2 and P3 - P2 which define the angle.
//
//    Output, double ANGLE_RAD_3D, the angle between the two rays,
//    in radians.  This value will always be between 0 and 2*PI.  If either ray has
//    zero length, then the angle is returned as zero.
//
{
# define DIM_NUM 2

  double p[DIM_NUM];
  double pi = 3.141592653589793;
  double value;

  p[0] = ( p3[0] - p2[0] ) * ( p1[0] - p2[0] )
       + ( p3[1] - p2[1] ) * ( p1[1] - p2[1] );


  p[1] = ( p3[0] - p2[0] ) * ( p1[1] - p2[1] )
       - ( p3[1] - p2[1] ) * ( p1[0] - p2[0] );

  if ( p[0] == 0.0 && p[1] == 0.0 )
  {
    value = 0.0;
    return value;
  }

  value = atan2 ( p[1], p[0] );

  if ( value < 0.0 )
  {
    value = value + 2.0 * pi;
  }

  return value;
# undef DIM_NUM
}
//****************************************************************************80

double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
# define PI 3.141592653589793

  double value;

  if ( c <= -1.0 )
  {
    value = PI;
  }
  else if ( 1.0 <= c )
  {
    value = 0.0;
  }
  else
  {
    value = acos ( c );
  }
  return value;
# undef PI
}
//****************************************************************************80

void area_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *area_min, double *area_max, double *area_ratio,
  double *area_ave, double *area_std )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_MEASURE determines the area ratio quality measure.
//
//  Discusion:
//
//    This measure computes the area of every triangle, and returns
//    the ratio of the minimum to the maximum triangle.  A value of
//    1 is "perfect", indicating that all triangles have the same area.
//    A value of 0 is the worst possible result.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double Z[2*N], the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the triangulation.
//
//    Output, double *AREA_MIN, *AREA_MAX, the minimum and maximum
//    areas.
//
//    Output, double *AREA_RATIO, the ratio of the minimum to the
//    maximum area.
//
//    Output, double *AREA_AVE, the average area.
//
//    Output, double *AREA_STD, the standard deviation of the areas.
//
{
  double area;
  int triangle;
  double value;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  *area_max = 0.0;
  *area_min = r8_huge ( );
  *area_ave = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    x1 = z[0+(triangle_node[0+triangle*triangle_order]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*triangle_order]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*triangle_order]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*triangle_order]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*triangle_order]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*triangle_order]-1)*2];

    area = 0.5 * r8_abs ( x1 * ( y2 - y3 )
                        + x2 * ( y3 - y1 )
                        + x3 * ( y1 - y2 ) );

    *area_min = r8_min ( *area_min, area );
    *area_max = r8_max ( *area_max, area );
    *area_ave = *area_ave + area;
  }

  *area_ave = *area_ave / ( double ) ( triangle_num );
  *area_std = 0.0;
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    x1 = z[0+(triangle_node[0+triangle*triangle_order]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*triangle_order]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*triangle_order]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*triangle_order]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*triangle_order]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*triangle_order]-1)*2];

    area = 0.5 * r8_abs ( x1 * ( y2 - y3 )
                        + x2 * ( y3 - y1 )
                        + x3 * ( y1 - y2 ) );

    *area_std = *area_std + pow ( area - *area_ave, 2 );
  }
  *area_std = sqrt ( *area_std / ( double ) ( triangle_num ) );

  if ( 0.0 < *area_max )
  {
    *area_ratio = *area_min / *area_max;
  }
  else
  {
    *area_ratio = 0.0;
  }

  return;
}
//****************************************************************************80

void bandwidth ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth associated with the finite element mesh.
//
//  Discussion:
//
//    The quantity computed here is the "geometric" bandwidth determined
//    by the finite element mesh alone.
//
//    If a single finite element variable is associated with each node
//    of the mesh, and if the nodes and variables are numbered in the
//    same way, then the geometric bandwidth is the same as the bandwidth
//    of a typical finite element matrix.
//
//    The bandwidth M is defined in terms of the lower and upper bandwidths:
//
//      M = ML + 1 + MU
//
//    where
//
//      ML = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but earlier column,
//
//      MU = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but later column.
//
//    Because the finite element node adjacency relationship is symmetric,
//    we are guaranteed that ML = MU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
//
//    Output, int *M, the bandwidth of the matrix.
//
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
//****************************************************************************80

bool delaunay_swap_test ( double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    DELAUNAY_SWAP_TEST performs the Delaunay swap test.
//
//  Discussion:
//
//    The current triangles are formed by nodes [0+2,3) and [0+3,4).
//    if a swap is recommended, the new triangles will be [0+2,4) and [1+3,4).
//
//      4     ?     4
//     / \         /|\
//    1---3  ==>  1 | 3
//     \ /         \|/
//      2           2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Graham Carey,
//    Computational Grids:
//    Generation, Adaptation and Solution Strategies,
//    Taylor and Francis, 1997,
//    ISBN13: 978-1560326359,
//    LC: QA377.C32.
//
//  Parameters:
//
//    Input, double XY[2*4], the coordinates of four points.
//
//    Output, bool SWAP, is TRUE if the triangles [0+2,4) and [1+3,4)
//    are to replace triangles [0+2,3) and [0+3,4).
//
{
  double a;
  double b;
  double c;
  double d;
  bool swap;
  double x13;
  double x14;
  double x23;
  double x24;
  double y13;
  double y14;
  double y23;
  double y24;

  x13 = xy[0+0*2] - xy[0+2*2];
  x14 = xy[0+0*2] - xy[0+3*2];
  x23 = xy[0+1*2] - xy[0+2*2];
  x24 = xy[0+1*2] - xy[0+3*2];

  y13 = xy[1+0*2] - xy[1+2*2];
  y14 = xy[1+0*2] - xy[1+3*2];
  y23 = xy[1+1*2] - xy[1+2*2];
  y24 = xy[1+1*2] - xy[1+3*2];

  a = x13 * x23 + y13 * y23;
  b = x24 * y14 - x14 * y24;
  c = x23 * y13 - x13 * y23;
  d = x24 * x14 + y14 * y24;
//
//  The reference gives two initial tests before the
//  main one.  However, there seems to be an error
//  in at least one of these tests.  Since they are
//  intended to avoid error in borderline cases, but
//  instead cause real error in common cases, they are
//  omitted for now.
//
// if ( 0.0 <= a && 0.0 <= d )
// {
//   swap = true;
// }
// else if ( a < d && d < 0.0 )
// {
//   swap = true;
//  }
//  else if ...

  if ( a * b < c * d )
  {
    swap = true;
  }
  else
  {
    swap = false;
  }

  return swap;
}
//****************************************************************************80

int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 )

//****************************************************************************80
//
//  Purpose:
//
//    DIAEDG chooses a diagonal edge.
//
//  Discussion:
//
//    The routine determines whether 0--2 or 1--3 is the diagonal edge
//    that should be chosen, based on the circumcircle criterion, where
//    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
//    quadrilateral in counterclockwise order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
//    vertices of a quadrilateral, given in counter clockwise order.
//
//    Output, int DIAEDG, chooses a diagonal:
//    +1, if diagonal edge 02 is chosen;
//    -1, if diagonal edge 13 is chosen;
//     0, if the four vertices are cocircular.
//
{
  double ca;
  double cb;
  double dx10;
  double dx12;
  double dx30;
  double dx32;
  double dy10;
  double dy12;
  double dy30;
  double dy32;
  double s;
  double tol;
  double tola;
  double tolb;
  int value;

  tol = 100.0 * r8_epsilon ( );

  dx10 = x1 - x0;
  dy10 = y1 - y0;
  dx12 = x1 - x2;
  dy12 = y1 - y2;
  dx30 = x3 - x0;
  dy30 = y3 - y0;
  dx32 = x3 - x2;
  dy32 = y3 - y2;

  tola = tol * r8_max ( fabs ( dx10 ),
               r8_max ( fabs ( dy10 ),
               r8_max ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

  tolb = tol * r8_max ( fabs ( dx12 ),
               r8_max ( fabs ( dy12 ),
               r8_max ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

  ca = dx10 * dx30 + dy10 * dy30;
  cb = dx12 * dx32 + dy12 * dy32;

  if ( tola < ca && tolb < cb )
  {
    value = -1;
  }
  else if ( ca < -tola && cb < -tolb )
  {
    value = 1;
  }
  else
  {
    tola = r8_max ( tola, tolb );
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb
      + ( dx32 * dy12 - dx12 * dy32 ) * ca;

    if ( tola < s )
    {
      value = -1;
    }
    else if ( s < -tola )
    {
      value = 1;
    }
    else
    {
      value = 0;
    }

  }

  return value;
}
//****************************************************************************80

unsigned long get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define UNSIGNED_LONG_MAX 4294967295UL

  time_t clock;
  int i;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,UNSIGNED_LONG_MAX].
//
    seed = ( unsigned long )
      ( ( ( double ) seed )
      * ( ( double ) UNSIGNED_LONG_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

# undef UNSIGNED_LONG_MAX
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Formula:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Discussion:
//
//    The sign of 0 and all positive integers is taken to be +1.
//    The sign of all negative integers is -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  if ( i < 0 )
  {
    return (-1);
  }
  else
  {
    return 1;
  }

}
//****************************************************************************80

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;

  return;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I  I4_WRAP
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

int i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
//
//  Check.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void i4col_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts an I4COL.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

int i4col_sorted_unique_count ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
//
//  Discussion:
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        unique_num = unique_num + 1;
        j1 = j2;
        break;
      }
    }
  }

  return unique_num;
}
//****************************************************************************80

void i4col_swap ( int m, int n, int a[], int icol1, int icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based!  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << "  ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void i4vec_heap_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
//
//  Discussion:
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//                /      \
//            A(1)         A(2)
//          /     \        /  \
//      A(3)       A(4)  A(5) A(6)
//      /  \       /   \
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhui, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, int A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ;; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }

      }

    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;

  }

  return;
}
//****************************************************************************80

int *i4vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR_NEW sets an I4VEC to the indicator vector.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR_NEW[N], the array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}
//****************************************************************************80

int i4vec_min ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the minimum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] < value )
    {
      value = a[i];
    }
  }
  return value;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(8) << a[i]  << "\n";
  }
  return;
}
//****************************************************************************80

void i4vec_reverse ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_REVERSE reverses the elements of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      N = 5,
//      A = ( 11, 12, 13, 14, 15 ).
//
//    Output:
//
//      A = ( 15, 14, 13, 12, 11 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A(N), the array to be reversed.
//
{
  int i;
  int j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}
//****************************************************************************80

void i4vec_sort_heap_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  i4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    i4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;

  }

  return;
}
//****************************************************************************80

int i4vec_sorted_unique ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORTED_UNIQUE finds unique elements in a sorted I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements in A.
//
//    Input/output, int A[N].  On input, the sorted
//    integer array.  On output, the unique elements in A.
//
//    Output, int I4VEC_SORTED_UNIQUE, the number of unique elements in A.
//
{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[unique_num-1] )
    {
      unique_num = unique_num + 1;
      a[unique_num-1] = a[i];
    }
  }

  return unique_num;
}
//****************************************************************************80

int i4vec2_compare ( int n, int a1[], int a2[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC2_COMPARE compares pairs of integers stored in two vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data items.
//
//    Input, int A1[N], A2[N], contain the two components of each item.
//
//    Input, int I, J, the items to be compared.  These values will be
//    1-based indices for the arrays A1 and A2.
//
//    Output, int I4VEC2_COMPARE, the results of the comparison:
//    -1, item I < item J,
//     0, item I = item J,
//    +1, item J < item I.
//
{
  int isgn;

  isgn = 0;

       if ( a1[i-1] < a1[j-1] )
  {
    isgn = -1;
  }
  else if ( a1[i-1] == a1[j-1] )
  {
         if ( a2[i-1] < a2[j-1] )
    {
      isgn = -1;
    }
    else if ( a2[i-1] < a2[j-1] )
    {
      isgn = 0;
    }
    else if ( a2[j-1] < a2[i-1] )
    {
      isgn = +1;
    }
  }
  else if ( a1[j-1] < a1[i-1] )
  {
    isgn = +1;
  }

  return isgn;
}
//****************************************************************************80

void i4vec2_sort_a ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
//
//  Discussion:
//
//    Each item to be sorted is a pair of integers (I,J), with the I
//    and J values stored in separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items of data.
//
//    Input/output, int A1[N], A2[N], the data to be sorted..
//
{
  int i;
  int indx;
  int isgn;
  int j;
  int temp;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

int i4vec2_sorted_unique ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC2_SORTED_UNIQUE keeps the unique elements in a array of pairs of integers.
//
//  Discussion:
//
//    Item I is stored as the pair A1(I), A2(I).
//
//    The items must have been sorted, or at least it must be the
//    case that equal items are stored in adjacent vector locations.
//
//    If the items were not sorted, then this routine will only
//    replace a string of equal values by a single representative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items.
//
//    Input/output, int A1[N], A2[N].
//    On input, the array of N items.
//    On output, an array of UNIQUE_NUM unique items.
//
//    Output, int I4VEC2_SORTED_UNIQUE, the number of unique items.
//
{
  int itest;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( itest = 1; itest < n; itest++ )
  {
    if ( a1[itest] != a1[unique_num-1] ||
         a2[itest] != a2[unique_num-1] )
    {
      a1[unique_num] = a1[itest];
      a2[unique_num] = a2[itest];
      unique_num = unique_num + 1;
    }
  }

  return unique_num;
}
//****************************************************************************80

int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv )

//****************************************************************************80
//
//  Purpose:
//
//    LRLINE determines where a point lies in relation to a directed line.
//
//  Discussion:
//
//    LRLINE determines whether a point is to the left of, right of,
//    or on a directed line parallel to a line through given points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
//    directed line is parallel to and at signed distance DV to the left of
//    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
//    which the position relative to the directed line is to be determined.
//
//    Input, double DV, the signed distance, positive for left.
//
//    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
//    to the right of, on, or left of the directed line.  LRLINE is 0 if
//    the line degenerates to a point.
//
{
  double dx;
  double dxu;
  double dy;
  double dyu;
  double t;
  double tol = 0.0000001;
  double tolabs;
  int value;

  dx = xv2 - xv1;
  dy = yv2 - yv1;
  dxu = xu - xv1;
  dyu = yu - yv1;

  tolabs = tol * r8_max ( fabs ( dx ),
                 r8_max ( fabs ( dy ),
                 r8_max ( fabs ( dxu ),
                 r8_max ( fabs ( dyu ), fabs ( dv ) ) ) ) );

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

  if ( tolabs < t )
  {
    value = 1;
  }
  else if ( -tolabs <= t )
  {
    value = 0;
  }
  else if ( t < -tolabs )
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

void lvec_print ( int n, bool a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    LVEC_PRINT prints a logical vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, bool A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8) << i
         << "  " << setw(1) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void mesh_base_one ( int node_num, int element_order, int element_num,
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ONE ensures that the element definition is 1-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//
//    If this function detects 0-based indexing, it converts to 1-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }
  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  This will be converted to 1-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] + 1;
      }
    }
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ONE - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

void mesh_base_zero ( int node_num, int element_order, int element_num,
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ZERO ensures that the element definition is zero-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//    Since a C++ program will naturally assume a 0-based indexing, it is
//    necessary to check a given element definition and, if it is actually
//    1-based, to convert it.
//
//    This function attempts to detect 1-based node indexing and correct it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  This will be converted to 0-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

void node_merge ( int dim_num, int node_num, double node_xy[],
  double tolerance, int node_rep[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_MERGE detects nodes that should be merged.
//
//  Discussion:
//
//    Two nodes "should" be merged if they are within TOLERANCE distance
//    of each other.
//
//    With a tolerance of 0, only exactly equal nodes are counted.
//
//    With a positive tolerance, a pair of nodes inside a circle of
//    radius TOLERANCE result in a count of 1 duplicate.
//
//    However, what do we do if nodes A, B and C are arranged in a line,!
//    with A and B just within TOLERANCE of each other, and B and C just
//    within tolerance of each other?  What we do here is make a choice
//    that can be defended consistently.  A and B define an equivalence
//    class because they are closer than TOLERANCE.  C is then added to
//    this equivalence class, because it is within TOLERANCE of at least
//    on thing in that equivalence class.
//
//    Thus, if 100 nodes are separated pairwise by slightly less
//    than TOLERANCE, a total of 99 duplicates will be counted.
//
//    The program starts out by giving each node its own label.
//    If it finds that two nodes should be merged, then the index of
//    one node is used as the label for both.  This process continues
//    until all nodes have been considered.  The number of unique nodes
//    is the number of unique values in the output quantity NODE_REP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[DIM_NUM*NODE_NUM], the nodes.
//
//    Input, double TOLERANCE, the maximum distance between
//    two nodes regarded as duplicate.
//
//    Output, int NODE_REP[NODE_NUM], the "representative" of each node.
//    NODE_REP(NODE) is the index of a node which is within TOLERANCE of node
//    NODE, or for which a chain of nodes can be found, all having the
//    same representative, and all of which are pairwise closer than TOLERANCE.
//
{
  double dist;
  int i;
  int j;
  int node1;
  int node2;
  int rep;
  double *rep_dist;

  rep_dist = new double[node_num];

  for ( node1 = 0; node1 < node_num; node1++ )
  {
    node_rep[node1] = node1;
  }

  for ( node1 = 0; node1 < node_num; node1++ )
  {
    for ( j = 0; j < node_num; j++ )
    {
      rep_dist[j] = r8_huge ( );
    }

    for ( node2 = 0; node2 < node_num; node2++ )
    {
      dist = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        dist = dist
          + pow ( node_xy[i+node1*dim_num] - node_xy[i+node2*dim_num], 2 );
      }
      dist = sqrt ( dist );

      rep = node_rep[node2];

      if ( dist < rep_dist[rep] )
      {
        rep_dist[rep] = dist;
      }
    }

    for ( node2 = 0; node2 < node_num; node2++ )
    {
      rep = node_rep[node2];
      if ( rep_dist[rep] <= tolerance )
      {
        node_rep[node2] = node1;
      }
    }

  }

  delete [] rep_dist;

  return;
}
//****************************************************************************80

int ns_adj_col_set ( int node_num, int triangle_num, int variable_num,
  int triangle_node[], int triangle_neighbor[], int node_u_variable[],
  int node_v_variable[], int node_p_variable[], int adj_col[] )

//****************************************************************************80
//
//  Purpose:
//
//    NS_ADJ_COL_SET sets the COL array in a Navier Stokes triangulation.
//
//  Discussion:
//
//    This routine also counts the the value and returns the value of
//    ADJ_NUM, the number of Navier-Stokes variable adjacencies, which
//    should be identical to the value that would have been computed
//    by calling NS_ADJ_COUNT.
//
//    This routine is called to set up the ADJ_COL array, which indicates
//    the number of entries needed to store each column in the sparse
//    compressed column storage to be used for the adjancency matrix.
//
//    The triangulation is assumed to involve 6-node triangles.
//
//    Variables for the horizontal and vertical velocities are associated
//    with every node.  Variables for the pressure are associated only with
//    the vertex nodes.
//
//    We are interested in determining the number of nonzero entries in the
//    stiffness matrix of the Stokes equations, or the jacobian matrix of
//    the Navier Stokes equations.  To this end, we will say, somewhat
//    too broadly, that two variables are "adjacent" if their associated
//    nodes both occur in some common element.  This adjacency of variables
//    I and J is taken to be equivalent to the possible nonzeroness of
//    matrix entries A(I,J) and A(J,I).
//
//    A sparse compressed column format is used to store the counts for
//    the nonzeroes.  In other words, while the value ADJ_NUM reports the
//    number of adjacencies, the vector ADJ_COL is sufficient to allow us
//    to properly set up a sparse compressed matrix for the actual storage
//    of the sparse matrix, if we desire to proceed.
//
//  Local Node Numbering:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//  Variable Diagram:
//
//      UVP
//       |\
//       | \
//       |  \
//      UV   UV
//       |    \
//       |     \
//       |      \
//      UVP--UV--UVP
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
//    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
//    vertical velocity and pressure variables associated with a node,
//    or -1 if no such variable is associated with the node.
//
//    Output, int ADJ_COL[VARIABLE_NUM+1].  Information about variable J
//    is stored in entries ADJ_COL[J] through ADJ_COL[J+1]-1 of ADJ.
//
//    Output, int NS_ADJ_COL_SET, the number of Navier Stokes variable
//    adjacencies.
//
{
  int adj_num;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int node;
  int p1;
  int p2;
  int p3;
  int triangle;
  int triangle_order = 6;
  int triangle2;
  int u1;
  int u2;
  int u3;
  int u4;
  int u5;
  int u6;
  int v1;
  int v2;
  int v3;
  int v4;
  int v5;
  int v6;
  int variable;

  adj_num = 0;
//
//  Set every variable to be adjacent to itself.
//
  for ( variable = 0; variable < variable_num; variable++ )
  {
    adj_col[variable] = 1;
  }
//
//  Set every variable to be adjacent to the other variables associated with
//  that node.
//
//  U <=> V
//  U <=> P (if there is a P variable)
//  V <=> P (if there is a P variable)
//
  for ( node = 0; node < node_num; node++ )
  {
    u1 = node_u_variable[node] - 1;
    v1 = node_v_variable[node] - 1;
    p1 = node_p_variable[node] - 1;

    adj_col[u1] = adj_col[u1] + 1;
    adj_col[v1] = adj_col[v1] + 1 ;

    if ( 0 <= p1 )
    {
      adj_col[u1] = adj_col[u1] + 1;
      adj_col[v1] = adj_col[v1] + 1;
      adj_col[p1] = adj_col[p1] + 2;
    }
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order] - 1;
    n2 = triangle_node[1+triangle*triangle_order] - 1;
    n3 = triangle_node[2+triangle*triangle_order] - 1;
    n4 = triangle_node[3+triangle*triangle_order] - 1;
    n5 = triangle_node[4+triangle*triangle_order] - 1;
    n6 = triangle_node[5+triangle*triangle_order] - 1;

    u1 = node_u_variable[n1] - 1;
    v1 = node_v_variable[n1] - 1;
    p1 = node_p_variable[n1] - 1;

    u2 = node_u_variable[n2] - 1;
    v2 = node_v_variable[n2] - 1;
    p2 = node_p_variable[n2] - 1;

    u3 = node_u_variable[n3] - 1;
    v3 = node_v_variable[n3] - 1;
    p3 = node_p_variable[n3] - 1;

    u4 = node_u_variable[n4] - 1;
    v4 = node_v_variable[n4] - 1;

    u5 = node_u_variable[n5] - 1;
    v5 = node_v_variable[n5] - 1;

    u6 = node_u_variable[n6] - 1;
    v6 = node_v_variable[n6] - 1;
//
//  For sure, we add the new adjacencies:
//
//    U5 V5 <=> U1 V1 P1
//    U6 V6 <=> U2 V2 P2
//    U4 V4 <=> U3 V3 P3
//    U5 V5 <=> U4 V4
//    U6 V6 <=> U4 V4
//    U6 V6 <=> U5 V5
//
    adj_col[u1] = adj_col[u1] + 2;
    adj_col[v1] = adj_col[v1] + 2;
    adj_col[p1] = adj_col[p1] + 2;

    adj_col[u2] = adj_col[u2] + 2;
    adj_col[v2] = adj_col[v2] + 2;
    adj_col[p2] = adj_col[p2] + 2;

    adj_col[u3] = adj_col[u3] + 2;
    adj_col[v3] = adj_col[v3] + 2;
    adj_col[p3] = adj_col[p3] + 2;

    adj_col[u4] = adj_col[u4] + 7;
    adj_col[v4] = adj_col[v4] + 7;

    adj_col[u5] = adj_col[u5] + 7;
    adj_col[v5] = adj_col[v5] + 7;

    adj_col[u6] = adj_col[u6] + 7;
    adj_col[v6] = adj_col[v6] + 7;
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//
//    U1 V1 P1 <=> U2 V2 P2
//    U1 V1 P1 <=> U4 V4
//    U2 V2 P2 <=> U4 V4
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[u1] = adj_col[u1] + 5;
      adj_col[v1] = adj_col[v1] + 5;
      adj_col[p1] = adj_col[p1] + 5;

      adj_col[u2] = adj_col[u2] + 5;
      adj_col[v2] = adj_col[v2] + 5;
      adj_col[p2] = adj_col[p2] + 5;

      adj_col[u4] = adj_col[u4] + 6;
      adj_col[v4] = adj_col[v4] + 6;
    }
//
//  Maybe add
//
//    U2 V2 P2 <=> U3 V3 P3
//    U2 V2 P2 <=> U5 V5
//    U3 V3 P3 <=> U5 V5
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[u2] = adj_col[u2] + 5;
      adj_col[v2] = adj_col[v2] + 5;
      adj_col[p2] = adj_col[p2] + 5;

      adj_col[u3] = adj_col[u3] + 5;
      adj_col[v3] = adj_col[v3] + 5;
      adj_col[p3] = adj_col[p3] + 5;

      adj_col[u5] = adj_col[u5] + 6;
      adj_col[v5] = adj_col[v5] + 6;
    }
//
//  Maybe add
//
//    U1 V1 P1 <=> U3 V3 P3
//    U1 V1 P1 <=> U6 V6
//    U3 V3 P3 <=> U6 V6
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[u1] = adj_col[u1] + 5;
      adj_col[v1] = adj_col[v1] + 5;
      adj_col[p1] = adj_col[p1] + 5;

      adj_col[u3] = adj_col[u3] + 5;
      adj_col[v3] = adj_col[v3] + 5;
      adj_col[p3] = adj_col[p3] + 5;

      adj_col[u6] = adj_col[u6] + 6;
      adj_col[v6] = adj_col[v6] + 6;
    }

  }
//
//  We used ADJ_COL to count the number of entries in each column.
//  Convert it to pointers into the ADJ array.
//
  for ( variable = variable_num; 0 < variable; variable-- )
  {
    adj_col[variable] = adj_col[variable-1];
  }

  adj_col[0] = 1;
  for ( variable = 1; variable <= variable_num; variable++ )
  {
    adj_col[variable] = adj_col[variable-1] + adj_col[variable];
  }

  adj_num = adj_col[variable_num] - 1;

  return adj_num;
}
//****************************************************************************80

int ns_adj_count ( int node_num, int triangle_num, int variable_num,
  int triangle_node[], int triangle_neighbor[], int node_u_variable[],
  int node_v_variable[], int node_p_variable[] )

//****************************************************************************80
//
//  Purpose:
//
//    NS_ADJ_COUNT counts adjacencies in a Navier Stokes triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The value of ADJ_NUM computed and returned by this routine should
//    be identical to the value computed by NS_ADJ_COL_SET.
//
//    The triangulation is assumed to involve 6-node triangles.
//
//    Variables for the horizontal and vertical velocities are associated
//    with every node.  Variables for the pressure are associated only with
//    the vertex nodes.
//
//    We are interested in determining the number of nonzero entries in the
//    stiffness matrix of the Stokes equations, or the jacobian matrix of
//    the Navier Stokes equations.  To this end, we will say, somewhat
//    too broadly, that two variables are "adjacent" if their associated
//    nodes both occur in some common element.  This adjacency of variables
//    I and J is taken to be equivalent to the possible nonzeroness of
//    matrix entries A(I,J) and A(J,I).
//
//    A sparse compressed column format is used to store the counts for
//    the nonzeroes.  In other words, while the value ADJ_NUM reports the
//    number of adjacencies, the vector ADJ_COL is sufficient to allow us
//    to properly set up a sparse compressed matrix for the actual storage
//    of the sparse matrix, if we desire to proceed.
//
//  Local Node Numbering:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//  Variable Diagram:
//
//      UVP
//       |\
//       | \
//       |  \
//      UV   UV
//       |    \
//       |     \
//       |      \
//      UVP--UV--UVP
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
//    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
//    vertical velocity and pressure variables associated with a node,
//    or -1 if no such variable is associated with the node.
//
//    Output, int NS_ADJ_COUNT, the value of ADJ_NUM, the number of
//    Navier Stokes variable adjacencies.
//
{
  int adj_num;
  int node;
  int p1;
  int triangle;
  int triangle_order = 6;
  int triangle2;
  int variable;

  adj_num = 0;
//
//  Set every variable to be adjacent to itself.
//
  adj_num = variable_num;
//
//  Set every variable to be adjacent to the other variables associated with
//  that node.
//
//  U <=> V
//  U <=> P (if there is a P variable)
//  V <=> P (if there is a P variable)
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_num = adj_num + 2;

    p1 = node_p_variable[node] - 1;

    if ( 0 <= p1 )
    {
      adj_num = adj_num + 4;
    }
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {

//
//  For sure, we add the new adjacencies:
//
//    U5 V5 <=> U1 V1 P1
//    U6 V6 <=> U2 V2 P2
//    U4 V4 <=> U3 V3 P3
//    U5 V5 <=> U4 V4
//    U6 V6 <=> U4 V4
//    U6 V6 <=> U5 V5
//
    adj_num = adj_num + 60;
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//
//    U1 V1 P1 <=> U2 V2 P2
//    U1 V1 P1 <=> U4 V4
//    U2 V2 P2 <=> U4 V4
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_num = adj_num + 42;
    }
//
//  Maybe add
//
//    U2 V2 P2 <=> U3 V3 P3
//    U2 V2 P2 <=> U5 V5
//    U3 V3 P3 <=> U5 V5
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_num = adj_num + 42;
    }
//
//  Maybe add
//
//    U1 V1 P1 <=> U3 V3 P3
//    U1 V1 P1 <=> U6 V6
//    U3 V3 P3 <=> U6 V6
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_num = adj_num + 42;
    }

  }

  return adj_num;
}
//****************************************************************************80

void ns_adj_insert ( int v1, int v2, int variable_num, int adj_num,
  int adj_col_free[], int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    NS_ADJ_INSERT inserts an adjacency into a compressed column adjacency matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int V1, V2, the indices of two items which are adjacent.
//
//    Input, int VARIABLE_NUM, the number of items.
//
//    Input, int ADJ_NUM, the number of entries available in ADJ_ROW.
//
//    Input/output, int ADJ_COL_FREE[VARIABLE_NUM], contains the next free
//    location in which an entry for a given column can be stored.  On output,
//    two pointers have been updated.
//
//    Input/output, int ADJ_ROW[ADJ_NUM], the row indices of the Navier Stokes
//    variable adjacency matrix.  On output, two new entries have been added.
//
{
  adj_row[adj_col_free[v1-1]-1] = v2;
  adj_col_free[v1-1] = adj_col_free[v1-1] + 1;

  if ( v1 == v2 )
  {
    return;
  }

  adj_row[adj_col_free[v2-1]-1] = v1;
  adj_col_free[v2-1] = adj_col_free[v2-1] + 1;

  return;
}
//****************************************************************************80

void ns_adj_row_set ( int node_num, int triangle_num, int variable_num,
  int triangle_node[], int triangle_neighbor[], int node_u_variable[],
  int node_v_variable[], int node_p_variable[], int adj_num, int adj_col[],
  int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    NS_ADJ_ROW_SET sets the Navier Stokes sparse compressed column row indices.
//
//  Discussion:
//
//    After NS_ADJ_COUNT has been called to count ADJ_NUM, the number of
//    variable adjacencies and to set up ADJ_COL, the compressed column pointer,
//    this routine can be called to assign values to ADJ_ROW, the row
//    indices of the sparse compressed column adjacency matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int VARIABLE_NUM, the number of variables.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, int NODE_U_VARIABLE[NODE_NUM], NODE_V_VARIABLE[NODE_NUM],
//    NODE_P_VARIABLE[NODE_NUM], the index of the horizontal velocity,
//    vertical velocity and pressure variables associated with a node,
//    or -1 if no such variable is associated with the node.
//
//    Input, int ADJ_NUM, the number of Navier Stokes variable adjacencies.
//
//    Input, int ADJ_COL[VARIABLE_NUM+1].  Information about variable J
//    is stored in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
//    Output, int ADJ_ROW[ADJ_NUM], the row indices of the Navier Stokes
//    variable adjacency matrix.
//
//  Local Parameters:
//
//    Local, int ADJ_COL_FREE[VARIABLE_NUM], for each column,
//    the location in ADJ_ROW which can store the next row index.
//
{
  int *adj_col_free;
  int k1;
  int k2;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int node;
  int p1;
  int p2;
  int p3;
  int triangle;
  int triangle_order = 6;
  int triangle2;
  int u1;
  int u2;
  int u3;
  int u4;
  int u5;
  int u6;
  int v;
  int v1;
  int v2;
  int v3;
  int v4;
  int v5;
  int v6;

  for ( v = 0; v < adj_num; v++ )
  {
    adj_row[v] = -1;
  }

  adj_col_free = new int[variable_num];

  for ( v = 0; v < variable_num; v++ )
  {
    adj_col_free[v] = adj_col[v];
  }
//
//  Set every variable to be adjacent to itself.
//  Here, we have to be careful to start at index 1.
//
  for ( v = 1; v <= variable_num; v++ )
  {
    ns_adj_insert ( v, v, variable_num, adj_num, adj_col_free, adj_row );
  }
//
//  Set every variable to be adjacent to the other variables associated with
//  that node.
//
//  U <=> V
//  U <=> P (if there is a P variable)
//  V <=> P (if there is a P variable)
//
  for ( node = 0; node < node_num; node++ )
  {
    u1 = node_u_variable[node];
    v1 = node_v_variable[node];
    p1 = node_p_variable[node];

    ns_adj_insert ( u1, v1, variable_num, adj_num, adj_col_free, adj_row );

    if ( 0 < p1 )
    {
      ns_adj_insert ( u1, p1, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, p1, variable_num, adj_num, adj_col_free, adj_row );
    }
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*6];
    n2 = triangle_node[1+triangle*6];
    n3 = triangle_node[2+triangle*6];
    n4 = triangle_node[3+triangle*6];
    n5 = triangle_node[4+triangle*6];
    n6 = triangle_node[5+triangle*6];

    u1 = node_u_variable[n1-1];
    v1 = node_v_variable[n1-1];
    p1 = node_p_variable[n1-1];

    u2 = node_u_variable[n2-1];
    v2 = node_v_variable[n2-1];
    p2 = node_p_variable[n2-1];

    u3 = node_u_variable[n3-1];
    v3 = node_v_variable[n3-1];
    p3 = node_p_variable[n3-1];

    u4 = node_u_variable[n4-1];
    v4 = node_v_variable[n4-1];

    u5 = node_u_variable[n5-1];
    v5 = node_v_variable[n5-1];

    u6 = node_u_variable[n6-1];
    v6 = node_v_variable[n6-1];
//
//  For sure, we add the new adjacencies:
//
//    U5 V5 <=> U1 V1 P1
//    U6 V6 <=> U2 V2 P2
//    U4 V4 <=> U3 V3 P3
//    U5 V5 <=> U4 V4
//    U6 V6 <=> U4 V4
//    U6 V6 <=> U5 V5
//
    ns_adj_insert ( u5, u1, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u5, v1, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u5, p1, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v5, u1, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v5, v1, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v5, p1, variable_num, adj_num, adj_col_free, adj_row );

    ns_adj_insert ( u6, u2, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u6, v2, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u6, p2, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, u2, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, v2, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, p2, variable_num, adj_num, adj_col_free, adj_row );

    ns_adj_insert ( u4, u3, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u4, v3, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u4, p3, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v4, u3, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v4, v3, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v4, p3, variable_num, adj_num, adj_col_free, adj_row );

    ns_adj_insert ( u5, u4, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u5, v4, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v5, u4, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v5, v4, variable_num, adj_num, adj_col_free, adj_row );

    ns_adj_insert ( u6, u4, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u6, v4, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, u4, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, v4, variable_num, adj_num, adj_col_free, adj_row );

    ns_adj_insert ( u6, u5, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( u6, v5, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, u5, variable_num, adj_num, adj_col_free, adj_row );
    ns_adj_insert ( v6, v5, variable_num, adj_num, adj_col_free, adj_row );
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//
//    U1 V1 P1 <=> U2 V2 P2
//    U1 V1 P1 <=> U4 V4
//    U2 V2 P2 <=> U4 V4
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      ns_adj_insert ( u1, u2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u1, v2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u1, p2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, u2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, v2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, p2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, u2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, v2, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, p2, variable_num, adj_num, adj_col_free, adj_row );

      ns_adj_insert ( u1, u4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u1, v4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, u4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, v4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, u4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, v4, variable_num, adj_num, adj_col_free, adj_row );

      ns_adj_insert ( u2, u4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u2, v4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, u4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, v4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, u4, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, v4, variable_num, adj_num, adj_col_free, adj_row );
    }
//
//  Maybe add
//
//    U2 V2 P2 <=> U3 V3 P3
//    U2 V2 P2 <=> U5 V5
//    U3 V3 P3 <=> U5 V5
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      ns_adj_insert ( u2, u3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u2, v3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u2, p3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, u3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, v3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, p3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, u3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, v3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, p3, variable_num, adj_num, adj_col_free, adj_row );

      ns_adj_insert ( u2, u5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u2, v5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, u5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v2, v5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, u5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p2, v5, variable_num, adj_num, adj_col_free, adj_row );

      ns_adj_insert ( u3, u5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u3, v5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v3, u5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v3, v5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p3, u5, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p3, v5, variable_num, adj_num, adj_col_free, adj_row );
    }
//
//  Maybe add
//
//    U1 V1 P1 <=> U3 V3 P3
//    U1 V1 P1 <=> U6 V6
//    U3 V3 P3 <=> U6 V6
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      ns_adj_insert ( u1, u3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u1, v3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u1, p3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, u3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, v3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, p3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, u3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, v3, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, p3, variable_num, adj_num, adj_col_free, adj_row );

      ns_adj_insert ( u1, u6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u1, v6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, u6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v1, v6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, u6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p1, v6, variable_num, adj_num, adj_col_free, adj_row );

      ns_adj_insert ( u3, u6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( u3, v6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v3, u6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( v3, v6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p3, u6, variable_num, adj_num, adj_col_free, adj_row );
      ns_adj_insert ( p3, v6, variable_num, adj_num, adj_col_free, adj_row );
    }
  }
//
//  Ascending sort the entries for each variable.
//
  for ( v = 0; v < variable_num; v++ )
  {
    k1 = adj_col[v];
    k2 = adj_col[v+1]-1;

    i4vec_sort_heap_a ( k2+1-k1, adj_row+k1-1 );
  }

  return;
}
//****************************************************************************80

bool perm_check2 ( int n, int p[], int base )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK2 checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK2, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = base; seek < base + n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      return false;
    }

  }

  return true;
}
//****************************************************************************80

void perm_inverse ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE inverts a permutation "in place".
//
//  Discussion:
//
//    This algorithm assumes that the entries in the permutation vector are
//    strictly positive.  In particular, the value 0 must not occur.
//
//    When necessary, this function shifts the data temporarily so that
//    this requirement is satisfied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
  int base;
  int i;
  int i0;
  int i1;
  int i2;
  int is;
  int p_min;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "PERM_INVERSE - Fatal error!\n";
    cout << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }
//
//  Find the least value, and shift data so it begins at 1.
//
  p_min = i4vec_min ( n, p );
  base = 1;

  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - p_min + base;
  }
//
//  Now we can safely check the permutation.
//
  if ( !perm_check2 ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "PERM_INVERSE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  Now we can invert the permutation.
//
  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = - i4_sign ( p[i-1] );
    p[i-1] = i4_sign ( is ) * abs ( p[i-1] );
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = - p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }
        i0 = i1;
        i1 = i2;
      }
    }
  }
//
//  Now we can restore the permutation.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + p_min - base;
  }

  return;
}
//****************************************************************************80

int *points_delaunay_naive_2d ( int node_num, double node_xy[],
  int *triangle_num )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_DELAUNAY_NAIVE_2D computes the Delaunay triangulation in 2D.
//
//  Discussion:
//
//    A naive and inefficient (but extremely simple) method is used.
//
//    This routine is only suitable as a demonstration code for small
//    problems.  Its running time is of order NODE_NUM^4.  Much faster
//    algorithms are available.
//
//    Given a set of nodes in the plane, a triangulation is a set of
//    triples of distinct nodes, forming triangles, so that every
//    point with the convex hull of the set of  nodes is either one
//    of the nodes, or lies on an edge of one or more triangles,
//    or lies within exactly one triangle.
//
//    The number of nodes must be at least 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Joseph ORourke,
//    Computational Geometry,
//    Cambridge University Press,
//    Second Edition, 1998, page 187.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, int *TRIANGLE_NUM, the number of triangles.
//
//    Output, int POINTS_DELAUNAY_NAIVE_2D[3*TRIANGLE_NUM], the indices of the
//    nodes making each triangle.
//
{
  int count;
  int flag;
  int i;
  int j;
  int k;
  int m;
  int pass;
  int *tri;
  double xn;
  double yn;
  double zn;
  double *z;

  count = 0;

  z = new double [ node_num ];

  for ( i = 0; i < node_num; i++ )
  {
    z[i] = node_xy[0+i*2] * node_xy[0+i*2] + node_xy[1+i*2] * node_xy[1+i*2];
  }
//
//  First pass counts triangles,
//  Second pass allocates triangles and sets them.
//
  for ( pass = 1; pass <= 2; pass++ )
  {
    if ( pass == 2 )
    {
      tri = new int[3*count];
    }
    count = 0;
//
//  For each triple (I,J,K):
//
    for ( i = 0; i < node_num - 2; i++ )
    {
      for ( j = i+1; j < node_num; j++ )
      {
        for ( k = i+1; k < node_num; k++ )
        {
          if ( j != k )
          {
            xn = ( node_xy[1+j*2] - node_xy[1+i*2] ) * ( z[k] - z[i] )
               - ( node_xy[1+k*2] - node_xy[1+i*2] ) * ( z[j] - z[i] );
            yn = ( node_xy[0+k*2] - node_xy[0+i*2] ) * ( z[j] - z[i] )
               - ( node_xy[0+j*2] - node_xy[0+i*2] ) * ( z[k] - z[i] );
            zn = ( node_xy[0+j*2] - node_xy[0+i*2] )
               * ( node_xy[1+k*2] - node_xy[1+i*2] )
               - ( node_xy[0+k*2] - node_xy[0+i*2] )
               * ( node_xy[1+j*2] - node_xy[1+i*2] );

            flag = ( zn < 0 );

            if ( flag )
            {
              for ( m = 0; m < node_num; m++ )
              {
                flag = flag && ( ( node_xy[0+m*2] - node_xy[0+i*2] ) * xn
                               + ( node_xy[1+m*2] - node_xy[1+i*2] ) * yn
                               + ( z[m] - z[i] ) * zn <= 0 );
              }
            }

            if ( flag )
            {
              if ( pass == 2 )
              {
                tri[0+count*3] = i + 1;
                tri[1+count*3] = j + 1;
                tri[2+count*3] = k + 1;
              }
              count = count + 1;
            }
          }
        }
      }
    }
  }

  *triangle_num = count;
  delete [] z;

  return tri;
}
//****************************************************************************80

void points_hull_2d ( int node_num, double node_xy[], int *hull_num,
  int hull[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_HULL_2D computes the convex hull of a set of nodes in 2D.
//
//  Discussion:
//
//    The work involved is N*log(H), where N is the number of points, and H is
//    the number of points that are on the hull.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, int *HULL_NUM, the number of nodes that lie on the convex hull.
//
//    Output, int HULL[NODE_NUM].  The first HULL_NUM entries contain
//    the indices of the nodes that form the convex hull, in order.
//    These indices are 1-based, not 0-based!
//
{
  double angle;
  double angle_max;
  double di;
  double dr;
  int first;
  int i;
  double p_xy[2];
  int q;
  double q_xy[2];
  int r;
  double r_xy[2];

  *hull_num = 0;

  if ( node_num < 1 )
  {
    return;
  }
//
//  If NODE_NUM = 1, the hull is the node.
//
  if ( node_num == 1 )
  {
    hull[*hull_num] = 1;
    *hull_num = *hull_num + 1;
    return;
  }
//
//  If NODE_NUM = 2, then the convex hull is either the two distinct nodes,
//  or possibly a single (repeated) node.
//
  if ( node_num == 2 )
  {
    hull[*hull_num] = 1;
    *hull_num = *hull_num + 1;

    if ( node_xy[0+0*2] != node_xy[0+1*2] || node_xy[1+0*2] != node_xy[1+1*2] )
    {
      hull[*hull_num] = 2;
      *hull_num = *hull_num + 1;
    }

    return;
  }
//
//  Find the leftmost point, and take the bottom-most in a tie.
//  Call it "Q".
//
  q = 1;
  for ( i = 2; i <= node_num; i++ )
  {
    if ( node_xy[0+(i-1)*2] < node_xy[0+(q-1)*2] ||
      ( node_xy[0+(i-1)*2] == node_xy[0+(q-1)*2] &&
        node_xy[1+(i-1)*2] < node_xy[1+(q-1)*2] ) )
    {
      q = i;
    }
  }

  q_xy[0] = node_xy[0+(q-1)*2];
  q_xy[1] = node_xy[1+(q-1)*2];
//
//  Remember the starting point.
//
  first = q;
  hull[*hull_num] = q;
  *hull_num = *hull_num + 1;
//
//  For the first point, make a dummy previous point, 1 unit south,
//  and call it "P".
//
  p_xy[0] = q_xy[0];
  p_xy[1] = q_xy[1] - 1.0;
//
//  Now, having old point P, and current point Q, find the new point R
//  so the angle PQR is maximal.
//
//  Watch out for the possibility that the two nodes are identical.
//
  for ( ; ; )
  {
    r = 0;
    angle_max = 0.0;

    for ( i = 1; i <= node_num; i++ )
    {
      if ( i != q && ( node_xy[0+(i-1)*2] != q_xy[0] || node_xy[1+(i-1)*2] != q_xy[1] ) )
      {
        angle = angle_rad_2d ( p_xy, q_xy, node_xy+(i-1)*2 );

        if ( r == 0 || angle_max < angle )
        {
          r = i;
          r_xy[0] = node_xy[0+(r-1)*2];
          r_xy[1] = node_xy[1+(r-1)*2];
          angle_max = angle;
        }
//
//  In case of ties, choose the nearer point.
//
        else if ( r != 0 && angle == angle_max )
        {
          di = sqrt ( pow ( node_xy[0+(i-1)*2] - q_xy[0], 2 )
                    + pow ( node_xy[1+(i-1)*2] - q_xy[1], 2 ) );

          dr = sqrt ( pow ( r_xy[0] - q_xy[0], 2 )
                    + pow ( r_xy[1] - q_xy[1], 2 ) );

          if ( di < dr )
          {
            r = i;
            r_xy[0] = node_xy[0+(r-1)*2];
            r_xy[1] = node_xy[1+(r-1)*2];
            angle_max = angle;
          }
        }
      }
    }
//
//  If we've returned to our starting node, exit.
//
    if ( r == first )
    {
      break;
    }

    if ( node_num < *hull_num + 1 )
    {
      cout << "\n";
      cout << "POINTS_HULL_2D - Fatal error!\n";
      cout << "  The algorithm failed.\n";
      exit ( 1 );
    }
//
//  Add point R to the convex hull.
//
    hull[*hull_num] = r;
    *hull_num = *hull_num + 1;
//
//  Set Q := P, P := R, and repeat.
//
    q = r;

    p_xy[0] = q_xy[0];
    p_xy[1] = q_xy[1];

    q_xy[0] = r_xy[0];
    q_xy[1] = r_xy[1];
  }

  return;
}
//****************************************************************************80

int points_point_near_naive_nd ( int dim_num, int nset, double pset[],
  double ptest[], double *d_min )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
//
//  Discussion:
//
//    A naive algorithm is used.  The distance to every point is calculated,
//    in order to determine the smallest.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[DIM_NUM*NSET], the coordinates of the points
//    in the set.
//
//    Input, double PTEST[DIM_NUM], the point whose nearest neighbor is sought.
//
//    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
//
//    Output, int POINTS_POINT_NEAR_NAIVE_ND, I_MIN, the index of the nearest
//    point in PSET to P.
//
{
  double d;
  int i;
  int j;
  int p_min;

  *d_min = r8_huge ( );
  p_min = -1;

  for ( j = 0; j < nset; j++ )
  {
    d = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      d = d + ( ptest[i] - pset[i+j*dim_num] )
            * ( ptest[i] - pset[i+j*dim_num] );
    }
    if ( d < *d_min )
    {
      *d_min = d;
      p_min = j + 1;
    }
  }

  *d_min = sqrt ( *d_min );

  return p_min;
}
//****************************************************************************80

void q_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *q_min, double *q_max, double *q_ave,
  double *q_area )

//****************************************************************************80
//
//  Purpose:
//
//    Q_MEASURE determines the triangulated pointset quality measure Q.
//
//  Discussion:
//
//    The Q measure evaluates the uniformity of the shapes of the triangles
//    defined by a triangulated pointset.
//
//    For a single triangle T, the value of Q(T) is defined as follows:
//
//      TAU_IN = radius of the inscribed circle,
//      TAU_OUT = radius of the circumscribed circle,
//
//      Q(T) = 2 * TAU_IN / TAU_OUT
//        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
//
//    where A, B and C are the lengths of the sides of the triangle T.
//
//    The Q measure computes the value of Q(T) for every triangle T in the
//    triangulation, and then computes the minimum of this
//    set of values:
//
//      Q_MEASURE = min ( all T in triangulation ) Q(T)
//
//    In an ideally regular mesh, all triangles would have the same
//    equilateral shape, for which Q = 1.  A good mesh would have
//    0.5 < Q.
//
//    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
//    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
//    triangles.  Generally, a maximal triangulation is expected, namely,
//    a triangulation whose image is a planar graph, but for which the
//    addition of any new triangle would mean the graph was no longer planar.
//    A Delaunay triangulation is a maximal triangulation which maximizes
//    the minimum angle that occurs in any triangle.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//    Per-Olof Persson and Gilbert Strang,
//    A Simple Mesh Generator in MATLAB,
//    SIAM Review,
//    Volume 46, Number 2, pages 329-345, June 2004.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double Z[2*N], the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the triangulation.
//
//    Output, double *Q_MIN, *Q_MAX, the minimum and maximum values
//    of Q over all triangles.
//
//    Output, double *Q_AVE, the average value of Q.
//
//    Output, double *Q_AREA, the average value of Q, weighted by
//    the area of each triangle.
//
{
  int a_index;
  double ab_length;
  double area;
  double area_total;
  int b_index;
  double bc_length;
  int c_index;
  double ca_length;
  double q;
  int triangle;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  *q_min = r8_huge ( );
  *q_max = - r8_huge ( );
  *q_ave = 0.0;
  *q_area = 0.0;
  area_total = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    a_index = triangle_node[0+triangle*triangle_order];
    b_index = triangle_node[1+triangle*triangle_order];
    c_index = triangle_node[2+triangle*triangle_order];

    ab_length = sqrt (
        pow ( z[0+(a_index-1)*2] - z[0+(b_index-1)*2], 2 )
      + pow ( z[1+(a_index-1)*2] - z[1+(b_index-1)*2], 2 ) );

    bc_length = sqrt (
        pow ( z[0+(b_index-1)*2] - z[0+(c_index-1)*2], 2 )
      + pow ( z[1+(b_index-1)*2] - z[1+(c_index-1)*2], 2 ) );

    ca_length = sqrt (
        pow ( z[0+(c_index-1)*2] - z[0+(a_index-1)*2], 2 )
      + pow ( z[1+(c_index-1)*2] - z[1+(a_index-1)*2], 2 ) );

    q = ( bc_length + ca_length - ab_length )
      * ( ca_length + ab_length - bc_length )
      * ( ab_length + bc_length - ca_length )
      / ( ab_length * bc_length * ca_length );

    x1 = z[0+(triangle_node[0+triangle*triangle_order]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*triangle_order]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*triangle_order]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*triangle_order]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*triangle_order]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*triangle_order]-1)*2];

    area = 0.5 * r8_abs ( x1 * ( y2 - y3 )
                        + x2 * ( y3 - y1 )
                        + x3 * ( y1 - y2 ) );

    *q_min = r8_min ( *q_min, q );
    *q_max = r8_max ( *q_max, q );
    *q_ave = *q_ave + q;
    *q_area = *q_area + q * area;

    area_total = area_total + area;
  }

  *q_ave = *q_ave / ( double ) ( triangle_num );

  if ( 0.0 < area_total )
  {
    *q_area = *q_area / area_total;
  }
  else
  {
    *q_area = 0.0;
  }

  return;
}
//****************************************************************************80

void quad_convex_random ( int *seed, double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
//
//  Description:
//
//    The quadrilateral is constrained in that the vertices must all lie
//    with the unit square.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number
//    generator.
//
//    Output, double XY[2*NODE_NUM], the coordinates of the
//    nodes of the quadrilateral, given in counterclockwise order.
//
{
  int hull[4];
  int hull_num;
  int i;
  int j;
  double xy_random[2*4];

  for ( ; ; )
  {
//
//  Generate 4 random points.
//
    r8mat_uniform_01 ( 2, 4, seed, xy_random );
//
//  Determine the convex hull.
//
    points_hull_2d ( 4, xy_random, &hull_num, hull );
//
//  If HULL_NUM < 4, then our convex hull is a triangle.
//  Try again.
//
    if ( hull_num == 4 )
    {
      break;
    }
  }
//
//  Make an ordered copy of the random points.
//
  for ( j = 0; j < 4; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      xy[i+j*2] = xy_random[i+(hull[j]-1)*2];
    }
  }
  return;
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
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

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r8_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r8_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate,
//    strictly between 0 and 1.
//
{
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r82vec_permute ( int n, int p[], int base, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PERMUTE permutes an R82VEC in place.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type double precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
//
//    Input/output, double A[2*N], the array to be permuted.
//
{
  double a_temp[2];
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check2 ( n, p, base ) )
  {
    cerr << "\n";
    cerr << "R82VEC_PERMUTE - Fatal error!\n";
    cerr << "  PERM_CHECK rejects this permutation.\n";
    exit ( 1 );
  }
//
//  In order for the sign negation trick to work, we need to assume that the
//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
//  So temporarily add 1-BASE to each entry to force positivity.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] + 1 - base;
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = - p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = - p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cout << "\n";
          cout << "R82VEC_PERMUTE - Fatal error!\n";
          cout << "  Entry IPUT = " << iput << " of the permutation has\n";
          cout << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = - p[i];
  }
//
//  Restore the base of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = p[i] - 1 + base;
  }
  return;
}
//****************************************************************************80

int *r82vec_sort_heap_index_a ( int n, int base, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type double precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(*,indx(*))
//
//    or explicitly, by the call
//
//      r82vec_permute ( n, indx, base, a )
//
//    after which a(*,*) is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int BASE, the desired indexing for the sort index:
//    0 for 0-based indexing,
//    1 for 1-based indexing.
//
//    Input, double A[2*N], an array to be index-sorted.
//
//    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
//    I-th element of the sorted array is A(0:1,R8VEC_SORT_HEAP_INDEX_A(I)).
//
{
  double aval[2];
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  if ( n < 1 )
  {
    return NULL;
  }

  indx = new int[n];

  for ( i = 0; i < n; i++ )
  {
    indx[i] = i;
  }

  if ( n == 1 )
  {
    indx[0] = indx[0] + base;
    return indx;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
    }
    else
    {
      indxt = indx[ir-1];
      aval[0] = a[0+indxt*2];
      aval[1] = a[1+indxt*2];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }
    }
    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if (   a[0+indx[j-1]*2] <  a[0+indx[j]*2] ||
             ( a[0+indx[j-1]*2] == a[0+indx[j]*2] &&
               a[1+indx[j-1]*2] <  a[1+indx[j]*2] ) )
        {
          j = j + 1;
        }
      }

      if (   aval[0] <  a[0+indx[j-1]*2] ||
           ( aval[0] == a[0+indx[j-1]*2] &&
             aval[1] <  a[1+indx[j-1]*2] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
//
//  Take care of the base.
//
  for ( i = 0; i < n; i++ )
  {
    indx[i] = indx[i] + base;
  }

  return indx;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
}
//****************************************************************************80

int r8tris2 ( int node_num, double node_xy[], int *triangle_num,
  int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8TRIS2 constructs a Delaunay triangulation of 2D vertices.
//
//  Discussion:
//
//    The routine constructs the Delaunay triangulation of a set of 2D vertices
//    using an incremental approach and diagonal edge swaps.  Vertices are
//    first sorted in lexicographically increasing (X,Y) order, and
//    then are inserted one at a time from outside the convex hull.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input/output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//    On output, the coordinates have been sorted into dictionary order.
//
//    Output, int *TRIANGLE_NUM, the number of triangles in the triangulation;
//    TRIANGLE_NUM is equal to 2*node_num - NB - 2, where NB is the number
//    of boundary vertices.
//
//    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up each
//    triangle.  The elements are indices of NODE_XY.  The vertices of the
//    triangles are in counterclockwise order.
//
//    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
//    Positive elements are indices of TIL; negative elements are used for links
//    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
//    where I, J = triangle, edge index; TRIANGLE_NEIGHBOR[I,J] refers to
//    the neighbor along edge from vertex J to J+1 (mod 3).
//
//    Output, int R8TRIS2, is 0 for no error.
{
  int base;
  double cmax;
  int e;
  int error;
  int i;
  int *indx;
  int j;
  int k;
  int l;
  int ledg;
  int lr;
  int ltri;
  int m;
  int m1;
  int m2;
  int n;
  int redg;
  int rtri;
  int *stack;
  int t;
  double tol;
  int top;

  stack = new int[node_num];

  tol = 100.0 * r8_epsilon ( );
//
//  Sort the vertices by increasing (x,y).
//
  base = 0;

  indx = r82vec_sort_heap_index_a ( node_num, base, node_xy );

  r82vec_permute ( node_num, indx, base, node_xy );
//
//  Make sure that the nodes are "reasonably" distinct.
//
  m1 = 1;

  for ( i = 2; i <= node_num; i++ )
  {
    m = m1;
    m1 = i;

    k = -1;

    for ( j = 0; j <= 1; j++ )
    {
      cmax = r8_max ( fabs ( node_xy[2*(m-1)+j] ),
                     fabs ( node_xy[2*(m1-1)+j] ) );

      if ( tol * ( cmax + 1.0 )
           < fabs ( node_xy[2*(m-1)+j] - node_xy[2*(m1-1)+j] ) )
      {
        k = j;
        break;
      }

    }

    if ( k == -1 )
    {
      cout << "\n";
      cout << "R8TRIS2 - Fatal error!\n";
      cout << "  Fails for point number I = " << i << "\n";
      cout << "  M =  " << m  << "\n";
      cout << "  M1 = " << m1 << "\n";
      cout << "  X,Y(M)  = " << node_xy[2*(m-1)+0] << "  "
                             << node_xy[2*(m-1)+1] << "\n";
      cout << "  X,Y(M1) = " << node_xy[2*(m1-1)+0] << "  "
                             << node_xy[2*(m1-1)+1] << "\n";
      exit ( 1 );
    }

  }
//
//  Starting from nodes M1 and M2, search for a third point M that
//  makes a "healthy" triangle (M1,M2,M)
//
  m1 = 1;
  m2 = 2;
  j = 3;

  for ( ; ; )
  {
    if ( node_num < j )
    {
      cout << "\n";
      cout << "R8TRIS2 - Fatal error!\n";
      delete [] stack;
      return 225;
    }

    m = j;

    lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
      node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
      node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

    if ( lr != 0 )
    {
      break;
    }

    j = j + 1;

  }
//
//  Set up the triangle information for (M1,M2,M), and for any other
//  triangles you created because nodes were collinear with M1, M2.
//
  *triangle_num = j - 2;

  if ( lr == -1 )
  {
    triangle_node[3*0+0] = m1;
    triangle_node[3*0+1] = m2;
    triangle_node[3*0+2] = m;
    triangle_neighbor[3*0+2] = -3;

    for ( i = 2; i <= *triangle_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      triangle_node[3*(i-1)+0] = m1;
      triangle_node[3*(i-1)+1] = m2;
      triangle_node[3*(i-1)+2] = m;
      triangle_neighbor[3*(i-1)+0] = -3 * i;
      triangle_neighbor[3*(i-1)+1] = i;
      triangle_neighbor[3*(i-1)+2] = i - 1;

    }

    triangle_neighbor[3*(*triangle_num-1)+0] = -3 * (*triangle_num) - 1;
    triangle_neighbor[3*(*triangle_num-1)+1] = -5;
    ledg = 2;
    ltri = *triangle_num;
  }
  else
  {
    triangle_node[3*0+0] = m2;
    triangle_node[3*0+1] = m1;
    triangle_node[3*0+2] = m;
    triangle_neighbor[3*0+0] = -4;

    for ( i = 2; i <= *triangle_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      triangle_node[3*(i-1)+0] = m2;
      triangle_node[3*(i-1)+1] = m1;
      triangle_node[3*(i-1)+2] = m;
      triangle_neighbor[3*(i-2)+2] = i;
      triangle_neighbor[3*(i-1)+0] = -3 * i - 3;
      triangle_neighbor[3*(i-1)+1] = i - 1;
    }

    triangle_neighbor[3*(*triangle_num-1)+2] = -3 * (*triangle_num);
    triangle_neighbor[3*0+1] = -3 * (*triangle_num) - 2;
    ledg = 2;
    ltri = 1;

  }
//
//  Insert the vertices one at a time from outside the convex hull,
//  determine visible boundary edges, and apply diagonal edge swaps until
//  Delaunay triangulation of vertices (so far) is obtained.
//
  top = 0;

  for ( i = j+1; i <= node_num; i++ )
  {
    m = i;
    m1 = triangle_node[3*(ltri-1)+ledg-1];

    if ( ledg <= 2 )
    {
      m2 = triangle_node[3*(ltri-1)+ledg];
    }
    else
    {
      m2 = triangle_node[3*(ltri-1)+0];
    }

    lr = lrline ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1],
      node_xy[2*(m1-1)+0], node_xy[2*(m1-1)+1],
      node_xy[2*(m2-1)+0], node_xy[2*(m2-1)+1], 0.0 );

    if ( 0 < lr )
    {
      rtri = ltri;
      redg = ledg;
      ltri = 0;
    }
    else
    {
      l = -triangle_neighbor[3*(ltri-1)+ledg-1];
      rtri = l / 3;
      redg = (l % 3) + 1;
    }

    vbedg ( node_xy[2*(m-1)+0], node_xy[2*(m-1)+1], node_num,
      node_xy, *triangle_num, triangle_node, triangle_neighbor,
      &ltri, &ledg, &rtri, &redg );

    n = *triangle_num + 1;
    l = -triangle_neighbor[3*(ltri-1)+ledg-1];

    for ( ; ; )
    {
      t = l / 3;
      e = ( l % 3 ) + 1;
      l = -triangle_neighbor[3*(t-1)+e-1];
      m2 = triangle_node[3*(t-1)+e-1];

      if ( e <= 2 )
      {
        m1 = triangle_node[3*(t-1)+e];
      }
      else
      {
        m1 = triangle_node[3*(t-1)+0];
      }

      *triangle_num = *triangle_num + 1;
      triangle_neighbor[3*(t-1)+e-1] = *triangle_num;
      triangle_node[3*(*triangle_num-1)+0] = m1;
      triangle_node[3*(*triangle_num-1)+1] = m2;
      triangle_node[3*(*triangle_num-1)+2] = m;
      triangle_neighbor[3*(*triangle_num-1)+0] = t;
      triangle_neighbor[3*(*triangle_num-1)+1] = *triangle_num - 1;
      triangle_neighbor[3*(*triangle_num-1)+2] = *triangle_num + 1;
      top = top + 1;

      if ( node_num < top )
      {
        cout << "\n";
        cout << "R8TRIS2 - Fatal error!\n";
        cout << "  Stack overflow.\n";
        delete [] stack;
        return 8;
      }

      stack[top-1] = *triangle_num;

      if ( t == rtri && e == redg )
      {
        break;
      }

    }

    triangle_neighbor[3*(ltri-1)+ledg-1] = -3 * n - 1;
    triangle_neighbor[3*(n-1)+1] = -3 * (*triangle_num) - 2;
    triangle_neighbor[3*(*triangle_num-1)+2] = -l;
    ltri = n;
    ledg = 2;

    error = swapec ( m, &top, &ltri, &ledg, node_num, node_xy, *triangle_num,
      triangle_node, triangle_neighbor, stack );

    if ( error != 0 )
    {
      cout << "\n";
      cout << "R8TRIS2 - Fatal error!\n";
      cout << "  Error return from SWAPEC.\n";
      delete [] stack;
      return error;
    }

  }
//
//  Now account for the sorting that we did.
//
  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < *triangle_num; j++ )
    {
      triangle_node[i+j*3] = indx [ triangle_node[i+j*3] - 1 ];
    }
  }

  perm_inverse ( node_num, indx );

  r82vec_permute ( node_num, indx, base, node_xy );

  delete [] indx;
  delete [] stack;

  return 0;
}
//****************************************************************************80

void r8vec_bracket ( int n, double x[], double xval, int *left,
  int *right )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1];
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, double X[N], an array that has been sorted into ascending order.
//
//    Input, double XVAL, a value to be bracketed.
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
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
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = - r8_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
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

  value = r8_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Nijenhuis and Wilf,
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }

    k = k - 1;
    k1 = k;

  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ;; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

int swapec ( int i, int *top, int *btri, int *bedg, int node_num,
  double node_xy[], int triangle_num, int triangle_node[],
  int triangle_neighbor[], int stack[] )

//****************************************************************************80
//
//  Purpose:
//
//    SWAPEC swaps diagonal edges until all triangles are Delaunay.
//
//  Discussion:
//
//    The routine swaps diagonal edges in a 2D triangulation, based on
//    the empty circumcircle criterion, until all triangles are Delaunay,
//    given that I is the index of the new vertex added to the triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int I, the index of the new vertex.
//
//    Input/output, int *TOP, the index of the top of the stack.
//    On output, TOP is zero.
//
//    Input/output, int *BTRI, *BEDG; on input, if positive, are the
//    triangle and edge indices of a boundary edge whose updated indices
//    must be recorded.  On output, these may be updated because of swaps.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence
//    list.  May be updated on output because of swaps.
//
//    Input/output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor
//    list; negative values are used for links of the counter-clockwise linked
//    list of boundary edges;  May be updated on output because of swaps.
//
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
//    contain the indices of initial triangles (involving vertex I)
//    put in stack; the edges opposite I should be in interior;  entries
//    TOP+1 through MAXST are used as a stack.
//
//    Output, int SWAPEC, is set to 8 for abnormal return.
//
{
  int a;
  int b;
  int c;
  int e;
  int ee;
  int em1;
  int ep1;
  int f;
  int fm1;
  int fp1;
  int l;
  int r;
  int s;
  int swap;
  int t;
  int tt;
  int u;
  double x;
  double y;
//
//  Determine whether triangles in stack are Delaunay, and swap
//  diagonal edge of convex quadrilateral if not.
//
  x = node_xy[2*(i-1)+0];
  y = node_xy[2*(i-1)+1];

  for ( ; ; )
  {
    if ( *top <= 0 )
    {
      break;
    }

    t = stack[(*top)-1];
    *top = *top - 1;

    if ( triangle_node[3*(t-1)+0] == i )
    {
      e = 2;
      b = triangle_node[3*(t-1)+2];
    }
    else if ( triangle_node[3*(t-1)+1] == i )
    {
      e = 3;
      b = triangle_node[3*(t-1)+0];
    }
    else
    {
      e = 1;
      b = triangle_node[3*(t-1)+1];
    }

    a = triangle_node[3*(t-1)+e-1];
    u = triangle_neighbor[3*(t-1)+e-1];

    if ( triangle_neighbor[3*(u-1)+0] == t )
    {
      f = 1;
      c = triangle_node[3*(u-1)+2];
    }
    else if ( triangle_neighbor[3*(u-1)+1] == t )
    {
      f = 2;
      c = triangle_node[3*(u-1)+0];
    }
    else
    {
      f = 3;
      c = triangle_node[3*(u-1)+1];
    }

    swap = diaedg ( x, y,
      node_xy[2*(a-1)+0], node_xy[2*(a-1)+1],
      node_xy[2*(c-1)+0], node_xy[2*(c-1)+1],
      node_xy[2*(b-1)+0], node_xy[2*(b-1)+1] );

    if ( swap == 1 )
    {
      em1 = i4_wrap ( e - 1, 1, 3 );
      ep1 = i4_wrap ( e + 1, 1, 3 );
      fm1 = i4_wrap ( f - 1, 1, 3 );
      fp1 = i4_wrap ( f + 1, 1, 3 );

      triangle_node[3*(t-1)+ep1-1] = c;
      triangle_node[3*(u-1)+fp1-1] = i;
      r = triangle_neighbor[3*(t-1)+ep1-1];
      s = triangle_neighbor[3*(u-1)+fp1-1];
      triangle_neighbor[3*(t-1)+ep1-1] = u;
      triangle_neighbor[3*(u-1)+fp1-1] = t;
      triangle_neighbor[3*(t-1)+e-1] = s;
      triangle_neighbor[3*(u-1)+f-1] = r;

      if ( 0 < triangle_neighbor[3*(u-1)+fm1-1] )
      {
        *top = *top + 1;
        stack[(*top)-1] = u;
      }

      if ( 0 < s )
      {
        if ( triangle_neighbor[3*(s-1)+0] == u )
        {
          triangle_neighbor[3*(s-1)+0] = t;
        }
        else if ( triangle_neighbor[3*(s-1)+1] == u )
        {
          triangle_neighbor[3*(s-1)+1] = t;
        }
        else
        {
          triangle_neighbor[3*(s-1)+2] = t;
        }

        *top = *top + 1;

        if ( node_num < *top )
        {
          return 8;
        }

        stack[(*top)-1] = t;
      }
      else
      {
        if ( u == *btri && fp1 == *bedg )
        {
          *btri = t;
          *bedg = e;
        }

        l = - ( 3 * t + e - 1 );
        tt = t;
        ee = em1;

        while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
        {
          tt = triangle_neighbor[3*(tt-1)+ee-1];

          if ( triangle_node[3*(tt-1)+0] == a )
          {
            ee = 3;
          }
          else if ( triangle_node[3*(tt-1)+1] == a )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        triangle_neighbor[3*(tt-1)+ee-1] = l;

      }

      if ( 0 < r )
      {
        if ( triangle_neighbor[3*(r-1)+0] == t )
        {
          triangle_neighbor[3*(r-1)+0] = u;
        }
        else if ( triangle_neighbor[3*(r-1)+1] == t )
        {
          triangle_neighbor[3*(r-1)+1] = u;
        }
        else
        {
          triangle_neighbor[3*(r-1)+2] = u;
        }
      }
      else
      {
        if ( t == *btri && ep1 == *bedg )
        {
          *btri = u;
          *bedg = f;
        }

        l = - ( 3 * u + f - 1 );
        tt = u;
        ee = fm1;

        while ( 0 < triangle_neighbor[3*(tt-1)+ee-1] )
        {
          tt = triangle_neighbor[3*(tt-1)+ee-1];

          if ( triangle_node[3*(tt-1)+0] == b )
          {
            ee = 3;
          }
          else if ( triangle_node[3*(tt-1)+1] == b )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        triangle_neighbor[3*(tt-1)+ee-1] = l;

      }

    }

  }

  return 0;
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
//****************************************************************************80

double *triangle_angles_2d_new ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ANGLES_2D_NEW computes the angles of a triangle in 2D.
//
//  Discussion:
//
//    The law of cosines is used:
//
//      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
//
//    where GAMMA is the angle opposite side C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double ANGLE[3], the angles opposite
//    sides P1-P2, P2-P3 and P3-P1, in radians.
//
{
  double *angle;
  double a;
  double b;
  double c;
  double pi = 3.141592653589793;

  angle = new double[3];

  a = sqrt ( pow ( t[0+1*2] - t[0+0*2], 2 )
           + pow ( t[1+1*2] - t[1+0*2], 2 ) );

  b = sqrt ( pow ( t[0+2*2] - t[0+1*2], 2 )
           + pow ( t[1+2*2] - t[1+1*2], 2 ) );

  c = sqrt ( pow ( t[0+0*2] - t[0+2*2], 2 )
           + pow ( t[1+0*2] - t[1+2*2], 2 ) );
//
//  Take care of a ridiculous special case.
//
  if ( a == 0.0 && b == 0.0 && c == 0.0 )
  {
    angle[0] = 2.0 * pi / 3.0;
    angle[1] = 2.0 * pi / 3.0;
    angle[2] = 2.0 * pi / 3.0;
    return angle;
  }

  if ( c == 0.0 || a == 0.0 )
  {
    angle[0] = pi;
  }
  else
  {
    angle[0] = arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) );
  }

  if ( a == 0.0 || b == 0.0 )
  {
    angle[1] = pi;
  }
  else
  {
    angle[1] = arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) );
  }

  if ( b == 0.0 || c == 0.0 )
  {
    angle[2] = pi;
  }
  else
  {
    angle[2] = arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
  }

  return angle;
}
//****************************************************************************80

double triangle_area_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed
//    area of a triangle, it is possible to easily compute the area
//    of a nonconvex polygon as the sum of the (possibly negative)
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.
//
{
  double area;

  area = 0.5 * (
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) +
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) +
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );

  return area;
}
//****************************************************************************80

double *triangle_circumcenter_2d ( double t[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
//
//  Discussion:
//
//    The circumcenter of a triangle is the center of the circumcircle, the
//    circle that passes through the three vertices of the triangle.
//
//    The circumcircle contains the triangle, but it is not necessarily the
//    smallest triangle to do so.
//
//    If all angles of the triangle are no greater than 90 degrees, then
//    the center of the circumscribed circle will lie inside the triangle.
//    Otherwise, the center will lie outside the triangle.
//
//    The circumcenter is the intersection of the perpendicular bisectors
//    of the sides of the triangle.
//
//    In geometry, the circumcenter of a triangle is often symbolized by "O".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double *TRIANGLE_CIRCUMCENTER_2D[2], the circumcenter of
//    the triangle.
//
{
# define DIM_NUM 2

  double asq;
  double bot;
  double *center;
  double csq;
  double top1;
  double top2;

  center = new double[DIM_NUM];

  asq = ( t[0+1*2] - t[0+0*2] ) * ( t[0+1*2] - t[0+0*2] )
      + ( t[1+1*2] - t[1+0*2] ) * ( t[1+1*2] - t[1+0*2] );

  csq = ( t[0+2*2] - t[0+0*2] ) * ( t[0+2*2] - t[0+0*2] )
      + ( t[1+2*2] - t[1+0*2] ) * ( t[1+2*2] - t[1+0*2] );

  top1 =   ( t[1+1*2] - t[1+0*2] ) * csq - ( t[1+2*2] - t[1+0*2] ) * asq;
  top2 = - ( t[0+1*2] - t[0+0*2] ) * csq + ( t[0+2*2] - t[0+0*2] ) * asq;

  bot  =  ( t[1+1*2] - t[1+0*2] ) * ( t[0+2*2] - t[0+0*2] )
        - ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] );

  center[0] = t[0+0*2] + 0.5 * top1 / bot;
  center[1] = t[1+0*2] + 0.5 * top2 / bot;

  return center;

# undef DIM_NUM
}
//****************************************************************************80

void triangle_order3_physical_to_reference ( double t[], int n,
  double phy[], double ref[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE maps physical points to reference points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (X,Y) in the physical triangle, the routine computes the value
//    of the corresponding image point (XSI,ETA) in reference space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the X and Y coordinates
//    of the vertices.  The vertices are assumed to be the images of
//    (0,0), (1,0) and (0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double PHY[2*N], the coordinates of physical points
//    to be transformed.
//
//    Output, double REF[2*N], the coordinates of the corresponding
//    points in the reference space.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {

    ref[0+j*2] = ( ( t[1+2*2] - t[1+0*2] ) * ( phy[0+j*2] - t[0+0*2] )
                 - ( t[0+2*2] - t[0+0*2] ) * ( phy[1+j*2] - t[1+0*2] ) )
               / ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2]   - t[0+0*2] )
                 - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2]   - t[1+0*2] ) );

    ref[1+j*2] = ( ( t[0+1*2] - t[0+0*2] ) * ( phy[1+j*2] - t[1+0*2] )
                 - ( t[1+1*2] - t[1+0*2] ) * ( phy[0+j*2] - t[0+0*2] ) )
               / ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2]   - t[0+0*2] )
                 - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2]   - t[1+0*2] ) );
  }
  return;
}
//****************************************************************************80

void triangle_order3_reference_to_physical ( double t[], int n,
  double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL maps reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] )
                 + t[i+1*2] *       + ref[0+j*2]
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
//****************************************************************************80

void triangle_order6_physical_to_reference ( double t[2*6], int n,
  double phy[], double ref[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE maps a physical point to a reference point.
//
//  Discussion:
//
//    Given the vertices of an order 6 physical triangle and a point
//    (X,Y) in the physical triangle, the routine computes the value
//    of the corresponding image point (R,S) in reference space.
//
//    The mapping from (R,S) to (X,Y) has the form:
//
//      X(R,S) = A1 * R * R + B1 * R * S + C1 * S * S
//             + D1 * R     + E1 * S     + F1
//
//      Y(R,S) = A2 * R * R + B2 * R * S + C2 * S * S
//             + D2 * R     + E2 * S     + F2
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T(2,6), the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0), (0,1),
//    (1/2,0), (1/2,1/2) and (0,1/2), in that order.
//
//    Input, int N, the number of points to transform.
//
//    Input, double PHY(2,N), the coordinates of points in the
//    physical space.
//
//    Output, double REF(2,N), the coordinates of the corresponding
//    points in the reference space.
//
{
  double a[2];
  double b[2];
  double c[2];
  double d[2];
  double det;
  double dx[2];
  double e[2];
  double f[2];
  double fun[2];
  double fun_norm;
  int i;
  int it;
  int j;
  double jac[2*2];
  int it_max = 10;
  double it_tol = 0.000001;
//
//  Set iteration parameters.
//
  for ( i = 0; i < 2; i++ )
  {
    a[i] =   2.0 * t[i+0*2] + 2.0 * t[i+1*2] - 4.0 * t[i+3*2];
    b[i] =   4.0 * t[i+0*2] - 4.0 * t[i+3*2] + 4.0 * t[i+4*2] - 4.0 * t[i+5*2];
    c[i] =   2.0 * t[i+0*2] + 2.0 * t[i+2*2] - 4.0 * t[i+5*2];

    d[i] = - 3.0 * t[i+0*2] - t[i+1*2] + 4.0 * t[i+3*2];
    e[i] = - 3.0 * t[i+0*2] - t[i+2*2] + 4.0 * t[i+5*2];

    f[i] =   t[i+0*2];
  }
//
//  Initialize the points by inverting the linear map.
//
  triangle_order3_physical_to_reference ( t, n, phy, ref );
//
//  Carry out the Newton iteration.
//
  for ( j = 0; j < n; j++ )
  {
    for ( it = 0; it < it_max; it++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        fun[i] = a[i] * ref[0+j*2] * ref[0+j*2]
               + b[i] * ref[0+j*2] * ref[1+j*2]
               + c[i] * ref[1+j*2] * ref[1+j*2]
               + d[i] * ref[0+j*2]
               + e[i] * ref[1+j*2]
               + f[i]
               - phy[i+j*2];
      }

      fun_norm = sqrt ( pow ( fun[0], 2 ) + pow ( fun[1], 2 ) );

      if ( fun_norm <= it_tol )
      {
        break;
      }

      jac[0+0*2] = 2.0 * a[0] * ref[0+j*2] + b[0] * ref[1+j*2] + d[0];
      jac[1+0*2] = 2.0 * a[1] * ref[0+j*2] + b[1] * ref[1+j*2] + d[1];
      jac[0+1*2] = b[0] * ref[0+j*2] + 2.0 * c[0] * ref[1+j*2] + e[0];
      jac[1+1*2] = b[1] * ref[0+j*2] + 2.0 * c[1] * ref[1+j*2] + e[1];

      det = jac[0+0*2] * jac[1+1*2] - jac[0+1*2] * jac[1+0*2];

      if ( det == 0.0 )
      {
        cout << "\n";
        cout << "TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE - Fatal error!\n";
        cout << "  The jacobian of the mapping is singular.\n";
      }

      dx[0] = (  jac[1+1*2] * fun[0] - jac[0+1*2] * fun[1] ) / det;
      dx[1] = ( -jac[1+0*2] * fun[0] + jac[0+0*2] * fun[1] ) / det;

      ref[0+j*2] = ref[0+j*2] - dx[0];
      ref[1+j*2] = ref[1+j*2] - dx[1];
    }
  }

  return;
}
//****************************************************************************80

void triangle_order6_reference_to_physical ( double t[], int n,
  double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL maps reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 6 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    The mapping from (XSI,ETA) to (X,Y) has the form:
//
//      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
//                 + D1 * XSI    + E1 * ETA     + F1
//
//      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
//                 + D2 * XSI    + E2 * ETA     + F2
//
//  Reference Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0),
//    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
//
//    Input, integer N, the number of points to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  double a[2];
  double b[2];
  double c[2];
  double d[2];
  double e[2];
  double f[2];
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    a[i] =   2.0 * t[i+0*2] + 2.0 * t[i+1*2]
           - 4.0 * t[i+3*2];

    b[i] =   4.0 * t[i+0*2]
           - 4.0 * t[i+3*2] + 4.0 * t[i+4*2] - 4.0 * t[i+5*2];

    c[i] =   2.0 * t[i+0*2]                  + 2.0 * t[i+2*2]
                                             - 4.0 * t[i+5*2];

    d[i] = - 3.0 * t[i+0*2] -       t[i+1*2]
           + 4.0 * t[i+3*2];

    e[i] = - 3.0 * t[i+0*2]                  -       t[i+2*2]
                                             + 4.0 * t[i+5*2];
    f[i] =         t[i+0*2];

  }

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = a[i] * ref[0+j*2] * ref[0+j*2]
                 + b[i] * ref[0+j*2] * ref[1+j*2]
                 + c[i] * ref[1+j*2] * ref[1+j*2]
                 + d[i] * ref[0+j*2]
                 + e[i] * ref[1+j*2]
                 + f[i];
    }
  }

  return;
}
//****************************************************************************80

void triangle_reference_sample ( int n, int *seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_REFERENCE_SAMPLE returns random points in the reference triangle.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  |   \  side 2
//       |    \
//    3  |     \
//       |      \
//       1-------2
//
//         side 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[2*N], a random point in the triangle.
//
{
# define DIM_NUM 2

  double alpha;
  double beta;
  int j;
  double r;

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    alpha = sqrt ( r );
//
//  Now choose, uniformly at random, a point on the line L.
//
    beta = r8_uniform_01 ( seed );

    p[0+j*2] = ( 1.0 - beta ) * alpha;
    p[1+j*2] =         beta   * alpha;
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void triangle_sample ( double t[2*3], int n, int *seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE returns random points in a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[2*N], a random point in the triangle.
//
{
# define DIM_NUM 2

  double alpha;
  double beta;
  int j;
  double r;
  double p12[DIM_NUM];
  double p13[DIM_NUM];

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    alpha = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
    p12[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+1*2];
    p12[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+1*2];

    p13[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+2*2];;
    p13[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+2*2];;
//
//  Now choose, uniformly at random, a point on the line L.
//
    beta = r8_uniform_01 ( seed );

    p[0+j*2] = ( 1.0 - beta ) * p12[0] + beta * p13[0];
    p[1+j*2] = ( 1.0 - beta ) * p12[1] + beta * p13[1];
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

double triangulation_area ( int node_num, double node_xy[], int element_order,
  int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_AREA computes the area of a triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int ELEMENT_ORDER, the order of the triangles.
//
//    Input, int ELEMENT_NUM, the number of triangles.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM],
//    the nodes making up each triangle.
//
//    Output, double TRIANGULATION_AREA, the area.
//
{
  int element;
  double element_area;
  double element_xy[2*3];
  int j;
  int nj;
  double value;

  value = 0.0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      nj = element_node[j+element*element_order];
      element_xy[0+j*2] = node_xy[0+nj*2];
      element_xy[1+j*2] = node_xy[1+nj*2];
    }

    element_area = 0.5 * (
      element_xy[0+0*2] * ( element_xy[1+1*2] - element_xy[1+2*2] ) +
      element_xy[0+1*2] * ( element_xy[1+2*2] - element_xy[1+0*2] ) +
      element_xy[0+2*2] * ( element_xy[1+0*2] - element_xy[1+1*2] ) );

    value = value + element_area;
  }
  return value;
}
//****************************************************************************80

double triangulation_areas ( int node_num, double node_xy[], int triangle_order,
  int triangle_num, int triangle_node[], double triangle_area[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_AREAS computes triangle and triangulation areas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes in the
//    triangulation.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_ORDER, the order of triangles in
//    the triangulation.
//
//    Input, int TRIANGLE_NUM, the number of triangles in
//    the triangulation.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the nodes making up each triangle.
//
//    Output, double TRIANGLE_AREA[TRIANGLE_NUM], the area of
//    the triangles.
//
//    Output, double TRIANGULATION_AREAS, the area of the triangulation.
//
{
  int j;
  int nj;
  double t_xy[2*3];
  int triangle;
  double triangulation_area;

  triangulation_area = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      nj = triangle_node[j+triangle*triangle_order];
      t_xy[0+j*2] = node_xy[0+nj*2];
      t_xy[1+j*2] = node_xy[1+nj*2];
    }

    triangle_area[triangle] = 0.5 * (
      t_xy[0+0*2] * ( t_xy[1+1*2] - t_xy[1+2*2] ) +
      t_xy[0+1*2] * ( t_xy[1+2*2] - t_xy[1+0*2] ) +
      t_xy[0+2*2] * ( t_xy[1+0*2] - t_xy[1+1*2] ) );

    triangulation_area = triangulation_area + triangle_area[triangle];
  }

  return triangulation_area;
}
//****************************************************************************80

double triangulation_delaunay_discrepancy_compute ( int node_num,
  double node_xy[], int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double *angle_min, int *angle_min_triangle,
  double *angle_max, int *angle_max_triangle )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE reports if a triangulation is Delaunay.
//
//  Discussion:
//
//    A (maximal) triangulation is Delaunay if and only if it is locally
//    Delaunay.
//
//    A triangulation is Delaunay if the minimum angle over all triangles
//    in the triangulation is maximized.  That is, there is no other
//    triangulation of the points which has a larger minimum angle.
//
//    A triangulation is locally Delaunay if, for every pair of triangles
//    that share an edge, the minimum angle in the two triangles is larger
//    than the minimum angle in the two triangles formed by removing the
//    common edge and joining the two opposing vertices.
//
//    This function examines the question of whether a given triangulation
//    is locally Delaunay.  It does this by looking at every pair of
//    neighboring triangles and comparing the minimum angle attained
//    for the current triangle pair and the alternative triangle pair.
//
//    Let A(i,j) be the minimum angle formed by triangles T(i) and T(j),
//    which are two triangles in the triangulation which share a common edge.
//    Let B(I,J) be the minimum angle formed by triangles S(i) and S(j),
//    where S(i) and S(j) are formed by removing the common edge of T(i)
//    and T(j), and joining the opposing vertices.
//
//    Then the triangulation is Delaunay if B(i,j) <= A(i,j) for every
//    pair of neighbors T(i) and T(j).
//
//    If A(i,j) < B(i,j) for at least one pair of neighbors, the triangulation
//    is not a Delaunay triangulation.
//
//    This program returns VALUE = min ( A(i,j) - B(i,j) ) over all
//    triangle neighbors.  VALUE is scaled to be in degrees, for
//    comprehensibility.  If VALUE is negative, then at least one pair
//    of triangles violates the Delaunay condition, and so the entire
//    triangulation is not a Delaunay triangulation.  If VALUE is nonnegative,
//    then the triangulation is a Delaunay triangulation.
//
//    It is useful to return VALUE, rather than a simple True/False value,
//    because there can be cases where the Delaunay condition is only
//    "slightly" violated.  A simple example is a triangulation formed
//    by starting with a mesh of squares and dividing each square into
//    two triangles by choosing one of the diagonals of the square.
//    The Delaunay discrepancy for this mesh, if computed exactly, is 0,
//    but roundoff could easily result in discrepancies that were very
//    slightly negative.
//
//    If VALUE is positive, and not very small in magnitude, then every
//    pair of triangles in the triangulation satisfies the local Delaunay
//    condition, and so the triangulation is a Delaunay triangulation.
//
//    If VALUE is negative, and not very small in magnitude, then at least
//    one pair of triangles violates the Delaunay condition, and to a
//    significant degree.  The triangulation is not a Delaunay triangulation.
//
//    If the magnitude of VALUE is very close to zero, then the triangulation
//    is numerically ambiguous.  At least one pair of triangles violates
//    or almost violates the condition, but no triangle violates the
//    condition to a great extent.  The user must judge whether the
//    violation is significant or not.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles in
//    the triangulation.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the nodes that make up each triangle.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the
//    triangle neighbor list.
//
//    Output, double *ANGLE_MIN, the minimum angle that occurred in
//    the triangulation.
//
//    Output, int *ANGLE_MIN_TRIANGLE, the triangle in which
//    the minimum angle occurred.
//
//    Output, double *ANGLE_MAX, the maximum angle that occurred in
//    the triangulation.
//
//    Output, int *ANGLE_MAX_TRIANGLE, the triangle in which
//    the maximum angle occurred.
//
//    Output, double TRIANGULATION_DELAUNAY_DISCREPANCY,
//    the minimum value of ( A(i,j) - B(i,j) ).
//    POSITIVE indicates the triangulation is Delaunay.
//    VERY NEAR ZERO is a numerically ambiguous case.
//    NEGATIVE indicates the triangulation is not Delaunay.
//
{
  double angle_min1;
  double angle_min2;
  double *angles1;
  double *angles2;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int n1;
  int n2;
  int n3;
  int n4;
  int neighbor;
  double pi = 3.141592653589793;
  double t[2*3];
  int triangle_index;
  int triangle1;
  int triangle2;
  double value;

  *angle_max = 0.0;
  *angle_max_triangle = - 1;
  *angle_min = pi;
  *angle_min_triangle = -1;
  value = 0.0;
//
//  Consider triangle TRIANGLE1
//
  for ( triangle1 = 0; triangle1 < triangle_num; triangle1++ )
  {
//
//  Consider the side opposite vertex NEIGHBOR.
//
    for ( neighbor = 0; neighbor < 3; neighbor++ )
    {
      triangle2 = triangle_neighbor[neighbor+triangle1*3];
//
//  There might be no neighbor on side NEIGHBOR.
//
      if ( triangle2 < 0 )
      {
        continue;
      }
//
//  We only need to check a pair of triangles once.
//
      if ( triangle2 < triangle1 )
      {
        continue;
      }
//
//  List the vertices of the quadrilateral in such a way
//  that the nodes of triangle 1 come first.
//
//  We rely on a property of the TRIANGLE_NEIGHBOR array, namely, that
//  neighbor #1 is on the side opposite to vertex #1, and so on.
//
      i1 = i4_wrap ( neighbor + 2, 0, 2 );
      i2 = i4_wrap ( neighbor,     0, 2 );
      i3 = i4_wrap ( neighbor + 1, 0, 2 );

      n1 = triangle_node[i1+triangle1*triangle_order];
      n2 = triangle_node[i2+triangle1*triangle_order];
      n3 = triangle_node[i3+triangle1*triangle_order];
//
//  The "odd" or "opposing" node of the neighboring triangle
//  is the one which follows common node I3.
//
      n4 = -1;
      for ( i = 0; i < 3; i++ )
      {
        if ( triangle_node[i+triangle2*triangle_order] == n3 )
        {
          i4 = i + 1;
          i4 = i4_wrap ( i4, 0, 2 );
          n4 = triangle_node[i4+triangle2*triangle_order];
          break;
        }
      }

      if ( n4 == -1 )
      {
        cout << "\n";
        cout << "TRIANGULATION_DELAUNAY_DISCREPANCY_COMPUTE - Fatal error/!\n";
        cout << "  Could not identify the fourth node.\n";
        cout << "\n";
        cout << "  Triangle1 = " << triangle1 << "\n";
        cout << "  Nodes =     ";
        for ( i = 0; i < 3; i++ )
        {
          cout << "  " << triangle_node[i+triangle1*triangle_order];
        }
        cout << "\n";
        cout << "  Neighbors =     ";
        for ( i = 0; i < 3; i++ )
        {
          cout << "  " << triangle_neighbor[i+triangle1*3];
        }
        cout << "\n";
        cout << "\n";
        cout << "  Neighbor index = " << neighbor << "\n";
        cout << "\n";
        cout << "  Triangle2 = " << triangle2 << "\n";
        cout << "  Nodes =     ";
        for ( i = 0; i < 3; i++ )
        {
          cout << "  " << triangle_node[i+triangle2*triangle_order];
        }
        cout << "\n";
        cout << "  Neighbors =     ";
        for ( i = 0; i < 3; i++ )
        {
          cout << "  " << triangle_neighbor[i+triangle2*3];
        }
        cout << "\n";
        exit ( 1 );
      }
//
//  Compute the minimum angle for (I1,I2,I3) and (I1,I3,I4).
//
      t[0+0*2] = node_xy[0+n1*2];
      t[1+0*2] = node_xy[1+n1*2];
      t[0+1*2] = node_xy[0+n2*2];
      t[1+1*2] = node_xy[1+n2*2];
      t[0+2*2] = node_xy[0+n3*2];
      t[1+2*2] = node_xy[1+n3*2];
      angles1 = triangle_angles_2d_new ( t );

      t[0+0*2] = node_xy[0+n1*2];
      t[1+0*2] = node_xy[1+n1*2];
      t[0+1*2] = node_xy[0+n3*2];
      t[1+1*2] = node_xy[1+n3*2];
      t[0+2*2] = node_xy[0+n4*2];
      t[1+2*2] = node_xy[1+n4*2];
      angles2 = triangle_angles_2d_new ( t );

      angle_min1 =
        r8_min ( r8vec_min ( 3, angles1 ), r8vec_min ( 3, angles2 ) );

      if ( *angle_max < r8vec_max ( 3, angles1 ) )
      {
        *angle_max = r8vec_max ( 3, angles1 );
        *angle_max_triangle = triangle1;
      }

      if ( *angle_max < r8vec_max ( 3, angles2 ) )
      {
        *angle_max = r8vec_max ( 3, angles2 );
        *angle_max_triangle = triangle2;
      }

      if ( r8vec_min ( 3, angles1 ) < *angle_min )
      {
        *angle_min = r8vec_min ( 3, angles1 );
        *angle_min_triangle = triangle1;
      }

      if ( r8vec_min ( 3, angles2 ) < *angle_min )
      {
        *angle_min = r8vec_min ( 3, angles2 );
        *angle_min_triangle = triangle2;
      }

     delete [] angles1;
     delete [] angles2;
//
//  Compute the minimum angle for (I1,I2,I4) and (I2,I3,I4).
//
      t[0+0*2] = node_xy[0+n1*2];
      t[1+0*2] = node_xy[1+n1*2];
      t[0+1*2] = node_xy[0+n2*2];
      t[1+1*2] = node_xy[1+n2*2];
      t[0+2*2] = node_xy[0+n4*2];
      t[1+2*2] = node_xy[1+n4*2];
      angles1 = triangle_angles_2d_new ( t );

      t[0+0*2] = node_xy[0+n3*2];
      t[1+0*2] = node_xy[1+n3*2];
      t[0+1*2] = node_xy[0+n3*2];
      t[1+1*2] = node_xy[1+n3*2];
      t[0+2*2] = node_xy[0+n4*2];
      t[1+2*2] = node_xy[1+n4*2];
      angles2 = triangle_angles_2d_new ( t );

      angle_min2 =
        r8_min ( r8vec_min ( 3, angles1 ), r8vec_min ( 3, angles2 ) );

     delete [] angles1;
     delete [] angles2;
//
//  Compare this value to the current minimum.
//
      value = r8_min ( value, angle_min1 - angle_min2 );
    }
  }
//
//  Scale the results to degrees.
//
  value = value * 180.0 / pi;
  *angle_max = *angle_max * 180.0 / pi;
  *angle_min = *angle_min * 180.0 / pi;

  return value;
}
//****************************************************************************80

int *triangulation_neighbor_elements ( int triangle_order, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_NEIGHBOR_ELEMENTS determines element neighbors.
//
//  Discussion:
//
//    A triangulation of a set of nodes can be completely described by
//    the coordinates of the nodes, and the list of nodes that make up
//    each triangle.  However, in some cases, it is necessary to know
//    triangle adjacency information, that is, which triangle, if any,
//    is adjacent to a given triangle on a particular side.
//
//    This routine creates a data structure recording this information.
//
//    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
//    data items.
//
//    This routine was modified to work with columns rather than rows.
//
//  Example:
//
//    The input information from TRIANGLE_NODE:
//
//    Triangle   Nodes
//    --------   ---------------
//     1         3      4      1
//     2         3      1      2
//     3         3      2      8
//     4         2      1      5
//     5         8      2     13
//     6         8     13      9
//     7         3      8      9
//     8        13      2      5
//     9         9     13      7
//    10         7     13      5
//    11         6      7      5
//    12         9      7      6
//    13        10      9      6
//    14         6      5     12
//    15        11      6     12
//    16        10      6     11
//
//    The output information in TRIANGLE_NEIGHBOR:
//
//    Triangle  Neighboring Triangles
//    --------  ---------------------
//
//     1        -1     -1      2
//     2         1      4      3
//     3         2      5      7
//     4         2     -1      8
//     5         3      8      6
//     6         5      9      7
//     7         3      6     -1
//     8         5      4     10
//     9         6     10     12
//    10         9      8     11
//    11        12     10     14
//    12         9     11     13
//    13        -1     12     16
//    14        11     -1     15
//    15        16     14     -1
//    16        13     15     -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes that
//    make up each triangle.
//
//    Output, int TRIANGLE_NEIGHBOR_TRIANGLES[3*TRIANGLE_NUM],
//    the three triangles
//    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I)
//    is the index of the triangle which touches side 1, defined by nodes 2
//    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no
//    neighbor on that side.  In this case, that side of the triangle lies
//    on the boundary of the triangulation.
//
{
  int *col;
  int i;
  int icol;
  int j;
  int k;
  int side1;
  int side2;
  int tri;
  int tri1;
  int tri2;
  int *triangle_neighbor;

  triangle_neighbor = new int[3*triangle_num];
  col = new int[4*(3*triangle_num)];
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the three neighbor relations:
//
//    (I,J,3,T) or (J,I,3,T),
//    (J,K,1,T) or (K,J,1,T),
//    (K,I,2,T) or (I,K,2,T)
//
//  where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
//
  for ( tri = 0; tri < triangle_num; tri++ )
  {
    i = triangle_node[0+tri*triangle_order];
    j = triangle_node[1+tri*triangle_order];
    k = triangle_node[2+tri*triangle_order];

    if ( i < j )
    {
      col[0+(3*tri+0)*4] = i;
      col[1+(3*tri+0)*4] = j;
      col[2+(3*tri+0)*4] = 3;
      col[3+(3*tri+0)*4] = tri + 1;
    }
    else
    {
      col[0+(3*tri+0)*4] = j;
      col[1+(3*tri+0)*4] = i;
      col[2+(3*tri+0)*4] = 3;
      col[3+(3*tri+0)*4] = tri + 1;
    }

    if ( j < k )
    {
      col[0+(3*tri+1)*4] = j;
      col[1+(3*tri+1)*4] = k;
      col[2+(3*tri+1)*4] = 1;
      col[3+(3*tri+1)*4] = tri + 1;
    }
    else
    {
      col[0+(3*tri+1)*4] = k;
      col[1+(3*tri+1)*4] = j;
      col[2+(3*tri+1)*4] = 1;
      col[3+(3*tri+1)*4] = tri + 1;
    }

    if ( k < i )
    {
      col[0+(3*tri+2)*4] = k;
      col[1+(3*tri+2)*4] = i;
      col[2+(3*tri+2)*4] = 2;
      col[3+(3*tri+2)*4] = tri + 1;
    }
    else
    {
      col[0+(3*tri+2)*4] = i;
      col[1+(3*tri+2)*4] = k;
      col[2+(3*tri+2)*4] = 2;
      col[3+(3*tri+2)*4] = tri + 1;
    }
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1 and 2; the routine we call here
//  sorts on rows 1 through 4 but that won't hurt us.
//
//  What we need is to find cases where two triangles share an edge.
//  Say they share an edge defined by the nodes I and J.  Then there are
//  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
//  we make sure that these two columns occur consecutively.  That will
//  make it easy to notice that the triangles are neighbors.
//
  i4col_sort_a ( 4, 3*triangle_num, col );
//
//  Step 3. Neighboring triangles show up as consecutive columns with
//  identical first two entries.  Whenever you spot this happening,
//  make the appropriate entries in TRIANGLE_NEIGHBOR.
//
  for ( j = 0; j < triangle_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_neighbor[i+j*3] = -1;
    }
  }

  icol = 1;

  for ( ; ; )
  {
    if ( 3 * triangle_num <= icol )
    {
      break;
    }

    if ( col[0+(icol-1)*4] != col[0+icol*4] ||
         col[1+(icol-1)*4] != col[1+icol*4] )
    {
      icol = icol + 1;
      continue;
    }

    side1 = col[2+(icol-1)*4];
    tri1 =  col[3+(icol-1)*4];
    side2 = col[2+ icol   *4];
    tri2 =  col[3+ icol   *4];

    triangle_neighbor[side1-1+(tri1-1)*3] = tri2;
    triangle_neighbor[side2-1+(tri2-1)*3] = tri1;

    icol = icol + 2;
  }

  delete [] col;

  return triangle_neighbor;
}
//****************************************************************************80

int *triangulation_node_order ( int triangle_order, int triangle_num,
  int triangle_node[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_NODE_ORDER determines the order of nodes in a triangulation.
//
//  Discussion:
//
//    The order of a node is the number of triangles that use that node
//    as a vertex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer TRIANGLE_ORDER, the order of the triangulation.
//
//    Input, integer TRIANGLE_NUM, the number of triangles.
//
//    Input, integer TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes
//    that make up the triangles.
//
//    Input, integer NODE_NUM, the number of nodes.
//
//    Output, integer TRIANGULATION_NODE_ORDER[NODE_NUM], the order of
//    each node.
//
{
  int i;
  int node;
  int *node_order;
  int triangle;

  node_order = new int[node_num];

  for ( node = 0; node < node_num; node++ )
  {
    node_order[node] = 0;
  }

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    for ( i = 0; i < triangle_order; i++ )
    {
      node = triangle_node[i+triangle*triangle_order];

      if ( node < 1 || node_num < node )
      {
        cout << "\n";
        cout << "TRIANGULATION_NODE_ORDER - Fatal error!\n";
        cout << "  Illegal entry in TRIANGLE_NODE.\n";
        node_order = NULL;
        exit ( 1 );
      }
      else
      {
        node_order[node-1] = node_order[node-1] + 1;
      }
    }
  }

  return node_order;
}
//****************************************************************************80

int triangulation_order3_adj_count ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_col[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 3-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  |   \  side 2
//       |    \
//    3  |     \
//       |      \
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//    A sample grid.
//
//
//    Below, we have a chart that summarizes the adjacency relationships
//    in the sample grid.  On the left, we list the node, and its neighbors,
//    with an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).  On the right, we list the number of adjancencies to
//    lower-indexed nodes, to the node itself, to higher-indexed nodes,
//    the total number of adjacencies for this node, and the location
//    of the first and last entries required to list this set of adjacencies
//    in a single list of all the adjacencies.
//
//    N   Adjacencies                Below  Self   Above   Total First  Last
//
//   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
//    1:  *  2  6                        0     1       2       3     1     3
//    2:  1  *  3  6  7                  1     1       3       5     4     8
//    3:  2  *  4  7  8                  1     1       3       5     9    13
//    4:  3  *  5  8  9                  1     1       3       5    14    18
//    5:  4  *  9 10                     1     1       2       4    19    22
//    6:  1  2  *  7 11                  2     1       2       5    23    27
//    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
//    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
//    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
//   10:  5  9  * 14 15                  2     1       2       5    49    53
//   11:  6  7  * 12 16                  2     1       2       5    54    58
//   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
//   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
//   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
//   15: 10 14  * 19 20                  2     1       2       5    80    84
//   16: 11 12  * 17 21                  2     1       2       5    85    89
//   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
//   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
//   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
//   20: 15 19  * 24 25                  2     1       2       5   111   115
//   21: 16 17  * 22                     2     1       1       4   116   119
//   22: 17 18 21  * 23                  3     1       1       5   120   124
//   23: 18 19 22  * 24                  3     1       1       5   125   129
//   24: 19 20 23  * 25                  3     1       1       5   130   134
//   25: 20 24  *                        2     1       0       3   135   137
//   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
//    make up each triangle, in counterclockwise order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Output, TRIANGULATION_ORDER3_ADJ_COUNT, the number of adjacencies.
//
//    Output, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
{
  int adj_num;
  int i;
  int n1;
  int n2;
  int n3;
  int node;
  int triangle;
  int triangle_order = 3;
  int triangle2;

  adj_num = 0;
//
//  Set every node to be adjacent to itself.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_col[node] = 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
    }
//
//  Add edge (2,3).
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
    }
//
//  Add edge (3,1).
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
    }
  }
//
//  We used ADJ_COL to count the number of entries in each column.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num; 1 <= node; node-- )
  {
    adj_col[node] = adj_col[node-1];
  }
  adj_col[0] = 1;
  for ( i = 1; i <= node_num; i++ )
  {
    adj_col[i]= adj_col[i-1] + adj_col[i];
  }

  adj_num = adj_col[node_num] - 1;

  return adj_num;
}
//****************************************************************************80

int *triangulation_order3_adj_set ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_num, int adj_col[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to set the adjacencies, after the
//    appropriate amount of memory has been set aside for storage.
//
//    The triangulation is assumed to involve 3-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//    This routine can be used to create the compressed column storage
//    for a linear triangle finite element discretization of
//    Poisson's equation in two dimensions.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  |   \  side 2
//       |    \
//    3  |     \
//       |      \
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//    A sample grid
//
//
//    Below, we have a chart that summarizes the adjacency relationships
//    in the sample grid.  On the left, we list the node, and its neighbors,
//    with an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).  On the right, we list the number of adjancencies to
//    lower-indexed nodes, to the node itself, to higher-indexed nodes,
//    the total number of adjacencies for this node, and the location
//    of the first and last entries required to list this set of adjacencies
//    in a single list of all the adjacencies.
//
//    N   Adjacencies                Below  Self    Above  Total First  Last
//
//   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
//    1:  *  2  6                        0     1       2       3     1     3
//    2:  1  *  3  6  7                  1     1       3       5     4     8
//    3:  2  *  4  7  8                  1     1       3       5     9    13
//    4:  3  *  5  8  9                  1     1       3       5    14    18
//    5:  4  *  9 10                     1     1       2       4    19    22
//    6:  1  2  *  7 11                  2     1       2       5    23    27
//    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
//    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
//    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
//   10:  5  9  * 14 15                  2     1       2       5    49    53
//   11:  6  7  * 12 16                  2     1       2       5    54    58
//   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
//   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
//   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
//   15: 10 14  * 19 20                  2     1       2       5    80    84
//   16: 11 12  * 17 21                  2     1       2       5    85    89
//   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
//   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
//   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
//   20: 15 19  * 24 25                  2     1       2       5   111   115
//   21: 16 17  * 22                     2     1       1       4   116   119
//   22: 17 18 21  * 23                  3     1       1       5   120   124
//   23: 18 19 22  * 24                  3     1       1       5   125   129
//   24: 19 20 23  * 25                  3     1       1       5   130   134
//   25: 20 24  *                        2     1       0       3   135   137
//   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
//    make up each triangle in counterclockwise order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, int ADJ_NUM, the number of adjacencies.
//
//    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
//    Output, int TRIANGULATION_ORDER3_ADJ_SET[ADJ_NUM], the adjacency
//    information.
//
{
  int *adj;
  int *adj_copy;
  int k;
  int k1;
  int k2;
  int n1;
  int n2;
  int n3;
  int node;
  int triangle;
  int triangle2;
  int triangle_order = 3;

  adj = new int[adj_num];
  for ( k = 0; k < adj_num; k++ )
  {
    adj[k] = -1;
  }

  adj_copy = new int[node_num];
  for ( node = 0; node < node_num; node++ )
  {
    adj_copy[node] = adj_col[node];
  }
//
//  Set every node to be adjacent to itself.
//
  for ( node = 1; node <= node_num; node++ )
  {
    adj[adj_copy[node-1]-1] = node;
    adj_copy[node-1] = adj_copy[node-1] + 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n2;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n2-1]-1] = n1;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
    }
//
//  Add edge (2,3).
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n2-1]-1] = n3;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n3-1]-1] = n2;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
    }
//
//  Add edge (3,1).
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n3;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n3-1]-1] = n1;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
    }
  }
//
//  Ascending sort the entries for each node.
//
  for ( node = 1; node <= node_num; node++ )
  {
    k1 = adj_col[node-1];
    k2 = adj_col[node]-1;
    i4vec_sort_heap_a ( k2+1-k1, adj+k1-1 );
  }

  delete [] adj_copy;

  return adj;
}
//****************************************************************************80

void triangulation_order3_adj_set2 ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_num, int adj_col[],
  int ia[], int ja[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to set up the arrays IA and JA that
//    record which nodes are adjacent in a triangulation.
//
//    The triangulation is assumed to involve 3-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//    This routine can be used to create the compressed column storage
//    for a linear triangle finite element discretization of
//    Poisson's equation in two dimensions.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  |   \  side 2
//       |    \
//    3  |     \
//       |      \
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//    A sample grid
//
//
//    Below, we have a chart that summarizes the adjacency relationships
//    in the sample grid.  On the left, we list the node, and its neighbors,
//    with an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).  On the right, we list the number of adjancencies to
//    lower-indexed nodes, to the node itself, to higher-indexed nodes,
//    the total number of adjacencies for this node, and the location
//    of the first and last entries required to list this set of adjacencies
//    in a single list of all the adjacencies.
//
//    N   Adjacencies                Below  Self    Above  Total First  Last
//
//   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
//    1:  *  2  6                        0     1       2       3     1     3
//    2:  1  *  3  6  7                  1     1       3       5     4     8
//    3:  2  *  4  7  8                  1     1       3       5     9    13
//    4:  3  *  5  8  9                  1     1       3       5    14    18
//    5:  4  *  9 10                     1     1       2       4    19    22
//    6:  1  2  *  7 11                  2     1       2       5    23    27
//    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
//    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
//    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
//   10:  5  9  * 14 15                  2     1       2       5    49    53
//   11:  6  7  * 12 16                  2     1       2       5    54    58
//   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
//   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
//   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
//   15: 10 14  * 19 20                  2     1       2       5    80    84
//   16: 11 12  * 17 21                  2     1       2       5    85    89
//   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
//   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
//   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
//   20: 15 19  * 24 25                  2     1       2       5   111   115
//   21: 16 17  * 22                     2     1       1       4   116   119
//   22: 17 18 21  * 23                  3     1       1       5   120   124
//   23: 18 19 22  * 24                  3     1       1       5   125   129
//   24: 19 20 23  * 25                  3     1       1       5   130   134
//   25: 20 24  *                        2     1       0       3   135   137
//   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
//
//    For this example, the initial portion of the IA and JA arrays will be:
//
//      (1,1), (1,2), (1,6),
//      (2,1), (2,2), (2,3), (2,6), (2,7),
//      (3,2), (3,3), (3,4), (3,7), (3,8),
//     ...
//      (25,20), (25,24), (25,25)
//
//    for a total of 137 pairs of values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
//    make up each triangle in counterclockwise order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, int ADJ_NUM, the number of adjacencies.
//
//    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
//    Output, int IA[ADJ_NUM], JA[ADJ_NUM], the adjacency information.
//
{
  int adj;
  int *adj_copy;
  int k;
  int k1;
  int k2;
  int n1;
  int n2;
  int n3;
  int node;
  int triangle;
  int triangle2;
  int triangle_order = 3;

  for ( adj = 0; adj < adj_num; adj++ )
  {
    ia[adj] = -1;
  }

  for ( adj = 0; adj < adj_num; adj++ )
  {
    ja[adj] = -1;
  }

  adj_copy = new int[node_num];
  for ( node = 0; node < node_num; node++ )
  {
    adj_copy[node] = adj_col[node];
  }
//
//  Set every node to be adjacent to itself.
//
  for ( node = 1; node <= node_num; node++ )
  {
    ia[adj_copy[node-1]-1] = node;
    ja[adj_copy[node-1]-1] = node;
    adj_copy[node-1] = adj_copy[node-1] + 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      ia[adj_copy[n1-1]-1] = n1;
      ja[adj_copy[n1-1]-1] = n2;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;

      ia[adj_copy[n2-1]-1] = n2;
      ja[adj_copy[n2-1]-1] = n1;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
    }
//
//  Add edge (2,3).
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      ia[adj_copy[n2-1]-1] = n2;
      ja[adj_copy[n2-1]-1] = n3;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;

      ia[adj_copy[n3-1]-1] = n3;
      ja[adj_copy[n3-1]-1] = n2;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
    }
//
//  Add edge (3,1).
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      ia[adj_copy[n1-1]-1] = n1;
      ja[adj_copy[n1-1]-1] = n3;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;

      ia[adj_copy[n3-1]-1] = n3;
      ja[adj_copy[n3-1]-1] = n1;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
    }
  }
//
//  Lexically sort the IA, JA values.
//
  i4vec2_sort_a ( adj_num, ia, ja );

  delete [] adj_copy;

  return;
}
//****************************************************************************80

int *triangulation_order3_adjacency ( int node_num, int element_num, 
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_ADJACENCY computes the full adjacency matrix
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes in the
//    triangulation.
//
//    Input, int ELEMENT_NUM, the number of triangles in
//    the triangulation.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM],
//    the nodes making up each triangle.
//
//    Output, int TRIANGULATION_ORDER3_ADJACENCY[NODE_NUM*NODE_NUM], the adjacency
//    matrix.  ADJ(I,J) is 1 if nodes I and J are adjacent, that is,
//    they are immediate neighbors on an edge of the triangulation.
//
{
  int *adj;
  int element;
  int i;
  int j;
  int k;

  adj = new int[node_num*node_num];

  for ( j = 0; j < node_num; j++ )
  {
    for ( i = 0; i < node_num; i++ )
    {
      adj[i+j*node_num] = 0;
    }
  }

  for ( element = 0; element < element_num; element++ )
  {
    i = element_node[0+element*3];
    j = element_node[1+element*3];
    k = element_node[2+element*3];

    adj[i+j*node_num] = 1;
    adj[i+k*node_num] = 1;
    adj[j+i*node_num] = 1;
    adj[j+k*node_num] = 1;
    adj[k+i*node_num] = 1;
    adj[k+j*node_num] = 1;
  }

  return adj;
}
//****************************************************************************80

int triangulation_order3_boundary_edge_count ( int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the boundary edges.
//
//  Discussion:
//
//    This routine is given a triangulation, an abstract list of triples
//    of nodes.  It is assumed that the nodes in each triangle are listed
//    in a counterclockwise order, although the routine should work
//    if the nodes are consistently listed in a clockwise order as well.
//
//    It is assumed that each edge of the triangulation is either
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes - as long
//    as the boundary of the hole comprises more than 3 edges!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Output, integer TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT, the number
//    of boundary edges.
//
{
  int boundary_edge_num;
  int e1;
  int e2;
  int *edge;
  int i;
  int interior_edge_num;
  int j;
  int m;
  int n;
  int unique_num;

  m = 2;
  n = 3 * triangle_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < triangle_num; j++ )
  {
    edge[0+(j               )*m] = triangle_node[0+j*3];
    edge[1+(j               )*m] = triangle_node[1+j*3];
    edge[0+(j+  triangle_num)*m] = triangle_node[1+j*3];
    edge[1+(j+  triangle_num)*m] = triangle_node[2+j*3];
    edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*3];
    edge[1+(j+2*triangle_num)*m] = triangle_node[0+j*3];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+j*m], edge[1+j*m] );
    e2 = i4_max ( edge[0+j*m], edge[1+j*m] );
    edge[0+j*m] = e1;
    edge[1+j*m] = e2;
  }
//
//  Ascending sort the column array.
//
  i4col_sort_a ( m, n, edge );
//
//  Get the number of unique columns in EDGE.
//
  unique_num = i4col_sorted_unique_count ( m, n, edge );

  interior_edge_num = 3 * triangle_num - unique_num;

  boundary_edge_num = 3 * triangle_num - 2 * interior_edge_num;

  delete [] edge;

  return boundary_edge_num;
}
//****************************************************************************80

int triangulation_order3_boundary_edge_count_euler ( int node_num,
  int triangle_num, int hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
//
//  Discussion:
//
//    We assume we are given information about a triangulation
//    of a set of nodes in the plane.
//
//    Given the number of nodes and triangles, we are going to apply
//    Euler's formula to determine the number of edges that lie on the
//    boundary of the set of nodes.
//
//    The number of faces, including the infinite face and internal holes,
//    is TRIANGLE_NUM + HOLE_NUM + 1.
//
//    Let BOUNDARY_NUM denote the number of edges on the boundary.
//    Each of the TRIANGLE_NUM triangles uses three edges.  Every edge
//    occurs in two different faces, so the number of edges must be
//    ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2.
//
//    The number of nodes used in the triangulation is NODE_NUM.
//
//    Euler's formula asserts that, for a simple connected figure in the
//    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
//    FACE_NUM faces:
//
//      NODE_NUM - EDGE_NUM + FACE_NUM = 2
//
//    In our context, this becomes
//
//      NODE_NUM - ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2
//      + TRIANGLE_NUM + HOLE_NUM + 1 = 2
//
//    or
//
//      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - TRIANGLE_NUM - 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Marc deBerg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
//    Computational Geometry,
//    Springer, 2000,
//    ISBN: 3-540-65620-0.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int HOLE_NUM, the number of holes.
//
//    Output, int TRIANGULATION_BOUNDARY_COUNT, the number of edges that
//    lie on the convex hull of the triangulation.
//
{
  return ( 2 * node_num + 2 * hole_num - triangle_num - 2 );
}
//****************************************************************************80

bool *triangulation_order3_boundary_node ( int node_num, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_BOUNDARY_NODE indicates nodes on the boundary.
//
//  Discussion:
//
//    This routine is given a triangulation, an abstract list of triples
//    of nodes.  It is assumed that the nodes in each triangle are listed
//    in a counterclockwise order, although the routine should work
//    if the nodes are consistently listed in a clockwise order as well.
//
//    It is assumed that each edge of the triangulation is either
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes - as long
//    as the boundary of the hole comprises more than 3 edges!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Output, bool TRIANGULATION_ORDER3_BOUNDARY_NODE[NODE_NUM],
//    is TRUE if the node is on a boundary edge.
//
{
  int e1;
  int e2;
  int *edge;
  bool equal;
  int i;
  int j;
  int m;
  int n;
  bool *node_boundary;

  m = 2;
  n = 3 * triangle_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < triangle_num; j++ )
  {
    edge[0+(j               )*m] = triangle_node[0+j*3];
    edge[1+(j               )*m] = triangle_node[1+j*3];
    edge[0+(j+  triangle_num)*m] = triangle_node[1+j*3];
    edge[1+(j+  triangle_num)*m] = triangle_node[2+j*3];
    edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*3];
    edge[1+(j+2*triangle_num)*m] = triangle_node[0+j*3];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+j*m], edge[1+j*m] );
    e2 = i4_max ( edge[0+j*m], edge[1+j*m] );
    edge[0+j*m] = e1;
    edge[1+j*m] = e2;
  }
//
//  Ascending sort the column array.
//
  i4col_sort_a ( m, n, edge );
//
//  Records which appear twice are internal edges and can be ignored.
//
  node_boundary = new bool[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    node_boundary[i] = false;
  }

  j = 0;

  while ( j < 3 * triangle_num )
  {
    j = j + 1;

    if ( j == 3 * triangle_num )
    {
      for ( i = 0; i < m; i++ )
      {
        node_boundary[edge[i+(j-1)*m]-1] = true;
      }
      break;
    }

    equal = true;

    for ( i = 0; i < m; i++ )
    {
      if ( edge[i+(j-1)*m] != edge[i+j*m] )
      {
        equal = false;
      }
    }

    if ( equal )
    {
      j = j + 1;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        node_boundary[edge[i+(j-1)*m]-1] = true;
      }
    }

  }

  delete [] edge;

  return node_boundary;
}
//****************************************************************************80

int triangulation_order3_check ( int node_num, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_CHECK makes some simple checks on a triangulation.
//
//  Discussion:
//
//    Because this routine does not receive the physical coordinates of
//    the nodes, it cannot ensure that the triangulation is maximal,
//    that is, that no more triangles can be created.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Output, int TRIANGULATION_CHECK, error flag.
//    0, no error occurred.
//    nonzero, an error occurred, the triangulation is not valid.
//
{
  int boundary_num;
  int error;
  int euler;
  int i;
  int j;
  int *used;
//
//  Checks 1 and 2:
//  node_num must be at least 3.
//  TRIANGLE_NUM must be at least 1.
//
  if ( node_num < 3 )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
    cout << "  The number of nodes is less than 3!\n";
    return 1;
  }

  if ( triangle_num < 1 )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
    cout << "  The number of triangles is less than 1!\n";
    return 2;
  }
//
//  Checks 3 and 4:
//  Verify that all node values are greater than or equal to 1
//  and less than or equal to node_num.
//
  for ( j = 0; j < triangle_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      if ( triangle_node[i+j*3] < 1 )
      {
        cout << "\n";
        cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
        cout << "  Some vertices are less than 1!\n";
        return 3;
      }
    }
  }

  for ( j = 0; j < triangle_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      if ( node_num < triangle_node[i+j*3] )
      {
        cout << "\n";
        cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
        cout << "  Some vertices are greater than node_num!\n";
        return 4;
      }
    }
  }
//
//  Check 5:
//  Verify that every node is used at least once.
//
  used = new int[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    used[i] = 0;
  }

  for ( j = 0; j < triangle_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      used[triangle_node[i+j*3]-1] = used[triangle_node[i+j*3]-1] + 1;
    }
  }

  for ( i = 0; i < node_num; i++ )
  {
    if ( used[i] == 0 )
    {
      cout << "\n";
      cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
      cout << "  Some nodes are never used as triangle vertices!\n";
      cout << "  First example is node " << i+1 << "\n";
      delete [] used;
      return 5;
    }
  }
  delete [] used;
//
//  Check 6:
//  Verify that no node is repeated in a triangle.
//
  for ( j = 0; j < triangle_num; j++ )
  {
    if ( triangle_node[0+j*3] == triangle_node[1+j*3] ||
         triangle_node[1+j*3] == triangle_node[2+j*3] ||
         triangle_node[2+j*3] == triangle_node[0+j*3] )
    {
      cout << "\n";
      cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
      cout << "  A triangle contains a null edge!\n";
      return 6;
    }
  }
//
//  Check 7:
//  Verify that no edge is repeated, and that repeated edges occur in
//  negated pairs.
//
  boundary_num = triangulation_order3_edge_check ( triangle_num,
    triangle_node );

  if ( boundary_num < 0 )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
    cout << "  Some edges are repeated or given in the wrong direction!\n";
    return 7;
  }
//
//  Check 8:
//  Does the triangulation satisfy Euler's criterion?
//  If not, then the triangulation is not proper.  (For instance, there
//  might be a hole in the interior.)
//
  euler = boundary_num + triangle_num + 2 - 2 * node_num;

  if ( euler != 0 )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER3_CHECK - Fatal error!\n";
    cout << "  The triangulation does not satisfy Euler's criterion!\n";
    return 8;
  }

  return 0;
}
//****************************************************************************80

int triangulation_order3_edge_check ( int triangle_num, int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EDGE_CHECK checks the edges of a triangulation.
//
//  Discussion:
//
//    Converted from a row-based to a column-based calculation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
//    each triangle.
//
//    Output, int TRIANGULATION_EDGE_CHECK is negative if an error was
//    detected; otherwise, it is the number of edges that lie on the boundary.
//
{
  int boundary_num;
  int i;
  int j;
  int k;
  int *col;
  int tri;
  int triangle_order = 3;
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the three neighbor relations:
//
//    (I,J,+1) or (J,I,-1),
//    (J,K,+1) or (K,J,-1),
//    (K,I,+1) or (I,K,-1)
//
//  where we choose (I,J,+1) if I < J, or else (J,I,-1) and so on.
//
  col = new int[3*(3*triangle_num)];

  for ( tri = 0; tri < triangle_num; tri++ )
  {
    i = triangle_node[0+tri*triangle_order];
    j = triangle_node[1+tri*triangle_order];
    k = triangle_node[2+tri*triangle_order];

    if ( i < j )
    {
      col[0+(3*tri+0)*3] =  i;
      col[1+(3*tri+0)*3] =  j;
      col[2+(3*tri+0)*3] = +1;
    }
    else
    {
      col[0+(3*tri+0)*3] =  j;
      col[1+(3*tri+0)*3] =  i;
      col[2+(3*tri+0)*3] = -1;
    }

    if ( j < k )
    {
      col[0+(3*tri+1)*3] =  j;
      col[1+(3*tri+1)*3] =  k;
      col[2+(3*tri+1)*3] = +1;
    }
    else
    {
      col[0+(3*tri+1)*3] =  k;
      col[1+(3*tri+1)*3] =  j;
      col[2+(3*tri+1)*3] = -1;
    }

    if ( k < i )
    {
      col[0+(3*tri+2)*3] =  k;
      col[1+(3*tri+2)*3] =  i;
      col[2+(3*tri+2)*3] = +1;
    }
    else
    {
      col[0+(3*tri+2)*3] =  i;
      col[1+(3*tri+2)*3] =  k;
      col[2+(3*tri+2)*3] = -1;
    }
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//
  i4col_sort_a ( 3, 3*triangle_num, col );
//
//  Step 3.
//
//  If any record occurs twice, we have an error.
//  Unpaired records lie on the convex hull.
//
  i = 0;
  boundary_num = 0;

  while ( i < 3 * triangle_num )
  {
    i = i + 1;

    if ( i == 3 * triangle_num )
    {
      boundary_num = boundary_num + 1;
    }
    else
    {
      if ( col[0+(i-1)*3] == col[0+i*3] &&
           col[1+(i-1)*3] == col[1+i*3] )
      {
        if ( col[2+(i-1)*3] == col[2+i*3] )
        {
          cout << "\n";
          cout << "TRIANGULATION_ORDER3_EDGE_CHECK - Warning!\n";
          cout << "  An edge occurs twice.\n";
          delete [] col;
          boundary_num = -1;
          return boundary_num;
        }
        else
        {
          i = i + 1;
        }
      }
      else
      {
        boundary_num = boundary_num + 1;
      }
    }
  }

  delete [] col;

  return boundary_num;
}
//****************************************************************************80

void triangulation_order3_example1 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EXAMPLE1 sets up a sample triangulation.
//
//  Discussion:
//
//    This triangulation is actually a Delaunay triangulation.
//
//    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
//    determined by calling TRIANGULATION_ORDER3_EXAMPLE1_SIZE first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  Negative values indicate edges that lie on the exterior.
//
{
# define DIM_NUM 2
# define NODE_NUM 13
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int i;
  static int triangle_neighbor_save[3*TRIANGLE_NUM] = {
       -4,  -13,    2,
        1,    4,    3,
        2,    5,    7,
        2,  -43,    8,
        3,    8,    6,
        5,    9,    7,
        3,    6,   -3,
        5,    4,   10,
        6,   10,   12,
        9,    8,   11,
       12,   10,   14,
        9,   11,   13,
      -23,   12,   16,
       11,  -47,   15,
       16,   14,  -50,
       13,   15,  -39 };
  static int triangle_node_save[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     3,   4,   1,
     3,   1,   2,
     3,   2,   8,
     2,   1,   5,
     8,   2,  13,
     8,  13,   9,
     3,   8,   9,
    13,   2,   5,
     9,  13,   7,
     7,  13,   5,
     6,   7,   5,
     9,   7,   6,
    10,   9,   6,
     6,   5,  12,
    11,   6,  12,
    10,   6,  11 };
  static double node_xy_save[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       2.0, 2.0,
      -1.0, 3.0,
      -2.0, 2.0,
       8.0, 2.0,
       9.0, 5.0,
       7.0, 4.0,
       5.0, 6.0,
       6.0, 7.0,
       8.0, 8.0,
      11.0, 7.0,
      10.0, 4.0,
       6.0, 4.0 };

  for ( i = 0; i < 3 * TRIANGLE_NUM; i++ )
  {
    triangle_neighbor[i] = triangle_neighbor_save[i];
  }

  for ( i = 0; i < TRIANGLE_ORDER * TRIANGLE_NUM; i++ )
  {
    triangle_node[i] = triangle_node_save[i];
  }

  for ( i = 0; i < DIM_NUM * NODE_NUM; i++ )
  {
    node_xy[i] = node_xy_save[i];
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void triangulation_order3_example1_size ( int *node_num, int *triangle_num,
  int *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EXAMPLE1_SIZE sets sizes for a sample triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *TRIANGLE_NUM, the number of triangles.
//
//    Output, int *HOLE_NUM, the number of holes.
//
{
  *node_num = 13;
  *triangle_num = 16;
  *hole_num = 0;

  return;
}
//****************************************************************************80

void triangulation_order3_example2 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EXAMPLE2 sets up a sample triangulation.
//
//  Discussion:
//
//    This triangulation is actually a Delaunay triangulation.
//
//    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
//    determined by calling TRIANGULATION_ORDER3_EXAMPLE2_SIZE first.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
//    triangles.
//
//    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  Negative values indicate edges that lie on the exterior.
//
{
# define DIM_NUM 2
# define NODE_NUM 25
# define TRIANGLE_NUM 32
# define TRIANGLE_ORDER 3

  int i;
  static int triangle_neighbor_save[3*TRIANGLE_NUM] = {
    -1,  2, -1,
     9,  1,  3,
    -1,  4,  2,
    11,  3,  5,
    -1,  6,  4,
    13,  5,  7,
    -1,  8,  6,
    15,  7, -1,
     2, 10, -1,
    17,  9, 11,
     4, 12, 10,
    19, 11, 13,
     6, 14, 12,
    21, 13, 15,
     8, 16, 14,
    23, 15, -1,
    10, 18, -1,
    25, 17, 19,
    12, 20, 18,
    27, 19, 21,
    14, 22, 20,
    29, 21, 23,
    16, 24, 22,
    31, 23, -1,
    18, 26, -1,
    -1, 25, 27,
    20, 28, 26,
    -1, 27, 29,
    22, 30, 28,
    -1, 29, 31,
    24, 32, 30,
    -1, 31, -1 };
  static int triangle_node_save[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  2,  6,
     7,  6,  2,
     2,  3,  7,
     8,  7,  3,
     3,  4,  8,
     9,  8,  4,
     4,  5,  9,
    10,  9,  5,
     6,  7, 11,
    12, 11,  7,
     7,  8, 12,
    13, 12,  8,
     8,  9, 13,
    14, 13,  9,
     9, 10, 14,
    15, 14, 10,
    11, 12, 16,
    17, 16, 12,
    12, 13, 17,
    18, 17, 13,
    13, 14, 18,
    19, 18, 14,
    14, 15, 19,
    20, 19, 15,
    16, 17, 21,
    22, 21, 17,
    17, 18, 22,
    23, 22, 18,
    18, 19, 23,
    24, 23, 19,
    19, 20, 24,
    25, 24, 20 };
  static double node_xy_save[DIM_NUM*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    3.0, 3.0,
    4.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0,
    3.0, 4.0,
    4.0, 4.0  };

  for ( i = 0; i < 3 * TRIANGLE_NUM; i++ )
  {
    triangle_neighbor[i] = triangle_neighbor_save[i];
  }

  for ( i = 0; i < TRIANGLE_ORDER * TRIANGLE_NUM; i++ )
  {
    triangle_node[i] = triangle_node_save[i];
  }

  for ( i = 0; i < DIM_NUM * NODE_NUM; i++ )
  {
    node_xy[i] = node_xy_save[i];
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void triangulation_order3_example2_size ( int *node_num, int *triangle_num,
  int *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_EXAMPLE2_SIZE sets sizes for a sample triangulation.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   16-17-18-19-20
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//   11-12-13-14-15
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    6--7--8--9-10
//    |\ |\ |\ |\ |
//    | \| \| \| \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *TRIANGLE_NUM, the number of triangles.
//
//    Output, int *HOLE_NUM, the number of holes.
//
{
  *node_num = 25;
  *triangle_num = 32;
  *hole_num = 0;

  return;
}
//****************************************************************************80

void triangulation_order3_neighbor ( int triangle_num, int triangle_node[],
  int t1, int s1, int  *t2, int *s2 )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_NEIGHBOR determines a neighbor of a given triangle.
//
//  Discussion:
//
//    A set of nodes is given.  A triangulation of the nodes has been
//    defined and recorded in TRIANGLE_NODE.  The TRIANGLE_NODE data structure
//    records triangles as sets of three nodes, N1, N2, N3, that implicitly
//    define three sides, being the line segments N1-N2, N2-N3, and N3-N1.
//
//    The nodes of the triangle are listed in counterclockwise order.
//    This means that if two triangles share a side, then the nodes
//    defining that side occur in the order (N1,N2) for one triangle,
//    and (N2,N1) for the other.
//
//    The routine is given a triangle and a side, and asked to find
//    another triangle (if any) that shares that side.  The routine
//    simply searches the TRIANGLE_NODE structure for an occurrence of the
//    nodes in the opposite order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that define
//    each triangle.
//
//    Input, int T1, the index of the triangle.
//
//    Input, int S1, the index of the triangle side.
//
//    Output, int *T2, the index of the triangle which is the neighbor
//    to T1 on side S1, or -1 if there is no such neighbor.
//
//    Output, int *S2, the index of the side of triangle T2 which
//    is shared with triangle T1, or -1 if there is no such neighbor.
//
{
  int n1;
  int n2;
  int s;
  int ss;
  int t;

  n1 = triangle_node[s1-1+(t1-1)*3];
  ss = i4_wrap ( s1+1, 1, 3 );
  n2 = triangle_node[ss-1+(t1-1)*3];

  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      if ( triangle_node[s+t*3] == n1 )
      {
        ss = i4_wrap ( s-1, 0, 2 );
        if ( triangle_node[ss+t*3] == n2 )
        {
          *t2 = t + 1;
          *s2 = ss + 1;
          return;
        }
      }
    }
  }

  *t2 = -1;
  *s2 = -1;

  return;
}
//****************************************************************************80

void triangulation_order3_neighbor_nodes ( int node_num, int triangle_num,
  int triangle_node[], int nabes_first[], int nabes_num[], int nabes_max,
  int *nabes_dim, int nabes[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_NEIGHBOR_NODES determines node neighbors.
//
//  Example:
//
//    On input, the triangle data structure is:
//
//    Triangle  Nodes
//    --------  ----------
//     1        3,   4,   1
//     2        3,   1,   2
//     3        3,   2,   6
//     4        2,   1,   5
//     5        6,   2,   5
//
//  On output, the auxilliary neighbor arrays are:
//
//    Node  Num  First
//    ----  ---  -----
//     1     4     1
//     2     4     5
//     3     4     9
//     4     2    13
//     5     3    15
//     6     3    18
//
//  and the neighbor array is:
//
//    Position  Node
//    --------  ----
//
//     1        2
//     2        3
//     3        4
//     4        5
//    -----------
//     5        1
//     6        3
//     7        5
//     8        6
//    -----------
//     9        1
//    10        2
//    11        4
//    12        6
//    -----------
//    13        1
//    14        3
//    -----------
//    15        1
//    16        2
//    17        6
//    -----------
//    18        2
//    19        3
//    20        5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
//    each triangle.
//
//    Output, int NABES_FIRST[NODE_NUM], the index in NABES of the first
//    neighbor in the list for each node.
//
//    Output, int NABES_NUM[NODE_NUM], the number of neighbors of each node.
//
//    Input, int NABES_MAX, the maximum dimension of NABES.
//
//    Output, int *NABES_DIM, the dimension of NABES.
//
//    Output, int NABES[*NABES_DIM], a list of the neighbors of all the nodes.
//    Neighbors of node 1 are listed first, and so on.
//
{
  int i;
  int i_current;
  int j;
  int k;
  int n;
  int nabe;
  int *nabes1;
  int tri;

  nabes = new int[nabes_max];
//
//  Step 1.  From the triangle list (I,J,K)
//  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
//
  n = 0;

  for ( tri = 0; tri < triangle_num; tri++ )
  {
    i = triangle_node[0+tri*3];
    j = triangle_node[1+tri*3];
    k = triangle_node[2+tri*3];
    nabes1[n]   = i;
    nabes1[n+1] = i;
    nabes1[n+2] = j;
    nabes1[n+3] = j;
    nabes1[n+4] = k;
    nabes1[n+5] = k;
    nabes[n]    = j;
    nabes[n+1]  = k;
    nabes[n+2]  = i;
    nabes[n+3]  = k;
    nabes[n+4]  = i;
    nabes[n+5]  = j;

    n = n + 6;
  }
//
//  Step 2. Dictionary sort the neighbor relations.
//
  i4vec2_sort_a ( n, nabes1, nabes );
//
//  Step 3. Remove duplicate entries.
//
  n = i4vec2_sorted_unique ( n, nabes1, nabes );
//
//  Step 4. Construct the NABES_NUM and NABES_FIRST data.
//
  for ( i = 0; i < node_num; i++ )
  {
    nabes_num[i] = 0;
  }
  for ( i = 0; i < node_num; i++ )
  {
    nabes_first[i] = 0;
  }

  i_current = 0;

  for ( nabe = 1; nabe <= n; nabe++ )
  {
    i = nabes1[nabe-1];
    if ( i == i_current )
    {
      nabes_num[i-1] = nabes_num[i-1] + 1;
    }
    else
    {
      i_current = i;
      nabes_first[i-1] = nabe;
      nabes_num[i-1] = 1;
    }
  }

  *nabes_dim = n;

  delete [] nabes1;

  return;
}
//****************************************************************************80

void triangulation_order3_neighbor_nodes_print ( int node_num,
  int nabes_first[], int nabes_num[], int nabes_dim, int nabes[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_NEIGHBOR_NODES_PRINT prints a node neighbor array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 July 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NABES_FIRST[NODE_NUM], the index in NABES of the first
//    neighbor in the list for each node.
//
//    Input, int NABES_NUM[NODE_NUM], the number of neighbors of each node.
//
//    Input, int NABES_DIM, the dimension of NABES.
//
//    Input, int NABES[NABES_DIM], a list of the neighbors of all the nodes.
//    Neighbors of node 1 are listed first, and so on.
//
{
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "  Node Nabes Index  List\n";
  cout << "\n";

  for ( i = 0; i < node_num; i++ )
  {
    cout << setw(4) << i              << "  "
         << setw(4) << nabes_num[i]   << "  "
         << setw(4) << nabes_first[i] << "  ";

    k = 0;
    for ( j = nabes_first[i] - 1; j < nabes_first[i] + nabes_num[i]; j++ )
    {
      if ( k == 10 )
      {
        cout << "\n";
        cout << "                  ";
        k = 0;
      }
      cout << setw(4) << nabes[j] << "  ";
      k = k + 1;
    }
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void triangulation_order3_plot ( string file_name, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_PLOT plots a triangulation of a set of nodes.
//
//  Discussion:
//
//    The triangulation is most usually a Delaunay triangulation,
//    but this is not necessary.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the output file.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists, for each triangle,
//    the indices of the nodes that form the vertices of the triangle.
//
//    Input, int NODE_SHOW:
//    0, do not show nodes;
//    1, show nodes;
//    2, show nodes and label them.
//
//    Input, int TRIANGLE_SHOW:
//    0, do not show triangles;
//    1, show triangles;
//    2, show triangles and label them.
//
{
  double ave_x;
  double ave_y;
  int circle_size;
  int delta;
  int e;
  ofstream file_unit;
  int i;
  int node;
  int triangle;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( file_name.c_str ( ) );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER3_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: triangulation_order3_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";

  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "%  Increase line width from default 0.\n";
  file_unit << "%\n";
  file_unit << "2 setlinewidth\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Draw filled dots at each node:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750  setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps << "  "
        << y_ps << "  "
        << circle_size << " 0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to darker blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850  setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";

    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps     << "  "
        << y_ps + 5 << "  moveto ("
        << node+1   << ") show\n";
    }
  }
//
//  Draw the triangles.
//
  if ( 1 <= triangle_show )
  {
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to red.\n";
    file_unit << "%\n";
    file_unit << "0.900  0.200  0.100 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Draw the triangles.\n";
    file_unit << "%\n";

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      file_unit << "newpath\n";

      for ( i = 1; i <= 4; i++ )
      {
        e = i4_wrap ( i, 1, 3 );

        node = triangle_node[e-1+triangle*3] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
          / ( y_max                     - y_min ) );

        if ( i == 1 )
        {
          file_unit << x_ps << "  " << y_ps << "  moveto\n";
        }
        else
        {
          file_unit << x_ps << "  " << y_ps << "  lineto\n";
        }
      }
      file_unit << "stroke\n";
    }
  }
//
//  Label the triangles.
//
  if ( 2 <= triangle_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the triangles.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker red.\n";
    file_unit << "%\n";
    file_unit << "0.950  0.250  0.150 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 1; i <= 3; i++ )
      {
        node = triangle_node[i-1+triangle*3] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }
      ave_x = ave_x / 3.0;
      ave_y = ave_y / 3.0;

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max         - y_min ) );

      file_unit << x_ps << "  "
                << y_ps << "  moveto ("
                << triangle+1 << ") show\n";
    }
  }

  file_unit << "%\n";
  file_unit << "restore  showpage\n";
  file_unit << "%\n";
  file_unit << "%  End of page.\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
//****************************************************************************80

void triangulation_order3_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_PRINT prints information defining a triangulation.
//
//  Discussion:
//
//    Triangulations created by R8TRIS2 include extra information encoded
//    in the negative values of TRIANGLE_NEIGHBOR.
//
//    Because some of the nodes counted in NODE_NUM may not actually be
//    used in the triangulation, I needed to compute the true number
//    of vertices.  I added this calculation on 13 October 2001.
//
//    Ernest Fasse pointed out an error in the indexing of VERTEX_LIST,
//    which was corrected on 19 February 2004.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  If there is no triangle neighbor on a particular side,
//    the value of TRIANGLE_NEIGHBOR should be negative.  If the
//    triangulation data was created by R8TRIS2, then there is more
//    information encoded in the negative values.
//
{
# define DIM_NUM 2

  int boundary_num;
  int i;
  int j;
  int k;
  int n1;
  int n2;
  int s;
  int s1;
  int s2;
  bool skip;
  int t;
  int *vertex_list;
  int vertex_num;

  cout << "\n";
  cout << "TRIANGULATION_ORDER3_PRINT\n";
  cout << "  Information defining a triangulation.\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num << "\n";

  r8mat_transpose_print ( DIM_NUM, node_num, node_xy, "  Node coordinates" );

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";
  cout << "\n";
  cout << "  Sets of three nodes are used as vertices of\n";
  cout << "  the triangles.  For each triangle, the nodes\n";
  cout << "  are listed in counterclockwise order.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_node, "  Triangle nodes" );

  cout << "\n";
  cout << "  On each side of a given triangle, there is either\n";
  cout << "  another triangle, or a piece of the convex hull.\n";
  cout << "  For each triangle, we list the indices of the three\n";
  cout << "  neighbors, or (if negative) the codes of the\n";
  cout << "  segments of the convex hull.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_neighbor,
    "  Triangle neighbors" );
//
//  Determine VERTEX_NUM, the number of vertices.
//
  vertex_list = new int[3*triangle_num];

  k = 0;
  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      vertex_list[k] = triangle_node[s+t*3];
      k = k + 1;
    }
  }

  i4vec_sort_heap_a ( 3*triangle_num, vertex_list );

  vertex_num = i4vec_sorted_unique ( 3*triangle_num, vertex_list );

  delete [] vertex_list;
//
//  Determine the number of boundary points.
//
  boundary_num = 2 * vertex_num - triangle_num - 2;

  cout << "\n";
  cout << "  The number of boundary points is " << boundary_num << "\n";
  cout << "\n";
  cout << "  The segments that make up the convex hull can be\n";
  cout << "  determined from the negative entries of the triangle\n";
  cout << "  neighbor list.\n";
  cout << "\n";
  cout << "     #   Tri  Side    N1    N2\n";
  cout << "\n";

  skip = false;

  k = 0;

  for ( i = 0; i < triangle_num; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      if ( triangle_neighbor[j+i*3] < 0 )
      {
        s = -triangle_neighbor[j+i*3];
        t = s / 3;

        if ( t < 1 || triangle_num < t )
        {
          cout << "\n";
          cout << "  Sorry, this data does not use the R8TRIS2\n";
          cout << "  convention for convex hull segments.\n";
          skip = true;
          break;
        }

        s1 = ( s % 3 ) + 1;
        s2 = i4_wrap ( s1+1, 1, 3 );
        k = k + 1;
        n1 = triangle_node[s1-1+(t-1)*3];
        n2 = triangle_node[s2-1+(t-1)*3];
        cout                  << "  "
             << setw(4) << k  << "  "
             << setw(4) << t  << "  "
             << setw(4) << s1 << "  "
             << setw(4) << n1 << "  "
             << setw(4) << n2 << "\n";
      }
    }

    if ( skip )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void triangulation_order3_quad ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  void quad_fun ( int n, double xy_vec[], double f_vec[] ), int quad_num,
  double quad_xy[], double quad_w[], double *quad_value, double *region_area )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_QUAD approximates an integral over a triangulation.
//
//  Discussion:
//
//    The routine will accept triangulations of order higher than 3.
//    However, only the first three nodes (the vertices) of each
//    triangle will be used.  This will still produce correct results
//    for higher order triangulations, as long as the sides of the
//    triangle are straight.
//
//    We assume that the vertices of each triangle are listed first
//    in the description of higher order triangles, and we assume that
//    the vertices are listed in counterclockwise order.
//
//    The approximation of the integral is made using a quadrature rule
//    defined on the unit triangle, and supplied by the user.
//
//    The user also supplies the name of a subroutine, here called "QUAD_FUN",
//    which evaluates the integrand at a set of points.  The form is:
//
//      void quad_fun ( int n, double xy_vec[], double f_vec[] )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes in the triangulation.
//
//    Input, double NODE_XY(2,NODE_NUM), the coordinates of the nodes.
//
//    Input, int TRIANGLE_ORDER, the order of triangles in the triangulation.
//
//    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the nodes making up each triangle.
//
//    Input, void QUAD_FUN ( int N, double XY_VEC[], double F_VEC[] ),
//    the name of the function that evaluates the integrand.
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_XY(2,QUAD_NUM), the abscissas of the
//    quadrature rule, in the unit triangle.
//
//    Input, double QUAD_W(QUAD_NUM), the weights of the
//    quadrature rule.
//
//    Output, double *QUAD_VALUE, the estimate of the integral
//    of F(X,Y) over the region covered by the triangulation.
//
//    Output, double *REGION_AREA, the area of the region.
//
{
  int i;
  int j;
  int quad;
  double *quad_f;
  double *quad2_xy;
  double temp;
  int triangle;
  double triangle_area;
  double triangle_xy[2*3];

  quad_f = new double[quad_num];
  quad2_xy = new double[2*quad_num];

  *quad_value = 0.0;
  *region_area = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      for ( i = 0; i < 2; i++ )
      {
        triangle_xy[i+j*2] = node_xy[i+(triangle_node[j+triangle*3]-1)*2];
      }
    }
    triangle_area = triangle_area_2d ( triangle_xy );

    triangle_order3_reference_to_physical ( triangle_xy,
      quad_num, quad_xy, quad2_xy );

    quad_fun ( quad_num, quad2_xy, quad_f );

    temp = 0.0;
    for ( quad = 0; quad < quad_num; quad++ )
    {
      temp = temp + quad_w[quad] * quad_f[quad];
    }

    *quad_value = *quad_value + triangle_area * temp;

    *region_area = *region_area + triangle_area;
  }

  delete [] quad_f;
  delete [] quad2_xy;

  return;
}
//****************************************************************************80

void triangulation_order3_refine_compute ( int node_num1, int triangle_num1,
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2,
  int edge_data[], double node_xy2[], int triangle_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
//
//  Discussion:
//
//    Given a triangle defined by nodes 1, 2, 3, we need to generate
//    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
//    and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.
//
//          3
//         / \
//        /T3 \
//      13----23
//      / \T4 / \
//     /T1 \ /T2 \
//    1----12-----2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes.
//
//    Input, int TRIANGLE_NUM1, the number of triangles.
//
//    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
//
//    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Input, int NODE_NUM2, the number of nodes in the refined mesh.
//
//    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
//
//    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
//    by TRIANGULATION_ORDER3_REFINE_SIZE.
//
//    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
//
//    Output, int TRIANGLE_NODE2[3*TRIANGLE_NUM2], the nodes that make up the
//    triangles in the refined mesh.
//
{
  int edge;
  int i;
  int j;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int node;
  int triangle1;
  int v1;
  int v2;
//
//  Copy the old nodes.
//
  for ( j = 0; j < node_num1; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      node_xy2[i+j*2] = node_xy1[i+j*2];
    }
  }
  for ( j = 0; j < triangle_num2; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_node2[i+j*3] = -1;
    }
  }
//
//  We can assign the existing nodes to the new triangles.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    triangle_node2[0+(triangle1*4+0)*3] = triangle_node1[0+triangle1*3];
    triangle_node2[1+(triangle1*4+1)*3] = triangle_node1[1+triangle1*3];
    triangle_node2[2+(triangle1*4+2)*3] = triangle_node1[2+triangle1*3];
  }

  node = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5] - 1;
    n2 = edge_data[1+edge*5] - 1;
//
//  If this edge is new, create the coordinates and index for this node.
//
    if ( n1 != n1_old || n2 != n2_old )
    {

      if ( node_num2 < node )
      {
        cout << "\n";
        cout << "TRIANGLE_MESH_ORDER3_REFINE - Fatal error!\n";
        cout << "  Node index exceeds NODE_NUM2.\n";
        exit ( 1 );
      }

      for ( i = 0; i < 2; i++ )
      {
        node_xy2[i+node*2] = ( node_xy2[i+n1*2] + node_xy2[i+n2*2] ) / 2.0;
      }

      node = node + 1;

      n1_old = n1;
      n2_old = n2;
    }
//
//  Assign the node to triangles.
//
    v1 = edge_data[2+edge*5];
    v2 = edge_data[3+edge*5];
    triangle1 = edge_data[4+edge*5];

    if ( v1 == 1 && v2 == 2 )
    {
      triangle_node2[0+(triangle1*4+1)*3] = node;
      triangle_node2[1+(triangle1*4+0)*3] = node;
      triangle_node2[2+(triangle1*4+3)*3] = node;
    }
    else if ( v1 == 1 && v2 == 3 )
    {
      triangle_node2[0+(triangle1*4+2)*3] = node;
      triangle_node2[1+(triangle1*4+3)*3] = node;
      triangle_node2[2+(triangle1*4+0)*3] = node;
    }
    else if ( v1 == 2 && v2 == 3 )
    {
      triangle_node2[0+(triangle1*4+3)*3] = node;
      triangle_node2[1+(triangle1*4+2)*3] = node;
      triangle_node2[2+(triangle1*4+1)*3] = node;
    }
  }
  return;
}
//****************************************************************************80

void triangulation_order3_refine_size ( int node_num1, int triangle_num1,
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
//
//  Discussion:
//
//    Given a triangle defined by nodes 1, 2, 3, we need to generate
//    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
//    and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.
//
//          3
//         / \
//        /T3 \
//      13----23
//      / \T4 / \
//     /T1 \ /T2 \
//    1----12-----2
//
//    This routine simply determines the sizes of the resulting node
//    and triangle arrays.
//
//    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
//    data items, one item for every edge of every triangle.  Each
//    data item records, for a given edge, the global indices
//    of the two endpoints, the local indices of the two endpoints,
//    and the index of the triangle.
//
//    Through careful sorting, it is possible to arrange this data in
//    a way that allows the proper generation of the interpolated nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the original mesh.
//
//    Input, int  TRIANGLE_NUM1, the number of triangles in the
//    original mesh.
//
//    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the indices of the nodes
//    that form the triangles in the input mesh.
//
//    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
//
//    Output, int *TRIANGLE_NUM2, the number of triangles in the
//    refined mesh.
//
//    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
//    be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
//
{
  int a;
  int b;
  int edge;
  int i;
  int j;
  int k;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int triangle;
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the edge relations:
//
//    (I,J,1,2,T)
//    (I,K,1,3,T)
//    (J,K,2,3,T)
//
//  In order to make matching easier, we reorder each pair of nodes
//  into ascending order.
//
  for ( triangle = 0; triangle < triangle_num1; triangle++ )
  {
    i = triangle_node1[0+triangle*3];
    j = triangle_node1[1+triangle*3];
    k = triangle_node1[2+triangle*3];

    a = i4_min ( i, j );
    b = i4_max ( i, j );

    edge_data[0+5*(3*triangle+0)] = a;
    edge_data[1+5*(3*triangle+0)] = b;
    edge_data[2+5*(3*triangle+0)] = 1;
    edge_data[3+5*(3*triangle+0)] = 2;
    edge_data[4+5*(3*triangle+0)] = triangle;

    a = i4_min ( i, k );
    b = i4_max ( i, k );

    edge_data[0+5*(3*triangle+1)] = a;
    edge_data[1+5*(3*triangle+1)] = b;
    edge_data[2+5*(3*triangle+1)] = 1;
    edge_data[3+5*(3*triangle+1)] = 3;
    edge_data[4+5*(3*triangle+1)] = triangle;

    a = i4_min ( j, k );
    b = i4_max ( j, k );

    edge_data[0+5*(3*triangle+2)] = a;
    edge_data[1+5*(3*triangle+2)] = b;
    edge_data[2+5*(3*triangle+2)] = 2;
    edge_data[3+5*(3*triangle+2)] = 3;
    edge_data[4+5*(3*triangle+2)] = triangle;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:2; the routine we call here
//  sorts on the full column but that won't hurt us.
//
//  What we need is to find all cases where triangles share an edge.
//  By sorting the columns of the EDGE_DATA array, we will put shared edges
//  next to each other.
//
  i4col_sort_a ( 5, 3*triangle_num1, edge_data );
//
//  Step 3. All the triangles which share an edge show up as consecutive
//  columns with identical first two entries.  Figure out how many new
//  nodes there are, and allocate space for their coordinates.
//
  *node_num2 = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
    if ( n1 != n1_old || n2 != n2_old )
    {
      *node_num2 = *node_num2 + 1;
      n1_old = n1;
      n2_old = n2;
    }
  }

  *triangle_num2 = 4 * triangle_num1;

  return;
}
//****************************************************************************80

void triangulation_order3_sample ( int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int num_ran, int *seed,
  double xd[], int td[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER3_SAMPLE returns random points in a triangulation.
//
//  Discussion:
//
//    It is assumed that the triangulation consists of a set of non-overlapping
//    triangles.
//
//    The point is chosen uniformly in the area covered by the triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
//    triangles.
//
//    Input, int NUM_RAN, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double XD[2*NUM_RAN], the sample points.
//
//    Output, int TD[NUM_RAN], the triangle to which each sample point
//    belongs.
//
{
  double area;
  double *area_cum;
  double area_total;
  int i;
  int i1;
  int i2;
  int i3;
  int left;
  double r;
  int right;
  double t[2*3];
//
//  Compute the areas of the triangles.
//  Build a cumulative area vector.
//  Convert it to a relative cumulative area vector.
//
  area_cum = new double[triangle_num+1];
  area_cum[0] = 0.0;

  for ( i = 0; i < triangle_num; i++ )
  {
    i1 = triangle_node[0+i*3];
    t[0+0*2] = node_xy[0+i1*2];
    t[1+0*2] = node_xy[1+i1*2];

    i2 = triangle_node[1+i*3];
    t[0+1*2] = node_xy[0+i2*2];
    t[1+1*2] = node_xy[1+i2*2];

    i3 = triangle_node[2+i*3];
    t[0+2*2] = node_xy[0+i3*2];
    t[1+2*2] = node_xy[1+i3*2];

    area_cum[i+1] = area_cum[i] + triangle_area_2d ( t );
  }

  area_total = area_cum[triangle_num];

  for ( i = 0; i <= triangle_num; i++ )
  {
    area_cum[i] = area_cum[i] / area_total;
  }
//
//  Pick random values.  A random value R indicates the corresponding triangle
//  whose cumulative relative area contains R.
//
//  Bracket the random value in the cumulative relative areas,
//  indicating a triangle.
//
//  Pick a random point in the triangle.
//
  for ( i = 0; i < num_ran; i++ )
  {
    r = r8_uniform_01 ( seed );

    r8vec_bracket ( triangle_num+1, area_cum, r, &left, &right );

    td[i] = right - 1;

    i1 = triangle_node[0+(td[i]-1)*3];
    t[0+0*2] = node_xy[0+i1*2];
    t[1+0*2] = node_xy[1+i1*2];

    i2 = triangle_node[1+(td[i]-1)*3];
    t[0+1*2] = node_xy[0+i2*2];
    t[1+1*2] = node_xy[1+i2*2];

    i3 = triangle_node[2+(td[i]-1)*3];
    t[0+2*2] = node_xy[0+i3*2];
    t[1+2*2] = node_xy[1+i3*2];

    triangle_sample ( t, 1, seed, xd+i*2 );
  }

  delete [] area_cum;

  return;
}
//****************************************************************************80

void triangulation_order4_plot ( string plot_filename, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER4_PLOT plots a 4-node triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string PLOT_FILENAME, the name of the output file.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[4*TRIANGLE_NUM], lists, for each triangle,
//    the indices of the nodes that form the vertices of the triangle,
//    and the centroid.
//
//    Input, int NODE_SHOW:
//    0, do not show nodes;
//    1, show nodes;
//    2, show nodes and label them.
//
//    Input, int TRIANGLE_SHOW:
//    0, do not show triangles;
//    1, show triangles;
//    2, show triangles and label them.
//
{
  double ave_x;
  double ave_y;
  int circle_size;
  int delta;
  int e;
  ofstream plot_unit;
  int i;
  int node;
  int triangle;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  plot_unit.open ( plot_filename.c_str ( ) );

  if ( !plot_unit )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER4_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  plot_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  plot_unit << "%%Creator: triangulation_order4_plot.C\n";
  plot_unit << "%%Title: " << plot_filename << "\n";

  plot_unit << "%%Pages: 1\n";
  plot_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  plot_unit << "%%Document-Fonts: Times-Roman\n";
  plot_unit << "%%LanguageLevel: 1\n";
  plot_unit << "%%EndComments\n";
  plot_unit << "%%BeginProlog\n";
  plot_unit << "/inch {72 mul} def\n";
  plot_unit << "%%EndProlog\n";
  plot_unit << "%%Page:      1     1\n";
  plot_unit << "save\n";
  plot_unit << "%\n";
  plot_unit << "%  Increase line width from default 0.\n";
  plot_unit << "%\n";
  plot_unit << "2 setlinewidth\n";
  plot_unit << "%\n";
  plot_unit << "% Set the RGB line color to very light gray.\n";
  plot_unit << "%\n";
  plot_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  plot_unit << "%\n";
  plot_unit << "% Draw a gray border around the page.\n";
  plot_unit << "%\n";
  plot_unit << "newpath\n";
  plot_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  plot_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  plot_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  plot_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  plot_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  plot_unit << "stroke\n";
  plot_unit << "%\n";
  plot_unit << "% Set RGB line color to black.\n";
  plot_unit << "%\n";
  plot_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  plot_unit << "%\n";
  plot_unit << "%  Set the font and its size:\n";
  plot_unit << "%\n";
  plot_unit << "/Times-Roman findfont\n";
  plot_unit << "0.50 inch scalefont\n";
  plot_unit << "setfont\n";
  plot_unit << "%\n";
  plot_unit << "%  Print a title:\n";
  plot_unit << "%\n";
  plot_unit << "%  210  702 moveto\n";
  plot_unit << "%(Pointset) show\n";
  plot_unit << "%\n";
  plot_unit << "% Define a clipping polygon\n";
  plot_unit << "%\n";
  plot_unit << "newpath\n";
  plot_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  plot_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  plot_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  plot_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  plot_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  plot_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= node_show )
  {
    plot_unit << "%\n";
    plot_unit << "%  Draw filled dots at each node:\n";
    plot_unit << "%\n";
    plot_unit << "%  Set the color to blue:\n";
    plot_unit << "%\n";
    plot_unit << "0.000  0.150  0.750  setrgbcolor\n";
    plot_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      plot_unit << "newpath  "
        << x_ps << "  "
        << y_ps << "  "
        << circle_size << " 0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    plot_unit << "%\n";
    plot_unit << "%  Label the nodes:\n";
    plot_unit << "%\n";
    plot_unit << "%  Set the color to darker blue:\n";
    plot_unit << "%\n";
    plot_unit << "0.000  0.250  0.850  setrgbcolor\n";
    plot_unit << "/Times-Roman findfont\n";
    plot_unit << "0.20 inch scalefont\n";
    plot_unit << "setfont\n";

    plot_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      plot_unit << "newpath  "
        << x_ps     << "  "
        << y_ps + 5 << "  moveto ("
        << node+1   << ") show\n";
    }
  }
//
//  Draw the triangles.
//
  if ( 1 <= triangle_show )
  {
    plot_unit << "%\n";
    plot_unit << "%  Set the RGB color to red.\n";
    plot_unit << "%\n";
    plot_unit << "0.900  0.200  0.100 setrgbcolor\n";
    plot_unit << "%\n";
    plot_unit << "%  Draw the triangles.\n";
    plot_unit << "%\n";

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      plot_unit << "newpath\n";

      for ( i = 1; i <= 4; i++ )
      {
        e = i4_wrap ( i, 1, 3 );

        node = triangle_node[e-1+triangle*4] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
          / ( y_max                     - y_min ) );

        if ( i == 1 )
        {
          plot_unit << x_ps << "  " << y_ps << "  moveto\n";
        }
        else
        {
          plot_unit << x_ps << "  " << y_ps << "  lineto\n";
        }
      }
      plot_unit << "stroke\n";
    }
  }
//
//  Label the triangles.
//
  if ( 2 <= triangle_show )
  {
    plot_unit << "%\n";
    plot_unit << "%  Label the triangles.\n";
    plot_unit << "%\n";
    plot_unit << "%  Set the RGB color to darker red.\n";
    plot_unit << "%\n";
    plot_unit << "0.950  0.250  0.150 setrgbcolor\n";
    plot_unit << "/Times-Roman findfont\n";
    plot_unit << "0.20 inch scalefont\n";
    plot_unit << "setfont\n";
    plot_unit << "%\n";

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 1; i <= 3; i++ )
      {
        node = triangle_node[i-1+triangle*4] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }
      ave_x = ave_x / 3.0;
      ave_y = ave_y / 3.0;

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max         - y_min ) );

      plot_unit << x_ps << "  "
                << y_ps << "  moveto ("
                << triangle+1 << ") show\n";
    }
  }

  plot_unit << "%\n";
  plot_unit << "restore  showpage\n";
  plot_unit << "%\n";
  plot_unit << "%  End of page.\n";
  plot_unit << "%\n";
  plot_unit << "%%Trailer\n";
  plot_unit << "%%EOF\n";

  plot_unit.close ( );

  return;
}
//****************************************************************************80

int triangulation_order6_adj_count ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_col[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 6-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\    |\    |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    |    \|    \|
//   11-12-13-14-15
//    |\    |\    |
//    | \   | \   |
//    6  7  8  9 10
//    |   \ |   \ |
//    |    \|    \|
//    1--2--3--4--5
//
//    A sample grid.
//
//
//    Below, we have a chart that lists the nodes adjacent to each node, with
//    an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).
//
//    N   Adjacencies
//
//    1:  *  2  3  6  7 11
//    2:  1  *  3  6  7 11
//    3:  1  2  *  4  5  6  7  8  9 11 12 13
//    4:  3  *  5  8  9 13
//    5:  3  4  *  8  9 10 13 14 15
//    6:  1  2  3  *  7 11
//    7:  1  2  3  6  *  8 11 12 13
//    8:  3  4  5  7  *  9 11 12 13
//    9:  3  4  5  8  * 10 13 14 15
//   10:  5  9  * 13 14 15
//   11:  1  2  3  6  7  8  * 12 13 16 17 21
//   12:  3  7  8 11  * 13 16 17 21
//   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
//   14:  5  9 10 13  * 15 18 19 23
//   15:  5  9 10 13 14  * 18 19 20 23 24 25
//   16: 11 12 13  * 17 21
//   17: 11 12 13 16  * 18 21 22 23
//   18: 13 14 15 17  * 19 21 22 23
//   19: 13 14 15 18  * 20 23 24 25
//   20: 15 19  * 23 24 25
//   21: 11 12 13 16 17 18  * 22 23
//   22: 13 17 18 21  * 23
//   23: 13 14 15 17 18 19 20 21 22  * 24 25
//   24: 15 19 20 23  * 25
//   25: 15 19 20 23 24  *
//
//    Below, we list the number of adjancencies to lower-indexed nodes, to
//    the node itself, to higher-indexed nodes, the total number of
//    adjacencies for this node, and the location of the first and last
//    entries required to list this set of adjacencies in a single list
//    of all the adjacencies.
//
//    N   Below  Self   Above   Total First  Last
//
//   --      --    --      --      --   ---     0
//    1:      0     1       5       6     1     6
//    2:      1     1       4       6     7    12
//    3:      2     1       9      12    13    24
//    4:      1     1       4       6    25    30
//    5:      2     1       6       9    31    39
//    6:      3     1       2       6    40    45
//    7:      4     1       4       9    46    54
//    8:      4     1       4       9    55    63
//    9:      4     1       4       9    62    72
//   10:      2     1       3       6    73    78
//   11:      6     1       5      12    79    90
//   12:      4     1       4       9    91    99
//   13:      9     1       9      19   100   118
//   14:      4     1       4       9   119   127
//   15:      5     1       6      12   128   139
//   16:      3     1       2       6   140   145
//   17:      4     1       4       9   146   154
//   18:      4     1       4       9   155   163
//   19:      4     1       4       9   164   172
//   20:      2     1       3       6   173   178
//   21:      6     1       2       9   179   187
//   22:      4     1       1       6   188   193
//   23:      9     1       2      12   194   205
//   24:      4     1       1       6   206   211
//   25:      5     1       0       6   212   217
//   --      --    --      --      --   218   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Output, int TRIANGULATION_ORDER6_ADJ_COUNT, the number of adjacencies.
//
//    Output, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
{
  int adj_num;
  int i;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int node;
  int triangle;
  int triangle_order = 6;
  int triangle2;

  adj_num = 0;
//
//  Set every node to be adjacent to itself.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_col[node] = 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
    n4 = triangle_node[3+triangle*triangle_order];
    n5 = triangle_node[4+triangle*triangle_order];
    n6 = triangle_node[5+triangle*triangle_order];
//
//  For sure, we add the adjacencies:
//    43 / (34)
//    51 / (15)
//    54 / (45)
//    62 / (26)
//    64 / (46)
//    65 / (56)
//
    adj_col[n3-1] = adj_col[n3-1] + 1;
    adj_col[n4-1] = adj_col[n4-1] + 1;
    adj_col[n1-1] = adj_col[n1-1] + 1;
    adj_col[n5-1] = adj_col[n5-1] + 1;
    adj_col[n4-1] = adj_col[n4-1] + 1;
    adj_col[n5-1] = adj_col[n5-1] + 1;
    adj_col[n2-1] = adj_col[n2-1] + 1;
    adj_col[n6-1] = adj_col[n6-1] + 1;
    adj_col[n4-1] = adj_col[n4-1] + 1;
    adj_col[n6-1] = adj_col[n6-1] + 1;
    adj_col[n5-1] = adj_col[n5-1] + 1;
    adj_col[n6-1] = adj_col[n6-1] + 1;
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//    21 / 12
//    41 / 14
//    42 / 24
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n4-1] = adj_col[n4-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n4-1] = adj_col[n4-1] + 1;
    }
//
//  Maybe add
//    32 / 23
//    52 / 25
//    53 / 35
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n2-1] = adj_col[n2-1] + 1;
      adj_col[n5-1] = adj_col[n5-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n5-1] = adj_col[n5-1] + 1;
    }
//
//  Maybe add
//    31 / 13
//    61 / 16
//    63 / 36
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n1-1] = adj_col[n1-1] + 1;
      adj_col[n6-1] = adj_col[n6-1] + 1;
      adj_col[n3-1] = adj_col[n3-1] + 1;
      adj_col[n6-1] = adj_col[n6-1] + 1;
    }
  }
//
//  We used ADJ_COL to count the number of entries in each column.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num; 1 <= node; node-- )
  {
    adj_col[node] = adj_col[node-1];
  }
  adj_col[0] = 1;
  for ( i = 1; i <= node_num; i++ )
  {
    adj_col[i]= adj_col[i-1] + adj_col[i];
  }

  adj_num = adj_col[node_num] - 1;

  return adj_num;
}
//****************************************************************************80

int *triangulation_order6_adj_set ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_num, int adj_col[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The triangulation is assumed to involve 6-node triangles.
//
//    Two nodes are "adjacent" if they are both nodes in some triangle.
//    Also, a node is considered to be adjacent to itself.
//
//    This routine can be used to create the compressed column storage
//    for a quadratic triangle finite element discretization of
//    Poisson's equation in two dimensions.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//    The local node numbering
//
//
//   21-22-23-24-25
//    |\    |\    |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    |    \|    \|
//   11-12-13-14-15
//    |\    |\    |
//    | \   | \   |
//    6  7  8  9 10
//    |   \ |   \ |
//    |    \|    \|
//    1--2--3--4--5
//
//    A sample grid.
//
//
//    Below, we have a chart that lists the nodes adjacent to each node, with
//    an asterisk to indicate the adjacency of the node to itself
//    (in some cases, you want to count this self adjacency and in some
//    you don't).
//
//    N   Adjacencies
//
//    1:  *  2  3  6  7 11
//    2:  1  *  3  6  7 11
//    3:  1  2  *  4  5  6  7  8  9 11 12 13
//    4:  3  *  5  8  9 13
//    5:  3  4  *  8  9 10 13 14 15
//    6:  1  2  3  *  7 11
//    7:  1  2  3  6  *  8 11 12 13
//    8:  3  4  5  7  *  9 11 12 13
//    9:  3  4  5  8  * 10 13 14 15
//   10:  5  9  * 13 14 15
//   11:  1  2  3  6  7  8  * 12 13 16 17 21
//   12:  3  7  8 11  * 13 16 17 21
//   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
//   14:  5  9 10 13  * 15 18 19 23
//   15:  5  9 10 13 14  * 18 19 20 23 24 25
//   16: 11 12 13  * 17 21
//   17: 11 12 13 16  * 18 21 22 23
//   18: 13 14 15 17  * 19 21 22 23
//   19: 13 14 15 18  * 20 23 24 25
//   20: 15 19  * 23 24 25
//   21: 11 12 13 16 17 18  * 22 23
//   22: 13 17 18 21  * 23
//   23: 13 14 15 17 18 19 20 21 22  * 24 25
//   24: 15 19 20 23  * 25
//   25: 15 19 20 23 24  *
//
//    Below, we list the number of adjancencies to lower-indexed nodes, to
//    the node itself, to higher-indexed nodes, the total number of
//    adjacencies for this node, and the location of the first and last
//    entries required to list this set of adjacencies in a single list
//    of all the adjacencies.
//
//    N   Below  Self   Above   Total First  Last
//
//   --      --    --      --      --   ---     0
//    1:      0     1       5       6     1     6
//    2:      1     1       4       6     7    12
//    3:      2     1       9      12    13    24
//    4:      1     1       4       6    25    30
//    5:      2     1       6       9    31    39
//    6:      3     1       2       6    40    45
//    7:      4     1       4       9    46    54
//    8:      4     1       4       9    55    63
//    9:      4     1       4       9    62    72
//   10:      2     1       3       6    73    78
//   11:      6     1       5      12    79    90
//   12:      4     1       4       9    91    99
//   13:      9     1       9      19   100   118
//   14:      4     1       4       9   119   127
//   15:      5     1       6      12   128   139
//   16:      3     1       2       6   140   145
//   17:      4     1       4       9   146   154
//   18:      4     1       4       9   155   163
//   19:      4     1       4       9   164   172
//   20:      2     1       3       6   173   178
//   21:      6     1       2       9   179   187
//   22:      4     1       1       6   188   193
//   23:      9     1       2      12   194   205
//   24:      4     1       1       6   206   211
//   25:      5     1       0       6   212   217
//   --      --    --      --      --   218   ---
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
//    a triangle, lists the neighboring triangle, or -1 if there is
//    no neighbor.
//
//    Input, int ADJ_NUM, the number of adjacencies.
//
//    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
//
//    Output, int TRIANGULATION_ORDER6_ADJ_SET[ADJ_NUM], the adjacency
//    information.
//
{
  int *adj;
  int *adj_copy;
  int k;
  int k1;
  int k2;
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int node;
  int triangle;
  int triangle2;
  int triangle_order = 6;

  adj = new int[adj_num];
  for ( k = 0; k < adj_num; k++ )
  {
    adj[k] = -1;
  }

  adj_copy = new int[node_num];
  for ( node = 0; node < node_num; node++ )
  {
    adj_copy[node] = adj_col[node];
  }
//
//  Set every node to be adjacent to itself.
//
  for ( node = 1; node <= node_num; node++ )
  {
    adj[adj_copy[node-1]-1] = node;
    adj_copy[node-1] = adj_copy[node-1] + 1;
  }
//
//  Examine each triangle.
//
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    n1 = triangle_node[0+triangle*triangle_order];
    n2 = triangle_node[1+triangle*triangle_order];
    n3 = triangle_node[2+triangle*triangle_order];
    n4 = triangle_node[3+triangle*triangle_order];
    n5 = triangle_node[4+triangle*triangle_order];
    n6 = triangle_node[5+triangle*triangle_order];
//
//  For sure, we add the adjacencies:
//    43 / (34)
//    51 / (15)
//    54 / (45)
//    62 / (26)
//    64 / (46)
//    65 / (56)
//
    adj[adj_copy[n3-1]-1] = n4;
    adj_copy[n3-1] = adj_copy[n3-1] + 1;
    adj[adj_copy[n4-1]-1] = n3;
    adj_copy[n4-1] = adj_copy[n4-1] + 1;

    adj[adj_copy[n1-1]-1] = n5;
    adj_copy[n1-1] = adj_copy[n1-1] + 1;
    adj[adj_copy[n5-1]-1] = n1;
    adj_copy[n5-1] = adj_copy[n5-1] + 1;

    adj[adj_copy[n4-1]-1] = n5;
    adj_copy[n4-1] = adj_copy[n4-1] + 1;
    adj[adj_copy[n5-1]-1] = n4;
    adj_copy[n5-1] = adj_copy[n5-1] + 1;

    adj[adj_copy[n2-1]-1] = n6;
    adj_copy[n2-1] = adj_copy[n2-1] + 1;
    adj[adj_copy[n6-1]-1] = n2;
    adj_copy[n6-1] = adj_copy[n6-1] + 1;

    adj[adj_copy[n4-1]-1] = n6;
    adj_copy[n4-1] = adj_copy[n4-1] + 1;
    adj[adj_copy[n6-1]-1] = n4;
    adj_copy[n6-1] = adj_copy[n6-1] + 1;

    adj[adj_copy[n5-1]-1] = n6;
    adj_copy[n5-1] = adj_copy[n5-1] + 1;
    adj[adj_copy[n6-1]-1] = n5;
    adj_copy[n6-1] = adj_copy[n6-1] + 1;
//
//  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
//  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
//  or if this triangle is the first of the pair in which the edge
//  occurs (TRIANGLE < TRIANGLE2).
//
//  Maybe add
//    21 / 12
//    41 / 14
//    42 / 24
//
    triangle2 = triangle_neighbor[0+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n2;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n2-1]-1] = n1;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n1-1]-1] = n4;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n4-1]-1] = n1;
      adj_copy[n4-1] = adj_copy[n4-1] + 1;
      adj[adj_copy[n2-1]-1] = n4;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n4-1]-1] = n2;
      adj_copy[n4-1] = adj_copy[n4-1] + 1;
    }
//
//  Maybe add
//    32 / 23
//    52 / 25
//    53 / 35
//
    triangle2 = triangle_neighbor[1+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n2-1]-1] = n3;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n3-1]-1] = n2;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n2-1]-1] = n5;
      adj_copy[n2-1] = adj_copy[n2-1] + 1;
      adj[adj_copy[n5-1]-1] = n2;
      adj_copy[n5-1] = adj_copy[n5-1] + 1;
      adj[adj_copy[n3-1]-1] = n5;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n5-1]-1] = n3;
      adj_copy[n5-1] = adj_copy[n5-1] + 1;
    }
//
//  Maybe add
//    31 / 13
//    61 / 16
//    63 / 36
//
    triangle2 = triangle_neighbor[2+triangle*3];

    if ( triangle2 < 0 || triangle < triangle2 )
    {
      adj[adj_copy[n1-1]-1] = n3;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n3-1]-1] = n1;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n1-1]-1] = n6;
      adj_copy[n1-1] = adj_copy[n1-1] + 1;
      adj[adj_copy[n6-1]-1] = n1;
      adj_copy[n6-1] = adj_copy[n6-1] + 1;
      adj[adj_copy[n3-1]-1] = n6;
      adj_copy[n3-1] = adj_copy[n3-1] + 1;
      adj[adj_copy[n6-1]-1] = n3;
      adj_copy[n6-1] = adj_copy[n6-1] + 1;
    }
  }
//
//  Ascending sort the entries for each node.
//
  for ( node = 1; node <= node_num; node++ )
  {
    k1 = adj_col[node-1];
    k2 = adj_col[node]-1;
    i4vec_sort_heap_a ( k2+1-k1, adj+k1-1 );
  }

  delete [] adj_copy;

  return adj;
}
//****************************************************************************80

int triangulation_order6_boundary_edge_count ( int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the boundary edges.
//
//  Discussion:
//
//    This routine is given a triangulation, a set of 6-node triangles.
//    It is assumed that, in each list of 6 nodes, the vertices are listed
//    first, in counterclockwise order, followed by the three midside nodes,
//    in counterclockwise order, starting with the node between vertices
//    1 and 2.
//
//    It is assumed that each edge of the triangulation is either
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes - as long
//    as the boundary of the hole comprises more than 3 edges!
//
//    Except for the dimension of TRIANGLE, this routine is identical
//    to the routine for the order 3 case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
//    triangles.  These should be listed in counterclockwise order.
//
//    Output, integer TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT, the number
//    of boundary edges.
//
{
  int boundary_edge_num;
  int e1;
  int e2;
  int *edge;
  int i;
  int interior_edge_num;
  int j;
  int m;
  int n;
  int unique_num;

  m = 2;
  n = 3 * triangle_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < triangle_num; j++ )
  {
    edge[0+(j               )*m] = triangle_node[0+j*6];
    edge[1+(j               )*m] = triangle_node[1+j*6];
    edge[0+(j+  triangle_num)*m] = triangle_node[1+j*6];
    edge[1+(j+  triangle_num)*m] = triangle_node[2+j*6];
    edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*6];
    edge[1+(j+2*triangle_num)*m] = triangle_node[0+j*6];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+j*m], edge[1+j*m] );
    e2 = i4_max ( edge[0+j*m], edge[1+j*m] );
    edge[0+j*m] = e1;
    edge[1+j*m] = e2;
  }
//
//  Ascending sort the column array.
//
  i4col_sort_a ( m, n, edge );
//
//  Get the number of unique columns in EDGE.
//
  unique_num = i4col_sorted_unique_count ( m, n, edge );

  interior_edge_num = 3 * triangle_num - unique_num;

  boundary_edge_num = 3 * triangle_num - 2 * interior_edge_num;

  delete [] edge;

  return boundary_edge_num;
}
//****************************************************************************80

int triangulation_order6_boundary_edge_count_euler ( int node_num,
  int triangle_num, int hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
//
//  Discussion:
//
//    We assume we are given information about an order 6 triangulation
//    of a set of nodes in the plane.
//
//    By ignoring the midside nodes, we can determine the corresponding
//    information for an order 3 triangulation, and apply
//    Euler's formula to determine the number of edges that lie on the
//    boundary of the set of nodes.
//
//    Thus, if we have TRIANGLE_NUM triangles, and NODE_NUM nodes, we
//    imagine that each triangle is replaced by 4 triangles, created
//    by adding the edges created by joining the midside nodes.
//
//    Thus, for 4 * TRIANGLE_NUM triangles, we can apply Euler's formula
//    to compute the number of boundary edges.
//
//    Now, to adjust the data to our order 6 triangles, we divide the
//    number of boundary edges by 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Marc deBerg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
//    Computational Geometry,
//    Springer, 2000,
//    ISBN: 3-540-65620-0.
//
//  Parameters:
//
//    Input, integer NODE_NUM, the number of nodes.
//
//    Input, integer TRIANGLE_NUM, the number of triangles.
//
//    Input, integer HOLE_NUM, the number of internal nodes.
//
//    Output, int TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT, the number of
//    edges that lie on the boundary of the triangulation.
//
{
  int boundary_num;

  boundary_num = ( 2 * node_num + 2 * hole_num - 4 * triangle_num - 2 ) / 2;

  return boundary_num;
}
//****************************************************************************80

bool *triangulation_order6_boundary_node ( int node_num, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_BOUNDARY_NODE indicates nodes on the boundary.
//
//  Discussion:
//
//    This routine is given an order 6 triangulation, an abstract list of
//    sets of six nodes.  The vertices are listed clockwise, then the
//    midside nodes.
//
//    It is assumed that each edge of the triangulation is either
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes - as long
//    as the boundary of the hole comprises more than 3 edges!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 January 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
//    triangles.
//
//    Output, bool TRIANGULATION_ORDER6_BOUNDARY_NODE[NODE_NUM],
//    is TRUE if the node is on a boundary edge.
//
{
  int e1;
  int e2;
  int *edge;
  bool equal;
  int i;
  int j;
  int m;
  int n;
  bool *node_boundary;

  m = 3;
  n = 3 * triangle_num;
//
//  Set up the edge array.
//
  edge = new int[m*n];

  for ( j = 0; j < triangle_num; j++ )
  {
    edge[0+(j               )*m] = triangle_node[0+j*6];
    edge[1+(j               )*m] = triangle_node[3+j*6];
    edge[2+(j               )*m] = triangle_node[1+j*6];

    edge[0+(j+  triangle_num)*m] = triangle_node[1+j*6];
    edge[1+(j+  triangle_num)*m] = triangle_node[4+j*6];
    edge[2+(j+  triangle_num)*m] = triangle_node[2+j*6];

    edge[0+(j+2*triangle_num)*m] = triangle_node[2+j*6];
    edge[1+(j+2*triangle_num)*m] = triangle_node[5+j*6];
    edge[2+(j+2*triangle_num)*m] = triangle_node[0+j*6];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+j*m], edge[2+j*m] );
    e2 = i4_max ( edge[0+j*m], edge[2+j*m] );
    edge[0+j*m] = e1;
    edge[2+j*m] = e2;
  }
//
//  Ascending sort the column array.
//
  i4col_sort_a ( m, n, edge );
//
//  Records which appear twice are internal edges and can be ignored.
//
  node_boundary = new bool[node_num];

  for ( i = 0; i < node_num; i++ )
  {
    node_boundary[i] = false;
  }

  j = 0;

  while ( j < 3 * triangle_num )
  {
    j = j + 1;

    if ( j == 3 * triangle_num )
    {
      for ( i = 0; i < m; i++ )
      {
        node_boundary[edge[i+(j-1)*m]-1] = true;
      }
      break;
    }

    equal = true;

    for ( i = 0; i < m; i++ )
    {
      if ( edge[i+(j-1)*m] != edge[i+j*m] )
      {
        equal = false;
      }
    }

    if ( equal )
    {
      j = j + 1;
    }
    else
    {
      for ( i = 0; i < m; i++ )
      {
        node_boundary[edge[i+(j-1)*m]-1] = true;
      }
    }

  }

  delete [] edge;

  return node_boundary;
}
//****************************************************************************80

void triangulation_order6_example1 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_EXAMPLE1 sets up a sample triangulation.
//
//  Discussion:
//
//    This triangulation is actually a Delaunay triangulation.
//
//    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
//    determined by calling TRIANGULATION_ORDER6_EXAMPLE1_SIZE first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, int TRIANGLE_ORDER[6*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  Negative values indicate edges that lie on the exterior.
//
{
# define DIM_NUM 2
# define NODE_NUM 48
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 6

  int i;
  int j;
  static double node_xy_save[DIM_NUM*NODE_NUM] = {
       0.0, 0.0,
       1.0, 0.0,
       2.0, 0.0,
       3.0, 0.0,
       4.0, 0.0,
       5.0, 0.0,
       6.0, 0.0,
       0.0, 1.0,
       1.0, 1.0,
       2.0, 1.0,
       3.0, 1.0,
       4.0, 1.0,
       5.0, 1.0,
       6.0, 1.0,
       0.0, 2.0,
       1.0, 2.0,
       2.0, 2.0,
       3.0, 2.0,
       4.0, 2.0,
       5.0, 2.0,
       6.0, 2.0,
       0.0, 3.0,
       1.0, 3.0,
       2.0, 3.0,
       3.0, 3.0,
       5.0, 3.0,
       6.0, 3.0,
       0.0, 4.0,
       1.0, 4.0,
       2.0, 4.0,
       3.0, 4.0,
       4.0, 4.0,
       5.0, 4.0,
       6.0, 4.0,
       0.0, 5.0,
       1.0, 5.0,
       2.0, 5.0,
       3.0, 5.0,
       4.0, 5.0,
       5.0, 5.0,
       6.0, 5.0,
       0.0, 6.0,
       1.0, 6.0,
       2.0, 6.0,
       3.0, 6.0,
       4.0, 6.0,
       5.0, 6.0,
       6.0, 6.0 };
  static int triangle_node_save[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  3, 15,  2,  9,  8,
    17, 15,  3, 16,  9, 10,
     5, 19,  3, 12, 11,  4,
    17,  3, 19, 10, 11, 18,
     7, 21,  5, 14, 13,  6,
    19,  5, 21, 12, 13, 20,
    17, 30, 15, 24, 23, 16,
    28, 15, 30, 22, 23, 29,
    30, 17, 32, 24, 25, 31,
    21, 34, 19, 27, 26, 20,
    30, 44, 28, 37, 36, 29,
    42, 28, 44, 35, 36, 43,
    32, 46, 30, 39, 38, 31,
    44, 30, 46, 37, 38, 45,
    32, 34, 46, 33, 40, 39,
    48, 46, 34, 47, 40, 41 };
  static int triangle_neighbor_save[3*TRIANGLE_NUM] = {
    -3,   2,  -5,
     7,   1,   4,
     6,   4, -11,
     2,   3, -14,
   -15,   6, -17,
     3,   5,  10,
     9,   8,   2,
   -24,   7,  11,
     7, -28,  13,
    27, -31,   6,
     8,  14,  12,
   -36,  11, -38,
    15,  14,   9,
    11,  13, -44,
   -45,  16,  13,
   -48,  15, -50 };

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      node_xy[i+j*DIM_NUM] = node_xy_save[i+j*DIM_NUM];
    }
  }

  for ( j = 0; j < TRIANGLE_NUM; j++ )
  {
    for ( i = 0; i < TRIANGLE_ORDER; i++ )
    {
      triangle_node[i+j*TRIANGLE_ORDER] =
        triangle_node_save[i+j*TRIANGLE_ORDER];
    }
  }

  for ( j = 0; j < TRIANGLE_NUM; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_neighbor[i+j*3] = triangle_neighbor_save[i+j*3];
    }
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void triangulation_order6_example1_size ( int *node_num, int *triangle_num,
  int *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_EXAMPLE1_SIZE sets sizes for a sample triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *TRIANGLE_NUM, the number of triangles.
//
//    Output, int *HOLE_NUM, the number of holes.
//
{
  *node_num = 48;
  *triangle_num = 16;
  *hole_num = 1;

  return;
}
//****************************************************************************80

void triangulation_order6_example2 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_EXAMPLE2 sets up a sample triangulation.
//
//  Discussion:
//
//    This triangulation is actually a Delaunay triangulation.
//
//    The appropriate input values of NODE_NUM and TRIANGLE_NUM can be
//    determined by calling TRIANGULATION_ORDER6_EXAMPLE2_SIZE first.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\  6 |\  8 |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    | 5  \| 7  \|
//   11-12-13-14-15
//    |\  2 |\  4 |
//    | \   | \   |
//    6  7  8  9 10
//    | 1 \ | 3 \ |
//    |    \|    \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Output, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, int TRIANGLE_ORDER[6*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  Negative values indicate edges that lie on the exterior.
//
{
# define DIM_NUM 2
# define NODE_NUM 48
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 6

  int i;
  int j;
  static double node_xy_save[DIM_NUM*NODE_NUM] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    3.0, 3.0,
    4.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0,
    3.0, 4.0,
    4.0, 4.0 };
  static int triangle_node_save[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  3, 11,  2,  7,  6,
    13, 11,  3, 12,  7,  8,
     3,  5, 13,  4,  9,  8,
    15, 13,  5, 14,  9, 10,
    11, 13, 21, 12, 17, 16,
    23, 21, 13, 22, 17, 18,
    13, 15, 23, 14, 19, 18,
    25, 23, 15, 24, 19, 20  };
  static int triangle_neighbor_save[3*TRIANGLE_NUM] = {
    -1,  2, -1,
     5,  1,  3,
    -1,  4,  2,
     7,  3, -1,
     2,  6, -1,
    -1,  5,  7,
     4,  8,  6,
    -1,  7, -1 };

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      node_xy[i+j*DIM_NUM] = node_xy_save[i+j*DIM_NUM];
    }
  }

  for ( j = 0; j < TRIANGLE_NUM; j++ )
  {
    for ( i = 0; i < TRIANGLE_ORDER; i++ )
    {
      triangle_node[i+j*TRIANGLE_ORDER] =
        triangle_node_save[i+j*TRIANGLE_ORDER];
    }
  }

  for ( j = 0; j < TRIANGLE_NUM; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      triangle_neighbor[i+j*3] = triangle_neighbor_save[i+j*3];
    }
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void triangulation_order6_example2_size ( int *node_num, int *triangle_num,
  int *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_EXAMPLE2_SIZE sets sizes for a sample triangulation.
//
//  Diagram:
//
//   21-22-23-24-25
//    |\  6 |\  8 |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    | 5  \| 7  \|
//   11-12-13-14-15
//    |\  2 |\  4 |
//    | \   | \   |
//    6  7  8  9 10
//    | 1 \ | 3 \ |
//    |    \|    \|
//    1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *TRIANGLE_NUM, the number of triangles.
//
//    Output, int *HOLE_NUM, the number of holes.
//
{
  *node_num = 25;
  *triangle_num = 8;
  *hole_num = 0;

  return;
}
//****************************************************************************80

void triangulation_order6_neighbor ( int triangle_num, int triangle_node[],
  int t1, int s1, int  *t2, int *s2 )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_NEIGHBOR determines a neighbor of a given triangle.
//
//  Discussion:
//
//    A set of nodes is given.  A triangulation of the nodes has been
//    defined and recorded in TRIANGLE_NODE.  The TRIANGLE_NODE data
//    structure records triangles as sets of six nodes, with the first three
//    being the vertices, in counterclockwise order.  The fourth node is the
//    midside node between nodes 1 and 2, and the other two are listed
//    logically.
//
//    The nodes of the triangle are listed in counterclockwise order.
//    This means that if two triangles share a side, then the nodes
//    defining that side occur in the order (N1,N2,N3) for one triangle,
//    and (N3,N2,N1) for the other.
//
//    The routine is given a triangle and a side, and asked to find
//    another triangle (if any) that shares that side.  The routine
//    simply searches the TRIANGLE_NODE structure for an occurrence of the
//    nodes in the opposite order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input/output, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that define
//    each triangle.
//
//    Input, int T1, the index of the triangle.
//
//    Input, int S1, the index of the triangle side.
//
//    Output, int *T2, the index of the triangle which is the neighbor
//    to T1 on side S1, or -1 if there is no such neighbor.
//
//    Output, int *S2, the index of the side of triangle T2 which
//    is shared with triangle T1, or -1 if there is no such neighbor.
//
{
  int n1;
  int n2;
  int s;
  int ss;
  int t;

  n1 = triangle_node[s1-1+(t1-1)*6];
  ss = i4_wrap ( s1+1, 1, 3 );
  n2 = triangle_node[ss-1+(t1-1)*6];

  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      if ( triangle_node[s+t*6] == n1 )
      {
        ss = i4_wrap ( s-1, 0, 2 );
        if ( triangle_node[ss+t*6] == n2 )
        {
          *t2 = t + 1;
          *s2 = ss + 1;
          return;
        }
      }
    }
  }

  *t2 = -1;
  *s2 = -1;

  return;
}
//****************************************************************************80

void triangulation_order6_plot ( string file_name, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_PLOT plots a 6-node triangulation of a set of nodes.
//
//  Discussion:
//
//    The triangulation is most usually a Delaunay triangulation,
//    but this is not necessary.
//
//    This routine has been specialized to deal correctly ONLY with
//    a mesh of 6 node elements, with the property that starting
//    from local node 1 and traversing the edges of the element will
//    result in encountering local nodes 1, 4, 2, 5, 3, 6 in that order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to create.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists, for each triangle,
//    the indices of the nodes that form the vertices and midsides
//    of the triangle.
//
//    Input, int NODE_SHOW:
//    0, do not show nodes;
//    1, show nodes;
//    2, show nodes and label them.
//
//    Input, int TRIANGLE_SHOW:
//    0, do not show triangles;
//    1, show triangles;
//    2, show triangles and label them.
//
{
  double ave_x;
  double ave_y;
  int circle_size;
  int delta;
  int e;
  ofstream file_unit;
  int i;
  int ip1;
  int node;
  int order[6] = { 1, 4, 2, 5, 3, 6 };
  int triangle;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = - r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min )
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  file_unit.open ( file_name.c_str ( ) );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "TRIANGULATION_ORDER6_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: triangulation_order6_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";

  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "%  Increase line width from default 0.\n";
  file_unit << "%\n";
  file_unit << "2 setlinewidth\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Draw filled dots at each node:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750  setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps << "  "
        << y_ps << "  "
        << circle_size << " 0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to darker blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850  setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";

    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  "
        << x_ps     << "  "
        << y_ps + 5 << "  moveto ("
        << node+1   << ") show\n";
    }
  }
//
//  Draw the triangles.
//
  if ( 1 <= triangle_show )
  {
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to red.\n";
    file_unit << "%\n";
    file_unit << "0.900  0.200  0.100 setrgbcolor\n";
    file_unit << "%\n";
    file_unit << "%  Draw the triangles.\n";
    file_unit << "%\n";

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      node = triangle_node[order[0]-1+triangle*6] - 1;

      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max                     - y_min ) );

      file_unit << "newpath  " << x_ps << "  " << y_ps << "  moveto\n";

      for ( i = 1; i <= 6; i++ )
      {
        ip1 = ( i % 6 ) + 1;
        node = triangle_node[order[ip1-1]-1+triangle*6] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) )
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) )
          / ( y_max                     - y_min ) );

        file_unit << x_ps << "  " << y_ps << "  lineto\n";
      }
      file_unit << "stroke\n";
    }
  }
//
//  Label the triangles.
//
  if ( 2 <= triangle_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the triangles.\n";
    file_unit << "%\n";
    file_unit << "%  Set the RGB color to darker red.\n";
    file_unit << "%\n";
    file_unit << "0.950  0.250  0.150 setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";
    file_unit << "%\n";

    for ( triangle = 0; triangle < triangle_num; triangle++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 0; i < 6; i++ )
      {
        node = triangle_node[i+triangle*6] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }

      ave_x = ave_x / 6.0;
      ave_y = ave_y / 6.0;

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) )
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) )
        / ( y_max         - y_min ) );

      file_unit << setw(4) << x_ps << "  "
                << setw(4) << y_ps << "  "
                << "moveto (" << triangle+1 << ") show\n";
    }
  }

  file_unit << "%\n";
  file_unit << "restore showpage\n";
  file_unit << "%\n";
  file_unit << "% End of page\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
//****************************************************************************80

void triangulation_order6_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_PRINT prints information defining a triangulation.
//
//  Discussion:
//
//    Triangulations created by R8TRIS2 include extra information encoded
//    in the negative values of TRIANGLE_NEIGHBOR.
//
//    Because some of the nodes counted in node_num may not actually be
//    used in the triangulation, I needed to compute the true number
//    of vertices.  I added this calculation on 13 October 2001.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
//    triangles.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors on each side.
//    If there is no triangle neighbor on a particular side, the value of
//    TRIANGLE_NEIGHBOR should be negative.  If the triangulation data was created by
//    R8TRIS2, then there is more information encoded in the negative values.
//
{
# define DIM_NUM 2

  int boundary_num;
  int i;
  int j;
  int k;
  int n1;
  int n2;
  int n3;
  int s;
  int s1;
  int s2;
  bool skip;
  int t;
  int *vertex_list;
  int vertex_num;

  cout << "\n";
  cout << "TRIANGULATION_ORDER6_PRINT\n";
  cout << "  Information defining a triangulation.\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num << "\n";

  r8mat_transpose_print ( DIM_NUM, node_num, node_xy, "  Node coordinates" );

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";
  cout << "\n";
  cout << "  Sets of six nodes are used as vertices of\n";
  cout << "  the triangles.  For each triangle, the vertices are listed\n";
  cout << "  in counterclockwise order, followed by the midside nodes.\n";

  i4mat_transpose_print ( 6, triangle_num, triangle_node, "  Triangle nodes" );

  cout << "\n";
  cout << "  On each side of a given triangle, there is either\n";
  cout << "  another triangle, or a piece of the convex hull.\n";
  cout << "  For each triangle, we list the indices of the three\n";
  cout << "  neighbors, or (if negative) the codes of the\n";
  cout << "  segments of the convex hull.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_neighbor,
    "  Triangle neighbors" );
//
//  Determine VERTEX_NUM, the number of vertices.
//
  vertex_list = new int[3*triangle_num];

  k = 0;
  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      vertex_list[k] = triangle_node[s+t*6];
      k = k + 1;
    }
  }

  i4vec_sort_heap_a ( 3*triangle_num, vertex_list );

  vertex_num = i4vec_sorted_unique ( 3*triangle_num, vertex_list );

  delete [] vertex_list;
//
//  Determine the number of boundary points.
//
  boundary_num = 2 * vertex_num - triangle_num - 2;

  cout << "\n";
  cout << "  The number of boundary points is " << boundary_num << "\n";
  cout << "\n";
  cout << "  The segments that make up the convex hull can be\n";
  cout << "  determined from the negative entries of the triangle\n";
  cout << "  neighbor list.\n";
  cout << "\n";
  cout << "     #   Tri  Side    N1    N2    N3\n";
  cout << "\n";

  skip = false;

  k = 0;

  for ( i = 0; i < triangle_num; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      if ( triangle_neighbor[j+i*3] < 0 )
      {
        s = -triangle_neighbor[j+i*3];
        t = s / 3;

        if ( t < 1 || triangle_num < t )
        {
          cout << "\n";
          cout << "  Sorry, this data does not use the R8TRIS2\n";
          cout << "  convention for convex hull segments.\n";
          skip = true;
          break;
        }

        s1 = ( s % 3 ) + 1;
        s2 = i4_wrap ( s1+1, 1, 3 );
        k = k + 1;
        n1 = triangle_node[s1-1+(t-1)*6];
        n2 = triangle_node[s1+3-1+(t-1)*6];
        n3 = triangle_node[s2-1+(t-1)*6];
        cout                  << "  "
             << setw(4) << k  << "  "
             << setw(4) << t  << "  "
             << setw(4) << s1 << "  "
             << setw(4) << n1 << "  "
             << setw(4) << n2 << "  "
             << setw(4) << n3 << "\n";
      }
    }

    if ( skip )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void triangulation_order6_refine_compute ( int node_num1, int triangle_num1,
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2,
  int edge_data[], double node_xy2[], int triangle_node2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
//
//  Discussion:
//
//    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
//    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
//    quadratic subtriangles T1, T2, T3 and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.  (In fact, only
//    the new nodes on the edges can be shared, and then only by at most
//    one other triangle.)
//
//            3
//           / \
//          36 35
//         / T3  \
//        6--56---5
//       / \ T4  / \
//      16 46  45  25
//     / T1  \ / T2  \
//    1--14---4--24---2
//
//    This routine is given sorted information defining the edges, and uses
//    it to build the new node and triangle arrays.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes.
//
//    Input, int TRIANGLE_NUM1, the number of triangles.
//
//    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
//
//    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the nodes that make up the
//    triangles.
//
//    Input, int NODE_NUM2, the number of nodes in the refined mesh.
//
//    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
//
//    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
//    by TRIANGULATION_ORDER6_REFINE_SIZE.
//
//    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
//
//    Output, int TRIANGLE_NODE2[6*TRIANGLE_NUM2], the nodes that make up the
//    triangles in the refined mesh.
//
{
  int edge;
  int i;
  int j;
  int l1;
  int l2;
  int l3;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int node;
  int t1;
  int t2;
  int t3;
  int t4;
  int triangle1;
  int v1;
  int v2;
  int v3;
  int v4;
  int v5;
  int v6;
//
//  Step 1:
//  Copy the old nodes.
//
  for ( j = 0; j < node_num1; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      node_xy2[i+j*2] = node_xy1[i+j*2];
    }
  }
  for ( j = 0; j < triangle_num2; j++ )
  {
    for ( i = 0; i < 6; i++ )
    {
      triangle_node2[i+j*6] = -1;
    }
  }
//
//  We can assign the existing nodes to the new triangles.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    t1 = triangle1 * 4 + 0;
    t2 = triangle1 * 4 + 1;
    t3 = triangle1 * 4 + 2;
    t4 = triangle1 * 4 + 3;

    triangle_node2[0+t1*6] = triangle_node1[0+triangle1*6];
    triangle_node2[1+t1*6] = triangle_node1[3+triangle1*6];
    triangle_node2[2+t1*6] = triangle_node1[5+triangle1*6];

    triangle_node2[0+t2*6] = triangle_node1[3+triangle1*6];
    triangle_node2[1+t2*6] = triangle_node1[1+triangle1*6];
    triangle_node2[2+t2*6] = triangle_node1[4+triangle1*6];

    triangle_node2[0+t3*6] = triangle_node1[5+triangle1*6];
    triangle_node2[1+t3*6] = triangle_node1[4+triangle1*6];
    triangle_node2[2+t3*6] = triangle_node1[2+triangle1*6];

    triangle_node2[0+t4*6] = triangle_node1[4+triangle1*6];
    triangle_node2[1+t4*6] = triangle_node1[5+triangle1*6];
    triangle_node2[2+t4*6] = triangle_node1[3+triangle1*6];
  }
//
//  Step 2.
//  Examine sorted edge information.  The first time an edge is encountered,
//  generate two new nodes, then assign them (usually) to the four subtriangles
//  of the two triangles that share that edge.
//
  node = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5] - 1;
    n2 = edge_data[1+edge*5] - 1;

    l1 = edge_data[2+edge*5];
    l3 = edge_data[3+edge*5];

    if ( l1 == 1 && l3 == 2 )
    {
      l2 = 4;
    }
    else if ( l1 == 1 && l3 == 3 )
    {
      l2 = 6;
    }
    else if ( l1 == 2 && l3 == 3 )
    {
      l2 = 5;
    }
    triangle1 = edge_data[4+edge*5];
//
//  If this is the first time we've encountered this edge,
//  create the new nodes.
//
    if ( n1 != n1_old || n2 != n2_old )
    {
      n1_old = n1;
      n2_old = n2;

      v1 = triangle_node1[l1-1+triangle1*6];
      v2 = triangle_node1[l2-1+triangle1*6];
      v3 = triangle_node1[l3-1+triangle1*6];

      for ( i = 0; i < 2; i++ )
      {
        node_xy2[i+node*2] = ( node_xy2[i+(v1-1)*2]
                             + node_xy2[i+(v2-1)*2] ) / 2.0;
      }
      node = node + 1;
      v4 = node;

      for ( i = 0; i < 2; i++ )
      {
        node_xy2[i+node*2] = ( node_xy2[i+(v2-1)*2]
                             + node_xy2[i+(v3-1)*2] ) / 2.0;
      }
      node = node + 1;
      v5 = node;
    }
    t1 = triangle1 * 4 + 0;
    t2 = triangle1 * 4 + 1;
    t3 = triangle1 * 4 + 2;

    if ( l1 == 1 && l3 == 2 )
    {
      if ( triangle_node1[0+triangle1*6] == v1 + 1 )
      {
        triangle_node2[3+t1*6] = v4;
        triangle_node2[3+t2*6] = v5;
      }
      else
      {
        triangle_node2[3+t1*6] = v5;
        triangle_node2[3+t2*6] = v4;
      }
    }
    else if ( l1 == 1 && l3 == 3 )
    {
      if ( triangle_node1[0+triangle1*6] == v1 + 1 )
      {
        triangle_node2[5+t1*6] = v4;
        triangle_node2[5+t3*6] = v5;
      }
      else
      {
        triangle_node2[5+t1*6] = v5;
        triangle_node2[5+t3*6] = v4;
      }
    }
    else if ( l1 == 2 && l3 == 3 )
    {
      if ( triangle_node1[1+triangle1*6] == v1 + 1 )
      {
        triangle_node2[4+t3*6] = v4;
        triangle_node2[4+t2*6] = v5;
      }
      else
      {
        triangle_node2[4+t3*6] = v5;
        triangle_node2[4+t2*6] = v4;
      }
    }
  }
//
//  Step 3.
//  Each old triangle has a single central subtriangle, for which we now
//  need to generate three new "interior" nodes.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    v4 = triangle_node1[3+triangle1*6];
    v5 = triangle_node1[4+triangle1*6];
    v6 = triangle_node1[5+triangle1*6];

    t1 = triangle1 * 4 + 0;
    t2 = triangle1 * 4 + 1;
    t3 = triangle1 * 4 + 2;
    t4 = triangle1 * 4 + 3;

    node_xy2[0+node*2] = 0.5 * ( node_xy1[0+(v5-1)*2] + node_xy1[0+(v6-1)*2] );
    node_xy2[1+node*2] = 0.5 * ( node_xy1[1+(v5-1)*2] + node_xy1[1+(v6-1)*2] );
    node = node + 1;
    triangle_node2[3+t4*6] = node;
    triangle_node2[3+t3*6] = node;

    node_xy2[0+node*2] = 0.5 * ( node_xy1[0+(v6-1)*2] + node_xy1[0+(v4-1)*2] );
    node_xy2[1+node*2] = 0.5 * ( node_xy1[1+(v6-1)*2] + node_xy1[1+(v4-1)*2] );
    node = node + 1;
    triangle_node2[4+t4*6] = node;
    triangle_node2[4+t1*6] = node;

    node_xy2[0+node*2] = 0.5 * ( node_xy1[0+(v4-1)*2] + node_xy1[0+(v5-1)*2] );
    node_xy2[1+node*2] = 0.5 * ( node_xy1[1+(v4-1)*2] + node_xy1[1+(v5-1)*2] );
    node = node + 1;
    triangle_node2[5+t4*6] = node;
    triangle_node2[5+t2*6] = node;
  }

  return;
}
//****************************************************************************80

void triangulation_order6_refine_size ( int node_num1, int triangle_num1,
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
//
//  Discussion:
//
//    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
//    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
//    quadratic subtriangles T1, T2, T3 and T4.
//
//    The task is more complicated by the fact that we are working with
//    a mesh of triangles, so that we want to create a node only once,
//    even though it may be shared by other triangles.  (In fact, only
//    the new nodes on the edges can be shared, and then only by at most
//    one other triangle.)
//
//            3
//           / \
//          36 35
//         / T3  \
//        6--56---5
//       / \ T4  / \
//      16 46  45  25
//     / T1  \ / T2  \
//    1--14---4--24---2
//
//    This routine determines the sizes of the resulting node and
//    triangles, and constructs an edge array that can be used to
//    properly number the new nodes.
//
//    The primary work occurs in sorting a list related to the edges.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM1, the number of nodes in the original mesh.
//
//    Input, int  TRIANGLE_NUM1, the number of triangles in the
//    original mesh.
//
//    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the indices of the nodes
//    that form the triangles in the input mesh.
//
//    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
//
//    Output, int *TRIANGLE_NUM2, the number of triangles in the
//    refined mesh.
//
//    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
//    be needed by TRIANGULATION_ORDER6_REFINE_COMPUTE.
//
{
  int a;
  int b;
  int edge;
  int i;
  int j;
  int k;
  int n1;
  int n1_old;
  int n2;
  int n2_old;
  int triangle1;
//
//  Step 1.
//  From the list of nodes for triangle T, of the form: (I,J,K)
//  construct the edge relations:
//
//    (I,J,1,2,T)
//    (I,K,1,3,T)
//    (J,K,2,3,T)
//
//  In order to make matching easier, we reorder each pair of nodes
//  into ascending order.
//
  for ( triangle1 = 0; triangle1 < triangle_num1; triangle1++ )
  {
    i = triangle_node1[0+triangle1*6];
    j = triangle_node1[1+triangle1*6];
    k = triangle_node1[2+triangle1*6];

    a = i4_min ( i, j );
    b = i4_max ( i, j );

    edge_data[0+5*(3*triangle1+0)] = a;
    edge_data[1+5*(3*triangle1+0)] = b;
    edge_data[2+5*(3*triangle1+0)] = 1;
    edge_data[3+5*(3*triangle1+0)] = 2;
    edge_data[4+5*(3*triangle1+0)] = triangle1;

    a = i4_min ( i, k );
    b = i4_max ( i, k );

    edge_data[0+5*(3*triangle1+1)] = a;
    edge_data[1+5*(3*triangle1+1)] = b;
    edge_data[2+5*(3*triangle1+1)] = 1;
    edge_data[3+5*(3*triangle1+1)] = 3;
    edge_data[4+5*(3*triangle1+1)] = triangle1;

    a = i4_min ( j, k );
    b = i4_max ( j, k );

    edge_data[0+5*(3*triangle1+2)] = a;
    edge_data[1+5*(3*triangle1+2)] = b;
    edge_data[2+5*(3*triangle1+2)] = 2;
    edge_data[3+5*(3*triangle1+2)] = 3;
    edge_data[4+5*(3*triangle1+2)] = triangle1;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1:2; the routine we call here
//  sorts on the full column but that won't hurt us.
//
//  What we need is to find all cases where triangles share an edge.
//  By sorting the columns of the EDGE_DATA array, we will put shared edges
//  next to each other.
//
  i4col_sort_a ( 5, 3*triangle_num1, edge_data );
//
//  Step 3. All the triangles which share an edge show up as consecutive
//  columns with identical first two entries.  Figure out how many new
//  nodes there are, and allocate space for their coordinates.
//
  *node_num2 = node_num1;

  n1_old = -1;
  n2_old = -1;

  for ( edge = 0; edge < 3 * triangle_num1; edge++ )
  {
    n1 = edge_data[0+edge*5];
    n2 = edge_data[1+edge*5];
    if ( n1 != n1_old || n2 != n2_old )
    {
      *node_num2 = *node_num2 + 2;
      n1_old = n1;
      n2_old = n2;
    }
  }

  *node_num2 = *node_num2 + 3 * triangle_num1;

  *triangle_num2 = 4 * triangle_num1;

  return;
}
//****************************************************************************80

int *triangulation_order6_to_order3 ( int triangle_num1, int triangle_node1[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_TO_ORDER3 linearizes a quadratic triangulation.
//
//  Discussion:
//
//    A quadratic triangulation is assumed to consist of 6-node triangles,
//    as in the following:
//
//    11-12-13-14-15
//     |\    |\    |
//     | \   | \   |
//     6  7  8  9 10
//     |   \ |   \ |
//     |    \|    \|
//     1--2--3--4--5
//
//   This routine rearranges information so as to define the 3-node
//   triangulation:
//
//    11-12-13-14-15
//     |\ |\ |\ |\ |
//     | \| \| \| \|
//     6--7--8--9-10
//     |\ |\ |\ |\ |
//     | \| \| \| \|
//     1--2--3--4--5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int TRIANGLE_NUM1, the number of triangles in the quadratic
//    triangulation.
//
//    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the quadratic triangulation.
//
//    Output, int TRIANGULATION_ORDER6_TO_ORDER3[3*TRIANGLE_NUM2], the linear
//    triangulation.  Here, TRIANGLE_NUM2 = 4 * TRIANGLE_NUM1.
//
{
  int n1;
  int n2;
  int n3;
  int n4;
  int n5;
  int n6;
  int triangle_num2;
  int tri1;
  int tri2;
  int *triangle_node2;

  triangle_num2 = 4 * triangle_num1;
  triangle_node2 = new int[3*triangle_num2];

  tri2 = 0;

  for ( tri1 = 0; tri1 < triangle_num1; tri1++ )
  {
    n1 = triangle_node1[0+tri1*6];
    n2 = triangle_node1[1+tri1*6];
    n3 = triangle_node1[2+tri1*6];
    n4 = triangle_node1[3+tri1*6];
    n5 = triangle_node1[4+tri1*6];
    n6 = triangle_node1[5+tri1*6];

    triangle_node2[0+tri2*3] = n1;
    triangle_node2[1+tri2*3] = n4;
    triangle_node2[2+tri2*3] = n6;
    tri2 = tri2 + 1;

    triangle_node2[0+tri2*3] = n2;
    triangle_node2[1+tri2*3] = n5;
    triangle_node2[2+tri2*3] = n4;
    tri2 = tri2 + 1;

    triangle_node2[0+tri2*3] = n3;
    triangle_node2[1+tri2*3] = n6;
    triangle_node2[2+tri2*3] = n5;
    tri2 = tri2 + 1;

    triangle_node2[0+tri2*3] = n4;
    triangle_node2[1+tri2*3] = n5;
    triangle_node2[2+tri2*3] = n6;
    tri2 = tri2 + 1;
  }

  return triangle_node2;
}
//****************************************************************************80

int triangulation_order6_vertex_count ( int tri_num, int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_ORDER6_VERTEX_COUNT counts vertex nodes in a triangulation.
//
//  Discussion:
//
//    In a triangulation of order 6, some nodes are midside nodes and some
//    nodes are vertex nodes.
//
//    Especially when an order 6 triangulation is used to handle the
//    Navier Stokes equations, it is useful to know the number of
//    vertex and midside nodes.
//
//    Note that the number of midside nodes is simple NODE_NUM - VERTEX_NUM.
//
//  Diagram:
//
//       3
//    s  |\
//    i  | \
//    d  |  \
//    e  6   5  side 2
//       |    \
//    3  |     \
//       |      \
//       1---4---2
//
//         side 1
//
//    The local node numbering.  Local nodes 1, 2 and 3 are "vertex nodes",
//    while nodes 4, 5 and 6 are "midside nodes".
//
//
//   21-22-23-24-25
//    |\    |\    |
//    | \   | \   |
//   16 17 18 19 20
//    |   \ |   \ |
//    |    \|    \|
//   11-12-13-14-15
//    |\    |\    |
//    | \   | \   |
//    6  7  8  9 10
//    |   \ |   \ |
//    |    \|    \|
//    1--2--3--4--5
//
//    A sample grid, which contains 9 vertex nodes and 16 midside nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[6*TRI_NUM], lists the nodes that
//    make up each triangle.  The first three nodes are the vertices,
//    in counterclockwise order.  The fourth value is the midside
//    node between nodes 1 and 2; the fifth and sixth values are
//    the other midside nodes in the logical order.
//
//    Output, int TRIANGULATION_ORDER6_VERTEX_COUNT, the number of nodes
//    which are vertices.
//
{
  int j;
  int vertex_num;
  int *vertices;

  vertices = new int[3*tri_num];

  for ( j = 0; j < tri_num; j++ )
  {
    vertices[j] = triangle_node[0+j*6];
  }
  for ( j = 0; j < tri_num; j++ )
  {
    vertices[tri_num+j] = triangle_node[1+j*6];
  }
  for ( j = 0; j < tri_num; j++ )
  {
    vertices[2*tri_num+j] = triangle_node[2+j*6];
  }

  i4vec_sort_heap_a ( 3*tri_num, vertices );

  vertex_num = i4vec_sorted_unique ( 3*tri_num, vertices );

  delete [] vertices;

  return vertex_num;
}
//****************************************************************************80

void triangulation_search_delaunay ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double p[2], int *triangle_index, 
  double *alpha, double *beta, double *gamma, int *edge,
  int *step_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_SEARCH_DELAUNAY searches a triangulation for a point.
//
//  Discussion:
//
//    The algorithm "walks" from one triangle to its neighboring triangle,
//    and so on, until a triangle is found containing point P, or P is found
//    to be outside the convex hull.
//
//    The algorithm computes the barycentric coordinates of the point with
//    respect to the current triangle.  If all three quantities are positive,
//    the point is contained in the triangle.  If the I-th coordinate is
//    negative, then (X,Y) lies on the far side of edge I, which is opposite
//    from vertex I.  This gives a hint as to where to search next.
//
//    For a Delaunay triangulation, the search is guaranteed to terminate.
//    For other triangulations, a cycle may occur.
//
//    Note the surprising fact that, even for a Delaunay triangulation of
//    a set of nodes, the nearest point to (X,Y) need not be one of the
//    vertices of the triangle containing (X,Y).
//
//    The code can be called for triangulations of any order, but only
//    the first three nodes in each triangle are considered.  Thus, if
//    higher order triangles are used, and the extra nodes are intended
//    to give the triangle a polygonal shape, these will have no effect,
//    and the results obtained here might be misleading.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 October 2012
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the nodes of each triangle.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
//
//    Input, double P[2], the coordinates of a point.
//
//    Output, int *TRIANGLE_INDEX, the index of the triangle where the search ended.
//    If a cycle occurred, then TRIANGLE_INDEX = -1.
//
//    Output, double *ALPHA, *BETA, *GAMMA, the barycentric coordinates
//    of the point with respect to triangle *TRIANGLE_INDEX.
//
//    Output, int *EDGE, indicates the position of the point (X,Y) in
//    triangle TRIANGLE:
//    0, the interior or boundary of the triangle;
//    -1, outside the convex hull of the triangulation, past edge 1;
//    -2, outside the convex hull of the triangulation, past edge 2;
//    -3, outside the convex hull of the triangulation, past edge 3.
//
//    Output, int *STEP_NUM, the number of steps.
{
  int a;
  int b;
  int c;
  double det;
  double dxp;
  double dxa;
  double dxb;
  double dyp;
  double dya;
  double dyb;
  static int triangle_index_save = -1;

  *step_num = - 1;
  *edge = 0;

  if ( triangle_index_save < 0 || triangle_num <= triangle_index_save )
  {
    *triangle_index = ( triangle_num + 1 ) / 2;
  }
  else
  {
    *triangle_index = triangle_index_save;
  }

  for ( ; ; )
  {
    *step_num = *step_num + 1;

    if ( triangle_num < *step_num )
    {
      cout << "\n";
      cout << "TRIANGULATION_SEARCH_DELAUNAY - Fatal error!\n";
      cout << "  The algorithm seems to be cycling.\n";
      cout << "  Current triangle is " << *triangle_index << "\n";
      *triangle_index = -1;
      *alpha = -1.0;
      *beta = -1.0;
      *gamma = -1.0;
      *edge = -1;
      return;
    }
//
//  Get the vertices of triangle TRIANGLE.
//
    a = triangle_node[0+(*triangle_index-1)*triangle_order];
    b = triangle_node[1+(*triangle_index-1)*triangle_order];
    c = triangle_node[2+(*triangle_index-1)*triangle_order];
//
//  Using vertex C as a base, compute the distances to vertices A and B,
//  and the point (X,Y).
//
    dxa = node_xy[0+a*2] - node_xy[0+c*2];
    dya = node_xy[1+a*2] - node_xy[1+c*2];

    dxb = node_xy[0+b*2] - node_xy[0+c*2];
    dyb = node_xy[1+b*2] - node_xy[1+c*2];

    dxp = p[0]           - node_xy[0+c*2];
    dyp = p[1]           - node_xy[1+c*2];

    det = dxa * dyb - dya * dxb;
//
//  Compute the barycentric coordinates of the point (X,Y) with respect
//  to this triangle.
//
    *alpha = ( dxp * dyb - dyp * dxb ) / det;
    *beta =  ( dxa * dyp - dya * dxp ) / det;
    *gamma = 1.0 - *alpha - *beta;
//
//  If the barycentric coordinates are all positive, then the point
//  is inside the triangle and we're done.
//
    if ( 0.0 <= *alpha &&
         0.0 <= *beta  &&
         0.0 <= *gamma )
    {
      break;
    }
//
//  At least one barycentric coordinate is negative.
//
//  If there is a negative barycentric coordinate for which there exists
//  an opposing triangle neighbor closer to the point, move to that triangle.
//
//  (Two coordinates could be negative, in which case we could go for the
//  most negative one, or the most negative one normalized by the actual
//  distance it represents).
//
    if ( *alpha < 0.0 && 0 <= triangle_neighbor[1+(*triangle_index-1)*3] )
    {
      *triangle_index = triangle_neighbor[1+(*triangle_index-1)*3];
      continue;
    }
    else if ( *beta < 0.0 && 0 <= triangle_neighbor[2+(*triangle_index-1)*3] )
    {
      *triangle_index = triangle_neighbor[2+(*triangle_index-1)*3];
      continue;
    }
    else if ( *gamma < 0.0 && 0 <= triangle_neighbor[0+(*triangle_index-1)*3] )
    {
      *triangle_index = triangle_neighbor[0+(*triangle_index-1)*3];
      continue;
    }
//
//  All negative barycentric coordinates correspond to vertices opposite
//  sides on the convex hull.
//
//  Note the edge and exit.
//
    if ( *alpha < 0.0 )
    {
      *edge = -2;
      break;
    }
    else if ( *beta < 0.0 )
    {
      *edge = -3;
      break;
    }
    else if ( *gamma < 0.0 )
    {
      *edge = -1;
      break;
    }
    else
    {
      cout << "\n";
      cout << "TRIANGULATION_ORDER3_SEARCH - Fatal error!\n";
      cout << "  The algorithm seems to have reached a dead end\n";
      cout << "  after " << *step_num << " steps.\n";
      *triangle_index = -1;
      *edge = -1;
      return;
    }
  }
  triangle_index_save = *triangle_index;

  return;
}
//****************************************************************************80

int triangulation_search_naive ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[], double p[2] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_SEARCH_NAIVE naively searches a triangulation for a point.
//
//  Discussion:
//
//    The algorithm simply checks each triangle to see if point P is
//    contained in it.  Surprisingly, this is not the fastest way to
//    do the check, at least if the triangulation is Delaunay.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the nodes of each triangle.
//
//    Input, double P[2], the coordinates of a point.
//
//    Output, int TRIANGULATION_SEARCH_NAIVE, the index of the triangle
//    containing the point, or -1 if no triangle was found containing
//    the point.
//
{
  int a;
  double alpha;
  int b;
  double beta;
  int c;
  double det;
  double dxp;
  double dxa;
  double dxb;
  double dyp;
  double dya;
  double dyb;
  double gamma;
  int triangle;
  int triangle_index;

  triangle_index = -1;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
//
//  Get the vertices of triangle TRIANGLE.
//
    a = triangle_node[0+triangle*triangle_order];
    b = triangle_node[1+triangle*triangle_order];
    c = triangle_node[2+triangle*triangle_order];
//
//  Using vertex C as a base, compute the distances to vertices A and B,
//  and the point (X,Y).
//
    dxa = node_xy[0+a*2] - node_xy[0+c*2];
    dya = node_xy[1+a*2] - node_xy[1+c*2];

    dxb = node_xy[0+b*2] - node_xy[0+c*2];
    dyb = node_xy[1+b*2] - node_xy[1+c*2];

    dxp = p[0]           - node_xy[0+c*2];
    dyp = p[1]           - node_xy[1+c*2];

    det = dxa * dyb - dya * dxb;
//
//  Compute the barycentric coordinates of the point (X,Y) with respect
//  to this triangle.
//
    alpha = ( dxp * dyb - dyp * dxb ) / det;
    beta =  ( dxa * dyp - dya * dxp ) / det;
    gamma = 1.0 - alpha - beta;
//
//  If the barycentric coordinates are all positive, then the point
//  is inside the triangle and we're done.
//
    if ( 0.0 <= alpha &&
         0.0 <= beta  &&
         0.0 <= gamma )
    {
      triangle_index = triangle + 1;
      break;
    }
  }

  return triangle_index;
}
//****************************************************************************80

void vbedg ( double x, double y, int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int triangle_neighbor[], int *ltri,
  int *ledg, int *rtri, int *redg )

//****************************************************************************80
//
//  Purpose:
//
//    VBEDG determines which boundary edges are visible to a point.
//
//  Discussion:
//
//    The point (X,Y) is assumed to be outside the convex hull of the
//    region covered by the 2D triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2008
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point outside the convex hull
//    of the current triangulation.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence list.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list;
//    negative values are used for links of a counter clockwise linked list
//    of boundary edges;
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
//    assumed to be already computed and are not changed, else they are updated.
//    On output, LTRI is the index of boundary triangle to the left of the
//    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
//    edge of triangle LTRI to the left of the leftmost boundary edge visible
//    from (X,Y).  1 <= LEDG <= 3.
//
//    Input/output, int *RTRI.  On input, the index of the boundary triangle
//    to begin the search at.  On output, the index of the rightmost boundary
//    triangle visible from (X,Y).
//
//    Input/output, int *REDG, the edge of triangle RTRI that is visible
//    from (X,Y).  1 <= REDG <= 3.
//
{
  int a;
  double ax;
  double ay;
  int b;
  double bx;
  double by;
  bool done;
  int e;
  int l;
  int lr;
  int t;
//
//  Find the rightmost visible boundary edge using links, then possibly
//  leftmost visible boundary edge using triangle neighbor information.
//
  if ( *ltri == 0 )
  {
    done = false;
    *ltri = *rtri;
    *ledg = *redg;
  }
  else
  {
    done = true;
  }

  for ( ; ; )
  {
    l = -triangle_neighbor[3*((*rtri)-1)+(*redg)-1];
    t = l / 3;
    e = 1 + l % 3;
    a = triangle_node[3*(t-1)+e-1];

    if ( e <= 2 )
    {
      b = triangle_node[3*(t-1)+e];
    }
    else
    {
      b = triangle_node[3*(t-1)+0];
    }

    ax = node_xy[2*(a-1)+0];
    ay = node_xy[2*(a-1)+1];

    bx = node_xy[2*(b-1)+0];
    by = node_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

    *rtri = t;
    *redg = e;

  }

  if ( done )
  {
    return;
  }

  t = *ltri;
  e = *ledg;

  for ( ; ; )
  {
    b = triangle_node[3*(t-1)+e-1];
    e = i4_wrap ( e-1, 1, 3 );

    while ( 0 < triangle_neighbor[3*(t-1)+e-1] )
    {
      t = triangle_neighbor[3*(t-1)+e-1];

      if ( triangle_node[3*(t-1)+0] == b )
      {
        e = 3;
      }
      else if ( triangle_node[3*(t-1)+1] == b )
      {
        e = 1;
      }
      else
      {
        e = 2;
      }

    }

    a = triangle_node[3*(t-1)+e-1];
    ax = node_xy[2*(a-1)+0];
    ay = node_xy[2*(a-1)+1];

    bx = node_xy[2*(b-1)+0];
    by = node_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

  }

  *ltri = t;
  *ledg = e;

  return;
}
//****************************************************************************80

double voronoi_polygon_area ( int node, int neighbor_num,
  int neighbor_index[], int node_num, double node_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    VORONOI_POLYGON_AREA computes the area of a Voronoi polygon.
//
//  Formula:
//
//    It is assumed that the Voronoi polygon is finite!  Every Voronoi
//    diagram includes some regions which are infinite, and for those,
//    this formula is not appropriate.
//
//    The routine is given the indices of the nodes that are
//    Voronoi "neighbors" of a given node.  These are also the nodes
//    that are paired to form edges of Delaunay triangles.
//
//    The assumption that the polygon is a Voronoi polygon is
//    used to determine the location of the boundaries of the polygon,
//    which are the perpendicular bisectors of the lines connecting
//    the center point to each of its neighbors.
//
//    The finiteness assumption is employed in part in the
//    assumption that the polygon is bounded by the finite
//    line segments from point 1 to 2, 2 to 3, ...,
//    M-1 to M, and M to 1, where M is the number of neighbors.
//
//    It is assumed that this routine is being called by a
//    process which has computed the Voronoi diagram of a large
//    set of nodes, so the arrays X and Y are dimensioned by
//    NODE_NUM, which may be much greater than the number of neighbor
//    nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
//    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
//    Second Edition,
//    Wiley, 2000, page 485.
//
//  Parameters:
//
//    Input, int NODE, the index of the node whose Voronoi
//    polygon is to be measured. 0 <= NODE < NODE_NUM.
//
//    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
//    the given node.
//
//    Input, int NEIGHBOR_INDEX[NEIGHBOR_NUM], the indices
//    of the neighbor nodes (used to access X and Y).  The neighbor
//    nodes should be listed in the (counter-clockwise) order in
//    which they occur as one circles the center node.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, double VORONOI_POLYGON_AREA, the area of the Voronoi polygon.
//
{
  double a;
  double area;
  double b;
  double c;
  int i;
  int ip1;
  double ui;
  double uip1;
  double vi;
  double vip1;
  double xc;
  double xi;
  double xip1;
  double yc;
  double yi;
  double yip1;

  area = 0.0;

  if ( node < 0 || node_num <= node )
  {
    cout << "\n";
    cout << "  VORONOI_POLYGON_AREA - Fatal error!\n";
    cout << "  Illegal value of input parameter NODE.\n";
    exit ( 1 );
  }

  xc = node_xy[0+node*2];
  yc = node_xy[1+node*2];

  i = neighbor_num - 1;
  i = neighbor_index[i];

  xi = node_xy[0+i*2];
  yi = node_xy[1+i*2];

  ip1 = 0;
  ip1 = neighbor_index[ip1];

  xip1 = node_xy[0+ip1*2];
  yip1 = node_xy[1+ip1*2];
  a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
  b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
  c = 2.0 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
  uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
  vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

  for ( i = 0; i < neighbor_num; i++ )
  {
    xi = xip1;
    yi = yip1;
    ui = uip1;
    vi = vip1;

    ip1 = i4_wrap ( i+1, 0, neighbor_num-1 );
    ip1 = neighbor_index[ip1];

    xip1 = node_xy[0+ip1*2];
    yip1 = node_xy[1+ip1*2];
    a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
    b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
    c = 2.0 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
    uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
    vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

    area = area + uip1 * vi - ui * vip1;
  }

  area = 0.5 * area;

  return area;
}
//****************************************************************************80

double *voronoi_polygon_centroid ( int node, int neighbor_num,
  int neighbor_index[], int node_num, double node_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    VORONOI_POLYGON_CENTROID_2D computes the centroid of a Voronoi polygon.
//
//  Formula:
//
//    It is assumed that the Voronoi polygon is finite!  Every Voronoi
//    diagram includes some regions which are infinite, and for those,
//    this formula is not appropriate.
//
//    The routine is given the indices of the nodes that are
//    Voronoi "neighbors" of a given node.  These are also the nodes
//    that are paired to form edges of Delaunay triangles.
//
//    The assumption that the polygon is a Voronoi polygon is
//    used to determine the location of the boundaries of the polygon,
//    which are the perpendicular bisectors of the lines connecting
//    the center point to each of its neighbors.
//
//    The finiteness assumption is employed in part in the
//    assumption that the polygon is bounded by the finite
//    line segments from point 1 to 2, 2 to 3, ...,
//    M-1 to M, and M to 1, where M is the number of neighbors.
//
//    It is assumed that this routine is being called by a
//    process which has computed the Voronoi diagram of a large
//    set of nodes, so the arrays X and Y are dimensioned by
//    NODE_NUM, which may be much greater than the number of neighbor
//    nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
//    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
//    Second Edition,
//    Wiley, 2000, page 490.
//
//  Parameters:
//
//    Input, int NODE, the index of the node whose Voronoi
//    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
//
//    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
//    the given node.
//
//    Input, int NEIGHBOR_INDEX[NEIGHBOR_NUM], the indices
//    of the neighbor nodes.  These indices are used to access the
//    X and Y arrays.  The neighbor nodes should be listed in the
//    (counter-clockwise) order in which they occur as one circles
//    the center node.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, double *VORONOI_POLYGON_CENTROID_2D, a pointer to a 2D array
//    containing the coordinates of the centroid of the Voronoi polygon
//    of node NODE.
//
{
  double a;
  double area;
  double b;
  double c;
  double *centroid;
  int i;
  int ip1;
  double ui;
  double uip1;
  double vi;
  double vip1;
  double xc;
  double xi;
  double xip1;
  double yc;
  double yi;
  double yip1;

  centroid = new double[2];

  centroid[0] = 0.0;
  centroid[1] = 0.0;

  if ( node < 0 || node_num <= node )
  {
    cout << "\n";
    cout << "VORONOI_POLYGON_CENTROID - Fatal error!\n";
    cout << "  Illegal value of input parameter NODE.\n";
    exit ( 1 );
  }

  xc = node_xy[0+node*2];
  yc = node_xy[1+node*2];

  i = neighbor_num - 1;
  i = neighbor_index[i];

  xi = node_xy[0+i*2];
  yi = node_xy[1+i*2];

  ip1 = 0;
  ip1 = neighbor_index[ip1];

  xip1 = node_xy[0+ip1*2];
  yip1 = node_xy[1+ip1*2];
  a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
  b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
  c = 2.0 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
  uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
  vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

  for ( i = 0; i < neighbor_num; i++ )
  {
    xi = xip1;
    yi = yip1;
    ui = uip1;
    vi = vip1;

    ip1 = i4_wrap ( i+1, 0, neighbor_num-1 );
    ip1 = neighbor_index[ip1];

    xip1 = node_xy[0+ip1*2];
    yip1 = node_xy[1+ip1*2];
    a = ( xi   * xi   + yi   * yi   - xc * xc - yc * yc );
    b = ( xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc );
    c = 2.0 * ( ( xi - xc ) * ( yip1 - yc ) - ( xip1 - xc ) * ( yi - yc ) );
    uip1 = ( a * ( yip1 - yc ) - b * ( yi - yc )  ) / c;
    vip1 = ( a * ( xip1 - xc ) - b * ( xi - xc )  ) / c;

    centroid[0] = centroid[0] + ( vi - vip1 )
      * ( ( uip1 + ui ) * ( uip1 + ui ) - uip1 * ui );
    centroid[1] = centroid[1] + ( ui - uip1 )
      * ( ( vip1 + vi ) * ( vip1 + vi ) - vip1 * vi );
  }

  area = voronoi_polygon_area ( node, neighbor_num, neighbor_index,
    node_num, node_xy );

  centroid[0] = centroid[0] / ( 6.0 * area );
  centroid[1] = centroid[1] / ( 6.0 * area );

  return centroid;
}
//****************************************************************************80

void voronoi_polygon_vertices ( int node, int neighbor_num,
  int neighbor_index[], int node_num, double node_xy[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    VORONOI_POLYGON_VERTICES_2D computes the vertices of a Voronoi polygon.
//
//  Formula:
//
//    This routine is only appropriate for Voronoi polygons that are finite.
//
//    The routine is given the indices of the nodes that are neighbors of a
//    given "center" node.  A node is a neighbor of the center node if the
//    Voronoi polygons of the two nodes share an edge.  The triangles of the
//    Delaunay triangulation are formed from successive pairs of these neighbor
//    nodes along with the center node.
//
//    Given only the neighbor node information, it is possible to determine
//    the location of the vertices of the polygonal Voronoi region by computing
//    the circumcenters of the Delaunay triangles.
//
//    It is assumed that this routine is being called by a process which has
//    computed the Voronoi diagram of a large set of nodes, so the arrays X and
//    Y are dimensioned by NODE_NUM, which may be much greater than the number
//    of neighbor nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
//    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
//    Second Edition,
//    Wiley, 2000.
//
//  Parameters:
//
//    Input, int NODE, the index of the node whose Voronoi
//    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
//
//    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
//    the given node.
//
//    Input, int NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
//    of the neighbor nodes.  These indices are used to access the
//    X and Y arrays.  The neighbor nodes should be listed in the
//    (counter-clockwise) order in which they occur as one circles
//    the center node.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Output, double V[2*NEIGHBOR_NUM], the vertices of the Voronoi polygon
//    around node NODE.
//
{
# define DIM_NUM 2

  double *center;
  int i;
  int ip1;
  double t[DIM_NUM*3];

  if ( node < 0 || node_num <= node )
  {
    cout << "\n";
    cout << "VORONOI_POLYGON_VERTICES - Fatal error!\n";
    cout << "  Illegal value of input parameter NODE.\n";
    exit ( 1 );
  }

  t[0+0*2] = node_xy[0+node*2];
  t[1+0*2] = node_xy[1+node*2];

  ip1 = neighbor_index[0];
  t[0+2*2] = node_xy[0+ip1*2];
  t[1+2*2] = node_xy[1+ip1*2];

  for ( i = 0; i < neighbor_num; i++ )
  {
    t[0+1*2] = t[0+2*2];
    t[1+1*2] = t[1+2*2];

    ip1 = i4_wrap ( i+1, 0, neighbor_num-1 );
    ip1 = neighbor_index[ip1];

    t[0+2*2] = node_xy[0+ip1*2];
    t[1+2*2] = node_xy[1+ip1*2];

    center = triangle_circumcenter_2d ( t );

    v[0+i*2] = center[0];
    v[1+i*2] = center[1];

    delete [] center;
  }

  return;
# undef DIM_NUM
}
