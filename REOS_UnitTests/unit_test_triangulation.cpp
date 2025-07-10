# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "triangulation.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test025 ( );
void test026 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

void test10 ( );
void test11 ( );
void test12 ( );
void test125 ( );
void test127 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );

void test20 ( );
void test21 ( );
void test213 ( );
void quad_fun ( int n, double xy_vec[], double f_vec[] );
void test215 ( );
void test217 ( );
void test219 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test265 ( );
void test27 ( );

void test31 ( );
void test32 ( );
void test33 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGULATION_PRB.
//
//  Discussion:
//
//    TRIANGULATION_PRB tests the TRIANGULATION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRIANGULATION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGULATION library.\n";

  test01 ( );
  test02 ( );
  test025 ( );
  test026 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test125 ( );
  test127 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test213 ( );
  test215 ( );
  test217 ( );
  test219 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test265 ( );
  test27 ( );

  test31 ( );
  test32 ( );
  test33 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGULATION_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests ALPHA_MEASURE.
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
{
  double alpha_ave;
  double alpha_area;
  double alpha_min;
  int hole_num;
  int node_num;
  double *node_xy;
  double quality;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ALPHA_MEASURE returns the ALPHA measure of\n";
  cout << "  quality of a triangulation.\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );
//
//  Allocate space.
//
  node_xy = new double[2*node_num];
  triangle_node = new int [triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];
//
//  Get the triangulation data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Compute the triangulation quality.
//
  alpha_measure ( node_num, node_xy, triangle_order, triangle_num,
    triangle_node, &alpha_min, &alpha_ave, &alpha_area );

  cout << "\n";
  cout << "  ALPHA_MIN  = " << alpha_min << "\n";
  cout << "  ALPHA_AVE  = " << alpha_ave << "\n";
  cout << "  ALPHA_AREA = " << alpha_area << "\n";
//
//  Free the memory.
//
  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests AREA_MEASURE.
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
{
  double area_ave;
  double area_max;
  double area_min;
  double area_ratio;
  double area_std;
  int hole_num;
  int node_num;
  double *node_xy;
  double quality;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  AREA_MEASURE returns the AREA measure of\n";
  cout << "  quality of a triangulation.\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );
//
//  Allocate space.
//
  node_xy = new double[2*node_num];
  triangle_node = new int [triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];
//
//  Get the triangulation data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Compute the triangulation quality.
//
  area_measure ( node_num, node_xy, triangle_order, triangle_num,
    triangle_node, &area_min, &area_max, &area_ratio, &area_ave, &area_std );

  cout << "\n";
  cout << "  AREA_MIN   = " << area_min << "\n";
  cout << "  AREA_MAX   = " << area_max << "\n";
  cout << "  AREA_RATIO = " << area_ratio << "\n";
  cout << "  AREA_AVE   = " << area_ave << "\n";
  cout << "  AREA_STD   = " << area_std << "\n";
//
//  Free the memory.
//
  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests DELAUNAY_SWAP_TEST.
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
{
  int node_num = 4;
  int triangle_num = 2;
  int triangle_order = 3;

  double alpha_area;
  double alpha_ave;
  double alpha_min_swapped;
  double alpha_min_unswapped;
  double node_xy[2*4];
  int seed = 123456789;
  bool swap;
  int test;
  int test_num = 10;
  int triangle_node[3*2];
  int value;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  DELAUNAY_SWAP_TEST determines whether two triangles\n";
  cout << "  with a common edge need to \"swap\" diagonals.\n";
  cout << "  If swapping is indicated, then ALPHA_MIN should increase.\n";
  cout << "\n";
  cout << "  Swap   ALPHA_MIN   ALPHA_MIN\n";
  cout << "         Unswapped   Swapped\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Generate a random quadrilateral (1,2,3,4).
//
    quad_convex_random ( &seed, node_xy );
//
//  Does it need swapping?
//
    swap = delaunay_swap_test ( node_xy );
//
//  Compute ALPHA_MIN unswapped.
//
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 3;
    triangle_node[0+1*3] = 1;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_unswapped, &alpha_ave, &alpha_area );
//
//  Compute ALPHA_MIN swapped.
//
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 4;
    triangle_node[0+1*3] = 2;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_swapped, &alpha_ave, &alpha_area );

    if ( false )
    {
      r8mat_transpose_print ( 2, node_num, node_xy, "  Quadrilateral" );
    }

    cout << "     " << swap
         << "  " << setw(10) << alpha_min_unswapped
         << "  " << setw(10) << alpha_min_swapped << "\n";
  }

  return;
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests DIAEDG.
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
{
  int node_num = 4;
  int triangle_num = 2;
  int triangle_order = 3;

  double alpha_area;
  double alpha_ave;
  double alpha_min_swapped;
  double alpha_min_unswapped;
  double node_xy[2*4];
  int seed = 123456789;
  bool swap;
  int test;
  int test_num = 10;
  int triangle_node[3*2];
  int value;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  DIAEDG determines whether two triangles\n";
  cout << "  with a common edge need to \"swap\" diagonals.\n";
  cout << "  If swapping is indicated, then ALPHA_MIN should increase.\n";
  cout << "\n";
  cout << "  Swap   ALPHA_MIN   ALPHA_MIN\n";
  cout << "         Unswapped   Swapped\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Generate a random quadrilateral (1,2,3,4).
//
    quad_convex_random ( &seed, node_xy );
//
//  Does it need swapping?
//
    value = diaedg ( 
      node_xy[0+0*2], node_xy[1+0*2], 
      node_xy[0+1*2], node_xy[1+1*2], 
      node_xy[0+2*2], node_xy[1+2*2], 
      node_xy[0+3*2], node_xy[1+3*2] );

    if ( value == 1 )
    {
      swap = false;
    }
    else
    {
      swap = true;
    }
//
//  Compute ALPHA_MIN unswapped.
//
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 3;
    triangle_node[0+1*3] = 1;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_unswapped, &alpha_ave, &alpha_area );
//
//  Compute ALPHA_MIN swapped.
//
    triangle_node[0+0*3] = 1;
    triangle_node[1+0*3] = 2;
    triangle_node[2+0*3] = 4;
    triangle_node[0+1*3] = 2;
    triangle_node[1+1*3] = 3;
    triangle_node[2+1*3] = 4;

    alpha_measure ( node_num, node_xy, triangle_order, triangle_num, 
      triangle_node, &alpha_min_swapped, &alpha_ave, &alpha_area );

    if ( false )
    {
      r8mat_transpose_print ( 2, node_num, node_xy, "  Quadrilateral" );
    }

    cout << "     " << swap
         << "  " << setw(10) << alpha_min_unswapped
         << "  " << setw(10) << alpha_min_swapped << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests NODE_MERGE.
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
{
# define DIM_NUM 2
# define NODE_NUM 15
# define TEST_NUM 4

  int node;
  int node_rep[NODE_NUM];
  double node_xy[DIM_NUM*NODE_NUM] = {
       0.0, 0.0, 
       1.0, 0.0,
       3.0, 0.0,
       4.0, 0.0,
       1.0, 1.0,
       4.0, 1.0,
       2.0, 2.0,
       3.0, 3.0,
       2.0, 3.5,
       0.5, 4.0,
       1.0, 4.0,
       1.5, 4.0,
       4.0, 4.0,
       1.0, 4.5,
       1.0, 4.5 };
  int rep;
  int rep_num;
  int test;
  double tolerance;
  double tolerance_test[TEST_NUM] = { 0.01, 0.75, 1.2, 1.5 };

  cout << "\n";
  cout << "TEST03\n";
  cout << "  NODE_MERGE identifies groups of nodes\n";
  cout << "  that can be merged, with a given tolerance.\n";

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_xy, "  Node coordinates:" );

  for ( test = 0; test < TEST_NUM; test++ )
  {
    tolerance = tolerance_test[test];

    node_merge ( DIM_NUM, NODE_NUM, node_xy, tolerance, node_rep );

    cout << "\n";
    cout << "  TOLERANCE = " << tolerance << "\n";
    cout << "\n";
    cout << "      Node  Representatives:\n";
    cout << "\n";

    for ( node = 0; node < NODE_NUM; node++ )
    {
      cout << "  " << setw(8) << node
           << "  " << setw(8) << node_rep[node] << "\n";
    }
//
//  Make a list of the node representatives.
//
    cout << "\n";
    cout << "      Rep   Coordinates:\n";
    cout << "\n";

    i4vec_sort_heap_a ( NODE_NUM, node_rep );

    rep_num = 0;

    for ( node = 0; node < NODE_NUM; node++ )
    {
      if ( 1 <= node )
      {
        if ( node_rep[node-1] == node_rep[node] )
        {
          continue;
        }
      }

      rep = node_rep[node];

      cout << "  " << setw(8)  << rep_num
           << "  " << setw(12) << node_xy[0+rep*2]
           << "  " << setw(12) << node_xy[1+rep*2] << "\n";

      rep_num = rep_num + 1;
    }
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TEST_NUM
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests NS_ADJ_COL_SET, NS_ADJ_COUNT and NS_ADJ_ROW_SET.
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
{
# define NODE_NUM 15
# define TRIANGLE_NUM 4
# define TRIANGLE_ORDER 6
# define VARIABLE_NUM 36

  int adj_col[VARIABLE_NUM+1];
  int adj_num;
  int *adj_row;
  string file_name = "ns_triangulation.eps";
  int node;
  int node_show;
  int node_p_variable[NODE_NUM] = {
    3, -1,  8, -1, 13, 
   -1, -1, -1, -1, 
   24, -1, 29, 
   -1, -1, 
   36 };
  int node_u_variable[NODE_NUM] = {
    1,  4,  6,  9, 11, 
   14, 16, 18, 20, 
   22, 25, 27, 
   30, 32, 
   34 };
  int node_v_variable[NODE_NUM] = {
    2,  5,  7, 10, 12, 
   15, 17, 19, 21, 
   23, 26, 28, 
   31, 33, 
   35 };
  double node_xy[2*NODE_NUM] = {
   0.0, 0.0, 
   0.0, 1.0, 
   0.0, 2.0, 
   0.0, 3.0, 
   0.0, 4.0, 
   1.0, 0.0, 
   1.0, 1.0, 
   1.0, 2.0, 
   1.0, 3.0, 
   2.0, 0.0, 
   2.0, 1.0, 
   2.0, 2.0, 
   3.0, 0.0, 
   3.0, 1.0, 
   4.0, 0.0 };
  int num;
  int r;
  int rhi;
  int rlo;
  int triangle_neighbor[3*TRIANGLE_NUM] = {
    -1,  2, -1, 
     3,  1,  4, 
     2, -1, -1, 
    -1, -1,  2 };
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1, 10,  3,  6,  7,  2, 
    12,  3, 10,  8,  7, 11, 
     3, 12,  5,  8,  9,  4, 
    10, 15, 12, 13, 14, 11 };
  int triangle_show;
  int variable;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For an order 3/order 6 Taylor Hood triangulation\n";
  cout << "  for Navier Stokes velocity and pressure,\n";
  cout << "  NS_ADJ_COUNT counts variable adjacencies\n";
  cout << "  and sets up the sparse compressed column\n";
  cout << "  column pointer array.\n";
  cout << "  NS_ADJ_COL_SET sets up the sparse compressed column\n";
  cout << "  COL vector.\n";
  cout << "  NS_ADJ_ROW_SET sets up the sparse compressed column\n";
  cout << "  ROW vector.\n";
//
//  Plot the example.
//
  node_show = 2;
  triangle_show = 2;

  triangulation_order6_plot ( file_name, NODE_NUM, node_xy, 
    TRIANGLE_NUM, triangle_node, node_show, triangle_show );
//
//  Get the count of the variable adjacencies.
//  We don't really need to make this call, since the next
//  call does the calculation as part of getting ADJ_COL.
//
  cout << "\n";
  cout << "  Number of variables is " << VARIABLE_NUM << "\n";

  adj_num = ns_adj_count ( NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
    triangle_neighbor, node_u_variable, node_v_variable, node_p_variable );

  cout << "\n";
  cout << "  As computed by NS_ADJ_COUNT,\n";
  cout << "  Number of variable adjacency entries is " << adj_num << "\n";
//
//  Get the count of the variable adjacencies and the COL vector.
//
  adj_num = ns_adj_col_set ( NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node,
    triangle_neighbor, node_u_variable, node_v_variable, node_p_variable, 
    adj_col );

  cout << "\n";
  cout << "  As computed by NS_ADJ_COL_SET,\n";
  cout << "  Number of variable adjacency entries is " << adj_num << "\n";

  cout << "\n";
  cout << "  Variable adjacency pointers:\n";
  cout << "\n";
  cout << "  Variable     First      Last    Number\n";
  cout << "\n";

  for ( variable = 0; variable < VARIABLE_NUM; variable++ )
  {
    num = adj_col[variable+1] - adj_col[variable];

    cout << "  " << setw(8) << variable+1
         << "  " << setw(8) << adj_col[variable]
         << "  " << setw(8) << adj_col[variable+1]-1
         << "  " << setw(8) << num << "\n";
  }
//
//  Get the variable adjacencies.
//
  adj_row = new int[adj_num];

  ns_adj_row_set ( NODE_NUM, TRIANGLE_NUM, VARIABLE_NUM, triangle_node, 
    triangle_neighbor, node_u_variable, node_v_variable, node_p_variable, 
    adj_num, adj_col, adj_row );
//
//  This is a huge array.  We only print out the beginning and end.
//
  cout << "\n";
  cout << "  Variable adjacency row entries:\n";
  cout << "  (Partial printout only)\n";
  cout << "\n";
  cout << "     Entry     Row       Col\n";
  cout << "\n";

  for ( variable = 0; variable < VARIABLE_NUM; variable++ )
  {
    rlo = adj_col[variable]-1;
    rhi = adj_col[variable+1]-2;

    if ( variable <= 2 || VARIABLE_NUM - 4 <= variable )
    {
      cout << "\n";

      for ( r = rlo; r <= rhi; r++ )
      {
        cout << "  " << setw(8) << r+1
             << "  " << setw(8) << adj_row[r]
             << "  " << setw(8) << variable+1 << "\n";
      }
    }

    if ( variable == 2 )
    {
      cout << "\n";
      cout << "  (SKIPPING MANY MANY ENTRIES...)\n";
      cout << "\n";
    }

  }

  delete [] adj_row;

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
# undef VARIABLE_NUM
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:  
//
//    TEST04 tests POINTS_DELAUNAY_NAIVE_2D.
//
//  Diagram:
//
//    !....3&11....
//    !............
//    !............
//    X..9.........
//    !.....5......
//    !...........6
//    !.4.2...10...
//    !.....8...12.
//    V............
//    !..7.........
//    !......1.....
//    !............
//    !............
//    !----V----X--
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
{
# define NODE_NUM 12
# define DIM_NUM 2

  int i;
  int triangle_num;
  int *triangle_node;
  double node_xy[DIM_NUM*NODE_NUM] = {
     7.0, 3.0,
     4.0,  7.0,
     5.0, 13.0,
     2.0,  7.0,
     6.0,  9.0,
    12.0, 10.0,
     3.0,  4.0,
     6.0,  6.0,
     3.0, 10.0,
     8.0,  7.0,
     5.0, 13.0,
    10.0,  6.0 };

  cout << "\n";
  cout << "TEST05\n";
  cout << "  POINTS_DELAUNAY_NAIVE_2D computes the Delaunay\n";
  cout << "  triangulation of a set of nodes.\n";

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_xy, "  The nodes:" );

  triangle_node = points_delaunay_naive_2d ( NODE_NUM, node_xy, &triangle_num );

  cout << "\n";
  cout << "  Number of triangles is TRIANGLE_NUM = " << triangle_num << "\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_node, 
    "  The Delaunay triangles:" );

  delete [] triangle_node;

  return;
# undef NODE_NUM
# undef DIM_NUM
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests POINTS_HULL_2D.
//
//  Diagram:
//
//    !....3.......
//    !............
//    !..9.........
//    !.....5......
//    !...........6
//    !.4.2...10...
//    !.....8......
//    !.........12.
//    !..7.........
//    !......1.....
//    !............
//    !............
//    !-----------
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
{
# define NODE_NUM 12
# define DIM_NUM 2

  int i;
  int j;
  int ival[NODE_NUM];
  int nval;
  double node_xy[DIM_NUM*NODE_NUM] = {
       7.0,  3.0, 
       4.0,  7.0, 
       5.0, 13.0, 
       2.0,  7.0, 
       6.0,  9.0, 
      12.0,  8.0, 
       3.0,  4.0, 
       6.0,  6.0, 
       3.0, 10.0, 
       8.0,  7.0, 
       5.0, 13.0, 
      10.0,  6.0 };

  cout << "\n";
  cout << "TEST06\n";
  cout << "  POINTS_HULL_2D computes the convex hull\n";
  cout << "  of a set of nodes.\n";

  r8mat_transpose_print ( DIM_NUM, NODE_NUM, node_xy, "  The nodes:" );

  points_hull_2d ( NODE_NUM, node_xy, &nval, ival );

  cout << "\n";
  cout << "  The convex hull is formed by connecting:\n";
  cout << "\n";
  for ( j = 0; j < nval; j++ )
  {
    cout << "  "
         << setw(3) << j << "  "
         << setw(3) << ival[j] << "  ";
    for ( i = 0; i < DIM_NUM; i++ )
    {
      cout << "  " << setw(14) << node_xy[i+(ival[j]-1)*DIM_NUM];
    }
    cout << "\n";    
  }

  cout << "\n";
  cout << "  The correct sequence of nodes is:\n";
  cout << "  4, 9, 3, 6, 12, 1, 7, (4).\n";

  return;
# undef NODE_NUM
# undef DIM_NUM
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests Q_MEASURE.
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
{
  int hole_num;
  int node_num;
  double *node_xy;
  double q_area;
  double q_ave;
  double q_max;
  double q_min;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Q_MEASURE returns the Q measure of\n";
  cout << "  quality of a triangulation.\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );
//
//  Allocate space.
//
  node_xy = new double[2*node_num];
  triangle_node = new int [triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];
//
//  Get the triangulation data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Compute the triangulation quality.
//
  q_measure ( node_num, node_xy, triangle_order, triangle_num,
    triangle_node, &q_min, &q_max, &q_ave, &q_area );

  cout << "\n";
  cout << "  Q_MIN  = " << q_min << "\n";
  cout << "  Q_MAX  = " << q_max << "\n";
  cout << "  Q_AVE  = " << q_ave << "\n";
  cout << "  Q_AREA = " << q_area << "\n";
//
//  Free the memory.
//
  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests R8TRIS2.
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
{
# define NODE_NUM 9

  int error;
  double node_xy[NODE_NUM*2] = {
       0.0, 0.0,
       0.0, 1.0,
       0.2, 0.5,
       0.3, 0.6,
       0.4, 0.5,
       0.6, 0.4,
       0.6, 0.5,
       1.0, 0.0,
       1.0, 1.0 };
  int triangle_node[2*NODE_NUM*3];
  int triangle_neighbor[2*NODE_NUM*3];
  int triangle_num;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  R8TRIS2 computes the Delaunay triangulation of a\n";
  cout << "  pointset in 2D.\n";
//
//  Set up the Delaunay triangulation.
//
  error = r8tris2 ( NODE_NUM, node_xy, &triangle_num, triangle_node,
    triangle_neighbor );

  if ( error == 0 ) 
  {
    cout << "\n";
    cout << "  R8TRIS2 computed the Delaunay triangulation with no\n";
    cout << "  errors detected.\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8TRIS2 detected an error condition of index " << error << "\n";
    return;
  }

  triangulation_order3_print ( NODE_NUM, triangle_num, node_xy,
    triangle_node, triangle_neighbor );

  return;
# undef NODE_NUM
}

//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE.
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
{
# define N 10

  int i;
  int j;
  double phy[2*N];
  double ref[2*N];
  double ref2[2*N];
  int seed;
  double t[2*3] = {
    1.0, 1.0, 
    3.0, 1.0, 
    2.0, 5.0 };

  seed = 123456789;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For an order 3 triangle,\n";
  cout << "  TRIANGLE_ORDER3_PHYSICAL_TO_REFERENCE\n";
  cout << "  maps a physical point to a reference point.\n";
  cout << "  TRIANGLE_ORDER3_REFERENCE_TO_PHYSICAL \n";
  cout << "  maps a reference point to a physical point.\n";
  cout << "\n";
  cout << "   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )\n";
  cout << "\n";

  triangle_reference_sample ( N, &seed, ref );

  triangle_order3_reference_to_physical ( t, N, ref, phy );

  triangle_order3_physical_to_reference ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(10) << ref[0+j*2]
         << "  " << setw(10) << ref[1+j*2]
         << "  "
         << "  " << setw(10) << phy[0+j*2]
         << "  " << setw(10) << phy[1+j*2]
         << "  "
         << "  " << setw(10) << ref2[0+j*2]
         << "  " << setw(10) << ref2[1+j*2] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE.
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
{
# define N 10

  int i;
  int j;
  double phy[2*N];
  double ref[2*N];
  double ref2[2*N];
  int seed;
  double t[2*6] = {
    7.0, 2.0, 
    9.0, 2.0, 
    7.0, 3.0, 
    8.0, 2.0, 
    8.0, 2.5, 
    7.0, 2.5 };

  seed = 123456789;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For an order 6 triangle,\n";
  cout << "  TRIANGLE_ORDER6_PHYSICAL_TO_REFERENCE\n";
  cout << "  maps a physical point to a reference point\n";
  cout << "  TRIANGLE_ORDER6_REFERENCE_TO_PHYSICAL\n";
  cout << "  maps a reference point to a physical point.\n";
  cout << "\n";
  cout << "   (  XSI    ETA ) ==> ( X      Y  ) ==> ( XSI2    ETA2 )\n";
  cout << "\n";

  triangle_reference_sample ( N, &seed, ref );

  triangle_order6_reference_to_physical ( t, N, ref, phy );

  triangle_order6_physical_to_reference ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(10) << ref[0+j*2]
         << "  " << setw(10) << ref[1+j*2]
         << "  "
         << "  " << setw(10) << phy[0+j*2]
         << "  " << setw(10) << phy[1+j*2]
         << "  "
         << "  " << setw(10) << ref2[0+j*2]
         << "  " << setw(10) << ref2[1+j*2] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests TRIANGULATION_NODE_ORDER.
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
{
# define NODE_NUM 36
# define TRIANGLE_NUM 41
# define TRIANGLE_ORDER 3

  int *node_order;
  int triangle_node[3*TRIANGLE_NUM] = {
     1,  8,  7, 
     1,  2,  8, 
     2,  9,  8, 
     2,  3,  9, 
     3, 10,  9, 
     3,  4, 10, 
     4, 11, 10, 
     4,  5, 11, 
     5, 12, 11, 
     5,  6, 12, 
     7, 14, 13, 
     7,  8, 14, 
     8, 15, 14, 
     8,  9, 15, 
    11, 18, 17, 
    11, 12, 18, 
    13, 20, 19, 
    13, 14, 20, 
    14, 21, 20, 
    14, 15, 21, 
    15, 22, 21, 
    15, 16, 22, 
    16, 23, 22, 
    16, 17, 23, 
    17, 24, 23, 
    17, 18, 24, 
    19, 26, 25, 
    19, 20, 26, 
    21, 28, 27, 
    21, 22, 28, 
    25, 30, 29, 
    25, 26, 30, 
    26, 31, 30, 
    27, 32, 31, 
    27, 28, 32, 
    29, 34, 33, 
    29, 30, 34, 
    30, 35, 34, 
    30, 31, 35, 
    31, 36, 35, 
    31, 32, 36 };

  cout << "\n";
  cout << "TEST11\n";
  cout << "  TRIANGULATION_NODE_ORDER computes the order\n";
  cout << "  of the nodes in a triangulation.\n";

  node_order = triangulation_node_order ( TRIANGLE_ORDER, TRIANGLE_NUM,
    triangle_node, NODE_NUM );

  i4vec_print ( NODE_NUM, node_order, "  NODE ORDER:" );

  delete [] node_order;

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests TRIANGULATION_ORDER3_ADJ_SET.
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
{
  int *adj;
  int adj_num;
  int *adj_col;
  int hole_num;
  int k;
  int node;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_order = 3;
  int triangle_num;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies\n";
  cout << "  TRIANGULATION_ORDER3_ADJ_SET sets adjacencies.\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  adj_col = new int[node_num+1];
  node_xy = new double[2*node_num];
  triangle_neighbor = new int[3*triangle_num];
  triangle_node = new int[triangle_order*triangle_num];
//
//  Get the example data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Get the count of the adjacencies.
//
  adj_num = triangulation_order3_adj_count ( node_num, triangle_num, 
    triangle_node, triangle_neighbor, adj_col );

  cout << "\n";
  cout << "  Number of adjacency entries is " << adj_num << "\n";

  cout << "\n";
  cout << "  Adjacency pointers:\n";
  cout << "\n";
  for ( node = 1; node <= node_num; node++ )
  {
    cout << "  " << setw(8) << node
         << "  " << setw(8) << adj_col[node-1]
         << "  " << setw(8) << adj_col[node]-1 << "\n";
  }
//
//  Get the adjacencies.
//
  adj = triangulation_order3_adj_set ( node_num, triangle_num, triangle_node,
    triangle_neighbor, adj_num, adj_col );
//
//  Print the adjacencies.
//
  for ( node = 1; node <= node_num; node++ )
  {
    cout << "\n";
    cout << "  Nodes adjacent to node " << node << "\n";
    cout << "\n";

    for ( k = adj_col[node-1]; k <= adj_col[node]-1; k++ )
    {
      cout << "  " << setw(8) << adj[k-1] << "\n";
    }
  }

  delete [] adj;
  delete [] adj_col;
  delete [] node_xy;
  delete [] triangle_neighbor;
  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST125 tests TRIANGULATION_ORDER3_ADJ_SET2.
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
{
  int adj;
  int adj_num;
  int *adj_col;
  int hole_num;
  int *ia;
  int *ja;
  int k;
  int node;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_order = 3;
  int triangle_num;

  cout << "\n";
  cout << "TEST125\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies\n";
  cout << "  TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies\n";
  cout << "  as a pair of vectors IA(*), JA(*).\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  adj_col = new int[node_num+1];
  node_xy = new double[2*node_num];
  triangle_neighbor = new int[3*triangle_num];
  triangle_node = new int[triangle_order*triangle_num];
//
//  Get the example data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Get the count of the adjacencies.
//
  adj_num = triangulation_order3_adj_count ( node_num, triangle_num, 
    triangle_node, triangle_neighbor, adj_col );

  cout << "\n";
  cout << "  Number of adjacency entries is " << adj_num << "\n";

  cout << "\n";
  cout << "  Adjacency pointers:\n";
  cout << "\n";
  for ( node = 1; node <= node_num; node++ )
  {
    cout << "  " << setw(8) << node
         << "  " << setw(8) << adj_col[node-1]
         << "  " << setw(8) << adj_col[node]-1 << "\n";
  }
//
//  Get the adjacencies.
//
  ia = new int[adj_num];
  ja = new int[adj_num];

  triangulation_order3_adj_set2 ( node_num, triangle_num, triangle_node,
    triangle_neighbor, adj_num, adj_col, ia, ja );
//
//  Print the adjacencies.
//
  cout << "\n";
  cout << "  Adjacency list:\n";
  cout << "\n";
  for ( adj = 0; adj < adj_num; adj++ )
  {
    cout << "  " << setw(8) << adj+1
         << "  (" << setw(2) << ia[adj]
         << "," << setw(2) << ja[adj] << ")\n";
  }

  delete [] adj_col;
  delete [] ia;
  delete [] ja;
  delete [] node_xy;
  delete [] triangle_neighbor;
  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test127 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST127 tests TRIANGULATION_ORDER3_ADJACENCY.
//
//  Discussion:
//
//    41--42--43--44  45--46--47--48
//     | \ | \ | \ |   | \ | \ | \ |
//    33--34--35--36  37--38--39--40
//     | \ |                   | \ |
//    29--30                  31--32
//     | \ |                   | \ |
//    25--26                  27--28
//     | \ |                   | \ |
//    21--22                  23--24
//     | \ |                   | \ |
//    17--18                  19--20
//     | \ |                   | \ |
//     9--10--11--12--13--14--15--16
//     | \ | \ | \ | \ | \ | \ | \ |
//     1---2---3---4---5---6---7---8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_NUM 46
# define NODE_NUM 48
  
  int *adj;
  int element_node[3*ELEMENT_NUM] = {
     1,  2,  9,
     2, 10,  9,
     2,  3, 10,
     3, 11, 10,
     3,  4, 11,
     4, 12, 11,
     4,  5, 12,
     5, 13, 12,
     5,  6, 13,
     6, 14, 13,
     6,  7, 14,
     7, 15, 14,
     7,  8, 15,
     8, 16, 15,
     9, 10, 17,
    10, 18, 17,
    15, 16, 19,
    16, 20, 19,
    17, 18, 21,
    18, 22, 21,
    19, 20, 23,
    20, 24, 23,
    21, 22, 25,
    22, 26, 25,
    23, 24, 27,
    24, 28, 27,
    25, 26, 29,
    26, 30, 29,
    27, 28, 31,
    28, 32, 31,
    29, 30, 33,
    30, 34, 33,
    31, 32, 39,
    32, 40, 39,
    33, 34, 41,
    34, 42, 41,
    34, 35, 42,
    35, 43, 42,
    35, 36, 43,
    36, 44, 43,
    37, 38, 45,
    38, 46, 45,
    38, 39, 46,
    39, 47, 46,
    39, 40, 47,
    40, 48, 47 };
  int element_num = ELEMENT_NUM;
  int i;
  int j;
  int node_num = NODE_NUM;

  cout << "\n";
  cout << "TEST127\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_ADJACENCY sets the full\n";
  cout << "  adjacency matrix.\n";

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      element_node[i+j*3] = element_node[i+j*3] - 1;
    }
  }

  adj = triangulation_order3_adjacency ( node_num, element_num, element_node );

  cout << "\n";
  cout << "  Adjacency matrix:\n";
  cout << "\n";
  cout << "                1         2         3         4       \n";
  cout << "      012345678901234567890123456789012345678901234567\n";
  cout << "\n";
  for ( i = 0; i < node_num; i++ )
  {
    cout << "  " << setw(2) << i << "  ";
    for ( j = 0; j < node_num; j++ )
    {
      cout << adj[i+j*node_num];
    }
    cout << "\n";
  }

  delete [] adj;

  return;
# undef ELEMENT_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT.
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
{
# define DIM_NUM 2
# define NODE_NUM 36
# define TRIANGLE_NUM 41
# define TRIANGLE_ORDER 3

  int boundary_edge_num;
  string file_name = "triangulation_order3_plot2.eps";
  int node_show = 2;
  double node_xy[DIM_NUM*NODE_NUM] = {
    0.0, 0.0, 
    1.0, 0.0, 
    2.0, 0.0, 
    3.0, 0.0, 
    4.0, 0.0, 
    5.0, 0.0, 
    0.0, 1.0, 
    1.0, 1.0, 
    2.0, 1.0, 
    3.0, 1.0, 
    4.0, 1.0, 
    5.0, 1.0, 
    0.0, 2.0, 
    1.0, 2.0, 
    2.0, 2.0, 
    3.0, 2.0, 
    4.0, 2.0, 
    5.0, 2.0, 
    0.0, 3.0, 
    1.0, 3.0, 
    2.0, 3.0, 
    3.0, 3.0, 
    4.0, 3.0, 
    5.0, 3.0, 
    0.0, 4.0, 
    1.0, 4.0, 
    2.0, 4.0, 
    3.0, 4.0, 
    0.0, 5.0, 
    1.0, 5.0, 
    2.0, 5.0, 
    3.0, 5.0, 
    0.0, 6.0, 
    1.0, 6.0, 
    2.0, 6.0, 
    3.0, 6.0 };
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  8,  7, 
     1,  2,  8, 
     2,  9,  8, 
     2,  3,  9, 
     3, 10,  9, 
     3,  4, 10, 
     4, 11, 10, 
     4,  5, 11, 
     5, 12, 11, 
     5,  6, 12, 
     7, 14, 13, 
     7,  8, 14, 
     8, 15, 14, 
     8,  9, 15, 
    11, 18, 17, 
    11, 12, 18, 
    13, 20, 19, 
    13, 14, 20, 
    14, 21, 20, 
    14, 15, 21, 
    15, 22, 21, 
    15, 16, 22, 
    16, 23, 22, 
    16, 17, 23, 
    17, 24, 23, 
    17, 18, 24, 
    19, 26, 25, 
    19, 20, 26, 
    21, 28, 27, 
    21, 22, 28, 
    25, 30, 29, 
    25, 26, 30, 
    26, 31, 30, 
    27, 32, 31, 
    27, 28, 32, 
    29, 34, 33, 
    29, 30, 34, 
    30, 35, 34, 
    30, 31, 35, 
    31, 36, 35, 
    31, 32, 36 };
  int triangle_show = 2;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the\n";
  cout << "  boundary edges;\n";
  cout << "  TRIANGULATION_ORDER3_PLOT plots the triangulation.\n";

  triangulation_order3_plot ( file_name, NODE_NUM, node_xy, 
    TRIANGLE_NUM, triangle_node, node_show, triangle_show );

  cout << "\n";
  cout << "  An Encapsulated PostScript image of this\n";
  cout << "  triangulation is in \"" << file_name << "\".\n";

  boundary_edge_num = triangulation_order3_boundary_edge_count ( 
    TRIANGLE_NUM, triangle_node );

  cout << "\n";
  cout << "  Number of boundary edges = " << boundary_edge_num << "\n";
  cout << "  Correct number =           " << 33 << "\n";

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER.
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
{
  int boundary_num;
  int hole_num = 2;
  int node_num = 36;
  int triangle_num = 41;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER\n";
  cout << "  determines the number of edges that lie on the\n";
  cout << "  boundary of a region that has been triangulated.\n";
  cout << "\n";
  cout << "  Number of points =         " << node_num << "\n";
  cout << "  Number of triangles =      " << triangle_num << "\n";
  cout << "  Number of holes =          " << hole_num << "\n";

  boundary_num = triangulation_order3_boundary_edge_count_euler ( node_num, 
    triangle_num, hole_num );

  cout << "  Number of boundary edges = " << boundary_num << "\n";

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests TRIANGULATION_ORDER3_BOUNDARY_NODE.
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
{
# define DIM_NUM 2
# define NODE_NUM 36
# define TRIANGLE_NUM 41
# define TRIANGLE_ORDER 3

  int i;
  bool *node_boundary;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
     1,  8,  7, 
     1,  2,  8, 
     2,  9,  8, 
     2,  3,  9, 
     3, 10,  9, 
     3,  4, 10, 
     4, 11, 10, 
     4,  5, 11, 
     5, 12, 11, 
     5,  6, 12, 
     7, 14, 13, 
     7,  8, 14, 
     8, 15, 14, 
     8,  9, 15, 
    11, 18, 17, 
    11, 12, 18, 
    13, 20, 19, 
    13, 14, 20, 
    14, 21, 20, 
    14, 15, 21, 
    15, 22, 21, 
    15, 16, 22, 
    16, 23, 22, 
    16, 17, 23, 
    17, 24, 23,
    17, 18, 24, 
    19, 26, 25, 
    19, 20, 26, 
    21, 28, 27, 
    21, 22, 28, 
    25, 30, 29, 
    25, 26, 30, 
    26, 31, 30, 
    27, 32, 31, 
    27, 28, 32, 
    29, 34, 33, 
    29, 30, 34, 
    30, 35, 34, 
    30, 31, 35, 
    31, 36, 35, 
    31, 32, 36 };

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_BOUNDARY_NODE determines which\n";
  cout << "  nodes lie on the boundary of a triangulation.\n";

  node_boundary = triangulation_order3_boundary_node ( NODE_NUM, TRIANGLE_NUM, 
    triangle_node );

  lvec_print ( NODE_NUM, node_boundary, "    Node  BN?" );

  delete [] node_boundary;

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests TRIANGULATION_ORDER3_CHECK.
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
{
# define NODE_NUM 13
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int i;
  int ierror;
  int isave;
  int node_num2;
  int triangle_num2;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
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

  cout << "\n";
  cout << "TEST16\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_CHECK checks the triangulation.\n";

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node, 
    "  Triangles:" );
//
//  Pass all tests.
//
  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  cout << "  Error code = " << ierror << "\n";
//
//  Fail test 1.
//
  node_num2 = 2;

  ierror = triangulation_order3_check ( node_num2, TRIANGLE_NUM, 
    triangle_node );

  cout << "  Error code = " << ierror << "\n";
//
//  Fail test 2.
//
  triangle_num2 = 0;

  ierror = triangulation_order3_check ( NODE_NUM, triangle_num2, 
    triangle_node );

  cout << "  Error code = " << ierror << "\n";
//
//  Fail test 3.
//
  isave = triangle_node[1+4*3];
  triangle_node[1+4*3] = 0;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  cout << "  Error code = " << ierror << "\n";
  triangle_node[1+4*3] = isave;
//
//  Fail test 4.
//
  isave = triangle_node[2+9*3];
  triangle_node[2+9*3] = 2 * NODE_NUM + 1;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  cout << "  Error code = " << ierror << "\n";
  triangle_node[2+9*3] = isave;
//
//  Fail test 5.
//
  triangle_node[2+3*3] = 3;
  triangle_node[2+7*3] = 3;
  triangle_node[2+9*3] = 3;
  triangle_node[2+10*3] = 3;
  triangle_node[1+13*3] = 3;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  cout << "  Error code = " << ierror << "\n";

  triangle_node[2+3*3] = 5;
  triangle_node[2+7*3] = 5;
  triangle_node[2+9*3] = 5;
  triangle_node[2+10*3] = 5;
  triangle_node[1+13*3] = 5;
//
//  Fail test 6.
//
  triangle_node[0+8*3] = 7;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  cout << "  Error code = " << ierror << "\n";
  triangle_node[0+8*3] = 9;
//
//  Fail test 7.
//
  triangle_node[2+6*3] = 2;

  ierror = triangulation_order3_check ( NODE_NUM, TRIANGLE_NUM, triangle_node );

  cout << "  Error code = " << ierror << "\n";

  triangle_node[2+6*3] = 9;

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests TRIANGULATION_ORDER3_EXAMPLE1.
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
{
  int hole_num;
  int node_num;
  double *node_xy;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 3;
  int triangle_num;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_EXAMPLE1_SIZE gives the sizes\n";
  cout << "  for an example triangulation;\n";
  cout << "  TRIANGULATION_ORDER3_EXAMPLE1 returns the information\n";
  cout << "  for an example triangulation;\n";
  cout << "  TRIANGULATION_ORDER3_PRINT prints a triangulation.\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = new double[2*node_num];
  triangle_node = new int[triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];
//
//  Get the data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  triangulation_order3_print ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests TRIANGULATION_ORDER3_NEIGHBOR.
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
{
# define NODE_NUM 13
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int s1;
  int s2;
  int t1;
  int t2;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
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

  cout << "\n";
  cout << "TEST18\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_NEIGHBOR determines the\n";
  cout << "  triangle neighbors.\n";
  cout << "\n";
  cout << "    T1    S1    T2    S2\n";
  cout << "\n";

  for ( t1 = 1; t1 <= TRIANGLE_NUM; t1++ )
  {
    for ( s1 = 1; s1 <= 3; s1++ )
    {
      triangulation_order3_neighbor ( TRIANGLE_NUM, triangle_node, 
        t1, s1, &t2, &s2 );

      cout << "  " << setw(4) << t1
           << "  " << setw(4) << s1
           << "  " << setw(4) << t2
           << "  " << setw(4) << s2 << "\n";
    }
  }

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests TRIANGULATION_NEIGHBOR_ELEMENTS.
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
{
# define TRIANGLE_NUM 16
# define TRIANGLE_ORDER 3

  int i;
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
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
  int *triangle_neighbor;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  TRIANGULATION_NEIGHBOR_ELEMENTS determines the\n";
  cout << "  adjacency relationships between elements.\n";

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM, triangle_node, 
    "  Elements:" );

  triangle_neighbor = triangulation_neighbor_elements ( TRIANGLE_ORDER, 
    TRIANGLE_NUM, triangle_node );

  i4mat_transpose_print ( 3, TRIANGLE_NUM, triangle_neighbor, 
    "  Element neighbors:" );

  delete [] triangle_neighbor;

  return;
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests TRIANGULATION_ORDER3_PLOT.
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
{
  string file_name = "triangulation_order3_plot.eps";
  int hole_num;
  int node_show = 0;
  int node_num;
  double *node_xy;
  int triangle_show = 2;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_PLOT can plot a triangulation.\n";
//
//  Get the sizes.
//
  triangulation_order3_example1_size ( &node_num, &triangle_num, &hole_num );

  cout << "\n";
  cout << "  Example data has " << node_num << " points,\n";
  cout << "  organized into " << triangle_num << " triangles.\n";
//
//  Allocate space.
//
  node_xy = new double[2*node_num];
  triangle_node = new int[triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];
//
//  Get the example data.
//
  triangulation_order3_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Make the plot.
//
  triangulation_order3_plot ( file_name, node_num, node_xy, triangle_num, 
    triangle_node, node_show, triangle_show );

  cout << "\n";
  cout << "  TRIANGULATION_ORDER3_PLOT has created an\n";
  cout << "  Encapsulated PostScript file (EPS) containing\n";
  cout << "  an image of the triangulation.\n";
  cout << "\n";
  cout << "  This file is called \"" << file_name << "\".\n";

  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests TRIANGULATION_ORDER3_PRINT.
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
{
# define NODE_NUM 9
# define TRIANGLE_NUM 12
# define TRIANGLE_ORDER 3

  double node_xy[2*NODE_NUM] = {
       0.0, 0.0, 
       0.0, 1.0, 
       0.2, 0.5,
       0.3, 0.6, 
       0.4, 0.5, 
       0.6, 0.4, 
       0.6, 0.5, 
       1.0, 0.0, 
       1.0, 1.0 };
  int triangle_node[TRIANGLE_ORDER*TRIANGLE_NUM] = {
       2, 1, 3, 
       3, 1, 6, 
       2, 3, 4, 
       4, 3, 5, 
       7, 4, 5, 
       5, 3, 6, 
       7, 5, 6, 
       9, 4, 7, 
       6, 1, 8, 
       7, 6, 8, 
       7, 8, 9, 
       2, 4, 9 };
  int triangle_neighbor[3*TRIANGLE_NUM] = {
       -28,   2,  3, 
         1,   9,  6, 
         1,   4, 12, 
         3,   6,  5, 
         8,   4,  7, 
         4,   2,  7, 
         5,   6, 10, 
        12,   5, 11, 
         2, -34, 10, 
         7,   9, 11, 
        10, -38,  8, 
         3,   8, -3 };

  cout << "\n";
  cout << "TEST21\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_PRINT prints out a triangulation.\n";

  triangulation_order3_print ( NODE_NUM, TRIANGLE_NUM, node_xy,
    triangle_node, triangle_neighbor );

  return;
# undef NODE_NUM
# undef TRIANGLE_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test213 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST213 tests TRIANGULATION_ORDER3_QUAD.
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
{
# define QUAD_NUM 6

  int i;
  int j;
  int k;
  int n;
  int n11;
  int n12;
  int n21;
  int n22;
  int node_num;
  double *node_xy;
  double quad_value;
  double quad_w[QUAD_NUM] = {
    0.1666666666666666, 
    0.1666666666666666, 
    0.1666666666666666, 
    0.1666666666666666, 
    0.1666666666666666, 
    0.16666666666666660 };
  double quad_xy[2*QUAD_NUM] = {
    0.659027622374092,  0.231933368553031, 
    0.659027622374092,  0.109039009072877, 
    0.231933368553031,  0.659027622374092, 
    0.231933368553031,  0.109039009072877, 
    0.109039009072877,  0.659027622374092, 
    0.109039009072877,  0.231933368553031 };
  double region_area;
  int test;
  int *triangle_node;
  int test_num = 4;
  int triangle_order = 3;
  int triangle_num;
  double x;
  double y;

  cout << "\n";
  cout << "TEST213\n";
  cout << "  TRIANGULATION_ORDER3_QUAD can apply a quadrature rule\n";
  cout << "  to every triangle in a triangulated region,\n";
  cout << "  and estimate the integral of a function over\n";
  cout << "  that region.\n";
  cout << "\n";
  cout << "  NODE_NUM   TRI_NUM  Integral estim  Area of Region\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Set up the grid.
//
    n = i4_power ( 2, test - 1 );
    node_num = ( n + 1 ) * ( n + 1 );

    node_xy = new double[2*node_num];

    k = 0;
    for ( j = 1; j <= n + 1; j++ )
    {
      y = ( double ) ( j - 1 ) / ( double ) ( n + 1 - 1 );
      for ( i = 1; i <= n + 1; i++ )
      {
        x = ( double ) ( i - 1 ) / ( double ) ( n + 1 - 1 );
        node_xy[0+k*2] = x;
        node_xy[1+k*2] = y;
        k = k + 1;
      }
    }
//
//  Set up the triangulation.
//
    triangle_num = 2 * n * n;

    triangle_node = new int[triangle_order*triangle_num];

    k = 0;
    for ( j = 1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        n11 = i     + ( j     - 1 ) * ( n + 1 );
        n12 = i     + ( j + 1 - 1 ) * ( n + 1 );
        n21 = i + 1 + ( j     - 1 ) * ( n + 1 );
        n22 = i + 1 + ( j + 1 - 1 ) * ( n + 1 );

        triangle_node[0+k*3] = n11;
        triangle_node[1+k*3] = n21;
        triangle_node[2+k*3] = n12;
        k = k + 1;

        triangle_node[0+k*3] = n22;
        triangle_node[1+k*3] = n12;
        triangle_node[2+k*3] = n21;
        k = k + 1;
      }
    }
//
//  Estimate the integral.
//
    triangulation_order3_quad ( node_num, node_xy, triangle_order, 
      triangle_num, triangle_node, &quad_fun, QUAD_NUM, quad_xy, quad_w, 
      &quad_value, &region_area );

    cout << "  " << setw(8)  << node_num
         << "  " << setw(8)  << triangle_num
         << "  " << setw(14) << quad_value
         << "  " << setw(14) << region_area << "\n";
//
//  Delete allocatables.
//
    delete [] node_xy;
    delete [] triangle_node;
  }

  return;;
# undef QUAD_NUM
}
//****************************************************************************80

void quad_fun ( int n, double xy_vec[], double f_vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_FUN is a sample integrand function for TRIANGULATION_QUAD.
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
//    Input, int N, the number of evaluation points.
//
//    Input, double XY_VEC[2*N], the evaluation points.
//
//    Output, double F_VEC[N], the value of the integrand
//    function at the evaluation points.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    f_vec[i] = exp ( pow ( xy_vec[0+i*2], 2 ) 
                   + pow ( xy_vec[1+i*2], 2 ) );
  }
  return;
}
//****************************************************************************80

void test215 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST215 tests TRIANGULATION_ORDER3_REFINE_COMPUTE.
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
{
# define DIM_NUM 2
# define NODE_NUM1 5
# define TRIANGLE_NUM1 3
# define TRIANGLE_ORDER 3

  int *edge_data;
  int node_num2;
  double node_xy1[DIM_NUM*NODE_NUM1] = {
       0.0, 0.0, 
       1.0, 0.0, 
       0.0, 1.0, 
       1.0, 1.0, 
       0.5, 1.5 };
  double *node_xy2;
  int triangle_node1[TRIANGLE_ORDER*TRIANGLE_NUM1] = {
       1, 2, 3, 
       4, 3, 2, 
       3, 4, 5 };
  int *triangle_node2;
  int triangle_num2;

  cout << "\n";
  cout << "TEST215\n";
  cout << "  For an order3 triangulation:\n";
  cout << "  TRIANGULATION_ORDER3_REFINE_SIZE determines the\n";
  cout << "  size of a refined triangulation.\n";
  cout << "  TRIANGULATION_ORDER3_REFINE_COMPUTES computes the\n";
  cout << "  refined triangulation.\n";

  cout << "\n";
  cout << "  The number of nodes is " << NODE_NUM1 << "\n";
  cout << "  The number of triangles is " << TRIANGLE_NUM1 << "\n";

  r8mat_transpose_print ( DIM_NUM, NODE_NUM1, node_xy1, 
    "  The nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM1, triangle_node1, 
    "  The triangles:" );

  edge_data = new int[5*(3*TRIANGLE_NUM1)];

  cout << "\n";
  cout << "  Sizing the refined mesh:\n";

  triangulation_order3_refine_size ( NODE_NUM1, TRIANGLE_NUM1, 
    triangle_node1, &node_num2, &triangle_num2, edge_data );

  cout << "\n";
  cout << "  Information about the refined mesh:\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num2 << "\n";
  cout << "  The number of triangles is " << triangle_num2 << "\n";

  cout << "\n";
  cout << "  Computing the refined mesh:\n";

  node_xy2 = new double[DIM_NUM*node_num2];
  triangle_node2 = new int[TRIANGLE_ORDER*triangle_num2];

  triangulation_order3_refine_compute ( NODE_NUM1, TRIANGLE_NUM1, 
    node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, 
    triangle_node2 );

  r8mat_transpose_print ( DIM_NUM, node_num2, node_xy2, 
    "  The refined nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, triangle_num2, triangle_node2, 
    "  The refined triangles:" );

  delete [] edge_data;
  delete [] node_xy2;
  delete [] triangle_node2;

  return;
# undef DIM_NUM
# undef NODE_NUM1
# undef TRIANGLE_NUM1
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test217 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST217 tests TRIANGULATION_SEARCH_DELAUNAY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define NODE_NUM 13
# define TEST_NUM 10
# define TRIANGLE_ORDER 3

  double alpha;
  double beta;
  double d1;
  double d2;
  double d3;
  double dist;
  double dnear;
  int edge;
  int error;
  double gamma;
  int i;
  int i1;
  int i2;
  int i3;
  int nnear;
  double node_xy[DIM_NUM*NODE_NUM] = {
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
  double p[DIM_NUM];
  int seed;
  int step_num;
  int td[TEST_NUM];
  int test;
  int triangle_index;
  int triangle_neighbor[3*2*NODE_NUM];
  int triangle_node[TRIANGLE_ORDER*2*NODE_NUM];
  int triangle_num;
  double xd[DIM_NUM*TEST_NUM];

  cout << "\n";
  cout << "TEST217\n";
  cout << "  Given a set of nodes NODE_XY, and a single point XD,\n";
  cout << "  find the nearest node in NODE_XY to XD.\n";
  cout << "\n";
  cout << "  POINTS_POINT_NEAR_NAIVE_ND uses a naive method.\n";
  cout << "  TRIANGULATION_SEARCH_DELAUNAY finds a triangle\n";
  cout << "  containing the point.  Often, one of these vertices\n";
  cout << "  is the closest point.\n";
//
//  Set up the Delaunay triangulation.
//
  error = r8tris2 ( NODE_NUM, node_xy, &triangle_num, triangle_node,
    triangle_neighbor );

  if ( error == 0 )
  {
    cout << "\n";
    cout << "  R8TRIS2 computed the Delaunay triangulation.\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8TRIS2 returned an error condition.\n";
    exit ( 1 );
  }
//
//  Get the test points.
//
  seed = 123456789;

  triangulation_order3_sample ( NODE_NUM, node_xy, triangle_num,
    triangle_node, TEST_NUM, &seed, xd, td );

  cout << "\n";
  cout << "              X         Y     Distance    Index     Steps\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    for ( i = 0; i < DIM_NUM; i++ )
    {
      p[i] = xd[i+test*DIM_NUM];
    }

    nnear = points_point_near_naive_nd ( DIM_NUM, NODE_NUM, node_xy, 
      p, &dnear );

    cout << "\n";
    cout << "  XD       " 
         << setw(8) << p[0] << "  "
         << setw(8) << p[1] << "\n";
    cout << "  Naive    " 
         << setw(8) << node_xy[0+(nnear-1)*DIM_NUM] << "  "
         << setw(8) << node_xy[1+(nnear-1)*DIM_NUM] << "  "
         << setw(8) << dnear << "  "
         << setw(6) << nnear << "\n";

    triangulation_search_delaunay ( NODE_NUM, node_xy, TRIANGLE_ORDER, triangle_num,
      triangle_node, triangle_neighbor, p, &triangle_index, 
      &alpha, &beta, &gamma, &edge, &step_num );
   
    if ( triangle_index < 1 )
    {
      cout << "\n";
      cout << "  Error: the search failed.\n";
      continue;
    }

    i1 = triangle_node[0+(triangle_index-1)*TRIANGLE_ORDER];
    d1 = sqrt ( pow ( p[0] - node_xy[0+i1*2], 2 ) 
              + pow ( p[1] - node_xy[1+i1*2], 2 ) );

    dist = d1;
    nnear = i1;

    i2 = triangle_node[1+(triangle_index-1)*TRIANGLE_ORDER];
    d2 = sqrt ( pow ( p[0] - node_xy[0+i2*2], 2 ) 
              + pow ( p[1] - node_xy[1+i2*2], 2 ) );

    if ( d2 < dist )
    {
      dnear = d2;
      nnear = i2;
    }

    i3 = triangle_node[2+(triangle_index-1)*TRIANGLE_ORDER];
    d3 = sqrt ( pow ( p[0] - node_xy[0+i3*2], 2 ) 
              + pow ( p[1] - node_xy[1+i3*2], 2 ) );

    if ( d3 < dist )
    {
      dnear = d3;
      nnear = i3;
    }

    cout << "  Delaunay "
         << setw(8) << node_xy[0+nnear*2] << "  "
         << setw(8) << node_xy[1+nnear*2] << "  "
         << setw(8) << dnear << "  "
         << setw(8) << nnear + 1 
         << setw(8) << step_num << "\n";
  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TEST_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test219 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST219 tests TRIANGULATION_SEARCH_DELAUNAY, TRIANGULATION_SEARCH_NAIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define NODE_NUM 13
# define TEST_NUM 10
# define TRIANGLE_ORDER 3

  double alpha;
  double beta;
  int edge;
  int error;
  double gamma;
  int i;
  int nnear;
  double node_xy[DIM_NUM*NODE_NUM] = {
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
  double p_test[DIM_NUM*TEST_NUM];
  int seed = 123456789;
  int step_num;
  int t_test[TEST_NUM];
  int test;
  int triangle_index1;
  int triangle_index2;
  int triangle_neighbor[3*2*NODE_NUM];
  int triangle_node[TRIANGLE_ORDER*2*NODE_NUM];
  int triangle_num;

  cout << "\n";
  cout << "TEST219\n";
  cout << "  Given a triangulation, and a point P,\n";
  cout << "  find the triangle T containing to P.\n";
  cout << "\n";
  cout << "  TRIANGULATION_SEARCH_NAIVE uses a naive method.\n";
  cout << "  TRIANGULATION_SEARCH_DELAUNAY uses a method that will work\n";
  cout << "  fast if the triangulation is Delaunay.\n";
//
//  Set up the Delaunay triangulation.
//
  error = r8tris2 ( NODE_NUM, node_xy, &triangle_num, triangle_node,
    triangle_neighbor );

  if ( error == 0 )
  {
    cout << "\n";
    cout << "  R8TRIS2 computed the Delaunay triangulation.\n";
  }
  else
  {
    cout << "\n";
    cout << "  R8TRIS2 returned an error condition.\n";
    exit ( 1 );
  }
//
//  Get the test points.
//
  triangulation_order3_sample ( NODE_NUM, node_xy, triangle_num,
    triangle_node, TEST_NUM, &seed, p_test, t_test );

  cout << "\n";
  cout << "         X           Y     Naive   Delaunay  Steps\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    triangle_index1 = triangulation_search_naive ( NODE_NUM, node_xy, 
      TRIANGLE_ORDER, triangle_num, triangle_node, p_test+DIM_NUM*test );

    triangulation_search_delaunay ( NODE_NUM, node_xy, TRIANGLE_ORDER, 
      triangle_num, triangle_node, triangle_neighbor, p_test+DIM_NUM*test, 
      &triangle_index2, &alpha, &beta, &gamma, &edge, &step_num );

    cout << "  " << setw(10) << p_test[0+test*DIM_NUM]
         << "  " << setw(10) << p_test[1+test*DIM_NUM]
         << "  " << setw(8) << triangle_index1
         << "  " << setw(8) << triangle_index2
         << "  " << setw(8) << step_num << "\n";

  }

  return;
# undef DIM_NUM
# undef NODE_NUM
# undef TEST_NUM
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests TRIANGULATION_ORDER6_ADJ_SET.
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
{
  int *adj;
  int adj_num;
  int *adj_col;
  int hole_num;
  int k;
  int node;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_order = 6;
  int triangle_num;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies\n";
  cout << "  TRIANGULATION_ORDER6_ADJ_SET sets adjacencies.\n";
//
//  Get the sizes.
//
  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  adj_col = new int[node_num+1];
  node_xy = new double[2*node_num];
  triangle_neighbor = new int[3*triangle_num];
  triangle_node = new int[triangle_order*triangle_num];
//
//  Get the example data.
//
  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Get the count of the adjacencies.
//
  adj_num = triangulation_order6_adj_count ( node_num, triangle_num, 
    triangle_node, triangle_neighbor, adj_col );

  cout << "\n";
  cout << "  Number of adjacency entries is " << adj_num << "\n";

  cout << "\n";
  cout << "  Adjacency pointers:\n";
  cout << "\n";
  for ( node = 1; node <= node_num; node++ )
  {
    cout << "  " << setw(8) << node
         << "  " << setw(8) << adj_col[node-1]
         << "  " << setw(8) << adj_col[node]-1 << "\n";
  }
//
//  Get the adjacencies.
//
adj = triangulation_order6_adj_set ( node_num, triangle_num, triangle_node,
  triangle_neighbor, adj_num, adj_col );
//
//  Print the adjacencies.
//
  for ( node = 1; node <= node_num; node++ )
  {
    cout << "\n";
    cout << "  Nodes adjacent to node " << node << "\n";
    cout << "\n";

    for ( k = adj_col[node-1]; k <= adj_col[node]-1; k++ )
    {
      cout << "  " << setw(8) << adj[k-1] << "\n";
    }
  }

  delete [] adj;
  delete [] adj_col;
  delete [] node_xy;
  delete [] triangle_neighbor;
  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT.
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
{
  int boundary_edge_num;
  int dim_num = 2;
  int hole_num;
  int node_num;
  double *node_xy;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_num;
  int triangle_order = 6;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the\n";
  cout << "  boundary edges.\n";

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = new double[dim_num*node_num];
  triangle_node = new int[triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  boundary_edge_num = triangulation_order6_boundary_edge_count ( triangle_num, 
    triangle_node );

  cout << "\n";
  cout << "  Number of boundary edges = " << boundary_edge_num << "\n";
  cout << "  Correct number =           " << 16 << "\n";

  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER.
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
{
  int boundary_num;
  int hole_num;
  int node_num;
  int triangle_num;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER\n";
  cout << "  determines the number of edges that lie on the\n";
  cout << "  boundary of a region that has been triangulated.\n";

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  cout << "\n";
  cout << "  Number of nodes =          " << node_num << "\n";
  cout << "  Number of triangles =      " << triangle_num << "\n";
  cout << "  Number of holes =          " << hole_num << "\n";

  boundary_num = triangulation_order6_boundary_edge_count_euler ( node_num, 
    triangle_num, hole_num );

  cout << "  Number of boundary edges = " << boundary_num << "\n";
  cout << "  Correct number =           " << 16 << "\n";

  return;
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests TRIANGULATION_ORDER6_BOUNDARY_NODE.
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
{
  string file_name = "triangulation_order6_plot.eps";
  int i;
  int dim_num = 2;
  int hole_num;
  bool *node_boundary;
  int node_num;
  int node_show = 2;
  double *node_xy;
  int triangle_num;
  int *triangle_node;
  int *triangle_neighbor;
  int triangle_order = 6;
  int triangle_show = 2;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_BOUNDARY_COUNT counts the boundary\n";
  cout << "  edges.\n";
  cout << "  TRIANGULATION_ORDER6_PLOT plots the triangulation.\n";

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = new double[dim_num*node_num];
  triangle_node = new int[triangle_order*triangle_num];
  triangle_neighbor = new int[3*triangle_num];

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );
//
//  Make the plot.
//
  triangulation_order6_plot ( file_name, node_num, node_xy, triangle_num, 
    triangle_node, node_show, triangle_show );

  cout << "\n";
  cout << "  An Encapsulated PostScript image of this\n";
  cout << "  triangulation is in \"" << file_name << "\"\n";;

  node_boundary = triangulation_order6_boundary_node ( node_num, triangle_num, 
    triangle_node );

  cout << "\n";
  cout << "    Node  BN?\n";
  cout << "\n";

  for ( i = 1; i <= node_num; i++ )
  {
    cout << "  "
         << setw(6) << i << "  "
         << node_boundary[i-1] << "\n";
  }

  delete [] node_boundary;
  delete [] node_xy;
  delete [] triangle_node;
  delete [] triangle_neighbor;

  return;
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests TRIANGULATION_ORDER6_PRINT.
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
{
  int hole_num;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  int triangle_order = 6;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_PRINT prints the data.\n";

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = new double[2*node_num];
  triangle_neighbor = new int[3*triangle_num];
  triangle_node = new int[triangle_order*triangle_num];

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  triangulation_order6_print ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  delete [] node_xy;
  delete [] triangle_neighbor;
  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test265 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST265 tests TRIANGULATION_ORDER6_REFINE_COMPUTE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define NODE_NUM1 12
# define TRIANGLE_NUM1 3
# define TRIANGLE_ORDER 6

  int *edge_data;
  int node_num2;
  double node_xy1[DIM_NUM*NODE_NUM1] = {
       0.0, 0.0, 
       2.0, 0.0, 
       0.0, 2.0, 
       2.0, 2.0, 
       1.0, 3.0, 
       1.0, 0.0, 
       0.0, 1.0, 
       1.0, 1.0, 
       2.0, 1.0, 
       1.0, 2.0, 
       0.5, 2.5, 
       1.5, 2.5 };
  double *node_xy2;
  int triangle_node1[TRIANGLE_ORDER*TRIANGLE_NUM1] = {
       1,  2,  3,  6,  8,  7, 
       4,  3,  2,  9, 10,  8, 
       3,  4,  5, 10, 12, 11 };
  int *triangle_node2;
  int triangle_num2;

  cout << "\n";
  cout << "TEST265\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_REFINE_SIZE determines the\n";
  cout << "  size of a refined triangulation.\n";
  cout << "  TRIANGULATION_ORDER6_REFINE_COMPUTES computes the\n";
  cout << "  refined triangulation.\n";

  cout << "\n";
  cout << "  The number of nodes is " << NODE_NUM1 << "\n";
  cout << "  The number of triangles is " << TRIANGLE_NUM1 << "\n";

  r8mat_transpose_print ( DIM_NUM, NODE_NUM1, node_xy1, 
    "  The nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, TRIANGLE_NUM1, triangle_node1, 
    "  The triangles:" );

  edge_data = new int[5*(3*TRIANGLE_NUM1)];

  cout << "\n";
  cout << "  Sizing the refined mesh:\n";

  triangulation_order6_refine_size ( NODE_NUM1, TRIANGLE_NUM1, 
    triangle_node1, &node_num2, &triangle_num2, edge_data );

  cout << "\n";
  cout << "  Information about the refined mesh:\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num2 << "\n";
  cout << "  The number of triangles is " << triangle_num2 << "\n";

  cout << "\n";
  cout << "  Computing the refined mesh:\n";

  node_xy2 = new double[DIM_NUM*node_num2];
  triangle_node2 = new int[TRIANGLE_ORDER*triangle_num2];

  triangulation_order6_refine_compute ( NODE_NUM1, TRIANGLE_NUM1, 
    node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, node_xy2, 
    triangle_node2 );

  r8mat_transpose_print ( DIM_NUM, node_num2, node_xy2, 
    "  The refined nodes" );

  i4mat_transpose_print ( TRIANGLE_ORDER, triangle_num2, triangle_node2, 
    "  The refined triangles:" );

  delete [] edge_data;
  delete [] node_xy2;
  delete [] triangle_node2;

  return;
# undef DIM_NUM
# undef NODE_NUM1
# undef TRIANGLE_NUM1
# undef TRIANGLE_ORDER
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests TRIANGULATION_ORDER6_VERTEX_COUNT.
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
{
  int hole_num;
  int midside_num;
  int node_num;
  double *node_xy;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  int triangle_order = 6;
  int vertex_num;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  For an order6 triangulation:\n";
  cout << "  TRIANGULATION_ORDER6_VERTEX_COUNT counts the \n";
  cout << "  vertex nodes and midside nodes.\n";

  triangulation_order6_example1_size ( &node_num, &triangle_num, &hole_num );

  node_xy = new double[2*node_num];
  triangle_neighbor = new int[3*triangle_num];
  triangle_node = new int[triangle_order*triangle_num];

  triangulation_order6_example1 ( node_num, triangle_num, node_xy, 
    triangle_node, triangle_neighbor );

  vertex_num = triangulation_order6_vertex_count ( triangle_num, 
    triangle_node );

  midside_num = node_num - vertex_num;

  cout << "\n";
  cout << "  Number of nodes =         " << node_num << "\n";
  cout << "  Number of vertex nodes =  " << vertex_num << "\n";
  cout << "  Number of midside nodes = " << midside_num << "\n";

  delete [] node_xy;
  delete [] triangle_neighbor;
  delete [] triangle_node;

  return;
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests VORONOI_POLYGON_AREA.
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
{
# define DIM_NUM 2
# define NEIGHBOR_NUM 4
# define NODE_NUM 5

  double area;
  double area_correct = 0.5;
  int neighbor_index[NEIGHBOR_NUM] = { 0, 1, 2, 3 };
  int node = 4;
  double node_xy[DIM_NUM*NODE_NUM] = { 
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    0.5, 0.5 };

  cout << "\n";
  cout << "TEST31\n";
  cout << "  VORONOI_POLYGON_AREA computes the area of\n";
  cout << "  a finite Voronoi polygon.\n";

  area = voronoi_polygon_area ( node, NEIGHBOR_NUM, neighbor_index,
    NODE_NUM, node_xy );

  cout << "\n";
  cout << "  The computed area is " << area         << "\n";
  cout << "  The correct area is  " << area_correct << "\n";

  return;
# undef DIM_NUM
# undef NEIGHBOR_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests VORONOI_POLYGON_CENTROID.
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
{
# define DIM_NUM 2
# define NEIGHBOR_NUM 4
# define NODE_NUM 5

  double *centroid;
  double centroid_exact[2] = { 0.5, 0.5 };
  int neighbor_index[NEIGHBOR_NUM] = { 0, 1, 2, 3 };
  int node = 4;
  double node_xy[DIM_NUM*NODE_NUM] = { 
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    0.5, 0.5 };

  cout << "\n";
  cout << "TEST32\n";
  cout << "  VORONOI_POLYGON_CENTROID computes the centroid of\n";
  cout << "  a finite Voronoi polygon.\n";
  cout << flush;

  centroid = voronoi_polygon_centroid ( node, NEIGHBOR_NUM, 
    neighbor_index, NODE_NUM, node_xy );

  cout << "\n";
  cout << "  The computed centroid is " 
    << setw(10) << centroid[0] << "  "
    << setw(10) << centroid[1] << "\n";
  cout << "  The correct centroid is  " 
    << setw(10) << centroid_exact[0] << "  "
    << setw(10) << centroid_exact[1] << "\n";

  return;
# undef DIM_NUM
# undef NEIGHBOR_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests VORONOI_POLYGON_VERTICES.
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
{
# define DIM_NUM 2
# define NEIGHBOR_NUM 4
# define NODE_NUM 5

  int neighbor_index[NEIGHBOR_NUM] = { 0, 1, 2, 3 };
  int node = 4;
  double v[DIM_NUM*NEIGHBOR_NUM];
  double v_y[NEIGHBOR_NUM];
  double node_xy[DIM_NUM*NODE_NUM] = { 
    0.0, 0.0,
    1.0, 0.0,
    1.0, 1.0,
    0.0, 1.0,
    0.5, 0.5 };
 
  cout << "\n";
  cout << "TEST33\n";
  cout << "  VORONOI_POLYGON_VERTICES computes the vertices of\n";
  cout << "  a finite Voronoi polygon.\n";
  cout << flush;

  voronoi_polygon_vertices ( node, NEIGHBOR_NUM, neighbor_index, 
    NODE_NUM, node_xy, v );

  r8mat_transpose_print ( DIM_NUM, NEIGHBOR_NUM, v, "  Vertices:" );

  return;
# undef DIM_NUM
# undef NEIGHBOR_NUM
# undef NODE_NUM
}
