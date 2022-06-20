# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "gnuplot_i.h"

# define SLEEP_LGTH  1
# define NPOINTS     50

int main ( int argc, char *argv[] );
void slow_print ( char *s );

/******************************************************************************/

int main ( int argc, char *argv[] ) 

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for ANIM.

  Modified:

    24 June 2011
*/
{
  double d[NPOINTS];
  dpoint dp[NPOINTS];
  gnuplot_ctrl *h1;
  gnuplot_ctrl *h2;
  gnuplot_ctrl *h3;
  gnuplot_ctrl *h4;
  int i;
  int j;
  double x[NPOINTS];
  double y[NPOINTS];

  printf ( "\n" );
  printf ( "EXAMPLE:\n" );
  printf ( "  C++ version.\n" );
  printf ( "  Demonstrate how a running C++ program can produce plots\n" );
  printf ( "  through GNUPLOT, by invoking the GNUPLOT_I interface library.\n" );
/*
  Initialize the gnuplot handle
*/
  h1 = gnuplot_init();

  if ( h1 == NULL )
  {
    printf ( "\n" );
    printf ( "EXAMPLE - Fatal error!\n" );
    printf ( "  GNUPLOT is not available in your path.\n" );
    exit ( 1 );
  }
/* 
  Slopes 
*/    
  gnuplot_setstyle(h1, "lines");
    
  slow_print("*** plotting slopes\n");
  slow_print("y = x\n");
  gnuplot_plot_slope(h1, 1.0, 0.0, "unity slope");
  sleep(SLEEP_LGTH);

  slow_print("y = 2*x\n");
  gnuplot_plot_slope(h1, 2.0, 0.0, "y=2x");
  sleep(SLEEP_LGTH);

  slow_print("y = -x\n");
  gnuplot_plot_slope(h1, -1.0, 0.0, "y=-x");
  sleep(SLEEP_LGTH);
/* 
  Equations 
*/
  gnuplot_resetplot(h1);
  printf("\n\n");
  slow_print("*** various equations\n");
  slow_print("y = sin(x)\n");
  gnuplot_plot_equation(h1, "sin(x)", "sine");
  sleep(SLEEP_LGTH);

  slow_print("y = log(x)\n");
  gnuplot_plot_equation(h1, "log(x)", "logarithm");
  sleep(SLEEP_LGTH);

  slow_print("y = sin(x)*cos(2*x)\n");
  gnuplot_plot_equation(h1, "sin(x)*cos(2*x)", "sine product");
  sleep(SLEEP_LGTH);
/* 
  Styles 
*/
  gnuplot_resetplot(h1);
  printf("\n\n");
  slow_print("*** Showing plot style options:\n");

  slow_print("sine(x) in points\n");
  gnuplot_setstyle(h1, "points");
  gnuplot_plot_equation(h1, "sin(x)", "sine");
  sleep(SLEEP_LGTH);
    
  slow_print("sine(x) in impulses\n");
  gnuplot_setstyle(h1, "impulses");
  gnuplot_plot_equation(h1, "sin(x)", "sine");
  sleep(SLEEP_LGTH);
    
  slow_print("sine(x) in steps\n");
  gnuplot_setstyle(h1, "steps");
  gnuplot_plot_equation(h1, "sin(x)", "sine");
  sleep(SLEEP_LGTH);
 /*
   User defined 1d and 2d point sets
 */
  gnuplot_resetplot(h1);
  gnuplot_setstyle(h1, "impulses");
  printf("\n\n");
  slow_print("*** user defined lists of points\n");
  slow_print("random doubles\n");

  srand48 ( getpid ( ) );
  for ( i = 0; i < NPOINTS; i++ ) 
  {
    d[i] = drand48();
  }
  gnuplot_plot1d_var1(h1, d, NPOINTS, "random doubles");
  sleep(SLEEP_LGTH);

  gnuplot_resetplot(h1);
  gnuplot_setstyle(h1, "points");
  slow_print("random points\n");
  for ( i = 0; i < NPOINTS; i++ ) 
  {
    dp[i].x = drand48();
    dp[i].y = drand48();
  }
  gnuplot_plot1d_var2(h1, dp, NPOINTS, "random points");
  sleep(SLEEP_LGTH);

  gnuplot_resetplot(h1);
  gnuplot_setstyle(h1, "points");
  slow_print("cosine points with var2v\n");
  for ( j = 0; j < NPOINTS; j++ ) 
  {
    x[j] = ( double ) j / 10.0;
    y[j] = cos ( x[j] );
  }
  gnuplot_plot1d_var2v(h1, x, y, NPOINTS, "cosine points");
  sleep(SLEEP_LGTH);
/* 
  Multiple output screens 
*/
  printf("\n\n");
  slow_print("*** multiple output windows\n");
  gnuplot_resetplot(h1);
  gnuplot_setstyle(h1, "lines");
  h2 = gnuplot_init();
  gnuplot_setstyle(h2, "lines");
  h3 = gnuplot_init();
  gnuplot_setstyle(h3, "lines");
  h4 = gnuplot_init();
  gnuplot_setstyle(h4, "lines");

  slow_print("window 1: sin(x)\n");
  gnuplot_plot_equation(h1, "sin(x)", "sin(x)");
  sleep(SLEEP_LGTH);
  slow_print("window 2: cos(x)\n");
  gnuplot_plot_equation(h2, "cos(x)", "cos(x)");
  sleep(SLEEP_LGTH);
  slow_print("window 3: asin(x)\n");
  gnuplot_plot_equation(h3, "asin(x)", "arcsin(x)");
  sleep(SLEEP_LGTH);
  slow_print("window 4: acos(x)\n");
  gnuplot_plot_equation(h4, "acos(x)", "arccos(x)");
  sleep(SLEEP_LGTH);
/*  
  Close gnuplot handles. 
*/
  printf ( "\n\n" );
  slow_print ( "*** end of gnuplot example\n" );
  gnuplot_close ( h1 );
  gnuplot_close ( h2 );
  gnuplot_close ( h3 );
  gnuplot_close ( h4 );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "EXAMPLE:\n" );
  printf ( "  Normal end of execution.\n" );

  return 0;
}
/******************************************************************************/

void slow_print ( char *s )

/******************************************************************************/
{
  int i;

  for ( i = 0; i < ( int ) strlen ( s ); i++ )
  {
    printf ( "%c", s[i] );
    fflush ( stdout );
    usleep ( 100000 );
  }
  return;
}

