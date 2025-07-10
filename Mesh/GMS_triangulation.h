void alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *alpha_min, double *alpha_ave,
  double *alpha_area );
double angle_rad_2d ( double p1[2], double p2[2], double p3[2] );
double arc_cosine ( double c );
void area_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *area_min, double *area_max, double *area_ratio,
  double *area_ave, double *area_std );
void bandwidth ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
bool delaunay_swap_test ( double xy[] );
int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 );
unsigned long get_seed ( );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_power ( int i, int j );
int i4_sign ( int i );
void i4_swap ( int *i, int *j );
int i4_uniform ( int a, int b, int *seed );
int i4_wrap ( int ival, int ilo, int ihi );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_sort_a ( int m, int n, int a[] );
int i4col_sorted_unique_count ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void i4vec_heap_d ( int n, int a[] );
int *i4vec_indicator_new ( int n );
int i4vec_min ( int n, int a[] );
void i4vec_print ( int n, int a[], string title );
void i4vec_reverse ( int n, int a[] );
void i4vec_sort_heap_a ( int n, int a[] );
int i4vec_sorted_unique ( int n, int a[] );
int i4vec2_compare ( int n, int a1[], int a2[], int i, int j );
void i4vec2_sort_a ( int n, int a1[], int a2[] );
int i4vec2_sorted_unique ( int n, int a1[], int a2[] );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv );
void lvec_print ( int n, bool a[], string title );
void mesh_base_one ( int node_num, int element_order,
  int element_num, int element_node[] );
void mesh_base_zero ( int node_num, int element_order,
  int element_num, int element_node[] );
void node_merge ( int dim_num, int node_num, double node_xy[],
  double tolerance, int node_rep[] );
int ns_adj_col_set ( int node_num, int triangle_num, int variable_num,
  int triangle_node[], int triangle_neighbor[], int node_u_variable[],
  int node_v_variable[], int node_p_variable[], int adj_col[] );
int ns_adj_count ( int node_num, int triangle_num, int variable_num,
  int triangle_node[], int triangle_neighbor[], int node_u_variable[],
  int node_v_variable[], int node_p_variable[] );
void ns_adj_insert ( int v1, int v2, int variable_num, int adj_num,
  int adj_col_free[], int adj_row[] );
void ns_adj_row_set ( int node_num, int triangle_num, int variable_num,
  int triangle_node[], int triangle_neighbor[], int node_u_variable[],
  int node_v_variable[], int node_p_variable[], int adj_num, int adj_col[],
  int adj_row[] );
bool perm_check2 ( int n, int p[], int base );
void perm_inverse ( int n, int p[] );
int *points_delaunay_naive_2d ( int node_num, double node_xy[],
  int *triangle_num );
void points_hull_2d ( int node_num, double node_xy[], int *nval, int ival[] );
int points_point_near_naive_nd ( int dim_num, int nset, double pset[],
  double ptest[], double *d_min );
void q_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *q_min, double *q_max, double *q_ave,
  double *q_area );
void quad_convex_random ( int *seed, double xy[] );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_abs ( double x );
double r8_epsilon ( );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
double r8_uniform_01 ( int *seed );
void r82vec_permute ( int n, int p[], int base, double a[] );
int *r82vec_sort_heap_index_a ( int n, int base, double a[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, string title );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
int r8tris2 ( int node_num, double node_xy[], int *triangle_num,
  int triangle_node[], int triangle_neighbor[] );
void r8vec_bracket ( int n, double x[], double xval, int *left, int *right );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
int s_len_trim ( string s );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
int swapec ( int i, int *top, int *btri, int *bedg, int node_num,
  double node_xy[], int triangle_num, int triangle_node[],
  int triangle_neighbor[], int stack[] );
void timestamp ( );
//
//  Triangle utilities.
//
double *triangle_angles_2d_new ( double t[2*3] );
double triangle_area_2d ( double t[2*3] );
double *triangle_circumcenter_2d ( double t[] );
void triangle_order3_physical_to_reference ( double t[], int n,
  double phy[], double ref[] );
void triangle_order3_reference_to_physical ( double t[], int n,
  double ref[], double phy[] );
void triangle_order6_physical_to_reference ( double t[2*6], int n,
  double phy[], double ref[] );
void triangle_order6_reference_to_physical ( double t[], int n,
  double ref[], double phy[] );
void triangle_reference_sample ( int n, int *seed, double p[] );
void triangle_sample ( double t[2*3], int n, int *seed, double p[] );
//
//  Triangulation routines that don't depend on order.
//
double triangulation_area ( int node_num, double node_xy[], int element_order,
  int element_num, int element_node[] );
double triangulation_areas ( int node_num, double node_xy[], int triangle_order,
  int triangle_num, int triangle_node[], double triangle_area[] );
double triangulation_delaunay_discrepancy_compute ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double *angle_min, int *angle_min_triangle,
  double *angle_max, int *angle_max_triangle );
int *triangulation_neighbor_elements ( int triangle_order, int triangle_num,
  int triangle_node[] );
int *triangulation_node_order ( int triangle_order, int triangle_num,
  int triangle_node[], int node_num );
//
//  Order3 Triangulations.
//
int triangulation_order3_adj_count ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_col[] );
int *triangulation_order3_adj_set ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_num, int adj_col[] );
void triangulation_order3_adj_set2 ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_num, int adj_col[],
  int ia[], int ja[] );
int *triangulation_order3_adjacency ( int node_num, int element_num, 
  int element_node[] );
int triangulation_order3_boundary_edge_count ( int triangle_num,
  int triangle_node[] );
int triangulation_order3_boundary_edge_count_euler ( int node_num,
  int triangle_num, int hole_num );
bool *triangulation_order3_boundary_node ( int node_num, int triangle_num,
  int triangle_node[] );
int triangulation_order3_check ( int node_num, int triangle_num,
  int triangle_node[] );
int triangulation_order3_edge_check ( int triangle_num, int triangle_node[] );
void triangulation_order3_example1 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_order3_example1_size ( int *node_num, int *triangle_num,
  int *hole_num );
void triangulation_order3_example2 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_order3_example2_size ( int *node_num, int *triangle_num,
  int *hole_num );
void triangulation_order3_neighbor ( int triangle_num, int triangle_node[],
  int t1, int s1, int *t2, int *s2 );
void triangulation_order3_neighbor_nodes ( int node_num, int triangle_num,
  int triangle_node[], int nabes_first[], int nabes_num[], int nabes_max,
  int *nabes_dim, int nabes[] );
void triangulation_order3_neighbor_nodes_print ( int node_num,
  int nabes_first[], int nabes_num[], int nabes_dim, int nabes[] );
void triangulation_order3_plot ( string file_name, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show );
void triangulation_order3_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_order3_quad ( int node_num, double node_xy[], int triangle_order,
  int triangle_num, int triangle_node[],
  void f ( int n, double xy_vec[], double f_vec[] ), int quad_num,
  double quad_xy[], double quad_w[], double *quad_value, double *region_area );
void triangulation_order3_refine_compute ( int node_num1, int triangle_num1,
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2,
  int edge_data[], double node_xy2[], int triangle_node2[] );
void triangulation_order3_refine_size ( int node_num1, int triangle_num1,
  int triangle_node1[], int *node_num2, int *triangle_num2, int edge_data[] );
void triangulation_order3_sample ( int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int num_ran, int *seed, double xd[],
  int td[] );
//
//  Order4 Triangulations.
//
void triangulation_order4_plot ( string plot_filename, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show );
//
//  Order6 Triangulations.
//
int triangulation_order6_adj_count ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_col[] );
int *triangulation_order6_adj_set ( int node_num, int triangle_num,
  int triangle_node[], int triangle_neighbor[], int adj_num, int adj_col[] );
int triangulation_order6_boundary_edge_count ( int triangle_num,
  int triangle_node[] );
int triangulation_order6_boundary_edge_count_euler ( int node_num,
  int triangle_num, int hole_num );
bool *triangulation_order6_boundary_node ( int node_num, int triangle_num,
  int triangle_node[] );
void triangulation_order6_example1 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_order6_example1_size ( int *node_num, int *triangle_num,
  int *hole_num );
void triangulation_order6_example2 ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_order6_example2_size ( int *node_num, int *triangle_num,
  int *hole_num );
void triangulation_order6_neighbor ( int triangle_num, int triangle_node[],
  int t1, int s1, int  *t2, int *s2 );
void triangulation_order6_plot ( string file_name, int node_num,
  double node_xy[], int triangle_num, int triangle_node[], int node_show,
  int triangle_show );
void triangulation_order6_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] );
void triangulation_order6_refine_compute ( int node_num1, int triangle_num1,
  double node_xy1[], int triangle_node1[], int node_num2, int triangle_num2,
  int edge_data[], double node_xy2[], int triangle_node2[] );
void triangulation_order6_refine_size ( int node_num1, int triangle_num1,
  int triangle_node1[], int *node_num2, int *triangle_num2,
  int edge_data[] );
int *triangulation_order6_to_order3 ( int triangle_num1, int triangle_node1[] );
int triangulation_order6_vertex_count ( int tri_num, int triangle_node[] );
//
//  Triangulation routines that don't depend on order.
//
void triangulation_search_delaunay ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[],
  int triangle_neighbor[], double p[2], int *triangle_index, 
  double *alpha, double *beta, double *gamma, int *edge,
  int *step_num );
int triangulation_search_naive ( int node_num, double node_xy[],
  int triangle_order, int triangle_num, int triangle_node[], double p[2] );
//
//  More utilities.
//
void vbedg ( double x, double y, int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int triangle_neighbor[],
  int *ltri, int *ledg, int *rtri, int *redg );
//
//  Voronoi polygons..
//
double voronoi_polygon_area ( int node, int neighbor_num,
  int neighbor_index[], int node_num, double node_xy[] );
double *voronoi_polygon_centroid ( int node, int neighbor_num,
  int neighbor_index[], int node_num, double node_xy[] );
void voronoi_polygon_vertices ( int node, int neighbor_num,
  int neighbor_index[], int node_num, double node_xy[], double v[] );
