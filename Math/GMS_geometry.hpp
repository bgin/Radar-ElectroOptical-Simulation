
#ifndef __GMS_GEOMETRY_HPP__
#define __GMS_GEOMETRY_HPP__

void angle_box_2d ( double dist, double p1[2], double p2[2], double p3[2], 
  double p4[2], double p5[2] );
bool angle_contains_ray_2d ( double p1[2], double p2[2], double p3[2], 
  double p[2] );
double angle_deg_2d ( double p1[2], double p2[2], double p3[2] );
double *angle_half_2d ( double p1[2], double p2[2], double p3[2] );
double angle_rad_2d ( double p1[2], double p2[2], double p3[2] );
double angle_rad_3d ( double p1[3], double p2[3], double p3[3] );
double angle_rad_nd ( int n, double vec1[], double vec2[] );
double angle_turn_2d ( double p1[2], double p2[2], double p3[2] );
double anglei_deg_2d ( double p1[2], double p2[2], double p3[2] );
double anglei_rad_2d ( double p1[2], double p2[2], double p3[2] );
double annulus_area_2d ( double r1, double r2 );
double annulus_sector_area_2d ( double r1, double r2, double theta1, 
  double theta2 );
double *annulus_sector_centroid_2d ( double pc[2], double r1, double r2, 
  double theta1, double theta2 );
double *ball_unit_sample_2d ( int &seed );
double *ball_unit_sample_3d ( int &seed );
double *ball_unit_sample_nd ( int n, int &seed );
double *basis_map_3d ( double u[3*3], double v[3*3] );
bool box_01_contains_point_2d ( double p[2] );
bool box_01_contains_point_nd ( int ndim, double p[] );
bool box_contains_point_2d ( double p1[2], double p2[2], double p[2] );
bool box_contains_point_nd ( int ndim, double p1[], double p2[], double p[] );
void box_ray_int_2d ( double p1[2], double p2[2], double pa[2], 
  double pb[2], double pint[2] );
int box_segment_clip_2d ( double p1[2], double p2[2], double pa[2], 
  double pb[2] );
void circle_arc_point_near_2d ( double r, double center[2], double theta1, 
  double theta2, double p[2], double pn[2], double *dist );
double circle_area_2d ( double r );
void circle_dia2imp_2d ( double p1[2], double p2[2], double *r, 
  double center[2] );
int circle_exp_contains_point_2d ( double p1[2], double p2[2], double p3[2], 
  double p[2] );
void circle_exp2imp_2d ( double p1[2], double p2[2], double p3[2], double *r, 
  double pc[2] );
bool circle_imp_contains_point_2d ( double r, double pc[2], double p[2] );
void circle_imp_line_par_int_2d ( double r, double pc[2], double x0, double y0, 
  double f, double g, int *int_num, double p[] );
double circle_imp_point_dist_2d ( double r, double pc[2], double p[2]  );
double circle_imp_point_dist_signed_2d ( double r, double pc[2], double p[2] );
double circle_imp_point_near_2d ( double r, double pc[2], double p[2], 
  double pn[2] );
double *circle_imp_points_2d ( double r, double pc[2], int n );
double *circle_imp_points_3d ( double r, double pc[3], double nc[3], int n );
void circle_imp_points_arc_2d ( double r, double pc[2], double theta1, 
  double theta2, int n, double p[] );
void circle_imp_print_2d ( double r, double center[2], string title );
void circle_imp_print_3d ( double r, double pc[3], double nc[3], string title );
void circle_imp2exp_2d ( double r, double pc[2], double p1[2], double p2[2], 
  double p3[2] );
double *circle_llr2imp_2d ( double p1[], double p2[], double q1[], double q2[], 
  double r );
double circle_lune_area_2d ( double r, double pc[2], double theta1, 
  double theta2 );
double *circle_lune_centroid_2d ( double r, double pc[2], double theta1,
  double theta2 );
void circle_pppr2imp_3d ( double p1[], double p2[], double p3[], double r, 
  double pc[], double normal[] );
double *circle_ppr2imp_2d ( double p1[], double p2[], double r );
double circle_sector_area_2d ( double r, double pc[2], double theta1, 
  double theta2 );
double *circle_sector_centroid_2d ( double r, double pc[2], double theta1, 
  double theta2 );
bool circle_sector_contains_point_2d ( double r, double pc[2], double theta1, 
  double theta2, double p[2] );
void circle_sector_print_2d ( double r, double pc[2], double theta1, 
  double theta2 );
double circle_triangle_area_2d ( double r, double pc[2], double theta1, 
  double theta2 );
void circle_triple_angles_2d ( double r1, double r2, double r3, double *angle1, 
  double *angle2, double *angle3 );
void circles_imp_int_2d ( double r1, double pc1[2], double r2, double pc2[2],
  int *num_int, double p[] );
double cone_area_3d ( double h, double r );
double *cone_centroid_3d ( double r, double pc[3], double pt[3] );
double cone_volume_3d ( double h, double r );
void conv3d ( char axis, double theta, int n, double cor3[], double cor2[] );
double cot_rad ( double angle );
void cube_shape_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[] );
void cube_size_3d ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
double cylinder_point_dist_3d ( double p1[3], double p2[3], double r, 
  double p[3] );
double cylinder_point_dist_signed_3d ( double p1[3], double p2[3], double r, 
  double p[3] );
bool cylinder_point_inside_3d ( double p1[3], double p2[3], double r, 
  double p[3] );
double *cylinder_point_near_3d ( double p1[3], double p2[3], double r, 
  double p[3] );
double *cylinder_sample_3d ( double p1[3], double p2[3], double r, int n, 
  int &seed );
double cylinder_volume_3d ( double p1[3], double p2[3], double r );
double degrees_to_radians ( double angle );
double dge_det ( int n, double a[], int pivot[] );
int dge_fa ( int n, double a[], int pivot[] );
void dge_sl ( int n, double a[], int pivot[], double b[], int job );
double *direction_pert_3d ( double sigma, double vbase[3], int &seed );
double *direction_uniform_2d ( int &seed );
double *direction_uniform_3d ( int &seed );
double *direction_uniform_nd ( int n, int &seed );
double disk_point_dist_3d ( double pc[3], double r, double axis[3], 
  double p[3] );
double dms_to_radians ( int degrees, int minutes, int seconds );
void dodec_shape_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[] );
void dodec_size_3d ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
void dual_shape_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[], int point_num2, 
  int face_num2, int face_order_max2, double point_coord2[], int face_order2[], 
  int face_point2[] );
void dual_size_3d ( int point_num, int edge_num, int face_num, 
  int face_order_max, double point_coord[], int face_order[], int face_point[], 
  int *point_num2, int *edge_num2, int *face_num2, int *face_order_max2 );
double ellipse_area_2d ( double r1, double r2 );
double ellipse_point_dist_2d ( double r1, double r2, double p[2] );
double *ellipse_point_near_2d ( double r1, double r2, double p[2] );
void ellipse_points_2d ( double pc[2], double r1, double r2, double psi, 
  int n, double p[] );
void ellipse_points_arc_2d ( double pc[2], double r1, double r2, double psi,
  double theta1, double theta2, int n, double p[] );
double enorm0_nd ( int n, double x[], double y[] );
int get_seed ( );
void glob2loc_3d ( double cospitch, double cosroll, double cosyaw, 
  double sinpitch, double sinroll, double sinyaw, double globas[3], 
  double glopts[3], double locpts[3] );
bool halfplane_contains_point_2d ( double pa[2], double pb[2], double p[2] );
int halfspace_imp_triangle_int_3d ( double a, double b, double c, double d, 
  double t[3*3], double p[3*4] );
int halfspace_norm_triangle_int_3d ( double pp[3], double pn[3], double t[3*3], 
  double p[3*4] );
int halfspace_triangle_int_3d ( double dist1, double dist2, double dist3, 
  double t[3*3], double p[3*4] );
double haversine ( double a );
void helix_shape_3d ( double a, int n, double r, double theta1, double theta2, 
  double p[] );
double hexagon_area_2d ( double r );
bool hexagon_contains_point_2d ( double v[2*6], double p[2] );
void hexagon_shape_2d ( double angle, double p[2] );
double hexagon_unit_area_2d ( );
void hexagon_vertices_2d ( double h[2*6] );
double i4_dedekind_factor ( int p, int q );
double i4_dedekind_sum ( int p, int q );
int i4_factorial2 ( int n );
int i4_gcd ( int i, int j );
int i4_lcm ( int i, int j );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_sign ( int i );
void i4_swap ( int *i, int *j );
int i4_uniform ( int a, int b, int &seed );
int i4_wrap ( int ival, int ilo, int ihi );
int i4col_compare ( int m, int n, int a[], int i, int j );
void i4col_find_item ( int m, int n, int a[], int item, int *row, int *col );
void i4col_find_pair_wrap ( int m, int n, int a[], int item1, int item2,
  int *row, int *col );
void i4col_sort_a ( int m, int n, int a[] );
int i4col_sorted_unique_count ( int m, int n, int a[] );
void i4col_swap ( int m, int n, int a[], int icol1, int icol2 );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int i4row_compare ( int m, int n, int a[], int i, int j );
void i4row_sort_a ( int m, int n, int a[] );
void i4row_swap ( int m, int n, int a[], int irow1, int irow2 );
void i4vec_copy ( int n, int a1[], int a2[] );
void i4vec_heap_d ( int n, int a[] );
int *i4vec_indicator_new ( int n );
int i4vec_lcm ( int n, int v[] );
void i4vec_print ( int n, int a[], string title );
int i4vec_product ( int n, int a[] );
void i4vec_reverse ( int n, int a[] );
void i4vec_sort_heap_a ( int n, int a[] );
void i4vec_sorted_unique ( int n, int a[], int *nuniq );
int *i4vec_uniform_new ( int n, int a, int b, int &seed );
void i4vec_zero ( int n, int a[] );
int i4vec2_compare ( int n, int a1[], int a2[], int i, int j );
void i4vec2_sort_a ( int n, int a1[], int a2[] );
void i4vec2_sorted_unique ( int n, int a1[], int a2[], int *nuniq );
void icos_shape ( int point_num, int edge_num, int face_num, 
  int face_order_max, double point_coord[], int edge_point[], int face_order[],
  int face_point[] );
void icos_size ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
bool line_exp_is_degenerate_nd ( int dim_num, double p1[], double p2[] );
double *line_exp_normal_2d ( double p1[2], double p2[2] );
double *line_exp_perp_2d ( double p1[2], double p2[2], double p3[2], bool *flag );
double line_exp_point_dist_2d ( double p1[2], double p2[2], double p[2] );
double line_exp_point_dist_3d ( double p1[3], double p2[3], double p[3] );
double line_exp_point_dist_signed_2d ( double p1[2], double p2[2], double p[2] );
void line_exp_point_near_2d ( double p1[2], double p2[2], double p[2], 
  double pn[2], double *dist, double *t );
void line_exp_point_near_3d ( double p1[3], double p2[3], double p[3], 
  double pn[3], double *dist, double *t );
void line_exp2imp_2d ( double p1[2], double p2[2], double *a, double *b, 
  double *c );
void line_exp2par_2d ( double p1[2], double p2[2], double *f, double *g, 
  double *x0, double *y0 );
void line_exp2par_3d ( double p1[3], double p2[3], double *f, double *g, 
  double *h, double *x0, double *y0, double *z0 );
bool line_imp_is_degenerate_2d ( double a, double b, double c );
double line_imp_point_dist_2d ( double a, double b, double c, double p[2] );
double line_imp_point_dist_signed_2d ( double a, double b, double c, 
  double p[2] );
void line_imp2exp_2d ( double a, double b, double c, double p1[2], 
  double p2[2] );
void line_imp2par_2d ( double a, double b, double c, double *f, double *g, 
  double *x0, double *y0 );
double line_par_point_dist_2d ( double f, double g, double x0, double y0, 
  double p[2] );
double line_par_point_dist_3d ( double f, double g, double h, double x0, 
  double y0, double z0, double p[3] );
double *line_par_point_near_2d ( double f, double g, double x0, double y0, 
  double p[2] );
double *line_par_point_near_3d ( double f, double g, double h, double x0, 
  double y0, double z0, double p[3] );
void line_par2exp_2d ( double f, double g, double x0, double y0, 
  double p1[2], double p2[2] );
void line_par2exp_3d ( double f, double g, double h, double x0, double y0, 
  double z0, double p1[3], double p2[3] );
void line_par2imp_2d ( double f, double g, double x0, double y0, double *a, 
  double *b, double *c );
double lines_exp_angle_3d ( double p1[3], double p2[3], double p3[3], 
  double p4[3] );
double lines_exp_angle_nd ( double p1[], double p2[], double q1[], double q2[], 
  int n );
double lines_exp_dist_3d ( double p1[3], double p2[3], double q1[3], 
  double q2[3] );
double lines_exp_dist_3d_2 ( double p1[3], double p2[3], double q1[3], 
  double q2[3] );
bool lines_exp_equal_2d ( double p1[2], double p2[2], double q1[2], 
  double q2[2] );
void lines_exp_int_2d ( double p1[2], double p2[2], double p3[2], 
  double p4[2], int *ival, double p[2] );
void lines_exp_near_3d ( double p1[3], double p2[3], double q1[3], 
  double q2[3], double pn[3], double qn[3] );
bool lines_exp_parallel_2d ( double p1[2], double p2[2], double q1[2], 
  double q2[2] );
bool lines_exp_parallel_3d ( double p1[3], double p2[3], double q1[3], 
  double q2[3] );
double lines_imp_angle_2d ( double a1, double b1, double c1, 
  double a2, double b2, double c2 );
double lines_imp_dist_2d ( double a1, double b1, double c1, double a2, 
  double b2, double c2 );
void lines_imp_int_2d ( double a1, double b1, double c1, double a2, double b2, 
  double c2, int *ival, double p[2] );
double lines_par_angle_2d ( double f1, double g1, double x01, double y01, 
  double f2, double g2, double x02, double y02 );
double lines_par_angle_3d ( double f1, double g1, double h1, double x01, 
  double y01, double z01, double f2, double g2, double h2, double x02, 
  double y02, double z02 );
double lines_par_dist_3d ( double f1, double g1, double h1, double x01, 
  double y01, double z01, double f2, double g2, double h2, double x02, 
  double y02, double z02 );
void lines_par_int_2d ( double f1, double g1, double x1, double y1, double f2, 
  double g2, double x2, double y2, double *t1, double *t2, double pint[2] );
void loc2glob_3d ( double cospitch, double cosroll, double cosyaw, 
  double sinpitch, double sinroll, double sinyaw, double locpts[3],
  double globas[3], double glopts[3] );
int lrline ( double xu, double yu, double xv1, double yv1, double xv2, 
  double yv2, double dv );
void lvec_print ( int n, bool a[], string title );
void minabs ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double *xmin, double *ymin );
bool minquad ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double *xmin, double *ymin );
void octahedron_shape_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[] );
void octahedron_size_3d ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
int parabola_ex ( double x1, double y1, double x2, double y2, double x3, 
  double y3, double *x, double *y );
int parabola_ex2 ( double x1, double y1, double x2, double y2, double x3, 
  double y3, double *x, double *y, double *a, double *b, double *c );
double parallelogram_area_2d ( double p[] );
double parallelogram_area_3d ( double p[] );
bool parallelogram_contains_point_2d ( double p1[2], double p2[2], double p3[2], 
  double p[2] );
bool parallelogram_contains_point_3d ( double p1[3], double p2[3], double p3[3], 
  double p[3] );
double parallelogram_point_dist_3d ( double p1[3], double p2[3], double p3[3], 
  double p[3] );
bool parallelepiped_contains_point_3d ( double p1[3], double p2[3], 
  double p3[3], double p4[3], double p[3] );
double parallelepiped_point_dist_3d ( double p1[3], double p2[3], double p3[3], 
  double p4[3], double p[3] );
bool perm_check ( int n, int p[] );
void perm_inv ( int n, int p[] );
void plane_exp_grid_3d ( double p1[3], double p2[3], double p3[3], int *ncor3, 
  int *line_num, double cor3[], int lines[], int maxcor3, int line_max, 
  int *ierror );
double plane_exp_point_dist_3d ( double p1[3], double p2[3], double p3[3], 
  double p[3] );
void plane_exp_normal_3d ( double p1[3], double p2[3], double p3[3], 
  double pn[3] );
void plane_exp_pro2 ( double p1[3], double p2[3], double p3[3], int npnt, 
  double pp[], double alpha[], double beta[] );
void plane_exp_pro3 ( double p1[3], double p2[3], double p3[3], int n, 
  double po[], double pp[] );
void plane_exp_project_3d ( double p1[3], double p2[3], double p3[3], 
  double pf[3], int n, double po[], double pp[], int ivis[] );
void plane_exp2imp_3d ( double p1[3], double p2[3], double p3[3], double *a, 
  double *b, double *c, double *d );
void plane_exp2normal_3d ( double p1[3], double p2[3], double p3[3],
  double pp[3], double pn[3] );
bool plane_imp_is_degenerate_3d ( double a, double b, double c );
bool plane_imp_line_par_int_3d ( double a, double b, double c, double d, 
  double x0, double y0, double z0, double f, double g, double h, double p[3] );
double plane_imp_point_dist_3d ( double a, double b, double c, double d, 
  double p[3] );
double plane_imp_point_dist_signed_3d ( double a, double b, double c, double d, 
  double p[3] );
void plane_imp_point_near_3d ( double a, double b, double c, double d, 
  double p[3], double pn[3] );
void plane_imp_segment_near_3d ( double p1[3], double p2[3], double a, double b, 
  double c, double d, double *dist, double pnp[3], double pnl[3] );
void plane_imp_triangle_int_3d ( double a, double b, double c, double d, 
  double t[3*3], int *int_num, double p[] );
void plane_imp_triangle_int_add_3d ( double p1[3], double p2[3], double dist1, 
  double dist2, int *int_num, double p[] );
int plane_imp_triangle_near_3d ( double t[3*3], double a, double b, double c, 
  double d, double *dist, double pn[] );
void plane_imp2exp_3d ( double a, double b, double c, double d, double p1[3], 
  double p2[3], double p3[3] );
void plane_imp2normal_3d ( double a, double b, double c, double d, 
  double pp[3], double pn[3] );
void plane_normal_basis_3d ( double pp[3], double pn[3], double pq[3], 
  double pr[3] );
int plane_normal_line_exp_int_3d ( double pp[3], double normal[3], 
  double p1[3], double p2[3], double pint[3] );
double *plane_normal_qr_to_xyz ( double pp[], double normal[], double pq[], 
  double pr[], int n, double qr[] );
void plane_normal_tetrahedron_intersect ( double pp[3], 
  double normal[3], double t[3*4], int *int_num, double pint[3*4] );
int plane_normal_triangle_int_3d ( double pp[3], double pn[3], double t[3*3], 
  double p[3*3] );
void plane_normal_uniform_3d ( int &seed, double pp[3], double normal[3] );
void plane_normal_uniform_nd ( int dim_num, int &seed, double pp[], 
  double normal[] );
double *plane_normal_xyz_to_qr ( double pp[], double normal[], double pq[], 
  double pr[], int n, double xyz[] );
void plane_normal2exp_3d ( double pp[3], double pn[3], double p1[3], 
  double p2[3], double p3[3] );
void plane_normal2imp_3d ( double pp[3], double pn[3], double *a, double *b, 
  double *c, double *d );
double planes_imp_angle_3d ( double a1, double b1, double c1, double d1, 
  double a2, double b2, double c2, double d2 );
bool points_avoid_point_naive_2d ( int n, double pset[], double p[2] );
void points_bisect_line_imp_2d ( double p1[2], double p2[2], double *a, 
  double *b, double *c );
void points_bisect_line_par_2d ( double p1[2], double p2[2], double *f, 
  double *g, double *x, double *y );
int points_centroid_2d ( int n, double p[] );
double points_colin_2d ( double p1[2], double p2[2], double p3[2] );
double points_colin_3d ( double p1[3], double p2[3], double p3[3] );
int *points_delaunay_naive_2d ( int n, double p[], int *ntri );
double points_dist_2d ( double p1[2], double p2[2] );
double points_dist_3d ( double p1[3], double p2[3] );
double points_dist_nd ( int n, double p1[], double p2[] );
void points_hull_2d ( int node_num, double node_xy[], int *hull_num, 
  int hull[] );
void points_plot ( string file_name, int node_num, double node_xy[], 
  bool node_label );
int points_point_near_naive_2d ( int nset, double pset[], double ptest[], 
  double *d_min );
int points_point_near_naive_3d ( int nset, double pset[], double ptest[], 
  double *d_min );
int points_point_near_naive_nd ( int ndim, int nset, double pset[], 
  double ptest[], double *d_min );
int *points_points_near_naive_2d ( int nset, double pset[], int ntest, 
  double ptest[] );
int *points_points_near_naive_3d ( int nset, double pset[], int ntest, 
  double ptest[] );
void polar_to_xy ( double r, double t, double xy[2] );
double polygon_1_2d ( int n, double v[] );
double *polygon_angles_2d ( int n, double v[] );
double polygon_area_2d ( int n, double v[] );
double polygon_area_2d_2 ( int n, double v[] );
double polygon_area_3d ( int n, double v[], double normal[] );
double polygon_area_3d_2 ( int n, double v[] );
double *polygon_centroid_2d ( int n, double v[] );
double *polygon_centroid_2d_2 ( int n, double v[] );
double *polygon_centroid_3d ( int n, double v[] );
bool polygon_contains_point_2d ( int n, double v[], double p[2] );
bool polygon_contains_point_2d_2 ( int n, double v[], double p[2] );
double polygon_diameter_2d ( int n, double v[] );
double *polygon_expand_2d ( int n, double v[], double h );
void polygon_inrad_data_2d ( int n, double radin, double *area, double *radout, 
  double *side );
int polygon_is_convex ( int n, double v[] );
double polygon_lattice_area_2d ( int i, int b );
double *polygon_normal_3d ( int n, double v[] );
void polygon_outrad_data_2d ( int n, double radout, double *area, double *radin, 
  double *side );
void polygon_side_data_2d ( int n, double side, double *area, double *radin, 
  double *radout );
double polygon_solid_angle_3d ( int n, double v[], double p[3] );
double polygon_x_2d ( int n, double v[] );
double polygon_y_2d ( int n, double v[] );
double polygon_xx_2d ( int n, double v[] );
double polygon_xy_2d ( int n, double v[] );
double polygon_yy_2d ( int n, double v[] );
double polyhedron_area_3d ( double coord[], int maxorder, int face_num, 
  int node[], int node_num, int order[] );
double *polyhedron_centroid_3d ( double coord[], int maxorder, int face_num, 
  int node[], int node_num, int order[] );
bool polyhedron_contains_point_3d ( int node_num, int face_num, 
  int face_order_max, double v[], int face_order[], int face_point[],
  double p[3] );
double polyhedron_volume_3d ( double coord[], int maxorder, int face_num, 
  int node[], int node_num, int order[] );
double polyhedron_volume_3d_2 ( double coord[], int maxorder, int face_num, 
  int node[], int node_num, int order[] );
double *polyline_arclength_nd ( int ndim, int n, double p[] );
double *polyline_index_point_nd ( int ndim, int n, double p[], double t );
double polyline_length_nd ( int ndim, int n, double p[] );
double *polyline_points_nd ( int ndim, int n, double p[], int nt );
double *polyloop_arclength_nd ( int dim_num, int nk, double pk[] );
double *polyloop_points_nd ( int dim_num, int nk, double pk[], int nt );
void provec ( int m, int n, double base[], double vecm[], double vecn[], 
  double vecnm[] );
double pyramid_volume_3d ( double h, double s );
double quad_area_2d ( double q[2*4] );
double quad_area2_2d ( double q[] );
double quad_area_3d ( double q[] );
bool quad_contains_point_2d ( double q[2*4], double p[2] );
void quad_convex_random ( int &seed, double xy[] );
double quad_point_dist_2d ( double q[2*4], double p[2] );
double quad_point_dist_signed_2d ( double q[2*4], double p[2] );
void quad_point_near_2d ( double q[2*4], double p[2], double pn[2], 
  double *dist );
double *quat_conj ( double q[] );
double *quat_inv ( double q[] );
double *quat_mul ( double q1[], double q2[] );
double quat_norm ( double q[] );
float r4_abs ( float x );
int r4_nint ( float x );
double r8_abs ( double x );
double r8_acos ( double c );
double r8_asin ( double s );
double r8_atan ( double y, double x );
double r8_cosd ( double angle );
double r8_cotd ( double angle );
double r8_cscd ( double angle );
double r8_epsilon ( );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_modp ( double x, double y );
int r8_nint ( double x );
double r8_normal_01 ( int &seed );
double r8_pi ( );
double r8_secd ( double angle );
double r8_sign ( double x );
bool r8_sign_opposite_strict ( double r1, double r2 );
double r8_sind ( double angle );
void r8_swap ( double *x, double *y );
double r8_tand ( double angle );
double r8_uniform ( double a, double b, int &seed );
double r8_uniform_01 ( int &seed );
void r82vec_part_quick_a ( int n, double a[], int *l, int *r );
void r82vec_permute ( int n, double a[], int p[] );
void r82vec_print ( int n, double a[], string title );
int *r82vec_sort_heap_index_a ( int n, double a[] );
void r82vec_sort_quick_a ( int n, double a[] );
void r8mat_copy ( int m, int n, double a1[], double a2[] );
double r8mat_det_2d ( double a[] );
double r8mat_det_3d ( double a[] );
double r8mat_det_4d ( double a[] );
double r8mat_det_5d ( double a[] );
double *r8mat_inverse_2d ( double a[] );
double *r8mat_inverse_3d ( double a[] );
double *r8mat_mv ( int m, int n, double a[], double x[] );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
int r8mat_solve ( int n, int rhs_num, double a[] );
double *r8mat_solve_2d ( double a[], double b[], double *det );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
double *r8mat_uniform_new ( int m, int n, double alo, double ahi, int &seed );
void r8mat_uniform_01 ( int m, int n, int &seed, double r[] );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
double r8vec_angle_3d ( double u[], double v[] );
double *r8vec_any_normal ( int dim_num, double v1[] );
void r8vec_bracket ( int n, double x[], double xval, int *left, int *right );
void r8vec_copy ( int n, double a1[], double a2[] );
double r8vec_cross_product_2d ( double v1[2], double v2[2] );
double r8vec_cross_product_affine_2d ( double v0[2], double v1[2], double v2[2] );
double *r8vec_cross_product_3d ( double v1[3], double v2[3] );
double *r8vec_cross_product_affine_3d ( double v0[3], double v1[3], double v2[3] );
double r8vec_distance ( int dim_num, double v1[], double v2[] );
double r8vec_dot_product ( int n, double a1[], double a2[] );
double r8vec_dot_product_affine ( int n, double v0[], double v1[], double v2[] );
bool r8vec_eq ( int n, double a1[], double a2[] );
bool r8vec_gt ( int n, double a1[], double a2[] );
bool r8vec_lt ( int n, double a1[], double a2[] );
double r8vec_max ( int n, double *rvec );
double r8vec_mean ( int n, double x[] );
double r8vec_min ( int n, double *rvec );
bool r8vec_negative_strict ( int n, double a[] );
double r8vec_norm ( int dim_num, double x[] );
double r8vec_norm_affine ( int n, double v0[], double v1[] );
double *r8vec_normal_01_new ( int n, int &seed );
double r8vec_normsq ( int n, double a[] );
double r8vec_normsq_affine ( int n, double v0[], double v1[] );
bool r8vec_positive_strict ( int n, double a[] );
void r8vec_print ( int n, double a[], string title );
void r8vec_print_2d ( double x, double y, string title );
void r8vec_print_3d ( double x, double y, double z, string title );
double r8vec_scalar_triple_product ( double v1[3], double v2[3], double v3[3] );
void r8vec_swap ( int n, double a1[], double a2[] );
double *r8vec_uniform_new ( int n, double alo, double ahi, int &seed );
double *r8vec_uniform_01_new ( int n, int &seed );
double *r8vec_uniform_unit_new ( int dim_num, int &seed );
double r8vec_variance ( int n, double x[] );
void r8vec_zero ( int n, double a1[] );
double radec_distance_3d ( double ra1, double dec1, double ra2, double dec2 );
double *radec_to_xyz ( double ra, double dec );
double radians_to_degrees ( double angle );
void radians_to_dms ( double radians, int *degrees, int *minutes, 
  int *seconds );
unsigned long random_initialize ( unsigned long seed );
void rotation_axis_vector_3d ( double axis[3], double angle, double v[3], 
  double w[3] );
void rotation_axis2mat_3d ( double axis[3], double angle, double a[3*3] );
void rotation_axis2quat_3d ( double axis[3], double angle, double q[4] );
void rotation_mat_vector_3d ( double a[3*3], double v[3], double w[3] );
void rotation_mat2axis_3d ( double a[3*3], double axis[3], double *angle );
void rotation_mat2quat_3d ( double a[3*3], double q[4] );
void rotation_quat_vector_3d ( double q[4], double v[3], double w[3] );
void rotation_quat2axis_3d ( double q[4], double axis[3], double *angle );
void rotation_quat2mat_3d ( double q[4], double a[3*3] );
void rtp_to_xyz ( double r, double theta, double phi, double xyz[3] );
int s_len_trim ( string s );
void segment_contains_point_1d ( double p1, double p2, double p3, double *u );
void segment_contains_point_2d ( double p1[2], double p2[2], double p3[2], 
  double u[2] );
void segment_point_coords_2d ( double p1[2], double p2[2], double p[2], 
  double *s, double *t );
void segment_point_coords_3d ( double p1[3], double p2[3], double p[3], 
  double *s, double *t );
double segment_point_dist_2d ( double p1[2], double p2[2], double p[2] );
double segment_point_dist_3d ( double p1[3], double p2[3], double p[3] );
void segment_point_near_2d ( double p1[2], double p2[2], double p[2], 
  double pn[2], double *dist, double *t );
void segment_point_near_3d ( double p1[3], double p2[3], double p[3],
  double pn[3], double *dist, double *t );
double segments_curvature_2d ( double p1[2], double p2[2], double p3[2] );
double segments_dist_2d ( double p1[2], double p2[2], double q1[2], 
  double q2[2] );
double segments_dist_3d ( double p1[3], double p2[3], double q1[3], 
  double q2[3] );
double segments_dist_3d_old ( double p1[3], double p2[3], double q1[3], 
  double q2[3] );
double segments_int_1d ( double p1, double p2, double q1, double q2, 
  double *r1, double *r2 );
void segments_int_2d ( double p1[2], double p2[2], double p3[2], 
  double p4[2], int *flag, double p5[2] );
double shape_point_dist_2d ( double pc[2], double p1[2], int nside, 
  double p[2] );
void shape_point_near_2d ( double pc[2], double p1[2], int nside, double p[2], 
  double pn[2], double *dist );
void shape_print_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[] );
void shape_ray_int_2d ( double pc[2], double p1[2], int nside, double pa[2], 
  double pb[2], double pi[2] );
void simplex_lattice_layer_point_next ( int n, int c[], int v[], bool *more );
void simplex_lattice_point_next ( int n, int c[], int v[], bool *more );
int simplex_unit_lattice_point_num_nd ( int d, int s );
double simplex_unit_volume_nd ( int ndim );
double simplex_volume_nd ( int ndim, double a[] );
double sin_power_int ( double a, double b, int n );
void soccer_shape_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[] );
void soccer_size_3d ( int *point_num, int *edge_num, int *face_num,
  int *face_order_max );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
double sphere_cap_area_2d ( double r, double h );
double sphere_cap_area_3d ( double r, double h );
double sphere_cap_area_nd ( int ndim, double r, double h );
double sphere_cap_volume_2d ( double r, double h );
double sphere_cap_volume_3d ( double r, double h );
double sphere_cap_volume_nd ( int ndim, double r, double h );
void sphere_dia2imp_3d ( double p1[3], double p2[3], double *r, double pc[3] );
double sphere_distance_xyz ( double xyz1[3], double xyz2[3] );
double sphere_distance1 ( double lat1, double long1, double lat2, 
  double long2, double radius );
double sphere_distance2 ( double lat1, double long1, double lat2, 
  double long2, double radius );
double sphere_distance3 ( double lat1, double long1, double lat2, 
  double long2, double radius );
bool sphere_exp_contains_point_3d ( double p1[3], double p2[3], double p3[3], 
  double p4[3], double p[3] );
void sphere_exp_point_near_3d ( double p1[3], double p2[3], double p3[3], 
  double p4[3], double p[3], double pn[3] );
void sphere_exp2imp_3d ( double p1[3], double p2[3], double p3[3], double p4[3], 
  double *r, double pc[3] );
void sphere_exp2imp_nd ( int n, double p[], double &r, double pc[] );
double sphere_imp_area_nd ( int n, double r );
bool sphere_imp_contains_point_3d ( double r, double pc[3], double p[3] );
void sphere_imp_grid_icos_size ( int factor, int *node_num,  int *edge_num,
  int *triangle_num );
void sphere_imp_gridfaces_3d ( int maxtri, int nlat, int nlong, int *ntri, 
  int tri[] );
int sphere_imp_line_project_3d ( double r, double pc[3], int n, double p[], 
  int maxpnt2, double pp[], double thetamin, double thetamax );
void sphere_imp_local2xyz_3d ( double r, double pc[3], double theta, 
  double phi, double p[3] );
void sphere_imp_point_near_3d ( double r, double pc[3], double p[3], 
  double pn[3] );
void sphere_imp_point_project_3d ( double r, double pc[3], double p[3], 
  double pp[3] );
double sphere_imp_volume_3d ( double r );
double sphere_imp_volume_nd ( int n, double r );
double sphere_imp_zone_area_3d ( double r, double h1, double h2 );
double sphere_imp_zone_volume_3d ( double r, double h1, double h2 );
void sphere_imp2exp_3d ( double r, double pc[3], double p1[3], double p2[3], 
  double p3[3], double p4[3] );
double sphere_k ( int n );
double sphere_triangle_angles_to_area ( double r, double a, double b, double c );
void sphere_triangle_sides_to_angles ( double r, double as, double bs, double cs, 
  double &a, double &b, double &c );
void sphere_triangle_vertices_to_angles ( double r, double v1[3], double v2[3], 
  double v3[3], double &a, double &b, double &c );
double sphere_triangle_vertices_to_area ( double r, double v1[3], double v2[3], 
  double v3[3] );
void sphere_triangle_vertices_to_centroid ( double r, double v1[3], double v2[3], 
  double v3[3], double vs[] );
int sphere_triangle_vertices_to_orientation ( double v1[3], double v2[3], double v3[3] );
void sphere_triangle_vertices_to_sides ( double r, double v1[3], double v2[3], 
  double v3[3], double &as, double &bs, double &cs );
double sphere_unit_area_nd ( int n );
void sphere_unit_area_values ( int &n_data, int &n, double &area );
double *sphere_unit_sample_2d ( int &seed );
double *sphere_unit_sample_3d ( int &seed );
double *sphere_unit_sample_3d_2 ( int &seed );
double *sphere_unit_sample_nd ( int n, int &seed );
double *sphere_unit_sample_nd_2 ( int n, int &seed );
double *sphere_unit_sample_nd_3 ( int n, int &seed );
double sphere_unit_volume_nd ( int n );
void sphere_unit_volume_values ( int &n_data, int &n, double &volume );
double sphere01_distance_xyz ( double xyz1[3], double xyz2[3] );
double sphere01_polygon_area ( int n, double lat[], double lon[] );
double sphere01_polygon_area_karney ( int n, double lat[], double lon[] );
double sphere01_triangle_angles_to_area ( double a, double b, double c );
void sphere01_triangle_sides_to_angles ( double as, double bs, double cs, 
  double &a, double &b, double &c );
void sphere01_triangle_vertices_to_angles ( double v1[3], double v2[3], 
  double v3[3], double &a, double &b, double &c );
double sphere01_triangle_vertices_to_area ( double v1[3], double v2[3], 
  double v3[3] );
void sphere01_triangle_vertices_to_midpoints ( double v1[3], double v2[3], 
  double v3[3], double m1[3], double m2[3], double m3[3] );
void sphere01_triangle_vertices_to_centroid ( double v1[3], double v2[3], 
  double v3[3], double vs[] );
void sphere01_triangle_vertices_to_sides ( double v1[3], double v2[3], 
  double v3[3], double &as, double &bs, double &cs );
void string_2d ( int vec_num, double p1[], double p2[], int *string_num, 
  int order[], int string[] );
void super_ellipse_points_2d ( double pc[2], double r1, double r2, 
  double expo, double psi, int n, double p[] );
double *tetrahedron_barycentric_3d ( double tetra[3*4], double p[3] );
double *tetrahedron_centroid_3d ( double tetra[3*4] );
void tetrahedron_circumsphere_3d ( double tetra[3*4], double *r, double pc[3] );
bool tetrahedron_contains_point_3d ( double tetra[3*4], double p[3] );
double *tetrahedron_dihedral_angles_3d ( double tetra[] );
double *tetrahedron_edge_length_3d ( double tetra[3*4] );
void tetrahedron_face_angles_3d ( double tetra[], double angles[] );
void tetrahedron_face_areas_3d ( double tetra[], double areas[] );
void tetrahedron_insphere_3d ( double tetra[3*4], double *r, double pc[3] );
void tetrahedron_lattice_layer_point_next ( int c[], int v[], bool *more );
void tetrahedron_lattice_point_next ( int c[], int v[], bool *more );
double tetrahedron_quality1_3d ( double tetra[3*4] );
double tetrahedron_quality2_3d ( double tetra[3*4] );
double tetrahedron_quality3_3d ( double tetra[3*4] );
double tetrahedron_quality4_3d ( double tetra[3*4] );
void tetrahedron_rhombic_shape_3d ( int point_num, int face_num, 
  int face_order_max, double point_coord[], int face_order[], 
  int face_point[] );
void tetrahedron_rhombic_size_3d ( int *point_num, int *edge_num, 
  int *face_num, int *face_order_max );
void tetrahedron_sample_3d ( double tetra[3*4], int n, int &seed, double p[] );
void tetrahedron_shape_3d ( int point_num, int face_num, int face_order_max, 
  double point_coord[], int face_order[], int face_point[] );
void tetrahedron_size_3d ( int *point_num, int *edge_num, int *face_num, 
  int *face_order_max );
double *tetrahedron_solid_angles_3d ( double tetra[] );
int tetrahedron_unit_lattice_point_num_3d ( int s );
double tetrahedron_volume_3d ( double tetra[3*4] );
void timestamp ( );
void tmat_init ( double a[4*4] );
void tmat_mxm ( double a[4*4], double b[4*4], double c[4*4] );
void tmat_mxp ( double a[4*4], double x[4], double y[4] );
void tmat_mxp2 ( double a[4*4], double p1[], double p2[], int n );
void tmat_mxv ( double a[4*4], double x[4], double y[4] );
void tmat_rot_axis ( double a[4*4], double b[4*4], double angle, 
  char axis );
void tmat_rot_vector ( double a[4*4], double b[4*4], double angle, 
  double v[3] );
void tmat_scale ( double a[4*4], double b[4*4], double s[3] );
void tmat_shear ( double a[4*4], double b[4*4], string axis, double s );
void tmat_trans ( double a[4*4], double b[4*4], double v[3] );
double torus_volume_3d ( double r1, double r2 );
double *tp_to_xyz ( double theta, double phi );
void triangle_angles_2d ( double t[2*3], double angle[3] );
double *triangle_angles_2d_new ( double t[2*3] );
void triangle_angles_3d ( double t[3*3], double angle[3] );
double *triangle_angles_3d_new ( double t[3*3] );
double triangle_area_2d ( double t[2*3] );
double triangle_area_3d ( double t[3*3] );
double triangle_area_3d_2 ( double t[3*3] );
double triangle_area_3d_3 ( double t[3*3] );
double triangle_area_heron ( double s[3] );
double *triangle_area_vector_3d ( double t[3*3] );
double *triangle_barycentric_2d ( double t[2*3], double p[2] );
double *triangle_centroid_2d ( double t[2*3] );
double *triangle_centroid_3d ( double t[3*3] );
double *triangle_circumcenter_2d ( double t[2*3] );
double *triangle_circumcenter_2d_2 ( double t[2*3] );
double *triangle_circumcenter ( int n, double t[] );
void triangle_circumcircle_2d ( double t[2*3], double *r, double pc[2] );
void triangle_circumcircle_2d_2 ( double t[2*3], double *r, double pc[2] );
double triangle_circumradius_2d ( double t[2*3] );
void triangle_contains_line_exp_3d ( double t[3*3], double p1[3], 
  double p2[3], bool *inside, double pint[3] );
void triangle_contains_line_par_3d ( double t[], double p0[], double pd[], 
  bool *inside, double p[] );
bool triangle_contains_point_2d_1 ( double t[2*3], double p[2] );
bool triangle_contains_point_2d_2 ( double t[2*3], double p[2] );
bool triangle_contains_point_2d_3 ( double t[2*3], double p[2] );
double triangle_diameter_2d ( double t[2*3] );
double *triangle_edge_length_2d ( double t[2*3] );
void triangle_gridpoints_2d ( double t[2*3], int sub_num, int grid_max, 
  int *grid_num, double p[] );
void triangle_incenter_2d ( double t[2*3], double pc[2] );
void triangle_incircle_2d ( double t[2*3], double pc[2], double *r );
double triangle_inradius_2d ( double t[2*3] );
bool triangle_is_degenerate_nd ( int dim_num, double t[] );
void triangle_lattice_layer_point_next ( int c[], int v[], bool *more );
void triangle_lattice_point_next ( int c[], int v[], bool *more );
void triangle_line_imp_int_2d ( double t[2*3], double a, double b, double c, 
  int *int_num, double pint[] );
int triangle_orientation_2d ( double t[2*3] );
void triangle_orthocenter_2d ( double t[2*3], double p[2], bool *flag );
double triangle_point_dist_2d ( double t[2*3], double p[2] );
double triangle_point_dist_3d ( double t[3*3], double p[3] );
double triangle_point_dist_signed_2d ( double t[2*3], double p[2] );
void triangle_point_near_2d ( double t[2*3], double p[2], double pn[2], 
  double *dist );
double triangle_quality_2d ( double t[2*3] );
int triangle_right_lattice_point_num_2d ( int a, int b );
void triangle_sample ( double t[2*3], int n, int &seed, double p[] );
int triangle_unit_lattice_point_num_2d ( int s );
void triangle_xsi_to_xy_2d ( double t[2*3], double xsi[3], double p[2] );
void triangle_xy_to_xsi_2d ( double t[2*3], double p[2], double xsi[3] );
void truncated_octahedron_shape_3d ( int point_num, int face_num, 
  int face_order_max, double point_coord[], int face_order[], 
  int face_point[] );
void truncated_octahedron_size_3d ( int *point_num, int *edge_num,
  int *face_num, int *face_order_max );
void tube_2d ( double dist, int n, double p[], double p1[], double p2[] );
void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank );
void vector_directions_nd ( int dim_num, double v[], double angle[] );
void vector_rotate_2d ( double v1[2], double angle, double v2[2] );
void vector_rotate_3d ( double p1[3], double pa[3], double angle, 
  double p2[3] );
void vector_rotate_base_2d ( double p1[2], double pb[2], double angle, 
  double p2[2] );
double vector_separation_2d ( double v1[], double v2[] );
double vector_separation_3d ( double v1[], double v2[] );
double vector_separation_nd ( int n, double v1[], double v2[] );
void vector_unit_nd ( int n, double p[] );
int voxels_dist_l1_3d ( int v1[3], int v2[3] );
int voxels_dist_l1_nd ( int dim_num, int v1[], int v2[] );
void voxels_line_3d ( int p1[3], int p2[3], int n, int p[] );
void voxels_region_3d ( int maxlist, int nx, int ny, int nz, int ishow[], 
  int *list_num, int list[], int *nregion );
void voxels_step_3d ( int v1[3], int v2[3], int inc, int jnc, int knc, 
  int v3[3] );
void xy_to_polar ( double xy[2], double *r, double *t );
void xyz_to_radec ( double p[3], double *ra, double *dec );
void xyz_to_rtp ( double xyz[3], double *r, double *theta, double *phi );
void xyz_to_tp ( double xyz[3], double *theta, double *phi );

#endif /*__GMS_GEOMETRY_HPP__*/
