/* This file is part of RVPM, a library for the Reformulated Vortex
 * Particle Method of Alvarez and Ning (see README.md and
 * documentation for references)
 *
 * Copyright (C) 2024 Michael Carley
 *
 * RVPM is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version. RVPM is distributed in the
 * hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RVPM.  If not, see <https://www.gnu.org/licenses/>.
 */

#define wbfmm_tree_point_index(_t,_i)			\
  ((RVPM_REAL *)(&((_t)->points[(_i)*((_t)->pstr)])))

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <wbfmm.h>

#include "rvpm.h"
#include "rvpm-private.h"

rvpm_tree_t *RVPM_FUNCTION_NAME(rvpm_tree_new)(rvpm_distribution_t *d,
					       gint depth,
					       gint order_max, RVPM_REAL *work)

{
  rvpm_tree_t *t ;
  RVPM_REAL xt[3], del, origin[3], D ;
  gsize pstr ;
  gint i ;

  if ( d->size != sizeof(RVPM_REAL) ) {
    g_error("%s: distribution precision (single or double) must match "
	    "tree precision", __FUNCTION__) ;
  }
  
  t = (rvpm_tree_t *)g_malloc0(sizeof(rvpm_tree_t)) ;
  t->d = d ;
  t->depth = depth ;
  t->order_max = order_max ;
  
  del = 1e-2 ;

  /*find bounding box for particles in d*/
  RVPM_FUNCTION_NAME(wbfmm_points_origin_width)((RVPM_REAL *)rvpm_distribution_particle(d,0),
						RVPM_DISTRIBUTION_PARTICLE_SIZE,
						rvpm_distribution_particle_number(d),
						origin, xt, &D, TRUE) ;
  origin[0] -= del ; origin[1] -= del ; origin[2] -= del ; D += 2.0*del ;
  t->origin[0] = origin[0] ;
  t->origin[1] = origin[1] ;
  t->origin[2] = origin[2] ;
  t->D = D ;

  /*allocate tree*/
  t->t = RVPM_FUNCTION_NAME(wbfmm_tree_new)(origin, D,
					    rvpm_distribution_particle_number_max(d)) ;

  /* /\*initialise the various operators*\/ */
  /* sizew = wbfmm_element_number_rotation(2*order_max) ; */
  /* sizew = MAX(sizew, (order_max+1)*(order_max+1)*3*16) ; */
  /* work = (RVPM_REAL *)g_malloc0(2*sizew*sizeof(RVPM_REAL)) ; */

  RVPM_FUNCTION_NAME(wbfmm_shift_angle_table_init)() ;
  t->shifts = RVPM_FUNCTION_NAME(wbfmm_shift_operators_new)(order_max,
							    TRUE, work) ;
  RVPM_FUNCTION_NAME(wbfmm_laplace_coaxial_translate_init)(order_max+2) ;

  /*add the points to the tree and refine*/
  pstr = RVPM_DISTRIBUTION_PARTICLE_SIZE*sizeof(RVPM_REAL) ;
  RVPM_FUNCTION_NAME(wbfmm_tree_sort_points)(t->t,
					     rvpm_distribution_particle(d,0), pstr,
					     rvpm_distribution_particle_number(d)) ;

  RVPM_FUNCTION_NAME(wbfmm_tree_add_points)(t->t, (gpointer)
					    (rvpm_distribution_particle(d,0)),
					    pstr, NULL, 0,
					    rvpm_distribution_particle_number(d), FALSE) ;
  for ( i = 0 ; i < depth ; i ++ ) RVPM_FUNCTION_NAME(wbfmm_tree_refine)(t->t) ;

  wbfmm_tree_problem(t->t) = WBFMM_PROBLEM_LAPLACE ;
  wbfmm_tree_source_size(t->t) = 3 ;
  for ( i = 1 ; i <= depth ; i ++ ) {
    RVPM_FUNCTION_NAME(wbfmm_tree_laplace_coefficient_init)(t->t, i, order_max, order_max) ;
  }  
  
  return t ;
}

/* static void check_indices(wbfmm_tree_t *t) */

/* { */
/*   gint i ; */
/*   RVPM_REAL *xi, *xim1 ; */
/*   guint64 idxi, idxim1 ; */
  
/*   for ( i = 1 ; i < wbfmm_tree_point_number(t) ; i ++ ) { */
/*     xim1 = wbfmm_tree_point_index(t, (i-1)) ; */
/*     xi   = wbfmm_tree_point_index(t, i  ) ; */
/*     idxi = wbfmm_point_index_3d(xi, wbfmm_tree_origin(t), wbfmm_tree_width(t)) ; */
/*     idxim1 = wbfmm_point_index_3d(xim1, wbfmm_tree_origin(t), wbfmm_tree_width(t)) ; */
/*     g_assert(idxi >= idxim1) ; */
/*   } */
  
/*   return ; */
/* }   */

/* static void check_boxes(wbfmm_tree_t *t, guint level) */

/* { */
/*   wbfmm_box_t *parents, *children ; */
/*   guint np, j, n ; */
/*   guint64 idx, child, xi, box ; */
/*   RVPM_REAL *x ; */

/*   /\* g_assert(t->problem != 0) ; *\/ */
/*   /\* wbfmm_tree_add_level(t) ; *\/ */

/*   /\*number of parent boxes to refine*\/ */
/*   np = 1 << 3*(level) ; */

/*   parents  = t->boxes[level  ] ; */
/*   /\* children = t->boxes[level+1] ; *\/ */

/*   /\*this could probably be done with binary searches*\/ */
/*   for ( idx = 0 ; idx < np ; idx ++ ) { */
/*     for ( j = parents[idx].i ; j < parents[idx].i + parents[idx].n ; j ++ ) { */
/*       x = wbfmm_tree_point_index(t, t->ip[j]) ; */
/*       xi = wbfmm_point_index_3d(x, wbfmm_tree_origin(t),  */
/* 				wbfmm_tree_width(t)) ; */
/*       box = wbfmm_point_locate_box(xi, level) ; */
/*       if ( box != idx ) { */
/* 	g_error("%uth point %u (%lg, %lg, %lg) in box %lu assigned to box %lu" */
/* 		" at level %u", */
/* 		j, t->ip[j], x[0], x[1], x[2], box, idx, level) ; */
/*       } */

/*     } */
/*   } */

/*   return ; */
/* } */

gint RVPM_FUNCTION_NAME(rvpm_tree_update)(rvpm_tree_t *t, gint nthreads,
					  gboolean init_limits, RVPM_REAL *work)

{
  RVPM_REAL xt[3], del, *x, origin[3], D ;
  rvpm_distribution_t *d ;
  gsize pstr, qstr ;
  gint i ;
  
  d = t->d ;
  /*find bounding box for particles in d*/
  del = 1e-2 ;
  if ( !init_limits ) {
    origin[0] = t->origin[0] ; xt[0] = origin[0] + t->D ;
    origin[1] = t->origin[1] ; xt[1] = origin[1] + t->D ;
    origin[2] = t->origin[2] ; xt[2] = origin[2] + t->D ;
  }
  RVPM_FUNCTION_NAME(wbfmm_points_origin_width)((RVPM_REAL *)rvpm_distribution_particle(d,0),
						RVPM_DISTRIBUTION_PARTICLE_SIZE,
						rvpm_distribution_particle_number(d),
						origin, xt, &D, init_limits) ;
  origin[0] -= del ; origin[1] -= del ; origin[2] -= del ; D += 2.0*del ;
  t->origin[0] = origin[0] ;
  t->origin[1] = origin[1] ;
  t->origin[2] = origin[2] ;
  t->D = D ;

  pstr = RVPM_DISTRIBUTION_PARTICLE_SIZE*sizeof(RVPM_REAL) ;
  qstr = RVPM_DISTRIBUTION_PARTICLE_SIZE ;

  x = wbfmm_tree_origin(t->t) ;
  x[0] = t->origin[0] ; x[1] = t->origin[1] ; x[2] = t->origin[2] ;
  wbfmm_tree_width(t->t) = t->D ;
  
  RVPM_FUNCTION_NAME(wbfmm_tree_sort_points)(t->t,
					     rvpm_distribution_particle(d,0), pstr,
					     rvpm_distribution_particle_number(d)) ;

  RVPM_FUNCTION_NAME(wbfmm_tree_add_points)(t->t, (gpointer)(rvpm_distribution_particle(d,0)),
			pstr, NULL, 0,
			rvpm_distribution_particle_number(d), FALSE) ;
  /* check_indices(t->t) ; */
  
  t->t->depth = 0 ;
  for ( i = 0 ; i < t->depth ; i ++ ) {
    /* check_boxes(t->t, i) ; */
    RVPM_FUNCTION_NAME(wbfmm_tree_refine)(t->t) ;
  }

  for ( i = 2 ; i <= t->t->depth ; i ++ ) {
    wbfmm_tree_coefficients_zero(t->t, i) ;
  }

  RVPM_FUNCTION_NAME(wbfmm_tree_laplace_leaf_expansions)(t->t,
							 (RVPM_REAL *)rvpm_distribution_vorticity(d,0), qstr,
				     NULL, 0, TRUE, work) ;

  for ( i = t->depth ; i >= 3 ; i -- ) {
    RVPM_FUNCTION_NAME(wbfmm_laplace_upward_pass)(t->t, t->shifts, i, work) ;
  }
  
  for ( i = 2 ; i <= t->depth ; i ++ ) {
    RVPM_FUNCTION_NAME(wbfmm_laplace_downward_pass)(t->t, t->shifts, i,
						    work, nthreads) ;
  }
  
  return 0 ;
}

static void box_curl_correct_WL(wbfmm_tree_t *t,gint i0, gint i1,
				RVPM_REAL *src, gint sstr,
				RVPM_REAL *sig, gint sigstr,
				RVPM_REAL *x, RVPM_REAL *f)

{
  gint idx, j ;
  RVPM_REAL *xs, K[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, xs, s, K, NULL) ;
    rvpm_vector_cross(u,K,&(src[idx*sstr])) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
  }

  return ;
}

static void box_curl_gradient_correct_WL(wbfmm_tree_t *t,
					 gint i0, gint i1,
					 RVPM_REAL *src, gint sstr,
					 RVPM_REAL *sig, gint sigstr,
					 RVPM_REAL *x,
					 RVPM_REAL *f, RVPM_REAL *df)

{
  gint idx, j ;
  RVPM_REAL *xs, K[12], u[12], *w, s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    w = &(src[idx*sstr]) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, xs, s, K, &(K[3])) ;
    rvpm_vector_cross(u,K,w) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
    rvpm_vector_cross_gradient(u,&(K[3]),w) ;
    df[0] += u[0] ; df[1] += u[1] ; df[2] += u[2] ; 
    df[3] += u[3] ; df[4] += u[4] ; df[5] += u[5] ; 
    df[6] += u[6] ; df[7] += u[7] ; df[8] += u[8] ; 
  }

  return ;
}

static void box_curl_correct_MR(wbfmm_tree_t *t,gint i0, gint i1,
				RVPM_REAL *src, gint sstr,
				RVPM_REAL *sig, gint sigstr,
				RVPM_REAL *x, RVPM_REAL *f)

{
  gint idx, j ;
  RVPM_REAL *xs, K[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    RVPM_FUNCTION_NAME(rvpm_kernel_MR)(x, xs, s, K, NULL) ;
    rvpm_vector_cross(u,K,&(src[idx*sstr])) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
  }

  return ;
}

static void box_curl_gradient_correct_MR(wbfmm_tree_t *t,
					 gint i0, gint i1,
					 RVPM_REAL *src, gint sstr,
					 RVPM_REAL *sig, gint sigstr,
					 RVPM_REAL *x,
					 RVPM_REAL *f, RVPM_REAL *df)

{
  gint idx, j ;
  RVPM_REAL *xs, K[12], u[12], *w, s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    w = &(src[idx*sstr]) ;
    s = sig[idx*sigstr] ;
    RVPM_FUNCTION_NAME(rvpm_kernel_MR)(x, xs, s, K, &(K[3])) ;
    rvpm_vector_cross(u,K,w) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
    rvpm_vector_cross_gradient(u,&(K[3]),w) ;
    df[0] += u[0] ; df[1] += u[1] ; df[2] += u[2] ; 
    df[3] += u[3] ; df[4] += u[4] ; df[5] += u[5] ; 
    df[6] += u[6] ; df[7] += u[7] ; df[8] += u[8] ; 
  }

  return ;
}

static void box_curl_correct_GS(wbfmm_tree_t *t,gint i0, gint i1,
				RVPM_REAL *src, gint sstr,
				RVPM_REAL *sig, gint sigstr,
				RVPM_REAL *x, RVPM_REAL *f)

{
  gint idx, j ;
  RVPM_REAL *xs, K[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    RVPM_FUNCTION_NAME(rvpm_kernel_GS)(x, xs, s, K, NULL) ;
    rvpm_vector_cross(u,K,&(src[idx*sstr])) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
  }

  return ;
}

static void box_curl_gradient_correct_GS(wbfmm_tree_t *t,gint i0, gint i1,
					 RVPM_REAL *src, gint sstr,
					 RVPM_REAL *sig, gint sigstr,
					 RVPM_REAL *x,
					 RVPM_REAL *f, RVPM_REAL *df)

{
  gint idx, j ;
  RVPM_REAL *xs, K[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    RVPM_FUNCTION_NAME(rvpm_kernel_GS)(x, xs, s, K, NULL) ;
    rvpm_vector_cross(u,K,&(src[idx*sstr])) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
  }

  return ;
}

/* static void box_curl_correct_sorted(wbfmm_tree_t *t, */
/* 				    gint i0, gint i1, */
/* 				    RVPM_REAL *src, gint sstr, */
/* 				    RVPM_REAL s, */
/* 				    RVPM_REAL *x, RVPM_REAL *f) */

/* { */
/*   gint idx ; */
/*   RVPM_REAL *xs, K[3], u[3] ; */

/*   for ( idx = i0 ; idx < i1 ; idx ++ ) { */
/*     xs = wbfmm_tree_point_index(t, idx) ; */
/*     RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, xs, s, K, NULL) ; */
/*     rvpm_vector_cross(u,K,&(src[idx*sstr])) ; */
/*     f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ;  */
/*   } */

/*   return ; */
/* } */

gint RVPM_FUNCTION_NAME(rvpm_tree_velocity_gradient)(rvpm_tree_t *t,
						     RVPM_REAL reg,
						     rvpm_kernel_t kernel,
						     RVPM_REAL *x, RVPM_REAL *u,
						     RVPM_REAL *du,
						     RVPM_REAL *work)

{
  rvpm_distribution_t *d ;
  wbfmm_box_t *boxes, *box ;
  guint64 b, neighbours[27] ;
  gint nnbr, i ;
  RVPM_REAL *src, utmp[12] = {0}, *sig ;
  gint sstr = RVPM_DISTRIBUTION_PARTICLE_SIZE, sigstr ;
  void (*correction_func)(wbfmm_tree_t *t,gint i0, gint i1,
			  RVPM_REAL *src, gint sstr,
			  RVPM_REAL *sig, gint sigstr,
			  RVPM_REAL *x, RVPM_REAL *f, RVPM_REAL *df) ;
  
  d = t->d ;
  b = RVPM_FUNCTION_NAME(wbfmm_point_box)(t->t, t->depth, x) ;
  g_assert(du != NULL) ;
  
  RVPM_FUNCTION_NAME(wbfmm_laplace_box_field)(t->t, t->depth, b,
					      (RVPM_REAL *)
					      rvpm_distribution_vorticity(d,0),
					      RVPM_DISTRIBUTION_PARTICLE_SIZE,
					      NULL, 0,
					      WBFMM_FIELD_CURL |
					      WBFMM_FIELD_GRADIENT,
					      FALSE, x, utmp, 1, work) ;
  u[0] = utmp[0] ; u[1] = utmp[1] ; u[2] = utmp[2] ;
  for ( i = 0 ; i < 9 ; i ++ ) du[i] = utmp[3+i] ;

  /*adjust for neighbours*/
  nnbr = wbfmm_box_neighbours(t->depth, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;
  boxes = t->t->boxes[t->depth] ;

  src = (RVPM_REAL *)rvpm_distribution_vorticity(d,0) ;
  g_assert(!(t->t->sorted)) ;
  switch ( kernel ) {
  default: g_assert_not_reached() ; break ;
  case RVPM_KERNEL_WINCKELMANS_LEONARD:
    correction_func = box_curl_gradient_correct_WL ;
    sig = &reg ; sigstr = 0 ; 
    break ;
  case RVPM_KERNEL_MOORE_ROSENHEAD:
    correction_func = box_curl_gradient_correct_MR ;
    sig = &reg ; sigstr = 0 ; 
    break ;
  case RVPM_KERNEL_GAUSSIAN:
    correction_func = box_curl_gradient_correct_GS ;
    sig = (RVPM_REAL *)rvpm_distribution_particle_radius(d,0) ;
    sigstr = sstr ; 
    break ;
  }
  
  for ( i = 0 ; i < nnbr ; i ++ ) {
    box = &(boxes[neighbours[i]]) ;
    correction_func(t->t, box->i, box->i+box->n, src, sstr, sig, sigstr, x,
		    u, du) ;
  }
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_tree_velocity)(rvpm_tree_t *t, RVPM_REAL reg,
					    rvpm_kernel_t kernel,
					    RVPM_REAL *x, RVPM_REAL *u,
					    RVPM_REAL *work)
 
{
  rvpm_distribution_t *d ;
  wbfmm_box_t *boxes, *box ;
  guint64 b, neighbours[27] ;
  gint nnbr, i ;
  RVPM_REAL *src, utmp[12] = {0}, *sig ;
  gint sstr = RVPM_DISTRIBUTION_PARTICLE_SIZE, sigstr ;
  void (*correction_func)(wbfmm_tree_t *t,gint i0, gint i1,
			  RVPM_REAL *src, gint sstr,
			  RVPM_REAL *sig, gint sigstr,
			  RVPM_REAL *x, RVPM_REAL *f) ;
  
  d = t->d ;
  b = RVPM_FUNCTION_NAME(wbfmm_point_box)(t->t, t->depth, x) ;
  
  RVPM_FUNCTION_NAME(wbfmm_laplace_box_field)(t->t, t->depth, b,
					      (RVPM_REAL *)rvpm_distribution_vorticity(d,0),
					      RVPM_DISTRIBUTION_PARTICLE_SIZE,
					      NULL, 0,
					      WBFMM_FIELD_CURL |
					      WBFMM_FIELD_GRADIENT,
					      FALSE, x, utmp, 1, work) ;
  u[0] += utmp[0] ; u[1] += utmp[1] ; u[2] += utmp[2] ;

  /*adjust for neighbours*/
  nnbr = wbfmm_box_neighbours(t->depth, b, neighbours) ;
  g_assert(nnbr >= 0 && nnbr < 28) ;
  boxes = t->t->boxes[t->depth] ;

  src = (RVPM_REAL *)rvpm_distribution_vorticity(d,0) ;
  g_assert(!(t->t->sorted)) ;
  
  switch ( kernel ) {
  default: g_assert_not_reached() ; break ;
  case RVPM_KERNEL_WINCKELMANS_LEONARD:
    correction_func = box_curl_correct_WL ;
    sig = &reg ; sigstr = 0 ; 
    break ;
  case RVPM_KERNEL_MOORE_ROSENHEAD:
    correction_func = box_curl_correct_MR ;
    sig = &reg ; sigstr = 0 ; 
    break ;
  case RVPM_KERNEL_GAUSSIAN:
    correction_func = box_curl_correct_GS ;
    sig = (RVPM_REAL *)rvpm_distribution_particle_radius(d,0) ;
    sigstr = sstr ; 
    break ;
  }

  for ( i = 0 ; i < nnbr ; i ++ ) {
    box = &(boxes[neighbours[i]]) ;
    correction_func(t->t, box->i, box->i+box->n, src, sstr, sig, sigstr, x, u) ;
  }
  
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_tree_velocity_self)(rvpm_tree_t *t, RVPM_REAL reg,
						 rvpm_kernel_t kernel,
						 RVPM_REAL *u, gint ustr,
						 RVPM_REAL *du, gint dustr,
						 RVPM_REAL *work)

{
  gint i ;

  if ( du == NULL ) {
    for ( i = 0 ; i < rvpm_distribution_particle_number(t->d) ; i ++ ) {
      RVPM_FUNCTION_NAME(rvpm_tree_velocity)(t, reg, kernel,
					     (RVPM_REAL *)
					     rvpm_distribution_particle(t->d,i),
					     &(u[i*ustr]), work) ;
    }

    return 0 ;
  }

  for ( i = 0 ; i < rvpm_distribution_particle_number(t->d) ; i ++ ) {
    RVPM_FUNCTION_NAME(rvpm_tree_velocity_gradient)(t, reg, kernel,
						    (RVPM_REAL *)rvpm_distribution_particle(t->d,i),
						    &(u[i*ustr]), &(du[i*dustr]), work) ;
  }

  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_tree_derivatives)(rvpm_tree_t *t,
					       rvpm_solver_t *solver,
					       RVPM_REAL *du, gint str,
					       RVPM_REAL al, 
					       RVPM_REAL *work)

/*
 * increment du by al*f(x) (compatible with low-storage Runge-Kutta)
 */
  
{
  gint i ;
  RVPM_REAL utmp[12], *p, *dG, *G, *s, f, g ;
  RVPM_REAL dG0[3], *ds, reg, dsigt ;
  rvpm_kernel_t kernel ;
  
  g_assert(str >= RVPM_DISTRIBUTION_PARTICLE_SIZE) ;

  kernel = rvpm_solver_kernel(solver) ;

  f   = rvpm_solver_model_parameter_f(solver) ;
  g   = rvpm_solver_model_parameter_g(solver) ;
  reg = rvpm_solver_regularisation(solver) ;
  
  for ( i = 0 ; i < rvpm_distribution_particle_number(t->d) ; i ++ ) {
    p = (RVPM_REAL *)rvpm_distribution_particle(t->d,i) ;
    s = (RVPM_REAL *)rvpm_distribution_particle_radius(t->d,i) ;
    RVPM_FUNCTION_NAME(rvpm_tree_velocity_gradient)(t, reg, kernel, p,
						    &(utmp[0]), &(utmp[3]),
						    work) ;
    /*particle velocity*/
    du[i*str+RVPM_DISTRIBUTION_PARTICLE_X+0] += al*utmp[0] ;
    du[i*str+RVPM_DISTRIBUTION_PARTICLE_X+1] += al*utmp[1] ;
    du[i*str+RVPM_DISTRIBUTION_PARTICLE_X+2] += al*utmp[2] ;
    /*vortex stretching*/
    G = (RVPM_REAL *)rvpm_distribution_vorticity(t->d,i) ;
    dG = &(du[i*str+RVPM_DISTRIBUTION_PARTICLE_VORTICITY]) ;
    ds = &(du[i*str+RVPM_DISTRIBUTION_PARTICLE_RADIUS]) ;

    RVPM_FUNCTION_NAME(rvpm_vorticity_derivatives)(G, *s, f, g,
						   &(utmp[3]), dG0, &dsigt) ;
    
    dG[0] += al*dG0[0] ;
    dG[1] += al*dG0[1] ;
    dG[2] += al*dG0[2] ;

    *ds += al*dsigt ;
  }

  return 0 ;
}
