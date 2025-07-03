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
    RVPM_FUNCTION_NAME(wbfmm_tree_laplace_coefficient_init)(t->t, i, order_max,
							    order_max) ;
    order_max -= 2 ;
    if ( order_max <= 0 ) {
      g_error("%s: expansion order (%d) out of range at level %d\n",
	      __FUNCTION__, order_max, i) ;
    }
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

static void kernel_bare(RVPM_REAL *x, RVPM_REAL *y, RVPM_REAL *K)

{
  RVPM_REAL R, R2, R3, r[3] ;
  
  memset( K, 0, 3*sizeof(RVPM_REAL)) ;
  
  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) ;
  R  = SQRT(R2) ;

  if ( R < WBFMM_LOCAL_CUTOFF_RADIUS ) return ;
  /* if ( R < 1e-9 ) return ; */

  R3 = R2*R ;

  K[0] = -r[0]/R3*0.25*M_1_PI ;
  K[1] = -r[1]/R3*0.25*M_1_PI ;
  K[2] = -r[2]/R3*0.25*M_1_PI ;
    
  return ;
}

static void correct_kernel(RVPM_REAL *K, RVPM_REAL *K0)

{
  K[0] -= K0[0] ; K[1] -= K0[1] ; K[2] -= K0[2] ; 

  return ;
}

static void kernel_gradient_bare(RVPM_REAL *x, RVPM_REAL *y,
				 RVPM_REAL *K, RVPM_REAL *dK)

{
  RVPM_REAL R, R2, R3, R5, r[3] ;
  
  memset( K, 0, 3*sizeof(RVPM_REAL)) ;
  memset(dK, 0, 9*sizeof(RVPM_REAL)) ;
  
  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) ;
  R  = SQRT(R2) ;

  if ( R < WBFMM_LOCAL_CUTOFF_RADIUS ) return ;
  /* if ( R < 1e-9 ) return ; */

  R3 = R2*R ;

  K[0] = -r[0]/R3*0.25*M_1_PI ;
  K[1] = -r[1]/R3*0.25*M_1_PI ;
  K[2] = -r[2]/R3*0.25*M_1_PI ;
  
  R5 = R3*R2 ;
  
  /*dK/dx*/
  dK[0] = -3.0*r[0]*r[0]/R5 + 1.0/R3 ;
  dK[1] = -3.0*r[1]*r[0]/R5 ;
  dK[2] = -3.0*r[2]*r[0]/R5 ;
  /*dK/dy*/
  /* dK[3] = -3.0*r[0]*r[1]/R5 ; */
  dK[3] = dK[1] ;
  dK[4] = -3.0*r[1]*r[1]/R5 + 1.0/R3 ;
  dK[5] = -3.0*r[2]*r[1]/R5 ;
  /*dK/dz*/
  /* dK[6] = -3.0*r[0]*r[2]/R5 ; */
  dK[6] = dK[2] ;
  /* dK[7] = -3.0*r[1]*r[2]/R5 ; */
  dK[7] = dK[5] ;
  dK[8] = -3.0*r[2]*r[2]/R5 + 1.0/R3 ;

  for ( gint i = 0 ; i < 9 ; i ++ ) {
    dK[i] *= -0.25*M_1_PI ;
  }
    
  return ;
}

static void correct_gradient_kernel(RVPM_REAL *K, RVPM_REAL *K0,
				    RVPM_REAL *dK, RVPM_REAL *dK0)

{
  K[0] -= K0[0] ; K[1] -= K0[1] ; K[2] -= K0[2] ; 

  dK[0] -= dK0[0] ; dK[1] -= dK0[1] ; dK[2] -= dK0[2] ; 
  dK[3] -= dK0[3] ; dK[4] -= dK0[4] ; dK[5] -= dK0[5] ; 
  dK[6] -= dK0[6] ; dK[7] -= dK0[7] ; dK[8] -= dK0[8] ; 

  return ;
}

static void box_curl_correct_WL(wbfmm_tree_t *t,gint i0, gint i1,
				RVPM_REAL *src, gint sstr,
				RVPM_REAL *sig, gint sigstr,
				RVPM_REAL *x, RVPM_REAL *f)

{
  gint idx, j ;
  RVPM_REAL *xs, K[3], K0[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    kernel_bare(x, xs, K0) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, xs, s, K, NULL) ;
    correct_kernel(K, K0) ;
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
  RVPM_REAL *xs, K0[12], K[12], u[12], *w, s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    w = &(src[idx*sstr]) ;
    kernel_gradient_bare(x, xs, K0, &(K0[3])) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, xs, s, K, &(K[3])) ;
    correct_gradient_kernel(K, K0, &(K[3]), &(K0[3])) ;
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
  RVPM_REAL *xs, K[3], K0[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    kernel_bare(x, xs, K0) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_MR)(x, xs, s, K, NULL) ;
    correct_kernel(K, K0) ;
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
  RVPM_REAL *xs, K0[12], K[12], u[12], *w, s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    w = &(src[idx*sstr]) ;
    s = sig[idx*sigstr] ;
    kernel_gradient_bare(x, xs, K0, &(K0[3])) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_MR)(x, xs, s, K, &(K[3])) ;
    correct_gradient_kernel(K, K0, &(K[3]), &(K0[3])) ;
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
  RVPM_REAL *xs, K[3], K0[3], u[3], s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    s = sig[idx*sigstr] ;
    kernel_bare(x, xs, K0) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_GS)(x, xs, s, K, NULL) ;
    correct_kernel(K, K0) ;
    rvpm_vector_cross(u,K,&(src[idx*sstr])) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
  }

  return ;
}

#if 0
static void kernel_correction_GS(RVPM_REAL *x, RVPM_REAL *y,
				 RVPM_REAL s,
				 RVPM_REAL *K, RVPM_REAL *dK)

{
  RVPM_REAL R, R2, R3, R5, r[3], k[3], dk[9], dg[3], g, cutoff, E, errfunc ;
  RVPM_REAL s3, Kb[3], dKb[9] ;
  const RVPM_REAL m_1_sqrtpi_3 = 0.5*M_1_PI*M_2_SQRTPI ;
  const RVPM_REAL m_1_4pi = 0.25*M_1_PI ;
  gint i, j ;

  /*argument at which g(R/\sigma) reaches machine precision*/
#ifdef RVPM_SINGLE_PRECISION
  cutoff = 4.2 ;
#else /*RVPM_SINGLE_PRECISION*/
  cutoff = 6.2 ;
#endif /*RVPM_SINGLE_PRECISION*/
  
  memset( K, 0, 3*sizeof(RVPM_REAL)) ;
  memset(dK, 0, 9*sizeof(RVPM_REAL)) ;
  
  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) ;
  R  = SQRT(R2) ;

  if ( R > cutoff*s ) return  ;
  
  s3 = s*s*s ;
  if ( R2 < 1e-12 ) {
    /*use series expansion at small R*/
    E = EXP(-R2/s/s) ;
    K[0] = -r[0]*E*m_1_sqrtpi_3/s3/3.0 ;
    K[1] = -r[1]*E*m_1_sqrtpi_3/s3/3.0 ;
    K[2] = -r[2]*E*m_1_sqrtpi_3/s3/3.0 ;

    dK[0] = -m_1_sqrtpi_3/3.0/s3 + 2.0*E/3.0*m_1_sqrtpi_3/s3/s/s*r[0]*r[0] ;
    dK[4] = -m_1_sqrtpi_3/3.0/s3 + 2.0*E/3.0*m_1_sqrtpi_3/s3/s/s*r[1]*r[1] ;
    dK[8] = -m_1_sqrtpi_3/3.0/s3 + 2.0*E/3.0*m_1_sqrtpi_3/s3/s/s*r[2]*r[2] ;

    return ;
  }
  
  R3 = R*R2 ;

  Kb[0] = -r[0]/R3*0.25*M_1_PI ;
  Kb[1] = -r[1]/R3*0.25*M_1_PI ;
  Kb[2] = -r[2]/R3*0.25*M_1_PI ;
  
  E = EXP(-R2/s/s) ;
  errfunc = ERF(R/s) ;
  g = errfunc - M_2_SQRTPI*R/s*E ;

  k[0] = -g*r[0]/R3*m_1_4pi ;
  k[1] = -g*r[1]/R3*m_1_4pi ;
  k[2] = -g*r[2]/R3*m_1_4pi ;

  K[0] = g*k[0] - Kb[0] ;
  K[1] = g*k[1] - Kb[1] ;
  K[2] = g*k[2] - Kb[2] ;
  
  R5 = R3*R2 ;
  
  /*dK/dx*/
  dKb[0] = -3.0*r[0]*r[0]/R5 + 1.0/R3 ;
  dKb[1] = -3.0*r[1]*r[0]/R5 ;
  dKb[2] = -3.0*r[2]*r[0]/R5 ;
  /*dK/dy*/
  /* dK[3] = -3.0*r[0]*r[1]/R5 ; */
  dKb[3] = dKb[1] ;
  dKb[4] = -3.0*r[1]*r[1]/R5 + 1.0/R3 ;
  dKb[5] = -3.0*r[2]*r[1]/R5 ;
  /*dK/dz*/
  /* dK[6] = -3.0*r[0]*r[2]/R5 ; */
  dKb[6] = dKb[2] ;
  /* dK[7] = -3.0*r[1]*r[2]/R5 ; */
  dKb[7] = dKb[5] ;
  dKb[8] = -3.0*r[2]*r[2]/R5 + 1.0/R3 ;

  for ( gint i = 0 ; i < 9 ; i ++ ) {
    dKb[i] *= -0.25*M_1_PI ;
  }

  /*dK/dx*/
  dk[0] = -m_1_4pi*(-3.0*r[0]*r[0]/R5 + 1.0/R3) ;
  dk[1] = -m_1_4pi*(-3.0*r[1]*r[0]/R5) ;
  dk[2] = -m_1_4pi*(-3.0*r[2]*r[0]/R5) ;
  /*dk/dy*/
  /* dk[3] = -3.0*r[0]*r[1]/R5 ; */
  dk[3] = dk[1] ;
  dk[4] = -m_1_4pi*(-3.0*r[1]*r[1]/R5 + 1.0/R3) ;
  dk[5] = -m_1_4pi*(-3.0*r[2]*r[1]/R5) ;
  /*dk/dz*/
  /* dk[6] = -3.0*r[0]*r[2]/R5 ; */
  dk[6] = dk[2] ;
  /* dk[7] = -3.0*r[1]*r[2]/R5 ; */
  dk[7] = dk[5] ;
  dk[8] = -m_1_4pi*(-3.0*r[2]*r[2]/R5 + 1.0/R3) ;

  dg[0] = 2.0*M_2_SQRTPI*R*E*r[0]/s3 ;
  dg[1] = 2.0*M_2_SQRTPI*R*E*r[1]/s3 ;
  dg[2] = 2.0*M_2_SQRTPI*R*E*r[2]/s3 ;

  for ( i = 0 ; i < 3 ; i ++ ) {
    for ( j = 0 ; j < 3 ; j ++ ) {
      dK[3*i+j] = g*dk[3*i+j] + k[j]*dg[i] - dKb[3*i+j] ;
    }
  }

  return ;
}
#endif

static void box_curl_gradient_correct_GS(wbfmm_tree_t *t,gint i0, gint i1,
					 RVPM_REAL *src, gint sstr,
					 RVPM_REAL *sig, gint sigstr,
					 RVPM_REAL *x,
					 RVPM_REAL *f, RVPM_REAL *df)

{
  gint idx, j ;
  RVPM_REAL *xs, K0[12], K[12], u[12], *w, s ;

  for ( j = i0 ; j < i1 ; j ++ ) {
    idx = t->ip[j] ;
    xs = wbfmm_tree_point_index(t, idx) ;
    w = &(src[idx*sstr]) ;
    s = sig[idx*sigstr] ;
    kernel_gradient_bare(x, xs, K0, &(K0[3])) ;
    RVPM_FUNCTION_NAME(rvpm_kernel_GS)(x, xs, s, K, &(K[3])) ;
    correct_gradient_kernel(K, K0, &(K[3]), &(K0[3])) ;
    /* kernel_correction_GS(x, xs, s, K, &(K[3])) ; */
    rvpm_vector_cross(u,K,w) ;
    f[0] += u[0] ; f[1] += u[1] ; f[2] += u[2] ; 
    rvpm_vector_cross_gradient(u,&(K[3]),w) ;
    df[0] += u[0] ; df[1] += u[1] ; df[2] += u[2] ;
    df[3] += u[3] ; df[4] += u[4] ; df[5] += u[5] ;
    df[6] += u[6] ; df[7] += u[7] ; df[8] += u[8] ;
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
					      TRUE, x, utmp, 1, work) ;
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
    /* g_assert_not_reached() ; */
    correction_func = box_curl_gradient_correct_WL ;
    sig = &reg ; sigstr = 0 ; 
    break ;
  case RVPM_KERNEL_MOORE_ROSENHEAD:
    correction_func = box_curl_gradient_correct_MR ;
    sig = &reg ; sigstr = 0 ; 
    break ;
  case RVPM_KERNEL_GAUSSIAN:
    /* correction_func = box_curl_gradient_correct_GS ; */
    sig = (RVPM_REAL *)rvpm_distribution_particle_radius(d,0) ;
    sigstr = sstr ;
    for ( i = 0 ; i < nnbr ; i ++ ) {
      box = &(boxes[neighbours[i]]) ;
      box_curl_gradient_correct_GS(t->t, box->i, box->i+box->n, src, sstr,
				   sig, sigstr, x, u, du) ;
    }
    return 0 ;
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
					      TRUE, x, utmp, 1, work) ;
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
    /* g_assert_not_reached() ; */
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
  RVPM_REAL dG0[3], *ds, reg, dsigt, nu ;
  rvpm_kernel_t kernel ;
  
  g_assert(str >= RVPM_DISTRIBUTION_PARTICLE_SIZE) ;

  kernel = rvpm_solver_kernel(solver) ;

  f   = rvpm_solver_model_parameter_f(solver) ;
  g   = rvpm_solver_model_parameter_g(solver) ;
  reg = rvpm_solver_regularisation(solver) ;
  nu  = rvpm_solver_viscosity(solver) ;
  
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
    ds = &(du[i*str+RVPM_DISTRIBUTION_PARTICLE_RADIUS   ]) ;

    RVPM_FUNCTION_NAME(rvpm_vorticity_derivatives)(G, *s, f, g, nu,
						   &(utmp[3]), dG0, &dsigt) ;
    
    dG[0] += al*dG0[0] ;
    dG[1] += al*dG0[1] ;
    dG[2] += al*dG0[2] ;

    *ds += al*dsigt ;
  }

  return 0 ;
}
