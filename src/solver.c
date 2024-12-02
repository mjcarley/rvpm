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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <blaswrap.h>

#include "rvpm.h"
#include "rvpm-private.h"

extern GTimer *timer ;

static gint time_step(rvpm_tree_t *tree, rvpm_solver_t *s,
		      RVPM_REAL *a, RVPM_REAL *b, gint nc,
		      RVPM_REAL t, RVPM_REAL dt,
		      RVPM_REAL *u, gint ustr,
		      RVPM_REAL *work)

{
  rvpm_distribution_t *d ;
  gint np, str, i, j ;
  RVPM_REAL *x ;

  g_assert(ustr >= RVPM_DISTRIBUTION_PARTICLE_SIZE) ;
  
  d = tree->d ;
  np = rvpm_distribution_particle_number(d) ;
  
  fprintf(stderr,
	  "%s: tree dimensions: %lg %lg %lg %lg\n", __FUNCTION__,
	  tree->origin[0], tree->origin[1], tree->origin[2], tree->D) ;

  for ( j = 0 ; j < nc ; j ++ ) {
    RVPM_FUNCTION_NAME(rvpm_tree_update)(tree, rvpm_solver_thread_number(s),
					 TRUE, work) ;
    for ( i = 0 ; i < RVPM_DISTRIBUTION_PARTICLE_SIZE ; i ++ ) {
#ifdef RVPM_SINGLE_PRECISION
      blaswrap_sscal(np, a[j], &(u[i]), ustr) ;
#else /*RVPM_SINGLE_PRECISION*/
      blaswrap_dscal(np, a[j], &(u[i]), ustr) ;
#endif /*RVPM_SINGLE_PRECISION*/
    }

    /*call out to auxiliary velocity and gradient calculation goes
      here*/
    RVPM_FUNCTION_NAME(rvpm_tree_derivatives)(tree, s, u, ustr, dt, work) ;
      
    x = (RVPM_REAL *)rvpm_distribution_particle(d,0) ;
    str = RVPM_DISTRIBUTION_PARTICLE_SIZE ;
    for ( i = 0 ; i < RVPM_DISTRIBUTION_PARTICLE_SIZE ; i ++ ) {
#ifdef RVPM_SINGLE_PRECISION
      blaswrap_saxpy(np, b[j], &(u[i]), ustr, &(x[i]), str) ;
#else /*RVPM_SINGLE_PRECISION*/
      blaswrap_daxpy(np, b[j], &(u[i]), ustr, &(x[i]), str) ;
#endif /*RVPM_SINGLE_PRECISION*/
    }

    /*update time with dt/dt = 1*/
    t += b[j]*(a[j]*t + dt*1.0) ;
  }
    
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_solver_solve)(rvpm_tree_t *tree,
					   rvpm_solver_t *s,
					   RVPM_REAL t, RVPM_REAL dt,
					   RVPM_REAL *u, gint ustr,
					   RVPM_REAL *work)

{
  RVPM_REAL a[4], b[4] ;
  gint nc, p ;

  RVPM_FUNCTION_NAME(rvpm_solver_coefficients)(rvpm_solver_time_step(s),
					       a, b, &nc, &p) ;
  if ( ustr < RVPM_DISTRIBUTION_PARTICLE_SIZE) {
    g_error("%s: derivative buffer stride (%d) must be at least %d",
	    __FUNCTION__, ustr, RVPM_DISTRIBUTION_PARTICLE_SIZE) ;
  }

  time_step(tree, s, a, b, nc, t, dt, u, ustr, work) ;
  
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_solver_coefficients)(rvpm_time_step_t s,
						  RVPM_REAL *a,
						  RVPM_REAL *b,
						  gint *n, gint *p)

{
  RVPM_REAL al[4], wt[4], bt ;
  
  switch ( s ) {
  default: g_error("%s: unhandled time step %u", __FUNCTION__, s) ;
  case RVPM_TIME_STEP_EULER:
    a[0] = 0.0 ; b[0] = 1.0 ;
    *n = 1 ; *p = 1 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_2:
    al[1] = 2.0/3.0 ;
    wt[0] = 1.0/4.0; wt[1] = 3.0/4.0 ;
    *n = 2 ; *p = 2 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_3:
    al[1] = 1.0/2.0 ;
    wt[0] = 0.0 ; wt[1] = 1.0 ;
    *n = 2 ; *p = 2 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_4:
    al[1] = 2.0/3.0 ; al[2] = 0.0 ; 
    bt = -3.0/4.0 ;
    wt[0] = 7.0/12.0 ; wt[1] = 3.0/4.0 ; wt[2] = -1.0/3.0 ;  
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_5:
    al[1] = 1.0/4.0 ; al[2] = 5.0/12.0 ;
    bt = 2.0/9.0 ;
    wt[0] = 1.0 ; wt[1] = -3.0 ; wt[2] = 3.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_6:
    al[1] = 1.0/4.0 ; al[2] = 2.0/3.0 ;
    bt = 8.0/9.0 ;
    wt[0] = 1.0/4.0 ; wt[1] = 0.0 ; wt[2] = 3.0/4.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_7:
    al[1] = 1.0/3.0 ; al[2] = 3.0/4.0 ;
    bt = 15.0/16.0 ;
    wt[0] = 1.0/6.0 ; wt[1] = 3.0/10.0 ; wt[2] = 8/15.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_8:
    al[1] = 3.0/5.0 - SQRT(6.0)/10.0 ;
    al[2] = 3.0/5.0 + SQRT(6.0)/15.0 ;
    bt = 2.0/3.0 + SQRT(6.0)/9.0 ;
    wt[0] = 1.0/6.0 ; wt[1] = 1.0/3.0 ; wt[2] = 1.0/2.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_12:
    al[1] = 2.0/3.0 ; al[2] = 2.0/3.0 ;
    bt = 3.0/4.0 ;
    wt[0] = 1.0/4.0 ; wt[1] = 5.0/12.0 ; wt[2] = 1.0/3.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_13:
    al[1] = 3.0/5.0 + SQRT(6.0)/10 ;
    al[2] = 3.0/5.0 - SQRT(6.0)/15 ;
    bt = 2.0/3.0 - SQRT(6.0)/9.0 ;
    wt[0] = 1.0/6.0 ; wt[1] = 1.0/3.0 ; wt[2] = 1.0/2.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_14:
    al[1] = 1.0 ; al[2] = 1.0/3.0 ;
    bt = 2.0/9.0 ;
    wt[0] = 0.0 ; wt[1] = 1.0/4.0 ; wt[2] = 3.0/4.0 ;
    *n = 3 ; *p = 3 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_17:
    al[1] = 1.0/2.0 ; al[2] = 1.0 ;
    bt = 1.0 ;
    wt[0] = 1.0/6.0 ; wt[1] = 2.0/3.0 ; wt[2] = 1.0/6.0 ;
    *n = 3 ; *p = 4 ;
    break ;
  case RVPM_TIME_STEP_WILLIAMSON_19:
    al[1] = 1.0 ; al[2] = 1.0/2.0 ;
    bt = 1.0/8.0 ;
    wt[0] = 1.0/6.0 ; wt[1] = 1.0/6.0 ; wt[2] = 2.0/3.0 ;
    *n = 3 ; *p = 4 ;
    break ;
  }

  if ( *n == 1 ) {
    return 0 ;
  }
  
  if ( *n == 2 ) {
    b[0] = al[1] ; b[1] = wt[1] ;
    a[0] = 0.0 ;
    a[1] = (wt[0] - b[0])/wt[1] ;

    return 0 ;
  }

  if ( *n == 3 ) {
    b[0] = al[1] ; b[1] = bt ; b[2] = wt[2] ;
    a[0] = 0.0 ;
    if ( wt[1] != 0 ) {
      a[1] = (wt[0] - b[0])/wt[1] ;
    } else {
      a[1] = (al[2] - bt - al[1])/b[1] ;
    }

    a[2] = (wt[1] - b[1])/wt[2] ;

    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

