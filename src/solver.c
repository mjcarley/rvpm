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

#include "rvpm.h"
#include "rvpm-private.h"

extern GTimer *timer ;

static gint time_step_classical_euler(rvpm_tree_t *t,
				      rvpm_solver_t *s,
				      RVPM_REAL dt, RVPM_REAL *u, gint ustr,
				      RVPM_REAL *work)

{
  rvpm_distribution_t *d ;
  gint np, i ;
  RVPM_REAL *x, *G, dG[3], *du, reg ;
  
  d = t->d ;
  np = rvpm_distribution_particle_number(d) ;
  reg = rvpm_solver_regularisation(s) ;
  
  memset(u, 0, ustr*np*sizeof(RVPM_REAL)) ;
  RVPM_FUNCTION_NAME(rvpm_tree_update)(t, rvpm_solver_thread_number(s), work) ;

  fprintf(stderr,
	  "%s: tree dimensions: %lg %lg %lg %lg\n", __FUNCTION__,
	  t->origin[0], t->origin[1], t->origin[2], t->D) ;

  /*velocity and velocity gradient*/
  RVPM_FUNCTION_NAME(rvpm_tree_velocity_self)(t, reg, rvpm_solver_kernel(s),
					      u, ustr, &(u[3]), ustr, work) ;

  /*equation 34 of Alvarez and Ning 2024 (cVPM)*/
  for ( i = 0 ; i < np ; i ++ ) {
    x = (RVPM_REAL *)rvpm_distribution_particle(d,i) ;
    G = (RVPM_REAL *)rvpm_distribution_vorticity(d,i) ;
    /*dx/dt*/
    x[0] += u[ustr*i+0]*dt ;
    x[1] += u[ustr*i+1]*dt ;
    x[2] += u[ustr*i+2]*dt ;
    /*convenience pointer for \nabla u*/
    du = &(u[ustr*i+3]) ;
    /*d\Gamma/dt = \Gamma_{x}\partial u/\partial x + ...*/
    dG[0] = G[0]*du[0] + G[1]*du[3] + G[2]*du[6] ;
    dG[1] = G[0]*du[1] + G[1]*du[4] + G[2]*du[7] ;
    dG[2] = G[0]*du[2] + G[1]*du[5] + G[2]*du[8] ;
    G[0] += dG[0]*dt ; 
    G[1] += dG[1]*dt ; 
    G[2] += dG[2]*dt ; 
  }

  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_solver_solve)(rvpm_tree_t *t, rvpm_solver_t *s,
					   RVPM_REAL dt,
					   RVPM_REAL *u, gint ustr,
					   RVPM_REAL *work)

{
  if ( ustr < 12 ) {
    g_error("%s: velocity stride (%d) must be at least 12",
	    __FUNCTION__, ustr) ;
  }

  if ( rvpm_solver_method(s) == RVPM_METHOD_CLASSICAL ) {
    if ( rvpm_solver_time_step(s) == RVPM_TIME_STEP_EULER ) {
      time_step_classical_euler(t, s, dt, u, ustr, work) ;

      return 0 ;
    }
  }

  g_assert_not_reached() ;
  
  return 0 ;
}
