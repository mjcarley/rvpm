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

#include <gqr.h>

#include "rvpm.h"
#include "rvpm-private.h"

static void stream_func(RVPM_REAL r, RVPM_REAL z,
			RVPM_REAL r1, RVPM_REAL z1,
			RVPM_REAL *G, RVPM_REAL *Gr, RVPM_REAL *Gz) 

{
  RVPM_REAL R1, R2, lm, lm_r, lm_z, K, E, dK, dE ;

  R1 = SQRT((r - r1)*(r - r1) + (z - z1)*(z - z1)) ;
  R2 = SQRT((r + r1)*(r + r1) + (z - z1)*(z - z1)) ;
  lm = (R2 - R1)/(R2 + R1) ;

  lm_r = 4.0*r1/R1/R2/(R1+R2)/(R1+R2)*(r1*r1 - r*r + (z-z1)*(z-z1)) ;
  lm_z = -2.0*lm*(z-z1)/R1/R2 ;

  RVPM_FUNCTION_NAME(rvpm_elliptic_KE)(lm, &K, &E) ;  

  dK = E/lm/(1.0-lm*lm) - K/lm ;
  dE = (E - K)/lm ;

  *G = (R1 + R2)*(K - E) ;

  *Gz = ((z-z1)/R1 + (z-z1)/R2)*(K - E) + (R1 + R2)*(dK - dE)*lm_z ;
  *Gr = ((r-r1)/R1 + (r+r1)/R2)*(K - E) + (R1 + R2)*(dK - dE)*lm_r ;

  *G /= 2.0*M_PI ;
  *Gz /= 2.0*M_PI ;
  *Gr /= 2.0*M_PI ;
  
  return ;
}
  
static void quad_limit_r(rvpm_stream_func_vorticity_t func,
			 gpointer data,
			 RVPM_REAL r, RVPM_REAL z, RVPM_REAL rmax,
			 RVPM_REAL *ur, RVPM_REAL *uz,
			 RVPM_REAL th0, RVPM_REAL th1,
			 gqr_rule_t *gs, gqr_rule_t *gth)

{
  gint i, j ;
  RVPM_REAL dth, thbar, s, s1, ds, sbar, th, w, G, Gr, Gz ;
  RVPM_REAL r1, z1 ;
  
  dth = 0.5*(th1 - th0) ; thbar = 0.5*(th1 + th0) ;

  for ( i = 0 ; i < gqr_rule_length(gth) ; i ++ ) {
    th = thbar + dth*gqr_rule_abscissa(gth, i) ;
    s1 = (r - rmax)/cos(th) ;
    ds = 0.5*(s1 - 0) ; sbar = 0.5*(s1 + 0) ;
    for ( j = 0 ; j < gqr_rule_length(gs) ; j ++ ) {
      s = sbar + ds*gqr_rule_abscissa(gs, j) ;
      r1 = r - s*cos(th) ; z1 = z + s*sin(th) ;
      w = func(r1, z1, data) ;
      stream_func(r, z, r1, z1, &G, &Gr, &Gz) ;
      *ur += -Gz*w*gqr_rule_weight(gs,j)*gqr_rule_weight(gth,i)*dth*s*ds ;
      *uz +=  Gr*w*gqr_rule_weight(gs,j)*gqr_rule_weight(gth,i)*dth*s*ds ;
      
    }
  }
  
  return ;
}

static void quad_limit_z(rvpm_stream_func_vorticity_t func,
			 gpointer data,
			 RVPM_REAL r, RVPM_REAL z, RVPM_REAL zmax,
			 RVPM_REAL *ur, RVPM_REAL *uz,
			 RVPM_REAL th0, RVPM_REAL th1,
			 gqr_rule_t *gs, gqr_rule_t *gth)

{
  gint i, j ;
  RVPM_REAL dth, thbar, s, s1, ds, sbar, th, w ;
  RVPM_REAL r1, z1, G, Gr, Gz ;

  /* fprintf(stderr, "%lg %lg\n", K, Kr) ; */
  
  dth = 0.5*(th1 - th0) ; thbar = 0.5*(th1 + th0) ;

  for ( i = 0 ; i < gqr_rule_length(gth) ; i ++ ) {
    th = thbar + dth*gqr_rule_abscissa(gth, i) ;
    s1 = (zmax - z)/sin(th) ;
    ds = 0.5*(s1 - 0) ; sbar = 0.5*(s1 + 0) ;
    for ( j = 0 ; j < gqr_rule_length(gs) ; j ++ ) {
      s = sbar + ds*gqr_rule_abscissa(gs, j) ;
      r1 = r - s*cos(th) ; z1 = z + s*sin(th) ;
      w = func(r1, z1, data) ;
      stream_func(r, z, r1, z1, &G, &Gr, &Gz) ;
      *ur += -Gz*w*gqr_rule_weight(gs,j)*gqr_rule_weight(gth,i)*dth*s*ds ;
      *uz +=  Gr*w*gqr_rule_weight(gs,j)*gqr_rule_weight(gth,i)*dth*s*ds ;
      /* /\* w = func(r - s*cos(th), z + s*sin(th), data) ; *\/ */

      /* *ur += gqr_rule_weight(gs,j)*gqr_rule_weight(gth,i)*dth*s*ds ; */
      
    }
  }
  
  return ;
}

gint RVPM_FUNCTION_NAME(rvpm_stream_func_velocity)(rvpm_stream_func_vorticity_t func,

						   gpointer data,
						   RVPM_REAL r, RVPM_REAL z,
						   RVPM_REAL rmax,
						   RVPM_REAL zmin,
						   RVPM_REAL zmax, 
						   RVPM_REAL *ur, RVPM_REAL *uz,
						   gqr_rule_t *gs,
						   gqr_rule_t *gth)
  
{
  gdouble th0, th1 ;

  *ur = *uz = 0.0 ;

  /*\rho_1 = \rho - s\cos\theta, z_1 = z + s\sin\theta*/
  th0 = 0 ; th1 = atan2(zmax - z, r - 0) ;
  quad_limit_r(func, data, r, z, 0, ur, uz, th0, th1, gs, gth) ;

  th0 = atan2(zmax - z, r - 0) ; th1 = atan2(zmax - z, r - rmax) ;
  quad_limit_z(func, data, r, z, zmax, ur, uz, th0, th1, gs, gth) ;

  th0 = atan2(zmax - z, r - rmax) ; th1 = M_PI ;
  quad_limit_r(func, data, r, z, rmax, ur, uz, th0, th1, gs, gth) ;

  th0 = M_PI ; th1 = atan2(zmin - z, r - rmax) + 2.0*M_PI ; 
  quad_limit_r(func, data, r, z, rmax, ur, uz, th0, th1, gs, gth) ;

  th0 = atan2(zmin - z, r - rmax) + 2.0*M_PI ;
  th1 = atan2(zmin - z, r - 0) + 2.0*M_PI ;
  quad_limit_z(func, data, r, z, zmin, ur, uz, th0, th1, gs, gth) ;

  th0 = atan2(zmin - z, r - 0) + 2.0*M_PI ; th1 = 2.0*M_PI ; 
  quad_limit_r(func, data, r, z, 0, ur, uz, th0, th1, gs, gth) ;

  *ur /= r ; *uz /= r ;
  
  return 0 ;
}
