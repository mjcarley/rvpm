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
#include "config.h"
#endif /*HAVE_CONFIG_H*/

#include <stdio.h>
#include <math.h>

#include <glib.h>

#include <rvpm.h>

#include "rvpm-private.h"

#ifdef RVPM_SINGLE_PRECISION
#define ELLIP_TOL 1e-7
#define EPS 1.0e-7
#define TINY 1.0e-32
#define NPRE 300
#else /*RVPM_SINGLE_PRECISION*/
#define ELLIP_TOL 1e-14
#define EPS 1.0e-16
#define TINY 1.0e-280
#define NPRE 300
#endif /*RVPM_SINGLE_PRECISION*/

static void carlson_R(RVPM_REAL x, RVPM_REAL y, RVPM_REAL *RF, RVPM_REAL *RG,
		      RVPM_REAL tol)

  /* 
   * Carlson's method for complete elliptic integrals
   * 
   * https://arxiv.org/pdf/math/9409227.pdf 
   */
  
{
  gint m ;
  RVPM_REAL c, tmp, two_pm ;
  
  tol = 2.7*SQRT(tol) ;
  x = SQRT(x) ; y = SQRT(y) ;
  
  two_pm = 0.25 ;
  c = 0.25*(x+y)*(x+y) ;
  for ( m = 0 ; m < 16 ; m ++ ) {
    tmp = x ;
    x = 0.5*(x+y) ;
    y = SQRT(y*tmp) ;
    two_pm *= 2.0 ;
    c -= two_pm*(x-y)*(x-y) ;
    if ( ABS(x-y) < tol*ABS(x) ) break ;
  }

  *RF = M_PI/(x+y) ;
  *RG = 0.5*c*(*RF) ;
  
  return ;
}

gint RVPM_FUNCTION_NAME(rvpm_elliptic_KE)(RVPM_REAL k,
					  RVPM_REAL *K, RVPM_REAL *E)

{
#ifdef RVPM_SINGLE_PRECISION
  RVPM_REAL tol = 1e-7 ;
#else /*RVPM_SINGLE_PRECISION*/
  RVPM_REAL tol = 1e-14 ;
#endif /*RVPM_SINGLE_PRECISION*/

  carlson_R(1.0 - k*k, 1.0, K, E, tol) ;
  *E *= 2.0 ;
  
  return 0 ;
}
