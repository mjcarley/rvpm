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

rvpm_distribution_t *RVPM_FUNCTION_NAME(rvpm_distribution_alloc)(gint np)

{
  rvpm_distribution_t *d ;

  d = (rvpm_distribution_t *)g_malloc0(sizeof(rvpm_distribution_t)) ;

  d->size = sizeof(RVPM_REAL) ;
  
  d->data = (char *)g_malloc0(np*RVPM_DISTRIBUTION_PARTICLE_SIZE*
			      sizeof(RVPM_REAL)) ;
  rvpm_distribution_particle_number(d) = 0 ;
  rvpm_distribution_particle_number_max(d) = np ;
  
  return d ;
}

gint RVPM_FUNCTION_NAME(rvpm_distribution_interleaved_to_contiguous)(rvpm_distribution_t *d)

/*
 * in-place permutation using method and code from:
 https://stackoverflow.com/questions/16501424/algorithm-to-apply-permutation-in-constant-memory-space/41472796#41472796
*/
  
{
  gint i, idx, np, ne ;
  RVPM_REAL tmp, *data ;

  np = rvpm_distribution_particle_number(d) ;
  ne = RVPM_DISTRIBUTION_PARTICLE_SIZE ;
  data = (RVPM_REAL *)(d->data) ;
  for ( i = 0 ; i < np*ne - 1 ; i ++ ) {
    /*convert from contiguous index to interleaved*/
    idx = (i - (i/np)*np)*ne + i/np ;
    while ( idx < i ) {
      idx = (idx - (idx/np)*np)*ne + idx/np ;
    }

    tmp = data[i] ;
    data[i] = data[idx] ;
    data[idx] = tmp ;
  }

  rvpm_distribution_interleaved(d) = FALSE ;
  
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_distribution_contiguous_to_interleaved)(rvpm_distribution_t *d)

/*
 * in-place permutation using method and code from:
 https://stackoverflow.com/questions/16501424/algorithm-to-apply-permutation-in-constant-memory-space/41472796#41472796
*/
  
{
  gint i, j, k, idx, np, ne ;
  RVPM_REAL tmp, *data ;

  np = rvpm_distribution_particle_number(d) ;
  ne = RVPM_DISTRIBUTION_PARTICLE_SIZE ;
  data = (RVPM_REAL *)d->data ;
  for ( i = 0 ; i < np*ne - 1 ; i ++ ) {
    /*convert from interleaved index to contiguous*/
    /* idx = (i - (i/np)*np)*ne + i/np ; */

    /*point index*/
    j = i/ne ;
    /*element index*/
    k = i-j*ne ;
    idx = k*np + j ;
    
    while ( idx < i ) {
      j = idx/ne ;
      /*element index*/
      k = idx-j*ne ;
      idx = k*np + j ;
    }

    tmp = data[i] ;
    data[i] = data[idx] ;
    data[idx] = tmp ;
  }

  rvpm_distribution_interleaved(d) = FALSE ;
  
  return 0 ;
}

rvpm_distribution_t *RVPM_FUNCTION_NAME(rvpm_distribution_read_alloc)(FILE *f,
								      gint n)

{
  rvpm_distribution_t *d ;
  RVPM_REAL limits[6], smax, *data ;
  gdouble *c ;
  gint np, i ;

  g_assert(n >= 0) ;
  
  fscanf(f, "%d", &np) ;
#ifdef RVPM_SINGLE_PRECISION
  for ( i = 0 ; i < 6 ; i ++ ) fscanf(f, "%g", &(limits[i])) ;
  fscanf(f, "%g", &smax) ;
#else /*RVPM_SINGLE_PRECISION*/
  for ( i = 0 ; i < 6 ; i ++ ) fscanf(f, "%lg", &(limits[i])) ;
  fscanf(f, "%lg", &smax) ;
#endif /*RVPM_SINGLE_PRECISION*/

  d = RVPM_FUNCTION_NAME(rvpm_distribution_alloc)(np+n) ;
  data = (RVPM_REAL *)rvpm_distribution_particle(d, 0) ;
  
  for ( i = 0 ; i < np*RVPM_DISTRIBUTION_PARTICLE_SIZE ; i ++ ) {
#ifdef RVPM_SINGLE_PRECISION
    fscanf(f, "%g", &(data[i])) ;
#else /*RVPM_SINGLE_PRECISION*/
    fscanf(f, "%lg", &(data[i])) ;
#endif /*RVPM_SINGLE_PRECISION*/
  }    

  rvpm_distribution_particle_number(d) = np ;

  c = rvpm_distribution_origin(d) ;
  c[0] = limits[0] ; c[1] = limits[2] ; c[2] = limits[4] ;
  rvpm_distribution_width(d) = MAX(limits[1] - limits[0],
				   MAX(limits[3] - limits[2],
				       limits[5] - limits[4])) ;
  
  return d ;
}

gint RVPM_FUNCTION_NAME(rvpm_distribution_read)(FILE *f, rvpm_distribution_t *d)

{
  RVPM_REAL limits[6], smax, *data ;
  gdouble *c ;
  gint np, i ;
  
  fscanf(f, "%d", &np) ;
  if ( np > rvpm_distribution_particle_number_max(d) ) {
    g_error("%s: not enough memory allocated (%d) for %d particles",
	    __FUNCTION__, rvpm_distribution_particle_number_max(d), np) ;
  }
    
#ifdef RVPM_SINGLE_PRECISION
  for ( i = 0 ; i < 6 ; i ++ ) fscanf(f, "%g", &(limits[i])) ;
  fscanf(f, "%g", &smax) ;
#else /*RVPM_SINGLE_PRECISION*/
  for ( i = 0 ; i < 6 ; i ++ ) fscanf(f, "%lg", &(limits[i])) ;
  fscanf(f, "%lg", &smax) ;
#endif /*RVPM_SINGLE_PRECISION*/

  data = (RVPM_REAL *)rvpm_distribution_particle(d, 0) ;
  
  for ( i = 0 ; i < np*RVPM_DISTRIBUTION_PARTICLE_SIZE ; i ++ ) {
#ifdef RVPM_SINGLE_PRECISION
    fscanf(f, "%g", &(data[i])) ;
#else /*RVPM_SINGLE_PRECISION*/
    fscanf(f, "%lg", &(data[i])) ;
#endif /*RVPM_SINGLE_PRECISION*/
  }    

  rvpm_distribution_particle_number(d) = np ;

  c = rvpm_distribution_origin(d) ;
  c[0] = limits[0] ; c[1] = limits[2] ; c[2] = limits[4] ;
  rvpm_distribution_width(d) = MAX(limits[1] - limits[0],
				   MAX(limits[3] - limits[2],
				       limits[5] - limits[4])) ;
  
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_distribution_particle_add)(rvpm_distribution_t *d,
							RVPM_REAL *x,
							RVPM_REAL *w,
							RVPM_REAL  s,
							RVPM_REAL G)
{
  gint np ;
  RVPM_REAL *p ;
  
  if ( rvpm_distribution_particle_number(d) >=
       rvpm_distribution_particle_number_max(d) ) {
    g_error("%s: not enough space allocated (%d) for %d particles",
	    __FUNCTION__,
	    rvpm_distribution_particle_number_max(d),
	    rvpm_distribution_particle_number(d)) ;
  }

  np = rvpm_distribution_particle_number(d) ;

  p = (RVPM_REAL *)rvpm_distribution_particle(d,np) ;

  p[RVPM_DISTRIBUTION_PARTICLE_X+0] = x[0] ;
  p[RVPM_DISTRIBUTION_PARTICLE_X+1] = x[1] ;
  p[RVPM_DISTRIBUTION_PARTICLE_X+2] = x[2] ; 
  p[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+0] = w[0]*G ;
  p[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+1] = w[1]*G ;
  p[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+2] = w[2]*G ; 
  p[RVPM_DISTRIBUTION_PARTICLE_RADIUS] = s ;

  rvpm_distribution_particle_number(d) ++ ;
  
  return 0 ;
}
			
gint RVPM_FUNCTION_NAME(rvpm_distribution_limits_crop)(rvpm_distribution_t *d,
						       RVPM_REAL gcrop,
						       RVPM_REAL *limits)

{
  RVPM_REAL bounds[6], *x, *w, absw ;  
  gint i ;

  bounds[0] = bounds[2] = bounds[4] =  G_MAXFLOAT ;
  bounds[1] = bounds[3] = bounds[5] = -G_MAXFLOAT ;
  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    w = (RVPM_REAL *)rvpm_distribution_vorticity(d, i) ;
    absw = rvpm_vector_length(w) ;
    if ( absw > gcrop ) {
      x = (RVPM_REAL *)rvpm_distribution_particle(d, i) ;
      bounds[0] = MIN(bounds[0], x[0]) ;
      bounds[1] = MAX(bounds[1], x[0]) ;
      bounds[2] = MIN(bounds[2], x[1]) ;
      bounds[3] = MAX(bounds[3], x[1]) ;
      bounds[4] = MIN(bounds[4], x[2]) ;
      bounds[5] = MAX(bounds[5], x[2]) ;
    }
  }

  limits[0] = bounds[0] ; limits[1] = bounds[1] ; 
  limits[2] = bounds[2] ; limits[3] = bounds[3] ; 
  limits[4] = bounds[4] ; limits[5] = bounds[5] ; 
  
  return 0 ;
}
