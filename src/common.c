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

static char *kernel_list[] = {
  "undefined",
  "moore_rosenhead",
  "winckelmans_leonard",
  "gaussian",
  NULL} ;

static gint distribution_write(FILE *f, rvpm_distribution_t *d)

{
  gint i, j ;
  gdouble smax, limits[6], *x ;

  limits[0] = limits[2] = limits[4] =  G_MAXDOUBLE ;
  limits[1] = limits[3] = limits[5] = -G_MAXDOUBLE ;
  smax = 0.0 ;

  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    smax = MAX(smax, (*(gdouble *)(rvpm_distribution_particle_radius(d,i)))) ;
    /* smax = MAX(smax, rvpm_distribution_particle_radius(d,i)) ; */
    x = (gdouble *)rvpm_distribution_particle(d,i) ;
    limits[0] = MIN(x[0],limits[0]) ;
    limits[1] = MAX(x[0],limits[1]) ;
    limits[2] = MIN(x[1],limits[2]) ;
    limits[3] = MAX(x[1],limits[3]) ;
    limits[4] = MIN(x[2],limits[4]) ;
    limits[5] = MAX(x[2],limits[5]) ;
  }

  fprintf(f, "%d %e %e %e %e %e %e %e\n",
	  rvpm_distribution_particle_number(d),
	  limits[0], limits[1], limits[2], limits[3], limits[4], limits[5],
	  smax) ;
  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    x = (gdouble *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      fprintf(f, "%1.16e ", x[j]) ;
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

static gint distribution_write_f(FILE *f, rvpm_distribution_t *d)

{
  gint i, j ;
  gfloat smax, limits[6], *x ;

  limits[0] = limits[2] = limits[4] =  G_MAXDOUBLE ;
  limits[1] = limits[3] = limits[5] = -G_MAXDOUBLE ;
  smax = 0.0 ;

  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    smax = MAX(smax, (*(gfloat *)(rvpm_distribution_particle_radius(d,i)))) ;
    x = (gfloat *)rvpm_distribution_particle(d,i) ;
    limits[0] = MIN(x[0],limits[0]) ;
    limits[1] = MAX(x[0],limits[1]) ;
    limits[2] = MIN(x[1],limits[2]) ;
    limits[3] = MAX(x[1],limits[3]) ;
    limits[4] = MIN(x[2],limits[4]) ;
    limits[5] = MAX(x[2],limits[5]) ;
  }

  fprintf(f, "%d %e %e %e %e %e %e %e\n",
	  rvpm_distribution_particle_number(d),
	  limits[0], limits[1], limits[2], limits[3], limits[4], limits[5],
	  smax) ;
  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    x = (gfloat *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      fprintf(f, "%1.16e ", x[j]) ;
    }
    fprintf(f, "\n") ;
  }
  
  return 0 ;
}

gint rvpm_distribution_write(FILE *f, rvpm_distribution_t *d)

{
  if ( d->size == sizeof(gfloat) ) {
    distribution_write_f(f, d) ;
    return 0 ;
  }

  if ( d->size == sizeof(gdouble) ) {
    distribution_write(f, d) ;
    return 0 ;
  }

  g_assert_not_reached() ;
  
  return 0 ;
}

gdouble rvpm_stream_func_vorticity_gaussian(gdouble r, gdouble z, gpointer data)

{
  gdouble w, z0, r0, sigma, G, s2 ;
  gdouble *p = data ;

  r0 = p[0] ; z0 = p[1] ; sigma = p[2] ; G = p[3] ;

  s2 = (r - r0)*(r - r0) + (z - z0)*(z - z0) ;

  w = G*exp(-s2/sigma/sigma)/(M_PI*sigma*sigma) ;
  
  return w ;
}

rvpm_kernel_t rvpm_kernel_parse(char *str)

{
  gint i ;

  for ( i = 1 ; kernel_list[i] != NULL ; i ++ ) {
    if ( strcmp(str, kernel_list[i]) == 0 ) return i ;
  }
  
  return RVPM_KERNEL_UNDEFINED ;
}

gchar *rvpm_kernel_name(rvpm_kernel_t kernel)

{
  if ( kernel > RVPM_KERNEL_GAUSSIAN ) return NULL ;
  
  return kernel_list[kernel] ;
}

gint rvpm_kernels_list(FILE *f)

{
  gint i ;

  for ( i = 0 ; kernel_list[i] != NULL ; i ++ ) {
    fprintf(f, "%s\n", kernel_list[i]) ;
  }
  
  return 0 ;
}

rvpm_kernel_func_t rvpm_kernel_func(rvpm_kernel_t kernel)

{
  switch ( kernel ) {
  default: return NULL ; break ;
  case RVPM_KERNEL_MOORE_ROSENHEAD: return rvpm_kernel_MR ; break ;
  case RVPM_KERNEL_WINCKELMANS_LEONARD: return rvpm_kernel_WL ; break ;
  case RVPM_KERNEL_GAUSSIAN: return rvpm_kernel_GS ; break ;
  }
  
  return NULL ;
}

rvpm_kernel_func_f_t rvpm_kernel_func_f(rvpm_kernel_t kernel)

{
  switch ( kernel ) {
  default: return NULL ; break ;
  case RVPM_KERNEL_MOORE_ROSENHEAD: return rvpm_kernel_MR_f ; break ;
  case RVPM_KERNEL_WINCKELMANS_LEONARD: return rvpm_kernel_WL_f ; break ;
  case RVPM_KERNEL_GAUSSIAN: return rvpm_kernel_GS_f ; break ;
  }
  
  return NULL ;
}
