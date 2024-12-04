/* This file is part of RVPM, a library for the Reformulated Vortex
 * Particle Method of Alvarez and Ning (see README.md and
 * documentation for references)
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
#include <unistd.h>
#include <math.h>

#include <glib.h>

#include <grbf.h>

#include <gqr.h>

#include "rvpm.h"
#include "rvpm-private.h"

GTimer *timer ;

char *progname ;

/* gint vortex_ring(gdouble *x, gdouble *p, gint np, gdouble *w, gboolean) ; */

/* gpointer source_functions[] = { */
/*   "vortex_ring", vortex_ring, */
/*   NULL, NULL */
/* } ; */

static FILE *file_open(char *file, char *mode,
		       char *file_default, FILE *fdefault)

{
  FILE *f ;

  if ( file_default != NULL ) {
    if ( strcmp(file, file_default) == 0 ) return fdefault ;
  }
  
  f = fopen(file, mode) ;

  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open file \"%s\"\n", progname, file) ;
    exit(1) ;
  }
  
  return f ;
}

static void file_close(FILE *f)

{
  if ( f == stdin || f == stdout || f == stderr ) return ;

  fclose(f) ;
  
  return ;
}

/* gint vortex_ring(gdouble *x, gdouble *p, gint np, gdouble *w, gboolean limits) */

/* { */
/*   gdouble r0, w0, z0, G, r, r2, th, z, s, ee ; */

/*   if ( np < 2 ) { */
/*     fprintf(stderr, "%s: vortex ring requires at least two parameters\n", */
/* 	    progname) ; */
/*     return 1 ; */
/*   } */

/*   r0 = p[0] ; s = p[1] ; */
/*   if ( np < 3 ) z0 = 0.0 ; else z0 = p[2] ; */
/*   if ( np < 4 ) G = 1.0 ; else G = p[3] ; */

/*   w0 = G/(M_PI*s*s) ; */
  
/*   if ( limits ) { */
/*     /\*treat the first entry of w as a tolerance to be used in finding */
/*       the limits of the grid box*\/ */
/*     ee = w[0] ; */
/*     r = sqrt(-s*s*log(ee/w0)) ; */
/*     x[0] = -r0 - r ; x[1] =  r0 + r ; */
/*     x[2] = -r0 - r ; x[3] =  r0 + r ; */
/*     x[4] =  z0 - r ; x[5] =  z0 + r ; */

/*     return 0 ; */
/*   } */

/*   r = sqrt(x[0]*x[0] + x[1]*x[1]) ; */
/*   th = atan2(x[1], x[0]) ; */
/*   z = x[2] ; */

/*   r2 = (r - r0)*(r - r0) + (z-z0)*(z-z0) ; */
/*   w0 *= exp(-r2/s/s) ; */

/*   w[0] += -w0*sin(th) ; */
/*   w[1] +=  w0*cos(th) ; */
/*   w[2] +=  0.0 ; */
  
/*   return 0 ; */
/* } */

/* static gpointer parse_func(char *str) */

/* { */
/*   gint i ; */

/*   for ( i = 0 ; source_functions[i] != NULL ; i += 2 ) { */
/*     if ( strcmp((char *)source_functions[i], str) == 0 ) { */
/*       return source_functions[i+1] ; */
/*     } */
/*   } */
  
/*   return NULL ; */
/* } */

/* static gint grid_limits(gdouble *wt, gdouble xmin, gdouble ymin, gdouble zmin, */
/* 			gint nx, gint ny, gint nz, gdouble h, */
/* 			gdouble *limits, gdouble tol) */

/* { */
/*   gint i, j, k, np ; */
/*   gdouble x[3], *w, norm ; */
  
/*   np = 0 ; */

/*   limits[0] = limits[2] = limits[4] =  G_MAXDOUBLE ; */
/*   limits[1] = limits[3] = limits[5] = -G_MAXDOUBLE ; */
  
/*   for ( i = 0 ; i < nx ; i ++ ) { */
/*     x[0] = xmin + h*i ; */
/*     for ( j = 0 ; j < ny ; j ++ ) { */
/*       x[1] = ymin + h*j ; */
/*       for ( k = 0 ; k < nz ; k ++ ) { */
/* 	x[2] = zmin + h*k ; */
/* 	w = &(wt[3*(i*ny*nz + j*nz + k)]) ; */
/* 	norm = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ; */
/* 	if ( norm > tol*tol ) { */
/* 	  limits[0] = MIN(x[0],limits[0]) ; */
/* 	  limits[1] = MAX(x[0],limits[1]) ; */
/* 	  limits[2] = MIN(x[1],limits[2]) ; */
/* 	  limits[3] = MAX(x[1],limits[3]) ; */
/* 	  limits[4] = MIN(x[2],limits[4]) ; */
/* 	  limits[5] = MAX(x[2],limits[5]) ; */
/* 	  np ++ ; */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   return np ; */
/* } */

/* static void write_active_nodes(FILE *f, gdouble *wt, */
/* 			       gdouble xmin, gdouble ymin, gdouble zmin, */
/* 			       gint nx, gint ny, gint nz, gdouble h, gdouble s, */
/* 			       gdouble tol) */

/* { */
/*   gint i, j, k ; */
/*   gdouble x[3], *w, norm, sc ; */

/*   sc = M_PI*sqrt(M_PI)*s*s*s ; */

/*   for ( i = 0 ; i < nx ; i ++ ) { */
/*     x[0] = xmin + h*i ; */
/*     for ( j = 0 ; j < ny ; j ++ ) { */
/*       x[1] = ymin + h*j ; */
/*       for ( k = 0 ; k < nz ; k ++ ) { */
/* 	x[2] = zmin + h*k ; */
/* 	w = &(wt[3*(i*ny*nz + j*nz + k)]) ; */
/* 	norm = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ; */
/* 	if ( norm > tol*tol ) { */
/* 	  fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", */
/* 		  x[0], x[1], x[2], w[0]*sc, w[1]*sc, w[2]*sc, s) ; */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   return ; */
/* } */

gint main(gint argc, char **argv)

{
  char ch, *vfile, *dfile ;
  gqr_rule_t *gs, *gt ;
  rvpm_distribution_t *d ;
  rvpm_solver_t solver ;
  gint ngs, ngt ;
  gdouble data[16], r, z, th, ur, uz, x[3], u[3], w[3], r0, z0, sigma, G ;
  gdouble rmax, zmin, zmax ;
  FILE *output, *input ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;
  output = stdout ;
  
  vfile = NULL ; dfile = NULL ;
  gs = gt = NULL ;
  ngs = ngt = 32 ;

  r0 = 0.9 ; z0 = 0.0 ; sigma = 0.1 ; G = 1.0 ;

  rvpm_solver_kernel(&solver)            = RVPM_KERNEL_GAUSSIAN ;
  /* rvpm_solver_time_step(&solver)         = RVPM_TIME_STEP_EULER ; */
  rvpm_solver_method(&solver)            = RVPM_METHOD_CLASSICAL ;
  rvpm_solver_thread_number(&solver)     = 1 ;
  rvpm_solver_regularisation(&solver)    = 1e-6 ;
  /* rvpm_solver_viscosity(&solver)         = 0.0 ; */
  /* rvpm_solver_model_parameter_f(&solver) = 0.0 ; */
  /* rvpm_solver_model_parameter_g(&solver) = 1/5.0 ; */
  
  while ( (ch = getopt(argc, argv, "d:M:N:v:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'd': dfile = g_strdup(optarg) ; break ;
    case 'M': ngs = atoi(optarg) ; break ;
    case 'N': ngt = atoi(optarg) ; break ;
    case 'v': vfile = g_strdup(optarg) ; break ;
    }
  }
  
  gs = gqr_rule_alloc(ngs) ;
  gqr_rule_select(gs, GQR_GAUSS_LEGENDRE, ngs, NULL) ;
  gt = gqr_rule_alloc(ngt) ;
  gqr_rule_select(gt, GQR_GAUSS_LEGENDRE, ngt, NULL) ;

  /* rmax = r0 + 10.0*sigma ;  */
  /* zmin = -3.0*r0 ; zmax = 3.0*r0 ;  */
  rmax = 5.0*r0 ;
  zmin = -5.0*r0 ; zmax = 5.0*r0 ; 
  
  if ( dfile != NULL ) {
    input = file_open(dfile, "r", "-", stdin) ;

    d = rvpm_distribution_read_alloc(input) ;

    file_close(input) ;
  }

  if ( vfile != NULL && dfile == NULL ) {
    data[0] = r0 ; data[1] = z0 ; data[2] = sigma ; data[3] = G ;

    input = file_open(vfile, "r", "-", stdin) ;

    while ( fscanf(input, "%lg %lg %lg", &(x[0]), &(x[1]), &(x[2]))
	    != EOF ) {
      r = sqrt(x[0]*x[0] + x[1]*x[1]) ;
      th = atan2(x[1], x[0]) ;
      z = x[2] ;
      ur = uz = 0.0 ;
      rvpm_stream_func_velocity(rvpm_stream_func_vorticity_gaussian, data,
				r, z, rmax, zmin, zmax, &ur, &uz,
				gs, gt) ;
      u[0] = ur*cos(th) ; u[1] = ur*sin(th) ; u[2] = uz ;
      fprintf(output, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
	      x[0], x[1], x[2], u[0], u[1], u[2]) ;
    }

    file_close(input) ;

    return 0 ;
  }

  if ( vfile != NULL && dfile != NULL ) {
    input = file_open(vfile, "r", "-", stdin) ;

    while ( fscanf(input, "%lg %lg %lg", &(x[0]), &(x[1]), &(x[2])) != EOF ) {
      u[0] = u[1] = u[2] = 0.0 ;
      w[0] = w[1] = w[2] = 0.0 ;
      grbf_gaussian_eval_3d((gdouble *)rvpm_distribution_particle(d,0),
			    RVPM_DISTRIBUTION_PARTICLE_SIZE,
			    rvpm_distribution_particle_number(d),
			    (gdouble *)
			    rvpm_distribution_particle_radius(d,0),
			    RVPM_DISTRIBUTION_PARTICLE_SIZE,
			    (gdouble *)rvpm_distribution_vorticity(d,0),
			    RVPM_DISTRIBUTION_PARTICLE_SIZE,
			    3, x, w) ;
      
      rvpm_vorticity_velocity_gradient(d,
				       rvpm_solver_regularisation
				       (&solver),
				       rvpm_solver_kernel(&solver),
				       x, u, NULL) ;
      fprintf(stdout,
	      "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e "
	      "%1.16e %1.16e %1.16e\n",
	      x[0], x[1], x[2], w[0], w[1], w[2], u[0], u[1], u[2]) ;
    }
    
    file_close(input) ;

    return 0 ;
  }
  
  return 0 ;
}
