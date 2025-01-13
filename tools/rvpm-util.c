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

static void print_help_text(FILE *f, gdouble gcrop)

{
  fprintf(f,
	  "%s: postprocessing and other utilities\n\n"
	  "Usage: %s [options]\n\n", progname, progname) ;
  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -d (vorticity distribution filename)\n"
	  "  -g # vorticity cutoff for grid cropping (%lg)\n"
	  "  -M # number of radial Gauss quadrature points in streamfunction\n"
	  "       evaluation\n"
	  "  -N # number of angular Gauss quadrature points in streamfunction\n"
	  "       evaluation\n"
	  "  -v (file of field points for evaluation of velocity and "
	  "vorticity)\n",
	  gcrop) ;
  
  return ;
}

gint main(gint argc, char **argv)

{
  char ch, *vfile, *dfile ;
  gqr_rule_t *gs, *gt ;
  rvpm_distribution_t *d ;
  rvpm_solver_t solver ;
  gint ngs, ngt, i, j ;
  gdouble data[16], r, z, th, ur, uz, x[3], u[3], w[3], r0, z0, sigma, G ;
  gdouble rmax, zmin, zmax, limits[6], s, t, *c, D, L[4], xc[12], gcrop ;
  gboolean calc_velocity ;
  FILE *output, *input ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;
  output = stdout ;
  
  vfile = NULL ; dfile = NULL ;
  calc_velocity = FALSE ;
  gs = gt = NULL ;
  ngs = ngt = 32 ;
  d = NULL ;
  gcrop = 1e-6 ;
  
  r0 = 0.9 ; z0 = 0.0 ; sigma = 0.1 ; G = 1.0 ;

  rvpm_solver_kernel(&solver)            = RVPM_KERNEL_GAUSSIAN ;
  /* rvpm_solver_time_step(&solver)         = RVPM_TIME_STEP_EULER ; */
  rvpm_solver_method(&solver)            = RVPM_METHOD_CLASSICAL ;
  rvpm_solver_thread_number(&solver)     = 1 ;
  rvpm_solver_regularisation(&solver)    = 1e-6 ;
  /* rvpm_solver_viscosity(&solver)         = 0.0 ; */
  /* rvpm_solver_model_parameter_f(&solver) = 0.0 ; */
  /* rvpm_solver_model_parameter_g(&solver) = 1/5.0 ; */
  
  while ( (ch = getopt(argc, argv, "hd:g:M:N:v:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help_text(stderr, gcrop) ;
      return 0 ;
      break ;
    case 'd': dfile = g_strdup(optarg) ; break ;
    case 'g': gcrop = atof(optarg) ; break ;
    case 'M': ngs = atoi(optarg) ; break ;
    case 'N': ngt = atoi(optarg) ; break ;
    case 'v': vfile = g_strdup(optarg) ; break ;
    }
  }
  
  /* rmax = r0 + 10.0*sigma ;  */
  /* zmin = -3.0*r0 ; zmax = 3.0*r0 ;  */
  rmax = 5.0*r0 ;
  zmin = -5.0*r0 ; zmax = 5.0*r0 ; 
  
  if ( dfile != NULL ) {
    input = file_open(dfile, "r", "-", stdin) ;

    d = rvpm_distribution_read_alloc(input, 0) ;

    file_close(input) ;
  }

  if ( vfile != NULL && dfile == NULL ) {
    gs = gqr_rule_alloc(ngs) ;
    gqr_rule_select(gs, GQR_GAUSS_LEGENDRE, ngs, NULL) ;
    gt = gqr_rule_alloc(ngt) ;
    gqr_rule_select(gt, GQR_GAUSS_LEGENDRE, ngt, NULL) ;

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
      
      if ( calc_velocity ) {
	rvpm_vorticity_velocity_gradient(d,
					 rvpm_solver_regularisation
					 (&solver),
					 rvpm_solver_kernel(&solver),
					 x, u, NULL) ;
      }
      fprintf(stdout,
	      "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e ",
	      x[0], x[1], x[2], w[0], w[1], w[2]) ;
      if ( calc_velocity ) {
	fprintf(stdout, "%1.16e %1.16e %1.16e", u[0], u[1], u[2]) ;
      }
      fprintf(stdout, "\n") ;
      /* rvpm_vorticity_velocity_gradient(d, */
      /* 				       rvpm_solver_regularisation */
      /* 				       (&solver), */
      /* 				       rvpm_solver_kernel(&solver), */
      /* 				       x, u, NULL) ; */
      /* fprintf(stdout, */
      /* 	      "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e " */
      /* 	      "%1.16e %1.16e %1.16e\n", */
      /* 	      x[0], x[1], x[2], w[0], w[1], w[2], u[0], u[1], u[2]) ; */
    }
    
    file_close(input) ;

    return 0 ;
  }

  /* ngs = 64 ; ngt = 64 ; */

  c = rvpm_distribution_origin(d) ;
  D = rvpm_distribution_width(d) ;
  limits[0] = c[0] ; limits[1] = c[0] + D ;
  limits[2] = c[1] ; limits[3] = c[1] + D ;
  limits[4] = c[2] ; limits[5] = c[2] + D ;

  rvpm_distribution_limits_crop(d, gcrop, limits) ;
  
  xc[3*0+0] = limits[0] ; xc[3*0+1] = 0.0 ; xc[3*0+2] = limits[4] ;
  xc[3*1+0] = limits[1] ; xc[3*1+1] = 0.0 ; xc[3*1+2] = limits[4] ;
  xc[3*2+0] = limits[1] ; xc[3*2+1] = 0.0 ; xc[3*2+2] = limits[5] ;
  xc[3*3+0] = limits[0] ; xc[3*3+1] = 0.0 ; xc[3*3+2] = limits[5] ;
  
  for ( i = 0 ; i < ngs ; i ++ ) {
    s = ((gdouble)i)/(ngs - 1) ;
    for ( j = 0 ; j < ngt ; j ++ ) {
      t = ((gdouble)j)/(ngt - 1) ;
      L[0] = (1.0 - s)*(1.0 - t) ;
      L[1] =        s *(1.0 - t) ;
      L[2] =        s *       t ;
      L[3] = (1.0 - s)*       t ;
      x[0] = L[0]*xc[3*0+0] + L[1]*xc[3*1+0] + L[2]*xc[3*2+0] + L[3]*xc[3*3+0] ;
      x[1] = L[0]*xc[3*0+1] + L[1]*xc[3*1+1] + L[2]*xc[3*2+1] + L[3]*xc[3*3+1] ;
      x[2] = L[0]*xc[3*0+2] + L[1]*xc[3*1+2] + L[2]*xc[3*2+2] + L[3]*xc[3*3+2] ;
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

      if ( calc_velocity ) {
	rvpm_vorticity_velocity_gradient(d,
					 rvpm_solver_regularisation
					 (&solver),
					 rvpm_solver_kernel(&solver),
					 x, u, NULL) ;
      }
      fprintf(stdout,
	      "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e ",
	      x[0], x[1], x[2], w[0], w[1], w[2]) ;
      if ( calc_velocity ) {
	fprintf(stdout, "%1.16e %1.16e %1.16e", u[0], u[1], u[2]) ;
      }
      fprintf(stdout, "\n") ;
    }
    
  }
  
  return 0 ;
}
