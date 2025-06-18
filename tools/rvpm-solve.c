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

#include "rvpm.h"
#include "rvpm-private.h"

#include <grbf.h>

GTimer *timer ;

char *progname ;

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

static void print_help_text(FILE *f, gint depth, gdouble dt,
			    rvpm_kernel_t kernel, gint order_max,
			    gint ns, gdouble reg, gint nthreads,
			    gdouble nu)

{
  fprintf(f,
	  "%s: advance a vorticity distribution in time\n\n"
	  "Usage: %s <options> < (input file) > (output file)\n\n",
	  progname, progname) ;
  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -D # depth of FMM tree for velocity evaluation (%d)\n"
	  "  -d # time step (%lg)\n"
	  "  -G select Gaussian interpolation kernel\n"
	  "  -L # maximum expansion order in FMM tree (%d)\n"
	  "  -M select Moore-Rosenhead kernel (deprecated)\n"
	  "  -n # number of time steps (%d)\n"
	  "  -r # regularisation parameter (%lg)\n"
	  "  -T # number of threads in velocity evaluation (%d)\n"
	  "  -v # viscosity (%lg)\n"
	  "  -W select Winckelmans-Leonard kernel (deprecated)\n",
	  depth, dt, order_max, ns, reg, nthreads, nu) ;

  return ;
}

gint main(gint argc, char **argv)

{
  char ch, *vfile ;
  rvpm_distribution_t *d ;
  rvpm_tree_t *tree ;
  rvpm_solver_t solver ;
  gdouble dt, *u, t0, *work ;
  gint ns, i, ustr, order_max, depth, wsize ;
  gboolean velocity_direct, velocity_gradient ;
  FILE *f ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  dt = 0.01 ; ns = 10 ; order_max = 10 ; depth = 4 ;
  vfile = NULL ; velocity_direct = FALSE ; velocity_gradient = FALSE ;
  
  rvpm_solver_kernel(&solver)            = RVPM_KERNEL_GAUSSIAN ;
  rvpm_solver_time_step(&solver)         = RVPM_TIME_STEP_EULER ;
  rvpm_solver_method(&solver)            = RVPM_METHOD_CLASSICAL ;
  rvpm_solver_thread_number(&solver)     = 1 ;
  rvpm_solver_regularisation(&solver)    = 1e-6 ;
  rvpm_solver_viscosity(&solver)         = 0.0 ;
  rvpm_solver_model_parameter_f(&solver) = 0.0 ;
  rvpm_solver_model_parameter_g(&solver) = 1/5.0 ;
  
  rvpm_solver_time_step(&solver)     = RVPM_TIME_STEP_WILLIAMSON_19 ;

  while ( (ch = getopt(argc, argv, "hD:d:gGL:Mn:r:sT:V:v:W")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'h':
      print_help_text(stderr, depth, dt,
		      rvpm_solver_kernel(&solver),
		      order_max, ns,
		      rvpm_solver_regularisation(&solver),
		      rvpm_solver_thread_number(&solver),
		      rvpm_solver_viscosity(&solver)) ;
      return 0 ;
      break ;
    case 'D': depth = atoi(optarg) ; break ;
    case 'd': dt = atof(optarg) ; break ;
    case 'g': velocity_gradient = TRUE ; break ;
    case 'L': order_max = atoi(optarg) ; break ;
    case 'n': ns = atoi(optarg) ; break ;
    case 's': velocity_direct = TRUE ; break ;
    case 'V': vfile = g_strdup(optarg) ; break ;
    case 'r': rvpm_solver_regularisation(&solver) = atof(optarg) ; break ;
    case 'T': rvpm_solver_thread_number(&solver)  = atoi(optarg) ; break ;
    case 'v': rvpm_solver_viscosity(&solver)      = atof(optarg) ; break ;
      /*kernel options*/
    case 'G': rvpm_solver_kernel(&solver) =
	RVPM_KERNEL_GAUSSIAN ; break ;
    case 'M': rvpm_solver_kernel(&solver) =
	RVPM_KERNEL_MOORE_ROSENHEAD ; break ;
    case 'W': rvpm_solver_kernel(&solver) =
	RVPM_KERNEL_WINCKELMANS_LEONARD ; break ;
    }
  }
  
  d = rvpm_distribution_read_alloc(stdin, 0) ;

  fprintf(stderr, "%s: particle distribution read\n", progname) ;
  fprintf(stderr, "%s: %d particles\n",
	  progname, rvpm_distribution_particle_number(d)) ;

  fprintf(stderr, "setting up FMM tree\n") ;
  wsize = wbfmm_element_number_rotation(2*order_max) ;
  wsize = MAX(wsize, (order_max+1)*(order_max+1)*3*16) ;
  work = (gdouble *)g_malloc0(wsize*sizeof(gdouble)) ;
  tree = rvpm_tree_new(d, depth, order_max, work) ;
  
  if ( vfile != NULL ) {
    gdouble x[3], u[12] ;
    
    f = file_open(vfile, "r", "-", stdin) ;

    if ( f == NULL ) {
      fprintf(stderr, "%s: cannot open file \"%s\"", progname, vfile) ;
      return 1 ;
    }

    rvpm_tree_update(tree,
			   rvpm_solver_thread_number(&solver), TRUE, work) ;

    
    while ( fscanf(f, "%lg %lg %lg",
		   &(x[0]), &(x[1]), &(x[2])) != EOF ) {
      memset(u, 0, 12*sizeof(gdouble)) ;
      if ( velocity_direct ) {
	  rvpm_vorticity_velocity_gradient(tree->d,
						 rvpm_solver_regularisation(&solver),
						 rvpm_solver_kernel(&solver),
						 x, &(u[0]), &(u[3])) ;
      } else {
	if ( velocity_gradient ) {
	  rvpm_tree_velocity_gradient(tree,
					    rvpm_solver_regularisation(&solver),
					    rvpm_solver_kernel(&solver),
					    x, &(u[0]), &(u[3]), work) ;
	} else {
	  rvpm_tree_velocity(tree,
				   rvpm_solver_regularisation(&solver),
				   rvpm_solver_kernel(&solver),
				   x, &(u[0]), work) ;
	}
      }
      fprintf(stdout, "%1.16e %1.16e %1.16e ", x[0], x[1], x[2]) ;
      fprintf(stdout, "%1.16e %1.16e %1.16e ", u[0], u[1], u[2]) ;
      if ( velocity_gradient ) {
	fprintf(stdout, "%1.16e %1.16e %1.16e ", u[3], u[4], u[5]) ;
	fprintf(stdout, "%1.16e %1.16e %1.16e ", u[6], u[7], u[8]) ;
	fprintf(stdout, "%1.16e %1.16e %1.16e ", u[9], u[10], u[11]) ;
      }
      fprintf(stdout, "\n") ;
    }
    
    file_close(f) ;
    
    return 0 ;
  }
  
  /*velocity buffer for time stepping*/
  ustr = 12 ;
  u = (gdouble *)g_malloc0(ustr*rvpm_distribution_particle_number(d)*
			   sizeof(gdouble)) ;

  fprintf(stderr, "%s: starting solve for %d steps\n", progname, ns) ;
  t0 = g_timer_elapsed(timer, NULL) ;   
  for ( i = 0 ; i < ns ; i ++ ) {
    fprintf(stderr, "%s: time step %d [%lg]\n",
	    progname, i, g_timer_elapsed(timer, NULL)-t0) ;
    rvpm_solver_solve(tree, &solver, i*dt, dt, u, ustr, work) ;
  }

  fprintf(stderr, "%s: solve completed [%lg]\n",
	  progname, g_timer_elapsed(timer, NULL)-t0) ;

  rvpm_distribution_write(stdout, d) ;
  
  return 0 ;
}
