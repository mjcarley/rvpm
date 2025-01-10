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

/* static FILE *file_open(char *file, char *mode, */
/* 		       char *file_default, FILE *fdefault) */

/* { */
/*   FILE *f ; */

/*   if ( file_default != NULL ) { */
/*     if ( strcmp(file, file_default) == 0 ) return fdefault ; */
/*   } */
  
/*   f = fopen(file, mode) ; */

/*   if ( f == NULL ) { */
/*     fprintf(stderr, "%s: cannot open file \"%s\"\n", progname, file) ; */
/*     exit(1) ; */
/*   } */
  
/*   return f ; */
/* } */

/* static void file_close(FILE *f) */

/* { */
/*   if ( f == stdin || f == stdout || f == stderr ) return ; */

/*   fclose(f) ; */
  
/*   return ; */
/* } */

static void print_help_text(FILE *f, gint depth, gfloat dt,
			    rvpm_kernel_t kernel, gint order_max,
			    gint ns, gfloat reg, gint nthreads,
			    gfloat nu)

{
  fprintf(f,
	  "%s: advance a vorticity distribution in time\n\n"
	  "Usage: %s <options> < (input file) > (output file)\n\n",
	  progname, progname) ;
  fprintf(f,
	  "Options:\n\n"
	  "  -h print this message and exit\n"
	  "  -D # depth of FMM tree for velocity evaluation (%d)\n"
	  "  -d # time step (%g)\n"
	  "  -G select Gaussian interpolation kernel\n"
	  "  -L # maximum expansion order in FMM tree (%d)\n"
	  "  -M select Moore-Rosenhead kernel (deprecated)\n"
	  "  -n # number of time steps (%d)\n"
	  "  -r # regularisation parameter (%g)\n"
	  "  -T # number of threads in velocity evaluation (%d)\n"
	  "  -v # viscosity (%g)\n"
	  "  -W select Winckelmans-Leonard kernel (deprecated)\n",
	  depth, dt, order_max, ns, reg, nthreads, nu) ;

  return ;
}

gint main(gint argc, char **argv)

{
  char ch ;
  rvpm_distribution_t *d ;
  rvpm_tree_t *tree ;
  rvpm_solver_t solver ;
  gfloat dt, *u, t0, *work ;
  gint ns, i, ustr, order_max, depth, wsize ;
  /* FILE *input ; */
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;

  dt = 0.01 ; ns = 10 ; order_max = 10 ; depth = 4 ;

  rvpm_solver_kernel(&solver)            = RVPM_KERNEL_GAUSSIAN ;
  rvpm_solver_time_step(&solver)         = RVPM_TIME_STEP_EULER ;
  rvpm_solver_method(&solver)            = RVPM_METHOD_CLASSICAL ;
  rvpm_solver_thread_number(&solver)     = 1 ;
  rvpm_solver_regularisation(&solver)    = 1e-6 ;
  rvpm_solver_viscosity(&solver)         = 0.0 ;
  rvpm_solver_model_parameter_f(&solver) = 0.0 ;
  rvpm_solver_model_parameter_g(&solver) = 1/5.0 ;
  
  rvpm_solver_time_step(&solver)     = RVPM_TIME_STEP_WILLIAMSON_19 ;

  while ( (ch = getopt(argc, argv, "hD:d:GL:Mn:r:T:v:W")) != EOF ) {
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
    case 'G': solver.kernel = RVPM_KERNEL_GAUSSIAN ; break ;
    case 'L': order_max = atoi(optarg) ; break ;
    case 'M': solver.kernel = RVPM_KERNEL_MOORE_ROSENHEAD ; break ;
    case 'n': ns = atoi(optarg) ; break ;
    case 'r': solver.reg = atof(optarg) ; break ;
    case 'T': solver.nthreads = atoi(optarg) ; break ;
    case 'v': rvpm_solver_viscosity(&solver) = atof(optarg) ; break ;
    case 'W': solver.kernel = RVPM_KERNEL_WINCKELMANS_LEONARD ; break ;
    }
  }
  
  d = rvpm_distribution_read_alloc_f(stdin, 0) ;

  fprintf(stderr, "%s: particle distribution read\n", progname) ;
  fprintf(stderr, "%s: %d particles\n",
	  progname, rvpm_distribution_particle_number(d)) ;

  fprintf(stderr, "setting up FMM tree\n") ;
  wsize = wbfmm_element_number_rotation(2*order_max) ;
  wsize = MAX(wsize, (order_max+1)*(order_max+1)*3*16) ;
  work = (gfloat *)g_malloc0(wsize*sizeof(gfloat)) ;
  tree = rvpm_tree_new_f(d, depth, order_max, work) ;
  
  /*velocity buffer for time stepping*/
  ustr = 12 ;
  u = (gfloat *)g_malloc0(ustr*rvpm_distribution_particle_number(d)*
			   sizeof(gfloat)) ;

  fprintf(stderr, "%s: starting solve for %d steps\n", progname, ns) ;
  t0 = g_timer_elapsed(timer, NULL) ;   
  for ( i = 0 ; i < ns ; i ++ ) {
    fprintf(stderr, "%s: time step %d [%lg]\n",
	    progname, i, g_timer_elapsed(timer, NULL)-t0) ;
    rvpm_solver_solve_f(tree, &solver, i*dt, dt, u, ustr, work) ;
  }

  fprintf(stderr, "%s: solve completed [%lg]\n",
	  progname, g_timer_elapsed(timer, NULL)-t0) ;

  rvpm_distribution_write(stdout, d) ;
  
  return 0 ;
}
