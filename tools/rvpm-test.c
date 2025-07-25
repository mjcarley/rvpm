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
  char *tests[] = {
    "gauss_fit",
    "permutation",
    "tree",
    "solver",
    "kernel_gradient",
    "stream_function",
    "stream_function_velocity",
    "time_step",
    "gauss_kernel",
    NULL} ;

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

static gint parse_test(char *str)

{
  gint i ;

  for ( i = 0 ; tests[i] != NULL ; i ++ ) {
    if ( strcmp(str, tests[i]) == 0 ) return i ;
  }
  
  return -1 ;
}

static void vortex_ring_func(gdouble r0, gdouble s, gdouble G,
			     gdouble *x, gdouble *f)

{
  gdouble r, z, th, r2, w, w0 ;

  w0 = G/(M_PI*s*s) ;
  
  r = sqrt(x[0]*x[0] + x[1]*x[1]) ;
  th = atan2(x[1], x[0]) ;
  z = x[2] ;

  r2 = (r - r0)*(r - r0) + z*z ;
  w = exp(-r2/s/s)*w0 ;

  f[0] = -w*sin(th) ;
  f[1] =  w*cos(th) ;
  f[2] = 0.0 ;
  
  return ;
}

static void distribution_random(rvpm_distribution_t *d, gint np, gdouble s)

{
  gdouble D, *x ;
  gdouble *c ;
  gint i ;
  
  g_assert(np <= rvpm_distribution_particle_number_max(d)) ;
  
  c = rvpm_distribution_origin(d) ;
  D = rvpm_distribution_width(d) ;

  for ( i = 0 ; i < np ; i ++ ) {
    x = (gdouble *)rvpm_distribution_particle(d,i) ;
    x[0] = c[0] + D*g_random_double() ;
    x[1] = c[1] + D*g_random_double() ;
    x[2] = c[2] + D*g_random_double() ;
    *((gdouble *)rvpm_distribution_particle_radius(d,i)) = s ;
    rvpm_distribution_particle_number(d) ++ ;
  }

  return ;
}

static void distribution_write(FILE *f, rvpm_distribution_t *d)

{
  gint i, j ;
  gdouble *y ;
  
  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    y = (gdouble *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      fprintf(f, "%e ", y[j]) ;
    }
    fprintf(f, "\n") ;
  }
  

  return ;
}

static gdouble kernel(gdouble r)

{
  gdouble f ;

  f = exp(-r*r)/M_PI/sqrt(M_PI) ;

  return f ;
}  

static void gauss_fit_test(gint np)

{
  gdouble *f, s, *work, *wt, x[3], omi[3], *om, *y, zeta ;
  gdouble r0, s0, sc, limits[6], h, tol, t0, al, r, gtol, dt ;
  char fname[64] ;
  grbf_workspace_t *ws ;
  rvpm_distribution_t *d ;
  gint i, j, k, n, N, wsize, nx, ny, nz ;
  FILE *fid ;
  
  s = 0.03 ; tol = 1e-3 ; gtol = 1e-2 ;
  r0 = 0.8 ; s0 = 0.1 ;
  /*fit parameters*/
  sc = 16.0*s0 ; al = 0.8 ;

  dt = 0.0125 ;
  
  fprintf(stderr, "Gauss transform coefficient fitting test\n") ;
  fprintf(stderr, "========================================\n") ;

  N = grbf_cardinal_function_length(al, tol) ;
  fprintf(stderr, "Gaussian fit parameter: alpha = %lg\n", al) ;
  fprintf(stderr, "Gaussian fit parameter: N     = %d\n", N) ;
  /*workspace for setting up FFTs*/
  wsize = (N+1)*(N+1) + 3*(N+1) ;
  
  /*spacing parameters*/
  h = al*s ;
  limits[0] = -r0 - sc ; limits[1] = r0 + sc ; 
  limits[2] = -r0 - sc ; limits[3] = r0 + sc ; 
  limits[4] =      -sc ; limits[5] =      sc ; 
  grbf_grid_adjust_3d(&(limits[0]), &(limits[1]), &nx,
			    &(limits[2]), &(limits[3]), &ny,
			    &(limits[4]), &(limits[5]), &nz, h) ;
  work = (gdouble *)g_malloc0((nx*ny*nz*3 + wsize)*sizeof(gdouble)) ;
  wt   = (gdouble *)g_malloc0(nx*ny*nz*3*sizeof(gdouble)) ;

  f = &(work[wsize]) ;
  
  fprintf(stderr, "ring parameters: s = %lg\n", s) ;
  fprintf(stderr, "ring parameters: r = %lg\n", r0) ;

  fprintf(stderr, "domain limits x: %lg %lg (%d)\n",
	  limits[0], limits[1], nx) ;
  fprintf(stderr, "domain limits y: %lg %lg (%d)\n",
	  limits[2], limits[3], ny) ;
  fprintf(stderr, "domain limits z: %lg %lg (%d)\n",
	  limits[4], limits[5], nz) ;  

  ws = grbf_workspace_alloc(3, MAX(nx,MAX(ny,nz))) ;
  grbf_workspace_init_3d(ws, f, 3, nx, 3*ny*nz, ny, 3*nz, nz,
			       wt, 3, 3*ny*nz, 3*nz, al, N, work) ;
  
  for ( i = 0 ; i < nx ; i ++ ) {
    x[0] = limits[0] + h*i ;
    for ( j = 0 ; j < ny ; j ++ ) {
      x[1] = limits[2] + h*j ;
      for ( k = 0 ; k < nz ; k ++ ) {
	x[2] = limits[4] + h*k ;
	vortex_ring_func(r0, s0, 1, x, &(f[3*(i*ny*nz + j*nz + k)])) ;
      }
    }
  }
  
  fprintf(stderr, "evaluating Gaussian weights (FFT) ") ;
  t0 = g_timer_elapsed(timer, NULL) ;
  grbf_interpolation_weights_fft_3d(f, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
					  ws, wt, 3, 3*ny*nz, 3*nz) ;
  fprintf(stderr, "[%lg]\n", g_timer_elapsed(timer, NULL)-t0) ;

  np = grbf_grid_count_points_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
				       gtol) ;

  fprintf(stderr, "computational nodes: %d\n", np) ;
  d = rvpm_distribution_alloc(np) ;

  n = 0 ;
  for ( i = 0 ; i < nx ; i ++ ) {
    for ( j = 0 ; j < ny ; j ++ ) {
      for ( k = 0 ; k < nz ; k ++ ) {
	om = &(wt[3*(i*ny*nz + j*nz + k)]) ;
	if ( rvpm_vector_scalar(om,om) > gtol*gtol ) {
	  y = (gdouble *)rvpm_distribution_particle(d,n) ;
	  y[0] = limits[0] + h*i ;
	  y[1] = limits[2] + h*j ;
	  y[2] = limits[4] + h*k ;

	  zeta = M_PI*sqrt(M_PI)*s*s*s ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+0] = om[0]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+1] = om[1]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+2] = om[2]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_RADIUS] = s ;
	  n ++ ;
	}	  
      }
    }
  }
  rvpm_distribution_particle_number(d) = n ;
  g_assert(n == np) ;

  /*check rescaled coefficients give correct vorticity*/
  for ( i = 0 ; i < 100 ; i ++ ) {
    /* xi = rvpm_distribution_particle(d,i) ; */

    x[0] = r0 ; x[2] = s0*0.5 ;
    x[1] = limits[2] + (limits[3] - limits[2])*i/100 ;
    
    omi[0] = omi[1] = omi[2] = 0.0 ;
    for ( j = 0 ; j < np ; j ++ ) {
      y = (gdouble *)rvpm_distribution_particle(d,j) ;
      /* r = sqrt(rvpm_vector_distance2(xi,y)) ; */
      r = sqrt(rvpm_vector_distance2(x,y)) ;
      om = (gdouble *)rvpm_distribution_vorticity(d,j) ;
      s = *((gdouble *)rvpm_distribution_particle_radius(d,j)) ;
      zeta = kernel(r/s)/s/s/s ;
      /* zeta = kernel(r/s) ; */
      omi[0] += om[0]*zeta ;
      omi[1] += om[1]*zeta ;
      omi[2] += om[2]*zeta ;
    }
    /* fprintf(stdout, "%e %e %e %e %e %e", */
    /* 	    xi[0], xi[1], xi[2], omi[0], omi[1], omi[2]) ; */
    /* vortex_ring_func(r0, s0, 1, xi, omi) ; */
    fprintf(stdout, "%e %e %e %1.16e %1.16e %1.16e",
	    x[0], x[1], x[2], omi[0], omi[1], omi[2]) ;
    vortex_ring_func(r0, s0, 1, x, omi) ;
    fprintf(stdout, " %1.16e %1.16e %1.16e\n", omi[0], omi[1], omi[2]) ;
  }
  
  fid = file_open("points.dat", "w", "-", stdout) ;

  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    y = (gdouble *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      fprintf(fid, "%e ", y[j]) ;
    }
    fprintf(fid, "\n") ;
  }
  
  file_close(fid) ;

  g_assert(np*3 < nx*ny*nz*3) ;
  /*time step for the crack*/
  for ( i = 0 ; i < 200 ; i ++ ) {
    fprintf(stderr, "time %d\n", i) ;
    memset(work, 0, 3*np*sizeof(gdouble)) ;
    for ( j = 0 ; j < np ; j ++ ) {
      rvpm_vorticity_velocity_gradient(d, s,
					     RVPM_KERNEL_WINCKELMANS_LEONARD,
					     (gdouble *)rvpm_distribution_particle(d,j),
					     &(work[3*j]), NULL) ;
    }
    for ( j = 0 ; j < np ; j ++ ) {
      y = (gdouble *)rvpm_distribution_particle(d,j) ;
      y[0] += work[3*j+0]*dt ;
      y[1] += work[3*j+1]*dt ;
      y[2] += work[3*j+2]*dt ;
    }

    if ( 10*(i/10) == i ) {
      sprintf(fname, "points-%d.dat", i) ;
    }

    fid = file_open(fname, "w", "-", stdout) ;

    distribution_write(fid, d) ;

    file_close(fid) ;    
  }

  fid = file_open("points-end.dat", "w", "-", stdout) ;

  distribution_write(fid, d) ;
  
  file_close(fid) ;
  
  return ;
}

static void tree_test(rvpm_kernel_t kernel, gint np)

{
  gdouble *f, s, *work, *wt, x[3], *om, *y, zeta, reg ;
  gdouble r0, s0, sc, limits[6], h, tol, t0, al, gtol, uv[12], ut[12] ;
  gdouble err[12] = {0} ;
  grbf_workspace_t *ws ;
  rvpm_distribution_t *d ;
  rvpm_tree_t *tree ;
  gint i, j, k, n, N, wsize, nx, ny, nz ;
  
  s = 0.03 ; tol = 1e-6 ; gtol = 1e-3 ; reg = 1e-3 ;
  r0 = 0.8 ; s0 = 0.1 ;
  /*fit parameters*/
  sc = 16.0*s0 ; al = 1.0 ;

  fprintf(stderr, "FMM tree velocity evaluation test\n") ;
  fprintf(stderr, "=================================\n") ;

  if ( kernel == RVPM_KERNEL_UNDEFINED ) {
    fprintf(stderr, "unrecognised kernel\n") ;

    return ;
  }

  fprintf(stderr, "kernel: %s\n", rvpm_kernel_name(kernel)) ;
  
  N = grbf_cardinal_function_length(al, tol) ;
  fprintf(stderr, "Gaussian fit parameter: alpha = %lg\n", al) ;
  fprintf(stderr, "Gaussian fit parameter: N     = %d\n", N) ;
  /*workspace for setting up FFTs*/
  wsize = (N+1)*(N+1) + 3*(N+1) ;
  
  /*spacing parameters*/
  h = al*s ;
  limits[0] = -r0 - sc ; limits[1] = r0 + sc ; 
  limits[2] = -r0 - sc ; limits[3] = r0 + sc ; 
  limits[4] =      -sc ; limits[5] =      sc ; 
  grbf_grid_adjust_3d(&(limits[0]), &(limits[1]), &nx,
			    &(limits[2]), &(limits[3]), &ny,
			    &(limits[4]), &(limits[5]), &nz, h) ;
  work = (gdouble *)g_malloc0((nx*ny*nz*3 + wsize)*sizeof(gdouble)) ;
  wt   = (gdouble *)g_malloc0(nx*ny*nz*3*sizeof(gdouble)) ;

  f = &(work[wsize]) ;
  
  fprintf(stderr, "ring parameters: s = %lg\n", s0) ;
  fprintf(stderr, "ring parameters: r = %lg\n", r0) ;

  fprintf(stderr, "domain limits x: %lg %lg (%d)\n",
	  limits[0], limits[1], nx) ;
  fprintf(stderr, "domain limits y: %lg %lg (%d)\n",
	  limits[2], limits[3], ny) ;
  fprintf(stderr, "domain limits z: %lg %lg (%d)\n",
	  limits[4], limits[5], nz) ;  

  ws = grbf_workspace_alloc(3, MAX(nx,MAX(ny,nz))) ;
  grbf_workspace_init_3d(ws, f, 3, nx, 3*ny*nz, ny, 3*nz, nz,
			       wt, 3, 3*ny*nz, 3*nz, al, N, work) ;
  
  for ( i = 0 ; i < nx ; i ++ ) {
    x[0] = limits[0] + h*i ;
    for ( j = 0 ; j < ny ; j ++ ) {
      x[1] = limits[2] + h*j ;
      for ( k = 0 ; k < nz ; k ++ ) {
	x[2] = limits[4] + h*k ;
	vortex_ring_func(r0, s0, 1, x, &(f[3*(i*ny*nz + j*nz + k)])) ;
      }
    }
  }
  
  fprintf(stderr, "evaluating Gaussian weights (FFT) ") ;
  t0 = g_timer_elapsed(timer, NULL) ;
  grbf_interpolation_weights_fft_3d(f, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
					  ws, wt, 3, 3*ny*nz, 3*nz) ;
  fprintf(stderr, "[%lg]\n", g_timer_elapsed(timer, NULL)-t0) ;

  np = grbf_grid_count_points_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
				       gtol) ;

  fprintf(stderr, "computational nodes: %d\n", np) ;
  d = rvpm_distribution_alloc(np) ;

  n = 0 ;
  for ( i = 0 ; i < nx ; i ++ ) {
    for ( j = 0 ; j < ny ; j ++ ) {
      for ( k = 0 ; k < nz ; k ++ ) {
	om = &(wt[3*(i*ny*nz + j*nz + k)]) ;
	if ( rvpm_vector_scalar(om,om) > gtol*gtol ) {
	  y = (gdouble *)rvpm_distribution_particle(d,n) ;
	  y[0] = limits[0] + h*i ;
	  y[1] = limits[2] + h*j ;
	  y[2] = limits[4] + h*k ;

	  zeta = M_PI*sqrt(M_PI)*s*s*s ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+0] = om[0]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+1] = om[1]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+2] = om[2]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_RADIUS] = s ;
	  n ++ ;
	}	  
      }
    }
  }
  rvpm_distribution_particle_number(d) = n ;
  g_assert(n == np) ;

  fprintf(stderr, "setting up FMM tree\n") ;
  tree = rvpm_tree_new(d, 4, 16, work) ;
  fprintf(stderr, "updating FMM tree\n") ;
  rvpm_tree_update(tree, 1, TRUE, work) ;
  
  /*velocity at sample points*/
  for ( i = 0 ; i < 1000 ; i ++ ) {
    x[0] = r0 ; x[2] = s0*0.5 ;
    x[1] = tree->origin[1] + 1e-2 + (tree->D-2.0*1e-2)*i/1000 ;
    
    memset(uv, 0, 12*sizeof(gdouble)) ;
    rvpm_vorticity_velocity_gradient(d, reg,
					   kernel,
					   x, uv, &(uv[3])) ;
    fprintf(stdout, "%1.16e %1.16e %1.16e", x[0], x[1], x[2]) ;
    for ( j = 0 ; j < 12 ; j ++ ) fprintf(stdout, " %1.16e", uv[j]) ;
    memset(ut, 0, 12*sizeof(gdouble)) ;
    rvpm_tree_velocity_gradient(tree, reg,
				      kernel,
				      x, ut, &(ut[3]), work) ;
    for ( j = 0 ; j < 12 ; j ++ ) fprintf(stdout, " %1.16e", ut[j]) ;
    fprintf(stdout, "\n") ;

    for ( j = 0 ; j < 12 ; j ++ ) {
      err[j] = MAX(err[j], fabs(ut[j]-uv[j])) ;
    }    
  }

  fprintf(stderr, "errors:") ;
  for ( j = 0 ; j < 12 ; j ++ ) {
    fprintf(stderr, " %e", err[j]) ;
  }

  fprintf(stderr, "\n") ;
  
  return ;
}

static void solver_test(gint np)

{
  gdouble s, *work, *wt, x[3], *om, *y, zeta, reg ;
  gdouble r0, s0, sc, limits[6], h, tol, t0, al, gtol, dt, *u ;
  char fname[64] ;
  grbf_workspace_t *ws ;
  rvpm_distribution_t *d ;
  rvpm_tree_t *tree ;
  gint i, j, k, n, N, wsize, sizew, nx, ny, nz, ustr, order_max ;
  gint nthreads = -1 ;
  FILE *fid ;
  
  al = 1.0 ; order_max = 10 ;
  
  tol = 1e-3 ; gtol = 1e-2 ; reg = 1e-1 ;
  r0 = 0.8 ; s0 = 0.1 ; ustr = 12 ;
  /*fit parameters*/
  sc = 16.0*s0 ;

  h = (2*r0 + 2*sc)/np ;

  s = h/al ;

  dt = 0.0125 ;
  
  fprintf(stderr, "basic solver test\n") ;
  fprintf(stderr, "=================\n") ;

  N = grbf_cardinal_function_length(al, tol) ;
  fprintf(stderr, "Gaussian fit parameter: alpha = %lg\n", al) ;
  fprintf(stderr, "Gaussian fit parameter: sigma = %lg\n", s) ;
  fprintf(stderr, "Gaussian fit parameter: N     = %d\n", N) ;
  /*workspace for setting up FFTs*/
  wsize = (N+1)*(N+1) + 3*(N+1) ;
  
  /*spacing parameters*/
  limits[0] = -r0 - sc ; limits[1] = r0 + sc ; 
  limits[2] = -r0 - sc ; limits[3] = r0 + sc ; 
  limits[4] =      -sc ; limits[5] =      sc ; 
  grbf_grid_adjust_3d(&(limits[0]), &(limits[1]), &nx,
			    &(limits[2]), &(limits[3]), &ny,
			    &(limits[4]), &(limits[5]), &nz, h) ;  
  fprintf(stderr, "ring parameters: s = %lg\n", s0) ;
  fprintf(stderr, "ring parameters: r = %lg\n", r0) ;

  fprintf(stderr, "domain limits x: %lg %lg (%d)\n",
	  limits[0], limits[1], nx) ;
  fprintf(stderr, "domain limits y: %lg %lg (%d)\n",
	  limits[2], limits[3], ny) ;
  fprintf(stderr, "domain limits z: %lg %lg (%d)\n",
	  limits[4], limits[5], nz) ;

  /*initialise the various operators*/
  sizew = wbfmm_element_number_rotation(2*order_max) ;
  sizew = MAX(sizew, (order_max+1)*(order_max+1)*3*16) ;
  wsize = MAX(sizew, wsize) ;
  
  fprintf(stderr, "allocating workspace of %lu bytes\n",
	  (wsize)*sizeof(gdouble)) ;
  work = (gdouble *)g_malloc0(wsize*sizeof(gdouble)) ;
  fprintf(stderr, "allocating workspace of %lu bytes\n",
	  nx*ny*nz*3*sizeof(gdouble)) ;
  wt   = (gdouble *)g_malloc0(nx*ny*nz*3*sizeof(gdouble)) ;

  ws = grbf_workspace_alloc(3, MAX(nx,MAX(ny,nz))) ;
  fprintf(stderr, "initialising workspace\n") ;
  grbf_workspace_init_3d(ws, wt, 3, nx, 3*ny*nz, ny, 3*nz, nz,
			       wt, 3, 3*ny*nz, 3*nz, al, N, work) ;
  
  for ( i = 0 ; i < nx ; i ++ ) {
    x[0] = limits[0] + h*i ;
    for ( j = 0 ; j < ny ; j ++ ) {
      x[1] = limits[2] + h*j ;
      for ( k = 0 ; k < nz ; k ++ ) {
	x[2] = limits[4] + h*k ;
	/* vortex_ring_func(r0, s0, 1, x, &(f[3*(i*ny*nz + j*nz + k)])) ; */
	vortex_ring_func(r0, s0, 1, x, &(wt[3*(i*ny*nz + j*nz + k)])) ;
      }
    }
  }
  
  fprintf(stderr, "evaluating Gaussian weights (FFT) ") ;
  t0 = g_timer_elapsed(timer, NULL) ;
  grbf_interpolation_weights_fft_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
					  ws, wt, 3, 3*ny*nz, 3*nz) ;
  fprintf(stderr, "[%lg]\n", g_timer_elapsed(timer, NULL)-t0) ;

  np = grbf_grid_count_points_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
				       gtol) ;

  fprintf(stderr, "computational nodes: %d\n", np) ;
  d = rvpm_distribution_alloc(np) ;

  n = 0 ;
  for ( i = 0 ; i < nx ; i ++ ) {
    for ( j = 0 ; j < ny ; j ++ ) {
      for ( k = 0 ; k < nz ; k ++ ) {
	om = &(wt[3*(i*ny*nz + j*nz + k)]) ;
	if ( rvpm_vector_scalar(om,om) > gtol*gtol ) {
	  y = (gdouble *)rvpm_distribution_particle(d,n) ;
	  y[0] = limits[0] + h*i ;
	  y[1] = limits[2] + h*j ;
	  y[2] = limits[4] + h*k ;

	  zeta = M_PI*sqrt(M_PI)*s*s*s ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+0] = om[0]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+1] = om[1]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_VORTICITY+2] = om[2]*zeta ;
	  y[RVPM_DISTRIBUTION_PARTICLE_RADIUS] = s ;
	  n ++ ;
	}	  
      }
    }
  }
  rvpm_distribution_particle_number(d) = n ;
  g_assert(n == np) ;

  /*set up the FMM solver for velocity evaluation*/
  /*should check the workspace size here to make sure it is large enough*/
  /* /\*initialise the various operators*\/ */
  /* sizew = wbfmm_element_number_rotation(2*order_max) ; */
  /* sizew = MAX(sizew, (order_max+1)*(order_max+1)*3*16) ; */
  /* work = (gdouble *)g_malloc0(2*sizew*sizeof(gdouble)) ; */

  fprintf(stderr, "setting up FMM tree\n") ;
  tree = rvpm_tree_new(d, 4, order_max, work) ;

  fid = file_open("points.dat", "w", "-", stdout) ;

  for ( i = 0 ; i < rvpm_distribution_particle_number(d) ; i ++ ) {
    y = (gdouble *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      fprintf(fid, "%e ", y[j]) ;
    }
    fprintf(fid, "\n") ;
  }
  
  file_close(fid) ;
  g_assert(np*3 < nx*ny*nz*3) ;
  
  /*time step for the crack*/
  u = (gdouble *)g_malloc0(ustr*np*sizeof(gdouble)) ;
  t0 = g_timer_elapsed(timer, NULL) ;
  for ( i = 0 ; i < 200 ; i ++ ) {
    fprintf(stderr, "time %d [%lg]\n", i,  g_timer_elapsed(timer, NULL)-t0) ;
    rvpm_tree_update(tree, nthreads, TRUE, work) ;
    memset(u, 0, ustr*np*sizeof(gdouble)) ;
    rvpm_tree_velocity_self(tree, reg, RVPM_KERNEL_WINCKELMANS_LEONARD,
				  u, ustr, &(u[3]), ustr, work) ;
    for ( j = 0 ; j < np ; j ++ ) {
      y = (gdouble *)rvpm_distribution_particle(d,j) ;
      y[0] += u[ustr*j+0]*dt ;
      y[1] += u[ustr*j+1]*dt ;
      y[2] += u[ustr*j+2]*dt ;
    }

    if ( 10*(i/10) == i ) {
      sprintf(fname, "points-%d.dat", i) ;
      fid = file_open(fname, "w", "-", stdout) ;

      distribution_write(fid, d) ;

      file_close(fid) ;
    }
  }

  return ;
}

static void permutation_test(gint np)

{
  rvpm_distribution_t *d, *e ;
  gdouble s = 0.1, *x ;
  gdouble *c ;
  gint i, j ;
  
  fprintf(stderr, "particle permutation test\n") ;
  fprintf(stderr, "========================================\n") ;

  d = rvpm_distribution_alloc(np) ;
  c = rvpm_distribution_origin(d) ;
  c[0] = -1 ; c[1] = -1 ; c[2] = -1 ; 
  rvpm_distribution_width(d) = 2 ;  

  distribution_random(d, np, s) ;
  fprintf(stderr, "particle number = %d\n",
	  np = rvpm_distribution_particle_number(d)) ;
  e = rvpm_distribution_alloc(np) ;

  rvpm_distribution_interleaved(d) = TRUE ;  
  rvpm_distribution_interleaved(e) = TRUE ;  
  memcpy(e->data, d->data, np*RVPM_DISTRIBUTION_PARTICLE_SIZE*sizeof(gdouble)) ;
  rvpm_distribution_particle_number(e) = np ;
  
  rvpm_distribution_interleaved_to_contiguous(e) ;  

  for ( i = 0 ; i < np ; i ++ ) {
    /* x = &(d->data[i*RVPM_DISTRIBUTION_PARTICLE_SIZE]) ; */
    x = (gdouble *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      if ( x[j] != e->data[j*np+i] ) {
	fprintf(stderr, "FAIL at point %d, element %d\n", i, j) ;
	return ;
      }
    }
  }

  fprintf(stderr, "interleaved to contiguous: PASS\n") ;
  
  rvpm_distribution_contiguous_to_interleaved(e) ;
  for ( i = 0 ; i < np ; i ++ ) {
    /* x = &(d->data[i*RVPM_DISTRIBUTION_PARTICLE_SIZE]) ; */
    x = (gdouble *)rvpm_distribution_particle(d,i) ;
    for ( j = 0 ; j < RVPM_DISTRIBUTION_PARTICLE_SIZE ; j ++ ) {
      if ( x[j] != e->data[i*RVPM_DISTRIBUTION_PARTICLE_SIZE+j] ) {
	fprintf(stderr, "FAIL at point %d, element %d\n", i, j) ;
	return ;
      }
    }
  }

  fprintf(stderr, "contiguous to interleaved: PASS\n") ;
  
  return ;
}

static void kernel_gradient_test(rvpm_kernel_t kernel)

{
  gdouble x[3], y[3], dy[3], x1[3], x2[3], K[12]={0}, K1[3], K2[3], ee, s ;
  gdouble err, w[3], du[9], u1[3], u2[3], cdu[9] ;
  rvpm_kernel_func_t kfunc ;
  gint i ;
  
  ee = 1e-9 ; s = 5e-1 ;

  fprintf(stderr, "kernel gradient test\n") ;
  fprintf(stderr, "====================\n") ;
  fprintf(stderr, "kernel: %s\n", rvpm_kernel_name(kernel)) ;

  /*this is to ensure the evaluation point is close enough to the
    source point to engage all terms in the analytical evaluation of
    the gradient (see the Gaussian smoothed kernel for example)*/
  dy[0] = 2.7*s ; dy[1] = -0.3*s ; dy[2] = 1.9*s ;
  dy[0] = 2.7e-3*s ; dy[1] = -0.3e-3*s ; dy[2] = 1.9e-3*s ;
  /* dy[2] = 1*s ; dy[1] = 1*s ; dy[0] = 1*s ; */
  dy[0] = dy[1] = dy[2] = 0.0 ;
  
  x[0] = 1.7 ; x[1] = -0.3 ; x[2] = 0.9 ;

  y[0] = x[0] + dy[0] ; 
  y[1] = x[1] + dy[1] ; 
  y[2] = x[2] + dy[2] ; 
  
  /* y[0] = 0.1 ; y[1] = 0.2 ; y[2] = -0.4 ; */
  w[0] = 0.5 ; w[1] = -0.1 ; w[2] = 0.3 ;

  kfunc = rvpm_kernel_func(kernel) ;
  
  kfunc(x, y, s, K, &(K[3])) ;
  rvpm_vector_cross(u1,K,w) ;

  fprintf(stderr, "u:    %lg %lg %lg\n", u1[0], u1[1], u1[2]) ;
  
  /*velocity gradient*/
  memset(du, 0, 9*sizeof(gdouble)) ;
  x1[0] = x[0] - ee/2 ; x1[1] = x[1] ; x1[2] = x[2] ;
  x2[0] = x[0] + ee/2 ; x2[1] = x[1] ; x2[2] = x[2] ;
  kfunc(x1, y, s, K1, NULL) ; kfunc(x2, y, s, K2, NULL) ;
  rvpm_vector_cross(u1,K1,w) ; rvpm_vector_cross(u2,K2,w) ;
  du[0] += (u2[0] - u1[0])/ee ;
  du[1] += (u2[1] - u1[1])/ee ;
  du[2] += (u2[2] - u1[2])/ee ;

  x1[0] = x[0] ; x1[1] = x[1] - ee/2 ; x1[2] = x[2] ;
  x2[0] = x[0] ; x2[1] = x[1] + ee/2 ; x2[2] = x[2] ;
  kfunc(x1, y, s, K1, NULL) ; kfunc(x2, y, s, K2, NULL) ;
  rvpm_vector_cross(u1,K1,w) ; rvpm_vector_cross(u2,K2,w) ;
  du[3] += (u2[0] - u1[0])/ee ;
  du[4] += (u2[1] - u1[1])/ee ;
  du[5] += (u2[2] - u1[2])/ee ;

  x1[0] = x[0] ; x1[1] = x[1] ; x1[2] = x[2] - ee/2 ;
  x2[0] = x[0] ; x2[1] = x[1] ; x2[2] = x[2] + ee/2 ;
  kfunc(x1, y, s, K1, NULL) ; kfunc(x2, y, s, K2, NULL) ;
  rvpm_vector_cross(u1,K1,w) ; rvpm_vector_cross(u2,K2,w) ;
  du[6] += (u2[0] - u1[0])/ee ;
  du[7] += (u2[1] - u1[1])/ee ;
  du[8] += (u2[2] - u1[2])/ee ;

  fprintf(stderr, "du/dx: %lg %lg %lg\n", du[0], du[1], du[2]) ;
  fprintf(stderr, "du/dy: %lg %lg %lg\n", du[3], du[4], du[5]) ;
  fprintf(stderr, "du/dz: %lg %lg %lg\n", du[6], du[7], du[8]) ;

  memset(cdu, 0, 9*sizeof(gdouble)) ;
  rvpm_vector_cross_gradient(cdu,&(K[3]),w) ;
  fprintf(stderr, "du/dx: %lg %lg %lg\n", cdu[0], cdu[1], cdu[2]) ;
  fprintf(stderr, "du/dy: %lg %lg %lg\n", cdu[3], cdu[4], cdu[5]) ;
  fprintf(stderr, "du/dz: %lg %lg %lg\n", cdu[6], cdu[7], cdu[8]) ;

  err = 0.0 ;
  for ( i = 0 ; i < 9 ; i ++ ) err = MAX(err,fabs(du[i]-cdu[i])) ;

  fprintf(stderr, "error:                 %lg\n", err) ;
  
  return ;
}

static void stream_function_test(void)

{
  gdouble r, z, ur, uz, rmax, zmin, zmax ;
  gdouble data[4] ;
  gqr_rule_t *gs, *gth ;
  gint ngs, ngth, i, j, nr, nz ;
  
  r = 0.6 ; z = 0.3 ;
  rmax = 2.5 ; zmin = -1.5 ; zmax = 1.5 ;

  ngs = 32 ; ngth = 16 ;

  gs = gqr_rule_alloc(ngs) ;
  gqr_rule_select(gs, GQR_GAUSS_LEGENDRE, ngs, NULL) ;
  gth = gqr_rule_alloc(ngth) ;
  gqr_rule_select(gth, GQR_GAUSS_LEGENDRE, ngth, NULL) ;

  data[0] = 0.9 ; data[1] = 0.1 ; data[2] = 0.3 ; data[3] = 1.0 ;

  nr = 32 ; nz = 64 ;
  for ( i = 0 ; i <= nr ; i ++ ) {
    r = 2.5*(gdouble)i/nr ;
    for ( j = 0 ; j <= nz ; j ++ ) {
      z = 3*(gdouble)j/nz - 1.5 ;
      rvpm_stream_func_velocity(rvpm_stream_func_vorticity_gaussian, data,
				      r, z, rmax, zmin, zmax, &ur, &uz,
				      gs, gth) ;
      fprintf(stdout, "%lg %lg %lg %lg\n", r, z, ur, uz) ;
    }
  }
  

  /* fprintf(stderr, "ur = %lg\n", ur) ; */
  
  return ;
}

static void stream_func(gdouble r, gdouble z,
			gdouble r1, gdouble z1,
			gdouble *G, gdouble *Gr, gdouble *Gz) 

{
  gdouble R1, R2, lm, lm_r, lm_z, K, E, dK, dE ;

  R1 = SQRT((r - r1)*(r - r1) + (z - z1)*(z - z1)) ;
  R2 = SQRT((r + r1)*(r + r1) + (z - z1)*(z - z1)) ;
  lm = (R2 - R1)/(R2 + R1) ;

  lm_r = 4.0*r1/R1/R2/(R1+R2)/(R1+R2)*(r1*r1 - r*r + (z-z1)*(z-z1)) ;
  lm_z = -2.0*lm*(z-z1)/R1/R2 ;
  
  rvpm_elliptic_KE(lm, &K, &E) ;  

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

static void stream_function_velocity_test(void)

{
  gdouble r, z, r1, z1, ee, G, Gr, Gz, G1, G2, gr, gz ;

  ee = 1e-6 ;
  
  r = 0.3 ; z = 0.9 ; z1 = 1.1 ; r1 = 0.35 ;
  
  stream_func(r-ee/2, z, r1, z1, &G1, &Gr, &Gz) ;
  stream_func(r+ee/2, z, r1, z1, &G2, &Gr, &Gz) ;
  gr = (G2 - G1)/ee ;
  stream_func(r, z-ee/2, r1, z1, &G1, &Gr, &Gz) ;
  stream_func(r, z+ee/2, r1, z1, &G2, &Gr, &Gz) ;
  gz = (G2 - G1)/ee ;
  
  stream_func(r, z, r1, z1, &G, &Gr, &Gz) ;

  fprintf(stderr, "radial: %lg %lg (%lg)\n",
	  Gr, gr, fabs(Gr-gr)) ;
  fprintf(stderr, " axial: %lg %lg (%lg)\n",
	  Gz, gz, fabs(Gz-gz)) ;
  
  return ;
}

static void solver_time_step_test(rvpm_time_step_t s, gdouble h)

{
  gdouble a[4], b[4], x[2], q[2], t ;
  gint n, ns, i, j, p ;

  rvpm_solver_coefficients(s, a, b, &n, &p) ;

  t = 0 ; x[0] = x[1] = 0 ; q[0] = q[1] = 0 ;
  ns = 5 ;
  
  for ( i = 0 ; i < ns ; i ++ ) {
    for ( j = 0 ; j < p ; j ++ ) {
      q[0] *= a[j] ; q[1] *= a[j] ;
      q[0] += h*pow(x[1],p-1) ; q[1] += h ;
      x[0] += b[j]*q[0] ; x[1] += b[j]*q[1] ;
    }
    t += h ;
  }    

  fprintf(stderr, "%lg %lg %lg %lg (%lg)\n",
	  t, x[1], x[0], pow(t,p)/(p), fabs(x[0] - pow(t,p)/p)) ;
  
  return ;
}

static void solver_step_test(void)

{
  gdouble h = 0.1 ;
  
  solver_time_step_test(RVPM_TIME_STEP_EULER, h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_2 , h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_3 , h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_4 , h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_5 , h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_6 , h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_7 , h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_8 , h) ;
  /* solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_10, h) ; */
  /* solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_11, h) ; */
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_12, h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_13, h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_14, h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_17, h) ;
  solver_time_step_test(RVPM_TIME_STEP_WILLIAMSON_19, h) ;
  
  return ;
}

static void gauss_fast_kernel_test(void)

{
  gdouble g, f, x[3], y[3], err, errd, Kr[12], Kf[12], R, s ;
  gint i ;
  
  fprintf(stderr, "interpolated Gaussian smoothing kernel test\n") ;
  fprintf(stderr, "===========================================\n") ;

  y[0] = 0.1 ; y[1] = 0.7 ; y[2] = -0.4 ;

  s = 0.079 ;
  
  err = errd = 0.0 ;
  for ( R = 0 ; R < 7.5*s ; R += 0.01*s ) {
    x[0] = y[0]+R/sqrt(3) ; x[1] = y[1]+R/sqrt(3) ; x[2] = y[2]-R/sqrt(3) ;

    rvpm_kernel_GS(x, y, s, Kr, &(Kr[3])) ;
    rvpm_kernel_fast_GS(x, y, s, Kf, &(Kf[3])) ;

    fprintf(stdout, "%1.16e", R) ;
    for ( i = 0 ; i < 12 ; i ++ ) {
      fprintf(stdout, " %1.16e", Kf[i]) ;
    }
    fprintf(stdout, "\n") ;
    err = MAX(err, ABS(Kf[0]-Kr[0])) ;
    err = MAX(err, ABS(Kf[1]-Kr[1])) ;
    err = MAX(err, ABS(Kf[2]-Kr[2])) ;
    for ( i = 0 ; i < 9 ; i ++ ) {
      errd = MAX(errd, ABS(Kr[3+i]-Kf[3+i])) ;
    }      
  }

  fprintf(stderr, "maximum error: %lg\n", err) ;
  fprintf(stderr, "maximum error: %lg\n", errd) ;
  
  return ;
}

gint main(gint argc, char **argv)

{
  char ch ;
  gint np, test ;
  rvpm_kernel_t kernel ;

  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;
  
  test = -1 ;
  np = 128 ;
  kernel = RVPM_KERNEL_GAUSSIAN ;
  
  while ( (ch = getopt(argc, argv, "k:n:t:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'k': kernel = rvpm_kernel_parse(optarg) ; break ;
    case 'n': np = atoi(optarg) ; break ;
    case 't': if ( (test = parse_test(optarg)) == -1 ) {
	fprintf(stderr, "unrecognised test \"%s\"\n", optarg) ;
	return 1 ;
      }
      break ;
    }
  }

  if ( test == -1 ) {
    fprintf(stderr, "%s: no test specified\n", progname) ;
    return 1 ;
  }

  if ( test == 0 ) {
    gauss_fit_test(np) ;

    return 0 ;
  }

  if ( test == 1 ) {
    permutation_test(np) ;

    return 0 ;
  }

  if ( test == 2 ) {
    tree_test(kernel, np) ;

    return 0 ;
  }

  if ( test == 3 ) {
    solver_test(np) ;

    return 0 ;
  }

  if ( test == 4 ) {
    kernel_gradient_test(kernel) ;

    return 0 ;
  }

  if ( test == 5 ) {
    stream_function_test() ;
    
    return 0 ;
  }

  if ( test == 6 ) {
    stream_function_velocity_test() ;
    
    return 0 ;
  }
  
  if ( test == 7 ) {
    solver_step_test() ;
    
    return 0 ;
  }

  if ( test == 8 ) {
    gauss_fast_kernel_test() ;

    return 0 ;
  }
  
  return 0 ;
}
