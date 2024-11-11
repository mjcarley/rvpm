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

gint vortex_ring(gdouble *x, gdouble *p, gint np, gdouble *w, gboolean) ;

gpointer source_functions[] = {
  "vortex_ring", vortex_ring,
  NULL, NULL
} ;

static FILE *file_open(char *file, char *mode,
		       char *file_default, FILE *fdefault)

{
  FILE *f ;

  if ( file_default != NULL ) {
    if ( strcmp(file, file_default) == 0 ) return fdefault ;
  }
  
  f = fopen(file, mode) ;

  if ( f == NULL ) {
    fprintf(stderr, "%s: cannot open file \"%s\"", progname, file) ;
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

gint vortex_ring(gdouble *x, gdouble *p, gint np, gdouble *w, gboolean limits)

{
  gdouble r0, w0, z0, G, r, r2, th, z, s, ee ;

  if ( np < 2 ) {
    fprintf(stderr, "%s: vortex ring requires at least two parameters\n",
	    progname) ;
    return 1 ;
  }

  r0 = p[0] ; s = p[1] ;
  if ( np < 3 ) z0 = 0.0 ; else z0 = p[2] ;
  if ( np < 4 ) G = 1.0 ; else G = p[3] ;

  w0 = G/(M_PI*s*s) ;
  
  if ( limits ) {
    /*treat the first entry of w as a tolerance to be used in finding
      the limits of the grid box*/
    ee = w[0] ;
    r = sqrt(-s*s*log(ee/w0)) ;
    x[0] = -r0 - r ; x[1] =  r0 + r ;
    x[2] = -r0 - r ; x[3] =  r0 + r ;
    x[4] =  z0 - r ; x[5] =  z0 + r ;

    return 0 ;
  }

  r = sqrt(x[0]*x[0] + x[1]*x[1]) ;
  th = atan2(x[1], x[0]) ;
  z = x[2] ;

  r2 = (r - r0)*(r - r0) + (z-z0)*(z-z0) ;
  w0 *= exp(-r2/s/s) ;

  w[0] += -w0*sin(th) ;
  w[1] +=  w0*cos(th) ;
  w[2] +=  0.0 ;
  
  return 0 ;
}

static gpointer parse_func(char *str)

{
  gint i ;

  for ( i = 0 ; source_functions[i] != NULL ; i += 2 ) {
    if ( strcmp((char *)source_functions[i], str) == 0 ) {
      return source_functions[i+1] ;
    }
  }
  
  return NULL ;
}

static gint grid_limits(gdouble *wt, gdouble xmin, gdouble ymin, gdouble zmin,
			gint nx, gint ny, gint nz, gdouble h,
			gdouble *limits, gdouble tol)

{
  gint i, j, k, np ;
  gdouble x[3], *w, norm ;
  
  np = 0 ;

  limits[0] = limits[2] = limits[4] =  G_MAXDOUBLE ;
  limits[1] = limits[3] = limits[5] = -G_MAXDOUBLE ;
  
  for ( i = 0 ; i < nx ; i ++ ) {
    x[0] = xmin + h*i ;
    for ( j = 0 ; j < ny ; j ++ ) {
      x[1] = ymin + h*j ;
      for ( k = 0 ; k < nz ; k ++ ) {
	x[2] = zmin + h*k ;
	w = &(wt[3*(i*ny*nz + j*nz + k)]) ;
	norm = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ;
	if ( norm > tol*tol ) {
	  limits[0] = MIN(x[0],limits[0]) ;
	  limits[1] = MAX(x[0],limits[1]) ;
	  limits[2] = MIN(x[1],limits[2]) ;
	  limits[3] = MAX(x[1],limits[3]) ;
	  limits[4] = MIN(x[2],limits[4]) ;
	  limits[5] = MAX(x[2],limits[5]) ;
	  np ++ ;
	}
      }
    }
  }

  return np ;
}

static void write_active_nodes(FILE *f, gdouble *wt,
			       gdouble xmin, gdouble ymin, gdouble zmin,
			       gint nx, gint ny, gint nz, gdouble h, gdouble s,
			       gdouble tol)

{
  gint i, j, k ;
  gdouble x[3], *w, norm, sc ;

  sc = M_PI*sqrt(M_PI)*s*s*s ;

  for ( i = 0 ; i < nx ; i ++ ) {
    x[0] = xmin + h*i ;
    for ( j = 0 ; j < ny ; j ++ ) {
      x[1] = ymin + h*j ;
      for ( k = 0 ; k < nz ; k ++ ) {
	x[2] = zmin + h*k ;
	w = &(wt[3*(i*ny*nz + j*nz + k)]) ;
	norm = w[0]*w[0] + w[1]*w[1] + w[2]*w[2] ;
	if ( norm > tol*tol ) {
	  fprintf(f, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n",
		  x[0], x[1], x[2], w[0]*sc, w[1]*sc, w[2]*sc, s) ;
	}
      }
    }
  }

  return ;
}

static void expand_limits(gdouble *limits, gint n, gint N)

{
  gdouble xmin, xmax, d ;
  
  g_assert(N > n) ;

  d = (limits[1] - limits[0])/(n-1) ;
  /* fprintf(stderr, "%lg\n", d) ; */
  
  xmin = limits[0] + n/2*d - N/2*d ;
  xmax = xmin + (N-1)*d ;

  limits[0] = xmin ; limits[1] = xmax ;
  /* d = (limits[1] - limits[0])/(N-1) ; */
  /* fprintf(stderr, "%lg\n", d) ; */

  return ;
}

gint main(gint argc, char **argv)

{
  char ch, *rfile ;
  gint nparams[32] = {0}, nx, ny, nz, N, i, j, k, m, np, ngfunc ;
  gdouble al, sg, params[32][16], limits[6], flimits[6], tol, h, *wt, *work ;
  gdouble gtol, x[3], q[3], t0, smax, len, origin[3], s ;
  gpointer gfuncs[32] ;
  gint (*gfunc)(gdouble *, gdouble *, gint, gdouble *, gboolean) ;
  grbf_workspace_t *ws ;
  FILE *output, *input ;
  
  progname = g_strdup(g_path_get_basename(argv[0])) ;
  timer = g_timer_new() ;
  output = stdout ;
  
  al = 0.5 ; sg = 0.1 ; gfunc = NULL ; tol = 1e-6 ; gtol = 1e-3 ; h = 0 ;
  rfile = NULL ;
  
  ngfunc = 0 ;
  while ( (ch = getopt(argc, argv, "a:f:g:p:r:s:t:")) != EOF ) {
    switch ( ch ) {
    default: g_assert_not_reached() ; break ;
    case 'a': al = atof(optarg) ; break ;
    case 'g': gtol = atof(optarg) ; break ;
    case 'f':
      gfuncs[ngfunc] = parse_func(optarg) ;
      if ( gfuncs[ngfunc] == NULL ) {
	fprintf(stderr, "%s: unrecognised vorticity function \"%s\"\n",
		progname, optarg) ;
	return 1 ;
      }
      ngfunc ++ ;
      break ;
    case 'p':
      params[ngfunc-1][nparams[ngfunc-1]] = atof(optarg) ;
      nparams[ngfunc-1] ++ ;
      break ;
    case 'r': rfile = g_strdup(optarg) ; break ;
    case 's': sg = atof(optarg) ; break ; 
    case 't': tol = atof(optarg) ; break ;
    }
  }

  if ( rfile != NULL && gfunc != NULL ) {
    fprintf(stderr, "%s: only select one of -f (fitted source) "
	    "or -r (regridding)\n", progname) ;
    return 1 ;
  }

  N = grbf_cardinal_function_length(al, tol) ;
  
  h = al*sg ;
  if ( ngfunc != 0 ) {
    /*generating a new source: find the domain limits*/
    limits[0] = limits[2] = limits[4] =  G_MAXDOUBLE ;
    limits[1] = limits[3] = limits[5] = -G_MAXDOUBLE ;
    for ( k = 0 ; k < ngfunc ; k ++ ) {
      gfunc = gfuncs[k] ;
      if ( gfunc(flimits, params[k], nparams[k], &tol, TRUE) != 0 ) {
	return 1 ;
      }
      limits[0] = MIN(limits[0], flimits[0]) ;
      limits[1] = MAX(limits[1], flimits[1]) ;
      limits[2] = MIN(limits[2], flimits[2]) ;
      limits[3] = MAX(limits[3], flimits[3]) ;
      limits[4] = MIN(limits[4], flimits[4]) ;
      limits[5] = MAX(limits[5], flimits[5]) ;
    }
    fprintf(stderr,
	    "%s: initial domain:\n"
	    "    %lg %lg\n"
	    "    %lg %lg\n"
	    "    %lg %lg\n",
	    progname,
	    limits[0], limits[1], limits[2],
	    limits[3], limits[4], limits[5]) ;
    grbf_grid_adjust_3d(&(limits[0]), &(limits[1]), &nx,
			&(limits[2]), &(limits[3]), &ny,
			&(limits[4]), &(limits[5]), &nz, h) ;
    origin[0] = limits[0] ; origin[1] = limits[2] ; origin[2] = limits[4] ;
  }

  if ( rfile != NULL ) {
    input = file_open(rfile, "r", "-", stdin) ;

    fscanf(input, "%d", &np) ;
    for ( i = 0 ; i < 6 ; i ++ ) fscanf(input, "%lg", &(limits[i])) ;
    fscanf(input, "%lg", &smax) ;

    /*adjust the limits to ensure zero on the boundaries*/
    len = sqrt(-smax*smax*log(tol/1000)) + N*h ;
    limits[0] -= len ; limits[1] += len ; 
    limits[2] -= len ; limits[3] += len ; 
    limits[4] -= len ; limits[5] += len ; 
    fprintf(stderr,
	    "%s: initial domain:\n"
	    "    %lg %lg\n"
	    "    %lg %lg\n"
	    "    %lg %lg\n",
	    progname,
	    limits[0], limits[1], limits[2],
	    limits[3], limits[4], limits[5]) ;
    origin[0] = limits[0] ; origin[1] = limits[2] ; origin[2] = limits[4] ;
    grbf_grid_adjust_3d(&(limits[0]), &(limits[1]), &nx,
			&(limits[2]), &(limits[3]), &ny,
			&(limits[4]), &(limits[5]), &nz, h) ;
  }

  /*check that domain is large enough for number of weights in Gaussian,
    expanding about current midpoint if necessary*/
  if ( nx < 2*(N+1) ) {
    expand_limits(&(limits[0]), nx, 2*(N+1)) ;
    nx = 2*(N+1) ;
  }

  if ( ny < 2*(N+1) ) {
    expand_limits(&(limits[2]), ny, 2*(N+1)) ;
    ny = 2*(N+1) ;
  }
  
  if ( nz < 2*(N+1) ) {
    expand_limits(&(limits[4]), nz, 2*(N+1)) ;
    nz = 2*(N+1) ;
  }
  
  fprintf(stderr,
	  "%s: adjusted domain:\n"
	  "    %lg %lg (%d)\n"
	  "    %lg %lg (%d)\n"
	  "    %lg %lg (%d)\n",
	  progname,
	  limits[0], limits[1], nx,
	  limits[2], limits[3], ny,
	  limits[4], limits[5], nz) ;
  
  wt = (gdouble *)g_malloc0(3*nx*ny*nz*sizeof(gdouble)) ;
  ws = grbf_workspace_alloc(3, MAX(nx,MAX(ny,nz))) ;
  work = (gdouble *)g_malloc0(((N+1)*(N+1) + 3*(N+1))*sizeof(gdouble)) ;
  
  grbf_workspace_init_3d(ws, wt, 3, nx, 3*ny*nz, ny, 3*nz, nz,
			 wt, 3, 3*ny*nz, 3*nz, al, N, work) ;

  if ( ngfunc != 0 ) {
    fprintf(stderr, "%s: evaluating vorticity on grid ", progname) ;
    t0 = g_timer_elapsed(timer, NULL) ;
    for ( m = 0 ; m < ngfunc ; m ++ ) {
      gfunc = gfuncs[m] ;
      for ( i = 0 ; i < nx ; i ++ ) {
	x[0] = origin[0] + h*i ;
	for ( j = 0 ; j < ny ; j ++ ) {
	  x[1] = origin[1] + h*j ;
	  for ( k = 0 ; k < nz ; k ++ ) {
	    x[2] = origin[2] + h*k ;
	    gfunc(x, params[m], nparams[m],
		  &(wt[3*(i*ny*nz + j*nz + k)]), FALSE) ;
	  }
	}
      }
    }
    fprintf(stderr, "[%lg]\n", g_timer_elapsed(timer, NULL)-t0) ;    
  }

  if ( rfile != NULL ) {
    memset(wt, 0, 3*nx*ny*nz*sizeof(gdouble)) ;
    
    for ( i = 0 ; i < np ; i ++ ) {
      fscanf(input, "%lg %lg %lg", &(x[0]), &(x[1]), &(x[2])) ;
      fscanf(input, "%lg %lg %lg", &(q[0]), &(q[1]), &(q[2])) ;
      fscanf(input, "%lg", &s) ;
      q[0] /= M_PI*sqrt(M_PI)*s*s*s ;
      q[1] /= M_PI*sqrt(M_PI)*s*s*s ;
      q[2] /= M_PI*sqrt(M_PI)*s*s*s ;
      grbf_grid_increment_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
			     origin, h, x, s, q, tol) ;
    }
    
    file_close(input) ;
  }
  
  fprintf(stderr, "%s: fitting Gaussian RBFs ", progname) ;
  t0 = g_timer_elapsed(timer, NULL) ;   
  grbf_interpolation_weights_fft_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3,
				    ws, wt, 3, 3*ny*nz, 3*nz) ;
  fprintf(stderr, "[%lg]\n", g_timer_elapsed(timer, NULL)-t0) ;
  
  np = grbf_grid_count_points_3d(wt, 3, nx, 3*ny*nz, ny, 3*nz, nz, 3, gtol) ;

  fprintf(stderr, "%s: %d source nodes above threshold (%lg)\n",
	  progname, np, gtol) ;

  np = grid_limits(wt, origin[0], origin[1], origin[2],
		   nx, ny, nz, h, limits, gtol) ;
  
  fprintf(stderr,
	  "%s: active domain:\n"
	  "    %lg %lg\n"
	  "    %lg %lg\n"
	  "    %lg %lg\n",
	  progname,
	  limits[0], limits[1], limits[2],
	  limits[3], limits[4], limits[5]) ;

  fprintf(output, "%d %e %e %e %e %e %e %e\n", np,
	  limits[0], limits[1], limits[2], limits[3], limits[4], limits[5],
	  sg) ;

  write_active_nodes(output, wt,
		     origin[0], origin[1], origin[2],
		     nx, ny, nz, h, sg, gtol) ;
		     
  return 0 ;
}
