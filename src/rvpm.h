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

#ifndef __RVPM_H_INCLUDED__
#define __RVPM_H_INCLUDED__

#include <glib.h>

#include <gqr.h>

#include <wbfmm.h>

typedef enum {
  RVPM_KERNEL_UNDEFINED           = 0,
  RVPM_KERNEL_MOORE_ROSENHEAD     = 1,
  RVPM_KERNEL_WINCKELMANS_LEONARD = 2,
  RVPM_KERNEL_GAUSSIAN            = 3
} rvpm_kernel_t ;

typedef gint (*rvpm_kernel_func_t)(gdouble *, gdouble *, gdouble,
				   gdouble *, gdouble *) ;
typedef gint (*rvpm_kernel_func_f_t)(gfloat *, gfloat *, gfloat,
				     gfloat *, gfloat *) ;

typedef enum {
  RVPM_TIME_STEP_UNDEFINED      =  0,
  RVPM_TIME_STEP_EULER          =  1,
  RVPM_TIME_STEP_WILLIAMSON_1   =  2,
  RVPM_TIME_STEP_WILLIAMSON_2   =  3,
  RVPM_TIME_STEP_WILLIAMSON_3   =  4,
  RVPM_TIME_STEP_WILLIAMSON_4   =  5,
  RVPM_TIME_STEP_WILLIAMSON_5   =  6,
  RVPM_TIME_STEP_WILLIAMSON_6   =  7,
  RVPM_TIME_STEP_WILLIAMSON_7   =  8,
  RVPM_TIME_STEP_WILLIAMSON_8   =  9,
  RVPM_TIME_STEP_WILLIAMSON_12  = 10,
  RVPM_TIME_STEP_WILLIAMSON_13  = 11,
  RVPM_TIME_STEP_WILLIAMSON_14  = 12,
  RVPM_TIME_STEP_WILLIAMSON_17  = 13,
  RVPM_TIME_STEP_WILLIAMSON_19  = 14
} rvpm_time_step_t ;

typedef enum {
  RVPM_METHOD_CLASSICAL    = 0,
  RVPM_METHOD_REFORMULATED = 1,
} rvpm_method_t ;

typedef struct _rvpm_distribution_t rvpm_distribution_t ;
struct _rvpm_distribution_t {
  gsize size ;
  gboolean interleaved ;
  gint np, npmax ;
  gdouble c[3], D ;
  char *data ;
} ;

#define RVPM_DISTRIBUTION_PARTICLE_SIZE          7
#define RVPM_DISTRIBUTION_PARTICLE_X             0
#define RVPM_DISTRIBUTION_PARTICLE_VORTICITY     3
#define RVPM_DISTRIBUTION_PARTICLE_RADIUS        6

#define rvpm_distribution_particle(_d,_i)		\
  (&((_d)->data[((_i)*RVPM_DISTRIBUTION_PARTICLE_SIZE+	\
		 RVPM_DISTRIBUTION_PARTICLE_X)*((_d)->size)]))
#define rvpm_distribution_vorticity(_d,_i)		\
  (&((_d)->data[((_i)*RVPM_DISTRIBUTION_PARTICLE_SIZE+	\
		 RVPM_DISTRIBUTION_PARTICLE_VORTICITY)*((_d)->size)]))
#define rvpm_distribution_particle_radius(_d,_i)		\
  (&((_d)->data[((_i)*RVPM_DISTRIBUTION_PARTICLE_SIZE+	\
		 RVPM_DISTRIBUTION_PARTICLE_RADIUS)*((_d)->size)]))
#define rvpm_distribution_particle_number(_d)     ((_d)->np)
#define rvpm_distribution_particle_number_max(_d) ((_d)->npmax)
#define rvpm_distribution_origin(_d)     ((_d)->c)
#define rvpm_distribution_width(_d)      ((_d)->D)
#define rvpm_distribution_interleaved(_d) ((_d)->interleaved)

typedef struct _rvpm_tree_t rvpm_tree_t ;
struct _rvpm_tree_t {
  rvpm_distribution_t *d ;
  wbfmm_tree_t *t ;
  wbfmm_shift_operators_t *shifts ;
  gdouble origin[3], D ;
  gint depth, order_max ;
} ;

typedef struct _rvpm_solver_t rvpm_solver_t ;
struct _rvpm_solver_t {
  rvpm_kernel_t kernel ;
  rvpm_method_t method ;
  rvpm_time_step_t step ;
  gint nthreads ;
  gdouble reg, nu, f, g ;
} ;

#define rvpm_solver_kernel(_s)                ((_s)->kernel)
#define rvpm_solver_method(_s)                ((_s)->method)
#define rvpm_solver_time_step(_s)             ((_s)->step)
#define rvpm_solver_thread_number(_s)         ((_s)->nthreads)
#define rvpm_solver_regularisation(_s)        ((_s)->reg)
#define rvpm_solver_viscosity(_s)             ((_s)->nu)
#define rvpm_solver_model_parameter_f(_s)     ((_s)->f)
#define rvpm_solver_model_parameter_g(_s)     ((_s)->g)

typedef gdouble (*rvpm_stream_func_vorticity_t)(gdouble r, gdouble z,
						gpointer data) ;

rvpm_distribution_t *rvpm_distribution_alloc(gint np) ;
rvpm_distribution_t *rvpm_distribution_alloc_f(gint np) ;
gint rvpm_distribution_interleaved_to_contiguous(rvpm_distribution_t *d) ;
gint rvpm_distribution_interleaved_to_contiguous_f(rvpm_distribution_t *d) ;
gint rvpm_distribution_contiguous_to_interleaved(rvpm_distribution_t *d) ;
gint rvpm_distribution_contiguous_to_interleaved_f(rvpm_distribution_t *d) ;
gint rvpm_distribution_write(FILE *f, rvpm_distribution_t *d) ;
gint rvpm_distribution_write_f(FILE *f, rvpm_distribution_t *d) ;
rvpm_distribution_t *rvpm_distribution_read_alloc(FILE *f) ;
rvpm_distribution_t *rvpm_distribution_read_alloc_f(FILE *f) ;

gint rvpm_distribution_particle_add(rvpm_distribution_t *d,
				    gdouble *x, gdouble *w,
				    gdouble  s, gdouble G) ;
gint rvpm_distribution_particle_add_f(rvpm_distribution_t *d,
				      gfloat *x, gfloat *w,
				      gfloat  s, gfloat G) ;

gint rvpm_kernel_WL(gdouble *x, gdouble *y, gdouble e,
		    gdouble *K, gdouble *dK) ;
gint rvpm_kernel_WL_f(gfloat *x, gfloat *y, gfloat e, gfloat *K, gfloat *dK) ;
gint rvpm_kernel_MR(gdouble *x, gdouble *y, gdouble e,
		    gdouble *K, gdouble *dK) ;
gint rvpm_kernel_MR_f(gfloat *x, gfloat *y, gfloat e, gfloat *K, gfloat *dK) ;
gint rvpm_kernel_GS(gdouble *x, gdouble *y, gdouble s,
		    gdouble *K, gdouble *dK) ;
gint rvpm_kernel_GS_f(gfloat *x, gfloat *y, gfloat s, gfloat *K, gfloat *dK) ;

gint rvpm_vorticity_velocity_gradient(rvpm_distribution_t *v, gdouble e,
				      rvpm_kernel_t kernel,
				      gdouble *x, gdouble *u,
				      gdouble *du) ;
gint rvpm_vorticity_velocity_gradient_f(rvpm_distribution_t *v, gfloat e,
					rvpm_kernel_t kernel,
					gfloat *x, gfloat *u,
					gfloat *du) ;
gint rvpm_vorticity_derivatives(gdouble *G, gdouble s,
				gdouble f, gdouble g, gdouble nu,
				gdouble *du, gdouble *dG, gdouble *ds) ;
gint rvpm_vorticity_derivatives_f(gfloat *G, gfloat s,
				  gfloat f, gfloat g, gfloat nu,
				  gfloat *du, gfloat *dG, gfloat *ds) ;

rvpm_tree_t *rvpm_tree_new(rvpm_distribution_t *d, gint depth,
			   gint order_max, gdouble *work) ;
rvpm_tree_t *rvpm_tree_new_f(rvpm_distribution_t *d, gint depth,
			     gint order_max, gfloat *work) ;
gint rvpm_tree_update(rvpm_tree_t *t, gint nthreads, gboolean init_limits,
		      gdouble *work) ;
gint rvpm_tree_update_f(rvpm_tree_t *t, gint nthreads, gboolean init_limits,
			gfloat *work) ;

gint rvpm_tree_velocity_gradient(rvpm_tree_t *t, gdouble reg,
				 rvpm_kernel_t kernel,
				 gdouble *x, gdouble *u,
				 gdouble *du, gdouble *work) ;
gint rvpm_tree_velocity_gradient_f(rvpm_tree_t *t, gfloat reg,
				   rvpm_kernel_t kernel,
				   gfloat *x, gfloat *u,
				   gfloat *du, gfloat *work) ;
gint rvpm_tree_velocity(rvpm_tree_t *t, gdouble reg,
			rvpm_kernel_t kernel,
			gdouble *x, gdouble *u, gdouble *work) ;
gint rvpm_tree_velocity_f(rvpm_tree_t *t, gfloat reg,
			  rvpm_kernel_t kernel,
			  gfloat *x, gfloat *u, gfloat *work) ;
gint rvpm_tree_velocity_self(rvpm_tree_t *t, gdouble e,
			     rvpm_kernel_t kernel,
			     gdouble *u, gint ustr,
			     gdouble *du, gint dustr, gdouble *work) ;
gint rvpm_tree_velocity_self_f(rvpm_tree_t *t, gfloat e,
			       rvpm_kernel_t kernel,
			       gfloat *u, gint ustr,
			       gfloat *du, gint dustr, gfloat *work) ;
gint rvpm_tree_derivatives(rvpm_tree_t *t,
			   rvpm_solver_t *solver,
			   gdouble *du, gint str,
			   gdouble al,
			   gdouble *work) ;
gint rvpm_tree_derivatives_f(rvpm_tree_t *t,
			     rvpm_solver_t *solver,
			     gfloat *du, gint str,
			     gfloat al,
			     gfloat *work) ;

gint rvpm_solver_solve(rvpm_tree_t *tree, rvpm_solver_t *s,
		       gdouble t, gdouble dt, gdouble *u, gint ustr,
		       gdouble *work) ;
gint rvpm_solver_solve_f(rvpm_tree_t *tree, rvpm_solver_t *s,
			 gfloat t, gfloat dt, gfloat *u, gint ustr,
			 gfloat *work) ;
gint rvpm_solver_coefficients(rvpm_time_step_t s, gdouble *a, gdouble *b,
			      gint *n, gint *p) ;
gint rvpm_solver_coefficients_f(rvpm_time_step_t s, gfloat *a, gfloat *b,
				gint *n, gint *p) ;

rvpm_kernel_t rvpm_kernel_parse(char *str) ;
gint rvpm_kernels_list(FILE *f) ;
gchar *rvpm_kernel_name(rvpm_kernel_t kernel) ;
rvpm_kernel_func_t rvpm_kernel_func(rvpm_kernel_t kernel) ;
rvpm_kernel_func_f_t rvpm_kernel_func_f(rvpm_kernel_t kernel) ;

gint rvpm_stream_func_velocity(rvpm_stream_func_vorticity_t func,
			       gpointer data,
			       gdouble r, gdouble z,
			       gdouble rmax,
			       gdouble zmin,
			       gdouble zmax, 
			       gdouble *ur, gdouble *uz,
			       gqr_rule_t *gs,
			       gqr_rule_t *gth) ;
gint rvpm_stream_func_velocity_f(rvpm_stream_func_vorticity_t func,
				 gpointer data,
				 gfloat r, gfloat z,
				 gfloat rmax,
				 gfloat zmin,
				 gfloat zmax, 
				 gfloat *ur, gfloat *uz,
				 gqr_rule_t *gs,
				 gqr_rule_t *gth) ;
gdouble rvpm_stream_func_vorticity_gaussian(gdouble r, gdouble z,
					    gpointer data) ;

#endif /*__RVPM_H_INCLUDED__*/
