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

#ifndef __RVPM_PRIVATE_H_INCLUDED__
#define __RVPM_PRIVATE_H_INCLUDED__

#include <glib.h>

#ifdef RVPM_SINGLE_PRECISION

#define RVPM_REAL gfloat

#define RVPM_FUNCTION_NAME(_func) _func##_f

#define SQRT(_x) sqrtf((_x))
#define CBRT(_x) cbrtf((_x))
#define SIN(_x) sinf((_x))
#define COS(_x) cosf((_x))
#define ACOS(_x) acosf((_x))
#define ATAN(_x) atanf((_x))
#define ATAN2(_y,_x) atan2f((_y),(_x))
#define LOG(_x) logf((_x))
#define EXP(_x) expf((_x))
#define ERF(_x) erff((_x))

#else

#define RVPM_REAL gdouble

#define RVPM_FUNCTION_NAME(_func) _func

#define SQRT(_x) sqrt((_x))
#define CBRT(_x) cbrt((_x))
#define SIN(_x) sin((_x))
#define COS(_x) cos((_x))
#define ACOS(_x) acos((_x))
#define ATAN(_x) atan((_x))
#define ATAN2(_y,_x) atan2((_y),(_x))
#define LOG(_x) log((_x))
#define EXP(_x) exp((_x))
#define ERF(_x) erf((_x))

#endif /*RVPM_SINGLE_PRECISION*/

/*r = x - y*/
#define rvpm_vector_diff(_r,_x,_y)		\
  do {						\
    (_r)[0] = (_x)[0] - (_y)[0] ;		\
    (_r)[1] = (_x)[1] - (_y)[1] ;		\
    (_r)[2] = (_x)[2] - (_y)[2] ;		\
  } while (0)

#define rvpm_vector_scalar(_x,_y)		\
  ((_x)[0]*(_y)[0]+(_x)[1]*(_y)[1]+(_x)[2]*(_y)[2])

#define rvpm_vector_cross(_u,_x,_y)					\
  do {									\
    (_u)[0] = (_x)[1]*(_y)[2] - (_x)[2]*(_y)[1] ;			\
    (_u)[1] = (_x)[2]*(_y)[0] - (_x)[0]*(_y)[2] ;			\
    (_u)[2] = (_x)[0]*(_y)[1] - (_x)[1]*(_y)[0] ;			\
  } while (0)

/*
 * gradient of a cross product: du = \nabla K(x)\times y
 *
 * input: dx = K_{xx} K_{xy} K_{xz} K_{yx} K_{yy} K_{yz} K_{zx} K_{zy} K_{zz}
 * 
 * ouput: 
 * 
 */

#define rvpm_vector_cross_gradient(_du,_dx,_y)				\
  do {									\
  (_du)[0] = (_dx)[1]*(_y)[2] - (_dx)[2]*(_y)[1] ;			\
  (_du)[1] = (_dx)[2]*(_y)[0] - (_dx)[0]*(_y)[2] ;			\
  (_du)[2] = (_dx)[0]*(_y)[1] - (_dx)[1]*(_y)[0] ;			\
  (_du)[3] = (_dx)[4]*(_y)[2] - (_dx)[5]*(_y)[1] ;			\
  (_du)[4] = (_dx)[5]*(_y)[0] - (_dx)[3]*(_y)[2] ;			\
  (_du)[5] = (_dx)[3]*(_y)[1] - (_dx)[4]*(_y)[0] ;			\
  (_du)[6] = (_dx)[7]*(_y)[2] - (_dx)[8]*(_y)[1] ;			\
  (_du)[7] = (_dx)[8]*(_y)[0] - (_dx)[6]*(_y)[2] ;			\
  (_du)[8] = (_dx)[6]*(_y)[1] - (_dx)[7]*(_y)[0] ;			\
  } while (0)


#define rvpm_vector_length(_x) sqrt(rvpm_vector_scalar(_x,_x))

#define rvpm_vector_distance2(_x,_y)		\
  (((_x)[0] - (_y)[0])*((_x)[0] - (_y)[0])	\
   + ((_x)[1] - (_y)[1])*((_x)[1] - (_y)[1])	\
   + ((_x)[2] - (_y)[2])*((_x)[2] - (_y)[2]))

gint rvpm_elliptic_KE(gdouble k, gdouble *K, gdouble *E) ;
gint rvpm_elliptic_KE_f(gfloat k, gfloat *K, gfloat *E) ;

#endif /*__RVPM_PRIVATE_H_INCLUDED__*/

