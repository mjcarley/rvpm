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
#include <stdlib.h>
#include <math.h>

#include <glib.h>

#include "rvpm.h"
#include "rvpm-private.h"

/**
 *
 * @ingroup kernels
 * 
 * @{
 *
 */

/** 
 * Evaluate Winckelmans-Leonard kernel, equation 13 of Dufour, Pinon,
 * et al, https://dx.doi.org/10.1002/we.2905,
 * 
 * \f$K(\mathbf{x},\mathbf{y}) =  \left(K_{x}, K_{y}, K_{z}\right) =
 * -\frac{1}{4\pi}\frac{R^{2}+5\epsilon^{2}/2}{(R^{2}+\epsilon^{2})^{5/2}}
 * (\mathbf{x}-\mathbf{y})\f$
 * 
 * and its gradient
 * 
 * \f$\left(
 * \frac{\partial K_{x}}{\partial x}, 
 * \frac{\partial K_{y}}{\partial x}, 
 * \frac{\partial K_{z}}{\partial x}, 
 * \frac{\partial K_{x}}{\partial y}, 
 * \frac{\partial K_{y}}{\partial y}, 
 * \frac{\partial K_{z}}{\partial y}, 
 * \frac{\partial K_{x}}{\partial z}, 
 * \frac{\partial K_{y}}{\partial z}, 
 * \frac{\partial K_{z}}{\partial z}
 * \right)\f$
 * 
 * @param x field point;
 * @param y source point;
 * @param reg regularisation parameter \f$\epsilon\f$;
 * @param K on exit, contains value of kernel;
 * @param dK if not NULL, on exit, contains gradient of kernel.
 * 
 * @return 0 on success.
 */

gint RVPM_FUNCTION_NAME(rvpm_kernel_WL)(RVPM_REAL *x, RVPM_REAL *y,
					RVPM_REAL reg,
					RVPM_REAL *K, RVPM_REAL *dK)

{
  RVPM_REAL R2, R3, R5, R7, s, t, r[3] ;
  gint i ;
  
  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) + reg*reg ;
  
  R3 = R2*SQRT(R2) ; R5 = R3*R2 ;
  t = (R2 + 1.5*reg*reg)/R5 ;
  K[0] = -r[0]*t*0.25*M_1_PI ;
  K[1] = -r[1]*t*0.25*M_1_PI ;
  K[2] = -r[2]*t*0.25*M_1_PI ;
   
  if ( dK == NULL ) return 0 ;

  R7 = R5*R2 ;
  s = (3*R2 + 7.5*reg*reg)/R7 ;
  
  /*dK/dx*/
  dK[0] = -r[0]*r[0]*s + t ;
  dK[1] = -r[1]*r[0]*s ;
  dK[2] = -r[2]*r[0]*s ;
  /*dK/dy*/
  /* dK[3] = -r[0]*r[1]*s ; */
  dK[3] = dK[1] ;
  dK[4] = -r[1]*r[1]*s + t ;
  dK[5] = -r[2]*r[1]*s ;
  /*dK/dz*/
  /* dK[6] = -r[0]*r[2]*s ; */
  dK[6] = dK[2] ;
  /* dK[7] = -r[1]*r[2]*s ; */
  dK[7] = dK[5] ;
  dK[8] = -r[2]*r[2]*s + t ;
  
  for ( i = 0 ; i < 9 ; i ++ ) dK[i] *= -0.25*M_1_PI ;

  return 0 ;
}

/** 
 * Evaluate Moore-Rosenhead kernel, equation 12 of Dufour, Pinon,
 * et al, https://dx.doi.org/10.1002/we.2905,
 * 
 * \f$K(\mathbf{x},\mathbf{y}) = \left(K_{x}, K_{y}, K_{z}\right) =
 * -\frac{1}{4\pi}\frac{\mathbf{x}-\mathbf{y}}{(R^{2}+\epsilon^{2})^{3/2}}\f$
 * 
 * and its gradient
 * 
 * \f$\left(
 * \frac{\partial K_{x}}{\partial x}, 
 * \frac{\partial K_{y}}{\partial x}, 
 * \frac{\partial K_{z}}{\partial x}, 
 * \frac{\partial K_{x}}{\partial y}, 
 * \frac{\partial K_{y}}{\partial y}, 
 * \frac{\partial K_{z}}{\partial y}, 
 * \frac{\partial K_{x}}{\partial z}, 
 * \frac{\partial K_{y}}{\partial z}, 
 * \frac{\partial K_{z}}{\partial z}
 * \right)\f$
 * 
 * @param x field point;
 * @param y source point;
 * @param reg regularisation parameter \f$\epsilon\f$;
 * @param K on exit, contains value of kernel;
 * @param dK if not NULL, on exit, contains gradient of kernel.
 * 
 * @return 0 on success.
 */

gint RVPM_FUNCTION_NAME(rvpm_kernel_MR)(RVPM_REAL *x, RVPM_REAL *y,
					RVPM_REAL reg,
					RVPM_REAL *K, RVPM_REAL *dK)

{
  RVPM_REAL R2, R3, R5, r[3] ;
  gint i ;

  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) + reg*reg ;

  R3 = SQRT(R2)*R2 ;
  K[0] = -r[0]/R3*0.25*M_1_PI ;
  K[1] = -r[1]/R3*0.25*M_1_PI ;
  K[2] = -r[2]/R3*0.25*M_1_PI ;

  if ( dK == NULL ) return 0 ;

  R5 = R3*R2 ;
  
  /*dK/dx*/
  dK[0] = -3.0*r[0]*r[0]/R5 + 1.0/R3 ;
  dK[1] = -3.0*r[1]*r[0]/R5 ;
  dK[2] = -3.0*r[2]*r[0]/R5 ;
  /*dK/dy*/
  /* dK[3] = -3.0*r[0]*r[1]/R5 ; */
  dK[3] = dK[1] ;
  dK[4] = -3.0*r[1]*r[1]/R5 + 1.0/R3 ;
  dK[5] = -3.0*r[2]*r[1]/R5 ;
  /*dK/dz*/
  /* dK[6] = -3.0*r[0]*r[2]/R5 ; */
  dK[6] = dK[2] ;
  /* dK[7] = -3.0*r[1]*r[2]/R5 ; */
  dK[7] = dK[5] ;
  dK[8] = -3.0*r[2]*r[2]/R5 + 1.0/R3 ;

  for ( i = 0 ; i < 9 ; i ++ ) dK[i] *= -0.25*M_1_PI ;
  
  return 0 ;
}

static void kernel_GS(RVPM_REAL *x, RVPM_REAL *y,
		      RVPM_REAL s, RVPM_REAL *K)

{
  RVPM_REAL R, R2, R3, r[3], g, cutoff ;

  cutoff = 6.3*s ;
  /* cutoff = 3 ; */
  
  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) ;
  R  = SQRT(R2) ;
  R3 = R*R2 ;

  K[0] = -r[0]/R3*0.25*M_1_PI ;
  K[1] = -r[1]/R3*0.25*M_1_PI ;
  K[2] = -r[2]/R3*0.25*M_1_PI ;

  if ( R > cutoff ) return ;

  g = erf(R/s) - M_2_SQRTPI*R/s*exp(-R2/s/s) ;

  if ( R > 1e-9 ) {
    K[0] *= g ; K[1] *= g ; K[2] *= g ; 

    return ;
  }

  K[0] = -r[0]*0.25*M_1_PI*M_2_SQRTPI/s/s/s ;
  K[1] = -r[1]*0.25*M_1_PI*M_2_SQRTPI/s/s/s ;
  K[2] = -r[2]*0.25*M_1_PI*M_2_SQRTPI/s/s/s ;

  return ;
}

/** 
 * Evaluate Gaussian smoothed kernel
 * 
 * \f$K(\mathbf{x},\mathbf{y}) = \left(K_{x}, K_{y}, K_{z}\right) =
 * -\frac{g(R/\sigma)}{4\pi}
 * \frac{\mathbf{x}-\mathbf{y}}{|\mathbf{x}-\mathbf{y}|)^{3/2}}\f$ with 
 * \f$G(R/\sigma=\mathrm{erf}(R/\sigma) - 
 * 2/\sqrt{pi}(R/\sigma)\exp[-R^{2}/\sigma^{2}]\f$
 * 
 * and its gradient
 * 
 * \f$\left(
 * \frac{\partial K_{x}}{\partial x}, 
 * \frac{\partial K_{y}}{\partial x}, 
 * \frac{\partial K_{z}}{\partial x}, 
 * \frac{\partial K_{x}}{\partial y}, 
 * \frac{\partial K_{y}}{\partial y}, 
 * \frac{\partial K_{z}}{\partial y}, 
 * \frac{\partial K_{x}}{\partial z}, 
 * \frac{\partial K_{y}}{\partial z}, 
 * \frac{\partial K_{z}}{\partial z}
 * \right)\f$
 * 
 * @param x field point;
 * @param y source point;
 * @param s Gaussian parameter \f$\sigma\f$;
 * @param K on exit, contains value of kernel;
 * @param dK if not NULL, on exit, contains gradient of kernel.
 * 
 * @return 0 on success.
 */

gint RVPM_FUNCTION_NAME(rvpm_kernel_GS)(RVPM_REAL *x, RVPM_REAL *y,
					RVPM_REAL s,
					RVPM_REAL *K, RVPM_REAL *dK)

{
  RVPM_REAL R, R2, R3, R5, r[3], dk[9], dg[3], g, cutoff, E, errfunc ;
  gint i, j ;

  /*argument at which g(R/\sigma) reaches machine precision*/
#ifdef RVPM_SINGLE_PRECISION
  cutoff = 4.2 ;
#else /*RVPM_SINGLE_PRECISION*/
  cutoff = 6.2 ;
#endif /*RVPM_SINGLE_PRECISION*/
  
  K[0] = K[1] = K[2] = 0.0 ;
  
  rvpm_vector_diff(r,x,y) ;
  R2 = rvpm_vector_scalar(r,r) ;
  R  = SQRT(R2) ;

  if ( R < 1e-9 ) return 0 ;
  
  R3 = R*R2 ;

  if ( R/s > cutoff ) {
    g = 1 ;
    E = 0.0 ;
    errfunc = 1.0 ;
  } else {
    E = EXP(-R2/s/s) ;
    errfunc = ERF(R/s) ;
    g = errfunc - M_2_SQRTPI*R/s*E ;
  }

  /* if ( R > 1e-9 ) { */
  K[0] = -g*r[0]/R3*0.25*M_1_PI ;
  K[1] = -g*r[1]/R3*0.25*M_1_PI ;
  K[2] = -g*r[2]/R3*0.25*M_1_PI ;
  /* } */
  
  if ( dK == NULL ) return 0 ;

  R5 = R3*R2 ;
  
  /*dK/dx*/
  dk[0] = -3.0*r[0]*r[0]/R5 + 1.0/R3 ;
  dk[1] = -3.0*r[1]*r[0]/R5 ;
  dk[2] = -3.0*r[2]*r[0]/R5 ;
  /*dk/dy*/
  /* dk[3] = -3.0*r[0]*r[1]/R5 ; */
  dk[3] = dk[1] ;
  dk[4] = -3.0*r[1]*r[1]/R5 + 1.0/R3 ;
  dk[5] = -3.0*r[2]*r[1]/R5 ;
  /*dk/dz*/
  /* dk[6] = -3.0*r[0]*r[2]/R5 ; */
  dk[6] = dk[2] ;
  /* dk[7] = -3.0*r[1]*r[2]/R5 ; */
  dk[7] = dk[5] ;
  dk[8] = -3.0*r[2]*r[2]/R5 + 1.0/R3 ;

  /* M_2_SQRTPI*R/s*E ; */
  /* dg[0] = 4.0*r[0]/SQRT(M_PI)/s/s/s*R*E ; */
  /* dg[1] = 4.0*r[1]/SQRT(M_PI)/s/s/s*R*E ; */
  /* dg[2] = 4.0*r[2]/SQRT(M_PI)/s/s/s*R*E ; */
  dg[0] = 2.0*M_2_SQRTPI*R/s*E*r[0]/s/s ;
  dg[1] = 2.0*M_2_SQRTPI*R/s*E*r[1]/s/s ;
  dg[2] = 2.0*M_2_SQRTPI*R/s*E*r[2]/s/s ;

  for ( i = 0 ; i < 3 ; i ++ ) {
    for ( j = 0 ; j < 3 ; j ++ ) {
      dK[3*i+j]  = -0.25*M_1_PI*g*dk[3*i+j] ;
      dK[3*i+j] += K[j]*dg[i] ;
    }
  }
  
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_vorticity_velocity_gradient)(rvpm_distribution_t *v,
							  RVPM_REAL e,
							  rvpm_kernel_t kernel,
							  RVPM_REAL *x,
							  RVPM_REAL *u,
							  RVPM_REAL *du)

{
  gint i ;
  RVPM_REAL *y, *w, K[3], dK[9], KxW[9], s ;
  
  /*note that this function *increments* u and du*/
  if ( kernel == RVPM_KERNEL_UNDEFINED ) {
    g_error("%s: kernel must be defined on entry to velocity calculation",
	    __FUNCTION__) ;
  }

  if ( kernel != RVPM_KERNEL_MOORE_ROSENHEAD &&
       kernel != RVPM_KERNEL_WINCKELMANS_LEONARD &&
       kernel != RVPM_KERNEL_GAUSSIAN
       ) {
    g_error("%s: unknown kernel (%d)", __FUNCTION__, kernel) ;
  }

  if ( kernel == RVPM_KERNEL_MOORE_ROSENHEAD && du == NULL ) {
    for ( i = 0 ; i < rvpm_distribution_particle_number(v) ; i ++ ) {
      y = (RVPM_REAL *)rvpm_distribution_particle(v,i) ;
      w = (RVPM_REAL *)rvpm_distribution_vorticity(v,i) ;
      RVPM_FUNCTION_NAME(rvpm_kernel_MR)(x, y, e, K, NULL) ;
      rvpm_vector_cross(KxW,K,w) ;
      u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ; 
    }

    return 0 ;
  }

  if ( kernel == RVPM_KERNEL_MOORE_ROSENHEAD ) {
    for ( i = 0 ; i < rvpm_distribution_particle_number(v) ; i ++ ) {
      y = (RVPM_REAL *)rvpm_distribution_particle(v,i) ;
      w = (RVPM_REAL *)rvpm_distribution_vorticity(v,i) ;
      RVPM_FUNCTION_NAME(rvpm_kernel_MR)(x, y, e, K, dK) ;
      rvpm_vector_cross(KxW,K,w) ;
      u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ; 
      rvpm_vector_cross_gradient(KxW,dK,w) ;
      du[0] += KxW[0] ; du[1] += KxW[1] ; du[2] += KxW[2] ; 
      du[3] += KxW[3] ; du[4] += KxW[4] ; du[5] += KxW[5] ; 
      du[6] += KxW[6] ; du[7] += KxW[7] ; du[8] += KxW[8] ; 
    }

    return 0 ;
  }
  
  if ( kernel == RVPM_KERNEL_WINCKELMANS_LEONARD && du == NULL ) {
    for ( i = 0 ; i < rvpm_distribution_particle_number(v) ; i ++ ) {
      y = (RVPM_REAL *)rvpm_distribution_particle(v,i) ;
      w = (RVPM_REAL *)rvpm_distribution_vorticity(v,i) ;
      RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, y, e, K, NULL) ;
      rvpm_vector_cross(KxW,K,w) ;
      u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ; 
    }

    return 0 ;
  }

  if ( kernel == RVPM_KERNEL_WINCKELMANS_LEONARD ) {
    for ( i = 0 ; i < rvpm_distribution_particle_number(v) ; i ++ ) {
      y = (RVPM_REAL *)rvpm_distribution_particle(v,i) ;
      w = (RVPM_REAL *)rvpm_distribution_vorticity(v,i) ;
      RVPM_FUNCTION_NAME(rvpm_kernel_WL)(x, y, e, K, dK) ;
      rvpm_vector_cross(KxW,K,w) ;
      u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ; 
      rvpm_vector_cross_gradient(KxW,dK,w) ;
      du[0] += KxW[0] ; du[1] += KxW[1] ; du[2] += KxW[2] ; 
      du[3] += KxW[3] ; du[4] += KxW[4] ; du[5] += KxW[5] ; 
      du[6] += KxW[6] ; du[7] += KxW[7] ; du[8] += KxW[8] ; 
    }

    return 0 ;
  }

  if ( kernel == RVPM_KERNEL_GAUSSIAN && du == NULL ) {
    for ( i = 0 ; i < rvpm_distribution_particle_number(v) ; i ++ ) {
      y = (RVPM_REAL *)rvpm_distribution_particle(v,i) ;
      w = (RVPM_REAL *)rvpm_distribution_vorticity(v,i) ;
      s = *((RVPM_REAL *)rvpm_distribution_particle_radius(v,i)) ;
      kernel_GS(x, y, s, K) ;
      rvpm_vector_cross(KxW,K,w) ;
      u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ;
    }

    return 0 ;
  }

  if ( kernel == RVPM_KERNEL_GAUSSIAN ) {
    /* g_assert_not_reached() ; */
    for ( i = 0 ; i < rvpm_distribution_particle_number(v) ; i ++ ) {
      y = (RVPM_REAL *)rvpm_distribution_particle(v,i) ;
      w = (RVPM_REAL *)rvpm_distribution_vorticity(v,i) ;
      s = *((RVPM_REAL *)rvpm_distribution_particle_radius(v,i)) ;
      RVPM_FUNCTION_NAME(rvpm_kernel_GS)(x, y, e, K, dK) ;
      rvpm_vector_cross(KxW,K,w) ;
      u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ; 
      rvpm_vector_cross_gradient(KxW,dK,w) ;
      du[0] += KxW[0] ; du[1] += KxW[1] ; du[2] += KxW[2] ; 
      du[3] += KxW[3] ; du[4] += KxW[4] ; du[5] += KxW[5] ; 
      du[6] += KxW[6] ; du[7] += KxW[7] ; du[8] += KxW[8] ; 
      /* kernel_GS(x, y, s, K) ; */
      /* rvpm_vector_cross(KxW,K,w) ; */
      /* u[0] += KxW[0] ; u[1] += KxW[1] ; u[2] += KxW[2] ; */
    }

    return 0 ;
  }
  
  return 0 ;
}

gint RVPM_FUNCTION_NAME(rvpm_vorticity_derivatives)(RVPM_REAL *G, RVPM_REAL s,
						    RVPM_REAL f, RVPM_REAL g,
						    RVPM_REAL nu,
						    RVPM_REAL *du,
						    RVPM_REAL *dG,
						    RVPM_REAL *ds)

{
  RVPM_REAL dG0[3], dG1[3] ;
  RVPM_REAL absG, Gh[3], tmp0 ;

  *ds = 0.0 ;  
  absG = rvpm_vector_length(G) ;

  if ( absG < 1e-12 ) return 0 ;

  /*classical VPM term \Gamma.(\nabla u)*/
  dG0[0] =
    G[0]*du[RVPM_GRADIENT_U_X] +
    G[1]*du[RVPM_GRADIENT_U_Y] +
    G[2]*du[RVPM_GRADIENT_U_Z] ;
  dG0[1] =
    G[0]*du[RVPM_GRADIENT_V_X] +
    G[1]*du[RVPM_GRADIENT_V_Y] +
    G[2]*du[RVPM_GRADIENT_V_Z] ;
  dG0[2] =
    G[0]*du[RVPM_GRADIENT_W_X] +
    G[1]*du[RVPM_GRADIENT_W_Y] +
    G[2]*du[RVPM_GRADIENT_W_Z] ;
  
  /*Alvarez and Ning 2024, supplementary material, equations 15 and 16*/
  Gh[0] = G[0]/absG ; Gh[1] = G[1]/absG ; Gh[2] = G[2]/absG ;
  tmp0 = rvpm_vector_scalar(dG0,Gh) ;
  dG1[0] = tmp0*Gh[0] ;
  dG1[1] = tmp0*Gh[1] ; 
  dG1[2] = tmp0*Gh[2] ; 

  *ds = -(g + f)/(1.0 + 3*f)*tmp0*s/absG ;

  *ds += 2.0*nu/s ;
  
  dG[0] = dG0[0] - (g + f)/(1/3.0+f)*dG1[0] ;
  dG[1] = dG0[1] - (g + f)/(1/3.0+f)*dG1[1] ;
  dG[2] = dG0[2] - (g + f)/(1/3.0+f)*dG1[2] ;

  return 0 ;
}

/**
 *
 * @}
 *
 */
