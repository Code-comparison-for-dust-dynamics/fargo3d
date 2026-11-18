//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"

/* Model B damping zones, as in Athena++ (de Val-Borro et al. 2006): the gas
   relaxes toward its initial state and the dust toward the steady drift
   solution. Each step a cell moves the fraction 1 - exp(-DAMPING_RATE*Omega_K
   *R(R_cyl)*dt) of the way, where R is a parabolic profile in cylindrical
   radius across each radial zone. */
#ifndef DAMPING_RATE
#define DAMPING_RATE 1.0
#endif
#ifndef DUST_HRATIO
#define DUST_HRATIO 1.0
#endif
//<\INCLUDES>

void StockholmBoundary_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Density);
  INPUT2D(Density0);
  OUTPUT(Density);
#ifdef ADIABATIC
  INPUT(Energy);
  INPUT2D(Energy0);
  OUTPUT(Energy);
#endif
#ifdef X
  INPUT(Vx);
  INPUT2D(Vx0);
  OUTPUT(Vx);
#endif
#ifdef Y
  INPUT(Vy);
  INPUT2D(Vy0);
  OUTPUT(Vy);
#endif
#ifdef Z
  INPUT(Vz);
  INPUT2D(Vz0);
  OUTPUT(Vz);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho  = Density->field_cpu;
  real* rho0 = Density0->field_cpu;
#ifdef X
  real* vx  = Vx->field_cpu;
  real* vx0 = Vx0->field_cpu;
#endif
#ifdef Y
  real* vy  = Vy->field_cpu;
  real* vy0 = Vy0->field_cpu;
#endif
#ifdef Z
  real* vz  = Vz->field_cpu;
  real* vz0 = Vz0->field_cpu;
#endif
#ifdef ADIABATIC
  real* e    = Energy->field_cpu;
  real* e0   = Energy0->field_cpu;
#endif
  int pitch   = Pitch_cpu;
  int stride  = Stride_cpu;
  int size_x  = Nx+2*NGHX;
  int size_y  = Ny+2*NGHY;
  int size_z  = Nz+2*NGHZ;
  int pitch2d = Pitch2D;
  real y_min = YMIN;
  real y_max = YMAX;
  real z_min = ZMIN;
  real z_max = ZMAX;
  real dampingzone = DAMPINGZONE;
  real kbcol = KILLINGBCCOLATITUDE;
  real of    = OMEGAFRAME;
  real of0   = OMEGAFRAME0;
  real r0 = R0;
  int periodic_z = PERIODICZ;
  int fluidtype = Fluidtype;
  real aspect_ratio = ASPECTRATIO;
  real flaring_index = FLARINGINDEX;
  real damping_rate_val = DAMPING_RATE;
  real epsilon = EPSILON;
  real sigma0 = SIGMA0;
  real sigmaslope = SIGMASLOPE;
  real invstokes1 = INVSTOKES1;
  real dust_hratio = DUST_HRATIO;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  //  Similar to Benitez-Llambay et al. (2016), Eq. 7.
  real Y_inf = y_min*pow(dampingzone, 2.0/3.0);
  real Y_sup = y_max*pow(dampingzone,-2.0/3.0);
  real Z_inf = z_min - (z_max-z_min); // Here we push Z_inf & Z_sup
  real Z_sup = z_max + (z_max-z_min); // out of the mesh
#ifdef CYLINDRICAL
  Z_inf = z_min + (z_max-z_min)*0.1;
  Z_sup = z_max - (z_max-z_min)*0.1;
  if (periodic_z) { // Push Z_inf & Z_sup out of mesh if periodic in Z
    Z_inf = z_min-r0;
    Z_sup = z_max+r0;
  }
#endif
#ifdef SPHERICAL
  Z_inf = M_PI/2.0-(M_PI/2.0-z_min)*(1.0-kbcol);
  Z_sup = M_PI/2.0+(M_PI/2.0-z_min)*(1.0-kbcol); // Avoid damping in ghost zones
  // if only half upper disk is covered by the mesh
#endif
  real radius;
  real vx0_target;
  real rampy;
  real rampyi;
  real rampz;
  real rampzz;
  real rampi;
  real ramp;
  real omega_dyn_y;
//<\INTERNAL>

//<CONSTANT>
// real xmin(Nx+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for (k=0; k<size_z; k++) {
#endif
#ifdef Y
    for (j=0; j<size_y; j++) {
#endif
#ifdef X
      for (i=0; i<size_x; i++) {
#endif
//<#>
	rampy = 0.0;
	rampyi = 0.0;
	rampz = 0.0;
	rampzz = 0.0;
	omega_dyn_y = 0.0;
#ifdef Y
	/* Radial profile in cylindrical radius at the cell centre (rampy) and at
	   the colatitude face of v_theta (rampyi) */
	real R_cyl_r = ymed(j);
	real R_cyl_i = ymed(j);
#ifdef SPHERICAL
	R_cyl_r *= sin(zmed(k));
	R_cyl_i *= sin(zmin(k));
#endif
	if (R_cyl_r < 1e-12) R_cyl_r = 1e-12;
	omega_dyn_y = sqrt(G * MSTAR / (R_cyl_r * R_cyl_r * R_cyl_r));

	if(R_cyl_r > Y_sup) rampy  = pow((R_cyl_r - Y_sup) / (y_max - Y_sup), 2.0);
	if(R_cyl_r < Y_inf) rampy  = pow((Y_inf - R_cyl_r) / (Y_inf - y_min), 2.0);
	if(R_cyl_i > Y_sup) rampyi = pow((R_cyl_i - Y_sup) / (y_max - Y_sup), 2.0);
	if(R_cyl_i < Y_inf) rampyi = pow((Y_inf - R_cyl_i) / (Y_inf - y_min), 2.0);
#endif
#ifdef Z
	if(zmed(k) > Z_sup) {
	  rampz   = (zmed(k)-Z_sup)/(z_max-Z_sup);
	}
	if(zmed(k) < Z_inf) {
	  rampz   = (Z_inf-zmed(k))/(Z_inf-z_min);
	}
	rampz = rampz * rampz;		/* vertical ramp in X^2 */
	if(zmin(k) > Z_sup) {
	  rampzz  = (zmin(k)-Z_sup)/(z_max-Z_sup);
	}
	if(zmin(k) < Z_inf) {
	  rampzz  = (Z_inf-zmin(k))/(Z_inf-z_min);
	}
	rampzz= rampzz * rampzz;		/* vertical ramp in X^2 */
#endif
	if (periodic_z) {
	  rampz = 0.0;
	  rampzz = 0.0;
	}
	ramp = rampy+rampz;
	rampi= rampyi+rampzz;
	if (ramp > 0.0) {
	  real w = 1.0 - exp(-damping_rate_val * omega_dyn_y * ramp * dt);

	  if (fluidtype == GAS) {
	    rho[l] += (rho0[l2D] - rho[l]) * w;

#ifdef ADIABATIC
	    e[l]   += (e0[l2D]   - e[l])   * w;
#endif

#ifdef X
	    vx0_target = vx0[l2D];
	    radius = ymed(j);
#ifdef SPHERICAL
	    radius *= sin(zmed(k));
#endif
	    vx0_target -= (of - of0) * radius;
	    vx[l] += (vx0_target - vx[l]) * w;
#endif
#ifdef Y
	    vy[l] += (vy0[l2D] - vy[l]) * w;
#endif
	  } else {
	    /* Dust: toward the analytic profile and the steady drift velocities
	       of Nakagawa et al. (1986) */
	    real r_sph = ymed(j);
	    real theta = 0.5*M_PI;
	    real R_cyl = r_sph;
#ifdef SPHERICAL
	    theta = zmed(k);
	    R_cyl = r_sph * sin(theta);
#endif
	    if (R_cyl < 1.0e-12) R_cyl = 1.0e-12;

	    {
	      real omegaK = sqrt(G * MSTAR / (R_cyl * R_cyl * R_cyl));
	      real vK = omegaK * R_cyl;
	      real H_over_R = aspect_ratio * pow(R_cyl / r0, flaring_index);
#ifdef Z
	      /* 3D: midplane pressure slope and the vertical dust profile */
	      const real dlnP_dlnR = -11.0/4.0;
	      real cs2 = aspect_ratio * aspect_ratio *
	                 pow(R_cyl / r0, 2.0 * flaring_index) *
	                 (G * MSTAR / R_cyl);
	      real dust_rho_target = (epsilon * pow(R_cyl / r0, -2.25) / dust_hratio) *
	                             exp((G * MSTAR) / (dust_hratio * dust_hratio * cs2) *
	                                 (1.0 / r_sph - 1.0 / R_cyl));
#else
	      /* 2D: slope of P = Sigma c_s^2, and EPSILON times the initial Sigma_g */
	      const real dlnP_dlnR = -sigmaslope + 2.0 * flaring_index - 1.0;
	      real dust_rho_target = epsilon * sigma0 * pow(R_cyl / r0, -sigmaslope);
#endif
	      if (dust_rho_target < 1.0e-20) dust_rho_target = 1.0e-20;

	      {
	        real eps_tot = 1.0 + epsilon;
	        real eta = 0.5 * H_over_R * H_over_R * dlnP_dlnR;
	        real St = 1.0 / invstokes1;
	        real vR_d = (2.0 / (St + eps_tot * eps_tot / St)) * eta * vK;
	        real vphi_d = (1.0 + (eps_tot * eta) / (eps_tot * eps_tot + St * St)) * vK;

	        rho[l] += (dust_rho_target - rho[l]) * w;
#ifdef X
	        vx[l] += ((vphi_d - of * R_cyl) - vx[l]) * w;
#endif
#ifdef Y
#ifdef SPHERICAL
	        vy[l] += (vR_d * sin(theta) - vy[l]) * w;
#else
	        vy[l] += (vR_d - vy[l]) * w;
#endif
#endif
#ifdef Z
	        vz[l] += (vR_d * cos(theta) - vz[l]) * w;
#endif
	      }
	    }
	  }
	}
#ifdef Z
	if (rampi > 0.0 && fluidtype == GAS) {
	  real w = 1.0 - exp(-damping_rate_val * omega_dyn_y * rampi * dt);
	  vz[l] += (vz0[l2D] - vz[l]) * w;
	}
#endif
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
