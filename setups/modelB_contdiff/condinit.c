#include "fargo3d.h"

/* Model B of Paper II, initial state (as in the Athena++ problem generator).
   Spherical mesh: Y = r, Z = colatitude; R = r sin(theta), z = r cos(theta).

   Gas:  rho_g = (R/R0)^p exp[GM/c_s^2 (1/r - 1/R)],  c_s^2 = h0^2 GM (R/R0)^q,
         v_R = -3/2 alpha c_s^2/(Omega_K R),
         v_phi = sqrt(GM/R) sqrt((p+q) c_s^2 R/GM + 1 + q - q R/r),
   Dust: rho_d = EPSILON rho_g, with the steady drift of Nakagawa et al. (1986):
         v_R = 2 St/(St^2 + (1+eps)^2) eta v_K,
         v_phi = v_K + (1+eps) eta/((1+eps)^2 + St^2) v_K,
         eta = (h^2/2)(p + q + 2), St = 1/INVSTOKES1,
   with p = -2.25 and q = -0.5. Azimuthal velocities are in the frame rotating
   at OMEGAFRAME, and v_R is projected on r and theta. */

void _CondInit() {
  int i, j, k;
  real *rho = Density->field_cpu;
  real *e = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr = Vy->field_cpu;
  real *vtheta = Vz->field_cpu;

  real r_sph, theta, R_cyl, z_cyl;
  real omegaK, vK, H_over_R, cs2_loc;
  real rhog, rhod;
  real vR_g, vphi_g, vR_d, vphi_d;
  real St, eta, eps_tot;
  const real pvalue = -2.25;   // midplane density slope
  const real qvalue = -0.5;    // c_s^2 slope
  real safe_rho_floor = 1.0e-12;

  i = j = k = 0;
  for (k = 0; k < Nz + 2 * NGHZ; k++) {
    for (j = 0; j < Ny + 2 * NGHY; j++) {
      for (i = 0; i < Nx + 2 * NGHX; i++) {

        r_sph = Ymed(j);
        theta = Zmed(k);

        R_cyl = r_sph * sin(theta);
        z_cyl = r_sph * cos(theta);

        if (R_cyl < 1.0e-6) R_cyl = 1.0e-6;

        omegaK = sqrt(G * MSTAR / (R_cyl * R_cyl * R_cyl));
        vK = omegaK * R_cyl;

        H_over_R = ASPECTRATIO * pow(R_cyl / R0, FLARINGINDEX);

        real cs2_0_ref = ASPECTRATIO * ASPECTRATIO * G * MSTAR;
        cs2_loc = cs2_0_ref * pow(R_cyl / R0, qvalue);

        rhog = 1.0 * pow(R_cyl / R0, pvalue) *
               exp((G * MSTAR / cs2_loc) * (1.0 / r_sph - 1.0 / R_cyl));

        if (rhog < safe_rho_floor) rhog = safe_rho_floor;
        rhod = EPSILON * rhog;

        vR_g = -1.5 * ALPHA * cs2_loc * R_cyl / sqrt(G * MSTAR * R_cyl);

        {
          real term1 = (pvalue + qvalue) * cs2_loc / (G * MSTAR / R_cyl);
          real term2 = 1.0 + qvalue;
          real term3 = qvalue * R_cyl / r_sph;
          real bracket = term1 + term2 - term3;
          if (bracket < 0.0) bracket = 0.0;
          vphi_g = sqrt(G * MSTAR / R_cyl) * sqrt(bracket) - R_cyl * OMEGAFRAME;
        }

        eps_tot = 1.0 + EPSILON;
        eta = 0.5 * H_over_R * H_over_R * (pvalue + qvalue + 2.0);
        St = 1.0 / INVSTOKES1;
        vR_d = (2.0 * St / (St * St + eps_tot * eps_tot)) * eta * vK;
        vphi_d = vK + (eps_tot * eta) / (eps_tot * eps_tot + St * St) * vK;

        real vr_final, vtheta_final, vphi_final;

        if (Fluidtype == GAS) {
           rho[l] = rhog;
#ifdef ADIABATIC
           e[l] = rhog * cs2_loc / (GAMMA - 1.0);
#else
           e[l] = sqrt(cs2_loc);
#endif

           vr_final     = vR_g * sin(theta);
           vtheta_final = vR_g * cos(theta);
           vphi_final   = vphi_g;
        }
        else { // DUST
           rho[l] = rhod;
           e[l] = 0.0;

           vr_final     = vR_d * sin(theta);
           vtheta_final = vR_d * cos(theta);
           vphi_final   = vphi_d - R_cyl * OMEGAFRAME;
        }

        vr[l] = vr_final;
        vtheta[l] = vtheta_final;
        vphi[l] = vphi_final;
      }
    }
  }
}

void CondInit()
{

  int id_gas = 0;
  int feedback = YES;

  /* Gas fluid */
  Fluids[id_gas] = CreateFluid("gas", GAS);
  SelectFluid(id_gas);
  _CondInit();

  /* Dust fluid(s) */
  char dust_name[MAXNAMELENGTH];
  int id_dust;

  for (id_dust = 1; id_dust < NFLUIDS; id_dust++)
  {
    sprintf(dust_name, "dust%d", id_dust);
    Fluids[id_dust] = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    _CondInit();
  }

  /* Dust-gas drag with back-reaction */
  ColRate(INVSTOKES1, id_gas, 1, feedback);
}
