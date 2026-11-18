#include "fargo3d.h"

/* Model B of Paper II in two dimensions (cylindrical: X = phi, Y = R), as in
   the 2D runs of Huang et al. (2025).

   Gas:  Sigma_g = SIGMA0 (R/R0)^-SIGMASLOPE,  c_s^2 = h^2 GM/R,  h = ASPECTRATIO (R/R0)^FLARINGINDEX,
         e = Sigma_g c_s^2/(gamma - 1),
         v_R = -3/2 alpha gamma c_s^2/(Omega_K R), the steady viscous inflow
         for nu = alpha gamma c_s^2/Omega_K (adiabatic sound speed),
         v_phi = v_K sqrt(1 + h^2 dlnP/dlnR), with P = Sigma_g c_s^2,
   Dust: Sigma_d = EPSILON Sigma_g, with the steady drift of Nakagawa et al. (1986):
         v_R = 2 St/(St^2 + (1+eps)^2) eta v_K,
         v_phi = v_K + (1+eps) eta/((1+eps)^2 + St^2) v_K,
         eta = (h^2/2) dlnP/dlnR, St = 1/INVSTOKES1,
   where dlnP/dlnR = -SIGMASLOPE + 2 FLARINGINDEX - 1 (= -1.5). Azimuthal
   velocities are in the frame rotating at OMEGAFRAME. */

void _CondInit() {
  int i, j, k;
  real *rho = Density->field_cpu;
  real *e = Energy->field_cpu;
  real *vphi = Vx->field_cpu;
  real *vr = Vy->field_cpu;

  real R, omegaK, vK, h, cs2;
  real sigmag, sigmad;
  real vR_g, vphi_g, vR_d, vphi_d;
  real St, eta, eps_tot;
  const real dlnP_dlnR = -SIGMASLOPE + 2.0 * FLARINGINDEX - 1.0;

  i = j = k = 0;
  for (j = 0; j < Ny + 2 * NGHY; j++) {
    for (i = 0; i < Nx + 2 * NGHX; i++) {

      R = Ymed(j);
      omegaK = sqrt(G * MSTAR / (R * R * R));
      vK = omegaK * R;
      h = ASPECTRATIO * pow(R / R0, FLARINGINDEX);
      cs2 = h * h * vK * vK;

      sigmag = SIGMA0 * pow(R / R0, -SIGMASLOPE);
      sigmad = EPSILON * sigmag;

      vR_g = -1.5 * ALPHA * GAMMA * cs2 / vK;
      vphi_g = vK * sqrt(1.0 + h * h * dlnP_dlnR) - R * OMEGAFRAME;

      eps_tot = 1.0 + EPSILON;
      eta = 0.5 * h * h * dlnP_dlnR;
      St = 1.0 / INVSTOKES1;
      vR_d = (2.0 * St / (St * St + eps_tot * eps_tot)) * eta * vK;
      vphi_d = vK + (eps_tot * eta) / (eps_tot * eps_tot + St * St) * vK - R * OMEGAFRAME;

      if (Fluidtype == GAS) {
        rho[l] = sigmag;
#ifdef ADIABATIC
        e[l] = sigmag * cs2 / (GAMMA - 1.0);
#else
        e[l] = sqrt(cs2);
#endif
        vr[l] = vR_g;
        vphi[l] = vphi_g;
      }
      else { // DUST
        rho[l] = sigmad;
        e[l] = 0.0;
        vr[l] = vR_d;
        vphi[l] = vphi_d;
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
