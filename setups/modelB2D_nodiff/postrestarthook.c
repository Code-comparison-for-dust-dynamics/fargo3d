#include "fargo3d.h"
extern int NbRestart;

/* Dust release in 2D (run with -p): restarting from output 20, the end of the
   gas stage at 100 orbits, sets Sigma_d to EPSILON times the local Sigma_g,
   so the gap and the ring are kept, with the steady drift velocities of
   Nakagawa et al. (1986) for P = Sigma_g c_s^2. Any other restart leaves the
   dust alone. */

void PostRestartHook() {
    if (NFLUIDS < 2) return;
    if (NbRestart != 20) return;

    masterprint("Releasing the dust at output 20\n");

    real *sigma_g = Fluids[0]->Density->field_cpu;
    real *sigma_d = Fluids[1]->Density->field_cpu;
    real *vx_d    = Fluids[1]->Vx->field_cpu;
    real *vy_d    = Fluids[1]->Vy->field_cpu;

#ifdef ADIABATIC
    real *e_d = NULL;
    if (Fluids[1]->Energy != NULL) {
        e_d = Fluids[1]->Energy->field_cpu;
    }
#endif

    int i, j;
    real R, omegaK, vK, h, eta, St, eps_tot, vR_d, vphi_d, sg;
    const real dlnP_dlnR = -SIGMASLOPE + 2.0 * FLARINGINDEX - 1.0;
    int pitch = Pitch_cpu;

    for (j = 0; j < Ny + 2 * NGHY; j++) {
        for (i = 0; i < Nx + 2 * NGHX; i++) {
            int idx = i + j * pitch;

            R = Ymed(j);
            omegaK = sqrt(G * MSTAR / (R * R * R));
            vK = omegaK * R;
            h = ASPECTRATIO * pow(R / R0, FLARINGINDEX);

            sg = sigma_g[idx];
            if (!(sg >= 0.0 && sg <= 1e30)) sg = 1e-12;
            sigma_d[idx] = (EPSILON * sg > 1e-12) ? EPSILON * sg : 1e-12;

            eps_tot = 1.0 + EPSILON;
            eta = 0.5 * h * h * dlnP_dlnR;
            St = 1.0 / INVSTOKES1;
            vR_d = (2.0 / (St + eps_tot * eps_tot / St)) * eta * vK;
            vphi_d = (1.0 + (eps_tot * eta) / (eps_tot * eps_tot + St * St)) * vK;

            vx_d[idx] = vphi_d - OMEGAFRAME * R;
            vy_d[idx] = vR_d;
#ifdef ADIABATIC
            if (e_d != NULL) e_d[idx] = 1.0e-20;
#endif
        }
    }

    /* The fields were written on the host */
    OUTPUT(Fluids[1]->Density);
    OUTPUT(Fluids[1]->Vx);
    OUTPUT(Fluids[1]->Vy);
#ifdef ADIABATIC
    if (Fluids[1]->Energy != NULL) {
        OUTPUT(Fluids[1]->Energy);
    }
#endif
}
