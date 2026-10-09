//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

/* Beta cooling (BETACOOLING): the gas internal energy relaxes toward the
   initial, locally isothermal temperature on t_cool = BETA/Omega_K, as in the
   Athena++ problem generator. Applied by Sources() to the gas only. */

void Edamp_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Energy);
  INPUT(Density);
  OUTPUT(Energy);
//<\USER_DEFINED>

//<EXTERNAL>
  real* e    = Energy->field_cpu;
  real* rho  = Density->field_cpu;
  real  beta = BETA;
  real  gam  = GAMMA;
  real  r0   = R0;
  real  asp  = ASPECTRATIO;
  real  flar = FLARINGINDEX;
  real  bigg = G;
  real  mstar= MSTAR;
  int   pitch     = Pitch_cpu;
  int   stride    = Stride_cpu;
  int   size_x    = Nx;
  int   size_y    = Ny + 2*NGHY;
  int   size_z    = Nz + 2*NGHZ;
  int   fluidtype = Fluidtype;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real bigrad;
  real Hg;
  real R3;
  real omk;
  real cs;
  real cs2;
  real rho_cell;
  real e_old;
  real e_eq;
  real inv_t_cool;
  real gas_pre;
  real delta_erg;
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

  for (k = 0; k < size_z; k++) {
    for (j = 0; j < size_y; j++) {
      for (i = 0; i < size_x; i++) {
//<#>
        ll = l;

#ifdef SPHERICAL
        real r_sph = ymed(j);
        real theta = zmed(k);
        bigrad = r_sph * sin(theta);
#endif
#ifdef CYLINDRICAL
        bigrad = ymed(j);
#endif

        if (bigrad <= 0.0) {
        } else {
          Hg  = asp * pow(bigrad / r0, flar) * bigrad;
          R3  = bigrad * bigrad * bigrad;
          omk = sqrt(bigg * mstar / R3);
          cs  = omk * Hg;
          cs2 = cs * cs;

#ifdef BETACOOLING
          if (fluidtype == GAS && beta > 0.0 && cs2 > 0.0 && omk > 0.0) {
            inv_t_cool = omk / beta;
            rho_cell = rho[ll];
            if (rho_cell > 0.0) {
              e_old = e[ll];
              if (e_old <= 0.0) e_old = 1e-30;
              e_eq = rho_cell * cs2 / (gam - 1.0);
              if (e_eq <= 0.0) e_eq = 1e-30;
              /* delta_e = (P - rho c_s,0^2)/(gamma-1) * dt/t_cool */
              gas_pre = (gam - 1.0) * e_old;
              delta_erg = (gas_pre - rho_cell * cs2) / (gam - 1.0) * inv_t_cool * dt;
              e[ll] = e_old - delta_erg;
              if (!(e[ll] > 0.0) || !isfinite(e[ll])) e[ll] = e_eq;
            }
          }
#endif
        }
//<\#>
      }
    }
  }

//<\MAIN_LOOP>
}
