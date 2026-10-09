//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

#define __NOPROTO
#define __GPU

//<INCLUDES>
#include "fargo3d.h"
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
//<\INCLUDES>

/* Floors and dust velocity limiters of Model B, after the safeguards of the
   Athena++ problem generator. GPU only: a CPU build of this setup would take
   the generic src/floor.c instead. One thread per cell, ghosts included;
   called for each fluid after the transport step. */

extern int Pitch_gpu;

__device__ __forceinline__ real clamp_real(real x, real xmin, real xmax) {
  return (x < xmin) ? xmin : ((x > xmax) ? xmax : x);
}

__device__ __forceinline__ real abs_real(real x) {
  return (x < (real)0.0) ? -x : x;
}

__device__ __forceinline__ real max_real(real a, real b) {
  return (a > b) ? a : b;
}

__device__ __forceinline__ real sqrt_real(real x) {
  return (x > (real)0.0) ? sqrt(x) : (real)0.0;
}

__global__ void Apply_Floor_Kernel(real *dens,
#ifdef ADIABATIC
                                  real *energy,
#else
                                  real * /*energy*/,
#endif
                                  real *vx, real *vy, real *vz,
                                  const real *ymin, const real *zmin, const real *xmin,
                                  real *gas_dens, real *gas_vx, real *gas_vy, real *gas_vz,
#ifdef ADIABATIC
                                  real *gas_energy,
#endif
                                  int size_x, int size_y, int size_z, int pitch,
                                  int fluidtype,
                                  real gas_abs_floor, real dust_abs_floor, real energy_floor,
                                  int tv_cap_flag,
                                  real tv_cap_factor, real tv_cap_d2g_max,
                                  real tv_cap_rho_floor_factor, real tv_cap_low_gas_floor_factor,
                                  real tv_cap_vratio_override, real tv_cap_mach_override,
                                  real tv_cap_z_hg_ratio,
                                  real aspectratio, real flaringindex, real r0,
                                  real gm, real invstokes, real epsilon,
#ifdef ADIABATIC
                                  real gamma,
#endif
                                  int relmedian_cap_flag, real relmedian_cap_factor) {

  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  int total = size_x * size_y * size_z;
  if (idx >= total) return;

  int i = idx % size_x;
  int t = idx / size_x;
  int j = t % size_y;
  int k = t / size_y;
  int mem = i + pitch * (j + k * size_y);

  real y_val = ymin[j];
#ifdef Z
  real z_val = (real)0.5 * (zmin[k] + zmin[k+1]);
#else
  real z_val = (real)0.0;
#endif

  // Density floor
  real limit_floor = (fluidtype == GAS) ? gas_abs_floor : dust_abs_floor;
  real d = dens[mem];
  if (!(d >= limit_floor)) dens[mem] = limit_floor;

  // Non-finite velocities are zeroed
  real v_x = vx[mem];
  if (!(v_x >= (real)-1.0e30 && v_x <= (real)1.0e30)) v_x = (real)0.0;

  real v_y = vy[mem];
  if (!(v_y >= (real)-1.0e30 && v_y <= (real)1.0e30)) v_y = (real)0.0;

#ifdef Z
  real v_z = (real)0.0;
  if (vz != NULL) {
    v_z = vz[mem];
    if (!(v_z >= (real)-1.0e30 && v_z <= (real)1.0e30)) v_z = (real)0.0;
  }
#endif

  real R_cyl = y_val;
#ifdef Z
  real r_sph = y_val;
  real theta = z_val;
  R_cyl = r_sph * sin(z_val);
#endif
  if (R_cyl < (real)1.0e-6) R_cyl = (real)1.0e-6;

  real vK = sqrt(gm / R_cyl);
  real h  = aspectratio * pow(R_cyl / r0, flaringindex);
  real cs_ana = h * vK;

  // Gas: |v_theta| <= c_s in tenuous cells more than 3.5 H off the midplane
#ifdef Z
#ifdef SPHERICAL
  if (fluidtype == GAS && vz != NULL) {
    const real NH_CEIL       = (real)3.5;
    const real MACH_CEIL     = (real)1.0;
    const real LOW_DENS_MULT = (real)150.0;

    if (dens[mem] < LOW_DENS_MULT * gas_abs_floor) {
      real dth = abs_real(z_val - (real)(0.5 * M_PI));
      if (dth > NH_CEIL * h) {
        real vlim = MACH_CEIL * cs_ana;
        v_z = clamp_real(v_z, -vlim, vlim);
      }
    }
  }
#endif
#endif

  if (fluidtype != GAS) {
    real rhod = dens[mem];
    real rhog = gas_dens[mem];
    if (!(rhog > (real)0.0)) rhog = gas_abs_floor;

    real vg_x = gas_vx[mem];
    real vg_y = gas_vy[mem];
#ifdef Z
    real vg_z = (gas_vz != NULL) ? gas_vz[mem] : (real)0.0;
#else
    real vg_z = (real)0.0;
#endif

    // Gas stage (EPSILON ~ 0): the passive dust moves with the gas
    if (epsilon <= (real)1.0e-15) {
      vx[mem] = vg_x;
      vy[mem] = vg_y;
#ifdef Z
      if (vz != NULL) vz[mem] = vg_z;
#endif
#ifdef ADIABATIC
      if (energy != NULL) energy[mem] = (real)0.0;
#endif
      return;
    }

    real dust_to_gas = rhod / rhog;

    // Terminal-velocity cap on the dust-gas drift, as in Athena++
    if (tv_cap_flag == 1) {
      real dv_x = v_x - vg_x;
      real dv_y = v_y - vg_y;
#ifdef Z
      real dv_z = v_z - vg_z;
#else
      real dv_z = (real)0.0;
#endif

      real dust_vel_abs2 = dv_x*dv_x + dv_y*dv_y + dv_z*dv_z;
      real gas_vel_abs2 = vg_x*vg_x + vg_y*vg_y + vg_z*vg_z;
      real gas_vel_abs2_safe = (gas_vel_abs2 > (real)1.0e-24) ? gas_vel_abs2 : (real)1.0e-24;
      real gas_cs2 = cs_ana * cs_ana;

      int apply_cap = 0;
      int force_override = 0;
      int force_low_gas_cap = 0;

      // Forced when |dv|/|v_g| > VRATIO_OVERRIDE, |dv|/c_s > MACH_OVERRIDE, or
      // |dv|^2 > 20 |v_g|^2 above Z_HG_RATIO scale heights (low-gas form)
      if (dust_vel_abs2 > tv_cap_vratio_override * tv_cap_vratio_override * gas_vel_abs2_safe) {
        force_override = 1;
      }
      if (tv_cap_mach_override > (real)0.0 &&
          dust_vel_abs2 > tv_cap_mach_override * tv_cap_mach_override * gas_cs2) {
        force_override = 1;
      }
#ifdef Z
      if (dust_vel_abs2 > (real)20.0 * gas_vel_abs2_safe) {
        real z_cyl = r_sph * cos(theta);
        real omega_dyn = sqrt(gm / (R_cyl * R_cyl * R_cyl));
        real hg = cs_ana / omega_dyn;
        if (abs_real(z_cyl) > tv_cap_z_hg_ratio * hg) {
          force_override = 1;
          force_low_gas_cap = 1;
        }
      }
#endif

      // Otherwise only in dust-poor cells (rho_d and d2g both small), and
      // only where |dv|^2 > 4 |v_g|^2
      bool d2g_gated = (rhod <= tv_cap_rho_floor_factor * dust_abs_floor && dust_to_gas < tv_cap_d2g_max);
      if (!d2g_gated && force_override == 0) {
      } else {
        if (dust_vel_abs2 > (real)4.0 * gas_vel_abs2_safe) {
          apply_cap = 1;
        }
      }

      if (apply_cap) {
        bool use_low_gas_cap = (rhog < tv_cap_low_gas_floor_factor * gas_abs_floor || force_low_gas_cap == 1);

        real vrel_cap1, vrel_cap2, vrel_cap3;

#ifdef Z
        const real dlnP_dlnR_tv = (real)-2.75;
#else
        const real dlnP_dlnR_tv = (real)-1.5;
#endif
        real St  = (invstokes > (real)0.0) ? ((real)1.0 / invstokes) : (real)0.0;
        real eta = (real)0.5 * h * h * dlnP_dlnR_tv;
        real St2 = St * St;

        if (use_low_gas_cap) {
          // Tenuous gas: cap from the unperturbed drift of the disk profile
          real vrel_r = tv_cap_factor * abs_real((real)2.0 * St / ((real)1.0 + St2)) * eta * vK;
          real vrel_phi = (real)0.0;
          real vrel_theta = tv_cap_factor * abs_real(St2 / ((real)1.0 + St2)) * eta * vK;

          vrel_cap1 = vrel_r;
          vrel_cap2 = vrel_phi;
          vrel_cap3 = vrel_theta;
        } else {
#ifdef ADIABATIC
          // Otherwise from the local pressure gradient, factor*t_s*|grad P|/rho_g,
          // falling back to the profile cap where that is unreliable

          int ll_ip = mem + 1;
          int ll_im = mem - 1;
          int ll_jp = mem + pitch;
          int ll_jm = mem - pitch;
#ifdef Z
          int ll_kp = mem + pitch * size_y;
          int ll_km = mem - pitch * size_y;
#endif

          int i_plus = (i < size_x - 1) ? 1 : 0;
          int i_minus = (i > 0) ? 1 : 0;
          int j_plus = (j < size_y - 1) ? 1 : 0;
          int j_minus = (j > 0) ? 1 : 0;
#ifdef Z
          int k_plus = (k < size_z - 1) ? 1 : 0;
          int k_minus = (k > 0) ? 1 : 0;
#endif

          real P_center = (gamma - (real)1.0) * gas_energy[mem];
          real P_ip = (i_plus) ? (gamma - (real)1.0) * gas_energy[ll_ip] : P_center;
          real P_im = (i_minus) ? (gamma - (real)1.0) * gas_energy[ll_im] : P_center;
          real P_jp = (j_plus) ? (gamma - (real)1.0) * gas_energy[ll_jp] : P_center;
          real P_jm = (j_minus) ? (gamma - (real)1.0) * gas_energy[ll_jm] : P_center;
#ifdef Z
          real P_kp = (k_plus) ? (gamma - (real)1.0) * gas_energy[ll_kp] : P_center;
          real P_km = (k_minus) ? (gamma - (real)1.0) * gas_energy[ll_km] : P_center;
#endif

          real dy = ymin[j+1] - ymin[j];
          real inv_dy = (dy > (real)1e-12) ? ((real)0.5 / dy) : (real)0.0;
#ifdef Z
          real dz = zmin[k+1] - zmin[k];
          real inv_dz_theta = (dz > (real)1e-12) ? ((real)0.5 / (r_sph * dz)) : (real)0.0;
#endif
          real dx = xmin[i+1] - xmin[i];
          real inv_dx = (dx > (real)1e-12) ? ((real)0.5 / (r_sph * sin(theta) * dx)) : (real)0.0;

          real grad_P_r = (P_jp - P_jm) * inv_dy;
          real grad_P_phi = (P_ip - P_im) * inv_dx;
#ifdef Z
          real grad_P_theta = (P_kp - P_km) * inv_dz_theta;
#endif

          real St_val = (invstokes > (real)0.0) ? ((real)1.0 / invstokes) : (real)0.0;
          real omega_K_local = sqrt_real(gm / (R_cyl * R_cyl * R_cyl));
          real tau_s = (omega_K_local > (real)0.0) ? (St_val / omega_K_local) : (real)0.0;
          real vrel_coef = tv_cap_factor * abs_real(tau_s) / rhog;

          vrel_cap2 = vrel_coef * abs_real(grad_P_r);
          vrel_cap1 = vrel_coef * abs_real(grad_P_phi);
#ifdef Z
          vrel_cap3 = vrel_coef * abs_real(grad_P_theta);
#else
          vrel_cap3 = (real)0.0;
#endif

          real vrel_cap_profile = tv_cap_factor * abs_real((real)2.0 * St / ((real)1.0 + St2)) * eta * vK;
          if (vrel_cap_profile < (real)0.1 * cs_ana) {
            vrel_cap_profile = tv_cap_factor * cs_ana;
          }

          real vrel_cap_max = vrel_cap1;
          if (vrel_cap2 > vrel_cap_max) vrel_cap_max = vrel_cap2;
          if (vrel_cap3 > vrel_cap_max) vrel_cap_max = vrel_cap3;

          if (vrel_cap_max > (real)10.0 * vrel_cap_profile || vrel_cap_max < (real)1e-10) {
            vrel_cap1 = vrel_cap_profile;
            vrel_cap2 = vrel_cap_profile;
            vrel_cap3 = vrel_cap_profile;
          }
#else
          real vrel_cap = tv_cap_factor * abs_real((real)2.0 * St / ((real)1.0 + St2)) * eta * vK;
          if (vrel_cap < (real)0.1 * cs_ana) {
            vrel_cap = tv_cap_factor * cs_ana;
          }
          vrel_cap1 = vrel_cap;
          vrel_cap2 = vrel_cap;
          vrel_cap3 = vrel_cap;
#endif
        }

        // Scale the drift down, keeping its direction
        real rel_vel2_tot = dv_x*dv_x + dv_y*dv_y + dv_z*dv_z;
        real vrel_cap2_tot = vrel_cap1*vrel_cap1 + vrel_cap2*vrel_cap2 + vrel_cap3*vrel_cap3;

        if (rel_vel2_tot > vrel_cap2_tot) {
          real damping = sqrt_real(vrel_cap2_tot / rel_vel2_tot);
          v_x = vg_x + damping * dv_x;
          v_y = vg_y + damping * dv_y;
#ifdef Z
          if (vz != NULL) v_z = vg_z + damping * dv_z;
#endif
        }
      }
    }

    // Repair of dust-poor cells near the floor. Not in Athena++: without its
    // mesh refinement, the thin dust at high altitude needs it to stay stable.
    real near_floor  = (rhod < (real)30.0 * dust_abs_floor) ? (real)1.0 : (real)0.0;

#ifdef Z
    const real dlnP_dlnR = (real)-2.75;
#else
    const real dlnP_dlnR = (real)-1.5;
#endif
    real St  = (invstokes > (real)0.0) ? ((real)1.0 / invstokes) : (real)0.0;
    real eta = (real)0.5 * h * h * dlnP_dlnR;
    real dv_r_theo   = ((real)2.0 * St / ((real)1.0 + St * St)) * eta * vK;
    real dv_phi_theo = (-(St * St) / ((real)1.0 + St * St)) * eta * vK;

    real dv_x = v_x - vg_x;
    real dv_y = v_y - vg_y;
#ifdef Z
    real dv_z = v_z - vg_z;
#else
    real dv_z = (real)0.0;
#endif

    real drift_cap = (real)8.0 * cs_ana;
    real bad_drift =
      (abs_real(dv_x) > drift_cap) ||
      (abs_real(dv_y) > drift_cap) ||
#ifdef Z
      (abs_real(dv_z) > drift_cap) ||
#endif
      (!(v_x >= (real)-1.0e30 && v_x <= (real)1.0e30)) ||
      (!(v_y >= (real)-1.0e30 && v_y <= (real)1.0e30));

#ifdef Z
    if (vz != NULL) {
      bad_drift = bad_drift || (!(v_z >= (real)-1.0e30 && v_z <= (real)1.0e30));
    }
#endif

    // Near the floor and below d2g = 1e-4: drift in R and phi reset to the
    // steady value, drift in theta clamped to 2 c_s
    if ((near_floor > (real)0.5) && (dust_to_gas < (real)1.0e-4)) {
      v_y = vg_y + dv_r_theo;
      v_x = vg_x + dv_phi_theo;
#ifdef Z
      if (vz != NULL) {
        real dvz_cap = (real)2.0 * cs_ana;
        v_z = vg_z + clamp_real(dv_z, -dvz_cap, dvz_cap);
      }
#endif
    }
    // Elsewhere: drift components clamped where they exceed 8 c_s
    else if (bad_drift) {
      real cap_r   = max_real((real)6.0 * abs_real(dv_r_theo),   (real)8.0 * cs_ana);
      real cap_phi = max_real((real)6.0 * abs_real(dv_phi_theo), (real)8.0 * cs_ana);
      real cap_th  = (real)8.0 * cs_ana;

      dv_y = clamp_real(dv_y, -cap_r,   cap_r);
      dv_x = clamp_real(dv_x, -cap_phi, cap_phi);
#ifdef Z
      dv_z = clamp_real(dv_z, -cap_th,  cap_th);
#endif

      v_y = vg_y + dv_y;
      v_x = vg_x + dv_x;
#ifdef Z
      if (vz != NULL) v_z = vg_z + dv_z;
#endif
    }
  }

  // Relative-median cap: |v_d|^2 against the median of its six neighbours
  if (relmedian_cap_flag == 1 && fluidtype != GAS) {
    real vg_x = gas_vx[mem];
    real vg_y = gas_vy[mem];
#ifdef Z
    real vg_z = (gas_vz != NULL) ? gas_vz[mem] : (real)0.0;
#else
    real vg_z = (real)0.0;
#endif

    real neighbor_vel[6];
    int n_idx = 0;
    int pitch_y = pitch;
    int stride_z = pitch * size_y;

    int ll_ip = mem + 1;
    int ll_im = mem - 1;
    int ll_jp = mem + pitch_y;
    int ll_jm = mem - pitch_y;
#ifdef Z
    int ll_kp = mem + stride_z;
    int ll_km = mem - stride_z;
#endif

    if (i < size_x - 1) neighbor_vel[n_idx++] = vx[ll_ip]*vx[ll_ip] + vy[ll_ip]*vy[ll_ip] + ((vz != NULL) ? vz[ll_ip]*vz[ll_ip] : (real)0.0);
    if (i > 0) neighbor_vel[n_idx++] = vx[ll_im]*vx[ll_im] + vy[ll_im]*vy[ll_im] + ((vz != NULL) ? vz[ll_im]*vz[ll_im] : (real)0.0);
    if (j < size_y - 1) neighbor_vel[n_idx++] = vx[ll_jp]*vx[ll_jp] + vy[ll_jp]*vy[ll_jp] + ((vz != NULL) ? vz[ll_jp]*vz[ll_jp] : (real)0.0);
    if (j > 0) neighbor_vel[n_idx++] = vx[ll_jm]*vx[ll_jm] + vy[ll_jm]*vy[ll_jm] + ((vz != NULL) ? vz[ll_jm]*vz[ll_jm] : (real)0.0);
#ifdef Z
    if (k < size_z - 1 && vz != NULL) neighbor_vel[n_idx++] = vx[ll_kp]*vx[ll_kp] + vy[ll_kp]*vy[ll_kp] + vz[ll_kp]*vz[ll_kp];
    if (k > 1 && vz != NULL) neighbor_vel[n_idx++] = vx[ll_km]*vx[ll_km] + vy[ll_km]*vy[ll_km] + vz[ll_km]*vz[ll_km];
#endif

    if (n_idx >= 3) {
      for (int mi = 0; mi < n_idx-1; mi++) {
        for (int mj = mi+1; mj < n_idx; mj++) {
          if (neighbor_vel[mi] > neighbor_vel[mj]) {
            real temp = neighbor_vel[mi];
            neighbor_vel[mi] = neighbor_vel[mj];
            neighbor_vel[mj] = temp;
          }
        }
      }
      real median_vel2 = neighbor_vel[n_idx/2];

      real dust_vel2 = v_x*v_x + v_y*v_y + v_z*v_z;
      real median_cap_sq = relmedian_cap_factor * relmedian_cap_factor * median_vel2;

      if (dust_vel2 > median_cap_sq && median_vel2 > (real)1e-20) {
        real damping = sqrt_real(median_cap_sq / dust_vel2);
        v_x = vg_x + damping * (v_x - vg_x);
        v_y = vg_y + damping * (v_y - vg_y);
#ifdef Z
        if (vz != NULL) v_z = vg_z + damping * (v_z - vg_z);
#endif
      }
    }
  }

  vx[mem] = v_x;
  vy[mem] = v_y;
#ifdef Z
  if (vz != NULL) vz[mem] = v_z;
#endif

  // Gas energy floor; the pressureless dust has no energy
#ifdef ADIABATIC
  if (energy != NULL) {
    if (fluidtype == GAS) {
      real e = energy[mem];
      if (!(e >= energy_floor)) energy[mem] = energy_floor;
    } else {
      energy[mem] = (real)0.0;
    }
  }
#endif
}

extern "C" void Floor_gpu() {
  int sx = Nx + 2 * NGHX;
  int sy = Ny + 2 * NGHY;
  int sz = Nz + 2 * NGHZ;
  int total_cells = sx * sy * sz;

  real *d_dens = Density->field_gpu;
  real *d_vx   = Vx->field_gpu;
  real *d_vy   = Vy->field_gpu;
#ifdef Z
  real *d_vz   = (Vz != NULL) ? Vz->field_gpu : NULL;
  real *d_zmin = Zmin_d;
#else
  real *d_vz   = NULL;
  real *d_zmin = NULL;
#endif

#ifdef ADIABATIC
  real *d_energy = (Energy != NULL) ? Energy->field_gpu : NULL;
#else
  real *d_energy = NULL;
#endif

  real *d_gas_dens = Fluids[0]->Density->field_gpu;
  real *d_gas_vx   = Fluids[0]->Vx->field_gpu;
  real *d_gas_vy   = Fluids[0]->Vy->field_gpu;
#ifdef Z
  real *d_gas_vz   = (Fluids[0]->Vz != NULL) ? Fluids[0]->Vz->field_gpu : NULL;
#else
  real *d_gas_vz   = NULL;
#endif
#ifdef ADIABATIC
  real *d_gas_energy = Fluids[0]->Energy->field_gpu;
#else
  real *d_gas_energy = NULL;
#endif

  const real gas_abs_floor = (real)GASFLOOR;
  const real dust_abs_floor = (real)DUSTFLOOR;
  const real energy_floor  = (real)ENERGYFLOOR;
  const int  tv_cap_flag = (int)DUST_TVCAP_FLAG;
  const real tv_cap_factor = (real)DUST_TVCAP_FACTOR;
  const real tv_cap_d2g_max = (real)DUST_TVCAP_D2G_MAX;
  const real tv_cap_rho_floor_factor = (real)DUST_TVCAP_RHO_FLOOR_FACTOR;
  const real tv_cap_low_gas_floor_factor = (real)DUST_TVCAP_LOW_GAS_FLOOR_FACTOR;
  const real tv_cap_vratio_override = (real)DUST_TVCAP_VRATIO_OVERRIDE;
  const real tv_cap_mach_override = (real)DUST_TVCAP_MACH_OVERRIDE;
  const real tv_cap_z_hg_ratio = (real)DUST_TVCAP_Z_HG_RATIO;

  const real aspectratio  = ASPECTRATIO;
  const real flaringindex = FLARINGINDEX;
  const real r0           = R0;
  const real gm           = G * MSTAR;
  const real invstokes    = INVSTOKES1;
  const real epsilon      = EPSILON;

  int blockSize = 256;
  int gridSize  = (total_cells + blockSize - 1) / blockSize;

  Apply_Floor_Kernel<<<gridSize, blockSize>>>(
      d_dens,
      d_energy,
      d_vx, d_vy, d_vz,
      Ymin_d, d_zmin, Xmin_d,
      d_gas_dens, d_gas_vx, d_gas_vy, d_gas_vz,
      d_gas_energy,
      sx, sy, sz, Pitch_gpu,
      Fluidtype,
      gas_abs_floor, dust_abs_floor, energy_floor,
      tv_cap_flag, tv_cap_factor, tv_cap_d2g_max,
      tv_cap_rho_floor_factor, tv_cap_low_gas_floor_factor,
      tv_cap_vratio_override, tv_cap_mach_override, tv_cap_z_hg_ratio,
      aspectratio, flaringindex, r0,
      gm, invstokes, epsilon,
#ifdef ADIABATIC
      (real)GAMMA,
#endif
      (int)RELATIVE_MEDIAN_CAP_FLAG, (real)RELATIVE_MEDIAN_CAP_FACTOR
  );

  cudaError_t err = cudaPeekAtLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "Floor launch failed on rank %d: %s\n",
            CPU_Rank, cudaGetErrorString(err));
    fflush(stderr);
    exit(1);
  }

  err = cudaDeviceSynchronize();
  if (err != cudaSuccess) {
    fprintf(stderr, "Floor sync failed on rank %d: %s\n",
            CPU_Rank, cudaGetErrorString(err));
    fflush(stderr);
    exit(1);
  }
}
