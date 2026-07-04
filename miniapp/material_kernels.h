/**
 * @file material_kernels.h
 *
 * @brief Declares `mat_model_kernel`, the miniapp's thin per-time-step call into the
 * ExaCMech material model.
 *
 * This is the middle stage of the miniapp's three-stage per-time-step pipeline (see
 * `orientation_evolution.cxx`'s main loop): `setup_data()` (in `setup_kernels.h`)
 * populates the arrays below from the persistent `state_vars`/`cauchy_stress_array`
 * buffers, `mat_model_kernel()` advances them one time step, and `retrieve_data()` (in
 * `retrieve_kernels.h`) writes the results back into the persistent buffers.
 */

#pragma once

namespace ecmech
{
    class matModelBase;
}

/**
 * @brief Advance all quadrature points one time step by calling into the material
 * model's `getResponseECM`.
 *
 * All parallelization (CPU/OpenMP/GPU dispatch across the `nqpts` points) happens inside
 * `matModelBase::getResponseECM` itself, based on the execution strategy the model was
 * configured with (`matModelBase::setExecutionStrategy`) during setup in
 * `orientation_evolution.cxx`; this function is just a direct passthrough.
 *
 * @param[in] mat_model_base Fully configured material model (see `ecmech::makeMatModel`
 * and `matModelBase::initFromParams`/`complete` in `orientation_evolution.cxx`).
 * @param[in] nqpts Number of quadrature points (independent material points) to advance.
 * @param[in] dt Time-step size.
 * @param[in,out] state_vars_array History (state) variables, `nqpts` × the model's
 * history stride (see `setup_kernels.h`'s `init_data` for this miniapp's specific
 * history-array layout).
 * @param[in,out] cauchy_stress_d6p_array Cauchy stress (deviatoric 6-vector + pressure),
 * `nqpts` × `ecmech::nsvp`; populated by `setup_data`, consumed/updated here.
 * @param[in] def_rate_d6v_array Deformation rate (deviatoric 6-vector + volumetric
 * rate), `nqpts` × `ecmech::nsvp`; populated by `setup_data`.
 * @param[in] spin_vec_array Spin (vorticity) axial vector, `nqpts` × `ecmech::nwvec`;
 * populated by `setup_data`.
 * @param[out] ddsdde_array Tangent stiffness matrix, `nqpts` × `ecmech::nsvec` ×
 * `ecmech::nsvec`; not consumed elsewhere in the miniapp, but computed since
 * `getResponseECM` always requires a (possibly unused) tangent output buffer.
 * @param[in,out] rel_vol_ratios_array Relative-volume bookkeeping, `nqpts` ×
 * `ecmech::nvr`; populated by `setup_data`.
 * @param[in,out] internal_energy_array Internal energy, `nqpts` × `ecmech::ne`;
 * populated by `setup_data`.
 * @param[in,out] tkelv_array Temperature in Kelvin, `nqpts` values; populated by
 * `setup_data` (held fixed at 300 K in this miniapp).
 * @param[out] sdd_array Auxiliary derived quantities (e.g. shear modulus), `nqpts` ×
 * `ecmech::nsdd`.
 */
void mat_model_kernel(const ecmech::matModelBase* mat_model_base,
                      const int nqpts, const double dt, double* state_vars_array,
                      double* cauchy_stress_d6p_array, double* def_rate_d6v_array,
                      double* spin_vec_array, double* ddsdde_array,
                      double* rel_vol_ratios_array, double* internal_energy_array,
                      double* tkelv_array, double* sdd_array);

