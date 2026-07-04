/**
 * @file setup_kernels.h
 *
 * @brief Declares the miniapp's setup-stage kernels: one-time state initialization
 * (`init_data`, `setup_velocity_grad`) and the per-time-step data prep that runs before
 * every call into the material model (`setup_data`).
 *
 * See `orientation_evolution.cxx` for how these are sequenced: `init_data` and
 * `setup_velocity_grad` run once before the time-step loop; `setup_data` runs once per
 * step, immediately before `mat_model_kernel` (`material_kernels.h`).
 */

#pragma once

#include <vector>

namespace ecmech
{
    class matModelBase;
}

/**
 * @brief One-time initialization of the persistent `state_vars` array for every
 * quadrature point, from the model's initial history values and the per-point crystal
 * orientations.
 *
 * The `state_vars` array this miniapp carries between time steps is laid out as the
 * model's own evptn history block (effective shear rate, effective shear, flow
 * strength, solver function-eval count, elastic strain, orientation quaternion,
 * hardening state, slip rates -- see `ecmech::evptn`'s `iHistA_*`/`iHistLb*` index
 * constants in `ECMech_base_classes.h`) immediately followed by two miniapp-specific
 * bookkeeping slots this function does *not* initialize from the model (relative volume
 * and internal energy; see `retrieve_data` in `retrieve_kernels.cxx`, which is what
 * updates them every step). Every other history slot is copied from the model's own
 * `getHistInfo()` initial values, except the orientation quaternion (`ori_vec`, unique
 * per point) and the relative volume (forced to `1.0`, undeformed reference state).
 *
 * @param[in] ori_vec Flat array of initial orientation quaternions, `nqpts` ×
 * `ecmech::qdim`.
 * @param[in] mat_model_base Fully configured material model; used only to query its
 * initial history values via `getHistInfo()`.
 * @param[in] nqpts Number of quadrature points.
 * @param[in] num_hardness Number of hardening state variables for this model
 * (`index_map["num_hardening"]` in `orientation_evolution.cxx`).
 * @param[in] num_slip Number of slip systems for this model
 * (`index_map["num_slip_system"]`).
 * @param[in] ind_gdot Starting index of the slip-rate block within the model's own
 * history layout (`index_map["index_slip_rates"]`, i.e. `iHistLbGdot`).
 * @param[in] state_var_vdim Per-point stride (width) of `state_vars`, i.e. the model's
 * history stride plus the two appended volume/energy slots.
 * @param[out] state_vars Persistent history array to initialize, `nqpts` ×
 * `state_var_vdim`.
 */
void init_data(const std::vector<double>& ori_vec, const ecmech::matModelBase* mat_model_base,
               const int nqpts, const int num_hardness,
               const int num_slip, const int ind_gdot,
               const int state_var_vdim, double* state_vars);

/**
 * @brief Broadcast a single 3x3 macroscopic velocity gradient to every quadrature point.
 *
 * This miniapp models `nqpts` independent single-crystal orientations, all subjected to
 * the identical macroscopic velocity gradient (a Taylor-type/full-constraint
 * polycrystal averaging setup) -- so every point simply gets a copy of the same
 * `velocity_grad_input` matrix. More interesting per-point deformation histories (with
 * spin, or point-to-point variation) could be built the same way just by varying what's
 * written per point.
 *
 * @param[in] velocity_grad_input 3x3 velocity gradient, row-major, flattened to 9
 * values.
 * @param[out] velocity_grad Output buffer, `nqpts` copies of the 3x3 matrix, stored in
 * the column-major-permuted RAJA layout `setup_data` (`setup_kernels.cxx`) expects.
 * @param[in] nqpts Number of quadrature points.
 */
void setup_velocity_grad(const std::vector<double>& velocity_grad_input, double* const velocity_grad, const int nqpts);

/**
 * @brief Populate one time step's worth of `getResponseECM` input arrays from the
 * persistent `state_vars`/`cauchy_stress_array`/velocity-gradient buffers.
 *
 * For every quadrature point, this: fixes temperature at 300 K; zeros the tangent
 * stiffness output buffer; copies internal energy from `state_vars`; extracts the spin
 * (skew part of the velocity gradient, as an axial vector) and the deviatoric +
 * volumetric-rate deformation rate (`ecmech::nsvp`-wide, matching
 * `matModelBase::getResponseECM`'s `def_rate_d6vV` convention: indices 0-2/3-5 are the
 * symmetrized, trace-removed normal/shear rates, index `ecmech::iSvecP` is `trace(L)`,
 * the volumetric rate) from the (symmetric-part-only) velocity gradient; integrates the
 * previous step's relative volume forward by that volumetric rate to get this step's
 * `rel_vol_ratios` bookkeeping (`[rel_vol_n, rel_vol_{n+1}, rate, delta]`, see
 * `ECMech_const.h`'s `nvr` doc); and converts the previous step's Cauchy stress from
 * plain Voigt form to ExaCMech's deviatoric + pressure (`svecp`) form (see
 * `ECMech_util.h`), the inverse of what `retrieve_data` (`retrieve_kernels.cxx`) does
 * going back out.
 *
 * @param[in] nqpts Number of quadrature points.
 * @param[in] nstatev Per-point stride (width) of `state_vars_array` (see `init_data`).
 * @param[in] dt Time-step size.
 * @param[in] vel_grad_array Macroscopic velocity gradient per point (see
 * `setup_velocity_grad`), `nqpts` × 3 × 3.
 * @param[in] cauchy_stress_array Beginning-of-step Cauchy stress in plain Voigt form,
 * `nqpts` × `ecmech::nsvec`.
 * @param[in] state_vars_array Persistent history array, `nqpts` × `nstatev` (only the
 * internal-energy and relative-volume slots are read here).
 * @param[out] cauchy_stress_d6p_array Beginning-of-step Cauchy stress in deviatoric +
 * pressure form, `nqpts` × `ecmech::nsvp`.
 * @param[out] def_rate_d6v_array Deformation rate in deviatoric + volumetric-rate form,
 * `nqpts` × `ecmech::nsvp`.
 * @param[out] spin_vec_array Spin (vorticity) axial vector, `nqpts` × `ecmech::nwvec`.
 * @param[out] ddsdde_array Tangent stiffness matrix, zeroed, `nqpts` × `ecmech::nsvec` ×
 * `ecmech::nsvec`.
 * @param[out] rel_vol_ratios_array Relative-volume bookkeeping, `nqpts` × `ecmech::nvr`.
 * @param[out] internal_energy_array Beginning-of-step internal energy, `nqpts` ×
 * `ecmech::ne`.
 * @param[out] tkelv_array Temperature in Kelvin, `nqpts` values (fixed at 300 K).
 */
void setup_data(const int nqpts, const int nstatev,
                const double dt, const double* vel_grad_array,
                const double* cauchy_stress_array, const double* state_vars_array,
                double* cauchy_stress_d6p_array, double* def_rate_d6v_array,
                double* spin_vec_array, double* ddsdde_array,
                double* rel_vol_ratios_array, double* internal_energy_array,
                double* tkelv_array);
