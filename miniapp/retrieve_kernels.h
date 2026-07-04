/**
 * @file retrieve_kernels.h
 *
 * @brief Declares `retrieve_data`, the miniapp's final per-time-step stage: copies the
 * material model's end-of-step outputs back into the miniapp's persistent
 * `state_vars`/`cauchy_stress` buffers, converting representations where needed.
 *
 * Runs after `setup_data()` (`setup_kernels.h`) and `mat_model_kernel()`
 * (`material_kernels.h`) in the main loop in `orientation_evolution.cxx`.
 */

#pragma once

namespace ecmech
{
    class matModelBase;
}

/**
 * @brief Write the material model's end-of-step relative volume, internal energy, and
 * stress back into the miniapp's persistent per-point arrays.
 *
 * Two conversions happen here:
 * -# The updated relative volume (`rel_vol_ratios_array[1]`, i.e. `rel_vol_{n+1}`) and
 *    internal energy are copied into the two miniapp-specific bookkeeping slots that
 *    `init_data` (`setup_kernels.cxx`) appends after the model's own history layout (see
 *    that function's doc for the full `state_vars` layout).
 * -# The Cauchy stress is converted from ExaCMech's deviatoric-6 + pressure
 *    representation (`cauchy_stress_d6p_array`, `ecmech::nsvp` wide, with the pressure
 *    term at `ecmech::iSvecP`) back to a plain Voigt 6-vector (`cauchy_stress_array`) by
 *    adding the (negated) pressure back onto the three normal components -- the inverse
 *    of the conversion `setup_data` (`setup_kernels.cxx`) performs going in.
 *
 * @param[in] nqpts Number of quadrature points (independent material points).
 * @param[in] nstatev Width of `state_vars_array` per point (the model's own history
 * stride plus the miniapp's appended volume/energy slots -- see `init_data`).
 * @param[in] cauchy_stress_d6p_array End-of-step Cauchy stress (deviatoric 6 +
 * pressure), `nqpts` × `ecmech::nsvp`, as written by `mat_model_kernel`.
 * @param[in] rel_vol_ratios_array End-of-step relative-volume bookkeeping, `nqpts` ×
 * `ecmech::nvr`, as written by `mat_model_kernel`.
 * @param[in] internal_energy_array End-of-step internal energy, `nqpts` × `ecmech::ne`,
 * as written by `mat_model_kernel`.
 * @param[in,out] state_vars_array Persistent history array, `nqpts` × `nstatev`; only
 * the trailing volume-ratio and internal-energy slots are updated here (the rest of the
 * model's own history state is written directly by `mat_model_kernel`/`getResponseECM`).
 * @param[out] cauchy_stress_array Persistent Cauchy stress in plain Voigt form, `nqpts`
 * × `ecmech::nsvec`.
 */
void retrieve_data(const int nqpts, const int nstatev,
                   const double* cauchy_stress_d6p_array, const double* rel_vol_ratios_array,
                   const double* internal_energy_array, double* state_vars_array,
                   double* cauchy_stress_array);

