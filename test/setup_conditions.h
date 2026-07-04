/**
 * @file setup_conditions.h
 *
 * @brief Fragment (see `setup_base.h` for the general `setup_*.h` inclusion pattern)
 * declaring the shared kinematic "boundary condition" test inputs -- deformation rate,
 * spin, time step, and initial relative-volume bookkeeping -- reused across the
 * `evptn`/`updst` test cases. Unlike most of the other `setup_*.h` fragments, this one
 * doesn't branch on `STACK_PARAMS`: it declares plain local values, not model
 * parameters.
 */

/**
 * @brief Deformation rate in ExaCMech's deviatoric-6 + volumetric-rate form
 * (`ecmech::nsvp` wide; see `matModelBase::getResponseECM`'s `def_rate_d6vV` doc).
 * The three normal components (`-0.5, -0.5, 1.0`) already sum to zero (a traceless,
 * z-direction-tension-like deformation) and the volumetric-rate slot is `0.0`, so this
 * is a purely deviatoric test case; it is then rescaled below.
 */
double def_rate_d6v_sample[ecmech::nsvp] = { -0.5, -0.5, 1.0,
                                      0.0, 0.0, 0.0,
                                      0.0 };
// Rescale the whole 7-vector by sqrt(2/3) so that, once converted to the 5-component
// deviatoric representation below, its Euclidean norm comes out to exactly 1 -- a
// convenient unit-magnitude reference deformation rate for the tests that use it.
vecsVsa<ecmech::nsvp>(def_rate_d6v_sample, ecmech::sqr2b3); // scale so that def_rate_d5_sample comes out to unit magnitude
//
/** @brief `def_rate_d6v_sample` converted to the 5-component deviatoric (`ntvec`-wide) representation via `svecToVecd` (see `ECMech_util.h`); has unit magnitude after the rescale above. */
double def_rate_d5_sample[ecmech::ntvec];
svecToVecd(def_rate_d5_sample, def_rate_d6v_sample);

/** @brief Time-step size used by the tests that include this fragment. */
double dt = 1e-1;

/** @brief Spin (vorticity) axial vector -- a nonzero rigid-body rotation rate about z, so the tests also exercise lattice-rotation update, not just pure stretch. */
double spin_vec_sample[ecmech::nwvec] = { 0.0, 0.0, 0.5 };

/** @brief Relative-volume bookkeeping `[rel_vol_n, rel_vol_{n+1}, rate, delta]` (see `ECMech_const.h`'s `nvr` doc), initialized to the undeformed reference state (both volumes 1.0, zero rate/delta). */
double rel_vol_ratios[ecmech::nvr] = { 1.0, 1.0, 0.0, 0.0 };
