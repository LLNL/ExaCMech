/**
 * @file vadl_vars_set.h
 * @brief Deviatoric vector component declarations for lattice frame tensors.
 * 
 * This header provides named variable declarations for the 5 independent components of a
 * symmetric deviatoric (traceless) second-order tensor in the lattice (crystal) coordinate
 * frame. It is designed to be #included within functions that perform tensor operations,
 * replacing indexed array accesses with meaningful component names.
 * 
 * Usage pattern:
 * Functions in ECMech_util.h that compute rotation operations or derivatives include this
 * header to extract deviatoric tensor components. The header expects a pointer `vec_lat`
 * to be defined in the including scope, pointing to a 5-component deviatoric vector.
 * 
 * Deviatoric tensor representation:
 * A symmetric 3×3 traceless tensor T (with T_11 + T_22 + T_33 = 0) has only 5 independent
 * components. These are stored in a compressed "deviatoric 5-vector" representation:
 * 
 * Full tensor:
 *     T = | T_11   T_12   T_13 |
 *         | T_12   T_22   T_23 |  with T_11 + T_22 + T_33 = 0
 *         | T_13   T_23   T_33 |
 * 
 * Deviatoric 5-vector components:
 *     vec_lat = [va1, va2, va3, va4, va5]
 * where each component relates to the full tensor through specific linear combinations
 * involving factors of √2 and √3 to preserve tensor norms under the mapping.
 * 
 * Component naming:
 * - va1, va2: Combinations of diagonal components (with √3 factors for normalization)
 * - va3, va4, va5: Off-diagonal shear components (with √2 factors)
 * - "va" prefix: Vector (reused naming from vad_vars_set.h for consistency)
 * - Source array: vec_lat (lattice frame, not vec_sm sample frame)
 * 
 * Lattice frame:
 * The "lattice frame" (also called crystal or material frame) is the coordinate system
 * attached to the crystal lattice that rotates with the material during deformation:
 * - x₁, x₂, x₃ aligned with crystal axes
 * - Rotates relative to sample frame via orientation quaternion
 * - Elastic constitutive law applies directly in this frame
 * - Slip system geometry defined in this frame
 * 
 * This contrasts with the sample/global frame which is fixed in space.
 * 
 * Applications:
 * This header is included in functions computing:
 * - d_rot_mat_vecd_smop(): Derivatives for sample→crystal tensor rotation
 * - d_rot_mat_vecd_latop(): Derivatives for crystal→sample tensor rotation
 * - Frame transformations for stress tensors
 * - Frame transformations for strain rate tensors
 * - Jacobian contributions for crystal plasticity updates
 * 
 * Physical quantities represented:
 * Common deviatoric tensors in lattice frame include:
 * - Elastic strain deviator ε'_e (used in elastic law)
 * - Plastic deformation rate D^p = Σ γ̇^α P^α
 * - Kirchhoff stress τ' in crystal frame
 * - Schmid tensor projections P^α
 * 
 * Crystal plasticity workflow:
 * 1. Sample frame deformation rate D_sample → rotate to D_lattice
 * 2. Decompose: D_lattice = D^e_lattice + D^p_lattice
 * 3. Compute stress in lattice frame from elastic strain
 * 4. Rotate stress back to sample frame for equilibrium
 * 
 * Design rationale:
 * Named components (va1-va5) improve readability of rotation derivative expressions.
 * Separate file from vad_vars_set.h makes clear distinction between sample and
 * lattice frame quantities, preventing frame confusion errors.
 * 
 * @note This file contains only variable declarations. It relies on `vec_lat` being
 *       defined in the including scope as a pointer to a 5-component array.
 * @note Uses same component names (va1-va5) as vad_vars_set.h but extracts from
 *       different source array (vec_lat vs vec_sm).
 * @note Components are in ECMech's deviatoric representation, not standard Voigt notation.
 * 
 * @ingroup ECMech_utilities
 * 
 * @see ECMech_util.h for functions using this header
 * @see vad_vars_set.h for deviatoric components in sample frame
 * @see mc_vars_set.h for rotation matrix components
 * @see ECMech_elastic.h for elastic law in lattice frame
 */

double va1 = vec_lat[0];
double va2 = vec_lat[1];
double va3 = vec_lat[2];
double va4 = vec_lat[3];
double va5 = vec_lat[4];
