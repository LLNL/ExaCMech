/**
 * @file vad_vars_set.h
 * @brief Deviatoric vector component declarations for sample frame tensors.
 * 
 * This header provides named variable declarations for the 5 independent components of a
 * symmetric deviatoric (traceless) second-order tensor in the sample (global) coordinate
 * frame. It is designed to be #included within functions that perform tensor operations,
 * replacing indexed array accesses with meaningful component names.
 * 
 * Usage pattern:
 * Functions in ECMech_util.h that compute rotation operations or derivatives include this
 * header to extract deviatoric tensor components. The header expects a pointer `vec_sm`
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
 *     vec_sm = [va1, va2, va3, va4, va5]
 * where each component relates to the full tensor through specific linear combinations
 * involving factors of √2 and √3 to preserve tensor norms under the mapping.
 * 
 * Component naming:
 * - va1, va2: Combinations of diagonal components (with √3 factors for normalization)
 * - va3, va4, va5: Off-diagonal shear components (with √2 factors)
 * - "va" prefix: Vector in sample ("a" for "ambient" or "applied") frame
 * - Index 1-5: The five degrees of freedom of a deviatoric tensor
 * 
 * Sample frame:
 * The "sample frame" (also called global, laboratory, or spatial frame) is the fixed
 * external coordinate system in which:
 * - Boundary conditions are applied
 * - Deformation gradients are measured
 * - Stress components are reported
 * 
 * This contrasts with the lattice/crystal frame which rotates with the material.
 * 
 * Applications:
 * This header is included in functions computing:
 * - d_rot_mat_vecd_smop(): Derivatives for sample→crystal tensor rotation
 * - d_rot_mat_vecd_latop(): Derivatives for crystal→sample tensor rotation  
 * - Frame transformations for deformation rate tensors
 * - Frame transformations for stress tensors
 * - Jacobian contributions for implicit Newton solves
 * 
 * Physical quantities represented:
 * Common deviatoric tensors in sample frame include:
 * - Deformation rate D (symmetric part of velocity gradient)
 * - Deviatoric stress σ' = σ - (1/3)tr(σ)I
 * - Elastic strain deviator ε'_e
 * 
 * Design rationale:
 * Named components (va1-va5) are more readable than vec_sm[0]-vec_sm[4] in the
 * complex algebraic expressions for rotation derivatives and tensor transformations.
 * 
 * @note This file contains only variable declarations. It relies on `vec_sm` being
 *       defined in the including scope as a pointer to a 5-component array.
 * @note Components are in the specific deviatoric representation used throughout ECMech,
 *       not standard Voigt notation.
 * 
 * @ingroup ECMech_utilities
 * 
 * @see ECMech_util.h for functions using this header
 * @see vadl_vars_set.h for deviatoric components in lattice frame
 * @see mc_vars_set.h for rotation matrix components
 * @see vw_vars_set.h for spin vector components
 */

double va1 = vec_sm[0];
double va2 = vec_sm[1];
double va3 = vec_sm[2];
double va4 = vec_sm[3];
double va5 = vec_sm[4];
