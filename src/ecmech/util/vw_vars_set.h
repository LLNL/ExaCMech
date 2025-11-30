/**
 * @file vw_vars_set.h
 * @brief Spin vector component declarations for skew-symmetric tensors.
 * 
 * This header provides named variable declarations for the 3 independent components of a
 * skew-symmetric (anti-symmetric) second-order tensor. It is designed to be #included
 * within functions that perform tensor operations on spin/vorticity quantities, replacing
 * indexed array accesses with meaningful component names.
 * 
 * Usage pattern:
 * Functions in ECMech_util.h that compute rotation operations or derivatives include this
 * header to extract spin tensor components. The header expects a pointer `cmv3w` to be
 * defined in the including scope, pointing to a 3-component spin (vorticity) vector.
 * 
 * Skew-symmetric tensor representation:
 * A skew-symmetric 3×3 tensor W (with W^T = -W and W_ii = 0) has only 3 independent
 * components, corresponding to the 3 degrees of freedom of a rotation. These are stored
 * as an axial vector:
 * 
 * Full tensor:
 *     W = |  0    -w₃    w₂ |
 *         |  w₃    0    -w₁ |
 *         | -w₂    w₁    0  |
 * 
 * Axial vector representation:
 *     cmv3w = [vw1, vw2, vw3]
 * 
 * Component mapping:
 * - vw1 corresponds to W₃₂ = -W₂₃ (rotation about x₁-axis)
 * - vw2 corresponds to W₁₃ = -W₃₁ (rotation about x₂-axis)
 * - vw3 corresponds to W₂₁ = -W₁₂ (rotation about x₃-axis)
 * 
 * Component naming:
 * - vw1, vw2, vw3: Three components of axial spin vector
 * - "vw" prefix: Vector for spin/vorticity ("w" for ω, vorticity)
 * - Index 1-3: The three rotational degrees of freedom
 * 
 * Physical interpretation:
 * The spin vector represents the rotational velocity of material elements:
 * - Magnitude: Angular velocity of rotation
 * - Direction: Axis of rotation (right-hand rule)
 * - Relation to vorticity: W = curl(v)/2 for velocity field v
 * 
 * Applications in crystal plasticity:
 * This header is included in functions computing:
 * - d_rot_mat_wveccp_latop(): Derivatives for spin tensor rotations
 * - d_rot_mat_wveccp_smop(): Derivatives for spin tensor transformations
 * - Plastic spin W^p = Σ γ̇^α Q^α from slip system contributions
 * - Lattice spin (rigid body rotation of crystal axes)
 * - Total spin decomposition: W = W^e + W^p
 * 
 * Physical quantities represented:
 * Common spin tensors include:
 * - Continuum spin: W = (∇v - (∇v)^T) / 2
 * - Plastic spin: W^p from asymmetric crystal slip
 * - Elastic spin: W^e from lattice rotation
 * - Corotational derivatives for objective stress rates
 * 
 * Spin vs. strain rate:
 * - Spin (skew W): Represents rotation without stretching (3 DOF)
 * - Strain rate (symmetric D): Represents stretching without rotation (6 DOF, 5 deviatoric)
 * - Velocity gradient: L = D + W (9 total DOF)
 * 
 * Frame transformations:
 * Unlike deviatoric tensors (5 components), spin vectors have special transformation
 * properties due to their pseudo-vector nature. Under rotation R:
 *   w_new = R w_old (transforms as regular vector)
 *   W_new = R W_old R^T (transforms as tensor)
 * 
 * Design rationale:
 * Named components (vw1-vw3) improve readability in rotation derivative calculations
 * and make the physical meaning of spin contributions more apparent in the code.
 * 
 * @note This file contains only variable declarations. It relies on `cmv3w` being
 *       defined in the including scope as a pointer to a 3-component array.
 * @note The array name `cmv3w` may refer to "crystal material vector 3-component w"
 *       or similar internal naming convention.
 * @note Components represent axial vector form, not full skew-symmetric tensor.
 * 
 * @ingroup ECMech_utilities
 * 
 * @see ECMech_util.h for functions using this header
 * @see vad_vars_set.h for deviatoric strain rate components
 * @see vadl_vars_set.h for deviatoric components in lattice frame
 * @see mc_vars_set.h for rotation matrix components
 */

double vw1 = cmv3w[0];
double vw2 = cmv3w[1];
double vw3 = cmv3w[2];
