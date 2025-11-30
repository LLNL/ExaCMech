/**
 * @file vb_d_vars_set.h
 * @brief Deviatoric 6-vector component declarations with explicit storage.
 * 
 * This header provides named variable declarations for 5 components extracted from a
 * 6-component deviatoric vector representation. It is designed to be #included within
 * functions that work with deviatoric tensors stored in an extended 6-component format,
 * extracting only the first 5 independent components.
 * 
 * Usage pattern:
 * Functions that receive deviatoric tensors in 6-vector format (with potential redundancy
 * or additional data in the 6th component) include this header to extract the meaningful
 * 5 components. The header expects a pointer `cmv6b` to be defined in the including scope.
 * 
 * Deviatoric tensor representations:
 * ECMech uses multiple representations for deviatoric (traceless) symmetric tensors:
 * 
 * 1. **Deviatoric 5-vector** (minimal representation):
 *    - 5 independent components of traceless symmetric tensor
 *    - Used in vad_vars_set.h and vadl_vars_set.h
 * 
 * 2. **Deviatoric 6-vector** (extended representation):
 *    - 6 components where first 5 are deviatoric
 *    - 6th component may contain:
 *      * Zero (explicit traceless enforcement)
 *      * Volumetric/pressure information
 *      * Numerical tolerance tracking
 *    - Used in some intermediate calculations
 * 
 * 3. **Full 6-vector** (Voigt notation):
 *    - All 6 independent components of general symmetric tensor
 *    - Not necessarily traceless
 * 
 * This header extracts the first 5 components from format (2).
 * 
 * Component naming:
 * - vb1, vb2, vb3, vb4, vb5: Five deviatoric tensor components
 * - "vb" prefix: Vector "b" variant (distinguishes from "va" in vad_vars_set.h)
 * - Source array: cmv6b (6-component array)
 * - Only indices 0-4 extracted (6th component cmv6b[5] not used here)
 * 
 * Array name interpretation:
 * `cmv6b` likely refers to:
 * - "cm": Crystal or continuum mechanics
 * - "v6": 6-component vector
 * - "b": Variant identifier (distinguishing from other vector types)
 * 
 * Applications:
 * This header may be included in functions that:
 * - Convert between 5-component and 6-component deviatoric representations
 * - Process stress/strain data from external sources using 6-vector format
 * - Handle intermediate results in tensor rotation operations
 * - Interface with codes using different deviatoric representations
 * 
 * Comparison with related headers:
 * 
 * | Header          | Source array | Components  | Frame   | Use case                   |
 * |-----------------|--------------|-------------|---------|----------------------------|
 * | vad_vars_set.h  | vec_sm[5]    | 5           | Sample  | Sample frame deviatoric    |
 * | vadl_vars_set.h | vec_lat[5]   | 5           | Lattice | Lattice frame deviatoric   |
 * | vb_d_vars_set.h | cmv6b[6]     | 5 (first)   | Either  | Extended format extraction |
 * | vw_vars_set.h   | cmv3w[3]     | 3           | Either  | Spin/vorticity             |
 * 
 * Design rationale:
 * Provides consistent interface for accessing deviatoric components regardless of
 * whether source data is in minimal 5-component or extended 6-component format.
 * Named variables improve code readability when working with mathematical expressions
 * involving deviatoric tensor components.
 * 
 * @note This file contains only variable declarations. It relies on `cmv6b` being
 *       defined in the including scope as a pointer to (at least) a 6-component array.
 * @note Only the first 5 components are extracted; cmv6b[5] is ignored.
 * @note The relationship between `cmv6b` array and physical tensor components depends
 *       on the specific deviatoric representation convention used in the including code.
 * 
 * @ingroup ECMech_utilities
 * 
 * @see vad_vars_set.h for 5-component sample frame deviatoric
 * @see vadl_vars_set.h for 5-component lattice frame deviatoric
 * @see mc_vars_set.h for rotation matrix components
 */

double vb1 = cmv6b[0];
double vb2 = cmv6b[1];
double vb3 = cmv6b[2];
double vb4 = cmv6b[3];
double vb5 = cmv6b[4];
