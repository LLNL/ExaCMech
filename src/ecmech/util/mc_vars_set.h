/**
 * @file mc_vars_set.h
 * @brief Matrix component variable declarations for rotation tensor operations.
 * 
 * This header provides named variable declarations for the 9 components of a 3×3 rotation
 * matrix C. It is designed to be #included within functions that perform tensor rotation
 * operations, improving code readability by replacing indexed array accesses with meaningful
 * component names.
 * 
 * Usage pattern:
 * Functions in ECMech_util.h that compute derivatives of rotation operations include this
 * header to extract rotation matrix components. The header expects a pointer `c` to be
 * defined in the including scope, pointing to a 3×3 matrix stored in row-major order.
 * 
 * Matrix structure:
 * The rotation matrix C is a 3×3 proper orthogonal matrix (C^T C = I, det(C) = 1) that
 * transforms vectors between coordinate frames:
 * 
 *     C = | c11  c12  c13 |
 *         | c21  c22  c23 |
 *         | c31  c32  c33 |
 * 
 * Component naming convention:
 * - c_ij: Component in row i, column j (1-indexed in notation, 0-indexed in array)
 * - First index: row (1=x, 2=y, 3=z in crystal frame)
 * - Second index: column (1=x, 2=y, 3=z in sample frame)
 * 
 * Array indexing:
 * Uses ECMECH_NN_INDX(row, col, ndim) macro to compute linear array index for
 * the (row, col) element of a ndim×ndim matrix stored in row-major order.
 * 
 * Coordinate transformations:
 * In crystal plasticity context, C often represents:
 * - Crystal-to-sample frame rotation
 * - Lattice orientation transformation
 * - Derived from quaternion representation of orientation
 * 
 * Applications:
 * This header is included in functions computing:
 * - d_rot_mat_vecd_smop(): Derivative of sample→crystal rotation for deviatoric tensors
 * - d_rot_mat_vecd_latop(): Derivative of crystal→sample rotation for deviatoric tensors
 * - rot_mat_symm(): 5×5 rotation matrix for symmetric deviatoric tensors
 * - Various Jacobian computations for implicit integration
 * 
 * Design rationale:
 * Using named variables (c11, c12, etc.) instead of c[ECMECH_NN_INDX(0,0,3)] makes
 * the extensive algebraic expressions in rotation derivative calculations much more
 * readable and maintainable, at the cost of slightly increased stack usage for local
 * variable storage.
 * 
 * @note This file contains only variable declarations, not definitions. It relies on
 *       a variable `c` (pointer to double array) being defined in the including scope.
 * @note The variable `ndim` must equal 3 for proper indexing.
 * 
 * @ingroup ECMech_utilities
 * 
 * @see ECMech_util.h for functions that include this header
 * @see vad_vars_set.h for deviatoric vector components (sample frame)
 * @see vadl_vars_set.h for deviatoric vector components (lattice frame)
 */

double  c11 = c[ECMECH_NN_INDX(0, 0, ndim)];
double  c21 = c[ECMECH_NN_INDX(1, 0, ndim)];
double  c31 = c[ECMECH_NN_INDX(2, 0, ndim)];
double  c12 = c[ECMECH_NN_INDX(0, 1, ndim)];
double  c22 = c[ECMECH_NN_INDX(1, 1, ndim)];
double  c32 = c[ECMECH_NN_INDX(2, 1, ndim)];
double  c13 = c[ECMECH_NN_INDX(0, 2, ndim)];
double  c23 = c[ECMECH_NN_INDX(1, 2, ndim)];
double  c33 = c[ECMECH_NN_INDX(2, 2, ndim)];
