/**
 * @file setup_base.h
 *
 * @brief Foundational fragment of the `setup_*.h` family: a `DUMPVEC` debug-printing
 * macro plus the baseline scalars (`density0`, `cvav`, `tolerance`) most test cases
 * start from.
 *
 * **The `setup_*.h` inclusion pattern** (documented in full once here; every other
 * `setup_*.h`/`setup_kin_*.h` file in this directory follows it and just links back to
 * this explanation): these are *not* standalone, independently compilable headers.
 * Each one is a raw block of C++ (often just a bare `{ ... }` scope) meant to be
 * `#include`d directly inside the body of a `TEST(...)` (or similar function scope) in
 * one of the `test_*.cxx` files, at a point where the surrounding code has already
 * declared whatever names that fragment needs (e.g. an `elastN`/`kinetics`/`slipGeom`
 * object to call `.setParams()` on, or a `params`/`opts`/`strs` triple to build up a
 * flat parameter vector for `matModelBase::initFromParams()`). The same fragment is
 * often `#include`d multiple times in the same translation unit (see `test_evptn.cxx`
 * and `test_updst.cxx`), once per `#if`/`#elif` branch selecting which concrete model
 * to test -- which is why these are text-substitution fragments rather than functions:
 * each inclusion site can be inside a different `#if KIN_TYPE == ...` branch, closing
 * over different template types.
 *
 * Most of the `setup_kin_*.h`/`setup_elastn*.h`/`setup_eos.h`/`setup_slipGeom*.h`
 * fragments branch on whether the including file has `#define`d `STACK_PARAMS`:
 * - **Without `STACK_PARAMS`** (the "white-box" path, e.g. `test_evptn.cxx`): the
 *   fragment calls `.setParams(paramsThese)` directly on an already-declared component
 *   object (`elastN`, `eos`, `slipGeom`, or `kinetics`), exercising that component's
 *   own public API in isolation.
 * - **With `STACK_PARAMS`** (the "black-box" path, e.g. `test_updst.cxx`,
 *   `test_orowan_px.cxx`): the fragment instead appends `paramsThese` onto an
 *   including-scope `std::vector<double> params`, building up the single flat
 *   parameter array that `matModelBase::initFromParams()` expects (see
 *   `ECMech_evptnWrap.h::matModel::initFromParams` for the concatenation order:
 *   density0/cvav/tolerance, then slip-geometry, elastic, kinetics, and finally the
 *   remaining EOS parameters) -- exercising the full model's public, string-driven
 *   construction path instead of the individual components.
 *
 * This file itself declares the shared `density0`/`cvav`/`tolerance` scalars (consumed
 * directly by `setup_eos.h`, and, in `STACK_PARAMS` mode, pushed as the first three
 * entries of `params` by the including test file) and the `DUMPVEC` macro used to print
 * the assembled `opts`/`params`/`strs` vectors for debugging.
 */

/** @brief Print a `std::vector`-like container `a` as a single comma-separated line prefixed with `aname`, for debugging the assembled `opts`/`params`/`strs` vectors. */
#define DUMPVEC(aname, a) std::cout << "# " << aname << " : "; for (unsigned int iThing = 0; \
                                                                    iThing<a.size(); \
                                                                           ++iThing) { if (iThing) { std::cout << ","; \
                                                                                       } std::cout << a[iThing]; \
   } std::cout << std::endl;

/** @brief Reference density (ρ₀) and specific heat (cᵥ), the first two entries `matModel::initFromParams` expects; consumed directly by `setup_eos.h`. */
double density0 = 3.0, cvav = 2.0e-5;
/** @brief Nonlinear-solver convergence tolerance, the third entry `matModel::initFromParams` expects. */
double tolerance = 1e-10;

