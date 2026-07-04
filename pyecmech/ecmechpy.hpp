/**
 * @file ecmechpy.hpp
 *
 * @brief Declares `pyECMech`, the thin C++ wrapper class that the pybind11 module
 * defined in `ecmech_pybind11.cpp` exposes to Python as `pyecmech.pyECMech`.
 *
 * The rest of ExaCMech is a templated C++ library: a concrete material model is a
 * `matModelBase`-derived type built from a `SlipGeom`/`Kinetics`/`ThermoElastN`/`EosModel`
 * combination (see `ECMech_evptnWrap.h`) and looked up by a model-name string via
 * `ecmech::makeMatModel()` (see `ECMech_cases.h`). `pyECMech` hides all of that behind a
 * small, numpy-friendly surface:
 * -# construct it with a friendly model-name string (see `ecmechpy.cpp` for the
 *    supported names) and a flat `numpy.ndarray` of model parameters,
 * -# call `getHistoryInfo()` once to get the names/initial-values/plot/state metadata
 *    for every history (state) variable the model carries,
 * -# call `solve()` once per time step, passing batches of material points (deformation
 *    rate, spin, volume ratio, ...) and mutating the stress/history/temperature arrays
 *    in place.
 *
 * The implementation lives in `ecmechpy.cpp`; the actual pybind11 bindings (the Python
 * docstrings users see via `help()`) live in `ecmech_pybind11.cpp`.
 */

#pragma once

#include "ECMech_core.h"
#include "ECMech_matModelBase.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include<pybind11/numpy.h>

/** @brief Convenience alias for a numpy array of `double`, used for all floating-point I/O. */
typedef typename pybind11::array_t<double> py_darray;
/** @brief Convenience alias for a numpy array of `int32_t`. Currently unused by `pyECMech` itself, but kept available for future bindings. */
typedef typename pybind11::array_t<int32_t> py_iarray;

/**
 * @brief Python-facing wrapper around a single `ecmech::matModelBase` instance.
 *
 * Owns exactly one heap-allocated material model, selected by name at construction time
 * (see `ecmechpy.cpp` for the recognized `model_name` strings and which
 * `ecmech::makeMatModel()` string each maps to). The model is fully configured
 * (`initFromParams`, `setExecutionStrategy(CPU)`, `complete()`) inside the constructor,
 * so a `pyECMech` object is ready to use as soon as it's built.
 *
 * @note Not copyable/movable at the C++ level (it owns `model` via a raw pointer with no
 * copy/move constructor defined) -- pybind11 only ever constructs it via the `py::init`
 * factory in `ecmech_pybind11.cpp`, so this isn't reachable from Python either.
 */
class pyECMech
{
   private:
      /** @brief The underlying material model, heap-allocated by `ecmech::makeMatModel()` and owned by this object. */
      ecmech::matModelBase* model = nullptr;
   public:
      /**
       * @brief Build and fully configure a material model by name.
       * @param model_name Friendly model identifier (e.g. `"voce_fcc_norm"`); see the
       * dispatch table in the constructor body in `ecmechpy.cpp` for the complete list
       * and which `ecmech::makeMatModel()` string (e.g. `"evptn_FCC_A"`) each one maps to.
       * @param params Flat 1-D array of model parameters, in the order the underlying
       * model's `initFromParams()` expects (density, specific heat, solver tolerance,
       * slip-geometry params, elastic constants, slip-kinetics params, remaining EOS
       * params -- see `ECMech_evptnWrap.h::matModel::initFromParams` for the general
       * concatenation order, and the individual `kinetics/`, `ECMech_elastic.h`, and
       * `ECMech_eosSimple.h` headers for what each model's own slice of the array means).
       * @throws std::runtime_error if `params` isn't 1-D, or if `model_name` isn't one of
       * the recognized strings.
       */
      pyECMech(std::string model_name, py_darray &params);

      /**
       * @brief Report the name, initial value, and plot/state flags for every history
       * (state) variable the model carries.
       * @return Tuple of `(names, vals, plot, state)`:
       * - `names` -- one string per history variable
       * - `vals` -- initial value of each history variable, as a numpy array
       * - `plot` -- whether each variable is generally interesting to plot
       * - `state` -- whether each variable is a true solved-for state variable (as
       *   opposed to a derived/diagnostic quantity)
       *
       * Mirrors `ecmech::matModelBase::getHistInfo()`; call this once up front to learn
       * the per-point history-array layout and initial values expected by `solve()`.
       */
      std::tuple<std::vector<std::string>, py_darray, std::vector<bool>, std::vector<bool>>
      getHistoryInfo();

      /** @brief Number of history (state) variables per material point; matches `getHistoryInfo()`'s array lengths and the required trailing dimension of `solve()`'s `hist` argument. */
      int getNumberHistory() { return model->getNumHist(); }

      /**
       * @brief Advance a batch of `nPassed` material points through one time step.
       *
       * Thin, shape-checked passthrough to `ecmech::matModelBase::getResponseECM()` for
       * all `nPassed` points at once (no tangent-stiffness matrix is computed or
       * returned). All array arguments are laid out as `nPassed` rows by the
       * model-independent widths documented on each parameter below (see
       * `ECMech_const.h` for the named width constants `nsvp`/`nwvec`/`nvr`/`ne`/`nsdd`).
       *
       * @param[in] dt Time-step size.
       * @param[in] def_rate_dev6_vol_sample Deviatoric (6) + volumetric (1) deformation
       * rate in the sample frame, shape `(nPassed, ecmech::nsvp)`.
       * @param[in] spin_vec_sample Sample-frame spin (skew-symmetric part of the velocity
       * gradient, axial-vector form), shape `(nPassed, ecmech::nwvec)`.
       * @param[in,out] volRatio Relative-volume bookkeeping `[rel_vol_n, rel_vol_{n+1},
       * rate, delta]` per point, shape `(nPassed, ecmech::nvr)` -- see `ECMech_const.h`'s
       * `nvr` doc for the exact per-slot meaning.
       * @param[in,out] internal_energy Internal energy per point, shape `(nPassed,
       * ecmech::ne)`; input is beginning-of-step, output is end-of-step.
       * @param[in,out] cauchy_stress_dev6_pressure Cauchy stress (deviatoric 6-vector +
       * pressure) per point, shape `(nPassed, ecmech::nsvp)`; input is beginning-of-step,
       * output is end-of-step.
       * @param[in,out] hist History (state) variables per point, shape `(nPassed,
       * getNumberHistory())`; input is beginning-of-step, output is end-of-step. Layout
       * and initial values come from `getHistoryInfo()`.
       * @param[in,out] temp_k Temperature in Kelvin per point, shape `(nPassed, 1)`;
       * input is beginning-of-step, output is end-of-step (computed from the EOS model).
       * @param[out] sddv Auxiliary derived quantities per point (e.g. shear modulus),
       * shape `(nPassed, ecmech::nsdd)`; see `ECMech_const.h`'s `i_sdd_*` indices.
       * @param[in] nPassed Number of material points contained in every array above.
       * @throws std::runtime_error if any array's shape doesn't match the dimensions
       * above (checked via `check2D_dim()` in `ecmechpy.cpp`).
       */
      void solve(double dt,
                 py_darray &def_rate_dev6_vol_sample,
                 py_darray &spin_vec_sample,
                 py_darray &volRatio,
                 py_darray &internal_energy,
                 py_darray &cauchy_stress_dev6_pressure,
                 py_darray &hist,
                 py_darray &temp_k,
                 py_darray &sddv,
                 const int nPassed);

      /** @brief Frees the owned material model. */
      ~pyECMech()
      {
         delete model;
      }
};
