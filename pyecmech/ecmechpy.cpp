/**
 * @file ecmechpy.cpp
 *
 * @brief Implementation of `pyECMech` (declared in `ecmechpy.hpp`).
 *
 * Three things live here:
 * -# `check2D_dim()` -- a small shape-validation helper used by `solve()` to turn a
 *    numpy shape mismatch into a clear Python exception instead of an out-of-bounds
 *    read once the raw buffer pointer is dereferenced.
 * -# The constructor's `model_name` dispatch table -- a fixed set of short, memorable
 *    Python-facing names (e.g. `"voce_fcc_norm"`) mapped onto the corresponding
 *    `ecmech::makeMatModel()` string (e.g. `"evptn_FCC_A"`; see
 *    `cases/ECMech_cases_fcc_defs.h` / `cases/ECMech_cases_bcc_defs.h` for what each one
 *    physically represents).
 * -# `getHistoryInfo()` / `solve()` -- thin passthroughs to the underlying
 *    `matModelBase` that convert between `std::vector`/raw pointers and numpy arrays.
 */

#include "ECMech_cases.h"

#include <random>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>

#include "ecmechpy.hpp"

namespace
{
   /**
    * @brief Validate that a numpy array is 2-D with the expected shape, throwing a
    * descriptive exception otherwise.
    * @param var Array to check.
    * @param dim1 Required size of dimension 0 (typically `nPassed`).
    * @param dim2 Required size of dimension 1 (typically a width constant like
    * `ecmech::nsvp`).
    * @param str Human-readable variable name, used only to make the thrown message
    * identify which argument failed.
    * @return `true` if the shape matches (the function never returns `false`; a mismatch
    * throws instead).
    * @throws std::runtime_error if `var` isn't 2-D or its shape doesn't match `(dim1,
    * dim2)`.
    */
   bool check2D_dim(py_darray &var, const int dim1, const int dim2, const std::string str)
   {
      pybind11::buffer_info buf1 = var.request();

      if (buf1.ndim != 2)
      {
         std::string err = "Number of dimensions must be two for variable: " + str;
         throw std::runtime_error(err.c_str());
      }

      if (buf1.shape[0] != dim1 || buf1.shape[1] != dim2)
      {
         const int bs0 = buf1.shape[0];
         const int bs1 = buf1.shape[1];
         std::string err = "Dimensions for variable: " + str + " should be ("
                                  + std::to_string(dim1) + " , " + std::to_string(dim2) + ") but got: ("
                                  + std::to_string(bs0) + " , " + std::to_string(bs1) + ")";
         throw std::runtime_error(err.c_str());
      }
      return true;
   }
}

pyECMech::pyECMech(std::string model_name, py_darray &params)
{

   pybind11::buffer_info buf1 = params.request();

   if (buf1.ndim != 1)
   {
      throw std::runtime_error("Number of dimensions must be one");
   }
   const int lparam = buf1.shape[0];
   double *param_data = (double*) buf1.ptr;

   std::vector<double> cparams;
   cparams.resize(lparam);

   for (int i = 0; i < lparam; i++)
   {
      cparams[i] = param_data[i];
   }

   // ecmech::matModelBase::initFromParams() also accepts integer "options" and string
   // arguments for models that need them; no currently-supported model here does, so
   // these stay empty. `opts`/`strs` are kept (rather than removed) so that hooking up a
   // future model needing them is a small, local change.
   std::vector<int> opts;
   std::vector<std::string> strs;

   // Friendly Python-facing model name -> ecmech::makeMatModel() string. Keep this list
   // in sync with the one in the pyECMech binding's pbdoc string in
   // ecmech_pybind11.cpp, which is what Python's help() actually shows users.
   if(std::string(model_name) == "voce_fcc_norm") {
      model = ecmech::makeMatModel("evptn_FCC_A");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "voce_nl_fcc_norm") {
      model = ecmech::makeMatModel("evptn_FCC_AH");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "voce_bcc_norm") {
      model = ecmech::makeMatModel("evptn_BCC_A");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "voce_nl_bcc_norm") {
      model = ecmech::makeMatModel("evptn_BCC_AH");    
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "km_bal_dd_fcc_norm") {
      model = ecmech::makeMatModel("evptn_FCC_B");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "km_bal_dd_bcc_norm") {
      model = ecmech::makeMatModel("evptn_BCC_B");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "km_bal_dd_hcp_norm") {
      model = ecmech::makeMatModel("evptn_HCP_A");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "oro_dd_bcc_iso_norm") {
      model = ecmech::makeMatModel("evptn_BCC_C");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "oro_dd_bcc_aniso_norm") {
      model = ecmech::makeMatModel("evptn_BCC_D");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "oro_dd_bcc_24_iso_norm") {
      model = ecmech::makeMatModel("evptn_BCC_C_24");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "oro_dd_bcc_24_aniso_norm") {
      model = ecmech::makeMatModel("evptn_BCC_D_24");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "oro_dd_bcc_aniso_non_schmid") {
      model = ecmech::makeMatModel("evptn_BCC_E");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "bcc_md") {
      model = ecmech::makeMatModel("evptn_BCC_MD");
      model->initFromParams(opts, cparams, strs);
   } else {
      throw std::runtime_error("Provided an unknown model name");
   }

   // Python runs should always be on the CPU for now unless we end up doing some form of memory
   // management on our side of things.
   model->setExecutionStrategy(ecmech::ExecutionStrategy::CPU);
   model->complete();
}

// See ecmechpy.hpp for the full @brief/@return doc -- this just forwards to
// matModelBase::getHistInfo() and converts the resulting std::vector<double> of initial
// values into a numpy array.
std::tuple<std::vector<std::string>, py_darray, std::vector<bool>, std::vector<bool>>
pyECMech::getHistoryInfo()
{
   std::vector<std::string> names;
   std::vector<bool> state;
   std::vector<bool> plot;
   std::vector<double> vals;

   model->getHistInfo(names, vals, plot, state);

   py_darray py_vals = pybind11::array(vals.size(), vals.data());

   return std::make_tuple(names, py_vals, plot, state);
}

// See ecmechpy.hpp for the full @brief/@param doc. Every array is validated up front via
// check2D_dim() (so a shape mistake raises a clear Python exception naming the offending
// argument) and then handed to matModelBase::getResponseECM() as raw pointers for all
// nPassed points at once; the trailing `nullptr` is the (unused here) tangent-stiffness
// output, which this binding doesn't expose.
void pyECMech::solve(double dt,
                     py_darray &def_rate_dev6_vol_sample, // defRate,
                     py_darray &spin_vec_sample, // spin
                     py_darray &volRatio,
                     py_darray &internal_energy,
                     py_darray &cauchy_stress_dev6_pressure,
                     py_darray &hist,
                     py_darray &temp_k,
                     py_darray &sdd,
                     const int nPassed)
{
   // Check that the dimensions for everything is correct
   check2D_dim(def_rate_dev6_vol_sample, nPassed, ecmech::nsvp, "def_rate_dev6_vol_sample");
   check2D_dim(spin_vec_sample, nPassed, ecmech::nwvec, "spin_vec_sample");
   check2D_dim(volRatio, nPassed, ecmech::nvr, "volRatio");
   check2D_dim(internal_energy, nPassed, ecmech::ne, "internal_energy");
   check2D_dim(cauchy_stress_dev6_pressure, nPassed, ecmech::nsvp, "cauchy_stress_dev6_pressure");
   check2D_dim(hist, nPassed, model->getNumHist(), "hist");
   check2D_dim(temp_k, nPassed, 1, "temp_k");
   check2D_dim(sdd, nPassed, ecmech::nsdd, "sdd");
   
   model->getResponseECM(dt,
                         (double*) def_rate_dev6_vol_sample.request().ptr,
                         (double*) spin_vec_sample.request().ptr,
                         (double*) volRatio.request().ptr,
                         (double*) internal_energy.request().ptr,
                         (double*) cauchy_stress_dev6_pressure.request().ptr,
                         (double*) hist.request().ptr,
                         (double*) temp_k.request().ptr,
                         (double*) sdd.request().ptr,
                         nullptr,
                         nPassed);
}