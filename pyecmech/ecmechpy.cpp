
#include "ECMech_cases.h"

#include <random>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>

#include "ecmechpy.hpp"

typedef ecmech::evptn::matModel<ecmech::SlipGeom_BCC_A, ecmech::Kin_FCC_A, 
               ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> >
               matModelEvptn_BCC_B;

typedef ecmech::evptn::matModel<ecmech::SlipGeom_BCC_A, ecmech::Kin_FCC_AH, 
               ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> >
               matModelEvptn_BCC_BH;

#if defined(ECMECH_PYDEV)

typedef pyEvptn_norm<ecmech::SlipGeomFCC, ecmech::Kin_FCC_A, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_FCC_A;
typedef pyEvptn_norm<ecmech::SlipGeomFCC, ecmech::Kin_FCC_AH, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_FCC_AH;
typedef pyEvptn_norm<ecmech::SlipGeomFCC, ecmech::Kin_FCC_B, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_FCC_B;
typedef pyEvptn_norm<ecmech::SlipGeom_BCC_A, ecmech::Kin_BCC_A, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_BCC_A;
typedef pyEvptn_norm<ecmech::SlipGeom_BCC_A, ecmech::Kin_FCC_A, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_BCC_B;
typedef pyEvptn_norm<ecmech::SlipGeom_BCC_A, ecmech::Kin_FCC_AH, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_BCC_BH;
typedef pyEvptn_norm<ecmech::SlipGeom_HCP_A, ecmech::Kin_HCP_A, ecmech::evptn::ThermoElastNHexag, ecmech::EosModelConst<false> > pyMatModelEvptn_HCP_A;

#endif

namespace
{
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

   std::vector<int> opts;
   std::vector<std::string> strs;

   // Only care about these for now. We can add ones later on.
   if(std::string(model_name) == "voce_fcc_norm") {
      model = ecmech::makeMatModel("evptn_FCC_A");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "voce_nl_fcc_norm") {
      model = ecmech::makeMatModel("evptn_FCC_AH");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "voce_bcc_norm") {
      matModelEvptn_BCC_B* mmECMEvptn = new matModelEvptn_BCC_B();
      model = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "voce_nl_bcc_norm") {
      matModelEvptn_BCC_BH* mmECMEvptn = new matModelEvptn_BCC_BH();
      model = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn);     
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "km_bal_dd_fcc_norm") {
      model = ecmech::makeMatModel("evptn_FCC_B");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "km_bal_dd_bcc_norm") {
      model = ecmech::makeMatModel("evptn_BCC_A");
      model->initFromParams(opts, cparams, strs);
   } else if (std::string(model_name) == "km_bal_dd_hcp_norm") {
      model = ecmech::makeMatModel("evptn_FCC_A");
      model->initFromParams(opts, cparams, strs);
   } else {
      throw std::runtime_error("Provided an unknown model name");
   }

   // Python runs should always be on the CPU for now unless we end up doing some form of memory
   // management on our side of things.
   model->setExecutionStrategy(ecmech::ExecutionStrategy::CPU);
   model->complete();
}

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

void pyECMech::solve(double dt,
                     py_darray &d_svec_kk_sm, // defRate,
                     py_darray &w_veccp_sm, // spin
                     py_darray &volRatio,
                     py_darray &eInt,
                     py_darray &stressSvecP,
                     py_darray &hist,
                     py_darray &tkelv,
                     py_darray &sdd,
                     const int nPassed)
{
   // Check that the dimensions for everything is correct
   check2D_dim(d_svec_kk_sm, nPassed, ecmech::nsvp, "d_svec_kk_sm");
   check2D_dim(w_veccp_sm, nPassed, ecmech::nwvec, "w_veccp_sm");
   check2D_dim(volRatio, nPassed, ecmech::nvr, "volRatio");
   check2D_dim(eInt, nPassed, ecmech::ne, "eInt");
   check2D_dim(stressSvecP, nPassed, ecmech::nsvp, "stressSvecP");
   check2D_dim(hist, nPassed, model->getNumHist(), "hist");
   check2D_dim(tkelv, nPassed, 1, "tkelv");
   check2D_dim(sdd, nPassed, ecmech::nsdd, "sdd");
   
   model->getResponseECM(dt,
                         (double*) d_svec_kk_sm.request().ptr,
                         (double*) w_veccp_sm.request().ptr,
                         (double*) volRatio.request().ptr,
                         (double*) eInt.request().ptr,
                         (double*) stressSvecP.request().ptr,
                         (double*) hist.request().ptr,
                         (double*) tkelv.request().ptr,
                         (double*) sdd.request().ptr,
                         nullptr,
                         nPassed);
}

#if defined(ECMECH_PYDEV)

pyECMechDev::pyECMechDev(std::string model_name, py_darray &params)
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

   // Only care about these for now. We can add ones later on.
   if(std::string(model_name) == "voce_fcc_norm") {
      model = new pyMatModelEvptn_FCC_A(cparams);
   } else if (std::string(model_name) == "voce_nl_fcc_norm") {
      model = new pyMatModelEvptn_FCC_AH(cparams);
   } else if (std::string(model_name) == "voce_bcc_norm") {
      model = new pyMatModelEvptn_BCC_B(cparams);
   } else if (std::string(model_name) == "voce_nl_bcc_norm") {
      model = new pyMatModelEvptn_BCC_BH(cparams);
   } else if (std::string(model_name) == "km_bal_dd_fcc_norm") {
      model = new pyMatModelEvptn_FCC_B(cparams);
   } else if (std::string(model_name) == "km_bal_dd_bcc_norm") {
      model = new pyMatModelEvptn_BCC_A(cparams);
   } else if (std::string(model_name) == "km_bal_dd_hcp_norm") {
      model = new pyMatModelEvptn_HCP_A(cparams);
   } else {
      throw std::runtime_error("Provided an unknown model name");
   }
}

std::tuple<std::vector<std::string>, py_darray, std::vector<bool>, std::vector<bool>>
pyECMechDev::getHistoryInfo()
{
   std::vector<std::string> names;
   std::vector<bool> state;
   std::vector<bool> plot;
   std::vector<double> vals;

   model->getHistoryInfo(names, vals, plot, state);

   py_darray py_vals = pybind11::array(vals.size(), vals.data());

   return std::make_tuple(names, py_vals, plot, state);
}

void pyECMechDev::setup(double dt,
                     double tolerance,
                     py_darray &d_svec_kk_sm, // defRate,
                     py_darray &w_veccp_sm, // spin
                     py_darray &volRatio,
                     py_darray &eInt,
                     py_darray &stressSvecP,
                     py_darray &hist,
                     double& tkelv)
{

   model->setup(dt, tolerance,
               (double*) d_svec_kk_sm.request().ptr,
               (double*) w_veccp_sm.request().ptr,
               (double*) volRatio.request().ptr,
               (double*) eInt.request().ptr,
               (double*) stressSvecP.request().ptr,
               (double*) hist.request().ptr,
               tkelv);
}

void pyECMechDev::computeRJ(py_darray &resid,
                    py_darray &J,
                    py_darray &x)
{
   model->computeRJ((double*) resid.request().ptr,
                    (double*) J.request().ptr,
                    (double*) x.request().ptr);
}

void pyECMechDev::getState(const py_darray &x,
                        py_darray &eInt,
                        py_darray &stressSvecP,
                        py_darray &hist,
                        double& tkelv,
                        py_darray &sdd)
{
   model->getState((double*) x.request().ptr,
                   (double*) eInt.request().ptr,
                   (double*) stressSvecP.request().ptr,
                   (double*) hist.request().ptr,
                   tkelv,
                   (double*) sdd.request().ptr);
}

#endif
