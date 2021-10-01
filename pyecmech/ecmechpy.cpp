#include "ecmechpy.hpp"

#include <random>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>

ECMechPy::ECMechPy(std::string model_name, py_darray &params)
{

   pybind11::buffer_info buf1 = params.request();

   if (buf1.ndim != 1)
   {
      throw std::runtime_error("Number of dimensions must be two");
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
      model = new matModelEvptn_BCC_B(cparams);
   } else if (std::string(model_name) == "voce_nl_bcc_norm") {
      model = new matModelEvptn_BCC_BH(cparams);
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
ECMechPy::getHistInfo()
{
   std::vector<std::string> names;
   std::vector<bool> state;
   std::vector<bool> plot;
   std::vector<double> vals;

   model->getHistInfo(names, vals, plot, state);

   py_darray py_vals = pybind11::array(vals.size(), vals.pointer());

   return std::make_tuple(names, py_vals, plot, state);
}

void ECMechPy::setup(double dt,
                     double tolerance,
                     py_darray d_svec_kk_sm, // defRate,
                     py_darray w_veccp_sm, // spin
                     py_darray volRatio,
                     py_darray eInt,
                     py_darray stressSvecP,
                     py_darray hist,
                     double& tkelv,
                     py_darray sdd,
                     py_darray mtanSD)
{

   model->setup(dt, tolerance,
               (double*) d_svec_kk_sm.request().ptr,
               (double*) w_veccp_sm.request().ptr,
               (double*) volRatio.request().ptr,
               (double*) eInt.request().ptr,
               (double*) stressSvecP.request().ptr,
               (double*) hist.request().ptr,
               tkelv,
               (double*) sdd.request().ptr,
               (double*) mtanSD.request().ptr);
}

void ECMechPy::computeRJ(py_darray &resid,
                    py_darray &J,
                    py_darray &x)
{
   model->computeRJ((double*) resid.request().ptr,
                    (double*) J.request().ptr,
                    (double*) x.request().ptr);
}

void ECMechPy::getState(const py_darray &x,
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