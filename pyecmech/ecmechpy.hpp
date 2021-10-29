#pragma once

#include "ECMech_core.h"
#include "ECMech_matModelBase.h"
#include "ecmech_classes.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include<pybind11/numpy.h>

typedef typename pybind11::array_t<double> py_darray;
typedef typename pybind11::array_t<int32_t> py_iarray;

class pyECMech
{
   private:
      ecmech::matModelBase* model = nullptr;
   public:
      pyECMech(std::string model_name, py_darray &params);

      std::tuple<std::vector<std::string>, py_darray, std::vector<bool>, std::vector<bool>>
      getHistoryInfo();

      int getNumberHistory() { return model->getNumHist(); }

      void solve(double dt,
                 py_darray &d_svec_kk_sm,
                 py_darray &w_veccp_sm,
                 py_darray &volRatio,
                 py_darray &eInt,
                 py_darray &stressSvecP,
                 py_darray &hist,
                 py_darray &tkelv,
                 py_darray &sddv,
                 const int nPassed);

      ~pyECMech()
      {
         delete model;
      }
};

#if defined(ECMECH_PYDEV)
class pyECMechDev
{
   private:
      pyevptn_base* model = nullptr;
   public:
      pyECMechDev(std::string model_name, py_darray &params);

      std::tuple<std::vector<std::string>, py_darray, std::vector<bool>, std::vector<bool>>
      getHistoryInfo();

      void setup(double dt,
                 double tolerance,
                 py_darray &d_svec_kk_sm, // defRate,
                 py_darray &w_veccp_sm, // spin
                 py_darray &volRatio,
                 py_darray &eInt,
                 py_darray &stressSvecP,
                 py_darray &hist,
                 double &tkelv);
      void computeRJ(py_darray &resid,
                     py_darray &J,
                     py_darray &x);
      
      void getState(const py_darray &x,
                    py_darray &eInt,
                    py_darray &stressSvecP,
                    py_darray &hist,
                    double &tkelv,
                    py_darray &sdd);

      ~pyECMechDev()
      {
         delete model;
      }
};

#endif
