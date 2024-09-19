#pragma once

#include "ECMech_core.h"
#include "ECMech_matModelBase.h"

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
                 py_darray &def_rate_dev6_vol_sample,
                 py_darray &spin_vec_sample,
                 py_darray &volRatio,
                 py_darray &internal_energy,
                 py_darray &cauchy_stress_dev6_pressure,
                 py_darray &hist,
                 py_darray &temp_k,
                 py_darray &sddv,
                 const int nPassed);

      ~pyECMech()
      {
         delete model;
      }
};
