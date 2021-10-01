#ifndef ECMECHPY_HPP
#define ECMECHPY_HPP

#include "ecmech_classes.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include<pybind11/numpy.h>

typedef typename pybind11::array_t<double> py_darray;
typedef typename pybind11::array_t<int32_t> py_iarray;

class ECMechPy
{
   private:
      pyevptn_base* model = nullptr;
   public:
      ECMechPy(std::string model_name, py_darray &params);

      std::tuple<std::vector<std::string>, py_darray, std::vector<bool>, std::vector<bool>>
      getHistInfo();

      void setup(double dt,
                 double tolerance,
                 py_darray d_svec_kk_sm, // defRate,
                 py_darray w_veccp_sm, // spin
                 py_darray volRatio,
                 py_darray eInt,
                 py_darray stressSvecP,
                 py_darray hist,
                 double& tkelv,
                 py_darray sdd,
                 py_darray mtanSD);
      void computeRJ(py_darray &resid,
                     py_darray &J,
                     py_darray &x);
      
      void getState(const py_darray &x,
                    py_darray &eInt,
                    py_darray &stressSvecP,
                    py_darray &hist,
                    double& tkelv,
                    py_darray &sdd);

      ~ECMechPy()
      {
         delete model;
      }
};


#endif