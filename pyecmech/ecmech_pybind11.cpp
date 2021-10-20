#include <pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include "ECMech_const.h"
#include "ecmechpy.hpp"

namespace py = pybind11;

PYBIND11_MODULE(ecmech, m) {
    m.doc() = R"pbdoc(
        ECMech Python Bindings
        -----------------------
        .. currentmodule:: ecmech
        .. autosummary::
           :toctree: _generate
           ECMech
    )pbdoc";
    py::module constants = m.def_submodule("constants"); // create namespace
    constants.attr("nsvec") = &ecmech::nsvec;
    constants.attr("nsvec2") = &ecmech::nsvec2;
    constants.attr("ntvec") = &ecmech::ntvec;
    constants.attr("nvr") = &ecmech::nvr;
    constants.attr("ne") = &ecmech::ne;
    constants.attr("nsvp") = &ecmech::nsvp;
    constants.attr("nwvec") = &ecmech::nwvec;
    constants.attr("nsdd") = &ecmech::nsdd;

    py::class_<pyECMech>(m, "pyECMech", "Provides pyECMech")
        .def(py::init([](std::string &model_name, py_darray &params) {
            return std::unique_ptr<pyECMech>(new pyECMech(model_name, params));
        }),
        R"pbdoc(
            std::string model_name - model name choices are:
                                     voce_fcc_norm,
                                     voce_nl_fcc_norm,
                                     voce_bcc_norm,
                                     voce_nl_bcc_norm,
                                     km_bal_dd_fcc_norm,
                                     km_bal_dd_bcc_norm,
                                     km_bal_dd_hcp_norm
                                     where voce refers to a Voce hardening law with power law slip kinetics,
                                     voce_nl refers to a nonlinear Voce hardening law with power law slip kinetics,
                                     km_bal_dd refers to a single Kocks-Mecking dislocation density hardening with
                                     balanced thermally activated MTS-like slip kinetics with phonon drag effects,
                                     and norm refers an implicit beginning of time step hardening state update and
                                     an implicit end of time step coupled elastic strain and lattice rotation update.
            py_darray params - model parameters for the provided model name.)pbdoc")
        .def("getHistoryInfo", &pyECMech::getHistoryInfo, py::return_value_policy::take_ownership,
             R"pbdoc(
                Output: names, vals, plot, state
                names: the history names as a list of strings, 
                vals: the initial values of the history name as numpy array, 
                plot: whether the history variable should be plotted as a list of booleans,
                state: whether the history variable is a state variable as a list of booleans.
                )pbdoc")
        .def("getNumberHistory", &pyECMech::getNumberHistory, py::return_value_policy::take_ownership,
            R"pbdoc(
                Output: numHist
                numHist: the number of history variables
            )pbdoc")
        .def("solve", &pyECMech::solve,
             R"pbdoc(
                 double dt, // delta time
                 py_darray& d_svec_kk_sm, // deformation rate in sample frame
                 py_darray& w_veccp_sm, // spin in sample rate
                 py_darray& volRatio, // volume ratio
                 py_darray& eInt, // internal energy
                 py_darray& stressSvecP, // stress deviatoric vector + pressure term
                 py_darray& hist, // history variable
                 py_darray& tkelv // current temperature in kelvin
                 py_darray& sdd // sdd array
             )pbdoc");

#if defined(ECMECH_PYDEV)
    py::class_<pyECMechDev>(m, "pyECMechDev", "Provides pyECMechDev")
        .def(py::init([](std::string &model_name, py_darray &params) {
            return std::unique_ptr<pyECMechDev>(new pyECMechDev(model_name, params));
        }),
        R"pbdoc(
            std::string model_name - model name choices are:
                                     voce_fcc_norm,
                                     voce_nl_fcc_norm,
                                     voce_bcc_norm,
                                     voce_nl_bcc_norm,
                                     km_bal_dd_fcc_norm,
                                     km_bal_dd_bcc_norm,
                                     km_bal_dd_hcp_norm
                                     where voce refers to a Voce hardening law with power law slip kinetics,
                                     voce_nl refers to a nonlinear Voce hardening law with power law slip kinetics,
                                     km_bal_dd refers to a single Kocks-Mecking dislocation density hardening with
                                     balanced thermally activated MTS-like slip kinetics with phonon drag effects,
                                     and norm refers an implicit beginning of time step hardening state update and
                                     an implicit end of time step coupled elastic strain and lattice rotation update.
            py_darray params - model parameters for the provided model name.)pbdoc")
        .def("getHistoryInfo", &pyECMechDev::getHistoryInfo, py::return_value_policy::take_ownership,
             R"pbdoc(
                Output: names, vals, plot, state
                names: the history names as a list of strings, 
                vals: the initial values of the history name as numpy array, 
                plot: whether the history variable should be plotted as a list of booleans,
                state: whether the history variable is a state variable as a list of booleans.
                )pbdoc")
        .def("setup", &pyECMechDev::setup,
             R"pbdoc(
                 double dt, // delta time
                 double tolerance, // solver tolerance - not used
                 py_darray d_svec_kk_sm, // deformation rate in sample frame
                 py_darray w_veccp_sm, // spin in sample rate
                 py_darray volRatio, // volume ratio
                 py_darray eInt, // internal energy
                 py_darray stressSvecP, // stress deviatoric vector + pressure term
                 py_darray hist, // history variable
                 double& tkelv // current temperature in kelvin
             )pbdoc")
        .def("computeRJ", &pyECMechDev::computeRJ,
             R"pbdoc(
                py_darray &resid, // residual of system
                py_darray &J, // Jacobian of system
                py_darray &x // input solution vector
             )pbdoc")
        .def("getState", &pyECMechDev::getState,
             R"pbdoc(
                const py_darray &x, // input solution vector
                py_darray &eInt, // internal energy
                py_darray &stressSvecP, // stress deviatoric vector + pressure term
                py_darray &hist, // history variable
                double& tkelv, // current temperature
                py_darray &sdd // sdd array
             )pbdoc");
#endif

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
