/**
 * @file ecmech_pybind11.cpp
 *
 * @brief `PYBIND11_MODULE` definition that builds the `pyecmech` Python extension
 * module: binds the `pyECMech` class (declared in `ecmechpy.hpp`) and exposes a
 * `pyecmech.constants` submodule mirroring the dimension/tolerance constants from
 * `ECMech_const.h`, so Python callers can correctly size numpy arrays for `solve()`
 * without hardcoding magic numbers.
 *
 * @note The `R"pbdoc(...)"` strings attached to each binding below *are* the Python
 * docstrings users see via `help(pyecmech.pyECMech)` -- keep them in sync with the
 * Doxygen docs on the corresponding declarations in `ecmechpy.hpp`.
 */

#include <pybind11/pybind11.h>
#include<pybind11/numpy.h>

#include "ECMech_core.h"
#include "ECMech_const.h"

#include "ecmechpy.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pyecmech, m) {
    m.doc() = R"pbdoc(
        ECMech Python Bindings
        -----------------------
        .. currentmodule:: pyecmech
        .. autosummary::
           :toctree: _generate
           pyECMech
    )pbdoc";
    // Each attribute below mirrors a same-named constant from ECMech_const.h -- see that
    // header for the full explanation of each. In short: nsvec/nsvec2/ntvec/nvr/ne/nsvp/
    // nwvec/nsdd/qdim are array-width constants needed to correctly shape the numpy
    // arrays passed to pyECMech.solve(), and dbl_tiny_sqrt/gam_ratio_*/ln_gam_ratio_min
    // are the slip-rate-ratio tolerances used internally by the kinetics models'
    // rate-independent/overflow-guard logic (exposed here mainly so Python-side
    // diagnostics can compare against the same thresholds the C++ core uses).
    py::module constants = m.def_submodule("constants"); // create namespace
    constants.attr("nsvec") = &ecmech::nsvec;
    constants.attr("nsvec2") = &ecmech::nsvec2;
    constants.attr("ntvec") = &ecmech::ntvec;
    constants.attr("nvr") = &ecmech::nvr;
    constants.attr("ne") = &ecmech::ne;
    constants.attr("nsvp") = &ecmech::nsvp;
    constants.attr("nwvec") = &ecmech::nwvec;
    constants.attr("nsdd") = &ecmech::nsdd;
    constants.attr("qdim") = &ecmech::qdim;
    constants.attr("dbl_tiny_sqrt") = &ecmech::idp_tiny_sqrt;
    constants.attr("gam_ratio_ovffx") = &ecmech::gam_ratio_ovffx;
    constants.attr("gam_ratio_ovf") = &ecmech::gam_ratio_ovf;
    constants.attr("gam_ratio_min") = &ecmech::gam_ratio_min;
    constants.attr("ln_gam_ratio_min") = &ecmech::ln_gam_ratio_min;

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
                                     km_bal_dd_hcp_norm,
                                     oro_dd_bcc_iso_norm,
                                     oro_dd_bcc_aniso_norm,
                                     oro_dd_bcc_24_iso_norm,
                                     oro_dd_bcc_24_aniso_norm,
                                     oro_dd_bcc_aniso_non_schmid,
                                     bcc_md,
                                     where voce refers to a Voce hardening law with power law slip kinetics,
                                     voce_nl refers to a nonlinear Voce hardening law with power law slip kinetics,
                                     km_bal_dd refers to a single Kocks-Mecking dislocation density hardening with
                                     balanced thermally activated MTS-like slip kinetics with phonon drag effects,
                                     oro_dd refers to a Orowanian slip kinetics-type model with a dislocation density
                                     hardening model with individual slip system DD evolution (iso and aniso options here
                                     refer to whether the hardening model is isotropic or anisotropic),
                                     non_schmid refers to a slip system construction based on non-schmid formulations popular
                                     with BCC materials,
                                     bcc_md refers to a BCC mobile-dislocation (pencil-glide) kinetics model with its own
                                     mobile/total dislocation density hardening evolution,
                                     and norm refers an implicit beginning of time step hardening state update and
                                     an implicit end of time step coupled elastic strain and lattice rotation update.
            py_darray params - model parameters for the provided model name, in the order the underlying model's
                               initFromParams() expects (see ecmechpy.hpp / ECMech_evptnWrap.h for the general
                               concatenation order).)pbdoc")
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
                 Advances nPassed material points through one time step of size dt, in place.
                 All array arguments below are shaped (nPassed, width); see pyecmech.constants
                 for the width constants (nsvp, nwvec, nvr, ne, nsdd).

                 double dt, // delta time
                 py_darray& def_rate_dev6_vol_sample, // deformation rate in sample frame [in]
                 py_darray& spin_vec_sample, // spin in sample frame [in]
                 py_darray& volRatio, // volume ratio bookkeeping [rel_vol_n, rel_vol_n+1, rate, delta] [in/out]
                 py_darray& internal_energy, // internal energy [in/out]
                 py_darray& cauchy_stress_dev6_pressure, // stress deviatoric vector + pressure term [in/out]
                 py_darray& hist, // history variable [in/out]
                 py_darray& temp_k, // current temperature in kelvin [in/out]
                 py_darray& sdd, // auxiliary derived quantities, e.g. shear modulus [out]
                 int nPassed // number of material points in every array above
             )pbdoc");
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
