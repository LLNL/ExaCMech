/**
 * @file ECMech_cases_util.h
 * @brief Helper for building the name-keyed parameter/history-index lookup map returned
 * by the `modelParamIndexMap*` functions declared in ECMech_cases.h.
 */

#pragma once

#include "ECMech_cases.h"

/**
 * @brief Builds a `std::string`-keyed lookup map of parameter counts and history-array
 * indices for a concrete model type `T`.
 *
 * `T` is expected to be an `evptn::matModel<SlipGeom, Kinetics, ThermoElastN, EosModel>`
 * instantiation (see ECMech_evptnWrap.h). The resulting map (`m_maps`) lets a calling
 * code (e.g. ExaConstit) look up how many parameters a model needs and where each
 * history quantity lives in the flattened history array, without needing compile-time
 * knowledge of `T`'s C++ type — this is exactly what the `modelParamIndexMap*` functions
 * in `cases/ECMech_cases*.cxx` return.
 *
 * @tparam T Concrete `evptn::matModel<...>` type to introspect.
 */
template<class T>
class NumParamIndexInfo {
public:
    /**
     * @brief Populate #m_maps with `T`'s parameter counts and history-array indices.
     */
    NumParamIndexInfo() {
        m_maps["num_params"] = T::nParams;
        m_maps["num_params_eos"] = T::nParamsEOS;
        m_maps["num_params_slip_geom"] = T::nParamsSlipGeom;
        m_maps["num_params_slip_kinetics"] = T::nParamsKinetics;
        m_maps["num_params_elasticity"] = T::nParamsThermoElastN;
        m_maps["num_hist"] = T::numHist;
        m_maps["num_hardening"] = T::nH;
        m_maps["num_slip_system"] = T::nslip;
        m_maps["index_effective_shear_rate"] = ecmech::evptn::iHistA_shrateEff;
        m_maps["index_effective_shear"] = ecmech::evptn::iHistA_shrEff;
        m_maps["index_flow_strength"] = ecmech::evptn::iHistA_flowStr;
        m_maps["index_num_func_evals"] = ecmech::evptn::iHistA_nFEval;
        m_maps["index_dev_elas_strain"] = ecmech::evptn::iHistLbE;
        m_maps["index_lattice_ori"] = ecmech::evptn::iHistLbQ;
        m_maps["index_hardness"] = ecmech::evptn::iHistLbH;
        m_maps["index_slip_rates"] = T::iHistLbGdot;
    }
public:
    /**
     * @brief Name -> value lookup table populated by the constructor; see the class
     * documentation for the full list of keys.
     */
    std::map<std::string, size_t> m_maps;
};
