#pragma once

#include "ECMech_cases.h"

template<class T>
class NumParamIndexInfo {
public:
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
    std::map<std::string, size_t> m_maps;
};