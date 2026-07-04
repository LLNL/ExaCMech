/**
 * @file ECMech_cases.h
 * @brief Model registry: convenience type aliases and string-name lookup functions for
 * building a concrete crystal-plasticity `matModelBase` from a model-name string.
 *
 * A "model" in ExaCMech is the combination of four independently-templated pieces (see
 * `evptn::matModel` in ECMech_evptnWrap.h):
 * 1. A slip geometry (crystal structure + slip system family, see ECMech_slipgeom.h)
 * 2. A kinetics model (slip-rate law + hardening evolution, see `kinetics/`)
 * 3. A thermoelastic model (cubic or hexagonal, see ECMech_elastic.h)
 * 4. An equation of state (see ECMech_eosSimple.h)
 *
 * This header declares a handful of commonly-used aliases for pieces (2)-(4) that are
 * shared across multiple concrete models, then declares the `makeMatModel*` factory
 * functions (implemented in `cases/ECMech_cases*.cxx`) that map a model-name string
 * (e.g. `"evptn_FCC_A"`) to a heap-allocated `matModelBase*`, and the
 * `modelParamIndexMap*` functions that map the same model-name strings to a lookup
 * table of parameter counts and history-array indices for that model.
 *
 * The concrete, crystal-structure-specific model type aliases and the string names that
 * select them live in `cases/ECMech_cases_bcc_defs.h`, `cases/ECMech_cases_fcc_defs.h`,
 * and `cases/ECMech_cases_hcp_defs.h`; those are the files to consult for "what does
 * model name X actually mean physically."
 *
 * @see cases/ECMech_cases_util.h for `NumParamIndexInfo`, which builds the map returned
 * by the `modelParamIndexMap*` functions
 */

// -*-c++-*-

#include "evptn/ECMech_evptn.h"
#include "kinetics/ECMech_kinetics.h"
#include "ECMech_slipgeom.h"
#include "ECMech_evptnWrap.h"

namespace ecmech {

   /**
    * @brief Constant-Gruneisen-parameter equation of state, non-isothermal (temperature
    * evolves with internal energy). Used by essentially all of the concrete models
    * declared in the `cases/` headers.
    * @see EosModelConst in ECMech_eosSimple.h
    */
   using EOS_const_model = EosModelConst<false>;

   // some common kinetic forms
   /**
    * @brief Power-law slip kinetics with linear Voce hardening (saturation-type
    * isotropic hardening with a single hardening state variable per slip system family).
    * @see KineticsVocePL in kinetics/ECMech_kinetics_VocePL.h
    */
   using Kin_Voce   = KineticsVocePL<false>;
   /**
    * @brief Power-law slip kinetics with the nonlinear ("NL") Voce hardening variant,
    * which adds one extra parameter controlling the shape of the hardening-rate curve
    * as it approaches saturation.
    * @see KineticsVocePL in kinetics/ECMech_kinetics_VocePL.h
    */
   using Kin_VoceNL = KineticsVocePL<true>;
   /**
    * @brief Kocks-Mecking balanced-thermally-activated (MTS-like) slip kinetics with a
    * single dislocation-density hardening variable, with the athermal/thermal-activation
    * split enabled (`withGAthermal = true`) and the MTS `p`/`q` exponents left general
    * (not pegged to 1). "TFF" names the first three template booleans in order:
    * `withGAthermal=true, pOne=false, qOne=false`.
    * @see KineticsKMBalD in kinetics/ECMech_kinetics_KMBalD.h
    */
   using Kin_KMBalD_TFF = KineticsKMBalD<true, false, false, false, 1>;
   /**
    * @brief Same as #Kin_KMBalD_TFF but with the athermal/thermal-activation split
    * disabled (`withGAthermal = false`). "FFF" names the first three template booleans:
    * `withGAthermal=false, pOne=false, qOne=false`.
    * @see KineticsKMBalD in kinetics/ECMech_kinetics_KMBalD.h
    */
   using Kin_KMBalD_FFF = KineticsKMBalD<false, false, false, false, 1>;

   /**
    * @brief Anisotropic thermoelastic model for cubic crystal symmetry (3 independent
    * elastic constants).
    * @see ThermoElastNCubic in ECMech_elastic.h
    */
   using EVPTN_cubic = evptn::ThermoElastNCubic;
   /**
    * @brief Anisotropic thermoelastic model for hexagonal crystal symmetry (5 independent
    * elastic constants plus a Gruneisen parameter).
    * @see ThermoElastNHexag in ECMech_elastic.h
    */
   using EVPTN_hex   = evptn::ThermoElastNHexag;

   /**
    * @brief Build a material model from its name, dispatching to the FCC/BCC/HCP
    * families based on a (case-sensitive) substring match on `modelName` — i.e. the name
    * must contain "FCC", "BCC", or "HCP".
    * @param modelName Model name, e.g. `"evptn_FCC_A"`; see the crystal-structure-specific
    * `cases/ECMech_cases_*_defs.h` headers for the full list of recognized names.
    * @return Newly heap-allocated `matModelBase*` for the requested model; ownership
    * transfers to the caller.
    * @note Calls `ECMECH_FAIL` (throwing on host) if `modelName` does not match any
    * recognized model.
    */
   __ecmech_host__
   matModelBase* makeMatModel(const std::string &modelName);
   /**
    * @brief Build an FCC material model from its name.
    * @param modelName Model name, e.g. `"evptn_FCC_A"`; see cases/ECMech_cases_fcc_defs.h.
    * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not
    * recognized.
    */
   __ecmech_host__
   matModelBase* makeMatModelFCC(const std::string &modelName);
   /**
    * @brief Build a BCC material model from its name.
    * @param modelName Model name, e.g. `"evptn_BCC_A"`; see cases/ECMech_cases_bcc_defs.h.
    * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not
    * recognized.
    */
   __ecmech_host__
   matModelBase* makeMatModelBCC(const std::string &modelName);
   /**
    * @brief Build an HCP material model from its name.
    * @param modelName Model name, e.g. `"evptn_HCP_A"`; see cases/ECMech_cases_hcp_defs.h.
    * @return Newly heap-allocated `matModelBase*`, or `nullptr` if `modelName` is not
    * recognized.
    */
   __ecmech_host__
   matModelBase* makeMatModelHCP(const std::string &modelName);

   /**
    * @brief Look up parameter-count and history-array-index information for a model by
    * name, dispatching to the FCC/BCC/HCP families based on a substring match on
    * `modelName` (same dispatch rule as `makeMatModel`).
    * @param modelName Model name, e.g. `"evptn_FCC_A"`.
    * @return Map from a fixed set of descriptive keys (`"num_params"`,
    * `"num_params_eos"`, `"num_params_slip_geom"`, `"num_params_slip_kinetics"`,
    * `"num_params_elasticity"`, `"num_hist"`, `"num_hardening"`, `"num_slip_system"`,
    * `"index_effective_shear_rate"`, `"index_effective_shear"`, `"index_flow_strength"`,
    * `"index_num_func_evals"`, `"index_dev_elas_strain"`, `"index_lattice_ori"`,
    * `"index_hardness"`, `"index_slip_rates"`) to the corresponding count or history
    * index for that model. Lets calling codes (e.g. ExaConstit) introspect a model's
    * layout without compile-time knowledge of its C++ type.
    * @note Calls `ECMECH_FAIL` if `modelName` does not match any recognized model.
    * @see NumParamIndexInfo in cases/ECMech_cases_util.h, which builds this map.
    */
   __ecmech_host__
   std::map<std::string, size_t>
   modelParamIndexMap(const std::string_view &modelName);
   /**
    * @brief Look up parameter-count and history-array-index information for an FCC model
    * by name; see modelParamIndexMap for the returned map's keys.
    * @param modelName Model name, e.g. `"evptn_FCC_A"`.
    * @return Lookup map as described in modelParamIndexMap, or an empty map if
    * `modelName` is not recognized.
    */
   __ecmech_host__
   std::map<std::string, size_t>
   modelParamIndexMapFCC(const std::string_view &modelName);
   /**
    * @brief Look up parameter-count and history-array-index information for a BCC model
    * by name; see modelParamIndexMap for the returned map's keys.
    * @param modelName Model name, e.g. `"evptn_BCC_A"`.
    * @return Lookup map as described in modelParamIndexMap, or an empty map if
    * `modelName` is not recognized.
    */
   __ecmech_host__
   std::map<std::string, size_t>
   modelParamIndexMapBCC(const std::string_view &modelName);
   /**
    * @brief Look up parameter-count and history-array-index information for an HCP model
    * by name; see modelParamIndexMap for the returned map's keys.
    * @param modelName Model name, e.g. `"evptn_HCP_A"`.
    * @return Lookup map as described in modelParamIndexMap, or an empty map if
    * `modelName` is not recognized.
    */
   __ecmech_host__
   std::map<std::string, size_t>
   modelParamIndexMapHCP(const std::string_view &modelName);

}
