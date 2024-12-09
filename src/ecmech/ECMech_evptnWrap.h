// -*-c++-*-

#ifndef ecmech_evptnWrap_include
#define ecmech_evptnWrap_include

#include "ECMech_core.h"
#include "SNLS_config.h"
#include "SNLS_device_forall.h"
#include "RAJA/RAJA.hpp"

#if defined(SNLS_RAJA_PORT_SUITE)
#include "SNLS_memory_manager.h"
#include <chai/managed_ptr.hpp>
#endif

#if defined(__ecmech_host_only__)
#include <sstream>
#endif

#include "ECMech_matModelBase.h"
#include "evptn/ECMech_evptnSngl.h"

#include "ECMech_unused.h"

namespace ecmech {
namespace internal {
#if !defined(SNLS_RAJA_PORT_SUITE)
   // Creating a pseudo-chai::managed_ptr<T> class that implements the same set of functions we need from Chai
   template<class T>
   class PseudoChaiManagedPtr {
      public:
         PseudoChaiManagedPtr() = default;
         template<typename ...Args>
         PseudoChaiManagedPtr(const double* const params, Args... args) : m_val(new T(params, std::forward<Args...>(args...))) {}
         PseudoChaiManagedPtr(const double* const params) : m_val(new T(params)) {}
         PseudoChaiManagedPtr(const PseudoChaiManagedPtr&) = default;
         ~PseudoChaiManagedPtr() = default;

         __ecmech_hdev__
         inline T& operator*() const { return *m_val; }
         __ecmech_hdev__
         inline void free() const { if (m_val) { delete m_val; } }

      public:
         T* m_val = nullptr;
   };

   template<class T>
   using pcmptr = PseudoChaiManagedPtr<T>;
#endif

   template<class T, typename ...Args>
   __ecmech_host__
   auto make_class_factory(const std::vector<double>& params, Args... args) {
#if defined(SNLS_RAJA_PORT_SUITE)
      auto mm = snls::memoryManager::getInstance();
      auto mvec = mm.allocManagedArray<double>(params.size());

      auto mvec_data = mvec.data(chai::ExecutionSpace::CPU);
      for (size_t i = 0; i < params.size(); i++ ) {
         mvec_data[i] = params[i];
      }

      return chai::make_managed<T>(chai::unpack(mvec), std::forward<Args...>(args...));
#else
      return pcmptr<T>(params.data(), std::forward<Args...>(args...));
#endif
   }

   template<class T>
   __ecmech_host__
   auto make_class_factory(const std::vector<double>& params) {
#if defined(SNLS_RAJA_PORT_SUITE)
      auto mm = snls::memoryManager::getInstance();
      auto mvec = mm.allocManagedArray<double>(params.size());

      auto mvec_data = mvec.data(chai::ExecutionSpace::CPU);
      for (size_t i = 0; i < params.size(); i++ ) {
         mvec_data[i] = params[i];
      }

      return chai::make_managed<T>(chai::unpack(mvec));
#else
      return pcmptr<T>(params.data());
#endif
   }

}
}

namespace ecmech {
   namespace evptn {
      //
      // template on the specifics of the crystal model ;
      // but have a base class so that the templating can stop here
      //
      template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
      class matModel : public matModelBase
      {
         public:

            static constexpr int iHistLbGdot = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;
            static constexpr int numHist = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
            static constexpr int nH = Kinetics::nH;
            static constexpr int nslip = SlipGeom::nslip;

            static constexpr int nParamsEOSHave = 3; // number that get from 'elsewhere' // these are assumed to go in first
            static constexpr int nParamsEOS = EosModel::nParams - nParamsEOSHave;
            static constexpr int nParamsSlipGeom = SlipGeom::nParams;
            static constexpr int nParamsKinetics = Kinetics::nParams;
            static constexpr int nParamsThermoElastN = ThermoElastN::nParams;
            static constexpr int nParams =
               2 + 1 + // density0, cvav, tolerance
               nParamsSlipGeom + nParamsKinetics + nParamsThermoElastN + nParamsEOS;

            // constructor
            __ecmech_host__
            matModel()
               : matModelBase()
            {
               // Should the tangent stiff matrix be included in these stride calculations?
               m_strides[istride_def_rate] = ecmech::nsvp;
               m_strides[istride_spin_v] = ecmech::ndim;
               m_strides[istride_vol_ratio] = ecmech::nvr;
               m_strides[istride_int_eng] = ecmech::ne;
               m_strides[istride_stress] = ecmech::nsvp;
               m_strides[istride_history] = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
               m_strides[istride_tkelv] = 1;
               m_strides[istride_sdd ] = ecmech::nsdd;
            }

            // constructor
            __ecmech_host__
            matModel(std::vector<size_t> strides)
               : matModelBase()
            {
               updateStrides(strides);
            }

            // deconstructor
            __ecmech_host__
            ~matModel()
            {
               m_slipGeom.free();
               m_kinetics.free();
               m_eosModel.free();
               m_elastN.free();
            }

            __ecmech_host__
            virtual void
            updateStrides(std::vector<size_t> strides) override final {
               if (m_complete) {
                  ECMECH_FAIL(__func__, "updateStrides should only be called before object is completed");
               }
               if (strides.size() != ecmech::nstride) {
                  // the order here needs to be consistent with ISTRIDE_* macros in ECMECH_const.h
                  std::ostringstream os;
                  os << "Stride vector needs to have a size of " << ecmech::nstride << " with strides of at least: " <<
                     ecmech::nsvp << ", " << ecmech::ndim << ", " << ecmech::nvr << ", " <<
                     ecmech::ne << ", " << ecmech::nsvp << ", " << numHist << ", 1, " << ecmech::nsdd
                  ;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               // Need to make sure all of the strides provided at least make sense
               if (strides[istride_def_rate] < ecmech::nsvp) {
                  std::ostringstream os;
                  os << "strides[istride_def_rate] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_spin_v] < ecmech::ndim) {
                  std::ostringstream os;
                  os << "strides[istride_spin_v] should have at least a length of: " << ecmech::ndim;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_vol_ratio] < ecmech::nvr) {
                  std::ostringstream os;
                  os << "strides[istride_int_eng] should have at least a length of: " << ecmech::nvr;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_int_eng] < ecmech::ne) {
                  std::ostringstream os;
                  os << "strides[istride_int_eng] should have at least a length of: " << ecmech::ne;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_stress] < ecmech::nsvp) {
                  std::ostringstream os;
                  os << "strides[istride_stress] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_history] < numHist) {
                  std::ostringstream os;
                  os << "strides[istride_history] should have at least a length of: " << numHist;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_tkelv] < 1) {
                  std::ostringstream os;
                  os << "strides[istride_tkelv] should have at least a length of: " << 1;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               if (strides[istride_sdd] < ecmech::nsdd) {
                  std::ostringstream os;
                  os << "strides[istride_sdd] should have at least a length of: " << ecmech::nsdd;
                  ECMECH_FAIL(__func__, os.str().c_str());
               }
               for (unsigned int i = 0; i < strides.size(); i++) {
                  m_strides[i] = strides[i];
               }
            }

            using matModelBase::initFromParams;
            __ecmech_host__
            void initFromParams(const std::vector<int>& opts,
                                const std::vector<double>& pars,
                                const std::vector<std::string>& strs,
                                void* /*callBackVoid*/ = nullptr
                                ) override final
            {
               // keep parameters for later
               m_opts = opts;
               m_pars = pars;
               m_strs = strs;

               if (pars.size() != (unsigned int) nParams) {
                  ECMECH_FAIL(__func__, "wrong number of pars");
               }
               if (opts.size() != 0) {
                  ECMECH_FAIL(__func__, "wrong number of opts");
               }
               if (strs.size() > 1) {
                  // strs[0] is optionally a name -- see makeMatModel
                  ECMECH_FAIL(__func__, "wrong number of strs");
               }

               // Want to make sure we free up any old memory before setting parameters just in-case we had a model around already...
               m_slipGeom.free();
               m_kinetics.free();
               m_eosModel.free();
               m_elastN.free();

               std::vector<double>::const_iterator parsIt = pars.begin();

               m_density0 = *parsIt; ++parsIt;
               m_cvav = *parsIt; ++parsIt;

               m_tolerance = *parsIt; ++parsIt;

               {
                  const std::vector<double> paramsThese(parsIt, parsIt + SlipGeom::nParams);
                  m_slipGeom = internal::make_class_factory<SlipGeom>(paramsThese);
                  parsIt += SlipGeom::nParams;
                  // m_slipGeom.setParams(paramsThese); parsIt += SlipGeom::nParams;
               }
               {
                  const std::vector<double> paramsThese(parsIt, parsIt + ThermoElastN::nParams);
                  m_elastN = internal::make_class_factory<ThermoElastN>(paramsThese);
                  parsIt += ThermoElastN::nParams;
                  // m_elastN.setParams(paramsThese); parsIt += ThermoElastN::nParams;
               }
               {
                  const std::vector<double> paramsThese(parsIt, parsIt + Kinetics::nParams);
                  m_kinetics = internal::make_class_factory<Kinetics>(paramsThese, SlipGeom::nslip);
                  parsIt += Kinetics::nParams;
                  // m_kinetics.setParams(paramsThese); parsIt += Kinetics::nParams;
               }
               {
                  double bulk_modulus = (*m_elastN).getBulkMod();
                  std::vector<double> paramsThese(EosModel::nParams);
                  paramsThese[0] = m_density0;
                  paramsThese[1] = bulk_modulus;
                  paramsThese[2] = m_cvav;
                  std::copy(parsIt, parsIt + nParamsEOS, paramsThese.begin() + nParamsEOSHave);

                  m_eosModel = internal::make_class_factory<EosModel>(paramsThese);
                  parsIt += nParamsEOS;
                  // m_eosModel.setParams(paramsThese); parsIt += nParamsEOS;

                  {
                     double rel_vol_min, rel_vol_max;
                     (*m_eosModel).getInfo(rel_vol_min, rel_vol_max, m_energy0, m_rel_vol0);
                  }
               }

               int iParam = parsIt - pars.begin();
               if (iParam != nParams) {
                  ECMECH_FAIL(__func__, "wrong number of params");
               }
               //////////////////////////////

               m_rhvNames.clear();
               m_rhvVals.clear();
               m_rhvPlot.clear();
               m_rhvState.clear();

#if defined(ECMECH_USE_DPEFF)
               m_rhvNames.push_back("dplas_eff"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrateEff
               m_rhvNames.push_back("eps"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrEff
#else
               m_rhvNames.push_back("shrate_eff"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrateEff
               m_rhvNames.push_back("shr_eff"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrEff
#endif
               m_rhvNames.push_back("flow_str"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(false); // iHistA_flowStr
               m_rhvNames.push_back("n_feval"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(false); // iHistA_nFEval
               // numHistAux
               //
               for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                  std::ostringstream os;
                  os << "t" << iTvec + 1;
                  m_rhvNames.push_back(os.str()); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true);
               }

               //
               {
                  double qVal = 1.0;
                  for (int iQ = 0; iQ < ecmech::qdim; ++iQ) {
                     std::ostringstream os;
                     os << "quat_" << iQ + 1;
                     m_rhvNames.push_back(os.str()); m_rhvVals.push_back(qVal); m_rhvPlot.push_back(true); m_rhvState.push_back(true);
                     qVal = 0.0;
                  }
               }
               //
               (*m_kinetics).getHistInfo(m_rhvNames, m_rhvVals, m_rhvPlot, m_rhvState);
               //
               for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                  std::ostringstream os;
                  os << "shrate_" << iSlip + 1;
                  m_rhvNames.push_back(os.str()); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true);
               }

               //
               if (m_rhvNames.size() != numHist) {
                  ECMECH_FAIL(__func__, "mismatch in numHist");
               }
            }

            using matModelBase::getParams;
            __ecmech_host__
            void getParams(std::vector<int>& opts,
                           std::vector<double>& pars,
                           std::vector<std::string>& strs) const override final
            {
               opts = m_opts;
               pars = m_pars;
               strs = m_strs;
            }

            using matModelBase::getResponseECM;
            __ecmech_host__
            void getResponseECM(const double & dt,
                                const double * def_rate_d6vV,
                                const double * spin_vecV,
                                const double * rel_vol_ratiosV,
                                double * internal_energyV,
                                double * cauchy_stress_d6pV,
                                double * histV,
                                double * tkelvV,
                                double * sddV,
                                double * mtanSDV,
                                const int& nPassed) const override final
            {
               if (!m_complete) {
                  ECMECH_FAIL(__func__, "not complete");
               }

               // All of the stride lengths are constant within this function
               const unsigned int def_rate_stride = m_strides[istride_def_rate];
               const unsigned int spin_v_stride = m_strides[istride_spin_v];
               const unsigned int vol_ratio_stride = m_strides[istride_vol_ratio];
               const unsigned int int_eng_stride = m_strides[istride_int_eng];
               const unsigned int stress_stride = m_strides[istride_stress];
               const unsigned int history_stride = m_strides[istride_history];
               const unsigned int tkelv_stride = m_strides[istride_tkelv];
               const unsigned int sdd_stride = m_strides[istride_sdd];

               const auto slipGeom = m_slipGeom;
               const auto kinetics = m_kinetics;
               const auto elastN = m_elastN;
               const auto eosModel = m_eosModel;
               const auto tolerance = m_tolerance;
               const auto outputLevel = m_outputLevel;

               snls::forall<ECMECH_GPU_THREADS>(0, nPassed, [=]
                  __ecmech_hdev__
                  (int i)
               {
                  double *mtanSDThis = (mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr);
                  auto get_response = [=] (const SlipGeom& slip_geom) -> bool {
                     return getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                     (slip_geom, *kinetics, *elastN, *eosModel,
                        dt,
                        tolerance,
                        &def_rate_d6vV[def_rate_stride * i],
                        &spin_vecV[spin_v_stride * i],
                        &rel_vol_ratiosV[vol_ratio_stride * i],
                        &internal_energyV[int_eng_stride * i],
                        &cauchy_stress_d6pV[stress_stride * i],
                        &histV[history_stride * i],
                        tkelvV[tkelv_stride * i],
                        &sddV[sdd_stride * i],
                        mtanSDThis,
                        outputLevel);
                  };
                  bool status;
                  // Thanks to NVCC being difficult we have to create an unnecessary temp variable just so we can use our
                  // lambda expression in an if constexpr...
                  auto slipGeom_tmp = slipGeom;
                  // If we have dynamic slip systems then we need to create
                  // a thread local slip geometery class or else we might run into
                  // race condition issues...
                  if constexpr (SlipGeom::dynamic) {
                     SlipGeom slip_geom = *slipGeom_tmp;
                     status = get_response(slip_geom);
                  } else {
                     status = get_response(*slipGeom_tmp);
                  }
                  if (!status) {
                     histV[history_stride * i + iHistA_nFEval] *= -1;
                  }
               });

               if (this->reduceStatus(histV, nPassed)) {
                  getResponseRetry(dt, def_rate_d6vV, spin_vecV, rel_vol_ratiosV, internal_energyV,
                                   cauchy_stress_d6pV, histV, tkelvV, sddV, mtanSDV, nPassed);
               }
            }// End of getResponse

            __ecmech_host__
            inline
            void getResponseRetry( const double & UNUSED_EXTRA(dt),
                                   const double * UNUSED_EXTRA(def_rate_d6vV),
                                   const double * UNUSED_EXTRA(spin_vecV),
                                   const double * UNUSED_EXTRA(rel_vol_ratiosV),
                                   double * UNUSED_EXTRA(internal_energyV),
                                   double * UNUSED_EXTRA(cauchy_stress_d6pV),
                                   double * UNUSED_EXTRA(histV),
                                   double * UNUSED_EXTRA(tkelvV),
                                   double * UNUSED_EXTRA(sddV),
                                   double * UNUSED_EXTRA(mtanSDV),
                                   const int& UNUSED_EXTRA(nPassed)
                                 ) const
            {
#if defined(ECMECH_EXTRA_SOLVERS)
               if (!m_complete) {
                  ECMECH_FAIL(__func__, "not complete");
               }

               // All of the stride lengths are constant within this function
               const unsigned int def_rate_stride = m_strides[istride_def_rate];
               const unsigned int spin_v_stride = m_strides[istride_spin_v];
               const unsigned int vol_ratio_stride = m_strides[istride_vol_ratio];
               const unsigned int int_eng_stride = m_strides[istride_int_eng];
               const unsigned int stress_stride = m_strides[istride_stress];
               const unsigned int history_stride = m_strides[istride_history];
               const unsigned int tkelv_stride = m_strides[istride_tkelv];
               const unsigned int sdd_stride = m_strides[istride_sdd];

               const auto slipGeom = m_slipGeom;
               const auto kinetics = m_kinetics;
               const auto elastN = m_elastN;
               const auto eosModel = m_eosModel;
               const auto tolerance = m_tolerance;
               const auto outputLevel = m_outputLevel;

               snls::forall(0, nPassed, [=]
                  __ecmech_hdev__
                  (int i)
               {
                  // skip elements that were successful
                  if (histV[history_stride * i + iHistA_nFEval] >= 0) return;

                  double *mtanSDThis = (mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr);
                  auto get_response = [=] (const SlipGeom& slip_geom) -> bool {
                     return getResponseNRSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                     (slip_geom, *kinetics, *elastN, *eosModel,
                        dt,
                        tolerance,
                        &def_rate_d6vV[def_rate_stride * i],
                        &spin_vecV[spin_v_stride * i],
                        &rel_vol_ratiosV[vol_ratio_stride * i],
                        &internal_energyV[int_eng_stride * i],
                        &cauchy_stress_d6pV[stress_stride * i],
                        &histV[history_stride * i],
                        tkelvV[tkelv_stride * i],
                        &sddV[sdd_stride * i],
                        mtanSDThis,
                        outputLevel);
                  };
                  bool status;
                  // Thanks to NVCC being difficult we have to create an unnecessary temp variable just so we can use our
                  // lambda expression in an if constexpr...
                  auto slipGeom_tmp = slipGeom;
                  // If we have dynamic slip systems then we need to create
                  // a thread local slip geometery class or else we might run into
                  // race condition issues...
                  if constexpr (SlipGeom::dynamic) {
                     SlipGeom slip_geom = *slipGeom_tmp;
                     status = get_response(slip_geom);
                  } else {
                     status = get_response(*slipGeom_tmp);
                  }
                  if (!status) {
                     histV[history_stride * i + iHistA_nFEval] *= -1;
                  }
               });

               if (this->reduceStatus(histV, nPassed)) {
                  ECMECH_FAIL(__func__, "Back-up solvers failed to converge for at least one point!");
               }
#else
               ECMECH_FAIL(__func__, "Solver failed to converge for at least one point!");
#endif
            } // End of getResponseRetry

            __ecmech_host__
            bool reduceStatus(const double* const histV,
                              const int& nPassed) const
            {
               RAJA::RangeSegment default_range(0, nPassed);
               const unsigned int history_stride = m_strides[istride_history];
               switch (m_accel) {
#if defined(RAJA_ENABLE_OPENMP)
                  case ECM_EXEC_STRAT_OPENMP:
                  {
                     RAJA::ReduceSum<RAJA::omp_reduce_ordered, int> status_all(0);
                     RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i) {
                        status_all += (histV[history_stride * i + iHistA_nFEval] < 0);
                     });
                     return status_all.get() > 0;
                     break;
                  }
#endif
#if defined(RAJA_ENABLE_CUDA) || defined(RAJA_ENABLE_HIP)
                  case ECM_EXEC_STRAT_GPU:
                  {
#if defined(RAJA_ENABLE_CUDA)
                     using gpu_reduce = RAJA::cuda_reduce;
                     using gpu_policy = RAJA::cuda_exec<ECMECH_GPU_THREADS>;
#else
                     using gpu_reduce = RAJA::hip_reduce;
                     using gpu_policy = RAJA::hip_exec<ECMECH_GPU_THREADS>;
#endif
                     RAJA::ReduceBitOr<gpu_reduce, bool> status_all(false);
                     RAJA::forall<gpu_policy>(default_range, [=] RAJA_DEVICE(int i) {
                        status_all |= (histV[history_stride * i + iHistA_nFEval] < 0);
                     });
                     return status_all.get();
                     break;
                  }
#endif
                  case ECM_EXEC_STRAT_CPU:
                  default: // fall through to CPU if other options are not available
                  {
                     RAJA::ReduceBitOr<RAJA::seq_reduce, bool> status_all(false);
                     RAJA::forall<RAJA::seq_exec>(default_range, [ = ] (int i) {
                        status_all |= (histV[history_stride * i + iHistA_nFEval] < 0);
                     });
                     return status_all.get();
                     break;
                  }
               } // switch _accel
            }

            using matModelBase::getHistInfo;
            __ecmech_host__
            void getHistInfo(std::vector<std::string> & names,
                             std::vector<double>      & vals,
                             std::vector<bool>        & plot,
                             std::vector<bool>        & state) const override final {
               if (m_rhvNames.size() != numHist) {
                  ECMECH_FAIL(__func__, "have not yet set up history information");
               }
               names.resize(numHist); std::copy(m_rhvNames.begin(), m_rhvNames.end(), names.begin() );
               vals.resize(numHist); std::copy(m_rhvVals.begin(), m_rhvVals.end(), vals.begin() );
               plot.resize(numHist); std::copy(m_rhvPlot.begin(), m_rhvPlot.end(), plot.begin() );
               state.resize(numHist); std::copy(m_rhvState.begin(), m_rhvState.end(), state.begin() );
            }

            __ecmech_host__
            int getNumHist( ) const override final {
               return numHist;
            }

            __ecmech_host__
            void complete( ) override final
            {
               m_bulkRef = (*m_eosModel).getBulkRef();
               m_complete = true;
            }

            using ExecStrat = ecmech::ExecutionStrategy;

            /**
             *  @brief
             *  Set the accelerator to be used for getResponse.
             */
            __ecmech_host__
            void setExecutionStrategy(ExecStrat accel) override final  {
               switch (accel) {
#ifdef __ecmech_gpu_active__
                  case (ECM_EXEC_STRAT_GPU): {
                     m_accel = ECM_EXEC_STRAT_GPU;
                     break;
                  }
#endif
#if defined(RAJA_ENABLE_OPENMP) && defined(OPENMP_ENABLE)
                  case (ECM_EXEC_STRAT_OPENMP): {
                     m_accel = ECM_EXEC_STRAT_OPENMP;
                     break;
                  }
#endif
                  case (ECM_EXEC_STRAT_CPU):
                  default: {
                     m_accel = ECM_EXEC_STRAT_CPU;
                     break;
                  }
               }
               this->setSNLSExecutionStrategy(m_accel);
            }

            __ecmech_host__
            void setSNLSExecutionStrategy(ExecStrat accel) const
            {
               snls::Device &device = snls::Device::GetInstance();
               switch (accel) {
#ifdef __ecmech_gpu_active__
                  case (ECM_EXEC_STRAT_GPU): {
                     device.SetBackend(snls::ExecutionStrategy::GPU);
                     break;
                  }
#endif
#if defined(RAJA_ENABLE_OPENMP) && defined(OPENMP_ENABLE)
                  case (ECM_EXEC_STRAT_OPENMP): {
                     device.SetBackend(snls::ExecutionStrategy::OPENMP);
                     break;
                  }
#endif
                  case (ECM_EXEC_STRAT_CPU):
                  default: {
                     device.SetBackend(snls::ExecutionStrategy::CPU);
                     break;
                  }
               }
            }

            // Constant getter functions to return the underlying templated classes.
            // Uses for these could be for example to compute the sample D^p tensor
            // using the symmetric schmid tensor from the SlipGeom class.
            //
            // Note: Stability of the underlying templated class API's is not
            // guaranteed, so breaking changes can occur from point release to
            // point release.
            const SlipGeom & getSlipGeom() const { return *m_slipGeom; }

            const Kinetics & getKinetics() const { return *m_kinetics; }

            const ThermoElastN & getThermoElastN() const { return *m_elastN; }

            const EosModel & getEosModel() const { return *m_eosModel; }

         private:

#if defined(SNLS_RAJA_PORT_SUITE)
            chai::managed_ptr<SlipGeom> m_slipGeom;
            chai::managed_ptr<Kinetics> m_kinetics;
            chai::managed_ptr<ThermoElastN> m_elastN;
            chai::managed_ptr<EosModel> m_eosModel;
#else
            internal::pcmptr<SlipGeom> m_slipGeom;
            internal::pcmptr<Kinetics> m_kinetics;
            internal::pcmptr<ThermoElastN> m_elastN;
            internal::pcmptr<EosModel> m_eosModel;
#endif

            double m_tolerance;
            unsigned int m_strides[ecmech::nstride];

            std::vector<std::string> m_rhvNames;
            std::vector<double>      m_rhvVals;
            std::vector<bool>        m_rhvPlot;
            std::vector<bool>        m_rhvState;

            // keep initFromParams vectors as a convenience
            std::vector<int>          m_opts;
            std::vector<double>       m_pars;
            std::vector<std::string>  m_strs;
      }; // class matModel
   } // namespace evptn
} // namespace ecmech

#endif // ecmech_evptnWrap_include
