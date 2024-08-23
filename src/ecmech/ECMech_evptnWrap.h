// -*-c++-*-

#ifndef ecmech_evptnWrap_include
#define ecmech_evptnWrap_include

#include "ECMech_core.h"
#include "RAJA/RAJA.hpp"

#if defined(__ecmech_host_only__)
#include <sstream>
#endif

#include "ECMech_matModelBase.h"
#include "evptn/ECMech_evptnSngl.h"

#include "ECMech_unused.h"

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
            static constexpr int nParams =
               2 + 1 + // density0, cvav, tolerance
               SlipGeom::nParams + Kinetics::nParams + ThermoElastN::nParams + nParamsEOS;

            // constructor
            __ecmech_host__
            matModel()
               : matModelBase(),
               m_kinetics(SlipGeom::nslip)
            {
               // Should the tangent stiff matrix be included in these stride calculations?
               m_strides[istride_def_rate] = ecmech::nsvp;
               m_strides[istride_spin_v] = ecmech::ndim;
               m_strides[istride_vol_ratio] = ecmech::nvr;
               m_strides[istride_int_eng] = ecmech::ne;
               m_strides[istride_stress] = ecmech::nsvp;
               m_strides[istride_history] = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
               m_strides[istride_temp_k] = 1;
               m_strides[istride_sdd ] = ecmech::nsdd;
            };

            // constructor
            __ecmech_host__
            matModel(const unsigned int* strides, const unsigned int stride_len)
               : matModelBase(),
               m_kinetics(SlipGeom::nslip)
            {
               unsigned int nhist = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;

               if (stride_len != ecmech::nstride) {
#if defined(__ecmech_host_only__)
                  // the order here needs to be consistent with ISTRIDE_* macros in ECMECH_const.h
                  std::ostringstream os;
                  os << "Stride vector needs to have a size of " << ecmech::nstride << " with strides of at least: " <<
                     ecmech::nsvp << ", " << ecmech::ndim << ", " << ecmech::nvr << ", " <<
                     ecmech::ne << ", " << ecmech::nsvp << ", " << nhist << ", 1, " << ecmech::nsdd
                  ;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "Stride vector is the wrong size");
#endif
               }
               // Need to make sure all of the strides provided at least make sense
               if (strides[istride_def_rate] < ecmech::nsvp) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_def_rate] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_spin_v] < ecmech::ndim) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_spin_v] should have at least a length of: " << ecmech::ndim;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_vol_ratio] < ecmech::nvr) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_int_eng] should have at least a length of: " << ecmech::nvr;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_int_eng] < ecmech::ne) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_int_eng] should have at least a length of: " << ecmech::ne;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_stress] < ecmech::nsvp) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_stress] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_history] < nhist) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_history] should have at least a length of: " << nhist;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_temp_k] < 1) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_temp_k] should have at least a length of: " << 1;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_sdd] < ecmech::nsdd) {
#if defined(__ecmech_host_only__)
                  std::ostringstream os;
                  os << "strides[istride_sdd] should have at least a length of: " << ecmech::nsdd;
                  ECMECH_FAIL(__func__, os.str().c_str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               for (unsigned int i = 0; i < stride_len; i++) {
                  m_strides[i] = strides[i];
               }
            };

            // deconstructor
            __ecmech_host__
            ~matModel(){}

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

               std::vector<double>::const_iterator parsIt = pars.begin();

               m_density0 = *parsIt; ++parsIt;
               m_cvav = *parsIt; ++parsIt;

               m_tolerance = *parsIt; ++parsIt;

               {
                  const std::vector<double> paramsThese(parsIt, parsIt + SlipGeom::nParams);
                  m_slipGeom.setParams(paramsThese); parsIt += SlipGeom::nParams;
               }
               {
                  const std::vector<double> paramsThese(parsIt, parsIt + ThermoElastN::nParams);
                  m_elastN.setParams(paramsThese); parsIt += ThermoElastN::nParams;
               }
               {
                  const std::vector<double> paramsThese(parsIt, parsIt + Kinetics::nParams);
                  m_kinetics.setParams(paramsThese); parsIt += Kinetics::nParams;
               }
               {
                  double bulk_modulus = m_elastN.getBulkMod();
                  std::vector<double> paramsThese(EosModel::nParams);
                  paramsThese[0] = m_density0;
                  paramsThese[1] = bulk_modulus;
                  paramsThese[2] = m_cvav;
                  std::copy(parsIt, parsIt + nParamsEOS, paramsThese.begin() + nParamsEOSHave);

                  m_eosModel.setParams(paramsThese); parsIt += nParamsEOS;

                  {
                     double rel_vol_min, rel_vol_max;
                     m_eosModel.getInfo(rel_vol_min, rel_vol_max, m_energy0, m_rel_vol0);
                  }
               }

               int iParam = parsIt - pars.begin();
               if (iParam != nParams) {
                  ECMECH_FAIL(__func__, "wrong number of params");
               }

#if defined(__ecmech_host_only__)
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
               m_kinetics.getHistInfo(m_rhvNames, m_rhvVals, m_rhvPlot, m_rhvState);
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
#endif
            };

            using matModelBase::getParams;
            __ecmech_host__
            void getParams(std::vector<int>& opts,
                           std::vector<double>& pars,
                           std::vector<std::string>& strs) const override final
            {
               opts = m_opts;
               pars = m_pars;
               strs = m_strs;
            };

            using matModelBase::getResponseECM;
            __ecmech_host__
            void getResponseECM(const double & dt,
                                const double * defRateV,
                                const double * spinV,
                                const double * volRatioV,
                                double * internal_energyV,
                                double * cauchy_stress_dev6_pressureV,
                                double * histV,
                                double * temp_kV,
                                double * sddV,
                                double * mtanSDV,
                                const int& nPassed) const override final
            {
               if (!m_complete) {
                  ECMECH_FAIL(__func__, "not complete");
               }

               RAJA::RangeSegment default_range(0, nPassed);
               // All of the stride lengths are constant within this function
               const unsigned int def_rate_stride = m_strides[istride_def_rate];
               const unsigned int spin_v_stride = m_strides[istride_spin_v];
               const unsigned int vol_ratio_stride = m_strides[istride_vol_ratio];
               const unsigned int int_eng_stride = m_strides[istride_int_eng];
               const unsigned int stress_stride = m_strides[istride_stress];
               const unsigned int history_stride = m_strides[istride_history];
               const unsigned int temp_k_stride = m_strides[istride_temp_k];
               const unsigned int sdd_stride = m_strides[istride_sdd];

               switch (m_accel) {
#if defined(RAJA_ENABLE_OPENMP)
                  case ECM_EXEC_STRAT_OPENMP:
                  {
                     RAJA::ReduceSum<RAJA::omp_reduce_ordered, int> status_all(0);
                     RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i) {
                        double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                        const bool status = 
                        getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                        (m_slipGeom, m_kinetics, m_elastN, m_eosModel,
                        dt,
                        m_tolerance,
                        &defRateV[def_rate_stride * i],
                        &spinV[spin_v_stride * i],
                        &volRatioV[vol_ratio_stride * i],
                        &internal_energyV[int_eng_stride * i],
                        &cauchy_stress_dev6_pressureV[stress_stride * i],
                        &histV[history_stride * i],
                        temp_kV[temp_k_stride * i],
                        &sddV[sdd_stride * i],
                        mtanSDThis,
                        m_outputLevel);

                        status_all += (int) (!status);
                        if (!status) {
                         histV[history_stride * i + iHistA_nFEval] *= -1;
                        }
                     });

                     if (status_all.get() > 0) {
                        getResponseRetry(dt, defRateV, spinV, volRatioV, internal_energyV,
                                         cauchy_stress_dev6_pressureV, histV, temp_kV, sddV, mtanSDV, nPassed);
                     }

                     break;
                  }
#endif
#if defined(RAJA_ENABLE_CUDA) || defined(RAJA_ENABLE_HIP)
                  case ECM_EXEC_STRAT_GPU:
                  {
#if defined(RAJA_ENABLE_CUDA)
                     using gpu_reduce = RAJA::cuda_reduce;
                     using gpu_policy = RAJA::cuda_exec<RAJA_CUDA_THREADS>;
#else
                     using gpu_reduce = RAJA::hip_reduce;
                     using gpu_policy = RAJA::hip_exec<RAJA_HIP_THREADS>;
#endif
                     RAJA::ReduceSum<gpu_reduce, int> status_all(0);
                     RAJA::forall<gpu_policy>(default_range, [ =
#if defined(ECMECH_NON_CORAL1_MACHINE)|| defined(RAJA_ENABLE_HIP)
                      , m_slipGeom=this->m_slipGeom, m_kinetics=this->m_kinetics, m_elastN=this->m_elastN, m_eosModel=this->m_eosModel, m_tolerance=this->m_tolerance, m_outputLevel=this->m_outputLevel
#endif
                      ] RAJA_DEVICE(int i) {
                        double *mtanSDThis = (mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr);
                        const bool status =
		                  getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                        (m_slipGeom, m_kinetics, m_elastN, m_eosModel,
                           dt,
                           m_tolerance,
                           &defRateV[def_rate_stride * i],
                           &spinV[spin_v_stride * i],
                           &volRatioV[vol_ratio_stride * i],
                           &internal_energyV[int_eng_stride * i],
                           &cauchy_stress_dev6_pressureV[stress_stride * i],
                           &histV[history_stride * i],
                           temp_kV[temp_k_stride * i],
                           &sddV[sdd_stride * i],
                           mtanSDThis,
                           m_outputLevel);

                        status_all += (int) (!status);
                        if (!status) {
                           histV[history_stride * i + iHistA_nFEval] *= -1;
                        }
                     });

                     if (status_all.get() > 0) {
                        getResponseRetry(dt, defRateV, spinV, volRatioV, internal_energyV,
                                         cauchy_stress_dev6_pressureV, histV, temp_kV, sddV, mtanSDV, nPassed);
                     }

                     break;
                  }
#endif
                  case ECM_EXEC_STRAT_CPU:
                  default: // fall through to CPU if other options are not available
                  {
                     RAJA::ReduceSum<RAJA::seq_reduce, int> status_all(0);
                     RAJA::forall<RAJA::seq_exec>(default_range, [ = ] (int i) {
                        double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                        const bool status = 
                        getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                           (m_slipGeom, m_kinetics, m_elastN, m_eosModel,
                              dt,
                              m_tolerance,
                              &defRateV[def_rate_stride * i],
                              &spinV[spin_v_stride * i],
                              &volRatioV[vol_ratio_stride * i],
                              &internal_energyV[int_eng_stride * i],
                              &cauchy_stress_dev6_pressureV[stress_stride * i],
                              &histV[history_stride * i],
                              temp_kV[temp_k_stride * i],
                              &sddV[sdd_stride * i],
                              mtanSDThis,
                              m_outputLevel);

                        status_all += (int) (!status);
                        if (!status) {
                           histV[history_stride * i + iHistA_nFEval] *= -1;
                        }
                     });

                     if (status_all.get() > 0) {
                        getResponseRetry(dt, defRateV, spinV, volRatioV, internal_energyV,
                                         cauchy_stress_dev6_pressureV, histV, temp_kV, sddV, mtanSDV, nPassed);
                     }

                     break;
                  }
               } // switch _accel
            }; // End of getResponse

            __ecmech_host__
            inline
            void getResponseRetry( const double & UNUSED_EXTRA(dt),
                                   const double * UNUSED_EXTRA(defRateV),
                                   const double * UNUSED_EXTRA(spinV),
                                   const double * UNUSED_EXTRA(volRatioV),
                                   double * UNUSED_EXTRA(internal_energyV),
                                   double * UNUSED_EXTRA(cauchy_stress_dev6_pressureV),
                                   double * UNUSED_EXTRA(histV),
                                   double * UNUSED_EXTRA(temp_kV),
                                   double * UNUSED_EXTRA(sddV),
                                   double * UNUSED_EXTRA(mtanSDV),
                                   const int& UNUSED_EXTRA(nPassed)
                                 ) const
            {
#if defined(ECMECH_EXTRA_SOLVERS)
               if (!m_complete) {
                  ECMECH_FAIL(__func__, "not complete");
               }

               RAJA::RangeSegment default_range(0, nPassed);
               // All of the stride lengths are constant within this function
               const unsigned int def_rate_stride = m_strides[istride_def_rate];
               const unsigned int spin_v_stride = m_strides[istride_spin_v];
               const unsigned int vol_ratio_stride = m_strides[istride_vol_ratio];
               const unsigned int int_eng_stride = m_strides[istride_int_eng];
               const unsigned int stress_stride = m_strides[istride_stress];
               const unsigned int history_stride = m_strides[istride_history];
               const unsigned int temp_k_stride = m_strides[istride_temp_k];
               const unsigned int sdd_stride = m_strides[istride_sdd];

               switch (m_accel) {
#if defined(RAJA_ENABLE_OPENMP)
                  case ECM_EXEC_STRAT_OPENMP:
                  {
                     RAJA::ReduceSum<RAJA::omp_reduce_ordered, int> status_all(0);
                     RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i) {
                        if (histV[history_stride * i + iHistA_nFEval] < 0) { // skip elements that were successful
                        double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                        bool status = getResponseNRSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                           (m_slipGeom, m_kinetics, m_elastN, m_eosModel,
                           dt,
                           m_tolerance,
                           &defRateV[def_rate_stride * i],
                           &spinV[spin_v_stride * i],
                           &volRatioV[vol_ratio_stride * i],
                           &internal_energyV[int_eng_stride * i],
                           &cauchy_stress_dev6_pressureV[stress_stride * i],
                           &histV[history_stride * i],
                           temp_kV[temp_k_stride * i],
                           &sddV[sdd_stride * i],
                           mtanSDThis,
                           m_outputLevel);

                        status_all += (int) (!status);
                        if (!status) {
                           histV[history_stride * i + iHistA_nFEval] *= -1;
                        }
                        }
                     });

                     if (status_all.get() > 0) {
                         ECMECH_FAIL(__func__, "Back-up solvers failed to converge for at least one point!");
                     }

                     break;
                  }
#endif
#if defined(RAJA_ENABLE_CUDA) || defined(RAJA_ENABLE_HIP)
                  case ECM_EXEC_STRAT_GPU:
                  {
#if defined(RAJA_ENABLE_CUDA)
                     using gpu_reduce = RAJA::cuda_reduce;
                     using gpu_policy = RAJA::cuda_exec<RAJA_CUDA_THREADS>;
#else
                     using gpu_reduce = RAJA::hip_reduce;
                     using gpu_policy = RAJA::hip_exec<RAJA_HIP_THREADS>;
#endif
                     RAJA::ReduceSum<gpu_reduce, int> status_all(0);
                     RAJA::forall<gpu_policy>(default_range, [ =
#if defined(ECMECH_NON_CORAL1_MACHINE)|| defined(RAJA_ENABLE_HIP)
                      , m_slipGeom=this->m_slipGeom, m_kinetics=this->m_kinetics, m_elastN=this->m_elastN, m_eosModel=this->m_eosModel, m_tolerance=this->m_tolerance, m_outputLevel=this->m_outputLevel
#endif
                      ] RAJA_DEVICE(int i) {
                        if (histV[history_stride * i + iHistA_nFEval] < 0) { // skip elements that were successful
                        double *mtanSDThis = (mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr);
                        bool status = getResponseNRSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                           (m_slipGeom, m_kinetics, m_elastN, m_eosModel,
                           dt,
                           m_tolerance,
                           &defRateV[def_rate_stride * i],
                           &spinV[spin_v_stride * i],
                           &volRatioV[vol_ratio_stride * i],
                           &internal_energyV[int_eng_stride * i],
                           &cauchy_stress_dev6_pressureV[stress_stride * i],
                           &histV[history_stride * i],
                           temp_kV[temp_k_stride * i],
                           &sddV[sdd_stride * i],
                           mtanSDThis,
                           m_outputLevel);

                        status_all += (int) (!status);
                        if (!status) {
                           histV[history_stride * i + iHistA_nFEval] *= -1;
                        }
                        }
                     });

                     if (status_all.get() > 0) {
                         ECMECH_FAIL(__func__, "Back-up solvers failed to converge for at least one point!");
                     }

                     break;
                  }
#endif
                  case ECM_EXEC_STRAT_CPU:
                  default: // fall through to CPU if other options are not available
                  {
                     RAJA::ReduceSum<RAJA::seq_reduce, int> status_all(0);
                     RAJA::forall<RAJA::seq_exec>(default_range, [ = ] (int i) {
                        if (histV[history_stride * i + iHistA_nFEval] < 0) { // skip elements that were successful
                        double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                        bool status = getResponseNRSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                           (m_slipGeom, m_kinetics, m_elastN, m_eosModel,
                           dt,
                           m_tolerance,
                           &defRateV[def_rate_stride * i],
                           &spinV[spin_v_stride * i],
                           &volRatioV[vol_ratio_stride * i],
                           &internal_energyV[int_eng_stride * i],
                           &cauchy_stress_dev6_pressureV[stress_stride * i],
                           &histV[history_stride * i],
                           temp_kV[temp_k_stride * i],
                           &sddV[sdd_stride * i],
                           mtanSDThis,
                           m_outputLevel);

                        status_all += (int) (!status);
                        if (!status) {
                           histV[history_stride * i + iHistA_nFEval] *= -1;
                        }

                        }
                     });

                     if (status_all.get() > 0) {
                         ECMECH_FAIL(__func__, "Back-up solvers failed to converge for at least one point!");
                     }

                     break;
                  }
               } // switch _accel
#else
               ECMECH_FAIL(__func__, "Solver failed to converge for at least one point!");
#endif
            } // End of getResponseRetry

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
            };

            __ecmech_host__
            int getNumHist( ) const override final {
               return numHist;
            };

            __ecmech_host__
            void complete( ) override final
            {
               m_bulkRef = m_eosModel.getBulkRef();
               m_complete = true;
            };

            // Constant getter functions to return the underlying templated classes.
            // Uses for these could be for example to compute the sample D^p tensor
            // using the symmetric schmid tensor from the SlipGeom class.
            //
            // Note: Stability of the underlying templated class API's is not
            // guaranteed, so breaking changes can occur from point release to
            // point release.
            const SlipGeom & getSlipGeom() const { return m_slipGeom; }

            const Kinetics & getemp_kinetics() const { return m_kinetics; }

            const ThermoElastN & getThermoElastN() const { return m_elastN; }

            const EosModel & getEosModel() const { return m_eosModel; }

         private:

            SlipGeom m_slipGeom;
            Kinetics m_kinetics;
            ThermoElastN m_elastN;
            EosModel m_eosModel;

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
