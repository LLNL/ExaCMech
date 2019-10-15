// -*-c++-*-

#ifndef ecmech_evptnWrap_include
#define ecmech_evptnWrap_include

#include "ECMech_core.h"
#include "RAJA/RAJA.hpp"

#ifdef __cuda_host_only__
#include <sstream>
#endif

#include "ECMech_matModelBase.h"
#include "ECMech_evptn.h"

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

            // this are assumed to go in first
            static constexpr int nParamsEOSHave = 3; // number that get from 'elsewhere'

            // constructor
            __ecmech_host__
            matModel()
               : matModelBase(),
               _kinetics(SlipGeom::nslip),
               _outputLevel(0)
            {
               // Should the tangent stiff matrix be included in these stride calculations?
               _strides[istride_def_rate] = ecmech::nsvp;
               _strides[istride_spin_v] = ecmech::ndim;
               _strides[istride_vol_ratio] = ecmech::nvr;
               _strides[istride_int_eng] = ecmech::ne;
               _strides[istride_stress] = ecmech::nsvp;
               _strides[istride_history] = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
               _strides[istride_tkelv] = 1;
               _strides[istride_sdd ] = ecmech::nsdd;
            };

            // constructor
            __ecmech_host__
            matModel(const uint* strides, const uint stride_len)
               : matModelBase(),
               _kinetics(SlipGeom::nslip),
               _outputLevel(0)
            {
               uint nhist = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;

               if (stride_len != 8) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "Stride vector needs to have a size of 8 with strides of at least: " <<
                     ecmech::nsvp << ", " << ecmech::ndim << ", " << ecmech::nvr << ", " <<
                     ecmech::ne << ", " << ecmech::nsvp << ", " << nhist << ", 1, " << ecmech::nsdd;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "Stride vector needs to have a size of 8 with strides");
#endif
               }
               // Need to make sure all of the strides provided at least make sense
               if (strides[istride_def_rate] < ecmech::nsvp) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_def_rate] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_spin_v] < ecmech::ndim) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_spin_v] should have at least a length of: " << ecmech::ndim;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_vol_ratio] < ecmech::nvr) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_int_eng] should have at least a length of: " << ecmech::nvr;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_int_eng] < ecmech::ne) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_int_eng] should have at least a length of: " << ecmech::ne;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_stress] < ecmech::nsvp) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_stress] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_history] < nhist) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_history] should have at least a length of: " << nhist;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_tkelv] < 1) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_tkelv] should have at least a length of: " << 1;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               if (strides[istride_sdd] < ecmech::nsdd) {
#ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[istride_sdd] should have at least a length of: " << ecmech::nsdd;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
#endif
               }
               for (uint i = 0; i < stride_len; i++) {
                  _strides[i] = strides[i];
               }
            };

            // deconstructor
            __ecmech_host__
            ~matModel(){}

            __ecmech_hdev__
            void setOutputLevel(int outputLevel) { _outputLevel = outputLevel; };

            using matModelBase::initFromParams;
            __ecmech_host__
            void initFromParams(const std::vector<int>& opts,
                                const std::vector<double>& pars,
                                const std::vector<std::string>& strs,
                                const ecmech::Accelerator& accel = ecmech::Accelerator::CPU) override final {
               const int nParamsEOS = EosModel::nParams - nParamsEOSHave;
               int nParams =
                  2 + 1 + // rho0, cvav, tolerance
                  Kinetics::nParams +
                  ThermoElastN::nParams +
                  nParamsEOS;

               if (pars.size() != (unsigned int) nParams) {
                  ECMECH_FAIL(__func__, "wrong number of pars");
               }
               if (opts.size() != 0) {
                  ECMECH_FAIL(__func__, "wrong number of opts");
               }
               if (strs.size() != 0) {
                  ECMECH_FAIL(__func__, "wrong number of strs");
               }

               _accel = accel;

               std::vector<double>::const_iterator parsIt = pars.begin();

               _rho0 = *parsIt; ++parsIt;
               _cvav = *parsIt; ++parsIt;

               _tolerance = *parsIt; ++parsIt;

               {
                  const std::vector<double> paramsThese(parsIt, parsIt + ThermoElastN::nParams);
                  _elastN.setParams(paramsThese); parsIt += ThermoElastN::nParams;
               }
               {
                  const std::vector<double> paramsThese(parsIt, parsIt + Kinetics::nParams);
                  _kinetics.setParams(paramsThese); parsIt += Kinetics::nParams;
               }
               {
                  double bulkMod = _elastN.getBulkMod();
                  std::vector<double> paramsThese(EosModel::nParams);
                  paramsThese[0] = _rho0;
                  paramsThese[1] = bulkMod;
                  paramsThese[2] = _cvav;
                  std::copy(parsIt, parsIt + nParamsEOS, paramsThese.begin() + nParamsEOSHave);

                  _eosModel.setParams(paramsThese); parsIt += nParamsEOS;

                  _eosModel.getEV0(_e0, _v0);
               }

               int iParam = parsIt - pars.begin();
               if (iParam != nParams) {
                  ECMECH_FAIL(__func__, "wrong number of params");
               }

#ifdef __cuda_host_only__
               //////////////////////////////

               _rhvNames.clear();
               _rhvVals.clear();
               _rhvPlot.clear();
               _rhvState.clear();

               _rhvNames.push_back("shrate_eff"); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(true); // iHistA_shrateEff
               _rhvNames.push_back("shr_eff"); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(true); // iHistA_shrEff
               _rhvNames.push_back("n_feval"); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(true); // iHistA_nFEval
               // numHistAux
               //
               for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                  std::ostringstream os;
                  os << "t" << iTvec + 1;
                  _rhvNames.push_back(os.str()); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(true);
               }

               //
               {
                  double qVal = 1.0;
                  for (int iQ = 0; iQ < ecmech::qdim; ++iQ) {
                     std::ostringstream os;
                     os << "quat_" << iQ + 1;
                     _rhvNames.push_back(os.str()); _rhvVals.push_back(qVal); _rhvPlot.push_back(true); _rhvState.push_back(true);
                     qVal = 0.0;
                  }
               }
               //
               _kinetics.getHistInfo(_rhvNames, _rhvVals, _rhvPlot, _rhvState);
               //
               for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                  std::ostringstream os;
                  os << "shrate_" << iSlip + 1;
                  _rhvNames.push_back(os.str()); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(true);
               }

               //
               if (_rhvNames.size() != numHist) {
                  ECMECH_FAIL(__func__, "mismatch in numHist");
               }
#endif
            };

            __ecmech_host__
            void getParams(std::vector<int>& opts,
                           std::vector<double>& pars,
                           std::vector<std::string>& strs) const override {
               // ...*** ;
               opts.clear();
               strs.clear();
               pars.clear();
               ECMECH_FAIL(__func__, "getParams not yet implemented");
            };

            __ecmech_host__
            void logParameters(std::ostringstream& oss) const override {
               // ...*** ;
               oss << "evptn constitutive model" << std::endl;
               ECMECH_FAIL(__func__, "logParameters not yet implemented");
            };

            using matModelBase::getResponse;
            __ecmech_host__
            void getResponse(const double & dt,
                             const double * defRateV,
                             const double * spinV,
                             const double * volRatioV,
                             double * eIntV,
                             double * stressSvecPV,
                             double * histV,
                             double * tkelvV,
                             double * sddV,
                             double * mtanSDV,
                             const int& nPassed) const override final {
               RAJA::RangeSegment default_range(0, nPassed);
               //All of the stride lengths are constant within this function
               const uint def_rate_stride = _strides[istride_def_rate];
               const uint spin_v_stride = _strides[istride_spin_v];
               const uint vol_ratio_stride = _strides[istride_vol_ratio];
               const uint int_eng_stride = _strides[istride_int_eng];
               const uint stress_stride = _strides[istride_stress];
               const uint history_stride = _strides[istride_history];
               const uint tkelv_stride = _strides[istride_tkelv];
               const uint sdd_stride = _strides[istride_sdd];

#if defined(RAJA_ENABLE_OPENMP)
               if (_accel == ecmech::Accelerator::OPENMP) {
                  RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i) {
                  double *mtanSDThis = nullptr;
                  if (mtanSDV != nullptr) {
                     mtanSDThis = &mtanSDV[ecmech::nsvec2 * i];
                  }
                  getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>(_slipGeom, _kinetics, _elastN, _eosModel,
                                                                              dt,
                                                                              _tolerance,
                                                                              &defRateV[def_rate_stride * i],
                                                                              &spinV[spin_v_stride * i],
                                                                              &volRatioV[vol_ratio_stride * i],
                                                                              &eIntV[int_eng_stride * i],
                                                                              &stressSvecPV[stress_stride * i],
                                                                              &histV[history_stride * i],
                                                                              tkelvV[tkelv_stride * i],
                                                                              &sddV[sdd_stride * i],
                                                                              mtanSDThis,
                                                                              _outputLevel);
               });
               }
#endif
#if defined(RAJA_ENABLE_CUDA)
               if (_accel == ecmech::Accelerator::CUDA) {
                  RAJA::forall<RAJA::cuda_exec<RAJA_CUDA_THREADS> >(default_range, [ = ] RAJA_DEVICE(int i) {
                  double *mtanSDThis = nullptr;
                  if (mtanSDV != nullptr) {
                     mtanSDThis = &mtanSDV[ecmech::nsvec2 * i];
                  }
                  getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>(_slipGeom, _kinetics, _elastN, _eosModel,
                                                                              dt,
                                                                              _tolerance,
                                                                              &defRateV[def_rate_stride * i],
                                                                              &spinV[spin_v_stride * i],
                                                                              &volRatioV[vol_ratio_stride * i],
                                                                              &eIntV[int_eng_stride * i],
                                                                              &stressSvecPV[stress_stride * i],
                                                                              &histV[history_stride * i],
                                                                              tkelvV[tkelv_stride * i],
                                                                              &sddV[sdd_stride * i],
                                                                              mtanSDThis,
                                                                              _outputLevel);
               });
               }
#endif
               if (_accel == ecmech::Accelerator::CPU) {
                  RAJA::forall<RAJA::loop_exec>(default_range, [ = ] (int i) {
                  double *mtanSDThis = nullptr;
                  if (mtanSDV != nullptr) {
                     mtanSDThis = &mtanSDV[ecmech::nsvec2 * i];
                  }
                  getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>(_slipGeom, _kinetics, _elastN, _eosModel,
                                                                              dt,
                                                                              _tolerance,
                                                                              &defRateV[def_rate_stride * i],
                                                                              &spinV[spin_v_stride * i],
                                                                              &volRatioV[vol_ratio_stride * i],
                                                                              &eIntV[int_eng_stride * i],
                                                                              &stressSvecPV[stress_stride * i],
                                                                              &histV[history_stride * i],
                                                                              tkelvV[tkelv_stride * i],
                                                                              &sddV[sdd_stride * i],
                                                                              mtanSDThis,
                                                                              _outputLevel);
               });
               }
            }; // End of getResponse

            using matModelBase::getHistInfo;
            __ecmech_host__
            void getHistInfo(std::vector<std::string> & names,
                             std::vector<double>       & vals,
                             std::vector<bool>        & plot,
                             std::vector<bool>        & state) const override{
               if (_rhvNames.size() != numHist) {
                  ECMECH_FAIL(__func__, "have not yet set up history information");
               }
               names.resize(numHist); std::copy(_rhvNames.begin(), _rhvNames.end(), names.begin() );
               vals.resize(numHist); std::copy(_rhvVals.begin(), _rhvVals.end(), vals.begin() );
               plot.resize(numHist); std::copy(_rhvPlot.begin(), _rhvPlot.end(), plot.begin() );
               state.resize(numHist); std::copy(_rhvState.begin(), _rhvState.end(), state.begin() );
            };

            __ecmech_hdev__
            int getNumHist( ) const override{
               return numHist;
            };

            __ecmech_host__
            void setAccelerator(ecmech::Accelerator accel) override final {
               _accel = accel;
            }

         private:

            SlipGeom _slipGeom;
            Kinetics _kinetics;
            ThermoElastN _elastN;
            EosModel _eosModel;

            double _tolerance;
            int _outputLevel;
            uint _strides[8];
            ecmech::Accelerator _accel;

            std::vector<std::string> _rhvNames;
            std::vector<double>       _rhvVals;
            std::vector<bool>        _rhvPlot;
            std::vector<bool>        _rhvState;
      }; // class matModel
   } // namespace evptn
} // namespace ecmech

#endif // ecmech_evptnWrap_include
