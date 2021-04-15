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

            static constexpr int nParamsEOSHave = 3; // number that get from 'elsewhere' // these are assumed to go in first
            static constexpr int nParamsEOS = EosModel::nParams - nParamsEOSHave;
            static constexpr int nParams = 
               2 + 1 + // rho0, cvav, tolerance
               Kinetics::nParams + ThermoElastN::nParams + nParamsEOS;
         
            // constructor
            __ecmech_host__
            matModel()
               : matModelBase(),
               _kinetics(SlipGeom::nslip)
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
            matModel(const unsigned int* strides, const unsigned int stride_len)
               : matModelBase(),
               _kinetics(SlipGeom::nslip)
            {
               unsigned int nhist = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;

               if (stride_len != ecmech::nstride) {
#ifdef __cuda_host_only__
                  // the order here needs to be consistent with ISTRIDE_* macros in ECMECH_const.h
                  std::ostringstream os;
                  os << "Stride vector needs to have a size of " << ecmech::nstride << " with strides of at least: " <<
                     ecmech::nsvp << ", " << ecmech::ndim << ", " << ecmech::nvr << ", " <<
                     ecmech::ne << ", " << ecmech::nsvp << ", " << nhist << ", 1, " << ecmech::nsdd
                     ;
                  ECMECH_FAIL(__func__, os.str());
#else
                  ECMECH_FAIL(__func__, "Stride vector is the wrong size");
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
               for (unsigned int i = 0; i < stride_len; i++) {
                  _strides[i] = strides[i];
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
               _opts = opts ;
               _pars = pars ;
               _strs = strs ;
               
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

               _rho0 = *parsIt; ++parsIt;
               _cvav = *parsIt; ++parsIt;

               _tolerance = *parsIt; ++parsIt;

               {
                  const std::vector<double> paramsThese(parsIt, parsIt + SlipGeom::nParams);
                  _slipGeom.setParams(paramsThese); parsIt += SlipGeom::nParams;
               }
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

                  {
                     double vMin, vMax;
                     _eosModel.getInfo(vMin, vMax, _e0, _v0);
                  }
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
               _rhvNames.push_back("flow_str"); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(false); // iHistA_flowStr
               _rhvNames.push_back("n_feval"); _rhvVals.push_back(0.); _rhvPlot.push_back(true); _rhvState.push_back(false); // iHistA_nFEval
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

            using matModelBase::getParams;
            __ecmech_host__
            void getParams(std::vector<int>& opts,
                           std::vector<double>& pars,
                           std::vector<std::string>& strs) const override final
            {
               opts = _opts;
               pars = _pars;
               strs = _strs;
            };

            using matModelBase::getResponseECM;
            __ecmech_host__
            void getResponseECM(const double & dt,
                                const double * defRateV,
                                const double * spinV,
                                const double * volRatioV,
                                double * eIntV,
                                double * stressSvecPV,
                                double * histV,
                                double * tkelvV,
                                double * sddV,
                                double * mtanSDV,
                                const int& nPassed) const override final
            {
               if ( !_complete ) {
                  ECMECH_FAIL(__func__,"not complete");
               }

               RAJA::RangeSegment default_range(0, nPassed);
               // All of the stride lengths are constant within this function
               const unsigned int def_rate_stride = _strides[istride_def_rate];
               const unsigned int spin_v_stride = _strides[istride_spin_v];
               const unsigned int vol_ratio_stride = _strides[istride_vol_ratio];
               const unsigned int int_eng_stride = _strides[istride_int_eng];
               const unsigned int stress_stride = _strides[istride_stress];
               const unsigned int history_stride = _strides[istride_history];
               const unsigned int tkelv_stride = _strides[istride_tkelv];
               const unsigned int sdd_stride = _strides[istride_sdd];

               switch ( _accel ) {
#if defined(RAJA_ENABLE_OPENMP)
                  case ECM_EXEC_STRAT_OPENMP :
                     RAJA::forall<RAJA::omp_parallel_for_exec>(default_range, [ = ] (int i) {
                           double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                           getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                              (_slipGeom, _kinetics, _elastN, _eosModel,
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
                     break;
#endif
#if defined(RAJA_ENABLE_CUDA)
                  case ECM_EXEC_STRAT_CUDA :
                     RAJA::forall<RAJA::cuda_exec<RAJA_CUDA_THREADS> >(default_range, [ = ] RAJA_DEVICE(int i) {
                           double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                           getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                              (_slipGeom, _kinetics, _elastN, _eosModel,
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
                     break;
#endif
                  case ECM_EXEC_STRAT_CPU :
                  default : // fall through to CPU if other options are not available
                     RAJA::forall<RAJA::loop_exec>(default_range, [ = ] (int i) {
                           double *mtanSDThis       = ( mtanSDV ? &mtanSDV[ecmech::nsvec2 * i] : nullptr );
                           getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>
                              (_slipGeom, _kinetics, _elastN, _eosModel,
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
                     break;
               } // switch _accel

            }; // End of getResponse

            using matModelBase::getHistInfo;
            __ecmech_host__
            void getHistInfo(std::vector<std::string> & names,
                             std::vector<double>      & vals,
                             std::vector<bool>        & plot,
                             std::vector<bool>        & state) const override final {
               if (_rhvNames.size() != numHist) {
                  ECMECH_FAIL(__func__, "have not yet set up history information");
               }
               names.resize(numHist); std::copy(_rhvNames.begin(), _rhvNames.end(), names.begin() );
               vals.resize(numHist); std::copy(_rhvVals.begin(), _rhvVals.end(), vals.begin() );
               plot.resize(numHist); std::copy(_rhvPlot.begin(), _rhvPlot.end(), plot.begin() );
               state.resize(numHist); std::copy(_rhvState.begin(), _rhvState.end(), state.begin() );
            };

            __ecmech_host__
            int getNumHist( ) const override final {
               return numHist;
            };

            __ecmech_host__
            void complete( ) override final
            {
               _bulkRef = _eosModel.getBulkRef();
               _complete = true;
            };

            // Constant getter functions to return the underlying templated classes.
            // Uses for these could be for example to compute the sample D^p tensor
            // using the symmetric schmid tensor from the SlipGeom class.
            const SlipGeom & getSlipGeom() const { return _slipGeom; }
            const Kinetics & getKinetics() const { return _kinetics; }
            const ThermoElastN & getThermoElastN() const { return _elastN; }
            const EosModel & getEosModel() const { return _eosModel; }

         private:

            SlipGeom _slipGeom;
            Kinetics _kinetics;
            ThermoElastN _elastN;
            EosModel _eosModel;

            double _tolerance;
            unsigned int _strides[ecmech::nstride];

            std::vector<std::string> _rhvNames;
            std::vector<double>      _rhvVals;
            std::vector<bool>        _rhvPlot;
            std::vector<bool>        _rhvState;

            // keep initFromParams vectors as a convenience
            std::vector<int>          _opts;
            std::vector<double>       _pars;
            std::vector<std::string>  _strs;
         
      }; // class matModel
   } // namespace evptn
} // namespace ecmech

#endif // ecmech_evptnWrap_include
