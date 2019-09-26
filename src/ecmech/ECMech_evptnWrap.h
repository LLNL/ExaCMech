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
            __ecmech_hdev__
            matModel()
               : matModelBase(),
               _kinetics(SlipGeom::nslip),
               _outputLevel(0),
               _rhvNames(nullptr),
               _rhvNamesLen(0),
               _rhvNamesPtrEndLoc(nullptr),
               _rhvVals(nullptr),
               _rhvValsLen(0),
               _rhvPlot(nullptr),
               _rhvPlotLen(0),
               _rhvState(nullptr),
               _rhvStateLen(0)
            {
               //Should the tangent stiff matrix be included in these stride calculations?
               _strides[0] = ecmech::nsvp;
               _strides[1] = ecmech::ndim;
               _strides[2] = ecmech::nvr;
               _strides[3] = ecmech::ne;
               _strides[4] = ecmech::nsvp;
               _strides[5] = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
               _strides[6] = 1;
               _strides[7] = ecmech::nsdd;
            };

            // constructor
            __ecmech_hdev__
            matModel(const uint* strides, const uint stride_len)
               : matModelBase(),
               _kinetics(SlipGeom::nslip),
               _outputLevel(0),
               _rhvNames(nullptr),
               _rhvNamesLen(0),
               _rhvNamesPtrEndLoc(nullptr),
               _rhvVals(nullptr),
               _rhvValsLen(0),
               _rhvPlot(nullptr),
               _rhvPlotLen(0),
               _rhvState(nullptr),
               _rhvStateLen(0)
            {
               uint nhist = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;

               if (stride_len != 8) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "Stride vector needs to have a size of 8 with strides of: " <<
                     ecmech::nsvp << ", " << ecmech::ndim << ", " << ecmech::nvr << ", " <<
                     ecmech::ne << ", " << ecmech::nsvp << ", " << nhist << ", 1, " << ecmech::nsdd;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "Stride vector needs to have a size of 8 with strides");
                  #endif
               }

               //Need to make sure all of the strides provided at least make sense
               if (strides[0] < ecmech::nsvp) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[0] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }
               if (strides[1] < ecmech::ndim) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[1] should have at least a length of: " << ecmech::ndim;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
                  
               }
               if (strides[2] < ecmech::nvr) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[2] should have at least a length of: " << ecmech::nvr;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }
               if (strides[3] < ecmech::ne) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[3] should have at least a length of: " << ecmech::ne;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }
               if (strides[4] < ecmech::nsvp) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[4] should have at least a length of: " << ecmech::nsvp;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }
               if (strides[5] < nhist) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[5] should have at least a length of: " << nhist;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }
               if (strides[6] < 1) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[6] should have at least a length of: " << 1;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }
               if (strides[7] < ecmech::nsdd) {
                  #ifdef __cuda_host_only__
                  std::ostringstream os;
                  os << "strides[7] should have at least a length of: " << ecmech::nsdd;
                  ECMECH_FAIL(__func__, os.str());
                  #else
                  ECMECH_FAIL(__func__, "One of the stride lengths was not long enough");
                  #endif
               }

               _strides[0] = strides[0];
               _strides[1] = strides[1];
               _strides[2] = strides[2];
               _strides[3] = strides[3];
               _strides[4] = strides[4];
               _strides[5] = strides[5];
               _strides[6] = strides[6];
               _strides[7] = strides[7];
            };

            // deconstructor
            __ecmech_hdev__
            ~matModel()
            {
               delete _rhvNames;
               delete _rhvNamesPtrEndLoc;
               delete _rhvVals;
               delete _rhvPlot;
               delete _rhvState;
            }

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

               //In the new method we would need to delete all data for the pointers
               //We would also set all of the lengths back to 0
               //Here we would declare temporary std::vectors here we can just add temp to the names
               //fix me...

               std::vector<std::string> _rhvNamesTemp;
               std::vector<double> _rhvValsTemp;
               std::vector<bool> _rhvPlotTemp;
               std::vector<bool> _rhvStateTemp;

               delete _rhvNames; _rhvNamesLen = 0; delete _rhvNamesPtrEndLoc;
               delete _rhvVals; _rhvValsLen = 0;
               delete _rhvPlot; _rhvPlot = 0;
               delete _rhvState; _rhvState = 0;

               _rhvNamesTemp.push_back("shrate_eff"); _rhvValsTemp.push_back(0.); _rhvPlotTemp.push_back(true);
               _rhvStateTemp.push_back(true);                                                                                                  // iHistA_shrateEff
               _rhvNamesTemp.push_back("shr_eff"); _rhvValsTemp.push_back(0.); _rhvPlotTemp.push_back(true);
               _rhvStateTemp.push_back(true);                                                                                               // iHistA_shrEff
               _rhvNamesTemp.push_back("n_feval"); _rhvValsTemp.push_back(0.); _rhvPlotTemp.push_back(true);
               _rhvStateTemp.push_back(true);                                                                                               // iHistA_nFEval
               // numHistAux
               //
               for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                  std::ostringstream os;
                  os << "t" << iTvec + 1;
                  _rhvNamesTemp.push_back(os.str()); _rhvValsTemp.push_back(0.); _rhvPlotTemp.push_back(true);
                  _rhvStateTemp.push_back(true);
               }

               //
               {
                  double qVal = 1.0;
                  for (int iQ = 0; iQ < ecmech::qdim; ++iQ) {
                     std::ostringstream os;
                     os << "quat_" << iQ + 1;
                     _rhvNamesTemp.push_back(os.str()); _rhvValsTemp.push_back(qVal); _rhvPlotTemp.push_back(true);
                     _rhvStateTemp.push_back(true);
                     qVal = 0.0;
                  }
               }
               //
               _kinetics.getHistInfo(_rhvNamesTemp, _rhvValsTemp, _rhvPlotTemp, _rhvStateTemp);
               //
               for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                  std::ostringstream os;
                  os << "shrate_" << iSlip + 1;
                  _rhvNamesTemp.push_back(os.str()); _rhvValsTemp.push_back(0.); _rhvPlotTemp.push_back(true);
                  _rhvStateTemp.push_back(true);
               }

               //
               if (_rhvNamesTemp.size() != numHist) {
                  ECMECH_FAIL(__func__, "mismatch in numHist");
               }

               //Now we can take all of the temp data stored in the vector and initialize all
               //the data using the new construct. The double and bool data will be rather trivial.
               //However, the string data will take a bit more effort then previously.

               _rhvVals = new double [_rhvValsTemp.size()];
               _rhvValsLen = _rhvValsTemp.size();
               _rhvState = new bool [_rhvStateTemp.size()];
               _rhvStateLen = _rhvStateTemp.size();
               _rhvPlot = new bool [_rhvPlotTemp.size()];
               _rhvPlotLen = _rhvPlotTemp.size();
               _rhvNamesPtrEndLoc = new uint [_rhvNamesTemp.size()];

               std::copy(_rhvValsTemp.begin(), _rhvValsTemp.end(), _rhvVals);
               std::copy(_rhvStateTemp.begin(), _rhvStateTemp.end(), _rhvState);
               std::copy(_rhvPlotTemp.begin(), _rhvPlotTemp.end(), _rhvPlot);
               //The below several lines are used to construct our raw char array
               //from the std::vector<std::string>. We store all the strings
               //contiguously in the raw char array with the null character included.
               uint offset = 0;
               //the _rhvNamesPtrEndLoc stores the starting loc for the next string
               for (uint i = 0; i < _rhvNamesTemp.size(); i++) {
                  offset += _rhvNamesTemp[i].size() + 1;
                  _rhvNamesPtrEndLoc[i] = offset;
               }

               _rhvNamesLen = offset;

               _rhvNames = new char [offset];
               uint init = 0;
               for (uint i = 0; i < _rhvNamesTemp.size(); i++) {
                  char* ca = _rhvNames + init;
                  std::copy(_rhvNamesTemp[i].begin(), _rhvNamesTemp[i].end(), ca);
                  ca[_rhvNamesTemp[i].size()] = '\0';
                  init = _rhvNamesPtrEndLoc[i];
               }

#endif
            };

            __ecmech_host__
            void getParams(std::vector<int>& opts,
                           std::vector<double>& pars,
                           std::vector<std::string>& strs) const {
               // ...*** ;
               opts.clear();
               strs.clear();
               pars.clear();
               ECMECH_FAIL(__func__, "getParams not yet implemented");
            };

            __ecmech_host__
            void logParameters(std::ostringstream& oss) const {
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
                                                                              &defRateV[_strides[0] * i],
                                                                              &spinV[_strides[1] * i],
                                                                              &volRatioV[_strides[2] * i],
                                                                              &eIntV[_strides[3] * i],
                                                                              &stressSvecPV[_strides[4] * i],
                                                                              &histV[_strides[5] * i],
                                                                              tkelvV[_strides[6] * i],
                                                                              &sddV[_strides[7] * i],
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
                                                                              &defRateV[_strides[0] * i],
                                                                              &spinV[_strides[1] * i],
                                                                              &volRatioV[_strides[2] * i],
                                                                              &eIntV[_strides[3] * i],
                                                                              &stressSvecPV[_strides[4] * i],
                                                                              &histV[_strides[5] * i],
                                                                              tkelvV[_strides[6] * i],
                                                                              &sddV[_strides[7] * i],
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
                                                                              &defRateV[_strides[0] * i],
                                                                              &spinV[_strides[1] * i],
                                                                              &volRatioV[_strides[2] * i],
                                                                              &eIntV[_strides[3] * i],
                                                                              &stressSvecPV[_strides[4] * i],
                                                                              &histV[_strides[5] * i],
                                                                              tkelvV[_strides[6] * i],
                                                                              &sddV[_strides[7] * i],
                                                                              mtanSDThis,
                                                                              _outputLevel);
               });
               }
            };//End of getResponse

            using matModelBase::getHistInfo;
            __ecmech_host__
            virtual void getHistInfo(std::vector<std::string> & names,
                                     std::vector<double>       & vals,
                                     std::vector<bool>        & plot,
                                     std::vector<bool>        & state) const override {
               if (_rhvValsLen != numHist) {
                  ECMECH_FAIL(__func__, "have not yet set up history information");
               }

               //The below replaces all of the previous std::vector usages
               //of _rhv* with the raw pointers. We still make use of
               //std::vectors to pass things back to the users.
               std::vector<std::string> _rhvNamesTemp;

               uint offset = 0;

               for (uint i = 0; i < numHist; i++) {
                  std::string strs(_rhvNames + offset);
                  _rhvNamesTemp.push_back(strs);
                  offset = _rhvNamesPtrEndLoc[i];
               }

               std::vector<double> _rhvValsTemp(_rhvVals, _rhvVals + _rhvValsLen);
               std::vector<bool> _rhvPlotTemp(_rhvPlot, _rhvPlot + _rhvPlotLen);
               std::vector<bool> _rhvStateTemp(_rhvState, _rhvState + _rhvStateLen);

               names.resize(numHist); std::copy(_rhvNamesTemp.begin(), _rhvNamesTemp.end(), names.begin() );
               vals.resize(numHist); std::copy(_rhvValsTemp.begin(), _rhvValsTemp.end(), vals.begin() );
               plot.resize(numHist); std::copy(_rhvPlotTemp.begin(), _rhvPlotTemp.end(), plot.begin() );
               state.resize(numHist); std::copy(_rhvStateTemp.begin(), _rhvStateTemp.end(), state.begin() );
            };

            __ecmech_hdev__
            int getNumHist( ) const {
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

            //I'm not sure we can even have these in here when dealing
            //with host and device code. So, I've replaced them with
            //raw ptr arrays instead along with the lengths of each array.
            // #ifdef __cuda_host_only__
            // // TO_DO : it is a good idea to have these host-only?
            // std::vector<std::string> _rhvNames;
            // std::vector<double>       _rhvVals;
            // std::vector<bool>        _rhvPlot;
            // std::vector<bool>        _rhvState;
            // #else//Not used at all but it keeps the compiler happy
            //If we want to have the equivalent of the std::vector for code that's
            //works on both the Device and Host we might need to make due with something
            //like the below where we hold raw arrays to everything. Since, we keep track
            //of the number of elements and in the string case the various pointer locations
            //we should be able to reconstruct the behavior of the vectors all the same...
            char* _rhvNames; uint _rhvNamesLen; uint* _rhvNamesPtrEndLoc;
            double* _rhvVals; uint _rhvValsLen;
            bool* _rhvPlot; uint _rhvPlotLen;
            bool* _rhvState; uint _rhvStateLen;
            // #endif
      }; // class matModel
   } //  namespace evptn
} // namespace ecmech

#endif // ecmech_evptnWrap_include
