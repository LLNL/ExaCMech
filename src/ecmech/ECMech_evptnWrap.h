// -*-c++-*-

#ifndef ecmech_evptnWrap_include
#define ecmech_evptnWrap_include

#include "ECMech_core.h"

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

            static const int iHistLbGdot = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;
            static const int numHist = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
            static const int nH = Kinetics::nH;
            static const int nslip = SlipGeom::nslip;

            // this are assumed to go in first
            static const int nParamsEOSHave = 3; // number that get from 'elsewhere'

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
            {};

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
            void initFromParams(const std::vector<int>         & opts,
                                const std::vector<double>       & pars,
                                const std::vector<std::string> & strs) override {
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
            void getParams(std::vector<int>         & opts,
                           std::vector<double>       & pars,
                           std::vector<std::string> & strs) const {
               // ...*** ;
               opts.clear();
               strs.clear();
               pars.clear();
               ECMECH_FAIL(__func__, "getParams not yet implemented");
            };

            __ecmech_host__
            void logParameters(std::ostringstream & oss) const {
               // ...*** ;
               oss << "evptn constitutive model" << std::endl;
               ECMECH_FAIL(__func__, "logParameters not yet implemented");
            };

            using matModelBase::getResponse;
            __ecmech_hdev__
            void getResponse(const double  & dt,
                             const double  * defRateV,
                             const double  * spinV,
                             const double  * volRatioV,
                             double  * eIntV,
                             double  * stressSvecPV,
                             double  * histV,
                             double  * tkelvV,
                             double  * sddV,
                             double  * mtanSDV,
                             const int    & nPassed) const override final {
               double *mtanSDThis = nullptr;
               //
               for (int i = 0; i<nPassed; ++i) {
                  if (mtanSDV != nullptr) {
                     mtanSDThis = &mtanSDV[ecmech::nsvec2 * i];
                  }
                  getResponseSngl<SlipGeom, Kinetics, ThermoElastN, EosModel>(_slipGeom, _kinetics, _elastN, _eosModel,
                                                                              dt,
                                                                              _tolerance,
                                                                              &defRateV[ecmech::nsvp*i],
                                                                              &spinV[ecmech::ndim*i],
                                                                              &volRatioV[ecmech::nvr*i],
                                                                              &eIntV[ecmech::ne*i],
                                                                              &stressSvecPV[ecmech::nsvp*i],
                                                                              &histV[numHist * i],
                                                                              tkelvV[i],
                                                                              &sddV[ecmech::nsdd*i],
                                                                              mtanSDThis,
                                                                              _outputLevel);
               }
            };

            #ifdef __cuda_host_only__
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
               std::vector<bool> _rhvPlotTemp(_rhvPlot, _rhvPlot + _rhvPlotLen);;
               std::vector<bool> _rhvStateTemp(_rhvState, _rhvState + _rhvStateLen);;

               names.resize(numHist); std::copy(_rhvNamesTemp.begin(), _rhvNamesTemp.end(), names.begin() );
               vals.resize(numHist); std::copy(_rhvValsTemp.begin(), _rhvValsTemp.end(), vals.begin() );
               plot.resize(numHist); std::copy(_rhvPlotTemp.begin(), _rhvPlotTemp.end(), plot.begin() );
               state.resize(numHist); std::copy(_rhvStateTemp.begin(), _rhvStateTemp.end(), state.begin() );
            };
            #endif

            int getNumHist( ) const {
               return numHist;
            };

         private:

            SlipGeom _slipGeom;
            Kinetics _kinetics;
            ThermoElastN _elastN;
            EosModel _eosModel;

            double _tolerance;
            int _outputLevel;

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
