#ifndef ECMECH_CLASSES_HPP
#define ECMECH_CLASSES_HPP

#include "ECMech_const.h"
#include "ECMech_core.h"
#include "ECMech_cases.h"

class pyevptn_base
{
    public:
    pyevptn_base(std::vector<double> &/*params*/) {}
    virtual ~pyevptn_base() {}
    /// Provides history names, initial values, 
    /// whether the param is meant to be plotted,
    /// and whether its a state variable
    void getHistoryInfo(std::vector<std::string>& names,
                        std::vector<double>& vals,
                        std::vector<bool>& plot,
                        std::vector<bool>& state)
    {
        const int numHist = m_rhvNames.size();

        names.resize(numHist); std::copy(m_rhvNames.begin(), m_rhvNames.end(), names.begin() );
        vals.resize(numHist); std::copy(m_rhvVals.begin(), m_rhvVals.end(), vals.begin() );
        plot.resize(numHist); std::copy(m_rhvPlot.begin(), m_rhvPlot.end(), plot.begin() );
        state.resize(numHist); std::copy(m_rhvState.begin(), m_rhvState.end(), state.begin() );
    };
    /// Should be called each time step and point to set-up problem
    virtual void setup( double    dt,
                        double    tolerance,
                        const double  * d_svec_kk_sm, // defRate,
                        const double  * w_veccp_sm, // spin
                        const double  * volRatio,
                        double  * eInt,
                        double  * stressSvecP,
                        double  * hist,
                        double  & tkelv) = 0;
    /// Should be handed to nonlinear solver for each time step and point
    /// for the Jacobian/gradient and function calculations
    virtual void computeRJ(double * const resid,
                           double * const Jacobian,
                           const double * const x) = 0;
    /// Should be called after a nonlinear solve for each time step and point
    /// to get out the state at end of each time step for a material point
    virtual void getState(const double * const x,
                          double * const eInt,
                          double * const stressSvecP,
                          double * const hist,
                          double& tkelv,
                          double * const sdd) = 0;
    protected:

    // Variables related to history names and their initial values
    std::vector<std::string> m_rhvNames;
    std::vector<double> m_rhvVals;
    std::vector<bool> m_rhvPlot;
    std::vector<bool> m_rhvState;

    // Various variables needed in set-up phases
    double m_dt;
    double m_tolerance;
    double m_d_svec_kk_sm[ecmech::nsvp]; // defRate,
    double m_w_veccp_sm[ecmech::nwvec]; // spin
    double m_volRatio[ecmech::nvr];
    double m_eInt[ecmech::ne];
    double m_stressSvecP[ecmech::nsvp];
    double m_quat_n[ecmech::qdim];
    double m_tkelv;
    double m_vNew;
    double m_pEOS;
    double m_d_vecd_sm[ecmech::ntvec];
    double m_e_vecd_n[ecmech::ntvec];
    double m_eNew;
    double m_eDevTot;
    double m_halfVMidDt;
    double m_rho0;
    double m_bulkNew;
    double m_cvav;
    double m_e0;
    double m_v0;

};

// As we add more variations of evptn classes here this won't be only one available
template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
class pyEvptn_norm : public pyevptn_base
{
    public:
    pyEvptn_norm(std::vector<double> &params) : pyevptn_base(params), 
                                                kinetics(SlipGeom::nslip)
    {

        if (params.size() != (unsigned int) nParams) {
            ECMECH_FAIL(__func__, "wrong number of params");
        }

        std::vector<double>::const_iterator parsIt = params.begin();

        m_rho0 = *parsIt; ++parsIt;
        m_cvav = *parsIt; ++parsIt;

        m_tolerance = *parsIt; ++parsIt;

        {
            const std::vector<double> paramsThese(parsIt, parsIt + SlipGeom::nParams);
            slipGeom.setParams(paramsThese); parsIt += SlipGeom::nParams;
        }
        {
            const std::vector<double> paramsThese(parsIt, parsIt + ThermoElastN::nParams);
            elastN.setParams(paramsThese); parsIt += ThermoElastN::nParams;
        }
        {
            const std::vector<double> paramsThese(parsIt, parsIt + Kinetics::nParams);
            kinetics.setParams(paramsThese); parsIt += Kinetics::nParams;
        }
        {
            double bulkMod = elastN.getBulkMod();
            std::vector<double> paramsThese(EosModel::nParams);
            paramsThese[0] = m_rho0;
            paramsThese[1] = bulkMod;
            paramsThese[2] = m_cvav;
            std::copy(parsIt, parsIt + nParamsEOS, paramsThese.begin() + nParamsEOSHave);

            eos.setParams(paramsThese); parsIt += nParamsEOS;

            {
                double vMin, vMax;
                eos.getInfo(vMin, vMax, m_e0, m_v0);
            }
        }

        int iParam = parsIt - params.begin();
        if (iParam != nParams) {
            ECMECH_FAIL(__func__, "wrong number of params");
        }
        //////////////////////////////

        m_rhvNames.clear();
        m_rhvVals.clear();
        m_rhvPlot.clear();
        m_rhvState.clear();

#if defined(ECMECH_USE_DPEFF)
        m_rhvNames.push_back("effective_plastic_deformation_rate"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrateEff
        m_rhvNames.push_back("effective_plastic_strain"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrEff
#else
        m_rhvNames.push_back("effective_shearing_rate"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrateEff
        m_rhvNames.push_back("effective_shearing"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true); // iHistA_shrEff
#endif
        m_rhvNames.push_back("flow_strength"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(false); // iHistA_flowStr
        m_rhvNames.push_back("n_function_evals"); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(false); // iHistA_nFEval
        // numHistAux
        //
        for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            std::ostringstream os;
            os << "deviatoric_elastic_strain_" << iTvec + 1;
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
        kinetics.getHistInfo(m_rhvNames, m_rhvVals, m_rhvPlot, m_rhvState);
        //
        for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
            std::ostringstream os;
            os << "shearing_rate_" << iSlip + 1;
            m_rhvNames.push_back(os.str()); m_rhvVals.push_back(0.); m_rhvPlot.push_back(true); m_rhvState.push_back(true);
        }

        //
        if (m_rhvNames.size() != numHist) {
            ECMECH_FAIL(__func__, "mismatch in numHist");
        }
    }

    virtual ~pyEvptn_norm() {}

    void setup( double    dt,
                double    tolerance,
                const double  * d_svec_kk_sm, // defRate,
                const double  * w_veccp_sm, // spin
                const double  * volRatio,
                double  * eInt,
                double  * stressSvecP,
                double  * hist,
                double  & tkelv) override final
    {
        static const int iHistLbGdot = ecmech::evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;

        m_dt = dt;
        m_tolerance = tolerance;
        m_tkelv = tkelv;

        for (int i = 0; i < ecmech::nsvp; i++)
        {
            m_d_svec_kk_sm[i] = d_svec_kk_sm[i];
            m_stressSvecP[i] = stressSvecP[i];
        }

        for (int i = 0; i < ecmech::nwvec; i++)
        {
            m_w_veccp_sm[i] = w_veccp_sm[i];
        }

        for (int i = 0; i < ecmech::ne; i++)
        {
            m_eInt[i] = eInt[i];
        }

        for (int i = 0; i < ecmech::nvr; i++)
        {
            m_volRatio[i] = volRatio[i];
        }

        for (int i = 0; i < numHist; i++)
        {
            m_hist[i] = hist[i];
        }

        // convert deformation rate convention
        //
        ecmech::svecToVecd(m_d_vecd_sm, m_d_svec_kk_sm);

        // pointers to state
        //
        double* h_state = &(hist[ecmech::evptn::iHistLbH]);
        double* gdot = &(hist[iHistLbGdot]);

        for (int i_hist = 0; i_hist < ecmech::ntvec; i_hist++) {
            m_e_vecd_n[i_hist] = hist[ecmech::evptn::iHistLbE + i_hist];
        }

        for (int i_hist = 0; i_hist < ecmech::qdim; i_hist++) {
            m_quat_n[i_hist] = hist[ecmech::evptn::iHistLbQ + i_hist];
        }

        //
        // normalize quat just in case
        ecmech::vecsVNormalize<ecmech::qdim>(m_quat_n);

        // EOS
        //
        double eOld = eInt[ecmech::i_ne_total];
        double pOld = m_stressSvecP[6];
        //
        // get tkelv from beginning-of-step to avoid tangent stiffness contributions
        {
            double pBOS;
            double vOld = m_volRatio[0];
            eos.evalPT(pBOS, m_tkelv, vOld, eOld);
        }
        //
        double tkelvNew ;
        {
        double dpde, dpdv, dtde;
        ecmech::updateSimple<EosModel>(eos, m_pEOS, tkelvNew, m_eNew, m_bulkNew,
                                       dpde, dpdv, dtde,
                                       m_volRatio[1], m_volRatio[3],
                                       eOld, pOld);
        }

        // update hardness state to the end of the step
        // gdot is still at beginning-of-step
        //
        kinetics.updateH(m_hard_u, h_state, dt, gdot);
        m_vNew = m_volRatio[1];

        if (prob != nullptr)
        {
            delete prob;
        }

        prob = new ecmech::evptn::EvptnUpdstProblem<SlipGeom, Kinetics, ThermoElastN>
                    (slipGeom, kinetics, elastN,
                    m_dt,
                    m_vNew, m_eNew, m_pEOS, m_tkelv,
                    m_hard_u, m_e_vecd_n, m_quat_n,
                    m_d_vecd_sm, m_w_veccp_sm);

    }

    void computeRJ(double * const resid,
                   double * const Jacobian,
                   const double * const x) override final
    {
        prob->computeRJ(resid, Jacobian, x);
    };

    void getState(const double * const x,
                          double * const eInt,
                          double * const stressSvecP,
                          double * const hist,
                          double& tkelv,
                          double * const sdd) override final 
    {
        double* h_state = &(hist[ecmech::evptn::iHistLbH]);
        double* gdot = &(hist[iHistLbGdot]);
        double* e_vecd_u = &(hist[ecmech::evptn::iHistLbE]);
        double* quat_u = &(hist[ecmech::evptn::iHistLbQ]);
        // store updated state
        //
        prob->stateFromX(e_vecd_u, quat_u, x);
        for (int i_hstate = 0; i_hstate < Kinetics::nH; i_hstate++) {
            h_state[i_hstate] = m_hard_u[i_hstate];
        }
        //
        {
            const double* gdot_u = prob->getGdot();
            for (int i_gdot = 0; i_gdot < SlipGeom::nslip; i_gdot++) {
                gdot[i_gdot] = gdot_u[i_gdot];
            }
        }
        //
        hist[ecmech::evptn::iHistA_shrateEff] = prob->getShrateEff();
        hist[ecmech::evptn::iHistA_shrEff] += hist[ecmech::evptn::iHistA_shrateEff] * m_dt;
        //
        {
            double dEff = ecmech::vecd_Deff(m_d_vecd_sm);
            double flow_strength = prob->getHdnScale();
            if (dEff > ecmech::idp_tiny_sqrt) {
                flow_strength = prob->getDisRate() / dEff;
            }
            hist[ecmech::evptn::iHistA_flowStr] = flow_strength;
        }
        //
        hist[ecmech::evptn::iHistA_nFEval] = 0; // does _not_ include updateH iterations

        // get Cauchy stress
        //
        double Cstr_vecds_lat[ecmech::nsvec];
        prob->elastNEtoC(Cstr_vecds_lat, e_vecd_u);

        double C_matx[ecmech::ndim * ecmech::ndim];
        ecmech::quat_to_tensor(C_matx, quat_u);
        //
        double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
        ecmech::get_rot_mat_vecd(qr5x5_ls, C_matx);
        //
        double Cstr_vecds_sm[ecmech::nsvec];
        ecmech::vecsVMa<ecmech::ntvec>(Cstr_vecds_sm, qr5x5_ls, Cstr_vecds_lat);
        Cstr_vecds_sm[ecmech::iSvecS] = Cstr_vecds_lat[ecmech::iSvecS];
        //
        // put end-of-step stress in stressSvecP
        ecmech::vecdsToSvecP(stressSvecP, Cstr_vecds_sm);
        //
        // and now the second half of the trapezoidal integration
        //
        m_eDevTot += m_halfVMidDt * ecmech::vecsInnerSvecDev(stressSvecP, m_d_svec_kk_sm);

        // adjust sign on quat so that as close as possible to quat_o;
        // more likely to keep orientations clustered this way;
        // this flip through the origin is equivalent under antipodal symmetry
        //
        if (ecmech::vecsyadotb<ecmech::qdim>(quat_u, m_quat_n) < ecmech::zero) {
        for (int iQ = 0; iQ < ecmech::qdim; ++iQ) {
            quat_u[iQ] = -quat_u[iQ];
        }
        }

        {
        double gmod = elastN.getGmod(m_tkelv, m_pEOS, m_eNew);
        sdd[ecmech::i_sdd_bulk] = m_bulkNew;
        sdd[ecmech::i_sdd_gmod] = gmod;
        tkelv = m_tkelv;

        }
#ifdef ECMECH_DEBUG
        assert(ecmech::nsdd == 2);
#endif

        m_eNew += m_eDevTot;
        //
        // could update pressure and temperature again, but do not bother

        eInt[ecmech::i_ne_total] = m_eNew;
#ifdef ECMECH_DEBUG
        assert(ecmech::ne == 1);
#endif 
    }
    private:

    SlipGeom slipGeom;
    Kinetics kinetics;
    ThermoElastN elastN;
    EosModel eos;
    ecmech::evptn::EvptnUpdstProblem<SlipGeom, Kinetics, ThermoElastN>* prob = nullptr;


    static constexpr int iHistLbGdot = ecmech::evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;
    static constexpr int numHist = ecmech::evptn::NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::numHist;
    static constexpr int nH = Kinetics::nH;
    static constexpr int nslip = SlipGeom::nslip;

    static constexpr int nParamsEOSHave = 3; // number that get from 'elsewhere' // these are assumed to go in first
    static constexpr int nParamsEOS = EosModel::nParams - nParamsEOSHave;
    static constexpr int nParams = 2 + 1 + // rho0, cvav, tolerance
                         Kinetics::nParams + ThermoElastN::nParams + nParamsEOS;

    //
    double m_hist[numHist];
    double m_hard_u[nH];

};


// Should probably figure out a better way to do this sometime in the future...
typedef pyEvptn_norm<ecmech::SlipGeomFCC, ecmech::Kin_FCC_A, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_FCC_A;
typedef pyEvptn_norm<ecmech::SlipGeomFCC, ecmech::Kin_FCC_AH, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_FCC_AH;
typedef pyEvptn_norm<ecmech::SlipGeomFCC, ecmech::Kin_FCC_B, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_FCC_B;
typedef pyEvptn_norm<ecmech::SlipGeom_BCC_A, ecmech::Kin_BCC_A, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_BCC_A;
typedef pyEvptn_norm<ecmech::SlipGeom_BCC_A, ecmech::Kin_FCC_A, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_BCC_B;
typedef pyEvptn_norm<ecmech::SlipGeom_BCC_A, ecmech::Kin_FCC_AH, ecmech::evptn::ThermoElastNCubic, ecmech::EosModelConst<false> > pyMatModelEvptn_BCC_BH;
typedef pyEvptn_norm<ecmech::SlipGeom_HCP_A, ecmech::Kin_HCP_A, ecmech::evptn::ThermoElastNHexag, ecmech::EosModelConst<false> > pyMatModelEvptn_HCP_A;

#endif