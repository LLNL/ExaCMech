#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_elastic.h"
#include "ECMech_eosSimple.h"
#include "ECMech_base_classes.h"
#include "ECMech_base_fcns.h"

#include "ECMech_evptn.h"

#include "SNLS_TrDLDenseG.h"
#include "SNLS_HybrdTrDLDenseG.h"

namespace ecmech {
namespace evptn {

__ecmech_hdev__
template<class SlipGeom, class Kinetics, class EosModel, class ThermoElastN, class ProbState, bool RStarSolve=false>
inline
void preprocess(const SlipGeom& slipGeom,
                const Kinetics& kinetics,
                const EosModel& eos,
                const ThermoElastN& thermoElastN,
                const double* const volRatio,
                const double* const eInt,
                const double* const d_svec_kk_sm,
                ProbState& prob_state,
                double& halfVMidDt,
                double& eDevTot)
{
    // total increment in the deviatoric part of the strain energy using
    // trapezoidal rule integration
    //
    // just beginning-of-step stress part so far
    //
    halfVMidDt = oneqrtr * (volRatio[0] + volRatio[1]) * prob_state.dt;
    eDevTot = halfVMidDt * vecsInnerSvecDev(prob_state.stressSvecP, d_svec_kk_sm);

    // EOS
    //
    const double eOld = eInt[ecmech::i_ne_total];
    //
    // get tkelv from beginning-of-step to avoid tangent stiffness contributions
    {
        double pBOS;
        const double vOld = volRatio[0];
        eos.evalPT(pBOS, prob_state.tkelv, vOld, eOld);
    }

    {
        const double pOld = prob_state.stressSvecP[6];
        double tkelvNew, dpde, dpdv, dtde;
        updateSimple(eos, prob_state.pEOS, tkelvNew, prob_state.eNew, prob_state.bulkNew,
                     dpde, dpdv, dtde,
                     volRatio[1], volRatio[3],
                     eOld, pOld);
    }

    // update hardness state to the end of the step
    // gdot is still at beginning-of-step
    //
    constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip) : 1;
    double hvals[nslip_dyn] = {}; // additional values needed to update the hardening state
    if constexpr(SlipGeom::dynamic) {

        double T_vecds[ecmech::nsvec] = {};
        double elas_dev_vol_vec[ecmech::nsvec] = {};
        double SvecP[ecmech::nsvec+1] = {};

       double m_a_vol = pow(prob_state.vNew, onethird);
       double m_inv_a_vol = 1.0 / prob_state.vNew;

        vecsVxa<ntvec>(elas_dev_vol_vec, m_inv_a_vol, prob_state.e_vecd_n);
        elas_dev_vol_vec[iSvecS] = sqr3 * log(m_a_vol);
        thermoElastN.eval(T_vecds, elas_dev_vol_vec, prob_state.tkelv, prob_state.pEOS, prob_state.eNew);
        vecdsToSvecP(SvecP, T_vecds);

        // For dynamic slip systems we need the chi angle
        double P[ecmech::ntvec * SlipGeom::nslip];
        double Q[ecmech::nwvec * SlipGeom::nslip];
        // still need to rotate stress state back to original value
        slipGeom.getPQ(hvals, P, Q, SvecP);
    }
    kinetics.updateH(prob_state.h_state_u, prob_state.h_state, prob_state.dt, prob_state.gdot, hvals, prob_state.tkelv);
#if defined(ECMECH_EXTRA_SOLVERS)
    if constexpr (RStarSolve) {
        auto prob = RotUpdProblem(slipGeom, thermoElastN, prob_state);
         // update Rstar aka Rdot * dt
         // gdot is still at beginning-of-step
        snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
        const bool status = main_problem(1e-8, solver, 0);
        if (!status) {
            return;
        }

        prob.stateFromX(prob_state.quat_u, solver._x);
    }

#endif
}

__ecmech_hdev__
template<class SNLS_Solver>
inline
bool main_problem(const double tolerance,
                  SNLS_Solver& solver,
                  const int outputLevel)
{
    snls::TrDeltaControl deltaControl;
    deltaControl._deltaInit = 1e0;
    {
        static constexpr int maxIter = 200;
        solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
    }

    // set initial guess
    //
    for (int iX = 0; iX < solver.getNDim(); ++iX) {
        solver._x[iX] = 0e0;
    }

    snls::SNLSStatus_t status = solver.solve( );
    if (status < snls::converged ) {
#if defined(__ecmech_host_only__)
        std::cout << "trust region solver residual " << solver.getRes() << " exit status " << status << std::endl;
        ECMECH_WARN(__func__, "Solver(s) failed to converge -- will try again with implicit elastic strain solve only");
#endif
        return false;
    }
    return true;
}

__ecmech_hdev__
template<class Problem, class Solver, class ProblemState>
inline
void computeTangentStiffness(Problem& prob,
                             Solver& solver,
                             ProblemState& prob_state,
                             double* const mtanSD)
{
    double mtanSD_vecds[ ecmech::nsvec2 ] = {};
    prob.provideMTan(mtanSD_vecds);
    {
        double residual[Problem::nDimSys] = {};
        double Jacobian[Problem::nDimSys * Problem::nDimSys] = {};
        solver.computeRJ(&residual[0], &Jacobian[0]);
    }
    prob.clearMTan();
    // currently have derivative with-respsect-to deformation rate;
    // to get derivative with-respsect-to strain increment,
    // multiply by 1/dt
    //
    double dt_ri = prob.getDtRi();
    for (int i = 0; i < ecmech::nsvec2; ++i) {
        mtanSD_vecds[i] = mtanSD_vecds[i] * dt_ri;
    }

    // contribution to stiffness from EOS
    // this is a bit crude, but should do the trick for now;
    // neglects effect of pEOS and vNew on workings of evptn
    //
    mtanSD_vecds[ECMECH_NN_INDX(iSvecS, iSvecS, ecmech::nsvec)] = three * prob_state.bulkNew;

    // convert from vecds notation to svec notation
    //
    mtan_conv_sd_svec<true>(mtanSD, mtanSD_vecds);
}

__ecmech_hdev__
template<int kinNH, class Problem, class ProblemState>
inline
void postprocess_prob(Problem& prob,
                      ProblemState& prob_state,
                      double* const Cstr_vecds_lat
                     )
{
    for (int i_hstate = 0; i_hstate < kinNH; i_hstate++) {
        prob_state.h_state[i_hstate] = prob_state.h_state_u[i_hstate];
    }

    double pl_disipation_rate = 0.0;
    double effective_shear_rate = 0.0;

    prob.get_slip_contribution(pl_disipation_rate, effective_shear_rate,
                               prob_state.gdot, prob_state.e_vecd_u);

    prob_state.eps_dot = effective_shear_rate;
    prob_state.eps += prob_state.eps_dot * prob_state.dt;
    //
    {
        double dEff = vecd_Deff(prob_state.d_vecd_sm);
        double flow_strength = prob.getHdnScale();
        if (dEff > idp_tiny_sqrt) {
            flow_strength = pl_disipation_rate / dEff;
        }
        prob_state.flow_strength = flow_strength;
    }
        // get Cauchy stress
        //
        prob.elastNEtoC(Cstr_vecds_lat, prob_state.e_vecd_u);
}


__ecmech_hdev__
template<class ProblemState, class ThermoElastN>
inline
void postprocess(ProblemState& prob_state,
                 const ThermoElastN& elastN,
                 const double* const d_svec_kk_sm,
                 double* const sdd,
                 double* const eInt,
                 double* const Cstr_vecds_lat,
                 double eDevTot,
                 double halfVMidDt
                )
{
    double C_matx[ecmech::ndim * ecmech::ndim];
    quat_to_tensor(C_matx, prob_state.quat_u);
    //
    double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
    get_rot_mat_vecd(qr5x5_ls, C_matx);
    //
    double Cstr_vecds_sm[ecmech::nsvec];
    vecsVMa<ntvec>(Cstr_vecds_sm, qr5x5_ls, Cstr_vecds_lat);
    Cstr_vecds_sm[iSvecS] = Cstr_vecds_lat[iSvecS];
    //
    // put end-of-step stress in stressSvecP
    vecdsToSvecP(prob_state.stressSvecP, Cstr_vecds_sm);
    //
    // and now the second half of the trapezoidal integration
    //
    eDevTot += halfVMidDt * vecsInnerSvecDev(prob_state.stressSvecP, d_svec_kk_sm);

    // adjust sign on quat so that as close as possible to quat_o;
    // more likely to keep orientations clustered this way;
    // this flip through the origin is equivalent under antipodal symmetry
    //
    if (vecsyadotb<qdim>(prob_state.quat_u, prob_state.quat_n) < zero) {
        for (int iQ = 0; iQ < ecmech::qdim; ++iQ) {
            prob_state.quat_u[iQ] = -prob_state.quat_u[iQ];
        }
    }

    {
        double gmod = elastN.getGmod(prob_state.tkelv, prob_state.pEOS, prob_state.eNew);
        sdd[i_sdd_bulk] = prob_state.bulkNew;
        sdd[i_sdd_gmod] = gmod;
    }
#ifdef ECMECH_DEBUG
    assert(ecmech::nsdd == 2);
#endif

    prob_state.eNew = prob_state.eNew + eDevTot;
    //
    // could update pressure and temperature again, but do not bother

    eInt[ecmech::i_ne_total] = prob_state.eNew;
#ifdef ECMECH_DEBUG
    assert(ecmech::ne == 1);
#endif
}

/*
* for steady-flow capability, might want to check out Dlsmm_getEnabled() stuff in EvpC.c
*
* convention for spin coming in should be consistent with w_veccp_sm convention
*/
template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
__ecmech_hdev__
inline
bool getResponseSngl(const SlipGeom& slipGeom,
                     const Kinetics& kinetics,
                     const ThermoElastN& elastN,
                     const EosModel& eos,
                     const double dt,
                     const double tolerance,
                     const double* const d_svec_kk_sm, // defRate,
                     const double* const w_veccp_sm, // spin
                     const double* const volRatio,
                     double* const eInt,
                     double* const stressSvecP,
                     double* const hist,
                     double& tkelv,
                     double* const sdd,
                     double* const mtanSD,
                     int outputLevel = 0)
{
    auto prob_state = ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModel>(hist, stressSvecP, tkelv, d_svec_kk_sm, w_veccp_sm, volRatio, dt);

    double halfVMidDt, eDevTot;
    preprocess(slipGeom, kinetics, eos, elastN, volRatio, eInt, d_svec_kk_sm, prob_state, halfVMidDt, eDevTot);

    double Cstr_vecds_lat[ecmech::nsvec];
    {
        EvptnUpdstProblem prob(slipGeom, kinetics, elastN, prob_state);

        // Solver update of things
        {
            snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
            bool status = main_problem(tolerance, solver, outputLevel);

            if (!status) {
                return false;
            }

            if (mtanSD != nullptr) {
                computeTangentStiffness(prob, solver, prob_state, mtanSD);
            }
            // store updated state
            //
            prob.stateFromX(prob_state.e_vecd_u, prob_state.quat_u, solver._x);
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations
        }
        postprocess_prob<Kinetics::nH>(prob, prob_state, Cstr_vecds_lat);
    }
    postprocess(prob_state, elastN, d_svec_kk_sm, sdd, eInt, Cstr_vecds_lat, eDevTot, halfVMidDt);
    return true;
} // getResponseSngl

#if defined(ECMECH_EXTRA_SOLVERS)
/*
* for steady-flow capability, might want to check out Dlsmm_getEnabled() stuff in EvpC.c
*
* convention for spin coming in should be consistent with w_veccp_sm convention
*/
template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
__ecmech_hdev__
inline
bool getResponseNRSngl(
                     const SlipGeom& slipGeom,
                     const Kinetics& kinetics,
                     const ThermoElastN& elastN,
                     const EosModel& eos,
                     const double dt,
                     const double tolerance,
                     const double* const d_svec_kk_sm, // defRate,
                     const double* const w_veccp_sm, // spin
                     const double* const volRatio,
                     double* const eInt,
                     double* const stressSvecP,
                     double* const hist,
                     double& tkelv,
                     double* const sdd,
                     double* const mtanSD,
                     int outputLevel = 0)
{
    auto prob_state = ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModel>(hist, stressSvecP, tkelv, d_svec_kk_sm, w_veccp_sm, volRatio, dt);

    double halfVMidDt, eDevTot;
    preprocess<SlipGeom, Kinetics, EosModel, ThermoElastN, decltype(prob_state), true>(slipGeom, kinetics, eos, elastN, volRatio, eInt, d_svec_kk_sm, prob_state, halfVMidDt, eDevTot);

    double Cstr_vecds_lat[ecmech::nsvec];
    {
        EvptnNRUpdstProblem prob(slipGeom, kinetics, elastN, prob_state);

        // Solver update of things
        {
            snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
            bool status = main_problem(tolerance, solver, outputLevel);

            if (!status) {
                return false;
            }

            if (mtanSD != nullptr) {
                computeTangentStiffness(prob, solver, prob_state, mtanSD);
            }
            // store updated state
            //
            prob.stateFromX(prob_state.e_vecd_u, solver._x);
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations
        }
        postprocess_prob<Kinetics::nH>(prob, prob_state, Cstr_vecds_lat);
    }
    postprocess(prob_state, elastN, d_svec_kk_sm, sdd, eInt, Cstr_vecds_lat, eDevTot, halfVMidDt);
    return true;
} // getResponseSngl
#endif

}
}