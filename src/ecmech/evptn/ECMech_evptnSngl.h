#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_elastic.h"
#include "ECMech_eosSimple.h"
#include "evptn/ECMech_base_classes.h"
#include "evptn/ECMech_base_fcns.h"
#include "evptn/ECMech_evptn.h"

#include "SNLS_TrDLDenseG.h"
#include "SNLS_HybrdTrDLDenseG.h"

namespace ecmech {
namespace evptn {

template<class SlipGeom, class Kinetics, class EosModel, class ThermoElastN, class ProbState, bool RStarSolve=false>
__ecmech_hdev__
inline
bool preprocess(const SlipGeom& slipGeom,
                const Kinetics& kinetics,
                const EosModel& eos,
                const ThermoElastN& thermoElastN,
                const double* const rel_vol_ratios,
                const double* const internal_energy,
                const double* const def_rate_d6v_sample,
                ProbState& prob_state,
                double& halfVMidDt,
                double& dev_strain_energy_total)
{
    // total increment in the deviatoric part of the strain energy using
    // trapezoidal rule integration
    //
    // just beginning-of-step stress part so far
    //
    halfVMidDt = oneqrtr * (rel_vol_ratios[0] + rel_vol_ratios[1]) * prob_state.dt;
    dev_strain_energy_total = halfVMidDt * vecsInnerSvecDev(prob_state.cauchy_stress_d6p, def_rate_d6v_sample);

    // EOS
    //
    const double energy_old = internal_energy[ecmech::i_ne_total];
    //
    // get tkelv from beginning-of-step to avoid tangent stiffness contributions
    {
        double pressure_BOS;
        const double rel_vol_old = rel_vol_ratios[0];
        eos.evalPT(pressure_BOS, prob_state.tkelv, rel_vol_old, energy_old);
    }

    {
        const double pressure_old = prob_state.cauchy_stress_d6p[6];
        double tkelv_new, dpde, dpdv, dtde;
        updateSimple(eos, prob_state.pressure_EOS, tkelv_new, prob_state.energy_new, prob_state.bulk_modulus_new,
                     dpde, dpdv, dtde,
                     rel_vol_ratios[1], rel_vol_ratios[3],
                     energy_old, pressure_old);
    }

    // update hardness state to the end of the step
    // gdot is still at beginning-of-step
    //
    constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip) : 1;
    double hvals[nslip_dyn] = {}; // additional values needed to update the hardening state
    if constexpr(SlipGeom::dynamic) {

        double kirchoff[ecmech::nsvec] = {};
        double elast_d5v[ecmech::nsvec] = {};
        double stress_dev6_press[ecmech::nsvec+1] = {};

       double a_vol = pow(prob_state.rel_vol_new, onethird);
       double inv_a_vol = 1.0 / a_vol;

        vecsVxa<ntvec>(elast_d5v, inv_a_vol, prob_state.elast_d5_n);
        elast_d5v[iSvecS] = sqr3 * log(a_vol);
        thermoElastN.eval(kirchoff, elast_d5v, prob_state.tkelv, prob_state.pressure_EOS, prob_state.energy_new);
        vecdsToSvecP(stress_dev6_press, kirchoff);

        // For dynamic slip systems we need the chi angle
        double P[ecmech::ntvec * SlipGeom::nslip];
        double Q[ecmech::nwvec * SlipGeom::nslip];
        // still need to rotate stress state back to original value
        slipGeom.getPQ(hvals, P, Q, stress_dev6_press);
    }
    const int nfevals = kinetics.updateH(prob_state.h_state_u, prob_state.h_state, prob_state.dt, prob_state.gdot, hvals, prob_state.tkelv);
    if (nfevals < 0) {
        ECMECH_WARN(__func__, "Hardening failed to converge");
        return false;
    }
#if defined(ECMECH_EXTRA_SOLVERS)
    if constexpr (RStarSolve) {
        auto prob = RotUpdProblem(slipGeom, thermoElastN, prob_state);
         // update Rstar aka Rdot * dt
         // gdot is still at beginning-of-step
        snls::SNLSTrDlDenseG<decltype(prob)> solver(prob);
        const bool status = main_problem(1e-8, solver, 0);
        if (!status) {
            ECMECH_WARN(__func__, "RStar solver failed to converge");
            return false;
        }
        prob.stateFromX(prob_state.quat_u, solver._x);
    }
#endif
    return true;
}

template<class SNLS_Solver>
__ecmech_hdev__
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

template<class Problem, class Solver, class ProblemState>
__ecmech_hdev__
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
    // neglects effect of pressure_EOS and rel_vol_new on workings of evptn
    //
    mtanSD_vecds[ECMECH_NN_INDX(iSvecS, iSvecS, ecmech::nsvec)] = three * prob_state.bulk_modulus_new;

    // convert from vecds notation to svec notation
    //
    mtan_conv_sd_svec<true>(mtanSD, mtanSD_vecds);
}

template<int kinNH, class Problem, class ProblemState>
__ecmech_hdev__
inline
void postprocess_prob(Problem& prob,
                      ProblemState& prob_state,
                      double* const cauchy_stress_d5p_xtal
                     )
{
    for (int i_hstate = 0; i_hstate < kinNH; i_hstate++) {
        prob_state.h_state[i_hstate] = prob_state.h_state_u[i_hstate];
    }

    double pl_disipation_rate = 0.0;
    double effective_shear_rate = 0.0;

    prob.get_slip_contribution(pl_disipation_rate, effective_shear_rate,
                               prob_state.gdot, prob_state.elast_d5_u);

    prob_state.eps_dot = effective_shear_rate;
    prob_state.eps += prob_state.eps_dot * prob_state.dt;
    //
    {
        double dEff = vecd_Deff(prob_state.def_rate_d5_sample);
        double flow_strength = prob.getHdnScale();
        if (dEff > idp_tiny_sqrt) {
            flow_strength = pl_disipation_rate / dEff;
        }
        prob_state.flow_strength = flow_strength;
    }
        // get Cauchy stress
        //
        prob.elastNEtoC(cauchy_stress_d5p_xtal, prob_state.elast_d5_u);
}


template<class ProblemState, class ThermoElastN>
__ecmech_hdev__
inline
void postprocess(ProblemState& prob_state,
                 const ThermoElastN& elastN,
                 const double* const def_rate_d6v_sample,
                 double* const sdd,
                 double* const internal_energy,
                 double* const cauchy_stress_d5p_xtal,
                 double dev_strain_energy_total,
                 double halfVMidDt
                )
{
    double xtal_rmat[ecmech::ndim * ecmech::ndim];
    quat_to_tensor(xtal_rmat, prob_state.quat_u);
    //
    double rmat_5x5_sample2xtal[ecmech::ntvec * ecmech::ntvec];
    get_rot_mat_vecd(rmat_5x5_sample2xtal, xtal_rmat);
    //
    double cauchy_stress_dev_press_sample[ecmech::nsvec];
    vecsVMa<ntvec>(cauchy_stress_dev_press_sample, rmat_5x5_sample2xtal, cauchy_stress_d5p_xtal);
    cauchy_stress_dev_press_sample[iSvecS] = cauchy_stress_d5p_xtal[iSvecS];
    //
    // put end-of-step stress in cauchy_stress_d6p
    vecdsToSvecP(prob_state.cauchy_stress_d6p, cauchy_stress_dev_press_sample);
    //
    // and now the second half of the trapezoidal integration
    //
    dev_strain_energy_total += halfVMidDt * vecsInnerSvecDev(prob_state.cauchy_stress_d6p, def_rate_d6v_sample);

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
        double shear_modulus = elastN.getGmod(prob_state.tkelv, prob_state.pressure_EOS, prob_state.energy_new);
        sdd[i_sdd_bulk] = prob_state.bulk_modulus_new;
        sdd[i_sdd_gmod] = shear_modulus;
    }
#ifdef ECMECH_DEBUG
    assert(ecmech::nsdd == 2);
#endif

    prob_state.energy_new = prob_state.energy_new + dev_strain_energy_total;
    //
    // could update pressure and temperature again, but do not bother

    internal_energy[ecmech::i_ne_total] = prob_state.energy_new;
#ifdef ECMECH_DEBUG
    assert(ecmech::ne == 1);
#endif
}

/*
* for steady-flow capability, might want to check out Dlsmm_getEnabled() stuff in EvpC.c
*
* convention for spin coming in should be consistent with spin_vec_sample convention
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
                     const double* const def_rate_d6v_sample, // defRate,
                     const double* const spin_vec_sample, // spin
                     const double* const rel_vol_ratios,
                     double* const internal_energy,
                     double* const cauchy_stress_d6p,
                     double* const hist,
                     double& tkelv,
                     double* const sdd,
                     double* const mtanSD,
                     int outputLevel = 0)
{
    auto prob_state = ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModel>(hist, cauchy_stress_d6p, tkelv, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios, dt);

    double halfVMidDt, dev_strain_energy_total;
    const bool pre_status = preprocess(slipGeom, kinetics, eos, elastN, rel_vol_ratios, internal_energy, def_rate_d6v_sample, prob_state, halfVMidDt, dev_strain_energy_total);

    if (!pre_status) { return false; }

    double cauchy_stress_d5p_xtal[ecmech::nsvec];
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
            prob.stateFromX(prob_state.elast_d5_u, prob_state.quat_u, solver._x);
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations
        }
        postprocess_prob<Kinetics::nH>(prob, prob_state, cauchy_stress_d5p_xtal);
    }
    postprocess(prob_state, elastN, def_rate_d6v_sample, sdd, internal_energy, cauchy_stress_d5p_xtal, dev_strain_energy_total, halfVMidDt);
    return true;
} // getResponseSngl

#if defined(ECMECH_EXTRA_SOLVERS)
/*
* for steady-flow capability, might want to check out Dlsmm_getEnabled() stuff in EvpC.c
*
* convention for spin coming in should be consistent with spin_vec_sample convention
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
                     const double* const def_rate_d6v_sample, // defRate,
                     const double* const spin_vec_sample, // spin
                     const double* const rel_vol_ratios,
                     double* const internal_energy,
                     double* const cauchy_stress_d6p,
                     double* const hist,
                     double& tkelv,
                     double* const sdd,
                     double* const mtanSD,
                     int outputLevel = 0)
{
    auto prob_state = ProblemState<SlipGeom, Kinetics, ThermoElastN, EosModel>(hist, cauchy_stress_d6p, tkelv, def_rate_d6v_sample, spin_vec_sample, rel_vol_ratios, dt);

    double halfVMidDt, dev_strain_energy_total;
    preprocess<SlipGeom, Kinetics, EosModel, ThermoElastN, decltype(prob_state), true>(slipGeom, kinetics, eos, elastN, rel_vol_ratios, internal_energy, def_rate_d6v_sample, prob_state, halfVMidDt, dev_strain_energy_total);

    double cauchy_stress_d5p_xtal[ecmech::nsvec];
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
            prob.stateFromX(prob_state.elast_d5_u, solver._x);
            //
            hist[iHistA_nFEval] = solver.getNFEvals(); // does _not_ include updateH iterations
        }
        postprocess_prob<Kinetics::nH>(prob, prob_state, cauchy_stress_d5p_xtal);
    }
    postprocess(prob_state, elastN, def_rate_d6v_sample, sdd, internal_energy, cauchy_stress_d5p_xtal, dev_strain_energy_total, halfVMidDt);
    return true;
} // getResponseSngl
#endif

}
}