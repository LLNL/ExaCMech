#pragma once

#include "ECMech_core.h"
#include "ECMech_util.h"

#include "SNLS_lup_solve.h"

namespace ecmech {
namespace evptn {

    __ecmech_hdev__
    inline
    void get_xtal_frame_vel_grad_terms(double* const def_rate_d5_xtal,
                                        double* const spin_vec_xtal,
                                        const double* const def_rate_d5_sample,
                                        const double* const spin_vec_sample,
                                        const double* const xtal_rmat,
                                        const double* const xtal_rot_mat5)
    {
        vecsVMTa<ecmech::ntvec>(def_rate_d5_xtal, xtal_rot_mat5, def_rate_d5_sample);
        vecsVMTa<ecmech::ndim>(spin_vec_xtal, xtal_rmat, spin_vec_sample);
    }

    template<class SlipGeom, class SlipKinetics>
    __ecmech_hdev__
    inline
    void get_slip_rate_terms(double* const dgdot_dtau,
                            double* const plastic_def_rate_d5,
                            double* const plastic_spin_vec,
                            const double* const kirchoff,
                            const double* const kinetic_values,
                            const SlipGeom& slip_geom,
                            const SlipKinetics& slip_kinetics
                            )
    {        
        if constexpr (SlipGeom::nslip > 0) {
        // default initialize everything to 0.0
        constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip + SlipGeom::nSlipExtra) : SlipGeom::nslip;
        double abs_resolved_shear_stress[nslip_dyn] = {};
        double gdot[SlipGeom::nslip] = {};
        // resolve stress onto slip systems
        // CALL resolve_tau_a_n(crys%tmp4_slp, s_meas%kirchoff, crys)
        //vecsVaTM<ntvec, SlipGeom::nslip>(taua, kirchoff, slipP);
        slip_geom.evalRSS(abs_resolved_shear_stress, kirchoff, slip_geom.getP());
        if constexpr (SlipGeom::dynamic) {
            slip_geom.getExtras(&abs_resolved_shear_stress[SlipGeom::nslip]);
        }
        //
        // CALL plaw_eval(plastic_def_rate_d5, plastic_spin_vec, gss, crys, tkelv, ierr)
        // chi values are passed within extended taua array
        slip_kinetics.evalGdots(gdot, dgdot_dtau, abs_resolved_shear_stress, kinetic_values);
        
        //
        // CALL sum_slip_def(plastic_def_rate_d5, plastic_spin_vec, crys%tmp1_slp, crys) ;
        vecsVMa<ntvec, SlipGeom::nslip>(plastic_def_rate_d5, slip_geom.getP(), gdot);
        vecsVMa<nwvec, SlipGeom::nslip>(plastic_spin_vec, slip_geom.getQ(), gdot);
        }
    }

    // This function performs the necessary chain rules to go from the:
    // dgammadot_dRSS -> dDp_hat_delast_strain
    // dgammadot_dRSS -> dWp_hat_delast_strain
    // terms used typically in either the Jacobian or material tangent stiffness matrix
    template<class SlipGeom, class ThermoElastN>
    __ecmech_hdev__
    inline
    void get_slip_rate_deriv_terms(double* const dDp_hat_delast_strain,
                                    double* const dWp_hat_delast_strain,
                                    const double* const dgdot_dtau,
                                    const double inv_a_vol,
                                    const SlipGeom& slip_geom,
                                    const ThermoElastN& thermoElastN
                                )
    {
        if constexpr (SlipGeom::nslip > 0) {
        double dtaua_deps[ ecmech::ntvec * SlipGeom::nslip ];
        thermoElastN.multDTDepsT(dtaua_deps, slip_geom.getP(), inv_a_vol, SlipGeom::nslip);

        double dgdot_deps[ ecmech::ntvec * SlipGeom::nslip ];
        for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            for (int iSlip = 0; iSlip < SlipGeom::nslip; ++iSlip) {
                int ijThis = ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, SlipGeom::nslip);
                dgdot_deps[ijThis] = dgdot_dtau[iSlip] * dtaua_deps[ijThis];
            }
        }
        vecsMABT<ntvec, SlipGeom::nslip>(dDp_hat_delast_strain, slip_geom.getP(), dgdot_deps);
        vecsMABT<nwvec, ntvec, SlipGeom::nslip>(dWp_hat_delast_strain, slip_geom.getQ(), dgdot_deps);
        }
    }

    // The values returned here are usually useful for post-processing and might have application in
    // application codes. However, they are not really state variables in that everything can be
    // calculated post-state variable update. 
    template<class SlipGeom, class SlipKinetics, class Elasticty>
    __ecmech_hdev__
    inline
    void get_slip_contributions(double& pl_disipation_rate,
                                double& effective_shear_rate,
                                double* const gdot,
                                const double inv_det_v_e,
                                const double* const elast_strain,
                                const double* const kinetic_values,
                                const SlipGeom& slip_geom,
                                const SlipKinetics& slip_kinetics,
                                const Elasticty& elasticity
                                )
    {
        pl_disipation_rate = 0.0;
        effective_shear_rate = 0.0;
        if constexpr (SlipGeom::nslip > 0) {
        // default initialize everything to 0.0 
        constexpr size_t nslip_dyn = (SlipGeom::dynamic) ? (SlipGeom::nslip + SlipGeom::nSlipExtra) : SlipGeom::nslip;
        double abs_resolved_shear_stress[nslip_dyn] = {};
        double junk[SlipGeom::nslip] = {};
        double kirchoff[ecmech::nsvec] = {};
        elasticity.elast_strain_to_kirchoff_stress(kirchoff, elast_strain);
        // resolve stress onto slip systems
        if constexpr (SlipGeom::dynamic) {
            slip_geom.getExtras(&abs_resolved_shear_stress[SlipGeom::nslip]);
        }
        slip_geom.evalRSS(abs_resolved_shear_stress, kirchoff, slip_geom.getP());
        slip_kinetics.evalGdots(gdot, junk, abs_resolved_shear_stress, kinetic_values);
#if defined(ECMECH_USE_DPEFF)
        double plastic_def_rate_d5[ntvec] = {};
        vecsVMa<ntvec, SlipGeom::nslip>(plastic_def_rate_d5, slip_geom.getP(), gdot);
        effective_shear_rate = vecd_Deff(plastic_def_rate_d5);
#else
        effective_shear_rate = vecsssumabs<SlipGeom::nslip>(gdot);
#endif
        pl_disipation_rate = inv_det_v_e * vecsyadotb<SlipGeom::nslip>(abs_resolved_shear_stress, gdot);
        }
    }

    // These terms are commonly used in the residual and jacobian terms related to omega
    // We should probably try to create better names for these outputs...
    __ecmech_hdev__
    inline
    void elasticity_higher_order_terms(double* const A_e_M35,
                                        double* const ee_spin_vec,
                                        double& ee_fac,
                                        const double inv_a_vol,
                                        const double* const elast_d5,
                                        const double* const elast_dt_d5
                                    )
    {
        // from e edot product term in spin (formerly neglected)
        M35_d_AAoB_dA(A_e_M35, elast_d5);
        vecsVMa<nwvec, ntvec>(ee_spin_vec, A_e_M35, elast_dt_d5);
        ee_fac = onehalf * inv_a_vol * inv_a_vol;
    }

    // Quite a few of these variables could be calculated on the fly
    // However, we already calculate most of them as part of the computeRJ portion of things
    // so just do it again here...
    // Might be able to rework this in a better way at some point...
    template<class ThermoElastN, size_t JAC_SIZE, size_t ind_sub_omega>
    __ecmech_hdev__
    inline
    void get_material_tangent_stiffness(double* const material_tangent,
                                        const double* const jacobian,
                                        const double* const dquat_domega_t, 
                                        const double* const rmat_5x5_sample2xtal,
                                        const double* const quat,
                                        const double* const rmat,
                                        const double* const cauchy_stress,
                                        const double inv_det_v_e,
                                        const double inv_a_vol,
                                        const ThermoElastN& thermo_elast_n
                                        )
    {
        // Only concerned with the dCauchy/dDefRate at this point in time
        // and more specifically only considering the deviatoric portion 
        constexpr int nRHS = ecmech::ntvec;

        // mtan_sample_frame
        // = d/dDefRate_sample(Cauchy) = d/dDefRate(Q * alpha * (C_{elas} : lattice_strain)) 
        // Q is the 5x5 rotation oper from crystal to sample 
        // alpha is necessary scaling from Kirchoff to Cauchy
        // C_{elas} is the elasticity tensor (deviatoric contributions)
        // = d/dDefRate_s (Cauchy) = d(Q Cauchy) / dOmega_c * dOmega / dDefRate_s 
        //   + Q * d(Cauchy)/delast_strain * delast_strain / dDefRate_s
        // From our Jacobian we can calculate the dOmega / dDefRate_s and delast_strain / dDefRate_s terms
        // The other ones are either simple to calculate or require some math...
        //
        // dstrainomega_ddef_rate_t => [dlat_strain_ddef_rate_sample; domega_ddef_rate_sample];
        // Initially we set it to be our RHS
        double dstrainomega_ddef_rate_t[ nRHS * JAC_SIZE ] = {}; // transpose for use in SNLS_LUP_SolveX !
        {
        // RHS calculations
        //
        // dstrainomega_ddef_rate_t => dResidual_ddef_rate_sample
        {
            // negatives cancel
            // dstrainomega_ddef_rate_t[0:ind_omega_vec,:] = qr5x5_c2s
            for (int jE = 0; jE < nRHS; ++jE) {
                for (int iE = 0; iE < ntvec; ++iE) { // ntvec, _not_ nDimSys // iE is same as index 
                    dstrainomega_ddef_rate_t[ECMECH_NM_INDX(jE, iE, nRHS, JAC_SIZE)] = rmat_5x5_sample2xtal[ECMECH_NN_INDX(jE, iE, ntvec)];
                }
            }
            // If we had hardening terms then we'd add those here as well...
            // dRhard_ddef_rate but typically we can but for most problems safe to treat that
            // set of terms as being = 0
        }
        // Now solve for our dstrain_ddefrate and domega_ddefrate terms
        int err = SNLS_LUP_SolveX<JAC_SIZE>(const_cast<double*>(jacobian), dstrainomega_ddef_rate_t, nRHS);
        if (err != 0) {
            ECMECH_FAIL(__func__, "error from SNLS_LUP_SolveX");
        }
        }

        {
        double temp_M6[ ecmech::nsvec2 ];
        for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
            for (int jTvec = iTvec; jTvec<nRHS; ++jTvec) {
                const int offset1 = ECMECH_NM_INDX(iTvec, jTvec, nRHS, JAC_SIZE);
                const int offset2 = ECMECH_NM_INDX(jTvec, iTvec, nRHS, JAC_SIZE);
                const double tmp = dstrainomega_ddef_rate_t[offset1];
                dstrainomega_ddef_rate_t[offset1] = dstrainomega_ddef_rate_t[offset2];
                dstrainomega_ddef_rate_t[offset2] = tmp;
            }
        }
        thermo_elast_n.template multCauchyDif<nRHS, JAC_SIZE>(temp_M6, dstrainomega_ddef_rate_t, inv_det_v_e, inv_a_vol);
        // Apply final rotation
        qr6x6_pre_mul<ecmech::nsvec, false>(material_tangent, temp_M6, rmat_5x5_sample2xtal);
        }

        // Calculate the d(QCauchy) / dDefRate_s term now and add that to
        {
        double dcauchy_dquat[ ecmech::ntvec * ecmech::qdim ];
        {
            double drmat_dquat[ ecmech::ndim * ecmech::ndim * ecmech::qdim ];
            d_quat_to_tensor(drmat_dquat, quat);
            double dcauchy_drmat[ ecmech::ntvec * ecmech::ndim * ecmech::ndim ];
            d_rot_mat_vecd_smop(dcauchy_drmat, rmat, cauchy_stress);
            vecsMAB<ntvec, qdim, ndim*ndim>(dcauchy_dquat, dcauchy_drmat, drmat_dquat);
        }

        // We now need to be able to go from our domega_ddef_rate to dquat_ddef_rate
        // dquat_ddef_rate = dquat_domega_t * domega_ddef_rate_t
        double dquat_ddef_rate[ ecmech::qdim * nRHS ];
        for (int ii_I = 0; ii_I < nRHS; ++ii_I) {
            for (int ii_Q = 0; ii_Q < ecmech::qdim; ++ii_Q) {
                int iiQI = ECMECH_NM_INDX(ii_Q, ii_I, ecmech::qdim, nRHS);
                dquat_ddef_rate[iiQI] = 0.0;
                for (int ii_W = 0; ii_W < ecmech::nwvec; ++ii_W) {
                    dquat_ddef_rate[iiQI] +=
                    dquat_domega_t[ECMECH_NM_INDX(ii_W, ii_Q, ecmech::nwvec, ecmech::qdim)] *
                    dstrainomega_ddef_rate_t[ECMECH_NM_INDX(ii_I, ind_sub_omega + ii_W, nRHS, JAC_SIZE)];
                }
            }
        }

        // Now get the dcauchy_lattice_dI terms by doing ->
        // dcauchy_lattice_dI = d_cauchy_lattice_dquat * dRmat_quat_dI
        double dqcauchy_ddefrate[ ecmech::ntvec * nRHS ];
        vecsMAB<ecmech::ntvec, nRHS, ecmech::qdim>(dqcauchy_ddefrate, dcauchy_dquat, dquat_ddef_rate);
        for (int ii_T = 0; ii_T < ecmech::ntvec; ++ii_T) {
            for (int ii_I = 0; ii_I < nRHS; ++ii_I) {
                // NOTE : only looping over ntvec, but mtan_sI is nsvec in the first dimension
                material_tangent[ECMECH_NN_INDX(ii_T, ii_I, ecmech::nsvec)] += dqcauchy_ddefrate[ECMECH_NM_INDX(ii_T, ii_I, ecmech::ntvec, nRHS)];
            }
        }
        }
    }
}
}