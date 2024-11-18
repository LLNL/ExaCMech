#pragma once

#include <cassert>

#include "ECMech_core.h"
#include "ECMech_util.h"

#include "RAJA/RAJA.hpp"

namespace ecmech {
namespace evptn {

    const int numHistAux = 4; // effective shearing rate, accumulated shear, flow strength, nFEval
    //
    const int iHistLbA = 0;
    const int iHistA_shrateEff = iHistLbA + 0;
    const int iHistA_shrEff = iHistLbA + 1;
    const int iHistA_flowStr = iHistLbA + 2;
    const int iHistA_nFEval = iHistLbA + 3;
    const int iHistLbE = numHistAux;
    const int iHistLbQ = numHistAux + ecmech::ntvec;
    const int iHistLbH = numHistAux + ecmech::ntvec + ecmech::qdim;

    /*
    * just a container for a traits
    */
    template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
    class NumHist
    {
        public:
        // see n_rsv_matmod in F90 code
        static constexpr int iHistLbGdot = iHistLbH + Kinetics::nH;
        static constexpr int numHist = iHistLbH + Kinetics::nH + SlipGeom::nslip;
    }; // NumHist

    // These are largely things that we need to persist between function calls / what we want to pass around various function calls
    template<class SlipGeom, class Kinetics, class ThermoElastN, class EosModel>
    struct ProblemState
    {
        static constexpr int iHistLbGdot = NumHist<SlipGeom, Kinetics, ThermoElastN, EosModel>::iHistLbGdot;

        /// Hardening state at beg-of-time / end-of-time step
        double* const h_state;
        /// shearing rate at beg-of-time / end-of-time step
        double* const gdot;
        /// elastic strain as deviatoric 5 vector at end-of-time step in crystal frame
        double* const elast_d5_u;
        /// lattice orientation quaternion that maps xtal2sample frame at end-of-time step
        double* const quat_u;
        /// equivalent plastic strain rate
        double& eps_dot;
        /// equivalent plastic strain
        double& eps;
        double& flow_strength;
        /// cauchy stress as deviatoric 6 vector with pressure in sample frame
        double* const cauchy_stress_d6p;
        /// spin vector of the velocity gradient in sample frame
        const double* const spin_vec_sample;
        /// Depending on the context this will be seen as either the
        /// relative volume change when in use cases like EOS evaluations
        /// or when dealing with the elasticity equations this is the
        /// $det(\mathbf{V}^e)$
        const double rel_vol_new;
        /// delta time
        const double dt;
        /// temperature as typically given in Kelvins
        double& tkelv;

        /// deformation rate as the deviatoric 5 vector in the sample frame
        double def_rate_d5_sample[ecmech::ntvec];
        /// elastic strain as deviatoric 5 vector at beg-of-time step in crystal frame
        double elast_d5_n[ecmech::ntvec];
        /// lattice orientation quaternion that maps xtal2sample frame at beg-of-time step
        double quat_n[ecmech::qdim];
        /// Hardening state at end-of-time step
        double h_state_u[Kinetics::nH];
        /// Values determined by the EOS (equation of state)
        /// our new pressure, energy from the EOS, and bulk modulus
        double pressure_EOS, energy_new, bulk_modulus_new;

        __ecmech_hdev__
        ProblemState(double* const hist, double* const cauchy_stress_d6p,
                    double& tkelv,
                    const double* const def_rate_d6v_sample,
                    const double* const spin_vec_sample,
                    const double* const rel_vol_ratios,
                    const double dt) :
        h_state(&(hist[iHistLbH])),
        gdot(&(hist[iHistLbGdot])),
        elast_d5_u(&(hist[iHistLbE])),
        quat_u(&(hist[iHistLbQ])),
        eps_dot(hist[iHistA_shrateEff]),
        eps(hist[iHistA_shrEff]),
        flow_strength(hist[iHistA_flowStr]),
        cauchy_stress_d6p(cauchy_stress_d6p),
        spin_vec_sample(spin_vec_sample),
        rel_vol_new(rel_vol_ratios[1]),
        dt(dt),
        tkelv(tkelv)
        {
            // convert deformation rate convention
            //
            // double def_rate_d5_sample[ecmech::ntvec];
            svecToVecd(def_rate_d5_sample, def_rate_d6v_sample);
            //
            // copies, to keep beginning-of-step state safe
            //
            for (int i_hist = 0; i_hist < ecmech::ntvec; i_hist++) {
                elast_d5_n[i_hist] = hist[iHistLbE + i_hist];
            }

            for (int i_hist = 0; i_hist < ecmech::qdim; i_hist++) {
                quat_n[i_hist] = hist[iHistLbQ + i_hist];
            }
            //
            // normalize quat just in case
            vecsVNormalize<qdim>(quat_n);
        }

        ~ProblemState() = default;
    };

    template<class ThermoElastN>
    class EvptnLatticeStrainProblem
    {
        public:
        static constexpr int nDimSys = ecmech::ntvec;

        __ecmech_hdev__
        EvptnLatticeStrainProblem(const ThermoElastN& thermoElastN,
                                const double dt,
                                const double det_v_e, 
                                const double energy_vol_ref, 
                                const double pressure_EOS, 
                                const double tkelv,
                                const double* const elast_d5_n)
        : m_thermo_elast_n(thermoElastN),
        m_dt(dt), m_det_v_e(det_v_e), m_energy_vol_ref(energy_vol_ref),
        m_pressure_EOS(pressure_EOS), m_tkelv(tkelv),
        m_elast_d5_n(elast_d5_n),
        m_inv_dt(1.0 / dt),
        m_inv_det_v_e(1.0 / m_det_v_e),
        m_a_vol(pow(m_det_v_e, onethird)),
        m_inv_a_vol(1.0 / m_a_vol)
        {}

        ~EvptnLatticeStrainProblem() = default;

        // used to be elastNEtoT
        __ecmech_hdev__
        inline
        void elast_strain_to_kirchoff_stress(double* const kirchoff_stress, // nsvec
                                            const double* const elast_d5 // ntvec
                                        ) const
        {
        //// do not need to use elaw_T_BT here as T and BT are the same
        //
        // specialize to cem%l_lin_lnsd
        // CALL elawn_T(s_meas, elast_d5_f, crys%elas, tkelv, .TRUE., a_V, &
        // & pressure_EOS, energy_vol_ref, crys%i_eos_model, crys%eos_const &
        // &)
        double elast_d5v[ecmech::nsvec];
        vecsVxa<ntvec>(elast_d5v, m_inv_a_vol, elast_d5);
        //// tr_Ee = three * DLOG(a_V%r)
        //// CALL trace_to_vecds_s(s_meas%elast_dev_press_vec(SVEC), tr_Ee)
        elast_d5v[iSvecS] = sqr3 * log(m_a_vol); // could go into constructor
        //
        //// Kirchhoff stress from elast_d5v
        // CALL elawn_lin_op(s_meas%kirchoff, s_meas%elast_dev_press_vec, cem, tkelv, &
        // & pressure_EOS, energy_vol_ref, i_eos_model, eos_const)
        m_thermo_elast_n.eval(kirchoff_stress, elast_d5v, m_tkelv, m_pressure_EOS, m_energy_vol_ref);
        }

        // used to be elastNEtoC
        __ecmech_hdev__
        inline
        void elast_strain_to_cauchy_stress(double* const cauchy, // nsvec
                                        const double* const elast_d5_f // ntvec
                                        ) const
        {
        double kirchoff[ecmech::nsvec];
        this->elast_strain_to_kirchoff_stress(kirchoff, elast_d5_f);
        m_thermo_elast_n.getCauchy(cauchy, kirchoff, m_inv_det_v_e);
        }

        template<bool calc_strain_rate = false>
        __ecmech_hdev__
        inline
        void get_elast_strain_state(double* const elast_delta_d5,
                                double* const elast_dt_d5,
                                const double* const x) const
        {
        //////////////////////////////
        // PULL VALUES out of x, with scalings
        //
        // double elast_dt_d5[ecmech::ntvec];
        vecsVxa<ntvec>(elast_dt_d5, ecmech::e_scale, x); // elast_dt_d5 is now the delta, _not_ yet elast_dt_d5
        // elast_d5_f is end-of-step
        // double elast_d5_f[ntvec];
        vecsVapb<ntvec>(elast_delta_d5, elast_dt_d5, m_elast_d5_n);
        if constexpr(calc_strain_rate) {
            vecsVsa<ntvec>(elast_dt_d5, m_inv_dt); // _now_ elast_dt_d5 has dt contributions
        }
        }

        /*
        * NOTES :
        * () should be equivalent to what happens in get_elast_strain_state<false>
        * () not necessarily safe if elast_dev_press_vec is the same memory as _elast_d5_n or quat is the same as _xtal_ori_quat_n
        */
        __ecmech_hdev__
        inline
        void stateFromX(double* const elast_d5,
                    const double* const x) const
        {
        double elast_d5_delta[ecmech::ntvec] = {};
        this->get_elast_strain_state(elast_d5, elast_d5_delta, x);
        }

        __ecmech_hdev__
        inline
        void get_elast_strain_residual(double* const residual,
                                    const double epsdot_scale_inv,
                                    const double* const elast_dt_d5,
                                    const double* const plastic_def_rate_d5,
                                    const double* const def_rate_d5_xtal) const
        {
        for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            residual[m_ind_sub_elas + iTvec] = epsdot_scale_inv * ( // SCALING
                m_inv_a_vol * elast_dt_d5[iTvec] + plastic_def_rate_d5[iTvec] - def_rate_d5_xtal[iTvec]);
        }
        }

        template<size_t JAC_SIZE>
        __ecmech_hdev__
        inline
        void get_deriv_elast_strain_wrt_elast_strain(double* const jacobian,
                                                    const double* const dDp_hat_delast_strain) const
        {
        RAJA::View<double, RAJA::Layout<2>> jacob_ee(jacobian, JAC_SIZE, JAC_SIZE);
        // dislocation plasticity;
        // first contribution; overwrite
        for (size_t jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
            for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                jacob_ee(iTvec, jTvec) = dDp_hat_delast_strain[ECMECH_NN_INDX(iTvec, jTvec, ecmech::ntvec)];
            }
        }
        // elastic rate
        //
        {
            const double adti = m_inv_a_vol * m_inv_dt;
            for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                jacob_ee(iTvec, iTvec) += adti;
            }
        }
        } 

        template<size_t JAC_SIZE, size_t ind_sub_r>
        __ecmech_hdev__
        inline
        void get_deriv_omega_wrt_elast_strain(double* const jacobian,
                                            const double* const elast_dt_d5,
                                            const double elast_elast_factor,
                                            const double* const dWp_hat_delast_strain,
                                            const double* const A_e_M35) const
        {
        // d(B_xi)/d(elast_dev_press_vecs_f)
        //
        RAJA::View<double, RAJA::Layout<2>> jacob_re(jacobian, JAC_SIZE, JAC_SIZE);

        double A_edot_M35[ecmech::nwvec * ecmech::ntvec];
        M35_d_AAoB_dA(A_edot_M35, elast_dt_d5);

        double dt_ee_fac = m_dt * elast_elast_factor;

        for (size_t iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            for (size_t jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                int ijWT = ECMECH_NM_INDX(iWvec, jTvec, ecmech::nwvec, ecmech::ntvec);
                jacob_re(iWvec + ind_sub_r, jTvec) =
                    m_dt * dWp_hat_delast_strain[ijWT] - dt_ee_fac * (A_e_M35[ijWT] * m_inv_dt - A_edot_M35[ijWT]);
            }
        }         
        }

        template<size_t JAC_SIZE, size_t ind_sub_h, size_t num_hard, size_t num_slip>
        __ecmech_hdev__
        inline
        void get_deriv_hardening_wrt_elast_strain(double* const jacobian,
                                                const double* const dhard_dgdot,
                                                const double* const dgdot_delast_strain) const
        {
        // d(B_h) / d(e)
        // jacob_he = dt * (dhdot/dgdot)(dgdot/de)
        // nh x ntvec matrix
        double dhdot_delast_strain[num_hard * ecmech::ntvec];
        vecsMABT<num_hard, ecmech::ntvec, num_slip>(dhdot_delast_strain, dhard_dgdot, dgdot_delast_strain);
        RAJA::View<double, RAJA::Layout<2>> jacob_he(jacobian, JAC_SIZE, JAC_SIZE);
        for (size_t iH = 0; iH < num_hard; ++iH) {
            for (size_t jE = 0; jE < ecmech::ntvec; ++jE) {
                // could also make dhdot_delast_strain into a RAJA view, but not really needed
                jacob_he(iH + ind_sub_h, jE) = -m_dt * dhdot_delast_strain[ECMECH_NM_INDX(iH, jE, num_hard, ecmech::ntvec) ];
            }
        }       
        }

        public:
        static constexpr size_t m_ind_sub_elas = 0; // ntvec end_point
        const ThermoElastN& m_thermo_elast_n;
        const double m_dt, m_det_v_e, m_energy_vol_ref;
        const double m_pressure_EOS, m_tkelv;
        const double* const m_elast_d5_n;
        const double m_inv_dt, m_inv_det_v_e, m_a_vol, m_inv_a_vol;
    };

    template <size_t ind_sub_r=ecmech::ntvec>
    class EvptnLatticeRotationProblem {
        public:
        static constexpr size_t nDimSys = ecmech::nwvec;

        public:
        __ecmech_hdev__
        EvptnLatticeRotationProblem(const double dt,
                                const double* const xtal_ori_quat_n)
        : m_dt(dt), m_xtal_ori_quat_n(xtal_ori_quat_n) {}
        
        ~EvptnLatticeRotationProblem() = default;

        __ecmech_hdev__
        inline
        void get_rotation_state(double* const delta_omega,
                                double* const xtal_rmat,
                                double* const xtal_rot_mat5,
                                const double* const x) const
        {
        vecsVxa<ecmech::nwvec>(delta_omega, ecmech::r_scale, &(x[ind_sub_r]));
        //
        // not done in EvpC :
        // CALL exp_map_cpvec(A, xi_f)
        // CALL get_c(c, A, C_n)
        //
        double xtal_ori_quat_delta[ecmech::qdim];
        double xtal_ori_quat_n1[ecmech::qdim];

        emap_to_quat(xtal_ori_quat_delta, delta_omega);
        get_c_quat(xtal_ori_quat_n1, xtal_ori_quat_delta, m_xtal_ori_quat_n);
        quat_to_tensor(xtal_rmat, xtal_ori_quat_n1);
        get_rot_mat_vecd(xtal_rot_mat5, xtal_rmat);
        }

        /*
        * NOTES :
        * () should be equivalent to what happens in get_elast_strain_state<false>
        * () not necessarily safe if elast_dev_press_vec is the same memory as _elast_d5_n or quat is the same as _xtal_ori_quat_n
        */
        // Assume that x has is at the location we need it to be at... 
        __ecmech_hdev__
        inline
        void stateFromX(double* const xtal_ori_quat,
                    const double* const x) const
        {
        double delta_omega[ecmech::nwvec];
        double xtal_ori_quat_delta[ecmech::qdim];
        vecsVxa<ecmech::nwvec>(delta_omega, ecmech::r_scale, x);
        emap_to_quat(xtal_ori_quat_delta, delta_omega);
        get_c_quat(xtal_ori_quat, xtal_ori_quat_delta, m_xtal_ori_quat_n);
        }

        __ecmech_hdev__
        inline
        void deltaOmegaFromState(double* const delta_omega,
                                 const double* const xtal_ori_quat_n,
                                 const double* const xtal_ori_quat_n1) const
        {
            double xtal_ori_quat_delta[ecmech::qdim] = {};

            quat_rel_rotation(xtal_ori_quat_delta, xtal_ori_quat_n1, xtal_ori_quat_n);
            quat_to_emap(delta_omega, xtal_ori_quat_delta);
        }

        __ecmech_hdev__
        inline
        void get_omega_residual(double* const residual,
                                const double rot_incr_scale_inv,
                                const double ee_fac,
                                const double* const delta_omega,
                                const double* const spin_vec_lat,
                                const double* const plastic_spin_vec,
                                const double* const ee_spin_vec) const
        {
        // RESIDUAL B_omega
        for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            residual[ind_sub_r + iWvec] = rot_incr_scale_inv * // SCALING
                                        (delta_omega[iWvec] - m_dt * (spin_vec_lat[iWvec]
                                        - plastic_spin_vec[iWvec] + ee_fac * ee_spin_vec[iWvec]));
        }
        }

        template<size_t JAC_SIZE>
        __ecmech_hdev__
        inline
        void get_deriv_omega_wrt_omega(double* const jacobian,
                                    const double* const dspin_samp_domega) const
        {
        // d(B_xi)/d(xi_f)
        RAJA::View<double, RAJA::Layout<2>> jacob_rr(jacobian, JAC_SIZE, JAC_SIZE);
        for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            for (int jWvec = 0; jWvec < ecmech::nwvec; ++jWvec) {
                int ijWW = ECMECH_NN_INDX(iWvec, jWvec, ecmech::nwvec);
                jacob_rr(iWvec + ind_sub_r, jWvec + ind_sub_r) = -m_dt * dspin_samp_domega[ijWW];
            }
            jacob_rr(iWvec + ind_sub_r, iWvec + ind_sub_r) += one;
        }
        }

        template<size_t JAC_SIZE>
        __ecmech_hdev__
        inline
        void get_deriv_elast_strain_wrt_omega(double* const jacobian,
                                            const double* const ddef_rate_samp_domega) const
        {
        // d(B_S)/d(xi_f)
        //
        // jacob_er = -dDsm_dxi(:,:)
        RAJA::View<double, RAJA::Layout<2> > jacob_er(jacobian, JAC_SIZE, JAC_SIZE);
        for (int jWvec = 0; jWvec<ecmech::nwvec; ++jWvec) {
            for (int iTvec = 0; iTvec<ecmech::ntvec; ++iTvec) {
                // could also make dDsm_dxi into a RAJA view, but not really needed
                const size_t ind = ECMECH_NM_INDX(iTvec, jWvec, ecmech::ntvec, ecmech::nwvec);
                jacob_er(iTvec, jWvec + ind_sub_r) = -ddef_rate_samp_domega[ind];
            }
        }
        }

        /**
        This method does nothing as most crystal plasticity hardening models have no dependency directly or indirectly with the
        lattice rotation.
        */
        __ecmech_hdev__
        inline
        void get_deriv_hardening_wrt_omega(double* const /* jacobian */){}

        public:
        const double m_dt;
        const double* const m_xtal_ori_quat_n;
    };

}
}