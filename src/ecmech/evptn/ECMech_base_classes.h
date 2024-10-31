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

        double* const h_state;
        double* const gdot;
        double* const elast_d5_u;
        double* const quat_u;
        double& eps_dot;
        double& eps;
        double& flow_strength;
        double* const cauchy_stress_d6p;
        const double* const spin_vec_sample;
        const double rel_vol_new;
        const double dt;
        double& tkelv;


        double def_rate_d5_sample[ecmech::ntvec];
        double elast_d5_n[ecmech::ntvec];
        double quat_n[ecmech::qdim];
        double h_state_u[Kinetics::nH];
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
                                const double det_vol, 
                                const double det_v_e, 
                                const double pressure_EOS, 
                                const double tkelv,
                                const double* const elast_d5_n)
        : m_thermo_elast_n(thermoElastN),
        m_dt(dt), m_det_vol(det_vol), m_det_v_e(det_v_e),
        m_pressure_EOS(pressure_EOS), m_tkelv(tkelv),
        m_elast_d5_n(elast_d5_n),
        m_inv_dt(1.0 / dt),
        m_inv_det_vol(1.0 / m_det_vol),
        m_a_vol(pow(m_det_vol, onethird)),
        m_inv_a_vol(1.0 / m_a_vol)
        {}

        ~EvptnLatticeStrainProblem() = default;

        // used to be elastNEtoT
        __ecmech_hdev__
        inline
        void elas_strain_to_kirchoff_stress(double* const kirchoff_stress, // nsvec
                                            const double* const elas_dev_vec // ntvec
                                        ) const
        {
        //// do not need to use elaw_T_BT here as T and BT are the same
        //
        // specialize to cem%l_lin_lnsd
        // CALL elawn_T(s_meas, elast_dev_vec_f, crys%elas, tkelv, .TRUE., a_V, &
        // & pressure_EOS, det_v_e, crys%i_eos_model, crys%eos_const &
        // &)
        double elas_dev_vol_vec[ecmech::nsvec];
        vecsVxa<ntvec>(elas_dev_vol_vec, m_inv_a_vol, elas_dev_vec);
        //// tr_Ee = three * DLOG(a_V%r)
        //// CALL trace_to_vecds_s(s_meas%elast_dev_press_vec(SVEC), tr_Ee)
        elas_dev_vol_vec[iSvecS] = sqr3 * log(m_a_vol); // could go into constructor
        //
        //// Kirchhoff stress from elas_dev_vol_vec
        // CALL elawn_lin_op(s_meas%kirchoff, s_meas%elast_dev_press_vec, cem, tkelv, &
        // & pressure_EOS, det_v_e, i_eos_model, eos_const)
        m_thermo_elast_n.eval(kirchoff_stress, elas_dev_vol_vec, m_tkelv, m_pressure_EOS, m_det_v_e);
        }

        // used to be elastNEtoC
        __ecmech_hdev__
        inline
        void elas_strain_to_cauchy_stress(double* const cauchy, // nsvec
                                        const double* const elast_dev_vec_f // ntvec
                                        ) const
        {
        double kirchoff[ecmech::nsvec];
        this->elas_strain_to_kirchoff_stress(kirchoff, elast_dev_vec_f);
        m_thermo_elast_n.getCauchy(cauchy, kirchoff, m_inv_det_vol);
        }

        template<bool calc_strain_rate = false>
        __ecmech_hdev__
        inline
        void get_elas_strain_state(double* const elas_delta_dev_vec,
                                double* const elas_dt_dev_vec,
                                const double* const x) const
        {
        //////////////////////////////
        // PULL VALUES out of x, with scalings
        //
        // double elas_dt_dev_vec[ecmech::ntvec];
        vecsVxa<ntvec>(elas_dt_dev_vec, ecmech::e_scale, x); // elas_dt_dev_vec is now the delta, _not_ yet elas_dt_dev_vec
        // elast_dev_vec_f is end-of-step
        // double elast_dev_vec_f[ntvec];
        vecsVapb<ntvec>(elas_delta_dev_vec, elas_dt_dev_vec, m_elast_d5_n);
        if constexpr(calc_strain_rate) {
            vecsVsa<ntvec>(elas_dt_dev_vec, m_inv_dt); // _now_ elas_dt_dev_vec has dt contributions
        }
        }

        /*
        * NOTES :
        * () should be equivalent to what happens in get_elas_strain_state<false>
        * () not necessarily safe if elast_dev_press_vec is the same memory as _elast_d5_n or quat is the same as _xtal_ori_quat_n
        */
        __ecmech_hdev__
        inline
        void stateFromX(double* const elast_dev_vec,
                    const double* const x) const
        {
        double elast_dev_vec_delta[ecmech::ntvec] = {};
        this->get_elas_strain_state(elast_dev_vec, elast_dev_vec_delta, x);
        }

        __ecmech_hdev__
        inline
        void get_elas_strain_residual(double* const residual,
                                    const double epsdot_scale_inv,
                                    const double* const elas_dt_dev_vec,
                                    const double* const plastic_def_rate_dev_vec,
                                    const double* const def_rate_dev_vec_lattice) const
        {
        for (size_t iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            residual[m_ind_sub_elas + iTvec] = epsdot_scale_inv * ( // SCALING
                m_inv_a_vol * elas_dt_dev_vec[iTvec] + plastic_def_rate_dev_vec[iTvec] - def_rate_dev_vec_lattice[iTvec]);
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
                                            const double* const elas_dt_dev_vec,
                                            const double elast_elast_factor,
                                            const double* const dWp_hat_delast_strain,
                                            const double* const A_e_M35) const
        {
        // d(B_xi)/d(elast_dev_press_vecs_f)
        //
        RAJA::View<double, RAJA::Layout<2>> jacob_re(jacobian, JAC_SIZE, JAC_SIZE);

        double A_edot_M35[ecmech::nwvec * ecmech::ntvec];
        M35_d_AAoB_dA(A_edot_M35, elas_dt_dev_vec);

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
        double dhdot_delas_strain[num_hard * ecmech::ntvec];
        vecsMABT<num_hard, ecmech::ntvec, num_slip>(dhdot_delas_strain, dhard_dgdot, dgdot_delast_strain);
        RAJA::View<double, RAJA::Layout<2>> jacob_he(jacobian, JAC_SIZE, JAC_SIZE);
        for (size_t iH = 0; iH < num_hard; ++iH) {
            for (size_t jE = 0; jE < ecmech::ntvec; ++jE) {
                // could also make dhdot_delas_strain into a RAJA view, but not really needed
                jacob_he(iH + ind_sub_h, jE) = -m_dt * dhdot_delas_strain[ECMECH_NM_INDX(iH, jE, num_hard, ecmech::ntvec) ];
            }
        }       
        }

        public:
        static constexpr size_t m_ind_sub_elas = 0; // ntvec end_point
        const ThermoElastN& m_thermo_elast_n;
        const double m_dt, m_det_vol, m_det_v_e;
        const double m_pressure_EOS, m_tkelv;
        const double* const m_elast_d5_n;
        const double m_inv_dt, m_inv_det_vol, m_a_vol, m_inv_a_vol;
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
        * () should be equivalent to what happens in get_elas_strain_state<false>
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