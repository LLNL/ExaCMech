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
        static const int iHistLbGdot = iHistLbH + Kinetics::nH;
        static const int numHist = iHistLbH + Kinetics::nH + SlipGeom::nslip;
    }; // NumHist

    template<class ThermoElastN>
    class EvptnLatticeStrainProblem
    {
        public:
        static const int nDimSys = ecmech::ntvec;

        __ecmech_hdev__
        EvptnLatticeStrainProblem(const ThermoElastN& thermoElastN,
                                const double dt,
                                const double detV, 
                                const double eVref, 
                                const double p_EOS, 
                                const double tK,
                                const double* const e_vecd_n)
        : m_thermo_elast_n(thermoElastN),
        m_dt(dt), m_det_vol(detV), m_elast_vol_ref(eVref),
        m_pressure_eos(p_EOS), m_temp_k(tK),
        m_elast_dev_vec_n(e_vecd_n),
        m_inv_dt(1.0 / dt),
        m_inv_det_vol(1.0 / m_det_vol),
        m_a_vol(pow(m_det_vol, onethird)),
        m_inv_a_vol(1.0 / m_a_vol)
        {};

        __ecmech_hdev__
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
        // CALL elawn_T(s_meas, e_vecd_f, crys%elas, tK, .TRUE., a_V, &
        // & p_EOS, eVref, crys%i_eos_model, crys%eos_const &
        // &)
        double elas_dev_vol_vec[ecmech::nsvec];
        vecsVxa<ntvec>(elas_dev_vol_vec, m_inv_a_vol, elas_dev_vec);
        //// tr_Ee = three * DLOG(a_V%r)
        //// CALL trace_to_vecds_s(s_meas%Ee_vecds(SVEC), tr_Ee)
        elas_dev_vol_vec[iSvecS] = sqr3 * log(m_a_vol); // could go into constructor
        //
        //// Kirchhoff stress from elas_dev_vol_vec
        // CALL elawn_lin_op(s_meas%T_vecds, s_meas%Ee_vecds, cem, tK, &
        // & p_EOS, eVref, i_eos_model, eos_const)
        m_thermo_elast_n.eval(kirchoff_stress, elas_dev_vol_vec, m_temp_k, m_pressure_eos, m_elast_vol_ref);
        }

        // used to be elastNEtoC
        __ecmech_hdev__
        inline
        void elas_strain_to_cauchy_stress(double* const cauchy, // nsvec
                                        const double* const e_vecd_f // ntvec
                                        ) const
        {
        double kirchoff[ecmech::nsvec];
        this->elas_strain_to_kirchoff_stress(kirchoff, e_vecd_f);
        m_thermo_elast_n.getCauchy(cauchy, kirchoff, m_inv_det_vol);
        }

        __ecmech_hdev__
        template<bool calc_strain_rate = false>
        inline
        void get_elas_strain_state(double* const elas_delta_dev_vec,
                                double* const elas_dt_dev_vec,
                                const double* const x) const
        {
        //////////////////////////////
        // PULL VALUES out of x, with scalings
        //
        // double edot_vecd[ecmech::ntvec];
        vecsVxa<ntvec>(elas_dt_dev_vec, ecmech::e_scale, x); // elas_dt_dev_vec is now the delta, _not_ yet edot_vecd
        // e_vecd_f is end-of-step
        // double e_vecd_f[ntvec];
        vecsVapb<ntvec>(elas_delta_dev_vec, elas_dt_dev_vec, m_elast_dev_vec_n);
        if constexpr(calc_strain_rate) {
            vecsVsa<ntvec>(elas_dt_dev_vec, m_inv_dt); // _now_ elas_dt_dev_vec has dt contributions
        }
        }

        /*
        * NOTES :
        * () should be equivalent to what happens in get_elas_strain_state<false>
        * () not necessarily safe if e_vecd is the same memory as _e_vecd_n or quat is the same as _Cn_quat
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

        __ecmech_hdev__
        template<size_t JAC_SIZE>
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

        __ecmech_hdev__
        template<size_t JAC_SIZE, size_t ind_sub_r>
        inline
        void get_deriv_omega_wrt_elast_strain(double* const jacobian,
                                            const double* const elas_dt_dev_vec,
                                            const double elast_elast_factor,
                                            const double* const dWp_hat_delast_strain,
                                            const double* const A_e_M35) const
        {
        // d(B_xi)/d(e_vecds_f)
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

        __ecmech_hdev__
        template<size_t JAC_SIZE, size_t ind_sub_h, size_t num_hard, size_t num_slip>
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
        const double m_dt, m_det_vol, m_elast_vol_ref;
        const double m_pressure_eos, m_temp_k;
        const double* const m_elast_dev_vec_n;
        const double m_inv_dt, m_inv_det_vol, m_a_vol, m_inv_a_vol;
    };

    template <size_t ind_sub_r=ecmech::ntvec>
    class EvptnLatticeRotationProblem {
        public:
        static constexpr size_t nDimSys = ecmech::nwvec;

        public:
        __ecmech_hdev__
        EvptnLatticeRotationProblem(const double dt,
                                const double* const Cn_quat)
        : m_dt(dt), m_xtal_ori_quat_n(Cn_quat) {};
        
        __ecmech_hdev__
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
        * () not necessarily safe if e_vecd is the same memory as _e_vecd_n or quat is the same as _Cn_quat
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

        __ecmech_hdev__
        template<size_t JAC_SIZE>
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
        __ecmech_hdev__
        template<size_t JAC_SIZE>
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

#if defined(ECMECH_EXTRA_SOLVERS)
    template<class SlipGeom, class ThermoElastN>
    class RstarUpdProblem
    {
        public:
        static const int nDimSys = ecmech::nwvec;
        
        __ecmech_hdev__
        RstarUpdProblem(const SlipGeom& slipGeom,
                        const ThermoElastN& thermoElastN,
                        double dt,
                        double detV, double eVref, double p_EOS, double tK,
                        const double* const gdot,
                        const double* const e_vecd_n,
                        const double* const Cn_quat,
                        const double* const d_vecd_sm, // okay to pass d_vecds_sm, but d_vecd_sm[iSvecS] is not used
                        const double* const w_veccp_sm)
                        : 
            _slipGeom(slipGeom),
            _thermoElastN(thermoElastN),
            _lattice_strain_prob(thermoElastN, dt, detV, eVref, p_EOS, tK, e_vecd_n),
            _lattice_rot_prob(dt, Cn_quat),
            _eVref(eVref),
            _p_EOS(p_EOS),
            _tK(tK),
            _gdot(gdot),
            _e_vecd_n(e_vecd_n),
            _Cn_quat(Cn_quat),
            _d_vecd_sm(d_vecd_sm), // vel_grad_sm%d_vecds
            _w_veccp_sm(w_veccp_sm) // vel_grad_sm%w_veccp
        {
            _dt_ri = 1.0 / _dt;
            _detV_ri = 1.0 / _detV;
            _a_V = pow(detV, onethird);
            _a_V_ri = 1.0 / _a_V;

            double adots_ref = vecNorm<SlipGeom::nslip>(gdot);

            double eff = vecNorm<ecmech::ntvec>(_d_vecd_sm); // do not worry about factor of sqrt(twothird)
            if (eff < epsdot_scl_nzeff * adots_ref) {
                _epsdot_scale_inv = one / adots_ref;
            }
            else {
                _epsdot_scale_inv = fmin(one / eff, 1e6 * _lattice_strain_prob.m_dt);
            }
            //
            _rotincr_scale_inv = _lattice_strain_prob.m_inv_dt * _epsdot_scale_inv;
        }

        // deconstructor
        __ecmech_hdev__
        ~RstarUpdProblem() {}

                    __ecmech_hdev__
        inline
        void stateFromX(double* const quat,
                        const double* const x) {
            _lattice_rot_prob.stateFromX(quat, x);
        }
        
        __ecmech_hdev__
        bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            bool doComputeJ = (Jacobian != nullptr);

            if (doComputeJ) {
                // zero the Jacobian so that do not need to worry about zero
                // entries in the midst of other things later
                //
                for (int ijJ = 0; ijJ< _nXnDim; ++ijJ) {
                    Jacobian[ijJ] = 0.0;
                }
            }
            //
            for (int iR = 0; iR<nDimSys; ++iR) {
                resid[iR] = 0.0;
            }

            double edot_vecd[ecmech::ntvec];
            // Calculate what this edot_vecd term should be given the current
            // state information.
            for (int i = 0; i < ecmech::ntvec; i++)
            {
                edot_vecd[i] = _a_V_ri * (d_vecd_lat[i] - pl_vecd[i]);
            }

            double xi_f[nwvec];
            vecsVxa<nwvec>(xi_f, ecmech::r_scale, x);
            //
            // not done in EvpC :
            // CALL exp_map_cpvec(A, xi_f)
            // CALL get_c(c, A, C_n)
            //
            double A_quat[ecmech::qdim];
            emap_to_quat(A_quat, xi_f);
            //
            double C_quat[ecmech::qdim];
            get_c_quat(C_quat, A_quat, _lattice_rot_prob.m_xtal_ori_quat_n);
            //
            double C_matx[ecmech::ndim * ecmech::ndim];
            quat_to_tensor(C_matx, C_quat);
            //
            double qr5x5_ls[ecmech::ntvec * ecmech::ntvec];
            get_rot_mat_vecd(qr5x5_ls, C_matx);

            // CALL matt_x_vec_5(qr5x5_ls, vel_grad_sm%d_vecds(1:TVEC), d_vecd_lat)
            double d_vecd_lat[ecmech::ntvec];
            // vecsVMTa<ntvec>(d_vecd_lat, qr5x5_ls, _d_vecd_sm);
            // d_vecds_lat(SVEC) = vel_grad_sm%d_vecds(SVEC)
            //
            //// CALL rot_mat_vecd(A, qr5x5_A)
            //
            // CALL rot_mat_wveccp(C_matx, qr3x3_ls) // amounts to qr3x3_ls = C_matx
            // CALL matt_x_vec_3(qr3x3_ls, vel_grad_sm%w_veccp, w_vec_lat)
            double w_vec_lat[ecmech::nwvec]; // assumes nwvec = ndim
            // vecsVMTa<ndim>(w_vec_lat, C_matx, _w_veccp_sm);

            get_xtal_frame_vel_grad_terms(d_vecd_lat, w_vec_lat, _d_vecd_sm,  _w_veccp_sm, C_matx, qr5x5_ls);

            double pl_vecd[ecmech::ntvec] = { 0.0 };
            double pl_wvec[ecmech::nwvec] = { 0.0 }; // \pcDhat
            {
                // Yeah we're going to do a bit more work here but at the gdot and elasticity stuff will always be consistent
                double T_vecds[ecmech::nsvec];
                _lattice_strain_prob.elas_strain_to_kirchoff_stress(T_vecds, _lattice_strain_prob.m_elast_dev_vec_n);
                double dgdot_dtau[SlipGeom::nslip] = { 0.0 }; // crys%tmp2_slp
                get_slip_rate_terms(dgdot_dtau, pl_vecd, pl_wvec, T_vecds, _kin_vals, _slipGeom, _kinetics);
            }

            // Higher-order terms related to the elasticity stuff that's used in the residuals and jacobian calculation
            double A_e_M35[ecmech::nwvec * ecmech::ntvec];
            double ee_wvec[ecmech::nwvec];
            double ee_fac;
            elasticity_higher_order_terms(A_e_M35, ee_wvec, ee_fac, _lattice_strain_prob.m_inv_a_vol, _lattice_strain_prob.m_elast_dev_vec_n, edot_vecd);

            // Residual Calculations
            _lattice_rot_prob.get_omega_residual(resid, _rotincr_scale_inv, ee_fac, xi_f, w_vec_lat, pl_wvec, ee_wvec);

            //////////////////////////////////////////////////////////////////////
            // JACOBIAN, fixed hardness and temperature
            //
            if (doComputeJ) {

                //
                //
                // derivatives with respect to lattice orientation changes
                double dC_quat_dxi_T[ ecmech::nwvec * ecmech::qdim ];
                double dDsm_dxi[ ecmech::ntvec * ecmech::nwvec ];
                double dWsm_dxi[ ecmech::nwvec * ecmech::nwvec ];
                eval_d_dxi_impl_quat(dC_quat_dxi_T, dDsm_dxi, dWsm_dxi,
                                    _d_vecd_sm, _w_veccp_sm,
                                    xi_f, 
                                    _lattice_rot_prob.m_xtal_ori_quat_n,
                                    C_matx, C_quat);

                // d(B_xi)/d(xi_f)
                //
                _lattice_rot_prob.template get_deriv_omega_wrt_omega<nDimSys>(Jacobian, dWsm_dxi);

                const double scaleFactorJ = _rotincr_scale_inv * ecmech::r_scale;
                for (int iJ = 0; iJ<nDimSys; ++iJ) {
                    // Jacobian(i_sub_r:i_sup_r,i_sub_r:i_sup_r) = jacob_rr * rotincr_scale_inv * r_scale
                    for (int jJ = 0; jJ<nDimSys; ++jJ) { // <_i_sup_r
                        int ijJ = ECMECH_NN_INDX(iJ, jJ, nDimSys);
                        Jacobian[ ijJ ] *= scaleFactorJ;
                    }
                }
            }
            return true;
        }

        private:

        const SlipGeom &_slipGeom;
        const ThermoElastN &_thermoElastN;
        const EvptnLatticeStrainProblem<ThermoElastN> _lattice_strain_prob;
        const EvptnLatticeRotationProblem<0> _lattice_rot_prob;

        double _eVref, _p_EOS, _tK, _a_V;
        double _epsdot_scale_inv, _rotincr_scale_inv;

        const double* const _d_vecd_sm; // d_vecds_sm would be fine too -- but do not use _d_vecd_sm[iSvecS];
        const double* const _w_veccp_sm;

        static const int _nXnDim = nDimSys * nDimSys;
    };
#endif

}
}