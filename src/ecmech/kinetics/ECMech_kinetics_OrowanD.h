// -*-c++-*-

#ifndef ECMECH_KINETICS_OROWAND_H
#define ECMECH_KINETICS_OROWAND_H

#include <cassert>
#include <cmath>

#include <string>
#include <vector>

#include "RAJA/RAJA.hpp"

namespace ecmech {
   /**
    * slip and hardening kinetics
    * based on a single Kocks-Mecking dislocation density
    * balanced thermally activated MTS-like slip kinetics with phonon drag effects
    *
    * see \cite{hmx}
    *
    * if withGAthermal then
    *  see subroutine kinetics_mtspwr_d in mdef : (l_mts, l_mtsp, l_plwr)
    *    ! like kinetics_mtswr_d, but with pl%tau_a (possible associated
    *    ! with the Peierls barrier) being the thermally activated part and
    *    ! g being athermal
    *       ! use balanced and pegged MTS model;
    *       ! add to it a low rate sensitivity power law model to take over for high stresses;
    *       ! and combine with drag limited kinetics
    * else then see subroutine kinetics_mtswr_d in mdef
    *
    *   ! note: gdot_w, gdot_r are always positive by definition
    *   !
    *   ! tkelv should only be used for derivative calculations
    *
    * templated on p and q being 1 or not;
    * might eventually template on number of slip systems, but to not do so just yet
    */
   template<bool withGAthermal,
            bool pOne, // l_p_1
            bool qOne, // l_q_1
            bool isotropic, // H^{\alpha\beta} = 1 so isotropic interaction matrix
            bool perSS, // If varying params per SS usually used for non-cubic materials
            int nVPer, // If perSS then nVPer should equal nslip
            class SlipGeom,
            bool LOGFORM = false> // LOGFORM dictates whether or not we use a logrithmic form for our hardness update
   class KineticsOrowanD
   {
      public:
         static constexpr int nH = 2 * SlipGeom::nslip; // Number of mobile and total dislocation density
         static constexpr int nIH = isotropic ? 1 : (SlipGeom::nslip * SlipGeom::nslip); // Number of params in interaction matrix
         static constexpr int nParams = 12 + 4 * nVPer + nH + nIH + SlipGeom::nParams;
         static constexpr int nVals = 1 + nVPer + 2 * SlipGeom::nslip; //Our ref_slip_rate, CRSS, C1/T, and b*q_m params
         static constexpr int nEvolVals = SlipGeom::nslip; // We really don't need to evolve anything here
         // constructor
         __ecmech_hdev__
         KineticsOrowanD(int _nslip) : nslip(_nslip) {
            assert(nslip == SlipGeom::nslip);
            if (perSS) {
               assert(nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
         }
         // deconstructor
         ~KineticsOrowanD() = default;

         // constructor
         __ecmech_hdev__
         KineticsOrowanD(const double* const params, int _nslip) :
         nslip(_nslip)
         {
            assert(nslip == SlipGeom::nslip);
            if (perSS) {
               assert(nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
            setParams(params);
         }

         __ecmech_host__
         inline void setParams(const std::vector<double> & params)
         {
            setParams(params.data());
         }

         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            //////////////////////////////
            // power-law stuff

            m_mu_ref = *parsIt; ++parsIt;
            m_tkelv_ref = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               m_berg_mag[iVal] = *parsIt; ++parsIt;
            }

            m_lbar_b = *parsIt; ++parsIt;

            m_gam_ro = *parsIt; ++parsIt;
            m_wrD = *parsIt; ++parsIt;

            // thermal activation params
            m_fD = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               m_c_1[iVal] = *parsIt; ++parsIt;
            }

            m_tau_a = *parsIt; ++parsIt;
            m_p = *parsIt; ++parsIt;
            m_q = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               m_c_2[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal < nIH; ++iVal) {
               m_inter_mat[iVal] = *parsIt; ++parsIt;
            }

            if (withGAthermal) {
               assert(m_tau_a > zero);
            }
            if (pOne) {
               assert(m_p == one);
            }
            if (qOne) {
               assert(m_q == one);
            }

            // Figure out what the equivalent from the KMBalD is for the down below
            // plaw_from_elawRef
            //
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               // pl%xm = getMtsxmEffective(pl, mu_ref, T_ref)
               double xm = one / (two * ((m_c_1[iVal] / m_tkelv_ref) * m_mu_ref * m_p * m_q));
               //
               // CALL fill_power_law(pl)
               // xmm  = xm - one ;
               m_xnn[iVal] = one / xm;
               m_xn[iVal] = m_xnn[iVal] - one;
               // xMp1 = xnn + one
               //
               // CALL set_t_min_max(pl)
               // These factors are the same from the balanced-MTS kinetic mobility law
               // and if the ratio tau/crss is below m_t_min we aren't moving at all
               // and if the ratio is above m_t_max we're strictly in phonon drag mobility
               m_t_min[iVal] = pow(ecmech::gam_ratio_min, xm);
               m_t_max[iVal] = pow(ecmech::gam_ratio_ovf, xm);

            }

            //////////////////////////////
            // Dislocation evolution stuff

            m_c_ann = *parsIt; ++parsIt;
            m_d_ann = *parsIt; ++parsIt;
            m_c_trap = *parsIt; ++parsIt;
            m_c_mult = *parsIt; ++parsIt;

            //////////////////////////////
            // Dislocation Densities
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               m_qM[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               m_qT[iVal] = *parsIt; ++parsIt;
            }

            m_hdn_min = m_qM[0];
            // The mobile dd should be the smallest so find the smallest
            // here and base our m_hdn_min on that.
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               if (m_hdn_min > m_qM[iVal]) {
                  m_hdn_min = m_qM[iVal];
               }
            }
            // Might want to make this smaller if provided large initial DD value?
            m_hdn_min *= 1.0e-4;

            //////////////////////////////
            // Initialize slip system matrix
            {
               // Unfortunately, it looks like the simplest way to have this work for various
               // slip systems is by passing in the SlipGeom's params in twice...
               SlipGeom slipgeom(parsIt);
               parsIt += SlipGeom::nParams;
               const double* mref = slipgeom.getM();
               const double* sref = slipgeom.getS();
               // our forest interaction matrix has the following calculation:
               // A^{\alpha\beta} = 1/2 * (|m^alpha \cdot s^alpha| + |m^alpha \cdot (m^beta \cross s^beta)|)
               RAJA::View<const double, RAJA::Layout<2> > mView(mref, SlipGeom::nslip, ecmech::ndim);
               RAJA::View<const double, RAJA::Layout<2> > sView(sref, SlipGeom::nslip, ecmech::ndim);
               RAJA::View<double, RAJA::Layout<2> > aView(&m_a_mat[0], SlipGeom::nslip, SlipGeom::nslip);

               for (int alpha = 0; alpha < SlipGeom::nslip; alpha++) {
                  for (int beta = 0; beta < SlipGeom::nslip; beta++) {
#ifndef ORO_USE_INTERMAT
                     const double mds = mView(alpha, 0) * sView(beta, 0) +
                                        mView(alpha, 1) * sView(beta, 1) +
                                        mView(alpha, 2) * sView(beta, 2);
                     const double mdmxs = mView(alpha, 0) * (mView(beta, 1) * sView(beta, 2) - mView(beta, 2) * sView(beta, 1)) +
                                          mView(alpha, 1) * (mView(beta, 2) * sView(beta, 0) - mView(beta, 0) * sView(beta, 2)) +
                                          mView(alpha, 2) * (mView(beta, 0) * sView(beta, 1) - mView(beta, 1) * sView(beta, 0));
                     aView(alpha, beta) = 1.0 / 2.0 * (std::abs(mds) + std::abs(mdmxs));
#else
                     // use interaction matrix
                     aView(alpha, beta) = m_inter_mat[alpha * SlipGeom::nslip + beta];
#endif
                  }
               }
            }

            //////////////////////////////

#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
         }

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
#ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
#endif
            params.push_back(m_mu_ref);
            params.push_back(m_tkelv_ref);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(m_berg_mag[iVal]);
            }

            params.push_back(m_lbar_b);
            // phonon drag params
            params.push_back(m_gam_ro);
            params.push_back(m_wrD);

            // thermal activation params
            params.push_back(m_fD);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(m_c_1[iVal]);
            }

            params.push_back(m_tau_a);
            params.push_back(m_p);
            params.push_back(m_q);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(m_c_2[iVal]);
            }

            for (int iVal = 0; iVal < nIH; ++iVal) {
               params.push_back(m_inter_mat[iVal]);
            }

            //////////////////////////////
            // Dislocation evolution stuff

            params.push_back(m_c_ann);
            params.push_back(m_d_ann);
            params.push_back(m_c_trap);
            params.push_back(m_c_mult);

            //////////////////////////////
            // Dislocation Densities
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               params.push_back(m_qM[iVal]);
            }

            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               params.push_back(m_qT[iVal]);
            }

            //////////////////////////////
#ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
#endif
         }

         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_dd_mobile_" + std::to_string(iSlip));
               init.push_back(m_qM[iSlip]);
               plot.push_back(true);
               state.push_back(true);
            }

            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_dd_total_" + std::to_string(iSlip));
               init.push_back(m_qT[iSlip]);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         const int nslip; // could template on this if there were call to do so

         //////////////////////////////
         // MTS-like stuff

         // parameters
         double m_lbar_b; // We might need to make this per SS as well
         double m_mu_ref;
         double m_tkelv_ref;
         double m_fD;
         double m_berg_mag[nVPer];
         double m_c_1[nVPer];
         double m_tau_a;
         double m_c_2[nVPer];
         double m_p; // only used if pOne is false
         double m_q; // only used if qOne is false
         double m_inter_mat[nIH]; // symmetric matrix

         double m_gam_ro;
         double m_wrD;

         // derived from parameters
         double m_t_max[nVPer], m_t_min[nVPer], m_xn[nVPer], m_xnn[nVPer];

         //////////////////////////////
         // Dislocation evolution stuff

         double m_c_ann;
         double m_d_ann;
         double m_c_trap;
         double m_c_mult;
         // stored c-style
         double m_a_mat[SlipGeom::nslip * SlipGeom::nslip]; // Forest interaction matrix

         //////////////////////////////
         // Initial dislocation densities
         // so _hdn_init in other models

         double m_qM[SlipGeom::nslip];
         double m_qT[SlipGeom::nslip];
         double m_hdn_min;

      public:

         __ecmech_hdev__
         inline
         double
         getFixedRefRate(const double* const vals) const
         {
            return vals[0];
         }

         /**
          * @brief Akin to hs_to_gss, power_law_tdep_vals, and plaw_from_hs
          *
          * Could eventually bring in additional pressure and temperature dependence through the dependence of _mu on such ;
          * see use of mu_factors in Fortran code
          */
         __ecmech_hdev__
         inline
         double
         getVals(double* const vals, // [nVals]
                 double, // p, not used
                 double tkelv,
                 const double* const h_state
                 ) const
         {
            double const nVPerInv = 1.0 / nslip;

            double maxRefRate = 0.0;
            double hdnScale = 0.;
            for (int iVal = 0; iVal < nslip; ++iVal) {
               const double int_q = isotropic ? sqrt(m_inter_mat[0] * vecsssumabs<SlipGeom::nslip>(&h_state[nslip])) :
                                                sqrt(vecsyadotb<SlipGeom::nslip>(&m_inter_mat[iVal * nslip], &h_state[nslip]));
               const double hdnI = perSS ? (m_c_2[iVal] * int_q) : (m_c_2[0] * int_q);
               hdnScale += hdnI;
               vals[1 + iVal] = hdnI;
               vals[1 + nslip + iVal] = h_state[iVal];
               // Thermal activation + phonon ref slip rate = (1/(f_D * \bar{L}/b * sqrt(qM_0)/sqrt(qM)) + 1/(gammadot_r0 * qM))^-1
               const double isqrth = 1.0 / sqrt(vals[1 + nslip + iVal]);
               const double rate = 1.0 / ((1.0 / (m_lbar_b * m_fD * isqrth)) + (1.0 / (m_gam_ro * vals[1 + nslip + iVal])));
               if (rate > maxRefRate) {
                  maxRefRate = rate;
               }
            }

            // average flow strength across all slip systems
            hdnScale = hdnScale * nVPerInv;
            vals[0] = maxRefRate;

            for (int iVal = 0; iVal < nVPer; ++iVal) {
               vals[1 + 2 * nslip + iVal] = m_c_1[iVal] / tkelv; // _c_t
            }

            return hdnScale;
         }

         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   const double* const tau,
                   const double* const vals
                   ) const
         {
            for (int iSlip = 0; iSlip<this->nslip; ++iSlip) {
               bool l_act;
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip],
                              vals, iSlip,
                              tau[iSlip],
                              m_mu_ref // gss%ctrl%mu(islip)
                              );
            }
         }

         /**
          * like mts_dG, but with output args first
          */
         __ecmech_hdev__
         inline
         void
         get_mts_dG(double &exp_arg,
                    double &mts_dfac,
                    const double c_e, const double denom_i,
                    const double t_frac) const {
            mts_dfac = c_e * denom_i;

            double p_func;
            if (pOne) {
               p_func = t_frac;
            }
            else {
               if (fabs(t_frac) < idp_tiny_sqrt) {
                  // !! p_dfac is either zero or blows up
                  // !IF (pl%p > one) THEN ! no longer allowed
                  // !   mts_dfac = zero
                  // !ELSE
                  // ! blows up, but just set big
                  p_func = zero;
                  mts_dfac = mts_dfac * 1e10;
                  // !END IF
               }
               else {
                  p_func = pow(fabs(t_frac), m_p);
                  p_func = copysign(p_func, t_frac);
                  mts_dfac = mts_dfac *
                             m_p * p_func / t_frac; // always positive
               }
            }

            double q_arg = one - p_func;
            double pq_fac;
            if (q_arg < idp_tiny_sqrt) {
               // peg
               q_arg = zero;
               mts_dfac = zero;
               pq_fac = zero;
            }
            else {
               if (qOne) {
                  pq_fac = q_arg;
               }
               else {
                  double temp = pow(fabs(q_arg), m_q);
                  mts_dfac = mts_dfac *
                             m_q * temp / fabs(q_arg); // always positive
                  pq_fac = copysign(temp, q_arg);
               }
            }

            exp_arg = -c_e * pq_fac;
         }

         /**
          * see subroutine kinetics_mtspwr_d in mdef
          */
         __ecmech_hdev__
         inline
         void
         evalGdot(
            double & gdot,
            bool   & l_act,
            double & dgdot_dtau, // wrt resolved shear stress
            const double* const vals,
            int      iSlip,
            double   tau,
            double   mu
            ) const
         {
            static const double gdot_w_pl_scaling = 10.0;
            static const double one = 1.0, zero = 0.0;

            const double gIn = vals[1 + iSlip];
            const double qm = vals[1 + nslip + iSlip];
            const double xn = perSS ? m_xn[iSlip] : m_xn[0];
            const double xnn = perSS ? m_xnn[iSlip] : m_xnn[0];
            const double t_min = perSS ? m_t_min[iSlip] : m_t_min[0];
            const double t_max = perSS ? m_t_max[iSlip] : m_t_max[0];
            const double c_t = perSS ? vals[1 + 2 * nslip + iSlip] : vals[1 + 2 * nslip];
            const double gam_w = m_lbar_b * m_fD / sqrt(qm);
            const double gam_r = m_gam_ro * qm;

            // zero things so that can more easily just return if inactive
            gdot = zero;
            //
            dgdot_dtau = zero;
            l_act = false;

            double g_i;
            double gAth;
            if (withGAthermal) {
               gAth = gIn;
               g_i = one / m_tau_a;
            }
            else {
               gAth = m_tau_a;
               if (tau == zero) {
                  return;
               }
               g_i = one / gIn;
            }
            double at_0 = fmax(zero, fabs(tau) - gAth) * g_i;
            // calculate drag limited kinetics
            //
            double gdot_r, dgdot_r_dtau;
            {
               double exp_arg = (fabs(tau) - gAth) / m_wrD;
               double temp;
               if (exp_arg < gam_ratio_min) { // ! IF (gdot_r < gam_ratio_min) THEN
                  // note that this should catch tau <= g
                  return;
               }
               else if (exp_arg < idp_eps_sqrt) {
                  // linear expansion is cheaper and more accurate
                  gdot_r = gam_r * exp_arg;
                  temp = one - exp_arg; // still use temp below as approximation to exp(-fabs(tau)/m_wrD)
               }
               else {
                  temp = exp(-exp_arg);
                  gdot_r = gam_r * (one - temp);
               }
               dgdot_r_dtau = gam_r * temp / m_wrD;
            }
            //
            if (at_0 > t_max) {
               // have overflow of thermally activated kinetics, purely drag limited

               gdot = gdot_r;

               dgdot_dtau = dgdot_r_dtau;
               gdot = copysign(gdot, tau);

               l_act = true;
               return;
            }

            double gdot_w, dgdot_w_dtau;
            double dgdot_w_dg; // only used if !withGAthermal
            //
            // calculate thermally activated kinetics
            {
               double c_e = c_t * mu;
               //
               const double t_frac = (fabs(tau) - gAth) * g_i;
               double exp_arg, mts_dfac;
               get_mts_dG(exp_arg, mts_dfac, c_e, g_i, t_frac);
               //
               if (exp_arg < ln_gam_ratio_min) {
                  // effectively zero due to thermally activated kinetics
                  l_act = false;
                  return;
               }
               //
               //
               // !IF (exp_arg > ln_gam_ratio_ovf) THEN
               // !END IF
               // ! do not need to check the above condition because have pegged the MTS part of the kinetics
               //
               gdot_w = gam_w * exp(exp_arg);
               dgdot_w_dtau = mts_dfac * gdot_w;
               if (!withGAthermal) {
                  dgdot_w_dg = dgdot_w_dtau * t_frac;
               }
               //
               double t_frac_m = (-fabs(tau) - gAth) * g_i;
               double exp_arg_m, mts_dfac_m;
               get_mts_dG(exp_arg_m, mts_dfac_m, c_e, g_i, t_frac_m);
               //
               if (exp_arg_m > ln_gam_ratio_min) {
                  // non-vanishing contribution from balancing MTS-like kinetics
                  double gdot_w_m = gam_w * exp(exp_arg_m);
                  gdot_w = gdot_w - gdot_w_m;
                  double contrib = mts_dfac_m * gdot_w_m;
                  dgdot_w_dtau = dgdot_w_dtau - contrib; // sign used to be the other way, but suspect that was a bug
                  if (!withGAthermal) {
                     dgdot_w_dg = dgdot_w_dg - contrib * t_frac_m;
                  }
                  if (fabs(gdot_w / gam_w) < gam_ratio_min) {
                     // effectively zero from roundoff
                     l_act = false;
                     return;
                  }
               }
            }

            if (at_0 > t_min) {
               // need power-law part

               double abslog = log(at_0);
               double blog = xn * abslog;
               double temp = (gam_w * gdot_w_pl_scaling) * exp(blog);

               double gdot_w_pl = temp * at_0; // not signed ! copysign(at_0,tau)
               gdot_w = gdot_w + gdot_w_pl;

               double contrib = temp * xnn * g_i;
               dgdot_w_dtau = dgdot_w_dtau + contrib;
               if (!withGAthermal) {
                  dgdot_w_dg = dgdot_w_dg + contrib * at_0;
               }
            }

            l_act = true;
            //
            {
               // gdot = gdot_r;
               gdot = one / (one / gdot_w + one / gdot_r);
               const double gdrdiv2 = one / (gdot_r * gdot_r);
               const double gdwdiv2 = one / (gdot_w * gdot_w);

               // dgdot_dtau = dgdot_r_dtau;
               dgdot_dtau = (gdot * gdot) * (dgdot_w_dtau * gdwdiv2 + dgdot_r_dtau * gdrdiv2);
               //
            }

            gdot = copysign(gdot, tau);
         } // evalGdot

         __ecmech_hdev__
         inline
         int
         updateH(double* const hs_u,
                 const double* const hs_o,
                 double dt,
                 const double* const gdot,
                 const double* const hvals,
                 double tkelv,
                 int outputLevel = 0) const
         {

            double ihs_o[SlipGeom::nslip * 2];
            double nu[SlipGeom::nslip];
            for (int i = 0; i < nslip * 2; i++) {
               if (i < nslip) {
                  const double div = perSS ? fmax(hs_o[i], m_hdn_min) * m_berg_mag[i] :
                                     fmax(hs_o[i], m_hdn_min) * m_berg_mag[0];
                  nu[i] = abs(gdot[i]) / (div);
               }
               ihs_o[i] = fmax(hs_o[i], m_hdn_min);
               if (LOGFORM) {
                  ihs_o[i] = log(ihs_o[i]);
               }
            }

            int nFEvals = updateHN<KineticsOrowanD>(this,
                                                   &hs_u[0], &ihs_o[0], dt, nu, hvals, tkelv,
                                                   outputLevel);
            if (LOGFORM) {
               for (int i = 0; i < nslip * 2; i++) {
                  hs_u[i] = exp(hs_u[i]);
               }
            }
            else
            {
               // We need to check that none of our solutions became negative
               // If we did obtain something negative then we should abort
               // It means our time step was too large for this step.
               // If this is not desirable / possible then we should probably
               // do a terrible hack and cut the dt by some factor resolve things by
               // assuming a constant slip rate during the time step, and then
               // evolve the dd content. We would get a solution, but it wouldn't necessarily
               // be correct.
               bool flag = false;
               for (int i = 0; i < 2 * nslip; i++) {
                  if(hs_u[i] < zero) {
                     flag = true;
                     break;
                  }
               }
               if (flag)
               {
                  ECMECH_WARN(__func__, "Solver returned negative dislocation values trying again by substepping through the solution");
                  // This is pretty ad-hoc but it seems to work fairly well for a number of simple test cases.
                  // It's definitely not the best way to probably do things though...
                  const double dtnew = dt / 10.0;
                  double hs_temp[2 * SlipGeom::nslip];

                  for (int iSlip = 0; iSlip < 2 * SlipGeom::nslip; iSlip++) {
                     hs_u[iSlip] = fmax(hs_o[iSlip], m_hdn_min);
                  }

                  for (int i = 0; i < 10; i++)
                  {
                     for (int iSlip = 0; iSlip < 2 * SlipGeom::nslip; iSlip++) {
                        hs_temp[iSlip] = fmax(hs_u[iSlip], m_hdn_min);
                        if (iSlip < nslip)
                        {
                           const double div = perSS ? fmax(hs_temp[iSlip], m_hdn_min) * m_berg_mag[iSlip] :
                           fmax(hs_temp[iSlip], m_hdn_min) * m_berg_mag[0];
                           nu[iSlip] = abs(gdot[iSlip]) / (div);
                        }
                     }
                     nFEvals += updateHN<KineticsOrowanD>(this,
                                                         &hs_u[0], hs_temp, dtnew, nu, hvals, tkelv,
                                                         outputLevel);
                     flag = false;
                     for (int iSlip = 0; iSlip < 2 * nslip; iSlip++) {
                        if(hs_u[iSlip] < zero) {
                           flag = true;
                           break;
                        }
                     }
                  }

                  if (flag)
                  {
                     for (int iSlip = 0; iSlip < 2 * SlipGeom::nslip; iSlip++) {
                        printf("dd[%d]: %lf ", iSlip, hs_u[iSlip]);
                     }
                     printf("\n");
                     ECMECH_FAIL(__func__, "Solver returned negative dislocation values!");
                  }
               }
            }

            return nFEvals;
         }

         __ecmech_hdev__
         inline
         void
         getEvolVals(double* const evolVals,
                     const double* const nu
                     ) const
         {
            // We're not really evolving anything here at this point in time
            // so we can just return...
            // The Jacobian doesn't take in gdots / nu so we're left with this.
            for (int i = 0; i < SlipGeom::nslip; i++) {
               evolVals[i] = nu[i];
            }
         }

         __ecmech_hdev__
         inline
         void
         getSdotN( double* sdot,
                   double* dsdot_ds,
                   const double* const h_i,
                   const double* const evolVals,
                   const double* const /*hvals*/,
                   double /*tkelv*/                   ) const
         {
            // Hopefully, the compiler is pretty smart here and is able to optimize these
            // loops as if we're using the templated values. Since, this is essentially
            // a constexpr for these variables.
            const int nslip = SlipGeom::nslip;
            const int JDIM = 2;
            const int nDimSys = 2 * SlipGeom::nslip;

            double forest_dis[nslip];
            constexpr int h_content = (LOGFORM) ? 2 * SlipGeom::nslip : 1;
            double hexp[h_content];
            if (LOGFORM) {
               for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                  hexp[iDD] = exp(h_i[iDD]);
               }
            }
            const double* const h = (LOGFORM) ? &hexp[0] : h_i;
            vecsVMa<SlipGeom::nslip>(&forest_dis[0], &m_a_mat[0], &h[nslip]);

            for (int iM = 0; iM < nslip; iM++) {
               const double sqrt_fD = sqrt(forest_dis[iM]);
               const double q_dmult = m_c_mult * sqrt_fD * h[iM] * evolVals[iM];
               const double q_dtrap = m_c_trap * sqrt_fD * h[iM] * evolVals[iM];
               // This could become a very large number and could become problematic
               // later on. Do we want to cap it at some large value?
               // Although, it might be that this is only a problem if q and qM are defined
               // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
               const double q_dann = m_c_ann * m_d_ann * h[iM] * h[iM] * evolVals[iM];
               // mobile dislocation density rate of change
               sdot[iM] = q_dmult - q_dtrap - q_dann;
               // total dislocation density rate of change
               sdot[iM + nslip] = q_dmult - q_dann;
            }

            if (LOGFORM) {
               for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                  sdot[iDD] *= (1.0 / h[iDD]);
               }
            }
            // The dsdot_ds calculation for our nonlinear solve
            if (dsdot_ds) {
               // zero out dsdot_ds matrix
               for (int i = 0; i < nDimSys * nDimSys; i++) {
                  dsdot_ds[i] = ecmech::zero;
               }

               RAJA::View<double, RAJA::Layout<JDIM> > dsdot_ds_view(dsdot_ds, nDimSys, nDimSys);
               // dqM/dqM portion of dsdot_ds
               for (int iM = 0; iM < nslip; iM++) {
                  const double sqrt_fD = sqrt(forest_dis[iM]);
                  const double q_dmult_dtrap = (m_c_mult - m_c_trap) * sqrt_fD;
                  // Although, it might be that this is only a problem if q and qM are defined
                  // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
                  const double q_dann = 2 * m_c_ann * m_d_ann * h[iM];
                  dsdot_ds_view(iM, iM) = evolVals[iM] * (q_dmult_dtrap - q_dann);
               }

               // dq/dqM portion of dsdot_ds
               for (int iT = 0; iT < nslip; iT++) {
                  const double sqrt_fD = sqrt(forest_dis[iT]);
                  const double q_dmult = m_c_mult * sqrt_fD;
                  // This could become a very large number and could become problematic
                  // later on. Do we want to cap it at some large value?
                  // Although, it might be that this is only a problem if q and qM are defined
                  // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
                  const double q_dann = 2 * m_c_ann * m_d_ann * h[iT];
                  dsdot_ds_view(iT + nslip, iT) =  evolVals[iT] * (q_dmult - q_dann);
               }

               RAJA::View<const double, RAJA::Layout<JDIM> > amat(&m_a_mat[0], nslip, nslip);
               // dq/dq portion of dsdot_ds
               for (int iT = 0; iT < nslip; iT++) {
                  for (int jT = 0; jT < nslip; jT++) {
                     // First, terms found only on the diagonal of this submatrix
                     const double ifact = ecmech::onehalf / sqrt(forest_dis[iT]);
                     const double q_dmult = m_c_mult * amat(iT, jT) * ifact;

                     dsdot_ds_view(iT + nslip, jT + nslip) = h[iT] * evolVals[iT] * q_dmult;
                  }
               }

               // dqM/dq portion of dsdot_dt
               for (int iT = 0; iT < nslip; iT++) {
                  for (int jT = 0; jT < nslip; jT++) {
                     const double ifact = ecmech::onehalf / sqrt(forest_dis[iT]);
                     const double q_dmult_dtrap = (m_c_mult - m_c_trap) * amat(iT, jT) * ifact;

                     dsdot_ds_view(iT, jT + nslip) = h[iT] * evolVals[iT] * (q_dmult_dtrap);
                  }
               }

               if (LOGFORM) {
                  for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                     dsdot_ds_view(iDD, iDD) -= sdot[iDD];
                  }
                  for (int iDD = 0; iDD < 2 * nslip; iDD++) {
                     for (int jDD = 0; jDD < 2 * nslip; jDD++) {
                        dsdot_ds_view(iDD, jDD) *= h[jDD] / h[iDD];
                     }
                  }
               }
            } // if dsdot_ds
         }
   }; // class KineticsOrowanD
} // namespace ecmech

#endif // ECMECH_KINETICS_OROWAND_H
