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
    *   ! tK should only be used for derivative calculations
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
            class SlipGeom>
   class KineticsOrowanD
   {
      public:
         static const int nH = 2 * SlipGeom::nslip; // Number of mobile and total dislocation density
         static const int nIH = isotropic ? 1 : (SlipGeom::nslip * SlipGeom::nslip); // Number of params in interaction matrix
         static const int nParams = 13 + 4 * nVPer + nH + nIH + SlipGeom::nParams;
         static const int nVals = 1 + nVPer + 2 * SlipGeom::nslip; //Our ref_slip_rate, CRSS, C1/T, and b*q_m params
         static const int nEvolVals = 2 * SlipGeom::nslip; // We really don't need to evolve anything here

         // constructor
         __ecmech_hdev__
         KineticsOrowanD(int nslip) : _nslip(nslip) {
            assert(_nslip == SlipGeom::nslip);
            if (perSS) {
               assert(_nslip == nVPer);
            }
            else {
               assert(nVPer == 1);
            }
         };
         // deconstructor
         __ecmech_hdev__
         ~KineticsOrowanD() {}

         __ecmech_host__
         void setParams(const std::vector<double> & params // const double* const params
                        ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            //////////////////////////////
            // power-law stuff

            _mu_ref = *parsIt; ++parsIt;
            _tK_ref = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               _berg_mag[iVal] = *parsIt; ++parsIt;
            }

            _lbar = *parsIt; ++parsIt;
            // phonon drag params
            _shear_speed = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               _c_3[iVal] = *parsIt; ++parsIt;
            }

            // thermal activation params
            _fD = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               _c_1[iVal] = *parsIt; ++parsIt;
            }

            _tau_a = *parsIt; ++parsIt;
            _p = *parsIt; ++parsIt;
            _q = *parsIt; ++parsIt;
            _tau_0 = *parsIt; ++parsIt;
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               _c_2[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal < nIH; ++iVal) {
               _inter_mat[iVal] = *parsIt; ++parsIt;
            }

            if (withGAthermal) {
               assert(_tau_a > zero);
            }
            if (pOne) {
               assert(_p == one);
            }
            if (qOne) {
               assert(_q == one);
            }



            // Figure out what the equivalent from the KMBalD is for the down below
            // plaw_from_elawRef
            //
            for (int iVal = 0; iVal<nVPer; ++iVal) {
               // pl%xm = getMtsxmEffective(pl, mu_ref, T_ref)
               double xm = one / (two * ((_c_1[iVal] / _tK_ref) * _mu_ref * _p * _q));
               //
               // CALL fill_power_law(pl)
               // xmm  = xm - one ;
               _xnn[iVal] = one / xm;
               _xn[iVal] = _xnn[iVal] - one;
               // xMp1 = xnn + one
               //
               // CALL set_t_min_max(pl)
               // These factors are the same from the balanced-MTS kinetic mobility law
               // and if the ratio tau/crss is below _t_min we aren't moving at all
               // and if the ratio is above _t_max we're strictly in phonon drag mobility
               _t_min[iVal] = pow(ecmech::gam_ratio_min, xm);
               _t_max[iVal] = pow(ecmech::gam_ratio_ovf, xm);

               // Factors related to bounding the phonon drag term
               // If we're above _l_max for tau our phonon drag term, \nu_r, approaches shear wave speed
               // If we're below _l_min for tau our phonon drag term, \nu_r, approaches 0.
               _l_min[iVal] = _c_3[iVal] * ecmech::lorentz_min;
               _l_max[iVal] = _c_3[iVal] * ecmech::lorentz_max;
            }

            //////////////////////////////
            // Dislocation evolution stuff

            _c_ann = *parsIt; ++parsIt;
            _d_ann = *parsIt; ++parsIt;
            _c_trap = *parsIt; ++parsIt;
            _c_mult = *parsIt; ++parsIt;

            //////////////////////////////
            // Dislocation Densities
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               _qM[iVal] = *parsIt; ++parsIt;
            }

            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               _qT[iVal] = *parsIt; ++parsIt;
            }

            _hdn_min = _qM[0];
            // The mobile dd should be the smallest so find the smallest
            // here and base our _hdn_min on that.
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               if (_hdn_min > _qM[iVal]) {
                  _hdn_min = _qM[iVal];
               }
            }

            _hdn_min *= 1.0e-4;

            //////////////////////////////
            // Initialize slip system matrix
            {
               // Unfortunately, it looks like the simplest way to have this work for various
               // slip systems is by passing in the SlipGeom's params in twice...
               SlipGeom slipgeom;
               const std::vector<double> paramsThese(parsIt, parsIt + SlipGeom::nParams);
               slipgeom.setParams(paramsThese); parsIt += SlipGeom::nParams;
               const double* mref = slipgeom.getM();
               const double* sref = slipgeom.getS();
               // our forest interaction matrix has the following calculation:
               // A^{\alpha\beta} = 1/2 * (m^alpha \cdot s^alpha + m^alpha \cdot (m^beta \cross s^beta))
               RAJA::View<const double, RAJA::Layout<2> > mView(mref, SlipGeom::nslip, ecmech::ndim);
               RAJA::View<const double, RAJA::Layout<2> > sView(sref, SlipGeom::nslip, ecmech::ndim);
               RAJA::View<double, RAJA::Layout<2> > aView(&_a_mat[0], SlipGeom::nslip, SlipGeom::nslip);

               for (int alpha = 0; alpha < SlipGeom::nslip; alpha++) {
                  for (int beta = 0; beta < SlipGeom::nslip; beta++) {
                     const double mds = mView(alpha, 0) * sView(beta, 0) +
                                        mView(alpha, 1) * sView(beta, 1) +
                                        mView(alpha, 2) * sView(beta, 2);
                     const double mdmxs = mView(alpha, 0) * (mView(beta, 1) * sView(beta, 2) - mView(beta, 2) * sView(beta, 1)) +
                                          mView(alpha, 1) * (mView(beta, 2) * sView(beta, 0) - mView(beta, 0) * sView(beta, 2)) +
                                          mView(alpha, 2) * (mView(beta, 0) * sView(beta, 1) - mView(beta, 1) * sView(beta, 0));
                     aView(alpha, beta) = 1.0 / 2.0 * (mds + mdmxs);
                  }
               }
            }

            //////////////////////////////

            int iParam = parsIt - params.begin();
            assert(iParam == nParams);
         };

         __ecmech_host__
         void getParams(std::vector<double> & params
                        ) const {
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();

            params.push_back(_mu_ref);
            params.push_back(_tK_ref);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(_berg_mag[iVal]);
            }

            params.push_back(_lbar);
            // phonon drag params
            params.push_back(_shear_speed);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(_c_3[iVal]);
            }

            // thermal activation params
            params.push_back(_fD);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(_c_1[iVal]);
            }

            params.push_back(_tau_a);
            params.push_back(_p);
            params.push_back(_q);
            params.push_back(_tau_0);
            for (int iVal = 0; iVal < nVPer; ++iVal) {
               params.push_back(_c_2[iVal]);
            }

            for (int iVal = 0; iVal < nIH; ++iVal) {
               params.push_back(_inter_mat[iVal]);
            }

            //////////////////////////////
            // Dislocation evolution stuff

            params.push_back(_c_ann);
            params.push_back(_d_ann);
            params.push_back(_c_trap);
            params.push_back(_c_mult);

            //////////////////////////////
            // Dislocation Densities
            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               params.push_back(_qM[iVal]);
            }

            for (int iVal = 0; iVal < SlipGeom::nslip; iVal++) {
               params.push_back(_qT[iVal]);
            }

            //////////////////////////////

            int iParam = params.size() - paramsStart;
            assert(iParam == nParams);
         };

         __ecmech_host__
         void getHistInfo(std::vector<std::string> & names,
                          std::vector<double>       & init,
                          std::vector<bool>        & plot,
                          std::vector<bool>        & state) const {
            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_dd_mobile_" + std::to_string(iSlip));
               init.push_back(_qM[iSlip]);
               plot.push_back(true);
               state.push_back(true);
            }

            for (int iSlip = 0; iSlip < SlipGeom::nslip; iSlip++) {
               names.push_back("rho_dd_total_" + std::to_string(iSlip));
               init.push_back(_qT[iSlip]);
               plot.push_back(true);
               state.push_back(true);
            }
         }

      private:

         const int _nslip; // could template on this if there were call to do so

         //////////////////////////////
         // MTS-like stuff

         // parameters
         double _lbar; // We might need to make this per SS as well
         double _mu_ref;
         double _tK_ref;
         double _fD;
         double _shear_speed;
         double _c_3[nVPer];
         double _berg_mag[nVPer];
         double _c_1[nVPer];
         double _tau_a;
         double _tau_0;
         double _c_2[nVPer];
         double _p; // only used if pOne is false
         double _q; // only used if qOne is false
         double _inter_mat[nIH]; // symmetric matrix

         // derived from parameters
         double _t_max[nVPer], _t_min[nVPer], _xn[nVPer], _xnn[nVPer];
         double _l_max[nVPer], _l_min[nVPer];

         //////////////////////////////
         // Dislocation evolution stuff

         double _c_ann;
         double _d_ann;
         double _c_trap;
         double _c_mult;
         // stored c-style
         double _a_mat[SlipGeom::nslip * SlipGeom::nslip]; // Forest interaction matrix

         //////////////////////////////
         // Initial dislocation densities
         // so _hdn_init in other models

         double _qM[SlipGeom::nslip];
         double _qT[SlipGeom::nslip];
         double _hdn_min;

      public:

         __ecmech_hdev__
         inline
         double
         getFixedRefRate(const double* const vals) const
         {
            // Thermal activation term
            // As mobile dislocation densities increase this number can get quite high,
            // so we scale it down by a few order of magnitudes to make it more reasonable.
            return 1.0e-3 * vals[0];
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
                 double tK,
                 const double* const h_state
                 ) const
         {
            double const nVPerInv = 1.0 / _nslip;
            const double fD = _fD;

            double maxThermalRefRate = 0.0;
            double hdnScale = 0.;
            for (int iVal = 0; iVal < _nslip; ++iVal) {
               const double int_q = isotropic ? sqrt(_inter_mat[0] * vecsssumabs<SlipGeom::nslip>(&h_state[_nslip])) :
                                                sqrt(vecsyadotb<SlipGeom::nslip>(&_inter_mat[iVal * _nslip], &h_state[_nslip]));
               const double hdnI = perSS ? (_tau_0 + _c_2[iVal] * int_q) : (_tau_0 + _c_2[0] * int_q);
               hdnScale += hdnI;
               vals[1 + iVal] = hdnI;
               vals[1 + _nslip + iVal] = perSS ? (_berg_mag[iVal] * h_state[iVal]) : (_berg_mag[0] * h_state[iVal]);
               // Thermal activation ref slip rate = b * q_M * f_D * \bar{L}
               const double rate = perSS ? (_berg_mag[iVal] * h_state[iVal] * _lbar * fD)
                                   : (_berg_mag[0] * h_state[iVal] * _lbar * fD);
               if (rate > maxThermalRefRate) {
                  maxThermalRefRate = rate;
               }
            }

            // average flow strength across all slip systems
            hdnScale = hdnScale * nVPerInv;
            vals[0] = maxThermalRefRate;

            for (int iVal = 0; iVal < nVPer; ++iVal) {
               vals[1 + 2 * _nslip + iVal] = _c_1[iVal] / tK; // _c_t
            }

            return hdnScale;
         }

         __ecmech_hdev__
         inline
         void
         evalGdots(double* const gdot,
                   double* const dgdot_dtau,
                   double* const dgdot_dg,
                   const double* const tau,
                   const double* const vals
                   ) const
         {
            for (int iSlip = 0; iSlip<this->_nslip; ++iSlip) {
               bool l_act;
               this->evalGdot(gdot[iSlip], l_act, dgdot_dtau[iSlip], dgdot_dg[iSlip],
                              vals, iSlip,
                              tau[iSlip],
                              _mu_ref // gss%ctrl%mu(islip)
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
                  p_func = pow(fabs(t_frac), _p);
                  p_func = copysign(p_func, t_frac);
                  mts_dfac = mts_dfac *
                             _p * p_func / t_frac; // always positive
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
                  double temp = pow(fabs(q_arg), _q);
                  mts_dfac = mts_dfac *
                             _q * temp / fabs(q_arg); // always positive
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
            double & dgdot_dg, // wrt slip system strength
#if MORE_DERIVS
            double & dgdot_dmu, // wrt shear modulus, not through g
            double & dgdot_dgamo, // wrt reference rate for thermal part
            double & dgdot_dgamr, // wrt reference rate for drag limited part
            double & dgdot_dtK, // wrt temperature, with other arguments fixed
#endif
            const double* const vals,
            int      iSlip,
            double   tau,
            double   mu
#if MORE_DERIVS
            ,
            double   tK
#endif
            ) const
         {
            static const double gdot_w_pl_scaling = 10.0;
            static const double one = 1.0, zero = 0.0;

            const double gIn = vals[1 + iSlip];
            const double bqm = vals[1 + _nslip + iSlip];
            const double xn = perSS ? _xn[iSlip] : _xn[0];
            const double xnn = perSS ? _xnn[iSlip] : _xnn[0];
            const double t_min = perSS ? _t_min[iSlip] : _t_min[0];
            const double t_max = perSS ? _t_max[iSlip] : _t_max[0];
            const double l_min = perSS ? _l_min[iSlip] : _l_min[0];
            const double l_max = perSS ? _l_max[iSlip] : _l_max[0];
            const double c_3 = perSS ? _c_3[iSlip] : _c_3[0];
            const double c_t = perSS ? vals[1 + 2 * _nslip + iSlip] : vals[1 + 2 * _nslip];

            // zero things so that can more easily just return if inactive
            gdot = zero;
            //
            dgdot_dtau = zero;
            dgdot_dg = zero;
#if MORE_DERIVS
            dgdot_dmu = zero;
            dgdot_dgamo = zero;
            dgdot_dgamr = zero;
            dgdot_dtK = zero;
#endif
            l_act = false;

            double g_i;
            double gAth;
            if (withGAthermal) {
               gAth = gIn;
               g_i = one / _tau_a;
            }
            else {
               gAth = _tau_a;
               if (tau == zero) {
                  return;
               }
               g_i = one / gIn;
            }
            double at_0 = fmax(zero, fabs(tau) - gAth) * g_i;

            // calculate drag limited kinetics
            //
            double nu_r, dnu_r_dtau;
#if MORE_DERIVS
            double dnu_r_dtK = zero;
#endif

            {
               const double arg = (fabs(tau) - gAth);
               const double iarg = 1.0 / arg;
               double temp;
               if (arg < l_min) { // ! IF (gdot_r < gam_ratio_min) THEN
                  // note that this should catch tau <= g
                  return;
               }
               else if (arg > l_max) {
                  nu_r = _shear_speed;
                  dnu_r_dtau = zero;
               }
               else {
                  const double factor = c_3 * onehalf;
                  const double sqrterm = sqrt((factor * factor * iarg * iarg) + 1);
                  nu_r = _shear_speed * (sqrt((factor * factor * iarg * iarg) + 1)
                                         - (factor * iarg));
                  // dnu_r / dtau
                  // = -c_s * ( c_3^2 / 4 * 1 / tau^3 / sqrt(c_3^2/4 * 1/tau^2 + 1) - c_3 / 2 * 1 / tau^2)
                  dnu_r_dtau = -_shear_speed * (factor * factor * iarg * iarg * iarg / sqrterm - (factor * iarg * iarg));
               }
#if MORE_DERIVS
#endif
            }
            //
            if (at_0 > t_max) {
               // have overflow of thermally activated kinetics, purely drag limited

               gdot = bqm * nu_r;

               dgdot_dtau = bqm * dnu_r_dtau;
               if (withGAthermal) {
                  dgdot_dg = -copysign(dgdot_dtau, tau);
               }
               else {
                  dgdot_dg = zero;
               }
#if MORE_DERIVS
               dgdot_dmu = zero;
               dgdot_dgamo = zero;
               dgdot_dtK = copysign(dgdotr_dtK, tau);
               dgdot_dgamr = zero;
#endif
               gdot = copysign(gdot, tau);

               l_act = true;
               return;
            }

            double nu_w, dnu_w_dtau;
            double dnu_w_dg; // only used if !withGAthermal
#if MORE_DERIVS
            double dnu_w_dmu;
            double dnu_w_dtK;
#endif
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
#if MORE_DERIVS
               dnu_w_dmu = zero;
               dnu_w_dtK = zero;
#endif
               //
               // !IF (exp_arg > ln_gam_ratio_ovf) THEN
               // !END IF
               // ! do not need to check the above condition because have pegged the MTS part of the kinetics
               //
               nu_w = _lbar * _fD * exp(exp_arg);
               dnu_w_dtau = mts_dfac * nu_w;
               if (!withGAthermal) {
                  dnu_w_dg = dnu_w_dtau * t_frac;
               }
#if MORE_DERIVS
               dnu_w_dmu = nu_w * (-exp_arg / mu);
               dnu_w_dtK = nu_w * (exp_arg / tK); // negatives cancel
#endif
               //
               double t_frac_m = (-fabs(tau) - gAth) * g_i;
               double exp_arg_m, mts_dfac_m;
               get_mts_dG(exp_arg_m, mts_dfac_m, c_e, g_i, t_frac_m);
               //
               if (exp_arg_m > ln_gam_ratio_min) {
                  // non-vanishing contribution from balancing MTS-like kinetics
                  double nu_w_m = _lbar * _fD * exp(exp_arg_m);
                  nu_w = nu_w - nu_w_m;
                  double contrib = mts_dfac_m * nu_w_m;
                  dnu_w_dtau = dnu_w_dtau - contrib; // sign used to be the other way, but suspect that was a bug
                  if (!withGAthermal) {
                     dnu_w_dg = dnu_w_dg - contrib * t_frac_m;
                  }
                  if (fabs(nu_w / (_lbar * _fD)) < gam_ratio_min) {
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
               double temp = (_lbar * _fD * gdot_w_pl_scaling) * exp(blog);

               double nu_w_pl = temp * at_0; // not signed ! copysign(at_0,tau)
               nu_w = nu_w + nu_w_pl;

               double contrib = temp * xnn * g_i;
               dnu_w_dtau = dnu_w_dtau + contrib;
               if (!withGAthermal) {
                  dnu_w_dg = dnu_w_dg + contrib * at_0;
               }
            }

            l_act = true;
            //
            {
               const double nu = one / (one / nu_w + one / nu_r);
               const double ndrdiv2 = one / (nu_r * nu_r);
               const double ndwdiv2 = one / (nu_w * nu_w);

               gdot = bqm * nu;
               dgdot_dtau = bqm * (nu * nu) * (dnu_w_dtau * ndwdiv2 + dnu_r_dtau * ndrdiv2);
               //
               double temp = gdot * copysign(gdot, tau) * ndwdiv2;
               // neglect difference in at_0 versus t_frac for dgdot_dg evaluation
               if (withGAthermal) {
                  dgdot_dg = -temp * dnu_w_dtau; // opposite sign as signed gdot
               }
               else {
                  dgdot_dg = -temp * dnu_w_dg; // opposite sign as signed gdot
               }
#if MORE_DERIVS
               // The reference rate is a bit different for the orowonian
               // framework then the previous version
               dgdot_dgamo = temp * (nu_w / (l_bar * _fD));
               dgdot_dmu = temp * dnu_w_dmu;
               dgdot_dtK = temp * dnu_w_dtK;
#endif
               //
               temp = gdot * copysign(gdot, tau) * ndrdiv2;
               if (withGAthermal) {
                  dgdot_dg = -temp * dnu_r_dtau + dgdot_dg; // opposite sign as signed gdot
               }
#if MORE_DERIVS
               // There reference value here is just the shear speed...
               dgdot_dgamr = temp * (nu_r / _shear_speed);
               dgdot_dtK = dgdot_dtK + temp * dnu_r_dtK;
#endif
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
                 int outputLevel = 0) const
         {
            // do not yet both with l_overdriven and setting-to-saturation machinery as in Fortran coding

            // update is done on log(h) -- h treated as a nomralized (unitless) dislocation density
            // double log_hs_u[SlipGeom::nslip * 2];
            double ihs_o[SlipGeom::nslip * 2];
            double nu[SlipGeom::nslip];
            for (int i = 0; i < _nslip * 2; i++) {
               if (i < _nslip) {
                  const double div = perSS ? fmax(hs_o[i], _hdn_min) * _berg_mag[i] :
                                     fmax(hs_o[i], _hdn_min) * _berg_mag[0];
                  nu[i] = gdot[i] / (div);
               }
               ihs_o[i] = fmax(hs_o[i], _hdn_min);
            }

            int nFEvals = updateHN<KineticsOrowanD>(this,
                                                   &hs_u[0], &ihs_o[0], dt, nu,
                                                   outputLevel);
            // We need to check that none of our solutions became negative
            // If we did obtain something negative then we should abort
            // It means our time step was too large for this step.
            // If this is not desirable / possible then we should probably
            // do a terrible hack and cut the dt by some factor resolve things by
            // assuming a constant slip rate during the time step, and then
            // evolve the dd content. We would get a solution, but it wouldn't necessarily
            // be correct.
            for (int i = 0; i < 2 * _nslip; i++) {
               if(hs_u[i] < zero) {
                  ECMECH_FAIL(__func__, "Solver returned negative dislocation values!");
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
            for (int i = 0; i < _nslip * 2; i++) {
               evolVals[i] = nu[i];
            }
         }

         __ecmech_hdev__
         inline
         void
         getSdotN( double* sdot,
                   double* dsdot_ds,
                   const double* const h,
                   const double* const evolVals
                   ) const
         {
            // Hopefully, the compiler is pretty smart here and is able to optimize these
            // loops as if we're using the templated values. Since, this is essentially
            // a constexpr for these variables.
            const int nslip = SlipGeom::nslip;
            const int JDIM = 2;
            const int nDimSys = 2 * SlipGeom::nslip;

            double forest_dis[nslip];
            vecsVMa<SlipGeom::nslip>(&forest_dis[0], &_a_mat[0], &h[nslip]);

            for (int iM = 0; iM < nslip; iM++) {
               const double sqrt_fd = sqrt(forest_dis[iM]);
               const double q_dmult = _c_mult * sqrt_fd * h[iM] * evolVals[iM];
               const double q_dtrap = _c_trap * sqrt_fd * h[iM] * evolVals[iM];
               // This could become a very large number and could become problematic
               // later on. Do we want to cap it at some large value?
               // Although, it might be that this is only a problem if q and qM are defined
               // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
               const double q_dann = _c_ann * _d_ann * h[iM] * h[iM + nslip] * evolVals[iM];
               // mobile dislocation density rate of change
               sdot[iM] = q_dmult - q_dtrap - q_dann;
               // total dislocation density rate of change
               sdot[iM + nslip] = q_dmult - q_dann;
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
                  const double sqrt_fd = sqrt(forest_dis[iM]);
                  const double q_dmult_dtrap = (_c_mult - _c_trap) * sqrt_fd;
                  // Although, it might be that this is only a problem if q and qM are defined
                  // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
                  const double q_dann = _c_ann * _d_ann * h[iM + nslip];
                  dsdot_ds_view(iM, iM) = evolVals[iM] * (q_dmult_dtrap - q_dann);
               }

               // dq/dqM portion of dsdot_ds
               for (int iT = 0; iT < nslip; iT++) {
                  const double sqrt_fd = sqrt(forest_dis[iT]);
                  const double q_dmult = _c_mult * sqrt_fd;
                  // This could become a very large number and could become problematic
                  // later on. Do we want to cap it at some large value?
                  // Although, it might be that this is only a problem if q and qM are defined
                  // with units 1/m^2 rather than 1/mm^2 or 1/micron^2
                  const double q_dann = _c_ann * _d_ann * h[iT + nslip];
                  dsdot_ds_view(iT + nslip, iT) =  evolVals[iT] * (q_dmult - q_dann);
               }

               RAJA::View<const double, RAJA::Layout<JDIM> > amat(&_a_mat[0], nslip, nslip);
               // dq/dq portion of dsdot_ds
               for (int iT = 0; iT < nslip; iT++) {
                  for (int jT = 0; jT < nslip; jT++) {
                     // First, terms found only on the diagonal of this submatrix
                     const double q_dann = (iT == jT) ? (_c_ann * _d_ann) : ecmech::zero;
                     const double ifact = ecmech::onehalf / sqrt(forest_dis[iT]);
                     const double q_dmult = _c_mult * amat(iT, jT) * ifact;

                     dsdot_ds_view(iT + nslip, jT + nslip) = h[iT] * evolVals[iT] * (q_dmult - q_dann);
                  }
               }

               // dM/dq portion of dsdot_dt
               for (int iT = 0; iT < nslip; iT++) {
                  for (int jT = 0; jT < nslip; jT++) {
                     const double q_dann = (iT == jT) ? (_c_ann * _d_ann) : ecmech::zero;

                     const double ifact = ecmech::onehalf / sqrt(forest_dis[iT]);
                     const double q_dmult_dtrap = (_c_mult - _c_trap) * amat(iT, jT) * ifact;

                     dsdot_ds_view(iT, jT + nslip) = h[iT] * evolVals[iT] * (q_dmult_dtrap - q_dann);
                  }
               }
            } // if dsdot_ds
         }
   }; // class KineticsOrowanD
} // namespace ecmech

#endif // ECMECH_KINETICS_OROWAND_H
