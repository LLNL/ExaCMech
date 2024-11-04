// -*-c++-*-

#ifndef ECMECH_kinetics_H
#define ECMECH_kinetics_H

#include "SNLS_TrDLDenseG.h"

#include "ECMech_core.h"
#include "ECMech_util.h"

/*
 * Kinetics models include both the slip system kinetics and 'hardness'
 * state evolution kinetics. Those two could be separated to increase
 * modularity, but then it might be a pain to keep the interactions
 * among them both generally flexible and simple to manage. So the
 * choice has been made to group them into a single class.
 *
 * Specific model cases are down below, through include files. Cases are expected to provide traits:
 *
 *    nH, nParams, nVals, nEvolVals
 *
 * and member functions:
 *
 *  constructor(int nslip) ;
 *
 *  void setParams( const std::vector<double> & params ) ;
 *
 *  void getParams( std::vector<double> & params ) const ;
 *
 *  void getHistInfo(std::vector<std::string> & names,
 *                   std::vector<double>       & init,
 *                   std::vector<bool>        & plot,
 *                   std::vector<bool>        & state) const ;
 *
 *  double getVals( double* const vals,
 *                  double p,
 *                  double tkelv,
 *                  const double* const h_state
 *                  ) const ;
 *
 *  void
 *  evalGdots( double* const gdot,
 *             double* const dgdot_dtau,
 *             const double* const tau,
 *             const double* const vals
 *             ) const ;
 *
 * This returns the number of function evaluations with a negative
 * value signaling a failed solve
 *
 * hs_u and hs_o are the hardening states
 * gdot is the slip rate
 * hvals are any additional variables that might required by the internal solver
 *
 *  int
 *  updateH( double* const hs_u,
             const double* const hs_o,
             double dt,
             const double* const gdot,
             const double* const hvals,
             double tkelv,
             int outputLevel = 0) const ;
 *
 *  void
 *  getEvolVals( double* const evolVals,
 *               const double* const gdot
 *               ) const ;
 *
 * If using updateH1, then provide member function getSdot1,
 * else if using getSdotN, then provide member function getSdotN.
 *
 *  void
 *  getSdot1( double &sdot,
 *            double &dsdot_ds,
 *            double h,
 *            const double* const evolVals,
 *            double tkelv) const const;
 *
 * The incoming dsdot_ds ptr should be checked to make sure its not a nullptr,
 * and if so the calculations relevant to dsdot_ds should be skipped.
 *
 * hvals are any additional variables that might required by the internal solver
 *
 * void
 * getSdotN( double* sdot,
 *           double* dsdot_ds,
 *           const double* const h,
 *           const double* const evolVals,
 *           const double* const hvals,
 *           double tkelv) const;
 */

namespace ecmech {
   /*
    * State update solver for cases in which there is a single hardness state variable.
    */
   template<class Kinetics>
   class Kinetics_H1Problem
   {
      public:
         static constexpr int nDimSys = Kinetics::nH;

         // constructor
         __ecmech_hdev__
         Kinetics_H1Problem(const Kinetics* const kinetics,
                            double h_o,
                            double dt,
                            const double* const evolVals,
                            double tkelv) :
            m_kinetics(kinetics), m_h_o(h_o), m_dt(dt), m_evolVals(evolVals), m_tkelv(tkelv)
         {
            m_x_scale = fmax(m_h_o, 1.0); // TO_DO -- generalize this to not max with 1
            // NOTE : see comment below about changing Jacobian calculation if m_res_scale != one / s_scale
            m_res_scale = one / m_x_scale;
         }

         // deconstructor
         __ecmech_hdev__
         ~Kinetics_H1Problem() {}

         __ecmech_hdev__
         inline
         double getHn(const double* const x) const {
            return m_h_o + x[0] * m_x_scale;
         }

         __ecmech_hdev__
         inline
         bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            bool doComputeJ = (Jacobian != nullptr);

            double h_delta = x[0] * m_x_scale;
            double h = m_h_o + h_delta;

            double sdot, dsdot_ds;
            m_kinetics->getSdot1(sdot, dsdot_ds, h, m_evolVals, m_tkelv);

            resid[0] = (h_delta - sdot * m_dt) * m_res_scale;

            if (doComputeJ) {
               // The below is based on the assumption that m_res_scale = 1/m_x_scale
               // if this were to change in the future version than this would need to become
               // Jacobian[0] = (one - dsdot_ds * m_dt) * m_res_scale * m_x_scale;
               Jacobian[0] = (one - dsdot_ds * m_dt);
            }

            return true;
         } // computeRJ

      private:
         const Kinetics* m_kinetics;
         const double m_h_o, m_dt;
         const double* const m_evolVals;
         const double m_tkelv;
         double m_x_scale, m_res_scale;
   }; // class Kinetics_H1Problem

   /*
    * Helper function to run the state update solver for cases in which there is a single hardness state variable.
    If the solve did not fail then it returns the number of function evaluations
    the solver required. However, if it did fail then it returns a -1 value
    to signal failure rather than throwing failures. This is done as we
    throw exceptions on GPUs so we need a way to handle failures.
    */
   template<class Kinetics, bool relaxed_solver = false>
   __ecmech_hdev__
   inline
   int
   updateH1(const Kinetics* const kinetics,
            double &hs_n,
            double hs_o,
            double dt,
            const double* const gdot,
            double tkelv,
            int outputLevel = 0)
   {
      double evolVals[Kinetics::nEvolVals];
      kinetics->getEvolVals(evolVals, gdot);

      Kinetics_H1Problem<Kinetics> prob(kinetics, hs_o, dt, evolVals, tkelv);
      snls::SNLSTrDlDenseG<Kinetics_H1Problem<Kinetics> > solver(prob);

      snls::TrDeltaControl deltaControl;
      deltaControl._deltaInit = 1e0;
      {
         int maxIter = 100;
         double tolerance = 1e-10;
         solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
      }

      for (int iX = 0; iX < prob.nDimSys; ++iX) {
         solver._x[iX] = 0e0;
      }

      snls::SNLSStatus_t status = solver.solve( );

      int nFevals = solver.getNFEvals();
      if (status != snls::converged) {
         snls::SNLSStatus_t status2 = status;
         if constexpr(relaxed_solver) {
            {
               int maxIter = 100;
               double tolerance = 1e-9;
               solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
            }
            for (int iX = 0; iX < prob.nDimSys; ++iX) {
               solver._x[iX] = 0e0;
            }
            status2= solver.solve();
            nFevals = solver.getNFEvals();
         }
         if (status2 != snls::converged) {
            nFevals = -1;
         }
      }

      hs_n = prob.getHn(solver._x);

      return nFevals;
   } // updateH1

   /*
    * State update solver for cases in which there are multiple hardness state variable.
    */
   template<class Kinetics>
   class Kinetics_HNProblem
   {
      public:
         static constexpr int nDimSys = Kinetics::nH;

         // constructor
         __ecmech_hdev__
         Kinetics_HNProblem(const Kinetics* const kinetics,
                            const double* const h_o,
                            double dt,
                            const double* const evolVals,
                            const double* const hvals,
                            double tkelv) :
            m_kinetics(kinetics), m_h_o(h_o), m_dt(dt), m_evolVals(evolVals), m_hvals(hvals), m_tkelv(tkelv)
         {
            for (int i = 0; i < nDimSys; i++) {
               m_x_scale[i] = fmax(m_h_o[i], 1.0); // TO_DO -- generalize this to not max with 1
               // NOTE : see comment below about changing Jacobian calculation if m_res_scale != one / s_scale
               m_res_scale[i] = one / m_x_scale[i];
            }
         }

         // deconstructor
         __ecmech_hdev__
         ~Kinetics_HNProblem() {}

         __ecmech_hdev__
         inline
         void getHn(double* h, const double* const x) const {
            for (int i = 0; i < nDimSys; i++) {
               h[i] = m_h_o[i] + x[i] * m_x_scale[i];
            }
         }

         __ecmech_hdev__
         inline
         bool computeRJ(double* const resid,
                        double* const Jacobian,
                        const double* const x) {
            double h[nDimSys];
            for (int i = 0; i < nDimSys; i++) {
               h[i] = m_h_o[i] + x[i] * m_x_scale[i];
            }

            double sdot[nDimSys];
            // The dsdot_ds portion of the Jacobian is set in here if it was provided
            m_kinetics->getSdotN(sdot, Jacobian, h, m_evolVals, m_hvals, m_tkelv);

            for (int i = 0; i < nDimSys; i++) {
               resid[i] = (x[i] * m_x_scale[i] - sdot[i] * m_dt) * m_res_scale[i];
            }

            if (Jacobian) {
               // Multiply dsdot_ds terms by the negative outer product of x_scale and res_scale and dt
               for (int i = 0; i < nDimSys; i++) {
                  for (int j = 0; j < nDimSys; j++) {
                     Jacobian[ECMECH_NN_INDX(i, j, nDimSys)] *= -m_x_scale[j] * m_res_scale[i] * m_dt;
                  }
               }

               // Now add in the identity term
               // The below is based on the assumption that m_res_scale = 1/m_x_scale
               // if this were to change in the future version than this would need to become
               // Jacobian[ECMECH_NN_INDX(i, i, nDimSys)] += ecmech::one * m_x_scale[i] * _r_scale[i]
               for (int i = 0; i < nDimSys; i++) {
                  Jacobian[ECMECH_NN_INDX(i, i, nDimSys)] += ecmech::one;
               }
            } // if Jacobian

            return true;
         } // computeRJ

      private:
         const Kinetics* m_kinetics;
         const double* const m_h_o;
         const double m_dt;
         const double* const m_evolVals;
         const double* const m_hvals;
         const double m_tkelv;
         double m_x_scale[nDimSys], m_res_scale[nDimSys];
   }; // class Kinetics_HNProblem

   /*
    * Helper function to run the state update solver for cases in which there are multiple hardness state variables.
    If the solve did not fail then it returns the number of function evaluations
    the solver required. However, if it did fail then it returns a -1 value
    to signal failure rather than throwing failures. This is done as we
    throw exceptions on GPUs so we need a way to handle failures.
    */
   template<class Kinetics, bool relaxed_solver = false>
   __ecmech_hdev__
   inline
   int
   updateHN(const Kinetics* const kinetics,
            double* hs_n,
            const double* const hs_o,
            double dt,
            const double* const gdot,
            const double* const hvals,
            double tkelv,
            int outputLevel = 0)
   {
      double evolVals[Kinetics::nEvolVals];
      kinetics->getEvolVals(evolVals, gdot);

      Kinetics_HNProblem<Kinetics> prob(kinetics, hs_o, dt, evolVals, hvals, tkelv);
      snls::SNLSTrDlDenseG<Kinetics_HNProblem<Kinetics> > solver(prob);

      snls::TrDeltaControl deltaControl;
      deltaControl._deltaInit = 1e0;
      {
         int maxIter = 100;
         double tolerance = 1e-10;
         solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
      }

      for (int iX = 0; iX < prob.nDimSys; ++iX) {
         solver._x[iX] = 0e0;
      }

      snls::SNLSStatus_t status = solver.solve( );

      int nFevals = solver.getNFEvals();
      if (status != snls::converged) {
         snls::SNLSStatus_t status2 = status;
         if constexpr(relaxed_solver) {
            {
               int maxIter = 100;
               double tolerance = 1e-9;
               solver.setupSolver(maxIter, tolerance, &deltaControl, outputLevel);
            }
            for (int iX = 0; iX < prob.nDimSys; ++iX) {
               solver._x[iX] = 0e0;
            }
            status2= solver.solve();
            nFevals = solver.getNFEvals();
         }
         if (status2 != snls::converged) {
            nFevals = -1;
         }
      }

      prob.getHn(hs_n, solver._x);

      return nFevals;
   } // updateHN
} // namespace ecmech

/*
 * And now, specific model cases.
 */
#include "ECMech_kinetics_KMBalD.h"
#include "ECMech_kinetics_VocePL.h"
#include "ECMech_kinetics_OrowanD.h"
#include "ECMech_kinetics_BCCMD.h"

#endif // ECMECH_KINETICS_H
