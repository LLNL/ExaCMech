#include "ECMech_core.h"
#include "ECMech_util.h"
#include "ECMech_eosSimple.h"


namespace ecmech {
namespace evptn {
    /**
    * for cubic cyrstal symmetry
    *
    * in Fortran mdef coding, corresponds to cem%l_lin_lnsd
    *
    */
    class ThermoElastNCubic
    {
        public:
        static constexpr int nParams = 3;

        // constructor and destructor
        __ecmech_hdev__
        inline ThermoElastNCubic() : m_bulk_modulus(-1.0), m_shear_modulus(-1.0) {}

        ~ThermoElastNCubic() = default;

         __ecmech_hdev__
         ThermoElastNCubic(const double* const params) {
            setParams(params);
         }

        __ecmech_host__
        inline void setParams(const std::vector<double> & params) {
            setParams(params.data());
        }


         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            m_c11 = *parsIt; ++parsIt;
            m_c12 = *parsIt; ++parsIt;
            m_c44 = *parsIt; ++parsIt;
            //
#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
            m_K_diag[0] = m_c11 - m_c12;
            m_K_diag[1] = m_c11 - m_c12;
            m_K_diag[2] = two * m_c44;
            m_K_diag[3] = two * m_c44;
            m_K_diag[4] = two * m_c44;
            double K_vecds_s = m_c11 + two * m_c12;
            m_bulk_modulus = onethird * K_vecds_s;
            m_shear_modulus = (two * m_c11 - two * m_c12 + six * m_c44) * 0.2; // average of m_K_diag entries
        }

        __ecmech_host__
        inline void getParams(std::vector<double> & params
                                ) const {
    #ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
    #endif

            params.push_back(m_c11);
            params.push_back(m_c12);
            params.push_back(m_c44);

    #ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
    #endif
        }

        __ecmech_hdev__
        inline
        void eval(double* const kirchoff,
                    const double* const elast_dev_press_vec,
                    double, // tkelv
                    double pressure_EOS,
                    double // energy_vol_ref
                    ) const {
            double ln_J = sqr3 * elast_dev_press_vec[iSvecS]; // vecds_s_to_trace
            double J = exp(ln_J);
            double kirchoff_pressure = -sqr3 * J * pressure_EOS;

            vecsVAdiagB<ntvec>(kirchoff, m_K_diag, elast_dev_press_vec);
            kirchoff[iSvecS] = kirchoff_pressure; // _K_vecds_s * elast_dev_press_vec(SVEC)
        }

        /**
            * dT_deps[0:ntvec-1,:]^T * A, for non-square A[ntvec,p] (with p likely being nSlip)
            * so that even if kirchoff[iSvecS] depends on elast_dev_press_vec, that is not in the result
            *
            * combines calls to elawn_T_dif and eval_dtaua_deps_n
            *
            * for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
            */
        __ecmech_hdev__
        inline
        void multDTDepsT(double* const P, // ntvec*p
                            const double* const A, // ntvec*p
                            double inv_a_vol,
                            int p) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double dTdepsThis = m_K_diag[iTvec] * inv_a_vol;
                for (int iP = 0; iP < p; ++iP) {
                    int ii = ECMECH_NM_INDX(iTvec, iP, ecmech::ntvec, p);
                    P[ii] = dTdepsThis * A[ii];
                }
            }
        }

        __ecmech_hdev__
        inline
        void getCauchy(double* const cauchy_xtal,
                        const double* const kirchoff,
                        double inv_det_v_e) const
        {
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                cauchy_xtal[iSvec] = inv_det_v_e * kirchoff[iSvec];
            }
        }

        /**
            * like dsigC_de * A, with disgC_de[nsvec,ntvec] having come from elawn_Cauchy_dif
            * for A[ntvec,ntvec]
            *
            * NOTE : dsigC_de[nsvec,ntvec] with nsvec in the first dimension
            * because in general distorational deformation can produce
            * pressure -- for example in materials with hexagonal symmetry,
            * even if it does not happen in cubic symmetry
            *
            * NOTE : M6[nsvec,nsvec] with nsvec in the second dimension
            * (instead of ntvec) to make things easier elsewhere
            */
        template<size_t N=ecmech::ntvec, size_t M=ecmech::ntvec>
        __ecmech_hdev__
        inline
        void multCauchyDif(double* const M6,
                            const double* const A,
                            double inv_det_v_e,
                            double inv_a_vol
                            ) const {
            // CALL vecds_s_to_trace(tr_ln_V, s_meas%elast_dev_press_vec(SVEC))
            // det_v_e = DEXP(tr_ln_V)
            // inv_det_v_e = one / det_v_e

            // dsigC_de(:,:) = inv_det_v_e * s_meas%dT_deps(:,:)
            // for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
            // M65_ij = dd_ii A_ij
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double vFact = inv_det_v_e * inv_a_vol * m_K_diag[iTvec];
                for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                    M6[ECMECH_NN_INDX(iTvec, jTvec, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvec, jTvec, N, M)];
                }
            }

            for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                M6[ECMECH_NN_INDX(iSvecS, jTvec, ecmech::nsvec)] = 0.0;
            }

            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                M6[ECMECH_NN_INDX(iSvec, iSvecS, ecmech::nsvec)] = 0.0;
            }
        }

        __ecmech_hdev__
        inline
        double getBulkMod( ) const {
            if (m_bulk_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "bulk modulus negative -- not initialized?");
            }
            return m_bulk_modulus;
        }

        __ecmech_hdev__
        inline
        double getGmod(double, // tkelv
                        double, // pressure_EOS
                        double // energy_vol_ref
                        ) const {
            if (m_shear_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "effective shear modulus negative -- not initialized?");
            }
            return m_shear_modulus;
        }

        private:
        double m_c11, m_c12, m_c44;
        double m_K_diag[ecmech::ntvec];
        double m_bulk_modulus, m_shear_modulus;
    };

    /**
    * for hexagonal cyrstal symmetry
    *
    * in Fortran mdef coding, corresponds to cem%l_lin_lnsd, cem%l_h
    *
    * Gruneisen gamma is diag(g_a,g_a,g_b)
    * g_vecd is
    *    (g11-g22)/sqrt(2.) = 0
    *    (2. * g33 - g11 - g22)/sqrt(6.) = 2.0 * (g_b - g_a) / sqrt(6.)
    *    sqrt(2.) * g12 = 0
    *    sqrt(2.) * g13 = 0
    *    sqrt(2.) * g23 = 0
    * and just store the one non-zero as m_g_vecd2
    *
    */
    class ThermoElastNHexag
    {
        public:
        static constexpr int nParams = 6;

        // constructor and destructor
        __ecmech_hdev__
        inline ThermoElastNHexag() : m_bulk_modulus(-1.0), m_shear_modulus(-1.0) {}

        ~ThermoElastNHexag() = default;

         __ecmech_hdev__
         ThermoElastNHexag(const double* const params) {
            setParams(params);
         }

        __ecmech_host__
        inline void setParams(const std::vector<double> & params) {
            setParams(params.data());
        }


         __ecmech_hdev__
         inline
         void setParams(const double* const params) {
            const double* parsIt = params;

            m_c11 = *parsIt; ++parsIt;
            m_c12 = *parsIt; ++parsIt;
            m_c13 = *parsIt; ++parsIt;
            m_c33 = *parsIt; ++parsIt;
            m_c44 = *parsIt; ++parsIt;
            //
            m_g_vecd2 = *parsIt; ++parsIt;
            //
#if defined(ECMECH_DEBUG)
            int iParam = parsIt - params;
            if (iParam != nParams) {
               ECMECH_FAIL(__func__, "iParam != nParams");
            }
#endif
            m_K_diag[0] = m_c11 - m_c12;
            m_K_diag[1] = m_c11 * onethird + m_c12 * onethird - fourthirds * m_c13 + twothird * m_c33;
            m_K_diag[2] = m_c11 - m_c12;
            m_K_diag[3] = two * m_c44;
            m_K_diag[4] = two * m_c44;
            double K_vecds_s = twothird * m_c11 + twothird * m_c12 + fourthirds * m_c13 + m_c33 * onethird;
            m_K_sdax3 = sqr2 * (-m_c11 - m_c12 + m_c13 + m_c33) * onethird;
            m_bulk_modulus = onethird * K_vecds_s;
            //
            // m_shear_modulus below ignores the m_K_sdax3 contribution, but it is just meant to be approximate anyway
            m_shear_modulus = 0.5 * 0.2 * vecsssum<ecmech::ntvec>(m_K_diag); // 0.5 * (average of m_K_diag entries)
        }

        __ecmech_host__
        inline void getParams(std::vector<double> & params
                                ) const {
    #ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
    #endif

            params.push_back(m_c11);
            params.push_back(m_c12);
            params.push_back(m_c13);
            params.push_back(m_c33);
            params.push_back(m_c44);
            //
            params.push_back(m_g_vecd2);

    #ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
    #endif
        }

        __ecmech_hdev__
        inline
        void eval(double* const kirchoff,
                    const double* const elast_dev_press_vec,
                    double, // tkelv
                    double pressure_EOS,
                    double energy_vol_ref
                    ) const {
            double ln_J = sqr3 * elast_dev_press_vec[iSvecS]; // vecds_s_to_trace
            double J = exp(ln_J);
            double kirchoff_pressure = -sqr3 * J * pressure_EOS;

            vecsVAdiagB<ntvec>(kirchoff, m_K_diag, elast_dev_press_vec);
            kirchoff[iSvecS] = kirchoff_pressure; // _K_vecds_s * elast_dev_press_vec(SVEC)

            kirchoff[iTvecHex] += m_K_sdax3 * elast_dev_press_vec[iSvecS];
            kirchoff[iSvecS] += m_K_sdax3 * elast_dev_press_vec[iTvecHex];

            // anisotropic Gruneisen contribution; pressure part of Gruneisen tensor contribution should already be in pressure_EOS
            // CALL eos_eval_e_Csdev(Cauchy_eos_vecd, energy_vol_ref, J, &
            // & i_eos_model, eos_const)
            // -(Gamma' + a' * mu) * energy_vol_ref // but do not do a'*mu part
            // Cauchy_eos_vecd(:) = -eos_const(4:8) * energy_vol_ref
            // kirchoff(1:TVEC) = kirchoff(1:TVEC) + J * Cauchy_eos_vecd(:)
            kirchoff[iTvecHex] += J * (-m_g_vecd2 * energy_vol_ref);
        }

        /**
            * multDTDepsT ends up looking the same as in the cubic case because m_K_sdax3 does not enter
            */
        __ecmech_hdev__
        inline
        void multDTDepsT(double* const P, // ntvec*p
                            const double* const A, // ntvec*p
                            double inv_a_vol,
                            int p) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double dTdepsThis = m_K_diag[iTvec] * inv_a_vol;
                for (int iP = 0; iP < p; ++iP) {
                    int ii = ECMECH_NM_INDX(iTvec, iP, ecmech::ntvec, p);
                    P[ii] = dTdepsThis * A[ii];
                }
            }
        }

        __ecmech_hdev__
        inline
        void getCauchy(double* const cauchy_xtal,
                        const double* const kirchoff,
                        double inv_det_v_e) const
        {
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                cauchy_xtal[iSvec] = inv_det_v_e * kirchoff[iSvec];
            }
        }

        template<size_t N=ecmech::ntvec, size_t M=ecmech::ntvec>
        __ecmech_hdev__
        inline
        void multCauchyDif(double* const M6,
                            const double* const A,
                            double inv_det_v_e,
                            double inv_a_vol
                            ) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double vFact = inv_det_v_e * inv_a_vol * m_K_diag[iTvec];
                for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                    M6[ECMECH_NN_INDX(iTvec, jTvec, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvec, jTvec, N, M)];
                }
            }

            // M6[iSvecS,:] = dsigC_de[iSvecS, iTvecHex] * A[iTvecHex,:] // for hexagonal specifically
            // dsigC_de[iTvecHex, iSvecS] does not end up getting used
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                M6[ECMECH_NN_INDX(iSvec, iSvecS, ecmech::nsvec)] = 0.0;
                M6[ECMECH_NN_INDX(iSvecS, iSvec, ecmech::nsvec)] = 0.0;
            }

            {
                double vFact = inv_det_v_e * inv_a_vol * m_K_sdax3;
                for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                    M6[ECMECH_NN_INDX(iSvecS, jTvec, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvecHex, jTvec, N, M)];
                }
            }
            {
                double vFact = inv_det_v_e * inv_a_vol * m_K_sdax3;
                for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                    M6[ECMECH_NN_INDX(jTvec, iSvecS, ecmech::nsvec)] = vFact * A[ECMECH_NM_INDX(iTvecHex, jTvec, N, M)];
                }
            }
        }

        __ecmech_hdev__
        inline
        double getBulkMod( ) const {
            if (m_bulk_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "bulk modulus negative -- not initialized?");
            }
            return m_bulk_modulus;
        }

        __ecmech_hdev__
        inline
        double getGmod(double, // tkelv
                        double, // pressure_EOS
                        double // energy_vol_ref
                        ) const {
            if (m_shear_modulus <= 0.0) {
                ECMECH_FAIL(__func__, "effective shear modulus negative -- not initialized?");
            }
            return m_shear_modulus;
        }

        private:
        double m_c11, m_c12, m_c13, m_c33, m_c44;
        double m_K_sdax3;
        double m_g_vecd2;
        double m_K_diag[ecmech::ntvec];
        double m_bulk_modulus, m_shear_modulus;
        static constexpr int iTvecHex = 1;
    };

}
}