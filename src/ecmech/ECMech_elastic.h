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
        static const int nParams = 3;

        // constructor and destructor
        __ecmech_hdev__
        inline ThermoElastNCubic() : _K_bulkMod(-1.0), _K_gmod(-1.0) {};
        __ecmech_hdev__
        inline ~ThermoElastNCubic() {};

        __ecmech_host__
        inline void setParams(const std::vector<double> & params // const double* const params
                                ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            _c11 = *parsIt; ++parsIt;
            _c12 = *parsIt; ++parsIt;
            _c44 = *parsIt; ++parsIt;
            //
            assert((parsIt - params.begin()) == nParams);

            _K_diag[0] = _c11 - _c12;
            _K_diag[1] = _c11 - _c12;
            _K_diag[2] = two * _c44;
            _K_diag[3] = two * _c44;
            _K_diag[4] = two * _c44;
            double K_vecds_s = _c11 + two * _c12;
            _K_bulkMod = onethird * K_vecds_s;
            _K_gmod = (two * _c11 - two * _c12 + six * _c44) * 0.2; // average of _K_diag entries
        }

        __ecmech_host__
        inline void getParams(std::vector<double> & params
                                ) const {
    #ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
    #endif

            params.push_back(_c11);
            params.push_back(_c12);
            params.push_back(_c44);

    #ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
    #endif
        }

        __ecmech_hdev__
        inline
        void eval(double* const T_vecds,
                    const double* const Ee_vecds,
                    double, // tK
                    double p_EOS,
                    double // eVref
                    ) const {
            double ln_J = sqr3 * Ee_vecds[iSvecS]; // vecds_s_to_trace
            double J = exp(ln_J);
            double Ts_bulk = -sqr3 * J * p_EOS;

            vecsVAdiagB<ntvec>(T_vecds, _K_diag, Ee_vecds);
            T_vecds[iSvecS] = Ts_bulk; // _K_vecds_s * Ee_vecds(SVEC)
        }

        /**
            * dT_deps[0:ntvec-1,:]^T * A, for non-square A[ntvec,p] (with p likely being nSlip)
            * so that even if T_vecds[iSvecS] depends on Ee_vecds, that is not in the result
            *
            * combines calls to elawn_T_dif and eval_dtaua_deps_n
            *
            * for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
            */
        __ecmech_hdev__
        inline
        void multDTDepsT(double* const P, // ntvec*p
                            const double* const A, // ntvec*p
                            double a_V_ri,
                            int p) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double dTdepsThis = _K_diag[iTvec] * a_V_ri;
                for (int iP = 0; iP < p; ++iP) {
                    int ii = ECMECH_NM_INDX(iTvec, iP, ecmech::ntvec, p);
                    P[ii] = dTdepsThis * A[ii];
                }
            }
        }

        __ecmech_hdev__
        inline
        void getCauchy(double* const sigC_vecds_lat,
                        const double* const T_vecds,
                        double detVi) const
        {
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                sigC_vecds_lat[iSvec] = detVi * T_vecds[iSvec];
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
        __ecmech_hdev__
        template<size_t N=ecmech::ntvec, size_t M=ecmech::ntvec>
        inline
        void multCauchyDif(double* const M6,
                            const double* const A,
                            double detVi,
                            double a_V_ri
                            ) const {
            // CALL vecds_s_to_trace(tr_ln_V, s_meas%Ee_vecds(SVEC))
            // detV = DEXP(tr_ln_V)
            // detVi = one / detV

            // dsigC_de(:,:) = detVi * s_meas%dT_deps(:,:)
            // for cubic, dT_deps is diag(K_diag * a_V%ri) (symmetric) ; dT_deps[iSvecS,:] = 0
            // M65_ij = dd_ii A_ij
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double vFact = detVi * a_V_ri * _K_diag[iTvec];
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
            if (_K_bulkMod <= 0.0) {
                ECMECH_FAIL(__func__, "bulk modulus negative -- not initialized?");
            }
            return _K_bulkMod;
        }

        __ecmech_hdev__
        inline
        double getGmod(double, // tK
                        double, // p_EOS
                        double // eVref
                        ) const {
            if (_K_gmod <= 0.0) {
                ECMECH_FAIL(__func__, "effective shear modulus negative -- not initialized?");
            }
            return _K_gmod;
        }

        private:
        double _c11, _c12, _c44;
        double _K_diag[ecmech::ntvec];
        double _K_bulkMod, _K_gmod;
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
    * and just store the one non-zero as _g_vecd2
    *
    */
    class ThermoElastNHexag
    {
        public:
        static const int nParams = 6;

        // constructor and destructor
        __ecmech_hdev__
        inline ThermoElastNHexag() : _K_bulkMod(-1.0), _K_gmod(-1.0) {};
        __ecmech_hdev__
        inline ~ThermoElastNHexag() {};

        __ecmech_host__
        inline void setParams(const std::vector<double> & params // const double* const params
                                ) {
            std::vector<double>::const_iterator parsIt = params.begin();

            _c11 = *parsIt; ++parsIt;
            _c12 = *parsIt; ++parsIt;
            _c13 = *parsIt; ++parsIt;
            _c33 = *parsIt; ++parsIt;
            _c44 = *parsIt; ++parsIt;
            //
            _g_vecd2 = *parsIt; ++parsIt;
            //
            assert((parsIt - params.begin()) == nParams);

            _K_diag[0] = _c11 - _c12;
            _K_diag[1] = _c11 * onethird + _c12 * onethird - fourthirds * _c13 + twothird * _c33;
            _K_diag[2] = _c11 - _c12;
            _K_diag[3] = two * _c44;
            _K_diag[4] = two * _c44;
            double K_vecds_s = twothird * _c11 + twothird * _c12 + fourthirds * _c13 + _c33 * onethird;
            _K_sdax3 = sqr2 * (-_c11 - _c12 + _c13 + _c33) * onethird;
            _K_bulkMod = onethird * K_vecds_s;
            //
            // _K_gmod below ignores the _K_sdax3 contribution, but it is just meant to be approximate anyway
            _K_gmod = 0.5 * 0.2 * vecsssum<ecmech::ntvec>(_K_diag); // 0.5 * (average of _K_diag entries)
        }

        __ecmech_host__
        inline void getParams(std::vector<double> & params
                                ) const {
    #ifdef ECMECH_DEBUG
            // do not clear params in case adding to an existing set
            int paramsStart = params.size();
    #endif

            params.push_back(_c11);
            params.push_back(_c12);
            params.push_back(_c13);
            params.push_back(_c33);
            params.push_back(_c44);
            //
            params.push_back(_g_vecd2);

    #ifdef ECMECH_DEBUG
            assert((params.size() - paramsStart) == nParams);
    #endif
        }

        __ecmech_hdev__
        inline
        void eval(double* const T_vecds,
                    const double* const Ee_vecds,
                    double, // tK
                    double p_EOS,
                    double eVref
                    ) const {
            double ln_J = sqr3 * Ee_vecds[iSvecS]; // vecds_s_to_trace
            double J = exp(ln_J);
            double Ts_bulk = -sqr3 * J * p_EOS;

            vecsVAdiagB<ntvec>(T_vecds, _K_diag, Ee_vecds);
            T_vecds[iSvecS] = Ts_bulk; // _K_vecds_s * Ee_vecds(SVEC)

            T_vecds[iTvecHex] += _K_sdax3 * Ee_vecds[iSvecS];
            T_vecds[iSvecS] += _K_sdax3 * Ee_vecds[iTvecHex];

            // anisotropic Gruneisen contribution; pressure part of Gruneisen tensor contribution should already be in p_EOS
            // CALL eos_eval_e_Csdev(Cauchy_eos_vecd, eVref, J, &
            // & i_eos_model, eos_const)
            // -(Gamma' + a' * mu) * eVref // but do not do a'*mu part
            // Cauchy_eos_vecd(:) = -eos_const(4:8) * eVref
            // T_vecds(1:TVEC) = T_vecds(1:TVEC) + J * Cauchy_eos_vecd(:)
            T_vecds[iTvecHex] += J * (-_g_vecd2 * eVref);
        }

        /**
            * multDTDepsT ends up looking the same as in the cubic case because _K_sdax3 does not enter
            */
        __ecmech_hdev__
        inline
        void multDTDepsT(double* const P, // ntvec*p
                            const double* const A, // ntvec*p
                            double a_V_ri,
                            int p) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double dTdepsThis = _K_diag[iTvec] * a_V_ri;
                for (int iP = 0; iP < p; ++iP) {
                    int ii = ECMECH_NM_INDX(iTvec, iP, ecmech::ntvec, p);
                    P[ii] = dTdepsThis * A[ii];
                }
            }
        }

        __ecmech_hdev__
        inline
        void getCauchy(double* const sigC_vecds_lat,
                        const double* const T_vecds,
                        double detVi) const
        {
            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                sigC_vecds_lat[iSvec] = detVi * T_vecds[iSvec];
            }
        }

        __ecmech_hdev__
        template<size_t N=ecmech::ntvec, size_t M=ecmech::ntvec>
        inline
        void multCauchyDif(double* const M6,
                            const double* const A,
                            double detVi,
                            double a_V_ri
                            ) const {
            for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
                double vFact = detVi * a_V_ri * _K_diag[iTvec];
                for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                    M6[ECMECH_NN_INDX(iTvec, jTvec, ecmech::nsvec)] = vFact * A[ECMECH_NN_INDX(iTvec, jTvec, ecmech::ntvec)];
                }
            }

            // M6[iSvecS,:] = dsigC_de[iSvecS, iTvecHex] * A[iTvecHex,:] // for hexagonal specifically
            // dsigC_de[iTvecHex, iSvecS] does not end up getting used
            {
                double vFact = detVi * a_V_ri * _K_sdax3;
                for (int jTvec = 0; jTvec < ecmech::ntvec; ++jTvec) {
                    M6[ECMECH_NN_INDX(iSvecS, jTvec, ecmech::nsvec)] = vFact * A[ECMECH_NN_INDX(iTvecHex, jTvec, ecmech::ntvec)];
                }
            }

            for (int iSvec = 0; iSvec < ecmech::nsvec; ++iSvec) {
                M6[ECMECH_NN_INDX(iSvec, iSvecS, ecmech::nsvec)] = 0.0;
            }
        }

        __ecmech_hdev__
        inline
        double getBulkMod( ) const {
            if (_K_bulkMod <= 0.0) {
                ECMECH_FAIL(__func__, "bulk modulus negative -- not initialized?");
            }
            return _K_bulkMod;
        }

        __ecmech_hdev__
        inline
        double getGmod(double, // tK
                        double, // p_EOS
                        double // eVref
                        ) const {
            if (_K_gmod <= 0.0) {
                ECMECH_FAIL(__func__, "effective shear modulus negative -- not initialized?");
            }
            return _K_gmod;
        }

        private:
        double _c11, _c12, _c13, _c33, _c44;
        double _K_sdax3;
        double _g_vecd2;
        double _K_diag[ecmech::ntvec];
        double _K_bulkMod, _K_gmod;
        static const int iTvecHex = 1;
    };

}
}