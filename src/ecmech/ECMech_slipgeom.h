// -*-c++-*-

#ifndef ECMECH_SLIPGEOM_H
#define ECMECH_SLIPGEOM_H

#include "ECMech_core.h"
#include "ECMech_util.h"

namespace ecmech {
   // A version using constexpr could probably be created for the below
   __ecmech_hdev__
   static void
   fillFromMS(double* const P, // ntvec * nslip
              double* const Q, // nwvec * nslip
              const double* const mVecs, // nslip * ndim
              const double* const sVecs, // nslip * ndim
              int nslip)
   {
      for (int iSlip = 0; iSlip<nslip; ++iSlip) {

         const double* mVec = &(mVecs[iSlip * ecmech::ndim]);
         const double* sVec = &(sVecs[iSlip * ecmech::ndim]);
#ifndef NO_CHECKS
         if ( fabs(vecsyadotb<ecmech::ndim>(mVec,sVec)) > idp_eps_sqrt ) {
            ECMECH_FAIL(__func__, "internal error");
         }
#endif
         
         // CALL vec_x_vect_mn(crys%vecs(:,is),crys%vecm(:,is),crys%t_ref(:,:,is),DIMS,DIMS)
         double T_ref[ ecmech::ndim*ecmech::ndim ];
         vecsMaTb<ndim>(T_ref, sVec, mVec );

         double P_vecd[ ecmech::ntvec ];
         double Q_veccp[ ecmech::nwvec ];
         matToPQ(P_vecd, Q_veccp, T_ref);

         for (int iTvec = 0; iTvec < ecmech::ntvec; ++iTvec) {
            P[ECMECH_NM_INDX(iTvec, iSlip, ecmech::ntvec, nslip)] = P_vecd[iTvec];
         }

         for (int iWvec = 0; iWvec < ecmech::nwvec; ++iWvec) {
            Q[ECMECH_NM_INDX(iWvec, iSlip, ecmech::nwvec, nslip)] = Q_veccp[iWvec];
         }

         //
         // in some approaches, it is useful to form the outer product of P_vecd with itself, for tangent stiffness contributions
      }
   }

   class SlipGeomFCC
   {
      public:

         static const int nslip = 12;
         // constructor and destructor
         __ecmech_hdev__
         SlipGeomFCC()
         {
            // m = (/ sqr3i, sqr3i, sqr3i /)
            // s = (/ zero, sqr2i, -sqr2i /)
            //
            // do not yet bother with making slip systems from symmetry group -- just write them out
            const double mVecs[ nslip * ecmech::ndim ] = {
               sqr3i, sqr3i, sqr3i,
               sqr3i, sqr3i, sqr3i,
               sqr3i, sqr3i, sqr3i,
               -sqr3i, sqr3i, sqr3i,
               -sqr3i, sqr3i, sqr3i,
               -sqr3i, sqr3i, sqr3i,
               -sqr3i, -sqr3i, sqr3i,
               -sqr3i, -sqr3i, sqr3i,
               -sqr3i, -sqr3i, sqr3i,
               sqr3i, -sqr3i, sqr3i,
               sqr3i, -sqr3i, sqr3i,
               sqr3i, -sqr3i, sqr3i,
            };
            const double sVecs[ nslip * ecmech::ndim ] = {
               zero, sqr2i, -sqr2i,
               -sqr2i, zero, sqr2i,
               sqr2i, -sqr2i, zero,
               -sqr2i, zero, -sqr2i,
               zero, -sqr2i, sqr2i,
               sqr2i, sqr2i, zero,
               zero, -sqr2i, -sqr2i,
               sqr2i, zero, sqr2i,
               -sqr2i, sqr2i, zero,
               sqr2i, zero, -sqr2i,
               zero, sqr2i, sqr2i,
               -sqr2i, -sqr2i, zero,
            };

            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       mVecs, sVecs, this->nslip);
         };

         __ecmech_hdev__ ~SlipGeomFCC(){};

         __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; };
         __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; };

      private:
         double _P_ref_vec[ ecmech::ntvec * nslip ];
         double _Q_ref_vec[ ecmech::nwvec * nslip ];
   }; // SlipGeomFCC

   /**
    * NOTE : the coding below is a hack just to get things going
    *
    * it is not the best way of doing things, and modifications should be made with great care
    */
   class SlipGeomHCPaBRYcaY1 // EVP_HCP_a_BRY_ca_Y1 32 in Fortran coding
   {
      public:

         //    3  slip systems in basal <a> family
         //    3  slip systems in prismatic <a> family
         //    6  slip systems in pyramidal <a> family
         //   12  slip systems in pyramidal 1 <c+a> family
         static const int nslip = 3+3+6+12;

         // constructor and destructor
         __ecmech_hdev__
         SlipGeomHCPaBRYcaY1(double cOverA)
            : _cOverA(cOverA)
         {
            
            //  pyramidal 10-11 1-210 depends on c/a
            //
            double m_ya[ecmech::ndim], s_ya[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one,  one } ; // plane
               double ab[ecmech::nMiller] = { one, -two,  one,  zero } ; // direction
               //
               miller_to_orthog_sngl(an, ab, 
                                     m_ya, s_ya,
                                     _cOverA);
            }
            double m_ya_pp = sqrt(1.0-m_ya[2]*m_ya[2]);

            // pyramidal 10-11 -1-123 depends on c/a
            //
            double m_y1ca[ecmech::ndim], s_y1ca[ecmech::ndim];
            {
               double an[ecmech::nMiller] = { one, zero, -one,  one } ; // plane
               double ab[ecmech::nMiller] = { -one, -one,  two,  three } ; // direction
               //
               miller_to_orthog_sngl(an, ab, 
                                     m_y1ca, s_y1ca,
                                     _cOverA);
            }
            double m_y1ca_pp = sqrt(1.0-m_y1ca[2]*m_y1ca[2]);
            double s_y1ca_pp = sqrt(1.0-s_y1ca[2]*s_y1ca[2]);

            const double mVecs[ nslip * ecmech::ndim ] = {
               zero,  zero,  one,
               zero,  zero,  one,
               zero,  zero,  one,

               -halfsqr3,   onehalf,   zero,
               -halfsqr3,  -onehalf,   zero,
               zero,  -one,   zero,

               m_ya[0],   m_ya[1],   m_ya[2],
               m_ya[0],  -m_ya[1],  -m_ya[2],
               m_ya[0],   m_ya[1],  -m_ya[2],
              -m_ya[0],   m_ya[1],  -m_ya[2],
                  zero,   m_ya_pp,  -m_ya[2],
                  zero,  -m_ya_pp,  -m_ya[2],

               m_y1ca[0],   m_y1ca[1],   m_y1ca[2],
               m_y1ca[0],  -m_y1ca[1],  -m_y1ca[2],
               m_y1ca[0],   m_y1ca[1],  -m_y1ca[2],
                    zero,   m_y1ca_pp,  -m_y1ca[2],
              -m_y1ca[0],   m_y1ca[1],  -m_y1ca[2],
              -m_y1ca[0],  -m_y1ca[1],  -m_y1ca[2],
                    zero,  -m_y1ca_pp,  -m_y1ca[2],
                    zero,   m_y1ca_pp,   m_y1ca[2],
              -m_y1ca[0],   m_y1ca[1],   m_y1ca[2],
              -m_y1ca[0],  -m_y1ca[1],   m_y1ca[2],
               m_y1ca[0],  -m_y1ca[1],   m_y1ca[2],
                    zero,  -m_y1ca_pp,   m_y1ca[2]
            };
            const double sVecs[ nslip * ecmech::ndim ] = {
               onehalf,   halfsqr3,   zero,
               onehalf,  -halfsqr3,   zero,
               one,   zero,   zero,

               onehalf,   halfsqr3,   zero,
               onehalf,  -halfsqr3,   zero,
               one,   zero,   zero,
               
                s_ya[0],   s_ya[1],   zero,
                s_ya[0],  -s_ya[1],   zero,               
               -s_ya[0],  -s_ya[1],   zero,
               -s_ya[0],   s_ya[1],   zero,
               -one,  zero,  zero,
                one,  zero,  zero,

                s_y1ca[0],   s_y1ca[1],   s_y1ca[2], 
                s_y1ca[0],  -s_y1ca[1],  -s_y1ca[2],
               -s_y1ca_pp,        zero,  -s_y1ca[2],
                s_y1ca[0],   s_y1ca[1],  -s_y1ca[2],
               -s_y1ca[0],   s_y1ca[1],  -s_y1ca[2],
                s_y1ca_pp,        zero,  -s_y1ca[2],
               -s_y1ca[0],  -s_y1ca[1],  -s_y1ca[2],
               -s_y1ca[0],   s_y1ca[1],   s_y1ca[2],
                s_y1ca_pp,        zero,   s_y1ca[2],
               -s_y1ca[0],  -s_y1ca[1],   s_y1ca[2],
               -s_y1ca_pp,        zero,   s_y1ca[2],
                s_y1ca[0],  -s_y1ca[1],   s_y1ca[2]
            };

            fillFromMS(this->_P_ref_vec, this->_Q_ref_vec,
                       mVecs, sVecs, this->nslip);
         };

         __ecmech_hdev__ ~SlipGeomHCPaBRYcaY1(){};

         __ecmech_hdev__ inline const double* getP() const { return _P_ref_vec; };
         __ecmech_hdev__ inline const double* getQ() const { return _Q_ref_vec; };

      private:
         double _cOverA;
         double _P_ref_vec[ ecmech::ntvec * nslip ];
         double _Q_ref_vec[ ecmech::nwvec * nslip ];
   }; // SlipGeomHCPaBRYcaY1
   
} // namespace ecmech

#endif // ECMECH_SLIPGEOM_H
