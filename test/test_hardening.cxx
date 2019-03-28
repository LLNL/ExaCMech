#include "SNLS_TrDLDenseG.h"

#include "ECMech_kinetics.h"
#include "ECMech_cases.h"

int main(int , // argc,
         char ** // *argv[]
         )
{
   // some convenience stuff
   using namespace ecmech ;

   const int nslip = 12 ;
   real8 dt = 1e-1 ;
   real8 gdot[nslip] = { 1.0/nslip } ;
   
   {
      KineticsVocePL kinetics(nslip) ;
#include "setup_kin_VocePL.h"
   
      std::vector<real8>       init ;
      {
         std::vector<std::string> names ;
         std::vector<bool>        plot ;
         std::vector<bool>        state ;
         kinetics.getHistInfo( names, init, plot, state ) ;
      }
      real8 hs_n[kinetics.nH] ;
      kinetics.updateH( hs_n, &(init[0]), dt, gdot ) ;
   
      std::cout << "Updated hardness state : " ;
      printVec<kinetics.nH>(hs_n, std::cout) ;
   }
   
   {
      Kin_KMBalD_FFF kinetics(nslip) ;
#include "setup_kin_KMBalD_FFF.h"
   
      std::vector<real8>       init ;
      {
         std::vector<std::string> names ;
         std::vector<bool>        plot ;
         std::vector<bool>        state ;
         kinetics.getHistInfo( names, init, plot, state ) ;
      }
      real8 hs_n[kinetics.nH] ;
      kinetics.updateH( hs_n, &(init[0]), dt, gdot ) ;
   
      std::cout << "Updated hardness state : " ;
      printVec<kinetics.nH>(hs_n, std::cout) ;
   }
   
}
