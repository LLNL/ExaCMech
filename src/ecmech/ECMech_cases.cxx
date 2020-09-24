#include "ECMech_cases.h"

namespace ecmech {

   /**
    * @brief These are not the only possible cases -- they are here as a convenience
    */
   __ecmech_host__
   matModelBase* makeMatModel( const std::string &modelName ) {
      
      matModelBase* matModel = nullptr ;
      //
      if ( modelName == "evptn_FCC_A" ) {
         ecmech::matModelEvptn_FCC_A* mmECMEvptn = new ecmech::matModelEvptn_FCC_A() ;
         matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn) ;
      }
      else if ( modelName == "evptn_FCC_AH" ) {
         ecmech::matModelEvptn_FCC_AH* mmECMEvptn = new ecmech::matModelEvptn_FCC_AH() ;
         matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn) ;
      }
      else if ( modelName == "evptn_FCC_B" ) {
         ecmech::matModelEvptn_FCC_B* mmECMEvptn = new ecmech::matModelEvptn_FCC_B() ;
         matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn) ;
      }
      else if ( modelName == "evptn_BCC_A" ) {
         ecmech::matModelEvptn_BCC_A* mmECMEvptn = new ecmech::matModelEvptn_BCC_A() ;
         matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn) ;
      }
      else if ( modelName == "evptn_HCP_A" ) {
         ecmech::matModelEvptn_HCP_A* mmECMEvptn = new ecmech::matModelEvptn_HCP_A() ;
         matModel = dynamic_cast<ecmech::matModelBase*>(mmECMEvptn) ;
      }
      else {
         std::string msg = std::string("model name not recognized : ") + modelName ;
         ECMECH_FAIL(__func__, msg.c_str()) ;
      }
      
      return matModel ;
   }
   
} // namespace ecmech
