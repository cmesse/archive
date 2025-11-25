//
// Created by Christian Messe on 15.04.21.
//
#include "cl_EN_ConfigChamber.hpp"
#include "constants.hpp"
namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        ConfigChamber::ConfigChamber( XML & aXML )
        {
            // select the correct tree
            aXML.select_subtree( "belfem/turbopump/chamber") ;

            // read chamber pressure
            mPc = aXML.get_real("combustionpressure") * 1e5 ;

            // read oxidizer fuel ratio
            mOF = aXML.get_real( "oxidizerfuelratio") ;

            // read throat diameter
            mDt = aXML.get_real( "throatdiameter" ) * 0.001 ;

            // compute cross section at throat
            mAt = 0.25 * constant::pi * mDt * mDt ;
        }

//------------------------------------------------------------------------------

    }
}
