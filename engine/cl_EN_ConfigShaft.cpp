//
// Created by Christian Messe on 15.04.21.
//

#include "cl_EN_ConfigShaft.hpp"
#include "constants.hpp"
namespace belfem
{
    namespace engine
    {
        ConfigShaft::ConfigShaft( XML & aXML )
        {
            // select the correct tree
            aXML.select_subtree( "belfem/turbopump/shaft") ;

            // read the rotation
            mN = aXML.get_real("rotation" );

            // read the efficiency
            mEta = aXML.get_real( "eta" );
        }
    }
}