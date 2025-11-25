//
// Created by Christian Messe on 15.04.21.
//

#ifndef BELFEM_CL_EN_CONFIGSHAFT_HPP
#define BELFEM_CL_EN_CONFIGSHAFT_HPP

#include "cl_XML.hpp"
#include "typedefs.hpp"

namespace belfem
{
    namespace engine
    {
        class ConfigShaft
        {
            // rotation in RPM
            real mN = BELFEM_QUIET_NAN ;

            // mechanical shaft efficiency
            real mEta = BELFEM_QUIET_NAN ;
//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            ConfigShaft( XML & aXML );

//------------------------------------------------------------------------------

            ~ConfigShaft() = default ;

//------------------------------------------------------------------------------

            /**
             * returns the shaft speed in RPM
             */
            inline const real &
            n()
            {
                return mN ;
            }

//------------------------------------------------------------------------------

            /**
             * returns the shaft efficiency
             */
            inline const real &
            eta()
            {
                return mEta ;
            }

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_EN_CONFIGSHAFT_HPP
