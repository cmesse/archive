//
// Created by Christian Messe on 15.04.21.
//

#ifndef BELFEM_CL_EN_CONFIGCHAMBER_HPP
#define BELFEM_CL_EN_CONFIGCHAMBER_HPP

#include "cl_XML.hpp"
#include "typedefs.hpp"

namespace belfem
{
    namespace engine
    {
        class ConfigChamber
        {
            // combustion pressure
            real mPc = BELFEM_QUIET_NAN ;

            // oxidizer-fuel-ratio
            real mOF = BELFEM_QUIET_NAN ;

            // throat diameter
            ream mDt = BELFEM_QUIET_NAN ;

            // cross section at throat
            ream mAt = BELFEM_QUIET_NAN ;

        public:
//------------------------------------------------------------------------------

            ConfigChamber( XML & aXML );

//------------------------------------------------------------------------------

            ~ConfigChamber() = default ;

//------------------------------------------------------------------------------

            // return the combustion pressure
            inline const real &
            pc const()
            {
                return mPc;
            }

//------------------------------------------------------------------------------

            // return the oxidizer-fuel ratio
            inline const real &
            of const()
            {
                return mOF ;
            }

//------------------------------------------------------------------------------

            // return the throat diameter
            inline const real &
            Dt const()
            {
                return mDt ;
            }

//------------------------------------------------------------------------------

            // return the throat cross section
            inline const real &
            At const()
            {
                return mAt ;
            }

//------------------------------------------------------------------------------

        };
    }
}
#endif //BELFEM_CL_EN_CONFIGCHAMBER_HPP
