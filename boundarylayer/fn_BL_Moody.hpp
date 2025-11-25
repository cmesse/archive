//
// Created by Christian Messe on 04.01.21.
//

#ifndef BELFEM_FN_BL_MOODY_HPP
#define BELFEM_FN_BL_MOODY_HPP

#include "typedefs.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        /**
         * the moody chart function, see VDI Heat Atlas
         * Chapter Lab 2
         * fully developed turbulent flow is assumed
         *
         * @param aReDh : Reynods Number with respect to hydraulic diameter
         * @param aDh   : Hydraulic Diameter in m
         * @param aK    : absolute roughness in m
         * @return
         */
        real
        cf_moody( const real & aReDh, const real & aDh, const real aK=0.0 );

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_FN_BL_MOODY_HPP
