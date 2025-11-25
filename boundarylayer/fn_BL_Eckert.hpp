//
// Created by Christian Messe on 15.04.20.
//

#ifndef BELFEM_FN_BL_ECKERT_HPP
#define BELFEM_FN_BL_ECKERT_HPP

#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        // eckert's method, writes heat loads into state
        // and returns the reference temperature
        real
        eckert( State & aState,
                const real & aX,
                bool aIsTurbulent,
                const real aMangler = 1.0 );

    }
}
#endif //BELFEM_FN_BL_ECKERT_HPP
