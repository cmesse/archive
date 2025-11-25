//
// Created by Christian Messe on 15.04.20.
//

#ifndef BELFEM_FN_BL_VANDRIEST_HPP
#define BELFEM_FN_BL_VANDRIEST_HPP

#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        void
        vandriest(
                State & aState,
                const real & aX,         // x-coordinate for reynolds number
                const real aKs = 0.0,    // Roughness Factor in mm
                const real aMangler=1.0  // parameter for cone transformation
                        );
    }
}

#endif //BELFEM_FN_BL_VANDRIEST_HPP
