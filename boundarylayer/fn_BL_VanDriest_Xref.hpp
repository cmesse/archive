//
// Created by Christian Messe on 07.05.20.
//

#ifndef BELFEM_FN_BL_VANDRIEST_XREF_HPP
#define BELFEM_FN_BL_VANDRIEST_XREF_HPP

#include "typedefs.hpp"
#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        // eckert's method, writes heat loads into state
        // and returns the reference temperature
        real
        vandriest_xref( State & aState,
                        const real & aDotQ,
                        const real aKs = 0.0,
                        const real aMangler = 1.0 );

    }
}

#endif //BELFEM_FN_BL_VANDRIEST_XREF_HPP
