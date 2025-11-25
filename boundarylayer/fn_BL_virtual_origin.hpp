//
// Created by Christian Messe on 06.05.20.
//

#include "typedefs.hpp"
#include "cl_BL_State.hpp"

#ifndef BELFEM_FN_BL_VIRTUAL_ORIGIN_HPP
#define BELFEM_FN_BL_VIRTUAL_ORIGIN_HPP

namespace belfem
{
    namespace boundarylayer
    {
        /**
         * see Hirschel Chapter 10.4.4
         * @param aState
         * @param x
         * @param aMode
         *      1: from laminar to turbulent on flat plate / cylinder
         *      2: from laminar to turbulent on cone
         *      3: from cone to cylinder, laminar
         *      4: from cone to cylinder, turbulent
         * @return
         */
        real
        virtual_origin(
                State & aState,
                const real & aX,
                const uint aMode );

    }
}
#endif //BELFEM_FN_BL_VIRTUAL_ORIGIN_HPP
