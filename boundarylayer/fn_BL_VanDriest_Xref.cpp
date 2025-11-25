//
// Created by Christian Messe on 07.05.20.
//

#ifndef BELFEM_FN_BL_VANDRIEST_XREF_HPP
#define BELFEM_FN_BL_VANDRIEST_XREF_HPP

#include "typedefs.hpp"
#include "cl_BL_State.hpp"
#include "fn_BL_VanDriest_Xref.hpp"
#include "fn_BL_VanDriest.hpp"
namespace belfem
{
    namespace boundarylayer
    {
        real
        vandriest_xref( State & aState,
                     const real & aDotQ,
                     const real aKs,
                     const real aMangler )
        {
            real tX0 = 1e-9;
            real tX1 = 2000.0;

            vandriest( aState, tX0, aKs, aMangler );
            real tF0 = aState.dot_q() - aDotQ;

            if ( tF0 < 0 ) return tX0 ;

            vandriest( aState, tX1, aKs, aMangler );
            real tF = aState.dot_q() - aDotQ;

            real aX;
            uint tCount = 0;

            while ( std::abs( tF0 - tF ) > 0.0001 )
            {
                aX = 0.5 * ( tX0 + tX1 );
                vandriest( aState, aX, aKs, aMangler );
                tF = aState.dot_q() - aDotQ;

                if ( tF0 * tF > 0 )
                {
                    tX0 = aX;
                    tF0 = tF;
                }
                else
                {
                    tX1 = aX;
                }

                BELFEM_ERROR( tCount++ < 1000, "too many iterations " );
            }

            return aX;
        }
    }
}

#endif //BELFEM_FN_BL_VANDRIEST_XREF_HPP
