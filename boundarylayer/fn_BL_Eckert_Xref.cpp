//
// Created by Christian Messe on 06.05.20.
//

#include "fn_BL_Eckert_Xref.hpp"
#include "fn_BL_Eckert.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        // eckert's method, writes heat loads into state
        // and returns the reference temperature
        real
        eckert_xref( State & aState,
                     const real & aDotQ,
                     bool aIsTurbulent,
                     const real aMangler )
        {
            real tX0 = 1e-9;
            real tX1 = 2000.0;

            eckert( aState, tX0, aIsTurbulent, aMangler );
            real tF0 = aState.dot_q() - aDotQ;

            if ( tF0 < 0 ) return tX0 ;

            eckert( aState, tX1, aIsTurbulent, aMangler );
            real tF = aState.dot_q() - aDotQ;

            real aX;
            uint tCount = 0;

            while ( std::abs( tF0 - tF ) > 0.0001 )
            {
                aX = 0.5 * ( tX0 + tX1 );
                eckert( aState, aX, aIsTurbulent, aMangler );
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
