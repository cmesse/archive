//
// Created by Christian Messe on 04.01.21.
//

#include "fn_BL_Moody.hpp"
#include "assert.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        real
        cf_moody( const real & aReDh, const real & aDh, const real aK )
        {
            // initial guess using dittus-boelther
            real tCf = 0.046 / std::pow( aReDh, 0.2 );

            // help magnitude
            real tX = 1.0 / std::sqrt( 4.0 * tCf ) ;

            // residuum
            real tF = BELFEM_REAL_MAX ;

            // counter
            uint tCount = 0 ;

            // relaxaiton
            real tOmega = 0.9 ;

            // other help magnitudes
            real tA = 2.51 / aReDh ;
            real tB = aK / ( aDh * 3.71 );
            real tC = 2.0 / std::log( 10.0 );
            real tD ;
            real tdF ;

            while( std::abs( tF ) > 1e-9 )
            {
                tD = tA * tX + tB ;

                BELFEM_ERROR( std::abs( tD ) > BELFEM_EPSILON, "Algorithm fail." );

                tF  = tX + tC * std::log( tD );
                tdF = 1.0 + tA * tC / ( tD );

                tX -= tOmega * tF / tdF ;

                BELFEM_ERROR( tCount++ < 100, "Too many iterations" );
            }

            // reconstruct friction factor
            return 0.25 / ( tX * tX );
        }

//------------------------------------------------------------------------------
    }
}