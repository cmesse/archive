//
// Created by Christian Messe on 15.04.20.
//
#include "assert.hpp"
#include "fn_BL_cf_flatplate_inc_turbulent.hpp"
namespace belfem
{
//------------------------------------------------------------------------------

    namespace boundarylayer
    {
        real
        cf_flatplate_inc_turbulent(
                const real & aReX,
                const real aBplus,
                const real aKarman,
                const real aPiWake,
                const real aFe )
        {
            // relaxation factor
            const real tOmega0 = 0.9;

            // Help values
            const real tLogval = std::max( std::log10( aReX ), 1.0 );
            const real tPsi = ( 2.0 * aPiWake - std::log( aFe )) / aKarman + aBplus;
            const real tSq2 = std::sqrt( 2.0 );

            // initial values for B+ = 0.0 and B+ = 5.0
            real tCf0 = 0.8446 * std::pow( tLogval, -2.793 );
            real tCf5 = 0.2929 * std::pow( tLogval, -2.421 );

            // interpolate initial guess
            real aCf = std::max( tCf0 + 0.2 * ( tCf5 - tCf0 ) * aBplus, tCf5 );

            // Backup values if Newton Fails
            tCf0 = 0.00001;
            real tCf1 = 5.0 * aCf;

            // function to be tested
            real tF = 1.0;

            // derivative of function
            real tdF;

            // help value
            real tSqCf;

            // initiaize counter
            uint tCount = 0;

            real tOmega1 ;
            real tOmega ;

            // begin newton
            while ( std::abs( tF ) > 1.0E-9 && tCount++ < 200 && aCf > 0.0 )
            {
                // get square root of cf
                tSqCf = std::sqrt( aCf );

                // compute function
                tF = std::log( 0.5 * aCf * aReX ) / aKarman + tPsi - tSq2 / tSqCf;

                // compute derivative
                tdF = 1.0 / ( aCf * tSqCf * tSq2 ) - 1.0 / ( aKarman * aCf );

                tOmega1 = std::abs( 0.9 * ( aCf - 0.0000001 ) * tdF/tF );
                tOmega = tOmega1 < tOmega0 ? tOmega1 : tOmega0 ;

                // perform newton step
                aCf -= tOmega * tF / tdF;
            }

            // check if we succeeded
            if ( !( std::abs( tF ) < 1E-9 ) || !( std::abs( aCf ) < 0.01 ) )
            {
                // the newton step failed. We try a bisection

                // initial step
                real tF0 = std::log( 0.5 * tCf0 * aReX ) / aKarman + tPsi - tSq2 / std::sqrt( tCf0 );

                // reset counter
                tCount = 0;

                while ( std::abs( tCf0 - tCf1 ) > 1.0E-9 && tCount++ < 1000 )
                {
                    aCf = 0.5 * ( tCf0 + tCf1 );

                    tF = std::log( 0.5 * aCf * aReX ) / aKarman + tPsi - tSq2 / std::sqrt( aCf );

                    if ( tF0 * tF > 0 )
                    {
                        tCf0 = aCf;
                        tF0 = tF;
                    }
                    else
                    {
                        tCf1 = aCf;
                    }
                }

                // check for sanity of result
                BELFEM_ERROR( tCount < 1000 && abs( aCf ) > 0.0,
                             "Function cf_flatplate_inc_turbulent failed for aRex= %E, aBplus = %E, aKarman=%E, aPiWake = %E, aFe= %E",
                             aReX, aBplus, aKarman, aPiWake, aFe );

            }

            // return result
            return aCf;
        }
//------------------------------------------------------------------------------
    }

//------------------------------------------------------------------------------
}