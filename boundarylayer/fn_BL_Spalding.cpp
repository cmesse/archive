//
// Created by Christian Messe on 2019-03-30.
//

#include "fn_BL_Spalding.hpp"
#include "assert.hpp"
namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        real
        spalding(
                const real & aB,
                const real & aKappa,
                const real & aExp,
                const real & aYplus,
                const real   aF_guess )
        {
            BELFEM_ASSERT( aYplus>=0.0, "Invalid value for Y+ = %7.0f", aYplus );
            real aF;

            // initial guess for f
            if ( aF_guess == 0.0 )
            {
                if ( aYplus < 10 )
                {
                    aF = aYplus;
                }
                else
                {
                    aF = std::log( aYplus ) / aKappa + aB;
                }
             }
             else
            {
                 aF = aF_guess;
            }

            real tF = BELFEM_REAL_MAX;

            real tYplus;
            real tdYplus;

            uint tCount = 0 ;
            const real tOmega = 0.99 ;

            // start Newton iteration
            while( std::abs( tF - aF ) > 20 * 1e-12 )
            {
                // shift F
                tF = aF;

                // Spalding wall function
                tYplus = spalding_y( aB, aKappa, aExp, aF  );

                // Derivative of Spalding Wall Funciton
                tdYplus = spalding_dydf( aB, aKappa, aExp, aF );

                // correct f
                aF -= tOmega * ( tYplus - aYplus ) / tdYplus;

                BELFEM_ERROR( tCount++ < 100 , "too many iterations");
            }

            // return value
            return aF;
        }

//------------------------------------------------------------------------------

        real
        spalding_y(
                const real & aB,
                const real & aKappa,
                const real & aExp,
                const real & aF )
        {
            real tkF = aKappa * aF;

            real tG = std::exp( tkF );

            real  tH = 1.0 + tkF *
                   (  24.0 + tkF *
                   (  12.0 + tkF *
                   (   4.0 + tkF ) ) ) / 24.0;


            return aF + aExp * ( tG - tH );
        }

//------------------------------------------------------------------------------

        real
        spalding_dydf(
                const real & aB,
                const real & aKappa,
                const real & aExp,
                const real & aF )
        {
            real tkF = aKappa * aF;
            real tdG = aKappa * std::exp( tkF );


            real tdH =
                    aKappa *
                    ( 6.0 + tkF *
                    ( 6.0 + tkF *
                    ( 3.0 + tkF ))) / 6.0;

            return 1.0 + aExp * ( tdG - tdH );
        }

//------------------------------------------------------------------------------
    }
}