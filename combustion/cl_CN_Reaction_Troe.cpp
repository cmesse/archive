//
// Created by Christian Messe on 28.09.19.
//

#include  <cmath> // for std::isnan
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_CN_Reaction_Troe.hpp"
#include "fn_CN_arrhenius.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Reaction_Troe::Reaction_Troe(
                Scheme & aScheme,
                const Vector <uint> & aEductIndices,
                const Vector <real> & aEductNu,
                const Vector <uint> & aProductIndices,
                const Vector <real> & aProductNu,
                const bool aHasThirdBody,
                const Vector < real > & aCoeffsLow,
                const Vector < real > & aCoeffsHigh,
                const Vector < real > & aTroe ) :
                Reaction(
                        aScheme,
                        aEductIndices,
                        aEductNu,
                        aProductIndices,
                        aProductNu,
                        aHasThirdBody),
                mCoeffs1( aCoeffsLow ),
                mCoeffs2( aCoeffsHigh ),
                mA( aTroe( 0 ) ),
                mTau3( -1.0/aTroe( 1 ) ),
                mTau1( -1.0/aTroe( 2 ) ),
                mTau2( -aTroe( 3 ) ),
                mLog10( std::log( 10.0 ) )
        {
            // test which low pressure funciton to use
            if ( aCoeffsLow( 1 ) == 0.0 )
            {
                mFunctionK1    = & arrhenius_simple;
                mFunctiondK1dT = & darrhenius_simpledT;
            }
            else
            {
                mFunctionK1    = & arrhenius;
                mFunctiondK1dT = & darrheniusdT;
            }

            // test which high pressure funciton to use
            if ( aCoeffsHigh( 1 ) == 0.0 )
            {
                mFunctionK2    = & arrhenius_simple;
                mFunctiondK2dT = & darrhenius_simpledT;
            }
            else
            {
                mFunctionK2    = & arrhenius;
                mFunctiondK2dT = & darrheniusdT;
            }

            // test which Fcent to use
            if ( std::isnan( aTroe( 3 ) ) )
            {
                mFunctionCent  = & Reaction_Troe::F_cent_simple ;
            }
            else
            {
                mFunctionCent  = & Reaction_Troe::F_cent ;
            }
        }

//------------------------------------------------------------------------------

        void
        Reaction_Troe::eval_forward_reaction_speed( const real & aT )
        {
            // see also http://akrmys.com/public/chemkin/CKm_inp.html.en

            // LOW pressure part
            real tk1 = ( *mFunctionK1 ) ( mCoeffs1, aT );
            real tdk1dT = ( *mFunctiondK1dT ) ( mCoeffs1, aT, tk1 );

            // HIGH pressure part
            real tk2 = ( *mFunctionK2 ) ( mCoeffs2, aT );
            real tdk2dT = ( *mFunctiondK2dT ) ( mCoeffs2, aT, tk2 );

            // Helpers
            real tX    = tk1 * mCm / tk2;
            real tdXdT = ( tk1 * tk2 * mdCmdT + mCm * tk2 * tdk1dT - mCm * tk1 * tdk2dT ) / ( tk2 * tk2 );

            real tY = std::log10( tX );
            real tdYdT = tdXdT / ( tX * mLog10 );

            // Lindemann part of reaction
            real tL = tk2 * tX / ( 1.0 + tX );
            real tdLdT = ( tX * ( 1.0 + tX ) * tdk2dT + tk2 * tdXdT ) / ( ( 1.0 + tX ) * ( 1.0 + tX ) );

            // Troe Function Fc
            real tFc;
            real tdFcdT;
            ( this->*mFunctionCent )( aT , tFc, tdFcdT );

            // Log 10-ize function Fc
            real tGc    = std::log10( tFc );
            real tdGcdT = tdFcdT / ( tFc * mLog10 );

            // N-Parameter
            real tN = 0.75 - 1.27 * tGc;
            real tdNdT = -1.27 * tdGcdT;

            // c-Parameter
            real tC    = -0.4-0.67*tGc;
            real tdCdT = -0.67*tdGcdT;

            // d-Parameter
            real tD = 0.14;

            // Temporary Parameter
            real tH    = ( tY + tC ) / ( tN - tD*( tY + tC ) );
            real tdHdT = ( tN * ( tdCdT + tdYdT) - ( tC + tY )*tdNdT ) /
                    std::pow( tN - tD*(tC+tY), 2 );

            // Temporary Parameter
            real tI = ( 1.0 + tH * tH );
            real tJ    = tGc / tI ;
            real tdJdT = ( tI * tdGcdT - 2.0 * tGc * tH * tdHdT ) / ( tI * tI );

            // Troe Function
            real tF    = std::pow( 10, tJ );
            real tdFdT =   mLog10 * tF * tdJdT;

            // final result
            mk1   =  tL * tF ;
            mdk1dT = tdLdT * tF + tL * tdFdT;
        }

//------------------------------------------------------------------------------

        void
        Reaction_Troe::F_cent( const real & aT, real & aF, real & adFdT  )
        {
            real tExp3 = std::exp( mTau3 * aT );
            real tExp2 = std::exp( mTau2 / aT );
            real tExp1 = std::exp( mTau1 * aT );

            aF = ( 1.0 - mA ) * tExp3 +  mA * tExp1 + tExp2;

            adFdT =    ( 1.0 - mA ) * tExp3 * mTau3
                     +         mA   * tExp1 * mTau1
                     -                tExp2 * mTau2 / ( aT * aT );
        }

//------------------------------------------------------------------------------

        void
        Reaction_Troe::F_cent_simple( const real & aT, real & aF, real & adFdT  )
        {
            real tExp3 = std::exp( mTau3 * aT );
            real tExp1 = std::exp( mTau1 * aT );

            aF = ( 1.0 - mA ) * tExp3 +  mA * tExp1;

            adFdT =    ( 1.0 - mA ) * tExp3 * mTau3
                       +         mA * tExp1 * mTau1;
        }

//------------------------------------------------------------------------------
    } /* namespace combustion */
} /* namespace belfem */