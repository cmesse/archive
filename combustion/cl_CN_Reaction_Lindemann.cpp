//
// Created by Christian Messe on 28.09.19.
//

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_CN_Reaction_Lindemann.hpp"
#include "fn_CN_arrhenius.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Reaction_Lindemann::Reaction_Lindemann(
                Scheme & aScheme,
                const Vector <uint> & aEductIndices,
                const Vector <real> & aEductNu,
                const Vector <uint> & aProductIndices,
                const Vector <real> & aProductNu,
                const bool aHasThirdBody,
                const Vector < real > & aCoeffsLow,
                const Vector < real > & aCoeffsHigh ) :
                Reaction(
                        aScheme,
                        aEductIndices,
                        aEductNu,
                        aProductIndices,
                        aProductNu,
                        aHasThirdBody),
                mCoeffs1( aCoeffsLow ),
                mCoeffs2( aCoeffsHigh )
        {
            // test which funciton we use
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

            // test which funciton we use
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
        }

//------------------------------------------------------------------------------

        void
        Reaction_Lindemann::eval_forward_reaction_speed( const real & aT )
        {
            // LOW pressure part
            real tk1 = ( *mFunctionK1 ) ( mCoeffs1, aT );
            real tdk1dT = ( *mFunctiondK1dT ) ( mCoeffs1, aT, tk1 );

            // HIGH pressure part
            real tk2 = ( *mFunctionK2 ) ( mCoeffs2, aT );
            real tdk2dT = ( *mFunctiondK2dT ) ( mCoeffs2, aT, tk2 );

            // Helpers
            real tX    = tk1 * mCm / tk2;
            real tdXdT = ( tk1 * tk2 * mdCmdT + mCm * tk2 * tdk1dT - mCm * tk1 * tdk2dT ) / ( tk2 * tk2 );

            // Lindemann reaction
            mk1 = tk2 * tX / ( 1.0 + tX );
            mdk1dT = ( tX * ( 1.0 + tX ) * tdk2dT + tk2 * tdXdT ) / ( ( 1.0 + tX ) * ( 1.0 + tX ) );
        }

//------------------------------------------------------------------------------
    } /* namespace combustion */
} /* namespace belfem */