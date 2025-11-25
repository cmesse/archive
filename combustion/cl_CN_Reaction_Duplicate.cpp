//
// Created by Christian Messe on 28.09.19.
//


#include "cl_CN_Reaction_Duplicate.hpp"
#include "assert.hpp"
#include "fn_CN_arrhenius.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Reaction_Duplicate::Reaction_Duplicate(
                Scheme & aScheme,
                const Vector <uint> & aEductIndices,
                const Vector <real> & aEductNu,
                const Vector <uint> & aProductIndices,
                const Vector <real> & aProductNu,
                const bool aHasThirdBody,
                const Vector < real > & aCoeffs1,
                const Vector < real > & aCoeffs2 ) :
                Reaction(
                        aScheme,
                        aEductIndices,
                        aEductNu,
                        aProductIndices,
                        aProductNu,
                        aHasThirdBody),
                mCoeffs1( aCoeffs1 ),
                mCoeffs2( aCoeffs2 )
        {
            // test which funciton we use
            if ( aCoeffs1( 1 ) == 0.0 )
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
            if ( aCoeffs2( 1 ) == 0.0 )
            {
                mFunctionK2    = & arrhenius_simple;
                mFunctiondK2dT = & darrhenius_simpledT;
            }
            else
            {
                mFunctionK2   = & arrhenius;
                mFunctiondK2dT = & darrheniusdT;
            }

        }

//------------------------------------------------------------------------------

        void
        Reaction_Duplicate::eval_forward_reaction_speed( const real & aT )
        {

            real tk1 = ( *mFunctionK1 ) ( mCoeffs1, aT );
            real tk2 = ( *mFunctionK2 ) ( mCoeffs2, aT );
            mk1    =  tk1 + tk2;
            mdk1dT =   ( *mFunctiondK1dT ) ( mCoeffs1, aT, tk1 )
                     + ( *mFunctiondK2dT ) ( mCoeffs2, aT, tk2 );
        }

//------------------------------------------------------------------------------

    }
}