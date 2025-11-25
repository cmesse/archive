//
// Created by Christian Messe on 26.09.19.
//

#include "cl_CN_Reaction_Arrhenius.hpp"
#include "assert.hpp"
#include "fn_CN_arrhenius.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Reaction_Arrhenius::Reaction_Arrhenius(
                             Scheme & aScheme,
                const Vector <uint> & aEductIndices,
                const Vector <real> & aEductNu,
                const Vector <uint> & aProductIndices,
                const Vector <real> & aProductNu,
                const bool aHasThirdBody,
                const Vector < real > & aCoeffs ) :
                Reaction(
                        aScheme,
                        aEductIndices,
                        aEductNu,
                        aProductIndices,
                        aProductNu,
                        aHasThirdBody),
                mCoeffs( aCoeffs )
        {
            // test which funciton we use
            if ( aCoeffs( 1 ) == 0.0 )
            {
                mFunctionK    = & arrhenius_simple;
                mFunctiondKdT = & darrhenius_simpledT;
            }
            else
            {
                mFunctionK    = & arrhenius;
                mFunctiondKdT = & darrheniusdT;
            }
        }

//------------------------------------------------------------------------------

        void
        Reaction_Arrhenius::eval_forward_reaction_speed( const real & aT )
        {
            mk1    = ( *mFunctionK ) ( mCoeffs, aT );
            mdk1dT =  ( *mFunctiondKdT ) ( mCoeffs, aT, mk1 );
        }

//------------------------------------------------------------------------------

    }
}