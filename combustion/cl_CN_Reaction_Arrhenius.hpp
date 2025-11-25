//
// Created by Christian Messe on 26.09.19.
//

#ifndef BELFEM_CL_CN_REACTION_ARRHENIUS_HPP
#define BELFEM_CL_CN_REACTION_ARRHENIUS_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_CN_Reaction.hpp"
#include "fn_CN_arrhenius.hpp"


namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        class Reaction_Arrhenius : public Reaction
        {

            const Vector< real > mCoeffs;

            // pointer to k-function
            real ( *mFunctionK )
                ( const Vector< real > & aCoeffs, const real & aT );

            // pointer to dkdT-function
            real ( *mFunctiondKdT )
                ( const Vector< real > & aCoeffs, const real & aT, const real & aK );

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Reaction_Arrhenius(              Scheme & aScheme,
                                const Vector <uint> & aEductIndices,
                                const Vector <real> & aEductNu,
                                const Vector <uint> & aProductIndices,
                                const Vector <real> & aProductNu,
                                const bool aHasThirdBody,
                                const Vector < real > & aCoeffs );

//------------------------------------------------------------------------------

            void
            eval_forward_reaction_speed( const real & aT );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace combustion */
} /* namespace belfem */
#endif //BELFEM_CL_CN_REACTION_ARRHENIUS_HPP
