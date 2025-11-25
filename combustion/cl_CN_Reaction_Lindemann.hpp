//
// Created by Christian Messe on 28.09.19.
//

#ifndef BELFEM_CL_CN_REACTION_LINDEMANN_HPP
#define BELFEM_CL_CN_REACTION_LINDEMANN_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_CN_Reaction.hpp"
#include "fn_CN_arrhenius.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        class Reaction_Lindemann : public Reaction
        {

            const Vector< real > mCoeffs1;
            const Vector< real > mCoeffs2;

            // pointer to k-function
            real ( *mFunctionK1 )
                    ( const Vector< real > & aCoeffs, const real & aT );

            // pointer to dkdT-function
            real ( *mFunctiondK1dT )
                    ( const Vector< real > & aCoeffs, const real & aT, const real & aK );

            // pointer to k2-function
            real ( *mFunctionK2 )
                    ( const Vector< real > & aCoeffs, const real & aT );

            // pointer to dkdT-function
            real ( *mFunctiondK2dT )
                    ( const Vector< real > & aCoeffs, const real & aT, const real & aK );

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

            Reaction_Lindemann(          Scheme              & aScheme,
                                         const Vector <uint> & aEductIndices,
                                         const Vector <real> & aEductNu,
                                         const Vector <uint> & aProductIndices,
                                         const Vector <real> & aProductNu,
                                         const bool            aHasThirdBody,
                                         const Vector < real > & aCoeffsLow,
                                         const Vector < real > & aCoeffsHigh );

//------------------------------------------------------------------------------

            void
            eval_forward_reaction_speed( const real & aT );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace combustion */
} /* namespace belfem */
#endif //BELFEM_CL_CN_REACTION_LINDEMANN_HPP
