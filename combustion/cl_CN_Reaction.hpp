//
// Created by Christian Messe on 23.09.19.
//

#ifndef BELFEM_CL_CN_REACTION_HPP
#define BELFEM_CL_CN_REACTION_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_Gas.hpp"

namespace belfem
{
    namespace combustion
    {
        class Scheme;

        class Reaction
        {
//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------
            Scheme & mScheme;

            const Vector <uint> mEductIndices;
            const Vector <real> mEductNu;
            const Vector <uint> mProductIndices;
            const Vector <real> mProductNu;
            const bool mHaveThirdBody;

            const uint mN1;
            const uint mN2;

            const real mPhi1;
            const real mPhi2;
            const real mSumNu;

            Vector <uint> mThirdBodyIndices;
            Vector <real> mThirdBodyWeights;

            real mAlpha;

            // third body concentration in ccm / kg
            real mCm;
            real mdCmdT;

            real mPsi1;
            real mPsi2;
            Vector< real > mDeltaNu;

            Vector <real> mdPsi1dY;
            Vector <real> mdPsi2dY;
            real mdPsi1dT;
            real mdPsi2dT;

            real mk1;
            real mdk1dT;
            real mk2;
            real mdk2dT;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Reaction(
                    Scheme & aScheme,
                    const Vector <uint> & aEductIndices,
                    const Vector <real> & aEductNu,
                    const Vector <uint> & aProductIndices,
                    const Vector <real> & aProductNu,
                    const bool aHasThirdBody );

//------------------------------------------------------------------------------

            inline bool
            has_third_body() const;

//------------------------------------------------------------------------------

            inline bool
            has_third_body_weights() const;

//------------------------------------------------------------------------------

            virtual ~Reaction() = default;

//------------------------------------------------------------------------------

            void
            set_third_body(
                    const Vector <uint> & aIndices,
                    const Vector <real> & aWeights );
//------------------------------------------------------------------------------

            void
            eval(   const real & aT,
                    const real & aP,
                    Vector< real > & aS,
                    Matrix< real > & aJ );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            eval_cm();

            void
            eval_Psi();

            void
            eval_dPsidY();

            void
            eval_dPsidT();

            virtual void
            eval_forward_reaction_speed( const real & aT ) = 0;

            void
            eval_backward_reaction_speed( const real & aT );

//------------------------------------------------------------------------------
            //protected :
//------------------------------------------------------------------------------
        };

        bool
        Reaction::has_third_body() const
        {
            return mHaveThirdBody;
        }

//------------------------------------------------------------------------------

        bool
        Reaction::has_third_body_weights() const
        {
            return mThirdBodyIndices.length() > 0;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CN_REACTION_HPP
