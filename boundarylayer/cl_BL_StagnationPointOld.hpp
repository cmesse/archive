//
// Created by Christian Messe on 11.05.20.
//

#ifndef BELFEM_CL_BL_STAGNATIONPOINTOLD_HPP
#define BELFEM_CL_BL_STAGNATIONPOINTOLD_HPP

#include "typedefs.hpp"
#include "cl_Node.hpp"
#include "cl_Vector.hpp"
#include "cl_Gas.hpp"
#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {
        class StagnationPoint
        {
            Gas   & mGas ;

            State & mFreesteam ;

            State & mStagnation ;

            const real mRadius ;

            // equilibrium: 0.52, frozen: 0.63
            real mPhi = 0.52 ;
            real mLewis  ;

            // Equation ( 62 )
            real mF ;

            // Equation ( 62 ) divided by ( 58 )
            real mG ;

            // stagnation point heatflux
            real mDotQs;

            // stagnation point friction
            real mTauws;

            // Stanton number at stagnaiton point
            real mSts;

            // conditions after perpendicular shock
            real mT1 ;
            real mP1 ;
            real mU1 ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            StagnationPoint(
                    Gas        & aGas,
                    State      & aFreesteam,
                    State      & aStagnation,
                    const real & aRadius );

//------------------------------------------------------------------------------

            ~StagnationPoint() = default ;

//------------------------------------------------------------------------------

            void
            compute( const real & aTwnose );


//------------------------------------------------------------------------------

            real
            compute_heatload( State & aState );

//------------------------------------------------------------------------------

            void
            compute_heatload( State & aState, const real & aX, const uint aCount=0 );

//------------------------------------------------------------------------------

            // returns the nose radius
            const real &
            nose_radius() const;

            State &
            freestream();

            State &
            stagnation();

//------------------------------------------------------------------------------

            // expose gas object
            Gas &
            gas() ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_stagnation_conditions();

//------------------------------------------------------------------------------

            void
            compute_stagnation_heatload();

//------------------------------------------------------------------------------

        };

//------------------------------------------------------------------------------

        // returns the nose radius
        inline const real &
        StagnationPoint::nose_radius() const
        {
            return mRadius;
        }

//------------------------------------------------------------------------------

        inline State &
        StagnationPoint::freestream()
        {
            return mFreesteam;
        }

//------------------------------------------------------------------------------

        inline State &
        StagnationPoint::stagnation()
        {
            return mStagnation;

        }

//------------------------------------------------------------------------------

        inline Gas &
        StagnationPoint::gas()
        {
            return mGas ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_BL_STAGNATIONPOINTOLD_HPP
