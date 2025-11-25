//
// Created by Christian Messe on 14.06.20.
//

#ifndef BELFEM_CL_BL_STAGNATIONPOINT_HPP
#define BELFEM_CL_BL_STAGNATIONPOINT_HPP

#include "typedefs.hpp"
#include "cl_Gas.hpp"
#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {

        class StagnationPoint
        {
            Gas & mGas;

            State & mFreesteam;

            State & mStagnation;

            const real mNoseRadius;

            // equilibrium: 0.52, frozen: 0.63
            real mPhi = 0.52 ;

            // Lewis number
            real mLewis ;


            real mF ;
            real mG ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            StagnationPoint(
                    Gas & aGas,
                    State & aFreesteam,
                    State & aStagnation,
                    const real & aRadius );

//------------------------------------------------------------------------------

            ~StagnationPoint() = default;

//------------------------------------------------------------------------------

            void
            compute_flowstates( const real & aT, const real & aP, const real & aU );

//------------------------------------------------------------------------------

            // returns the nose radius
            const real &
            nose_radius() const;

//------------------------------------------------------------------------------

            /**
             * expost the freestream state
             * @return
             */
            State &
            freestream();

//------------------------------------------------------------------------------

            /**
              * expose the stagnation state
              * @return
              */
            State &
            stagnation();

//------------------------------------------------------------------------------

            // expose gas object
            Gas &
            gas() ;

//------------------------------------------------------------------------------

            real
            compute_stagnation_heatload( const real & aTw  ) ;

//------------------------------------------------------------------------------

            real
            compute_stagnation_heatload( const real & aTw, const real & aX ) ;

//------------------------------------------------------------------------------

        private:

//------------------------------------------------------------------------------

            void
            compute_wallstate( const real & aTw );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        // returns the nose radius
        inline const real &
        StagnationPoint::nose_radius() const
        {
            return mNoseRadius;
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


    }
}

#endif //BELFEM_CL_BL_STAGNATIONPOINT_HPP
