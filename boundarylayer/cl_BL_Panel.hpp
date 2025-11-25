//
// Created by Christian Messe on 05.05.20.
//

#ifndef BELFEM_CL_BL_PANEL_HPP
#define BELFEM_CL_BL_PANEL_HPP

#include "typedefs.hpp"
#include "cl_Node.hpp"
#include "cl_Vector.hpp"
#include "cl_Gas.hpp"
#include "cl_BL_State.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//----------------------------------------------------------------------------

        class Panel
        {

            // gas object
            Gas   & mGas ;

            State & mFreestream ;
            State & mStagnation ;

            // state object
            State mState ;

            // corresponding node on mesh
            mesh::Node & mNode ;

            // surface coordinate
            const real mS ;

            // direction vector
            const Vector< real > mR ;

            // normal vector
            const Vector< real > mN ;

            // angle of attack for this panel
            real mAoA = BELFEM_QUIET_NAN ;

            // coordinate for X-vector
            real mX = BELFEM_QUIET_NAN ;

//----------------------------------------------------------------------------
        public:
//----------------------------------------------------------------------------

            Panel(  Gas                  & aGas,
                    State                & aFreestream,
                    State                & aStagnation,
                    mesh::Node           * aNode,
                    const real           & aSurfaceCoordinate,
                    const Vector< real > & aFlowDirection,
                    const Vector< real > & aNormalDirection );

//----------------------------------------------------------------------------

            ~Panel() = default ;

//----------------------------------------------------------------------------

            // expose the node
            mesh::Node *
            node() ;

//----------------------------------------------------------------------------

            // expose the state
            State *
            state();

//----------------------------------------------------------------------------

            void
            flag() ;

//----------------------------------------------------------------------------

            void
            unflag() ;

//----------------------------------------------------------------------------

            bool
            is_flagged() const ;

//----------------------------------------------------------------------------

            // expose the tangent vector
            const Vector< real > &
            flow_direction();

//----------------------------------------------------------------------------

            void
            compute_aoa( const Vector< real > & aFreestreamDirection );

//----------------------------------------------------------------------------

            void
            compute_newton() ;

//----------------------------------------------------------------------------

            /**
             * the surface coordinate of this panel
             * @return
             */
            const real &
            s() const ;

//----------------------------------------------------------------------------

            /**
             * the angle of attack of this panel
             * @return
             */
            const real &
            aoa() const ;

//----------------------------------------------------------------------------

            void
            set_x( const real & aX ) ;

//----------------------------------------------------------------------------

            const real &
            x() const ;

//----------------------------------------------------------------------------

            void
            compute_prandtl_meyer( const real & aT, const real & aP, const real & aU, const real & aNu );

//----------------------------------------------------------------------------
        };
//----------------------------------------------------------------------------

        // expose the node
        inline mesh::Node *
        Panel::node()
        {
            return & mNode ;
        }

//----------------------------------------------------------------------------

        // expose the state
        inline State *
        Panel::state()
        {
            return & mState ;
        }

//----------------------------------------------------------------------------

        inline void
        Panel::flag()
        {
            mNode.flag();
        }

//----------------------------------------------------------------------------

        inline void
        Panel::unflag()
        {
            mNode.unflag();
        }

//----------------------------------------------------------------------------

        inline bool
        Panel::is_flagged() const
        {
            return mNode.is_flagged();
        }

//----------------------------------------------------------------------------

        // expose the tangent vector
        inline const Vector< real > &
        Panel::flow_direction()
        {
            return mR ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Panel::s() const
        {
            return mS ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Panel::aoa() const
        {
            return mAoA;
        }

//------------------------------------------------------------------------------

        inline void
        Panel::set_x( const real & aX )
        {
            mX = aX ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Panel::x() const
        {
            return mX ;
        }

//----------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_BL_PANEL_HPP
