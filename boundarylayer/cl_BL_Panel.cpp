//
// Created by Christian Messe on 05.05.20.
//

#include "cl_BL_Panel.hpp"
#include "cl_BL_StreamLineOld.hpp"
#include "fn_dot.hpp"
namespace belfem
{
    namespace boundarylayer
    {
//----------------------------------------------------------------------------

        Panel::Panel(  Gas                  & aGas,
                       State                & aFreestream,
                       State                & aStagnation,
                       mesh::Node           * aNode,
                       const real           & aSurfaceCoordinate,
                       const Vector< real > & aFlowDirection,
                       const Vector< real > & aNormalDirection ) :
            mGas( aGas ),
            mFreestream( aFreestream ),
            mStagnation( aStagnation ),
            mState( aGas ),
            mNode( *aNode ),
            mS( aSurfaceCoordinate ),
            mR( aFlowDirection ),
            mN( aNormalDirection )
        {
            // set a default temperature for the surface
            mState.set_wall_temperature( 800.0 );
        }

//----------------------------------------------------------------------------

        void
        Panel::compute_aoa( const Vector<real> &aFreestreamDirection )
        {
            mAoA = -std::asin( dot( aFreestreamDirection, mN  ) );
        }

//----------------------------------------------------------------------------

        void
        Panel::compute_newton()
        {
            mState.set_Cp(
                    mStagnation.Cp()
                    * std::pow( std::sin( mAoA ), 2.0 ) );

            // compute the pressure from the pressure coefficient
            real tP = 0.5 * mState.Cp()
                      * mFreestream.rho()
                      * mFreestream.u()
                      * mFreestream.u()
                      + mFreestream.p();
            // this must be an isotropic expansion along the surface
            real tT = mGas.isen_T(
                    mStagnation.T(),
                    mStagnation.p(),
                    tP );

            // compute enthalpy
            real tH = mGas.h( tT, tP );
            real tU = std::sqrt( 2.0 * ( mStagnation.h() - tH ) );

            // finish this state
            mState.compute( tT, tP, tU );
        }

//----------------------------------------------------------------------------

        void
        Panel::compute_prandtl_meyer( const real & aT, const real & aP, const real & aU, const real & aNu )
        {
            real tT;
            real tP;
            real tU;
            mGas.prandtl_meyer( aT, aP, aU, aNu, tT, tP, tU );

            mState.compute( tT, tP, tU );
        }

//----------------------------------------------------------------------------

    }
}