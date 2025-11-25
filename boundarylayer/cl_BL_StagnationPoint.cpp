//
// Created by Christian Messe on 14.06.20.
//

#include "cl_BL_StagnationPoint.hpp"

namespace belfem
{
    namespace boundarylayer
    {
//------------------------------------------------------------------------------

        StagnationPoint::StagnationPoint(
                Gas   & aGas,
                State & aFreesteam,
                State & aStagnation,
        const real & aRadius ) :
        mGas( aGas ),
        mFreesteam( aFreesteam ),
        mStagnation( aStagnation ),
        mNoseRadius( aRadius )
        {

        }

//------------------------------------------------------------------------------

        void
        StagnationPoint::compute_flowstates(
                const real & aT,
                const real & aP,
                const real & aU )
        {
            mFreesteam.compute( aT, aP, aU );

            // state after perpendicular shock
            real tT1 ;
            real tP1 ;
            real tU1 ;
            mGas.shock( aT, aP, aU, tT1, tP1, tU1 );

            // state at stagnation
            real tTs ;
            real tPs ;
            mGas.total( tT1, tP1, tU1, tTs, tPs );

            mStagnation.compute( tTs, tPs, 0.0 );

            // set the pressure coefficient for the stagnation point
            mStagnation.set_Cp( 2.0 * ( tPs - aP ) / ( mFreesteam.rho() * aU * aU ) );

        }

//------------------------------------------------------------------------------

        void
        StagnationPoint::compute_wallstate( const real & aTw )
        {
            mStagnation.set_wall_temperature( aTw );
            mStagnation.compute_wall_state() ;
        }

//------------------------------------------------------------------------------

        real
        StagnationPoint::compute_stagnation_heatload( const real & aTw  )
        {
            this->compute_wallstate( aTw );

            // velocity derivative at stagnation point according to newton, Eq. 64
            real tDudXs = std::sqrt(
                    2.0 * ( mStagnation.p() - mFreesteam.p() ) /
                    mStagnation.rho() ) / mNoseRadius ;


            // see hirschel, L is reciprocal to definition in Wikipedia
            mLewis = 2.0 * mStagnation.Pr_w();

            // from equation 62
            mF = 0.67 * std::pow( mStagnation.rho() * mStagnation.mu() /
                                  ( mStagnation.rho_w() * mStagnation.mu_w() ), 0.4 );

            // Equation 52
            mG =  1.0 + ( std::pow( mLewis, mPhi ) - 1.0 ) *
                        mGas.hd( mStagnation.T(), mStagnation.p() ) / mStagnation.h();

            // Equation 63
            real aDotQ =   1.1343 / std::pow( mStagnation.Pr_w(), 0.6 )  * mF * mG
                           * ( mStagnation.h() - mStagnation.hw() ) * std::sqrt( tDudXs );

            std::cout << "stag " << mLewis << " " << mF << " " << mG << " " << tDudXs << " " << aDotQ << std::endl ;

            // exit( 0 );

            return aDotQ ;
        }

//------------------------------------------------------------------------------

        real
        StagnationPoint::compute_stagnation_heatload( const real & aTw, const real & aX )
        {
            this->compute_wallstate( aTw );

            // see hirschel, L is reciprocal to definition in Wikipedia
            mLewis = 2.0 * mStagnation.Pr_w();

            // from equation 62
            mF = 0.67 * std::pow( mStagnation.rho() * mStagnation.mu() /
                                  ( mStagnation.rho_w() * mStagnation.mu_w() ), 0.4 );

            // Equation 52
            mG =  1.0 + ( std::pow( mLewis, mPhi ) - 1.0 ) *
                        mGas.hd( mStagnation.T(), mStagnation.p() ) / mStagnation.h();

            // Equation 44
            real tRe = mStagnation.rho_w() * mFreesteam.u() * aX / mStagnation.mu_w() ;

            // Equation 62
            real tNu = mF * mG * std::sqrt( tRe );

            real tSt = tNu / ( mStagnation.Pr_w() * tRe );

            real aDotQ = mStagnation.rho_w() * mFreesteam.u() * tSt * ( mStagnation.h() - mStagnation.hw() );

            std::cout << mLewis << " " << mF << " " << mG << " " << tRe << " " << aDotQ << std::endl ;
            exit( 0 );

            return aDotQ ;
        }

//------------------------------------------------------------------------------
    }
}