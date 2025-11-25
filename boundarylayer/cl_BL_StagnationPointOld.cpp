//
// Created by Christian Messe on 11.05.20.
//

#include "cl_BL_StagnationPointOld.hpp"
#include "fn_BL_VanDriest.hpp"

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
                mRadius( aRadius )
        {

        }

//------------------------------------------------------------------------------

        void
        StagnationPoint::compute( const real & aTwnose  )
        {
            mStagnation.set_wall_temperature( aTwnose );

            this->compute_stagnation_conditions();

            // see hirschel
            mLewis = 2.0 * mStagnation.Pr_w();

            // from equation 62
            mF = 0.67 * std::pow( mStagnation.rho() * mStagnation.mu() /
                         ( mStagnation.rho_w() * mStagnation.mu_w() ), 0.4 );

            // Equation 52
            mG =  1.0 + ( std::pow( mLewis, mPhi ) - 1.0 ) *
                        mGas.hd( mStagnation.T(), mStagnation.p() ) / mStagnation.h();

            this->compute_stagnation_heatload();

        }

//------------------------------------------------------------------------------

        void
        StagnationPoint::compute_stagnation_conditions()
        {
            // perpendicular shock
            mGas.shock( mFreesteam.T(), mFreesteam.p(), mFreesteam.u(), mT1, mP1, mU1 );

            // total conditions
            real tTs ;
            real tPs ;
            mGas.total( mT1, mP1, mU1, tTs, tPs );
            mStagnation.compute( tTs, tPs, 0.0 );

            std::cout << "Ts " << tTs << std::endl;
            std::cout << "ps " << tPs << std::endl;

            mStagnation.set_Cp( 2.0 * ( tPs - mFreesteam.p() )  / ( mFreesteam.rho() * mFreesteam.u()*  mFreesteam.u() ) );

            mStagnation.compute_wall_state() ;


        }

//------------------------------------------------------------------------------

        void
        StagnationPoint::compute_stagnation_heatload()
        {
            // Eq. 64
            real tdudx = std::sqrt( 2.0 * ( mStagnation.p() - mFreesteam.p() ) / mStagnation.rho() ) / mRadius ;

            // Eq. 63
            mDotQs =  0.76 / std::pow( mStagnation.Pr_w(), 0.6 )
                    * std::pow( mStagnation.rho_w() * mStagnation.mu_w(), 0.1 )
                    * std::pow( mStagnation.mu() * mStagnation.rho(), 0.4 )
                    * mG * ( mStagnation.h() - mStagnation.hw() ) * std::sqrt( tdudx ) ;

            mSts = mDotQs / ( mFreesteam.rho() * mFreesteam.u() * ( mStagnation.h() - mStagnation.hw()) ) ;

            real tSigma = std::pow( mStagnation.Pr_w(), 2.0/3.0 );

            mTauws = mDotQs * tSigma * mU1 / ( mStagnation.h() - mStagnation.hw() );

            mStagnation.set_loads( mTauws, mDotQs, mStagnation.h() );

            std::cout << "stag " << std::endl;
            std::cout << mGas.hd( mStagnation.T(), mStagnation.p() )  / mStagnation.h() << " " << mG << " " << mDotQs << std::endl;
        }

//------------------------------------------------------------------------------

        real
        StagnationPoint::compute_heatload( State & aState )
        {

            real aDotQ = mSts * mFreesteam.rho() * mFreesteam.u()
                    * ( mStagnation.h() - aState.hw() );

            aState.set_loads( mTauws, aDotQ, mStagnation.h() );

            return aDotQ ;
        }

//------------------------------------------------------------------------------

        void
        StagnationPoint::compute_heatload( State & aState, const real & aX, const uint aCount )
        {
            // see hirschel
            real tL = 2.0 * aState.Pr_w();

            // from equation 62
            real tF = 0.67 * std::pow( mStagnation.rho() * mStagnation.mu() /
                                  ( aState.rho_w() * aState.mu_w() ), 0.4 );

            // Equation 52
            real tG =  1.0 + ( std::pow( tL, mPhi ) - 1.0 ) *
                        mGas.hd( mStagnation.T(), mStagnation.p() ) / mStagnation.h();

            // Reynolds number with respect to wall state, equation 44
            real tRe = aState.rho_w() * aState.u() * aX / aState.mu_w();

            // compute the Stanton number with respect to the wall state
            real tSt = tF * tG / ( std::sqrt( tRe ) * aState.Pr_w() ) ;

            // heat load
            real tDotQ = tSt * aState.rho_w() * aState.u() * ( mStagnation.h() - aState.hw() );

            // friction factor
            real tCf = 2.0 * tSt * std::pow( aState.Pr_w(), 0.6 );

            // shear stress
            real tTauw = 0.5 * aState.rho_w() * aState.u() * aState.u() * tCf ;

            // set values
            aState.set_loads( tTauw, tDotQ, mStagnation.h() );

            std::cout << aCount << " " << aX << " " << tDotQ << " " << tG << std::endl;
        }

//------------------------------------------------------------------------------
    }
}