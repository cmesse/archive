//
// Created by Christian Messe on 26.09.19.
//

#include "cl_CN_Reaction.hpp"
#include "assert.hpp"
#include "fn_sum.hpp"
#include "constants.hpp"
#include "GT_globals.hpp"
#include "cl_CN_Scheme.hpp"
#include "fn_dot.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Reaction::Reaction(
                Scheme               & aScheme,
                const Vector< uint > & aEductIndices,
                const Vector< real > & aEductNu,
                const Vector< uint > & aProductIndices,
                const Vector< real > & aProductNu,
               const          bool     aHasThirdBody ) :
            mScheme( aScheme ),
            mEductIndices( aEductIndices ),
            mEductNu( aEductNu ),
            mProductIndices( aProductIndices ),
            mProductNu( aProductNu ),
            mHaveThirdBody( aHasThirdBody ),
            mN1( aEductIndices.length() ),
            mN2( aProductIndices.length() ),
            mPhi1( -sum( aEductNu ) ),
            mPhi2( -sum( aProductNu ) ),
            mSumNu( sum( aEductNu )-sum( aProductNu ) )
        {
            BELFEM_ASSERT( aEductIndices.length() == aEductNu.length(),
                    "Lengths of indices and nu for educts do not match" );

            BELFEM_ASSERT( aProductIndices.length() == aProductNu.length(),
                          "Lengths of indices and nu for educts do not match" );

            mdPsi1dY.set_size( aScheme.combgas()->number_of_components(), 0.0 );
            mdPsi2dY.set_size( aScheme.combgas()->number_of_components(), 0.0 );

            // populate delta nu vector
            mDeltaNu.set_size( aScheme.combgas()->number_of_components(), 0.0 );

            for( uint k=0; k<mN1; ++k )
            {
                mDeltaNu( mEductIndices( k ) ) -= mEductNu( k );
            }
            for( uint k=0; k<mN2; ++k )
            {
                mDeltaNu( mProductIndices( k ) ) += mProductNu( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Reaction::set_third_body(
                const Vector< uint > & aIndices,
                const Vector< real > & aWeights )
        {
            BELFEM_ASSERT( mHaveThirdBody,
                    "Can not call Reaction::set_third_body if mHaveThirdBody==false." );

            BELFEM_ASSERT( aIndices.length() == aWeights.length(),
                          "Lengths of third body indices third body weights do not match" );

            mThirdBodyIndices = aIndices;
            mThirdBodyWeights = aWeights;
        }

//------------------------------------------------------------------------------

        void
        Reaction::eval(   const real     & aT,
                          const real     & aP,
                          Vector< real > & aS,
                          Matrix< real > & aJ )
        {
            // calculate alpha parameter
            mAlpha = mScheme.combgas()->alpha( aT, aP );

            // calculate third body concentration
            this->eval_cm();

            // psi parameters
            this->eval_Psi();
            this->eval_dPsidY();
            this->eval_dPsidT();

            this->eval_forward_reaction_speed( aT );
            this->eval_backward_reaction_speed( aT );


            // eval S-Term
            aS.vector_data() += mDeltaNu * ( mk1 * mPsi1 - mk2 * mPsi2 );

            uint tN = mScheme.number_of_reacting_species();


            for( uint j=0; j<tN; ++j )
            {
                for( uint i=0; i<tN; ++i )
                {
                    aJ( i, j ) += mDeltaNu( i )
                                    * ( mk1 * mdPsi1dY( j ) - mk2 * mdPsi2dY( j ) );
                }
            }

            // last column
            for( uint i=0; i<tN; ++i )
            {
                // fixme: for testing purpose only
                aJ( i, tN ) += mDeltaNu( i ) * (
                          mdk1dT * mPsi1
                        + mk1 * mdPsi1dT
                        - mdk2dT * mPsi2
                        - mk2 * mdPsi2dT );

                //aJ( i, tN ) += mDeltaNu( i ) * ( mdk1dT * mPsi1 - mdk2dT * mPsi2 );

            }
        }
//------------------------------------------------------------------------------

        void
        Reaction::eval_cm()
        {


            if( mHaveThirdBody )
            {
                mCm = 0.0;
                // Gerlinger ( 2.44 )
                for( uint k=0; k<mThirdBodyIndices.length(); ++k )
                {
                    mCm += mThirdBodyWeights( k ) * mScheme.c( mThirdBodyIndices( k ) );
                }

                mdCmdT = -mCm * mAlpha;
            }
            else
            {
                mCm = 1.0 ;
                mdCmdT = 0.0;
            }
        }

//------------------------------------------------------------------------------

        void
        Reaction::eval_Psi()
        {
            // forward reaction
            mPsi1 = mCm; // <<-- is 1 if there is no inert partner

            for( uint k=0; k<mN1; ++k )
            {
                mPsi1 *= std::pow( mScheme.c( mEductIndices( k ) ), mEductNu( k ) );
            }

            // backward reaction
            mPsi2 =  mCm; // <<-- is 1 if there is no inert partner
            for( uint k=0; k<mN2; ++k )
            {
                mPsi2 *= std::pow( mScheme.c( mProductIndices( k ) ), mProductNu( k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        Reaction::eval_dPsidY()
        {
            // forward reaction
            mdPsi1dY.fill( 0.0 );

            for( uint k=0; k<mN1; ++k )
            {
                uint j = mEductIndices( k );
                if( mScheme.Y( j ) > BELFEM_EPSILON )
                {
                    mdPsi1dY( j ) = mPsi1 * mEductNu( k ) / mScheme.Y( j );
                }
            }

            // backward reaction
            mdPsi2dY.fill( 0.0 );
            for( uint k=0; k<mN2; ++k )
            {
                uint j = mProductIndices( k );
                if( mScheme.Y( j ) > BELFEM_EPSILON )
                {
                    mdPsi2dY( j ) = mPsi2 * mProductNu( k ) / mScheme.Y( j );
                }
            }
        }


//------------------------------------------------------------------------------

        void
        Reaction::eval_dPsidT()
        {
            mdPsi1dT = mAlpha * mPhi1 * mPsi1;
            mdPsi2dT = mAlpha * mPhi2 * mPsi2;
        }

//------------------------------------------------------------------------------

        void
        Reaction::eval_backward_reaction_speed( const real & aT)
        {
            // ( 3.9 )
            // unit of A: mol / ccm
            real tA = std::pow( gastables::gPref *1e-6/ ( constant::Rm * aT ), mSumNu );

            //real tA = std::pow( gastables::gPref / ( constant::Rm * aT ), mSumNu );
            real tdAdT = -mSumNu * tA / aT;

            real tG = 0.0;
            real tdGdT = 0.0;
            for ( uint k = 0; k < mN1; ++k )
            {
                tG -= mEductNu( k ) * mScheme.G( mEductIndices( k ));
                tdGdT -= mEductNu( k ) * mScheme.dGdT( mEductIndices( k ));
            }
            for ( uint k = 0; k < mN2; ++k )
            {
                tG += mProductNu( k ) * mScheme.G( mProductIndices( k ));
                tdGdT += mProductNu( k ) * mScheme.dGdT( mProductIndices( k ));
            }

            // equilibrium constandt with reference to pressure, Gerlinger (2.59)
            // note: sign is correct
            real tB = std::exp( tG / ( constant::Rm * aT ));

            real tdBdT = tB * ( aT * tdGdT - tG ) / ( constant::Rm * aT * aT );

            mk2 = mk1 * tA * tB;
            mdk2dT = mdk1dT * tA * tB + mk1 * ( tdAdT * tB + tA * tdBdT );
        }

//------------------------------------------------------------------------------

    }

}