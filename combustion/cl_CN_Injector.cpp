//
// Created by Christian Messe on 23.04.20.
//

#include "cl_CN_Injector.hpp"
#include "cl_GT_RefGas.hpp"
#include "assert.hpp"
#include "cl_Vector.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Injector::Injector(
                Scheme & aScheme,
                const real aXinj,
                const real aDuctHeight,
                const real aEtaMix,
                const real aPulsonetti ) :
                mScheme( aScheme ),
                mGas( *aScheme.combgas()),
                mXinj( aXinj ),
                mDuctHeight( aDuctHeight ),
                mEtaMix( aEtaMix ),
                mPulsonetti( aPulsonetti ),
                mOFst( this->compute_stoich_fuel_oxidizer_ratio( aScheme ) ),
                mLnGamma( std::log( mGamma ) )
        {
            mY.set_size( aScheme.combgas()->number_of_components() );
            mX.set_size( aScheme.combgas()->number_of_components() );

            //this->compute_mixing_length();


        }

//------------------------------------------------------------------------------

        void
        Injector::print()
        {
            std::cout << "Duct Height " << mDuctHeight << std::endl;
            std::cout << "Mixing Length " << mMixingLength << std::endl;
            std::cout << "OF st " << mOFst << std::endl ;

            // get mass fractions
            std::cout << "reacting fuel" << mGas.mass_fraction( mScheme.reacting_fuel_index() ) << std::endl ;
            std::cout << "inert fuel" << mGas.mass_fraction( mScheme.inert_fuel_index() ) << std::endl ;
        }


//------------------------------------------------------------------------------

        void
        Injector::compute_mixing_length()
        {
            // AIAA 88-3258 Eq. ( 16 )
            // and 10.18419/opus-9381 Eq. ( 3.18 )
            mMixingLength = 0.179 * mDuctHeight * mPulsonetti
                            * std::exp( 1.72 * mPhi );

        }

//------------------------------------------------------------------------------

        real
        Injector::compute_stoich_fuel_oxidizer_ratio( Scheme & aScheme )
        {
            // get components
            Cell< gastables::RefGas * > & tComponents = aScheme.combgas()->components() ;

            // get oxidizer
            gastables::RefGas * tOxidizer = tComponents( aScheme.oxidizer_index() );

            // get fuel
            gastables::RefGas * tFuel = tComponents( aScheme.reacting_fuel_index() );

            // make sure that oxidizer contains only of C, H2 and O2
            const Cell< string > & tOx = tOxidizer->data()->elements() ;
            uint tN = tOx.size() ;
            for( uint k=0; k<tN; ++k )
            {
                string tString = tOx( k );
                BELFEM_ERROR( tString == "C" || tString == "H" || tString == "O",
                             "Unsupportet component in Oxidizer: %s", tString.c_str() );
            }

            // make sure that fuel contains only of C, H2 and O2
            const Cell< string > & tFl = tFuel->data()->elements() ;
            tN = tFl.size() ;
            for( uint k=0; k<tN; ++k )
            {
                string tString = tFl( k );
                BELFEM_ERROR( tString == "C" || tString == "H" || tString == "O",
                             "Unsupportet component in Fuel: %s", tString.c_str() );
            }

            // let:
            // CxHyOz + alpha * CaHbOc = beta * CO2 + gamma * H2O

            // fuel multiplcities
            real tX = tFuel->data()->component_multiplicity("C");
            real tY = tFuel->data()->component_multiplicity("H");
            real tZ = tFuel->data()->component_multiplicity("O");

            // oxidizer multiplicities
            real tA = tOxidizer->data()->component_multiplicity("C");
            real tB = tOxidizer->data()->component_multiplicity("H");
            real tC = tOxidizer->data()->component_multiplicity("O");

            // make sure we don't divide by zero
            BELFEM_ERROR( std::abs( tC ) > BELFEM_EPSILON,
                    "Oxidizer %s does not seem to contain any oxygen",
                         tOxidizer->label().c_str() );

            real tDet = 4.0 * tA + tB - 2.0 * tC ;

            BELFEM_ERROR( std::abs( tDet ) > BELFEM_EPSILON,
                    "Error in fuel oxidizer composition" );

            // Oxidizer multiplicity
            real tAlpha = ( 2.0 * tZ - tY - 4.0 * tX ) / tDet ;

            // carbon dioxide multiplicity
            // real tBeta  = ( tB * tX - tA * tY + 2.0 * ( tA * tZ - tC * tX  ) )/ tDet ;

            // water multiplicity
            // real tGamma = ( 2.0 * ( tA * tY - tB * tX ) + tB * tZ - tC * tY ) / tDet ;

            // compute stoiciometic ratio
            return tAlpha * tOxidizer->M() /
                    ( aScheme.combgas()->mass_fraction( aScheme.oxidizer_index()  )
                    * tFuel->M() );
        }

//------------------------------------------------------------------------------

        // set stoichiometric ratio
        void
        Injector::set_phi( const real & aPhi )
        {
            BELFEM_ERROR( aPhi > BELFEM_EPSILON, "stiochiometric ratio must be > 0" );

            mPhi = aPhi ;

            mOF = mOFst / aPhi ;

            if( mAutoMixingLength )
            {
                this->compute_mixing_length() ;
            }
        }

//------------------------------------------------------------------------------

        // set oxidizer to fuel ratio
        void
        Injector::set_of( const real & aOF )
        {
            BELFEM_ERROR( aOF > BELFEM_EPSILON, "oxidizer to fuel ratio must be > 0" );

            mOF = aOF ;
            mPhi = mOFst / aOF ;

            if( mAutoMixingLength )
            {
                this->compute_mixing_length() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Injector::set_oxidizer_massflow( const real & aOxidizerMassflow )
        {
            mOxidizerMassflow = aOxidizerMassflow ;
            mFuelMassflow     = aOxidizerMassflow / mOF ;
        }

//------------------------------------------------------------------------------

        void
        Injector::set_massflow( const real & aMassflow )
        {
            mFuelMassflow     = aMassflow / ( 1.0 + mOF );
            mOxidizerMassflow = mFuelMassflow * mOF ;
        }

//------------------------------------------------------------------------------

        void
        Injector::set_mixing_length( const real & aMixingLength )
        {
            mAutoMixingLength = false ;

            mMixingLength = aMixingLength ;
        }

//------------------------------------------------------------------------------

        real
        Injector::mix( const real & aX ) const
        {
            // derivati10.18419/opus-9381 Eq. ( 3.17 )
            return mEtaMix * ( 1.0 - std::pow( mGamma,
               ( aX - mXinj ) / mMixingLength ) ) ;
        }

//------------------------------------------------------------------------------
        // derivative of mixing function
        real
        Injector::dmix( const real & aX ) const
        {
            // 10.18419/opus-9381 Eq. ( 3.17 )
            return - mEtaMix * mLnGamma * std::pow( mGamma,
                    ( aX - mXinj ) / mMixingLength ) / mMixingLength ;
        }

//------------------------------------------------------------------------------

        // move fuel from inert to reacting
        void
        Injector::inject( const real & aX )
        {
            real tMu = mFuelMassflow /
                       ( mOxidizerMassflow + mFuelMassflow ) ;

            // get mass fractions from gas
            mY = mGas.mass_fractions() ;

            // molar fraction of inert fuel
            real & tY_inert = mY( mScheme.inert_fuel_index()  );

            // molar fraction of reacting fuel
            real & tY_react = mY( mScheme.reacting_fuel_index() );

            real tY_Expect = ( 1.0 - this->mix( aX ) ) * tMu ;

            // change of mass fraction
            real tDeltaY = tY_inert - tY_Expect ;

            //if( tY_inert > tDeltaY )
            //{
                tY_inert -= tDeltaY;
                tY_react += tDeltaY;
            //}

            BELFEM_ASSERT( tY_inert > 0, "negative mass fraction while injecting fuel" );

            // remix, but leave fuel properties alone
            mGas.remix_mass( mY, false );
        }

//------------------------------------------------------------------------------

    }
}