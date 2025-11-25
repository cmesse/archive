//
// Created by Christian Messe on 25.09.19.
//
#include "cl_CN_Scheme.hpp"
#include "cl_CN_Chemkin.hpp"

#include "cl_Map.hpp"
#include "fn_dot.hpp"
#include "constants.hpp"
#include "fn_sum.hpp"
#include "cl_CN_ReactionFactory.hpp"
#include "cl_CN_Entry.hpp"
#include "cl_GT_RefGas.hpp"
#include "fn_gesv.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Scheme::Scheme(
                const string & aChemkinFilePath,
                const Fuel aFuel,
                const Oxidizer aOxidizer,
                const GasModel aGasModel ) :
                mChemkinFileName( filename( aChemkinFilePath )),
                mFuel( aFuel ),
                mOxidizer( aOxidizer ),
                mGasModel( aGasModel )
        {
            if ( aChemkinFilePath.length() > 0 )
            {
                // read the chemkin file
                Chemkin tChemkinFile( aChemkinFilePath );

                // create the combustion gas
                this->create_combgas( tChemkinFile );

                // create the factory
                ReactionFactory tFactory( this );

                // create reactions
                mReactions.set_size( tChemkinFile.number_of_reactions(), nullptr );

                // counter
                uint tCount = 0;

                for( uint k=0; k<tChemkinFile.number_of_entries(); ++k )
                {
                    // get entry
                    Entry * tEntry = tChemkinFile.entry( k );

                    // test if this entry is active
                    if( tEntry->is_active() )
                    {
                        mReactions( tCount++ ) = tFactory.create_reaction( tEntry );
                    }
                }

                BELFEM_ASSERT( tCount == tChemkinFile.number_of_reactions(),
                        "number of reactions does not match" );
            }
            else
            {
                // An empty scheme with air as combgas.
                // Needed for testing only
                mCombgas = new Gas();

                mNumberOfReactingSpecies = mCombgas->number_of_components();
            }

            // Allocate containers
            mNumberOfAllSpecies = mCombgas->number_of_components();

            mY.set_size( mNumberOfAllSpecies );
            mY0.set_size( mNumberOfAllSpecies );

            mC.set_size( mNumberOfReactingSpecies );

            mH.set_size( mNumberOfAllSpecies, 0.0 );
            mCp.set_size( mNumberOfAllSpecies, 0.0 );
            mdCpdT.set_size( mNumberOfAllSpecies, 0.0 );
            mH0.set_size( mNumberOfAllSpecies, 0.0 );

            for( uint k=0; k<mNumberOfAllSpecies; ++k )
            {
                mH0( k )
                    = mCombgas->component( k )->data()->Hf() / mCombgas->component( k )->M()
                            - mCombgas->component( k )->h( gastables::gTref );
            }

            mGibbs.set_size( mNumberOfAllSpecies );
            mdGibbsdT.set_size( mNumberOfAllSpecies );

            mM.set_size( mNumberOfAllSpecies );

            for( uint k=0; k<mNumberOfAllSpecies; ++k )
            {
                mM( k ) = mCombgas->data( k )->M();
            }

            mdYdt.set_size( mNumberOfAllSpecies );
            mRHS.set_size( mNumberOfReactingSpecies + 1 );
            mLHS.set_size( mNumberOfReactingSpecies + 1, 0.0 );
            mPivot.set_size( mNumberOfReactingSpecies + 1 );

            mJacobi.set_size( mNumberOfReactingSpecies + 1, mNumberOfReactingSpecies + 1 );

            mTemperatureIndex = mNumberOfReactingSpecies ;

        }

//------------------------------------------------------------------------------

        Scheme::~Scheme()
        {
            delete mCombgas;

            for( Reaction * tReaction: mReactions )
            {
                delete tReaction;
            }
        }

//------------------------------------------------------------------------------

        void
        Scheme::create_combgas( Chemkin & aChemkinFile )
        {
            Cell<string> tSpecies;
            aChemkinFile.get_species( tSpecies );

            mNumberOfReactingSpecies = tSpecies.size();

            this->add_inert_fuel( tSpecies );

            if ( mOxidizer == Oxidizer::AIR )
            {
                this->add_air( tSpecies, mInitialMolarFractions );
            }
            else
            {
                this->add_oxidizer( tSpecies, mInitialMolarFractions );
            }

            // create the combustion gas
            mCombgas = new Gas( tSpecies, mInitialMolarFractions, mGasModel );

            // find entry for oxidzer
            this->find_oxidizer_index() ;

        }

//------------------------------------------------------------------------------

        void
        Scheme::add_inert_fuel( Cell<string> & aSpecies )
        {
            // get string of fuel
            string tFuelLabel = fuel_to_string( mFuel );

            mInertFuelIndex = aSpecies.size();

            mReactingFuelIndex = aSpecies.size();

            // find index in species
            for ( uint k = 0; k < mInertFuelIndex; ++k )
            {
                if ( aSpecies( k ) == tFuelLabel )
                {
                    mReactingFuelIndex = k;
                    break;
                }
            }

            // make sure that fuel exists in reaction scheme
            BELFEM_ERROR( mReactingFuelIndex < mInertFuelIndex,
                         "could not find fuel %s in reaction scheme %s.",
                         tFuelLabel.c_str(),
                         mChemkinFileName.c_str());

            // add intert fuel to list
            aSpecies.push( tFuelLabel );
        }

//------------------------------------------------------------------------------

        void
        Scheme::add_air(
                Cell<string> & aSpecies,
                Vector<real> & aMolarFractions )
        {

            // create a temporary map
            Map<string, uint> tMap;

            // in case of airbreathing and air, k<aSpecies.size()-1 would put
            // CH4 of air into the active position, passive if k<aSpecies.size()
            // Former one is intended.
            for ( uint k = 0; k < aSpecies.size() - 1; ++k )
            {
                tMap[ aSpecies( k ) ] = k;
            }

            // create a tempoarary air gas
            Gas tAir;

            // counter for indices
            uint tCount = aSpecies.size();

            Vector<uint> tIndex( tAir.number_of_components());

            for ( uint k = 0; k < tAir.number_of_components(); ++k )
            {
                // get label of gas
                const string & tSpecie = tAir.data( k )->label();

                // test if species already exists
                if ( tMap.key_exists( tSpecie ))
                {
                    // only set index of specie
                    tIndex( k ) = tMap( tSpecie );
                }
                else
                {
                    // add species to list
                    aSpecies.push( tSpecie );

                    // add entry to map
                    tMap[ tSpecie ] = tCount;

                    // add counter to index
                    tIndex( k ) = tCount;

                    ++tCount;
                }
            }

            aMolarFractions.set_size( tCount, 0.0 );

            // populate molar fractions
            const Vector<real> & tMolarFractions = tAir.molar_fractions();

            for ( uint k = 0; k < tAir.number_of_components(); ++k )
            {
                aMolarFractions( tIndex( k )) = tMolarFractions( k );
            }
        }

//------------------------------------------------------------------------------

        void
        Scheme::add_oxidizer( Cell<string> & aSpecies,
                              Vector<real> & aMolarFractions )
        {
            // find oxidizer flag in gas
            uint tOxidizerFlag = aSpecies.size();
            string tOxidizerLabel = oxidizer_to_string( mOxidizer );

            for ( uint k = 0; k < aSpecies.size(); ++k )
            {
                if ( aSpecies( k ) == tOxidizerLabel )
                {
                    tOxidizerFlag = k;
                    break;
                }
            }

            // make sure that fuel exists in reaction scheme
            BELFEM_ERROR( tOxidizerFlag < aSpecies.size(),
                         "could not find oxidizer %s in reaction scheme %s.",
                         tOxidizerLabel.c_str(),
                         mChemkinFileName.c_str());

            // prepare molar fractions
            aMolarFractions.set_size( aSpecies.size(), 0.0 );
            aMolarFractions( tOxidizerFlag ) = 1.0;
        }

//------------------------------------------------------------------------------

        void
        Scheme::reset_combgas_mixture()
        {
            mCombgas->remix( mInitialMolarFractions );
            mCount = 0 ;
            mLHS.fill( 0.0 );
        }

//------------------------------------------------------------------------------

        void
        Scheme::preprocess( const real & aT, const real & aP )
        {
            // we take the mass fractions from the combustion gas ...
            mY = mCombgas->mass_fractions() ;

            // specific enthalpy and specific heat capacity
            for ( uint k = 0; k < mNumberOfReactingSpecies; ++k )
            {
                // compute concentration ( part 1 )
                mC( k ) = mY( k ) / mM( k );
            }

            // specific volume in ccm / kg
            mV = mCombgas->v( aT, aP ) * 1e6;

            // concentration ( part 2 )
            mC /= mV;

            // we only need specific enthalpies for reacting species
            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
                mH( k ) = mCombgas->h( k, aT, aP ) + mH0( k ) ;
            }

            for ( uint k = 0; k < mNumberOfAllSpecies; ++k )
            {
                mCp( k ) = mCombgas->cp( k, aT, aP );
                mdCpdT( k ) = mCombgas->dcpdT( k, aT, aP );
            }

            // calculate gibbs free energy
            mCombgas->Gibbs( aT, mGibbs );

            // calculate temperature derivative of gibbs free energy
            mCombgas->dGibbsdT( aT, mdGibbsdT );
        }

//------------------------------------------------------------------------------

       real
       Scheme::compute( const real & aT, const real & aP, const real & aU, const real & aDeltaX )
       {
            // sum of reacting species
            this->preprocess( aT, aP );

            this->compute_jacobi( aT, aP );

            this->compute_rhs( aT, aP, aU, aDeltaX );

            /*mJacobi.print("J");
            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
               std::cout << k << " " << mCombgas->component( k )->label() << " " << mdYdt( k ) << std::endl ;
            }
            std::cout << mNumberOfReactingSpecies << " " << "dTdt " << mdTdt << std::endl ;

            exit( 0 ); */

            mJacobi *= - mC1 ;
            for( uint k=0; k<=mNumberOfReactingSpecies; ++k )
            {
                mJacobi( k, k ) += 1. ;
            }

            mLHS = mRHS ;
            gesv( mJacobi, mLHS, mPivot );


            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
                // limit reducition
                if( mY( k ) + mLHS( k ) < 0 )
                {
                    mLHS( k ) = -mY( k );
                }
                mY( k ) += mLHS( k );
            }

            //mCombgas->remix_mass( mY , false );

            // temperature rise
            //std::cout << "heat " << mLHS( mTemperatureIndex ) << " " << mC1 << " " << mY( mReactingFuelIndex ) << " " << mY( mInertFuelIndex ) << std::endl ;
            //std::cout << "-----" << std::endl ;

            return mLHS( mTemperatureIndex ) ;
        }

//------------------------------------------------------------------------------

        void
        Scheme::compute_jacobi( const real & aT, const real & aP )
        {
            // reset matrices
            mdYdt.fill( 0.0 );

            mJacobi.fill( 0.0 );

            for ( Reaction * tReaction : mReactions )
            {
                tReaction->eval( aT, aP, mdYdt, mJacobi );
            }
            mdYdt %= mM;
            mdYdt *= mV;

            //Vector< real > & tS = mdYdt;

            for( uint j=0; j<mNumberOfReactingSpecies; ++j )
            {
                for( uint i=0; i<mNumberOfReactingSpecies; ++i )
                {
                    mJacobi( i, j ) *= mM( i );
                    //mJacobi( i, j ) += tS( i );
                }
            }

            // last column
            //real tAlpha = mCombgas->alpha( aT, aP );
            real tCp = dot( mCp, mY );
            real tdCpdT =  dot( mdCpdT, mY );
            real tdHdt = dot( mH, mdYdt );

            for( uint i=0; i<mNumberOfReactingSpecies; ++i )
            {
                mJacobi( i, mTemperatureIndex ) *= mM( i );
                //mJacobi( i, mTemperatureIndex ) += tAlpha * tS( i );
            }

            mJacobi *= mV;

            // help value for temperature derivative
            real tValue ;

            // final row
            for( uint j=0; j<mNumberOfReactingSpecies; ++j )
            {
                tValue = 0.0;

                for( uint i=0; i<mNumberOfReactingSpecies; ++i )
                {
                    tValue += mH( i ) * mJacobi( i, j );
                }

                mJacobi( mTemperatureIndex, j ) =
                        ( tdHdt * tdCpdT / tCp - tValue ) / tCp ;
            }

            // last entry
            mJacobi( mTemperatureIndex, mTemperatureIndex )
                = tdCpdT / ( tCp * tCp ) * tdHdt - dot( mCp, mdYdt ) / tCp;
        }

//------------------------------------------------------------------------------

        void
        Scheme::compute_rhs( const real & aT, const real & aP, const real & aU, const real & aDeltaX )
        {

            //if( mCount++ < 1 )
            //{
            //    mC1 = aDeltaX / aU ;
            //    mC2 = 0.0 ;
            //}
            //else
            //{*/
                mC1 = 2./3. * aDeltaX / aU ;
                mC2 = 1.0/3.0 ;
            //}

            // contribution from molar fractions
            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
                mRHS( k ) = mC1 * mdYdt( k ) + mC2 * mLHS( k );
            }

            // temperature change
            mdTdt = 0.0 ;
            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
                mdTdt += mH( k ) * mdYdt( k );
            }
            mdTdt /= - dot( mCp, mY );

            mRHS( mTemperatureIndex )
                = mC1 * mdTdt + mC2 * mLHS( mTemperatureIndex );
        }

//------------------------------------------------------------------------------

        void
        Scheme::find_oxidizer_index()
        {
            // get oxidizer label
            string tLabel ;

            switch( mOxidizer )
            {
                case( Oxidizer::AIR ) :
                case( Oxidizer::LOX ) :
                {
                    tLabel = "O2" ;
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Unknown Oxidizer");
                }
            }

            // get cell with reference gases
            Cell< gastables::RefGas * > & tComponents = mCombgas->components() ;

            mOxidizerIndex = tComponents.size() ;

            // initialize a counter
            uint tCount = 0 ;

            // find oxygen in combustion gas ( we can't use a map
            // because fuel is used twice here )
            for( gastables::RefGas * tComponent : tComponents )
            {
                if( tComponent->label() == "O2" )
                {
                    // copy the index
                    mOxidizerIndex = tCount ;

                    // break the loop
                    break ;
                }

                // increment counter
                ++tCount ;
            }

            BELFEM_ERROR( mOxidizerIndex < tComponents.size(),
                         "Could not find oxygen entry in gas mixture" );
        }

//------------------------------------------------------------------------------

        real
        Scheme::delta_w()
        {
            real aDeltaW = 0.0 ;
            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
                aDeltaW += ( mY( k ) - mY0( k ) ) * mH( k );
            }
            return aDeltaW ;
        }

//------------------------------------------------------------------------------

        real
        Scheme::delta_R()
        {
            real aDeltaR = 0.0 ;
            for( uint k=0; k<mNumberOfReactingSpecies; ++k )
            {
                aDeltaR += ( mY( k ) - mY0( k ) ) / mM( k );
            }
            return aDeltaR ;
        }

//------------------------------------------------------------------------------



    }
}