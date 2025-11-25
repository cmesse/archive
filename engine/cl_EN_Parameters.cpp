//
// Created by Christian Messe on 31.08.20.
//

#include "cl_EN_Parameters.hpp"
#include "assert.hpp"
#include "cl_GT_RefGas.hpp"
#include "fn_sum.hpp"
namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        Parameters::Parameters()
        {
            // set temperatures and species based on default data
            this->set_fuel_and_oxidizer( mFuel, mOxidizer );
        }

//------------------------------------------------------------------------------

        // special constructor if meshinfo file is given
        Parameters::Parameters( HDF5 & aDatabase )
        {
            // set temperatures and species based on default data
            this->set_fuel_and_oxidizer( mFuel, mOxidizer );

            // Chamber data
            aDatabase.select_group("Chamber");
            real tThroatDiameter ;
            aDatabase.load_data("ThroatDiameter", tThroatDiameter );
            this->set_throat_diameter( tThroatDiameter );

            aDatabase.close_active_group() ;

            // Nozzle data
            aDatabase.select_group("Nozzle");
            real tExpansionRatio;
            aDatabase.load_data( "ExpansionRatio", tExpansionRatio );
            this->set_expansion_ratio( tExpansionRatio );
            aDatabase.close_active_group() ;

        }

//------------------------------------------------------------------------------

        void
        Parameters::set_thrust( const real & aF )
        {
            mDesignThrust = aF ;
        }

//------------------------------------------------------------------------------
        void
        Parameters::set_chamber_pressure( const real & aP )
        {
            mChamberPressure = aP ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_exit_pressure( const real & aP )
        {
            mExitPressure = aP ;
            mExpansionRatio = BELFEM_QUIET_NAN ;
            mNozzleMode   = NozzleMode::ComputeCrossSection ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_expansion_ratio( const real & aRatio )
        {
            mExitPressure   = BELFEM_QUIET_NAN ;
            mExpansionRatio = aRatio ;
            mNozzleMode     = NozzleMode::ComputeExitPressure ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_fuel_and_oxidizer(
                const Fuel     & aFuel,
                const Oxidizer & aOxidizer )
        {
            mFuel     = aFuel ;
            mOxidizer = aOxidizer ;

            if( ! mUserTemperatures )
            {
                this->set_fluid_temperatures();
            }

            // set the species
            this->set_species() ;

            // compute the mass fractions
            this->compute_mass_fractions() ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_fuel_and_oxidizer_temperatures(
                const real & aFuelTemp,
                const real & aOxTemp )
        {
            mUserTemperatures    = true ;
            mFuelTemperature     = aFuelTemp ;
            mOxidizerTemperature = aOxTemp ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_mixture_ratio( const real & aOF )
        {
            mOF = aOF ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_throat_diameter( const real & aD )
        {
            mThroatDiameter = aD ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::set_gas_generator_conditions( const real & aT, const real & aP )
        {
            mGasgeneratorTemperature = aT ;
            mGasgeneratorPressure    = aP ;
        }

//------------------------------------------------------------------------------

        Gas *
        Parameters::create_gas() const
        {
            Vector< real > tX( this->combgas_species().size(),
                               0.0 );
            tX( 0 ) = 1.0 ;
            return new Gas( this->combgas_species(), tX ) ;
        }

//------------------------------------------------------------------------------

        void
        Parameters::remix_gas( Gas * aGas, const real & aOF, const bool aRemixTransport ) const
        {
            BELFEM_ASSERT( aGas->number_of_components()
                == this->combgas_species().size(),
                "Number of components in gas is %u but should be %u",
                          ( unsigned int ) aGas->number_of_components(),
                          ( unsigned int ) this->combgas_species().size() );

            // create mixture
            Vector< real > tY( mOxidizerMassFractions * aOF + mFuelMassFractions );
            tY /= sum( tY );

            // remix
            aGas->remix_mass( tY, true, aRemixTransport );
        }

//------------------------------------------------------------------------------
// private
//------------------------------------------------------------------------------

        void
        Parameters::set_fluid_temperatures()
        {
            // set fuel temperature
            switch ( mFuel )
            {
                case( Fuel::LH2 ) :
                {
                    mFuelTemperature = 20.0 ;
                    mFuelMolarFractions.set_size( 1, 1.0 );
                    break ;
                }
                case( Fuel::LCH4 ) :
                {
                    mFuelTemperature = 110.0 ;
                    mFuelMolarFractions.set_size( 1, 1.0 );
                    break ;
                }
                case( Fuel::LNG ) :
                {
                    mFuelTemperature = 110.0 ;

                    // CH4, C2H6, C3H8, CO2
                    mFuelMolarFractions = { 0.935, 0.046, 0.012, 0.007 };

                    break ;
                }
                case( Fuel::C2H5OH ) :
                {
                    mFuelTemperature = BELFEM_TREF ;
                    mFuelMolarFractions.set_size( 1, 1.0 );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Fuel not supported: %s",
                                 fuel_to_string( mFuel ).c_str() );
                }
            }

            // set oxidizer temperature
            switch( mOxidizer )
            {
                case( Oxidizer::LOX ):
                {
                    mOxidizerTemperature = 90.0 ;
                    mOxidizerMolarFractions.set_size( 1, 1.0 );
                    break ;
                }
                default:
                {
                    BELFEM_ERROR( false, "Oxidizer not supported: %s",
                                 oxidizer_to_string( mOxidizer ).c_str() );
                }
            }

        }

//------------------------------------------------------------------------------

        void
        Parameters::set_species()
        {
            if( mOxidizer == Oxidizer::LOX )
            {
                switch ( mFuel )
                {
                    case( Fuel::LH2 ) :
                    {
                        mSpecies = { "H2","O2","H","H2O","H2O2","HO2","O","OH" } ;
                        mFuelIndices.set_size( 1, 0 );
                        mOxidizerIndices.set_size( 1, 1 );
                        break ;
                    }
                    case( Fuel::LCH4 ) :
                    {
                        mSpecies = {
                                "CH4",
                                "O2",
                                "C(gr)",
                                "C2H4",
                                "C2H6",
                                "C3H8",
                                "CH3CHO,ethanal",
                                "CH3OH",
                                "CO",
                                "CO2",
                                "CH2O",
                                "COOH",
                                "H",
                                "H2",
                                "H2O",
                                "H2O2",
                                "HCHO,formaldehy",
                                "HCO",
                                "HCOOH",
                                "HO2",
                                "O",
                                "OH"};

                        mFuelIndices.set_size( 1, 0 );
                        mOxidizerIndices.set_size( 1, 1 );
                        break ;
                    }
                    case( Fuel::LNG ) :
                    {
                        mSpecies = {
                                "CH4",
                                "C2H6",
                                "C3H8",
                                "CO2",
                                "O2",
                                "CH2O",
                                "CO",
                                "COOH",
                                "H",
                                "H2",
                                "H2O",
                                "H2O2",
                                "HCHO,formaldehy",
                                "HCO",
                                "HCOOH",
                                "O",
                                "O3",
                                "OH"
                                };

                        mFuelIndices = { 0, 1, 2, 3 };
                        mOxidizerIndices.set_size( 1, 4 );

                        break ;
                    }
                    case( Fuel::C2H5OH ) :
                    {
                        mSpecies = {
                                "C2H5OH",
                                "O2",
                                "C(gr)",
                                "C2H2,acetylene",
                                "C2H4",
                                "C2H6",
                                "CH2CO,ketene",
                                "CH3",
                                "CH3OH",
                                "CH4",
                                "CO",
                                "CO2",
                                "H",
                                "H2",
                                "H2O",
                                "H2O2",
                                "HCO",
                                "HCHO,formaldehy",
                                "HCOOH",
                                "O3",
                                "OH"
                        };

                        mFuelIndices.set_size( 1, 0 );
                        mOxidizerIndices.set_size( 1, 1 );
                        break ;
                    }
                    default:
                    {
                        BELFEM_ERROR( false, "Species have not been defined for %s-%s mixture",
                                     oxidizer_to_string( mOxidizer ).c_str(),
                                     fuel_to_string( mFuel ).c_str() );
                    }
                }
            }
            else
            {
                BELFEM_ERROR( false, "Species for Oxidizer have not been defined: %s",
                             oxidizer_to_string( mOxidizer ).c_str() );
            }

            // set the species lists for fuel and oxidizer
            mFuelSpecies.set_size( mFuelIndices.length(), "" );
            for( uint k=0; k<mFuelIndices.length(); ++k )
            {
                mFuelSpecies( k ) = mSpecies( mFuelIndices( k ) );
            }

            mOxidizerSpecies.set_size( mOxidizerIndices.length() , "" );
            for( uint k=0; k<mOxidizerIndices.length(); ++k )
            {
                mOxidizerSpecies( k ) = mSpecies( mOxidizerIndices( k ) );
            }
        }

//------------------------------------------------------------------------------

        void
        Parameters::compute_mass_fractions()
        {
            // number of species
            uint tN = this->combgas_species().size() ;
            // temporary vector with molar fractions
            Vector< real > tX( tN, 0.0 );
            tX( 0 ) = 1.0 ;

            // create a temporary gas so that we can get the molar masses
            Gas * tGas = new Gas( this->combgas_species(), tX ) ;

            // compute mass fractions for fuel
            uint tCount = 0 ;
            mFuelMassFractions.set_size( tN, 0.0 );
            for( uint k : mFuelIndices )
            {
                mFuelMassFractions( k ) = mFuelMolarFractions( tCount++ )
                        * tGas->component( k )->M() ;
            }
            mFuelMassFractions /= sum( mFuelMassFractions );

            // compute mass fractions for oxidizer
            tCount = 0 ;
            mOxidizerMassFractions.set_size( tN, 0.0 );
            for( uint k : mOxidizerIndices )
            {
                mOxidizerMassFractions( k ) = mOxidizerMolarFractions( tCount++ )
                                          * tGas->component( k )->M() ;
            }
            mOxidizerMassFractions /= sum( mOxidizerMassFractions );

            delete tGas ;
        }

//------------------------------------------------------------------------------
    }
}