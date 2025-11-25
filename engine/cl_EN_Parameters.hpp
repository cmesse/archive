//
// Created by Christian Messe on 31.08.20.
//

#ifndef BELFEM_CL_EN_PARAMETERS_HPP
#define BELFEM_CL_EN_PARAMETERS_HPP

#include "typedefs.hpp"
#include "constants.hpp"
#include "CN_Enums.hpp"
#include "cl_Cell.hpp"
#include "cl_Vector.hpp"
#include "en_EN_Enums.hpp"
#include "cl_Gas.hpp"
#include "cl_HDF5.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        class Parameters
        {
            // desired temperature in gas generator
            real mGasgeneratorTemperature = 850 ;

            // pressure in gas generator
            real mGasgeneratorPressure =  85e5 ;

            // combustion pressure in Pa
            real mChamberPressure       = 70e5 ;

            // ambient perssure in Pa
            real mExitPressure       = 0.5e5 ;

            // expansion ratio, if set
            real mExpansionRatio     = BELFEM_QUIET_NAN ;

            // fuel oxidizer ratio
            real mOF                    = 3.2 ;

            // reference thrust in kN
            real mDesignThrust = 60e3 ;

            // nozzle cross section
            real mThroatDiameter = 1.0 / constant::pi ;

            // fuel for this engine
            Fuel         mFuel = Fuel::LCH4 ;

            // oxidizer for this engine
            Oxidizer mOxidizer = Oxidizer::LOX ;

            // flag that tells if temperatures have been set by user
            bool         mUserTemperatures = false ;

            // temperature for oxidizer
            real         mOxidizerTemperature = BELFEM_QUIET_NAN ;

            // temperature for fuel
            real         mFuelTemperature     = BELFEM_QUIET_NAN ;

            Cell< string > mSpecies ;

            Vector< index_t > mFuelIndices ;
            Vector< index_t > mOxidizerIndices ;

            // molar fractions for fuel, length of mFuelIndices
            Vector< real    > mFuelMolarFractions ;

            // molar fractions of oxidizer, length of mOxidierIndices
            Vector< real    > mOxidizerMolarFractions ;

            Cell< string > mFuelSpecies ;
            Cell< string > mOxidizerSpecies ;

            // mass fractions for fuel, length of mSpecies
            Vector< real > mFuelMassFractions ;

            // mass fractions of oxidizer, length of mSpecies
            Vector< real > mOxidizerMassFractions ;

            NozzleMode     mNozzleMode = NozzleMode::ComputeCrossSection ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Parameters() ;

//------------------------------------------------------------------------------

            // special constructor if meshinfo file is given
            Parameters( HDF5 & aDatabase ) ;

//------------------------------------------------------------------------------

            ~Parameters() = default ;

//------------------------------------------------------------------------------
// setters
//------------------------------------------------------------------------------


            /**
             * reference pressure in Pa
             */
            void
            set_thrust( const real & aF );

//------------------------------------------------------------------------------
            /**
             * total chamber pressure in Pa
             */
            void
            set_chamber_pressure( const real & aP );

//------------------------------------------------------------------------------

            /**
             * pressure for nozzle exit pressure in Pa
             *
             * overwrites set_expansion_ratio
             */
            void
            set_exit_pressure( const real & aP );

//------------------------------------------------------------------------------

            /**
             * sets the expansion ratio to a fixed value
             *
             * overwrites set_exit_pressure
             * @param aP
             */
            void
            set_expansion_ratio( const real & aRatio );

//------------------------------------------------------------------------------

            void
            set_fuel_and_oxidizer(
                    const Fuel & aFuel,
                    const Oxidizer & aOxidizer );

//------------------------------------------------------------------------------

            void
            set_fuel_and_oxidizer_temperatures(
                    const real & aFuelTemp,
                    const real & aOxTemp );

//------------------------------------------------------------------------------

            void
            set_mixture_ratio( const real & aOF );

//------------------------------------------------------------------------------

            void
            set_throat_diameter( const real & aD );

//------------------------------------------------------------------------------

            void
            set_gas_generator_conditions( const real & aT, const real & aP );

//------------------------------------------------------------------------------

            /**
             * creates a gas based on the selected mixture
             * @return
             */
            Gas *
            create_gas() const;

//------------------------------------------------------------------------------

            /**
             * remixes a gas, using the given OF ratio
             * @return
             */
             void
             remix_gas( Gas * aGas, const real & aOF,
                        const bool aRemixTransport=false ) const;

//------------------------------------------------------------------------------
// getters
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

            const real &
            chamber_pressure() const ;

//------------------------------------------------------------------------------

            const real &
            exit_pressure() const ;

//------------------------------------------------------------------------------

            const real &
            expansion_ratio() const ;

//------------------------------------------------------------------------------

            const NozzleMode &
            nozzle_mode() const ;

//------------------------------------------------------------------------------

            const Fuel &
            fuel() const ;

//------------------------------------------------------------------------------

            const Oxidizer &
            oxidizer() const ;

//------------------------------------------------------------------------------

            const real &
            fuel_temperature() const ;

//------------------------------------------------------------------------------

            const real &
            oxidizer_temperature() const ;

//------------------------------------------------------------------------------

            const Vector< real > &
            oxidizer_molar_fractions() const;

//------------------------------------------------------------------------------

            const Vector< real > &
            fuel_molar_fractions() const;

//------------------------------------------------------------------------------

            const real  &
            thrust() const;

//------------------------------------------------------------------------------

            const Cell< string > &
            combgas_species() const ;

//------------------------------------------------------------------------------

            const Cell< string > &
            fuel_species() const ;

//------------------------------------------------------------------------------

            const Cell< string > &
            oxidizer_species() const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            fuel_indices() const ;

//------------------------------------------------------------------------------

            const Vector< index_t > &
            oxidizer_indices() const ;

//------------------------------------------------------------------------------

            const real &
            mixture_ratio() const ;

//------------------------------------------------------------------------------

            const real &
            throat_diameter() const ;

//------------------------------------------------------------------------------

            uint
            number_of_species() const ;

//------------------------------------------------------------------------------

            const real &
            gasgenerator_temperature() const ;

//------------------------------------------------------------------------------

            const real &
            gasgenerator_pressure() const ;

//------------------------------------------------------------------------------

            const Vector< real > &
            fuel_mass_fractions() const ;

//------------------------------------------------------------------------------

            const Vector< real > &
            oxidizer_mass_fractions() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            set_fluid_temperatures();

//-----------------------------------------------------------------------------

            void
            set_species();

//-----------------------------------------------------------------------------

            void
            compute_mass_fractions() ;

//-----------------------------------------------------------------------------
        };

//-----------------------------------------------------------------------------

        inline const real  &
        Parameters::thrust() const
        {
            return mDesignThrust ;
        }

//-----------------------------------------------------------------------------

        inline const Fuel &
        Parameters::fuel() const
        {
            return mFuel ;
        }

//-----------------------------------------------------------------------------

        inline const Oxidizer &
        Parameters::oxidizer() const
        {
            return mOxidizer ;
        }

//-----------------------------------------------------------------------------

        inline const real &
        Parameters::fuel_temperature() const
        {
            return mFuelTemperature ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::oxidizer_temperature() const
        {
            return mOxidizerTemperature ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::chamber_pressure() const
        {
            return mChamberPressure ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::exit_pressure() const
        {
            return mExitPressure ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::expansion_ratio() const
        {
            return mExpansionRatio ;
        }

//------------------------------------------------------------------------------

        inline const NozzleMode &
        Parameters::nozzle_mode() const
        {
            return mNozzleMode ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Parameters::oxidizer_molar_fractions() const
        {
            return mOxidizerMolarFractions ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Parameters::fuel_molar_fractions() const
        {
            return mFuelMolarFractions ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        Parameters::combgas_species() const
        {
            return mSpecies ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        Parameters::fuel_species() const
        {
            return mFuelSpecies ;
        }

//------------------------------------------------------------------------------

        inline const Cell< string > &
        Parameters::oxidizer_species() const
        {
            return mOxidizerSpecies ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        Parameters::fuel_indices() const
        {
            return mFuelIndices ;
        }

//------------------------------------------------------------------------------

        inline const Vector< index_t > &
        Parameters::oxidizer_indices() const
        {
            return mOxidizerIndices ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::mixture_ratio() const
        {
            return mOF ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::throat_diameter() const
        {
            return mThroatDiameter ;
        }

//------------------------------------------------------------------------------

        uint
        inline Parameters::number_of_species() const
        {
            return this->combgas_species().size() ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::gasgenerator_temperature() const
        {
            return mGasgeneratorTemperature ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Parameters::gasgenerator_pressure() const
        {
            return mGasgeneratorPressure ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Parameters::fuel_mass_fractions() const
        {
            return mFuelMassFractions ;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Parameters::oxidizer_mass_fractions() const
        {
            return mOxidizerMassFractions ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EN_PARAMETERS_HPP
