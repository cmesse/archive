//
// Created by Christian Messe on 31.08.20.
//
#include "cl_EN_State.hpp"
#include "cl_EN_Analysis.hpp"
#include "cl_Gas.hpp"
#include "cl_Logger.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        State::State( Gas & aCombgas,
                      const string & aLabel,
                      const uint aNumSpecies ) :
                mCombgas( aCombgas ),
                mLabel( aLabel )
        {
            // allocate values object
            mValues.set_size( BELFEM_ENGINE_NUMSTATES, BELFEM_QUIET_NAN );

            if( aNumSpecies == 0 )
            {
                mMassFractions.set_size( aCombgas.number_of_components() );
                mMolarFractions.set_size( aCombgas.number_of_components() );
            }
            else
            {
                mMassFractions.set_size( aNumSpecies );
                mMolarFractions.set_size( aNumSpecies );
            }
        }


//------------------------------------------------------------------------------

        void
        State::compute_caloric( const real & aT, const real & aP, const real aU )
        {
            this->value( BELFEM_ENGINE_STATE_T )      = aT ;
            this->value( BELFEM_ENGINE_STATE_P )      = aP ;
            this->value( BELFEM_ENGINE_STATE_U )      = aU ;

            this->value( BELFEM_ENGINE_STATE_RHO )    = mCombgas.rho( aT, aP );

            this->value( BELFEM_ENGINE_STATE_R )      = mCombgas.R( aT, aP );
            this->value( BELFEM_ENGINE_STATE_M )      = mCombgas.M( aT, aP );

            this->value( BELFEM_ENGINE_STATE_H )      = mCombgas.h( aT, aP );
            this->value( BELFEM_ENGINE_STATE_CP )     = mCombgas.cp( aT, aP );
            this->value( BELFEM_ENGINE_STATE_S )      = mCombgas.s( aT, aP );
            this->value( BELFEM_ENGINE_STATE_GAMMA )  = mCombgas.gamma( aT, aP );

            this->value( BELFEM_ENGINE_STATE_W )      = mCombgas.c( aT, aP );

            this->value( BELFEM_ENGINE_STATE_MA )     = aU / this->value( BELFEM_ENGINE_STATE_W ) ;

            mMassFractions  = mCombgas.mass_fractions() ;
            mMolarFractions = mCombgas.molar_fractions() ;
        }

//------------------------------------------------------------------------------

        void
        State::compute_equilibrium(
                const real & aT,
                const real & aP,
                const real & aH )
        {
            // temperature of this state
            real & tT = this->value( BELFEM_ENGINE_STATE_T );

            // pressure of this state
            real & tP = this->value( BELFEM_ENGINE_STATE_P );

            // specific enthalpy
            real & tH = this->value( BELFEM_ENGINE_STATE_H );

            // specific heat capacity
            real & tCp = this->value( BELFEM_ENGINE_STATE_CP );

            // set initial values
            tT = aT ;
            tP = aP ;

            // relaxation factor
            real tOmega = 0.3 ;

            // epsilon criterion
            real tEpsilon = 1e-4 ;

            // counter for loop
            uint tCount = 0;

            real tDeltaT = BELFEM_REAL_MAX ;

            while ( std::abs( tDeltaT ) > tEpsilon )
            {

                // remix, but do not recompute the splines
                mCombgas.remix_to_equilibrium( tT, tP, true, false );

                tH = mCombgas.h( tT, tP );
                tCp = mCombgas.cp( tT, tP );

                // correction step
                tDeltaT = ( tH - aH ) / tCp;

                // relax step
                tT -= tOmega * tDeltaT;

                BELFEM_ERROR( tCount++ < 1000, "Failed to compute equilibrium temperature" );
            }

            // one final time for consisiency
            mCombgas.remix_to_equilibrium( tT, tP, true, true );

            // rember mass fractions
            mMassFractions  = mCombgas.mass_fractions() ;

            // remember molar fractions
            mMolarFractions = mCombgas.molar_fractions() ;
        }

//------------------------------------------------------------------------------

        void
        State::print() const
        {

            std::fprintf( stdout, "     State : %s\n\n", mLabel.c_str() );
            std::fprintf( stdout, "         Temperature    : %10.3f K\n\n",
                          this->value( BELFEM_ENGINE_STATE_T ) );

            std::fprintf( stdout, "         Pressure       : %10.3f bar\n\n",
                          this->value( BELFEM_ENGINE_STATE_P ) *1e-5 );

            if( ! std::isnan( this->value( BELFEM_ENGINE_STATE_RHO ) ) )
            {
                std::fprintf( stdout, "         Density        : %10.3f kg/m³\n\n",
                              this->value( BELFEM_ENGINE_STATE_RHO ));
            }

            if( ! std::isnan( this->value( BELFEM_ENGINE_STATE_M ) ) )
            {
                std::fprintf( stdout, "         Molar Mass     : %10.3f g/Mol\n\n",
                              this->value( BELFEM_ENGINE_STATE_M ) * 1000);
            }

            std::fprintf( stdout, "         Enthalpy       : %10.3f kJ/kg\n\n",
                          this->value( BELFEM_ENGINE_STATE_H ) *1e-3 );

            std::fprintf( stdout, "         Entropy        : %10.3f J/(kgK)\n\n",
                          this->value( BELFEM_ENGINE_STATE_S ) );

            if( ! std::isnan( this->value( BELFEM_ENGINE_STATE_GAMMA ) ) )
            {
                std::fprintf( stdout, "         Ratio of Heats : %10.3f -\n\n",
                              this->value( BELFEM_ENGINE_STATE_GAMMA ));
            }

            if( ! std::isnan( this->value( BELFEM_ENGINE_STATE_MA ) ) )
            {
                std::fprintf( stdout, "         Mach Number    : %10.3f -\n\n",
                              this->value( BELFEM_ENGINE_STATE_MA ));
            }

            if( ! std::isnan( this->value( BELFEM_ENGINE_STATE_U ) ) )
            {
                std::fprintf( stdout, "         Velocity      : %10.3f m/s\n\n",
                              this->value( BELFEM_ENGINE_STATE_U ));
            }

            if( ! std::isnan( this->value( BELFEM_ENGINE_STATE_A ) ) )
            {
                std::fprintf( stdout, "         Area Ratio     : %10.3f -\n\n",
                              this->value( BELFEM_ENGINE_STATE_A )  /
                                      this->value( BELFEM_ENGINE_STATE_A ) );
            }
        }

//------------------------------------------------------------------------------
    }
}
