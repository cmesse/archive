//
// Created by Christian Messe on 25.09.19.
//

#ifndef BELFEM_CL_CN_SCHEME_HPP
#define BELFEM_CL_CN_SCHEME_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Gas.hpp"

#include "CN_Enums.hpp"
#include "en_GM_GasModel.hpp"
#include "cl_CN_Reaction.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        class Chemkin;

//------------------------------------------------------------------------------

        class Scheme
        {
//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            const string   mChemkinFileName;
            const Fuel     mFuel;
            const Oxidizer mOxidizer;
            const GasModel mGasModel;

            const real     mEpsilonY = 1e-12;

            Gas * mCombgas = nullptr;

            Vector< real >  mInitialMolarFractions;

            uint mReactingFuelIndex;
            uint mInertFuelIndex;
            uint mOxidizerIndex ;

            uint mNumberOfReactingSpecies;
            uint mNumberOfAllSpecies;
            uint mTemperatureIndex ;

            // specific volume of gas
            real mV;

            // Y
            Vector< real > mY;
            Vector< real > mY0;

            Vector< real > mdYdt;
            real mdTdt ;

            // rhs of equation system
            Vector< real > mRHS ;

            // lhs of equation system
            Vector< real > mLHS ;

            // pivot for gesv
            Vector< int > mPivot ;

            // concentration
            Vector< real > mC;
            Vector< real > mM;

            //Vector< real > mSpecificVolumes;
            Vector< real > mH;
            Vector< real > mH0 ;
            Vector< real > mCp;
            Vector< real > mdCpdT;

            // Work vector with gibbs data
            Vector< real >  mGibbs;
            Vector< real >  mdGibbsdT;

            Matrix< real > mJacobi;

            Cell< Reaction * > mReactions;

            uint mCount = 0 ;

            // factors for Jacobian stabilization
            real mC1 = BELFEM_QUIET_NAN ;
            real mC2 = BELFEM_QUIET_NAN ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Scheme( const string   &     aChemkinFilePath,
                    const Fuel           aFuel,
                    const Oxidizer       aOxidizer=Oxidizer::LOX,
                    const GasModel       aGasModel=GasModel::IDGAS );

//------------------------------------------------------------------------------

            ~Scheme();

//------------------------------------------------------------------------------

            /**
             * expose combustion gas
             */
            inline Gas *
            combgas();

//------------------------------------------------------------------------------

            real
            compute( const real & aT, const real & aP, const real & aU, const real & aDeltaX  );

//------------------------------------------------------------------------------

            /**
             * reset the mixture of the combustion gas
             */
            void
            reset_combgas_mixture();

//------------------------------------------------------------------------------

            /**
             * tell how many species contribute to the reaction
             */
            inline const uint &
            number_of_reacting_species() const;

//------------------------------------------------------------------------------

            /**
             * mass fractions
             */
            const Vector< real > &
            Y() const ;

            void
            set_Y0( const Vector< real > & aY );

            inline const real &
            Y( const uint & aIndex ) const;

            inline const real &
            dYdt( const uint & aIndex ) const;

            inline const Vector< real > &
            dYdt() const;

//------------------------------------------------------------------------------

            /**
            *  concentration
            */
            inline const real &
            c( const uint & aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * specific volume in ccm/kg
             */
             inline const real &
             v( const uint & aIndex ) const;

             inline const real &
             v() const;

//------------------------------------------------------------------------------

            /**
             * gibbs
             */
            inline const real &
            G( const uint & aIndex ) const;

            inline const real &
            dGdT( const uint & aIndex ) const;

//------------------------------------------------------------------------------

            /**
             * for channel
             * @return
             */
             real
             delta_w() ;

//------------------------------------------------------------------------------

             real
             delta_R() ;

//------------------------------------------------------------------------------
// Compinent indices
// ------------------------------------------------------------------------------

            inline const uint &
            reacting_fuel_index() const ;

//------------------------------------------------------------------------------

            const uint &
            inert_fuel_index() const ;

//------------------------------------------------------------------------------

            const uint &
            oxidizer_index() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_combgas( Chemkin & aChemkinFile );

//------------------------------------------------------------------------------

            void
            add_inert_fuel( Cell< string > & aSpecies );

//------------------------------------------------------------------------------

            void
            add_air(
                    Cell< string > & aSpecies,
                    Vector< real > & aMolarFractions );

//------------------------------------------------------------------------------

            void
            add_oxidizer( Cell< string > & aSpecies,
                          Vector< real > & aMolarFractions );

//------------------------------------------------------------------------------

            /**
             * preprocess for jacobi calculation
             */
            void
            preprocess( const real & aT, const real & aP );

            void
            compute_jacobi( const real & aT, const real & aP );

            void
            compute_rhs( const real & aT, const real & aP, const real & aU, const real & aDeltaX  );

//------------------------------------------------------------------------------

            void
            find_oxidizer_index() ;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        Gas *
        Scheme::combgas()
        {
            return mCombgas;
        }

//------------------------------------------------------------------------------

        const uint &
        Scheme::number_of_reacting_species() const
        {
            return mNumberOfReactingSpecies;
        }

//------------------------------------------------------------------------------

        inline const Vector< real > &
        Scheme::Y() const
        {
            return mY ;
        }

//------------------------------------------------------------------------------

        inline void
        Scheme::set_Y0( const Vector< real > & aY )
        {
            mY0 = aY ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Scheme::Y( const uint & aIndex ) const
        {
            return mY( aIndex );
        }

//------------------------------------------------------------------------------

        inline const real &
        Scheme::dYdt( const uint & aIndex ) const
        {
            return mdYdt( aIndex );
        }

//------------------------------------------------------------------------------

        inline const  Vector< real > &
        Scheme::dYdt() const
        {
            return mdYdt;
        }

//------------------------------------------------------------------------------

        const real &
        Scheme::c( const uint & aIndex ) const
        {
            return mC( aIndex );
        }

//------------------------------------------------------------------------------

        /*const real &
        Scheme::v( const uint & aIndex ) const
        {
            return mSpecificVolumes( aIndex );
        }*/
//------------------------------------------------------------------------------

        const real &
        Scheme::v() const
        {
            return mV;
        }

//------------------------------------------------------------------------------

        inline const real &
        Scheme::G( const uint & aIndex ) const
        {
            return mGibbs( aIndex );
        }

//------------------------------------------------------------------------------

        inline const real &
        Scheme::dGdT( const uint & aIndex ) const
        {
            return mdGibbsdT( aIndex );
        }

//------------------------------------------------------------------------------

        inline const uint &
        Scheme::reacting_fuel_index() const
        {
            return mReactingFuelIndex ;
        }

//------------------------------------------------------------------------------

        inline const uint &
        Scheme::inert_fuel_index() const
        {
            return mInertFuelIndex ;
        }

//------------------------------------------------------------------------------

        inline const uint &
        Scheme::oxidizer_index() const
        {
            return mOxidizerIndex ;
        }

//------------------------------------------------------------------------------

    }
}

#endif //BELFEM_CL_CN_SCHEME_HPP
