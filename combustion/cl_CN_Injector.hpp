//
// Created by Christian Messe on 23.04.20.
//

#ifndef BELFEM_CL_CN_INJECTOR_HPP
#define BELFEM_CL_CN_INJECTOR_HPP
#include "typedefs.hpp"
#include "cl_Gas.hpp"
#include "cl_CN_Scheme.hpp"

namespace belfem
{
    namespace combustion
    {
        class Injector
        {
            // Combustion scheme for this injector
            Scheme & mScheme ;

            // gas connected to this injector
            Gas & mGas ;

            // position where injection begins
            const real mXinj ;

            // duct height
            const real mDuctHeight ;

            // mixing efficienty
            const real mEtaMix ;

            // pulsonetti factor
            const real mPulsonetti ;

            // stoichiometric oxidizer fuel ratio
            const real mOFst ;

            const real mGamma = 0.01 ; // value of 1 - eta_mix_hat

            // log ( 1 - eta_mix_hat )
            const real mLnGamma ;

            real mMixingLength = BELFEM_QUIET_NAN ;

            // flag that tells if mixing length is computed automatically
            bool mAutoMixingLength = true ;

            // stoiciometric ratio
            real mPhi = BELFEM_QUIET_NAN ;

            // oxidizer fuel ratio
            real mOF = BELFEM_QUIET_NAN ;

            // mass flows
            real mOxidizerMassflow = BELFEM_QUIET_NAN ;
            real mFuelMassflow     = BELFEM_QUIET_NAN ;

            // work vector for mass flows
            Vector< real > mY ;

            // work vector for mass fractions
            Vector< real > mX ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Injector(
                    Scheme & aScheme,
                    const real aXinj,
                    const real aDuctHeight,
                    const real aEtaMix=0.8,
                    const real aPulsonetti=25.0 );

            ~Injector() = default ;

//----------------------------------------------------------------------------

            /**
             * print some debug info
             */
             void
             print();

//------------------------------------------------------------------------------

            // set stoichiometric ratio ( overrides set_of )
            void
            set_phi( const real & aPhi );

//------------------------------------------------------------------------------

            // set oxidizer to fuel ratio ( overrides set_phi )
            void
            set_of( const real & aOF );

//------------------------------------------------------------------------------

            // set the massflow of the oxidizer. Must be called after
            // set_of or set_phi ( overrides set_massflow )
            void
            set_oxidizer_massflow( const real & aOxidizerMassflow );

//------------------------------------------------------------------------------

            // set the total massflow of the channel. Must be called after
            // set_of or set_phi ( overrides set_oxidizer_massflow )
            void
            set_massflow( const real & aMassflow );

//------------------------------------------------------------------------------

            // set mixing length ( overrides automatic estimation )
            void
            set_mixing_length( const real & aMixingLength );

//------------------------------------------------------------------------------

            // returns the oxidizer to fuel ratio
            const real &
            of() const ;

//------------------------------------------------------------------------------

            // returns the stoichiometric ratio
            const real &
            phi() const ;

//------------------------------------------------------------------------------

            // returns the mixing length
            const real &
            mixing_length() const ;

//------------------------------------------------------------------------------

            // move fuel from inert to reacting
            void
            inject( const real & aX );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            compute_mixing_length() ;

//------------------------------------------------------------------------------

            real
            compute_stoich_fuel_oxidizer_ratio( Scheme & aScheme );

//------------------------------------------------------------------------------

            // mixing function
            real
            mix( const real & aX ) const;

            // derivative of mixing function
            real
            dmix( const real & aX ) const;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline const real &
        Injector::of() const
        {
            return mOF ;
        }
//------------------------------------------------------------------------------

        inline const real &
        Injector::phi() const
        {
            return mPhi ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Injector::mixing_length() const
        {
            return mMixingLength ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CN_INJECTOR_HPP
