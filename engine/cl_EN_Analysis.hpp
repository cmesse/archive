//
// Created by Christian Messe on 31.08.20.
//

#ifndef BELFEM_CL_EN_ANALYSIS_HPP
#define BELFEM_CL_EN_ANALYSIS_HPP

#include "typedefs.hpp"

#include "cl_Gas.hpp"
#include "cl_EN_Parameters.hpp"
#include "cl_EN_State.hpp"
#include "en_EN_Enums.hpp"

namespace belfem
{
    namespace engine
    {

        enum class IspMode
        {
            Sealevel,
            OptimalExpansion,
            Vacuum,
            UNDEFINED
        };

//------------------------------------------------------------------------------

        class Analysis
        {
            const Parameters & mParams ;
            Gas * mCombgas  = nullptr ;

            Gas * mFuel     = nullptr ;
            Gas * mOxidizer = nullptr ;

            // state of injector
            State mInjector ;

            // total state
            State mTotal ;

            // throat state
            State mThroat ;

            // nozzle exit
            State mNozzle ;

            // work matrix for jacobian
            Matrix< real > mJ2 ;

            // work marix for RHS
            Vector< real > mF2 ;

            // work matrix for jacobian
            Matrix< real > mJ3 ;

            // work marix for RHS
            Vector< real > mF3 ;

            // pivot for gesv
            Vector< int > mPivot ;

            // performance parameters:

            // mixture ratio as passed to compute_mixture
            real         mOF                  = BELFEM_QUIET_NAN ;

            // ambient pressure
            real mPambient = BELFEM_QUIET_NAN ;

            // mass flow
            real mDotM = BELFEM_QUIET_NAN ;

            // thrust coeffcient
            real mCF = BELFEM_QUIET_NAN ;

            // ISP
            real mISPref = BELFEM_QUIET_NAN ;
            real mISPsl  = BELFEM_QUIET_NAN ;
            real mISPvac = BELFEM_QUIET_NAN ;

            // thrust
            real mF   = BELFEM_QUIET_NAN ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Analysis( const Parameters & aParams ) ;

//------------------------------------------------------------------------------

            ~Analysis() ;

//------------------------------------------------------------------------------

            real
            run( const real aOF = 0.0,
                 const IspMode aMode=IspMode::OptimalExpansion ) ;

//------------------------------------------------------------------------------

            real
            compute_gas_generator( const real & aT,
                                   const real & aP,
                                   const real aOFmin=0.1,
                                   const real aOFmax=0.4,
                                   const real aT0=600.0 );

//------------------------------------------------------------------------------

            /**
             * expose the combusition gas object
             */
             Gas *
             combgas() ;

//------------------------------------------------------------------------------

            /**
             * expose the oxidizer object
             */
             Gas *
             oxidizer() ;

//------------------------------------------------------------------------------

            /**
             * expose the fuel object
             */
            Gas *
            fuel() ;

//------------------------------------------------------------------------------

            /**
             * expose the parameter object
             */
            const Parameters *
            params() ;

//------------------------------------------------------------------------------

            /**
             * expose the injector state
             */
            const State *
            injector() const ;

//------------------------------------------------------------------------------

            /**
             * expose the total state
             */
             const State *
             total() const ;

//------------------------------------------------------------------------------

            /**
             * expose the nozzle throat
             */
             const State *
             throat() const ;

//------------------------------------------------------------------------------

            /**
             * expose the nozzle exit
             */
            const State *
            nozzle() const ;

//------------------------------------------------------------------------------

            void
            print_performance();

//------------------------------------------------------------------------------

            real
            find_best_mixture(
                    const real & aOFmin,
                    const real & aOFmax,
                    const IspMode aMode = IspMode::OptimalExpansion );

//------------------------------------------------------------------------------

            /**
             * return the massflow
             */
             const real &
             massflow() const ;

//------------------------------------------------------------------------------

            void
            compute_injector( const real & aOF, const real & aP );

//------------------------------------------------------------------------------

            void
            compute_total( const real aInitialTemperatureGuess=2000.0 );

//------------------------------------------------------------------------------

            const real &
            isp_sl_opt() const ;

//------------------------------------------------------------------------------

            const real &
            isp_vac_opt() const ;

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            create_gasmodels();

//------------------------------------------------------------------------------

            void
            compute_throat();

//------------------------------------------------------------------------------

            void
            compute_nozzle_pconst();

//------------------------------------------------------------------------------

            void
            compute_nozzle_aconst();

//------------------------------------------------------------------------------

            void
            compute_performance();

//------------------------------------------------------------------------------

            void
            create_initial_mixture( const real & aOF, const real & aP );

//------------------------------------------------------------------------------

            // help function for gas generator mode
            real
            compute_combustion_temperature(
                    const real & aP,
                    const real & aOF,
                    const real   aT0=600.0 );

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline Gas *
        Analysis::combgas()
        {
            return mCombgas ;
        }

//------------------------------------------------------------------------------

        inline Gas *
        Analysis::oxidizer()
        {
            return mOxidizer ;
        }

//------------------------------------------------------------------------------

        inline Gas *
        Analysis::fuel()
        {
            return mFuel ;
        }

//------------------------------------------------------------------------------

        inline const Parameters *
        Analysis::params()
        {
            return & mParams ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Analysis::injector() const
        {
            return & mInjector ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Analysis::total() const
        {
            return & mTotal ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Analysis::throat() const
        {
            return & mThroat ;
        }

//------------------------------------------------------------------------------

        inline const State *
        Analysis::nozzle() const
        {
            return & mNozzle ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Analysis::massflow() const
        {
            return mDotM ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Analysis::isp_sl_opt() const
        {
            return mISPsl ;
        }

//------------------------------------------------------------------------------

        inline const real &
        Analysis::isp_vac_opt() const
        {
            return mISPvac ;
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EN_ANALYSIS_HPP
