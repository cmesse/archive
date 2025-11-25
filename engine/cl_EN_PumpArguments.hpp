//
// Created by Christian Messe on 29.04.21.
//

#ifndef BELFEM_CL_EN_PUMPARGUMENTS_HPP
#define BELFEM_CL_EN_PUMPARGUMENTS_HPP

#include "typedefs.hpp"
#include "cl_Arguments.hpp"
#include "en_Helmholtz.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        enum class RunState
        {
            PrintHelp    = 0,
            PrintUsage   = 1,
            Compute      = 2,
            Undefined    = 3
        };

//------------------------------------------------------------------------------


        class PumpArguments : public Arguments
        {
            // run state of main
            RunState   mState           = RunState::Undefined ;

            // path to shared object
            string  mLibraryPath     = "" ;

            // name of symbol to load from library
            string mSymbolName = "" ;

            // fluid for pump
            HelmholtzModel mFluid    = HelmholtzModel::UNDEFINED ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            PumpArguments( int & argc, char * argv[] );

//------------------------------------------------------------------------------

            ~PumpArguments() = default;

//------------------------------------------------------------------------------

            /**
             * return the run state
             */
            const RunState &
            state() const;

//------------------------------------------------------------------------------

            /**
             * return the selected fluid
             */
            const HelmholtzModel &
            fluid() const;

//------------------------------------------------------------------------------

            /**
             * return the library path
             */
             const string &
             library_path() const ;

//------------------------------------------------------------------------------

            /**
             * return the name of the function to be loaded
             */
            const string &
            symbol_name() const ;

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------

        inline const RunState &
        PumpArguments::state() const
        {
            return mState ;
        }

//------------------------------------------------------------------------------

        inline const HelmholtzModel &
        PumpArguments::fluid() const
        {
            return mFluid ;
        }

//------------------------------------------------------------------------------

        inline const string &
        PumpArguments::symbol_name() const
        {
            return mSymbolName ;
        }

//------------------------------------------------------------------------------

        inline const string &
        PumpArguments::library_path() const
        {
            return mLibraryPath ;
        }

//------------------------------------------------------------------------------
    } /* end namespace engine */
} /* end namespace belfem */

#endif //BELFEM_CL_EN_PUMPARGUMENTS_HPP
