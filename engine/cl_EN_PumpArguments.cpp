//
// Created by Christian Messe on 29.04.21.
//
#include "assert.hpp"
#include "stringtools.hpp"

#include "cl_Communicator.hpp"
#include "cl_EN_PumpArguments.hpp"

extern belfem::Communicator gComm;

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        PumpArguments::PumpArguments( int & argc, char ** argv )  :
                belfem::Arguments( argc, argv )
        {
            if ( mArguments.size() <= 1 )
            {
                mState = RunState::PrintUsage ;
            }
            else
            {
                mState = RunState::Compute ;

                // initialize counter
                uint tCount = 0 ;

                for ( string tArg: mArguments )
                {

                    // increment counter
                    ++tCount ;

                    if ( tArg == "-h" || tArg == "--help" )
                    {
                        mState = RunState::PrintHelp;
                        break;
                    }
                    if ( tArg == "-l" || tArg == "--library" )
                    {
                        // Check that Library was given
                        BELFEM_ERROR( tCount < mArguments.size(),
                                     "You must specify a library!" );

                        // get argument
                        mLibraryPath = mArguments( tCount );

                    }
                    if ( tArg == "-f" || tArg == "--fluid" )
                    {
                        // Check that Library was given
                        BELFEM_ERROR( tCount < mArguments.size(),
                                     "You must specify a fluid!" );

                        string tString = string_to_lower( mArguments( tCount ) );

                        if ( tString == "lox" || tString == "lo2" || tString == "o2"  || tString == "oxygen" )
                        {
                            mFluid = HelmholtzModel::Oxygen;
                        }
                        else if ( tString == "lh2" || tString == "h2"  || tString == "hydrogen" )
                        {
                            mFluid = HelmholtzModel::NormalHydrogen ;
                        }
                        else if ( tString == "lch4" || tString == "ch4"  || tString == "methane" )
                        {
                            mFluid = HelmholtzModel::Methane ;
                        }
                        else
                        {
                            BELFEM_ERROR( false , "Unkown fluid type: %s \n must be either LOX, LCH4, LH2.", mArguments( tCount ).c_str() );
                        }
                    }
                    if ( tArg == "-s" || tArg == "--symbol" )
                    {
                        // Check that Library was given
                        BELFEM_ERROR( tCount < mArguments.size(),
                                     "You must specify a symbol!" );

                        // get argument
                        mSymbolName = mArguments( tCount );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
    } /* end namespace engine */
} /* end namespace belfem */