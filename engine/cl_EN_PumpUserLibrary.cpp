//
// Created by Christian Messe on 29.04.21.
//
// dynamic linker function
#include "dlfcn.h"

#include "assert.hpp"

#include "cl_EN_PumpUserLibrary.hpp"
#include "filetools.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        Library::Library( const string & aPath )
        {
            // check if path is relative
            if ( aPath.at( 0 ) != '/' )
            {
                if (aPath.at( 0 ) == '.' )
                {
                    mPath = gComm.workdir() + aPath.substr( 1 );
                }
                else
                {
                    mPath = gComm.workdir() + "/" +  aPath;
                }
            }
            else // path is absolute
            {
                mPath = aPath ;
            }

            // test if file exists
            BELFEM_ERROR( file_exists( aPath ), "File %s does not exist", aPath.c_str() );

            // try to open library file
            mLibraryHandle = dlopen( mPath.c_str(), RTLD_NOW );

            // test if loading succeeded
            if( ! mLibraryHandle )
            {
                // get error string
                std::string tError = dlerror();

                // throw error
                BELFEM_ERROR( mLibraryHandle, tError.c_str() );
            }

        }

//------------------------------------------------------------------------------

        Library::~Library()
        {
            // close handle to library
            dlclose( mLibraryHandle );
        }

//------------------------------------------------------------------------------

        BELFEM_PUMP_USER_FUNCTION
        Library::load_function( const string & aFunctionName )
        {
            // try to load function from library
            BELFEM_PUMP_USER_FUNCTION aUserFunction
                    = reinterpret_cast<BELFEM_PUMP_USER_FUNCTION>
                    ( dlsym( mLibraryHandle, aFunctionName.c_str() ) );

            // create error message in case we fail
            string tError =  "Could not find symbol " + aFunctionName
                                  + "  within file " + mPath;

            // throw error if if loading succeeded
            BELFEM_ERROR( aUserFunction, tError.c_str() );

            // return the function handle
            return aUserFunction ;
        }

//------------------------------------------------------------------------------
    }
}