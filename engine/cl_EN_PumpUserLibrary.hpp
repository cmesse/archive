//
// Created by Christian Messe on 29.04.21.
//

#ifndef BELFEM_CL_EN_PUMPUSERLIBRARY_HPP
#define BELFEM_CL_EN_PUMPUSERLIBRARY_HPP


#include "typedefs.hpp"
#include "cl_Gas.hpp"
#include "cl_EN_Pump.hpp"

namespace belfem
{
    namespace engine
    {
//------------------------------------------------------------------------------

        /**
          * Interface for user defined function
          */
        typedef
        void ( *BELFEM_PUMP_USER_FUNCTION )
        (
                    Pump & aPump
        );

//------------------------------------------------------------------------------

        /**
                * Wrapper class for shared object file
                */
        class Library
        {
            // path to library file
            string mPath;

            // handle to shared object
            void * mLibraryHandle;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /**
             * create the library
             */
            Library( const string & aPath ) ;

//------------------------------------------------------------------------------

            /**
             * destroy the library
             */
            ~Library();

//------------------------------------------------------------------------------

            /**
             * load the user defined function
             */
            BELFEM_PUMP_USER_FUNCTION
            load_function( const string & aFunctionName );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_EN_PUMPUSERLIBRARY_HPP
