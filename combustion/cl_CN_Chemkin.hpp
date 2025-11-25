//
// Created by Christian Messe on 23.09.19.
//

#ifndef BELFEM_CL_CN_CHEMKIN_HPP
#define BELFEM_CL_CN_CHEMKIN_HPP

#include "typedefs.hpp"
#include "cl_Map.hpp"
#include "cl_Ascii.hpp"
#include "cl_Vector.hpp"
#include "cl_CN_Entry.hpp"

namespace belfem
{
    namespace combustion
    {
        // http://akrmys.com/public/chemkin/CKm_inp.html.en
        class Chemkin : public Ascii
        {
            uint mStartTag;
            uint mEndTag;
            uint mNumberOfEntries;

            // an index map needed to identify positions in the file
            Vector< uint > mIndex;

            uint mNumberOfReactions;

            Cell< Entry * >  mEntries;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------


            Chemkin( const string & aPath );

//------------------------------------------------------------------------------

            ~Chemkin();

//------------------------------------------------------------------------------

            /**
             * create a Cell of all species that exist in the reactions
             * @param aSpecies
             */
            void
            get_species( Cell< string > & aSpecies );

//------------------------------------------------------------------------------

            /**
             * how many ( non-duplicate ) reactions exist in this file
             * @return
             */
            inline uint
            number_of_reactions() const;

//-----------------------------------------------------------------------------

            /**
             * how many ( including duplicate ) reactions exist in this file
             * @return
             */
            inline uint
            number_of_entries() const;

//-----------------------------------------------------------------------------

            // get an entry from the list
            inline Entry *
            entry( const uint & aIndex );

//------------------------------------------------------------------------------
        private:
//------------------------------------------------------------------------------

            void
            find_tags();

//------------------------------------------------------------------------------

            void
            count_entries();

//------------------------------------------------------------------------------

            void
            read_entries();

//------------------------------------------------------------------------------

            void
            find_duplicates();

//------------------------------------------------------------------------------

            void
            count_active_reactions();

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        uint
        Chemkin::number_of_reactions() const
        {
            return mNumberOfReactions;
        }

//------------------------------------------------------------------------------

        uint
        Chemkin::number_of_entries() const
        {
            return mNumberOfEntries;
        }

//------------------------------------------------------------------------------

        Entry *
        Chemkin::entry( const uint & aIndex )
        {
            return mEntries( aIndex );
        }

//------------------------------------------------------------------------------
    }
}
#endif //BELFEM_CL_CHEMKIN_HPP
