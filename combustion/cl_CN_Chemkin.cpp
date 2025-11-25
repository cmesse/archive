//
// Created by Christian Messe on 23.09.19.
//

#include "cl_CN_Chemkin.hpp"
#include "cl_GT_GasData.hpp"
#include "stringtools.hpp"

#include "fn_CN_molars_from_string.hpp"
#include "fn_GT_fix_label.hpp"

namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        Chemkin::Chemkin( const string & aPath ) :
            Ascii(aPath, FileMode::OPEN_RDONLY )
        {
            this->find_tags();
            this->count_entries();
            this->read_entries();
            this->find_duplicates();
            this->count_active_reactions();
        }

//------------------------------------------------------------------------------

        Chemkin::~Chemkin()
        {
            for ( Entry * tEntry: mEntries )
            {
                delete tEntry;
            }

            mEntries.clear();
        }

//------------------------------------------------------------------------------
        void
        Chemkin::find_tags()
        {
            uint tCount=0;

            mStartTag = 0;
            mEndTag = mBuffer.size();

            // find start tag
            for ( const string & tLine : mBuffer )
            {
                // get first word in line
                string tWord = string_to_upper(  first_word( tLine ) );

                if( tWord == "REACTIONS" )
                {
                    mStartTag = tCount + 1;
                }
                else if ( tWord == "END" && mStartTag > 0 )
                {
                    mEndTag = tCount;
                    break;
                }
                ++tCount;
            }
        }

//------------------------------------------------------------------------------

        void
        Chemkin::count_entries()
        {
            uint tCount = 0;

            // reserve memory for index
            mIndex.set_size( mEndTag - mStartTag, mEndTag );

            for( uint k=mStartTag; k<mEndTag; ++k )
            {
                // split the line into words
                Cell< string > tWords = string_to_words(  string_to_upper(  mBuffer( k ) ) );

                // assemble reaction
                string tReaction = "";

                // test if this line is commented out
                if ( tWords(0).at( 0 ) != '!' )
                {
                    if( tWords.size() > 3 )
                    {
                        for ( uint i = 0; i < tWords.size() - 3; ++i )
                        {
                            tReaction += tWords( i );
                        }

                        // test if this is a reaction
                        uint tPos = tReaction.find( "=" );

                        if ( tPos < tReaction.length() )
                        {
                            // remember position in file
                            mIndex( tCount ) = k;

                            // increment counter
                            ++tCount;
                        }
                    }
                }
            }

            // remember number of reactions
            mNumberOfEntries = tCount;
        }

//------------------------------------------------------------------------------

        void
        Chemkin::read_entries()
        {
            // crete an empty entry
            Entry tEmpty;

            // allocate memory
            mEntries.set_size( mNumberOfEntries, nullptr );

            // loop over all entries in list
            for( uint k=0; k<mNumberOfEntries; ++k )
            {
                // get entry
                Entry * tEntry = new Entry();

                // get line index
                uint tIndex = mIndex( k );

                // split the line into words
                Cell< string > tWords = string_to_words(  string_to_upper(  mBuffer( tIndex ) ) );

                BELFEM_ASSERT( tWords.size() > 3, "something went wring while reading chemkin file" );

                // assemble Reaction string
                string & tReaction = tEntry->reaction();

                for ( uint i = 0; i < tWords.size() - 3; ++i )
                {
                    tReaction += tWords( i );
                }

                // fix equal sign
                tReaction = search_and_replace( tReaction, "<", "" );
                tReaction = search_and_replace( tReaction, ">", "" );

                // fix third party
                bool tThirdBodyFlag = tReaction.find( "(+M)") < tReaction.length()
                  || tReaction.find( "+M")  < tReaction.length() ;

                if( tThirdBodyFlag )
                {
                    tReaction = search_and_replace( tReaction, "(+M)", "+M");
                }
                tReaction = search_and_replace( tReaction, "(+AR)", "+AR");
                tReaction = search_and_replace( tReaction, "(+O2)", "+O2");
                tReaction = search_and_replace( tReaction, "(+N2)", "+N2");
                tReaction = search_and_replace( tReaction, "(+H2O)", "+H2O");

                // read coefficients
                Vector< real > & tCoeffs = tEntry->coeffs();
                tCoeffs.set_size( 3 );
                tCoeffs( 0 ) = std::stod( tWords( tWords.size() - 3 ) );
                tCoeffs( 1 ) = std::stod( tWords( tWords.size() - 2 ) );
                tCoeffs( 2 ) = std::stod( tWords( tWords.size() - 1 ) );

                // test if this reaction has more lines
                for( uint j=tIndex+1; j<mIndex( k+1 ); ++j )
                {
                    // get line from buffer
                    string tLine = string_to_upper( mBuffer( j ) );

                    // test if this line is a LOW flag
                    if( tLine.find("LOW") < tLine.length() )
                    {
                        Vector< real > & tLow = tEntry->low();

                        // tidy up line
                        tLine = search_and_replace( tLine, "LOW", "" );
                        tLine = search_and_replace( tLine, "/", "" );

                        tWords = string_to_words( tLine );

                        // reserve memory
                        tLow.set_size( 3 );

                        for ( uint i = 0; i < 3; ++i )
                        {
                            tLow( i ) = std::stod( tWords( i ) );
                        }
                    }
                    else if ( tLine.find("TROE") < tLine.length() )
                    {
                        Vector< real > & tTroe = tEntry->troe();

                        // tidy up line
                        tLine = search_and_replace( tLine, "TROE", "" );
                        tLine = search_and_replace( tLine, "/", "" );

                        // create list of words
                        tWords = string_to_words( tLine );

                        // reserve memory
                        tTroe.set_size( 4, BELFEM_QUIET_NAN );

                        for ( uint i = 0; i < tWords.size(); ++i )
                        {
                            tTroe( i ) =  std::stod( tWords( i ) );
                        }
                    }
                    else if ( tLine.find("DUPLICATE") < tLine.length() )
                    {
                        tEntry->set_duplicate_flag();
                    }
                    else if ( tLine.length() > 0 )
                    {
                        if( tLine.at( 0 ) != '!' && tLine.find( "/") < tLine.length() )
                        {

                            Cell< string > & tThirdBodySpecies = tEntry->third_body_species();
                            Vector< real > & tThirdBodyWeights = tEntry->third_body_weights();

                            // tidy up line
                            tLine = search_and_replace( tLine, "/", " " );



                            tWords = string_to_words( tLine );

                            // get number of third body reactions
                            uint tN =  0.5 * tWords.size() ;

                            // allocate memory
                            tThirdBodySpecies.set_size(tN+1, "");
                            tThirdBodyWeights.set_size( tN+1 );

                            uint tCount = 0;
                            for( uint i=0; i<tN; ++i )
                            {
                                tThirdBodySpecies( i ) = gastables::fix_label( tWords( tCount ++ ) );
                                tThirdBodyWeights( i ) = std::stod( tWords( tCount ++ ));
                            }
                            tThirdBodySpecies( tN ) = "N2";
                            tThirdBodyWeights( tN ) = 1.0 ;
                        }
                    }
                }

                if( tThirdBodyFlag && tEntry->third_body_weights().length() == 0 )
                {
                    Cell< string > & tThirdBodySpecies = tEntry->third_body_species();
                    Vector< real > & tThirdBodyWeights = tEntry->third_body_weights();
                    tThirdBodySpecies.set_size( 1, "N2" );
                    tThirdBodyWeights.set_size( 1, 1.0 );

                }
                mEntries( k ) = tEntry;
            }
        }

//------------------------------------------------------------------------------

        void
        Chemkin::find_duplicates()
        {
            // Map containing entries
            Map< string, uint > tMap;

            uint tIndex = 0;

            for ( Entry * tEntry: mEntries )
            {
                // test if this is a duplicate
                if( tEntry->is_duplicate() )
                {
                    // test if entry exists
                    if ( tMap.key_exists( tEntry->reaction() ) )
                    {
                        // get ref to entry
                        Entry * tOriginal = mEntries( tMap( tEntry->reaction() ) );

                        // link original to duplicate
                        tOriginal->set_duplicate( tEntry->coeffs() );

                        // deactivate duplicate
                        tEntry->deactivate();

                        BELFEM_ERROR(
                                tOriginal->low().length() == 0 &&
                                tOriginal->troe().length() == 0 &&
                                tOriginal->third_body_weights().length() == 0,
                                "Unsupported type of duplicate entry: %s",
                                tOriginal->reaction().c_str() );

                        BELFEM_ERROR(
                                tEntry->low().length() == 0 &&
                                tEntry->troe().length() == 0 &&
                                tEntry->third_body_weights().length() == 0,
                                "Unsupported type of duplicate entry: %s",
                                tEntry->reaction().c_str() );
                    }
                    else
                    {
                        // create entry in map
                        tMap[ tEntry->reaction() ] = tIndex;
                    }
                }

                ++tIndex;
            }
        }

//------------------------------------------------------------------------------

        void
        Chemkin::count_active_reactions()
        {
            mNumberOfReactions = 0;
            for ( Entry * tEntry: mEntries )
            {
                if( tEntry->is_active() )
                {
                    ++mNumberOfReactions;
                }
            }
        }


//------------------------------------------------------------------------------

        void
        Chemkin::get_species( Cell< string > & aSpecies )
        {
            aSpecies.clear();

            Cell< string > tEductLabels;
            Vector< real > tEductMolars;
            Cell< string > tProductLabels;
            Vector< real > tProductMolars;
            bool           tThirdBodyFlag;

            for ( Entry * tEntry: mEntries )
            {
                // split string of entry into species
                molars_from_string(
                        tEntry->reaction(),
                        tEductLabels,
                        tEductMolars,
                        tProductLabels,
                        tProductMolars,
                        tThirdBodyFlag );

                // append educts to list
                for( string tEduct : tEductLabels )
                {
                    aSpecies.push( tEduct );
                }

                // append products to list
                for( string tProduct : tProductLabels )
                {
                    aSpecies.push( tProduct );
                }
            }

            // make list unique
            unique( aSpecies );
        }

//------------------------------------------------------------------------------
    }
}