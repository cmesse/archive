//
// Created by Christian Messe on 25.09.19.
//
#include <cctype>

#include "fn_CN_molars_from_string.hpp"
#include "fn_GT_fix_label.hpp"
#include "stringtools.hpp"
#include "cl_Map.hpp"
namespace belfem
{
    namespace combustion
    {
//------------------------------------------------------------------------------

        void
        molars_from_string(
                const string  & aString,
                Cell <string> & aEductLabels,
                Vector <real> & aEductMolars,
                Cell <string> & aProductLabels,
                Vector <real> & aProductMolars,
                bool          & aHasThirdBody )
        {
            Cell <string>  tEductLabels;
            Vector <real>  tEductMolars;
            Cell <string>  tProductLabels;
            Vector <real>  tProductMolars;

            // remove brackets
            string tString = search_and_replace( aString, "<", "" );
            tString = search_and_replace( tString, ">", "" );

            // find tag
            uint tFlag = tString.find_first_of( "=" );

            string tEductString =  search_and_replace(
                    tString.substr( 0, tFlag ), "+", " " );

            string tProductString = search_and_replace(
                    tString.substr( tFlag+1, tString.size() - tFlag-1 ), "+", " " );

            // read educts
            Cell< string > tEductWords = string_to_words( tEductString );
            uint tN = tEductWords.size();


            // read products
            Cell< string > tProductWords = string_to_words( tProductString );
            uint tM = tProductWords.size();

            // Check if this reaction has third body entries
            aHasThirdBody = false;
            for( string tSpecie : tEductWords )
            {
                if ( tSpecie == "M" )
                {
                    aHasThirdBody = true;
                    break;
                }
            }

            if( ! aHasThirdBody )
            {
                for( string tSpecie : tProductWords )
                {
                    if ( tSpecie == "M" )
                    {
                        aHasThirdBody = true;
                        break;
                    }
                }
            }

            uint tCount;

            if( aHasThirdBody )
            {
                tCount = tN-1;
            }
            else
            {
                tCount = tN;
            }

            tEductLabels.set_size( tCount, "" );
            tEductMolars.set_size( tCount );
            tCount = 0;
            for( uint k=0; k<tN; ++k )
            {
                if( tEductWords( k ) != "M" )
                {
                    split_to_molar_and_label(
                            tEductWords( k ),
                            tEductLabels( tCount ),
                            tEductMolars( tCount ) );
                    ++tCount;
                }
            }


            if( aHasThirdBody )
            {
                tCount = tM-1;
            }
            else
            {
                tCount = tM;
            }

            tProductLabels.set_size( tCount, "" );
            tProductMolars.set_size( tCount );
            tCount = 0;
            for( uint k=0; k<tM; ++k )
            {
                if( tProductWords( k ) != "M")
                {
                    split_to_molar_and_label(
                            tProductWords( k ),
                            tProductLabels( tCount ),
                            tProductMolars( tCount ));
                    ++tCount;
                }
            }

            // make products unique
            tCount = 0 ;
            Map< string, uint > tEductMap ;

            for( string tLabel : tEductLabels )
            {
                if( ! tEductMap.key_exists( tLabel ) )
                {
                    tEductMap[ tLabel ] = tCount++ ;
                }
            }

            aEductLabels.set_size( tCount, "" );
            aEductMolars.set_size( tCount, 0 );

            tCount = 0 ;
            for( string tLabel : tEductLabels )
            {
                uint k = tEductMap[ tLabel ];
                aEductLabels( k ) = tLabel ;
                aEductMolars( k ) += tEductMolars( tCount++ ) ;
            }

            // make products unique
            tCount = 0 ;
            Map< string, uint > tProductMap ;

            for( string tLabel : tProductLabels )
            {
                if( ! tProductMap.key_exists( tLabel ) )
                {
                    tProductMap[ tLabel ] = tCount++ ;
                }
            }

            aProductLabels.set_size( tCount, "" );
            aProductMolars.set_size( tCount, 0 );

            tCount = 0 ;
            for( string tLabel : tProductLabels )
            {
                uint k = tProductMap[ tLabel ];
                aProductLabels( k ) = tLabel ;
                aProductMolars( k ) += tProductMolars( tCount++ ) ;
            }
        }

//------------------------------------------------------------------------------

        void
        split_to_molar_and_label(
                const string & aString,
                string & aLabel,
                real   & aMolar )
        {
            uint tN = aString.size();

            string tMolar = "";

            for( uint k=0; k<tN; ++k )
            {
                if( std::isdigit( aString.at( k ) ) == 1 ||  aString.at( k ) == '.' )
                {
                    tMolar += aString.at( k );
                }
                else
                {
                    aLabel = gastables::fix_label( aString.substr( k, tN-k ) );
                    break;
                }
            }

            if( tMolar.size() > 0 )
            {
                aMolar = std::stod( tMolar );
            }
            else
            {
                aMolar = 1.0;
            }
        }

//------------------------------------------------------------------------------
    }
}