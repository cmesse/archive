//
// Created by Christian Messe on 21.11.20.
//

#include "cl_CH_Segment.hpp"
#include "cl_CH_Boundarylayer.hpp"
#include "fn_sum.hpp"
namespace belfem
{
    namespace channel
    {
//------------------------------------------------------------------------------

        Segment::Segment( const id_t & aID,
                          const real & aX,
                          const real & aA,
                          const real & aU,
                          const uint & aNumWalls ) :
                mID( aID ),
                mNumWalls( aNumWalls )
        {
            mData.set_size( 24, 0.0 );
            mData( BELFEM_CHANNEL_X ) = aX;
            mData( BELFEM_CHANNEL_A ) = aA;
            mData( BELFEM_CHANNEL_DH ) = 4.0 * aA / aU;
            mData( BELFEM_CHANNEL_TW1 ) = 300.0;
            //mData( BELFEM_CHANNEL_TW2 ) = 300.0;

            mWalls.set_size( mNumWalls, nullptr );
        }

//------------------------------------------------------------------------------

        Segment::~Segment()
        {
            // delete the elements
            for ( Wall * tWall : mWalls )
            {
                delete tWall;
            }
        }

//------------------------------------------------------------------------------

        void
        Segment::print()
        {
            //std::cout << mID << " " << this->x() << " " << this->cross_section() << " " << this->perimeter() << " "
            //          << mWalls( 0 )->average_surface_temperature() << std::endl;

            //const real & tT = mData( BELFEM_CHANNEL_TM ) ;
            //const real & tP = mData( BELFEM_CHANNEL_PM );
            //const real & tMa = mData( BELFEM_CHANNEL_MAM );
            //const real & tR  = mData( BELFEM_CHANNEL_RM ) ;

            //real tRho = tP / ( tR * tT );
            //real tDotM = tRho * mData( BELFEM_CHANNEL_UM ) * this->cross_section() ;

            std::cout << mID << " | " << this->x()
            << "  " << mData( BELFEM_CHANNEL_A ) *1e6
            << "  " << mData( BELFEM_CHANNEL_DH ) *1e3
            << " | " << mData( BELFEM_CHANNEL_TM )
            << " " << mData( BELFEM_CHANNEL_PM ) *1e-5
            << " " << mData( BELFEM_CHANNEL_MAM )
            << " | " << mData( BELFEM_CHANNEL_RM )
            << "  " << mData( BELFEM_CHANNEL_SM )
            << " | " << mData( BELFEM_CHANNEL_TAUW )
            << "   " << mData( BELFEM_CHANNEL_DOTQ ) *1e-6
            << "   " << mData( BELFEM_CHANNEL_ALPHA1 )
            << std::endl ;
        }

//------------------------------------------------------------------------------

        void
        Segment::add_wall( const uint & aWallIndex, Wall * aWall )
        {
            mWalls( aWallIndex ) = aWall;
        }

//------------------------------------------------------------------------------

        void
        Segment::pull_surface_temperatures()
        {
            // push first wall
            mData( BELFEM_CHANNEL_TW1 ) = mWalls( 0 )->average_surface_temperature() ;

            if( mNumWalls == 2 )
            {
                mData( BELFEM_CHANNEL_TW2 ) = mWalls( 1 )->average_surface_temperature() ;
            }
        }

//------------------------------------------------------------------------------

        void
        Segment::push_heatloads()
        {

            real tDotQ =  mWalls( 0 )->average_heatload(
                    mData( BELFEM_CHANNEL_ALPHA1 ),
                                        mData( BELFEM_CHANNEL_TREC ) );

            if( mNumWalls == 2 )
            {
                tDotQ *= mWalls( 0 )->segment_length() ;

                tDotQ += mWalls( 1 )->average_heatload(
                        mData( BELFEM_CHANNEL_ALPHA2 ),
                                            mData( BELFEM_CHANNEL_TREC ) )
                        * mWalls( 1 )->segment_length() ;

                // average heatload
                tDotQ /=  mWalls( 0 )->segment_length()
                        + mWalls( 1 )->segment_length() ;
            }

            // set heatloads in container
            mData( BELFEM_CHANNEL_DOTQ ) = tDotQ ;

        }

//------------------------------------------------------------------------------

        void
        Segment::push_flowdata()
        {
            mWalls( 0 )->set_flowdata(
                    mData( BELFEM_CHANNEL_TM ),
                    mData( BELFEM_CHANNEL_PM ),
                    mData( BELFEM_CHANNEL_MAM ) );

            if( mNumWalls == 2 )
            {
                mWalls( 1 )->set_flowdata( mData( BELFEM_CHANNEL_TM ),
                                           mData( BELFEM_CHANNEL_PM ),
                                           mData( BELFEM_CHANNEL_MAM ) );
            }
        }

//------------------------------------------------------------------------------
    }
}