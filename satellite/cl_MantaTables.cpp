//
// Created by christian on 9/14/24.
//
#include "commtools.hpp"
#include "cl_MantaTables.hpp"

#include <petsctools.hpp>

#include "cl_HDF5.hpp"
#include "cl_Progressbar.hpp"
#include "constants.hpp"
#include "cl_Solver.hpp"
#include "filetools.hpp"
#include "fn_unique.hpp"
#include "fn_preferred_matrix_format.hpp"
namespace belfem
{


    MantaTables::MantaTables( Mesh  & aMesh ) :
            mRank( comm_rank()) ,
            mCommSize( comm_size() ),
            mMesh( aMesh )
    {
        mSolver = new Solver( SolverType::MUMPS ) ;

        if ( mRank == 0 )
        {
            mMesh.create_field( "Infrared", EntityType::ELEMENT ) ;
            mMesh.create_field( "Visible", EntityType::ELEMENT ) ;
            mMesh.create_field( "Solar", EntityType::ELEMENT ) ;
            mMesh.create_field( "Wetness", EntityType::ELEMENT ) ;
        }
    }

    MantaTables::~MantaTables()
    {
        for( MantaSurface * tSurface : mSurfaces )
        {
            delete tSurface ;
        }

        for (SpMatrix *V: mViewFactors)
        {
            delete V;
        }

        if( mRadiationMatrix != nullptr )
        {
            delete mRadiationMatrix ;
        }
        if( mSolver != nullptr )
        {
            delete mSolver ;
        }
    }

    void
    MantaTables::interpolate_geometry_info( const real aTime )
    {
        if ( mRank == 0 )
        {
            index_t a ;
            index_t b ;
            real xi ;
            real eta ;
            this->compute_interpolation_factors( aTime, a, b, xi, eta );

            const Matrix< real > & tNodeCoordsI = mNodeCoords( a ) ;
            const Matrix< real > & tNodeCoordsJ = mNodeCoords( b ) ;

            Vector< real > tX( 3 );

            for ( mesh::Node * tNode : mMesh.nodes() )
            {
                // interpolate coordinates
                for ( uint i=0; i<3; ++i )
                {
                    tX(i) =   xi * tNodeCoordsI( tNode->index(), i )
                           + eta * tNodeCoordsJ( tNode->index(), i ) ;
                }

                // write coordinates to node
                tNode->set_coords( tX );
            }

            const Vector< real > & tFractionsI = mSolarAreaFractions( a );
            const Vector< real > & tFractionsJ = mSolarAreaFractions( b );

            const Vector< real > & tVelocitiesI = mVelocityAreaFractions( a );
            const Vector< real > & tVelocitiesJ = mVelocityAreaFractions( b );


            index_t tNumSurfaces = mNumSpacecraftElements + mNumPlanetElements ;
            Vector< real > & tSolar = mMesh.field_data("Solar");
            Vector< real > & tVelocity = mMesh.field_data("Wetness");

            for ( index_t e=0; e<mNumSpacecraftElements; ++e )
            {
                real tF = xi * tFractionsI( e ) + eta * tFractionsJ( e ) ;
                real tV = xi * tVelocitiesI( e ) + eta * tVelocitiesJ( e ) ;
                id_t tID = mElementIDs( e ) ;

                if ( std::abs( tF ) < 1e-6 )
                {
                    mSurfaceMap( tID )->set_solar_fraction( 0.0 );  // positive sides
                    mSurfaceMap( tID+tNumSurfaces )->set_solar_fraction( 0.0 ); //negative sides
                }
                else if ( mSurfaceMap( tID  )->is_spacecraft() )
                {
                    mSurfaceMap( tID )->set_solar_fraction( 0.0 );
                    mSurfaceMap( tID + tNumSurfaces )->set_solar_fraction( -tF );
                }
                else
                {
                    mSurfaceMap( tID )->set_solar_fraction( 0.0 );
                    if ( tID + tNumSurfaces < mMaxID )
                    {
                        mSurfaceMap( tID + tNumSurfaces )->set_solar_fraction( 0.0 );
                    }
                }

                if ( tID < mMaxSpacecraftElementID )
                {
                    tSolar( mMesh.element( tID )->index() ) = std::abs( tF );
                    tVelocity( mMesh.element( tID )->index() ) = std::abs( tV );
                }

            }

        }
    }

    void
    MantaTables::solve_infrared( const real aTime )
    {
        if( mRank == 0 )
        {
            // obtain the factors
            index_t n = mSurfaces.size() ;

            // reset the matrix
            mRadiationMatrix->fill( 0.0 );
            mRadiationRHS.set_size( n, -1.0 );
            mRadiationLHS.fill( 0.0 );

            SpMatrix & M = *mRadiationMatrix ;

            index_t a ;
            index_t b ;
            real xi ;
            real eta ;
            this->compute_interpolation_factors( aTime, a, b, xi, eta );

            SpMatrix * VFA = mViewFactors( a ) ;
            SpMatrix * VFB = mViewFactors( b ) ;


            for( index_t i=0; i<n; ++i )
            {
                M( i, i ) = 1.0 ;
            }

            index_t nnz = VFA->number_of_nonzeros() ;

            for ( index_t k=0; k<nnz; ++k )
            {
                index_t i = VFA->rows()[k];
                index_t j = VFA->cols()[k];

                real F_ij = xi * VFA->data()[k] + eta * VFB->data()[k] ;

                // write value into matrix
                M( i, j ) += F_ij * ( mSurfaces( j )->emmisivity() - 1.0 );

                // add contribution to RHS
                mRadiationRHS( j ) += F_ij ;
            }

            // finalize RHS
            for( MantaSurface * tSurface : mSurfaces )
            {
                mRadiationRHS( tSurface->index() ) *= tSurface->emmisivity() * constant::sigma *
                        tSurface->temperature() * tSurface->temperature()
                       * tSurface->temperature() * tSurface->temperature() ;
            }

        } // parallel

        comm_barrier();

        mSolver->solve( *mRadiationMatrix, mRadiationLHS, mRadiationRHS );

        if( mRank == 0 )
        {
            Vector< real > & tQ = mMesh.field_data( "Infrared" );

            tQ.fill( 0.0 );


            for( MantaSurface * tSurface : mSurfaces )
            {
                if ( tSurface->element() != nullptr )
                {
                    tQ( tSurface->element()->index() ) += mRadiationLHS( tSurface->index() );
                }
            }
        }
    }



    void
    MantaTables::solve_solar( const real aTime )
      {
        if( mRank == 0 )
        {
            for( MantaSurface * tSurface : mSurfaces )
            {
                tSurface->compute_solar( mSolarHeatflux );
            }

            // obtain the factors
            index_t n = mSurfaces.size() ;

            // reset the matrix
            mRadiationMatrix->fill( 0.0 );
            mRadiationRHS.set_size( n, -1.0 );
            mRadiationLHS.fill( 0.0 );

            SpMatrix & M = *mRadiationMatrix ;

            index_t a ;
            index_t b ;
            real xi ;
            real eta ;
            this->compute_interpolation_factors( aTime, a, b, xi, eta );

            SpMatrix * VFA = mViewFactors( a ) ;
            SpMatrix * VFB = mViewFactors( b ) ;


            for( index_t i=0; i<n; ++i )
            {
                M( i, i ) = 1.0 ;
            }

            index_t nnz = VFA->number_of_nonzeros() ;

            for ( index_t k=0; k<nnz; ++k )
            {
                index_t i = VFA->rows()[k];
                index_t j = VFA->cols()[k];

                real F_ij = xi * VFA->data()[k] + eta * VFB->data()[k] ;

                // write value into matrix
                M( i, j ) -= F_ij * mSurfaces( j )->absorbtivity();

                // add contribution to RHS
                mRadiationRHS( j ) += F_ij ;
            }

            // finalize RHS
            for( MantaSurface * tSurface : mSurfaces )
            {
                mRadiationRHS( tSurface->index() ) *= tSurface->solar_reflection() ;
            }

        } // end parralel

        comm_barrier();

        mSolver->solve( *mRadiationMatrix, mRadiationLHS, mRadiationRHS );

        if( mRank == 0 )
        {
            Vector< real > & tQ = mMesh.field_data( "Visible" );

            tQ.fill( 0.0 );


            for( MantaSurface * tSurface : mSurfaces )
            {
                if ( tSurface->element() != nullptr )
                {
                    tQ( tSurface->element()->index() ) += mRadiationLHS( tSurface->index() ) + tSurface->solar_absorption() ;
                }
            }
        }
    }

    void
    MantaTables::load_database( const string & aPath )
    {
        if ( mRank == 0 )
        {
            HDF5 tDatabase( aPath, FileMode::OPEN_RDONLY );

            // read header
            tDatabase.load_data( "maxTime", mMaxTime );
            tDatabase.load_data( "numKeyframes", mNumKeyframes );

            // allocate memory
            mKeyframeTimesteps.set_size( mNumKeyframes );
            mViewFactors.set_size( mNumKeyframes, nullptr );
            mSolarAreaFractions.set_size( mNumKeyframes, {} );
            mVelocityAreaFractions.set_size( mNumKeyframes, {} );
            mNodeCoords.set_size( mNumKeyframes, {} );


            // read surfaces
            tDatabase.select_group( "Mesh");
            tDatabase.load_data( "ElementAreas", mSpacecraftSurfaceAreas );
            tDatabase.load_data( "MaxSpacecraftID", mMaxSpacecraftElementID );

            Vector< unsigned int > tSpacecraftIDs ;
            tDatabase.load_data( "ElementIDs", tSpacecraftIDs );

            Vector< unsigned int > tPlanetIDs ;
            tDatabase.load_data( "PlanetIDs", tPlanetIDs );

            tDatabase.close_active_group();

            mElementIDs.set_size( tSpacecraftIDs.length() + tPlanetIDs.length() ) ;

            index_t tCount = 0 ;
            Map< id_t, index_t > tIndexMap ;

            for ( id_t tID : tSpacecraftIDs )
            {
                tIndexMap[ tID ] = tCount ;
                mElementIDs( tCount++ ) = tID ;
            }
            for ( id_t tID : tPlanetIDs )
            {
                tIndexMap[ tID ] = tCount ;
                mElementIDs( tCount++ ) = tID ;
            }

            for (uint k = 0; k < mNumKeyframes; ++k)
            {
                string tLabel = sprint( "keyframe_%02u", k );

                tDatabase.select_group( tLabel );

                // load timestep
                tDatabase.load_data( "time", mKeyframeTimesteps( k ) );

                // load the node coordinates
                tDatabase.load_data("nodeCoords", mNodeCoords( k ) );

                tDatabase.select_group( "Areas" );

                tDatabase.load_data( "SolarFactor", mSolarAreaFractions( k ) );

                // we don't need this now, this is needed for drag etc
                tDatabase.load_data( "VelocityFactor", mVelocityAreaFractions( k ) );

                //tDatabase.load_data( "Planet", mAllPlanetSurfaceAreas( k ) );

                // close areas
                tDatabase.close_active_group();

                // close keyframe
                tDatabase.close_active_group();
            }


            // load the environment data
            tDatabase.select_group( "Environment");
            tDatabase.load_data( "Time", mEnvironmentTimesteps );
            tDatabase.load_data( "SolarHeatFlux", mEnvironmentSolarHeatflux );
            tDatabase.load_data( "PlanetAlbedo", mEnvironmentPlanetAlbedo );
            tDatabase.load_data( "PlanetTemperature", mEnvironmentPlanetTemperature );
            tDatabase.close_active_group() ;

            mNumSpacecraftElements = tSpacecraftIDs.length() ;
            mNumPlanetElements = tPlanetIDs.length();


            this->create_surfaces() ;



            // load the view factors
            tDatabase.select_group( "Mesh");



            id_t tMaxSpacecraftID ;
            tDatabase.load_data( "MaxSpacecraftID", tMaxSpacecraftID );

            id_t tMaxPlanetID ;
            tDatabase.load_data( "MaxPlanetID", tMaxPlanetID );


            tDatabase.close_active_group() ;



            tCount = 0 ;

            // set IDs for positive sides on spacecraft
            for( id_t tID : tSpacecraftIDs )
            {
                mSurfaces( tCount++ )->set_id( tID );
            }

            // indices for planet
            for( id_t tID : tPlanetIDs )
            {
                mSurfaces( tCount++ )->set_id( tID );
            }

            // indices for negative sides on spacecraft
            for( id_t tID : tSpacecraftIDs )
            {
                mSurfaces( tCount++ )->set_id( tID + mNumSpacecraftElements + mNumPlanetElements );
            }

            BELFEM_ASSERT( tCount == mSurfaces.size(), "Error in reading view factors" );

            // create the map
            tCount = 0 ;

            mMaxID = 0 ;

            for ( MantaSurface * tSurface : mSurfaces )
            {
                // add the surface to the map
                mSurfaceMap[ tSurface->id() ] = tSurface ;

                // make sure that the index is correct (probably not neccessary)
                tSurface->set_index( tCount++ );

                if ( tSurface->id() > mMaxID )
                {
                    mMaxID = tSurface->id() ;
                }
            }


            // load the keyframe data
            Cell< Vector< unsigned int > > tSurfacesI( mNumKeyframes, {} );
            Cell< Vector< unsigned int > > tSurfacesJ( mNumKeyframes, {} );
            Cell< Vector< real > > tViewFactorsIJ( mNumKeyframes, {} );
            Cell< Vector< real > > tViewFactorsJI( mNumKeyframes, {} );

            for (uint k = 0; k < mNumKeyframes; ++k)
            {
                string tLabel = sprint( "keyframe_%02u", k );
                tDatabase.select_group( tLabel );
                tDatabase.select_group( "ViewFactors");


                tDatabase.load_data( "SurfacesI",tSurfacesI( k ) );
                tDatabase.load_data( "SurfacesJ",tSurfacesJ( k ) );
                tDatabase.load_data( "ViewFactorsIJ",tViewFactorsIJ( k ) );
                tDatabase.load_data( "ViewFactorsJI",tViewFactorsJI( k ) );

                // close view factors
                tDatabase.close_active_group();

                // close keyframe
                tDatabase.close_active_group();
            } // end loop over keyframes


            // count the surfaces
            Vector< index_t > tCounters( mSurfaces.size(), 1 );

            for ( uint k = 0; k < mNumKeyframes; ++k )
            {
                // populate data
                const Vector< unsigned int > & tSi  = tSurfacesI( k );
                const Vector< unsigned int > & tSj  = tSurfacesJ( k );
                const Vector< real > & tFij = tViewFactorsIJ( k );
                const Vector< real > & tFji = tViewFactorsJI( k );

                index_t tNumValues = tSi.length() ;

                for ( index_t l = 0; l < tNumValues; ++l )
                {
                    /*id_t tID = 714 ;
                    if ( tSi( l ) == tID || tSj( l ) == tID )
                    {
                        std::cout << "#CHECK " << tSi( l ) << " " << tSj( l ) << " " << tFij( l ) << " " << tFji( l ) << std::endl;
                    }*/
                    // skip negative sides of planet
                    if ( tSi( l ) > mMaxID || tSj( l ) > mMaxID )
                    {
                        continue;
                    }
                    // check counter for senders I
                    if ( tFij( l ) > mMinVF )
                    {
                        ++tCounters( mSurfaceMap( tSi( l ) )->index() );
                    }

                    // check counter for senders J
                    if ( tFji( l ) > mMinVF )
                    {
                        ++tCounters( mSurfaceMap( tSj( l ) )->index() );
                    }
                }
            }

            // with the counters set, we allocate temporary vectors with the ids
            Cell< Vector< index_t > > tTargets( tCount, {} );

            for ( index_t s=0; s < tCount; ++s )
            {
                tTargets( s ).set_size( tCounters( s ), 0 );
            }

            // add self
            for ( index_t s=0; s < tCount; ++s )
            {
                tTargets( s )( 0 ) = s ;
            }

            tCounters.fill( 1 );

            for ( uint k = 0; k < mNumKeyframes; ++k )
            {
                // populate data
                const Vector< unsigned int > & tSi  = tSurfacesI( k );
                const Vector< unsigned int > & tSj  = tSurfacesJ( k );
                const Vector< real > & tFij = tViewFactorsIJ( k );
                const Vector< real > & tFji = tViewFactorsJI( k );


                index_t tNumValues = tSi.length() ;

                for ( index_t l = 0; l < tNumValues; ++l )
                {
                    // skip negative sides of planet
                    if ( tSi( l ) > mMaxID || tSj( l ) > mMaxID )
                    {
                        continue;
                    }

                    index_t i = mSurfaceMap( tSi( l ) )->index() ;
                    index_t j = mSurfaceMap( tSj( l ) )->index() ;

                    // check counter for senders I
                    if ( tFij( l ) > mMinVF )
                    {
                        tTargets( i )( tCounters( i )++ ) = j ;
                    }

                    // check counter for senders J
                    if ( tFji( l ) > mMinVF )
                    {
                        tTargets( j )( tCounters( j )++ ) = i ;
                    }
                }
            }

            // create the graph
            Cell< graph::Vertex * > tGraph( tCount, nullptr );
            for ( index_t i=0; i < tCount; ++i )
            {
                Vector< index_t > & tTarget = tTargets( i );
                unique( tTarget );
                mSurfaces( i )->init_vertex_container( tTarget.length() );
                for ( index_t j : tTarget )
                {
                    mSurfaces( i )->insert_vertex( mSurfaces( j ) );
                }
                tGraph( i ) = mSurfaces( i );
            }

            // with the graph, we can create the matrices
            mRadiationMatrix = new SpMatrix( tGraph, preferred_matrix_format( mSolver->type() ), tCount, tCount );

            mViewFactors.set_size( mNumKeyframes, nullptr );

            for ( uint k = 0; k < mNumKeyframes; ++k )
            {
                // create matrix
                mViewFactors( k ) = new SpMatrix( tGraph, SpMatrixType::CSR, tCount, tCount );

                SpMatrix & tVF = *mViewFactors( k );
                tVF.create_coo_indices() ;

                // populate data
                const Vector< unsigned int > & tSi  = tSurfacesI( k );
                const Vector< unsigned int > & tSj  = tSurfacesJ( k );
                const Vector< real > & tFij = tViewFactorsIJ( k );
                const Vector< real > & tFji = tViewFactorsJI( k );


                index_t tNumValues = tSi.length() ;

                for ( index_t l = 0; l < tNumValues; ++l )
                {
                    // skip negative sides of planet
                    if ( tSi( l ) > mMaxID || tSj( l ) > mMaxID )
                    {
                        continue;
                    }

                    index_t i = mSurfaceMap( tSi( l ) )->index() ;
                    index_t j = mSurfaceMap( tSj( l ) )->index() ;

                    if ( tFij( l ) > mMinVF )
                    {
                        tVF( i, j ) = tFij( l ) ;
                    }

                    // check counter for senders J
                    if ( tFji( l ) > mMinVF )
                    {
                        tVF( j, i ) = tFji( l ) ;
                    }
                }
            }
        } // end rank = 0




    }

    void
    MantaTables::create_surfaces()
    {

        mSurfaces.set_size( 2 * mNumSpacecraftElements + mNumPlanetElements, nullptr );
        mPlanetSurfaces.set_size( mNumPlanetElements, nullptr );

		index_t tCount = 0  ;

        // create field on mesh
        Vector< real > & tT = mMesh.create_field( "Temperature", EntityType::ELEMENT ) ;

        tT.fill( mTinit );

        // element container on mesh
        Cell< mesh::Element * > & tElements = mMesh.elements();

        // create positive sides
        for( index_t k=0; k<mNumSpacecraftElements; ++k )
        {
            mSurfaces( tCount )= new MantaSurface( tCount, tElements( k ), tT( k ) );
            mSurfaces( tCount )->set_surface_area( mSpacecraftSurfaceAreas( k ) );
            ++tCount ;
        }

        // planet sides: in the future, we can link each planet surface to an individual temperature,
        // which would be useful for the moon
        for( index_t k=0; k<mNumPlanetElements; ++k )
        {
          	mSurfaces( tCount ) = new MantaSurface( tCount, nullptr, mPlanetTemperature );
            mPlanetSurfaces( k ) = mSurfaces( tCount );

            // for celestial bodies, we always refer to the black body temperature, henche epsilon=1
            mSurfaces( tCount )->set_emmisivity( 1.0 );
            ++tCount ;
        }

        // create the negative sides
        for( index_t k=0; k<mNumSpacecraftElements; ++k )
        {
          	mSurfaces( tCount ) = new MantaSurface( tCount, tElements( k ), tT( k ) );
            mSurfaces( tCount )->set_surface_area( mSpacecraftSurfaceAreas( k ) );
            ++tCount ;
        }
    }

    void
    MantaTables::create_communication_table()
    {
        if( mRank == 0 && mCommSize > 1 )
        {
            // create communication list
            uint c = 0 ;
            mCommTable.set_size( mCommSize-1 );
            for( proc_t k=0; k<mRank; ++k )
            {
                if( k != mRank )
                {
                    mCommTable( c++ ) = k ;
                }
            }
        }
    }

    void
    MantaTables::compute_environment( const real aTime )
    {
        // reduce time based on orbit time
        real tTime = aTime - std::floor( aTime / mMaxTime ) * mMaxTime ;
        real tDeltaTime =  mEnvironmentTimesteps(1)-mEnvironmentTimesteps(0) ;

        // identify interval
        index_t i =  static_cast<index_t>( std::floor( tTime / tDeltaTime ) );
        if ( i >= mEnvironmentTimesteps.length() )
        {
            i = mEnvironmentTimesteps.length() - 1 ;
        }

        // reduced time
        real tTau = ( tTime - mEnvironmentTimesteps( i ) ) / tDeltaTime ;

        // interpolate values
        mPlanetAlbedo      = mEnvironmentPlanetAlbedo( i ) * tTau + mEnvironmentPlanetAlbedo( i + 1 ) * ( 1.0 - tTau );
        mPlanetTemperature = mEnvironmentPlanetTemperature( i ) * tTau + mEnvironmentPlanetTemperature( i + 1 ) * ( 1.0 - tTau );
        mSolarHeatflux     = mEnvironmentSolarHeatflux( i ) * tTau + mEnvironmentSolarHeatflux( i + 1 ) * ( 1.0 - tTau );

        //mPlanetAlbedo = 0.0 ;
        //mPlanetTemperature = 0.0 ;

        // update albedo data
        for ( MantaSurface * tSurface : mPlanetSurfaces )
        {
            tSurface->set_absorbtivity( 1.0 - mPlanetAlbedo );
        }

    }

    void
    MantaTables::compute_interpolation_factors( const real aTime, index_t & aI, index_t & aJ, real & aXi, real & aEta )
    {
        real tTime = aTime - std::floor( aTime / mMaxTime ) * mMaxTime ;

        // find interval
        aI = 0 ;
        aJ = 0 ;
        real tTimeI  = 0 ;
        real tTimeJ  ;

        // look for left keyframe
        for( index_t k=0; k<mNumKeyframes; ++k )
        {
            if( mKeyframeTimesteps( k ) < tTime )
            {
                tTimeI = mKeyframeTimesteps( k );
                aI = k;
            }
            else
            {
                break ;
            }
        }

        // look for right keyframe
        if( aI == mNumKeyframes-1 )
        {
            aJ = 0 ;
            tTimeJ = mMaxTime ;
        }
        else
        {
            aJ = aI + 1 ;
            tTimeJ = mKeyframeTimesteps( aJ );
        }

        // compute interpolation factors
        aEta = ( tTime - tTimeI ) / ( tTimeJ - tTimeI );
        aXi = ( 1.0 - aEta ) ;
    }


}