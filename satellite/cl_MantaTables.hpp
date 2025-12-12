//
// Created by christian on 9/14/24.
//

#ifndef BELFEM_CL_MANTATABLES_HPP
#define BELFEM_CL_MANTATABLES_HPP
#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_SpMatrix.hpp"
#include "cl_Mesh.hpp"
#include "cl_MantaSurface.hpp"
#include "cl_Solver.hpp"

namespace belfem
{
    class MantaTables
    {
        const proc_t mCommRank;
        const proc_t mCommSize ;

        Mesh &mMesh;

        const string mMatrixFile = "radiation_matrix.hdf5";

        const real mMinVF = 1e-8 ;
        const real mTinit = 298.15 ;

        Vector< proc_t > mCommTable ;

        // Lookup table data
        Vector <real> mEnvironmentTimesteps;
        Vector <real> mEnvironmentPlanetAlbedo;
        Vector <real> mEnvironmentPlanetTemperature;
        Vector <real> mEnvironmentSolarHeatflux;

        Vector <real>   mSpacecraftSurfaceAreas ;

        // keyframe data
        Vector <real> mKeyframeTimesteps;
        Cell <Vector< real >>  mSolarAreaFractions;
        Cell <Vector< real > > mVelocityAreaFractions ;
		//Cell <Vector< real > > mAllPlanetSurfaceAreas ;
		Cell <Matrix< real > > mNodeCoords ;

        Cell< SpMatrix * > mViewFactors;

        // contains all surfaces, both spacecraft and planet
        Cell< MantaSurface * > mSurfaces ;
        Map< id_t, MantaSurface * > mSurfaceMap ;

        // contains only planet surfaces
        Cell< MantaSurface * > mPlanetSurfaces ;

        uint mNumKeyframes = 0;
        real mMaxTime = 0;

        real mPlanetAlbedo = 0.3;          // default value for earth
        real mPlanetTemperature = 261.15;  // default value for earth

        real mSolarHeatflux = 0;
        id_t mMaxSpacecraftElementID = 0 ;
        id_t mMaxID = 0 ;

        index_t mNumSpacecraftElements = 0 ;
        index_t mNumPlanetElements = 0 ;

        SpMatrix * mRadiationMatrix = nullptr  ;
        Vector< real > mRadiationLHS ;
        Vector< real > mRadiationRHS ;
        Vector< id_t > mElementIDs ;

        Solver * mSolver = nullptr ;
    public:

        MantaTables( Mesh & aMesh );

        ~MantaTables();

        void
        solve_infrared( const real aTime );

        void
        solve_solar( const real aTime );

        void
        load_database( const string &aPath );

        void
        create_surfaces();

        void
        compute_environment( const real aTime );

        void
        compute_interpolation_factors( const real aTime, index_t & aI, index_t & aJ, real & aXi, real & aEta );

        void
        interpolate_geometry_info( const real aTime );

    private :

        void
        create_communication_table();

        bool
        load_radiation_matrix();

    };
}

#endif //BELFEM_CL_MANTATABLES_HPP
