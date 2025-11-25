//
// Created by christian on 9/14/24.
//

#ifndef BELFEM_CL_MANTASURFACE_HPP
#define BELFEM_CL_MANTASURFACE_HPP

#include "typedefs.hpp"
#include "cl_Element.hpp"
#include "cl_Graph_Vertex.hpp"

namespace belfem
{
    class MantaSurface : public graph::Vertex
    {

        mesh::Element * mElement ;
        real & mTemperature;
        real mSurfaceArea = 0.0 ;
        real mSolarCellEfficiency = 0.0 ;
        real mEpsilon = 0.5 ;
        real mAlpha   = 0.10 ;
        real mSolarFraction = 0.0 ;
        real mSolarReflection = 0.0 ;
        real mSolarAbsorption = 0.0 ;
    public:

        MantaSurface( const index_t aIndex, mesh::Element * aElement, real & aTemperature );

        ~MantaSurface();


        void
        set_surface_area( const real aSurfaceArea );

        void
        set_emmisivity( const real aEpsilon );

        void
        set_absorbtivity( const real aAlpha ) ;

        void
        set_solar_cell_efficiency( const real aEta );

        real
        emmisivity() const ;

        real
        absorbtivity() const ;

        real
        temperature() const ;

        real
        solar_cell_efficiency() const ;

        mesh::Element *
        element() ;

        real
        solar_area() const ;

        real
        area() const ;

        void
        set_temperature( const real aTemperature );

        void
        set_solar_fraction( const real aSolarFraction );

        real
        solar_fraction() const ;

        bool
        is_spacecraft() const ;

        void
        compute_solar( const real aSolarHeatflux );

        real
        solar_reflection() const ;

        real
        solar_absorption() const ;

    };

    inline void
    MantaSurface::set_emmisivity( const real aEpsilon )
    {
        mEpsilon = aEpsilon ;
    }

    inline void
    MantaSurface::set_absorbtivity( const real aAlpha )
    {
        mAlpha = aAlpha ;
    }

    inline void
    MantaSurface::set_surface_area( const real aSurfaceArea )
    {
        mSurfaceArea = aSurfaceArea ;
    }

    inline real
    MantaSurface::emmisivity() const
    {
        return mEpsilon ;
    }

    inline real
    MantaSurface::absorbtivity() const
    {
        return mAlpha;
    }

    inline real
    MantaSurface::temperature() const
    {
        return mTemperature;
    }

    inline void
    MantaSurface::set_temperature( const real aTemperature )
    {
        mTemperature = aTemperature ;
    }

    inline mesh::Element *
    MantaSurface::element()
    {
        return mElement ;
    }

    inline real
    MantaSurface::solar_cell_efficiency() const
    {
        return mSolarCellEfficiency ;
    }

    inline real
    MantaSurface::area() const
    {
      return mSurfaceArea ;
    }

    inline void
    MantaSurface::set_solar_fraction( const real aSolarFraction )
    {
      mSolarFraction = aSolarFraction ;
    }

    inline real
    MantaSurface::solar_fraction() const
    {
        return mSolarFraction ;
    }

    inline bool
    MantaSurface::is_spacecraft() const
    {
      return ( mElement != nullptr ) ;
    }

    inline void
    MantaSurface::compute_solar( const real aSolarHeatflux )
    {
        mSolarAbsorption = this->absorbtivity() * aSolarHeatflux * this->solar_fraction() ;
        mSolarReflection = ( 1.0 - this->absorbtivity() ) * aSolarHeatflux * this->solar_fraction() ;
    }

    inline real
    MantaSurface::solar_reflection() const
    {
        return mSolarReflection ;
    }

    inline real
    MantaSurface::solar_absorption() const
    {
        return mSolarAbsorption ;
    }

}
#endif //BELFEM_CL_MANTASURFACE_HPP
