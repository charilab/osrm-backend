#ifndef GEOTIFF_SOURCE_HPP
#define GEOTIFF_SOURCE_HPP

#include "util/coordinate.hpp"
#include "util/exception.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <boost/assert.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_int.hpp>
#include <storage/io.hpp>

#include <iterator>
#include <unordered_map>

//#define _RASTER_GEOTIFF_FLOAT 1

#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>

namespace osrm
{
namespace extractor
{

/**
    \brief Small wrapper around raster source queries to optionally provide results
    gracefully, depending on source bounds
*/
struct GeoTiffDatum
{
    static std::int32_t get_invalid() { return std::numeric_limits<std::int32_t>::max(); }

    std::int32_t datum = get_invalid();

    GeoTiffDatum() = default;

    GeoTiffDatum(std::int32_t _datum) : datum(_datum) {}
};

class GeoTiffGrid
{
  public:
    GeoTiffGrid(GDALDataset* poDataset, std::size_t _xdim, std::size_t _ydim)
    {
        xdim = _xdim;
        ydim = _ydim;
        _data.reserve(ydim * xdim);

        GDALRasterBand  *poBand = poDataset->GetRasterBand(1);
#ifdef _RASTER_GEOTIFF_FLOAT
        float *buf = (float*) CPLMalloc(sizeof(float)*xdim);
        for(std::size_t y=0; y<ydim; y++) {
            poBand->RasterIO(GF_Read, 0, y, xdim, 1, buf, xdim, 1, GDT_Float32, 0, 0 );
            for(std::size_t x=0; x<xdim; x++) {
                _data.push_back((std::int32_t)buf[x]);
            }
        }
        CPLFree(buf);
#else
        std::uint8_t *buf = (std::uint8_t*) CPLMalloc(sizeof(std::uint8_t)*xdim);
        for(std::size_t y=0; y<ydim; y++) {
            poBand->RasterIO(GF_Read, 0, y, xdim, 1,
                            buf, xdim, 1, GDT_Byte, 0, 0 );
            for(unsigned int x=0; x<xdim; x++) {
                _data.push_back((std::int32_t)buf[x]);
            }
        }
        CPLFree(buf);
#endif
    }

    GeoTiffGrid(const GeoTiffGrid &) = default;
    GeoTiffGrid &operator=(const GeoTiffGrid &) = default;

    GeoTiffGrid(GeoTiffGrid &&) = default;
    GeoTiffGrid &operator=(GeoTiffGrid &&) = default;

    std::int32_t operator()(std::size_t x, std::size_t y) { return _data[y * xdim + x]; }
    std::int32_t operator()(std::size_t x, std::size_t y) const { return _data[(y)*xdim + (x)]; }

  private:
    std::vector<std::int32_t> _data;
    std::size_t xdim, ydim;
};

/**
    \brief Stores raster source data in memory and provides lookup functions.
*/
class GeoTiffSource
{
  private:
    const float xstep;
    const float ystep;

    float CalcSize(int min, int max, std::size_t count) const;

  public:
    GeoTiffGrid raster_data;

    const std::size_t width;
    const std::size_t height;
    const int xmin;
    const int xmax;
    const int ymin;
    const int ymax;

    GeoTiffDatum GetGeoTiffData(const int lon, const int lat) const;

    GeoTiffDatum GetGeoTiffInterpolate(const int lon, const int lat) const;

    GeoTiffSource(GeoTiffGrid _raster_data,
                 std::size_t width,
                 std::size_t height,
                 int _xmin,
                 int _xmax,
                 int _ymin,
                 int _ymax);
};

class GeoTiffContainer
{
  public:
    GeoTiffContainer() = default;

    int LoadGeoTiffSource(const std::string &path_string);

    GeoTiffDatum GetGeoTiffDataFromSource(unsigned int source_id, double lon, double lat);

    GeoTiffDatum GetGeoTiffInterpolateFromSource(unsigned int source_id, double lon, double lat);

  private:
    std::vector<GeoTiffSource> LoadedSources;
    std::unordered_map<std::string, int> LoadedSourcePaths;
};
}
}

#endif /* GEOTIFF_SOURCE_HPP */
