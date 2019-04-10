#ifndef RASTER_SOURCE_HPP
#define RASTER_SOURCE_HPP

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

#define _RASTER_GEOTIFF 1
//#define _RASTER_GEOTIFF_FLOAT 1

#ifdef _RASTER_GEOTIFF
#include <gdal/gdal_priv.h>
#include <gdal/cpl_conv.h>
#endif

namespace osrm
{
namespace extractor
{

/**
    \brief Small wrapper around raster source queries to optionally provide results
    gracefully, depending on source bounds
*/
struct RasterDatum
{
    static std::int32_t get_invalid() { return std::numeric_limits<std::int32_t>::max(); }

    std::int32_t datum = get_invalid();

    RasterDatum() = default;

    RasterDatum(std::int32_t _datum) : datum(_datum) {}
};

class RasterGrid
{
  public:
#ifdef _RASTER_GEOTIFF
    RasterGrid(GDALDataset* poDataset, std::size_t _xdim, std::size_t _ydim)
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
#else
    RasterGrid(const boost::filesystem::path &filepath, std::size_t _xdim, std::size_t _ydim)
    {
        xdim = _xdim;
        ydim = _ydim;
        _data.reserve(ydim * xdim);

        storage::io::FileReader file_reader(filepath, storage::io::FileReader::HasNoFingerprint);

        std::string buffer;
        buffer.resize(file_reader.GetSize());

        BOOST_ASSERT(buffer.size() > 1);

        file_reader.ReadInto(&buffer[0], buffer.size());

        boost::algorithm::trim(buffer);

        auto itr = buffer.begin();
        auto end = buffer.end();

        bool r = false;
        try
        {
            r = boost::spirit::qi::parse(
                itr, end, +boost::spirit::qi::int_ % +boost::spirit::qi::space, _data);
        }
        catch (std::exception const &ex)
        {
            throw util::exception("Failed to read from raster source " + filepath.string() + ": " +
                                  ex.what() + SOURCE_REF);
        }

        if (!r || itr != end)
        {
            throw util::exception("Failed to parse raster source: " + filepath.string() +
                                  SOURCE_REF);
        }
    }
#endif

    RasterGrid(const RasterGrid &) = default;
    RasterGrid &operator=(const RasterGrid &) = default;

    RasterGrid(RasterGrid &&) = default;
    RasterGrid &operator=(RasterGrid &&) = default;

    std::int32_t operator()(std::size_t x, std::size_t y) { return _data[y * xdim + x]; }
    std::int32_t operator()(std::size_t x, std::size_t y) const { return _data[(y)*xdim + (x)]; }

  private:
    std::vector<std::int32_t> _data;
    std::size_t xdim, ydim;
};

/**
    \brief Stores raster source data in memory and provides lookup functions.
*/
class RasterSource
{
  private:
    const float xstep;
    const float ystep;

    float CalcSize(int min, int max, std::size_t count) const;

  public:
    RasterGrid raster_data;

    const std::size_t width;
    const std::size_t height;
    const int xmin;
    const int xmax;
    const int ymin;
    const int ymax;

    RasterDatum GetRasterData(const int lon, const int lat) const;

    RasterDatum GetRasterInterpolate(const int lon, const int lat) const;

    RasterSource(RasterGrid _raster_data,
                 std::size_t width,
                 std::size_t height,
                 int _xmin,
                 int _xmax,
                 int _ymin,
                 int _ymax);
};

class RasterContainer
{
  public:
    RasterContainer() = default;

    int LoadRasterSource(const std::string &path_string,
                         double xmin,
                         double xmax,
                         double ymin,
                         double ymax,
                         std::size_t nrows,
                         std::size_t ncols);

    RasterDatum GetRasterDataFromSource(unsigned int source_id, double lon, double lat);

    RasterDatum GetRasterInterpolateFromSource(unsigned int source_id, double lon, double lat);

  private:
    std::vector<RasterSource> LoadedSources;
    std::unordered_map<std::string, int> LoadedSourcePaths;
};
}
}

#endif /* RASTER_SOURCE_HPP */
