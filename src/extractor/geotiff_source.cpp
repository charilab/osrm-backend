#include "extractor/geotiff_source.hpp"

#include "util/exception.hpp"
#include "util/exception_utils.hpp"
#include "util/log.hpp"
#include "util/timing_util.hpp"
#include "util/typedefs.hpp"

#include <cmath>

namespace osrm
{
namespace extractor
{

GeoTiffSource::GeoTiffSource(GeoTiffGrid _raster_data,
                           std::size_t _width,
                           std::size_t _height,
                           int _xmin,
                           int _xmax,
                           int _ymin,
                           int _ymax)
    : xstep(CalcSize(_xmin, _xmax, _width)), ystep(CalcSize(_ymin, _ymax, _height)),
      raster_data(std::move(_raster_data)), width(_width), height(_height), xmin(_xmin),
      xmax(_xmax), ymin(_ymin), ymax(_ymax)
{
    BOOST_ASSERT(xstep != 0);
    BOOST_ASSERT(ystep != 0);
}

float GeoTiffSource::CalcSize(int min, int max, std::size_t count) const
{
    BOOST_ASSERT(count > 0);
    return (max - min) / (static_cast<float>(count) - 1);
}

// Query raster source for nearest data point
GeoTiffDatum GeoTiffSource::GetGeoTiffData(const int lon, const int lat) const
{
    if (lon < xmin || lon > xmax || lat < ymin || lat > ymax)
    {
        return {};
    }

    const std::size_t xth = static_cast<std::size_t>(round((lon - xmin) / xstep));
    const std::size_t yth = static_cast<std::size_t>(round((ymax - lat) / ystep));

    return {raster_data(xth, yth)};
}

// Query raster source using bilinear interpolation
GeoTiffDatum GeoTiffSource::GetGeoTiffInterpolate(const int lon, const int lat) const
{
    if (lon < xmin || lon > xmax || lat < ymin || lat > ymax)
    {
        return {};
    }

    const auto xthP = (lon - xmin) / xstep;
    const auto ythP =
        (ymax - lat) /
        ystep; // the raster texture uses a different coordinate system with y pointing downwards

    const std::size_t top = static_cast<std::size_t>(fmax(floor(ythP), 0));
    const std::size_t bottom = static_cast<std::size_t>(fmin(ceil(ythP), height - 1));
    const std::size_t left = static_cast<std::size_t>(fmax(floor(xthP), 0));
    const std::size_t right = static_cast<std::size_t>(fmin(ceil(xthP), width - 1));

    // Calculate distances from corners for bilinear interpolation
    const float fromLeft = xthP - left; // this is the fraction part of xthP
    const float fromTop = ythP - top;   // this is the fraction part of ythP
    const float fromRight = 1 - fromLeft;
    const float fromBottom = 1 - fromTop;

    return {static_cast<std::int32_t>(raster_data(left, top) * (fromRight * fromBottom) +
                                      raster_data(right, top) * (fromLeft * fromBottom) +
                                      raster_data(left, bottom) * (fromRight * fromTop) +
                                      raster_data(right, bottom) * (fromLeft * fromTop))};
}

// Load raster source into memory
int GeoTiffContainer::LoadGeoTiffSource(const std::string &path_string)
{
    const auto itr = LoadedSourcePaths.find(path_string);
    if (itr != LoadedSourcePaths.end())
    {
        util::Log() << "[source loader] Already loaded source '" << path_string << "' at source_id "
                    << itr->second;
        return itr->second;
    }

    int source_id = static_cast<int>(LoadedSources.size());

    util::Log() << "[source loader] Loading from " << path_string << "  ... ";
    TIMER_START(loading_source);

    boost::filesystem::path filepath(path_string);
    if (!boost::filesystem::exists(filepath))
    {
        throw util::RuntimeError(
            path_string, ErrorCode::FileOpenError, SOURCE_REF, "File not found");
    }

    GDALDataset  *poDataset;
    double  adfGeoTransform[6];

    GDALAllRegister();
    poDataset = (GDALDataset *) GDALOpen(path_string.c_str(), GA_ReadOnly);
    if (!poDataset) {
        throw util::RuntimeError(
            path_string, ErrorCode::FileOpenError, SOURCE_REF, "File cannot opened");
    }

    // Read Data from GeoTIFF
    std::size_t ncols = poDataset->GetRasterXSize();
    std::size_t nrows = poDataset->GetRasterYSize();
    GeoTiffGrid rasterData{poDataset, ncols, nrows};

    if (poDataset->GetGeoTransform(adfGeoTransform) == CE_None) {
        auto _xmin = static_cast<int>(util::toFixed(util::FloatLongitude{adfGeoTransform[0]}));
        auto _ymax = static_cast<int>(util::toFixed(util::FloatLongitude{adfGeoTransform[3]}));
        auto _xmax = _xmin + static_cast<int>(util::toFixed(util::FloatLongitude{adfGeoTransform[1] * ncols})); 
        auto _ymin = _ymax + static_cast<int>(util::toFixed(util::FloatLongitude{adfGeoTransform[5] * nrows}));

        util::Log() << "[GeoTIFF loader] Size : " << ncols << " x " << nrows << " [" << _xmin << ", " << _ymin << "]-[" << _xmax << ", " << _ymax << "]";
        GeoTiffSource source{std::move(rasterData), ncols, nrows, _xmin, _xmax, _ymin, _ymax};
        TIMER_STOP(loading_source);
        LoadedSourcePaths.emplace(path_string, source_id);
        LoadedSources.push_back(std::move(source));

        util::Log() << "[source loader] ok, after " << TIMER_SEC(loading_source) << "s";
    }
    GDALClose(poDataset);

    return source_id;
}

// External function for looking up nearest data point from a specified source
GeoTiffDatum GeoTiffContainer::GetGeoTiffDataFromSource(unsigned int source_id, double lon, double lat)
{
    if (LoadedSources.size() < source_id + 1)
    {
        throw util::exception("Attempted to access source " + std::to_string(source_id) +
                              ", but there are only " + std::to_string(LoadedSources.size()) +
                              " loaded" + SOURCE_REF);
    }

    BOOST_ASSERT(lat < 90);
    BOOST_ASSERT(lat > -90);
    BOOST_ASSERT(lon < 180);
    BOOST_ASSERT(lon > -180);

    const auto &found = LoadedSources[source_id];
    return found.GetGeoTiffData(static_cast<std::int32_t>(util::toFixed(util::FloatLongitude{lon})),
                               static_cast<std::int32_t>(util::toFixed(util::FloatLatitude{lat})));
}

// External function for looking up interpolated data from a specified source
GeoTiffDatum
GeoTiffContainer::GetGeoTiffInterpolateFromSource(unsigned int source_id, double lon, double lat)
{
    if (LoadedSources.size() < source_id + 1)
    {
        throw util::exception("Attempted to access source " + std::to_string(source_id) +
                              ", but there are only " + std::to_string(LoadedSources.size()) +
                              " loaded" + SOURCE_REF);
    }

    BOOST_ASSERT(lat < 90);
    BOOST_ASSERT(lat > -90);
    BOOST_ASSERT(lon < 180);
    BOOST_ASSERT(lon > -180);

    const auto &found = LoadedSources[source_id];
    return found.GetGeoTiffInterpolate(
        static_cast<std::int32_t>(util::toFixed(util::FloatLongitude{lon})),
        static_cast<std::int32_t>(util::toFixed(util::FloatLatitude{lat})));
}
}
}
