#! /usr/bin/python

import sys
import os
import numpy
from flask import Flask, request, send_file
try:
    from osgeo import ogr, osr, gdal
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')


base_dir = os.path.dirname(os.path.abspath(__file__))
image_dir = base_dir + '/satellite_images/'

app = Flask(__name__)

# http://127.0.0.1:5000/date/20160501/crop_by_bounding_box?min_x=44.2&max_x=44.15&min_y=27.4&max_y=27.45


@app.route('/date/<file_name>/crop_by_bounding_box')
def crop_by_bounding_box(file_name):
    min_x = float(request.args.get('min_x'))
    max_x = float(request.args.get('max_x'))
    min_y = float(request.args.get('min_y'))
    max_y = float(request.args.get('max_y'))

    # Register Imagine driver and open file
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    dataset = gdal.Open('{}/{}.tif'.format(image_dir, file_name))
    if dataset is None:
        print 'Could not open {}.tif'.format(file_name)
        sys.exit(1)

    # Getting spatial reference and projection of input raster
    projection = dataset.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)

    # WGS84 projection reference
    OSR_WGS84_REF = osr.SpatialReference()
    OSR_WGS84_REF.ImportFromEPSG(4326)

    # OSR transformation
    wgs84_to_image_trasformation = osr.CoordinateTransformation(OSR_WGS84_REF, srs)

    # Convert input coordinates from lat/lon to UTM
    XYmin = wgs84_to_image_trasformation.TransformPoint(min_y, min_x)
    XYmax = wgs84_to_image_trasformation.TransformPoint(max_y, max_x)

    # Getting image dimensions
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    bands = dataset.RasterCount

    # Getting georeference info
    transform = dataset.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    # Computing Point1(i1,j1), Point2(i2,j2) and new raster size
    i1 = int((XYmin[0] - xOrigin) / pixelWidth)
    j1 = int((yOrigin - XYmin[1]) / pixelHeight)
    i2 = int((XYmax[0] - xOrigin) / pixelWidth)
    j2 = int((yOrigin - XYmax[1]) / pixelHeight)

    new_cols = i2 - i1 + 1
    new_rows = j2 - j1 + 1

    # Create  list to store band data in
    band_list = []

    # Read in bands and store all the data in bandList
    for i in range(bands):
        band = dataset.GetRasterBand(i + 1)  # 1-based index
        data = band.ReadAsArray(i1, j1, new_cols, new_rows)
        band_list.append(data)

    new_x = xOrigin + i1 * pixelWidth
    new_y = yOrigin - j1 * pixelHeight
    new_transform = (new_x, transform[1], transform[2], new_y, transform[4], transform[5])

    # Create gtif file
    output_file = '{}/{}_cutted_by_bounding_box.tif'.format(image_dir, file_name)
    driver = gdal.GetDriverByName("GTiff")
    dst_ds = driver.Create(output_file, new_cols, new_rows, bands, gdal.GDT_Byte)

    # Writting output raster
    for j in range(bands):
        data = band_list[j]
        dst_ds.GetRasterBand(j + 1).WriteArray(data)

    # Setting extension of output raster
    dst_ds.SetGeoTransform(new_transform)
    wkt = dataset.GetProjection()

    # Setting spatial reference of output raster
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dst_ds.SetProjection(srs.ExportToWkt())

    # Close output raster dataset
    dataset = None
    dst_ds = None
    return send_file(output_file, as_attachment=True)

# curl -H "Content-Type=application/json" --data """ GeoJSON """ -X POST http://127.0.0.1:5000/date/20160518/crop_by_geojson


@app.route('/date/<file_name>/crop_by_geojson', methods=['GET', 'POST'])
def crop_by_geojson(file_name):
    # geojson = jsonify(request.get_json(force=True))
    geojson = {
        "type": "Feature",
        "properties": {
        "id": 1
        },
        "geometry": {
        "type": "Multipolygon",
        "coordinates": [[[
          [27.397511483867362, 44.161117466010516],
          [27.393924221672666, 44.159751598403503],
          [27.393556666460618, 44.159252063395591],
          [27.393726740035870, 44.158373985750522],
          [27.392040835956994, 44.157378400690988],
          [27.390354358253163, 44.156239034941315],
          [27.390977924658255, 44.152849194060536],
          [27.391438333095618, 44.149298658002031],
          [27.386781918912796, 44.147461728155896],
          [27.384487250437232, 44.146859408403664],
          [27.382636468741264, 44.156671855578281],
          [27.383891699721374, 44.156645049015140],
          [27.384649769913505, 44.157388133683327],
          [27.385547083122507, 44.160232076255667],
          [27.387997850095061, 44.160722084482430],
          [27.390672446485077, 44.161638147279866],
          [27.395361188085396, 44.163429614137918],
          [27.396513835695238, 44.162325787855522],
          [27.397511483867362, 44.161117466010516]
    ]]]
    }
    }

    min_x, max_x, min_y, max_y = bbox(geojson)

    # Register Imagine driver and open file
    driver = gdal.GetDriverByName('GTiff')
    driver.Register()
    dataset = gdal.Open('{}/{}.tif'.format(image_dir, file_name))
    if dataset is None:
        print 'Could not open ' + file_name + '.tif'
        sys.exit(1)

    # Getting image dimensions
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    bands = dataset.RasterCount

    # Getting georeference info
    transform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = -transform[5]

    # Getting spatial reference of input raster
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)

    # WGS84 projection reference
    OSR_WGS84_REF = osr.SpatialReference()
    OSR_WGS84_REF.ImportFromEPSG(4326)

    # OSR transformation
    wgs84_to_image_trasformation = osr.CoordinateTransformation(OSR_WGS84_REF, srs)
    XYmin = wgs84_to_image_trasformation.TransformPoint(min_x, max_y)
    XYmax = wgs84_to_image_trasformation.TransformPoint(max_x, min_y)

    # Computing Point1(i1,j1), Point2(i2,j2)
    i1 = int((XYmin[0] - xOrigin) / pixelWidth)
    j1 = int((yOrigin - XYmin[1]) / pixelHeight)
    i2 = int((XYmax[0] - xOrigin) / pixelWidth)
    j2 = int((yOrigin - XYmax[1]) / pixelHeight)
    new_cols = i2 - i1 + 1
    new_rows = j2 - j1 + 1

    # New upper-left X,Y values
    new_x = xOrigin + i1 * pixelWidth
    new_y = yOrigin - j1 * pixelHeight
    new_transform = (new_x, transform[1], transform[2], new_y, transform[4], transform[5])

    wkt_geom = ogr.CreateGeometryFromJson("""{}""".format(geojson['geometry']))
    wkt_geom.Transform(wgs84_to_image_trasformation)

    target_ds = gdal.GetDriverByName('MEM').Create('', new_cols, new_rows, bands, gdal.GDT_Byte)
    target_ds.SetGeoTransform(new_transform)
    target_ds.SetProjection(projection)

    # Create a memory layer to rasterize from.
    driver = ogr.GetDriverByName('Memory')
    memds = driver.CreateDataSource('tmpmemds')

    lyr = memds.CreateLayer('poly', geom_type=ogr.wkbMultiPolygon)
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetGeometryDirectly(ogr.Geometry(wkt=wkt_geom.ExportToWkt()))
    lyr.CreateFeature(feat)

    gdal.RasterizeLayer(target_ds, [1, 2, 3], lyr, burn_values=[1, 1, 1])

    # Create output file
    output_file = '{}/{}_cutted_by_geojson.tif'.format(image_dir, file_name)
    driver = gdal.GetDriverByName('GTiff')
    outds = driver.Create(output_file, new_cols, new_rows, bands, gdal.GDT_Byte)

    # Read in bands and store all the data in bandList
    mask_array = []
    band_list = []
    for i in range(bands):
        band_list.append(dataset.GetRasterBand(i + 1).ReadAsArray(i1, j1, new_cols, new_rows))
        mask_array.append(target_ds.GetRasterBand(i + 1).ReadAsArray())

    for j in range(bands):
        data = numpy.where(mask_array[j] == 1, band_list[j], 0)
        outds.GetRasterBand(j + 1).WriteArray(data)

    outds.SetProjection(projection)
    outds.SetGeoTransform(new_transform)

    dataset = None
    outds = None
    return send_file(output_file, as_attachment=True)


def explode(coords):
    """Explode a GeoJSON geometry's coordinates object and yield coordinate tuples.
    As long as the input is conforming, the type of the geometry doesn't matter."""
    for e in coords:
        if isinstance(e, (float, int, long)):
            yield coords
            break
        else:
            for f in explode(e):
                yield f


def bbox(f):
    x, y = zip(*list(explode(f['geometry']['coordinates'])))
    return min(x), max(x), min(y), max(y)


if __name__ == '__main__':
    app.run(debug=True)
