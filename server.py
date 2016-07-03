#! /usr/bin/python

import sys, math
from flask import Flask, request, jsonify, send_file
from osgeo import gdal, osr, ogr


# Ellipsoid model constants (actual values here are for WGS84)
sm_a = 6378137.0
sm_b = 6356752.314
sm_EccSquared = 6.69437999013e-03
UTMScaleFactor = 0.9996

app = Flask(__name__)

#http://127.0.0.1:5000/date/20160501/crop_by_bounding_box?min_x=44.2&max_x=44.15&min_y=27.4&max_y=27.45
@app.route('/date/<file_name>/crop_by_bounding_box')
def crop_by_bounding_box(file_name):
	min_x = float(request.args.get('min_x'))
	max_x = float(request.args.get('max_x'))
	min_y = float(request.args.get('min_y'))
	max_y = float(request.args.get('max_y'))
	#
	# Convert input coordinates from lat/lon to UTM
	#
	XYmin = lat_lon_to_utm(min_x, min_y)
	XYmax = lat_lon_to_utm(max_x, max_y)
	#
	# Register Imagine driver and open file
	#
	driver = gdal.GetDriverByName('GTiff')
	driver.Register()
	dataset = gdal.Open('/home/sant/test/satellite_images/{}.tif'.format(file_name))
	if dataset is None:
		print 'Could not open {}.tif'.format(file_name)
		sys.exit(1)
	#
	# Getting image dimensions
	#
	cols = dataset.RasterXSize
	rows = dataset.RasterYSize
	bands = dataset.RasterCount
	#
	# Getting georeference info
	#
	transform = dataset.GetGeoTransform()
	xOrigin = transform[0]
	yOrigin = transform[3]
	pixelWidth = transform[1]
	pixelHeight = -transform[5]
	#
	# Computing Point1(i1,j1), Point2(i2,j2) and new raster size
	#
	i1 = int((XYmin[0] - xOrigin) / pixelWidth)
	j1 = int((yOrigin - XYmin[1]) / pixelHeight)
	i2 = int((XYmax[0] - xOrigin) / pixelWidth)
	j2 = int((yOrigin - XYmax[1]) / pixelHeight)
	new_cols = i2 - i1 + 1
	new_rows = j2 - j1 + 1
	#
	# Create  list to store band data in
	#
	band_list = []
	#
	# Read in bands and store all the data in bandList
	#
	for i in range(bands):
		band = dataset.GetRasterBand(i+1) # 1-based index
		data = band.ReadAsArray(i1, j1, new_cols, new_rows)
		band_list.append(data)

	new_x = xOrigin + i1 * pixelWidth
	new_y = yOrigin - j1 * pixelHeight
	new_transform = (new_x, transform[1], transform[2], new_y, transform[4], transform[5])
	#
	# Create gtif file
	#
	output_file = '/home/sant/test/{}_cutted_by_bounding_box.tif'.format(file_name)
	driver = gdal.GetDriverByName("GTiff")
	dst_ds = driver.Create(output_file, new_cols, new_rows, 3, gdal.GDT_Byte)
	#
	# Writting output raster
	#
	for j in range(bands):
		data = band_list[j]
		dst_ds.GetRasterBand(j+1).WriteArray(data)
	#
	# Setting extension of output raster
	#
	dst_ds.SetGeoTransform(new_transform)
	wkt = dataset.GetProjection()
	#
	# Setting spatial reference of output raster
	#
	srs = osr.SpatialReference()
	srs.ImportFromWkt(wkt)
	dst_ds.SetProjection( srs.ExportToWkt() )
	#
	# Close output raster dataset
	#
	dataset = None
	dst_ds = None
	return send_file(output_file,  as_attachment=True)

#curl -H "Content-Type=application/json" --data (input json) -X POST http://127.0.0.1:5000/date/20160518/crop_by_geojson
@app.route('/date/<file_name>/crop_by_geojson', methods=['GET', 'POST']) 
def crop_by_geojson(file_name):
	geojson = jsonify(request.get_json(force=True))
	# You can use this geojson directly, just uncomment
	"""
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
	"""
	min_x, max_x, min_y, max_y = bbox(geojson)

	XYmin = lat_lon_to_utm(max_y, min_x)
	XYmax = lat_lon_to_utm(min_y, max_x)
	xValues = [XYmin[0], XYmax[0]]
	yValues = [XYmin[1], XYmax[1]]
	# Convert GeoJSON to list of coordinates
	x, y = zip(*list(explode(geojson['geometry']['coordinates'])))
	#
	# Register Imagine driver and open file
	#
	driver = gdal.GetDriverByName('GTiff')
	driver.Register()
	dataset = gdal.Open('/home/sant/test/satellite_images/{}.tif'.format(file_name))
	if dataset is None:
		print 'Could not open ' + file_name + '.tif'
		sys.exit(1)
	#
	# Getting image dimensions
	#
	cols = dataset.RasterXSize
	rows = dataset.RasterYSize
	bands = dataset.RasterCount
	#
	# Getting georeference info
	#
	transform = dataset.GetGeoTransform()
	projection = dataset.GetProjection()
	xOrigin = transform[0]
	yOrigin = transform[3]
	pixelWidth = transform[1]
	pixelHeight = -transform[5]
	#
	# Computing Point1(i1,j1), Point2(i2,j2)
	#
	i1 = int((xValues[0] - xOrigin) / pixelWidth)
	j1 = int((yOrigin - yValues[0]) / pixelHeight)
	i2 = int((xValues[1] - xOrigin) / pixelWidth)
	j2 = int((yOrigin - yValues[1]) / pixelHeight)
	new_cols = i2 - i1 + 1
	new_rows = j2 - j1 + 1
	# New upper-left X,Y values
	new_x = xOrigin + i1 * pixelWidth
	new_y = yOrigin - j1 * pixelHeight
	new_transform = (new_x, transform[1], transform[2], new_y, transform[4], transform[5])
	# Create inner ring
	innerRing = ogr.Geometry(ogr.wkbLinearRing)
	coords = []
	for i in range(len(x)):
		coords.append(lat_lon_to_utm(y[i],x[i]))
		innerRing.AddPoint(coords[i][0], coords[i][1])
	# Create outer box
	box = ogr.Geometry(ogr.wkbLinearRing)
	box.AddPoint(new_x, new_y) #Upper Left
	box.AddPoint((new_x + new_cols * pixelWidth), new_y) #Upper right
	box.AddPoint((new_x + new_cols * pixelWidth), (new_y - new_rows * pixelHeight)) #Lower Right
	box.AddPoint(new_x, (new_y - new_rows * pixelHeight)) #Lower left
	# Create polygon
	wkt_geom  = ogr.Geometry(ogr.wkbPolygon)
	wkt_geom.AddGeometry(box)
	wkt_geom.AddGeometry(innerRing)
	# Create a memory layer to rasterize from.
	driver = ogr.GetDriverByName('Memory')
	memds = driver.CreateDataSource('tmpmemds')

	lyr = memds.CreateLayer('poly', geom_type=ogr.wkbPolygon)
	feat = ogr.Feature(lyr.GetLayerDefn())
	feat.SetGeometryDirectly( ogr.Geometry(wkt = wkt_geom.ExportToWkt()) )
	lyr.CreateFeature(feat)
	# Create output file
	output_file = '/home/sant/test/{}_cutted_by_geojson.tif'.format(file_name)
	driver = gdal.GetDriverByName('GTiff')
	outds = driver.Create(output_file, new_cols, new_rows, bands, gdal.GDT_Byte)
	#
	# Read in bands and store all the data in bandList
	#
	band_list = []
	for i in range(bands):
		band = dataset.GetRasterBand(i+1) # 1-based index
		data = band.ReadAsArray(i1, j1, new_cols, new_rows)
		band_list.append(data)
	for j in range(bands):
		data = band_list[j]
		outds.GetRasterBand(j+1).WriteArray(data)

	outds.SetProjection(projection)
	outds.SetGeoTransform(new_transform)
	gdal.RasterizeLayer(outds, [1,2,3], lyr, burn_values = [0,0,0])
	dataset = None
	outds = None
	return send_file(output_file,  as_attachment=True)

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

def lat_lon_to_utm(lat, lon):
	""" Converts lat lon to utm
	Inputs:
	lat - lattitude in degrees
	lon - longitude in degrees
	Outputs:
	xy - utm x(easting), y(northing)
	zone - utm zone
	hemi - 'N' or 'S' """
	if ((lon < -180.0) or (180.0 <= lon)):
		print 'The longitude you entered is out of range -', lon
		print 'Please enter a number in the range [-180, 180).'
		return 0

	if ((lat < -90.0) or (90.0 < lat)):
		print 'The latitude you entered is out of range -', lat
		print 'Please enter a number in the range [-90, 90].'

	# Compute the UTM zone.
	zone = math.floor ((lon + 180.0) / 6) + 1

	# Convert
	xy = lat_lon_to_utm_xy (deg_to_rad(lat), deg_to_rad(lon), zone)

	# Determine hemisphere
	hemi = 'N'
	if (lat < 0):
		hemi = 'S'

	return xy

def lat_lon_to_utm_xy(lat, lon, zone):
	'''
	Converts a latitude/longitude pair to x and y coordinates in the
	Universal Transverse Mercator projection.
	Inputs:
	lat - Latitude of the point, in radians.
	lon - Longitude of the point, in radians.
	zone - UTM zone to be used for calculating values for x and y.
	If zone is less than 1 or greater than 60, the routine
	will determine the appropriate zone from the value of lon.
	Outputs:
	xy - A 2-element array where the UTM x and y values will be stored.
	'''

	xy = map_lat_lon_to_xy(lat, lon, utm_central_meridian(zone))

	# Adjust easting and northing for UTM system.
	xy[0] = xy[0] * UTMScaleFactor + 500000.0
	xy[1] = xy[1] * UTMScaleFactor
	if (xy[1] < 0.0):
		xy[1] = xy[1] + 10000000.0

	return xy

def map_lat_lon_to_xy(phi, lambda_pt, lambda_ctr):
	'''
	Converts a latitude/longitude pair to x and y coordinates in the
	Transverse Mercator projection.  Note that Transverse Mercator is not
	the same as UTM a scale factor is required to convert between them.
	Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
	GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
	Inputs:
	phi - Latitude of the point, in radians.
	lambda_pt - Longitude of the point, in radians.
	lambda_ctr - Longitude of the central meridian to be used, in radians.
	Outputs:
	xy - A 2-element array containing the x and y coordinates
	of the computed point.
	'''

	# Precalculate ep2
	ep2 = (math.pow (sm_a, 2.0) - math.pow (sm_b, 2.0)) / math.pow (sm_b, 2.0)

	# Precalculate nu2
	nu2 = ep2 * math.pow (math.cos (phi), 2.0)

	# Precalculate N
	N = math.pow (sm_a, 2.0) / (sm_b * math.sqrt (1 + nu2))

	# Precalculate t
	t = math.tan (phi)
	t2 = t * t
	# tmp = (t2 * t2 * t2) - math.pow (t, 6.0)

	# Precalculate l
	l = lambda_pt - lambda_ctr

	# Precalculate coefficients for l**n in the equations below
	#   so a normal human being can read the expressions for easting
	#   and northing
	#   -- l**1 and l**2 have coefficients of 1.0
	l3coef = 1.0 - t2 + nu2

	l4coef = 5.0 - t2 + 9 * nu2 + 4.0 * (nu2 * nu2)

	l5coef = 5.0 - 18.0 * t2 + (t2 * t2) + 14.0 * nu2 \
		- 58.0 * t2 * nu2

	l6coef = 61.0 - 58.0 * t2 + (t2 * t2) + 270.0 * nu2 \
		- 330.0 * t2 * nu2

	l7coef = 61.0 - 479.0 * t2 + 179.0 * (t2 * t2) - (t2 * t2 * t2)

	l8coef = 1385.0 - 3111.0 * t2 + 543.0 * (t2 * t2) - (t2 * t2 * t2)

	# Calculate easting (x)
	xy = [0.0, 0.0]
	xy[0] = N * math.cos (phi) * l \
		+ (N / 6.0 * math.pow (math.cos (phi), 3.0) * l3coef * math.pow (l, 3.0)) \
		+ (N / 120.0 * math.pow (math.cos (phi), 5.0) * l5coef * math.pow (l, 5.0)) \
		+ (N / 5040.0 * math.pow (math.cos (phi), 7.0) * l7coef * math.pow (l, 7.0))

	# Calculate northing (y)
	xy[1] = arc_length_of_meridian (phi) \
		+ (t / 2.0 * N * math.pow (math.cos (phi), 2.0) * math.pow (l, 2.0)) \
		+ (t / 24.0 * N * math.pow (math.cos (phi), 4.0) * l4coef * math.pow (l, 4.0)) \
		+ (t / 720.0 * N * math.pow (math.cos (phi), 6.0) * l6coef * math.pow (l, 6.0)) \
		+ (t / 40320.0 * N * math.pow (math.cos (phi), 8.0) * l8coef * math.pow (l, 8.0))
	return xy

def deg_to_rad(deg):
	'''
	Converts degrees to radians.
	'''
	return (deg / 180.0 * math.pi)

def utm_central_meridian(zone):
	'''
	Determines the central meridian for the given UTM zone.
	Inputs:
	zone - An integer value designating the UTM zone, range [1,60].
	Outputs:
	The central meridian for the given UTM zone, in radians, or zero
	if the UTM zone parameter is outside the range [1,60].
	Range of the central meridian is the radian equivalent of [-177,+177].
	'''
	return deg_to_rad(-183.0 + (zone * 6.0))

def arc_length_of_meridian(phi):
	'''
	Computes the ellipsoidal distance from the equator to a point at a
	given latitude.
	Reference: Hoffmann-Wellenhof, B., Lichtenegger, H., and Collins, J.,
	GPS: Theory and Practice, 3rd ed.  New York: Springer-Verlag Wien, 1994.
	Inputs:
	phi - Latitude of the point, in radians.
	Globals:
	sm_a - Ellipsoid model major axis.
	sm_b - Ellipsoid model minor axis.
	Outputs:
	The ellipsoidal distance of the point from the equator, in meters.
	'''

	# Precalculate n
	n = (sm_a - sm_b) / (sm_a + sm_b)

	# Precalculate alpha
	alpha = ((sm_a + sm_b) / 2.0) \
		* (1.0 + (math.pow (n, 2.0) / 4.0) + (math.pow (n, 4.0) / 64.0))

	# Precalculate beta
	beta = (-3.0 * n / 2.0) + (9.0 * math.pow (n, 3.0) / 16.0) \
		+ (-3.0 * math.pow (n, 5.0) / 32.0)

	# Precalculate gamma
	gamma = (15.0 * math.pow (n, 2.0) / 16.0) \
		+ (-15.0 * math.pow (n, 4.0) / 32.0)

	# Precalculate delta
	delta = (-35.0 * math.pow (n, 3.0) / 48.0) \
		+ (105.0 * math.pow (n, 5.0) / 256.0)

	# Precalculate epsilon
	epsilon = (315.0 * math.pow (n, 4.0) / 512.0)

	# Now calculate the sum of the series and return
	result = alpha \
		* (phi + (beta * math.sin (2.0 * phi)) \
		   + (gamma * math.sin (4.0 * phi)) \
		   + (delta * math.sin (6.0 * phi)) \
		   + (epsilon * math.sin (8.0 * phi)))

	return result


if __name__ == '__main__':
	app.run(debug=True)