#!/usr/local/bin/python3.7

# preliminaries:
# createdb -U [user] [dbname]
# psql -U [user] -d [dbname] -c "CREATE EXTENSION postgis;"
# add functions contour_lines + helper function

import psycopg2
import osgeo.ogr
import shapely
import shapely.wkt
from shapely.geometry import Point, LineString, Polygon
import numpy as np

import sys
import subprocess
import os

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt

from postgis_functions import *


def are_perpendicular(lineOne, lineTwo):
    # https://stackoverflow.com/a/7572668
    x1, y1 = lineOne.coords[0][0], lineOne.coords[0][1]
    x2, y2 = lineOne.coords[1][0], lineOne.coords[1][1]
    x3, y3 = lineTwo.coords[0][0], lineTwo.coords[0][1]
    x4, y4 = lineTwo.coords[1][0], lineTwo.coords[1][1]

    a = np.array([x2-x1, y2-y1])
    b = np.array([x4-x3, y4-y3])
    scalarProduct = np.dot(a, b)
    orthoValue = abs(scalarProduct/(lineOne.length*lineTwo.length))

    return orthoValue


def progressBar(value, endvalue, bar_length=30):
    # https://stackoverflow.com/a/37630397
    percent = float(value) / endvalue
    arrow = '#' * int(round(percent * bar_length))  # + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rProgress: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()


def get_default_database(dftFile='default_db'):
    with open(dftFile) as fi:
        lines = [line.rstrip('\n') for line in fi]
    return lines[0], lines[1], lines[2]


def db_connection(dbName, dbUser, dbPass):
    # establishes the connection needed for psycopg2
    connection = psycopg2.connect(database=dbName, user=dbUser, password=dbPass)
    cursor = connection.cursor()
    return connection, cursor


def load_measurements(dbName, dbUser, dbPass, inputMeasurements):
    # inserts measurements from a shapefile into the database in a table called measurements
    # https://github.com/agaidus/PostgreSQL_PostGIS_Databases_Python/blob/master/Build_Query_Spatial_Database.ipynb
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    # create table
    cursor.execute('DROP TABLE IF EXISTS measurements')
    cursor.execute(
        'CREATE TABLE measurements (id SERIAL PRIMARY KEY, pos GEOMETRY, depth NUMERIC)')
    cursor.execute('CREATE INDEX measurements_index ON measurements USING GIST(pos)')

    # load points from shapefile
    shapefile = osgeo.ogr.Open(inputMeasurements)
    layer = shapefile.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        depth = feature.GetField('depth')
        geometry = feature.GetGeometryRef()
        wkt = geometry.ExportToWkt()
        # Insert data into database, converting WKT geometry to a PostGIS geometry
        cursor.execute(
            "INSERT INTO measurements (pos, depth) VALUES (ST_GeomFromText('{}'), {})".format(wkt, depth))

    connection.commit()

    return


def load_mask(dbName, dbUser, dbPass, inputMask):
    # inserts mask (polygon) from a shapefile into the database in a table called mask
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    # create table
    cursor.execute('DROP TABLE IF EXISTS mask')
    cursor.execute(
        'CREATE TABLE mask (id SERIAL PRIMARY KEY, geom GEOMETRY)')
    cursor.execute('CREATE INDEX mask_index ON mask USING GIST(geom)')

    # load points from shapefile
    shapefile = osgeo.ogr.Open(inputMask)
    layer = shapefile.GetLayer(0)
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        geometry = feature.GetGeometryRef()
        wkt = geometry.ExportToWkt()
        # Insert data into database, converting WKT geometry to a PostGIS geometry
        cursor.execute(
            "INSERT INTO mask (geom) VALUES (ST_GeomFromText('{}'))".format(wkt))

    connection.commit()

    return


def construct_delaunay(dbName, dbUser, dbPass):
    # creates delaunay triangles and stores the polygons into table delaunay_cells
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    # create table
    cursor.execute('DROP TABLE IF EXISTS delaunay_cells')
    cursor.execute(
        'CREATE TABLE delaunay_cells (id SERIAL PRIMARY KEY, geom GEOMETRY)')
    cursor.execute('CREATE INDEX delaunay_cells_index ON delaunay_cells USING GIST(geom)')

    cursor.execute(
        'SELECT ST_AsText((ST_Dump(geom)).geom) As triangle FROM(SELECT ST_DelaunayTriangles(ST_Collect(pos)) As geom FROM measurements) AS foo')

    for tri in cursor.fetchall():
        cursor.execute(
            "INSERT INTO delaunay_cells (geom) VALUES (ST_GeomFromText('{}'))".format(tri[0]))

    connection.commit()

    # select d.geom
    # from delaunay_cells d,
    # 	mask m
    # WHERE ST_Within(d.geom, m.geom);

    return


def construct_voronoi(dbName, dbUser, dbPass):
    # creates voronoi diagram  and stores the polygons into table voronoi_cells
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    # create table
    cursor.execute('DROP TABLE IF EXISTS voronoi_cells')
    cursor.execute(
        'CREATE TABLE voronoi_cells (id SERIAL PRIMARY KEY, geom GEOMETRY)')
    cursor.execute('CREATE INDEX voronoi_cells_index ON voronoi_cells USING GIST(geom)')

    cursor.execute(
        'SELECT ST_AsText((ST_Dump(geom)).geom) As cell FROM(SELECT ST_VoronoiPolygons(ST_Collect(pos)) As geom FROM measurements) AS foo')

    for cell in cursor.fetchall():
        cursor.execute(
            "INSERT INTO voronoi_cells (geom) VALUES (ST_GeomFromText('{}'))".format(cell[0]))

    connection.commit()

    # SELECT ST_Intersection(ST_Buffer(a.geom,50), b.geom) AS geom
    # FROM mask a
    # CROSS JOIN voronoi_cells b;

    return


def get_measurements(dbName, dbUser, dbPass):
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    cursor.execute('SELECT id, ST_AsText(pos), depth FROM measurements')
    rowsList = []
    for fid, pos, depth in cursor:
        data = {'FID': fid, 'geometry': shapely.wkt.loads(pos), 'depth': depth}
        rowsList.append(data)
    gdf = gpd.GeoDataFrame(rowsList, crs='epsg:28992').set_index('FID')

    return gdf


def status_measurements(indexedValues):
    print('Status\nNumber of measurements: {}\nMin depth: {}\nMax depth: {}\nAvg depth: {}'.format(
        len(indexedValues), indexedValues.depth.min(), indexedValues.depth.max(), indexedValues.depth.mean()))
    print(indexedValues.head())
    return


def visualize_measurements(indexedValues):
    indexedValues.plot(marker='x', column='depth', cmap='ocean', s=0.5)
    plt.grid(b=True, linewidth='0.5')

    plt.show()
    return


def get_delaunay_voronoi(dbName, dbUser, dbPass, dOrV):
    connection, cursor = db_connection(dbName, dbUser, dbPass)
    tableName = dOrV

    cursor.execute('SELECT id, ST_AsText(geom) FROM {}'.format(tableName))
    rowsList = []
    for fid, geom in cursor:
        data = {'FID': fid, 'geometry': shapely.wkt.loads(geom)}
        rowsList.append(data)
    gdf = gpd.GeoDataFrame(rowsList, crs='epsg:28992').set_index('FID')

    return gdf


def visualize_delaunay_voronoi(dbName, dbUser, dbPass, dOrV):
    dataFrame = get_delaunay_voronoi(dbName, dbUser, dbPass, dOrV)

    print(dataFrame.head())
    dataFrame.plot()
    plt.grid(b=True, linewidth='0.5')
    plt.show()

    return


def establish_worker_table(dbName, dbUser, dbPass):
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    # some generic info
    nrMeasurements = len(get_measurements(dbName, dbUser, dbPass))

    # create table
    cursor.execute('DROP TABLE IF EXISTS interpolator')
    cursor.execute(
        'CREATE TABLE interpolator (mid INT, neighbor_id INT, distance_d NUMERIC, distance_v NUMERIC)')
    connection.commit()

    # find the voronoi cell and neighboring triangles for each measurement
    cursor.execute('SELECT m.id mid, ARRAY_AGG(v.id) vids, ARRAY_AGG(d.id) dids FROM measurements m JOIN voronoi_cells v ON ST_Contains(v.geom, m.pos) JOIN delaunay_cells d ON ST_Intersects(m.pos, d.geom) GROUP BY m.id ORDER BY m.id')

    print('> found voronois and delaunays..')
    print('> iterating over every measurement..')

    nr = 0
    for measurement in cursor.fetchall():

        # print progress
        progressBar(nr, nrMeasurements)

        # select the measurement
        cursor.execute(
            'SELECT ST_AsText(pos) FROM measurements WHERE id = {}'.format(measurement[0]))
        currentPoint = cursor.fetchall()[0][0]
        currentPointGeom = shapely.wkt.loads(currentPoint)

        # select voronoi cell of the current measurement
        cursor.execute(
            'SELECT ST_AsText(geom) FROM voronoi_cells WHERE id = {}'.format(measurement[1][0]))
        voronoiPolygon = shapely.wkt.loads(cursor.fetchall()[0][0])

        # split the voronoi cell into line segments
        voronoiCoordinates = voronoiPolygon.exterior.coords
        voronoiSegments = []
        for i in range(len(voronoiCoordinates)-1):
            segment = LineString([Point(voronoiCoordinates[i][0], voronoiCoordinates[i][1]), Point(
                voronoiCoordinates[i+1][0], voronoiCoordinates[i+1][1])])
            voronoiSegments.append(segment)

        # find neighboring points through neighboring triangles
        triangleNeighbors = []
        for triangleId in measurement[2]:
            triangleNeighbors.append(triangleId)
        cursor.execute(
            "SELECT DISTINCT ST_AsText(geom) FROM (SELECT (ST_DumpPoints(g.geom)).* FROM (SELECT geom::geometry AS geom FROM delaunay_cells WHERE id IN ({})) AS g) j".format(str(triangleNeighbors)[1:-1]))

        measurementNeighbors = []
        for neighbor in cursor.fetchall():
            if neighbor[0] != currentPoint:  # wktPoint
                neighborPointGeom = shapely.wkt.loads(neighbor[0])

                # d_i
                delaunayLine = LineString([currentPointGeom, neighborPointGeom])
                distanceDelaunay = delaunayLine.length

                # find id of the neighbor
                cursor.execute(
                    "SELECT id FROM measurements WHERE pos = ST_GeomFromText('{}')".format(neighbor[0]))
                neighborId = cursor.fetchall()[0][0]
                measurementNeighbors.append(neighborId)

                # find voronoi twin segment for this neighbor
                smallestOrthovalue = 1
                voronoiTwin = None
                for voronoiSegment in voronoiSegments:
                    currentOrthovalue = are_perpendicular(voronoiSegment, delaunayLine)
                    if currentOrthovalue < smallestOrthovalue:
                        smallestOrthovalue = currentOrthovalue
                        voronoiTwin = voronoiSegment
                # v_i
                distanceVoronoi = voronoiTwin.length

                # print('\n--Measurement--\mid:\t\t{}\nnid:\t\t{}'.format(measurement[0], neighborId))
                cursor.execute(
                    "INSERT INTO interpolator (mid, neighbor_id, distance_d, distance_v) VALUES ({}, {}, {}, {})".format(measurement[0], neighborId, distanceDelaunay, distanceVoronoi))
                connection.commit()
        nr += 1

    connection.commit()

    return


def establish_contours(dbName, dbUser, dbPass):
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    measurementLimit = 500
    contourLimit = 500
    # SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT 100), ARRAY(SELECT depth FROM measurements LIMIT 100), array[2.0,5.0,8.0])).* FROM measurements LIMIT 100;

    cursor.execute(function_get_cell_intersects())
    cursor.execute(function_contour_lines())

    # create table
    cursor.execute('DROP TABLE IF EXISTS contour_lines_output')
    cursor.execute(
        "SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT {}), ARRAY(SELECT depth FROM measurements LIMIT {}), array[2.0])).* INTO contour_lines_output FROM measurements LIMIT {}".format(measurementLimit, measurementLimit, contourLimit))

    connection.commit()
    return


def export_contours(dbName, dbUser, dbPass, outputFile, breaks):
    connection, cursor = db_connection(dbName, dbUser, dbPass)
    serverLocation = 'localhost'

    measurementLimit = 100
    contourLimit = 10

    cursor.execute(function_get_cell_intersects())
    cursor.execute(function_contour_lines())

    breakList = str(breaks)

    # contour lines as SQL
    # create table
    cursor.execute('DROP TABLE IF EXISTS contour_lines_export')
    connection.commit()
    sqlString = 'SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT {}),ARRAY(SELECT depth FROM measurements LIMIT {}),ARRAY{})).* INTO contour_lines_export FROM measurements LIMIT {}'.format(
        measurementLimit, measurementLimit, breakList, contourLimit)
    # unlimited contours
    sqlString = 'SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT {}),ARRAY(SELECT depth FROM measurements LIMIT {}),ARRAY{})).* INTO contour_lines_export FROM measurements'.format(
        measurementLimit, measurementLimit, breakList)
    # print(sqlString)
    cursor.execute(sqlString)
    connection.commit()

    print('> contours saved to table')
    print('> exporting to shapefile')

    cwd = os.getcwd()
    cwd = ''
    filePath = os.path.join(cwd, outputFile)

    shpString = 'SELECT * FROM contour_lines_export'
    subprocess.call(['pgsql2shp', '-f', filePath, '-h', serverLocation, '-u',
                     dbUser, '-P', dbPass, dbName, shpString])

    connection.commit()


def manual_contours(dbName, dbUser, dbPass, outputFile, breaks):
    # creates delaunay triangles and stores the polygons into table delaunay_cells
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    cursor.execute(
        'SELECT ST_AsText((ST_Dump(geom)).geom) As triangle FROM(SELECT ST_DelaunayTriangles(ST_Collect(ST_MakePoint(ST_X(pos), ST_Y(pos), depth))) As geom FROM measurements) AS foo')

    contour_lines = []
    for tri in cursor.fetchall():
        # print('\n---NEW TRIANGLE---')
        # print(tri[0])
        triangle = shapely.wkt.loads(tri[0])
        triangleCoords = triangle.exterior.coords
        zValues = [triangleCoords[0][2], triangleCoords[1][2], triangleCoords[2][2]]
        zMin = min(zValues)
        zMax = max(zValues)
        # print(zValues)

        triangleSegments = []
        for i in range(len(triangleCoords)-1):
            segment = LineString([Point(triangleCoords[i]), Point(triangleCoords[i+1])])
            # segment = LineString([Point(voronoiCoordinates[i][0], voronoiCoordinates[i][1]), Point(
            #     voronoiCoordinates[i+1][0], voronoiCoordinates[i+1][1])])
            triangleSegments.append(segment)

        # could check if all breaks are above the max or all breaks below the min, later
        # TODO
        for value in breaks:  # iterate over every traingle just once, but check for every break
            # print('-- {}'.format(value))
            pointList = []
            if value >= zMin and value <= zMax:
                # break is inside triangle
                valueCount = zValues.count(value)
                if valueCount != 3:
                    # traingle not horizontal
                    if valueCount == 2:
                        # one edge is simply the contour segment
                        # print('edge is segment')
                        for segment in triangleSegments:
                            if segment.coords[0][2] == segment.coords[1][2]:
                                # this is the exact isobath
                                # print(segment)
                                contour_lines.append([segment, value])

                    elif valueCount == 1:
                        # may be tocuhing, or intersectecs to the opposite
                        if value != zMin and value != zMax:
                            # opposite intersect + vertex
                            # print('vertex + opposite')

                            for segment in triangleSegments:
                                xOne, yOne, zOne = segment.coords[0][0], segment.coords[0][1], segment.coords[0][2]
                                xTwo, yTwo, zTwo = segment.coords[1][0], segment.coords[1][1], segment.coords[1][2]
                                sMin, sMax = min(zOne, zTwo), max(zOne, zTwo)
                                if zOne == value:
                                    pointList.append(Point(xOne, yOne))
                                    #print(xOne, yOne)
                                if sMin < value < sMax:  # !!! not being the precise vertex
                                    x = (value*xTwo-value*xOne-zOne*xTwo+zTwo*xOne)/(zTwo-zOne)
                                    y = (x*(yTwo-yOne)-xOne*yTwo+xTwo*yOne)/(xTwo-xOne)
                                    #print(x, y)
                                    pointList.append(Point(x, y))

                    else:
                        # needs interpolation
                        # print('interpolated segment')
                        #pX, pY, qX, qY
                        for segment in triangleSegments:
                            xOne, yOne, zOne = segment.coords[0][0], segment.coords[0][1], segment.coords[0][2]
                            xTwo, yTwo, zTwo = segment.coords[1][0], segment.coords[1][1], segment.coords[1][2]
                            sMin, sMax = min(zOne, zTwo), max(zOne, zTwo)
                            if sMin <= value <= sMax:
                                x = (value*xTwo-value*xOne-zOne*xTwo+zTwo*xOne)/(zTwo-zOne)
                                y = (x*(yTwo-yOne)-xOne*yTwo+xTwo*yOne)/(xTwo-xOne)
                                #print(x, y)
                                pointList.append(Point(x, y))
                        # print(pointList)

            # add the actual isobath-segments
            if not len(pointList) == 0:
                # print(LineString(pointList).coords[0])
                # print(value)
                contour_lines.append([LineString(pointList), value])
                #save_shp(LineString(pointList), value)
    save_shp(contour_lines)

    return


def save_shp(contourList):

    # Here's an example Shapely geometry
    poly = Polygon([(0, 0), (0, 1), (1, 1), (0, 0)])

    # Now convert it to a shapefile with OGR
    driver = osgeo.ogr.GetDriverByName('Esri Shapefile')
    ds = driver.CreateDataSource('manual_contours.shp')
    layer = ds.CreateLayer('', None, osgeo.ogr.wkbLineString)
    # Add one attribute
    layer.CreateField(osgeo.ogr.FieldDefn('depth', osgeo.ogr.OFTReal))
    defn = layer.GetLayerDefn()

    # If there are multiple geometries, put the "for" loop here

    for entry in contourList:

        # Create a new feature (attribute and geometry)
        feat = osgeo.ogr.Feature(defn)
        feat.SetField('depth', entry[1])

        # Make a geometry, from Shapely object
        geom = osgeo.ogr.CreateGeometryFromWkb(entry[0].wkb)
        feat.SetGeometry(geom)

        layer.CreateFeature(feat)
        feat = geom = None  # destroy these

    # Save and close everything
    ds = layer = feat = geom = None

    return
