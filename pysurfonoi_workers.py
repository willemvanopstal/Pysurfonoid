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
    cursor.execute('CREATE INDEX measurements_index ON measuremets USING GIST(pos)')

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

    # create table
    cursor.execute('DROP TABLE IF EXISTS interpolator')
    cursor.execute(
        'CREATE TABLE interpolator (mid INT, neighbor_id INT, distance_d NUMERIC, distance_v NUMERIC)')
    connection.commit()

    # find the voronoi cell and neighboring triangles for each measurement
    cursor.execute('SELECT m.id mid, ARRAY_AGG(v.id) vids, ARRAY_AGG(d.id) dids FROM measurements m JOIN voronoi_cells v ON ST_Contains(v.geom, m.pos) JOIN delaunay_cells d ON ST_Intersects(m.pos, d.geom) GROUP BY m.id ORDER BY m.id')
    print('> found voronois and delaunays..')

    for measurement in cursor.fetchall():
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

                #print('\n--Measurement--\mid:\t\t{}\nnid:\t\t{}'.format(measurement[0], neighborId))
                cursor.execute(
                    "INSERT INTO interpolator (mid, neighbor_id, distance_d, distance_v) VALUES ({}, {}, {}, {})".format(measurement[0], neighborId, distanceDelaunay, distanceVoronoi))
                connection.commit()

    connection.commit()

    return


def establish_contours(dbName, dbUser, dbPass):
    connection, cursor = db_connection(dbName, dbUser, dbPass)

    measurementLimit = 10
    contourLimit = 10
    # SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT 100), ARRAY(SELECT depth FROM measurements LIMIT 100), array[2.0,5.0,8.0])).* FROM measurements LIMIT 100;

    cursor.execute(function_get_cell_intersects())
    cursor.execute(function_contour_lines())

    # create table
    cursor.execute('DROP TABLE IF EXISTS contour_lines_output')
    cursor.execute(
        "SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT {}), ARRAY(SELECT depth FROM measurements LIMIT {}), array[2.0,4,5,8,10])).* INTO contour_lines_output FROM measurements LIMIT {}".format(measurementLimit, measurementLimit, contourLimit))

    connection.commit()
