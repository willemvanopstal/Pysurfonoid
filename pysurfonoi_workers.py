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
from mpl_toolkits.mplot3d import Axes3D
from math import sqrt, pow

from postgis_functions import *


def are_perpendicular(lineOne, lineTwo):
    # https://stackoverflow.com/a/7572668
    x1, y1 = lineOne.coords[0][0], lineOne.coords[0][1]
    x2, y2 = lineOne.coords[1][0], lineOne.coords[1][1]
    x3, y3 = lineTwo.coords[0][0], lineTwo.coords[0][1]
    x4, y4 = lineTwo.coords[1][0], lineTwo.coords[1][1]

    # a = np.array([[x1, y1], [x2, y2]])
    # b = np.array([[x3, y3], [x4, y4]])
    a = np.array([x2-x1, y2-y1])
    b = np.array([x4-x3, y4-y3])
    scalarProduct = np.dot(a, b)
    orthoValue = abs(scalarProduct/(lineOne.length*lineTwo.length))

    return orthoValue


def xy_from_point(wktPoint):
    # POINT(x.dx y.dy)
    wktX, wktY = wktPoint[6:-1].split(" ")
    return (float(wktX), float(wktY))


def distance_between_points(xOne, yOne, xTwo, yTwo):
    dist = sqrt(pow(xTwo - xOne, 2) + pow(yTwo - yOne, 2))
    return dist


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
    # inserts measurements from a shapefile into the database in a table calles measurements
    # https://github.com/agaidus/PostgreSQL_PostGIS_Databases_Python/blob/master/Build_Query_Spatial_Database.ipynb
    connection, cursor = db_connection(dbName, dbUser, dbPass)
    tableName = 'measurements'

    # create table
    cursor.execute('DROP TABLE IF EXISTS {}'.format(tableName))
    cursor.execute(
        'CREATE TABLE {} (id SERIAL PRIMARY KEY, pos GEOMETRY, depth NUMERIC)'.format(tableName))
    cursor.execute('CREATE INDEX {}_index ON {} USING GIST(pos)'.format(tableName, tableName))

    # load points from shapefile
    shapefile = osgeo.ogr.Open(inputMeasurements)
    layer = shapefile.GetLayer(0)
    cursor.execute('DELETE FROM {}'.format(tableName))
    for i in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(i)
        depth = feature.GetField('depth')
        geometry = feature.GetGeometryRef()
        wkt = geometry.ExportToWkt()
        # Insert data into database, converting WKT geometry to a PostGIS geography
        cursor.execute("INSERT INTO {} (pos, depth) VALUES (ST_GeomFromText('{}'), {})".format(
            tableName, wkt, depth))

    connection.commit()

    return


def construct_delaunay(dbName, dbUser, dbPass):
    # inserts measurements from a shapefile into the database in a table calles measurements
    # https://github.com/agaidus/PostgreSQL_PostGIS_Databases_Python/blob/master/Build_Query_Spatial_Database.ipynb
    connection, cursor = db_connection(dbName, dbUser, dbPass)
    tableName = 'delaunay_cells'

    # create table
    cursor.execute('DROP TABLE IF EXISTS {}'.format(tableName))
    cursor.execute(
        'CREATE TABLE {} (id SERIAL PRIMARY KEY, geom GEOMETRY)'.format(tableName))
    cursor.execute('CREATE INDEX {}_index ON {} USING GIST(geom)'.format(tableName, tableName))
    # connection.commit()

    cursor.execute(
        'SELECT ST_AsText((ST_Dump(geom)).geom) As triangle FROM(SELECT ST_DelaunayTriangles(ST_Collect(pos)) As geom FROM measurements) AS foo')

    for tri in cursor.fetchall():
        # print(tri[0])
        cursor.execute(
            "INSERT INTO {} (geom) VALUES (ST_GeomFromText('{}'))".format(tableName, tri[0]))

    # SELECT(ST_Dump(geom)).geom As wkt
    # FROM(SELECT ST_DelaunayTriangles(ST_Collect(pos)) As geom
    #      FROM measurements) As foo
    connection.commit()

    return


def construct_voronoi(dbName, dbUser, dbPass):
    # inserts measurements from a shapefile into the database in a table calles measurements
    # https://github.com/agaidus/PostgreSQL_PostGIS_Databases_Python/blob/master/Build_Query_Spatial_Database.ipynb
    connection, cursor = db_connection(dbName, dbUser, dbPass)
    tableName = 'voronoi_cells'

    # create table
    cursor.execute('DROP TABLE IF EXISTS {}'.format(tableName))
    cursor.execute(
        'CREATE TABLE {} (id SERIAL PRIMARY KEY, geom GEOMETRY)'.format(tableName))
    cursor.execute('CREATE INDEX {}_index ON {} USING GIST(geom)'.format(tableName, tableName))
    # connection.commit()

    cursor.execute(
        'SELECT ST_AsText((ST_Dump(geom)).geom) As triangle FROM(SELECT ST_VoronoiPolygons(ST_Collect(pos)) As geom FROM measurements) AS foo')

    for cell in cursor.fetchall():
        # print(tri[0])
        cursor.execute(
            "INSERT INTO {} (geom) VALUES (ST_GeomFromText('{}'))".format(tableName, cell[0]))

    # SELECT(ST_Dump(geom)).geom As wkt
    # FROM(SELECT ST_DelaunayTriangles(ST_Collect(pos)) As geom
    #      FROM measurements) As foo
    connection.commit()

    return


def get_measurements(dbName, dbUser, dbPass):
    connection, cursor = db_connection(dbName, dbUser, dbPass)
    tableName = 'measurements'

    cursor.execute('SELECT id, ST_AsText(pos), depth FROM {}'.format(tableName))
    rowsList = []
    for fid, pos, depth in cursor:
        data = {'FID': fid, 'geometry': shapely.wkt.loads(pos), 'depth': depth}
        rowsList.append(data)
    gdf = gpd.GeoDataFrame(rowsList, crs='epsg:28992').set_index('FID')
    # print(gdf.head())

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
    # segmentation fault 11......
    dataFrame = get_delaunay_voronoi(dbName, dbUser, dbPass, dOrV)

    print(dataFrame.head())
    dataFrame.plot()
    plt.grid(b=True, linewidth='0.5')
    plt.show()

    return


def establish_worker_table(dbName, dbUser, dbPass):
    #     SELECT	m.id AS mId,
    # 		ARRAY_AGG(v.id) AS vId,
    # 		ARRAY_AGG(d.id) AS dIds
    # FROM 	measurements m
    # 	JOIN
    # 		voronoi_cells v
    # 	ON ST_Contains(v.geom, m.pos)
    # 	JOIN
    # 		delaunay_cells d
    # 	ON ST_Intersects(m.pos, d.geom)
    # GROUP BY m.id
    # ORDER BY m.id;

    # retrieve points from polygon
    #     SELECT path, ST_AsText(geom)
    # FROM (
    #   SELECT (ST_DumpPoints(g.geom)).*
    #   FROM
    #     (SELECT
    #        geom::geometry AS geom
    # 	 FROM delaunay_cells
    # 	 WHERE id = 29
    #     ) AS g
    #   ) j;

    connection, cursor = db_connection(dbName, dbUser, dbPass)

    # create table
    cursor.execute('DROP TABLE IF EXISTS interpolator')
    cursor.execute(
        'CREATE TABLE interpolator (mid INT, neighbor_id INT, distance_d NUMERIC, distance_v NUMERIC)')
    connection.commit()

    cursor.execute('SELECT m.id mid, ARRAY_AGG(v.id) vids, ARRAY_AGG(d.id) dids FROM measurements m JOIN voronoi_cells v ON ST_Contains(v.geom, m.pos) JOIN delaunay_cells d ON ST_Intersects(m.pos, d.geom) GROUP BY m.id ORDER BY m.id')
    print('> found voronois and delaunays..')
    for measurement in cursor.fetchall():
        print('>')
        # print('\n\n------MEASUREMENT------')
        # print(measurement)
        cursor.execute(
            'SELECT ST_AsText(pos) FROM measurements WHERE id = {}'.format(measurement[0]))
        # print(cursor.fetchall()[0][0])
        currentPoint = cursor.fetchall()[0][0]
        # print(currentPoint)
        currentPointGeom = shapely.wkt.loads(currentPoint)
        #currentPointX, currentPointY = xy_from_point(currentPoint)

        # voronoi cell
        # print('--VD--')
        # print(measurement[1][0])
        cursor.execute(
            'SELECT ST_AsText(geom) FROM voronoi_cells WHERE id = {}'.format(measurement[1][0]))
        voronoiPolygon = shapely.wkt.loads(cursor.fetchall()[0][0])
        # print(voronoiPolygon)
        voronoiCoordinates = voronoiPolygon.exterior.coords
        voronoiSegments = []
        for i in range(len(voronoiCoordinates)-1):
            # LineString([Point(0, 0), Point(1, 1)]).wkt
            segment = LineString([Point(voronoiCoordinates[i][0], voronoiCoordinates[i][1]), Point(
                voronoiCoordinates[i+1][0], voronoiCoordinates[i+1][1])])
            # print(segment.length)
            voronoiSegments.append(segment)

        # triangles
        # print('--DT neighbors--')
        triangleNeighbors = []
        for triangleId in measurement[2]:
            # print(triangleId)
            triangleNeighbors.append(triangleId)

            # cursor.execute("SELECT ST_AsText(geom) FROM (SELECT (ST_DumpPoints(g.geom)).* FROM (SELECT geom::geometry AS geom FROM delaunay_cells WHERE id = {}) AS g) j GROUP BY st_astext".format(triangleId))
            # cursor.execute(
            #     'SELECT ST_AsText(geom) FROM delaunay_cells WHERE id = {}'.format(triangleId))
            # print(cursor.fetchall())
        # print(str(triangleNeighbors)[1:-1])
        cursor.execute(
            "SELECT DISTINCT ST_AsText(geom) FROM (SELECT (ST_DumpPoints(g.geom)).* FROM (SELECT geom::geometry AS geom FROM delaunay_cells WHERE id IN ({})) AS g) j".format(str(triangleNeighbors)[1:-1]))
        # print(cursor.fetchall())
        measurementNeighbors = []
        for neighbor in cursor.fetchall():
            # print(neighbor[0])
            if neighbor[0] != currentPoint:
                # print(neighbor[0])
                neighborPointGeom = shapely.wkt.loads(neighbor[0])
                # print(neighborPointGeom)
                delaunayLine = LineString([currentPointGeom, neighborPointGeom])
                #neighborX, neighborY = xy_from_point(neighbor[0])
                distanceDelaunay = delaunayLine.length
                # distanceDelaunay = distance_between_points(
                #     currentPointX, currentPointY, neighborX, neighborY)
                cursor.execute(
                    "SELECT id FROM measurements WHERE pos = ST_GeomFromText('{}')".format(neighbor[0]))
                neighborId = cursor.fetchall()[0][0]
                measurementNeighbors.append(neighborId)
                # print(neighborId)
                # print(distanceDelaunay)

                smallestOrthovalue = 1
                voronoiTwin = None
                for voronoiSegment in voronoiSegments:
                    currentOrthovalue = are_perpendicular(voronoiSegment, delaunayLine)
                    # print(currentOrthovalue)
                    if currentOrthovalue < smallestOrthovalue:
                        smallestOrthovalue = currentOrthovalue
                        voronoiTwin = voronoiSegment
                distanceVoronoi = voronoiTwin.length
                # print(distanceVoronoi)
                # print(smallestOrthovalue)

                #print('\n--Measurement--\nid:\t\t{}\nnid:\t\t{}'.format(measurement[0], neighborId))
                cursor.execute(
                    "INSERT INTO interpolator (mid, neighbor_id, distance_d, distance_v) VALUES ({}, {}, {}, {})".format(measurement[0], neighborId, distanceDelaunay, distanceVoronoi))
                connection.commit()
                # print(len(voronoiSegments), len(measurementNeighbors))  # outer voronois dont have neighbors

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
    # cursor.execute(
    #     'CREATE TABLE contour_lines_output (geom GEOMETRY, break NUMERIC)')

    cursor.execute(
        "SELECT (contour_lines(ARRAY(SELECT pos FROM measurements LIMIT {}), ARRAY(SELECT depth FROM measurements LIMIT {}), array[2.0,4,5,8,10])).* INTO contour_lines_output FROM measurements LIMIT {}".format(measurementLimit, measurementLimit, contourLimit))
    connection.commit()
