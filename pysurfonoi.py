#!/usr/local/bin/python3.7

# preliminaries:
# createdb -U [user] [dbname]
# psql -U [user] -d [dbname] -c "CREATE EXTENSION postgis;"

import argparse
import textwrap
import matplotlib.pyplot as plt


from pysurfonoi_workers import *


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
                                 ______      _____             __                  _
                                 | ___ \    /  ___|           / _|                (_)
                                 | |_/ /   _\ `--. _   _ _ __| |_ ___  _ __   ___  _
                                 |  __/ | | |`--. \ | | | '__|  _/ _ \| '_ \ / _ \| |
                                 | |  | |_| /\__/ / |_| | |  | || (_) | | | | (_) | |
                                 \_|   \__, \____/ \__,_|_|  |_| \___/|_| |_|\___/|_|
                                        __/ |
                                       |___/

                                 for interacting with hydrographic depth measurements
                                 to create safe navigational isobaths.
                                 ----------------------------------------------------
                                 '''))
parser.add_argument('-dft', action='store_true',
                    help='get db credentials from default file')
parser.add_argument('-db', default=None, type=str, help='database name')
parser.add_argument('-u', default=None, type=str, help='database user')
parser.add_argument('-pw', default='', type=str, help='database password')
parser.add_argument('-load', metavar='input.shp', default=None, type=str,
                    help='load measurements from shapefile, inputfile.shp')
parser.add_argument('-status', default=None, type=str, help='retrieve status')
#parser.add_argument('-vis', default=None, type=str, help='visualize features')
parser.add_argument("-vis", choices=["m", "d", "v"], type=str, help="visualize features")
parser.add_argument('-tri', action='store_true', help='constructs the delaunay triangulation')
parser.add_argument('-vor', action='store_true', help='constructs the voronoi diagram')
parser.add_argument('-wt', action='store_true', help='establishes the worker table')
parser.add_argument('-cont', action='store_true', help='contours')


args = parser.parse_args()
dbName = args.db
dbUser = args.u
dbPass = args.pw
if args.dft:
    dbName, dbUser, dbPass = get_default_database()
    print('accessing database: {}\nuser: {}\n'.format(dbName, dbUser))

if args.load:
    print('> loading measurements')
    load_measurements(dbName, dbUser, dbPass, args.load)
if args.tri:
    print('> constructing Delaunay triangulation')
    construct_delaunay(dbName, dbUser, dbPass)
if args.vor:
    print('> constructing Voronoi diagram')
    construct_voronoi(dbName, dbUser, dbPass)
if args.wt:
    print('> establishing worker table')
    establish_worker_table(dbName, dbUser, dbPass)
if args.cont:
    print('> establishing contours')
    contouring(dbName, dbUser, dbPass)

if args.status == 'm':
    print('> retrieving status of measurements')
    status_measurements(get_measurements(dbName, dbUser, dbPass))
if args.vis == 'm':
    print('> visualizing measurements')
    visualize_measurements(get_measurements(dbName, dbUser, dbPass))
if args.vis in ['d', 'v']:
    if args.vis == 'd':
        print('> visualizing delaunay triangles')
        visualize_delaunay_voronoi(dbName, dbUser, dbPass, 'delaunay_cells')
    elif args.vis == 'v':
        print('> visualizing voronoi cells')
        visualize_delaunay_voronoi(dbName, dbUser, dbPass, 'voronoi_cells')