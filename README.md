usage: pysurfonoi.py [-h] [-dft] [-db DB] [-u U] [-pw PW] [-load input.shp]
                     [-loadmask input.shp] [-status STATUS] [-vis {m,d,v}]
                     [-tri] [-vor] [-wt] [-cont] [-inter INTER] [-resetinter]
                     [-cout contour_output.shp] [-c C [C ...]]

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

optional arguments:
  -h, --help            show this help message and exit
  -dft                  get db credentials from default file
  -db DB                database name
  -u U                  database user
  -pw PW                database password
  -load input.shp       load measurements from shapefile, inputfile.shp
  -loadmask input.shp   load mask polygon from shapefile, inputfile.shp
  -status STATUS        retrieve status
  -vis {m,d,v}          visualize features
  -tri                  constructs the delaunay triangulation
  -vor                  constructs the voronoi diagram
  -wt                   establishes the worker table
  -cont                 contours
  -inter INTER          interpolate
  -resetinter           reset interpolation
  -cout contour_output.shp
                        contour shapefile
  -c C [C ...]          isobaths values
  
  # quickstart
  - createdb
  - create extension postgis;
  - -load input.shp
  - -status m / -vis m
  - -tri -vor -wt
  - -cout contours_original.shp -c 99901 (shortcut for a set of isobaths)
  - inter [iterations]
  - -cout contours_generalised.shp -c 99901
  
