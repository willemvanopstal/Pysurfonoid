def function_get_cell_intersects():

    functionString = '''CREATE OR REPLACE FUNCTION _get_cell_intersects(
    IN vertex geometry[], -- vertices geometries
    IN vv numeric[], -- vertices values
    IN bu integer[], -- vertices buckets
    IN breaks numeric[], -- breaks
    IN i1 integer, -- first vertex index
    IN i2 integer -- last vertex index
    )
    RETURNS geometry[]  AS $$
    DECLARE
    result geometry[];
    BEGIN
    -- init the array of results
    result := array_fill(null::geometry, breaks::int[]);

    -- if the vertices buckets are different,
    -- there are might be intersections
    IF bu[i1] <> bu[i2] THEN
        SELECT
            array_agg(b.point) INTO result
        FROM
        (
            SELECT
                -- linear interpolation of vertices values
                -- to find the break point
                (t.x-vv[i1])/(vv[i2]-vv[i1]) AS p
            FROM unnest(breaks) AS t(x)
        ) AS a
        LEFT JOIN LATERAL
        (
        SELECT
            ST_LineInterpolatePoint(
                ST_MakeLine(vertex[i1], vertex[i2]),
                a.p
                )
            AS point
        WHERE
            -- we only want intersections between i1 and i2
            a.p BETWEEN 0 AND 1
        ) AS b
        ON 1=1;
    END IF;

    RETURN result;

    END;
    $$ LANGUAGE PLPGSQL IMMUTABLE'''

    return functionString


def function_contour_lines():
    functionString = '''CREATE OR REPLACE FUNCTION contour_lines(
    	IN geomin geometry[],
    	IN colin numeric[],
    	IN breaks numeric[]
    )
    RETURNS TABLE(geom geometry, break numeric)   AS $$
    DECLARE
    	bucketin integer[];
    	gs geometry[];
    	g geometry;
    	vertex geometry[];
    	vv numeric[];
    	bu integer[];
    	inter numeric[]            ;
    	interp12 geometry[];
    	interp23 geometry[];
    	interp31 geometry[];
    	segment geometry[];
    	running_merge geometry[];
    	i integer;
    BEGIN
    	WITH
    	a AS(
    		SELECT
    			width_bucket(t.x, breaks) AS bin
    		FROM unnest(colin) AS t(x)
    	)
    	SELECT array_agg(bin) INTO bucketin FROM a;

    	WITH
    	a AS (SELECT unnest(geomin) AS e),
    	b AS (SELECT ST_DelaunayTriangles(ST_Collect(a.e)) AS t FROM a),
    	c AS (SELECT (ST_Dump(t)).geom AS v FROM b)
    	SELECT array_agg(v) INTO gs FROM c;

    	i:= 0;

    	FOREACH g IN ARRAY gs
    	LOOP

    		SELECT
    			array_agg(a.v),
    			array_agg(b.c),
    			array_agg(b.bk)
    		INTO vertex, vv, bu
    		FROM
    		(
    			SELECT (ST_DumpPoints(g)).geom AS v limit 3
    		) as a
    		CROSS JOIN
    		LATERAL(
    			SELECT
    				t.*
    			FROM
    				unnest(geomin, colin, bucketin) AS t(geo, c, bk)
    			WHERE ST_Equals(geo, a.v)
    			LIMIT 1
    		) AS b;

    		CONTINUE WHEN bu[1] = bu[2] and bu[1] = bu[3];

    		interp12 := _get_cell_intersects(vertex, vv, bu, breaks,1,2);
    		interp23 := _get_cell_intersects(vertex, vv, bu, breaks,2,3);
    		interp31 := _get_cell_intersects(vertex, vv, bu, breaks,3,1);

    		WITH
    		a AS(
                SELECT
                    t.*
                FROM
                unnest(breaks, interp12, interp23, interp31) AS t(br, p12 , p23, p31)
    		),
    		b AS(
    		SELECT
                CASE
                WHEN
                (p12 IS NOT NULL AND p23 IS NOT NULL AND ST_equals(p12, p23)=false) OR
                (p23 IS NOT NULL AND p31 IS NOT NULL AND ST_equals(p23, p31)=false) OR
                (p31 IS NOT NULL AND p12 IS NOT NULL AND ST_equals(p31, p12)=false)
                THEN ST_MakeLine(ARRAY[p12, p23, p31]::geometry[])
                ELSE null::geometry END AS segm,
                br
    		FROM a
    		)
    		SELECT
    		    array_agg(b.segm) into segment
    		FROM unnest(breaks) AS c(x)
            LEFT JOIN b ON b.br = c.x;

    		IF i = 0 THEN
                running_merge = segment;
                i := 1;
    		ELSE
                WITH
                a AS(
                    SELECT
                        ST_CollectionExtract(x, 2) AS x,
                        y
                    FROM unnest(running_merge,segment) AS t(x,y)
                ),
                b AS(
                    SELECT
                    ST_collect(x,y) AS element
                    FROM a
                )
                SELECT
                    array_agg(element) INTO running_merge
                FROM b;
    		END IF;

    	END LOOP;

    	RETURN QUERY
    	WITH a AS(
                SELECT
                    br,
                    ST_CollectionExtract(geo, 2) AS geo
                FROM unnest(running_merge, breaks) AS t(geo, br)
            ),
            b AS(
                SELECT
                    ST_LineMerge(geo) AS geo,
                    br
                FROM a
            )
        SELECT
            geo AS geom,
            br AS break
    	FROM b;

    END;
    $$ LANGUAGE PLPGSQL IMMUTABLE'''

    return functionString
