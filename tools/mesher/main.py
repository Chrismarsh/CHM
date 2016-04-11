from osgeo import gdal,ogr,osr
import subprocess
import re
import json
import os
import numpy as np
import scipy.stats.mstats as sp
import sys


# Print iterations progress
# http://stackoverflow.com/a/34325723/410074
def printProgress(iteration, total, prefix='', suffix='', decimals=2, barLength=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iterations  - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
    """
    filledLength = int(round(barLength * iteration / float(total)))
    percents = round(100.00 * (iteration / float(total)), decimals)
    bar = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('%s [%s] %s%s %s\r' % (prefix, bar, percents, '%', suffix)),
    sys.stdout.flush()
    if iteration == total:
        print("\n")

#Zonal stats from here https://gist.github.com/perrygeo/5667173
def bbox_to_pixel_offsets(gt, bbox):
    originX = gt[0]
    originY = gt[3]
    pixel_width = gt[1]
    pixel_height = gt[5]
    x1 = int((bbox[0] - originX) / pixel_width)
    x2 = int((bbox[1] - originX) / pixel_width) + 1

    y1 = int((bbox[3] - originY) / pixel_height)
    y2 = int((bbox[2] - originY) / pixel_height) + 1

    xsize = x2 - x1
    ysize = y2 - y1
    return (x1, y1, xsize, ysize)



def extract_point(raster,mx,my):
    # Convert from map to pixel coordinates.
    # Only works for geotransforms with no rotation.

    rb = raster.GetRasterBand(1)
    gt = raster.GetGeoTransform()

    px = int((mx - gt[0]) / gt[1])  # x pixel
    py = int((my - gt[3]) / gt[5])  # y pixel

    #boundary verticies from Triangle can end up outside of the domain by 1 pixel.
    #if we adjusted back by  1 pixel to get the dz value, it's no problem and still gives a good boundary
    if px == raster.RasterXSize:
        px = px -1
    if py == raster.RasterYSize:
        py = py - 1
    mz = rb.ReadAsArray(px, py, 1, 1)
    mz = float(mz.flatten()[0])

    if mz == rb.GetNoDataValue():
        dx = 1
        dy = 1

        #attempt to pick points from the surrounding 8
        z1 = rb.ReadAsArray(px - dx if px > 0 else px, py, 1, 1)[0]
        z2 = rb.ReadAsArray(px + dx if px < raster.RasterXSize else px, py, 1, 1)[0]
        z3 = rb.ReadAsArray(px, py + dy if py < raster.RasterYSize else py, 1, 1)[0]
        z4 = rb.ReadAsArray(px, py - dy if py > 0 else py, 1, 1)[0]

        z5 = rb.ReadAsArray(px - dx if px > 0 else px,                  py - dy if py > 0 else py, 1, 1)[0]
        z6 = rb.ReadAsArray(px + dx if px < raster.RasterXSize else px, py - dy if py > 0 else py, 1, 1)[0]
        z7 = rb.ReadAsArray(px - dx if px > 0 else px,                  py + dy if py < raster.RasterYSize else py, 1, 1)[0]
        z8 = rb.ReadAsArray(px + dx if px < raster.RasterXSize else px, py - dy if py > 0 else py, 1, 1)[0]

        z=[z1, z2, z3, z4,z5,z6,z7,z8]
        z = [x for x in z if x != rb.GetNoDataValue()]

        if len(z) == 0:
            print 'Error: The point (%f,%f) and its 8-neighbours lies outside of the DEM domain' %(mx,my)
            exit(1)

        mz = np.mean(z)
    return mz

def main():

#######  user configurable paramters here    #######
    dem_filename = 'wolf_dem_30m.tif'

    parameter_files = {
        'landcover': { 'file' : 'eosd.tif',
                       'method':'mode'},  # mode, mean
        'svf':{'file':'wolf_svf1.tif',
               'method':'mean'}
    }

    #plgs_shp = None #plgs_shp constrains the triangulation. If this is None, the input raster's extent is used

    simplify     =   False
    simplify_tol =   5   #amount in meters to simplify the polygon by. Careful as too much will cause many lines to be outside of the bounds of the raster.

########################################################

    src_ds=gdal.Open(dem_filename)
    base_name = dem_filename[:dem_filename.rfind('.')]
    if src_ds is None:
        print 'Unable to open %s' % dem_filename
        exit(1)

    # paramater_fhandles={} #file handles for paramters
    # for key,file in parameter_files.iteritems():
    #     paramater_fhandles[key] = open(base_name + '_' + key + '.ele_data','w')

    #load up all the paramter files
    for key,data in parameter_files.iteritems():
        parameter_files[key]['file'] = gdal.Open( data['file'])
        if parameter_files[key]['file'] is None:
            print 'Error: Unable to open raster for: %s' % key
            exit(1)


    #if plgs_shp is None:

    plgs_shp = base_name + '.shp'

    dem = src_ds.GetRasterBand(1)


    # gt = src_ds.GetGeoTransform()

    #Step 1: Create a mask raster that has a uniform value for all cells
    # read all elevation data. Might be worth changing to the for-loop approach we use below so we don't have to read in all into ram.
#some raster
    Z = dem.ReadAsArray()


    #set all the non-nodata values to our mask value
    Z[ Z != dem.GetNoDataValue()] = 1

    #save to file
    (x, y) = Z.shape
    driver = gdal.GetDriverByName("GTiff")
    dst_datatype = gdal.GDT_Int16 #gdal.GDT_Float32   #<--- this is fine as we don't actually store elevation data.
    tmp_raster = 'mask_'+base_name+'.tif'
    dst_ds = driver.Create(tmp_raster, y, x, 1, dst_datatype)
    dst_ds.SetGeoTransform(src_ds.GetGeoTransform())
    dst_ds.SetProjection(src_ds.GetProjection())
    dst_ds.GetRasterBand(1).SetNoDataValue(dem.GetNoDataValue())
    dst_ds.GetRasterBand(1).WriteArray(Z)
    dst_ds.FlushCache()  #super key to get this to disk
    dst_ds = None #close file

    #raster -> polygon
    subprocess.check_call(['gdal_polygonize.py %s -b 1 -f "ESRI Shapefile" %s' % (tmp_raster,plgs_shp)], shell=True)

    print 'Converting polygon to linestring'
    exec_string = 'ogr2ogr -overwrite %s %s  -nlt LINESTRING' % ('line_'+plgs_shp, plgs_shp )
    if simplify:
        exec_string = exec_string + ' -simplify ' + simplify_tol

    subprocess.check_call(exec_string, shell=True)

    #convert to geoJSON because it's easy to parse
    poly_plgs = base_name + '.geojson'
    subprocess.check_call(['ogr2ogr -f GeoJSON %s %s' % (poly_plgs,'line_'+plgs_shp )], shell=True)

    with open(poly_plgs) as f:
        plgs = json.load(f)

    os.remove(poly_plgs)

    #assuming just the first feature is what we want. Need to add in more to support rivers and lakes
    if plgs['features'][0]['geometry']['type'] != 'LineString':
        print('Not linestring')
        exit(1)

    coords = plgs['features'][0]['geometry']['coordinates']

    #Create the PLGS to constrain the triangulation
    poly_file = 'PLGS'+base_name+'.poly'
    with open(poly_file,'w') as f:
        header = '%d 2 0 0\n' % (len(coords))
        f.write(header)
        vert = 1
        for c in coords:
            f.write( '%d %f %f\n' % (vert, c[0],c[1]))
            vert = vert + 1

        f.write('\n')
        header = '%d 0\n' % (len(coords))
        f.write(header)

        for i in range(len(coords)):
            if i+1 == len(coords): # last time has to loop back to start
                f.write('%d %d %d\n' % (i + 1, i + 1, 1))
            else:
                f.write('%d %d %d\n' % (i+1, i+1, i+2))

        f.write('0\n')

    subprocess.check_call(['./triangle -a1000 -p %s -n' % (poly_file)],shell=True)

    #read in the node, ele, and neigh from

    #holds our main mesh structure which we will write out to json to read into CHM
    mesh = {}
    mesh['mesh']={}
    mesh['mesh']['vertex']= []

    read_header=False

    print 'Reading nodes'
    with open('PLGS'+base_name+'.1.node') as f:
        for line in f:
            if '#' not in line:
                if not read_header:
                    header = re.findall(r"[+-]?\d+(?:\.\d+)?",line)
                    num_nodes = int(header[0])
                    read_header=True
                    mesh['mesh']['nvertex'] = num_nodes
                else:
                    items = re.findall(r"[+-]?\d+(?:\.\d+)?",line)
                    mx = float(items[1])
                    my = float(items[2])
                    mz = extract_point(src_ds,mx,my)

                    mesh['mesh']['vertex'].append( [mx,my,mz] )


    # read in the neighbour file
    print 'Reading in neighbour file'
    read_header = False
    mesh['mesh']['neigh'] = []
    with open('PLGS' + base_name + '.1.neigh') as elem:
        for line in elem:
            if '#' not in line:
                if not read_header:
                    read_header = True  # skip header, nothing to do here
                else:
                    items = re.findall(r"[+-]?\d+(?:\.\d+)?", line)
                    v0 = int(items[1]) - 1  # convert to zero indexing
                    v1 = int(items[2]) - 1
                    v2 = int(items[3]) - 1

                    mesh['mesh']['neigh'].append([v0, v1, v2])

    #Create the shape file to hold the triangulation. Could do all in memory but it is nice to see this in a GIS
    # set up the shapefile driver
    driver = ogr.GetDriverByName("ESRI Shapefile")
    # create the data source
    output_usm = driver.CreateDataSource(base_name + '_USM.shp')

    # create the spatial reference from the raster dataset
    wkt = src_ds.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)

    # create the layer
    layer = output_usm.CreateLayer(base_name, srs, ogr.wkbPolygon)
    layer.CreateField(ogr.FieldDefn("triangle", ogr.OFTInteger))  # holds the triangle id.

    for key,value in parameter_files.iteritems():
        layer.CreateField(ogr.FieldDefn(key,       ogr.OFTReal))

    read_header=False
    i=1

    print 'Parameterizing triangles'

    triangles_to_fix=[]
    mesh['mesh']['elem'] = []


    mesh['parameters']={}
    for key, data in parameter_files.iteritems():
        mesh['parameters'][key]=[]


    i=0
    with open('PLGS'+base_name+'.1.ele') as elem:
        for line in elem:
            if '#' not in line:
                if not read_header:
                    header = re.findall(r"[+-]?\d+(?:\.\d+)?",line)
                    nelem= int(header[0])
                    mesh['mesh']['nelem'] = nelem
                    read_header = True #skip header
                else:

                    printProgress(i,nelem)

                    items = re.findall(r"[+-]?\d+(?:\.\d+)?",line)
                    v0 = int(items[1])-1 #convert to zero indexing
                    v1 = int(items[2])-1
                    v2 = int(items[3])-1

                    mesh['mesh']['elem'].append( [v0, v1, v2] )

                    ring = ogr.Geometry(ogr.wkbLinearRing)
                    ring.AddPoint( mesh['mesh']['vertex'][v0][0], mesh['mesh']['vertex'][v0][1] )
                    ring.AddPoint( mesh['mesh']['vertex'][v1][0], mesh['mesh']['vertex'][v1][1] )
                    ring.AddPoint( mesh['mesh']['vertex'][v2][0], mesh['mesh']['vertex'][v2][1] )

                    tpoly = ogr.Geometry(ogr.wkbPolygon)
                    tpoly.AddGeometry(ring)

                    feature = ogr.Feature( layer.GetLayerDefn() )
                    feature.SetField('triangle',int(items[0])-1)
                    feature.SetGeometry(tpoly)

                    #get the value under each triangle from each paramter file
                    for key, data in parameter_files.iteritems():
                        raster = data['file']

                        wkt = raster.GetProjection()
                        srs = osr.SpatialReference()
                        srs.ImportFromWkt(wkt)

                        gt = raster.GetGeoTransform()
                        rb = raster.GetRasterBand(1)

                        src_offset = bbox_to_pixel_offsets(gt, feature.geometry().GetEnvelope())
                        src_array = rb.ReadAsArray(*src_offset)

                        # calculate new geotransform of the feature subset
                        new_gt = (
                            (gt[0] + (src_offset[0] * gt[1])),
                            gt[1],
                            0.0,
                            (gt[3] + (src_offset[1] * gt[5])),
                            0.0,
                            gt[5]
                        )
                        # Create a temporary vector layer in memory
                        mem_drv = ogr.GetDriverByName('Memory')
                        driver = gdal.GetDriverByName('MEM')
                        mem_ds = mem_drv.CreateDataSource('out')
                        mem_layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)
                        mem_layer.CreateFeature(feature.Clone())

                        # Rasterize it
                        rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
                        rvds.SetGeoTransform(new_gt)
                        gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1])
                        rv_array = rvds.ReadAsArray()


                        # Mask the source data array with our current feature
                        # we take the logical_not to flip 0<->1 to get the correct mask effect
                        # we also mask out nodata values explictly
                        masked = np.ma.MaskedArray(
                            src_array,
                            mask=np.logical_or(
                                src_array == rb.GetNoDataValue(),
                                np.logical_not(rv_array)
                            )
                        )

                        # feature_stats = {
                        #     'min': float(masked.min()),
                        #     'mean': float(masked.mean()),
                        #     'max': float(masked.max()),
                        #     'std': float(masked.std()),
                        #     'sum': float(masked.sum()),
                        #     'count': int(masked.count()),
                        #     'fid': int(feat.GetFID())
                        # }

                        output = -9999
                        if data['method'] == 'mode':
                            output = sp.mode(masked.flatten())[0][0]
                        elif data['method'] == 'mean':
                            if masked.count() > 0:
                                output = float(masked.mean())  #if it's entirely masked, then we get nan and a warning printed to stdout. would like to avoid showing this warning.
                            else:
                                output = 0
                        else:
                            print 'Error: unknown data aggregation method %s' % data['method']

                        #we have a triangle that only touches a no-value. So give up and use surrounding values.
                        if output == 0:
                            masked = np.ma.masked_where(
                                    src_array == rb.GetNoDataValue(), src_array )
                            output = sp.mode(masked.flatten())[0][0]


                        if output == 0:
                            #alright, we're in a bit of trouble now. We will save this triangle for later, and try to just duplciate any neighbouring triangle's data.
                            #this only seems to happen right at the edge of the triangulation.
                            triangles_to_fix.append( { 'key':key,'elem': int(items[0])-1 ,'nodata':rb.GetNoDataValue(),'method':data['method'] } )

                        feature.SetField(key, output)
                        mesh['parameters'][key].append( output )

                    layer.CreateFeature(feature)

                    i=i+1



    print 'There are %d triangles with no data' % len(triangles_to_fix)
    i=1
    for t in triangles_to_fix:
        f = layer.GetFeature( t['elem'])

        values = []
        #get each neighbour index to the problem triangle
        n0 = mesh['mesh']['elem'][t['elem']][0]
        n1 = mesh['mesh']['elem'][t['elem']][1]
        n2 = mesh['mesh']['elem'][t['elem']][2]

        nv = [ layer.GetFeature( n0 ).GetField( t['key'] ),
                        layer.GetFeature( n1 ).GetField( t['key'] ),
                        layer.GetFeature( n2 ).GetField( t['key'] ) ]

        masked = np.ma.masked_where(
            nv == t['nodata'], nv)

        output = -9999
        if t['method'] == 'mode':
            output = sp.mode(masked.flatten())[0][0]
        elif t['method'] == 'mean':
            if masked.count() > 0:
                output = float(masked.mean())  #if it's entirely masked, then we get nan and a warning printed to stdout. would like to avoid showing this warning.
            else:
                output = 0

        else:
            print 'Error: unknown data aggregation method %s' % data['method']

        print '[ %d / %d ] Changed nodata value to %f for triangle id: %d' % (i,len(triangles_to_fix),output,t['elem'])
        mesh['parameters'][t['key']][t['elem']] = output
        f.SetField( t['key'], output)
        layer.SetFeature(f)
        i=i+1

    print 'Saving mesh and parameters to file ' + base_name+'.mesh'
    with open(base_name+'.mesh', 'w') as outfile:
        json.dump(mesh, outfile,indent=4)


if __name__ == "__main__":

    main()