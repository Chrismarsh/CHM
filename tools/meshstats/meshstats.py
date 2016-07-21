#Given a mesh in a .shp file and a set of raster, this computes various error metrics of the mesh to the raster

from osgeo import gdal,ogr,osr
import os
import numpy as np
# import matplotlib.pyplot as plt
import math
import sys
import csv

def main():
    print 'Reading in files'
    mesher_output_dir = '../mesher/wolf1m_fill/'
    # raster_file ='/Users/chris/Documents/PhD/research/CHM/paper1/figures/meshes/granger/1mtol/wolf_lidar1_projected.tif'
    # shp_file='/Users/chris/Documents/PhD/research/CHM/paper1/figures/meshes/granger/1mtol/wolf_lidar1_1m_rmse.shp'

    #############
    if 'mesher_output_dir' in locals():
        base_name =os.path.basename(os.path.normpath(mesher_output_dir))
        raster_file = os.path.normpath(mesher_output_dir)+'/'+base_name+'_projected.tif'
        shp_file = os.path.normpath(mesher_output_dir)+'/'+base_name+'_USM.shp'
    elif 'raster_file' not in locals() or 'shp_file' not in locals():
        print 'If mesher output folder not given, must manually specify both shp and raster file'
        exit(1)
    else:
        base_name =os.path.basename(os.path.normpath(raster_file))
        base_shp_name =os.path.basename(os.path.normpath(shp_file))

    raster_ds = gdal.Open(raster_file)
    if raster_ds is None:
        print 'Could not open %s' % (raster_ds)
        exit(1)


    rb = raster_ds.GetRasterBand(1)
    src_array = rb.ReadAsArray(0,0,raster_ds.RasterXSize-1,raster_ds.RasterYSize-1)
    masked = np.ma.masked_where(src_array == rb.GetNoDataValue(),src_array )
    c = masked.count()

    src_array = None
    masked = None
    rb = None
    if c == 0:
        print "Only no data present in raster"
        exit(1)


    print "Number of raster cells = " + str(c)


    # driver = ogr.GetDriverByName('ESRI Shapefile')
    mesh = ogr.Open(shp_file, update=True)
    if mesh is None:
        print 'Could not open %s' % (mesh)
        exit(1)


    layer = mesh.GetLayer()
    num_elem = layer.GetFeatureCount()
    print "Number of triangles = %d" % (num_elem)

    print "# triangles = " + str( num_elem / float(c) * 100.) + '%'

    area = []
    angles = []
    rmse_value = []

    layerDefinition = layer.GetLayerDefn()
    field_names = []
    for i in range(layerDefinition.GetFieldCount()):
        field_names.append(layerDefinition.GetFieldDefn(i).GetName())

    if 'area' not in field_names:
        layer.CreateField(ogr.FieldDefn('area', ogr.OFTReal))
    if 'min_angle' not in field_names:
        layer.CreateField(ogr.FieldDefn('min_angle', ogr.OFTReal))
    if 'max_angle' not in field_names:
        layer.CreateField(ogr.FieldDefn('max_angle', ogr.OFTReal))
    if 'rmse' not in field_names:
        layer.CreateField(ogr.FieldDefn('rmse', ogr.OFTReal))

    i=0
    for feature in layer:
        printProgress(i,num_elem)

        geom = feature.GetGeometryRef()
        ring = geom.GetGeometryRef(0)

        # tri_id = feature.GetField('triangle')

        geom_area = geom.GetArea()
        area.append(geom_area)
        feature.SetField('area', geom_area)

        angle = tri_angles(raster_ds, ring)
        #if angles comes back none, we've hit an edge case where the triangle slightly sits outside of the domain.
        #this would have been fixed in the output mesh from mesher, but just ignore it here

        triangle=[] #this triangles angles
        if angle is not None:
            for a in angle:
                angles.append(a)
                triangle.append(a)

            feature.SetField('min_angle',min(triangle))
            feature.SetField('max_angle',max(triangle))

        rmse = tri_rmse(raster_ds, feature)
        if rmse is not None:
            rmse_value.append(rmse)
            feature.SetField('rmse',rmse)

        layer.SetFeature(feature)
        i+=1

    mesh= None

    #
    with open(base_shp_name+'_stats.csv','w') as f:
        writer = csv.writer(f)
        writer.writerow(["rmse", "area", "angle"])
        max_len = max(len(rmse_value), len(area), len(angles))
        for i in range(max_len):
            try:
                cur_rmse_value = rmse_value[i]
            except IndexError:
                cur_rmse_value = None
            try:
                cur_area = area[i]
            except IndexError:
                cur_area = None

            try:
                cur_angle = angles[i]
            except IndexError:
                cur_angle = None

            writer.writerow([cur_rmse_value, cur_area, cur_angle])

def tri_angles(raster_ds, ring):
    x1, y1, z = ring.GetPoint(0)  # this z is 0,
    z1 = extract_point(raster_ds, x1, y1);
    x2, y2, z = ring.GetPoint(1)
    z2 = extract_point(raster_ds, x2, y2);
    x3, y3, z = ring.GetPoint(2)
    z3 = extract_point(raster_ds, x3, y3);

    rb = raster_ds.GetRasterBand(1)
    if z1 == rb.GetNoDataValue() or z2 == rb.GetNoDataValue() or z3 == rb.GetNoDataValue():
        return None

    # vector length
    def vl(v):
        return math.sqrt(v[0] ** 2 + v[1] ** 2)

    # dot product
    def dp(u, v):
        return u[0] * v[0] + u[1] * v[1]

    u = (y3 - y1, x3 - x1)
    v = (y2 - y1, x2 - x1)
    angle32 = math.acos(dp(u, v) / (vl(u) * vl(v))) * 180. / math.pi
    u = (y1 - y3, x1 - x3)
    v = (y2 - y3, x2 - x3)
    angle12 = math.acos(dp(u, v) / (vl(u) * vl(v))) * 180. / math.pi
    u = (y1 - y2, x1 - x2)
    v = (y3 - y2, x3 - x2)
    angle13 = math.acos(dp(u, v) / (vl(u) * vl(v))) * 180. / math.pi

    return(angle32,angle12,angle13)


def tri_rmse(raster_ds,feature):
    geom = feature.GetGeometryRef()
    ring = geom.GetGeometryRef(0)

    new_gt, rtri,xsize,ysize = rasterize_elem(raster_ds, feature, '')

    #get triangle verticies
    x1, y1, z = ring.GetPoint(0)  # this z is 0,
    z1 = extract_point(raster_ds, x1, y1)
    x1, y1 = xyToPixel(new_gt,x1,y1,xsize,ysize)

    x2, y2, z = ring.GetPoint(1)
    z2 = extract_point(raster_ds, x2, y2)
    x2, y2 = xyToPixel(new_gt,x2,y2,xsize,ysize)

    x3, y3, z = ring.GetPoint(2)
    z3 = extract_point(raster_ds, x3, y3)
    x3, y3 = xyToPixel(new_gt,x3,y3,xsize,ysize)

    rb = raster_ds.GetRasterBand(1)
    if z1 == rb.GetNoDataValue() or z2 == rb.GetNoDataValue() or z3 == rb.GetNoDataValue():
        return None

    u1 = x2 - x1
    u2 = y2 - y1
    u3 = z2 - z1

    v1 = x3 - x1
    v2 = y3 - y1
    v3 = z3 - z1

    a = u2 * v3 - v2 * u3
    b = v1 * u3 - u1 * v3
    c = u1 * v2 - v1 * u2

    d =a*x1 + b*y1 + c*z1

    rmse = 0
    n = 0
    if c == 0:    #collinear
        return None

    max_diff = -9999;
    for y in xrange(rtri.shape[0]):
        for x in xrange(rtri.shape[1]):
            if not np.ma.is_masked(rtri[y,x]):
                pred = -(a*x+b*y-d)/c #plane eqn solved for z. allows us to predict z values via x,y coords
                value = rtri[y,x]
                diff = value - pred
                rmse += diff ** 2
                n+=1
                if abs(diff) > max_diff:
                    max_diff = abs(diff)

    if n==0:
        return None
    rmse /= n
    rmse = math.sqrt(rmse)

    # return max_diff
    # rmse = ((z1+z2+z3)/3)-float(rtri.mean())
    return abs(rmse)

#Zonal stats from here https://gist.github.com/perrygeo/5667173
def bbox_to_pixel_offsets(gt, bbox, rasterXsize, rasterYsize):
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

    #only apply this correction if we are touching the underlying raster.
    # if x1 < rasterXsize and y1 < rasterYsize:
        #deal with small out of bounds
    if x1 < 0:
        x1 = 0

    if y1 < 0:
        y1 = 0

    if x1 + xsize >= rasterXsize:
        xsize = rasterXsize-x1

    if y1 + ysize >= rasterYsize:
        ysize = rasterYsize-y1

    return (x1, y1, xsize, ysize)

def xyToPixel( gt, x,  y,  max_x,  max_y):
    # adfGeoTransform[0] /* top left x */
    # adfGeoTransform[1] /* w-e pixel resolution */
    # adfGeoTransform[2] /* 0 */
    # adfGeoTransform[3] /* top left y */
    # adfGeoTransform[4] /* 0 */
    # adfGeoTransform[5] /* n-s pixel resolution (negative value) */

    px = int((x - gt[0]) / gt[1])  # x pixel
    py = int((y - gt[3]) / gt[5])  # y pixel

    #out of bound issue when we are off by one
    if px >= max_x:
        px = max_x - 1;

    if py >= max_y:
        py = max_y - 1;

    if py < 0:
        py = 0;

    if px <0 :
        px = 0;

    return (px,py)

def rasterize_elem(raster, feature, key):

    wkt = raster.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    gt = raster.GetGeoTransform()
    rb = raster.GetRasterBand(1)
    src_offset = bbox_to_pixel_offsets(gt, feature.geometry().GetEnvelope(), raster.RasterXSize, raster.RasterYSize)

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
    mem_ds = mem_drv.CreateDataSource('out')
    # mem_drv = ogr.GetDriverByName('ESRI Shapefile')
    # mem_ds = mem_drv.CreateDataSource('rastertri.shp')

    mem_layer = mem_ds.CreateLayer('poly', srs, ogr.wkbPolygon)
    mem_layer.CreateFeature(feature.Clone())

    # Rasterize it
    driver = gdal.GetDriverByName('MEM')
    rvds = driver.Create('', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
    # driver = gdal.GetDriverByName('GTiff')
    # rvds = driver.Create('rasterized.tif', src_offset[2], src_offset[3], 1, gdal.GDT_Byte)
    rvds.SetGeoTransform(new_gt)

    gdal.RasterizeLayer(rvds, [1], mem_layer, burn_values=[1], options=['ALL_TOUCHED=TRUE'])
    rv_array = rvds.ReadAsArray()  # holds a mask of where the triangle is on the raster

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
    # output = rb.GetNoDataValue()
    #
    # if data['method'] == 'mode':
    #     output = sp.mode(masked.flatten())[0][0]
    # elif data['method'] == 'mean':
    #     if masked.count() > 0:
    #         output = float(
    #             masked.mean())  # if it's entirely masked, then we get nan and a warning printed to stdout. would like to avoid showing this warning.
    # else:
    #     print 'Error: unknown data aggregation method %s' % data['method']

    return new_gt,masked,src_offset[2],src_offset[3]
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
        z1 = rb.ReadAsArray(px - dx if px > 0 else px, py, 1, 1)
        z1 = z1[0] if z1 is not None else rb.GetNoDataValue()

        z2 = rb.ReadAsArray(px + dx if px < raster.RasterXSize else px, py, 1, 1)
        z2 = z2[0] if z2 is not None else rb.GetNoDataValue()

        z3 = rb.ReadAsArray(px, py + dy if py < raster.RasterYSize else py, 1, 1)
        z3 = z3[0] if z3 is not None else rb.GetNoDataValue()

        z4 = rb.ReadAsArray(px, py - dy if py > 0 else py, 1, 1)[0]
        z4 = z4[0] if z4 is not None else rb.GetNoDataValue()

        z5 = rb.ReadAsArray(px - dx if px > 0 else px,                  py - dy if py > 0 else py, 1, 1)
        z5 = z5[0] if z5 is not None else rb.GetNoDataValue()

        z6 = rb.ReadAsArray(px + dx if px < raster.RasterXSize else px, py - dy if py > 0 else py, 1, 1)
        z6 = z6[0] if z6 is not None else rb.GetNoDataValue()

        z7 = rb.ReadAsArray(px - dx if px > 0 else px,                  py + dy if py < raster.RasterYSize else py, 1, 1)
        z7 = z7[0] if z7 is not None else rb.GetNoDataValue()

        z8 = rb.ReadAsArray(px + dx if px < raster.RasterXSize else px, py - dy if py > 0 else py, 1, 1)
        z8 = z8[0] if z8 is not None else rb.GetNoDataValue()

        z=[z1, z2, z3, z4,z5,z6,z7,z8]
        z = [x for x in z if x != rb.GetNoDataValue()]

        if len(z) == 0:
            print 'Warning: The point (%f,%f) and its 8-neighbours lies outside of the DEM domain' %(mx,my)
            return rb.GetNoDataValue()
            #exit(1)

        mz = float(np.mean(z))
    return mz


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

if __name__ == "__main__":

    main()
