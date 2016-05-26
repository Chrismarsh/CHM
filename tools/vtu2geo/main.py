import vtk
import sys
from osgeo import gdal,ogr,osr
import glob

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


def main():
    gdal.UseExceptions()  # Enable errors
    base = "granger"
    # input_path = '/Users/chris/Documents/PhD/code/CHM/output_no_slope_30m/'
    input_path = '/home/chris/Documents/PhD/research/emergence/output/output_slope_1m/'
    # output_path = '/home/chris/Documents/PhD/code/CHM/output_tif/'
    EPSG=26908 #wolf 8N

    variables = ['total_inf','total_excess']  #set to None to dump all variables
    parameters = ['Aspect'] # paramters are one offs we want to extract from the vtu files
    pixel_size = 10 # (m)
    #####

    reader = vtk.vtkXMLUnstructuredGridReader()
    files =  glob.glob(input_path + base+'*.vtu')


    iter=1
    for path in files:
        printProgress(iter,len(files))
        reader.SetFileName(path)
        reader.Update()

        mesh = reader.GetOutput()

        driver = ogr.GetDriverByName('Memory')
        # driver = ogr.GetDriverByName("ESRI Shapefile")
        # output_usm = driver.CreateDataSource('lol.shp')
        output_usm = driver.CreateDataSource('out')

        srs = osr.SpatialReference()
        srs.ImportFromEPSG(EPSG)

        layer = output_usm.CreateLayer('poly', srs, ogr.wkbPolygon)

        cd = mesh.GetCellData()


        for i in range(0,cd.GetNumberOfArrays()):
            # print cd.GetArrayName(i)
            layer.CreateField(ogr.FieldDefn(cd.GetArrayName(i), ogr.OFTReal))

        #build the triangulation geometery
        for i in range(0, mesh.GetNumberOfCells()):
            v0 = mesh.GetCell(i).GetPointId(0)
            v1 = mesh.GetCell(i).GetPointId(1)
            v2 = mesh.GetCell(i).GetPointId(2)

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint( mesh.GetPoint(v0)[0], mesh.GetPoint(v0)[1] )
            ring.AddPoint (mesh.GetPoint(v1)[0], mesh.GetPoint(v1)[1] )
            ring.AddPoint( mesh.GetPoint(v2)[0], mesh.GetPoint(v2)[1] )
            ring.AddPoint( mesh.GetPoint(v0)[0], mesh.GetPoint(v0)[1] ) # add again to complete the ring.

            tpoly = ogr.Geometry(ogr.wkbPolygon)
            tpoly.AddGeometry(ring)

            feature = ogr.Feature( layer.GetLayerDefn() )
            feature.SetGeometry(tpoly)

            if variables is None:
                for j in range(0, cd.GetNumberOfArrays()):
                    name = cd.GetArrayName(j)
                    variables.append(name)

            for v in variables:
                data = cd.GetArray(v).GetTuple(i)
                feature.SetField(str(v), float(data[0]))

            if parameters is not None:
                for p in parameters:
                    data = cd.GetArray(p).GetTuple(i)
                    feature.SetField(str(p), float(data[0]))


            layer.CreateFeature(feature)


        x_min, x_max, y_min, y_max = layer.GetExtent()


        NoData_value = -9999
        x_res = int((x_max - x_min) / pixel_size)
        y_res = int((y_max - y_min) / pixel_size)

        for var in variables:
            target_ds = gdal.GetDriverByName('GTiff').Create(path[:-4]+'_'+var+'.tif', x_res, y_res, 1, gdal.GDT_Float32)
            target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
            # target_ds.SetProjection(srs)
            band = target_ds.GetRasterBand(1)
            band.SetNoDataValue(NoData_value)

            # Rasterize
            gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[0],options=['ALL_TOUCHED=TRUE',"ATTRIBUTE="+var])

        if parameters is not None:
            for p in parameters:
                target_ds = gdal.GetDriverByName('GTiff').Create(path[:-4] + '_' + p + '.tif', x_res, y_res, 1,
                                                                 gdal.GDT_Float32)
                target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
                band = target_ds.GetRasterBand(1)
                band.SetNoDataValue(NoData_value)

                # Rasterize
                gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[0],
                                    options=['ALL_TOUCHED=TRUE', "ATTRIBUTE=" + p])
            parameters = None

        iter += 1
# output_usm = None

if __name__ == "__main__":

    main()