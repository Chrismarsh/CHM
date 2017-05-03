import vtk
from osgeo import gdal,ogr,osr
import math
import imp
import sys
import xml.etree.ElementTree as ET
import subprocess
import numpy as np
from gdalconst import GA_ReadOnly
gdal.UseExceptions()  # Enable errors

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
    
    #####  load user configurable paramters here    #######
    # Check user defined configuraiton file
    if len(sys.argv) == 1:
        print('ERROR: main.py requires one argument [configuration file] (i.e. python main.py vtu2geo_config.py)')
        return

    # Get name of configuration file/module
    configfile = sys.argv[-1]

    # Load in configuration file as module
    X = imp.load_source('',configfile)

    # Grab variables
    input_path = X.input_path

    variables = X.variables
    parameters = []
    if hasattr(X,'parameters'):
        parameters = X.parameters

    # Check if we want to constrain output to a example geotif
    constrain_flag = False
    if hasattr(X,'constrain_tif_file'):
        constrain_tif_file = X.constrain_tif_file
        var_resample_method = X.var_resample_method
        param_resample_method = X.param_resample_method
        constrain_flag = True

    output_path = input_path[:input_path.rfind('/')+1]
    if hasattr(X,'output_path'):
        output_path = X.output_path

    pixel_size = 0
    if hasattr(X,'pixel_size'):
        pixel_size = X.pixel_size

    user_define_extent = False
    if hasattr(X,'user_define_extent'):
        user_define_extent = X.user_define_extent

    # Get size for first rasterization (less than triangle min area)
    pixel_size = X.pixel_size

    #####
    reader = vtk.vtkXMLUnstructuredGridReader()
    pvd = ET.parse(input_path)
    pvd = pvd.findall(".//*[@file]")

    # Get info for constrained output extent/resolution if selected
    if(constrain_flag):

        ex_ds = gdal.Open(constrain_tif_file,GA_ReadOnly)
        gt = ex_ds.GetGeoTransform()
        pixel_width = np.abs(gt[1])
        pixel_height = np.abs(gt[5])
        # Take extent from user input
        if user_define_extent:
            o_xmin = X.o_xmin
            o_xmax = X.o_xmax
            o_ymin = X.o_ymin
            o_ymax = X.o_ymax
        else: # Get extent for clipping from input tif
            o_xmin = gt[0]
            o_ymax = gt[3]
            o_xmax = o_xmin + gt[1] * ex_ds.RasterXSize
            o_ymin = o_ymax + gt[5] * ex_ds.RasterYSize

        print "Output pixel size is " + str(pixel_width) + " by " + str(pixel_height)
        ex_ds = None



    iter=1
    for vtu in pvd:

        vtu_file  = vtu.get('file')
        path = input_path[:input_path.rfind('/')+1] + vtu_file
        printProgress(iter,len(pvd))
        reader.SetFileName(path)
        reader.Update()

        mesh = reader.GetOutput()

        #default the pixel size to (min+max)/2
        if not pixel_size:
            area_range = mesh.GetCellData().GetArray('Area').GetRange()
            pixel_size = (math.sqrt(area_range[0]) + math.sqrt(area_range[1]))/2
            pixel_size = int( math.ceil(pixel_size) )

        driver = ogr.GetDriverByName('Memory')
        output_usm = driver.CreateDataSource('out')

        srsin = osr.SpatialReference()

        if not mesh.GetFieldData().HasArray("proj4"):
            print "VTU file does not contain a proj4 field"
            return -1

        vtu_proj4 = mesh.GetFieldData().GetAbstractArray("proj4").GetValue(0)
        srsin.ImportFromProj4(vtu_proj4)

        is_geographic = srsin.IsGeographic()

        #output conic equal area for geotiff
        srsout = osr.SpatialReference()

        if not is_geographic:
            srsout.ImportFromProj4(vtu_proj4)
        else:
            srsout.ImportFromWkt("PROJCS[\"North_America_Albers_Equal_Area_Conic\",     "
                                 "GEOGCS[\"GCS_North_American_1983\",         "
                                 "DATUM[\"North_American_Datum_1983\",            "
                                 " SPHEROID[\"GRS_1980\",6378137,298.257222101]],         "
                                 "PRIMEM[\"Greenwich\",0],        "
                                 " UNIT[\"Degree\",0.017453292519943295]],     "
                                 "PROJECTION[\"Albers_Conic_Equal_Area\"],     "
                                 "PARAMETER[\"False_Easting\",0],    "
                                 " PARAMETER[\"False_Northing\",0],     "
                                 "PARAMETER[\"longitude_of_center\",-96],     "
                                 "PARAMETER[\"Standard_Parallel_1\",20],     "
                                 "PARAMETER[\"Standard_Parallel_2\",60],     "
                                 "PARAMETER[\"latitude_of_center\",40],     "
                                 "UNIT[\"Meter\",1],     "
                                 "AUTHORITY[\"EPSG\",\"102008\"]]")

        trans = osr.CoordinateTransformation(srsin,srsout)

        layer = output_usm.CreateLayer('poly', srsout, ogr.wkbPolygon)

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
            if is_geographic:
                scale = 100000.  #undo the scaling that CHM does for the paraview output
                ring.AddPoint( mesh.GetPoint(v0)[0] / scale, mesh.GetPoint(v0)[1] / scale)
                ring.AddPoint (mesh.GetPoint(v1)[0]/ scale, mesh.GetPoint(v1)[1]/ scale )
                ring.AddPoint( mesh.GetPoint(v2)[0]/ scale, mesh.GetPoint(v2)[1]/ scale )
                ring.AddPoint( mesh.GetPoint(v0)[0]/ scale, mesh.GetPoint(v0)[1]/ scale ) # add again to complete the ring.

                ring.Transform(trans) #only transform if needed

            else:
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
            target_ds = gdal.GetDriverByName('GTiff').Create(output_path+'/'+vtu_file[:-4] + '_'+ var.replace(" ","_")+'.tif', x_res, y_res, 1, gdal.GDT_Float32)
            target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
            # target_ds.SetProjection(srs)
            band = target_ds.GetRasterBand(1)

            band.SetNoDataValue(NoData_value)

            # Rasterize
            gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[0],options=['ALL_TOUCHED=TRUE',"ATTRIBUTE="+var])
            target_ds.SetProjection(srsout.ExportToWkt())
            target_ds = None

            # Optional clip file
            if(constrain_flag):
                subprocess.check_call(['gdalwarp -overwrite -s_srs \"%s\" -t_srs \"%s\" -te %f %f %f %f -r \"%s\" -tr %f %f \"%s\" \"%s\"' %
                                       (srsout.ExportToProj4(),srsout.ExportToProj4(),o_xmin, o_ymin, o_xmax, o_ymax, var_resample_method[var], pixel_width, pixel_height, path[:-4]+'_'+var+'.tif',path[:-4]+'_'+var+'_clipped.tif')], shell=True)

            if parameters is not None:
                for p in parameters:
                    target_ds = gdal.GetDriverByName('GTiff').Create(output_path+'/'+vtu_file[:-4] + '_'+ p.replace(" ","_") + '.tif', x_res, y_res, 1, gdal.GDT_Float32)
                    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
                    target_ds.SetProjection(srsout.ExportToWkt())
                    band = target_ds.GetRasterBand(1)
                    band.SetNoDataValue(NoData_value)

                    # Rasterize
                    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[0],options=['ALL_TOUCHED=TRUE', "ATTRIBUTE=" + p])
                    target_ds = None

                    # Optional clip file
                    if(constrain_flag):
                        subprocess.check_call([
                          'gdalwarp -overwrite -s_srs \"%s\" -t_srs \"%s\" -te %f %f %f %f -r \"%s\" -tr %f %f \"%s\" \"%s\"' %
                          (srsout.ExportToProj4(), srsout.ExportToProj4(), o_xmin, o_ymin, o_xmax,
                           o_ymax, param_resample_method[p], pixel_width, pixel_height,
                           output_path + '/' + vtu_file[:-4] + '_' + p.replace(" ","_") + '.tif',
                           output_path + '/' + vtu_file[:-4] + '_' + p.replace(" ","_") + '_clipped.tif')],
                          shell=True)

        # we don't need to dump parameters for each timestep as they are currently assumed invariant with time.
        parameters = None

        #no parameters and no variables, just exit at this point
        if not variables and parameters is None:
            break

        iter += 1

if __name__ == "__main__":

    main()
