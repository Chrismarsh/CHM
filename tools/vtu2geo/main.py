import vtk
from osgeo import gdal,ogr,osr
from osgeo.gdalconst import GA_ReadOnly
import math
import imp
import sys
import xml.etree.ElementTree as ET
import subprocess
import numpy as np
import os

import resource

from rasterio.warp import transform
import xarray as xr

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

import contextlib


#suppress stderr from one of the vtk calls that has deprecated warning out put for something outside of our control with vtk8.1
#https://stackoverflow.com/a/17753620/410074
@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr

    e.g.:


    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """

    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()

# layer Input USM to be rasterized
# srsout Output CRS
# the filename of the raster
# pixel size, in srsout units. Constant dx,dy
# attribute The value from the USM we want to burn into the raster
def rasterize(layer, srsout, target_fname, pixel_size, attribute, user_define_extent=False, o_xmin=None, o_ymin=None, o_xmax=None, o_ymax=None, all_touched=True):
    x_min, x_max, y_min, y_max = layer.GetExtent()

    NoData_value = -9999
    # width/height of created raster in pixels.
    xsize = int((x_max - x_min) / pixel_size)
    ysize = int((y_max - y_min) / pixel_size)


    target_ds = gdal.GetDriverByName('GTiff').Create(target_fname, xsize, ysize, 1, gdal.GDT_Float32)
    target_ds.SetProjection(srsout.ExportToWkt())
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(NoData_value)

    # Rasterize. If all_touched is false, then only the cell centres that are within the triangle are used
    # if lower resolution rasterization is used, then this can miss cells. When all touched is true, as long as a triangle
    # touches a cell, then it's counted
    if all_touched:
        gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[0],
                            options=['ALL_TOUCHED=TRUE',"ATTRIBUTE=" + attribute])
    else:
        gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[0], options=[ "ATTRIBUTE=" + attribute]) ,

    target_ds = None

    # Optional clip file
    if user_define_extent: # or constrain_flag:
        subprocess.check_call(
            ['gdalwarp -overwrite -s_srs \"%s\" -t_srs \"%s\" -te %s %s %s %s -r \"%s\" -tr %s %s \"%s\" \"%s\"' %
             (srsout.ExportToProj4(),
              srsout.ExportToProj4(),
              o_xmin, o_ymin, o_xmax, o_ymax,
              'average', #var_resample_method[var]
              pixel_size, pixel_size,
              target_fname,
              target_fname[:-4] + '_clipped.tif')], shell=True)

def main():
    gdal.UseExceptions()  # Enable errors

    #####  load user configurable paramters here    #######
    # Check user defined configuraiton file
    if len(sys.argv) == 1:
        print('ERROR: main.py requires one argument [configuration file] (i.e. python main.py vtu2geo_config.py)')
        return

    # Get name of configuration file/module
    configfile = sys.argv[1]

    # Load in configuration file as module
    X = imp.load_source('',configfile)


    # if a 2nd command line argument is present, it is the input_path so use that, otherwise try to use the one from passed script
    input_path = ''
    if len(sys.argv) == 3: # we have a 2nd CLI arg
        input_path = sys.argv[2]
        if hasattr(X,'input_path'):
            print('Warning: Overwriting script defined input path with CL path')
    elif hasattr(X,'input_path'):
        input_path = X.input_path
    else:
        print('ERROR: No input path. A pvd or vtu file must be specified.')
        exit(-1)

    if os.path.isdir(input_path):
        print('ERROR: Either a pvd or vtu file must be specified.')
        exit(-1)

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

        if not os.path.isdir(output_path):
            os.makedirs(output_path)

    pixel_size = 10
    if hasattr(X,'pixel_size'):
        pixel_size = X.pixel_size
    else:
        print('Default pixel size of 10 mx10 m will be used.')

    user_define_extent = False
    if hasattr(X,'user_define_extent'):
        user_define_extent = X.user_define_extent

    #Produces a lat/long regular grid in CF netCDF format instead of the TIF files
    nc_archive = False
    if hasattr(X,'nc_archive'):
        nc_archive = X.nc_archive

    #user defined output EPSG to use instead of the proj4 as defined in the vtu
    out_EPSG = None
    if hasattr(X,'out_EPSG'):
        out_EPSG = X.out_EPSG

    if parameters is not None and nc_archive:
        print('Parameters are ignored when writing the nc archive.')
        parameters = []

    all_touched = True
    if hasattr(X,'all_touched'):
        all_touched = X.all_touched

    keep_tifs = False
    if hasattr(X,'keep_tifs'):
        keep_tifs = X.keep_tifs

    #####
    reader = vtk.vtkXMLUnstructuredGridReader()

    # see if we were given a single vtu file or a pvd xml file
    filename, file_extension = os.path.splitext(input_path)
    is_pvd = False
    pvd = [input_path] # if not given a pvd file, make this iterable for the below code
    timesteps=None
    if file_extension == '.pvd':
        print('Detected pvd file, processing all linked vtu files')
        is_pvd = True
        if not os.path.exists(input_path):
            print(f'ERROR: Path does not exists: {input_path}')
            exit(-1)

        parse = ET.parse(input_path)
        pvd = parse.findall(".//*[@file]")

        timesteps = parse.findall(".//*[@timestep]")

    o_xmin = o_xmax = o_ymin = o_ymax = None
    if user_define_extent:
        o_xmin = X.o_xmin
        o_xmax = X.o_xmax
        o_ymin = X.o_ymin
        o_ymax = X.o_ymax

    # Get info for constrained output extent/resolution if selected
    if constrain_flag :
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

        print(("Output pixel size is " + str(pixel_width) + " by " + str(pixel_height)))
        ex_ds = None

    if constrain_flag:
        print(" Constrain flag currently not supported!")
        return -1


    files_processed=1 # this really should be 1 for useful output

    #information to build up the nc meta data
    nc_rasters = {}
    if nc_archive:
        for v in variables:
            nc_rasters[v]=[]

    nc_time_counter=0
    tifs_to_remove = []
    epoch=np.datetime64(0,'s')

    # if we are loading a pvd, we have access to the timestep information if we want to build
    if is_pvd and timesteps is not None:
        epoch = np.datetime64(int(timesteps[0].get('timestep')),'s')

    dt=1 # model timestep, in seconds

    if timesteps is not None and len(timesteps) > 1:
        dt = int(timesteps[1].get('timestep')) - int(timesteps[0].get('timestep'))


    print(('Start epoch: %s, model dt = %i (s)' %(epoch,dt)))

    #because of how the netcdf is built we hold a file of file handles before converting. ensure we can do so
    if nc_archive:
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        total_output_files = len(pvd) * (len(variables)+len(parameters))
        if soft < total_output_files or hard < total_output_files:
          
            print('The users soft or hard file limit is less than the total number of tmp files to be created.')
            print('The system ulimit will be raised to at least ' + str(total_output_files))

            try:
                resource.setrlimit(resource.RLIMIT_NOFILE, (total_output_files,resource.RLIM_INFINITY ) )
            except ValueError as e:
                print("Failed to raise the limit")
                raise e


    for vtu in pvd:
        path = vtu
        vtu_file = ''
        if is_pvd:
            vtu_file  = vtu.get('file')
            path = input_path[:input_path.rfind('/')+1] + vtu_file
        else:
            base = os.path.basename(path) # since we have a full path to vtu, we need to get just the vtu filename
            vtu_file = os.path.splitext(base)[0] #we strip out vtu later so keep here

        if not os.path.exists(path):
            print(f'ERROR: Path does not exists: {path}')
            exit(-1)

        printProgress(files_processed,len(pvd),decimals=0)
        reader.SetFileName(path)
        reader.Update()
        #shut up a deprecated warning from vtk 8.1
        # with stdchannel_redirected(sys.stderr, os.devnull):
        #     reader.Update()

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
            print("VTU file does not contain a proj4 field")
            return -1

        vtu_proj4 = mesh.GetFieldData().GetAbstractArray("proj4").GetValue(0)
        srsin.ImportFromProj4(vtu_proj4)

        is_geographic = srsin.IsGeographic()

        #output conic equal area for geotiff
        srsout = osr.SpatialReference()
        srsout.ImportFromProj4(vtu_proj4)

        if out_EPSG:
            srsout.ImportFromEPSG(out_EPSG)

        trans = osr.CoordinateTransformation(srsin,srsout)
        layer = output_usm.CreateLayer('poly', srsout, ogr.wkbPolygon)
        cd = mesh.GetCellData()

        for i in range(0,cd.GetNumberOfArrays()):
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
            else:
                ring.AddPoint( mesh.GetPoint(v0)[0], mesh.GetPoint(v0)[1] )
                ring.AddPoint (mesh.GetPoint(v1)[0], mesh.GetPoint(v1)[1] )
                ring.AddPoint( mesh.GetPoint(v2)[0], mesh.GetPoint(v2)[1] )
                ring.AddPoint( mesh.GetPoint(v0)[0], mesh.GetPoint(v0)[1] ) # add again to complete the ring.

            ring.Transform(trans)
            tpoly = ogr.Geometry(ogr.wkbPolygon)
            tpoly.AddGeometry(ring)

            feature = ogr.Feature( layer.GetLayerDefn() )
            feature.SetGeometry(tpoly)

            if variables is None:
                variables = []
                for j in range(0, cd.GetNumberOfArrays()):
                    name = cd.GetArrayName(j)
                    if not '[param]' in name:
                        variables.append(name)

            for v in variables:
                try:
                    data = cd.GetArray(v).GetTuple(i)
                    feature.SetField(str(v), float(data[0]))
                except:
                    print(("Variable %s not present in mesh" % v))
                    return -1

            if parameters is not None:
                for p in parameters:
                    try:
                        data = cd.GetArray(p).GetTuple(i)
                        feature.SetField(str(p), float(data[0]))
                    except:
                        print(("Parameter %s not present in mesh" % v))
                        return -1
            layer.CreateFeature(feature)

        for var in variables:
            print(var)
            target_fname = os.path.join(output_path,
                                        vtu_file + '_' + var.replace(" ", "_") + str(pixel_size) + 'x' + str(
                                            pixel_size) + '.tif')
            rasterize(layer, srsout, target_fname, pixel_size, var, user_define_extent, o_xmin, o_ymin, o_xmax, o_ymax, all_touched)

            if nc_archive:
                df = xr.open_rasterio(target_fname).sel(band=1).drop('band')
                df = df.rename({'x':'lon','y':'lat'})
                df.coords['time']=epoch + nc_time_counter*np.timedelta64(dt,'s') # this will automatically get converted to min or hours in the output nc
                df.name=var
                nc_rasters[var].append(df) # these are lazy loaded at the to netcdf call

                # remove the tifs we produce
                tifs_to_remove.append(target_fname)


        if parameters is not None:
            print(parameters)
            for p in parameters:
                target_param_fname = os.path.join(output_path, vtu_file + '_'+ p.replace(" ","_") + str(pixel_size)+'x'+str(pixel_size)+'.tif')
                rasterize(layer, srsout, target_param_fname, pixel_size, p, all_touched = all_touched)

        nc_time_counter += 1
        # we don't need to dump parameters for each timestep as they are currently assumed invariant with time.
        parameters = None

        #no parameters and no variables, just exit at this point
        if not variables and parameters is None:
            break

        files_processed += 1

    if nc_archive:

        datasets = []

        for var, rasters in list(nc_rasters.items()):
            a = xr.concat(rasters,dim='time')
            datasets.append(a.to_dataset())

        arr = xr.merge(datasets)
        print('Writing netCDF file')
        fname=os.path.join(os.path.splitext(input_path)[0]+'.nc')
        arr.to_netcdf(fname,engine='netcdf4')

        for f in tifs_to_remove:
            try:
                os.remove(f)
            except:
                pass

if __name__ == "__main__":

    main()
