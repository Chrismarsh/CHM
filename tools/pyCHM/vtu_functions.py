# Functions to work with VHM .vtu files

# Get mesh (returns vtk file)
def get_mesh(vtu_file):
    import vtk
    # Read with vtk
    reader = vtk.vtkXMLUnstructuredGridReader()
    # Read in the mesh
    reader.SetFileName(vtu_file)
    reader.Update()
    mesh = reader.GetOutput()
    return mesh

# get datetime from vtu file name
def get_vtu_time(vtu_files,prefix):
    import datetime
    import re
    import pandas as pd 
    DT = []
    for cf in vtu_files:
        ctimestamp = float(re.split(prefix,cf)[1].split('.')[0])
        DT.append(datetime.datetime.utcfromtimestamp(ctimestamp))
    return pd.DatetimeIndex(DT)


# Get mesh coord info (returns dictionary)
def get_triangle(mesh):
    from vtk.util import numpy_support as VN
    # Get coords
    points =  VN.vtk_to_numpy(mesh.GetPoints().GetData())
    
    # Get location of nodes
    X = [cpt[0] for cpt in points] 
    Y = [cpt[1] for cpt in points]
    E = [cpt[2] for cpt in points]
    
    # Build face adjacency index
    triang=[]
    for i in range(0, mesh.GetNumberOfCells()):
        v0 = mesh.GetCell(i).GetPointId(0)
        v1 = mesh.GetCell(i).GetPointId(1)
        v2 = mesh.GetCell(i).GetPointId(2)

        triang.append( [v0,v1,v2])
        
    return {'X':X,'Y':Y,'E':E,'triang':triang}


# In[5]:

# Save triangle info
def save_triangle_info(c_mesh,tri_file):
    import numpy as np
    c_triang = get_triangle(get_mesh(c_mesh))
    np.save(tri_file, c_triang) 


# In[6]:

# Grab face array (returns numpy array)
def get_face_var(mesh,var_name):
    from vtk.util import numpy_support as VN
    cvar = VN.vtk_to_numpy(mesh.GetCellData().GetArray(var_name))
    return cvar


# In[7]:

# Get many mesh (returns vtk file)
def get_multi_mesh(vtu_files):
    meshes = []
    for cf in vtu_files:
        meshes.append(get_mesh(cf))
    return meshes


# In[8]:

# Get dataframe of variable (time x face)
def get_multi_mesh_var(meshes,var_name,time_stamp):
    import pandas as pd
    # Grab first var to get length
    temp = get_face_var(meshes[0],var_name)
    # Initilize dataframe
    df = pd.DataFrame(index=time_stamp,columns=np.arange(len(temp)))
    for cf,ct in zip(meshes,time_stamp):
        cnp = get_face_var(cf,var_name)
        df.loc[ct] = cnp
    return df


# In[9]:

# Get dataframe of variable (time x face)
def get_multi_mesh_var_2(vtu_files,var_name,time_stamp):
    # Grab first var to get length
    temp_mesh = get_mesh(vtu_files[0])
    temp_var  = get_face_var(temp_mesh,var_name)
    temp_mesh = None
    # Initilize dataframe
    df = pd.DataFrame(index=time_stamp,columns=np.arange(len(temp_var)))
    for cf,ct in zip(vtu_files,time_stamp):
        print cf
        c_mesh = get_mesh(cf) # current mesh
        df.loc[ct] = get_face_var(c_mesh,var_name) # current variable
        c_mesh = None
    return df


# In[10]:

# With dask
# Get dataframe of variable (time x face)
def get_multi_mesh_var_dask(vtu_files,var_name,time_stamp):
    from dask.delayed import delayed
    import time
    import pandas as pd
    # Grab first var to get length
    temp_mesh = get_mesh(vtu_files[0])
    temp_var  = get_face_var(temp_mesh,var_name)
    temp_mesh = None
    # Initilize dic of variable arrays
    ar = {}
    # Loop through each vtu file
    for cf,ct in zip(vtu_files,time_stamp):
        c_mesh = delayed(get_mesh)(cf) # current mesh
        cvar   = delayed(get_face_var)(c_mesh,var_name) # current variable
        ar[ct] = cvar
    #df.visualize() # not working, need to install dot.exe for windows path
    df = delayed(pd.DataFrame.from_dict)(ar) # ,orient='index'
    start_time = time.time()
    df_out = df.compute()
    df = None
    print("--- Took %s minutes ---" % ((time.time() - start_time)/60))
    return df_out


# In[11]:

# Import and save to hdf/netcdf multi variables from vtu files
def vtu_to_hdf(all_vars,vtu_files,time_stamp,hdf_dir):
    
    # For each variable
    for cvar in all_vars:
        
        # Get dataframe of current variable
        c_df = get_multi_mesh_var_dask(vtu_files,cvar,time_stamp).T # Transpose is needed to make time index
    
        # Save dataframe to hdf
        hdf_file = hdf_dir + '/' + cvar + '.hdf'
        c_df.to_hdf(hdf_file,'w')
        print 'Saved ' + hdf_file
        c_df = None


