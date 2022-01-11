# -*- coding: utf-8 -*-
# ----------------------------------
# Code to calculate ILCI on all BAS PASIN data 

# This code reads in all of BAS PASIN data published on BAS Polar Aerogeophysics
# Data Portal and calculates the ILCI for the chirp data.

# Details on ILCI methodology can be found in Karlsson et al. 
# (2012, EPSL, doi: 10.1016/j.epsl.2012.04.034)

# ILCI calculations are exported as two large files: 1x geopackage and 1x CSV
# file per survey, with 5x columns as follows:
# FltID (Flight ID), X (EPSG:3031), Y (EPSG:3031), ILCI, Normalised_ILCI

# As per Karlsson et al. (2012), we disregard the upper and lower fifths of the
# radar profile and smooth the data using a 1000-point average. 
# ----------------------------------

# Code written by J. Bodart (BAS/UoE) - 18/10/2021

# To cite the code + data output, refer to the associated manuscript
# Fremand and Bodart et al. (2022, ESSD, doi: xxxxxx)
# ----------------------------------

# import python libraries
# ----------------------------------
import numpy as np
import pandas as pd
import os
import netCDF4 as nc4
import glob
import geopandas as gpd
# ----------------------------------

# Get paths where NetCDF data lives on computer
#-----------------------------------------------
os.chdir('D:/British_Antarctic_Survey/data/')
surveys = sorted(glob.glob('D:/British_Antarctic_Survey/data/*/netcdf/')) # sort files in order
#-----------------------------------------------

# Loop starts here
#-----------------------------------------------
# Loop through all survey folders
for p in range(0, len(surveys)):
    
    # Initialise export of data (includes setting up output path, folderID, arrays, 
    # CSV where data will be stored, etc.)
    #------------------------------------------------
    Main_path = surveys[p]
    if 'FISS_1516' in Main_path:
        folderID = Main_path.split('/')[-6]
    else:
    	folderID = Main_path.split('\\')[1]
    
    print ('############## Importing NetCDF files from %s survey ##############' %folderID)
    output_path = "D:/British_Antarctic_Survey/Aerogeophysics_data_paper/ILCI_figure/ILCI_results/test/"
    nc_files = sorted(glob.glob(Main_path + '*.nc')) # sort files in order
    
    # initialise CSV file with 5x columns as set-up
    result = open('D:/British_Antarctic_Survey/Aerogeophysics_data_paper/ILCI_figure/ILCI_results/test/%s_ILCI.csv' %folderID, 'w+')
    result.write('FltID, X_EPSG3031, Y_EPSG3031, ILCI, Normalised_ILCI\n')
    
    # initialise array to store results into after loop
    # initialisation set-up: 1 row, 4 columns
    dfss = np.full((1,4), 0)
    #-----------------------------------------------
    
    # Loop through all NetCDF files for that specific survey
    #-----------------------------------------------
    for m in range(0,len(nc_files)):
        
        # import netcdf file
        #-----------------------------------------------
        ncin = nc4.Dataset(nc_files[m]) # open NetCDF file from path
        fltid = ncin.flight # get flight line name from netCDF
        print ('####### Importing and reading netcdf: %s #######' %fltid)
        #-----------------------------------------------
    
        # First, extract radar variable from nc files
        #-----------------------------------------------
        # AGAP, BBAS, GRADES-IMAGE, FISS 2015 & 2016, ICEGRAV, POLARGAP (non-polarised) WISE-ISODYN 
        # all have only 1x chirp
        if 'chirp_data' in list(ncin.variables):
            chirpData = ncin.variables['chirp_data'][:].data # read in chirp radar data array
 
        # Make sure you select correct chirp in NetCDF (in case of multiple)
        # or if chirp does not exist, avoid this NetCDF
        elif not 'chirp_data' in list(ncin.variables):
            if 'IMAFI' in folderID:
                chirpData = ncin.variables['chirp_DLRsar_data'][:].data # read in one of two chirp data array
            elif 'POLARGAP' in folderID and 'polarised_chirp_PPVV_data' in list(ncin.variables):
                chirpData = ncin.variables['polarised_chirp_PPVV_data'][:].data # read in in one of two chirp data array
            else:
                break # if not chirp, then break loop here and jump to next NetCDF

        # Make sure you use the correct trace variable if there are multiple (e.g. POLARGAP and AGAP)
        if 'POLARGAP' in folderID or 'AGAP' in folderID:
            traces_nc = ncin.variables['traces_chirp'][:].data # read in trace number from chirp variable
        else: # for all others
            traces_nc = ncin.variables['traces'][:].data # read in trace number (only one exists in NetCDF)
        #-----------------------------------------------
    
        # Second, import other NetCDF variables (X, Y, surf and bed picks)
        #-----------------------------------------------
        x_nc = ncin.variables['x_coordinates'][:].data # read in X positions (EPSG: 3031)
        y_nc = ncin.variables['y_coordinates'][:].data # read in Y positions (EPSG: 3031)
        
        # read in surface and bed pick from NetCDF file and convert -9999 to NaN
        surf_pick = ncin.variables['surface_pick_layerData'][:].data # read in surf pick array
        bed_pick = ncin.variables['bed_pick_layerData'][:].data # read in bed pick array
        surf_pick[surf_pick == -9999] = 'nan' # convert -9999 to NaNs
        bed_pick[bed_pick == -9999] = 'nan' # convert -9999 to NaNs
        #-----------------------------------------------
        
        # If a flightline only contains NaNs in surface or bed pick, ILCI cannot
        # be calculated, thus skip flightline entirely if this occurs
        #-----------------------------------------------
        if all(np.isnan(surf_pick)) or all(np.isnan(bed_pick)):
            print('Flightline %s only contains NaN - excluded from calculations' %fltid)
            pass
    
        # If there is valid values, then continue with the calculations
        else:
            
            # Remove data where there are NaNs in the surface pick variable
            #-----------------------------------------------
            nan_idx = np.where(~np.isnan(surf_pick)) # get index of non-NaN values
            x_nc = x_nc[nan_idx] # remove NaN in X variable
            y_nc = y_nc[nan_idx] # remove NaN in Y variable
            traces_nc = traces_nc[nan_idx] # remove NaN in trace variable
            chirpData = chirpData[:, nan_idx[0]]
            surf_pick = surf_pick[nan_idx] # remove NaN in surface pick variable
            bed_pick = bed_pick[nan_idx] # remove NaN in bed pick variable
            #-----------------------------------------------
            
            # split data into 8000-trace shunks
            #-----------------------------------------------
            nTrace = len(traces_nc) # read length of trace variable
            idx = np.arange(1,nTrace,8000) # get index of each 8000-trace segments
            dfs = np.full((1,4), 0) # initialisation set-up: 1 row, 4 columns
            new_df = []
            
            # loop over each index enveloppe
            print ('Calculating the indices for the data')
            for i in range(1, len(idx)):
                
                # calculate indices for all values that fall within 'idx'
                # variable above for last array in record
                #-----------------------------------------------
                # for last array in record
                if i == len(idx)-1:
                    boolean = (traces_nc>idx[i])
                    seg = np.where(boolean)[0]
                    indices = [(idx[i]), nTrace]
                    
                # for all arrays except last array in record
                else:
                    boolean = (traces_nc>=idx[i]) & (traces_nc < idx[i+1])
                    seg = np.where(boolean)[0]
                    indices = [(idx[i]), (idx[i+1])]
                    
                # select trace number, radar data, x/y positions,
                # and surf/bed picks based on 'idx' variable above
                traces_seg = traces_nc[seg]
                chirp_seg = chirpData[:,seg]
                surfPick_seg = surf_pick[seg]
                bedPick_seg = bed_pick[seg]
                x_seg = x_nc[seg]
                y_seg = y_nc[seg]
                
                # create array with name of flightline as string
                fltID_arr = [fltid for i in range(len(seg))]
                
                # smooth chirp_data using a window of 1000 traces (500 either side) 
                print ('Smoothing the chirp data for index: %s/%s' %(i, len(idx)-1))
                data = chirp_seg
                kk = np.arange(499,len(chirp_seg[1])-498)
                for j in np.arange(499,len(chirp_seg[1])-498):
                    data[:,j]  = chirp_seg[:,j-498:j+498].mean(axis=1)
                #-----------------------------------------------
            
                # Calculate ILCI
                #-----------------------------------------------
                print ('Calculating ILCI and normalising results')
                ilci_data = [] # initialize array
                
                for l in range(0, len(data[1])):
                               
                    # create log of data
                    data_log = np.log(data[:,l])
                    
                    # round surf/bed pick to nearest integer if not a nan
                    # otherwise leave as nan
                    if np.isnan(surfPick_seg[l]) or np.isnan(bedPick_seg[l]):
                        spick = surfPick_seg[l]
                        bpick = bedPick_seg[l]
                    else:
                        spick = int(np.round(surfPick_seg[l]))
                        bpick = int(np.round(bedPick_seg[l]))
                    #-----------------------------------------------
        
                    # exclude first and last fifths of radar data
                    #-----------------------------------------------
                    # if bed-surf pick is smaller than 20 (e.g. too close to call)
                    # or bed pick is larger than y-axis of radar data (happens...)
                    # then place nan. Otherwise, calculate ILCI
                    if bpick-spick < 20 or bpick > len(data) or np.isnan(bpick) or np.isnan(spick):
                        ilci = np.nan # place nan
                    else:
                        good = data_log[spick:bpick]
                        bad = int(np.round(len(good)/5)) # length of data to cut
                        data_cut = good[bad:-bad] # array without 1st and last fifth
                        #-----------------------------------------------
            
                        # calculate ILCI
                        #-----------------------------------------------
                        ilci = np.nanmean(abs(np.gradient(data_cut)))
                    
                    # append array after each iteration in the loop
                    ilci_data.append(ilci)
                    #-----------------------------------------------
                        
                # append data to large array
                new_df = np.column_stack([fltID_arr, x_seg,y_seg,ilci_data])
                new_df = new_df[499:-498] # cut first and last 500 to match radar smoothing window
                    
                # concatenate all results        
                dfs = np.concatenate((dfs,new_df),axis=0)
                #-----------------------------------------------
        
            # combine each flightline data into a large data array
            #-----------------------------------------------
            dfss = np.concatenate((dfss,dfs), axis=0)
            #-----------------------------------------------
    
    # Place data into panda dataframe and drop any consecutive zeros at the start
    # arising from initialisation of arrays in code above
    #----------------------------------------------- 
    # only put last 3x columns into the dataframe initially (i.e. not Flight ID yet)
    results = pd.DataFrame(dfss[:,1:], columns=['X_EPSG3031', 'Y_EPSG3031', 'ILCI'], dtype=float)
    
    # find amount of zeros to delete at start of record due to array initialisation
    # (ugly way)
    remove_zeros = len(results.loc[(results==0).all(axis=1)])
    
    # remove zeros at start of record in dataframe
    results = results.loc[~(results==0).all(axis=1)] # remove zeros at start of record
    
    # insert flight ID column into dataframe - make sure those zeros are excluded
    results.insert(0,"FltID", dfss[remove_zeros:,0])
    #-----------------------------------------------

    # Normalise ILCI between 0 and 1 for entire survey and attach new column to dataframe
    # (optional)
    #-----------------------------------------------
    normalised_ILCI = (results['ILCI'] - np.min(results['ILCI'])) / (np.max(results['ILCI']) - np.min(results['ILCI'])) # normalise ILCI
    results["Norm_ILCI"] = normalised_ILCI # attach to dataframe as last column
    #-----------------------------------------------
    
    # Write dataframe to point shapefile and CSV file
    #-----------------------------------------------
    # Convert dataframe to geopandas GeoDataFrame
    gdf = gpd.GeoDataFrame(results, geometry=gpd.points_from_xy(results['X_EPSG3031'], results['Y_EPSG3031']))#, dtype='float') 
    gdf = gdf.set_crs("EPSG:3031") #Set the coordinate system
    schema = gpd.io.file.infer_schema(gdf) # extract schema of dataframe

    # make sure all other variables after flight ID are floats
    vars_gdf = gdf.columns.tolist()
    vars_gdf.remove('geometry')
    for k in range(1, len(vars_gdf)):
        schema['properties'][vars_gdf[k]] = 'float'
    
    # Export to point shapefile
    print ('Writing data to shapefile and CSV file')
    #gdf.to_file(output_path + folderID + '_ILCI.shp', driver="ESRI Shapefile",  schema=schema) # To save as a shapefile
    gdf.to_file(output_path + folderID + '_ILCI.gpkg', driver="GPKG",  schema=schema) # To save as a geopackage

    # Export to csv file
    results.to_csv(output_path + folderID + '_ILCI.csv',index=False)
    print ('Done with %s survey - now onto next survey' %folderID)
    #-----------------------------------------------

    
        