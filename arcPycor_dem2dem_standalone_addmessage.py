'''
pycor_dem2dem.py

Description: Coregister slave DEM to master DEM using Nuth-Kaeaeb algorithm.
Reference: Nuth, C. and Kaeaeb, A. (2011): http://doi.org/10.5194/tc-5-271-2011

Inputs:
(1) master DEM
(2) slave DEM
(3) polygon shapefile of stable terrain

Outputs:
(1) x_bin (.csv; aspect)
(2) y_bin (.csv; dh/tan(slope))
(3) sigma_bin (.csv; standard deviation of y)
(4) shiftVec txt (shift vector)

History:
* 2014-11-8, Created
* 2019-5-7, path management
* 2019-5-11, delete unused variables, optimize the codes for converging
* 2023-5-26, set snap raster, test on ArcMap 10.8

Author: Yuzhe Wang
E-mail: yuzhe.wang@foxmail.com
Affiliation: Colledge of Geography and Environment, Shandong Normal University
'''

import os
import time
import shutil
import numpy as np
from scipy.optimize import curve_fit
import arcpy
from arcpy import env
from arcpy.sa import Raster, Slope, Aspect, ExtractByMask

startTime = time.clock()

if arcpy.CheckExtension("spatial")=="Available":
    arcpy.CheckOutExtension("spatial")
else:
    raise LicenseError

# Set environment workspace in current directory
path_script = os.path.dirname(os.path.abspath(__file__))
env.workspace = os.path.join(path_script, 'benchmark_data')

# Folder for outputs
dirOutputs = os.path.join(env.workspace, 'outputs')
if os.path.exists(dirOutputs):
    shutil.rmtree(dirOutputs)
os.makedirs(dirOutputs)

# Specify master DEM, slave DEM, stable terrain and corrected DEM
DEM_master = "DEM_master.tif"
DEM_slave = "DEM_slave4.tif"
stable_terrain = "OffGlacier.shp"

corrected_DEM = "DEM_slave_correct.tif"

# Initializations
iteration = 0
DEM_slave_before = DEM_slave
DEM_slave_after = DEM_slave
result_mean_dh = [0]
result_std_dh = [0]
ShiftX = [0]
ShiftY = [0]
file_shiftVec = os.path.join(dirOutputs, "shiftVec" + ".csv")

# Define CosineFitting function
def CosineFitting(x, a, b, c):
    return a*np.cos(b - np.pi/180*x) + c

while 1:
    iteration += 1
    print "--------------------------------------------------------------"
    print "Iteration {0} is running!".format(iteration)
    
    # DEM difference [m]
    dh = Raster(DEM_master) - Raster(DEM_slave_after)
        
    if iteration==1:
        dh.save(os.path.join(dirOutputs, "dh_init"))
    
    # slope of the slave DEM [degree]
    slp = Slope(DEM_slave_after, "DEGREE", "1")
    
    # aspect of the slave DEM [degree]
    asp = Aspect(DEM_slave_after)
    
    # Mask 'dh' using statale terrain polygon
    dh_mask = ExtractByMask(dh, stable_terrain)    
    
    # Mask 'slp' and 'asp' using 'dh_mask' in order to keep same georeference as 'dh_mask'
    env.snapRaster = dh_mask
    slp_mask = ExtractByMask(slp, dh_mask)
    asp_mask = ExtractByMask(asp, dh_mask)
    
    env.snapRaster = None
    # delete dh, slp, and asp after extraction
    del dh, slp, asp

    # Raster to Array
    dh_mask_arr = arcpy.RasterToNumPyArray(dh_mask, nodata_to_value=-32768)
    slp_mask_arr = arcpy.RasterToNumPyArray(slp_mask, nodata_to_value=-32768)
    asp_mask_arr = arcpy.RasterToNumPyArray(asp_mask, nodata_to_value=-32768)
    
    del dh_mask, slp_mask, asp_mask
       
    # Criteria: |dh| < 70 m and 5 < slope < 45.
    # slope>5 is cited from Purinton&Bookhagen, 2018, Earth Surface Dynamics.
    # slope<45 is cited from Berthier et al., 2019, Journal of Glaciology.
    index1 = np.where((dh_mask_arr>-70) & (dh_mask_arr<70) & (slp_mask_arr>5) & (slp_mask_arr<45))
    dh_mask_arr1 = dh_mask_arr[index1[0], index1[1]]
    slp_mask_arr1 = slp_mask_arr[index1[0], index1[1]]
    asp_mask_arr1 = asp_mask_arr[index1[0], index1[1]]
    
    del dh_mask_arr, slp_mask_arr, asp_mask_arr, index1
    
    # statistics of "dh"
    mean_dh_mask1 = np.mean(dh_mask_arr1)
    result_mean_dh.append(mean_dh_mask1)
    
    std_dh_mask1 = np.std(dh_mask_arr1)    
    result_std_dh.append(std_dh_mask1)
    
    print "MEAN(dh) of iteration {0}: {1:.1f}".format(iteration, mean_dh_mask1)
    print "STD(dh) of iteration {0}: {1:.1f}".format(iteration, std_dh_mask1)
    
    # Prepare the x and y values for curve fitting
    x = asp_mask_arr1
    y = dh_mask_arr1 / np.tan(np.pi*slp_mask_arr1/180)

    del dh_mask_arr1, slp_mask_arr1, asp_mask_arr1
    
    # Get the x, y values in bins
    range_asp = range(0, 370, 10)
    n = len(range_asp) - 1
    x_bin = np.zeros(n)
    y_bin = np.zeros(n)
    sigma_bin = np.zeros(n)
    for i in range(n):
        index2 = np.where((x >= range_asp[i]) & (x < range_asp[i+1]))
        x_bin[i] = range_asp[i] + 5
        y_bin[i] = np.median(y[index2])
        sigma_bin[i] = np.std(y[index2])

    # Save the results (x, y, x_bin, y_bin, sigma_bin) in csv format
    file_x = os.path.join(dirOutputs, "x" + str(iteration) + '.csv')
    file_y = os.path.join(dirOutputs, "y" + str(iteration) + '.csv')
    file_x_bin = os.path.join(dirOutputs, "x_bin" + str(iteration) + '.csv')
    file_y_bin = os.path.join(dirOutputs, "y_bin" + str(iteration) + '.csv')
    file_sigma_bin = os.path.join(dirOutputs, "sigma_bin" + str(iteration) + '.csv')
    np.savetxt(file_x, x, delimiter=',')
    np.savetxt(file_y, y, delimiter=',')
    np.savetxt(file_x_bin, x_bin, delimiter=',')
    np.savetxt(file_y_bin, y_bin, delimiter=',')
    np.savetxt(file_sigma_bin, sigma_bin, delimiter=',')

    del x, y
    
    # curve fitting
    p0 = [(np.max(y_bin) - np.min(y_bin))/2, 0.7, 0.4]
    popt, pcov = curve_fit(CosineFitting, x_bin, y_bin, p0=p0)
    ShiftX1 = popt[0]*np.sin(popt[1])
    ShiftY1 = popt[0]*np.cos(popt[1])
    ShiftX.append(ShiftX1)
    ShiftY.append(ShiftY1)

    del x_bin, y_bin

    print "Parameter a of iteration {0}: {1:.1f}".format(iteration, popt[0])
    print "Parameter b of iteration {0}: {1:.1f}".format(iteration, popt[1])
    print "Parameter c of iteration {0}: {1:.1f}".format(iteration, popt[2])
    print "Shift vector X of iteration {0}: {1:.1f}".format(iteration, ShiftX1)
    print "Shift vector Y of iteration {0}: {1:.1f}".format(iteration, ShiftY1)
    
    # Solve for parameters (a, b and c) iteratively until the improvement of std less than 2%
    if iteration>1:
        logic1 = abs(result_std_dh[iteration]) < 0.1
        logic2 = abs(result_std_dh[iteration]) <= abs(result_std_dh[iteration-1])
        logic3 = abs((result_std_dh[iteration-1] - result_std_dh[iteration])/(result_std_dh[iteration-1]+1e-4)) < 0.02       
        logic4 = logic2 and logic3
        logic5 = iteration>=7

        if logic1 or logic4 or logic5:
            #DEM_slave_final = "DEM_shift" + str(iteration-1) # just string
            sum_ShiftX = np.sum(ShiftX)
            sum_ShiftY = np.sum(ShiftY)

            # final shift vector [unit: m]
            shiftVec = [sum_ShiftX, sum_ShiftY, round(result_mean_dh[iteration-1],1)]
            np.savetxt(file_shiftVec, shiftVec, delimiter=',')
               
            print "********************Final Result********************"
            print "Note: results are saved in {0}".format(dirOutputs)
            print "The final shift X: {0:.1f}".format(sum_ShiftX)
            print "The final shift Y: {0:.1f}".format(sum_ShiftY)
            break        
    
    # Shift the slave DEM
    DEM_slave_after = "DEM_shift" + str(iteration)
    if arcpy.Exists(DEM_slave_after):
        arcpy.Delete_management(DEM_slave_after)
    
    arcpy.Shift_management(DEM_slave_before, DEM_slave_after, str(ShiftX1), str(ShiftY1))
    DEM_slave_before = DEM_slave_after

# correct the shifted slave DEM
DEM_slave_correct = Raster("DEM_shift" + str(iteration-1)) + round(result_mean_dh[iteration-1],1)
# DEM_slave_correct = Raster(DEM_slave_final) + round(result_mean_dh[iteration-1],1)

if arcpy.Exists(corrected_DEM):
    arcpy.Delete_management(corrected_DEM)

DEM_slave_correct.save(corrected_DEM)
print "The final shifted DEM: " + corrected_DEM

for i in range(1,iteration):
    arcpy.Delete_management("DEM_shift"+str(iteration))

endTime = time.clock()
print "Running time: {0} s".format(int(endTime-startTime))