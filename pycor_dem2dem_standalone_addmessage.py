'''
pycor_dem2dem.py

Description: Coregister subordinate DEM to main DEM using Nuth-Kaeaeb algorithm.
Reference: Nuth, C. and Kaeaeb, A. (2011): http://doi.org/10.5194/tc-5-271-2011

Inputs:
(1) main DEM
(2) subordinate DEM
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

Author: Yuzhe Wang
E-mail: wangyuzhe@ucas.ac.cn
Affiliation:
1. Colledge of Resources and Environment, University of Chinese Academy Sciences, Beijing
2. State Key Laboratory of Cryospheric Sciences, Chinese Academy of Sciences, Lanzhou
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
env.workspace = r"E:\Mix\WYZ_AcademicWriting\paper_pycor\demo_data"

# Folder for outputs
dirOutputs= os.path.join(env.workspace, 'outputs')
if os.path.exists(dirOutputs):
    shutil.rmtree(dirOutputs)
os.makedirs(dirOutputs)

# main DEM
DEM_main = "DEM_main.tif"

# subordinate DEM
DEM_subordinate = "DEM_subordinate1.tif"

# stable terrain
OffGlacier = "OffGlacier.shp"

# Initializations
iteration = 0
DEM_subordinate_before = DEM_subordinate
DEM_subordinate_after = DEM_subordinate
result_mean = [0]
result_std = [0]
ShiftX = [0]
ShiftY = [0]
file_shiftVec = os.path.join(dirOutputs, "shiftVec" + ".csv")

# delete unused variables
del DEM_subordinate

# Define CosineFitting function
def CosineFitting(x, a, b, c):
    return a*np.cos(b - np.pi/180*x) + c

while 1:
    iteration = iteration + 1
    arcpy.AddMessage("--------------------------------------------------------------")
    arcpy.AddMessage("Iteration {0} is running!".format(iteration))
    
    # DEM difference [m]
    dh = Raster(DEM_main) - Raster(DEM_subordinate_after)
    
    # slope of the subordinate DEM [degree]
    slp = Slope(DEM_subordinate_after, "DEGREE", "1")
    
    # aspect of the subordinate DEM [degree]
    asp = Aspect(DEM_subordinate_after)
    
    # Mask 'dh' using statale terrain polygon
    dh_mask = ExtractByMask(dh, OffGlacier)
    
    # Mask 'slp' and 'asp' using 'dh_mask' in order to keep same georeference as 'dh_mask'
    slp_mask = ExtractByMask(slp, dh_mask)
    asp_mask = ExtractByMask(asp, dh_mask)
    
    del dh, slp, asp

    # Raster to Array
    dh_mask_arr = arcpy.RasterToNumPyArray(dh_mask, nodata_to_value=-32768)
    slp_mask_arr = arcpy.RasterToNumPyArray(slp_mask, nodata_to_value=-32768)
    asp_mask_arr = arcpy.RasterToNumPyArray(asp_mask, nodata_to_value=-32768)
    
    del dh_mask, slp_mask, asp_mask
    
    # save "dh_mask" as csv file
    file_dh = os.path.join(dirOutputs, "dh" + str(iteration) + '.csv')
    np.savetxt(file_dh, dh_mask_arr, delimiter=',')
    
    # Criteria: |dh| < 70 m and slope > 5 degrees
    index1 = np.where((dh_mask_arr>-70) & (dh_mask_arr<70) & (slp_mask_arr>5))
    dh_mask_arr1 = dh_mask_arr[index1[0], index1[1]]
    slp_mask_arr1 = slp_mask_arr[index1[0], index1[1]]
    asp_mask_arr1 = asp_mask_arr[index1[0], index1[1]]
    
    del dh_mask_arr, slp_mask_arr, asp_mask_arr, index1
    
    # statistics of "dh"
    mean_dh_mask1 = np.mean(dh_mask_arr1)
    std_dh_mask1 = np.std(dh_mask_arr1)
    result_mean.append(mean_dh_mask1)
    result_std.append(std_dh_mask1)
    arcpy.AddMessage("Mean dh of iteration {0}: {1:.1f}".format(iteration, mean_dh_mask1))
    arcpy.AddMessage("Standard deviation of dh of iteration {0}: {1:.1f}".format(iteration, std_dh_mask1))
    
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

    arcpy.AddMessage("Parameter a of iteration {0}: {1:.1f}".format(iteration, popt[0]))
    arcpy.AddMessage("Parameter b of iteration {0}: {1:.1f}".format(iteration, popt[1]))
    arcpy.AddMessage("Parameter c of iteration {0}: {1:.1f}".format(iteration, popt[2]))
    arcpy.AddMessage("Shift vector X of iteration {0}: {1:.1f}".format(iteration, ShiftX1))
    arcpy.AddMessage("Shift vector Y of iteration {0}: {1:.1f}".format(iteration, ShiftY1))
    
    # Solve for parameters (a, b and c) iteratively until the improvement of std less than 2%
    if iteration>1:
        logic1 = result_std[iteration] < 1e-1
        logic2 = result_std[iteration] <= result_std[iteration-1]
        logic3 = (result_std[iteration-1] - result_std[iteration])/(result_std[iteration-1]+1e-4) < 0.02        
        logic4 = logic2 and logic3

        if logic1 or logic4:
            DEM_subordinate_final = "DEM_shift" + str(iteration-1) # just string
            sum_ShiftX = np.sum(ShiftX)
            sum_ShiftY = np.sum(ShiftY)

            # final shift vector [unit: m]
            shiftVec = [sum_ShiftX, sum_ShiftY]
            np.savetxt(file_shiftVec, shiftVec, delimiter=',')
               
            arcpy.AddMessage("********************Final Result********************")
            arcpy.AddMessage("Note: results are saved in {0}".format(dirOutputs))
            arcpy.AddMessage("The final shifted DEM: {0}".format(DEM_subordinate_final))
            arcpy.AddMessage("The final shift X: {0:.1f}".format(sum_ShiftX))
            arcpy.AddMessage("The final shift Y: {0:.1f}".format(sum_ShiftY))
            
            break
    
    # Shift the subordinate DEM
    DEM_subordinate_after = "DEM_shift" + str(iteration)
    arcpy.Shift_management(DEM_subordinate_before, DEM_subordinate_after, str(ShiftX1), str(ShiftY1))
    DEM_subordinate_before = DEM_subordinate_after

    
for i in range(iteration-1):
    iter_tag = i+1
    indata_del = "DEM_shift"+str(iter_tag)
    arcpy.Delete_management(indata_del)

endTime = time.clock()
arcpy.AddMessage("Running time: {0}".format(int(endTime-startTime)))
