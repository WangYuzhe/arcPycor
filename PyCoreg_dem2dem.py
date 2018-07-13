'''
PyCoR_dem2dem.py

Created on: 2014-11-8
Latest modify: 2017-6-5
Description: Coregister slave DEM (DEM_B) to master DEM (DEM_A) using Nuth's method.
Reference: Nuth, Kaab, 2011, the Cryosphere

Author: Yuzhe Wang
E-mail: yuzhe.wang@foxmail.com
Affiliation: State Key Laboratory of Cryospheric Sciences, Chinese Academy of Sciences, Lanzhou
'''

import arcpy
import numpy as np
import scipy as sp

#import matplotlib.pyplot as plt
#import matplotlib as mpl

from scipy.optimize import curve_fit
from arcpy import env
from arcpy.sa import *

# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Environment setting
env.workspace = "E:/cloud/AcademicWriting_WYZ/paper_PyCoR/data/run_sample_1arc/"

# File saving the results of x and y as txt
FileResults = "E:/cloud/AcademicWriting_WYZ/paper_PyCoR/results/run_sample_1arc/"

# Master DEM (template)
DEM_A = "dem_a.tif"

# Slave DEM
DEM_B = "dem_b.tif"

# Stable terrain assumed no elevation changes
OffGlacier = "OffGlacier.shp"

# Initializations
coregister = 1
iteration = 0
DEM_B_before = DEM_B
DEM_B_after = DEM_B
result_mean = [0]
result_std = [0]
ShiftX = [0]
ShiftY = [0]

while coregister == 1:
    iteration = iteration + 1
    print('--------------------------------------------------------------------')
    print('Iteration %d is ruuning!' %iteration)
    
    # DEM differencing [m]
    dh = Raster(DEM_A) - Raster(DEM_B_after)
    
    # Get the slope of the slave DEM [degree]
    slp = Slope(DEM_B, "DEGREE", "1")
    
    # Get the aspect of the slave DEM [degree]
    asp = Aspect(DEM_B)
    
    # Mask 'dh' using 'OffGlacier'
    dh_mask = ExtractByMask(dh, OffGlacier)
    
    # Mask 'slp' and 'asp' using 'dh_mask'
    # In order to keep the same georeference as 'dh_mask', here the 'slp' and 'asp' are masked by 'dh_mask'.
    slp_mask = ExtractByMask(slp, dh_mask)
    asp_mask = ExtractByMask(asp, dh_mask)
    
    # Raster to Array
    dh_mask = arcpy.RasterToNumPyArray(dh_mask)
    slp_mask = arcpy.RasterToNumPyArray(slp_mask)
    asp_mask = arcpy.RasterToNumPyArray(asp_mask)
    
    # Criteria: |dh| < 70 m and slope > 5 degrees
    index1 = np.where((dh_mask > -70) & (dh_mask < 70) & (slp_mask > 5))
    dh_mask1 = dh_mask[index1[0], index1[1]]
    slp_mask1 = slp_mask[index1[0], index1[1]]
    asp_mask1 = asp_mask[index1[0], index1[1]]
    
    # print the statistic results
    print('Mean dh of iteration %d is %0.1f' %(iteration, np.mean(dh_mask1)))
    print('Standard deviation of dh of iteration %d is %0.1f' %(iteration, np.std(dh_mask1)))
    
    # 待完善, 全部使用array
    result_mean.append(np.mean(dh_mask1))
    result_std.append(np.std(dh_mask1))
    
    # Prepare the x and y values for curve fitting
    x = asp_mask1
    y = dh_mask1/np.tan(np.pi*slp_mask1/180)
    
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

    # Save the results (x, y, x_bin, y_bin, sigma_bin) in texts
    file_x = FileResults + "x" + str(iteration) + '.txt'
    file_y = FileResults + "y" + str(iteration) + '.txt'
    file_x_bin = FileResults + "x_bin" + str(iteration) + '.txt'
    file_y_bin = FileResults + "y_bin" + str(iteration) + '.txt'
    file_sigma_bin = FileResults + "sigma_bin" + str(iteration) + '.txt'
    
    np.savetxt(file_x, x)
    np.savetxt(file_y, y)
    np.savetxt(file_x_bin, x_bin)
    np.savetxt(file_y_bin, y_bin)
    np.savetxt(file_sigma_bin, sigma_bin)

    # Define CosineFitting function
    def CosineFitting(x, a, b, c):
        return a*np.cos(b - np.pi/180*x) + c

    p0 = [(np.max(y_bin) - np.min(y_bin))/2, 0.7, 0.4]
    popt, pcov = curve_fit(CosineFitting, x_bin, y_bin, p0 = p0)
    ShiftX1 = popt[0]*np.sin(popt[1])
    ShiftY1 = popt[0]*np.cos(popt[1])
    ShiftX.append(ShiftX1)
    ShiftY.append(ShiftY1)

    print('Parameter a of iteration %d is %0.1f' %(iteration, popt[0]))
    print('Parameter b of iteration %d is %0.1f' %(iteration, popt[1]))
    print('Parameter c of iteration %d is %0.1f' %(iteration, popt[2]))
    print('Shift vector X of iteration %d is %0.1f' %(iteration, ShiftX1))
    print('Shift vector Y of iteration %d is %0.1f' %(iteration, ShiftY1))
    
    # Solve for parameters (a, b and c) iteratively until the improvement of std less than 2%
    if iteration == 1:
        coregister = 1
    else:
        logic1 = result_std[iteration] <= result_std[iteration-1]
        logic2 = (result_std[iteration-1] - result_std[iteration])/(result_std[iteration-1]+0.0001) < 0.02
        
        if logic1 and logic2:
            coregister = 0
            DEM_B_final = "DEM_B" + str(iteration-1)
            sum_ShiftX = np.sum(ShiftX)
            sum_ShiftY = np.sum(ShiftY)            
            
            print('********************Final Result********************')
            print('Note: all results (i.e., x, y, x_bin, y_bin, sigma_bin) are saved in the file %s' %FileResults)
            print('The final shifted DEM is %s' %DEM_B_final)
            print('The final shift X is %0.1f' %sum_ShiftX)
            print('The final shift Y is %0.1f' %sum_ShiftY)
            break
        else:
            coregister = 1
    
    # Shift the slave DEM
    DEM_B_after = "DEM_B" + str(iteration)    
    arcpy.Shift_management(DEM_B_before, DEM_B_after, str(ShiftX1), str(ShiftY1))
    DEM_B_before = DEM_B_after