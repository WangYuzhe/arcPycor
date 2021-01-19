'''
pycor_dem2point.py

Description: Coregister slave DEM to master DEM using Nuth-Kaeaeb algorithm.
Reference: Nuth, C. and Kaeaeb, A. (2011): http://doi.org/10.5194/tc-5-271-2011

Inputs:
(1) master elevation points
(2) slave DEM
(3) polygon shapefile of stable terrain

Outputs:
(1) x_bin (.csv; aspect)
(2) y_bin (.csv; dh/tan(slope))
(3) sigma_bin (.csv; standard deviation of y)
(4) shiftVec txt (shift vector)

History:
* 2014-11-17, Created
* 2019-5-11, delete unused variables, optimize the codes for converging

Author: Yuzhe Wang
E-mail: wangyuzhe@ucas.ac.cn
Affiliation:
1. Colledge of Resources and Environment, University of Chinese Academy Sciences, Beijing
2. State Key Laboratory of Cryospheric Sciences, Chinese Academy of Sciences, Lanzhou
'''

import os
import time
import numpy as np
from scipy.optimize import curve_fit
import arcpy
from arcpy import env
from arcpy.sa import *
import arcpy.da
import shutil
import matplotlib.pyplot as plt
# import matplotlib as mpl

startTime = time.clock()

if arcpy.CheckExtension("spatial")=="Available":
    arcpy.CheckOutExtension("spatial")
else:
    raise LicenseError

# Set environment workspace in current directory
# path_script = os.path.dirname(os.path.abspath(__file__))
# env.workspace = os.path.join(path_script, 'benchmark_data')
env.workspace = r"E:\Mix\WYZ_AcademicWriting\paper_pycor\demo_data\benchmark_data"
              
# Folder for outputs
dirOutputs = os.path.join(env.workspace, 'outputs_dem2point')
if os.path.exists(dirOutputs):
    shutil.rmtree(dirOutputs)
os.makedirs(dirOutputs)

# point shapefile
shp_points = "points_center.shp"

# the elevation field in point shapefile
elev_point_fieldName = "z"

# slave DEM
DEM = "DEM_slave1.tif"

# stable terrain
OffGlacier = "OffGlacier.shp"

# Determine the points located in the polygons
points_in_polygon = "points_in_polygon.shp"
if arcpy.Exists(points_in_polygon):
    arcpy.Delete_management(points_in_polygon)
arcpy.SpatialJoin_analysis(shp_points, OffGlacier, points_in_polygon)
arcpy.DeleteField_management(points_in_polygon, ["Join_Count", "TARGET_FID"])

# Initializations
iteration = 0
DEM_before = DEM
DEM_after = DEM
result_mean_dh = [0]
result_std_dh = [0]
ShiftX = [0]
ShiftY = [0]
file_shiftVec = os.path.join(dirOutputs, "shiftVec" + ".csv")

# Define CosineFitting function
def CosineFitting(x, a, b, c):
    return a*np.cos(b - np.pi/180*x) + c

while 1:
    iteration = iteration + 1
    arcpy.AddMessage("--------------------------------------------------------------")
    arcpy.AddMessage("Iteration {0} is running!".format(iteration))

    # Get the slope of the slave DEM [degree]
    slp = Slope(DEM_after, "DEGREE", "1")

    # Get the aspect of the slave DEM [degree]
    asp = Aspect(DEM_after)

    # Extract raster values using point shapefile
    dem_fieldName = "dem" + str(iteration)
    slp_fieldName = "slp" + str(iteration)
    asp_fieldName = "asp" + str(iteration)
    inRasterList = [[DEM_after, dem_fieldName], [slp, slp_fieldName], [asp, asp_fieldName]]
    ExtractMultiValuesToPoints(points_in_polygon, inRasterList, "BILINEAR")

    #del slp, asp
    arcpy.Delete_management(slp)
    arcpy.Delete_management(asp)

    # Read attribute table
    elev_point = [0]
    elev_dem = [0]
    slp_table = [0]
    asp_table = [0]
    dh_table = [0]
    with arcpy.da.SearchCursor(points_in_polygon, (elev_point_fieldName,dem_fieldName,slp_fieldName,asp_fieldName)) as cursor:
        for row in cursor:
            elev_point.append(row[0])
            elev_dem.append(row[1])
            slp_table.append(row[2])
            asp_table.append(row[3])

    del cursor, row

    elev_point = np.array(elev_point[1:])
    elev_dem = np.array(elev_dem[1:])
    slp_table = np.array(slp_table[1:])
    asp_table = np.array(asp_table[1:])
    dh_table = elev_point - elev_dem

    # Criteria: |dh| < 70 m and 5 < slope < 45.
    # slope>5 is cited from Purinton&Bookhagen, 2018, Earth Surface Dynamics.
    # slope<45 is cited from Berthier et al., 2019, Journal of Glaciology.
    index1 = np.where((dh_mask_arr>-70) & (dh_mask_arr<70) & (slp_mask_arr>5) & (slp_mask_arr<45))
    dh_table1 = dh_table[index1]
    slp_table1 = slp_table[index1]
    asp_table1 = asp_table[index1]

    # save "dh" as csv file
    file_dh = os.path.join(dirOutputs, "dh" + str(iteration) + '.csv')
    np.savetxt(file_dh, dh_table, delimiter=',')
    
    del dh_table, slp_table, asp_table, index1

    # statistic results of dh
    mean_dh = np.mean(dh_table1)
    result_mean_dh.append(mean_dh)
    
    std_dh = np.std(dh_table1)
    result_std_dh.append(std_dh)

    arcpy.AddMessage("Mean dh of iteration {0}: {1:.1f}".format(iteration, mean_dh))
    arcpy.AddMessage("Standard deviation of dh of iteration {0}: {1:.1f}".format(iteration, std_dh))

    # Prepare the x and y values for curve fitting
    x = asp_table1
    y = dh_table1 / np.tan(np.pi*slp_table1/180)    

    del dh_table1, slp_table1, asp_table1

    # Get the x, y values in bins
    range_asp = range(0, 370, 10)
    n = len(range_asp) - 1
    x_bin = np.zeros(n)
    y_bin = np.zeros(n)
    sigma_bin = np.zeros(n)
    for i in range(n):
        index2 = np.where( (x >= range_asp[i])&(x < range_asp[i+1]) )
        x_bin[i] = range_asp[i] + 5
        y_bin[i] = np.median(y[index2])
        sigma_bin[i] = np.std(y[index2])

    # Save the results (x, y, x_bin, y_bin, sigma_bin) in txt format
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
        logic1 = abs(result_std_dh[iteration]) < 0.1
        logic2 = abs(result_std_dh[iteration]) <= abs(result_std_dh[iteration-1])
        logic3 = abs((result_std_dh[iteration-1] - result_std_dh[iteration])/(result_std_dh[iteration-1]+1e-4)) < 0.02       
        logic4 = logic2 and logic3
        logic5 = iteration>=7

        if logic1 or logic4 or logic5:
            DEM_final = "DEMsh" + str(iteration-1)
            sum_ShiftX = np.sum(ShiftX)
            sum_ShiftY = np.sum(ShiftY)

            # final shift vector [unit: m]
            shiftVec = [sum_ShiftX, sum_ShiftY]
            np.savetxt(file_shiftVec, shiftVec, delimiter=',')

            arcpy.AddMessage("********************Final Result********************")
            arcpy.AddMessage("Note: results are saved in {0}".format(dirOutputs))
            arcpy.AddMessage("The final shift X: {0:.1f}".format(sum_ShiftX))
            arcpy.AddMessage("The final shift Y: {0:.1f}".format(sum_ShiftY))            
            break

    # Output shifted DEM
    DEM_after = "DEMsh" + str(iteration)
    arcpy.Shift_management(DEM_before, DEM_after, str(ShiftX1), str(ShiftY1))
    DEM_before = DEM_after

# correct the shifted slave DEM
DEM_slave_correct = Raster(DEM_final) + round(result_mean_dh[iteration-1],1)
DEM_slave_correct.save("DEM_slave_correct.tif")
arcpy.AddMessage("The final shifted DEM: " + "DEM_slave_correct.tif")

if iteration<2:
    for i in range(iteration-1):
        iter_tag = i+1
        indata_del = "DEMsh"+str(iter_tag)
        arcpy.Delete_management(indata_del)
    
endTime = time.clock()
arcpy.AddMessage("Running time: {0}".format(int(endTime-startTime)))