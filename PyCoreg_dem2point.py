# ---------------------------------------------------------------------------
# coregister_dem2point_v3.1.py
# Created on: 2014-11-17
# Description: Coregister DEM to point shapefile (ICESat, LiDAR, ...) using Nuth's method
# Reference: Nuth, Kaab, 2011, the Cryosphere
# Author: Yuzhe Wang
# E-mail: wgyuzhe@gmail.com
# Affiliation: State Key Laboratory of Cryospheric Sciences, Lanzhou
# History:
# 2014-11-17, 1st version;
# ---------------------------------------------------------------------------
# Import modules
import arcpy
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
# import matplotlib as mpl

from scipy.optimize import curve_fit
from arcpy import env
from arcpy.sa import *
# ---------------------------------------------------------------------------
# Check out any necessary licenses
arcpy.CheckOutExtension("spatial")

# Set environment settings
env.workspace = "E:/ICESat/CoregisterDEM/test/"

# File saving the results of x and y as txt
ResultFile = "E:/ICESat/CoregisterDEM/test/Result/"

# Local variables:
point = "point.shp"
    # point elevation dataset
DEM = "srtm.tif"
    # Need to be coregistered DEM
point_fieldName = "ELEV_WGS84"
	# Field name of point's elevation in point.shp 
# ---------------------------------------------------------------------------
coregister = 1
iteration = 0
DEM_before = DEM
DEM_after = DEM
result_mean=[0]
result_std=[0]
ShiftX=[0]
ShiftY=[0]
while coregister==1:
	iteration = iteration+1
	print('--------------------------------------------------------------------')
	print("iteration "+" : "),
	print(iteration)
	
	slp = Slope(DEM_after, "DEGREE", "1")
		# Process: Slope
	asp = Aspect(DEM_after)
		# Process: Aspect
	
	# Extract raster values using point shapefile
	dem_fieldName = "dem" + str(iteration-1)
	slp_fieldName = "slp" + str(iteration-1)
	asp_fieldName = "asp" + str(iteration-1)
	inRasterList = [[DEM_after, dem_fieldName], [slp, slp_fieldName], [asp, asp_fieldName]]
	ExtractMultiValuesToPoints(point, inRasterList, "BILINEAR")
	
	# Read attribute table	
	cursor = arcpy.SearchCursor(point)
	fields = [point_fieldName,dem_fieldName,slp_fieldName,asp_fieldName]
	point_table=[0]
	dem_table=[0]
	dh_table=[0]
	slp_table=[0]
	asp_table=[0]
	for row in cursor:
		point_table.append(row.getValue(fields[0]))
		dem_table.append(row.getValue(fields[1]))
		dh_table.append(row.getValue(fields[0])-row.getValue(fields[1]))
		slp_table.append(row.getValue(fields[2]))
		asp_table.append(row.getValue(fields[3]))
	
	point_table=np.array(point_table[1:])
	dem_table=np.array(dem_table[1:])
	dh_table=np.array(dh_table[1:])
	slp_table=np.array(slp_table[1:])
	asp_table=np.array(asp_table[1:])

	# Select |dh|<70m and slope > 5
	pos1 = np.where((dh_table>-70)&(dh_table<70)&(slp_table>5))
	dh_table1 = dh_table[pos1]
	slp_table1 = slp_table[pos1]
	asp_table1 = asp_table[pos1]

	# print("Mean dh of iteration "+str(iteration)+" : "),
	# print(np.mean(dh_table1))
	# print("Standard deviation dh of iteration "+str(iteration)+" : "),
	# print(np.std(dh_table1))
	
	print("Mean dh of iteration "+str(iteration)+" : "),
	print(np.mean(dh_table))
	print("Standard deviation dh of iteration "+str(iteration)+" : "),
	print(np.std(dh_table))
	
	result_mean.append(np.mean(dh_table1))
	result_std.append(np.std(dh_table1))
	# ---------------------------------------------------------------------------
	# Preparing the x and y value for curve fitting
	# x is the aspect of DEM, and y the dh/slope
	x = asp_table1
	y = dh_table1/sp.tan(np.pi*slp_table1/180)
		# degree to radian
	interval = range(0,370,10)
	LenInterval = len(interval)-1
	bin_x = np.zeros(LenInterval)
	bin_y = np.zeros(LenInterval)
	bin_sigma = np.zeros(LenInterval)
	for i in range(LenInterval):
		j = np.where( (x>=interval[i])&(x<interval[i+1]) )
		bin_x[i] = interval[i]+5
		bin_y[i] = np.median(y[j])
		bin_sigma[i] = np.std(y[j])
	# ---------------------------------------------------------------------------
	# Save the result x, y, bin_x, bin_y, bin_sigma as text file
	file_x = ResultFile + "x"+ str(iteration)+'.txt'
	file_y = ResultFile + "y"+ str(iteration)+'.txt'
	file_bin_x = ResultFile + "bin_x"+ str(iteration)+'.txt'
	file_bin_y = ResultFile + "bin_y"+ str(iteration)+'.txt'
	file_bin_sigma = ResultFile + "bin_sigma"+ str(iteration)+'.txt'
	
	np.savetxt(file_x, x)
	np.savetxt(file_y, y)
	np.savetxt(file_bin_x, bin_x)
	np.savetxt(file_bin_y, bin_y)
	np.savetxt(file_bin_sigma, bin_sigma)
	# ---------------------------------------------------------------------------
	# Curve fitting
	def CosineFitting(x,a,b,c):
		return a*np.cos(b - np.pi/180*x) + c
	
	p0 = [(np.max(bin_y)-np.min(bin_y))/2, 0.7, 0.4]
	popt, pcov = curve_fit(CosineFitting, bin_x, bin_y, p0 = p0)
	ShiftX1 = popt[0]*np.sin(popt[1])
	ShiftY1 = popt[0]*np.cos(popt[1])
	ShiftX.append(ShiftX1)
	ShiftY.append(ShiftY1)
		
	print('Parameter a is : '),
	print(popt[0])
	print('Parameter b is : '),
	print(popt[1])
	print('Parameter c is : '),
	print(popt[2])
	print("Shift vector X of iteration "+str(iteration) + " : "),
	print(ShiftX1)
	print("Shift vector Y of iteration "+str(iteration) + " : "),
	print(ShiftY1)
	print('\n')

	# Solve for parameters (a, b and c) iteratively until the improvement of std less than 2%
	if iteration==1:
		coregister=1
	else:
		logic_judge = (result_std[iteration]<=result_std[iteration-1]) and ((result_std[iteration-1] - result_std[iteration])/(result_std[iteration-1]+0.0001)<0.02)
		if logic_judge:
			coregister=0
			DEM_final = "DEM" + str(iteration-1)
			sum_ShiftX = np.sum(ShiftX)
			sum_ShiftY = np.sum(ShiftY)
			print("********************Final Result********************")
			print("The results of x and y value are saved in the file: " + ResultFile)
			print("The final shifted DEM is " + DEM_final)
			print("The final shift X is : "),
			print(sum_ShiftX)
			print("The final shift Y is : "),
			print(sum_ShiftY)
			break
		else:
			coregister=1

	# Output shifted DEM
	DEM_after = "DEM" + str(iteration)
	arcpy.Shift_management (DEM_before, DEM_after, str(ShiftX1), str(ShiftY1))
	DEM_before = DEM_after