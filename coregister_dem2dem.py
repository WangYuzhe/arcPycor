# ---------------------------------------------------------------------------
# coregister_dem2dem.py
# Created on: 2014-11-8
# Description: Coregister DEM_B to DEM_A using Nuth's method
# Reference: Nuth, Kaab, 2011, the Cryosphere
# Author: Yuzhe Wang
# E-mail: yuzhe.wang@foxmail.com
# Affiliation: State Key Laboratory of Cryospheric Sciences, Lanzhou
# History:
# 2014-11-8, 1st version;
# 2014-11-10, add the curve fitting function and while iteration;
# 2014-11-11, edit the threshold of dh and judgement for iteration;
# 2014-11-11, add statistics result of all dh of off glacier;
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
env.workspace = "E:/ICESat/CoregisterDEM/"

# File saving the results of x and y as txt
ResultFile = "E:/ICESat/CoregisterDEM/Result/"

# Local variables:
DEM_A = "dem_a"
    # Master DEM
DEM_B = "dem_b"
    # Need to be coregistered DEM
OffGlacier = "OffGlacier.shp"
    # Stable terrain assumed no elevation change
# ---------------------------------------------------------------------------
coregister = 1
iteration = 0
DEM_B_before = DEM_B
DEM_B_after = DEM_B
result_mean=[0]
result_std=[0]
ShiftX=[0]
ShiftY=[0]
while coregister==1:
	iteration = iteration+1
	print('--------------------------------------------------------------------')
	print("iteration "+" : "),
	print(iteration)
	# Process: Raster Calculator
	dh = Raster(DEM_A)-Raster(DEM_B_after)
	# Process: Slope
	slp = Slope(DEM_B, "DEGREE", "1")
	# Process: Aspect
	asp = Aspect(DEM_B)
	# Process: Extract by Mask. dh1 is masked by OffGlacier.shp
	dh_mask = ExtractByMask(dh, OffGlacier)
	# Process: Extract by Mask. slp1 and asp1 are masked by dh1_mask
	slp_mask = ExtractByMask(slp, dh_mask)
	asp_mask = ExtractByMask(asp, dh_mask)
	# ---------------------------------------------------------------------------
	dh_mask = arcpy.RasterToNumPyArray(dh_mask)
	slp_mask = arcpy.RasterToNumPyArray(slp_mask)
	asp_mask = arcpy.RasterToNumPyArray(asp_mask)
	
	# Select |dh|<70m and slope > 5
	pos1 = np.where((dh_mask>-70)&(dh_mask<70)&(slp_mask>5))
	dh_mask1 = dh_mask[pos1[0], pos1[1]]
	slp_mask1 = slp_mask[pos1[0], pos1[1]]
	asp_mask1 = asp_mask[pos1[0], pos1[1]]	

	# print("Mean dh of iteration "+str(iteration)+" : "),
	# print(np.mean(dh_mask1))
	# print("Standard deviation dh of iteration "+str(iteration)+" : "),
	# print(np.std(dh_mask1))
	
	print("Mean dh of iteration "+str(iteration)+" : "),
	print(np.mean(dh_mask))
	print("Standard deviation dh of iteration "+str(iteration)+" : "),
	print(np.std(dh_mask))
	
	result_mean.append(np.mean(dh_mask1))
	result_std.append(np.std(dh_mask1))
	# ---------------------------------------------------------------------------
	# Preparing the x and y value for curve fitting
	# x is the aspect of DEM_B, and y the dh/slope
	x = asp_mask1
	y = dh_mask1/sp.tan(np.pi*slp_mask1/180)
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
			DEM_B_final = "DEM_B" + str(iteration-1)
			sum_ShiftX = np.sum(ShiftX)
			sum_ShiftY = np.sum(ShiftY)
			print("********************Final Result********************")
			print("The results of x and y value are saved in the file: " + ResultFile)
			print("The final shifted DEM is " + DEM_B_final)
			print("The final shift X is : "),
			print(sum_ShiftX)
			print("The final shift Y is : "),
			print(sum_ShiftY)
			break
		else:
			coregister=1

	# Output shifted DEM_B
	DEM_B_after = "DEM_B" + str(iteration)
	arcpy.Shift_management (DEM_B_before, DEM_B_after, str(ShiftX1), str(ShiftY1))
	DEM_B_before = DEM_B_after
