# arcPycor
arcPycor can automatically align the geolocation of DEMs and elevation points (e.g. ICESat, LiDAR) by specifying stable terrain areas. We provide two modes of arcPycor, i.e. ArcToolbox mode and standalone mode.

arcPycor was developed by Dr. Yuzhe Wang (yuzhe.wang AT foxmail.com)  now at Shandong Normal University.

If arcPycor helped you, please 

# Implementation
**Dependencies**:

* python2.7
* arcpy
* scipy
* numpy

arcPycor is executed as a standard ArcGIS geoprocessing toolbox in ArcMap. To add the toolbox, open an ArcMap document, click the toolbox button, right click in the ArcToolbox window and select ‘Add Toolbox’. Navigate to the folder downloaded from this github site and select the toolbox file (‘Coregister_DEM.tbx’). Coregister_DEM will now appear in the list, containing two tools: **Coregister DEM to DEM** and **Coregister DEM to Elevation points**.

To co-register elevation datasets, open the **Coregister_DEM** tool and specify the parameters as described below.

**Workspace**: Set the workspace

**Master DEM**: Select the master DEM (any raster format)

**Slave DEM**: Select the slave DEM (any raster format)

**Stable Terrain**: Select the stable terrain areas (shapefile)

**Master Elevation points**: Select the elevation points (shapefile)

**Elevation field name**: Set the field name representing the elevation of 'Master Elevation points' (e.g. z)

![Implementation details](https://github.com/WangYuzhe/arcPycor/blob/master/fig_implementation.png)

# How to cite acrPycor?
Please reference the Wang Yuzhe and Ye Qinghua (2021) paper.
@Article{WangYuzhe2021_jms,
  author   = {Wang, Yuzhe and Ye, Qinghua},
  journal  = {Journal of Mountain Science},
  title    = {arcPycor: an open-source automated GIS tool to co-register elevation datasets},
  year     = {2021},
  number   = {4},
  pages    = {923--931},
  volume   = {18},
  doi      = {10.1007/s11629-020-6305-y}
}
