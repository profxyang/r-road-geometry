# 'r-road-geometry' Project

## Introduction
The r-road-geometry project is an open-source project to extract 3D geometry of roads and identify homogeneous alignment segments.
* horizontal curves
* horizontal tangents
* vertical curves
* vertical tangents
* combinations (i.e., horizontal curve on a vertical tangent)

## Input and Output
The 2D road geometry can be supplied by the user in shapefile or, more flexibly, extracted from public data sources such as OpenStreetMap. The output of the program is simple feature collections of homogeneous alignment segments. 

## Application of the Project
The original purpose of the program is to collect homogeneous alignment segments for traffic safety analysis. 

## Files
* The 'docs' folder contains helpful tutorials for each step of the analysis.
* The 'r' folder contains the up-to-date code of the functions.
* The 'demo' folder contains the full source code to run alignment analysis.

## Functions
The 'r' folder contains the up-to-date code of the functions. 

* get_routes(): This function performs network analysis to identify unique routes for the horizontal alignment analysis.
  - input: an sf or sfc object of 2D road geometry
  - output: an sfc of unique routes 
* get_elevation(): This function extracts elevation data from 3DEP at a specific resolution 'tres'.
  - input: an sfc of unique routes in XY format
  - output: an sfc of unique routes in XYZ format
* horizontal_analysis(): This function identifies and measures horizontal curves and tangents.
  - input: an sfc of unique routes in XY format
  - output: an sfc of horizontal curves and tangents with measurements
* vertical_analysis(): This function identifies and measures vertical curves and tangents.
  - input: an sfc of unique routes in XYZ format
  - output: an sfc of vertical curves and tangents with measurements

## Webpage
Articles related to the program can be found at the project website: https://profxyang.github.io/r-road-geometry/
