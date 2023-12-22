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

* get_routes(): This function takes input of sf or sfc.
  - sfds
  - sdfs
* get_elevation(): This function returns the 3D geometry of the road


