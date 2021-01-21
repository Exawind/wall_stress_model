# Complex Terrain 

## Overview
This directory contains the files and data for running the flow over a two-dimensional ridge and
and flow over a three-dimensional hill

## Basic settings

| Parameter       | Value                        |
|-----------------|------------------------------|
| Reynolds number | 10,500 based on ridge height |
| Hill height     | 0.04 m                       |
| Fluid density   | 1 kg/m^3                     |
| Fluid velocity  | 0.3 m/s                      |

**Mesh generation**  
Meshes for this case are generated using wind_utils abl_mesh and pystk based python script
