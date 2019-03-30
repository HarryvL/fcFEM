# fcFEM
Finite Element macro for FreeCAD

<img src="https://user-images.githubusercontent.com/35259498/55271886-bec9d500-52b4-11e9-936c-122f9072406d.png" width="425"/> <img src="https://user-images.githubusercontent.com/35259498/55271882-bd98a800-52b4-11e9-912f-b319b0c0694f.png" width="425"/>


### Description
fcFEM is a general finite element framework for 3D deformation analysis, currently supporting elasto-plastic collapse analysis and interface elements

### Installation and running
Install files fcFEM.FCMacro and femTools.py in the .FreeCAD/Macro folder on your machine  
Open with FreeCAD macro edititor

### Documentation
Please refer to source code for in-line comments and to the FreeCAD forum (https://forum.freecadweb.org/viewforum.php?f=18)

### TODO
Optimal use of LAPACK and BLAS libraries to increase speed of the incremental-iterative colver  
Addition of beam and shell elements  
Linear buckling and initial imperfections for non-liner buckling  
Loading stages  
Geometric non-linearity  
Advanced material modelling  

### Licence information
                                                                         
Copyright (c) 2019 - Harry van Langen <hvlanalysis@icloud.com>        
                                                                         
This program is free software; you can redistribute it and/or modify  
it under the terms of the GNU Lesser General Public License (LGPL)    
as published by the Free Software Foundation; either version 2 of     
the License, or (at your option) any later version.                   
for detail see the LICENCE text file.                                 
                                                                         
This program is distributed in the hope that it will be useful,       
but WITHOUT ANY WARRANTY; without even the implied warranty of        
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         
GNU Library General Public License for more details.                  
                                                                         
You should have received a copy of the GNU Library General Public     
License along with this program; if not, write to the Free Software   
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
USA                                                                   
