# AneuSI - v1.0

## Description

AneuSI (Aneurysm Surface Isolation) is a command line tool for isolating the aneurysm dome and its neighboring vessels from the artery tree in 3D surfaces models of intracranial aneurysms. This tool receives the .vtk files corresponding to the surface model of the aneurysm, its centerlines and its neck polygon as inputs, and gives the isolated model surface as a .vtk PolyData file. 

The algorithm selects reference points for each branch of the artery, and calculates the isolation distance from those points as the product of the inner diameter of the blood vessel in those reference points and a user-defined factor named "clipFactor", that allows the user to control the length of the remaining vasculature connected to the aneurysm dome. Reference points are selected depending of the location of the vessel: if it branches before reaching the aneurysm neck, the point is selected immediately after the bifurcation, and if the vessel reaches the neck, the point is selected where the vessel touches the ostium.

Examples of the output can be found in the "examples" directory in the repository. 

## Requirements

- OS: Linux
- VTK 8.2.0 or higher
- CMake 3.0 or higher

## Compilation and Installation

This tool is meant to be compiled and installed by using CMake. We encourage the use of CMake gui. Before installing AneuSI, remember to install VTK and to add local directory to your path.

Procedure: 
	
	1) Inside the building directory, use cmake gui to open the configuration. 

		ccmake ~/PATH_TO_FOLDER/AneuSI/src

	2) Set the variables:
		
		CMAKE_BUILD_TYPE                 build_type                                      
		CMAKE_INSTALL_PREFIX             desired/prefix

	and press c. After configuration, press g to generate. The default options are already defined as:

		CMAKE_BUILD_TYPE                 Release
		CMAKE_INSTALL_PREFIX             /home/username/local


	3) Compile the application, by executing "make" in the terminal.

	4) Install by executing "make install".
	
## Usage

This tool requires the use of a configuration file (.cfg), which can be found in the "Examples" folder in the repository. This configuration file contains the path of the .vtk input files, the name of the output file and the parameters used in the process, such as the aneurysm type (1 for laterals, 2 for terminals) and the clip factor, the multiple of the inner diameter of the vessel that will be used as distance for the clip. 

Once defined all paths and parameters in the .cfg file, the tool can be executed from the terminal as:

	AneuSI AneuConfigFile.cfg 

If the verbose parameter was set as "true", the tool will display some information about the process. AneuSI was developed as a purely command line application, so it does not have a GUI. Nevertheless, setting the rendering parameter "true" will open a windows after the process showing the final result, in order to perform a "quick" visual check of the isolated model. 

## Troubleshooting

Software has been tested in Linux Debian 10-12 with VTK 8.2-9.0. Issues with other distributions or problems with the algorithm will be shown here. 

Please, report issues or bugs to _______________.




