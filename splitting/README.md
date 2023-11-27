# vortices_statistics
Get the statistics of vortices in turbulent flows based on the geometric and physical analysis.

Installation instructions Windows 10.

Prerequisites: 
1. CMake >= 3.20. Get the installer from https://cmake.org/download/.
Make sure to select "Add CMake to the System Path for the Current User" 
when prompt.

2. Git. Get it from https://git-scm.com/downloads. Install with the default options. 

3. Visual Studio Community 2022. Get it from https://visualstudio.microsoft.com/vs/community/.
While installing, select the package "Desktop development with C++".

4. Download and install VTK.
	
	Open Command Prompt and cd to where you want to build and install VTK. 
	Execute the following commands.

	```
	git clone https://github.com/Kitware/VTK.git. 
	```

	This will download the VTK source code and make a folder VTK. Let Path/To/VTK is the VTK source directory.

	```
	cd VTK & git checkout v9.1.0
	```

	```
	cd Remote & git clone https://github.com/lorensen/PoissonReconstruction.git
	```

	Go to Path/To/VTK/Remote/PoissonReconstruction and replace the CMakeLists.txt file with the CMakeLists.txt file provided in the misc folder of this repository. Also delete the file PoissonReconstruction.remote.cmake in Path/To/VTK/Remote/ folder.

	```
	cd .. & mkdir build & cd build
	```

	```
	cmake -DCMAKE_BUILD_TYPE=Release -DVTK_MODULE_ENABLE_VTK_FiltersParallelDIY2=YES -DVTK_MODULE_ENABLE_VTK_PoissonReconstruction=YES -DVTK_MODULE_ENABLE_VTK_TestingCore=DONT_WANT -DVTK_MODULE_ENABLE_VTK_TestingRendering=DONT_WANT ..
	```

	You should see "Configure Done" and "Generating Done" after the configuration is finished.
	
	Run Visual Studio as administrator. Select "Open Project or Solution" and 
	open the VTK solution file in Path/To/VTK/build folder.
	
	Change the build type from "Debug" to "Release" on top of the visual studio.
	
	Right click ALL_BUILD in the solution explorer and press Build. Wait for the build process to finish. 
	
	Make sure you get "========== Build: 259 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========" at the end of the build process.
	
	Now in the Solution Explorer scroll down to INSTALL, right click and select Build. This will build the required headers and libraries in C:\Program Files (x86)\VTK or C:\Program Files\VTK.

	Add C:\Program Files (x86)\VTK\bin to the Path environment variable and restart the system.
			
5. Download and Configure CGAL. In command prompt, cd to the directory where you want to install Vcpkg.

	```
	git clone https://github.com/microsoft/vcpkg
	cd vcpkg & .\bootstrap-vcpkg.bat
	vcpkg.exe install yasm-tool:x86-windows
	vcpkg.exe install cgal:x64-windows eigen3:x64-windows
	```
	
7. Build this project.

	Clone this repository using git or download the code directly.
	In Command Prompt, cd to Path/To/Code.

	```
	mkdir build & cd build & cmake-gui ..
	```
	
	Press Configure. 
	
	Select the option "Specify toolchain file for cross-compiling" and click Next. 
	
	Select the file Path/To/vcpkg/scripts/buildsystems/vcpkg.cmake and click Finish. 
	
	Press Configure and let the cmake configure files.
	
	After the configuration is done and you see the message "Configuring done", press Generate. 
	
	After "Generating done", select Open Project. Right click "ALL_BUILD" and build.
	
8. Instructions to run the code. In the Command Prompt, cd to Path/To/Code/build/Release and run the following commands in order.

	```
	Convert_Dataset.exe "Path\To\RawFile" "Path\To\VTKfile.vtk" xDim yDim zDim
	Create_Full_Dataset.exe "Path\To\VTKfile.vtk" "Path\To\FullDataVTKfile.vtk"
	Compute_Histogram.exe "Path\To\FullDataVTKfile.vtk" "Path\To\Dataset_steps.csv"
	Region_Growing.exe "Path\To\FullDataVTKfile.vtk" "Path\To\Dataset_steps.csv" "Path\To\Dataset_Regions.vtk" "cutoff_threshold"
	Region_Splitting.exe "Path\To\Dataset_Regions.vtk" "Path\To\Dataset_steps.csv" "Path\To\Dataset_SplitRegions.vtk" "Path\To\Dataset_SplitRegions.json" "VSF"
	Region_Simplification.exe "Path\To\Dataset_SplitRegions.vtk" "Path\To\Dataset_SimpleRegions.vtk"
	Skeleton_Extraction.exe "Path\To\Dataset_SimpleRegions.vtk" "Path\To\Dataset_Skeleton.vtk" "Path\To\Dataset_Surface.vtk"
	Vortex_Statistics.exe "Path\To\Dataset_Skeleton.vtk" "Path\To\Dataset_Surface.vtk" 7 "Path\To\Dataset_VortSkeleton.vtk" "Path\To\Dataset_VortSurface.vtk"
	Profile_Construction.exe "Path\To\Dataset_SplitRegions.json" "Path\To\Dataset_SplitRegions.vtk" "Path\To\Dataset_SkelBlock.vtk" "Path\To\Dataset_SurfBlock.vtk" "Path\To\Dataset_NewSplitRegions.json"
	```
	
	
