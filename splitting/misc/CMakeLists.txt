message(STATUS "PoissonReconstruction: Building as a Remote VTK module")

set(Module_SRCS
    Factor.cpp
    MarchingCubes.cpp
    MultiGridOctreeData.cpp
    Geometry.cpp
    MultiGridOctest.cpp
  )
set_source_files_properties(
    Factor.cpp
    MarchingCubes.cpp
    MultiGridOctreeData.cpp
    Geometry.cpp
    MultiGridOctest.cpp
    PROPERTIES WRAP_EXCLUDE 1
)

set(classes
vtkPoissonReconstruction)

if (VTK_VERSION VERSION_LESS "8.90.0")
  vtk_module_library(PoissonReconstruction ${Module_SRCS})
else()
  vtk_module_add_module(VTK::PoissonReconstruction
    SOURCES ${Module_SRCS}
	CLASSES ${classes}
    )
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options (PoissonReconstruction PUBLIC -Wno-class-memaccess -Wno-unused-variable -Wno-unused-but-set-variable)
  endif()
endif()

