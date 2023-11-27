cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(vortexMain)

add_subdirectory("src/extract_skeleton")
add_subdirectory("src/create_full_dataset")
add_subdirectory("src/compute_histogram")
add_subdirectory("src/region_growing")
add_subdirectory("src/region_splitting")
add_subdirectory("src/region_simplification")
add_subdirectory("src/skeleton_extraction")
add_subdirectory("src/vortex_statistics")
add_subdirectory("src/profile_construction")