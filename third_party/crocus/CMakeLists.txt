cmake_minimum_required(VERSION 3.20)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    project(CHMCrocus LANGUAGES CXX Fortran)
    include(CTest)
endif()

add_library(chm_crocus STATIC
    sfc_options.f90
    crocus_bridge.f90
    upstream/ini_csts.F90
    upstream/ini_surf_csts.F90
    upstream/modd_blowsnw_surf.F90
    upstream/modd_const_atm.F90
    upstream/modd_const_tartes.F90
    upstream/modd_csts.F90
    upstream/modd_meb_par.F90
    upstream/modd_prep_snow.F90
    upstream/modd_reprod_oper.F90
    upstream/modd_snow_metamo.F90
    upstream/modd_snow_par.F90
    upstream/modd_surf_atm.F90
    upstream/modd_surf_par.F90
    upstream/modd_type_date_surf.F90
    upstream/mode_crodebug.F90
    upstream/mode_snow3l.F90
    upstream/mode_snowcro.F90
    upstream/mode_tartes.F90
    upstream/mode_thermos.F90
    upstream/modi_abor1_sfx.F90
    upstream/modi_ini_csts.F90
    upstream/modi_ini_surf_csts.F90
    upstream/modi_snowcro.F90
    upstream/modi_surface_aero_cond.F90
    upstream/modi_surface_cd.F90
    upstream/modi_surface_ri.F90
    upstream/modi_wind_threshold.F90
    upstream/parkind1.F90
    upstream/snowcro.F90
    upstream/surface_aero_cond.F90
    upstream/surface_cd.F90
    upstream/surface_ri.F90
    upstream/tridiag_ground_snowcro.F90
    upstream/wind_threshold.F90
    upstream/yomhook.F90
)
target_compile_definitions(chm_crocus PRIVATE CROCUS_EXT ECCCGEM)
target_include_directories(chm_crocus PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
set_target_properties(chm_crocus PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
# The CHM face loop can enter the native column concurrently from OpenMP threads.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(chm_crocus PRIVATE -frecursive -ffree-line-length-none)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(chm_crocus PRIVATE -recursive)
endif()

if(BUILD_TESTS OR (PROJECT_NAME STREQUAL "CHMCrocus" AND BUILD_TESTING))
    find_package(Threads REQUIRED)
    add_executable(test_crocus_column tests/test_column.cpp)
    target_link_libraries(test_crocus_column PRIVATE chm_crocus Threads::Threads)
    add_test(NAME crocus_column COMMAND test_crocus_column "${CMAKE_CURRENT_SOURCE_DIR}/tests")
endif()
