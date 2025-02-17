# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.5)
   message(FATAL_ERROR "CMake >= 2.6.0 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.6)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_targetsDefined)
set(_targetsNotDefined)
set(_expectedTargets)
foreach(_expectedTarget itkdouble-conversion itksys itkvcl itknetlib itkv3p_netlib itkvnl itkvnl_algo itktestlib ITKVNLInstantiation ITKCommon itkNetlibSlatec ITKStatistics ITKTransform ITKLabelMap ITKMesh itkzlib ITKMetaIO ITKSpatialObjects ITKPath ITKQuadEdgeMesh ITKIOImageBase ITKOptimizers ITKPolynomials ITKBiasCorrection ITKBioCell ITKDICOMParser ITKEXPAT ITKIOXML ITKIOSpatialObjects ITKFEM gdcmCommon gdcmDICT gdcmDSED gdcmIOD gdcmMSFF gdcmMEXD gdcmjpeg8 gdcmjpeg12 gdcmjpeg16 gdcmopenjpeg gdcmcharls gdcmsocketxx ITKznz ITKniftiio ITKgiftiio hdf5-static hdf5_cpp-static ITKIOBMP ITKIOBioRad ITKIOCSV ITKIOGDCM ITKIOIPL ITKIOGE ITKIOGIPL ITKIOHDF5 itkjpeg ITKIOJPEG itktiff ITKIOTIFF ITKIOLSM ITKIOMRC ITKIOMesh ITKIOMeta ITKIONIFTI ITKNrrdIO ITKIONRRD itkpng ITKIOPNG ITKIOSiemens ITKIOStimulate ITKTransformFactory ITKIOTransformBase ITKIOTransformHDF5 ITKIOTransformInsightLegacy ITKIOTransformMatlab ITKIOVTK ITKKLMRegionGrowing ITKOptimizersv4 itkTestDriver ITKVTK ITKVideoCore ITKVideoIO ITKWatersheds)
  list(APPEND _expectedTargets ${_expectedTarget})
  if(NOT TARGET ${_expectedTarget})
    list(APPEND _targetsNotDefined ${_expectedTarget})
  endif()
  if(TARGET ${_expectedTarget})
    list(APPEND _targetsDefined ${_expectedTarget})
  endif()
endforeach()
if("${_targetsDefined}" STREQUAL "${_expectedTargets}")
  unset(_targetsDefined)
  unset(_targetsNotDefined)
  unset(_expectedTargets)
  set(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT "${_targetsDefined}" STREQUAL "")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_targetsDefined}\nTargets not yet defined: ${_targetsNotDefined}\n")
endif()
unset(_targetsDefined)
unset(_targetsNotDefined)
unset(_expectedTargets)


# Compute the installation prefix relative to this file.
get_filename_component(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
if(_IMPORT_PREFIX STREQUAL "/")
  set(_IMPORT_PREFIX "")
endif()

# Create imported target itkdouble-conversion
add_library(itkdouble-conversion STATIC IMPORTED)

# Create imported target itksys
add_library(itksys STATIC IMPORTED)

# Create imported target itkvcl
add_library(itkvcl STATIC IMPORTED)

set_target_properties(itkvcl PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/ITK-4.12"
)

# Create imported target itknetlib
add_library(itknetlib STATIC IMPORTED)

set_target_properties(itknetlib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/ITK-4.12"
)

# Create imported target itkv3p_netlib
add_library(itkv3p_netlib STATIC IMPORTED)

set_target_properties(itkv3p_netlib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/ITK-4.12"
)

# Create imported target itkvnl
add_library(itkvnl STATIC IMPORTED)

set_target_properties(itkvnl PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/ITK-4.12"
)

# Create imported target itkvnl_algo
add_library(itkvnl_algo STATIC IMPORTED)

set_target_properties(itkvnl_algo PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/ITK-4.12"
)

# Create imported target itktestlib
add_library(itktestlib STATIC IMPORTED)

set_target_properties(itktestlib PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/ITK-4.12"
)

# Create imported target ITKVNLInstantiation
add_library(ITKVNLInstantiation STATIC IMPORTED)

# Create imported target ITKCommon
add_library(ITKCommon SHARED IMPORTED)

# Create imported target itkNetlibSlatec
add_library(itkNetlibSlatec STATIC IMPORTED)

# Create imported target ITKStatistics
add_library(ITKStatistics SHARED IMPORTED)

# Create imported target ITKTransform
add_library(ITKTransform SHARED IMPORTED)

# Create imported target ITKLabelMap
add_library(ITKLabelMap SHARED IMPORTED)

# Create imported target ITKMesh
add_library(ITKMesh SHARED IMPORTED)

# Create imported target itkzlib
add_library(itkzlib STATIC IMPORTED)

# Create imported target ITKMetaIO
add_library(ITKMetaIO STATIC IMPORTED)

# Create imported target ITKSpatialObjects
add_library(ITKSpatialObjects STATIC IMPORTED)

# Create imported target ITKPath
add_library(ITKPath STATIC IMPORTED)

# Create imported target ITKQuadEdgeMesh
add_library(ITKQuadEdgeMesh SHARED IMPORTED)

# Create imported target ITKIOImageBase
add_library(ITKIOImageBase SHARED IMPORTED)

# Create imported target ITKOptimizers
add_library(ITKOptimizers SHARED IMPORTED)

# Create imported target ITKPolynomials
add_library(ITKPolynomials SHARED IMPORTED)

# Create imported target ITKBiasCorrection
add_library(ITKBiasCorrection SHARED IMPORTED)

# Create imported target ITKBioCell
add_library(ITKBioCell SHARED IMPORTED)

# Create imported target ITKDICOMParser
add_library(ITKDICOMParser STATIC IMPORTED)

# Create imported target ITKEXPAT
add_library(ITKEXPAT STATIC IMPORTED)

# Create imported target ITKIOXML
add_library(ITKIOXML SHARED IMPORTED)

# Create imported target ITKIOSpatialObjects
add_library(ITKIOSpatialObjects STATIC IMPORTED)

# Create imported target ITKFEM
add_library(ITKFEM SHARED IMPORTED)

# Create imported target gdcmCommon
add_library(gdcmCommon STATIC IMPORTED)

# Create imported target gdcmDICT
add_library(gdcmDICT STATIC IMPORTED)

# Create imported target gdcmDSED
add_library(gdcmDSED STATIC IMPORTED)

# Create imported target gdcmIOD
add_library(gdcmIOD STATIC IMPORTED)

# Create imported target gdcmMSFF
add_library(gdcmMSFF STATIC IMPORTED)

# Create imported target gdcmMEXD
add_library(gdcmMEXD STATIC IMPORTED)

# Create imported target gdcmjpeg8
add_library(gdcmjpeg8 STATIC IMPORTED)

# Create imported target gdcmjpeg12
add_library(gdcmjpeg12 STATIC IMPORTED)

# Create imported target gdcmjpeg16
add_library(gdcmjpeg16 STATIC IMPORTED)

# Create imported target gdcmopenjpeg
add_library(gdcmopenjpeg STATIC IMPORTED)

# Create imported target gdcmcharls
add_library(gdcmcharls STATIC IMPORTED)

# Create imported target gdcmsocketxx
add_library(gdcmsocketxx STATIC IMPORTED)

# Create imported target ITKznz
add_library(ITKznz STATIC IMPORTED)

# Create imported target ITKniftiio
add_library(ITKniftiio STATIC IMPORTED)

# Create imported target ITKgiftiio
add_library(ITKgiftiio STATIC IMPORTED)

# Create imported target hdf5-static
add_library(hdf5-static STATIC IMPORTED)

set_target_properties(hdf5-static PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
)

# Create imported target hdf5_cpp-static
add_library(hdf5_cpp-static STATIC IMPORTED)

set_target_properties(hdf5_cpp-static PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include"
)

# Create imported target ITKIOBMP
add_library(ITKIOBMP SHARED IMPORTED)

# Create imported target ITKIOBioRad
add_library(ITKIOBioRad SHARED IMPORTED)

# Create imported target ITKIOCSV
add_library(ITKIOCSV SHARED IMPORTED)

# Create imported target ITKIOGDCM
add_library(ITKIOGDCM SHARED IMPORTED)

# Create imported target ITKIOIPL
add_library(ITKIOIPL SHARED IMPORTED)

# Create imported target ITKIOGE
add_library(ITKIOGE SHARED IMPORTED)

# Create imported target ITKIOGIPL
add_library(ITKIOGIPL SHARED IMPORTED)

# Create imported target ITKIOHDF5
add_library(ITKIOHDF5 SHARED IMPORTED)

# Create imported target itkjpeg
add_library(itkjpeg STATIC IMPORTED)

# Create imported target ITKIOJPEG
add_library(ITKIOJPEG SHARED IMPORTED)

# Create imported target itktiff
add_library(itktiff STATIC IMPORTED)

# Create imported target ITKIOTIFF
add_library(ITKIOTIFF SHARED IMPORTED)

# Create imported target ITKIOLSM
add_library(ITKIOLSM SHARED IMPORTED)

# Create imported target ITKIOMRC
add_library(ITKIOMRC SHARED IMPORTED)

# Create imported target ITKIOMesh
add_library(ITKIOMesh SHARED IMPORTED)

# Create imported target ITKIOMeta
add_library(ITKIOMeta SHARED IMPORTED)

# Create imported target ITKIONIFTI
add_library(ITKIONIFTI SHARED IMPORTED)

# Create imported target ITKNrrdIO
add_library(ITKNrrdIO STATIC IMPORTED)

# Create imported target ITKIONRRD
add_library(ITKIONRRD SHARED IMPORTED)

# Create imported target itkpng
add_library(itkpng STATIC IMPORTED)

# Create imported target ITKIOPNG
add_library(ITKIOPNG SHARED IMPORTED)

# Create imported target ITKIOSiemens
add_library(ITKIOSiemens SHARED IMPORTED)

# Create imported target ITKIOStimulate
add_library(ITKIOStimulate SHARED IMPORTED)

# Create imported target ITKTransformFactory
add_library(ITKTransformFactory STATIC IMPORTED)

# Create imported target ITKIOTransformBase
add_library(ITKIOTransformBase SHARED IMPORTED)

# Create imported target ITKIOTransformHDF5
add_library(ITKIOTransformHDF5 SHARED IMPORTED)

# Create imported target ITKIOTransformInsightLegacy
add_library(ITKIOTransformInsightLegacy SHARED IMPORTED)

# Create imported target ITKIOTransformMatlab
add_library(ITKIOTransformMatlab SHARED IMPORTED)

# Create imported target ITKIOVTK
add_library(ITKIOVTK SHARED IMPORTED)

# Create imported target ITKKLMRegionGrowing
add_library(ITKKLMRegionGrowing SHARED IMPORTED)

# Create imported target ITKOptimizersv4
add_library(ITKOptimizersv4 SHARED IMPORTED)

# Create imported target itkTestDriver
add_executable(itkTestDriver IMPORTED)

# Create imported target ITKVTK
add_library(ITKVTK SHARED IMPORTED)

# Create imported target ITKVideoCore
add_library(ITKVideoCore SHARED IMPORTED)

# Create imported target ITKVideoIO
add_library(ITKVideoIO SHARED IMPORTED)

# Create imported target ITKWatersheds
add_library(ITKWatersheds SHARED IMPORTED)

# Load information for each installed configuration.
get_filename_component(_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(GLOB CONFIG_FILES "${_DIR}/ITKTargets-*.cmake")
foreach(f ${CONFIG_FILES})
  include(${f})
endforeach()

# Cleanup temporary variables.
set(_IMPORT_PREFIX)

# Loop over all imported files and verify that they actually exist
foreach(target ${_IMPORT_CHECK_TARGETS} )
  foreach(file ${_IMPORT_CHECK_FILES_FOR_${target}} )
    if(NOT EXISTS "${file}" )
      message(FATAL_ERROR "The imported target \"${target}\" references the file
   \"${file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
    endif()
  endforeach()
  unset(_IMPORT_CHECK_FILES_FOR_${target})
endforeach()
unset(_IMPORT_CHECK_TARGETS)

# This file does not depend on other imported targets which have
# been exported from the same project but in a separate export set.

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
