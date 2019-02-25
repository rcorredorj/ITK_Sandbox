# This is the root ITKSandbox CMakeLists file. 
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT(ITKSandbox)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

# recurse into subdirectories
SUBDIRS(apps)
#-----------------------------------------------------------------------------



