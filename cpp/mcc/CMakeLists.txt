cmake_minimum_required(VERSION 3.10)
# Every project needs a name.  We call this the "examples" project.

# Project output alias
set(ProjectOutputName mcc)

# Src files
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# main component, src files
add_library(${ProjectOutputName} STATIC ${SOURCES})

# link libs and include dirs
#

# Properties
set(HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/mcc.h" "${CMAKE_CURRENT_SOURCE_DIR}/fingerprint_base.h")

set_target_properties(${ProjectOutputName} PROPERTIES PUBLIC_HEADER "${HEADERS}")
set_property(TARGET ${ProjectOutputName} PROPERTY CXX_STANDARD 20)
set_property(TARGET ${ProjectOutputName} PROPERTY CXX_STANDARD_REQUIRED ON)

# Installation targets
install(TARGETS ${ProjectOutputName} 
	CONFIGURATIONS Debug 
	DESTINATION ${CMAKE_SOURCE_DIR}/out/Debug
	PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/out/include
)

install(TARGETS ${ProjectOutputName} 
	CONFIGURATIONS Release 
	DESTINATION ${CMAKE_SOURCE_DIR}/out/Release
	PUBLIC_HEADER DESTINATION ${CMAKE_SOURCE_DIR}/out/include
)