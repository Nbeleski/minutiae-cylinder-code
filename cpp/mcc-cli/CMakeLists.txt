cmake_minimum_required(VERSION 3.10)
# Every project needs a name.  We call this the "examples" project.

# Project output alias
set(ProjectOutputName mcc-cli)

# Add source to this project's executable.
add_executable (${ProjectOutputName} main.cpp base64.cpp)

# Link dependencies
target_link_libraries(${ProjectOutputName} PUBLIC mcc)

# Add dependencies headers
target_include_directories(${ProjectOutputName} PRIVATE $<TARGET_PROPERTY:mcc,INCLUDE_DIRECTORIES>)

# Set properties
set_property(TARGET ${ProjectOutputName} PROPERTY CXX_STANDARD 20)
set_property(TARGET ${ProjectOutputName} PROPERTY CXX_STANDARD_REQUIRED ON)

# Installation targets
install(TARGETS ${ProjectOutputName} 
	CONFIGURATIONS Debug 
	DESTINATION ${CMAKE_SOURCE_DIR}/out/Debug
)

install(TARGETS ${ProjectOutputName} 
	CONFIGURATIONS Release 
	DESTINATION ${CMAKE_SOURCE_DIR}/out/Release
)