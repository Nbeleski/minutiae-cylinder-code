#
# CMake root CMakeLists.txt
#
cmake_minimum_required (VERSION 3.10)

project ("mcc-cpp")

#
# Force msvc to compile in MT
#
if (POLICY CMP0054)
   cmake_policy(SET CMP0054 NEW)
endif()

if (MSVC OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
      foreach(flag_var
         CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
         CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
         if(${flag_var} MATCHES "/MD")
            string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
         endif()
      endforeach(flag_var)
endif()

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Include sub-projects.
add_subdirectory ("mcc")       # build the static libs with opencv / openvino
add_subdirectory ("mcc-cli")	# build cli app, linked to the static lib

# Properties
# ---

# Install
# ---