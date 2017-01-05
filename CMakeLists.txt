cmake_minimum_required(VERSION 3.0)

project(mathematics_project)

###############################################################################
## file globbing ##############################################################
###############################################################################

# these instructions search the directory tree when cmake is
# invoked and put all files that match the pattern in the variables
# `sources` and `data`
file(GLOB_RECURSE sources                   src/*.cpp)
file(GLOB_RECURSE sources_geometry          src/Geometry/*.cpp)
file(GLOB_RECURSE sources_intersections     src/Intersections/*.cpp)
file(GLOB_RECURSE sources_matrix            src/Matrix/*.cpp)
file(GLOB_RECURSE sources_vector            src/Vector/*.cpp)
file(GLOB_RECURSE sources_helper            src/Helper/*.cpp)

file(GLOB_RECURSE includes                   include/*.h)
file(GLOB_RECURSE includes_geometry          include/Geometry/*.h)
file(GLOB_RECURSE includes_intersections     include/Intersections/*.h)
file(GLOB_RECURSE includes_matrix            include/Matrix/*.h)
file(GLOB_RECURSE includes_vector            include/Vector/*.h)
file(GLOB_RECURSE includes_helper            include/Helper/*.h)


file(GLOB_RECURSE sources_test src/test/*.cpp)
# file(GLOB_RECURSE data resources/*)
# you can use set(sources src/main.cpp) etc if you don't want to
# use globing to find files automatically


###############################################################################
## target definitions #########################################################
###############################################################################

# Allow the developer to select if Dynamic or Static libraries are built
OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
# Set the LIB_TYPE variable to STATIC
SET (LIB_TYPE STATIC)
IF (BUILD_SHARED_LIBS)
  # User wants to build Dynamic Libraries, so change the LIB_TYPE variable to CMake keyword 'SHARED'
  SET (LIB_TYPE SHARED)
ENDIF (BUILD_SHARED_LIBS)

# add the data to the target, so it becomes visible in some IDE
add_library(mathematics ${LIB_TYPE}
  ${includes}
  ${includes_geometry}
  ${includes_intersections}
  ${includes_matrix}
  ${includes_vector}
  ${includes_helper}
  ${sources}
  ${sources_geometry}
  ${sources_intersections}
  ${sources_matrix}
  ${sources_vector}
  ${sources_helper}
  ${data}
)

# just for example add some compiler flags
target_compile_options(mathematics PUBLIC -std=c++13 -Wall)

#MAKE IT EXPORT!!!
#add_definitions(-MATHEMATICS_LIBRARY_DYNAMIC)
add_definitions(-DMATHEMATICS_LIBRARY_STATIC)



# this lets me include files relative to the root src dir with a <> pair
target_include_directories(mathematics PUBLIC include/)

# this copies all resource files in the build directory
# we need this, because we want to work with paths relative to the executable
file(COPY ${data} DESTINATION resources)

###############################################################################
## dependencies ###############################################################
###############################################################################

# this defines the variables Boost_LIBRARIES that contain all library names
# that we need to link to
# find_package(Boost 1.36.0 COMPONENTS filesystem system REQUIRED)

# target_link_libraries(engine PUBLIC
  # ${Boost_LIBRARIES}
  # here you can add any library dependencies
# )

###############################################################################
## testing ####################################################################
###############################################################################

# this is for our testing framework
# we don't add REQUIRED because it's just for testing
find_package(GTest)

if(GTEST_FOUND)
  add_executable(unit_tests ${sources_test} ${sources})

  # we add this define to prevent collision with the main
  # this might be better solved by not adding the source with the main to the
  # testing target
  target_compile_definitions(unit_tests PUBLIC UNIT_TESTS)

  # this allows us to use our executable as a link library
  # therefore we can inherit all compiler options and library dependencies
  set_target_properties(engine PROPERTIES ENABLE_EXPORTS on)

  target_link_libraries(unit_tests PUBLIC
    ${GTEST_BOTH_LIBRARIES}
    engine
  )

  target_include_directories(unit_tests PUBLIC
    ${GTEST_INCLUDE_DIRS} # doesn't do anything on Linux
  )
endif()

###############################################################################
## packaging ##################################################################
###############################################################################

# all install commands get the same destination. this allows us to use paths
# relative to the executable.
install(TARGETS mathematics DESTINATION mathematics_destination)
# this is basically a repeat of the file copy instruction that copies the
# resources in the build directory, but here we tell cmake that we want it
# in the package
install(DIRECTORY resources DESTINATION mathematics_destination)

# now comes everything we need, to create a package
# there are a lot more variables you can set, and some
# you need to set for some package types, but we want to
# be minimal here
set(CPACK_PACKAGE_NAME "Mathematics")
set(CPACK_PACKAGE_VERSION "1.0.0")

# we don't want to split our program up into several things
set(CPACK_MONOLITHIC_INSTALL 1)

# This must be last
include(CPack)
