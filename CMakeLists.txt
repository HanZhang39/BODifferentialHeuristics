project (BiObjectiveDH)

cmake_minimum_required (VERSION 2.6)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
SET (CMAKE_LINKER  "/usr/bin/ld")
# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -Wall -fsanitize=address")

include_directories("inc")

file(GLOB SOURCES "src/BiCriteria/*.cpp" "src/Utils/*.cpp")

# Find Boost
find_package( Boost REQUIRED COMPONENTS program_options log thread filesystem system)
include_directories( ${Boost_INCLUDE_DIRS} )


message(${CMAKE_BINARY_DIR})

add_executable(solver src/driver.cpp ${SOURCES})
add_executable(bod src/single_to_all.cpp ${SOURCES})
add_executable(solver_dh src/dh_driver.cpp ${SOURCES})
add_executable(compress src/compress.cpp ${SOURCES})

target_link_libraries(solver ${Boost_LIBRARIES})
target_link_libraries(bod ${Boost_LIBRARIES})
target_link_libraries(solver_dh ${Boost_LIBRARIES})
target_link_libraries(compress ${Boost_LIBRARIES})
