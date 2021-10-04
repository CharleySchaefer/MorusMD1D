project (MorusMD)
cmake_minimum_required (VERSION 2.6)

set (MORUSMD_VERSION_MAJOR 0)
set (MORUSMD_VERSION_MINOR 1)

include_directories ("lib")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/src")


set(MORUSMD_SRC src/main.c src/printHelp.c)

add_executable (MorusMD ${MORUSMD_SRC})
target_link_libraries(MorusMD m)
