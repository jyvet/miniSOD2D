project(miniSOD2D)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.f90)

include_directories(${libfortran_INCLUDE_DIRS})

add_definitions(-D__BUILD_VERSION__="${TEST_VERSION}")

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${PROJECT_NAME} tmu_utils tmu_geometry tmu_fem tmu_kernels)

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "app_F90" VERSION ${TEST_VERSION} SOVERSION ${TEST_VERSION_MAJOR})

install(TARGETS ${PROJECT_NAME} DESTINATION bin)