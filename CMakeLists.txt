cmake_minimum_required(VERSION 3.14)
project(micro-c)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

aux_source_directory(${CMAKE_CURRENT_LIST_DIR}/src sources)

foreach(source ${sources})
  get_filename_component(binary ${source} NAME_WE)
  add_executable(${binary} ${source})
  target_include_directories(${binary} PRIVATE ${CMAKE_CURRENT_LIST_DIR}/include)
  target_link_libraries(${binary})
endforeach()
