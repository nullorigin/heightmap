project(libInterpolate)

set(LIB_INTERPOLATE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libInterpolate/src)

add_library(${PROJECT_NAME} INTERFACE)
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} INTERFACE ${LIB_INTERPOLATE_DIR})

target_precompile_headers(
  ${PROJECT_NAME}
  INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/libInterpolate/src/libInterpolate/Interpolate.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/libInterpolate/src/libInterpolate/AnyInterpolator.hpp
)
