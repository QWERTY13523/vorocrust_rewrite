cmake_minimum_required(VERSION 3.16)

set(CMAKE_TOOLCHAIN_FILE "/home/yiming/.vcpkg-clion/vcpkg/scripts/buildsystems/vcpkg.cmake"
    CACHE STRING "Vcpkg toolchain file")

# 只有设置了工具链后，才能定义 project
project(vorocrust LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -march=native")

option(USE_CGAL "Enable CGAL" ON)
option(USE_EMBREE "Enable Embree" ON)

find_package(Eigen3 CONFIG REQUIRED)
find_package(libigl CONFIG REQUIRED)
find_package(nanoflann CONFIG REQUIRED)

if (USE_CGAL)
  find_package(CGAL CONFIG REQUIRED)
endif()

if (USE_EMBREE)
  find_package(embree CONFIG REQUIRED)
endif()

set(VOROCRUST_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Generator.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/MeshingTree.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/Methods.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/offset.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/RandomSampler.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/VoronoiMesh.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/inout.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp"
)

add_executable(vorocrust ${VOROCRUST_SOURCES})

target_include_directories(vorocrust PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(vorocrust PRIVATE
  Eigen3::Eigen
  igl::igl_core
  nanoflann::nanoflann
)

# ============================================================
# Release optimization flags
# ============================================================

target_compile_definitions(vorocrust PRIVATE
  $<$<CONFIG:Release>:NDEBUG>
)

if (MSVC)
  target_compile_options(vorocrust PRIVATE
    $<$<CONFIG:Release>:/O2>
    $<$<CONFIG:Release>:/Ob2>
    $<$<CONFIG:Release>:/DNDEBUG>
    $<$<CONFIG:Release>:/arch:AVX2>
  )
else()
  target_compile_options(vorocrust PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Release>:-DNDEBUG>
    $<$<CONFIG:Release>:-march=native>
  )
endif()

# ============================================================
# LTO / IPO optimization
# ============================================================

include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)

if (IPO_SUPPORTED)
  set_property(TARGET vorocrust PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  message(STATUS "IPO/LTO enabled for Release build.")
else()
  message(WARNING "IPO/LTO is not supported: ${IPO_ERROR}")
endif()

# ============================================================
# Optional dependencies
# ============================================================

if (USE_CGAL)
  target_compile_definitions(vorocrust PRIVATE USE_CGAL)
  target_link_libraries(vorocrust PRIVATE CGAL::CGAL)
endif()

if (USE_EMBREE)
  target_compile_definitions(vorocrust PRIVATE USE_EMBREE)
  target_link_libraries(vorocrust PRIVATE embree)
endif()