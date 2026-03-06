
cmake_minimum_required(VERSION 3.16)
project(vorocrust LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_CGAL "Enable CGAL" ON)

# Configure vcpkg
set(CMAKE_TOOLCHAIN_FILE /home/yiming/.vcpkg-clion/vcpkg/scripts/buildsystems/vcpkg.cmake
    CACHE STRING "Vcpkg toolchain file")

find_package(Eigen3 REQUIRED)

if (USE_CGAL)
  find_package(CGAL CONFIG REQUIRED)
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

target_link_libraries(vorocrust PRIVATE Eigen3::Eigen)

if (USE_CGAL)
  target_compile_definitions(vorocrust PRIVATE USE_CGAL)
  target_link_libraries(vorocrust PRIVATE CGAL::CGAL)
endif()
