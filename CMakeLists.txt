cmake_minimum_required(VERSION 3.16)

# =========================================================
# 【核心修复】必须在 project() 之前定义工具链文件！
# 否则 CMake 会使用默认编译器设置，忽略 vcpkg
# =========================================================
set(CMAKE_TOOLCHAIN_FILE "D:/vcpkg/scripts/buildsystems/vcpkg.cmake" 
    CACHE STRING "Vcpkg toolchain file")

# 只有设置了工具链后，才能定义 project
project(vorocrust LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(USE_CGAL "Enable CGAL" ON)

find_package(Eigen3 CONFIG REQUIRED)
find_package(libigl CONFIG REQUIRED)

if (USE_CGAL)
  find_package(CGAL CONFIG REQUIRED)
endif()

file(GLOB VOROCRUST_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

list(FILTER VOROCRUST_SOURCES EXCLUDE REGEX "VoronoiMesh1\\.cpp$")

add_executable(vorocrust ${VOROCRUST_SOURCES})

target_include_directories(vorocrust PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(vorocrust PRIVATE 
    Eigen3::Eigen
    igl::igl_core
)

if (USE_CGAL)
  target_compile_definitions(vorocrust PRIVATE USE_CGAL)
  target_link_libraries(vorocrust PRIVATE CGAL::CGAL)
endif()