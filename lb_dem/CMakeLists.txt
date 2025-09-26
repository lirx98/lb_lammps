cmake_minimum_required(VERSION 3.10)
project(lb_dem_project VERSION 1.0 LANGUAGES CXX)



# 设置C++标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 关键：指定所有目标的输出目录为 lb_dem/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)  # 可执行文件输出目录
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)  # 共享库（插件）输出目录
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)  # 静态库输出目录（可选）

# 添加这两行：启用 RPATH 并设置运行时搜索路径
set(CMAKE_SKIP_BUILD_RPATH FALSE)           # 编译时使用 RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)     # 构建时使用安装后的 RPATH
set(CMAKE_INSTALL_RPATH /work/e283/e283/lirx/test/phase-field/lb_lammps/lb_dem/bin)           # 运行时在程序所在目录搜索库（$ORIGIN 表示程序自身目录）
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)  # 保留链接时的库路径

# 控制编译类型
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# 输出编译信息
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CXX flags: ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")

# 查找MPI库（关键：添加MPI支持）
find_package(MPI REQUIRED)
if(MPI_FOUND)
    message(STATUS "MPI found: ${MPI_CXX_LIBRARIES}")
else()
    message(FATAL_ERROR "MPI not found. Please install MPI (e.g., OpenMPI or MPICH)")
endif()

# 先编译插件，再编译主程序
add_subdirectory(plugins/lbforce)
add_subdirectory(src)

# 确保主程序依赖于插件
add_dependencies(lb_dem_main lbforceplugin)
