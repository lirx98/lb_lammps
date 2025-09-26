cmake_minimum_required(VERSION 3.10)
project(lb_dem_src LANGUAGES CXX)

include_directories(    
  /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/src       # 包含 json_fwd.h
)

# 查找LAMMPS库和头文件
find_path(LAMMPS_INCLUDE_DIR 
  NAMES lammps.h 
  PATHS /usr/local/include/lammps 
        /opt/lammps/include 
        $ENV{LAMMPS_DIR}/src
        /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps_install/include/lammps
        /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/STUBS
)

find_library(LAMMPS_LIBRARY 
  NAMES lammps 
  PATHS /usr/local/lib 
        /opt/lammps/lib 
        $ENV{LAMMPS_DIR}/build
        /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps_install/lib64
        /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/STUBS
)

if(NOT LAMMPS_INCLUDE_DIR OR NOT LAMMPS_LIBRARY)
  message(FATAL_ERROR "LAMMPS library or headers not found. Set LAMMPS_DIR environment variable.")
endif()

# 创建主程序可执行文件
add_executable(lb_dem_main main.cpp)

# 包含MPI和LAMMPS头文件
target_include_directories(lb_dem_main
  PRIVATE
    ${LAMMPS_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/../plugins/lbforce  # 插件头文件
)

# 链接MPI、LAMMPS库和插件
target_link_libraries(lb_dem_main
  PRIVATE
    ${LAMMPS_LIBRARY}
    ${MPI_CXX_LIBRARIES}
    lbforceplugin  # 链接插件库
)

# 设置MPI编译选项
target_compile_options(lb_dem_main
  PRIVATE
    ${MPI_CXX_COMPILE_FLAGS}
)

# 包含LAMMPS头文件目录
target_include_directories(lb_dem_main 
  PRIVATE 
    ${LAMMPS_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/plugins/lbforce  # 包含插件头文件以便在main中使用
)

# 链接LAMMPS库
target_link_libraries(lb_dem_main 
  PRIVATE 
    ${LAMMPS_LIBRARY}
)

# 设置运行时动态库搜索路径
set_target_properties(lb_dem_main PROPERTIES
  INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${LAMMPS_LIBRARY_DIR}"
  BUILD_WITH_INSTALL_RPATH TRUE
)

# 安装主程序
install(TARGETS lb_dem_main 
  RUNTIME DESTINATION bin
)
