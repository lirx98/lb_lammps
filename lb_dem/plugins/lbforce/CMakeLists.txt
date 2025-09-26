cmake_minimum_required(VERSION 3.10)
project(lbforceplugin LANGUAGES CXX)

# 设置C++标准，需与主程序和LAMMPS一致
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 对于共享库，必须启用位置无关代码
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# 查找MPI（插件可能需要MPI支持）
find_package(MPI REQUIRED)

include_directories(    
  /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/src       # 包含 json_fwd.h
)

# 查找LAMMPS
find_path(LAMMPS_INCLUDE_DIR 
  NAMES lammps.h 
  PATHS /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps_install/include/lammps
        # /usr/local/include/lammps 
        # /opt/lammps/include 
        # $ENV{LAMMPS_DIR}/include/lammps
        # /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/STUBS
        # /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/src
)

find_library(LAMMPS_LIBRARY 
  NAMES lammps 
  PATHS /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps_install/lib64
        # /usr/local/lib 
        # /opt/lammps/lib 
        # $ENV{LAMMPS_DIR}/lib64
        # /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/STUBS
        # /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/src
)

if(NOT LAMMPS_INCLUDE_DIR)
  message(FATAL_ERROR "LAMMPS headers not found. Set LAMMPS_DIR environment variable.")
endif()

if(NOT LAMMPS_LIBRARY)
  message(FATAL_ERROR "LAMMPS library not found")
endif()

# 创建插件共享库
add_library(lbforceplugin SHARED
  fix_lb_force.cpp
  fix_lb_force.h
)

# 包含必要的头文件（添加MPI头文件）
target_include_directories(lbforceplugin
  PRIVATE
    ${LAMMPS_INCLUDE_DIR}
    ${MPI_CXX_INCLUDE_DIRS}
)

# 链接LAMMPS和MPI库
target_link_libraries(lbforceplugin
  PRIVATE
    ${LAMMPS_LIBRARY}
    ${MPI_CXX_LIBRARIES}
)

# 添加MPI编译选项
target_compile_options(lbforceplugin
  PRIVATE
    ${MPI_CXX_COMPILE_FLAGS}
)

# 包含必要的头文件
target_include_directories(lbforceplugin
  PRIVATE
    ${LAMMPS_INCLUDE_DIR}
)

# 链接LAMMPS库
target_link_libraries(lbforceplugin
  PRIVATE
    ${LAMMPS_LIBRARY}
)

# # 插件编译选项 - 隐藏内部符号，避免与主程序冲突
# set_target_properties(lbforceplugin PROPERTIES
#   CXX_VISIBILITY_PRESET hidden
#   VISIBILITY_INLINES_HIDDEN ON
#   PREFIX ""  # 生成的库名为lbforceplugin.so而不是liblbforceplugin.so
# )

# 安装插件到系统库目录或LAMMPS插件目录
# set(LAMMPS_PLUGIN_DIR "${CMAKE_INSTALL_PREFIX}/lib/lammps/plugins" CACHE PATH "LAMMPS plugin installation directory")

install(TARGETS lbforceplugin
  LIBRARY DESTINATION ${LAMMPS_PLUGIN_DIR}
  RUNTIME DESTINATION ${LAMMPS_PLUGIN_DIR}  # 适用于Windows
)

# 确保安装目录存在
install(CODE "file(MAKE_DIRECTORY \${CMAKE_INSTALL_PREFIX}/bin)")
install(CODE "file(MAKE_DIRECTORY ${LAMMPS_PLUGIN_DIR})")
