cmake_minimum_required(VERSION 3.20)
project(lbforceplugin LANGUAGES CXX)

if(NOT DEFINED LAMMPS_SOURCE_DIR)
  message(FATAL_ERROR "Please call cmake with -DLAMMPS_SOURCE_DIR=/path/to/lammps/src or set LAMMPS_SOURCE_DIR")
endif()

# ====== 你的 LAMMPS 安装目录 ======
set(LAMMPS_DIR /home/e283/e283/lirx/lammps_install)

# ====== 头文件路径 ======
include_directories(
  ${LAMMPS_DIR}/include/lammps   # 已安装的 API 头文件 (lammps.h 等)
  ${LAMMPS_SOURCE_DIR}           # 源码目录 (lammpsplugin.h 在这里)
  ${LAMMPS_SOURCE_DIR}/STUBS     # STUBS 目录 (json_fwd.h 等)
)

# ====== 找 LAMMPS 库 ======
find_library(LAMMPS_LIB
  NAMES lammps lammps_mpi lammps_serial
  PATHS
    ${LAMMPS_DIR}/lib64
    ${LAMMPS_DIR}/lib
  NO_DEFAULT_PATH
)

# ====== 插件目标 ======
add_library(lbforceplugin MODULE lbforceplugin.cpp fix_lb_force.cpp)

if(UNIX)
  # 插件必须生成 xxx.so（去掉 lib 前缀）
  set_target_properties(lbforceplugin PROPERTIES PREFIX "" SUFFIX ".so" LINK_FLAGS "-rdynamic")
endif()

# ====== 链接 LAMMPS ======
if(LAMMPS_LIB)
  target_link_libraries(lbforceplugin PRIVATE ${LAMMPS_LIB})
else()
  message(FATAL_ERROR "Could not find liblammps; adjust LAMMPS_DIR or find_library PATH")
endif()

