cmake_minimum_required(VERSION 3.20)
project(dem_with_lbforce LANGUAGES CXX)

find_package(MPI REQUIRED)

# LAMMPS 安装路径
set(LAMMPS_DIR          "/work/e283/e283/lirx/test/phase-field/lb_lammps/lammps_install")
set(PLUGIN_DIR          "/mnt/lustre/a2fs-work3/work/e283/e283/lirx/test/phase-field/lb_lammps/my_lbplugin/build")

include_directories(
    ${LAMMPS_DIR}/include/lammps
    ${PLUGIN_DIR}/..
    /work/e283/e283/lirx/test/phase-field/lb_lammps/lammps/src       # 包含 json_fwd.h

)
link_directories(${LAMMPS_DIR}/lib64)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# 定义可执行目标
add_executable(dem_main main.cpp)
# target_link_libraries(dem_main PRIVATE lammps MPI::MPI_CXX dl)

# add_executable(dem_main main.cpp ../my_lbplugin/fix_lb_force.cpp)
target_link_libraries(dem_main PRIVATE lammps MPI::MPI_CXX dl)

# 确定rpath
# set_target_properties(dem_main PROPERTIES
#     BUILD_RPATH ${LAMMPS_DIR}/lib64
#     INSTALL_RPATH ${LAMMPS_DIR}/lib64
# )

# 定义复制文件的自定义命令（修正依赖路径）
add_custom_command(
    OUTPUT "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/lbforceplugin.so"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${PLUGIN_DIR}/lbforceplugin.so"  # 源文件（正确路径）
        "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/lbforceplugin.so"  # 目标文件
    MAIN_DEPENDENCY "${PLUGIN_DIR}/lbforceplugin.so"  # 依赖源文件，而非目标文件
    COMMENT "Copying lbforceplugin.so to output directory"
)

# 创建自定义目标
add_custom_target(
    copy_plugin ALL
    DEPENDS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/lbforceplugin.so"
)

# 关联依赖
add_dependencies(dem_main copy_plugin)

