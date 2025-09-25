#include <iostream>
#include "lammps.h"
#include "library.h"
#include "modify.h"
#include <string>
#include "mpi.h"
#include "fix_lb_force.h"
#include "domain.h"
#include "atom.h"

using namespace LAMMPS_NS;

void print_all_fixes(LAMMPS *lmp) {
    int nfix = lmp->modify->nfix;
    printf("Total %d fixes found:\n", nfix);

    for (int i = 0; i < nfix; i++) {
        Fix *f = lmp->modify->fix[i];
        if (f) {
            printf("Fix ID: %s, Style: %s\n", f->id, f->style);
        }
    }
}

int main(int argc, char** argv)
{
    // 初始化 LAMMPS（单进程）
    char* args[] = {
        (char*)"liblammps",
        (char*)"-log",(char*)"./log.txt",
        (char*)"-screen",(char*)"none"
    };
    LAMMPS* lmp = (LAMMPS*)lammps_open_no_mpi(5, args, nullptr);
    if (!lmp) {
        std::cerr << "Failed to init LAMMPS\n";
        return -1;
    }

    // 指定输入文件路径（绝对路径或相对 Release 目录）
    lammps_file(lmp, "./input/in.lammps");

    print_all_fixes(lmp);


    std::string fixName="fext";
    int ifix=lmp->modify->find_fix(fixName);
    // std::cout<<fixName.empty()<<std::endl;
    std::cout<<"ifix="<<ifix<<std::endl;

    FixLBForce*coupligFix=dynamic_cast<FixLBForce*>
    (lmp->modify->get_fix_by_id(fixName));    
    double**f_coupling=coupligFix->get_force_ptr();

    std::cout<<"Fix pointer="<<coupligFix<<std::endl;
    std::cout<<"f pointer="<<f_coupling<<std::endl;

    int nsteps = 10;

    for (int step = 0; step < nsteps; step++) {

        int nlocal = lammps_get_natoms(lmp);

        // 提取力数组
        // double** x=(double**)lammps_extract_atom(lmp,"x");
        // double** v=(double**)lammps_extract_atom(lmp,"v");
        // double** f = (double**)lammps_extract_atom(lmp, "f");

        // -------------------------
        // 每步持续叠加外力
        //-------------------------
        f_coupling=coupligFix->get_force_ptr();

        for (int i = 0; i < nlocal; i++) {
            f_coupling[i][0] = 0.0;   // Fx
            f_coupling[i][1] = 0.0;   // Fy
            f_coupling[i][2] = -0.1;  // Fz 向下
        }

        // 推进一步
        //lammps_command(lmp, "fix fext all addforce 0.0 0.0 0.000002");
        // lammps_command(lmp, "run 1 pre yes post yes");
        lammps_command(lmp,"run 1");

        // 每 100 步打印信息
        if (step % 2 == 0) {
            double** x = (double**)lammps_extract_atom(lmp, "x");
            double** v = (double**)lammps_extract_atom(lmp, "v");
            double** f = (double**)lammps_extract_atom(lmp, "f");

            std::cout << "Step: " << step << std::endl;
            for (int i = 0; i < nlocal; i++) {
                std::cout << "Atom " << i
                    << " z: " << x[i][2]
                    << " vz: " << v[i][2]
                    << " fz_ext: " << f[i][2] << std::endl;
            }
        }
    }

    lammps_close(lmp);
    return 0;
}



////#include "lammps_c.h"  // LAMMPS C API核心头文件
//#include <stdlib.h>
//#include "mpi.h"
//#include "library.h"
//#include "lammps.h"
//
//using namespace LAMMPS_NS;
//
//int main(int argc, char** argv) {
//    // 1. 初始化MPI（C API需手动初始化MPI，串行时传MPI_COMM_SELF）
//    MPI_Init(&argc, &argv);
//    MPI_Comm comm = MPI_COMM_WORLD;
//
//    // 2. 创建LAMMPS实例（C API入口函数）
//    void* lmp = lammps_create(0, NULL, comm);
//    if (lmp == NULL) {
//        fprintf(stderr, "Failed to create LAMMPS instance\n");
//        MPI_Abort(comm, 1);
//    }
//
//    // 3. 执行初始化命令（创建模拟框、粒子、力场，与C++ API逻辑一致）
//    lammps_command(lmp, "units metal");                  // 单位体系（需与流体力一致）
//    lammps_command(lmp, "create_box 1 -10 10 -10 10 -10 10");
//    lammps_command(lmp, "create_atoms 1 box");
//    lammps_command(lmp, "pair_style lj/cut 2.5");
//    lammps_command(lmp, "pair_coeff * * 1.0 1.0");
//    lammps_command(lmp, "neighbor 0.3 bin");
//    lammps_command(lmp, "neigh_modify every 1 delay 0 check no");
//    lammps_command(lmp, "fix 1 all nve");              // 时间积分器
//    lammps_command(lmp, "timestep 0.001");
//
//    // 4. 提取粒子核心数据（C API通过lammps_extract_atom获取）
//    int natoms;
//    tagint* tag;          // 粒子全局ID数组（tag[local_idx] = 全局ID）
//    double** f;           // 粒子受力数组（f[local_idx][0=x,1=y,2=z]）
//
//    // 提取粒子数（natoms为输出参数，需先传指针）
//    lammps_extract_atom(lmp, "natoms", &natoms, NULL, NULL);
//    // 提取全局ID数组（类型为"tag"，返回tagint*指针）
//    tag = (tagint*)lammps_extract_atom(lmp, "tag", NULL, NULL);
//    // 提取受力数组（类型为"f"，返回double**指针，三维数组[natoms][3]）
//    f = (double**)lammps_extract_atom(lmp, "f", NULL, NULL);
//
//    // 检查提取结果（避免空指针）
//    if (natoms == 0 || tag == NULL || f == NULL) {
//        fprintf(stderr, "Failed to extract atom data\n");
//        lammps_destroy(lmp);
//        MPI_Abort(comm, 1);
//    }
//    // 示例：流体力存储（模拟流体计算结果，按全局ID映射）
//    typedef struct {
//        tagint global_id;
//        double fx, fy, fz;
//    } FluidForce;
//    FluidForce* fluid_forces;  // 流体计算的所有粒子力（需提前分配内存并填充）
//    int n_fluid_forces;         // 流体力的粒子数量
//
//    // 5. 遍历本地粒子，叠加流体力
//    for (int i = 0; i < natoms; i++) {
//        tagint local_global_id = tag[i];  // 当前本地粒子的全局ID
//        // 查找该粒子对应的流体力（此处简化为线性查找，实际可优化为哈希表）
//        for (int j = 0; j < n_fluid_forces; j++) {
//            if (fluid_forces[j].global_id == local_global_id) {
//                // 叠加流体力到LAMMPS受力数组（总力 = 固体力 + 流体力）
//                f[i][0] += fluid_forces[j].fx;
//                f[i][1] += fluid_forces[j].fy;
//                f[i][2] += fluid_forces[j].fz;
//                break;
//            }
//        }
//    }
//}