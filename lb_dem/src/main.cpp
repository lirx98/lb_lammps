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
        std::cout<<"f_coupling="<<f_coupling[0][2]<<std::endl;

        // 推进一步
        //lammps_command(lmp, "fix fext all addforce 0.0 0.0 0.000002");
        // lammps_command(lmp, "run 1 pre yes post yes");
        lammps_command(lmp,"run 1 post yes");

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

