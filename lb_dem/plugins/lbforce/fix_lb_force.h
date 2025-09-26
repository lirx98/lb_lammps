#ifdef FIX_CLASS
FixStyle(lb/force,FixLBForce)
#else

#ifndef LMP_FIX_LB_FORCE_H
#define LMP_FIX_LB_FORCE_H

#include "fix.h"

// visibility macro define
#ifdef _WIN32
#define PLUGIN_EXPORT __declspec(dllexport)
#else
#define PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

namespace LAMMPS_NS {

class FixLBForce : public Fix {
 public:
  FixLBForce(class LAMMPS *, int, char **);
  ~FixLBForce();

  int setmask();
  void init();
  void post_force(int);
  int get_nlocal();

  // Return pointer of force and torque
  double** get_force_ptr();
  double** get_torque_ptr();

  // Allocate and deallocate 
  void allocate_arrays();
  void deallocate_arrays();

 private:
  int nlocal_saved;
  double **f_coupling;
  double **t_coupling;
};

}

#endif
#endif

