#include "fix_lb_force.h"
#include "atom.h"
#include "error.h"
#include "group.h"
#include "utils.h"
#include "comm.h"
#include "update.h"
#include "lammps.h"
#include "library.h"

using namespace LAMMPS_NS;
using namespace FixConst;

// double **f_coupling;
// double **t_coupling;

FixLBForce::FixLBForce(LAMMPS *lmp, int narg, char **arg) 
  : Fix(lmp, narg, arg),nlocal_saved(0),f_coupling(nullptr),t_coupling(nullptr) {
  // : Fix(lmp, narg, arg),nlocal_saved(0) {
  // if (narg < 6) error->all(FLERR,"Illegal fix lb/force command");
  // fx = utils::numeric(FLERR,arg[3],false,lmp);
  // fy = utils::numeric(FLERR,arg[4],false,lmp);
  // fz = utils::numeric(FLERR,arg[5],false,lmp);
  init();
}

FixLBForce::~FixLBForce(){
  deallocate_arrays();
};

int FixLBForce::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

void FixLBForce::init() {
  nlocal_saved=atom->nlocal;
  allocate_arrays();
  for(int i=0;i<nlocal_saved;i++){
    for(int j=0;j<3;j++){
      f_coupling[i][j]=0;
      t_coupling[i][j]=0;
    }
  }
}

void FixLBForce::post_force(int vflag) {
  // double **f = atom->f;
  double **f=(double**)lammps_extract_atom(lmp,"f");
  double **t=atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    // if (mask[i] & groupbit) {
      for(int j=0;j<3;j++){
        f[i][j]+=f_coupling[i][j];
        t[i][j]+=t_coupling[i][j];
      }
    // }
  }
  if (comm->me == 0)
  printf("[FixLBForce] Applied force at step %d, first atom fz=%g + %g\n",
         update->ntimestep, f[0][2], f_coupling[0][2]);

}

int FixLBForce::get_nlocal(){
  return nlocal_saved;
}

double** FixLBForce::get_force_ptr(){
  return f_coupling;
}

double** FixLBForce::get_torque_ptr(){
  return t_coupling;
}

void FixLBForce::allocate_arrays(){
  if(f_coupling||t_coupling){
    deallocate_arrays();
  }
  f_coupling=new double*[nlocal_saved];
  t_coupling=new double*[nlocal_saved];
  for(int i=0;i<nlocal_saved;i++){
    f_coupling[i]=new double[3];
    t_coupling[i]=new double[3];
    memset(f_coupling[i],0,3*sizeof(double));
    memset(t_coupling[i],0,3*sizeof(double));
  }
}

void FixLBForce::deallocate_arrays(){
  if (f_coupling) {
    for (int i = 0; i < nlocal_saved; i++) delete [] f_coupling[i];
    delete [] f_coupling;
    f_coupling = nullptr;
  }
  if (t_coupling) {
    for (int i = 0; i < nlocal_saved; i++) delete [] t_coupling[i];
    delete [] t_coupling;
    t_coupling = nullptr;
  }
}