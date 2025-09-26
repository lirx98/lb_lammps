#include "lammpsplugin.h"
#include "version.h"
#include "fix_lb_force.h"

using namespace LAMMPS_NS;

// factory: for fix-style plugin factory takes (LAMMPS*, int, char**)
static Fix *lbforce_creator(LAMMPS *lmp, int argc, char **argv) {
  return new FixLBForce(lmp, argc, argv);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc) {
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;
  lammpsplugin_t plugin;

  plugin.version = LAMMPS_VERSION;
  plugin.style   = "fix";
  plugin.name    = "lb/force";               // 你希望在输入脚本中使用的名字： fix ID lb/force ...
  plugin.info    = "LB force fix plugin";
  plugin.author  = "Your Name (youremail@host)";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &lbforce_creator;
  plugin.handle  = handle;

  (*register_plugin)(&plugin, (LAMMPS *) lmp);
}

