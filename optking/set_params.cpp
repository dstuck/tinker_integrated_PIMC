/*! \file set_params.cc
    \ingroup optking
    \brief set optimization parameters
*/

#define EXTERN
#include "globals.h"
//DES Temp:
#include <iostream>

//#if defined(OPTKING_PACKAGE_PSI)
//#include <psi4-dec.h>
//#elif defined(OPTKING_PACKAGE_QCHEM)
//#include <qchem.h>
//#endif

namespace opt {

void print_params(void);

//#if defined(OPTKING_PACKAGE_PSI)
//void set_params(psi::Options & options)
//#else
void set_params(void)
//#endif
{

  std::string s;

//#if defined(OPTKING_PACKAGE_PSI)
//
//// optimization type
//    s = options.get_str("OPT_TYPE");
//    if(s == "MIN")  Opt_params.opt_type = OPT_PARAMS::MIN;
//    else if (s == "TS")  Opt_params.opt_type = OPT_PARAMS::TS;
//    else if (s == "IRC")  Opt_params.opt_type = OPT_PARAMS::IRC;
//
//    if (options["OPT_TYPE"].has_changed()) {
//      s = options.get_str("STEP_TYPE");
//      if (s == "RFO") {
//        if (Opt_params.opt_type == OPT_PARAMS::MIN)
//          Opt_params.step_type = OPT_PARAMS::RFO;
//        else if (Opt_params.opt_type == OPT_PARAMS::TS)
//          Opt_params.step_type = OPT_PARAMS::P_RFO;
//      }
//      else if (s == "NR") Opt_params.step_type = OPT_PARAMS::NR;
//   }
//   else { // set defaults for step type
//     if (Opt_params.opt_type == OPT_PARAMS::MIN)
//       Opt_params.step_type = OPT_PARAMS::RFO;
//     else if (Opt_params.opt_type == OPT_PARAMS::TS)
//       Opt_params.step_type = OPT_PARAMS::P_RFO;
//     // else if (Opt_params.opt_type == OPT_PARAMS::IRC) options?
//   }
//
//// Maximum step size in bohr or radian along an internal coordinate {double}
//    Opt_params.intrafragment_step_limit = options.get_double("INTRAFRAGMENT_STEP_LIMIT");
//
//// Whether to 'follow' the initial RFO vector after the first step {true, false}
//    Opt_params.rfo_follow_root = options.get_bool("RFO_FOLLOW_ROOT");
//
//// Which RFO root to follow; internally 0 (externally 1) indicates minimum; {integer}
//    Opt_params.rfo_root = options.get_int("RFO_ROOT");
//
//// When determining connectivity, a bond is assigned if interatomic distance
//// is less than (this number) * sum of covalent radii {double}
//    Opt_params.scale_connectivity = options.get_double("SCALE_CONNECTIVITY");
//
//// Whether to treat multiple molecule fragments as a single bonded molecule;
//// or via interfragment coordinates ; a primary difference is that in MULTI mode,
//// the interfragment coordinates are not redundant. {SINGLE, MULTI}
//    s = options.get_str("FRAGMENT_MODE");
//    if (s == "SINGLE")     Opt_params.fragment_mode = OPT_PARAMS::SINGLE;
//    else if (s == "MULTI") Opt_params.fragment_mode = OPT_PARAMS::MULTI;
//
//// whether to use fixed linear combinations of atoms as reference points for
////interfragment coordinates or whether to use principal axes {FIXED, PRINCIPAL_AXES}
//    s = options.get_str("INTERFRAGMENT_MODE");
//    if (s == "FIXED")               Opt_params.interfragment_mode = OPT_PARAMS::FIXED;
//    else if (s == "PRINCIPAL_AXES") Opt_params.interfragment_mode = OPT_PARAMS::PRINCIPAL_AXES;
//
//// Whether to only generate the internal coordinates and then stop {true, false}
//    Opt_params.generate_intcos_only = options.get_bool("GENERATE_INTCOS_ONLY");
//
//// What model Hessian to use to guess intrafragment force constants {SCHLEGEL, FISCHER}
//    s = options.get_str("INTRAFRAGMENT_H");
//    if (s == "FISCHER")       Opt_params.intrafragment_H = OPT_PARAMS::FISCHER;
//    else if (s == "SCHLEGEL") Opt_params.intrafragment_H = OPT_PARAMS::SCHLEGEL;
//
//// Whether to use the default of FISCHER_LIKE force constants for the initial guess {DEFAULT, FISCHER_LIKE}
//    s = options.get_str("INTERFRAGMENT_H");
//    if (s == "DEFAULT")           Opt_params.interfragment_H = OPT_PARAMS::DEFAULT;
//    else if (s == "FISCHER_LIKE") Opt_params.interfragment_H = OPT_PARAMS::FISCHER_LIKE;
//
//// Whether to freeze all fragments rigid
//    Opt_params.freeze_intrafragment = options.get_bool("FREEZE_INTRAFRAGMENT");
//
//// By default, optking prints and saves the last (previous) geometry at the end of an
//// optimization, i.e., the one at which a gradient was computed.  If this keyword is
//// set to true, then the structure obtained from the last projected step is printed out and saved instead.
//    Opt_params.write_final_step_geometry = options.get_bool("WRITE_FINAL_STEP_GEOMETRY");
//
//// Choose from supported Hessian updates {NONE, BFGS, MS, POWELL, BOFILL}
//    s = options.get_str("H_UPDATE");
//    if (s == "NONE")        Opt_params.H_update = OPT_PARAMS::NONE;
//    else if (s == "BFGS")   Opt_params.H_update = OPT_PARAMS::BFGS;
//    else if (s == "MS")     Opt_params.H_update = OPT_PARAMS::MS;
//    else if (s == "POWELL") Opt_params.H_update = OPT_PARAMS::POWELL;
//    else if (s == "BOFILL") Opt_params.H_update = OPT_PARAMS::BOFILL;
//
////  How many previous steps' data to use in Hessian update; 0=use them all ; {integer}
//    Opt_params.H_update_use_last = options.get_int("H_UPDATE_USE_LAST");
//
//// Whether to limit the magnitutde of changes caused by the Hessian update {true, false}
//    Opt_params.H_update_limit = options.get_bool("H_UPDATE_LIMIT");
//
//// If the above is true, changes to the Hessian from the update are limited to the larger of
//// (H_update_limit_scale)*(the previous value) and H_update_limit_max (in au).
//    Opt_params.H_update_limit_scale = options.get_double("H_UPDATE_LIMIT_SCALE");
//    Opt_params.H_update_limit_max = options.get_double("H_UPDATE_LIMIT_MAX");
//
//// Whether to use 1/R(AB) for stretching coordinate between fragments (or just R(AB))
//    Opt_params.interfragment_distance_inverse = options.get_bool("INTERFRAGMENT_DISTANCE_INVERSE");
//
//// For now, this is a general maximum distance for the definition of H-bonds
//    Opt_params.maximum_H_bond_distance = options.get_double("MAXIMUM_H_BOND_DISTANCE");
//
//// QCHEM optimization criteria
//    Opt_params.conv_max_force = options.get_double("CONV_MAX_FORCE");
//    Opt_params.conv_max_DE = options.get_double("CONV_MAX_DE");
//    Opt_params.conv_max_disp = options.get_double("CONV_MAX_DISP");
//
//// Whether to test B matrix and derivative B matrix numerically
//    Opt_params.test_B = options.get_bool("TEST_B");
//    Opt_params.test_derivative_B = options.get_bool("TEST_DERIVATIVE_B");
//
//    Opt_params.print_lvl = options.get_int("PRINT");
//
//// read cartesian Hessian
//    Opt_params.read_cartesian_H = options.get_bool("READ_CARTESIAN_H");
//
//// only treating "dummy fragments"
//    // These are not found in psi4/read_options.cc
//    // Not sure if we need these.
//  Opt_params.efp_fragments = false;
//  Opt_params.efp_fragments_only = false;
//
//#elif defined(OPTKING_PACKAGE_QCHEM)

  int i;

// MIN = 0 ; TS = 1 ; IRC = 2      (default 0)
//  i = rem_read(REM_GEOM_OPT2_OPT_TYPE);
   i = 0;
  if (i == 0)      Opt_params.step_type = OPT_PARAMS::RFO;
  else if (i == 1) Opt_params.step_type = OPT_PARAMS::NR;
  else if (i == 2) Opt_params.step_type = OPT_PARAMS::P_RFO;

// RFO = 0 ; NR = 1 ; P_RFO = 2      (default 0)
  //i = rem_read(REM_GEOM_OPT2_STEP_TYPE);
   i=0;
  if (i == 0)      Opt_params.step_type = OPT_PARAMS::RFO;
  else if (i == 1) Opt_params.step_type = OPT_PARAMS::NR;
  else if (i == 2) Opt_params.step_type = OPT_PARAMS::P_RFO;

// max = i / 10  (default 4)
  //i = rem_read(REM_GEOM_OPT2_INTRAFRAGMENT_STEP_LIMIT);
   i=4;
  Opt_params.intrafragment_step_limit = i / 10.0;

// Whether to follow root   (default 0 means 'no')
  //Opt_params.rfo_follow_root = rem_read(REM_GEOM_OPT2_RFO_FOLLOW_ROOT);

// Which root is desired    (default 0 for minimum)
  //Opt_params.rfo_root = rem_read(REM_GEOM_OPT2_RFO_ROOT);

// scale = i / 10   (default 13)
  //i = rem_read(REM_GEOM_OPT2_SCALE_CONNECTIVITY);
//DES: Set below default!
   i=10;
  Opt_params.scale_connectivity = i / 10.0;

// multi-mode (0=single ; 1= multi) (default 0)
  //i = rem_read(REM_GEOM_OPT2_FRAGMENT_MODE);
   //i=0;
   i=1;
  if (i == 0) {
     //std::cout << "DES: setting fragment_mode=SINGLE" << std::endl;
     Opt_params.fragment_mode = OPT_PARAMS::SINGLE;
  }
  else if (i == 1) { 
     //std::cout << "DES: setting fragment_mode=MULTI" << std::endl;
     Opt_params.fragment_mode = OPT_PARAMS::MULTI;
  }

// interfragment mode type (0=fixed; 1=principal axes)
  //i = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_MODE);
   //i=1;
   i=0;
  if (i == 0)      Opt_params.interfragment_mode = OPT_PARAMS::FIXED;
  else if (i == 1) Opt_params.interfragment_mode = OPT_PARAMS::PRINCIPAL_AXES;

// only generate intcos
  //Opt_params.generate_intcos_only = rem_read(REM_GEOM_OPT2_GENERATE_INTCOS_ONLY);

// model 0=FISCHER ; 1 = SCHLEGEL (default 0)
  //i = rem_read(REM_GEOM_OPT2_INTRAFRAGMENT_H);
   i=0;
  if (i == 0)      Opt_params.intrafragment_H = OPT_PARAMS::FISCHER;
  else if (i == 1) Opt_params.intrafragment_H = OPT_PARAMS::SCHLEGEL;

// interfragment 0=DEFAULT ; 1=FISCHER_LIKE
  //i = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_H);
   i=0;
  if (i == 0)      Opt_params.interfragment_H = OPT_PARAMS::DEFAULT;
  else if (i == 1) Opt_params.interfragment_H = OPT_PARAMS::FISCHER_LIKE;

// Whether to freeze all fragments rigid (default 0);
  //Opt_params.freeze_intrafragment = rem_read(REM_GEOM_OPT2_FREEZE_INTRAFRAGMENT);
  Opt_params.freeze_intrafragment = 0;

// write final step ; default (false);
  //Opt_params.write_final_step_geometry = rem_read(REM_GEOM_OPT2_WRITE_FINAL_STEP_GEOMETRY);
  Opt_params.write_final_step_geometry = false;

// {NONE, BFGS, MS, POWELL, BOFILL} (default 1)
  //i = rem_read(REM_GEOM_OPT2_H_UPDATE);
   i=1;
  if (i == 0) Opt_params.H_update = OPT_PARAMS::NONE;
  else if (i == 1) Opt_params.H_update = OPT_PARAMS::BFGS;
  else if (i == 2) Opt_params.H_update = OPT_PARAMS::MS;
  else if (i == 3) Opt_params.H_update = OPT_PARAMS::POWELL;
  else if (i == 4) Opt_params.H_update = OPT_PARAMS::BOFILL;

// previous steps to use ; (0=all) ; default (6)
  //Opt_params.H_update_use_last = rem_read(REM_GEOM_OPT2_H_UPDATE_USE_LAST);
  Opt_params.H_update_use_last = 6;

// limit hessian changes (default true)
  //Opt_params.H_update_limit = rem_read(REM_GEOM_OPT2_H_UPDATE_LIMIT);

// scale is i / 10 (default 5 -> 0.5)  ; max is i / 10 (default 10 -> 1.0)
//  i = rem_read(REM_GEOM_OPT2_H_UPDATE_LIMIT_SCALE);
   i=5;
  Opt_params.H_update_limit_scale = i / 10;

  //i = rem_read(REM_GEOM_OPT2_H_UPDATE_LIMIT_MAX);
   i=10;
  Opt_params.H_update_limit_max  = i / 10;

// use 1/R(AB) ; (default 0)
  //Opt_params.interfragment_distance_inverse = rem_read(REM_GEOM_OPT2_INTERFRAGMENT_DISTANCE_INVERSE);
  Opt_params.interfragment_distance_inverse = 0;

// H-bond distance = i / 10 ; default (43 -> 4.3)
  //i = rem_read(REM_GEOM_OPT2_MAXIMUM_H_BOND_DISTANCE);
   i=43;
  Opt_params.maximum_H_bond_distance = i / 10;

// QCHEM optimization criteria ; names adapted from old QCHEM
  //i = rem_read(REM_GEOM_OPT2_TOL_GRADIENT);
   i=300;
  Opt_params.conv_max_force = i / 1.0e6; // default (300 -> 3e-4)

  //i = rem_read(REM_GEOM_OPT2_TOL_DISPLACEMENT);
   i = 1200;
  Opt_params.conv_max_disp  = i / 1.0e6; // default (1200 -> 1.2e-3)

  //i = rem_read(REM_GEOM_OPT2_TOL_ENERGY);
   i=100;
  Opt_params.conv_max_DE    = i / 1.0e8; // default (100 -> 1.0e-6)

// test B (default 0)
  //Opt_params.test_B = rem_read(REM_GEOM_OPT2_TEST_B);
  Opt_params.test_B = 0;
  //Opt_params.test_derivative_B = rem_read(REM_GEOM_OPT2_TEST_DERIVATIVE_B);
  Opt_params.test_derivative_B = 0;

// (default 1)
  //Opt_params.print_lvl = rem_read(REM_GEOM_OPT2_PRINT_LVL);
  Opt_params.print_lvl = 1;

// read Hessian (default 0)
  //Opt_params.read_cartesian_H = rem_read(REM_GEOM_OPT2_READ_CARTESIAN_H);
  Opt_params.read_cartesian_H = 0;

  Opt_params.efp_fragments = 0;

// are ONLY EFP fragments present
  if(Opt_params.efp_fragments)
    Opt_params.efp_fragments_only = false;
  else {
    Opt_params.efp_fragments_only = false;
  }

//#endif

// ** Items are below unlikely to need modified

// For RFO step, eigenvectors of augmented Hessian are divided by the last
// element unless it is smaller than this value {double}
  Opt_params.rfo_normalization_min = 1.0e-8;

// how close to pi should a torsion be to assume it may have passed through 180
  Opt_params.fix_tors_near_pi = 2.618; // ~150 degrees

// torsional angles will not be computed if the contained bond angles are within
// this fraction of pi from 0 or from pi
  Opt_params.tors_angle_lim = 0.005;

// only used for determining which atoms in a fragment are acceptable for use
// as reference atoms.  We avoid collinear sets.
// angle is 0/pi if the bond angle is within this fraction of pi from 0/pi
  Opt_params.interfrag_collinear_tol = 0.01;

// cos(torsional angle) must be this close to -1/+1 for angle to count as 0/pi
  Opt_params.tors_cos_tol = 1e-10;

// if bend exceeds this value, then also create linear bend complement
  Opt_params.linear_bend_threshold = 3.05; // about 175 degrees

// threshold for which entries in diagonalized redundant matrix are kept and inverted
// while computing a generalized inverse of a matrix
  Opt_params.redundant_eval_tol = 1.0e-10;

// Parameters that control the backtransformation to cartesian coordinates
  Opt_params.bt_max_iter = 25;
  Opt_params.bt_dx_conv = 1.0e-6;
  Opt_params.bt_dx_conv_rms_change = 1.0e-12;

// Hessian update is avoided if the denominators (Dq*Dq) or (Dq*Dg) are smaller than this
  Opt_params.H_update_den_tol = 1e-7;

// Hessian update is avoided if any internal coordinate has changed by more than this in radians/au
  Opt_params.H_update_dq_tol = 0.5;

// Some parameter error-checking / modification
  if (Opt_params.efp_fragments_only) {
    Opt_params.test_B = false;
    Opt_params.test_derivative_B = false;
  }

  if (Opt_params.print_lvl > 1) print_params();

}

void print_params(void) {

  if (Opt_params.opt_type == OPT_PARAMS::MIN)
  fprintf(outfile, "opt_type              = %18s\n", "Minimum-energy search");
  else if (Opt_params.opt_type == OPT_PARAMS::TS)
  fprintf(outfile, "opt_type              = %18s\n", "Transition-state search");
  else if (Opt_params.opt_type == OPT_PARAMS::IRC)
  fprintf(outfile, "opt_type              = %18s\n", "IRC path following");

  if (Opt_params.step_type == OPT_PARAMS::NR)
  fprintf(outfile, "step_type              = %18s\n", "N-R");
  else if (Opt_params.step_type == OPT_PARAMS::RFO)
  fprintf(outfile, "step_type              = %18s\n", "RFO");
  else if (Opt_params.step_type == OPT_PARAMS::P_RFO)
  fprintf(outfile, "step_type              = %18s\n", "P_RFO");

  fprintf(outfile, "conv_max_force         = %18.2e\n", Opt_params.conv_max_force);
  fprintf(outfile, "conv_max_DE            = %18.2e\n", Opt_params.conv_max_DE);
  fprintf(outfile, "conv_max_disp          = %18.2e\n", Opt_params.conv_max_disp);

  fprintf(outfile, "scale_connectivity     = %18.2e\n", Opt_params.scale_connectivity);

  if (Opt_params.fragment_mode == OPT_PARAMS::SINGLE)
  fprintf(outfile, "fragment_mode          = %18s\n", "single");
  else if (Opt_params.fragment_mode == OPT_PARAMS::MULTI)
  fprintf(outfile, "fragment_mode          = %18s\n", "multi");

  if (Opt_params.interfragment_mode == OPT_PARAMS::FIXED)
  fprintf(outfile, "interfragment_mode        = %18s\n", "fixed");
  else if (Opt_params.interfragment_mode == OPT_PARAMS::PRINCIPAL_AXES)
  fprintf(outfile, "interfragment_mode        = %18s\n", "principal axes");

  if (Opt_params.generate_intcos_only)
  fprintf(outfile, "generate_intcos_only   = %18s\n", "true");
  else
  fprintf(outfile, "generate_intcos_only   = %18s\n", "false");

  if (Opt_params.rfo_follow_root)
  fprintf(outfile, "rfo_follow_root        = %18s\n", "true");
  else
  fprintf(outfile, "rfo_follow_root        = %18s\n", "false");

  fprintf(outfile, "rfo_root               = %18d\n", Opt_params.rfo_root);

  fprintf(outfile, "rfo_normalization_min  = %18.2e\n", Opt_params.rfo_normalization_min);

  if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER)
  fprintf(outfile, "intrafragment_H        = %18s\n", "Fischer");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL)
  fprintf(outfile, "intrafragment_H        = %18s\n", "Schlegel");

  if (Opt_params.interfragment_H == OPT_PARAMS::DEFAULT)
  fprintf(outfile, "interfragment_H        = %18s\n", "Default");
  else if (Opt_params.interfragment_H == OPT_PARAMS::FISCHER_LIKE)
  fprintf(outfile, "interfragment_H        = %18s\n", "Fischer_like");

  if (Opt_params.H_update == OPT_PARAMS::NONE)
  fprintf(outfile, "H_update               = %18s\n", "None");
  else if (Opt_params.H_update == OPT_PARAMS::BFGS)
  fprintf(outfile, "H_update               = %18s\n", "BFGS");
  else if (Opt_params.H_update == OPT_PARAMS::MS)
  fprintf(outfile, "H_update               = %18s\n", "MS");
  else if (Opt_params.H_update == OPT_PARAMS::POWELL)
  fprintf(outfile, "H_update               = %18s\n", "Powell");
  else if (Opt_params.H_update == OPT_PARAMS::BOFILL)
  fprintf(outfile, "H_update               = %18s\n", "Bofill");

  fprintf(outfile, "H_update_use_last      = %18d\n", Opt_params.H_update_use_last);

  if (Opt_params.freeze_intrafragment)
  fprintf(outfile, "freeze_intrafragment   = %18s\n", "true");
  else
  fprintf(outfile, "freeze_intrafragment   = %18s\n", "false");

  fprintf(outfile, "intrafragment_step_limit=%18.2e\n", Opt_params.intrafragment_step_limit);

  if (Opt_params.H_update_limit == OPT_PARAMS::DEFAULT)
  fprintf(outfile, "H_update_limit         = %18s\n", "true");
  else if (Opt_params.H_update_limit == OPT_PARAMS::FISCHER_LIKE)
  fprintf(outfile, "H_update_limit         = %18s\n", "false");

  fprintf(outfile, "H_update_limit_scale   = %18.2e\n", Opt_params.H_update_limit_scale);
  fprintf(outfile, "H_update_limit_max     = %18.2e\n", Opt_params.H_update_limit_max);

  if (Opt_params.interfragment_distance_inverse)
  fprintf(outfile, "interfragment_distance_inverse=%12s\n", "true");
  else
  fprintf(outfile, "interfragment_distance_inverse=%12s\n", "false");

  if (Opt_params.write_final_step_geometry)
  fprintf(outfile, "write_final_step_geometry= %16s\n", "true");
  else
  fprintf(outfile, "write_final_step_geometry= %16s\n", "false");

  fprintf(outfile, "maximum_H_bond_distance= %18.2e\n", Opt_params.maximum_H_bond_distance);

  if (Opt_params.read_cartesian_H)
  fprintf(outfile, "read_cartesian_H       = %18s\n", "true");
  else
  fprintf(outfile, "read_cartesian_H       = %18s\n", "false");

  if (Opt_params.efp_fragments)
  fprintf(outfile, "efp_fragments          = %18s\n", "true");
  else
  fprintf(outfile, "efp_fragments          = %18s\n", "false");

  if (Opt_params.efp_fragments_only)
  fprintf(outfile, "efp_fragments_only     = %18s\n", "true");
  else
  fprintf(outfile, "efp_fragments_only     = %18s\n", "false");

}

}

