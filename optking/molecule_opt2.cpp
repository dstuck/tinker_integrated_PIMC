/*! \file molecule.cc
    \ingroup optking
    \brief molecule class (really, molecular system class)
*/

#include "molecule.h"

#include <iostream>
#include <sstream>

#include "linear_algebra.h"
#include "print.h"
#include "atom_data.h"
#include "physconst.h"

#define EXTERN
#include "globals.h"

//#if defined(OPTKING_PACKAGE_PSI)
//#include <cmath>
//#elif defined(OPTKING_PACKAGE_QCHEM)
#include "cmath"
//#include "qcmath.h"
//#include "EFP.h"
//#endif

namespace opt {

using namespace std;

// Adds one fragment with memory to accommodate the given number of atoms.
// Add EFP fragments which contains no atoms.
MOLECULE::MOLECULE(int num_atoms) {

  if (num_atoms > 0) {
    FRAG *one_frag = new FRAG(num_atoms);
    fragments.push_back(one_frag);
  }

#if defined(OPTKING_PACKAGE_QCHEM)
//  fprintf(outfile,"\tCreating molecule with %d real atoms.\n", num_atoms);
//  if (REM_EFP) {
//    int num_efp = EFP::GetInstance()->NFragments();
//    
//    for (int i=0; i<num_efp; ++i) {
//      EFP_FRAG *one_efp = new EFP_FRAG();
//      efp_fragments.push_back(one_efp);
//    }
//    fprintf(outfile,"\tAdding %d EFP fragments.\n", num_efp);
//  }
#endif

  return;
}


// compute forces in internal coordinates in au
// forces in internal coordinates, f_q = G_inv B u f_x
// if u is unit matrix, f_q = (BB^T)^(-1) * B f_x
void MOLECULE::forces(void) {
  double *f_x, *temp_arr, **B, **G, **G_inv;
  int Ncart = 3*g_natom();
  int Nintco = g_nintco();

  B = compute_B();
  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile, "B matrix\n");
    print_matrix(outfile, B, Nintco, Ncart);
  }

  f_x = g_grad_array(); // Hartree / bohr
  array_scm(f_x, -1, Ncart); // switch gradient -> forces

  temp_arr = init_array(Nintco);
  opt_matrix_mult(B, 0, &f_x, 1, &temp_arr, 1, Nintco, Ncart, 1, 0);
  free_array(f_x);

  G = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

  G_inv = symm_matrix_inv(G, Nintco, 1);
  free_matrix(G);

  double * f_q = p_Opt_data->g_forces_pointer();
  opt_matrix_mult(G_inv, 0, &temp_arr, 1, &f_q, 1, Nintco, Nintco, 1, 0);
  free_matrix(G_inv);
  free_array(temp_arr);

#if defined(OPTKING_PACKAGE_QCHEM)
//  // append exernally determiend efp forces
//  double * efp_force;
//  for (int f=0; f<efp_fragments.size(); ++f) {
//    efp_force = efp_fragments[f]->get_forces_pointer();
//    for (int i=0; i<efp_fragments[f]->g_nintco(); ++i)
//      f_q[ g_efp_fragment_intco_offset(f) + i ] =  efp_force[i] ;
//  }
#endif

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"Internal forces in au\n");
    print_matrix(outfile, &f_q, 1, g_nintco()); fflush(outfile);
  }

/*
    // test by transforming f_q back to cartesian forces and compare
  if (Opt_params.print_lvl > 3) {
    B = compute_B();
    temp_arr = init_array(Ncart);
    opt_matrix_mult(B, 1, &f_q, 1, &temp_arr, 1, Ncart, Nintco, 1, 0);
    fprintf(outfile,"Recomputed forces in cartesian coordinates\n");
    print_matrix(outfile, &temp_arr, 1, Ncart);
    free_array(temp_arr);
    free_matrix(B);
  }
*/

  fflush(outfile);
  return ;
}

// project redundancies (and constraints) out of forces and Hessian matrix
// add constraints here later
void MOLECULE::project_f_and_H(void) {
  int Nintco = g_nintco();
  int Ncart = 3*g_natom();

  // compute G = B B^t
  double **G = init_matrix(Nintco, Nintco);
  double **B = compute_B();
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(B);

#if defined (OPTKING_PACKAGE_QCHEM)
//  // Put 1's on diagonal for EFP coordinates
//  for (int I=0; I<efp_fragments.size(); ++I)
//    for (int i=0; i<efp_fragments[I]->g_nintco(); ++i)
//      G[g_efp_fragment_intco_offset(I) + i][g_efp_fragment_intco_offset(I) + i] = 1.0;
#endif

  // compute P = G G^-1
  double **G_inv = symm_matrix_inv(G, Nintco, 1);
  double **P = init_matrix(Nintco, Nintco);
  opt_matrix_mult(G, 0, G_inv, 0, P, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(G);
  free_matrix(G_inv);

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\tProjection matrix for redundancies.\n");
    print_matrix(outfile, P, Nintco, Nintco); fflush(outfile);
  }

  // add constraints to projection matrix
  double **C = compute_constraints();

  // P = P' - P' C (CPC)^-1 C P'
  double **T = init_matrix(Nintco,Nintco);
  double **T2 = init_matrix(Nintco,Nintco);
  opt_matrix_mult(P, 0, C, 0,  T, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult(C, 0, T, 0, T2, 0, Nintco, Nintco, Nintco, 0);
  double **T3 = symm_matrix_inv(T2, Nintco, 1);

  opt_matrix_mult( C, 0,  P, 0,  T, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult(T3, 0,  T, 0, T2, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult( C, 0, T2, 0, T3, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult( P, 0, T3, 0, T2, 0, Nintco, Nintco, Nintco, 0);
  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<Nintco; ++j)
      P[i][j] -= T2[i][j];
  free_matrix(T);
  free_matrix(T2);
  free_matrix(T3);
  free_matrix(C);

  // Project redundancies and contraints out of forces
  // Unconstrained forces will not need this unless we have interpolated forces,
  // or some other as yet unimplemented method.
  double *f_q = p_Opt_data->g_forces_pointer();
  // f_q~ = P f_q
  double * temp_arr = init_array(Nintco);
  opt_matrix_mult(P, 0, &f_q, 1, &temp_arr, 1, Nintco, Nintco, 1, 0);
  array_copy(temp_arr, f_q, g_nintco());
  free_array(temp_arr);

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile,"\tInternal forces in au, after projection of redundancies and constraints.\n");
    if (Opt_params.efp_fragments)
      fprintf(outfile,"\tEFP external coordinates are not projected.\n");
    print_matrix(outfile, &f_q, 1, g_nintco());
  }

  // Project redundances and constraints out of Hessian matrix
  // Peng, Ayala, Schlegel, JCC 1996 give H -> PHP + 1000(1-P)
  // The second term appears unnecessary and sometimes messes up Hessian updating.

  double **H = p_Opt_data->g_H_pointer();
  double **temp_mat = init_matrix(Nintco, Nintco);
  opt_matrix_mult(H, 0, P, 0, temp_mat, 0, Nintco, Nintco, Nintco, 0);
  opt_matrix_mult(P, 0, temp_mat, 0, H, 0, Nintco, Nintco, Nintco, 0);
  free_matrix(temp_mat);

/* 2nd term unnecessary
  for (int i=0; i<Nintco;++i)
    H[i][i] += 1000 * (1.0 - P[i][i]);

  double tval;
  for (int i=0; i<Nintco; ++i)
    for (int j=0; j<i; ++j) {
      tval = H[i][j] - 1000 * P[i][j];
      H[j][i] = H[i][j] = tval;
    } */

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Projected (PHP) Hessian matrix\n");
//    if (Opt_params.efp_fragments)
//      fprintf(outfile,"EFP external coordinates are not projected.\n");
    print_matrix(outfile, H, g_nintco(), g_nintco());
  }
  free_matrix(P);
  fflush(outfile);
}

void MOLECULE::print_geom(void) {
  fprintf(outfile,"\tCartesian Geometry (au)\n");
  fflush(outfile);
  for (int i=0; i<fragments.size(); ++i)
    fragments[i]->print_geom(outfile);
}

void MOLECULE::apply_intrafragment_step_limit(double * & dq) {
  int i, f;
  double scale = 1.0;
  double limit = Opt_params.intrafragment_step_limit;

  for (f=0; f<fragments.size(); ++f)
    for (i=0; i<fragments[f]->g_nintco(); ++i)
      if (scale * fabs(dq[g_intco_offset(f)+i]) > limit)
        scale = limit / fabs(dq[g_intco_offset(f)+i]);

  if (scale != 1.0) {
    fprintf(outfile,"\tChange in coordinate exceeds step limit of %10.5lf.\n", limit);
    fprintf(outfile,"\tScaling displacements by %10.5lf\n", scale);

    for (f=0; f<fragments.size(); ++f)
      for (i=0; i<fragments[f]->g_nintco(); ++i)
        dq[g_intco_offset(f)+i] *= scale;
  }
  fflush(outfile);
}

// don't let any angles get smaller than 0.0
void MOLECULE::check_intrafragment_zero_angles(double const * const dq) {
  for (int f=0; f<fragments.size(); ++f)
    fragments[f]->check_zero_angles(&(dq[g_intco_offset(f)]));
}

void MOLECULE::H_guess(void) const {
  double **H = p_Opt_data->g_H_pointer();

  if (Opt_params.intrafragment_H == OPT_PARAMS::SCHLEGEL)
    fprintf(outfile,"\tGenerating empirical Hessian (Schlegel '84) for each fragment.\n");
  else if (Opt_params.intrafragment_H == OPT_PARAMS::FISCHER)
    fprintf(outfile,"\tGenerating empirical Hessian (Fischer & Almlof '92) for each fragment.\n");

  for (int f=0; f<fragments.size(); ++f) {
    double **H_frag = fragments[f]->H_guess();

    for (int i=0; i<fragments[f]->g_nintco(); ++i)
      for (int j=0; j<fragments[f]->g_nintco(); ++j)
        H[g_intco_offset(f) + i][g_intco_offset(f) + j] = H_frag[i][j];

    free_matrix(H_frag);
  }

  for (int I=0; I<interfragments.size(); ++I) {
    double **H_interfrag = interfragments[I]->H_guess();

    for (int i=0; i<interfragments[I]->g_nintco(); ++i)
      for (int j=0; j<interfragments[I]->g_nintco(); ++j)
        H[g_interfragment_intco_offset(I) + i][g_interfragment_intco_offset(I) + j] =
          H_interfrag[i][j];

    free_matrix(H_interfrag);
  }

#if defined(OPTKING_PACKAGE_QCHEM)
//  for (int I=0; I<efp_fragments.size(); ++I) {
//    double **H_efp_frag = efp_fragments[I]->H_guess();
//
//    for (int i=0; i<efp_fragments[I]->g_nintco(); ++i)
//      for (int j=0; j<efp_fragments[I]->g_nintco(); ++j)
//        H[g_efp_fragment_intco_offset(I) + i][g_efp_fragment_intco_offset(I) + j] = H_efp_frag[i][j];
//
//    free_matrix(H_efp_frag);
//  }
#endif

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Initial Hessian guess\n");
    print_matrix(outfile, H, g_nintco(), g_nintco());
  }
  fflush(outfile);
  return;
}

bool MOLECULE::cartesian_H_to_internals(void) const {
  int Nintco = g_nintco();
  int Ncart = 3*g_natom();
  bool success = true; // to be dynamic later

  double **H_int = p_Opt_data->g_H_pointer();

  // compute A = u B^t (B u B^t)^-1 where u=unit matrix and -1 is generalized inverse
  double **B = compute_B();
  double **G = init_matrix(Nintco, Nintco);
  opt_matrix_mult(B, 0, B, 1, G, 0, Nintco, Ncart, Nintco, 0);

  double **G_inv = symm_matrix_inv(G, Nintco, true);
  free_matrix(G);

  double **A = init_matrix(Ncart, Nintco);
  opt_matrix_mult(B, 1, G_inv, 0, A, 0, Ncart, Nintco, Nintco, 0);
  free_matrix(G_inv);
  free_matrix(B);

  // compute gradient in internal coordinates, A^t g_x = g_q
  double *grad_x = g_grad_array();
  double *grad_q = init_array(Nintco);
  opt_matrix_mult(A, 1, &grad_x, 1, &grad_q, 1, Nintco, Ncart, 1, 0);
  free_array(grad_x);

  // read in cartesian H
  //double **H_cart = p_Opt_data->read_cartesian_H();
  double **H_cart;
   std::cout << "DES Error: Cant't read in hessian" << std::endl;
   exit(-1);

  // transform cartesian H to internals; A^t (H_x - K) A
  // K_ij = sum_q ( grad_q[q] d^2(q)/(dxi dxj) )
  double **dq2dx2;

  for (int q=0; q<Nintco; ++q) {
    dq2dx2 = compute_derivative_B(q); // d^2(q)/ dx_i dx_j

    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        H_cart[i][j] -= grad_q[q] * dq2dx2[i][j];

    free_matrix(dq2dx2);
  }
  free_array(grad_q);

  double **temp_mat = init_matrix(Ncart, Nintco);
  opt_matrix_mult(H_cart, 0, A, 0, temp_mat, 0, Ncart, Ncart, Nintco, 0);
  free_matrix(H_cart);

  //double **H_int = init_matrix(Nintco, Nintco);
  opt_matrix_mult(A, 1, temp_mat, 0, H_int, 0, Nintco, Ncart, Nintco, 0);
  free_matrix(temp_mat);

  free_matrix(A);

  if (Opt_params.print_lvl >= 3) {
    fprintf(outfile, "Hessian transformed to internal coordinates:\n");
    print_matrix(outfile, H_int, Nintco, Nintco); fflush(outfile);
  }

  // Check by transforming internal coordinate Hessian back into cartesian coordinates:
/*
  B = compute_B();
  temp_mat = init_matrix(Ncart, Nintco);
  opt_matrix_mult(B, 1, H_int, 0, temp_mat, 0, Ncart, Nintco, Nintco, 0);
  H_cart = init_matrix(Ncart, Ncart);
  opt_matrix_mult(temp_mat, 0, B, 0, H_cart, 0, Ncart, Nintco, Ncart, 0);
  free_matrix(temp_mat);

  for (int q=0; q<Nintco; ++q) {
    dq2dx2 = compute_derivative_B(q); // d^2(q)/ dx_i dx_j
    for (int i=0; i<Ncart; ++i)
      for (int j=0; j<Ncart; ++j)
        H_cart[i][j] += grad_q[q] * dq2dx2[i][j];
    free_matrix(dq2dx2);
  }
  free_array(grad_q);
  fprintf(outfile, "Hessian transformed back into Cartesian coordinates\n");
  print_matrix(outfile, H_cart, Ncart, Ncart); fflush(outfile);
  free_matrix(B);
*/
  return success;
}

// compute B matrix - leave rows for EFP coordinates empty
double ** MOLECULE::compute_B(void) const {
  double **B, **B_frag, **B_inter;
  int i, j;

  B = init_matrix(g_nintco(), 3*g_natom());

  for (int f=0; f<fragments.size(); ++f) {
    B_frag = fragments[f]->compute_B();

    for (i=0; i<fragments[f]->g_nintco(); ++i) {
      for (j=0; j<fragments[f]->g_natom(); ++j) {
         for(int k=0;k<3;k++) {
            B[g_intco_offset(f)+i][3*fragments[f]->globalIndex[j]+k] = B_frag[i][3*j+k];
         }
      }
   }

    free_matrix(B_frag);
  }

  for (int I=0; I<interfragments.size(); ++I) {
    B_inter = interfragments[I]->compute_B(); // ->g_nintco() X (3*atom A)+3(natom_B)

    int iA = interfragments[I]->g_A_index();
    int iB = interfragments[I]->g_B_index();
    int nA = interfragments[I]->g_natom_A();
    int nB = interfragments[I]->g_natom_B();

    for (i=0; i<interfragments[I]->g_nintco(); ++i) { // for each of up to 6 coordinates

      for (j=0; j<3*nA; ++j)
        B[g_interfragment_intco_offset(I)+i][3*g_atom_offset(iA)+j] = B_inter[i][j];

      for (j=0; j<3*nB; ++j)
        B[g_interfragment_intco_offset(I)+i][3*g_atom_offset(iB)+j] = B_inter[i][3*nA+j];

    }
    free_matrix(B_inter);
  }
  return B;
}

double ** MOLECULE::compute_derivative_B(int intco_index) const {
  int cnt_intcos = 0;
  int fragment_index = -1;
  int coordinate_index = 0;
  bool is_interfragment = true;
  int f;

  for (f=0; f<fragments.size(); ++f) {
    for (int i=0; i<fragments[f]->g_nintco(); ++i) {
      if (cnt_intcos++ == intco_index) {
        fragment_index = f;
        coordinate_index = i;
        is_interfragment = false;
        break;
      }
    }
  }

  if (is_interfragment) {  // intco_index not yet found
    for (f=0; f<interfragments.size(); ++f) {
      for (int i=0; i<interfragments[f]->g_nintco(); ++i) {
        if (cnt_intcos++ == intco_index) {
          fragment_index = f;
          coordinate_index = i;
          break;
        }
      }
    }
  }

  if (fragment_index == -1)
    throw(INTCO_EXCEPT("MOLECULE::compute_derivative_B() could not find intco_index"));

  double **dq2dx2_frag;

  double **dq2dx2 = init_matrix(3*g_natom(), 3*g_natom());

  if (!is_interfragment) {
    dq2dx2_frag = fragments[fragment_index]->compute_derivative_B(coordinate_index);
    // dimension is 2x2, 3x3 or 4x4 and given by natom_intco
    int natom_intco = fragments[fragment_index]->g_intco_natom(coordinate_index);
    int atom_a, atom_b;

    for (int a=0; a<natom_intco; ++a) {
      atom_a = g_atom_offset(fragment_index) + fragments[fragment_index]->g_intco_atom(coordinate_index, a);
      for (int b=0; b<natom_intco; ++b) {
        atom_b = g_atom_offset(fragment_index) + fragments[fragment_index]->g_intco_atom(coordinate_index, b);

      for (int xyz_a=0; xyz_a<3; ++xyz_a)
        for (int xyz_b=0; xyz_b<3; ++xyz_b)
          dq2dx2[3*atom_a + xyz_a][3*atom_b + xyz_b] = dq2dx2_frag[3*a+xyz_a][3*b+xyz_b];

      }
    }
    free_matrix(dq2dx2_frag);
  }
  else { // interfragment cordinate
    dq2dx2_frag = interfragments[fragment_index]->compute_derivative_B(coordinate_index);

    // dimension is (3*natomA + 3*natomB) x (3*natomA + 3*natomB) -> 3*natom X 3*natom
    int nA = interfragments[fragment_index]->g_natom_A();
    int nB = interfragments[fragment_index]->g_natom_B();
    int iA = 3*g_atom_offset(interfragments[fragment_index]->g_A_index());
    int iB = 3*g_atom_offset(interfragments[fragment_index]->g_B_index());

    //print_matrix(outfile,dq2dx2_frag, 3*(nA+nB), 3*(nA+nB));
    for (int a=0; a<3*nA; ++a) // A-A block
      for (int aa=0; aa<3*nA; ++aa)
        dq2dx2[iA + a][iA + aa] = dq2dx2_frag[a][aa];

    for (int a=0; a<3*nA; ++a) // A-B block
      for (int bb=0; bb<3*nB; ++bb)
        dq2dx2[iA + a][iB + bb] = dq2dx2_frag[a][3*nA + bb];

    for (int b=0; b<3*nB; ++b) // B-A block
      for (int aa=0; aa<3*nA; ++aa)
        dq2dx2[iB + b][iA + aa] = dq2dx2_frag[3*nA + b][aa];

    for (int b=0; b<3*nB; ++b) // B-B block
      for (int bb=0; bb<3*nB; ++bb)
        dq2dx2[iB + b][iB + bb] = dq2dx2_frag[3*nA + b][3*nA + bb];

    free_matrix(dq2dx2_frag);
  }
  return dq2dx2;
}

// print internal coordinates to text file
void MOLECULE::print_intco_dat(FILE *fp_intco) {
  for (int i=0; i<fragments.size(); ++i) {
    int first = g_atom_offset(i);
    fprintf(fp_intco,"F %d %d\n", first+1, first + fragments[i]->g_natom()); fflush(outfile);
    fragments[i]->print_intco_dat(fp_intco, g_atom_offset(i));
  }

  for (int I=0; I<interfragments.size(); ++I) {
    int frag_a = interfragments[I]->g_A_index();
    int frag_b = interfragments[I]->g_B_index();
    fprintf(fp_intco,"I %d %d\n", frag_a+1, frag_b+1); fflush(outfile);

    for (int i=0; i<6; ++i) 
      fprintf(fp_intco," %d", (int) interfragments[I]->coordinate_on(i));
    fprintf(fp_intco,"\n"); fflush(outfile);

    interfragments[I]->print_intco_dat(fp_intco, g_atom_offset(frag_a), g_atom_offset(frag_b));
  }
}

// tell whether internal coordinate is frozen
double ** MOLECULE::compute_constraints(void) {
  double **C, **C_frag, **C_inter;
  int i, j;

  C = init_matrix(g_nintco(), g_nintco());

  for (int f=0; f<fragments.size(); ++f) {
    C_frag = fragments[f]->compute_constraints();

    for (i=0; i<fragments[f]->g_nintco(); ++i)
      for (j=0; j<fragments[f]->g_nintco(); ++j)
        C[g_intco_offset(f)+i][g_intco_offset(f)+j] = C_frag[i][j];

    free_matrix(C_frag);
  }

  for (int I=0; I<interfragments.size(); ++I) {
    C_inter = interfragments[I]->compute_constraints(); // g_nintco() X g_nintco

    for (i=0; i<interfragments[I]->g_nintco(); ++i)
      for (j=0; j<interfragments[I]->g_nintco(); ++j)
        C[g_interfragment_intco_offset(I)+i][g_interfragment_intco_offset(I)+j] = C_inter[i][j];

    free_matrix(C_inter);
  }

  if (Opt_params.print_lvl >= 2) {
    fprintf(outfile,"Constraint matrix\n");
    print_matrix(outfile, C, g_nintco(), g_nintco()); fflush(outfile);
  }
  fflush(outfile);
  return C;
}

// The values of the EFP coordinates will be taken to be the total change
// since the beginning of the optimization.  These values are determined
// from the data in opt_data after that data is read.
//void MOLECULE::update_efp_values(void) {
//  for (int i=0; i<efp_fragments.size(); ++i) {
//    double *vals = init_array(6);
//
//    double *dq;
//    for (int iter=0; iter<p_Opt_data->g_iteration(); ++iter) {
//      dq = p_Opt_data->g_dq_pointer(iter);
//      for (int coord=0; coord < 6; ++coord)
//        vals[coord] += dq[g_efp_fragment_intco_offset(i)+coord] ;
//    }
//    efp_fragments[i]->set_values(vals);
//    free_array(vals);
//  }
//}

}

