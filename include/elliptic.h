



//=========================== Olliptic elliptic.h =========================//
//  
//  Elliptic problem. (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//










#ifndef ELLIPTIC_H

#define ELLIPTIC_H



//============================ Standar libraries =========================//


#include <time.h>


//============================ Olliptic libraries =========================//

#include "ffunction.h"

#include "interface.h"

#include "multigrid.h"




using namespace std;


#define T_ETA 1
#define T_PHI 1

#define BYT1 0.1
#define BYT2 0.5
#define BYT3 4.0


class elliptic
{


  //== Grid ffunctions ==//



  multigrid vars;


  vector < domain > D;



  vector < string > var_names;

  vector < size_t > var_index;


  vector < string > coeff_names;

  vector < size_t > coeff_index;

  
  vector < string > source_names;
    
  vector < string > aux_names;

  vector < string > total_names;



  ofstream stdout_file;

  string stdout_file_name;


  ofstream norm_file;

  string norm_file_name;

  size_t ps_case;

  size_t my_rank;

  size_t total_nodes;


  size_t levels;

  size_t mgrid_levels;

  size_t num_dom;

  size_t boxes;


  t_output output;

  size_t order;



  string dir_name;



  string tov_filename;


  t_method method;

  t_equation equation;

  t_boundary boundary;

  t_domain dom;


  bool print_1d;

  bool print_2d;

  bool print_3d;

  bool print_1dt;

  bool print_interpol;

  bool interpol_output;

  bool fit_approx;

  vector<double> u_at_p;

  double alpha_GB=0;

  double eps_coarse;

  double eps_solv;


  size_t nu1;

  size_t nu2;


  double ps_d1, ps_r0, ps_isigma, ps_phi0;

  double ps_e1, ps_e2, ps_de1;

  size_t it;

  size_t git;



  vector < size_t > falloff_n;

  vector < double >boundary_coeff;


  double PI;


  size_t NumBH;

    vector < double >bhmp;

    vector < double >bhqp;

    vector < double >ADM_mass;

    vector < double >bhx;

    vector < double >bhy;

    vector < double >bhz;


   vector < double >bhpx;

    vector < double >bhpy;

    vector < double >bhpz;


    vector < double >bhsx;

    vector < double >bhsy;

    vector < double >bhsz;






  void App_Ell_Operator (size_t l, size_t i);

  double Smooth (size_t l, size_t i, size_t nu = 1);

  void App_Boundary_Operator (size_t l, size_t i);



  double FASN (size_t k, size_t k_top, bool Wcycle = true, size_t cycle = 0);



  void Set_Method ();

  void Set_Grids (size_t nx, size_t ny, size_t nz, double dx, double dy,
		  double dz);

  void Set_Equation ();

  void Set_Coefficients ();

  void Print_Coeff ();

  void Print_Pre_Cycle (size_t k, size_t ktop, double norm_k, size_t cycle);

  void Print_Post_Cycle (size_t k, size_t ktop, double norm_k, size_t cycle);

  void Print_Cycle (size_t k, size_t ktop, double norm_k, size_t cycle);

 

  void (elliptic::*solve_ell) ();


  typedef double (elliptic::*pt2funct) (double x, double y, double z);

    vector < pt2funct > Boundary;


  typedef double (elliptic::*pt2rhs) (size_t l, size_t ibox);
    vector < pt2rhs > Operator;

    vector < pt2rhs > duOperator;



  void Info_Punctures ();

  void Compute_ADM_Mass();

  void Compute_Solution();


    public:~elliptic ()
  {
  };

  elliptic (interface * oll)
  {
    Make (oll);
  };




  void Make (interface * oll);





  void Set_Boundary ();





  void Print_Info ();

  void Print_Solution (double pd = 0.1,
		       double Lx = 1.0, double Ly = 1.0, double Lz = 1.0);



  void IO_bam (string file);

  void IO_Zcode (string file);




  void App_Bound (size_t l, size_t ibox);

  void App_Symmetry (ffunction * U);




  void Gauss_Seidel ();

  void ApproximateID ();

  void Vcycle ();

  void Wcycle ();

  void FullMG ()
  {
  };


  void Get_Puncture_App();


  inline void Solve ()
  {
    return ((this->*solve_ell) ());
  }


  inline double OPERATOR (size_t lvar, size_t l, size_t ibox)
  {
    return ((this->*Operator[lvar]) (l, ibox));
  }

  inline double duOPERATOR (size_t lvar, size_t l, size_t ibox)
  {
    return ((this->*duOperator[lvar]) (l, ibox));
  }






//======================== elliptic operator =======================//





//======================== Poisson =======================//

  inline double poisson_OP (size_t l, size_t ibox)
  {

    return (vars.Lap (l, ibox, var_index[0]) + vars.Get_val(l, ibox, var_index[0])*vars.Get_val(l, ibox, var_index[0]));

  }

  inline double poisson_duOP (size_t l, size_t ibox)
  {


    return (vars.duLap (l, ibox, var_index[0])+2.0*vars.Get_val(l, ibox, var_index[0]));

  }




//======================== Brill =======================//


  inline double brill_OP (size_t l, size_t ibox)
  {

    return (vars.Lap (l, ibox, var_index[0]) +
	    vars.Get_val (l, ibox, coeff_index[0]) * vars.Get_val (l, ibox,
								   var_index
								   [0]));

  }


  inline double brill_duOP (size_t l, size_t ibox)
  {
    return (vars.duLap (l, ibox, var_index[0]) +
	    vars.Get_val (l, ibox, coeff_index[0]));

  }



//======================== Punctures =======================//



  inline double puncture_OP (size_t l, size_t ibox)
  {



    return (vars.Lap (l, ibox, var_index[0]) +  vars.Get_val (l, ibox, coeff_index[0]) / pow (vars.Get_val (l, ibox,coeff_index[1])+ vars.Get_val (l, ibox,var_index[0]), 7));

  }


  inline double puncture_duOP (size_t l, size_t ibox)
  {

    return (vars.duLap (l, ibox, var_index[0]) - 7.0 * vars.Get_val (l, ibox, coeff_index[0]) / pow (vars.Get_val (l, ibox,coeff_index[1]) +vars.Get_val (l, ibox,var_index[0]),8));


  }



//=================== Punctures + scalar field ====================//



  inline double puncture_SF_OP (size_t l, size_t ibox)
  {



    return vars.Lap(l, ibox, var_index[0]) +  vars.Get_val (l, ibox, coeff_index[0]) / pow (vars.Get_val (l, ibox,coeff_index[1])+ vars.Get_val (l, ibox,var_index[0]), 7) + vars.Get_val (l, ibox, coeff_index[2])*(vars.Get_val(l, ibox,coeff_index[1])+ vars.Get_val (l, ibox,var_index[0]))+ vars.Get_val (l, ibox, coeff_index[3])*pow(vars.Get_val (l, ibox,coeff_index[1])+ vars.Get_val (l, ibox,var_index[0]),5);

  }


  inline double puncture_SF_duOP (size_t l, size_t ibox)
  {

    return vars.duLap (l, ibox, var_index[0]) - 7.0 * vars.Get_val (l, ibox, coeff_index[0]) / pow (vars.Get_val (l, ibox,coeff_index[1]) +vars.Get_val (l, ibox,var_index[0]),8) + vars.Get_val (l, ibox, coeff_index[2])+ 5*vars.Get_val (l, ibox, coeff_index[3])*pow(vars.Get_val (l, ibox,coeff_index[1])+ vars.Get_val (l, ibox,var_index[0]),4);


  }



//=================== Punctures in Gauss-Bonnet ====================//



  inline double puncture_GB_OP (size_t l, size_t ibox)
  {

    double psi = vars.Get_val (l, ibox,var_index[0]);
    double Lap_psi=vars.Lap(l, ibox, var_index[0]);
    
    double psi_dx = vars.Dx (l, ibox,var_index[0]);
    double psi_dy = vars.Dy (l, ibox,var_index[0]);
    double psi_dz = vars.Dz (l, ibox,var_index[0]);

    double psi_dxdx = vars.DDx (l, ibox,var_index[0]);
    double psi_dydx = vars.DDxy (l, ibox,var_index[0]);
    double psi_dzdx = vars.DDxz (l, ibox,var_index[0]);

    double psi_dxdy = vars.DDxy (l, ibox,var_index[0]);
    double psi_dydy = vars.DDy (l, ibox,var_index[0]);
    double psi_dzdy = vars.DDyz (l, ibox,var_index[0]);

    double psi_dxdz = vars.DDxz (l, ibox,var_index[0]);
    double psi_dydz = vars.DDyz (l, ibox,var_index[0]);
    double psi_dzdz = vars.DDz (l, ibox,var_index[0]);
    
    double phi    = vars.Get_val (l, ibox,coeff_index[0]);

    double phi_dx = vars.Get_val (l, ibox,coeff_index[1]);
    double phi_dy = vars.Get_val (l, ibox,coeff_index[2]);
    double phi_dz = vars.Get_val (l, ibox,coeff_index[3]);

    double phi_dxdx = vars.Get_val (l, ibox,coeff_index[4]);
    double phi_dydx = vars.Get_val (l, ibox,coeff_index[5]);
    double phi_dzdx = vars.Get_val (l, ibox,coeff_index[6]);

    double phi_dxdy = vars.Get_val (l, ibox,coeff_index[7]);
    double phi_dydy = vars.Get_val (l, ibox,coeff_index[8]);
    double phi_dzdy = vars.Get_val (l, ibox,coeff_index[9]);

    double phi_dxdz = vars.Get_val (l, ibox,coeff_index[10]);
    double phi_dydz = vars.Get_val (l, ibox,coeff_index[11]);
    double phi_dzdz = vars.Get_val (l, ibox,coeff_index[12]);

    
    return Lap_psi -4*alpha_GB*pow(phi_dx, 2)*psi_dxdx*exp(-phi)/pow(psi, 4) - 2*alpha_GB*pow(phi_dx, 2)*psi_dydy*exp(-phi)/pow(psi, 4) - 2*alpha_GB*pow(phi_dx, 2)*psi_dzdz*exp(-phi)/pow(psi, 4) - 4*alpha_GB*pow(phi_dx, 2)*pow(psi_dx, 2)*exp(-phi)/pow(psi, 5) + 2*alpha_GB*pow(phi_dx, 2)*pow(psi_dy, 2)*exp(-phi)/pow(psi, 5) + 2*alpha_GB*pow(phi_dx, 2)*pow(psi_dz, 2)*exp(-phi)/pow(psi, 5) - 4*alpha_GB*phi_dx*phi_dy*psi_dxdy*exp(-phi)/pow(psi, 4) - 12*alpha_GB*phi_dx*phi_dy*psi_dx*psi_dy*exp(-phi)/pow(psi, 5) - 4*alpha_GB*phi_dx*phi_dz*psi_dxdz*exp(-phi)/pow(psi, 4) - 12*alpha_GB*phi_dx*phi_dz*psi_dx*psi_dz*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dx*psi_dx*psi_dxdx*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dx*psi_dxdy*psi_dy*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dx*psi_dxdz*psi_dz*exp(-phi)/pow(psi, 5) - 16*alpha_GB*phi_dx*pow(psi_dx, 3)*exp(-phi)/pow(psi, 6) - 16*alpha_GB*phi_dx*psi_dx*pow(psi_dy, 2)*exp(-phi)/pow(psi, 6) - 16*alpha_GB*phi_dx*psi_dx*pow(psi_dz, 2)*exp(-phi)/pow(psi, 6) + 2*alpha_GB*phi_dxdx*psi_dydy*exp(-phi)/pow(psi, 4) + 2*alpha_GB*phi_dxdx*psi_dzdz*exp(-phi)/pow(psi, 4) - 2*alpha_GB*phi_dxdx*pow(psi_dx, 2)*exp(-phi)/pow(psi, 5) - 2*alpha_GB*phi_dxdx*pow(psi_dy, 2)*exp(-phi)/pow(psi, 5) - 2*alpha_GB*phi_dxdx*pow(psi_dz, 2)*exp(-phi)/pow(psi, 5) - 4*alpha_GB*phi_dxdy*psi_dxdy*exp(-phi)/pow(psi, 4) - 4*alpha_GB*phi_dxdz*psi_dxdz*exp(-phi)/pow(psi, 4) - 2*alpha_GB*pow(phi_dy, 2)*psi_dxdx*exp(-phi)/pow(psi, 4) - 4*alpha_GB*pow(phi_dy, 2)*psi_dydy*exp(-phi)/pow(psi, 4) - 2*alpha_GB*pow(phi_dy, 2)*psi_dzdz*exp(-phi)/pow(psi, 4) + 2*alpha_GB*pow(phi_dy, 2)*pow(psi_dx, 2)*exp(-phi)/pow(psi, 5) - 4*alpha_GB*pow(phi_dy, 2)*pow(psi_dy, 2)*exp(-phi)/pow(psi, 5) + 2*alpha_GB*pow(phi_dy, 2)*pow(psi_dz, 2)*exp(-phi)/pow(psi, 5) - 4*alpha_GB*phi_dy*phi_dz*psi_dydz*exp(-phi)/pow(psi, 4) - 12*alpha_GB*phi_dy*phi_dz*psi_dy*psi_dz*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dy*psi_dx*psi_dxdy*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dy*psi_dy*psi_dydy*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dy*psi_dydz*psi_dz*exp(-phi)/pow(psi, 5) - 16*alpha_GB*phi_dy*pow(psi_dx, 2)*psi_dy*exp(-phi)/pow(psi, 6) - 16*alpha_GB*phi_dy*pow(psi_dy, 3)*exp(-phi)/pow(psi, 6) - 16*alpha_GB*phi_dy*psi_dy*pow(psi_dz, 2)*exp(-phi)/pow(psi, 6) + 2*alpha_GB*phi_dydy*psi_dxdx*exp(-phi)/pow(psi, 4) + 2*alpha_GB*phi_dydy*psi_dzdz*exp(-phi)/pow(psi, 4) - 2*alpha_GB*phi_dydy*pow(psi_dx, 2)*exp(-phi)/pow(psi, 5) - 2*alpha_GB*phi_dydy*pow(psi_dy, 2)*exp(-phi)/pow(psi, 5) - 2*alpha_GB*phi_dydy*pow(psi_dz, 2)*exp(-phi)/pow(psi, 5) - 4*alpha_GB*phi_dydz*psi_dydz*exp(-phi)/pow(psi, 4) - 2*alpha_GB*pow(phi_dz, 2)*psi_dxdx*exp(-phi)/pow(psi, 4) - 2*alpha_GB*pow(phi_dz, 2)*psi_dydy*exp(-phi)/pow(psi, 4) - 4*alpha_GB*pow(phi_dz, 2)*psi_dzdz*exp(-phi)/pow(psi, 4) + 2*alpha_GB*pow(phi_dz, 2)*pow(psi_dx, 2)*exp(-phi)/pow(psi, 5) + 2*alpha_GB*pow(phi_dz, 2)*pow(psi_dy, 2)*exp(-phi)/pow(psi, 5) - 4*alpha_GB*pow(phi_dz, 2)*pow(psi_dz, 2)*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dz*psi_dx*psi_dxdz*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dz*psi_dy*psi_dydz*exp(-phi)/pow(psi, 5) - 8*alpha_GB*phi_dz*psi_dz*psi_dzdz*exp(-phi)/pow(psi, 5) - 16*alpha_GB*phi_dz*pow(psi_dx, 2)*psi_dz*exp(-phi)/pow(psi, 6) - 16*alpha_GB*phi_dz*pow(psi_dy, 2)*psi_dz*exp(-phi)/pow(psi, 6) - 16*alpha_GB*phi_dz*pow(psi_dz, 3)*exp(-phi)/pow(psi, 6) + 2*alpha_GB*phi_dzdz*psi_dxdx*exp(-phi)/pow(psi, 4) + 2*alpha_GB*phi_dzdz*psi_dydy*exp(-phi)/pow(psi, 4) - 2*alpha_GB*phi_dzdz*pow(psi_dx, 2)*exp(-phi)/pow(psi, 5) - 2*alpha_GB*phi_dzdz*pow(psi_dy, 2)*exp(-phi)/pow(psi, 5) - 2*alpha_GB*phi_dzdz*pow(psi_dz, 2)*exp(-phi)/pow(psi, 5) - 6*alpha_GB*pow(psi_dx, 2)*psi_dxdx*exp(-phi)/pow(psi, 5) - 12*alpha_GB*psi_dx*psi_dxdy*psi_dy*exp(-phi)/pow(psi, 5) - 12*alpha_GB*psi_dx*psi_dxdz*psi_dz*exp(-phi)/pow(psi, 5) - 6*alpha_GB*pow(psi_dy, 2)*psi_dydy*exp(-phi)/pow(psi, 5) - 12*alpha_GB*psi_dy*psi_dydz*psi_dz*exp(-phi)/pow(psi, 5) - 6*alpha_GB*pow(psi_dz, 2)*psi_dzdz*exp(-phi)/pow(psi, 5) + (1.0L/16.0L)*pow(phi_dx, 2)*psi + (1.0L/16.0L)*pow(phi_dy, 2)*psi + (1.0L/16.0L)*pow(phi_dz, 2)*psi;

  }


  inline double puncture_GB_duOP (size_t l, size_t ibox)
  {

    return 0;


  }

//===================== Electromagnetic BH =========================//


  inline double puncture_EM_OP (size_t l, size_t ibox)
  {
    //index ---> var 
    // 0  Ap
    // 1  Am
    // 2  B
    // 3 C
    // 4 D
    // 5 E
    // 6 B1z 
    // 7 B2x
    // 8 B2y
    // 9 B2z 
    // 10 C
    // 11 D 
    // 12 E 
   

    //Lap U 
    //+ (B (dx(U)^2+dy(U)^2+dz(U)^2)
    //+ (B1x U + B2x) Dx(U) 
    //+ (B1y U + B2y) Dy(U) 
    //+ (B1z U + B2z) Dz(U) 
    //+ C U^2 + D U + E)/(U^3 + A1 U^2 + A2 U + A3) 
    double uu = vars.Get_val (l, ibox,var_index[0]);
    double uu2 = uu*uu;

    double udx = vars.Dx (l, ibox,var_index[0]);
    double udy = vars.Dy (l, ibox,var_index[0]);
    double udz = vars.Dz (l, ibox,var_index[0]);

    return vars.Lap(l, ibox, var_index[0]) 
      + (vars.Get_val (l, ibox, coeff_index[3]) * ( udx*udx + udy*udy + udz*udz )
	 + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * udx 
	 + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * udy 
	 + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * udz 
	 + vars.Get_val (l, ibox, coeff_index[10]) * uu2 
	 + vars.Get_val (l, ibox, coeff_index[11]) * uu 
	 + vars.Get_val (l, ibox, coeff_index[12]))/(uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]) ); 



  }


  inline double puncture_EM_duOP (size_t l, size_t ibox)
  {

    //index ---> var 
    // 0  A1
    // 1  A2
    // 2  A3
    // 3 B
    // 4 B1x
    // 5 B1y
    // 6 B1z 
    // 7 B2x
    // 8 B2y
    // 9 B2z 
    // 10 C
    // 11 D 
    // 12 E


    //(U^3 + A1 U^2 + A2 U + A3) duLap U 
    //+ (3*U^2 + 2*A1 U + A2) Lap U 
    //+ 2*B (dx(U)*duDx(U) + dy(U)*duDy(U) + dz(U)*duDz(U))
    //+ (B1x U + B2x) duDx(U) 
    //+ B1x Dx(U) 
    //+ (B1y U + B2y) duDy(U) 
    //+ B1y Dy(U) 
    //+ (B1z U + B2z) duDz(U) 
    //+ B1z Dz(U) 
    //+ 2*C U + D 
    double uu = vars.Get_val (l, ibox,var_index[0]);
    double uu2 = uu*uu;

    double udx = vars.Dx (l, ibox,var_index[0]);
    double udy = vars.Dy (l, ibox,var_index[0]);
    double udz = vars.Dz (l, ibox,var_index[0]);

    double dudx = vars.duDx (l, ibox,var_index[0]);
    double dudy = vars.duDy (l, ibox,var_index[0]);
    double dudz = vars.duDz (l, ibox,var_index[0]);

    /*
    return vars.duLap(l, ibox, var_index[0]) 
      + ( 2*vars.Get_val (l, ibox, coeff_index[3]) * ( dudx*udx + dudy*udy + dudz*udz )
	 + (vars.Get_val (l, ibox, coeff_index[4]) ) * udx 
	 + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * dudx 
	 + (vars.Get_val (l, ibox, coeff_index[5]) ) * udy 
	 + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * dudy 
	 + (vars.Get_val (l, ibox, coeff_index[6]) ) * udz 
	 + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * dudz 
	 + 2 * vars.Get_val (l, ibox, coeff_index[10]) * uu 
	 + vars.Get_val (l, ibox, coeff_index[11]) );
    */
    
    return vars.duLap(l, ibox, var_index[0]) 
      +((2*vars.Get_val (l, ibox, coeff_index[3]) * ( udx*dudx + udy*dudy + udz*dudz )
	 + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * vars.duDx (l, ibox,var_index[0]) 
	 + vars.Get_val (l, ibox, coeff_index[4]) * vars.Dx (l, ibox,var_index[0]) 
	 + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * vars.duDy (l, ibox,var_index[0]) 
	 + vars.Get_val (l, ibox, coeff_index[5]) * vars.Dy (l, ibox,var_index[0]) 
	 + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * vars.duDz (l, ibox,var_index[0]) 
	 + vars.Get_val (l, ibox, coeff_index[6]) * vars.Dz (l, ibox,var_index[0]) 
	 + 2*vars.Get_val (l, ibox, coeff_index[10]) * uu 
	 + vars.Get_val (l, ibox, coeff_index[11]) )*(uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]) ) 
	- (vars.Get_val (l, ibox, coeff_index[3]) * ( udx*udx + udy*udy + udz*udz )
	   + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * udx 
	   + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * udy 
	   + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * udz 
	   + vars.Get_val (l, ibox, coeff_index[10]) * uu2 
	   + vars.Get_val (l, ibox, coeff_index[11]) * uu 
	   + vars.Get_val (l, ibox, coeff_index[12]))*(3*uu2 + 2*vars.Get_val(l, ibox, coeff_index[0])*uu + vars.Get_val(l, ibox, coeff_index[1])) )/pow(uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]),2); 

    

  }




//================== Neutron Star Oscillations ID =======================//



  inline double NSO_OP (size_t l, size_t ibox)
  {


    return (vars.Lap (l, ibox, var_index[0]) +
	    vars.Get_val (l, ibox,coeff_index[0])*pow (vars.Get_val(l, ibox, var_index[0]), -3));

  }


  inline double NSO_duOP (size_t l, size_t ibox)
  {

    return (vars.duLap (l, ibox, var_index[0]) - 
	    3.0 * vars.Get_val (l, ibox, coeff_index[0]) * pow (vars.Get_val (l, ibox, var_index[0]),-4));


  }



//======================== Elliptic system =======================//


  inline double system_OP_u (size_t l, size_t ibox)
  {
    return (vars.Lap (l, ibox, var_index[0]) +
	    vars.Get_val (l, ibox, coeff_index[0]) * vars.Get_val (l, ibox,
								   var_index
								   [1]) *
	    pow (vars.Get_val (l, ibox, var_index[0]), 2));
  }


  inline double system_duOP_u (size_t l, size_t ibox)
  {
    return (vars.duLap (l, ibox, var_index[0]) +
	    2.0 * vars.Get_val (l, ibox, coeff_index[0]) * vars.Get_val (l,
									 ibox,
									 var_index
									 [1])
	    * vars.Get_val (l, ibox, var_index[0]));
  }




  inline double system_OP_v (size_t l, size_t ibox)
  {
    return (vars.Lap (l, ibox, var_index[1]) +
	    vars.Get_val (l, ibox,
			  coeff_index[1]) * pow (vars.Get_val (l, ibox,
							       var_index[0]),
						 2) * vars.Get_val (l, ibox,
								    var_index
								    [1]));
  }


  inline double system_duOP_v (size_t l, size_t ibox)
  {
    return (vars.duLap (l, ibox, var_index[1]) +
	    vars.Get_val (l, ibox,
			  coeff_index[1]) * pow (vars.Get_val (l, ibox,
							       var_index[0]),
						 2));
  }




//======================== Trumpet =======================//

  inline double trumpet_OP (size_t l, size_t ibox)
  {




    return (vars.Lap (l, ibox, var_index[0]) -
	    1.25 * (pow (vars.Dx (l, ibox, var_index[0]), 2) +
		    pow (vars.Dy (l, ibox, var_index[0]),
			 2) + pow (vars.Dz (l, ibox, var_index[0]),
				   2)) / vars.Get_val (l, ibox,
						       var_index[0]));





  }

  inline double trumpet_duOP (size_t l, size_t ibox)
  {

    return (vars.duLap (l, ibox, var_index[0]) -
	    1.25 * (2.0 * vars.Get_val (l, ibox, var_index[0]) *
		    (vars.Dx (l, ibox, var_index[0]) *
		     vars.duDx (l, ibox, var_index[0]) + vars.Dy (l, ibox,
								  var_index
								  [0]) *
		     vars.duDy (l, ibox, var_index[0]) + vars.Dz (l, ibox,
								  var_index
								  [0]) *
		     vars.duDz (l, ibox,
				var_index[0])) - (pow (vars.Dx (l, ibox,
								var_index[0]),
						       2) + pow (vars.Dy (l,
									  ibox,
									  var_index
									  [0]),
								 2) +
						  pow (vars.
						       Dz (l, ibox,
							   var_index[0]),
						       2))) /
	    pow (vars.Get_val (l, ibox, var_index[0]), 2));






  }






//======================== Bowen-York =======================//




  inline double BowenYork_OP_X (size_t l, size_t ibox)
  {return vars.Lap (l, ibox, var_index[0]);}

  inline double BowenYork_duOP_X (size_t l, size_t ibox)
  { return vars.duLap (l, ibox, var_index[0]);}




  inline double BowenYork_OP_Y (size_t l, size_t ibox)
  {return vars.Lap (l, ibox, var_index[1]) ;}

  inline double BowenYork_duOP_Y (size_t l, size_t ibox)
  {return vars.duLap (l, ibox, var_index[1]) ;}



  inline double BowenYork_OP_Z (size_t l, size_t ibox)
  { return vars.Lap (l, ibox, var_index[2]);}

  inline double BowenYork_duOP_Z (size_t l, size_t ibox)
  {return vars.duLap (l, ibox, var_index[2]);}


  inline double BowenYork_OP_L (size_t l, size_t ibox)
  {return (vars.Lap (l, ibox, var_index[3]) -(
	   vars.Dx(l, ibox, var_index[0]) +
	   vars.Dy(l, ibox, var_index[1]) +
	   vars.Dz(l, ibox, var_index[2])) );}

  inline double BowenYork_duOP_L (size_t l, size_t ibox)
  {return vars.duLap (l, ibox, var_index[3]);}


//======================== BY-EM Operator =======================//


  //index ---> var 
  // 0 A1
  // 1 A2
  // 2 A3
  // 3 B
  // 4 B1x
  // 5 B1y
  // 6 B1z 
  // 7 B2x
  // 8 B2y
  // 9 B2z 
  // 10 C
  // 11 D 
  // 12 E
  // 13 BiBi
  //14 Jx
  //15 Jy
  //16 Jz
  //17 eta
  //18 phi2/4
  //19 Axx
  //20 Ayy
  //21 Azz
  //22 Axy
  //23 Axz
  //24 Ayz

  inline double BY_EM_OP_X (size_t l, size_t ibox)
  {
    double psi2 = pow(vars.Get_val (l, ibox,var_index[0])-vars.Get_val (l, ibox, coeff_index[17]),2)-vars.Get_val (l, ibox, coeff_index[18]);

    return vars.Lap (l, ibox, var_index[1])-2*vars.Get_val (l, ibox, coeff_index[14])/psi2;}

  inline double BY_EM_duOP_X (size_t l, size_t ibox)
  { return vars.duLap (l, ibox, var_index[1]);}



  inline double BY_EM_OP_Y (size_t l, size_t ibox)
  {
    double psi2 = pow(vars.Get_val (l, ibox,var_index[0])-vars.Get_val (l, ibox, coeff_index[17]),2)-vars.Get_val (l, ibox, coeff_index[18]);

    return vars.Lap (l, ibox, var_index[2])-2*vars.Get_val (l, ibox, coeff_index[15])/psi2;}

  inline double BY_EM_duOP_Y (size_t l, size_t ibox)
  {return vars.duLap (l, ibox, var_index[2]) ;}



  inline double BY_EM_OP_Z (size_t l, size_t ibox)
  {
    double psi2 = pow(vars.Get_val (l, ibox,var_index[0])-vars.Get_val (l, ibox, coeff_index[17]),2)-vars.Get_val (l, ibox, coeff_index[18]);

    return vars.Lap (l, ibox, var_index[3])-2*vars.Get_val (l, ibox, coeff_index[16])/psi2;}

  inline double BY_EM_duOP_Z (size_t l, size_t ibox)
  {return vars.duLap (l, ibox, var_index[3]);}


  inline double BY_EM_OP_L (size_t l, size_t ibox)
  {return vars.Lap (l, ibox, var_index[4]) -(
	   vars.Dx(l, ibox, var_index[1]) +
	   vars.Dy(l, ibox, var_index[2]) +
	   vars.Dz(l, ibox, var_index[3])) ;}

  inline double BY_EM_duOP_L (size_t l, size_t ibox)
  {return vars.duLap (l, ibox, var_index[4]); }



  inline double BY_EM_OP_U (size_t l, size_t ibox)
  {

    //(U^3 + A1 U^2 + A2 U + A3) Lap U 
    //+ (B (dx(U)^2+dy(U)^2+dz(U)^2)
    //+ (B1x U + B2x) Dx(U) 
    //+ (B1y U + B2y) Dy(U) 
    //+ (B1z U + B2z) Dz(U) 
    //+ C U^2 + D U + E + AijAij/(4 psi2^2) + BiBi)
    double uu = vars.Get_val (l, ibox,var_index[0]);
    double uu2 = uu*uu;

    double udx = vars.Dx (l, ibox,var_index[0]);
    double udy = vars.Dy (l, ibox,var_index[0]);
    double udz = vars.Dz (l, ibox,var_index[0]);

    double psi2 = pow(vars.Get_val (l, ibox,var_index[0])-vars.Get_val (l, ibox, coeff_index[17]),2)-vars.Get_val (l, ibox, coeff_index[18]);


    double DkXk = 2*(vars.Dx(l, ibox,var_index[1]) + vars.Dy(l, ibox,var_index[2])+vars.Dz(l, ibox,var_index[3]) -0.25*vars.Lap(l, ibox, var_index[4]) )/3;

   double Aij[3][3]; 

   Aij[0][0] =( 2*vars.Dx(l, ibox,var_index[1]) - DkXk -0.5*vars.DDx(l, ibox,var_index[4])   )
     + vars.Get_val (l, ibox, coeff_index[19]);


   Aij[1][1] =( 2*vars.Dy(l, ibox,var_index[2]) -DkXk -0.5*vars.DDy(l, ibox,var_index[4])   )
     + vars.Get_val (l, ibox, coeff_index[20]);

   Aij[2][2] =( 2*vars.Dz(l, ibox,var_index[3]) -DkXk -0.5*vars.DDz(l, ibox,var_index[4])  )
     + vars.Get_val (l, ibox, coeff_index[21]);


   Aij[0][1] =( vars.Dx(l, ibox,var_index[2])+vars.Dy(l, ibox,var_index[1]) -0.5*vars.DDxy(l, ibox,var_index[4]) )
     +vars.Get_val (l, ibox, coeff_index[22]);

   Aij[0][2] =( vars.Dx(l, ibox,var_index[3])+vars.Dz(l, ibox,var_index[1]) -0.5*vars.DDxz(l, ibox,var_index[4]) )
     + vars.Get_val (l, ibox, coeff_index[23]);
   
   Aij[1][2] =( vars.Dy(l, ibox,var_index[3])+vars.Dz(l, ibox,var_index[2]) -0.5*vars.DDyz(l, ibox,var_index[4]) )
     + vars.Get_val (l, ibox, coeff_index[24]);


    Aij[1][0] = Aij[0][1];
    Aij[2][0] = Aij[0][2];
    Aij[2][1] = Aij[1][2];


    double AijAij = 0;

    for(int i=0; i<3 ; i++)
      for(int j=0; j<3 ; j++)
	AijAij += Aij[i][j]*Aij[i][j];



    return (uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]) )*vars.Lap(l, ibox, var_index[0]) 
      + (vars.Get_val (l, ibox, coeff_index[3]) * ( udx*udx + udy*udy + udz*udz )
	 + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * udx 
	 + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * udy 
	 + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * udz 
	 + vars.Get_val (l, ibox, coeff_index[10]) * uu2 
	 + vars.Get_val (l, ibox, coeff_index[11]) * uu 
	 + vars.Get_val (l, ibox, coeff_index[12]) +
	 0.25*AijAij/(psi2*psi2) + vars.Get_val (l, ibox, coeff_index[13])
	 ); 



  }


  inline double BY_EM_duOP_U (size_t l, size_t ibox)
  {


    //( 3 U^2 +2 A1 U + A2) Lap U 
    //(U^3 + A1 U^2 + A2 U + A3) duLap U 
    //+ 2 B (dx(U) dudx(U) + dy(U) dudy(U) + dz(U) dudz(U) )
    //+ (B1x) Dx(U) 
    //+ (B1x U + B2x) duDx(U) 
    //+ (B1y) Dy(U) 
    //+ (B1y U + B2y) duDy(U) 
    //+ (B1z) Dz(U) 
    //+ (B1z U + B2z) duDz(U) 
    //+ 2 C U + D -  AijAij(U-eta)/(psi2^3) 
  

    //(U^3 + A1 U^2 + A2 U + A3) duLap U 
    //+ (3*U^2 + 2*A1 U + A2) Lap U 
    //+ 2*B (dx(U)*duDx(U) + dy(U)*duDy(U) + dz(U)*duDz(U))
    //+ (B1x U + B2x) duDx(U) 
    //+ B1x Dx(U) 
    //+ (B1y U + B2y) duDy(U) 
    //+ B1y Dy(U) 
    //+ (B1z U + B2z) duDz(U) 
    //+ B1z Dz(U) 
    //+ 2*C U + D 
    //- (u+eta) AijAij/( psi2^3)
    double uu = vars.Get_val (l, ibox,var_index[0]);
    double uu2 = uu*uu;

    double udx = vars.Dx (l, ibox,var_index[0]);
    double udy = vars.Dy (l, ibox,var_index[0]);
    double udz = vars.Dz (l, ibox,var_index[0]);

    double dudx = vars.duDx (l, ibox,var_index[0]);
    double dudy = vars.duDy (l, ibox,var_index[0]);
    double dudz = vars.duDz (l, ibox,var_index[0]);

    double psi2 = pow(uu-vars.Get_val (l, ibox, coeff_index[17]),2)-vars.Get_val (l, ibox, coeff_index[18]);


    double DkXk = 2*(vars.Dx(l, ibox,var_index[1]) + vars.Dy(l, ibox,var_index[2])+vars.Dz(l, ibox,var_index[3]) -0.25*vars.Lap(l, ibox, var_index[4]) )/3;

   double Aij[3][3]; 

   Aij[0][0] = 2*vars.Dx(l, ibox,var_index[1]) - DkXk -0.5*vars.DDx(l, ibox,var_index[4])
     + vars.Get_val (l, ibox, coeff_index[19]);


   Aij[1][1] = 2*vars.Dy(l, ibox,var_index[2]) -DkXk -0.5*vars.DDy(l, ibox,var_index[4])
     + vars.Get_val (l, ibox, coeff_index[20]);

   Aij[2][2] = 2*vars.Dz(l, ibox,var_index[3]) -DkXk -0.5*vars.DDz(l, ibox,var_index[4])
     + vars.Get_val (l, ibox, coeff_index[21]);


   Aij[0][1] = vars.Dx(l, ibox,var_index[2])+vars.Dy(l, ibox,var_index[1]) 
     -0.5*vars.DDxy(l, ibox,var_index[4]) 
     +vars.Get_val (l, ibox, coeff_index[22]);

   Aij[0][2] = vars.Dx(l, ibox,var_index[3])+vars.Dz(l, ibox,var_index[1]) 
     -0.5*vars.DDxz(l, ibox,var_index[4]) 
     + vars.Get_val (l, ibox, coeff_index[23]);
   
   Aij[1][2] = vars.Dy(l, ibox,var_index[3])+vars.Dz(l, ibox,var_index[2]) 
     -0.5*vars.DDyz(l, ibox,var_index[4]) 
     + vars.Get_val (l, ibox, coeff_index[24]);


    Aij[1][0] = Aij[0][1];
    Aij[2][0] = Aij[0][2];
    Aij[2][1] = Aij[1][2];


    double AijAij = 0;

    for(int i=0; i<3 ; i++)
      for(int j=0; j<3 ; j++)
	AijAij += Aij[i][j]*Aij[i][j];


    return (3*uu2 +2* vars.Get_val(l, ibox, coeff_index[0])*uu + vars.Get_val(l, ibox, coeff_index[1])  )*vars.Lap(l, ibox, var_index[0]) 
      +(uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]) )*vars.duLap(l, ibox, var_index[0]) 
      + 2*vars.Get_val (l, ibox, coeff_index[3]) * ( udx*vars.duDx (l, ibox,var_index[0]) + udy*vars.duDy (l, ibox,var_index[0]) + udz*vars.duDz (l, ibox,var_index[0]) )
      + (vars.Get_val (l, ibox, coeff_index[4]) ) * udx 
      + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * vars.duDx (l, ibox,var_index[0])
      + (vars.Get_val (l, ibox, coeff_index[5]) ) * udy 
      + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * vars.duDy (l, ibox,var_index[0]) 
      + (vars.Get_val (l, ibox, coeff_index[6]) ) * udz 
      + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * vars.duDz (l, ibox,var_index[0])
      + 2*vars.Get_val (l, ibox, coeff_index[10]) * uu
      + vars.Get_val (l, ibox, coeff_index[11])  
      - (uu-vars.Get_val (l, ibox, coeff_index[17]))*AijAij/pow(psi2,3); 


    /*    
    return vars.duLap(l, ibox, var_index[0]) 
      +( (2*vars.Get_val (l, ibox, coeff_index[3]) * ( udx*dudx + udy*dudy + udz*dudz )
	 + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * vars.duDx (l, ibox,var_index[0]) 
	 + vars.Get_val (l, ibox, coeff_index[4]) * vars.Dx (l, ibox,var_index[0]) 
	 + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * vars.duDy (l, ibox,var_index[0]) 
	 + vars.Get_val (l, ibox, coeff_index[5]) * vars.Dy (l, ibox,var_index[0]) 
	 + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * vars.duDz (l, ibox,var_index[0]) 
	 + vars.Get_val (l, ibox, coeff_index[6]) * vars.Dz (l, ibox,var_index[0]) 
	 + 2*vars.Get_val (l, ibox, coeff_index[10]) * uu 
	 + vars.Get_val (l, ibox, coeff_index[11])
	  -AijAij*(uu+vars.Get_val (l, ibox, coeff_index[17]))/pow(psi2,3)
 )*(uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]) ) 
	- (vars.Get_val (l, ibox, coeff_index[3]) * ( udx*udx + udy*udy + udz*udz )
	   + (vars.Get_val (l, ibox, coeff_index[4]) * uu + vars.Get_val (l, ibox, coeff_index[7])) * udx 
	   + (vars.Get_val (l, ibox, coeff_index[5]) * uu + vars.Get_val (l, ibox, coeff_index[8])) * udy 
	   + (vars.Get_val (l, ibox, coeff_index[6]) * uu + vars.Get_val (l, ibox, coeff_index[9])) * udz 
	   + vars.Get_val (l, ibox, coeff_index[10]) * uu2 
	   + vars.Get_val (l, ibox, coeff_index[11]) * uu 
	   + vars.Get_val (l, ibox, coeff_index[12]) +
	   0.25*AijAij/(psi2*psi2) + vars.Get_val (l, ibox, coeff_index[13])  )*(3*uu2 + 2*vars.Get_val(l, ibox, coeff_index[0])*uu + vars.Get_val(l, ibox, coeff_index[1])) 
	   )/pow(uu2*uu + vars.Get_val(l, ibox, coeff_index[0])*uu2 + vars.Get_val(l, ibox, coeff_index[1])*uu + vars.Get_val(l, ibox, coeff_index[2]),2); 
    */
    

  }





//======================== Variables and ffunctions =======================//





//============================== Brill waves =============================//


  double zero (double x, double y, double z)
  {
    return (0.0);
  }

  double one (double x, double y, double z)
  {
    return (1.0);
  }



  double scalar_phi_cI (double r)
  {
    return ps_phi0*tanh((r-ps_r0)*ps_isigma);
  }

  double d_scalar_phi_cI (double r)
  {
    return ps_d1/pow(cosh((r-ps_r0)*ps_isigma),4);
  }


  double scalar_phi_cII (double r)
  {
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);
  }

  double d_scalar_phi_cII (double r)
  {
    return ps_de1*(r-ps_r0)*(r-ps_r0)*exp(-2*ps_isigma*(r-ps_r0)*(r-ps_r0));
  }


  double scalar_phi_ground (double x){
    
    double l2=0.6;

    return (0.0481646*exp(-0.0581545*(x-1.8039e-08)*(x-1.8039e-08)/(l2))+0.298408*exp(-0.111412*(x+9.6741e-09)*(x+9.6741e-09)/l2)+ 0.42755*exp(-0.207156*(x-1.09822e-08)*(x-1.09822e-08)/l2)+0.204229*exp(-0.37742*(x+2.13778e-08)*(x+2.13778e-08)/l2)+0.021649*exp(-0.68406*(x-8.78608e-08)*(x-8.78608e-08)/l2))*l2;
    
  }


  double d_scalar_phi_ground (double x){
    
    double l2=10000;

    return M_PI*pow(-2*(0.0581545*(x-1.8039e-08)*0.0481646*exp(-0.0581545*(x-1.8039e-08)*(x-1.8039e-08)/(l2))+0.111412*(x+9.6741e-09)*0.298408*exp(-0.111412*(x+9.6741e-09)*(x+9.6741e-09)/l2)+0.207156*(x-1.09822e-08)* 0.42755*exp(-0.207156*(x-1.09822e-08)*(x-1.09822e-08)/l2)+0.37742*(x+2.13778e-08)*0.204229*exp(-0.37742*(x+2.13778e-08)*(x+2.13778e-08)/l2)+0.68406*(x-8.78608e-08)*0.021649*exp(-0.68406*(x-8.78608e-08)*(x-8.78608e-08)/l2))/(l2*l2),2);
    
  }


  double scalar_phi_cIV (double r)
  {
    return ps_phi0*0.5*(tanh((r+ps_r0)*ps_isigma)-tanh((r-ps_r0)*ps_isigma));
  }

  double d_scalar_phi_cIV (double r)
  {
    return ps_de1*0.25*pow(pow(1./cosh((r+ps_r0)*ps_isigma),2)-pow(1./cosh((r-ps_r0)*ps_isigma),2) ,2);
  }




  double brill (double x, double y, double z);


//============================= Gaussian test =============================//


  double Gaussian_source (double x, double y, double z);

  double Gaussian (double x, double y, double z);


//=============================== 1/r^2 test ==============================//


  double One_over_r2_source (double x, double y, double z);

  double One_over_r2 (double x, double y, double z);


//============================= Test ffunction =============================//


  double test_funct (double x, double y, double z);

  double test_funct_source (double x, double y, double z);


//=================== Test ffunction (for Robin boundary) ==================//


  double test_boundary (double x, double y, double z);

  double test_boundary_source (double x, double y, double z);


//===================== Test ffunction (for 2 boxes) =======================//



  double test_mfunct (double x, double y, double z);

  double test_mfunct_source (double x, double y, double z);


//===================== Test ffunction (non-linear) ========================//

  double Non_linear_source (double x, double y, double z);

  double Non_linear (double x, double y, double z);

  double Non_linear_b (double x, double y, double z);

  double Non_linear_c (double x, double y, double z);

  double Non_linear_d (double x, double y, double z);


//============================ Punctures ==================================//


  double punctures_b (double x, double y, double z);

  double punctures_c (double x, double y, double z);


//===================== Punctures_scalar_field =========================//


  double punctures_scalar_d (double x, double y, double z);

  double punctures_scalar_e (double x, double y, double z);

  double test_ps(double x, double y, double z);


//===================== Punctures_Gauss_Bonnet =========================//


  double phi (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}

  double phi_dx (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    if (r !=0)
      return -2*(r-ps_r0)*ps_isigma*x*phi(x,y,z)*x/r;
    else
      return 0;}

  
  double phi_dy (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    if (r !=0)
      return -2*(r-ps_r0)*ps_isigma*x*phi(x,y,z)*y/r;
    else
      return 0;}

  double phi_dz (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    if (r !=0)
      return -2*(r-ps_r0)*ps_isigma*x*phi(x,y,z)*z/r;
    else
      return 0;}

  double phi_dxdx (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    if (r !=0)
      return -2*(r-ps_r0)*ps_isigma*x*phi(x,y,z)*x/r;
    else
      return 0;}
  
  double phi_dydx (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}
  double phi_dzdx (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}

  double phi_dxdy (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}
  double phi_dydy (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}
  double phi_dzdy (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}

  double phi_dxdz (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}
  double phi_dydz (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}
  double phi_dzdz (double x, double y, double z)
  { double r = sqrt (x * x + y * y + z * z);
    return ps_phi0*exp(-(r-ps_r0)*(r-ps_r0)*ps_isigma);}
  
//============================ NOS_ID ==================================//


  int chebycoe_Extremes( valarray<double> &u, int n, int inv ); 

  double chebyeva( double a, double b, valarray<double> &c, int m, double x ); 

  int set_TOV();

  double  plegendre(int l, int m, double theta);

  void computeSphericalHarmonic(double theta, double phi, double *Yr, double *Yi, int l, int m);

  void add_perturbation();

  double ffact(double n);

  double fact(double n);

//============================== Test system ==============================//


  double system_u_a (double x, double y, double z);

  double system_v_a (double x, double y, double z);



  double system_u_source (double x, double y, double z);

  double system_v_source (double x, double y, double z);



  double system_u (double x, double y, double z);

  double system_v (double x, double y, double z);



//================================== Trumpet ===================================//

  double Trumpet_source (double x, double y, double z);

  double Trumpet (double x, double y, double z);



//================================== Bowen-York ===================================//


  double Bowen_York_X1 (double X, double Y, double Z)
  {

  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      double p_dot_x = bhpx[i]*x+ bhpy[i]*y+bhpz[i]*z;

      cval -= 0.25 * (7*bhpx[i]+p_dot_x * x/(r*r)) / r;

    }


  return (cval);


  }



  double Bowen_York_X2 (double X, double Y, double Z)
  {

  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      double p_dot_x = bhpx[i]*x+ bhpy[i]*y+bhpz[i]*z;

      cval -= 0.25 * (7*bhpy[i]+p_dot_x * x/(r*r)) / r;

    }


  return (cval);


  }


  double Bowen_York_X3 (double X, double Y, double Z)
  {

  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      double p_dot_x = bhpx[i]*x+ bhpy[i]*y+bhpz[i]*z;

      cval -= 0.25 * (7*bhpz[i]+p_dot_x * x/(r*r)) / r;

    }


  return (cval);


  }

  double Bowen_York_rho0(double x,double y,double z, double dx, double dy, double dz);

  double Bowen_York_rho1(double x,double y,double z, double dx, double dy, double dz);

  double Bowen_York_rho2(double x,double y,double z, double dx, double dy, double dz);


  double Bowen_York_rho0_test(double x,double y,double z);
  double Bowen_York_rho1_test(double x,double y,double z);
  double Bowen_York_rho2_test(double x,double y,double z);


//===================== Electromagnetic BH =========================//


  inline double EM_BH_A1 (double x, double y, double z) 
  { return 3*EM_BH_eta(x,y,z);}

  inline double EM_BH_A2 (double x, double y, double z) 
  { return (3*pow(EM_BH_eta(x,y,z),2)-0.25*pow(EM_BH_phi(x,y,z),2));}
  
  double EM_BH_A3 (double x, double y, double z) { 
    double eta = EM_BH_eta(x,y,z);
    double phi = EM_BH_phi(x,y,z);
    return (pow(eta,3)-0.25*eta*phi*phi);}

  inline double EM_BH_B (double x, double y, double z) 
  { return -0.25*pow(EM_BH_phi(x,y,z),2);}


  inline double EM_BH_B1x (double x, double y, double z) 
  { return 0.5*EM_BH_phi(x,y,z)*EM_BH_phi_dx(x,y,z);}

  inline double EM_BH_B1y (double x, double y, double z) 
  { return 0.5*EM_BH_phi(x,y,z)*EM_BH_phi_dy(x,y,z);}

  inline double EM_BH_B1z (double x, double y, double z)
  { return 0.5*EM_BH_phi(x,y,z)*EM_BH_phi_dz(x,y,z);}

  double EM_BH_B2x (double x, double y, double z)  { 
    double phi = EM_BH_phi(x,y,z);
    return (0.5*phi*(EM_BH_eta(x,y,z)*EM_BH_phi_dx(x,y,z) - phi * EM_BH_eta_dx(x,y,z)));}

  double EM_BH_B2y (double x, double y, double z) { 
    double phi = EM_BH_phi(x,y,z);
    return (0.5*phi*(EM_BH_eta(x,y,z)*EM_BH_phi_dy(x,y,z) - phi * EM_BH_eta_dy(x,y,z)));}

  double EM_BH_B2z (double x, double y, double z){ 
    double phi = EM_BH_phi(x,y,z);
    return (0.5*phi*(EM_BH_eta(x,y,z)*EM_BH_phi_dz(x,y,z) - phi * EM_BH_eta_dz(x,y,z)));}


  inline double EM_BH_C (double x, double y, double z) 
  { return -0.25*(pow(EM_BH_phi_dx(x,y,z),2)+pow(EM_BH_phi_dy(x,y,z),2)+pow(EM_BH_phi_dz(x,y,z),2));}

  double EM_BH_D (double x, double y, double z) 
  { return (0.5*EM_BH_phi(x,y,z)*(EM_BH_phi_dx(x,y,z)*EM_BH_eta_dx(x,y,z)+EM_BH_phi_dy(x,y,z)*EM_BH_eta_dy(x,y,z)+ EM_BH_phi_dz(x,y,z)*EM_BH_eta_dz(x,y,z))  -0.5*EM_BH_eta(x,y,z)*(pow(EM_BH_phi_dx(x,y,z),2)+pow(EM_BH_phi_dy(x,y,z),2)+pow(EM_BH_phi_dz(x,y,z),2)));}

  double EM_BH_E (double x, double y, double z){

    double eta = EM_BH_eta(x,y,z);
    double phi = EM_BH_phi(x,y,z);

    double eta_dx = EM_BH_eta_dx(x,y,z);
    double phi_dx = EM_BH_phi_dx(x,y,z);

    double eta_dy = EM_BH_eta_dy(x,y,z);
    double phi_dy = EM_BH_phi_dy(x,y,z);

    double eta_dz = EM_BH_eta_dz(x,y,z);
    double phi_dz = EM_BH_phi_dz(x,y,z);

    return -0.25*phi*phi*(eta_dx*eta_dx+eta_dy*eta_dy+eta_dz*eta_dz)
      +0.5*phi*eta*(phi_dx*eta_dx+phi_dy*eta_dy+phi_dz*eta_dz)
      +0.25*(1-eta*eta)*(phi_dx*phi_dx+phi_dy*phi_dy+phi_dz*phi_dz);
  }

  
  double EM_BH_eta(double x, double y, double z);

  double EM_BH_eta_dx(double x, double y, double z);
  double EM_BH_eta_dy(double x, double y, double z);
  double EM_BH_eta_dz(double x, double y, double z);

  double EM_BH_phi(double x, double y, double z);

  double EM_BH_phi_dx(double x, double y, double z);
  double EM_BH_phi_dy(double x, double y, double z);
  double EM_BH_phi_dz(double x, double y, double z);
  
  //Test 
  /*
  double EM_BH_eta(double x, double y, double z){ return exp(-T_ETA*(x*x+y*y+z*z));}
  double EM_BH_eta_dx(double x, double y, double z){ return -2*T_ETA*exp(-T_ETA*(x*x+y*y+z*z))*x;};
  double EM_BH_eta_dy(double x, double y, double z){ return -2*T_ETA*exp(-T_ETA*(x*x+y*y+z*z))*y;}
  double EM_BH_eta_dz(double x, double y, double z){ return -2*T_ETA*exp(-T_ETA*(x*x+y*y+z*z))*z;}

  double EM_BH_phi(double x, double y, double z){ return exp(-T_PHI*(x*x+y*y+z*z));}
  double EM_BH_phi_dx(double x, double y, double z){ return -2*T_PHI*exp(-T_PHI*(x*x+y*y+z*z))*x;}
  double EM_BH_phi_dy(double x, double y, double z){ return -2*T_PHI*exp(-T_PHI*(x*x+y*y+z*z))*y;}
  double EM_BH_phi_dz(double x, double y, double z){ return -2*T_PHI*exp(-T_PHI*(x*x+y*y+z*z))*z;}



  double EM_BH_rho_test(double x, double y, double z){

    double r2 = x*x + y*y + z*z;

    double t_eta = T_ETA;
    double t_u = 0.2;
    double t_phi = T_PHI;
    double u0 = 0.1;

    double eta=EM_BH_eta(x, y,z);
    double u=1-u0*exp(-t_u*r2);
    double phi= EM_BH_phi(x, y,z) ;

 
    
    return -2*exp(-r2*t_u)*u0*t_u*(2*r2*t_u-3) 
      -(r2*u0*(2*exp(-r2*(2+t_u))+exp(-2*r2*(1+t_u))*u0*(t_u-1))*(t_u-1)) / ( (eta+u)*( pow(eta+u,2)-0.25*phi*phi ) );
    


  }
  */




//===================== Bowen-York Electromagnetic BH =========================//


  inline double BY_EM_A1 (double x, double y, double z) 
  { return 3*BY_EM_eta(x,y,z);}

  inline double BY_EM_A2 (double x, double y, double z) 
  { return (3*pow(BY_EM_eta(x,y,z),2)-0.25*pow(BY_EM_phi(x,y,z),2));}
  
  double BY_EM_A3 (double x, double y, double z) { 
    double eta = BY_EM_eta(x,y,z);
    double phi = BY_EM_phi(x,y,z);
    return (pow(eta,3)-0.25*eta*phi*phi);}

  inline double BY_EM_B (double x, double y, double z) 
  { return -0.25*pow(BY_EM_phi(x,y,z),2);}


  inline double BY_EM_B1x (double x, double y, double z) 
  { return 0.5*BY_EM_phi(x,y,z)*BY_EM_phi_dx(x,y,z);}

  inline double BY_EM_B1y (double x, double y, double z) 
  { return 0.5*BY_EM_phi(x,y,z)*BY_EM_phi_dy(x,y,z);}

  inline double BY_EM_B1z (double x, double y, double z)
  { return 0.5*BY_EM_phi(x,y,z)*BY_EM_phi_dz(x,y,z);}

  double BY_EM_B2x (double x, double y, double z)  { 
    double phi = BY_EM_phi(x,y,z);
    return (0.5*phi*(BY_EM_eta(x,y,z)*BY_EM_phi_dx(x,y,z) - phi * BY_EM_eta_dx(x,y,z)));}

  double BY_EM_B2y (double x, double y, double z) { 
    double phi = BY_EM_phi(x,y,z);
    return (0.5*phi*(BY_EM_eta(x,y,z)*BY_EM_phi_dy(x,y,z) - phi * BY_EM_eta_dy(x,y,z)));}

  double BY_EM_B2z (double x, double y, double z){ 
    double phi = BY_EM_phi(x,y,z);
    return (0.5*phi*(BY_EM_eta(x,y,z)*BY_EM_phi_dz(x,y,z) - phi * BY_EM_eta_dz(x,y,z)));}


  inline double BY_EM_C (double x, double y, double z) 
  { return -0.25*(pow(BY_EM_phi_dx(x,y,z),2)+pow(BY_EM_phi_dy(x,y,z),2)+pow(BY_EM_phi_dz(x,y,z),2));}

  double BY_EM_D (double x, double y, double z) 
  { return (0.5*BY_EM_phi(x,y,z)*(BY_EM_phi_dx(x,y,z)*BY_EM_eta_dx(x,y,z)+BY_EM_phi_dy(x,y,z)*BY_EM_eta_dy(x,y,z)+ BY_EM_phi_dz(x,y,z)*BY_EM_eta_dz(x,y,z))  -0.5*BY_EM_eta(x,y,z)*(pow(BY_EM_phi_dx(x,y,z),2)+pow(BY_EM_phi_dy(x,y,z),2)+pow(BY_EM_phi_dz(x,y,z),2)));}

  double BY_EM_E (double x, double y, double z){

    double eta = BY_EM_eta(x,y,z);
    double phi = BY_EM_phi(x,y,z);

    double eta_dx = BY_EM_eta_dx(x,y,z);
    double phi_dx = BY_EM_phi_dx(x,y,z);

    double eta_dy = BY_EM_eta_dy(x,y,z);
    double phi_dy = BY_EM_phi_dy(x,y,z);

    double eta_dz = BY_EM_eta_dz(x,y,z);
    double phi_dz = BY_EM_phi_dz(x,y,z);

    return -0.25*phi*phi*(eta_dx*eta_dx+eta_dy*eta_dy+eta_dz*eta_dz)
      +0.5*phi*eta*(phi_dx*eta_dx+phi_dy*eta_dy+phi_dz*eta_dz)
      +0.25*(1-eta*eta)*(phi_dx*phi_dx+phi_dy*phi_dy+phi_dz*phi_dz);
  }

  
  double BY_EM_eta(double x, double y, double z)
  {return EM_BH_eta(x,y,z);  }

  double BY_EM_eta_dx(double x, double y, double z)
  {return EM_BH_eta_dx(x,y,z);  }

  double BY_EM_eta_dy(double x, double y, double z)
  {return EM_BH_eta_dy(x,y,z);  }

  double BY_EM_eta_dz(double x, double y, double z)
  {return EM_BH_eta_dz(x,y,z);  }

  double BY_EM_phi(double x, double y, double z);

  double BY_EM_phi_dx(double x, double y, double z);
  double BY_EM_phi_dy(double x, double y, double z);
  double BY_EM_phi_dz(double x, double y, double z);

  double BY_EM_BiBi(double x, double y, double z)
  { 

    double Bx = BY_EM_Bfield_x(x,y,z); 
    double By = BY_EM_Bfield_y(x,y,z); 
    double Bz = BY_EM_Bfield_z(x,y,z); 
    return  0.5*(Bx*Bx+By*By+Bz*Bz); }


  inline double BY_EM_Efield_x(double x,double y,double z)
  { return -BY_EM_phi_dx(x,y,z);} 

  inline double BY_EM_Efield_y(double x,double y,double z)
  { return -BY_EM_phi_dy(x,y,z);} 

  inline double BY_EM_Efield_z(double x,double y,double z)
  { return -BY_EM_phi_dz(x,y,z);} 


  double BY_EM_Bfield_x(double x,double y,double z);
  double BY_EM_Bfield_y(double x,double y,double z);
  double BY_EM_Bfield_z(double x,double y,double z);


  inline double BY_EM_Jx(double x,double y,double z)
  {return BY_EM_Efield_y(x,y,z)*BY_EM_Bfield_z(x,y,z)-BY_EM_Efield_z(x,y,z)*BY_EM_Bfield_y(x,y,z);}


  inline double BY_EM_Jy(double x,double y,double z)
  {return BY_EM_Efield_z(x,y,z)*BY_EM_Bfield_x(x,y,z)-BY_EM_Efield_x(x,y,z)*BY_EM_Bfield_z(x,y,z);}

  inline double BY_EM_Jz(double x,double y,double z)
  {return BY_EM_Efield_x(x,y,z)*BY_EM_Bfield_y(x,y,z)-BY_EM_Efield_y(x,y,z)*BY_EM_Bfield_x(x,y,z);}



//================================== End ===================================//

};





void Boundary_Face (double (ffunction::*dxu) (void),
		    double (ffunction::*du) (void),
		    double (ffunction::*cd) (void),
		    ffunction * U, t_iterator face, size_t pn, ffunction * S);




void Boundary_Edge (double (ffunction::*dxu1) (void),
		    double (ffunction::*du1) (void),
		    double (ffunction::*cd1) (void),
		    double (ffunction::*dxu2) (void),
		    double (ffunction::*du2) (void),
		    double (ffunction::*cd2) (void),
		    ffunction * U, t_iterator face, size_t pn, ffunction * S);


void Boundary_Corner (double (ffunction::*dxu1) (void),
		      double (ffunction::*du1) (void),
		      double (ffunction::*cd1) (void),
		      double (ffunction::*dxu2) (void),
		      double (ffunction::*du2) (void),
		      double (ffunction::*cd2) (void),
		      double (ffunction::*dxu3) (void),
		      double (ffunction::*du3) (void),
		      double (ffunction::*cd3) (void),
		      ffunction * U, t_iterator face, size_t pn, ffunction * S);











void App_Boundary_Face (double (ffunction::*dxu) (void),
			double (ffunction::*du) (void),
			double (ffunction::*cd) (void),
			ffunction * U,
			t_iterator face, size_t pn, ffunction * R);




void App_Boundary_Edge (double (ffunction::*dxu1) (void),
			double (ffunction::*du1) (void),
			double (ffunction::*cd1) (void),
			double (ffunction::*dxu2) (void),
			double (ffunction::*du2) (void),
			double (ffunction::*cd2) (void),
			ffunction * U,
			t_iterator face, size_t pn, ffunction * R);


void App_Boundary_Corner (double (ffunction::*dxu1) (void),
			  double (ffunction::*du1) (void),
			  double (ffunction::*cd1) (void),
			  double (ffunction::*dxu2) (void),
			  double (ffunction::*du2) (void),
			  double (ffunction::*cd2) (void),
			  double (ffunction::*dxu3) (void),
			  double (ffunction::*du3) (void),
			  double (ffunction::*cd3) (void),
			  ffunction * U,
			  t_iterator face, size_t pn, ffunction * R);




inline double
Khat (double a, double R, double M)
{
  double R2 = R * R;

  double a2 = a * a;

  return (0.25 *
	  (2 * a * (6.0 * M + 4.0 * R * (a2 - 1)) /
	   (2 * M + R * (a2 - 2 * a - 1)) + a2 - 1.0) / R2);
};



double a2R (double a, double M);


double r2a (double r, double M);


double rhs_arM (double a, double r, double M);


#endif
