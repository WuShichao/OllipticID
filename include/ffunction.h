



//=========================== Olliptic ffunction.h =========================//
//  
//  Grid fnctions (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//










#ifndef ffunction_H

#define ffunction_H



//============================ Standar libraries =========================//

#include <valarray>

#include <limits>



//============================ Olliptic librarie =========================//

#include "domain.h"




//== constants ==//




enum t_iterator
{ ALL, INSIDE, BOUNDARY,

  CORNER_XIYIZI, CORNER_XIYIZF,
  CORNER_XFYIZI, CORNER_XFYIZF,

  CORNER_XIYFZI, CORNER_XIYFZF,
  CORNER_XFYFZI, CORNER_XFYFZF,

  EDGE_XIYI, EDGE_XIYF,
  EDGE_XFYI, EDGE_XFYF,

  EDGE_XIZI, EDGE_XIZF,
  EDGE_XFZI, EDGE_XFZF,

  EDGE_YIZI, EDGE_YIZF,
  EDGE_YFZI, EDGE_YFZF,


  FACE_XI, FACE_XF,
  FACE_YI, FACE_YF,
  FACE_ZI, FACE_ZF
};


enum cut
{ CUT_X, CUT_Y, CUT_Z };

enum coord
{ CoordX, CoordY, CoordZ };


enum m_interpol
{ LINEAR, LAGRANGE, NEWTON, HERMITE, AVERAGE };


//typedef double (*pt2ffunction)(double x, double y, double z);


using namespace std;



class ffunction
{


  // total number of cpus and rank
  // !!! Cambiar esto, ya se definio en domain.

  int total_nodes;

  int my_rank;



    valarray < double >f;






  int init_inter;

  int init_corner[8];

  int init_edge[12];

  int init_face[6];

  int n;




  t_iterator iter;

  int index;

  int index_end;




  size_t order;

  size_t int_ord;


  size_t it_2d;


  double dx;

  double dy;

  double dz;


  double idx;

  double idy;

  double idz;


  double idxO;

  double idyO;

  double idzO;



  double idx2;

  double idy2;

  double idz2;


  double i180;

  double pi;

  double N_L1;

  double N_L2;

  double N_LInf;







  void Sync_Bound_Xi ();

  void Sync_Bound_Xf ();


  void Sync_Bound_Yi ();

  void Sync_Bound_Yf ();


  void Sync_Bound_Zi ();

  void Sync_Bound_Zf ();





  double (ffunction::*lap) (void);

  double (ffunction::*dulap) (void);

  double (ffunction::*lapB) (void);


  double (ffunction::*d_x) (void);

  double (ffunction::*dud_x) (void);


  double (ffunction::*d_y) (void);

  double (ffunction::*dud_y) (void);


  double (ffunction::*d_z) (void);

  double (ffunction::*dud_z) (void);



  double (ffunction::*dd_x) (void);

  double (ffunction::*dudd_x) (void);


  double (ffunction::*dd_y) (void);

  double (ffunction::*dudd_y) (void);


  double (ffunction::*dd_z) (void);

  double (ffunction::*dudd_z) (void);



  double (ffunction::*dd_xy) (void);

  double (ffunction::*dudd_xy) (void);


  double (ffunction::*dd_xz) (void);

  double (ffunction::*dudd_xz) (void);


  double (ffunction::*dd_yz) (void);

  double (ffunction::*dudd_yz) (void);




  double (ffunction::*interpol) (double x, double y, double z);


  double Linear_Interpol_In (double x, double y, double z);

  double Lagrange_Interpol_In (double x, double y, double z);

  double LagrangePoly (coord co, const double t, int i, int ni, int nf);


  double Newton_Interpol_In (double x, double y, double z);

  double NewtonPoly (const double t, valarray < double >F,
		     valarray < double >x);


  double Hermite_Interpol_In (double x, double y, double z);

  double Hermite_Poly (const double t,
		       const double ti,
		       const double tip1,
		       const double h,
		       const double fti,
		       const double ftip1,
		       const double dfti, const double dftip1);

  double Average_Interpol_In (double x, double y, double z);



    protected:domain * D;


public:




  //default constructor (using a default domain)

    ffunction ()
  {
    domain Dom;
      Make (&Dom);
  }



  //destructor 
   ~ffunction ()
  {
  }



  //constructor

  ffunction (domain * Dom)
  {
    Make (Dom);
  }





  void Make (domain * Dom);



//== Print methods ==//


  friend ostream & operator<< (ostream &, ffunction &);



  void Print_Info (ostream &);

  void Print_3D (string name,
		 double (*func_to_comp) (double x, double y, double z) =
		 NULL);


  void Print_2D (string name,
		 const cut c = CUT_X,
		 const double x1 = 0,
		 bool interpolar = true,
		 double (*func_to_comp) (double x, double y, double z) = NULL,
		 m_interpol minter = LAGRANGE, size_t order = 2);




  void Print_1D (string name,
		 const cut c = CUT_X,
		 const double x1 = 0,
		 const double x2 = 0,
		 bool interpolar = true,
		 double (*func_to_comp) (double x, double y, double z) = NULL,
		 m_interpol minter = LAGRANGE, size_t order = 2);




  void Print_t1D (string name, const double t);


  void Print_t2D (string name, size_t cyc);







  //== inline methods ==//

  inline int Get_rank ()
  {
    return (my_rank);
  }

  inline int Get_nodes ()
  {
    return (total_nodes);
  }


  inline void Set_Iterator (t_iterator it = INSIDE)
  {
    iter = it;
    Begin_End ();
  }


  inline const int Get_N () const
  {
    return (n);
  }


  inline const double Get_val (const int i) const
  {
    return (f[i]);
  }

  inline const double Get_val (const int i, const int j, const int k) const
  {
    return (f[Get_Index (i, j, k)]);
  }

  inline const double operator () (int i, int j, int k) const
  {
    return (Get_val (i, j, k));
  }

  inline const double operator () (int i) const
  {
    return (Get_val (i));
  }


  inline const size_t Get_order () const
  {
    return (order);
  }


  inline const double Get_idx () const
  {
    return (idx);
  }

  inline const double Get_idy () const
  {
    return (idy);
  }

  inline const double Get_idz () const
  {
    return (idz);
  }



  inline const double Get_idx2 () const
  {
    return (idx2);
  }

  inline const double Get_idy2 () const
  {
    return (idy2);
  }

  inline const double Get_idz2 () const
  {
    return (idz2);
  }



  inline const double Get_x () const
  {
    return (Get_x (index));
  }

  inline const double Get_y () const
  {
    return (Get_y (index));
  }

  inline const double Get_z () const
  {
    return (Get_z (index));
  }

  inline const double Get_val () const
  {
    return (Get_val (index));
  }

  inline const double operator () () const
  {
    return (Get_val (index));
  }

  inline const int Get_Index () const
  {
    return (index);
  }



  inline double Get_X ()
  {
    return (Get_x (index));
  }

  inline double Get_Y ()
  {
    return (Get_y (index));
  }

  inline double Get_Z ()
  {
    return (Get_z (index));
  }



  inline void Set (const int i, const double y)
  {
    f[i] = y;
  }


  inline void Set (const int i, const int j, const int k, const double y)
  {
    f[Get_Index (i, j, k)] = y;
  }


  inline void Set (const double y)
  {
    f[index] = y;
  }


  inline void Next ()
  {
    index++;
  };

  inline bool new_End ()
  {
    return (index < index_end);
  }

  inline bool End ()
  {
    index++;
    return (index < index_end);
  }


  inline ffunction & operator-= (const double u)
  {
    double old = Get_val ();

    f[index] = old - u;
    return (*this);
  }

  inline ffunction & operator+= (const double u)
  {
    double old = Get_val ();

    f[index] = old + u;
    return (*this);
  }


  inline void Norms (double &l1, double &l2, double &lInf)
  {
    l1 = Get_Norm_L1 ();
    l2 = Get_Norm_L2 ();
    lInf = Get_Norm_LInf ();
  }




  inline void Reset_Norms ()
  {
    N_L1 = 0;
    N_L2 = 0;
    N_LInf = 0;
  }



  inline const double Get_Norm_L1 () const
  {
    return (N_L1 * dx * dy * dz);
  }

  inline const double Get_Norm_L2 () const
  {
    return (sqrt (N_L2 * dx * dy * dz));
  }

  inline const double Get_Norm_LInf () const
  {
    return (N_LInf);
  }


  inline void Get_MIndex (int &i, int &j, int &k)
  {
    Get_MIndex (i, j, k, index);
  }


  inline void Get_MIndex (int &i, int &j, int &k, const int indx)
  {
    i = D->Get_Index_I (indx);
    j = D->Get_Index_J (indx);
    k = D->Get_Index_K (indx);
  }



  inline const int Fix_Index_i (const int ni, const int i, const int nf) const
  {

    if (ni >= 0 && nf < int (Get_nx ()))

      return (i);

    else
     if (ni < 0)

      return (i - ni);

    else

      return (i - nf + Get_nx () - 1);



  };



  inline const int Fix_Index_j (const int ni, const int i, const int nf) const
  {

    if (ni >= 0 && nf < int (Get_ny ()))

      return (i);

    else
     if (ni < 0)

      return (i - ni);

    else

      return (i - nf + Get_ny () - 1);


  };


  inline const int Fix_Index_k (const int ni, const int i, const int nf) const
  {

    if (ni >= 0 && nf < int (Get_nz ()))

      return (i);

    else
     if (ni < 0)

      return (i - ni);

    else

      return (i - nf + Get_nz () - 1);


  };


  void Set_Interpol (m_interpol minter = LINEAR, size_t ord = 1);

  inline double Interpol_In (double x, double y, double z)
  {
    return ((this->*interpol) (x, y, z));
  }



//== interface ==//
//!!! cambiar, usar un apuntador directo al dominio



  inline domain *Get_domain ()
  {
    return (D);
  }


  inline const size_t size () const
  {
    return (D->Get_N ());
  }



  inline const size_t Get_nx () const
  {
    return (D->Get_Nx ());
  }

  inline const size_t Get_ny () const
  {
    return (D->Get_Ny ());
  }

  inline const size_t Get_nz () const
  {
    return (D->Get_Nz ());
  }



  inline const double Get_x (const int i) const
  {
    return (D->Get_x (Get_i (i)));
  }

  inline const double Get_x (const int i, const int j, const int k) const
  {
    return (D->Get_x (i));
  }



  inline const double Get_y (const int i) const
  {
    return (D->Get_y (Get_j (i)));
  }

  inline const double Get_y (const int i, const int j, const int k) const
  {
    return (D->Get_y (j));
  }



  inline const double Get_z (const int i) const
  {
    return (D->Get_z (Get_k (i)));
  }

  inline const double Get_z (const int i, const int j, const int k) const
  {
    return (D->Get_z (k));
  }





  inline const double Get_dx () const
  {
    return (D->Get_dx ());
  }

  inline const double Get_dy () const
  {
    return (D->Get_dy ());
  }

  inline const double Get_dz () const
  {
    return (D->Get_dz ());
  }



  inline const int Get_Index (const int i, const int j, const int k) const
  {
    return (D->Get_Index (i, j, k));
  }




  inline const int Get_i (const int indx) const
  {
    return (D->Get_Index_I (indx));
  }

  inline const int Get_j (const int indx) const
  {
    return (D->Get_Index_J (indx));
  }

  inline const int Get_k (const int indx) const
  {
    return (D->Get_Index_K (indx));
  }


  inline const int Get_i (const double x) const
  {
    return (D->Get_i (x));
  }

  inline const int Get_j (const double y) const
  {
    return (D->Get_j (y));
  }

  inline const int Get_k (const double z) const
  {
    return (D->Get_k (z));
  }


  inline const bool In_Cell (int i, int j, int k) const
  {
    return (i - 1 <= Get_i (index) && Get_i (index) <= i &&
	    j - 1 <= Get_j (index) && Get_j (index) <= j &&
	    k - 1 <= Get_k (index) && Get_k (index) <= k);
  }


  inline const int Is_In_rank (const double x, const double y, const double z)
  {
    return (D->Is_In_rank (x, y, z));
  };

  inline const bool Is_In (const double x, const double y, const double z,
			   const size_t rank)
  {
    return (D->Is_In (x, y, z, rank));
  };

  inline const bool Is_In (const double x, const double y, const double z)
  {
    return (D->is_in (x, y, z));
  };

  inline const bool Is_In_Nb (const double x, const double y, const double z,
			      const size_t rank)
  {
    return (D->is_in_nb (x, y, z, rank));
  };

  inline const bool Is_In_Nb (const double x, const double y, const double z)
  {
    return (D->is_in_nb (x, y, z));
  };

  inline const bool Is_Near (const double x, const double y, const double z)
  {
    return (D->Is_Near (x, y, z));
  };

  inline const bool Is_In_Boundary (const double x, const double y,
				    const double z)
  {
    return (D->Is_In_Boundary (x, y, z));
  };

  inline const bool Is_Near (const double x, const double y, const double z,
			     const double Dx, const double Dy,
			     const double Dz)
  {
    return (D->Is_Near (x, y, z, Dx, Dy, Dz));
  };



//== Operators ==//


  inline ffunction & operator= (const double u)
  {
    f[index] = u;
    return (*this);
  }

  inline ffunction & operator* (const double u)
  {
    f *= u;
    return (*this);
  }

  inline ffunction & operator+= (const ffunction & u)
  {
    f += u.f;
    return (*this);
  }


  inline ffunction & operator-= (const ffunction & u)
  {
    f -= u.f;
    return (*this);
  }


  inline friend ffunction const operator+ (const ffunction & u,
					  const ffunction & v)
  {
    return (ffunction (u) += v);
  }


  inline friend ffunction const operator- (const ffunction & u,
					  const ffunction & v)
  {
    return (ffunction (u) -= v);
  }









  bool Begin_End ();



  double NormL2 ();

  double NormL1 ();

  double NormLInf ();



  void Add_To_Norms ();

  void Sync_Norms ();



  void Sync ();



  void Transferir (ffunction & s, t_iterator iter = ALL, m_interpol method =
		   LINEAR, size_t ord = 1, bool all = true, bool to_zero =
		   false);


  void Set_In (ffunction & s, domain * Dom, bool all = true);


  void Add (ffunction &);




  void Prolong (ffunction &);

  void Restrict (ffunction &);









//== differential operators   ==//





  inline double Lap ()
  {
    return ((this->*lap) ());
  }

  inline double duLap ()
  {
    return ((this->*dulap) ());
  }

  inline double LapB ()
  {
    return ((this->*lapB) ());
  }


  inline double Dx ()
  {
    return ((this->*d_x) ());
  }

  inline double duDx ()
  {
    return ((this->*dud_x) ());
  }


  inline double Dy ()
  {
    return ((this->*d_y) ());
  }

  inline double duDy ()
  {
    return ((this->*dud_y) ());
  }


  inline double Dz ()
  {
    return ((this->*d_z) ());
  }

  inline double duDz ()
  {
    return ((this->*dud_z) ());
  }




  inline double DDx ()
  {
    return ((this->*dd_x) ());
  }

  inline double duDDx ()
  {
    return ((this->*dudd_x) ());
  }


  inline double DDy ()
  {
    return ((this->*dd_y) ());
  }

  inline double duDDy ()
  {
    return ((this->*dudd_y) ());
  }


  inline double DDz ()
  {
    return ((this->*dd_z) ());
  }

  inline double duDDz ()
  {
    return ((this->*dudd_z) ());
  }



  inline double DDxy ()
  {
    return ((this->*dd_xy) ());
  }

  inline double duDDxy ()
  {
    return ((this->*dudd_xy) ());
  }


  inline double DDxz ()
  {
    return ((this->*dd_xz) ());
  }

  inline double duDDxz ()
  {
    return ((this->*dudd_xz) ());
  }


  inline double DDyz ()
  {
    return ((this->*dd_yz) ());
  }

  inline double duDDyz ()
  {
    return ((this->*dudd_yz) ());
  }



  double dxl2 ();

  inline double dudxl2 ()
  {
    return (3.0 * 0.5 * idx);
  }

  double dyl2 ();

  inline double dudyl2 ()
  {
    return (3.0 * 0.5 * idy);
  }

  double dzl2 ();

  inline double dudzl2 ()
  {
    return (3.0 * 0.5 * idz);
  }


  double dxc2 ();

  double dyc2 ();

  double dzc2 ();


  inline double dxc2 (double fim1, double fip1)
  {
    return ((fip1 - fim1) * 0.5 * idx);
  }

  inline double dyc2 (double fim1, double fip1)
  {
    return ((fip1 - fim1) * 0.5 * idy);
  }

  inline double dzc2 (double fim1, double fip1)
  {
    return ((fip1 - fim1) * 0.5 * idz);
  }


  inline double dx1 (double fim1, double fi)
  {
    return ((fi - fim1) * idx);
  }

  inline double dy1 (double fim1, double fi)
  {
    return ((fi - fim1) * idy);
  }

  inline double dz1 (double fim1, double fi)
  {
    return ((fi - fim1) * idz);
  }



  double dxr2 ();

  inline double dudxr2 ()
  {
    return (-3.0 * 0.5 * idx);
  }

  double dyr2 ();

  inline double dudyr2 ()
  {
    return (-3.0 * 0.5 * idy);
  }

  double dzr2 ();

  inline double dudzr2 ()
  {
    return (-3.0 * 0.5 * idz);
  }



  double dx2 ();

  double dudx2 ();

  double dy2 ();

  double dudy2 ();

  double dz2 ();

  double dudz2 ();





  double dx3 ();

  double dudx3 ();

  double dy3 ();

  double dudy3 ();

  double dz3 ();

  double dudz3 ();




  double dx4 ();

  double dudx4 ();

  double dy4 ();

  double dudy4 ();

  double dz4 ();

  double dudz4 ();


  double dx6 ();

  double dudx6 ();

  double dy6 ();

  double dudy6 ();

  double dz6 ();

  double dudz6 ();



  double dx8 ();

  double dudx8 ();

  double dy8 ();

  double dudy8 ();

  double dz8 ();

  double dudz8 ();





  double ddxc2 ();

  inline double duddxc2 ()
  {
    return (-2.0 * idx2);
  }

  double ddxl2 ();

  inline double duddxl2 ()
  {
    return (idx2);
  }

  double ddxr2 ();

  inline double duddxr2 ()
  {
    return (idx2);
  }


  double ddyc2 ();

  inline double duddyc2 ()
  {
    return (-2.0 * idy2);
  }

  double ddzc2 ();

  inline double duddzc2 ()
  {
    return (-2.0 * idz2);
  }


  inline double Lapc2 ()
  {

    int i,
      j,
      k;

    Get_MIndex (i, j, k);


    return ((Get_val (i + 1, j, k) - 2.0 * Get_val () +
	     Get_val (i - 1, j, k)) * idx2 + (Get_val (i, j + 1,
						       k) - 2.0 * Get_val () +
					      Get_val (i, j - 1,
						       k)) * idy2 +
	    (Get_val (i, j, k + 1) - 2.0 * Get_val () +
	     Get_val (i, j, k - 1)) * idz2);

  }

  inline double duLapc2 ()
  {
    return (-2.0 * (idx2 + idy2 + idz2));
  }

  double LapC2 ();

  double duLapC2 ();

  double Lapc2 (ffunction & u);

  double LapB2 ();



  double LapS (int i0, int j0, int k0);

  inline double duLapS ()
  {
    return (idx2 + idy2 + idz2);
  }




  double ddxyc2 ();

  inline double duddxyc2 ()
  {
    return (0.0);
  }

  double ddxzc2 ();

  inline double duddxzc2 ()
  {
    return (0.0);
  }

  double ddyzc2 ();

  inline double duddyzc2 ()
  {
    return (0.0);
  }



  double ddxc4 ();

  inline double duddxc4 ()
  {
    return (-2.5 * idx2);
  }

  double ddyc4 ();

  inline double duddyc4 ()
  {
    return (-2.5 * idy2);
  }

  double ddzc4 ();

  inline double duddzc4 ()
  {
    return (-2.5 * idz2);
  }


  double Lapc4 ();

  inline double duLapc4 ()
  {
    return (-2.5 * (idx2 + idy2 + idz2));
  }


  double ddxyc4 ();

  inline double duddxyc4 ()
  {
    return (0.0);
  }

  double ddxzc4 ();

  inline double duddxzc4 ()
  {
    return (0.0);
  }

  double ddyzc4 ();

  inline double duddyzc4 ()
  {
    return (0.0);
  }




  double ddxc6 ();

  inline double duddxc6 ()
  {
    return (-49 * idx2 / 18.0);
  }

  double ddyc6 ();

  inline double duddyc6 ()
  {
    return (-49 * idy2 / 18.0);
  }

  double ddzc6 ();

  inline double duddzc6 ()
  {
    return (-49 * idz2 / 18.0);
  }


//    double Lapc6();
  inline double Lapc6 ()
  {


    int i,
      j,
      k;


    Get_MIndex (i, j, k);


    return (2.0 * (Get_val (i + 3, j, k) + Get_val (i - 3, j, k))
	     - 27.0 * (Get_val (i + 2, j, k) + Get_val (i - 2, j, k)
		       - 10.0 * (Get_val (i + 1, j, k) + Get_val (i - 1, j,k))) -
	     490.0 * Get_val ()) * idx2 * i180 +
      
	    (2.0 * (Get_val (i, j + 3, k) + Get_val (i, j - 3, k)) -
	     27.0 * (Get_val (i, j + 2, k) + Get_val (i, j - 2, k) -
		     10.0 * (Get_val (i, j + 1, k) + Get_val (i, j - 1, k))) -
	     490.0 * Get_val ()) * idy2 * i180 +
      
	    (2.0 * (Get_val (i, j, k + 3) + Get_val (i, j, k - 3)) -
	     27.0 * (Get_val (i, j, k + 2) + Get_val (i, j, k - 2) -
		     10.0 * (Get_val (i, j, k + 1) + Get_val (i, j, k - 1))) -
	     490.0 * Get_val ()) * idz2 * i180;




  }
  inline double duLapc6 ()
  {
    return (-490 * (idx2 + idy2 + idz2) * i180);
  }



  double ddxyc6 ();

  inline double duddxyc6 ()
  {
    return (0.0);
  }

  double ddxzc6 ();

  inline double duddxzc6 ()
  {
    return (0.0);
  }

  double ddyzc6 ();

  inline double duddyzc6 ()
  {
    return (0.0);
  }




  double ddxc8 ();

  inline double duddxc8 ()
  {
    return (-205 * idx2 / 72.0);
  }

  double ddyc8 ();

  inline double duddyc8 ()
  {
    return (-205 * idy2 / 72.0);
  }

  double ddzc8 ();

  inline double duddzc8 ()
  {
    return (-205 * idz2 / 72.0);
  }


  double Lapc8 ();

  inline double duLapc8 ()
  {
    return (-205 * (idx2 + idy2 + idz2) / 72.0);
  }


  double ddxyc8 ();

  inline double duddxyc8 ()
  {
    return (0.0);
  }

  double ddxzc8 ();

  inline double duddxzc8 ()
  {
    return (0.0);
  }

  double ddyzc8 ();

  inline double duddyzc8 ()
  {
    return (0.0);
  }




  inline double max ()
  {
    return (f.max ());
  }

  inline double min ()
  {
    return (f.min ());
  }




};










#endif // ffunction_H
