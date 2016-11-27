







#ifndef DOMAIN_H

#define DOMAIN_H



//============================ Standar libraries =========================//


#include <iostream>

#include <iomanip>

#include <math.h>

#include <fstream>

#include <string>

#include <sstream>

#include <vector>

#include <valarray>

#include <limits>

#include <cstdlib>



#include <interface.h>

//=============================== MPI library ============================//


#ifdef OLLIN_MPI

#include <mpi.h>

#endif


//============================== End headers =============================//




//index order macro

#define INDEXORDER_IJK i + (j + k * ny[my_rank]) * nx[my_rank]


//== constants ==//




enum corner
{ XiYiZi, XiYiZf,
  XfYiZi, XfYiZf,
  XiYfZi, XiYfZf,
  XfYfZi, XfYfZf
};

enum edge
{ XiYi, XiYf,
  XfYi, XfYf,
  XiZi, XiZf,
  XfZi, XfZf,
  YiZi, YiZf,
  YfZi, YfZf
};

enum face
{ Xi, Xf, Yi, Yf, Zi, Zf };


//kind of symmetry

enum t_symmetry
{ NONE, REFLECT, ROTANT };





using namespace std;





class domain
{


  // number of grid points in each direction (for each cpu). 

  vector < size_t > nx;

  vector < size_t > ny;

  vector < size_t > nz;





  //number of buffers points (example inner points = nx - 2 * buff_nx)

  size_t buff_nx;

  size_t buff_ny;

  size_t buff_nz;


  size_t order;


    vector < int >ijk;


    vector < int >I;

    vector < int >J;

    vector < int >K;


  //limit coordinates (for each cpu)

    vector < double >xi;

    vector < double >xf;


    vector < double >yi;

    vector < double >yf;


    vector < double >zi;

    vector < double >zf;


  //Global domain limits


  double gxi;

  double gxf;


  double gyi;

  double gyf;


  double gzi;

  double gzf;


  //domain's center


  double cx;

  double cy;

  double cz;


  //Lenght

  double Lx;

  double Ly;

  double Lz;



  //grid size in each direction


  double dx;

  double dy;

  double dz;



  // total number of cpus and rank (node)


  size_t size;

  size_t my_rank;



  // kind of boundary 

  // n >= 0  : boundary with cpu n
  // n == -1 : physical boundary
  // n == -2 : symmetry, reflection
  // n == -3 : symmetry, rotation

  int bound_xi;

  int bound_xf;



  int bound_yi;

  int bound_yf;



  int bound_zi;

  int bound_zf;


  // kind of symmetry in each direction
  // example: quadrant domain should have
  // symm_x = REFLECT
  // symm_y = REFLECT
  // symm_z = REFLECT

  t_symmetry symm_x;

  t_symmetry symm_y;

  t_symmetry symm_z;



  //store an index for each corner,
  //edge, face and for the domain's interior.

  int corner[8];

  int edge[12];

  int face[6];

  int inter;



  size_t level;

  //set the grid size, center and kind of boundaries

  void Set (double dx = 0.01, double dy = 0.01, double dz = 0.01,
	    double cx = 0.0, double cy = 0.0, double cz = 0.0);


  //set the number of points in each direction.
  void Set_N (size_t nx = 10, size_t ny = 10, size_t nz = 10);


  //set the buffer points
  void Set_Buffer (size_t buff_nx = 2, size_t buff_ny = 2, size_t buff_nz =
		   2);







public:


  //constructor (call Make)
    domain (size_t nx = 13, size_t ny = 13, size_t nz = 13,
	    double dx = 0.1, double dy = 0.1, double dz = 0.1,
	    double cx = 0.0, double cy = 0.0, double cz = 0.0,
	    size_t ord = 2,
	    t_symmetry s_x = NONE, t_symmetry s_y = NONE, t_symmetry s_z =
	    NONE)
  {

    Make (nx, ny, nz, dx, dy, dz, cx, cy, cz, ord, s_x, s_y, s_z);


  }


  //destructor 
   ~domain ()
  {
  }


  //constructor 
  void Make (size_t nx, size_t ny, size_t nz,
	     double dx, double dy, double dz,
	     double cx, double cy, double cz,
	     size_t ord, t_symmetry s_x, t_symmetry s_y, t_symmetry s_z);



  void Fix_Iterators (int Corner[8], int Edge[12], int Face[6], int &Inter);



  void Fix_Points ();


  void Get_Partition (size_t & a, size_t & b, size_t & c);

  void Get_Index_Topology (const size_t l, const size_t m, const size_t n,
			   size_t & i, size_t & j, size_t & k, int rk = -1);



  void Set_level (size_t l = 0)
  {
    level = l;
  }

  const size_t Get_level ()
  {
    return (level);
  }

  const size_t Get_i (const double x) const;

  const size_t Get_j (const double x) const;

  const size_t Get_k (const double x) const;


  void Find_cell_xyz (const double x, const double y, const double z,
		      double &xi, double &yi, double &zi,
		      double &xf, double &yf, double &zf);



  const bool Is_1_to_2 (const domain & D) const;


  void Print_to_File (string name);


  const size_t Get_order ()
  {
    return (order);
  }




// Operators:



  friend ostream & operator<< (ostream &, const domain &);


  //homonomia de + para 2
  //dominios cubicos, el
  //resultado es el minimo dominio
  //cubico que cubre a ambos

  friend domain operator+ (const domain &, const domain &);


  //homonomia de - para 2
  //dominios cubicos, el
  //resultado es el minimo dominio
  //cubico contenido
  //entre ambos

  friend domain operator- (const domain &, const domain &);



  friend bool operator == (const domain & D1, const domain & D2);


  friend bool operator != (const domain & D1, const domain & D2)
  {
    return (!(D1 == D2));
  }





  //homonomia de A < B para 2
  //dominios cubicos, el
  //resultado es el valor logico
  //de evaluar si A esta
  //contenido completamente en B

  const bool operator< (const domain & D) const;




  //homonomia de A <= B para 2
  //dominios cubicos, el
  //resultado es el valor logico
  //de evaluar si A esta
  //contenido en B
  const bool operator<= (const domain & D) const;



  //homonomia de A > B para 2
  //dominios cubicos, el
  //resultado es el valor logico
  //de evaluar si B esta
  //contenido completamente en A

  const bool operator> (const domain & D) const;



  //homonomia de A >= B para 2
  //dominios cubicos, el
  //resultado es el valor logico
  //de evaluar si B esta
  //contenido en A

  const bool operator>= (const domain & D) const;


  //homonomia de A && B para 2
  //dominios cubicos, el
  //resultado es el valor logico
  //de evaluar si B intersecta A
  const bool operator&& (const domain & D) const;





  inline const double Get_dx () const
  {
    return (dx);
  }

  inline const double Get_dy () const
  {
    return (dy);
  }

  inline const double Get_dz () const
  {
    return (dz);
  }



  inline const int Get_Bound_Xi () const
  {
    return (bound_xi);
  }

  inline const int Get_Bound_Xf () const
  {
    return (bound_xf);
  }


  inline const int Get_Bound_Yi () const
  {
    return (bound_yi);
  }

  inline const int Get_Bound_Yf () const
  {
    return (bound_yf);
  }


  inline const int Get_Bound_Zi () const
  {
    return (bound_zi);
  }

  inline const int Get_Bound_Zf () const
  {
    return (bound_zf);
  }



  inline const t_symmetry Get_symmetry_x () const
  {
    return (symm_x);
  }

  inline const t_symmetry Get_symmetry_y () const
  {
    return (symm_y);
  }

  inline const t_symmetry Get_symmetry_z () const
  {
    return (symm_z);
  }





  inline size_t Get_N ()
  {
    return (nx[my_rank] * ny[my_rank] * nz[my_rank]);
  }


  inline const size_t Get_Nx () const
  {
    return (nx[my_rank]);
  }

  inline const size_t Get_Ny () const
  {
    return (ny[my_rank]);
  }

  inline const size_t Get_Nz () const
  {
    return (nz[my_rank]);
  }


  inline const size_t Get_Nx (size_t rank) const
  {
    return (nx[rank]);
  }

  inline const size_t Get_Ny (size_t rank) const
  {
    return (ny[rank]);
  }

  inline const size_t Get_Nz (size_t rank) const
  {
    return (nz[rank]);
  }


  inline const size_t Get_buffer_x () const
  {
    return (buff_nx);
  }

  inline const size_t Get_buffer_y () const
  {
    return (buff_ny);
  }

  inline const size_t Get_buffer_z () const
  {
    return (buff_nz);
  }


  inline const size_t Get_Index (const size_t i,
				 const size_t j, const size_t k) const
  {
    return (ijk[INDEXORDER_IJK]);

  };


  inline const size_t Get_Index_I (const size_t index) const
  {
    return (I[index]);
  }

  inline const size_t Get_Index_J (const size_t index) const
  {
    return (J[index]);
  }

  inline const size_t Get_Index_K (const size_t index) const
  {
    return (K[index]);
  }




  inline const double Get_Center_x () const
  {
    return (cx);
  }

  inline const double Get_Center_y () const
  {
    return (cy);
  }

  inline const double Get_Center_z () const
  {
    return (cz);
  }



  inline const double Get_Length_x () const
  {
    return (Lx);
  }

  inline const double Get_Length_y () const
  {
    return (Ly);
  }

  inline const double Get_Length_z () const
  {
    return (Lz);
  }





  inline const double Get_Xi () const
  {
    return (xi[my_rank]);
  }

  inline const double Get_Xf () const
  {
    return (xf[my_rank]);
  }


  inline const double Get_Yi () const
  {
    return (yi[my_rank]);
  }

  inline const double Get_Yf () const
  {
    return (yf[my_rank]);
  }


  inline const double Get_Zi () const
  {
    return (zi[my_rank]);
  }

  inline const double Get_Zf () const
  {
    return (zf[my_rank]);
  }




  inline const double Get_xi () const
  {
    return (xi[my_rank]);
  }

  inline const double Get_xf () const
  {
    return (xf[my_rank]);
  }


  inline const double Get_yi () const
  {
    return (yi[my_rank]);
  }

  inline const double Get_yf () const
  {
    return (yf[my_rank]);
  }


  inline const double Get_zi () const
  {
    return (zi[my_rank]);
  }

  inline const double Get_zf () const
  {
    return (zf[my_rank]);
  }





  inline const double Get_x (const int i) const
  {
    return (xi[my_rank] + double (i) * dx);
  }

  inline const double Get_y (const int i) const
  {
    return (yi[my_rank] + double (i) * dy);
  }

  inline const double Get_z (const int i) const
  {
    return (zi[my_rank] + double (i) * dz);
  }



  inline const double Get_x (const int i, const int j, const int k,
			     const size_t rank) const
  {
    return (xi[rank] + double (i) * dx);
  }

  inline const double Get_y (const int i, const int j, const int k,
			     const size_t rank) const
  {
    return (yi[rank] + double (j) * dy);
  }

  inline const double Get_z (const int i, const int j, const int k,
			     const size_t rank) const
  {
    return (zi[rank] + double (k) * dz);
  }



  inline const int Is_In_rank (const double x, const double y, const double z) const
  {

    int rank = -1;

    for (size_t i = 0; i < size; i++)

      if (Is_In (x, y, z, i))

	  rank = i;

      return (rank);

  }


  inline const bool Is_In (const double x, const double y, const double z,
			   const size_t rank) const
  {
    return (xi[rank] - dx * 0.1 <= x &&
	    x <= xf[rank] + dx * 0.1 &&
	    yi[rank] - dy * 0.1 <= y &&
	    y <= yf[rank] + dy * 0.1 &&
	    zi[rank] - dz * 0.1 <= z && z <= zf[rank] + dz * 0.1);
  }

  inline const bool is_in (const double x, const double y, const double z) const
  {
    return (Is_In (x, y, z, my_rank));
  }


  inline const bool Is_Near (const double x, const double y, const double z) const
  {
    return (xi[my_rank] - 2 * dx <= x && x <= xf[my_rank] + 2 * dx &&
	    yi[my_rank] - 2 * dy <= y && y <= yf[my_rank] + 2 * dy &&
	    zi[my_rank] - 2 * dz <= z && z <= zf[my_rank] + 2 * dz);
  }

  inline const bool Is_In_Boundary (const double x, const double y,
				    const double z) const
  {
    return ((xi[my_rank] - dx <= x && x <= xf[my_rank] + dx &&
	     yi[my_rank] - dy <= y && y <= yf[my_rank] + dy &&
	     zi[my_rank] - dz <= z && z <= zf[my_rank] + dz) &&
	    !(xi[my_rank] + dx <= x && x <= xf[my_rank] - dx &&
	      yi[my_rank] + dy <= y && y <= yf[my_rank] - dy &&
	      zi[my_rank] + dz <= z && z <= zf[my_rank] - dz));
  }


  inline const bool Is_Near (const double x, const double y, const double z,
			     const double Dx, const double Dy,
			     const double Dz) const
  {
    return (xi[my_rank] - Dx <= x && x <= xf[my_rank] + Dx &&
	    yi[my_rank] - Dy <= y && y <= yf[my_rank] + Dy &&
	    zi[my_rank] - Dz <= z && z <= zf[my_rank] + Dz);
  }


  inline const bool x_is_in (const double x) const
  {
    return (xi[my_rank] <= x && x <= xf[my_rank]);
  }

  const bool y_is_in (const double y) const
  {
    return (yi[my_rank] <= y && y <= yf[my_rank]);
  }

  const bool z_is_in (const double z) const
  {
    return (zi[my_rank] <= z && z <= zf[my_rank]);
  }



  const bool Is_In_Nb (const double x, const double y, const double z) const
  {
    return (is_in_nb (x, y, z, my_rank));
  }

  inline const bool is_in_nb (const double x, const double y, const double z) const
  {
    return (is_in_nb (x, y, z, my_rank));
  }

  inline const bool is_in_nb (const double x, const double y, const double z,
			      const size_t rank) const
  {
    return (xi[rank] + dx * (double (buff_nx) - 0.1) <x &&
	    x <
	    xf[rank] - dx * (double (buff_nx) - 0.1) &&yi[rank] +
	    dy * (double (buff_ny) - 0.1) <y
	    && y <
	    yf[rank] - dy * (double (buff_ny) - 0.1) &&zi[rank] +
	    dz * (double (buff_nz) - 0.1) <z
	    && z < zf[rank] - dz * (double (buff_nz) - 0.1));
  }


  inline size_t Get_rank (const int l,
			  const int m,
			  const int n, const int i, const int j, const int k)
  {
    return (n * (m * i + j) + k);
  }




};

//== utility ==//


void ordena (size_t & a, size_t & b, size_t & c);


inline bool
dcomp (const double x, const double y)
{
  return (fabs (x - y) <= 3 * numeric_limits < double >::epsilon ());
}


inline size_t
Min (size_t nx, size_t ny, size_t nz)
{
  size_t n = nx <= ny ? nx : ny;

  n = n <= nz ? n : nz;
  return (n);
}

inline double
Min (double nx, double ny, double nz)
{
  double n = nx <= ny ? nx : ny;

  n = n <= nz ? n : nz;
  return (n);
}

void Set_Boundary (double &ti,
		   double &tf,
		   size_t & Nt,
		   double gti,
		   double gtf,
		   double Lt,
		   double dt, size_t k, size_t c, size_t b_nt, bool symm);


#endif // DOMAIN_H
