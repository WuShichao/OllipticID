



//=========================== Olliptic elliptic.h =========================//
//  
//  Multigrids. (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//










#ifndef MULTIGRID_H

#define MULTIGRID_H



//============================ Standar libraries =========================//


#include <time.h>


//============================ Olliptic libraries =========================//

#include "scalarfields.h"

#include "interface.h"



using namespace std;




class multigrid
{



  vector < size_t > num_dom_in_level;


  vector < domain > D;


  size_t levels;

  size_t l0;

  size_t total_grids;


  t_domain dom;

  size_t l_it;

  size_t l_end;



    protected:vector < scalarfields > u;


    public:~multigrid ()
  {
  };

  multigrid ()
  {
    num_dom_in_level.clear ();
    u.clear ();
  };


  multigrid (vector < domain > D, t_domain d, vector < string > name)
  {
    Make (D, d, name);
  };

  void Make (vector < domain > D, t_domain d, vector < string > name);



  inline ffunction *Get_ffunction (size_t level, size_t ibox, string var_name)
  {
    return (u[Get_dom_index (level, ibox)].Get_ffunction (var_name));
  }




  inline void Set_Iterator (vector < string > var_name, const t_iterator it =
			    INSIDE, const size_t low = 0, const int up = -1)
  {
    l_it = low < levels ? Get_dom_index (low) : Get_dom_index (levels - 1);

    l_end = up < int (levels)
      && up >= 0 ? Get_dom_index (up,
				  num_dom_in_level[up] - 1) +
      1 : Get_dom_index (levels - 1, num_dom_in_level[levels - 1] - 1) + 1;

    for (size_t ibox = l_it; ibox < l_end; ibox++)

      u[ibox].Set_Iterator (it, var_name);

  }



  inline void Set_Iterator_ibox (vector < string > var_name,
				 const t_iterator it, const size_t l,
				 const size_t ibox)
  {


    l_it = Get_dom_index (l, ibox);

    u[l_it].Set_Iterator (it, var_name);

    l_end = l_it + 1;

  }



  inline bool End (vector < string > var_name)
  {

    if (!u[l_it].End (var_name))

      l_it++;

    return (l_it < l_end);
  }



  inline bool End ()
  {

    if (!u[l_it].End ())

      l_it++;

    return (l_it < l_end);
  }



  inline void Set_Iterator (string var_name, const t_iterator it =
			    INSIDE, const size_t low = 0, const int up = -1)
  {
    l_it = low < levels ? Get_dom_index (low) : Get_dom_index (levels - 1);

    l_end = up < int (levels)
      && up >= 0 ? Get_dom_index (up,
				  num_dom_in_level[up] - 1) +
      1 : Get_dom_index (levels - 1, num_dom_in_level[levels - 1] - 1) + 1;

    for (size_t ibox = l_it; ibox < l_end; ibox++)

      u[ibox].Set_Iterator (it, var_name);

  }


  inline bool End (string var_name)
  {
    if (!u[l_it].End (var_name))
      l_it++;
    return (l_it < l_end);
  }






  inline double Get_x (string var_name)
  {
    return (u[l_it].Get_x (var_name));
  }

  inline double Get_y (string var_name)
  {
    return (u[l_it].Get_y (var_name));
  }

  inline double Get_z (string var_name)
  {
    return (u[l_it].Get_z (var_name));
  }


  inline double Get_dx (string var_name)
  {
    return (u[l_it].Get_dx (var_name));
  }

  inline double Get_dy (string var_name)
  {
    return (u[l_it].Get_dy (var_name));
  }

  inline double Get_dz (string var_name)
  {
    return (u[l_it].Get_dz (var_name));
  }


  inline double Get_val (string var_name)
  {
    return (u[l_it].Get_val (var_name));
  }

  inline void Set_val (double f, string var_name)
  {
    return (u[l_it].Set_val (f, var_name));
  }

  inline void Set_val (double f, vector < string > var_name)
  {
    return (u[l_it].Set_val (f, var_name));
  }


  inline const size_t Get_levels () const
  {
    return (levels);
  }

  inline const size_t Get_num_dom (size_t l) const
  {
    return (num_dom_in_level[l]);
  }

  size_t Get_dom_index (const size_t lev, const size_t boxi = 0);


  void Transferir (size_t from, string from_var, size_t to, string to_var,
		   t_iterator iter = ALL, m_interpol method =
		   LINEAR, size_t ord = 1, bool all = true, bool to_zero =
		   false);




  inline void Print_3D (string dir, string name,
			size_t l,
			double (*func_to_comp) (double x, double y,
						double z) = NULL)
  {
    l = l < levels ? l : levels - 1;

    for (size_t ibox = 0; ibox < num_dom_in_level[l]; ibox++)
      {

	stringstream level_name;

	level_name << "-l" << l << "b" << ibox << "-";

	u[Get_dom_index (l, ibox)].Print_3D (dir, name, level_name.str (),
					     func_to_comp);
      }
  }





  inline void Print_2D (string dir, string name,
			size_t l,
			const cut c = CUT_Z,
			const double x1 = 0,
			bool interpolar = false,
			double (*func_to_comp) (double x, double y,
						double z) =
			NULL, m_interpol minter = LAGRANGE, size_t order = 2)
  {
    l = l < levels ? l : levels - 1;

    for (size_t ibox = 0; ibox < num_dom_in_level[l]; ibox++)
      {

	stringstream level_name;

	level_name << "-l" << l << "b" << ibox << "-";

	u[Get_dom_index (l, ibox)].Print_2D (dir, name, level_name.str (), c,
					     x1, interpolar, func_to_comp,
					     minter, order);
      }
  }




  inline void Print_1D (string dir, string name,
			size_t l,
			const cut c = CUT_X,
			const double x1 = 0,
			const double x2 = 0,
			bool interpolar = false,
			double (*func_to_comp) (double x, double y,
						double z) =
			NULL, m_interpol minter = LAGRANGE, size_t order = 2)
  {

    l = l < levels ? l : levels - 1;

    for (size_t ibox = 0; ibox < num_dom_in_level[l]; ibox++)
      {

	stringstream level_name;

	level_name << "-l" << l << "b" << ibox << "-";

	u[Get_dom_index (l, ibox)].Print_1D (dir, name, level_name.str (), c,
					     x1, x2, interpolar, func_to_comp,
					     minter, order);
      }
  }



  void Print_Interpol_1D (string dir, string name,
			  double LD,
			  double _dx = -1,
			  const cut c = CUT_X,
			  const double x1 = 0,
			  const double x2 = 0,
			  double (*func_to_comp) (double x, double y,
						  double z) =
			  NULL, m_interpol minter = LAGRANGE, size_t order =
			  9);

  void Print_Interpol_2D (string dir, string name,
			  double LD1,
			  double LD2,
			  double _dx1 = -1,
			  double _dx2 = -1,
			  const cut c = CUT_Z,
			  const double x1 = 0,
			  double (*funct_to_comp) (double x, double y,
						   double z) =
			  NULL, m_interpol minter = LAGRANGE, size_t order =
			  9);



  void Print_IO_BAM (string name, string var_name, m_interpol minter =
		     LAGRANGE, size_t order = 9);

  void Print_IO_Zcode (string name, string var_name, m_interpol minter =
		       LAGRANGE, size_t order = 9);


  inline void Print_t1D (string dir, string name, size_t l, const double t)
  {

    l = l < levels ? l : levels - 1;

    for (size_t ibox = 0; ibox < num_dom_in_level[l]; ibox++)
      {
	stringstream level_name;

	level_name << "-l" << l << "b" << ibox;

	u[Get_dom_index (l, ibox)].Print_t1D (dir, name, level_name.str (),
					      t);
      }
  }


  inline void Print_t2D (string dir, string name, size_t l, size_t cyc)
  {

    l = l < levels ? l : levels - 1;

    for (size_t ibox = 0; ibox < num_dom_in_level[l]; ibox++)
      {
	stringstream level_name;

	level_name << "-l" << l << "b" << ibox;

	u[Get_dom_index (l, ibox)].Print_t2D (dir, name, level_name.str (),
					      cyc);

      }

  }



  double Get_Norm_LInf (size_t lev, string var_name);



  inline double Lap (size_t l, size_t ibox, string var)
  {
    return (u[Get_dom_index (l, ibox)].Lap (var));
  }

  inline double duLap (size_t l, size_t ibox, string var)
  {
    return (u[Get_dom_index (l, ibox)].duLap (var));
  }

  inline double Get_val (size_t l, size_t ibox, string var)
  {
    return (u[Get_dom_index (l, ibox)].Get_val (var));
  }



  inline double Lap (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].Lap (iv));
  }

  inline double duLap (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duLap (iv));
  }

  inline double Get_val (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].Get_val (iv));
  }


  inline double Dx (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].Dx (iv));
  }

  inline double duDx (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDx (iv));
  }


  inline double Dy (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].Dy (iv));
  }

  inline double duDy (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDy (iv));
  }


  inline double Dz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].Dz (iv));
  }

  inline double duDz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDz (iv));
  }




  inline double DDxy (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].DDxy (iv));
  }

  inline double duDDxy (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDDxy (iv));
  }


  inline double DDyz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].DDyz (iv));
  }

  inline double duDDyz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDDyz (iv));
  }


  inline double DDxz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].DDxz (iv));
  }

  inline double duDDxz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDDxz (iv));
  }





  inline double DDx (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].DDx (iv));
  }

  inline double duDDx (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDDx (iv));
  }


  inline double DDy (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].DDy (iv));
  }

  inline double duDDy (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDDy (iv));
  }


  inline double DDz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].DDz (iv));
  }

  inline double duDDz (size_t l, size_t ibox, size_t iv)
  {
    return (u[Get_dom_index (l, ibox)].duDDz (iv));
  }



};




#endif
