



//=========================== Olliptic scalarfield.h =========================//
//  
//  scalarfields. (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//










#ifndef SCALARFIELDS_H

#define SCALARFIELDS_H



//============================ Standar libraries =========================//


#include <string>

#include <map>




//============================ Olliptic libraries =========================//

#include "ffunction.h"

#include "interface.h"



using namespace std;




class scalarfields
{


  domain *D;


  size_t NumSFs;



    map < string, size_t > NameFields;


    protected:vector < ffunction > u;



    public:~scalarfields ()
  {
  };

  scalarfields ()
  {
    NameFields.clear ();
    NumSFs = 0;
    u.clear ();
  };

  scalarfields (domain * D, vector < string > name)
  {
    Make (D, name);
  };

  void Make (domain * D, vector < string > name);


  scalarfields (domain * D, string name)
  {
    Make (D, name);
  };

  void Make (domain * D, string name);



  inline void New_Field (string name)
  {
    NameFields[name] = u.size ();
    u.push_back (ffunction (D));
    NumSFs = u.size ();
  }

  inline ffunction *Get_ffunction (size_t indx = 0)
  {
    if (check_index (indx))
      return (&u[indx]);
    else
      return (NULL);
  }

  inline ffunction *Get_ffunction (string var)
  {
    return (&u[NameFields[var]]);
  }


  inline void Set_Iterator (const t_iterator it = INSIDE, const size_t indx =
			    0)
  {
    if (check_index (indx))
      u[indx].Set_Iterator (it);
  }


  inline void Set_Iterator (const t_iterator it =
			    INSIDE, const size_t indx_down =
			    0, const size_t indx_up = 1)
  {
    for (size_t i = indx_down; i < indx_up; i++)
      if (check_index (i))
	u[i].Set_Iterator (it);
  }


  inline void Set_Iterator (const t_iterator it, string name)
  {
    u[NameFields[name]].Set_Iterator (it);
  }


  inline void Set_Iterator (const t_iterator it, vector < string > name)
  {
    for (size_t i = 0; i < name.size (); i++)
      u[NameFields[name[i]]].Set_Iterator (it);
  }


  inline void Set_val (const double f, vector < string > name)
  {
    for (size_t i = 0; i < name.size (); i++)
      u[NameFields[name[i]]].Set (f);
  }

  inline void Set_val (const double f, string name)
  {
    u[NameFields[name]].Set (f);
  }


  inline bool End (size_t indx)
  {
    if (check_index (indx))
      return (u[indx].End ());
    else
      return (true);
  }

  inline bool End (size_t indx_down, size_t indx_up)
  {
    bool end = true;

    for (size_t i = indx_down; i < indx_up; i++)
      if (check_index (i))
	end = (end && u[i].End ());

    return (end);
  }




  inline bool End ()
  {
    bool end = true;

    for (size_t i = 0; i < NumSFs; i++)
      end = (end && u[i].End ());

    return (end);
  }


  inline bool End (string var)
  {
    return (u[NameFields[var]].End ());
  }


  inline bool End (vector < string > var)
  {
    bool end = true;

    for (size_t i = 0; i < var.size (); i++)

      end = (end && u[NameFields[var[i]]].End ());


    return (end);
  }







  inline bool check_index (size_t indx)
  {
    return (indx >= 0 && indx < NumSFs);
  }







  inline const double Get_x (size_t indx) const
  {
    return (u[indx].Get_x ());
  }

  inline const double Get_dx (size_t indx) const
  {
    return (u[indx].Get_dx ());
  }


  inline const double Get_y (size_t indx) const
  {
    return (u[indx].Get_y ());
  }

  inline const double Get_dy (size_t indx) const
  {
    return (u[indx].Get_dy ());
  }


  inline const double Get_z (size_t indx) const
  {
    return (u[indx].Get_z ());
  }

  inline const double Get_dz (size_t indx) const
  {
    return (u[indx].Get_dz ());
  }


  inline const double Get_val (size_t indx) const
  {
    return (u[indx].Get_val ());
  }





  inline const double Get_x (size_t indx, size_t i) const
  {
    return (u[indx].Get_x (i));
  }

  inline const double Get_y (size_t indx, size_t i) const
  {
    return (u[indx].Get_y (i));
  }

  inline const double Get_z (size_t indx, size_t i) const
  {
    return (u[indx].Get_z (i));
  }

  inline const double Get_val (size_t indx, size_t i) const
  {
    return (u[indx].Get_val (i));
  }




  inline const double Get_x (string var)
  {
    return (u[NameFields[var]].Get_x ());
  }

  inline const double Get_dx (string var)
  {
    return (u[NameFields[var]].Get_dx ());
  }


  inline const double Get_y (string var)
  {
    return (u[NameFields[var]].Get_y ());
  }

  inline const double Get_dy (string var)
  {
    return (u[NameFields[var]].Get_dy ());
  }


  inline const double Get_z (string var)
  {
    return (u[NameFields[var]].Get_z ());
  }

  inline const double Get_dz (string var)
  {
    return (u[NameFields[var]].Get_dz ());
  }


  inline const double Get_val (string var)
  {
    return (u[NameFields[var]].Get_val ());
  }





  inline const double Get_x (string var, size_t i)
  {
    return (u[NameFields[var]].Get_x (i));
  }

  inline const double Get_y (string var, size_t i)
  {
    return (u[NameFields[var]].Get_y (i));
  }

  inline const double Get_z (string var, size_t i)
  {
    return (u[NameFields[var]].Get_z (i));
  }

  inline const double Get_val (string var, size_t i)
  {
    return (u[NameFields[var]].Get_val (i));
  }




  inline const int Get_Index (size_t indx) const
  {
    return (u[indx].Get_Index ());
  }

  inline const int Get_Index (string var)
  {
    return (u[NameFields[var]].Get_Index ());
  }



  inline const size_t Get_Num_Fields () const
  {
    return (NumSFs);
  }


  //const string Get_Name_Field(size_t i) const { return(NumSFs); }





  inline double Get_Norm_LInf (size_t l)
  {
    return (Get_ffunction (l)->Get_Norm_LInf ());
  }

  inline double Get_Norm_LInf (string var)
  {
    return (Get_ffunction (var)->Get_Norm_LInf ());
  }


  inline double Lap (size_t l)
  {
    return (Get_ffunction (l)->Lap ());
  }

  inline double duLap (size_t l)
  {
    return (Get_ffunction (l)->duLap ());
  }


  inline double Dx (size_t l)
  {
    return (Get_ffunction (l)->Dx ());
  }

  inline double duDx (size_t l)
  {
    return (Get_ffunction (l)->duDx ());
  }


  inline double Dy (size_t l)
  {
    return (Get_ffunction (l)->Dy ());
  }

  inline double duDy (size_t l)
  {
    return (Get_ffunction (l)->duDy ());
  }


  inline double Dz (size_t l)
  {
    return (Get_ffunction (l)->Dz ());
  }

  inline double duDz (size_t l)
  {
    return (Get_ffunction (l)->duDz ());
  }




  inline double DDxy (size_t l)
  {
    return (Get_ffunction (l)->DDxy ());
  }

  inline double duDDxy (size_t l)
  {
    return (Get_ffunction (l)->duDDxy ());
  }


  inline double DDyz (size_t l)
  {
    return (Get_ffunction (l)->DDyz ());
  }

  inline double duDDyz (size_t l)
  {
    return (Get_ffunction (l)->duDDyz ());
  }


  inline double DDxz (size_t l)
  {
    return (Get_ffunction (l)->DDxz ());
  }

  inline double duDDxz (size_t l)
  {
    return (Get_ffunction (l)->duDDxz ());
  }



  inline double DDx (size_t l)
  {
    return (Get_ffunction (l)->DDx ());
  }

  inline double duDDx (size_t l)
  {
    return (Get_ffunction (l)->duDDx ());
  }


  inline double DDy (size_t l)
  {
    return (Get_ffunction (l)->DDy ());
  }

  inline double duDDy (size_t l)
  {
    return (Get_ffunction (l)->duDDy ());
  }


  inline double DDz (size_t l)
  {
    return (Get_ffunction (l)->DDz ());
  }

  inline double duDDz (size_t l)
  {
    return (Get_ffunction (l)->duDDz ());
  }





  inline double Lap (string var)
  {
    return (Get_ffunction (var)->Lap ());
  }

  inline double duLap (string var)
  {
    return (Get_ffunction (var)->duLap ());
  }


  inline double Dx (string var)
  {
    return (Get_ffunction (var)->Dx ());
  }

  inline double duDx (string var)
  {
    return (Get_ffunction (var)->duDx ());
  }


  inline double Dy (string var)
  {
    return (Get_ffunction (var)->Dy ());
  }

  inline double duDy (string var)
  {
    return (Get_ffunction (var)->duDy ());
  }


  inline double Dz (string var)
  {
    return (Get_ffunction (var)->Dz ());
  }

  inline double duDz (string var)
  {
    return (Get_ffunction (var)->duDz ());
  }



  inline double DDxy (string var)
  {
    return (Get_ffunction (var)->DDxy ());
  }

  inline double duDDxy (string var)
  {
    return (Get_ffunction (var)->duDDxy ());
  }


  inline double DDyz (string var)
  {
    return (Get_ffunction (var)->DDyz ());
  }

  inline double duDDyz (string var)
  {
    return (Get_ffunction (var)->duDDyz ());
  }


  inline double DDxz (string var)
  {
    return (Get_ffunction (var)->DDxz ());
  }

  inline double duDDxz (string var)
  {
    return (Get_ffunction (var)->duDDxz ());
  }



  inline double DDx (string var)
  {
    return (Get_ffunction (var)->DDx ());
  }

  inline double duDDx (string var)
  {
    return (Get_ffunction (var)->duDDx ());
  }


  inline double DDy (string var)
  {
    return (Get_ffunction (var)->DDy ());
  }

  inline double duDDy (string var)
  {
    return (Get_ffunction (var)->duDDy ());
  }


  inline double DDz (string var)
  {
    return (Get_ffunction (var)->DDz ());
  }

  inline double duDDz (string var)
  {
    return (Get_ffunction (var)->duDDz ());
  }




  inline double max (size_t l)
  {
    return (u[l].max ());
  }

  inline double min (size_t l)
  {
    return (u[l].min ());
  }


  inline double max (string var)
  {
    return (u[NameFields[var]].max ());
  }

  inline double min (string var)
  {
    return (u[NameFields[var]].min ());
  }



  inline void Print_3D (string dir, string name, string postfix,
			double (*func_to_comp) (double x, double y,
						double z) = NULL)
  {
    u[NameFields[name]].Print_3D (dir + name + postfix, func_to_comp);
  }


  inline void Print_2D (string dir, string name, string postfix,
			const cut c = CUT_X,
			const double x1 = 0,
			bool interpolar = true,
			double (*func_to_comp) (double x, double y,
						double z) =
			NULL, m_interpol minter = LAGRANGE, size_t order = 2)
  {
    u[NameFields[name]].Print_2D (dir + name + postfix, c, x1, interpolar,
				  func_to_comp, minter, order);
  }





  inline void Print_1D (string dir, string name, string postfix,
			const cut c = CUT_X,
			const double x1 = 0,
			const double x2 = 0,
			bool interpolar = true,
			double (*func_to_comp) (double x, double y,
						double z) =
			NULL, m_interpol minter = LAGRANGE, size_t order = 2)
  {
    u[NameFields[name]].Print_1D (dir + name + postfix, c, x1, x2, interpolar,
				  func_to_comp, minter, order);
  }





  inline void Print_t1D (string dir, string name, string postfix,
			 const double t)
  {
    u[NameFields[name]].Print_t1D (dir + name + postfix, t);
  }


  inline void Print_t2D (string dir, string name, string postfix, size_t cyc)
  {
    u[NameFields[name]].Print_t2D (dir + name + postfix, cyc);
  }






};




#endif
