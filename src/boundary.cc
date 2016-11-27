#include "elliptic.h"





void
elliptic::Set_Boundary ()
{




  if (source_names.size () > 0)
    {

      vars.Set_Iterator (source_names, BOUNDARY);

      do
	{

	  double x = vars.Get_x (source_names[0]);

	  double y = vars.Get_y (source_names[0]);

	  double z = vars.Get_z (source_names[0]);


	  for (size_t iv = 0; iv < var_names.size (); iv++)
	    {


	      if (boundary == ROBIN)

		vars.Set_val (boundary_coeff[iv], source_names);

	      else
		{
		  stringstream stage_var;

		  stage_var << iv;

		  vars.Set_val ((this->*Boundary[iv]) (x, y, z),
				"rho" + stage_var.str ());


		}
	    }



	}
      while (vars.End (source_names));
    }



  for (size_t i = 0; i < var_names.size (); i++)
    {

      vars.Set_Iterator (var_names[i], ALL);

      do
	{


	  double x = vars.Get_x (var_names[i]);

	  double y = vars.Get_y (var_names[i]);

	  double z = vars.Get_z (var_names[i]);


	  vars.Set_val ((this->*Boundary[i]) (x, y, z), var_names[i]);
	  

	}
      while (vars.End (var_names[i]));

    }



}







void
elliptic::App_Symmetry (ffunction * U)
{




  if (U->Get_domain ()->Get_Bound_Xf () == -2)
    {


      int i,
        j,
        k;

      int buf = U->Get_domain ()->Get_buffer_x ();

      int nx = U->Get_domain ()->Get_Nx ();

      t_iterator bound[] = { FACE_XF,

	CORNER_XFYIZI,
	CORNER_XFYIZF,

	CORNER_XFYFZI,
	CORNER_XFYFZF,


	EDGE_XFYI,
	EDGE_XFYF,

	EDGE_XFZI,
	EDGE_XFZF
      };


      for (size_t ib = 0; ib < 9; ib++)
	{

	  U->Set_Iterator (bound[ib]);

	  do
	    {


	      U->Get_MIndex (i, j, k);

	      *U = U->Get_val (2 * (nx - 1 - buf) - i, j, k);


	    }
	  while (U->End ());

	}


    }


  else
   if (U->Get_domain ()->Get_Bound_Xf () == -3)
    {


      int i,
        j,
        k;

      int buf = U->Get_domain ()->Get_buffer_x ();

      int nx = U->Get_domain ()->Get_Nx ();

      int ny = U->Get_domain ()->Get_Nx ();

      t_iterator bound[] = { FACE_XF,

	CORNER_XFYIZI,
	CORNER_XFYIZF,

	CORNER_XFYFZI,
	CORNER_XFYFZF,


	EDGE_XFYI,
	EDGE_XFYF,

	EDGE_XFZI,
	EDGE_XFZF
      };


      for (size_t ib = 0; ib < 9; ib++)
	{

	  U->Set_Iterator (bound[ib]);

	  do
	    {


	      U->Get_MIndex (i, j, k);

	      *U = U->Get_val (2 * (nx - 1 - buf) - i, ny - 1 - j, k);


	    }
	  while (U->End ());

	}


    }








  if (U->Get_domain ()->Get_Bound_Yf () == -2)
    {


      int i,
        j,
        k;

      int buf = U->Get_domain ()->Get_buffer_y ();

      int ny = U->Get_domain ()->Get_Ny ();


      t_iterator bound[] = { FACE_YF,

	CORNER_XIYFZI,
	CORNER_XFYFZI,

	CORNER_XIYFZF,
	CORNER_XFYFZF,


	EDGE_XIYF,
	EDGE_XFYF,

	EDGE_YFZI,
	EDGE_YFZF
      };


      for (size_t ib = 0; ib < 9; ib++)
	{

	  U->Set_Iterator (bound[ib]);

	  do
	    {


	      U->Get_MIndex (i, j, k);

	      *U = U->Get_val (i, 2 * (ny - 1 - buf) - j, k);


	    }
	  while (U->End ());

	}


    }


  else

   if (U->Get_domain ()->Get_Bound_Yf () == -3)
    {


      int i,
        j,
        k;

      int buf = U->Get_domain ()->Get_buffer_y ();

      int ny = U->Get_domain ()->Get_Ny ();

      int nz = U->Get_domain ()->Get_Nz ();


      t_iterator bound[] = { FACE_YF,

	CORNER_XIYFZI,
	CORNER_XFYFZI,

	CORNER_XIYFZF,
	CORNER_XFYFZF,


	EDGE_XIYF,
	EDGE_XFYF,

	EDGE_YFZI,
	EDGE_YFZF
      };


      for (size_t ib = 0; ib < 9; ib++)
	{

	  U->Set_Iterator (bound[ib]);

	  do
	    {


	      U->Get_MIndex (i, j, k);

	      *U = U->Get_val (i, 2 * (ny - 1 - buf) - j, nz - 1 - k);


	    }
	  while (U->End ());

	}


    }








  if (U->Get_domain ()->Get_Bound_Zf () == -2)
    {


      int i,
        j,
        k;

      int buf = U->Get_domain ()->Get_buffer_z ();

      int nz = U->Get_domain ()->Get_Nz ();

      t_iterator bound[] = { FACE_ZF,

	CORNER_XIYIZF,
	CORNER_XFYIZF,

	CORNER_XIYFZF,
	CORNER_XFYFZF,


	EDGE_XIZF,
	EDGE_XFZF,

	EDGE_YIZF,
	EDGE_YFZF
      };


      for (size_t ib = 0; ib < 9; ib++)
	{

	  U->Set_Iterator (bound[ib]);

	  do
	    {


	      U->Get_MIndex (i, j, k);

	      *U = U->Get_val (i, j, 2 * (nz - 1 - buf) - k);


	    }
	  while (U->End ());

	}


    }


  else
   if (U->Get_domain ()->Get_Bound_Zf () == -3)
    {


      int i,
        j,
        k;

      int buf = U->Get_domain ()->Get_buffer_z ();

      int nx = U->Get_domain ()->Get_Nx ();

      int nz = U->Get_domain ()->Get_Nz ();

      t_iterator bound[] = { FACE_ZF,

	CORNER_XIYIZF,
	CORNER_XFYIZF,

	CORNER_XIYFZF,
	CORNER_XFYFZF,


	EDGE_XIZF,
	EDGE_XFZF,

	EDGE_YIZF,
	EDGE_YFZF
      };


      for (size_t ib = 0; ib < 9; ib++)
	{

	  U->Set_Iterator (bound[ib]);

	  do
	    {


	      U->Get_MIndex (i, j, k);

	      *U = U->Get_val (nx - 1 - i, j, 2 * (nz - 1 - buf) - k);


	    }
	  while (U->End ());

	}


    }




  U->Sync ();



}





void elliptic::App_Bound (size_t ll, size_t ibox)
{


  for (size_t iv = 0; iv < var_names.size (); iv++)
    {

      ffunction *U = vars.Get_ffunction (ll, ibox, var_names[iv]);

      ffunction *S = vars.Get_ffunction (ll, ibox, source_names[iv]);

      size_t pn = falloff_n[iv];

      if (boundary == ROBIN)

	{



	  if (U->Get_domain ()->Get_Bound_Xi () == -1)

	    Boundary_Face (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X, U, FACE_XI, pn, S);



	  if (U->Get_domain ()->Get_Bound_Xf () == -1)

	    Boundary_Face (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X, U, FACE_XF, pn, S);



	  if (U->Get_domain ()->Get_Bound_Yi () == -1)

	    Boundary_Face (&ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y, U, FACE_YI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Yf () == -1)

	    Boundary_Face (&ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y, U, FACE_YF, pn, S);





	  if (U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Face (&ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, FACE_ZI, pn, S);



	  if (U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Face (&ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, FACE_ZF, pn, S);





	  //  XY //

	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y, U, EDGE_XIYI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y, U, EDGE_XIYF, pn, S);




	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y, U, EDGE_XFYI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y, U, EDGE_XFYF, pn, S);




	  //  XZ //

	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_XIZI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_XIZF, pn, S);




	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Edge (&ffunction::Dx,
			   &ffunction::duDx,
			   &ffunction::Get_X,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_XFZI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)
	    {

	      Boundary_Edge (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, EDGE_XFZF, pn, S);


	    }



	  //  YZ //

	  if (U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Edge (&ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_YIZI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Edge (&ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_YIZF, pn, S);




	  if (U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Edge (&ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_YFZI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Edge (&ffunction::Dy,
			   &ffunction::duDy,
			   &ffunction::Get_Y,
			   &ffunction::Dz,
			   &ffunction::duDz,
			   &ffunction::Get_Z, U, EDGE_YFZF, pn, S);




	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XIYIZI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XIYIZF, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XFYIZI, pn, S);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XFYIZF, pn, S);





	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XIYFZI, pn, S);




	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XIYFZF, pn, S);




	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XFYFZI, pn, S);





	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    Boundary_Corner (&ffunction::Dx,
			     &ffunction::duDx,
			     &ffunction::Get_X,
			     &ffunction::Dy,
			     &ffunction::duDy,
			     &ffunction::Get_Y,
			     &ffunction::Dz,
			     &ffunction::duDz,
			     &ffunction::Get_Z, U, CORNER_XFYFZF, pn, S);


	}
      else
       if (boundary == DIRICHLET)

	{



	  U->Set_Iterator (BOUNDARY);
	  S->Set_Iterator (BOUNDARY);


	  do
	    {


	      *U = S->Get_val ();

	      S->Next ();


	    }
	  while (U->End ());




	}

      else
       if (boundary == ASYMPTOTIC)

	{


	  size_t i,j,k;

	  double x,y,z, r, coeffB;



	  i = U->Get_nx () / 2;

	  j = U->Get_ny () / 2;

	  k = U->Get_order () / 2;

	  x = U->Get_x (i, j, k);

	  y = U->Get_y (i, j, k);

	  z = U->Get_z (i, j, k);

	  r = sqrt (x * x + y * y + z * z);

	  coeffB = (U->Get_val (i, j, k) - S->Get_val (i, j, k)) * pow (r,double (pn));



	  U->Set_Iterator (BOUNDARY);
	  S->Set_Iterator (BOUNDARY);

	  do
	    {


	      x = U->Get_x ();

	      y = U->Get_y ();

	      z = U->Get_z ();


	      *U = coeffB / sqrt (x * x + y * y + z * z) + S->Get_val ();


	      S->Next ();



	    }
	  while (U->End ());










	}


      else
	{


	  cout << "Error: wrong boundary set" << endl;

	  exit (1);


	}
    }

}






void
Boundary_Face (double (ffunction::*dxu) (void),
	       double (ffunction::*du) (void),
	       double (ffunction::*cd) (void),
	       ffunction * U, t_iterator face, size_t pn, ffunction * S)
{




  double ipn = 1.0 / double (pn);



  U->Set_Iterator (face);

  S->Set_Iterator (face);

  do
    {


      double x2 = U->Get_x () * U->Get_x ();

      double y2 = U->Get_y () * U->Get_y ();

      double z2 = U->Get_z () * U->Get_z ();


      double r2 = x2 + y2 + z2;



      double F =
	r2 * ipn * (U->*dxu) () / (U->*cd) () + U->Get_val () - S->Get_val ();

      double dF = r2 * ipn * (U->*du) () / (U->*cd) () + 1.0;


      U->Set (U->Get_val () - F / dF);



      S->Next ();


    }
  while (U->End ());



}









void
Boundary_Edge (double (ffunction::*dxu1) (void),
	       double (ffunction::*du1) (void),
	       double (ffunction::*cd1) (void),
	       double (ffunction::*dxu2) (void),
	       double (ffunction::*du2) (void),
	       double (ffunction::*cd2) (void),
	       ffunction * U, t_iterator face, size_t pn, ffunction * S)
{




  double ipn = 1.0 / double (pn);




  U->Set_Iterator (face);

  S->Set_Iterator (face);

  do
    {


      double x2 = U->Get_x () * U->Get_x ();

      double y2 = U->Get_y () * U->Get_y ();

      double z2 = U->Get_z () * U->Get_z ();


      double r2 = x2 + y2 + z2;


      double F = r2 * ipn * ((U->*dxu1) () / (U->*cd1) () +
			     (U->*dxu2) () / (U->*cd2) ()) + U->Get_val () -
	S->Get_val ();

      double dF = r2 * ipn * ((U->*du1) () / (U->*cd1) () +
			      (U->*du2) () / (U->*cd2) ()) + 1.0;




      U->Set (U->Get_val () - F / dF);




      S->Next ();

    }
  while (U->End ());



}







void
Boundary_Corner (double (ffunction::*dxu1) (void),
		 double (ffunction::*du1) (void),
		 double (ffunction::*cd1) (void),
		 double (ffunction::*dxu2) (void),
		 double (ffunction::*du2) (void),
		 double (ffunction::*cd2) (void),
		 double (ffunction::*dxu3) (void),
		 double (ffunction::*du3) (void),
		 double (ffunction::*cd3) (void),
		 ffunction * U, t_iterator face, size_t pn, ffunction * S)
{



  double ipn = 1.0 / double (pn);



  U->Set_Iterator (face);

  S->Set_Iterator (face);


  do
    {


      double x2 = U->Get_x () * U->Get_x ();

      double y2 = U->Get_y () * U->Get_y ();

      double z2 = U->Get_z () * U->Get_z ();


      double r2 = x2 + y2 + z2;


      double F = r2 * ipn * ((U->*dxu1) () / (U->*cd1) () +
			     (U->*dxu2) () / (U->*cd2) () +
			     (U->*dxu3) () / (U->*cd3) ()) + U->Get_val () -
	S->Get_val ();

      double dF = r2 * ipn * ((U->*du1) () / (U->*cd1) () +
			      (U->*du2) () / (U->*cd2) () +
			      (U->*du3) () / (U->*cd3) ()) + 1.0;

      U->Set (U->Get_val () - F / dF);



      S->Next ();
    }
  while (U->End ());






}








//==  Solo se aplica el operador no se ==//
//==  actualiza el valor de la funcion ==//





void
elliptic::App_Boundary_Operator (size_t ll, size_t ibox)
{



  for (size_t iv = 0; iv < var_names.size (); iv++)
    {

      ffunction *U = vars.Get_ffunction (ll, ibox, var_names[iv]);

      stringstream stage_var;

      stage_var << iv;

      ffunction *R = vars.Get_ffunction (ll, ibox, "r" + stage_var.str ());


      size_t pn = falloff_n[iv];



      if (boundary == ROBIN)

	{



	  if (U->Get_domain ()->Get_Bound_Xi () == -1)

	    App_Boundary_Face (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X, U, FACE_XI, pn, R);



	  if (U->Get_domain ()->Get_Bound_Xf () == -1)

	    App_Boundary_Face (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X, U, FACE_XF, pn, R);



	  if (U->Get_domain ()->Get_Bound_Yi () == -1)

	    App_Boundary_Face (&ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y, U, FACE_YI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Yf () == -1)

	    App_Boundary_Face (&ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y, U, FACE_YF, pn, R);





	  if (U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Face (&ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, FACE_ZI, pn, R);



	  if (U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Face (&ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, FACE_ZF, pn, R);





	  //  XY //

	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y, U, EDGE_XIYI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y, U, EDGE_XIYF, pn, R);




	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y, U, EDGE_XFYI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y, U, EDGE_XFYF, pn, R);




	  //  XZ //

	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_XIZI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_XIZF, pn, R);




	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Edge (&ffunction::Dx,
			       &ffunction::duDx,
			       &ffunction::Get_X,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_XFZI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)
	    {

	      App_Boundary_Edge (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, EDGE_XFZF, pn, R);


	    }



	  //  YZ //

	  if (U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Edge (&ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_YIZI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Edge (&ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_YIZF, pn, R);




	  if (U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Edge (&ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_YFZI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Edge (&ffunction::Dy,
			       &ffunction::duDy,
			       &ffunction::Get_Y,
			       &ffunction::Dz,
			       &ffunction::duDz,
			       &ffunction::Get_Z, U, EDGE_YFZF, pn, R);




	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XIYIZI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XIYIZF, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XFYIZI, pn, R);


	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yi () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XFYIZF, pn, R);





	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XIYFZI, pn, R);




	  if (U->Get_domain ()->Get_Bound_Xi () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XIYFZF, pn, R);




	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zi () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XFYFZI, pn, R);





	  if (U->Get_domain ()->Get_Bound_Xf () == -1 &&
	      U->Get_domain ()->Get_Bound_Yf () == -1 &&
	      U->Get_domain ()->Get_Bound_Zf () == -1)

	    App_Boundary_Corner (&ffunction::Dx,
				 &ffunction::duDx,
				 &ffunction::Get_X,
				 &ffunction::Dy,
				 &ffunction::duDy,
				 &ffunction::Get_Y,
				 &ffunction::Dz,
				 &ffunction::duDz,
				 &ffunction::Get_Z, U, CORNER_XFYFZF, pn, R);


	}
      else
       if (boundary == DIRICHLET || boundary == ASYMPTOTIC)

	{



	  U->Set_Iterator (BOUNDARY);
	  R->Set_Iterator (BOUNDARY);

	  do
	    {



	      *R = U->Get_val ();




	      R->Next ();



	    }
	  while (U->End ());




	}



      else
	{


	  cout << "Error: wrong boundary set" << endl;

	  exit (1);


	}



      App_Symmetry (R);



    }



}







void
App_Boundary_Face (double (ffunction::*dxu) (void),
		   double (ffunction::*du) (void),
		   double (ffunction::*cd) (void),
		   ffunction * U, t_iterator face, size_t pn, ffunction * R)
{




  double ipn = 1.0 / double (pn);



  U->Set_Iterator (face);

  R->Set_Iterator (face);

  do
    {


      double x2 = U->Get_x () * U->Get_x ();

      double y2 = U->Get_y () * U->Get_y ();

      double z2 = U->Get_z () * U->Get_z ();


      double r2 = x2 + y2 + z2;

      *R = r2 * ipn * (U->*dxu) () / (U->*cd) () + U->Get_val ();



      R->Next ();


    }
  while (U->End ());



}









void
App_Boundary_Edge (double (ffunction::*dxu1) (void),
		   double (ffunction::*du1) (void),
		   double (ffunction::*cd1) (void),
		   double (ffunction::*dxu2) (void),
		   double (ffunction::*du2) (void),
		   double (ffunction::*cd2) (void),
		   ffunction * U, t_iterator face, size_t pn, ffunction * R)
{





  double ipn = 1.0 / double (pn);



  U->Set_Iterator (face);

  R->Set_Iterator (face);

  do
    {


      double x2 = U->Get_x () * U->Get_x ();

      double y2 = U->Get_y () * U->Get_y ();

      double z2 = U->Get_z () * U->Get_z ();


      double r2 = x2 + y2 + z2;


      *R = r2 * ipn * ((U->*dxu1) () / (U->*cd1) () +
		       (U->*dxu2) () / (U->*cd2) ()) + U->Get_val ();





      R->Next ();

    }
  while (U->End ());



}







void
App_Boundary_Corner (double (ffunction::*dxu1) (void),
		     double (ffunction::*du1) (void),
		     double (ffunction::*cd1) (void),
		     double (ffunction::*dxu2) (void),
		     double (ffunction::*du2) (void),
		     double (ffunction::*cd2) (void),
		     double (ffunction::*dxu3) (void),
		     double (ffunction::*du3) (void),
		     double (ffunction::*cd3) (void),
		     ffunction * U, t_iterator face, size_t pn, ffunction * R)
{





  double ipn = 1.0 / double (pn);


  U->Set_Iterator (face);

  R->Set_Iterator (face);


  do
    {


      double x2 = U->Get_x () * U->Get_x ();

      double y2 = U->Get_y () * U->Get_y ();

      double z2 = U->Get_z () * U->Get_z ();


      double r2 = x2 + y2 + z2;


      *R = r2 * ipn * ((U->*dxu1) () / (U->*cd1) () +
		       (U->*dxu2) () / (U->*cd2) () +
		       (U->*dxu3) () / (U->*cd3) ()) + U->Get_val ();





      R->Next ();
    }
  while (U->End ());






}
