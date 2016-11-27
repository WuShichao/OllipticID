



//============================ Olliptic multigrid.cc ========================//
//  
//     
//
//
//
//
//
//============================== Pablo Galaviz 2009 ========================//





#include "multigrid.h"



void
multigrid::Make (vector < domain > Dom, t_domain d, vector < string > name)
{


  D = Dom;

  dom = d;

  total_grids = D.size ();



  for (size_t d = 0; d < total_grids; d++)

    u.push_back (scalarfields (&D[d], name));

  l0 = 0;

  for (size_t d = 1; d < total_grids; d++)
    {

      if (D[l0] < D[d])

	l0 = d;

    }


  size_t l = 0;

  size_t ibox = 0;

  num_dom_in_level.clear ();

  for (size_t d = 0; d < total_grids; d++)

    if (l < D[d].Get_level ())

      l = D[d].Get_level ();

  levels = l + 1;




  for (l = 0; l < levels; l++)
    {

      ibox = 0;

      for (size_t d = 0; d < total_grids; d++)

	if (l == D[d].Get_level ())

	  ibox++;

      num_dom_in_level.push_back (ibox);


//    cout << "level:" << l << "\t"
//         << "boxes:" << ibox << endl;

    }



}





size_t
multigrid::Get_dom_index (const size_t lev, const size_t boxi)
{

  size_t indx = 0;

  if (lev < levels)
    {
      for (size_t l = 0; l < lev; l++)

	indx += num_dom_in_level[l];

      if (boxi < num_dom_in_level[lev])

	indx += boxi;

    }

  return (indx);

}




void
multigrid::Transferir (size_t from, string from_var, size_t to, string to_var,
		       t_iterator iter, m_interpol method, size_t ord,
		       bool all, bool to_zero)
{

  from = from < levels ? from : levels - 1;

  to = to < levels ? to : levels - 1;

  for (size_t ibox = 0; ibox < num_dom_in_level[from]; ibox++)

    for (size_t jbox = 0; jbox < num_dom_in_level[to]; jbox++)

      u[Get_dom_index (to, jbox)].Get_ffunction (to_var)->
	Transferir (*u[Get_dom_index (from, ibox)].Get_ffunction (from_var),
		    iter, method, ord, all, to_zero);



}








void
multigrid::Print_Interpol_1D (string dir, string name,
			      double LD,
			      double _dx,
			      const cut c,
			      const double x1,
			      const double x2,
			      double (*funct_to_comp) (double x, double y,
						       double z),
			      m_interpol minter, size_t order)
{


  for (size_t ibox = 0; ibox < total_grids; ibox++)

    u[ibox].Get_ffunction (name)->Set_Interpol (minter, order);


  double x = 0;

  double y = 0;

  double z = 0;

  double xi = -0.5 * LD;

  double xf = 0.5 * LD;

  double Dx = _dx > 0.0 ? _dx : LD / 1000;

  size_t in_level;


  switch (c)
    {




    case CUT_X:
      {

	string outfile = dir + name + ".X";

	ofstream File (outfile.c_str (), fstream::app);


	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (18);




	if (D[l0].y_is_in (x1))

	  y = x1;

	if (D[l0].z_is_in (x2))

	  z = x2;


	if (dom == OCTANT ||
	    dom == BITANT_X ||
	    dom == QUADRANT_XY || dom == QUADRANT_XZ || dom == ROTATE_Z)

	  xf = 0;





	for (x = xi; x <= xf; x += Dx)

	  if (D[l0].is_in (x, y, z))
	    {

	      in_level = Get_dom_index (l0, 0);

	      for (size_t l = 0; l < levels; l++)

		for (size_t i = 0; i < num_dom_in_level[l]; i++)
		  {

		    size_t ibox = Get_dom_index (l, i);

		    if (D[ibox].is_in_nb (x, y, z))

		      in_level = ibox;


		  }

/*
                    cout << x << "\t"
                         << y << "\t"
                         << z << endl;
*/
	      double fxyz =
		u[in_level].Get_ffunction (name)->Interpol_In (x, y, z);



	      File << x << "\t" << fxyz;


	      if (funct_to_comp != NULL)

		File << "\t" << fxyz - funct_to_comp (x, y, z);

	      File << endl << flush;





	    }








	File << endl << flush;


#ifdef MPI_OLLIN

	//== Finish MPI ==//    

	MPI::COMM_WORLD.Barrier ();


#endif

	if (File.is_open ())

	  File.close ();



      }

      break;





    case CUT_Y:
      {

	string outfile = dir + name + ".Y";

	ofstream File (outfile.c_str (), fstream::app);


	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (18);




	if (D[l0].x_is_in (x1))

	  x = x1;

	if (D[l0].z_is_in (x2))

	  z = x2;



	if (dom == OCTANT ||
	    dom == BITANT_Y ||
	    dom == QUADRANT_XY || dom == QUADRANT_YZ || dom == ROTATE_X)

	  xf = 0;




	if (funct_to_comp == NULL)

	  for (y = xi; y <= xf; y += Dx)
	    {


	      in_level = Get_dom_index (l0, 0);

	      for (size_t l = 0; l < levels; l++)

		for (size_t i = 0; i < num_dom_in_level[l]; i++)
		  {

		    size_t ibox = Get_dom_index (l, i);

		    if (D[ibox].is_in_nb (x, y, z))

		      in_level = ibox;


		  }


	      File << y << "\t" << u[in_level].Get_ffunction (name)->
		Interpol_In (x, y, z) << endl;


	    }

	else

	  for (y = xi; y <= xf; y += Dx)
	    {


	      in_level = Get_dom_index (l0, 0);

	      for (size_t l = 0; l < levels; l++)



		for (size_t i = 0; i < num_dom_in_level[l]; i++)
		  {

		    size_t ibox = Get_dom_index (l, i);

		    if (D[ibox].is_in_nb (x, y, z))

		      in_level = ibox;
		  }

	      double fxyz =
		u[in_level].Get_ffunction (name)->Interpol_In (x, y, z);

	      File << y << "\t" << fxyz << "\t"
		<< fxyz - funct_to_comp (x, y, z) << endl;


	    }



	File << endl << flush;


#ifdef MPI_OLLIN

	//== Finish MPI ==//    

	MPI::COMM_WORLD.Barrier ();


#endif

	if (File.is_open ())

	  File.close ();



      }

      break;







    case CUT_Z:
      {

	string outfile = dir + name + ".Z";

	ofstream File (outfile.c_str (), fstream::app);


	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (18);




	if (D[l0].x_is_in (x1))

	  x = x1;

	if (D[l0].y_is_in (x2))

	  y = x2;



	if (dom == OCTANT ||
	    dom == BITANT_Z ||
	    dom == QUADRANT_XZ || dom == QUADRANT_YZ || dom == ROTATE_Y)

	  xf = 0;




	if (funct_to_comp == NULL)

	  for (z = xi; z <= xf; z += Dx)
	    {


	      in_level = Get_dom_index (l0, 0);

	      for (size_t l = 0; l < levels; l++)

		for (size_t i = 0; i < num_dom_in_level[l]; i++)
		  {

		    size_t ibox = Get_dom_index (l, i);

		    if (D[ibox].is_in_nb (x, y, z))

		      in_level = ibox;


		  }


	      File << z << "\t" << u[in_level].Get_ffunction (name)->
		Interpol_In (x, y, z) << endl;


	    }

	else

	  for (z = xi; z <= xf; z += Dx)
	    {


	      in_level = Get_dom_index (l0, 0);

	      for (size_t l = 0; l < levels; l++)



		for (size_t i = 0; i < num_dom_in_level[l]; i++)
		  {

		    size_t ibox = Get_dom_index (l, i);

		    if (D[ibox].is_in_nb (x, y, z))

		      in_level = ibox;
		  }

	      double fxyz =
		u[in_level].Get_ffunction (name)->Interpol_In (x, y, z);

	      File << z << "\t" << fxyz << "\t"
		<< fxyz - funct_to_comp (x, y, z) << endl;


	    }



	File << endl << flush;


#ifdef MPI_OLLIN

	//== Finish MPI ==//    

	MPI::COMM_WORLD.Barrier ();


#endif

	if (File.is_open ())

	  File.close ();



      }

    }

}








void
multigrid::Print_Interpol_2D (string dir,
			      string name,
			      double LD1,
			      double LD2,
			      double _dx1,
			      double _dx2,
			      const cut c,
			      const double x1,
			      double (*funct_to_comp) (double x, double y,
						       double z),
			      m_interpol minter, size_t order)
{


  cout << "Print 2D interpol dir: " << dir << " name: " << name << " dir+name: "  << dir+name+".XY"<< endl;


  for (size_t ibox = 0; ibox < total_grids; ibox++)

    u[ibox].Get_ffunction (name)->Set_Interpol (minter, order);


  double x = 0;

  double y = 0;

  double z = 0;

  double x1i = -0.25 * LD1;

  double x2i = -0.25 * LD2;

  double x1f = 0.25 * LD1;

  double x2f = 0.25 * LD2;

  double Dx1 = _dx1 > 0.0 ? _dx1 : LD1 / 200;

  double Dx2 = _dx2 > 0.0 ? _dx2 : LD2 / 200;

  size_t in_level;
  string outfile;

  switch (c)
    {




    case CUT_X:
      {

	outfile = dir + name + ".YZ";

	ofstream File (outfile.c_str(), fstream::app);


	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (18);




	if (D[l0].x_is_in (x1))

	  x = x1;


	if (dom == OCTANT ||
	    dom == BITANT_Y ||
	    dom == QUADRANT_XY || dom == QUADRANT_YZ || dom == ROTATE_X)

	  x1f = 0;



	if (dom == OCTANT ||
	    dom == BITANT_Z ||
	    dom == QUADRANT_XZ || dom == QUADRANT_YZ || dom == ROTATE_Y)

	  x2f = 0;



	if (funct_to_comp == NULL)

	  for (y = x1i; y <= x1f; y += Dx1)
	    for (z = x2i; z <= x2f; z += Dx2)
	      {


		in_level = Get_dom_index (l0, 0);

		for (size_t l = 0; l < levels; l++)

		  for (size_t i = 0; i < num_dom_in_level[l]; i++)
		    {

		      size_t ibox = Get_dom_index (l, i);

		      if (D[ibox].is_in_nb (x, y, z))

			in_level = ibox;


		    }


		File << y << "\t" << z << "\t" << u[in_level].
		  Get_ffunction (name)->Interpol_In (x, y, z) << endl;


	      }

	else

	  for (y = x1i; y <= x1f; y += Dx1)
	    for (z = x2i; z <= x2f; z += Dx2)
	      {


		in_level = Get_dom_index (l0, 0);

		for (size_t l = 0; l < levels; l++)



		  for (size_t i = 0; i < num_dom_in_level[l]; i++)
		    {

		      size_t ibox = Get_dom_index (l, i);

		      if (D[ibox].is_in_nb (x, y, z))

			in_level = ibox;
		    }

		double fxyz =
		  u[in_level].Get_ffunction (name)->Interpol_In (x, y, z);

		File << y << "\t" << z << "\t" << fxyz << "\t"
		  << fxyz - funct_to_comp (x, y, z) << endl;


	      }



	File << endl;

	File.close ();



      }

      break;


    case CUT_Y:
      {

	outfile = dir + name + ".XZ";

	ofstream File (outfile.c_str (), fstream::app);


	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (18);




	if (D[l0].y_is_in (x1))

	  y = x1;


	if (dom == OCTANT ||
	    dom == BITANT_X ||
	    dom == QUADRANT_XY || dom == QUADRANT_XZ || dom == ROTATE_Z)

	  x1f = 0;



	if (dom == OCTANT ||
	    dom == BITANT_Z ||
	    dom == QUADRANT_XZ || dom == QUADRANT_YZ || dom == ROTATE_Y)

	  x2f = 0;



	if (funct_to_comp == NULL)

	  for (x = x1i; x <= x1f; x += Dx1)
	    for (z = x2i; z <= x2f; z += Dx2)
	      {


		in_level = Get_dom_index (l0, 0);

		for (size_t l = 0; l < levels; l++)

		  for (size_t i = 0; i < num_dom_in_level[l]; i++)
		    {

		      size_t ibox = Get_dom_index (l, i);

		      if (D[ibox].is_in_nb (x, y, z))

			in_level = ibox;


		    }


		File << x << "\t" << z << "\t" << u[in_level].
		  Get_ffunction (name)->Interpol_In (x, y, z) << endl;


	      }

	else

	  for (x = x1i; x <= x1f; x += Dx1)
	    for (z = x2i; z <= x2f; z += Dx2)
	      {


		in_level = Get_dom_index (l0, 0);

		for (size_t l = 0; l < levels; l++)



		  for (size_t i = 0; i < num_dom_in_level[l]; i++)
		    {

		      size_t ibox = Get_dom_index (l, i);

		      if (D[ibox].is_in_nb (x, y, z))

			in_level = ibox;
		    }

		double fxyz =
		  u[in_level].Get_ffunction (name)->Interpol_In (x, y, z);

		File << x << "\t" << z << "\t" << fxyz << "\t"
		  << fxyz - funct_to_comp (x, y, z) << endl;


	      }



	File << endl;

	File.close ();



      }

      break;





    case CUT_Z:
      {

	outfile = dir + name + ".XY";

	ofstream File (outfile.c_str (), fstream::app);


	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (18);




	if (D[l0].z_is_in (x1))

	  z = x1;


	if (dom == OCTANT ||
	    dom == BITANT_X ||
	    dom == QUADRANT_XY || dom == QUADRANT_XZ || dom == ROTATE_Z)

	  x1f = 0;



	if (dom == OCTANT ||
	    dom == BITANT_Y ||
	    dom == QUADRANT_XY || dom == QUADRANT_YZ || dom == ROTATE_X)

	  x2f = 0;



	if (funct_to_comp == NULL)

	  for (x = x1i; x <= x1f; x += Dx1){

	    for (y = x2i; y <= x2f; y += Dx2)
	      {


		in_level = Get_dom_index (l0, 0);

		for (size_t l = 0; l < levels; l++)

		  for (size_t i = 0; i < num_dom_in_level[l]; i++)
		    {

		      size_t ibox = Get_dom_index (l, i);

		      if (D[ibox].is_in_nb (x, y, z))

			in_level = ibox;


		    }


		File << x << "\t" << y << "\t" << u[in_level].
		  Get_ffunction (name)->Interpol_In (x, y, z) << endl;


	      }
	    File << endl;

	  }

	else

	  for (x = x1i; x <= x1f; x += Dx1)
	    for (y = x2i; y <= x2f; y += Dx2)
	      {


		in_level = Get_dom_index (l0, 0);

		for (size_t l = 0; l < levels; l++)



		  for (size_t i = 0; i < num_dom_in_level[l]; i++)
		    {

		      size_t ibox = Get_dom_index (l, i);

		      if (D[ibox].is_in_nb (x, y, z))

			in_level = ibox;
		    }

		double fxyz =
		  u[in_level].Get_ffunction (name)->Interpol_In (x, y, z);

		File << x << "\t" << y << "\t" << fxyz << "\t"
		  << fxyz - funct_to_comp (x, y, z) << endl;


	      }



	File << endl;

	File.close ();



      }

      break;





    }


      cout << "Out file: " << outfile << endl;


}






double
multigrid::Get_Norm_LInf (size_t k, string var_name)
{


  double g_norm = 0;

  if (k < levels)

    for (size_t i = 0; i < num_dom_in_level[k]; i++)
      {

	double norm = Get_ffunction (k, i, var_name)->Get_Norm_LInf ();

	g_norm = norm > g_norm ? norm : g_norm;


      }


  return (g_norm);

}








void
multigrid::Print_IO_BAM (string name, string var_name, m_interpol minter,
			 size_t order)
{




  fstream File_In (name.c_str (), fstream::in | fstream::out);


  double x,
    y,
    z,
    val;


  File_In.setf (ios_base::scientific, ios_base::floatfield);


  cout.setf (ios_base::scientific, ios_base::floatfield);

  cout.precision (16);


  cout << "reading from: " << name << endl;


  for (size_t ibox = 0; ibox < total_grids; ibox++)

    u[ibox].Get_ffunction (var_name)->Set_Interpol (minter, order);


  do
    {


      int pos_r1 = File_In.tellg ();


      File_In >> x >> y >> z;


      File_In >> val;

      int pos_r2 = File_In.tellg ();



      size_t in_level = Get_dom_index (l0, 0);

      if (!D[in_level].is_in (x, y, z))
	{

	  cout << "Error: bam grid bigger than Olliptic grid" << endl;

	  exit (1);

	}

      for (size_t l = 0; l < levels; l++)

	for (size_t i = 0; i < num_dom_in_level[l]; i++)
	  {

	    size_t ibox = Get_dom_index (l, i);

	    if (D[ibox].is_in_nb (x, y, z))

	      in_level = ibox;


	  }




      File_In.seekp (pos_r1);


      File_In.precision (9);

      File_In.width (16);
      File_In << x;

      File_In.width (17);
      File_In << y;

      File_In.width (17);
      File_In << z;


      File_In.width (24);
      File_In.precision (15);
      File_In << u[in_level].Get_ffunction (var_name)->Interpol_In (x, y,
								   z) <<
	flush;



      File_In.seekg (pos_r2 + 1);




    }
  while (!File_In.eof ());



  File_In.close ();



}

















void
multigrid::Print_IO_Zcode (string name, string var_name, m_interpol minter,
			   size_t order)
{




  fstream File_In (name.c_str (), fstream::in);

  string name_out = name+"_m";

  fstream File_Out (name_out.c_str (), fstream::out);

  double x,y,z,val;


  File_In.setf (ios_base::scientific, ios_base::floatfield);
  File_Out.setf (ios_base::scientific, ios_base::floatfield);



  cout << "reading from: " << name << endl;


  for (size_t ibox = 0; ibox < total_grids; ibox++)

    u[ibox].Get_ffunction (var_name)->Set_Interpol (minter, order);


  do
    {


      File_In >> x >> y >> z;


      File_In >> val;

      
      size_t in_level = Get_dom_index (l0, 0);

      if (!D[in_level].is_in (x, y, z))
	{

	  cout << "X: " << x << "\t" 
	       << "Y: " << y << "\t" 
	       << "Z: " << z << endl;

	  cout << "Error: bam grid bigger than Olliptic grid" << endl;

	  exit (1);

	}

      for (size_t l = 0; l < levels; l++)

	for (size_t i = 0; i < num_dom_in_level[l]; i++)
	  {

	    size_t ibox = Get_dom_index (l, i);

	    if (D[ibox].is_in_nb (x, y, z))

	      in_level = ibox;


	  }


           

      File_Out.precision (16);

      //File_In.width(16);
      File_Out << x << " ";

      //File_In.width(17);
      File_Out << y << " ";

      //File_In.width(17);
      File_Out << z << " ";


      File_Out << u[in_level].Get_ffunction (var_name)->Interpol_In (x, y,z) << endl;
      
 



    }
  while (!File_In.eof ());



  File_In.close ();
  File_Out.close ();



}
