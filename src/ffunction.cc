



//============================ Olliptic ffunction.cc ========================//
//  
//  ffunction : (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//






#include "ffunction.h"







void
ffunction::Make (domain * Dom)
{



  total_nodes = 1;

  my_rank = 0;

#ifdef OLLIN_MPI

  total_nodes = MPI::COMM_WORLD.Get_size ();

  my_rank = MPI::COMM_WORLD.Get_rank ();

#endif



  D = Dom;

  size_t ord_inter = D->Get_order ();

  order = ord_inter > 0 ? ord_inter : 2;


  switch (order)
    {

    case 2:

      lap = &ffunction::Lapc2;

      dulap = &ffunction::duLapc2;



      dd_xy = &ffunction::ddxyc2;

      dudd_xy = &ffunction::duddxyc2;


      dd_xz = &ffunction::ddxzc2;

      dudd_xz = &ffunction::duddxzc2;


      dd_yz = &ffunction::ddyzc2;

      dudd_yz = &ffunction::duddyzc2;



      d_x = &ffunction::dx2;

      dud_x = &ffunction::dudx2;


      d_y = &ffunction::dy2;

      dud_y = &ffunction::dudy2;


      d_z = &ffunction::dz2;

      dud_z = &ffunction::dudz2;



      dd_x = &ffunction::ddxc2;

      dudd_x = &ffunction::duddxc2;


      dd_y = &ffunction::ddyc2;

      dudd_y = &ffunction::duddyc2;


      dd_z = &ffunction::ddzc2;

      dudd_z = &ffunction::duddzc2;


      break;


    case 4:

      lap = &ffunction::Lapc4;

      dulap = &ffunction::duLapc4;



      dd_xy = &ffunction::ddxyc4;

      dudd_xy = &ffunction::duddxyc4;


      dd_xz = &ffunction::ddxzc4;

      dudd_xz = &ffunction::duddxzc4;


      dd_yz = &ffunction::ddyzc4;

      dudd_yz = &ffunction::duddyzc4;



      d_x = &ffunction::dx4;

      dud_x = &ffunction::dudx4;


      d_y = &ffunction::dy4;

      dud_y = &ffunction::dudy4;


      d_z = &ffunction::dz4;

      dud_z = &ffunction::dudz4;


      dd_x = &ffunction::ddxc4;

      dudd_x = &ffunction::duddxc4;


      dd_y = &ffunction::ddyc4;

      dudd_y = &ffunction::duddyc4;


      dd_z = &ffunction::ddzc4;

      dudd_z = &ffunction::duddzc4;

      break;


    case 6:

      lap = &ffunction::Lapc6;

      dulap = &ffunction::duLapc6;


/*            
            dd_xy = &ffunction::ddxyc6;

            dudd_xy = &ffunction::duddxyc6;


            dd_xz = &ffunction::ddxzc6;

            dudd_xz = &ffunction::duddxzc6;

            
            dd_yz = &ffunction::ddyzc6;

            dudd_yz = &ffunction::duddyzc6;
*/



      d_x = &ffunction::dx4;

      dud_x = &ffunction::dudx4;


      d_y = &ffunction::dy4;

      dud_y = &ffunction::dudy4;


      d_z = &ffunction::dz4;

      dud_z = &ffunction::dudz4;


      dd_x = &ffunction::ddxc4;

      dudd_x = &ffunction::duddxc4;


      dd_y = &ffunction::ddyc4;

      dudd_y = &ffunction::duddyc4;


      dd_z = &ffunction::ddzc4;

      dudd_z = &ffunction::duddzc4;


      break;

    case 8:

      lap = &ffunction::Lapc8;

      dulap = &ffunction::duLapc8;


/*
            dd_xy = &ffunction::ddxyc8;

            dudd_xy = &ffunction::duddxyc8;


            dd_xz = &ffunction::ddxzc8;

            dudd_xz = &ffunction::duddxzc8;

            
            dd_yz = &ffunction::ddyzc8;

            dudd_yz = &ffunction::duddyzc8;
*/

/*
      d_x = &ffunction::dx8;

      dud_x = &ffunction::dudx8;


      d_y = &ffunction::dy8;

      dud_y = &ffunction::dudy8;


      d_z = &ffunction::dz8;

      dud_z = &ffunction::dudz8;

*/

      d_x = &ffunction::dx4;

      dud_x = &ffunction::dudx4;


      d_y = &ffunction::dy4;

      dud_y = &ffunction::dudy4;


      d_z = &ffunction::dz4;

      dud_z = &ffunction::dudz4;


      
      dd_x = &ffunction::ddxc4;

      dudd_x = &ffunction::duddxc4;


      dd_y = &ffunction::ddyc4;

      dudd_y = &ffunction::duddyc4;


      dd_z = &ffunction::ddzc4;

      dudd_z = &ffunction::duddzc4;





      
      break;

    default:

      if (my_rank == 0)

	cout << "Seting Laplace operator error, wrong order " << endl;

      exit (1);

    }


  interpol = &ffunction::Linear_Interpol_In;





  it_2d = 0;





  n = D->Get_N ();


  f.resize (n);


  D->Fix_Iterators (init_corner, init_edge, init_face, init_inter);





  dx = D->Get_dx ();

  dy = D->Get_dy ();

  dz = D->Get_dz ();


  idx = 1.0 / dx;

  idy = 1.0 / dy;

  idz = 1.0 / dz;


  idx2 = idx * idx;

  idy2 = idy * idy;

  idz2 = idz * idz;



  i180 = 1.0 / 180;

  pi = 4 * atan (1);


}






ostream & operator<< (ostream & s, ffunction & F)
{

  s.setf (ios_base::scientific, ios_base::floatfield);

  s.precision (12);

  for (int i = 0; i < F.n; i++)

    s << F.Get_x (i) << "\t"
      << F.Get_y (i) << "\t" << F.Get_z (i) << "\t" << F (i) << endl;


  return (s);

}


void
ffunction::Print_Info (ostream & s)
{

  s << "Interpolation: " << order << endl
    << "  ------------------------------------------------------------------"
    << endl << *D << flush;



}







void
ffunction::Print_3D (string name,
		    double (*funct_to_comp) (double x, double y, double z))
{



  stringstream m_r;

  m_r << my_rank;



  double x,
    y,
    z;


  string outfile = name + m_r.str () + ".xyz";


  ofstream File (outfile.c_str (), fstream::app);

  File.setf (ios_base::scientific, ios_base::floatfield);

  File.precision (12);


  Set_Iterator (ALL);

  do
    {
      x = Get_x ();
      y = Get_y ();
      z = Get_z ();

      File << x << "\t" << y << "\t" << z << "\t" << Get_val () << "\t";

      if (funct_to_comp == NULL)

	File << endl;

      else

	File << Get_val () - funct_to_comp (x, y, z) << endl;


    }
  while (End ());


  File << flush;

  File.close ();


}








void
ffunction::Print_2D (string name,
		    const cut c,
		    const double x1,
		    bool interpolar,
		    double (*funct_to_comp) (double x, double y, double z),
		    m_interpol m_inter, size_t ord)
{


  Set_Interpol (m_inter, ord);

  stringstream m_r;

  m_r << my_rank;



  int nx = D->Get_Nx ();

  int ny = D->Get_Ny ();

  int nz = D->Get_Nz ();



  int i = (nx % 2 == 0) ? nx / 2 : (nx - 1) / 2;

  int j = (ny % 2 == 0) ? ny / 2 : (ny - 1) / 2;

  int k = (nz % 2 == 0) ? nz / 2 : (nz - 1) / 2;


  double x,
    y,
    z;

  switch (c)
    {




    case CUT_X:

      {
	string outfile = name + m_r.str () + ".yz";


	ofstream File (outfile.c_str (), fstream::app);

	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (12);

	if (D->x_is_in (x1) && interpolar)

	  x = x1;

	else

	  x = Get_x (i, j, k);


	if (interpolar)

	  for (j = 0; j < ny; j++)
	    {

	      for (k = 0; k < nz; k++)
		{


		  y = Get_y (i, j, k);
		  z = Get_z (i, j, k);

		  File << y << "\t"
		    << z << "\t" << Interpol_In (x, y, z) << "\t";

		  if (funct_to_comp == NULL)

		    File << endl;

		  else

		    File << Interpol_In (x, y, z) - funct_to_comp (x, y,
								   z) << endl;

		}

	      File << endl;
	    }


	else

	  for (j = 0; j < ny; j++)
	    {

	      for (k = 0; k < nz; k++)
		{


		  y = Get_y (i, j, k);
		  z = Get_z (i, j, k);


		  File << y << "\t" << z << "\t" << Get_val (i, j, k) << "\t";

		  if (funct_to_comp == NULL)

		    File << endl;

		  else

		    File << Get_val (i, j, k) - funct_to_comp (x, y,
							       z) << endl;


		}

	      File << endl;
	    }


	File << endl << flush;

	File.close ();


      }

      break;




    case CUT_Y:

      {
	string outfile = name + m_r.str () + ".xz";


	ofstream File (outfile.c_str (), fstream::app);

	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (12);

	if (D->y_is_in (x1) && interpolar)

	  y = x1;

	else

	  y = Get_y (i, j, k);



	if (interpolar)



	  for (i = 0; i < nx; i++)
	    {

	      for (k = 0; k < nz; k++)
		{

		  x = Get_x (i, j, k);
		  z = Get_z (i, j, k);

		  File << x << "\t"
		    << z << "\t" << Interpol_In (x, y, z) << "\t";

		  if (funct_to_comp == NULL)

		    File << endl;

		  else

		    File << Interpol_In (x, y, z) - funct_to_comp (x, y,
								   z) << endl;

		}

	      File << endl;
	    }



	else

	  for (i = 0; i < nx; i++)
	    {

	      for (k = 0; k < nz; k++)
		{

		  x = Get_x (i, j, k);
		  z = Get_z (i, j, k);

		  File << x << "\t" << z << "\t" << Get_val (i, j, k) << "\t";

		  if (funct_to_comp == NULL)

		    File << endl;

		  else

		    File << Get_val (i, j, k) - funct_to_comp (x, y,
							       z) << endl;

		}

	      File << endl;
	    }



	File << endl << flush;

	File.close ();


      }

      break;







    case CUT_Z:

      {
	string outfile = name + m_r.str () + ".xy";


	ofstream File (outfile.c_str (), fstream::app);

	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (12);

	if (D->z_is_in (x1) && interpolar)

	  z = x1;

	else

	  z = Get_z (i, j, k);




	if (interpolar)
	  for (i = 0; i < nx; i++)
	    {
	      for (j = 0; j < ny; j++)
		{


		  x = Get_x (i, j, k);
		  y = Get_y (i, j, k);

		  File << x << "\t"
		    << y << "\t" << Interpol_In (x, y, z) << "\t";
		  if (funct_to_comp == NULL)

		    File << endl;

		  else

		    File << Interpol_In (x, y, z) - funct_to_comp (x, y,
								   z) << endl;
		}

	      File << endl;
	    }

	else
	  for (i = 0; i < nx; i++)
	    {
	      for (j = 0; j < ny; j++)
		{


		  x = Get_x (i, j, k);
		  y = Get_y (i, j, k);

		  File << x << "\t" << y << "\t" << Get_val (i, j, k) << "\t";

		  if (funct_to_comp == NULL)

		    File << endl;

		  else

		    File << Get_val (i, j, k) - funct_to_comp (x, y,
							       z) << endl;
		}

	      File << endl;
	    }



	File << endl << flush;

	File.close ();


      }

      break;









    }





}








void
ffunction::Print_1D (string name,
		    const cut c,
		    const double x1,
		    const double x2,
		    bool interpolar,
		    double (*funct_to_comp) (double x, double y, double z),
		    m_interpol m_inter, size_t ord)
{


  Set_Interpol (m_inter, ord);

  stringstream m_r;

  m_r << my_rank;


  int nx = D->Get_Nx ();

  int ny = D->Get_Ny ();

  int nz = D->Get_Nz ();



  int i = (nx % 2 == 0) ? nx / 2 : (nx - 1) / 2;

  int j = (ny % 2 == 0) ? ny / 2 : (ny - 1) / 2;

  int k = (nz % 2 == 0) ? nz / 2 : (nz - 1) / 2;


  double x,
    y,
    z;

  switch (c)
    {




    case CUT_X:

      {
	string outfile = name + m_r.str () + ".x";


	ofstream File (outfile.c_str (), fstream::app);

	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (12);

	if (D->y_is_in (x1) && interpolar)

	  y = x1;

	else
	 if (D->y_is_in (x1))
	  {

	    j = Get_j (x1);

	    y = Get_y (i, j, k);


	  }
	else

	  y = Get_y (i, j, k);


	if (D->z_is_in (x2) && interpolar)

	  z = x2;

	else
	 if (D->z_is_in (x2))
	  {

	    k = Get_k (x2);

	    z = Get_z (i, j, k);


	  }
	else

	  z = Get_z (i, j, k);

	if (interpolar)

	  for (i = 0; i < nx; i++)
	    {


	      x = Get_x (i, j, k);

	      File << x << "\t" << Interpol_In (x, y, z) << "\t";

	      if (funct_to_comp == NULL)

		File << endl;

	      else

		File << Interpol_In (x, y, z) - funct_to_comp (x, y,
							       z) << endl;

	    }

	else

	  for (i = 0; i < nx; i++)
	    {


	      x = Get_x (i, j, k);

	      File << x << "\t" << Get_val (i, j, k) << "\t";

	      if (funct_to_comp == NULL)

		File << endl;

	      else

		File << Get_val (i, j, k) - funct_to_comp (x, y, z) << endl;


	    }




	File << endl << flush;


	File.close ();


      }

      break;








    case CUT_Y:

      {
	string outfile = name + m_r.str () + ".y";


	ofstream File (outfile.c_str (), fstream::app);

	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (12);

	if (D->x_is_in (x1) && interpolar)

	  x = x1;

	else
	 if (D->x_is_in (x1))
	  {

	    i = Get_i (x1);

	    x = Get_x (i, j, k);


	  }
	else

	  x = Get_x (i, j, k);


	if (D->z_is_in (x2) && interpolar)

	  z = x2;

	else
	 if (D->z_is_in (x2))
	  {

	    k = Get_k (x2);

	    z = Get_z (i, j, k);


	  }
	else

	  z = Get_z (i, j, k);



	if (interpolar)

	  for (j = 0; j < ny; j++)
	    {

	      y = Get_y (i, j, k);

	      File << y << "\t" << Interpol_In (x, y, z) << "\t";

	      if (funct_to_comp == NULL)

		File << endl;

	      else

		File << Interpol_In (x, y, z) - funct_to_comp (x, y,
							       z) << endl;

	    }

	else
	  for (j = 0; j < ny; j++)
	    {

	      y = Get_y (i, j, k);

	      File << y << "\t" << Get_val (i, j, k) << "\t";

	      if (funct_to_comp == NULL)

		File << endl;

	      else

		File << Get_val (i, j, k) - funct_to_comp (x, y, z) << endl;

	    }


	File << endl << flush;


	File.close ();


      }

      break;







    case CUT_Z:

      {
	string outfile = name + m_r.str () + ".z";


	ofstream File (outfile.c_str (), fstream::app);

	File.setf (ios_base::scientific, ios_base::floatfield);

	File.precision (12);

	if (D->x_is_in (x1) && interpolar)

	  x = x1;

	else
	 if (D->x_is_in (x1))
	  {

	    i = Get_i (x1);

	    x = Get_x (i, j, k);


	  }
	else

	  x = Get_x (i, j, k);


	if (D->y_is_in (x2) && interpolar)

	  y = x2;

	else
	 if (D->y_is_in (x2))
	  {

	    j = Get_j (x2);

	    y = Get_y (i, j, k);


	  }
	else

	  y = Get_y (i, j, k);


	if (interpolar)
	  for (k = 0; k < nz; k++)
	    {


	      z = Get_z (i, j, k);

	      File << z << "\t" << Interpol_In (x, y, z) << "\t";

	      if (funct_to_comp == NULL)

		File << endl;

	      else

		File << Interpol_In (x, y, z) - funct_to_comp (x, y,
							       z) << endl;
	    }

	else
	  for (k = 0; k < nz; k++)
	    {


	      z = Get_z (i, j, k);

	      File << z << "\t" << Get_val (i, j, k) << "\t";

	      if (funct_to_comp == NULL)

		File << endl;

	      else

		File << Get_val (i, j, k) - funct_to_comp (x, y, z) << endl;
	    }



	File << endl << flush;

	File.close ();


      }

      break;









    }





}






















void
ffunction::Print_t1D (string name, const double t)
{



  Set_Interpol (LINEAR);


  string outfile = name + ".xg";


  ofstream File (outfile.c_str (), fstream::app);

  File.setf (ios_base::scientific, ios_base::floatfield);

  File.precision (12);

  File << "\n\"Time = " << t << endl;


  int nx = D->Get_Nx ();

  int ny = D->Get_Ny ();

  int nz = D->Get_Nz ();


  int i = (nx % 2 == 0) ? nx / 2 : (nx - 1) / 2;

  int j = (ny % 2 == 0) ? ny / 2 : (ny - 1) / 2;

  int k = (nz % 2 == 0) ? nz / 2 : (nz - 1) / 2;

  for (i = 0; i < nx; i++)
    {

      double x = Get_x (i, j, k);

      File << x << "\t" << Interpol_In (x, 0, 0) << endl;
    }


  File << flush;

  File.close ();




}




void
ffunction::Print_t2D (string name, size_t cycle)
{



  string outfile = name + ".xyg";


  ofstream File (outfile.c_str (), fstream::app);


  File.setf (ios_base::scientific, ios_base::floatfield);

  File.precision (12);



  int nx = D->Get_Nx ();

  int ny = D->Get_Ny ();

  int nz = D->Get_Nz ();


  int i = (nx % 2 == 0) ? nx / 2 : (nx - 1) / 2;

  int j = (ny % 2 == 0) ? ny / 2 : (ny - 1) / 2;

  int k = (nz % 2 == 0) ? nz / 2 : (nz - 1) / 2;


  for (i = 0; i < nx; i++)
    {
      for (k = 0; k < nz; k++)
	{

	  double x = Get_x (i, j, k);

	  double z = Get_z (i, j, k);

	  File << x << "\t" << z << "\t" << Get_val (i, j, k) << endl;
	}

      File << endl;
    }
  File << endl << flush;

  File.close ();


  string scriptfile = name + ".gplot";


  ofstream FScript (scriptfile.c_str (), fstream::app);

  if (cycle == 0)

    FScript << "set term wxt" << endl
      << "set pm3d" << endl << "set hidden3d" << endl;

  FScript << "splot \"" << outfile.c_str ()
    << "\" index " << cycle << endl << "#pause -1" << endl;


  it_2d++;

}






bool
ffunction::Begin_End ()
{



  switch (iter)
    {


    case ALL:

      index = 0;

      index_end = n;

      break;

    case INSIDE:

      index = init_inter;

      index_end = n;



      break;



    case CORNER_XIYIZI:

      index = init_corner[XiYiZi];

      index_end = init_corner[XiYiZf];


      break;

    case CORNER_XIYIZF:

      index = init_corner[XiYiZf];

      index_end = init_corner[XfYiZi];


      break;


    case CORNER_XFYIZI:

      index = init_corner[XfYiZi];

      index_end = init_corner[XfYiZf];


      break;


    case CORNER_XFYIZF:


      index = init_corner[XfYiZf];

      index_end = init_corner[XiYfZi];

      break;



    case CORNER_XIYFZI:

      index = init_corner[XiYfZi];
      index_end = init_corner[XiYfZf];


      break;


    case CORNER_XIYFZF:

      index = init_corner[XiYfZf];
      index_end = init_corner[XfYfZi];


      break;


    case CORNER_XFYFZI:

      index = init_corner[XfYfZi];
      index_end = init_corner[XfYfZf];


      break;


    case CORNER_XFYFZF:

      index = init_corner[XfYfZf];
      index_end = init_edge[XiYi];


      break;



    case EDGE_XIYI:

      index = init_edge[XiYi];
      index_end = init_edge[XiYf];


      break;


    case EDGE_XIYF:

      index = init_edge[XiYf];
      index_end = init_edge[XfYi];


      break;


    case EDGE_XFYI:

      index = init_edge[XfYi];
      index_end = init_edge[XfYf];


      break;


    case EDGE_XFYF:

      index = init_edge[XfYf];
      index_end = init_edge[XiZi];


      break;



    case EDGE_XIZI:

      index = init_edge[XiZi];
      index_end = init_edge[XiZf];


      break;


    case EDGE_XIZF:

      index = init_edge[XiZf];
      index_end = init_edge[XfZi];


      break;


    case EDGE_XFZI:

      index = init_edge[XfZi];
      index_end = init_edge[XfZf];


      break;


    case EDGE_XFZF:

      index = init_edge[XfZf];
      index_end = init_edge[YiZi];


      break;



    case EDGE_YIZI:

      index = init_edge[YiZi];
      index_end = init_edge[YiZf];


      break;


    case EDGE_YIZF:

      index = init_edge[YiZf];
      index_end = init_edge[YfZi];


      break;


    case EDGE_YFZI:

      index = init_edge[YfZi];
      index_end = init_edge[YfZf];


      break;


    case EDGE_YFZF:

      index = init_edge[YfZf];
      index_end = init_face[Xi];


      break;





    case FACE_XI:

      index = init_face[Xi];
      index_end = init_face[Xf];


      break;


    case FACE_XF:

      index = init_face[Xf];
      index_end = init_face[Yi];


      break;


    case FACE_YI:

      index = init_face[Yi];
      index_end = init_face[Yf];


      break;


    case FACE_YF:

      index = init_face[Yf];
      index_end = init_face[Zi];


      break;



    case FACE_ZI:

      index = init_face[Zi];
      index_end = init_face[Zf];


      break;


    case FACE_ZF:

      index = init_face[Zf];
      index_end = init_inter;


      break;


    case BOUNDARY:

      index = 0;
      index_end = init_inter;


      break;



    default:


      index = 0;
      index_end = n;



    }


  return (true);

}






double
ffunction::NormL2 ()
{

  double l_norm = 0.0;

  double t_norm = 0.0;



  for (int i = 0; i < n; i++)

    l_norm += Get_val (i) * Get_val (i);



#ifdef OLLIN_MPI

  MPI_Allreduce (&l_norm, &t_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

  t_norm = l_norm;

#endif


  return (sqrt (t_norm * dx * dy * dz));

}




double
ffunction::NormL1 ()
{



  double l_norm = 0.0;

  double t_norm = 0.0;

  for (int i = 0; i < n; i++)

    l_norm += fabs (Get_val (i));




#ifdef OLLIN_MPI

  MPI_Allreduce (&l_norm, &t_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

  t_norm = l_norm;



#endif


  return (t_norm * dx * dy * dz);


}






double
ffunction::NormLInf ()
{



  double l_norm = 0.0;

  double t_norm = 0.0;



  for (int i = 0; i < n; i++)

    if (l_norm < fabs (Get_val (i)))

      l_norm = fabs (Get_val (i));





#ifdef OLLIN_MPI

  MPI_Allreduce (&l_norm, &t_norm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#else

  t_norm = l_norm;

#endif

  return (t_norm);



}






void
ffunction::Add_To_Norms ()
{


  double f = fabs (Get_val ());

  N_L1 += f;

  N_L2 += f * f;


  if (N_LInf < f)

    N_LInf = f;


}






void
ffunction::Sync_Norms ()
{




  double L1 = N_L1;

  double L2 = N_L2;

  double LInf = N_LInf;



  double global_L1;

  double global_L2;

  double global_LInf;



#ifdef OLLIN_MPI

  MPI_Allreduce (&L1, &global_L1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (&L2, &global_L2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (&LInf, &global_LInf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


#else

  global_L1 = L1;

  global_L2 = L2;

  global_LInf = LInf;


#endif



  N_L1 = global_L1;

  N_L2 = sqrt (global_L2);

  N_LInf = global_LInf;


}




void
ffunction::Sync ()
{



  if (D->Get_Bound_Zi () >= 0)

    Sync_Bound_Zi ();


  if (D->Get_Bound_Zf () >= 0)

    Sync_Bound_Zf ();




  if (D->Get_Bound_Yi () >= 0)

    Sync_Bound_Yi ();


  if (D->Get_Bound_Yf () >= 0)

    Sync_Bound_Yf ();




  if (D->Get_Bound_Xi () >= 0)

    Sync_Bound_Xi ();


  if (D->Get_Bound_Xf () >= 0)

    Sync_Bound_Xf ();




#ifdef OLLIN_MPI


  MPI::COMM_WORLD.Barrier ();


#endif


}




void
ffunction::Sync_Bound_Zi ()
{


#ifdef OLLIN_MPI



  const int dest = D->Get_Bound_Zi ();


  const int sendtag = Zi;

  const int recvtag = Zf;




  const int lnx = D->Get_Nx ();

  const int lny = D->Get_Ny ();

  const int bfz = D->Get_buffer_z ();

  int sr_size = lnx * lny * bfz;

  double data_buffer[sr_size];




  for (int k = bfz + 1; k < 2 * bfz + 1; k++)

    for (int j = 0; j < lny; j++)

      for (int i = 0; i < lnx; i++)

	data_buffer[i + (j + (k - (bfz + 1)) * lny) * lnx] =
	  Get_val (i, j, k);



  MPI_Status status;




  MPI_Sendrecv_replace (&data_buffer[0],
			sr_size, MPI_DOUBLE,
			dest, sendtag,
			dest, recvtag, MPI_COMM_WORLD, &status);







  for (int k = 0; k < bfz; k++)

    for (int j = 0; j < lny; j++)

      for (int i = 0; i < lnx; i++)

	Set (i, j, k, data_buffer[i + (j + k * lny) * lnx]);








#endif



}







void
ffunction::Sync_Bound_Zf ()
{




#ifdef OLLIN_MPI



  const int dest = D->Get_Bound_Zf ();



  const int sendtag = Zf;

  const int recvtag = Zi;




  const int lnx = D->Get_Nx ();

  const int lny = D->Get_Ny ();

  const int lnz = D->Get_Nz ();

  const int bfz = D->Get_buffer_z ();


  int sr_size = lnx * lny * bfz;

  double data_buffer[sr_size];




  for (int k = lnz - 2 * bfz - 1; k < lnz - bfz - 1; k++)

    for (int j = 0; j < lny; j++)

      for (int i = 0; i < lnx; i++)

	data_buffer[i + (j + (k - (lnz - 2 * bfz - 1)) * lny) * lnx] =
	  Get_val (i, j, k);



  MPI_Status status;





  MPI_Sendrecv_replace (&data_buffer[0],
			sr_size, MPI_DOUBLE,
			dest, sendtag,
			dest, recvtag, MPI_COMM_WORLD, &status);





  for (int k = lnz - bfz; k < lnz; k++)

    for (int j = 0; j < lny; j++)

      for (int i = 0; i < lnx; i++)

	Set (i, j, k, data_buffer[i + (j + (k - (lnz - bfz)) * lny) * lnx]);








#endif

}






void
ffunction::Sync_Bound_Yi ()
{



#ifdef OLLIN_MPI



  const int dest = D->Get_Bound_Yi ();


  const int sendtag = Yi;

  const int recvtag = Yf;



  const int lnx = D->Get_Nx ();

  const int bfy = D->Get_buffer_y ();

  const int lnz = D->Get_Nz ();



  int sr_size = lnx * bfy * lnz;

  double data_buffer[sr_size];



  for (int j = bfy + 1; j < 2 * bfy + 1; j++)

    for (int k = 0; k < lnz; k++)

      for (int i = 0; i < lnx; i++)

	data_buffer[i + (j - (bfy + 1) + k * bfy) * lnx] = Get_val (i, j, k);


  MPI_Status status;



  MPI_Sendrecv_replace (&data_buffer[0],
			sr_size, MPI_DOUBLE,
			dest, sendtag,
			dest, recvtag, MPI_COMM_WORLD, &status);






  for (int j = 0; j < bfy; j++)

    for (int k = 0; k < lnz; k++)

      for (int i = 0; i < lnx; i++)

	Set (i, j, k, data_buffer[i + (j + k * bfy) * lnx]);





#endif



}











void
ffunction::Sync_Bound_Yf ()
{




#ifdef OLLIN_MPI


  const int dest = D->Get_Bound_Yf ();



  const int sendtag = Yf;

  const int recvtag = Yi;



  const int lnx = D->Get_Nx ();

  const int lny = D->Get_Ny ();

  const int lnz = D->Get_Nz ();

  const int bfy = D->Get_buffer_y ();

  int sr_size = lnx * bfy * lnz;

  double data_buffer[sr_size];



  for (int j = lny - 2 * bfy - 1; j < lny - bfy - 1; j++)

    for (int k = 0; k < lnz; k++)

      for (int i = 0; i < lnx; i++)

	data_buffer[i + (j - (lny - 2 * bfy - 1) + k * bfy) * lnx] =
	  Get_val (i, j, k);



  MPI_Status status;





  MPI_Sendrecv_replace (&data_buffer[0],
			sr_size, MPI_DOUBLE,
			dest, sendtag,
			dest, recvtag, MPI_COMM_WORLD, &status);






  for (int j = lny - bfy; j < lny; j++)

    for (int k = 0; k < lnz; k++)

      for (int i = 0; i < lnx; i++)

	Set (i, j, k, data_buffer[i + (j - (lny - bfy) + k * bfy) * lnx]);









#endif



}








void
ffunction::Sync_Bound_Xi ()
{



#ifdef OLLIN_MPI


  const int dest = D->Get_Bound_Xi ();


  const int sendtag = Xi;

  const int recvtag = Xf;


  const int bfx = D->Get_buffer_x ();

  const int lny = D->Get_Ny ();

  const int lnz = D->Get_Nz ();



  int sr_size = bfx * lny * lnz;

  double data_buffer[sr_size];



  for (int i = bfx + 1; i < 2 * bfx + 1; i++)

    for (int k = 0; k < lnz; k++)

      for (int j = 0; j < lny; j++)

	data_buffer[i - (bfx + 1) + (j + k * lny) * bfx] = Get_val (i, j, k);


  MPI_Status status;



  MPI_Sendrecv_replace (&data_buffer[0],
			sr_size, MPI_DOUBLE,
			dest, sendtag,
			dest, recvtag, MPI_COMM_WORLD, &status);






  for (int i = 0; i < bfx; i++)

    for (int k = 0; k < lnz; k++)

      for (int j = 0; j < lny; j++)

	Set (i, j, k, data_buffer[i + (j + k * lny) * bfx]);








#endif



}



void
ffunction::Sync_Bound_Xf ()
{




#ifdef OLLIN_MPI


  const int dest = D->Get_Bound_Xf ();



  const int sendtag = Xf;

  const int recvtag = Xi;




  const int lnx = D->Get_Nx ();

  const int lny = D->Get_Ny ();

  const int lnz = D->Get_Nz ();

  const int bfx = D->Get_buffer_x ();


  int sr_size = bfx * lny * lnz;

  double data_buffer[sr_size];


  for (int i = lnx - 2 * bfx - 1; i < lnx - bfx - 1; i++)

    for (int k = 0; k < lnz; k++)

      for (int j = 0; j < lny; j++)

	data_buffer[i - (lnx - 2 * bfx - 1) + (j + k * lny) * bfx] =
	  Get_val (i, j, k);


  MPI_Status status;





  MPI_Sendrecv_replace (&data_buffer[0],
			sr_size, MPI_DOUBLE,
			dest, sendtag,
			dest, recvtag, MPI_COMM_WORLD, &status);







  for (int i = lnx - bfx; i < lnx; i++)

    for (int k = 0; k < lnz; k++)

      for (int j = 0; j < lny; j++)

	Set (i, j, k, data_buffer[i - (lnx - bfx) + (j + k * lny) * bfx]);



#endif



}







void
ffunction::Add (ffunction & u)
{

  Set_Iterator (ALL);

  do
    {

      double fxyz = Get_val () + u.Get_val ();

      Set (fxyz);


    }
  while (End ());


}






void
ffunction::Transferir (ffunction & u, t_iterator inter, m_interpol method,
		      size_t ord, bool all, bool to_zero)
{





  u.Set_Interpol (method, ord);


  Set_Iterator (inter);

  do
    {

      double x = Get_x ();

      double y = Get_y ();

      double z = Get_z ();


      if (u.Is_In_Nb (x, y, z))

	Set (u.Interpol_In (x, y, z));

      else
       if (all && u.Is_In (x, y, z))

	Set (u.Interpol_In (x, y, z));

      else
       if (to_zero)

	Set (0.0);


    }
  while (End ());




#ifdef OLLIN_MPI



  vector < double >send_buffer;

  vector < int >send_buffer_index_I;

  vector < int >send_buffer_index_J;

  vector < int >send_buffer_index_K;


  vector < int >send_counts (total_nodes, 0);

  int size_prev = 0;

  for (int r = 0; r < my_rank; r++)
    {
      size_t nx = Get_domain ()->Get_Nx (r);

      size_t ny = Get_domain ()->Get_Ny (r);

      size_t nz = Get_domain ()->Get_Nz (r);


      size_t bfx = Get_domain ()->Get_buffer_x ();

      size_t bfy = Get_domain ()->Get_buffer_y ();

      size_t bfz = Get_domain ()->Get_buffer_z ();


      size_prev = send_buffer.size ();

      for (size_t i = bfx; i < nx - bfx; i++)

	for (size_t j = bfy; j < ny - bfy; j++)

	  for (size_t k = bfz; k < nz - bfz; k++)
	    {

	      double x = Get_domain ()->Get_x (i, j, k, r);

	      double y = Get_domain ()->Get_y (i, j, k, r);

	      double z = Get_domain ()->Get_z (i, j, k, r);

	      if (u.Is_In_Nb (x, y, z) && !u.Is_In_Nb (x, y, z, r))
		{

		  send_buffer.push_back (u.Interpol_In (x, y, z));

		  send_buffer_index_I.push_back (i);

		  send_buffer_index_J.push_back (j);

		  send_buffer_index_K.push_back (k);

		}

	      else
	       if (all && u.Is_In (x, y, z) && !u.Is_In_Nb (x, y, z, r))
		{

		  send_buffer.push_back (u.Interpol_In (x, y, z));

		  send_buffer_index_I.push_back (i);

		  send_buffer_index_J.push_back (j);

		  send_buffer_index_K.push_back (k);

		}



	    }




      send_counts[r] = send_buffer.size () - size_prev;


    }






  for (int r = my_rank + 1; r < total_nodes; r++)
    {
      size_t nx = Get_domain ()->Get_Nx (r);

      size_t ny = Get_domain ()->Get_Ny (r);

      size_t nz = Get_domain ()->Get_Nz (r);


      size_t bfx = Get_domain ()->Get_buffer_x ();

      size_t bfy = Get_domain ()->Get_buffer_y ();

      size_t bfz = Get_domain ()->Get_buffer_z ();


      size_prev = send_buffer.size ();

      for (size_t i = bfx; i < nx - bfx; i++)

	for (size_t j = bfy; j < ny - bfy; j++)

	  for (size_t k = bfz; k < nz - bfz; k++)
	    {

	      double x = Get_domain ()->Get_x (i, j, k, r);

	      double y = Get_domain ()->Get_y (i, j, k, r);

	      double z = Get_domain ()->Get_z (i, j, k, r);

	      if (u.Is_In_Nb (x, y, z) && !u.Is_In_Nb (x, y, z, r))
		{


		  send_buffer.push_back (u.Interpol_In (x, y, z));

		  send_buffer_index_I.push_back (i);

		  send_buffer_index_J.push_back (j);

		  send_buffer_index_K.push_back (k);


		}

	      else if (all && u.Is_In (x, y, z) && !u.Is_In (x, y, z, r))
		{

		  send_buffer.push_back (u.Interpol_In (x, y, z));

		  send_buffer_index_I.push_back (i);

		  send_buffer_index_J.push_back (j);

		  send_buffer_index_K.push_back (k);

		}




	    }



      send_counts[r] = send_buffer.size () - size_prev;



    }





  size_t size_send_buffer = 0;

  vector < int >send_displs (total_nodes, 0);

  for (int r = 0; r < total_nodes; r++)
    {

      send_displs[r] = size_send_buffer;

      size_send_buffer += send_counts[r];

    }


/*

    for(int r = 0; r < total_nodes; r++)
        
            cout << "my_rank: " << my_rank << "\t"
                 << "to rank: " << r << "\t"
                 << "send counts: " << send_counts[r] << "\t"
                 << "send displs: " << send_displs[r] << endl;
    
    cout << endl << flush;
    
*/



  vector < int >total_counts (total_nodes * total_nodes);


  MPI_Allgather (&send_counts[0],
		 total_nodes,
		 MPI_INT,
		 &total_counts[0], total_nodes, MPI_INT, MPI_COMM_WORLD);




  size_t size_recv_buffer = 0;

  vector < int >recv_counts (total_nodes, 0);

  vector < int >recv_displs (total_nodes, 0);

  for (int r = 0; r < total_nodes; r++)
    {

      recv_displs[r] = size_recv_buffer;

      size_recv_buffer += total_counts[r * total_nodes + my_rank];

      recv_counts[r] = total_counts[r * total_nodes + my_rank];


    }




  /*  
     for(int r = 0; r < total_nodes; r++)

     cout << "from rank: " << r << "\t"
     << "my_rank: " << my_rank << "\t"
     << "recv counts: " << recv_counts[r] << "\t\t"
     << "recv displs: " << recv_displs[r] << endl
     << "size_recv_buffer: " << size_recv_buffer << endl;

     cout << endl << flush;

   */




  vector < double >recv_buffer (size_recv_buffer);



  MPI_Alltoallv (&send_buffer[0], &send_counts[0], &send_displs[0],
		 MPI_DOUBLE, &recv_buffer[0], &recv_counts[0],
		 &recv_displs[0], MPI_DOUBLE, MPI_COMM_WORLD);






  vector < int >recv_buffer_index_I (size_recv_buffer);



  MPI_Alltoallv (&send_buffer_index_I[0], &send_counts[0], &send_displs[0],
		 MPI_INT, &recv_buffer_index_I[0], &recv_counts[0],
		 &recv_displs[0], MPI_INT, MPI_COMM_WORLD);






  vector < int >recv_buffer_index_J (size_recv_buffer);



  MPI_Alltoallv (&send_buffer_index_J[0], &send_counts[0], &send_displs[0],
		 MPI_INT, &recv_buffer_index_J[0], &recv_counts[0],
		 &recv_displs[0], MPI_INT, MPI_COMM_WORLD);







  vector < int >recv_buffer_index_K (size_recv_buffer);



  MPI_Alltoallv (&send_buffer_index_K[0], &send_counts[0], &send_displs[0],
		 MPI_INT, &recv_buffer_index_K[0], &recv_counts[0],
		 &recv_displs[0], MPI_INT, MPI_COMM_WORLD);





  for (size_t indx = 0; indx < size_recv_buffer; indx++)

    Set (recv_buffer_index_I[indx],
	 recv_buffer_index_J[indx],
	 recv_buffer_index_K[indx], recv_buffer[indx]);



  MPI::COMM_WORLD.Barrier ();




#endif





  Sync ();



}






void
ffunction::Set_In (ffunction & u, domain * Dom, bool all)
{





  u.Set_Interpol (LINEAR, 1);


  Set_Iterator (ALL);

  do
    {

      double x = Get_x ();

      double y = Get_y ();

      double z = Get_z ();


      if (u.Is_In_Nb (x, y, z) && Dom->is_in_nb (x, y, z))

	Set (u.Interpol_In (x, y, z));

      else

       if (all && u.Is_In (x, y, z) && Dom->is_in (x, y, z))

	Set (u.Interpol_In (x, y, z));



    }
  while (End ());





}
