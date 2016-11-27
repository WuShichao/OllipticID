



//============================ Olliptic domain.cc ========================//
//  
//  domain : (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//




#include "domain.h"







void
domain::Make (size_t Nx, size_t Ny, size_t Nz,
	      double dx, double dy, double dz,
	      double cx, double cy, double cz,
	      size_t ord, t_symmetry s_x, t_symmetry s_y, t_symmetry s_z)
{


  symm_x = s_x;

  symm_y = s_y;

  symm_z = s_z;


  level = 0;

  order = ord;

  size = 1;

  my_rank = 0;


#ifdef OLLIN_MPI


  size = MPI::COMM_WORLD.Get_size ();

  my_rank = MPI::COMM_WORLD.Get_rank ();

#endif



  xi.resize (size);

  yi.resize (size);

  zi.resize (size);


  xf.resize (size);

  yf.resize (size);

  zf.resize (size);



  nx.resize (size);

  ny.resize (size);

  nz.resize (size);




  Set (dx, dy, dz, cx, cy, cz);


  Set_Buffer (order / 2, order / 2, order / 2);


  Set_N (Nx, Ny, Nz);



  Fix_Points ();



}










void
domain::Set (double Dx, double Dy, double Dz, double Cx, double Cy, double Cz)
{


  dx = Dx;

  dy = Dy;

  dz = Dz;


  cx = Cx;

  cy = Cy;

  cz = Cz;


  size_t a,
    b,
    c;

  size_t i = 0,
    j = 0,
    k = 0;

  Get_Partition (a, b, c);

  Get_Index_Topology (a, b, c, i, j, k);




  if (i < a - 1)

    bound_xf = Get_rank (a, b, c, i + 1, j, k);

  else
    {

      if (symm_x == ROTANT)
	{

	  bound_xf = -1;

	  bound_yf = -3;


	}
      else
       if (symm_x == REFLECT)

	bound_xf = -2;


      else
       if (symm_z != ROTANT)

	bound_xf = -1;





    }

  // para la cara "derecha"

  if (i > 0)

    bound_xi = Get_rank (a, b, c, i - 1, j, k);

  else

    bound_xi = -1;






  //Para la cara "izquierda"

  if (j < b - 1)

    bound_yf = Get_rank (a, b, c, i, j + 1, k);

  else
    {

      if (symm_y == ROTANT)
	{

	  bound_yf = -1;

	  bound_zf = -3;

	}
      else
       if (symm_y == REFLECT)

	bound_yf = -2;

      else if (symm_x != ROTANT)

	bound_yf = -1;



    }


  // para la cara "derecha"

  if (j > 0)

    bound_yi = Get_rank (a, b, c, i, j - 1, k);

  else

    bound_yi = -1;







  //Para la cara "izquierda"

  if (k < c - 1)

    bound_zf = Get_rank (a, b, c, i, j, k + 1);

  else
    {

      if (symm_z == ROTANT)
	{

	  bound_zf = -1;

	  bound_xf = -3;

	}
      else
       if (symm_z == REFLECT)

	bound_zf = -2;


      else
       if (symm_y != ROTANT)

	bound_zf = -1;





    }

  // para la cara "derecha"

  if (k > 0)

    bound_zi = Get_rank (a, b, c, i, j, k - 1);

  else

    bound_zi = -1;




}





void
domain::Set_N (size_t Nx, size_t Ny, size_t Nz)
{





  size_t a,
    b,
    c;

  size_t i = 0,
    j = 0,
    k = 0;

  Get_Partition (a, b, c);

  Get_Index_Topology (a, b, c, i, j, k);




  // Fija Frontera X global


  Valid_N (Nx);

  Valid_N (Ny);

  Valid_N (Nz);



  gxi = cx - 0.5 * double (Nx - 1) * dx - double (buff_nx - 1) * dx;

  gyi = cy - 0.5 * double (Ny - 1) * dy - double (buff_ny - 1) * dy;

  gzi = cz - 0.5 * double (Nz - 1) * dz - double (buff_nz - 1) * dz;

  Nx += buff_nx - 1;

  Ny += buff_ny - 1;

  Nz += buff_nz - 1;




  gxf = gxi + double (Nx - 1) * dx + double (buff_nx - 1) * dx;

  gyf = gyi + double (Ny - 1) * dy + double (buff_ny - 1) * dy;

  gzf = gzi + double (Nz - 1) * dz + double (buff_nz - 1) * dz;

  Nx += buff_nx - 1;

  Ny += buff_ny - 1;

  Nz += buff_nz - 1;


  double gLx = gxf - gxi;

  double gLy = gyf - gyi;

  double gLz = gzf - gzi;



  for (size_t ii = 0; ii < size; ii++)

    {


      nx[ii] = Nx;

      ny[ii] = Ny;

      nz[ii] = Nz;



      xi[ii] = gxi;

      xf[ii] = gxf;


      yi[ii] = gyi;

      yf[ii] = gyf;


      zi[ii] = gzi;

      zf[ii] = gzf;


      size_t I = 0,
	J = 0,
	K = 0;

      Get_Index_Topology (a, b, c, I, J, K, ii);

      if (symm_x == REFLECT || symm_z == ROTANT)

	Set_Boundary (xi[ii], xf[ii], nx[ii], gxi, gxf, gLx, dx, I, a,
		      buff_nx, true);

      else

	Set_Boundary (xi[ii], xf[ii], nx[ii], gxi, gxf, gLx, dx, I, a,
		      buff_nx, false);


      if (symm_y == REFLECT || symm_x == ROTANT)

	Set_Boundary (yi[ii], yf[ii], ny[ii], gyi, gyf, gLy, dy, J, b,
		      buff_ny, true);

      else

	Set_Boundary (yi[ii], yf[ii], ny[ii], gyi, gyf, gLy, dy, J, b,
		      buff_ny, false);



      if (symm_z == REFLECT || symm_y == ROTANT)

	Set_Boundary (zi[ii], zf[ii], nz[ii], gzi, gzf, gLz, dz, K, c,
		      buff_nz, true);

      else

	Set_Boundary (zi[ii], zf[ii], nz[ii], gzi, gzf, gLz, dz, K, c,
		      buff_nz, false);




    }



  Lx = xf[my_rank] - xi[my_rank];

  Ly = yf[my_rank] - yi[my_rank];

  Lz = zf[my_rank] - zi[my_rank];




}




void
Set_Boundary (double &ti,
	      double &tf,
	      size_t & Nt,
	      double gti,
	      double gtf,
	      double Lt,
	      double dt, size_t k, size_t c, size_t b_nt, bool symm)
{


  if (symm)
    {

      gtf = 0;

      Lt = 0.5 * Lt;
    }

  double lti = gti + k * Lt / double (c);

  double ltf = gtf - (c - 1 - k) * Lt / double (c);






  while (ti + double (Nt - 1) * dt > ltf)

    Nt--;

  if (ti + double (Nt - 1) * dt < ltf - 0.5 * dt)

    Nt++;



  tf = ti + double (Nt - 1) * dt;



  while (tf - double (Nt - 1) * dt < lti)

    Nt--;

  if (tf - double (Nt - 1) * dt > lti + 0.5 * dt)

    Nt++;

  ti = tf - double (Nt - 1) * dt;



  if (k > 0)
    {

      ti -= double (b_nt) * dt;


      Nt += b_nt;
    }

  if (k < c - 1 || symm)
    {

      tf += double (b_nt) * dt;


      Nt += b_nt;

    }



}





void
domain::Set_Buffer (size_t bf_nx, size_t bf_ny, size_t bf_nz)
{
  buff_nx = bf_nx >= 1 ? bf_nx : 1;

  buff_ny = bf_ny >= 1 ? bf_ny : 1;

  buff_nz = bf_nz >= 1 ? bf_nz : 1;
}








void
domain::Fix_Iterators (int Corner[8], int Edge[12], int Face[6], int &Inter)
{




  for (size_t i = 0; i < 8; i++)

    Corner[i] = corner[i];

  for (size_t i = 0; i < 12; i++)

    Edge[i] = edge[i];

  for (size_t i = 0; i < 6; i++)

    Face[i] = face[i];

  Inter = inter;


}



void
domain::Fix_Points ()
{


  I.clear ();

  J.clear ();

  K.clear ();

  ijk.clear ();



  I.reserve (nx[my_rank] * ny[my_rank] * nz[my_rank]);

  J.reserve (nx[my_rank] * ny[my_rank] * nz[my_rank]);

  K.reserve (nx[my_rank] * ny[my_rank] * nz[my_rank]);

  ijk.assign (nx[my_rank] * ny[my_rank] * nz[my_rank], 0);



  size_t index = 0;


  const size_t bfx = Get_buffer_x ();

  const size_t bfy = Get_buffer_y ();

  const size_t bfz = Get_buffer_z ();




  //Esquina (xi,yi,zi)


  corner[XiYiZi] = index;




  for (int k = bfz - 1; k >= 0; k--)

    for (int j = bfy - 1; j >= 0; j--)

      for (int i = bfx - 1; i >= 0; i--)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}




  //Esquina (xi,yi,zf)


  corner[XiYiZf] = index;



  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (int j = bfy - 1; j >= 0; j--)

      for (int i = bfx - 1; i >= 0; i--)
	{



	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}


  //Esquina (xf,yi,zi)


  corner[XfYiZi] = index;


  for (int k = bfz - 1; k >= 0; k--)

    for (int j = bfy - 1; j >= 0; j--)

      for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}





  //Esquina (xf,yi,zf)


  corner[XfYiZf] = index;

  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (int j = bfy - 1; j >= 0; j--)

      for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}





  //Esquina (xi,yf,zi)


  corner[XiYfZi] = index;

  for (int k = bfz - 1; k >= 0; k--)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (int i = bfx - 1; i >= 0; i--)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}




  //Esquina (xi,yf,zf)


  corner[XiYfZf] = index;

  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (int i = bfx - 1; i >= 0; i--)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}




  //Esquina (xf,yf,zi)

  corner[XfYfZi] = index;

  for (int k = bfz - 1; k >= 0; k--)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}




  //Esquina (xf,yf,zf)

  corner[XfYfZf] = index;

  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}



  // Esquina (xi,yi,z)

  edge[XiYi] = index;

  for (int i = bfx - 1; i >= 0; i--)

    for (int j = bfy - 1; j >= 0; j--)

      for (size_t k = bfz; k < nz[my_rank] - bfz; k++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Esquina (xi,yf,z)

  edge[XiYf] = index;

  for (int i = bfx - 1; i >= 0; i--)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (size_t k = bfz; k < nz[my_rank] - bfz; k++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Esquina (xf,yi,z)

  edge[XfYi] = index;

  for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)

    for (int j = bfy - 1; j >= 0; j--)

      for (size_t k = bfz; k < nz[my_rank] - bfz; k++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}




  // Esquina (xf,yf,z)

  edge[XfYf] = index;

  for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (size_t k = bfz; k < nz[my_rank] - bfz; k++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Esquina (xi,y,zi)

  edge[XiZi] = index;

  for (int i = bfx - 1; i >= 0; i--)

    for (int k = bfz - 1; k >= 0; k--)

      for (size_t j = bfy; j < ny[my_rank] - bfy; j++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}




  // Esquina (xi,y,zf)

  edge[XiZf] = index;

  for (int i = bfx - 1; i >= 0; i--)

    for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

      for (size_t j = bfy; j < ny[my_rank] - bfy; j++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Esquina (xf,y,zi)

  edge[XfZi] = index;

  for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)

    for (int k = bfz - 1; k >= 0; k--)

      for (size_t j = bfy; j < ny[my_rank] - bfy; j++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}




  // Esquina (xf,y,zf)

  edge[XfZf] = index;

  for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)

    for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

      for (size_t j = bfy; j < ny[my_rank] - bfy; j++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}




  // Esquina (x,yi,zi)

  edge[YiZi] = index;


  for (int k = bfz - 1; k >= 0; k--)

    for (int j = bfy - 1; j >= 0; j--)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}






  // Esquina (x,yi,zf)

  edge[YiZf] = index;


  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (int j = bfy - 1; j >= 0; j--)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{

	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Esquina (x,yf,zi)

  edge[YfZi] = index;


  for (int k = bfz - 1; k >= 0; k--)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}




  // Esquina (x,yf,zf)

  edge[YfZf] = index;


  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}






  // Cara (xi,y,z)

  face[Xi] = index;

  for (int i = bfx - 1; i >= 0; i--)

    for (size_t k = bfz; k < nz[my_rank] - bfz; k++)

      for (size_t j = bfy; j < ny[my_rank] - bfy; j++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Cara (xf,y,z)


  face[Xf] = index;

  for (size_t i = nx[my_rank] - bfx; i < nx[my_rank]; i++)

    for (size_t k = bfz; k < nz[my_rank] - bfz; k++)

      for (size_t j = bfy; j < ny[my_rank] - bfy; j++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}







  // Cara (x,yi,z)


  face[Yi] = index;

  for (int j = bfy - 1; j >= 0; j--)

    for (size_t k = bfz; k < nz[my_rank] - bfz; k++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}







  // Cara (x,yf,z)

  face[Yf] = index;

  for (size_t j = ny[my_rank] - bfy; j < ny[my_rank]; j++)

    for (size_t k = bfz; k < nz[my_rank] - bfz; k++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}




  // Cara (x,y,zi)


  face[Zi] = index;

  for (int k = bfz - 1; k >= 0; k--)

    for (size_t j = bfy; j < ny[my_rank] - bfy; j++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;


	}






  // Cara (x,y,zf) incluye fronteras


  face[Zf] = index;

  for (size_t k = nz[my_rank] - bfz; k < nz[my_rank]; k++)

    for (size_t j = bfy; j < ny[my_rank] - bfy; j++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





  // Puntos interiores

  inter = index;

  for (size_t k = bfz; k < nz[my_rank] - bfz; k++)

    for (size_t j = bfy; j < ny[my_rank] - bfy; j++)

      for (size_t i = bfx; i < nx[my_rank] - bfx; i++)
	{


	  I.push_back (i);

	  J.push_back (j);

	  K.push_back (k);

	  ijk[INDEXORDER_IJK] = index;

	  index++;

	}





}











void
domain::Get_Partition (size_t & a, size_t & b, size_t & c)
{





  bool cont = true;


  //prime numbers :s

  double P[20] = { 1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67
  };


  for (double l = 0; l < 4 && cont; l++)

    for (double m = 0; m < 4 && cont; m++)

      for (double n = 1; n < 4 && cont; n++)

	for (int i = 0; i < 20 && cont; i++)

	  for (int j = 0; j < 20 && cont; j++)

	    for (int k = 0; k < 20 && cont; k++)

	      if (size ==
		  size_t (pow (P[i], n) * pow (P[j], n + m) *
			  pow (P[k], n + m + l)))
		{



		  a = size_t (pow (P[i], n));

		  b = size_t (pow (P[j], n + m));

		  c = size_t (pow (P[k], n + m + l));


		  ordena (a, b, c);


		  cont = false;

		}



}




void
domain::Get_Index_Topology (const size_t l, const size_t m, const size_t n,
			    size_t & i, size_t & j, size_t & k, int rk)
{



  bool cont = true;

  size_t rank = rk < 0 ? my_rank : rk;

  for (i = 0; i < l && cont; i++)

    for (j = 0; j < m && cont; j++)

      for (k = 0; k < n && cont; k++)

	if (rank == Get_rank (l, m, n, i, j, k))

	  cont = false;


  i--;
  j--;
  k--;


}






const size_t
domain::Get_i (const double x) const
{


  //Binary search algorithm for  the index:
  //      o     o   x o     o     o
  //      i-2   i-1   i     i+1   i+2

  size_t i = 0;


  int sup = nx[my_rank];

  int inf = 0;

  int mid = nx[my_rank] % 2 == 0 ? nx[my_rank] / 2 : (nx[my_rank] - 1) / 2;

  if (xi[my_rank] - x < 0.0 && x - xf[my_rank] < 0.0)
    {

      i = mid;

      do
	{


	  if (xi[my_rank] + i * dx > x || dcomp (xi[my_rank] + i * dx, x))
	    {

	      sup = i;

	      i =
		(sup + inf) % 2 == 0 ? (sup + inf) / 2 : (sup + inf - 1) / 2;

	    }

	  else

	    {

	      inf = i;

	      i =
		(sup + inf) % 2 == 0 ? (sup + inf) / 2 : (sup + inf - 1) / 2;

	    }

	}
      while (inf + 1 < sup);


    }

  else

    i = (xi[my_rank] >= x) ? 0 : nx[my_rank] - 2;

  i++;

  return (i);

}



const size_t
domain::Get_j (const double y) const
{


  size_t i = 0;


  int sup = ny[my_rank];

  int inf = 0;

  int mid = ny[my_rank] % 2 == 0 ? ny[my_rank] / 2 : (ny[my_rank] - 1) / 2;

  if (yi[my_rank] < y && y < yf[my_rank])
    {

      i = mid;

      do
	{


	  if (yi[my_rank] + i * dy > y || dcomp (yi[my_rank] + i * dy, y))
	    {

	      sup = i;

	      i =
		(sup + inf) % 2 == 0 ? (sup + inf) / 2 : (sup + inf - 1) / 2;

	    }

	  else

	    {

	      inf = i;

	      i =
		(sup + inf) % 2 == 0 ? (sup + inf) / 2 : (sup + inf - 1) / 2;

	    }

	}
      while (inf + 1 < sup);


    }
  else

    i = (yi[my_rank] >= y) ? 0 : ny[my_rank] - 2;

  i++;

  return (i);

}


const size_t
domain::Get_k (const double z) const
{


  size_t i = 0;


  int sup = nz[my_rank];

  int inf = 0;

  int mid = nz[my_rank] % 2 == 0 ? nz[my_rank] / 2 : (nz[my_rank] - 1) / 2;

  if (zi[my_rank] < z && z < zf[my_rank])
    {

      i = mid;

      do
	{


	  if (zi[my_rank] + i * dz > z || dcomp (zi[my_rank] + i * dz, z))
	    {

	      sup = i;

	      i =
		(sup + inf) % 2 == 0 ? (sup + inf) / 2 : (sup + inf - 1) / 2;

	    }

	  else

	    {

	      inf = i;

	      i =
		(sup + inf) % 2 == 0 ? (sup + inf) / 2 : (sup + inf - 1) / 2;

	    }

	}
      while (inf + 1 < sup);


    }
  else

    i = (zi[my_rank] >= z) ? 0 : nz[my_rank] - 2;

  i++;

  return (i);

}





void
domain::Find_cell_xyz (const double x, const double y, const double z,
		       double &xi, double &yi, double &zi,
		       double &xf, double &yf, double &zf)
{



  size_t i = Get_i (x);

  size_t j = Get_j (y);

  size_t k = Get_k (z);


  xi = Get_x (i - 1);
  xf = Get_x (i);

  yi = Get_y (j - 1);
  yf = Get_y (j);

  zi = Get_z (k - 1);
  zf = Get_z (k);


}






const bool
domain::Is_1_to_2 (const domain & D) const
{
  return (dcomp (dx, 2 * D.dx) &&
	  dcomp (dy, 2 * D.dy) && dcomp (dz, 2 * D.dz));

}





void
domain::Print_to_File (string name)
{


  stringstream m_r;

  m_r << my_rank;

  string outfile = name + m_r.str () + ".xg";


  ofstream File (outfile.c_str (), fstream::out);


  for (double x = xi[my_rank]; x <= xf[my_rank] + 0.1 * dx; x += dx)

    for (double y = yi[my_rank]; y <= yf[my_rank] + 0.1 * dy; y += dy)

      for (double z = zi[my_rank]; z <= zf[my_rank] + 0.1 * dz; z += dz)

	File << x << "\t" << y << "\t" << z << endl;



  File.close ();

}






ostream & operator<< (ostream & s, const domain & D)
{


  s.setf (ios_base::scientific, ios_base::floatfield);

  s.precision (3);


  s << "  nx = " << D.Get_Nx () << ", "
    << "\t\tny = " << D.Get_Ny () << ", "
    << "\t\tnz = " << D.Get_Nz () << endl
    << "  dx = " << D.Get_dx () << ", "
    << "\tdy = " << D.Get_dy () << ", " << "\tdz = " << D.Get_dz () << endl;

  s.setf (ios_base::fixed, ios_base::floatfield);

  s.precision (2);


  s << "  lx = " << D.Get_Length_x () << ", "
    << "\t\tly = " << D.Get_Length_y () << ", "
    << "\t\tlz = " << D.Get_Length_z () << endl;



  s << "  x-->[" << D.Get_Xi () << "," << D.Get_Xf () << "]"
    << "\ty-->[" << D.Get_Yi () << "," << D.Get_Yf () << "]"
    << "\tz-->[" << D.Get_Zi () << "," << D.Get_Zf () << "]"

/*
      << ":: ["        
      << D.Get_xi() << "," << D.Get_xf() << "]X["
      << D.Get_yi() << "," << D.Get_yf() << "]X["
      << D.Get_zi() << "," << D.Get_zf() << "]"
*/
    << endl;


  return (s);

}









domain
operator+ (const domain & A, const domain & B)
{


  double dx = A.dx;

  double dy = A.dy;

  double dz = A.dz;



  size_t order = A.order;



  t_symmetry s_x = A.symm_x;

  t_symmetry s_y = A.symm_y;

  t_symmetry s_z = A.symm_z;



  double xi = min (A.gxi, B.gxi);

  double xf = max (A.gxf, B.gxf);


  double yi = min (A.gyi, B.gyi);

  double yf = max (A.gyf, B.gyf);


  double zi = min (A.gzi, B.gzi);

  double zf = max (A.gzf, B.gzf);


  double cx = 0.5 * (xi + xf);

  double cy = 0.5 * (yi + yf);

  double cz = 0.5 * (zi + zf);


  size_t nx = size_t ((xf - xi) / dx + 1);

  size_t ny = size_t ((yf - yi) / dy + 1);

  size_t nz = size_t ((zf - zi) / dz + 1);




  domain C (nx, ny, nz, dx, dy, dz, cx, cy, cz, order, s_x, s_y, s_z);


  return (C);

}





domain
operator- (const domain & A, const domain & B)
{



  domain C,
    D;


  C.xi = max (A.xi, B.xi);

  C.xf = min (A.xf, B.xf);


  D.xi = min (C.xi, C.xf);

  D.xf = max (C.xi, C.xf);



  C.yi = max (A.yi, B.yi);

  C.yf = min (A.yf, B.yf);


  D.yi = min (C.yi, C.yf);

  D.yf = max (C.yi, C.yf);



  C.zi = max (A.zi, B.zi);

  C.zf = min (A.zf, B.zf);


  D.zi = min (C.zi, C.zf);

  D.zf = max (C.zi, C.zf);



  return (D);

}






bool
operator == (const domain & D1, const domain & D2)
{

  return (D1.xi[D1.my_rank] == D2.xi[D2.my_rank]
	  && D1.xf[D1.my_rank] == D2.xf[D2.my_rank]
	  && D1.yi[D1.my_rank] == D2.yi[D2.my_rank]
	  && D1.yf[D1.my_rank] == D2.yf[D2.my_rank]
	  && D1.zi[D1.my_rank] == D2.zi[D2.my_rank]
	  && D1.zf[D1.my_rank] == D2.zf[D2.my_rank]);

}





const bool
domain::operator< (const domain & D) const
{

  return (D.xi[my_rank] < xi[my_rank] && xf[my_rank] < D.xf[my_rank] &&
	  D.yi[my_rank] < yi[my_rank] && yf[my_rank] < D.yf[my_rank] &&
	  D.zi[my_rank] < zi[my_rank] && zf[my_rank] < D.zf[my_rank]);

}


const bool
domain::operator<= (const domain & D) const
{

  return (D.xi[my_rank] <= xi[my_rank] && xf[my_rank] <= D.xf[my_rank] &&
	  D.yi[my_rank] <= yi[my_rank] && yf[my_rank] <= D.yf[my_rank] &&
	  D.zi[my_rank] <= zi[my_rank] && zf[my_rank] <= D.zf[my_rank]);

}







const bool
domain::operator> (const domain & D) const
{

  return (xi[my_rank] < D.xi[my_rank] && D.xf[my_rank] < xf[my_rank] &&
	  yi[my_rank] < D.yi[my_rank] && D.yf[my_rank] < yf[my_rank] &&
	  zi[my_rank] < D.zi[my_rank] && D.zf[my_rank] < zf[my_rank]);

}


const bool
domain::operator>= (const domain & D) const
{

  return (xi[my_rank] <= D.xi[my_rank] && D.xf[my_rank] <= xf[my_rank] &&
	  yi[my_rank] <= D.yi[my_rank] && D.yf[my_rank] <= yf[my_rank] &&
	  zi[my_rank] <= D.zi[my_rank] && D.zf[my_rank] <= zf[my_rank]);

}




const bool
domain::operator&& (const domain & D) const
{


  return (D >= D - *this);


}





// Auxiliar function out of the class


void
ordena (size_t & a, size_t & b, size_t & c)
{

  if (a > b)

    swap (a, b);

  if (a > c)

    swap (a, c);

  if (b > c)

    swap (b, c);

}
