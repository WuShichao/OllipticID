//
// interpolation.cc
//  
// Made by Pablo Galaviz
// Login   <pablo@NerV>
// 
// Started on  Wed May 27 12:21:47 2009 Pablo Galaviz
// Started on  Wed May 27 12:21:47 2009 Pablo Galaviz
//



//======================== Olliptic interpolation.cc =====================//
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
ffunction::Set_Interpol (m_interpol minter, size_t ord)
{




  switch (minter)
    {


    case LINEAR:

      interpol = &ffunction::Linear_Interpol_In;

      int_ord = 1;

//            cout << "Linear interpolation selected :p" << endl;

      break;



    case LAGRANGE:


      interpol = &ffunction::Lagrange_Interpol_In;


      int_ord = ord;


      idxO = pow (idx, int (int_ord));

      idyO = pow (idy, int (int_ord));

      idzO = pow (idz, int (int_ord));

      //          cout << "Lagrange interpolation selected :p" << endl;


      break;


    case NEWTON:


      interpol = &ffunction::Newton_Interpol_In;


      int_ord = 2;

      //cout << "Newton interpolation selected :p" << endl;



      break;


    case HERMITE:


      interpol = &ffunction::Hermite_Interpol_In;


      int_ord = 2;


      //cout << "Hermite interpolation selected :p" << endl;


      break;


    case AVERAGE:


      interpol = &ffunction::Average_Interpol_In;



      break;



    }


}




double
ffunction::Linear_Interpol_In (double x, double y, double z)
{




  size_t i = D->Get_i (x);

  size_t j = D->Get_j (y);

  size_t k = D->Get_k (z);


  if (dcomp (x, Get_x (i, j, k)) && dcomp (y, Get_y (i, j, k))
      && dcomp (z, Get_z (i, j, k)))

    {

      return (Get_val (i, j, k));
    }

  else
    {

      i = i == 0 ? 1 : i;

      j = j == 0 ? 1 : j;

      k = k == 0 ? 1 : k;


      double dxf = (x - Get_x (i - 1, j - 1, k - 1)) * idx;

      double dyf = (y - Get_y (i - 1, j - 1, k - 1)) * idy;

      double dzf = (z - Get_z (i - 1, j - 1, k - 1)) * idz;

      double dxi = (1 - dxf);

      double dyi = (1 - dyf);

      double dzi = (1 - dzf);


      return (((Get_val (i - 1, j - 1, k - 1) * dxi +
		Get_val (i, j - 1, k - 1) * dxf) * dyi +
	       (Get_val (i - 1, j, k - 1) * dxi +
		Get_val (i, j, k - 1) * dxf) * dyf) * dzi +
	      ((Get_val (i - 1, j - 1, k) * dxi +
		Get_val (i, j - 1, k) * dxf) * dyi +
	       (Get_val (i - 1, j, k) * dxi +
		Get_val (i, j, k) * dxf) * dyf) * dzf);


    }


}





double
ffunction::Lagrange_Interpol_In (double x, double y, double z)
{




  size_t i = D->Get_i (x);

  size_t j = D->Get_j (y);

  size_t k = D->Get_k (z);


  if (dcomp (x, Get_x (i, j, k)) && dcomp (y, Get_y (i, j, k))
      && dcomp (z, Get_z (i, j, k)))

    {

      return (Get_val (i, j, k));
    }

  else
    {


      int ii = int_ord % 2 == 0 ? i - int_ord / 2 : i - (int_ord + 1) / 2;

      int ji = int_ord % 2 == 0 ? j - int_ord / 2 : j - (int_ord + 1) / 2;

      int ki = int_ord % 2 == 0 ? k - int_ord / 2 : k - (int_ord + 1) / 2;



      if (ii < 0)

	ii = 0;

      else
       if (ii + int_ord >= D->Get_Nx ())

	ii = D->Get_Nx () - int_ord - 1;


      if (ji < 0)

	ji = 0;

      else
       if (ji + int_ord >= D->Get_Ny ())

	ji = D->Get_Ny () - int_ord - 1;



      if (ki < 0)

	ki = 0;

      else
       if (ki + int_ord >= D->Get_Nz ())

	ki = D->Get_Nz () - int_ord - 1;





      double uk[int_ord + 1][int_ord + 1];

      for (size_t I = 0; I <= int_ord; I++)

	for (size_t J = 0; J <= int_ord; J++)
	  {
	    uk[I][J] = 0;

	    for (size_t K = 0; K <= int_ord; K++)

	      uk[I][J] +=
		Get_val (ii + I, ji + J, ki + K) * LagrangePoly (CoordZ, z,
								 ki + K, ki,
								 ki +
								 int_ord);

	  }


      double ujk[int_ord + 1];

      for (size_t I = 0; I <= int_ord; I++)
	{

	  ujk[I] = 0;

	  for (size_t J = 0; J <= int_ord; J++)

	    ujk[I] +=
	      uk[I][J] * LagrangePoly (CoordY, y, ji + J, ji, ji + int_ord);

	}


      double fxyz = 0;

      for (size_t I = 0; I <= int_ord; I++)

	fxyz += ujk[I] * LagrangePoly (CoordX, x, ii + I, ii, ii + int_ord);


      return (fxyz);

    }



}




double
ffunction::LagrangePoly (coord co, const double t, int i, int ni, int nf)
{


  double lp = 1.0;

  switch (co)
    {

    case CoordX:

      {

	const double ti = D->Get_x (i);

	for (int j = ni; j < i; j++)

	  lp *= (t - D->Get_x (j)) / (ti - D->Get_x (j));

	for (int j = i + 1; j <= nf; j++)

	  lp *= (t - D->Get_x (j)) / (ti - D->Get_x (j));

      }

      break;

    case CoordY:

      {

	const double ti = D->Get_y (i);

	for (int j = ni; j < i; j++)

	  lp *= (t - D->Get_y (j)) / (ti - D->Get_y (j));

	for (int j = i + 1; j <= nf; j++)

	  lp *= (t - D->Get_y (j)) / (ti - D->Get_y (j));

      }

      break;

    case CoordZ:

      {

	const double ti = D->Get_z (i);

	for (int j = ni; j < i; j++)

	  lp *= (t - D->Get_z (j)) / (ti - D->Get_z (j));

	for (int j = i + 1; j <= nf; j++)

	  lp *= (t - D->Get_z (j)) / (ti - D->Get_z (j));

      }

    }

  return (lp);


}







double
ffunction::Newton_Interpol_In (double x, double y, double z)
{




  size_t i = D->Get_i (x);

  size_t j = D->Get_j (y);

  size_t k = D->Get_k (z);




  int ii = int_ord % 2 == 0 ? i - int_ord / 2 : i - (int_ord + 1) / 2;

  int ji = int_ord % 2 == 0 ? j - int_ord / 2 : j - (int_ord + 1) / 2;

  int ki = int_ord % 2 == 0 ? k - int_ord / 2 : k - (int_ord + 1) / 2;



  if (ii < 0)

    ii = 0;

  else
   if (ii + int_ord >= D->Get_Nx ())

    ii = D->Get_Nx () - int_ord - 1;


  if (ji < 0)

    ji = 0;

  else
   if (ji + int_ord >= D->Get_Ny ())

    ji = D->Get_Ny () - int_ord - 1;



  if (ki < 0)

    ki = 0;

  else
   if (ki + int_ord >= D->Get_Nz ())

    ki = D->Get_Nz () - int_ord - 1;






  valarray < double >f (int_ord + 1);

  valarray < double >t (int_ord + 1);


  double uk[int_ord + 1][int_ord + 1];

  for (size_t I = 0; I <= int_ord; I++)

    for (size_t J = 0; J <= int_ord; J++)
      {


	for (size_t K = 0; K <= int_ord; K++)
	  {

	    f[K] = Get_val (ii + I, ji + J, ki + K);

	    t[K] = Get_z (ii + I, ji + J, ki + K);

	  }

	uk[I][J] = NewtonPoly (z, f, t);

      }


  double ujk[int_ord + 1];

  for (size_t I = 0; I <= int_ord; I++)
    {


      for (size_t J = 0; J <= int_ord; J++)
	{

	  f[J] = uk[I][J];

	  t[J] = Get_y (ii + I, ji + J, ki);

	}

      ujk[I] = NewtonPoly (y, f, t);

    }




  for (size_t I = 0; I <= int_ord; I++)
    {

      f[I] = ujk[I];

      t[I] = Get_x (ii + I, ji, ki);

    }



  return (NewtonPoly (x, f, t));


}





double
ffunction::NewtonPoly (const double t, valarray < double >F,
		      valarray < double >x)
{



  size_t ord = F.size ();


  if (ord != 3)

    exit (1);



  double dx0 = (t - x[0]);

  double dx1 = (t - x[1]);



  double Fx0x1 = (F[1] - F[0]) / (x[1] - x[0]);

  double Fx1x2 = (F[2] - F[1]) / (x[2] - x[1]);

  double Fx0x1x2 = (Fx1x2 - Fx0x1) / (x[2] - x[0]);


  return (F[0] + Fx0x1 * dx0 + Fx0x1x2 * dx0 * dx1);


}






double
ffunction::Hermite_Interpol_In (double x, double y, double z)
{




  size_t i = D->Get_i (x);

  size_t j = D->Get_j (y);

  size_t k = D->Get_k (z);




  size_t ii = i == 0 ? 0 : i - 1;

  size_t ji = j == 0 ? 0 : j - 1;

  size_t ki = k == 0 ? 0 : k - 1;


  size_t nx = D->Get_Nx ();

  size_t ny = D->Get_Ny ();

  size_t nz = D->Get_Nz ();



  valarray < double >f (2);

  valarray < double >df (2);

  valarray < double >t (2);


  size_t Ii[4],
    Ji[4],
    Ki[4];

  if (ii > 0 && ii < nx - 2)
    {


      Ii[0] = ii - 1;

      Ii[1] = ii;

      Ii[2] = ii + 1;

      Ii[3] = ii + 2;


    }
  else
   if (ii == 0)
    {


      Ii[0] = ii;

      Ii[1] = ii + 1;

      Ii[2] = ii + 1;

      Ii[3] = ii + 2;


    }


  else
    {


      Ii[0] = ii - 1;

      Ii[1] = ii;

      Ii[2] = ii;

      Ii[3] = ii + 1;


    }








  if (ji > 0 && ji < ny - 2)
    {


      Ji[0] = ji - 1;

      Ji[1] = ji;

      Ji[2] = ji + 1;

      Ji[3] = ji + 2;


    }
  else
   if (ji == 0)
    {


      Ji[0] = ji;

      Ji[1] = ji + 1;

      Ji[2] = ji + 1;

      Ji[3] = ji + 2;


    }


  else
    {


      Ji[0] = ji - 1;

      Ji[1] = ji;

      Ji[2] = ji;

      Ji[3] = ji + 1;


    }





  if (ki > 0 && ki < nz - 2)
    {


      Ki[0] = ki - 1;

      Ki[1] = ki;

      Ki[2] = ki + 1;

      Ki[3] = ki + 2;


    }
  else
   if (ki == 0)
    {


      Ki[0] = ki;

      Ki[1] = ki + 1;

      Ki[2] = ki + 1;

      Ki[3] = ki + 2;


    }


  else
    {


      Ki[0] = ki - 1;

      Ki[1] = ki;

      Ki[2] = ki;

      Ki[3] = ki + 1;


    }





  double uk[4][4];

  for (size_t I = 0; I < 4; I++)

    for (size_t J = 0; J < 4; J++)
      {

	for (size_t K = 0; K < 2; K++)
	  {



	    f[K] = Get_val (Ii[I], Ji[J], ki + K);


	    if (ki + K - 1 >= 0 && ki + K + 1 < nz)

	      df[K] =
		dzc2 (Get_val (Ii[I], Ji[I], Ki[2 * K]),
		      Get_val (Ii[I], Ji[J], Ki[2 * K + 1]));

	    else

	      df[K] =
		dy1 (Get_val (Ii[I], Ji[J], Ki[2 * K]),
		     Get_val (Ii[I], Ji[J], Ki[2 * K + 1]));


	    t[K] = Get_z (ii, ji, ki + K);

	  }

	uk[I][J] = Hermite_Poly (z, t[0], t[1], dz, f[0], f[1], df[0], df[1]);


      }



  double ujk[4];

  for (size_t I = 0; I < 4; I++)
    {


      for (size_t J = 0; J < 2; J++)
	{

	  f[J] = uk[I][J + 1];

	  if (ji > 0 && ji < ny - 2)

	    df[J] = dyc2 (uk[I][J], uk[I][J + 2]);

	  else

	    df[J] = dy1 (uk[I][J], uk[I][J + 2]);



	  t[J] = Get_y (ii, ji + J, ki);

	}

      ujk[I] = Hermite_Poly (y, t[0], t[1], dy, f[0], f[1], df[0], df[1]);

    }




  for (size_t I = 0; I < 2; I++)
    {

      f[I] = ujk[I + 1];

      if (ii > 0 && ii < nx - 2)

	df[I] = dxc2 (ujk[I], ujk[I + 2]);

      else

	df[I] = dx1 (ujk[I], ujk[I + 2]);


      t[I] = Get_x (ii + I, ji, ki);

    }



  return (Hermite_Poly (x, t[0], t[1], dx, f[0], f[1], df[0], df[1]));





}







double
ffunction::Hermite_Poly (const double t,
			const double ti,
			const double tip1,
			const double h,
			const double p0,
			const double p1, const double m0, const double m1)
{


  return ((1 + 2 * (t - ti) / h) * pow ((t - tip1) / h, 2) * p0 +
	  (1 - 2 * (t - tip1) / h) * pow ((t - ti) / h, 2) * p1 +
	  (t - ti) * pow ((t - tip1) / h, 2) * m0 +
	  (t - tip1) * pow ((t - ti) / h, 2) * m1);


}











double
ffunction::Average_Interpol_In (double x, double y, double z)
{




  size_t i = D->Get_i (x);

  size_t j = D->Get_j (y);

  size_t k = D->Get_k (z);


  size_t nx = D->Get_Nx ();

  size_t ny = D->Get_Ny ();

  size_t nz = D->Get_Nz ();


  size_t im1 = i == 0 || i >= nx - 1 ? i : i - 1;

  size_t jm1 = j == 0 || j >= ny - 1 ? j : j - 1;

  size_t km1 = k == 0 || k >= nz - 1 ? k : k - 1;

  size_t ip1 = i == 0 || i >= nx - 1 ? i : i + 1;

  size_t jp1 = j == 0 || j >= ny - 1 ? j : j + 1;

  size_t kp1 = k == 0 || k >= nz - 1 ? k : k + 1;




  return (0.5 * Get_val (i, j, k) +
	  (Get_val (im1, j, k) +
	   Get_val (i, jm1, k) +
	   Get_val (i, j, km1) +
	   Get_val (ip1, j, k) +
	   Get_val (i, jp1, k) + Get_val (i, j, kp1)) / 12.0);








}
