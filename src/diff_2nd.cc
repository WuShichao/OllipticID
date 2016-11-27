#include "ffunction.h"




double
ffunction::dxl2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int im = Get_Index (i - 1, j, k);

  int imm = Get_Index (i - 2, j, k);

  return ((3.0 * Get_val () - 4.0 * Get_val (im) +
	   Get_val (imm)) * 0.5 * idx);

}



double
ffunction::dyl2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int im = Get_Index (i, j - 1, k);

  int imm = Get_Index (i, j - 2, k);

  return ((3.0 * Get_val () - 4.0 * Get_val (im) +
	   Get_val (imm)) * 0.5 * idy);

}




double
ffunction::dzl2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int im = Get_Index (i, j, k - 1);

  int imm = Get_Index (i, j, k - 2);

  return ((3.0 * Get_val () - 4.0 * Get_val (im) +
	   Get_val (imm)) * 0.5 * idz);

}





double
ffunction::dxc2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i + 1, j, k);

  int im = Get_Index (i - 1, j, k);

  return ((Get_val (ip) - Get_val (im)) * 0.5 * idx);

}


double
ffunction::dyc2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i, j + 1, k);

  int im = Get_Index (i, j - 1, k);

  return ((Get_val (ip) - Get_val (im)) * 0.5 * idy);

}






double
ffunction::dzc2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i, j, k + 1);

  int im = Get_Index (i, j, k - 1);



  return ((Get_val (ip) - Get_val (im)) * 0.5 * idz);


}






double
ffunction::dxr2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i + 1, j, k);

  int ipp = Get_Index (i + 2, j, k);

  return ((-3.0 * Get_val () + 4.0 * Get_val (ip) -
	   Get_val (ipp)) * 0.5 * idx);


}






double
ffunction::dyr2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i, j + 1, k);

  int ipp = Get_Index (i, j + 2, k);

  return ((-3.0 * Get_val () + 4.0 * Get_val (ip) -
	   Get_val (ipp)) * 0.5 * idy);


}






double
ffunction::dzr2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i, j, k + 1);

  int ipp = Get_Index (i, j, k + 2);

  return ((-3.0 * Get_val () + 4.0 * Get_val (ip) -
	   Get_val (ipp)) * 0.5 * idz);


}










double
ffunction::dx2 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i == 0)
    {

      int ip1 = Get_Index (i + 1, j, k);

      int ip2 = Get_Index (i + 2, j, k);


      return ((-3.0 * Get_val () + 4.0 * Get_val (ip1) -
	       Get_val (ip2)) * 0.5 * idx);

    }
  else
   if (i == nx - 1)
    {


      int im1 = Get_Index (i - 1, j, k);

      int im2 = Get_Index (i - 2, j, k);


      return ((3.0 * Get_val () - 4.0 * Get_val (im1) +
	       Get_val (im2)) * 0.5 * idx);


    }

  else
    {

      int ip1 = Get_Index (i + 1, j, k);

      int im1 = Get_Index (i - 1, j, k);

      return ((Get_val (ip1) - Get_val (im1)) * 0.5 * idx);


    }



}







double
ffunction::dudx2 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i == 0)

    return (-3.0 * 0.5 * idx);


  else
   if (i == nx - 1)

    return (3.0 * 0.5 * idx);

  else

    return (0.0);



}












double
ffunction::dy2 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j == 0)
    {

      int ip1 = Get_Index (i, j + 1, k);

      int ip2 = Get_Index (i, j + 2, k);


      return ((-3.0 * Get_val () + 4.0 * Get_val (ip1) -
	       Get_val (ip2)) * 0.5 * idy);

    }
  else
   if (j == ny - 1)
    {


      int im1 = Get_Index (i, j - 1, k);

      int im2 = Get_Index (i, j - 2, k);


      return ((3.0 * Get_val () - 4.0 * Get_val (im1) +
	       Get_val (im2)) * 0.5 * idy);


    }

  else
    {

      int ip1 = Get_Index (i, j + 1, k);

      int im1 = Get_Index (i, j - 1, k);

      return ((Get_val (ip1) - Get_val (im1)) * 0.5 * idy);


    }



}







double
ffunction::dudy2 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j == 0)

    return (-1.5 * idy);


  else
   if (j == ny - 1)

    return (1.5 * idy);

  else

    return (0.0);



}












double
ffunction::dz2 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k == 0)
    {

      int ip1 = Get_Index (i, j, k + 1);

      int ip2 = Get_Index (i, j, k + 2);


      return ((-3.0 * Get_val () + 4.0 * Get_val (ip1) -
	       Get_val (ip2)) * 0.5 * idz);

    }
  else
   if (k == nz - 1)
    {


      int im1 = Get_Index (i, j, k - 1);

      int im2 = Get_Index (i, j, k - 2);


      return ((3.0 * Get_val () - 4.0 * Get_val (im1) +
	       Get_val (im2)) * 0.5 * idz);


    }

  else
    {

      int ip1 = Get_Index (i, j, k + 1);

      int im1 = Get_Index (i, j, k - 1);

      return ((Get_val (ip1) - Get_val (im1)) * 0.5 * idz);


    }



}







double
ffunction::dudz2 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k == 0)

    return (-1.5 * idz);


  else
   if (k == nz - 1)

    return (1.5 * idz);

  else

    return (0.0);



}


















// *******************    Segunda derivada    *****************************************



double
ffunction::ddxc2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i + 1, j, k);

  int im = Get_Index (i - 1, j, k);

  return ((Get_val (ip) - 2.0 * Get_val () + Get_val (im)) * idx2);

}





double
ffunction::ddxl2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int im1 = Get_Index (i - 1, j, k);

  int im2 = Get_Index (i - 2, j, k);

  return ((Get_val () - 2.0 * Get_val (im1) + Get_val (im2)) * idx2);

}



double
ffunction::ddxr2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip1 = Get_Index (i + 1, j, k);

  int ip2 = Get_Index (i + 2, j, k);

  return ((Get_val () - 2.0 * Get_val (ip1) + Get_val (ip2)) * idx2);

}



double
ffunction::ddyc2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i, j + 1, k);

  int im = Get_Index (i, j - 1, k);


  return ((Get_val (ip) - 2.0 * Get_val () + Get_val (im)) * idy2);

}



double
ffunction::ddzc2 ()
{



  int i,
    j,
    k;


  Get_MIndex (i, j, k);

  int ip = Get_Index (i, j, k + 1);

  int im = Get_Index (i, j, k - 1);

  return ((Get_val (ip) - 2.0 * Get_val () + Get_val (im)) * idz2);

}



double
ffunction::LapS (int i0, int j0, int k0)
{


  int i,
    j,
    k;

  Get_MIndex (i, j, k);


  double lapi =
    i == i0 ? Get_val () - 2.0 * Get_val (i + 1, j, k) + Get_val (i + 2, j,
								  k) :
    Get_val (i - 2, j, k) - 2.0 * Get_val (i - 1, j, k) + Get_val ();

  double lapj =
    j == j0 ? Get_val () - 2.0 * Get_val (i, j + 1, k) + Get_val (i, j + 2,
								  k) :
    Get_val (i, j - 2, k) - 2.0 * Get_val (i, j - 1, k) + Get_val ();

  double lapk =
    k == k0 ? Get_val () - 2.0 * Get_val (i, j, k + 1) + Get_val (i, j,
								  k +
								  2) :
    Get_val (i, j, k - 2) - 2.0 * Get_val (i, j, k - 1) + Get_val ();


/*
    double x = Get_x();
    double y = Get_y();
    double z = Get_z();
    
    double lapi =  Interpol_In(x-0.5*dx,y,z) - 2.0 * Get_val() + Interpol_In(x+0.5*dx,y,z);

    double lapj =  Interpol_In(x,y-0.5*dy,z) - 2.0 * Get_val() + Interpol_In(x,y+0.5*dy,z);

    double lapk =  Interpol_In(x,y,z-0.5*dz) - 2.0 * Get_val() + Interpol_In(x,y,z+0.5*dz);

*/

  return (lapi * idx2 + lapj * idy2 + lapk * idz2);



}









double
ffunction::LapC2 ()
{




  int i,
    j,
    k;

  Get_MIndex (i, j, k);

  double x = Get_x ();

  double y = Get_y ();

  double z = Get_z ();

/*
    double iL2 = 1/(4);
    
    return(

        iL2*( (x*x-1) * ( (x*x-1) * ( Get_val(i+1,j,k) - 2.0 * Get_val() + Get_val(i-1,j,k) ) * idx2 -
                       2 * x * ( Get_val(i+1,j,k) - Get_val(i-1,j,k) ) * 0.5 * idx) +

          (y*y-1) * ( (y*y-1) * ( Get_val(i,j+1,k) - 2.0 * Get_val() + Get_val(i,j-1,k) ) * idy2 -
                       2 * y * ( Get_val(i,j+1,k) - Get_val(i,j-1,k) ) * 0.5 * idy) +

          (z*z-1) * ( (z*z-1) * ( Get_val(i,j,k+1) - 2.0 * Get_val() + Get_val(i,j,k-1) ) * idz2 -
                       2 * z * ( Get_val(i,j,k+1) - Get_val(i,j,k-1) ) * 0.5 * idz) ) 
        

           );

*/

/*    

    return(

        4.0 * ( pow(cos(0.5*pi*x),4) * ( ( Get_val(i+1,j,k) - 2.0 * Get_val() + Get_val(i-1,j,k) ) * idx2 -
                                         pi * tan(0.5 * pi * x) * ( Get_val(i+1,j,k) - Get_val(i-1,j,k) ) * 0.5 * idx) +
                
                pow(cos(0.5*pi*y),4) * ( ( Get_val(i,j+1,k) - 2.0 * Get_val() + Get_val(i,j-1,k) ) * idy2 -
                                         pi * tan(0.5 * pi * y) * ( Get_val(i,j+1,k) - Get_val(i,j-1,k) ) * 0.5 * idy) +

                pow(cos(0.5*pi*z),4) * ( ( Get_val(i,j,k+1) - 2.0 * Get_val() + Get_val(i,j,k-1) ) * idz2 -
                                         pi * tan(0.5 * pi * z) * ( Get_val(i,j,k+1) - Get_val(i,j,k-1) ) * 0.5 * idz)
                ) / (pi*pi)

           );
*/

  double du_x =
    0.125 * idx * (Get_val (i + 1, j + 1, k + 1) -
		   Get_val (i - 1, j + 1, k + 1) + Get_val (i + 1, j + 1,
							    k - 1) -
		   Get_val (i - 1, j + 1, k - 1) + Get_val (i + 1, j - 1,
							    k + 1) -
		   Get_val (i - 1, j - 1, k + 1) + Get_val (i + 1, j - 1,
							    k - 1) -
		   Get_val (i - 1, j - 1, k - 1));

  double du_y =
    0.125 * idy * (Get_val (i + 1, j + 1, k + 1) -
		   Get_val (i + 1, j - 1, k + 1) + Get_val (i + 1, j + 1,
							    k - 1) -
		   Get_val (i + 1, j - 1, k - 1) + Get_val (i - 1, j + 1,
							    k + 1) -
		   Get_val (i - 1, j - 1, k + 1) + Get_val (i - 1, j + 1,
							    k - 1) -
		   Get_val (i - 1, j - 1, k - 1));

  double du_z =
    0.125 * idz * (Get_val (i + 1, j + 1, k + 1) -
		   Get_val (i + 1, j + 1, k - 1) + Get_val (i - 1, j + 1,
							    k + 1) -
		   Get_val (i - 1, j + 1, k - 1) + Get_val (i + 1, j - 1,
							    k + 1) -
		   Get_val (i + 1, j - 1, k - 1) + Get_val (i - 1, j - 1,
							    k + 1) -
		   Get_val (i - 1, j - 1, k - 1));

  return (4.0 *
	  (pow (cos (0.5 * pi * x), 4) *
	   ((Get_val (i + 1, j, k) - 2.0 * Get_val () +
	     Get_val (i - 1, j,
		      k)) * idx2 - pi * tan (0.5 * pi * x) * du_x) +
	   pow (cos (0.5 * pi * y),
		4) * ((Get_val (i, j + 1, k) - 2.0 * Get_val () +
		       Get_val (i, j - 1,
				k)) * idy2 - pi * tan (0.5 * pi * y) * du_y) +
	   pow (cos (0.5 * pi * z),
		4) * ((Get_val (i, j, k + 1) - 2.0 * Get_val () +
		       Get_val (i, j,
				k - 1)) * idz2 -
		      pi * tan (0.5 * pi * z) * du_z)) / (pi * pi));



}



double
ffunction::duLapC2 ()
{




  int i,
    j,
    k;

  Get_MIndex (i, j, k);

  double x = Get_x ();

  double y = Get_y ();

  double z = Get_z ();


  return (-8.0 * (pow (cos (0.5 * pi * x), 4) * idx2 +
		  pow (cos (0.5 * pi * y), 4) * idy2 +
		  pow (cos (0.5 * pi * z), 4) * idz2) / (pi * pi));


/*
    double iL2 = 1/(4);

    return(

        iL2*( (x*x-1) * ( (x*x-1) * ( - 2.0 ) * idx2 ) +
          
          (y*y-1) * ( (y*y-1) * ( - 2.0 ) * idy2 ) +

          (z*z-1) * ( (z*z-1) * ( - 2.0 ) * idz2 ) ) 
        

           );

*/

}



double
ffunction::Lapc2 (ffunction & u)
{




  double x = Get_x ();

  double y = Get_y ();

  double z = Get_z ();



  int i,
    j,
    k;

  Get_MIndex (i, j, k);



  size_t nx = D->Get_Nx ();

  double uip =
    (size_t) i < nx - 1 ? Get_val (i + 1, j, k) : u.Interpol_In (x + dx, y,
								 z);

  double uim =
    (size_t) i > 0 ? Get_val (i - 1, j, k) : u.Interpol_In (x - dx, y, z);


  size_t ny = D->Get_Ny ();

  double ujp =
    (size_t) j < ny - 1 ? Get_val (i, j + 1, k) : u.Interpol_In (x, y + dy,
								 z);

  double ujm =
    (size_t) j > 0 ? Get_val (i, j - 1, k) : u.Interpol_In (x, y - dy, z);


  size_t nz = D->Get_Nz ();

  double ukp =
    (size_t) k < nz - 1 ? Get_val (i, j, k + 1) : u.Interpol_In (x, y,
								 z + dz);

  double ukm =
    (size_t) k > 0 ? Get_val (i, j, k - 1) : u.Interpol_In (x, y, z - dz);


  cout << ((uip - 2.0 * Get_val () + uim) * idx2
	  + (ujp - 2.0 * Get_val () + ujm) * idy2
	   + (ukp - 2.0 * Get_val () + ukm) * idz2);


  return ((uip - 2.0 * Get_val () + uim) * idx2
	  + (ujp - 2.0 * Get_val () + ujm) * idy2
	  + (ukp - 2.0 * Get_val () + ukm) * idz2);



}














double
ffunction::dx3 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i == 0)
    {

      int ip = Get_Index (i + 1, j, k);

      int ipp = Get_Index (i + 2, j, k);

      int ippp = Get_Index (i + 3, j, k);


      return (-(11 * Get_val ()
		- 18 * Get_val (ip)
		+ 9 * Get_val (ipp) - 2 * Get_val (ippp)) * idx / 6.0);
    }
  else
   if (i >= 1 && i <= nx / 2)
    {

      int im = Get_Index (i - 1, j, k);

      int ip = Get_Index (i + 1, j, k);

      int ipp = Get_Index (i + 2, j, k);



      return (-(2 * Get_val (im)
		+ 3 * Get_val ()
		- 6 * Get_val (ip) + Get_val (ipp)) * idx / 6.0);
    }
  else
   if (i < nx - 1)
    {

      int im = Get_Index (i - 1, j, k);

      int imm = Get_Index (i - 2, j, k);

      int ip = Get_Index (i + 1, j, k);

      return ((2 * Get_val (ip)
	       + 3 * Get_val ()
	       - 6 * Get_val (im) + Get_val (imm)) * idx / 6.0);

    }
  else
   if (i == nx - 1)
    {


      int im = Get_Index (i - 1, j, k);

      int imm = Get_Index (i - 2, j, k);

      int immm = Get_Index (i - 3, j, k);


      return ((11 * Get_val ()
	       - 18 * Get_val (im)
	       + 9 * Get_val (imm) - 2 * Get_val (immm)) * idx / 6.0);

    }


  return (0);

}







double
ffunction::dudx3 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i == 0)

    return (-11 * idx / 6.0);

  else
   if (i >= 1 && i <= nx / 2)

    return (-0.5 * idx);


  else
   if (i < nx - 1)

    return (0.5 * idx);

  else
   if (i == nx - 1)


    return (11 * idx / 6.0);


  else

    exit (0);


}







double
ffunction::dy3 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j == 0)
    {

      int ip = Get_Index (i, j + 1, k);

      int ipp = Get_Index (i, j + 2, k);

      int ippp = Get_Index (i, j + 3, k);


      return (-(11 * Get_val ()
		- 18 * Get_val (ip)
		+ 9 * Get_val (ipp) - 2 * Get_val (ippp)) * idy / 6.0);
    }
  else
   if (j >= 1 && j <= ny / 2)
    {

      int im = Get_Index (i, j - 1, k);

      int ip = Get_Index (i, j + 1, k);

      int ipp = Get_Index (i, j + 2, k);



      return (-(2 * Get_val (im)
		+ 3 * Get_val ()
		- 6 * Get_val (ip) + Get_val (ipp)) * idy / 6.0);
    }
  else
   if (j < ny - 1)
    {

      int im = Get_Index (i, j - 1, k);

      int imm = Get_Index (i, j - 2, k);

      int ip = Get_Index (i, j + 1, k);

      return ((2 * Get_val (ip)
	       + 3 * Get_val ()
	       - 6 * Get_val (im) + Get_val (imm)) * idy / 6.0);

    }
  else
   if (j == ny - 1)
    {


      int im = Get_Index (i, j - 1, k);

      int imm = Get_Index (i, j - 2, k);

      int immm = Get_Index (i, j - 3, k);


      return ((11 * Get_val ()
	       - 18 * Get_val (im)
	       + 9 * Get_val (imm) - 2 * Get_val (immm)) * idy / 6.0);

    }


  return (0);

}







double
ffunction::dudy3 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j == 0)

    return (-11 * idy / 6.0);

  else
   if (j >= 1 && j <= ny / 2)

    return (-0.5 * idy);


  else
   if (j < ny - 1)

    return (0.5 * idy);

  else
   if (j == ny - 1)


    return (11 * idy / 6.0);


  else

    exit (0);


}













double
ffunction::dz3 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k == 0)
    {

      int ip = Get_Index (i, j, k + 1);

      int ipp = Get_Index (i, j, k + 2);

      int ippp = Get_Index (i, j, k + 3);


      return (-(11 * Get_val ()
		- 18 * Get_val (ip)
		+ 9 * Get_val (ipp) - 2 * Get_val (ippp)) * idz / 6.0);
    }
  else
   if (k >= 1 && k <= nz / 2)
    {

      int im = Get_Index (i, j, k - 1);

      int ip = Get_Index (i, j, k + 1);

      int ipp = Get_Index (i, j, k + 2);



      return (-(2 * Get_val (im)
		+ 3 * Get_val ()
		- 6 * Get_val (ip) + Get_val (ipp)) * idz / 6.0);
    }
  else
   if (k < nz - 1)
    {

      int im = Get_Index (i, j, k - 1);

      int imm = Get_Index (i, j, k - 2);

      int ip = Get_Index (i, j, k + 1);

      return ((2 * Get_val (ip)
	       + 3 * Get_val ()
	       - 6 * Get_val (im) + Get_val (imm)) * idz / 6.0);

    }
  else
   if (k == nz - 1)
    {


      int im = Get_Index (i, j, k - 1);

      int imm = Get_Index (i, j, k - 2);

      int immm = Get_Index (i, j, k - 3);


      return ((11 * Get_val ()
	       - 18 * Get_val (im)
	       + 9 * Get_val (imm) - 2 * Get_val (immm)) * idz / 6.0);

    }


  return (0);

}







double
ffunction::dudz3 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k == 0)

    return (-11 * idz / 6.0);

  else
   if (k >= 1 && k <= nz / 2)

    return (-0.5 * idz);


  else
   if (k < nz - 1)

    return (0.5 * idz);

  else
   if (k == nz - 1)


    return (11 * idz / 6.0);


  else

    exit (0);


}






//== Mixed derivatives ==//











double
ffunction::ddxyc2 ()
{


  int i,j, k;


  Get_MIndex (i, j, k);


  return (0.25 *
	  (Get_val (i + 1, j + 1, k) - Get_val (i + 1, j - 1, k) -
	   Get_val (i - 1, j + 1, k) + Get_val (i - 1, j - 1, k) ) * idx * idy);

}





double
ffunction::ddxzc2 ()
{


 int i,
    j,
    k;


  Get_MIndex (i, j, k);


  return (0.25 *
	  (Get_val (i + 1, j, k + 1) - Get_val (i + 1, j, k - 1) -
	   Get_val (i - 1, j, k + 1) + Get_val (i - 1, j,
						k - 1)) * idx * idz);



}




double
ffunction::ddyzc2 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  return (0.25 *
	  (Get_val (i, j + 1, k + 1) - Get_val (i, j + 1, k - 1) -
	   Get_val (i, j - 1, k + 1) + Get_val (i, j - 1,
						k - 1)) * idy * idz);

}
