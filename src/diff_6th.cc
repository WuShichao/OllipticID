#include "ffunction.h"




// X


double
ffunction::dx6 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i <= nx / 2)
    {

      int ip1 = Get_Index (i + 1, j, k);

      int ip2 = Get_Index (i + 2, j, k);

      int ip3 = Get_Index (i + 3, j, k);

      int ip4 = Get_Index (i + 4, j, k);

      int ip5 = Get_Index (i + 5, j, k);

      int ip6 = Get_Index (i + 6, j, k);

      return (-(147 * Get_val ()
		- 360 * Get_val (ip1)
		+ 450 * Get_val (ip2)
		- 400 * Get_val (ip3)
		+ 225 * Get_val (ip4)
		- 72 * Get_val (ip5) + 10 * Get_val (ip6)) * idx / 60.0);
    }
  else

    {

      int im1 = Get_Index (i - 1, j, k);

      int im2 = Get_Index (i - 2, j, k);

      int im3 = Get_Index (i - 3, j, k);

      int im4 = Get_Index (i - 4, j, k);

      int im5 = Get_Index (i - 5, j, k);

      int im6 = Get_Index (i - 6, j, k);

      return ((147 * Get_val ()
	       - 360 * Get_val (im1)
	       + 450 * Get_val (im2)
	       - 400 * Get_val (im3)
	       + 225 * Get_val (im4)
	       - 72 * Get_val (im5) + 10 * Get_val (im6)) * idx / 60.0);
    }


}







double
ffunction::dudx6 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i <= nx / 2)

    return (-147 * idx / 60.0);

  else

    return (147 * idx / 60.0);


}







//  Y





double
ffunction::dy6 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j <= ny / 2)
    {

      int ip1 = Get_Index (i, j + 1, k);

      int ip2 = Get_Index (i, j + 2, k);

      int ip3 = Get_Index (i, j + 3, k);

      int ip4 = Get_Index (i, j + 4, k);

      int ip5 = Get_Index (i, j + 5, k);

      int ip6 = Get_Index (i, j + 6, k);

      return (-(147 * Get_val ()
		- 360 * Get_val (ip1)
		+ 450 * Get_val (ip2)
		- 400 * Get_val (ip3)
		+ 225 * Get_val (ip4)
		- 72 * Get_val (ip5) + 10 * Get_val (ip6)) * idy / 60.0);
    }
  else

    {

      int im1 = Get_Index (i, j - 1, k);

      int im2 = Get_Index (i, j - 2, k);

      int im3 = Get_Index (i, j - 3, k);

      int im4 = Get_Index (i, j - 4, k);

      int im5 = Get_Index (i, j - 5, k);

      int im6 = Get_Index (i, j - 6, k);

      return ((147 * Get_val ()
	       - 360 * Get_val (im1)
	       + 450 * Get_val (im2)
	       - 400 * Get_val (im3)
	       + 225 * Get_val (im4)
	       - 72 * Get_val (im5) + 10 * Get_val (im6)) * idy / 60.0);
    }


}







double
ffunction::dudy6 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j <= ny / 2)

    return (-147 * idy / 60.0);

  else

    return (147 * idy / 60.0);


}











//  Z





double
ffunction::dz6 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k <= nz / 2)
    {

      int ip1 = Get_Index (i, j, k + 1);

      int ip2 = Get_Index (i, j, k + 2);

      int ip3 = Get_Index (i, j, k + 3);

      int ip4 = Get_Index (i, j, k + 4);

      int ip5 = Get_Index (i, j, k + 5);

      int ip6 = Get_Index (i, j, k + 6);

      return (-(147 * Get_val ()
		- 360 * Get_val (ip1)
		+ 450 * Get_val (ip2)
		- 400 * Get_val (ip3)
		+ 225 * Get_val (ip4)
		- 72 * Get_val (ip5) + 10 * Get_val (ip6)) * idz / 60.0);
    }
  else

    {

      int im1 = Get_Index (i, j, k - 1);

      int im2 = Get_Index (i, j, k - 2);

      int im3 = Get_Index (i, j, k - 3);

      int im4 = Get_Index (i, j, k - 4);

      int im5 = Get_Index (i, j, k - 5);

      int im6 = Get_Index (i, j, k - 6);

      return ((147 * Get_val ()
	       - 360 * Get_val (im1)
	       + 450 * Get_val (im2)
	       - 400 * Get_val (im3)
	       + 225 * Get_val (im4)
	       - 72 * Get_val (im5) + 10 * Get_val (im6)) * idz / 60.0);
    }


}







double
ffunction::dudz6 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k <= nz / 2)

    return (-147 * idz / 60.0);

  else

    return (147 * idz / 60.0);


}


















// *******************    Segunda derivada    *****************************************

















double
ffunction::ddxc6 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip3 = Get_Index (i + 3, j, k);

  int ip2 = Get_Index (i + 2, j, k);

  int ip1 = Get_Index (i + 1, j, k);

  int im1 = Get_Index (i - 1, j, k);

  int im2 = Get_Index (i - 2, j, k);

  int im3 = Get_Index (i - 3, j, k);


  return ((2.0 * (Get_val (ip3) + Get_val (im3))
	   - 27.0 * (Get_val (ip2) + Get_val (im2)
		     - 10.0 * (Get_val (ip1) + Get_val (im1)))
	   - 490.0 * Get_val ()) * idx2 / 180.0);

}





double
ffunction::ddyc6 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip3 = Get_Index (i, j + 3, k);

  int ip2 = Get_Index (i, j + 2, k);

  int ip1 = Get_Index (i, j + 1, k);

  int im1 = Get_Index (i, j - 1, k);

  int im2 = Get_Index (i, j - 2, k);

  int im3 = Get_Index (i, j - 3, k);


  return ((2.0 * (Get_val (ip3) + Get_val (im3))
	   - 27.0 * (Get_val (ip2) + Get_val (im2)
		     - 10.0 * (Get_val (ip1) + Get_val (im1)))
	   - 490.0 * Get_val ()) * idy2 / 180.0);

}



double
ffunction::ddzc6 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip3 = Get_Index (i, j, k + 3);

  int ip2 = Get_Index (i, j, k + 2);

  int ip1 = Get_Index (i, j, k + 1);

  int im1 = Get_Index (i, j, k - 1);

  int im2 = Get_Index (i, j, k - 2);

  int im3 = Get_Index (i, j, k - 3);


  return ((2.0 * (Get_val (ip3) + Get_val (im3))
	   - 27.0 * (Get_val (ip2) + Get_val (im2)
		     - 10.0 * (Get_val (ip1) + Get_val (im1)))
	   - 490.0 * Get_val ()) * idz2 / 180.0);

}
