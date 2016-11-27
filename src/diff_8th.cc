#include "ffunction.h"






//  X




double
ffunction::dx8 ()
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

      int ip7 = Get_Index (i + 7, j, k);

      int ip8 = Get_Index (i + 8, j, k);


      return (-(2283 * Get_val ()
		- 6720 * Get_val (ip1)
		+ 11760 * Get_val (ip2)
		- 15680 * Get_val (ip3)
		+ 14700 * Get_val (ip4)
		- 9408 * Get_val (ip5)
		+ 3920 * Get_val (ip6)
		- 960 * Get_val (ip7) + 105 * Get_val (ip8)) * idx / 840.0);


    }
  else

    {

      int im1 = Get_Index (i - 1, j, k);

      int im2 = Get_Index (i - 2, j, k);

      int im3 = Get_Index (i - 3, j, k);

      int im4 = Get_Index (i - 4, j, k);

      int im5 = Get_Index (i - 5, j, k);

      int im6 = Get_Index (i - 6, j, k);

      int im7 = Get_Index (i - 7, j, k);

      int im8 = Get_Index (i - 8, j, k);


      return ((2283 * Get_val ()
	       - 6720 * Get_val (im1)
	       + 11760 * Get_val (im2)
	       - 15680 * Get_val (im3)
	       + 14700 * Get_val (im4)
	       - 9408 * Get_val (im5)
	       + 3920 * Get_val (im6)
	       - 960 * Get_val (im7) + 105 * Get_val (im8)) * idx / 840.0);



    }






}











double
ffunction::dudx8 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i <= nx / 2)

    return (-2283 * idx / 840.0);

  else

    return (2283 * idx / 840.0);





}











//  Y




double
ffunction::dy8 ()
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

      int ip7 = Get_Index (i, j + 7, k);

      int ip8 = Get_Index (i, j + 8, k);


      return (-(2283 * Get_val ()
		- 6720 * Get_val (ip1)
		+ 11760 * Get_val (ip2)
		- 15680 * Get_val (ip3)
		+ 14700 * Get_val (ip4)
		- 9408 * Get_val (ip5)
		+ 3920 * Get_val (ip6)
		- 960 * Get_val (ip7) + 105 * Get_val (ip8)) * idy / 840.0);


    }
  else

    {

      int im1 = Get_Index (i, j - 1, k);

      int im2 = Get_Index (i, j - 2, k);

      int im3 = Get_Index (i, j - 3, k);

      int im4 = Get_Index (i, j - 4, k);

      int im5 = Get_Index (i, j - 5, k);

      int im6 = Get_Index (i, j - 6, k);

      int im7 = Get_Index (i, j - 7, k);

      int im8 = Get_Index (i, j - 8, k);


      return ((2283 * Get_val ()
	       - 6720 * Get_val (im1)
	       + 11760 * Get_val (im2)
	       - 15680 * Get_val (im3)
	       + 14700 * Get_val (im4)
	       - 9408 * Get_val (im5)
	       + 3920 * Get_val (im6)
	       - 960 * Get_val (im7) + 105 * Get_val (im8)) * idy / 840.0);



    }






}











double
ffunction::dudy8 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j <= ny / 2)

    return (-2283 * idy / 840.0);

  else

    return (2283 * idy / 840.0);





}




















//  Z




double
ffunction::dz8 ()
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

      int ip7 = Get_Index (i, j, k + 7);

      int ip8 = Get_Index (i, j, k + 8);


      return (-(2283 * Get_val ()
		- 6720 * Get_val (ip1)
		+ 11760 * Get_val (ip2)
		- 15680 * Get_val (ip3)
		+ 14700 * Get_val (ip4)
		- 9408 * Get_val (ip5)
		+ 3920 * Get_val (ip6)
		- 960 * Get_val (ip7) + 105 * Get_val (ip8)) * idz / 840.0);


    }
  else

    {

      int im1 = Get_Index (i, j, k - 1);

      int im2 = Get_Index (i, j, k - 2);

      int im3 = Get_Index (i, j, k - 3);

      int im4 = Get_Index (i, j, k - 4);

      int im5 = Get_Index (i, j, k - 5);

      int im6 = Get_Index (i, j, k - 6);

      int im7 = Get_Index (i, j, k - 7);

      int im8 = Get_Index (i, j, k - 8);


      return ((2283 * Get_val ()
	       - 6720 * Get_val (im1)
	       + 11760 * Get_val (im2)
	       - 15680 * Get_val (im3)
	       + 14700 * Get_val (im4)
	       - 9408 * Get_val (im5)
	       + 3920 * Get_val (im6)
	       - 960 * Get_val (im7) + 105 * Get_val (im8)) * idz / 840.0);



    }






}











double
ffunction::dudz8 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k <= nz / 2)

    return (-2283 * idz / 840.0);

  else

    return (2283 * idz / 840.0);





}




























//  *******************    Segunda derivada    *************************








double
ffunction::ddxc8 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip4 = Get_Index (i + 4, j, k);

  int ip3 = Get_Index (i + 3, j, k);

  int ip2 = Get_Index (i + 2, j, k);

  int ip1 = Get_Index (i + 1, j, k);

  int im1 = Get_Index (i - 1, j, k);

  int im2 = Get_Index (i - 2, j, k);

  int im3 = Get_Index (i - 3, j, k);

  int im4 = Get_Index (i - 4, j, k);


  return ((-9.0 * (Get_val (ip4) + Get_val (im4))
	   + 128.0 * (Get_val (ip3) + Get_val (im3))
	   - 1008.0 * (Get_val (ip2) + Get_val (im2))
	   + 8064.0 * (Get_val (ip1) + Get_val (im1))
	   - 14350.0 * Get_val ()) * idx2 / 5040.0);

}




double
ffunction::ddyc8 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip4 = Get_Index (i, j + 4, k);

  int ip3 = Get_Index (i, j + 3, k);

  int ip2 = Get_Index (i, j + 2, k);

  int ip1 = Get_Index (i, j + 1, k);

  int im1 = Get_Index (i, j - 1, k);

  int im2 = Get_Index (i, j - 2, k);

  int im3 = Get_Index (i, j - 3, k);

  int im4 = Get_Index (i, j - 4, k);


  return ((-9.0 * (Get_val (ip4) + Get_val (im4))
	   + 128.0 * (Get_val (ip3) + Get_val (im3))
	   - 1008.0 * (Get_val (ip2) + Get_val (im2))
	   + 8064.0 * (Get_val (ip1) + Get_val (im1))
	   - 14350.0 * Get_val ()) * idy2 / 5040.0);

}




double
ffunction::ddzc8 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip4 = Get_Index (i, j, k + 4);

  int ip3 = Get_Index (i, j, k + 3);

  int ip2 = Get_Index (i, j, k + 2);

  int ip1 = Get_Index (i, j, k + 1);

  int im1 = Get_Index (i, j, k - 1);

  int im2 = Get_Index (i, j, k - 2);

  int im3 = Get_Index (i, j, k - 3);

  int im4 = Get_Index (i, j, k - 4);


  return ((-9.0 * (Get_val (ip4) + Get_val (im4))
	   + 128.0 * (Get_val (ip3) + Get_val (im3))
	   - 1008.0 * (Get_val (ip2) + Get_val (im2))
	   + 8064.0 * (Get_val (ip1) + Get_val (im1))
	   - 14350.0 * Get_val ()) * idz2 / 5040.0);

}





double
ffunction::Lapc8 ()
{


  int i,j,k;


  Get_MIndex (i, j, k);


  int ip4 = Get_Index (i + 4, j, k);

  int ip3 = Get_Index (i + 3, j, k);

  int ip2 = Get_Index (i + 2, j, k);

  int ip1 = Get_Index (i + 1, j, k);

  int im1 = Get_Index (i - 1, j, k);

  int im2 = Get_Index (i - 2, j, k);

  int im3 = Get_Index (i - 3, j, k);

  int im4 = Get_Index (i - 4, j, k);




  int jp4 = Get_Index (i, j + 4, k);

  int jp3 = Get_Index (i, j + 3, k);

  int jp2 = Get_Index (i, j + 2, k);

  int jp1 = Get_Index (i, j + 1, k);

  int jm1 = Get_Index (i, j - 1, k);

  int jm2 = Get_Index (i, j - 2, k);

  int jm3 = Get_Index (i, j - 3, k);

  int jm4 = Get_Index (i, j - 4, k);



  int kp4 = Get_Index (i, j, k + 4);

  int kp3 = Get_Index (i, j, k + 3);

  int kp2 = Get_Index (i, j, k + 2);

  int kp1 = Get_Index (i, j, k + 1);

  int km1 = Get_Index (i, j, k - 1);

  int km2 = Get_Index (i, j, k - 2);

  int km3 = Get_Index (i, j, k - 3);

  int km4 = Get_Index (i, j, k - 4);



  return (
	  (-9.0*(Get_val(ip4)+Get_val(im4))+128.0*(Get_val(ip3)+Get_val (im3))-1008.0*(Get_val(ip2)+Get_val(im2))+8064.0*(Get_val(ip1)+Get_val(im1))-14350.0*Get_val()) * idx2 / 5040.0
	  +
	  (-9.0 * (Get_val (jp4) + Get_val (jm4))
	     + 128.0 * (Get_val (jp3) + Get_val (jm3))
	     - 1008.0 * (Get_val (jp2) + Get_val (jm2))
	     + 8064.0 * (Get_val (jp1) + Get_val (jm1))
	     - 14350.0 * Get_val ()) * idy2 / 5040.0
	  + (-9.0 * (Get_val (kp4) + Get_val (km4))
	     + 128.0 * (Get_val (kp3) + Get_val (km3))
	     - 1008.0 * (Get_val (kp2) + Get_val (km2))
	     + 8064.0 * (Get_val (kp1) + Get_val (km1))
	     - 14350.0 * Get_val ()) * idz2 / 5040.0);

}
