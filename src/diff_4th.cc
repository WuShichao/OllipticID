#include "ffunction.h"




double ffunction::dx4()
{


    int i, j, k;

    const int nx = D->Get_Nx();

    const int nm = nx%2==0 ? nx/2 : (nx-1)/2;

    Get_MIndex(i,j,k);


    if (i <= nm )
    {

        int ip1 = Get_Index(i+1,j,k);

        int ip2 = Get_Index(i+2,j,k);

        int ip3 = Get_Index(i+3,j,k);

        int ip4 = Get_Index(i+4,j,k);

        return ( -(  25 * Get_val()
                    -48 * Get_val(ip1)
                    +36 * Get_val(ip2)
                    -16 * Get_val(ip3)
                     +3 * Get_val(ip4) ) * idx / 12.0  );
    }
    else

    {

        int im1 = Get_Index(i-1,j,k);

        int im2 = Get_Index(i-2,j,k);

        int im3 = Get_Index(i-3,j,k);

        int im4 = Get_Index(i-4,j,k);

        return ( ( 25  * Get_val()
                   -48 * Get_val(im1)
                   +36 * Get_val(im2)
                   -16 * Get_val(im3)
                   +3  * Get_val(im4) ) * idx / 12.0  );
    }
     


}







double ffunction::dudx4()
{


    int i, j, k;

    const int nx = D->Get_Nx();

    const int nm = nx%2==0 ? nx/2 : (nx-1)/2;


    Get_MIndex(i,j,k);


    if (i <= nm)

        return ( -25 * idx / 12.0  );

    else


        return ( 25  * idx / 12.0  );


}














double ffunction::dy4()
{


    int i, j, k;

    const int ny = D->Get_Ny();

    const int nm = ny%2==0 ? ny/2 : (ny-1)/2;


    Get_MIndex(i,j,k);


    if (j <= nm)
    {

        int ip1 = Get_Index(i,j+1,k);

        int ip2 = Get_Index(i,j+2,k);

        int ip3 = Get_Index(i,j+3,k);

        int ip4 = Get_Index(i,j+4,k);

        return ( -( 25  * Get_val()
                    -48 * Get_val(ip1)
                    +36 * Get_val(ip2)
                    -16 * Get_val(ip3)
                    +3  * Get_val(ip4) ) * idy / 12.0  );



    }
    else

    {

        int im1 = Get_Index(i,j-1,k);

        int im2 = Get_Index(i,j-2,k);

        int im3 = Get_Index(i,j-3,k);

        int im4 = Get_Index(i,j-4,k);

        return ( ( 25  * Get_val()
                   -48 * Get_val(im1)
                   +36 * Get_val(im2)
                   -16 * Get_val(im3)
                   +3  * Get_val(im4) ) * idy / 12.0  );
    }



}





double ffunction::dudy4()
{


    int i, j, k;

    const int ny = D->Get_Ny();

    const int nm = ny%2==0 ? ny/2 : (ny-1)/2;


    Get_MIndex(i,j,k);


    if (j <= nm)

        return ( -25 * idy / 12.0  );

    else
        
        return ( 25  * idy / 12.0  );


}





double ffunction::dz4()
{


    int i, j, k;

    const int nz = D->Get_Nz();

    const int nm = nz%2==0 ? nz/2 : (nz-1)/2;


    Get_MIndex(i,j,k);


    if (k <= nm)
    {

        int ip1 = Get_Index(i,j,k+1);

        int ip2 = Get_Index(i,j,k+2);

        int ip3 = Get_Index(i,j,k+3);

        int ip4 = Get_Index(i,j,k+4);

        return ( -( 25  * Get_val()
                    -48 * Get_val(ip1)
                    +36 * Get_val(ip2)
                    -16 * Get_val(ip3)
                    +3  * Get_val(ip4) ) * idz / 12.0  );
    }
    else


    {

        int im1 = Get_Index(i,j,k-1);
        
        int im2 = Get_Index(i,j,k-2);

        int im3 = Get_Index(i,j,k-3);

        int im4 = Get_Index(i,j,k-4);

        return ( ( 25  * Get_val()
                   -48 * Get_val(im1)
                   +36 * Get_val(im2)
                   -16 * Get_val(im3)
                   +3  * Get_val(im4) ) * idz / 12.0  );
    }



}







double ffunction::dudz4()
{


    int i, j, k;

    const int nz = D->Get_Nz();

    const int nm = nz%2==0 ? nz/2 : (nz-1)/2;


    Get_MIndex(i,j,k);


    if (k <= nm)

        return ( -25 * idz / 12.0  );

    else

        return ( 25  * idz / 12.0  );



}












// *******************    Segunda derivada    *****************************************








double
ffunction::ddxc4 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip2 = Get_Index (i + 2, j, k);

  int ip1 = Get_Index (i + 1, j, k);

  int im1 = Get_Index (i - 1, j, k);

  int im2 = Get_Index (i - 2, j, k);


  return ((-(Get_val (ip2) + Get_val (im2))
	   + 16.0 * (Get_val (ip1) + Get_val (im1))
	   - 30.0 * Get_val ()) * idx2 / 12.0);

}



double
ffunction::ddyc4 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip2 = Get_Index (i, j + 2, k);

  int ip1 = Get_Index (i, j + 1, k);

  int im1 = Get_Index (i, j - 1, k);

  int im2 = Get_Index (i, j - 2, k);


  return ((-(Get_val (ip2) + Get_val (im2))
	   + 16.0 * (Get_val (ip1) + Get_val (im1))
	   - 30.0 * Get_val ()) * idy2 / 12.0);

}


double
ffunction::ddzc4 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  int ip2 = Get_Index (i, j, k + 2);

  int ip1 = Get_Index (i, j, k + 1);

  int im1 = Get_Index (i, j, k - 1);

  int im2 = Get_Index (i, j, k - 2);


  return ((-(Get_val (ip2) + Get_val (im2))
	   + 16.0 * (Get_val (ip1) + Get_val (im1))
	   - 30.0 * Get_val ()) * idx2 / 12.0);

}





double
ffunction::Lapc4 ()
{

  int i,
    j,
    k;

  Get_MIndex (i, j, k);

  return (((-(Get_val (i + 2, j, k) + Get_val (i - 2, j, k))
	    + 16.0 * (Get_val (i + 1, j, k) + Get_val (i - 1, j, k))
	    - 30.0 * Get_val ()) * idx2
	   + (-(Get_val (i, j + 2, k) + Get_val (i, j - 2, k))
	      + 16.0 * (Get_val (i, j + 1, k) + Get_val (i, j - 1, k))
	      - 30.0 * Get_val ()) * idy2
	   + (-(Get_val (i, j, k + 2) + Get_val (i, j, k - 2))
	      + 16.0 * (Get_val (i, j, k + 1) + Get_val (i, j, k - 1))
	      - 30.0 * Get_val ()) * idz2) / 12.0);

}













/*




double
ffunction::dx4 ()
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

      int ip3 = Get_Index (i + 3, j, k);

      int ip4 = Get_Index (i + 4, j, k);

      return (-(25 * Get_val ()
		- 48 * Get_val (ip1)
		+ 36 * Get_val (ip2)
		- 16 * Get_val (ip3) + 3 * Get_val (ip4)) * idx / 12.0);
    }
  else
   if (i == 1)
    {

      int im1 = Get_Index (i - 1, j, k);

      int ip1 = Get_Index (i + 1, j, k);

      int ip2 = Get_Index (i + 2, j, k);

      int ip3 = Get_Index (i + 3, j, k);


      return ((-3 * Get_val (im1)
	       - 10 * Get_val ()
	       + 18 * Get_val (ip1)
	       - 6 * Get_val (ip2) + Get_val (ip3)) * idx / 12.0);
    }
  else
   if (i == nx - 1)
    {

      int im1 = Get_Index (i - 1, j, k);

      int im2 = Get_Index (i - 2, j, k);

      int im3 = Get_Index (i - 3, j, k);

      int im4 = Get_Index (i - 4, j, k);

      return ((25 * Get_val ()
	       - 48 * Get_val (im1)
	       + 36 * Get_val (im2)
	       - 16 * Get_val (im3) + 3 * Get_val (im4)) * idx / 12.0);
    }
  else
   if (i == nx - 2)
    {

      int ip1 = Get_Index (i + 1, j, k);

      int im1 = Get_Index (i - 1, j, k);

      int im2 = Get_Index (i - 2, j, k);

      int im3 = Get_Index (i - 3, j, k);


      return ((3 * Get_val (ip1)
	       + 10 * Get_val ()
	       - 18 * Get_val (im1)
	       + 6 * Get_val (im2) - Get_val (im3)) * idx / 12.0);

    }

  else
    {
      int ip2 = Get_Index (i + 2, j, k);

      int ip1 = Get_Index (i + 1, j, k);

      int im1 = Get_Index (i - 1, j, k);

      int im2 = Get_Index (i - 2, j, k);



      return ((Get_val (im2) - Get_val (ip2)
	       + 8.0 * (Get_val (ip1) - Get_val (im1))) * idx / 12.0);
    }



}







double
ffunction::dudx4 ()
{


  int i,
    j,
    k;

  const int nx = D->Get_Nx ();



  Get_MIndex (i, j, k);


  if (i == 0)

    return (-25 * idx / 12.0);

  else
   if (i == 1)

    return (-10 * idx / 12.0);


  else
   if (i == nx - 1)

    return (25 * idx / 12.0);

  else
   if (i == nx - 2)


    return (10 * idx / 12.0);



  else

    return (0);


}














double
ffunction::dy4 ()
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

      int ipppp = Get_Index (i, j + 4, k);

      return (-(25 * Get_val ()
		- 48 * Get_val (ip)
		+ 36 * Get_val (ipp)
		- 16 * Get_val (ippp) + 3 * Get_val (ipppp)) * idy / 12.0);
    }
  else
   if (j == 1)
    {

      int im = Get_Index (i, j - 1, k);

      int ip = Get_Index (i, j + 1, k);

      int ipp = Get_Index (i, j + 2, k);

      int ippp = Get_Index (i, j + 3, k);


      return ((-3 * Get_val (im)
	       - 10 * Get_val ()
	       + 18 * Get_val (ip)
	       - 6 * Get_val (ipp) + Get_val (ippp)) * idy / 12.0);
    }
  else
   if (j == ny - 1)
    {

      int im = Get_Index (i, j - 1, k);

      int imm = Get_Index (i, j - 2, k);

      int immm = Get_Index (i, j - 3, k);

      int immmm = Get_Index (i, j - 4, k);

      return ((25 * Get_val ()
	       - 48 * Get_val (im)
	       + 36 * Get_val (imm)
	       - 16 * Get_val (immm) + 3 * Get_val (immmm)) * idy / 12.0);
    }
  else
   if (j == ny - 2)
    {

      int ip = Get_Index (i, j + 1, k);

      int im = Get_Index (i, j - 1, k);

      int imm = Get_Index (i, j - 2, k);

      int immm = Get_Index (i, j - 3, k);


      return ((3 * Get_val (ip)
	       + 10 * Get_val ()
	       - 18 * Get_val (im)
	       + 6 * Get_val (imm) - Get_val (immm)) * idy / 12.0);

    }

  else
    {
      int ipp = Get_Index (i, j + 2, k);

      int ip = Get_Index (i, j + 1, k);

      int im = Get_Index (i, j - 1, k);

      int imm = Get_Index (i, j - 2, k);



      return ((Get_val (imm)
	       - 8 * Get_val (im)
	       + 8 * Get_val (ip) - Get_val (ipp)) * idy / 12.0);
    }



}





double
ffunction::dudy4 ()
{


  int i,
    j,
    k;

  const int ny = D->Get_Ny ();



  Get_MIndex (i, j, k);


  if (j == 0)

    return (-25 * idy / 12.0);

  else
   if (j == 1)

    return (-10 * idy / 12.0);


  else
   if (j == ny - 1)

    return (25 * idy / 12.0);

  else
   if (j == ny - 2)


    return (10 * idy / 12.0);

  else

    return (0);


}





double
ffunction::dz4 ()
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

      int ipppp = Get_Index (i, j, k + 4);

      return (-(25 * Get_val ()
		- 48 * Get_val (ip)
		+ 36 * Get_val (ipp)
		- 16 * Get_val (ippp) + 3 * Get_val (ipppp)) * idz / 12.0);
    }
  else
   if (k == 1)
    {

      int im = Get_Index (i, j, k - 1);

      int ip = Get_Index (i, j, k + 1);

      int ipp = Get_Index (i, j, k + 2);

      int ippp = Get_Index (i, j, k + 3);


      return ((-3 * Get_val (im)
	       - 10 * Get_val ()
	       + 18 * Get_val (ip)
	       - 6 * Get_val (ipp) + Get_val (ippp)) * idz / 12.0);
    }
  else
   if (k == nz - 1)
    {

      int im = Get_Index (i, j, k - 1);

      int imm = Get_Index (i, j, k - 2);

      int immm = Get_Index (i, j, k - 3);

      int immmm = Get_Index (i, j, k - 4);

      return ((25 * Get_val ()
	       - 48 * Get_val (im)
	       + 36 * Get_val (imm)
	       - 16 * Get_val (immm) + 3 * Get_val (immmm)) * idz / 12.0);
    }
  else
   if (k == nz - 2)
    {

      int ip = Get_Index (i, j, k + 1);

      int im = Get_Index (i, j, k - 1);

      int imm = Get_Index (i, j, k - 2);

      int immm = Get_Index (i, j, k - 3);


      return ((3 * Get_val (ip)
	       + 10 * Get_val ()
	       - 18 * Get_val (im)
	       + 6 * Get_val (imm) - Get_val (immm)) * idz / 12.0);

    }

  else
    {

      int ipp = Get_Index (i, j, k + 2);

      int ip = Get_Index (i, j, k + 1);

      int im = Get_Index (i, j, k - 1);

      int imm = Get_Index (i, j, k - 2);



      return ((Get_val (imm)
	       - 8 * Get_val (im)
	       + 8 * Get_val (ip) - Get_val (ipp)) * idz / 12.0);
    }



}







double
ffunction::dudz4 ()
{


  int i,
    j,
    k;

  const int nz = D->Get_Nz ();



  Get_MIndex (i, j, k);


  if (k == 0)

    return (-25 * idz / 12.0);

  else
   if (k == 1)

    return (-10 * idz / 12.0);


  else
   if (k == nz - 1)

    return (25 * idz / 12.0);

  else
   if (k == nz - 2)


    return (10 * idz / 12.0);

  else

    return (0);


}


*/

double
ffunction::ddxyc4 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);


  return ((Get_val (i - 2, j - 2, k) -
	   8.0 * Get_val (i - 2, j - 1, k) +
	   8.0 * Get_val (i - 2, j + 1, k) -
	   Get_val (i - 2, j + 2, k) -
	   8.0 * Get_val (i - 1, j - 2, k) +
	   64.0 * Get_val (i - 1, j - 1, k) -
	   64.0 * Get_val (i - 1, j + 1, k) +
	   8.0 * Get_val (i - 1, j + 2, k) +
	   8.0 * Get_val (i + 1, j - 2, k) -
	   64.0 * Get_val (i + 1, j - 1, k) +
	   64.0 * Get_val (i + 1, j + 1, k) -
	   8.0 * Get_val (i + 1, j + 2, k) -
	   Get_val (i + 2, j - 2, k) +
	   8.0 * Get_val (i + 2, j - 1, k) -
	   8.0 * Get_val (i + 2, j + 1, k) +
	   Get_val (i + 2, j + 2, k)) * idx * idy / 144.0);

}


double
ffunction::ddxzc4 ()
{


  int i,
    j,
    k;


  Get_MIndex (i, j, k);



  return ((Get_val (i - 2, j, k - 2) -
	   8.0 * Get_val (i - 2, j, k - 1) +
	   8.0 * Get_val (i - 2, j, k + 1) -
	   Get_val (i - 2, j, k + 2) -
	   8.0 * Get_val (i - 1, j, k - 2) +
	   64.0 * Get_val (i - 1, j, k - 1) -
	   64.0 * Get_val (i - 1, j, k + 1) +
	   8.0 * Get_val (i - 1, j, k + 2) +
	   8.0 * Get_val (i + 1, j, k - 2) -
	   64.0 * Get_val (i + 1, j, k - 1) +
	   64.0 * Get_val (i + 1, j, k + 1) -
	   8.0 * Get_val (i + 1, j, k + 2) -
	   Get_val (i + 2, j, k - 2) +
	   8.0 * Get_val (i + 2, j, k - 1) -
	   8.0 * Get_val (i + 2, j, k + 1) +
	   Get_val (i + 2, j, k + 2)) * idx * idz / 144.0);


}


double
ffunction::ddyzc4 ()
{



  int i,
    j,
    k;


  Get_MIndex (i, j, k);



  return ((Get_val (i, j - 2, k - 2) -
	   8.0 * Get_val (i, j - 2, k - 1) +
	   8.0 * Get_val (i, j - 2, k + 1) -
	   Get_val (i, j - 2, k + 2) -
	   8.0 * Get_val (i, j - 1, k - 2) +
	   64.0 * Get_val (i, j - 1, k - 1) -
	   64.0 * Get_val (i, j - 1, k + 1) +
	   8.0 * Get_val (i, j - 1, k + 2) +
	   8.0 * Get_val (i, j + 1, k - 2) -
	   64.0 * Get_val (i, j + 1, k - 1) +
	   64.0 * Get_val (i, j + 1, k + 1) -
	   8.0 * Get_val (i, j + 1, k + 2) -
	   Get_val (i, j + 2, k - 2) +
	   8.0 * Get_val (i, j + 2, k - 1) -
	   8.0 * Get_val (i, j + 2, k + 1) +
	   Get_val (i, j + 2, k + 2)) * idy * idz / 144.0);


}
