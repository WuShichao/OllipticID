



//============================ Olliptic elliptic.cc ========================//
//  
//     
//
//
//
//
//
//============================== Pablo Galaviz 2009 ========================//





#include "elliptic.h"










void
elliptic::Make (interface * Oll)
{



  //== Seting parameters (integers) ==//

  levels = Oll->Get_integer_par ("levels");

  order = Oll->Get_integer_par ("order");


  nu1 = Oll->Get_integer_par ("nu1");

  nu2 = Oll->Get_integer_par ("nu2");


  ps_case = Oll->Get_integer_par("ps_case");

 

  if(!(ps_case == 1 || ps_case == 2 || ps_case == 3 || ps_case == 4))
    error_exit("Ollitic can solve 4 cases for scalar field+puncture","set -psc to 1, 2, 3 or 4");


  //== Seting parameters (doubles) ==//


  eps_solv = Oll->Get_double_par ("tolerance");


  eps_coarse = 2 * numeric_limits < double >::epsilon ();

  if (eps_solv < eps_coarse)

    eps_solv = eps_coarse;




  ps_r0  = Oll->Get_double_par ("ps_r0");
  ps_isigma  = 1/Oll->Get_double_par ("ps_sigma");

  double a2 = Oll->Get_double_par ("ps_a2");
  ps_phi0 = (a2*a2/(1+a2*a2))*Oll->Get_double_par ("ps_phi0");

  ps_d1 = M_PI*ps_phi0*ps_phi0*ps_isigma*ps_isigma; 

  ps_e1 = 1/(16*Oll->Get_double_par ("ps_a2"));
  ps_e2 = 4*sqrt(M_PI/3);
  ps_de1 = 4*M_PI*ps_phi0*ps_phi0*ps_isigma*ps_isigma;




  //== Seting parameters (logics) ==//


  print_1dt = Oll->Get_logic_par ("print_1d_it");

  print_1d = Oll->Get_logic_par ("print_1d");

  print_2d = Oll->Get_logic_par ("print_2d");

  print_3d = Oll->Get_logic_par ("print_3d");

  print_interpol = Oll->Get_logic_par ("print_interpol");

  interpol_output = Oll->Get_logic_par ("interpol_output");

  fit_approx = Oll->Get_logic_par ("fit_approx");



  dir_name = Oll->Get_string_par ("dir_name");


  tov_filename = Oll->Get_string_par("file_TOV");


  dom = Oll->Get_Domain_Sym ();
  equation = Oll->Get_Equation ();
  method = Oll->Get_Method ();
  boundary = Oll->Get_Boundary ();
  output = Oll->Get_TOutput ();




  bhmp = Oll->Get_bhmp ();
  bhqp = Oll->Get_bhqp ();
 
  ADM_mass.resize(bhmp.size());
  

  bhx = Oll->Get_bhx ();
  bhy = Oll->Get_bhy ();
  bhz = Oll->Get_bhz ();

  bhpx = Oll->Get_bhpx ();
  bhpy = Oll->Get_bhpy ();
  bhpz = Oll->Get_bhpz ();


  bhsx = Oll->Get_bhsx ();
  bhsy = Oll->Get_bhsy ();
  bhsz = Oll->Get_bhsz ();


  NumBH = bhmp.size ();


  u_at_p.resize(NumBH);



  if (!(bhx.size () == bhy.size () && bhx.size () == bhz.size ()))
    {
      cout << "Error: making elliptic object, wrong positions." << endl;

      exit (1);

    }


  boxes = bhx.size ();

  num_dom = levels * boxes;


  PI = 4.0*atan(1.0);

  it = 0;

  git = 0;


  total_nodes = 1;

  my_rank = 0;

#ifdef OLLIN_MPI

  total_nodes = MPI::COMM_WORLD.Get_size ();

  my_rank = MPI::COMM_WORLD.Get_rank ();

#endif


  Set_Method ();


  Set_Grids (Oll->Get_integer_par ("nx"),
	     Oll->Get_integer_par ("ny"),
	     Oll->Get_integer_par ("nz"),
	     Oll->Get_double_par ("dx"),
	     Oll->Get_double_par ("dy"), Oll->Get_double_par ("dz"));


  Set_Equation ();


  Set_Coefficients ();


  if (Oll->Get_logic_par ("print_coeff"))

    Print_Coeff ();




  stringstream file_output;

  file_output << dir_name << "stdout" << my_rank << ".out";

  stdout_file_name = file_output.str ();


  stringstream file_norm;

  file_norm << dir_name << "Loo.gp";

  norm_file_name = file_norm.str ();



  

}







void
elliptic::Set_Method ()
{




  switch (method)
    {

    case GAUSS_SEIDEL:


      solve_ell = &elliptic::Gauss_Seidel;



      break;


    case FMG:


      solve_ell = &elliptic::FullMG;


      break;


    case VCYCLE:


      solve_ell = &elliptic::Vcycle;


      break;


    case WCYCLE:


      solve_ell = &elliptic::Wcycle;



      break;


    case APP:


      solve_ell = &elliptic::ApproximateID;


      break;



    }






}




void
elliptic::Set_Grids (size_t Nx, size_t Ny, size_t Nz, double Dx, double Dy,
		     double Dz)
{


  vector < size_t > nx (1, Nx);

  vector < size_t > ny (1, Ny);

  vector < size_t > nz (1, Nz);


  vector < double >dx (1, Dx);

  vector < double >dy (1, Dy);

  vector < double >dz (1, Dz);


  D.clear ();




  size_t subgridord = 2;


  t_symmetry sym_x;

  t_symmetry sym_y;

  t_symmetry sym_z;




  vector < size_t > num_dom_in_level;


  switch (dom)
    {


    case BITANT_X:


      sym_x = REFLECT;

      sym_y = NONE;

      sym_z = NONE;



      break;


    case BITANT_Y:


      sym_x = NONE;

      sym_y = REFLECT;

      sym_z = NONE;



      break;


    case BITANT_Z:


      sym_x = NONE;

      sym_y = NONE;

      sym_z = REFLECT;



      break;


    case QUADRANT_XY:


      sym_x = REFLECT;

      sym_y = REFLECT;

      sym_z = NONE;



      break;


    case QUADRANT_XZ:


      sym_x = REFLECT;

      sym_y = NONE;

      sym_z = REFLECT;



      break;


    case QUADRANT_YZ:


      sym_x = NONE;

      sym_y = REFLECT;

      sym_z = REFLECT;



      break;



    case OCTANT:


      sym_x = REFLECT;

      sym_y = REFLECT;

      sym_z = REFLECT;



      break;




    case ROTATE_X:


      sym_x = ROTANT;

      sym_y = NONE;

      sym_z = NONE;



      break;



    case ROTATE_Y:


      sym_x = NONE;

      sym_y = ROTANT;

      sym_z = NONE;



      break;


    case ROTATE_Z:


      sym_x = NONE;

      sym_y = NONE;

      sym_z = ROTANT;



      break;



    case FULL:

    default:

      sym_x = NONE;

      sym_y = NONE;

      sym_z = NONE;


    }


  size_t base_x = sym_x == NONE ? 12 : 16;

  size_t base_y = sym_y == NONE ? 12 : 16;

  size_t base_z = sym_z == NONE ? 12 : 16;



  if (method != GAUSS_SEIDEL && method != APP)
    {
      bool new_level = true;

      bool limit_x = false;

      bool limit_y = false;

      bool limit_z = false;


      size_t l = 0;

      size_t ntmp;

      do
	{

	  //sub-levels for multigrid direction X

	  ntmp = (nx[l] + 1) / 2;


	  if (ntmp < base_x)
	    {

	      limit_x = true;

	      ntmp = nx[l];

	    }

	  nx.push_back (ntmp);




	  if (nx[l + 1] == nx[l])

	    dx.push_back (dx[l]);

	  else

	    dx.push_back (dx[l] * 2.0);





	  //sub-levels for multigrid direction Y

	  ntmp = (ny[l] + 1) / 2;


	  if (ntmp < base_y)
	    {

	      limit_y = true;

	      ntmp = ny[l];

	    }

	  ny.push_back (ntmp);

	  if (ny[l + 1] == ny[l])

	    dy.push_back (dy[l]);

	  else

	    dy.push_back (dy[l] * 2.0);



	  //sub-levels for multigrid direction Z
	  ntmp = (nz[l] + 1) / 2;



	  if (ntmp < base_z)
	    {

	      limit_z = true;

	      ntmp = nz[l];

	    }

	  nz.push_back (ntmp);

	  if (nz[l + 1] == nz[l])

	    dz.push_back (dz[l]);

	  else

	    dz.push_back (dz[l] * 2.0);


	  if (limit_x && limit_y && limit_z)
	    {
	      new_level = false;

	      nx.pop_back ();

	      ny.pop_back ();

	      nz.pop_back ();


	      dx.pop_back ();

	      dy.pop_back ();

	      dz.pop_back ();

	    }


	  l++;

	}
      while (new_level);


      if (nx.size () <= 1)
	{

	  if (my_rank == 0)

	    cout << "Error: There are not enough grid-points for multigrid" <<
	      endl;

	  exit (1);

	}
      else

	mgrid_levels = nx.size () - 1;

    }

  else
    {

      mgrid_levels = 0;


      D.push_back (domain
		   (Nx, Ny, Nz, Dx, Dy, Dz, 0, 0, 0, order, sym_x, sym_y,
		    sym_z));

      num_dom_in_level.push_back (1);


    }











  if (method != GAUSS_SEIDEL && method != APP)
    {


      levels += mgrid_levels;



      for (int l = mgrid_levels; l >= 0; l--)
	{


	  size_t gz = l == 0 ? order : subgridord;

	  D.push_back (domain
		       (nx[l], ny[l], nz[l], dx[l], dy[l], dz[l], 0, 0, 0, gz,
			sym_x, sym_y, sym_z));

	  num_dom_in_level.push_back (1);




	}

    }








  for (int l = mgrid_levels + 1; (size_t) l < levels; l++)
    {

      for (size_t i = 0; i < boxes; i++)

	{


	  double pot = l - mgrid_levels;

	  double _dx = Dx / pow (2.0, pot);


	  double _dy = Dy / pow (2.0, pot);


	  double _dz = Dz / pow (2.0, pot);


	  double xi,
	    yi,
	    zi;

	  double xf,
	    yf,
	    zf;

	  double xc,
	    yc,
	    zc;


	  D[mgrid_levels].Find_cell_xyz (bhx[i], bhy[i], bhz[i], xi, yi, zi,
					 xf, yf, zf);




	  if (dcomp (xi, bhx[i]) || dcomp (xf, bhx[i]))

	    xc = bhx[i];

	  else
	    {

	      while (xi + _dx < bhx[i])

		xi += _dx;

	      while (xf - _dx > bhx[i])

		xf -= _dx;


	      xc = fabs (xi - bhx[i]) < fabs (xf - bhx[i]) ? xi : xf;

	    }



	  if (dcomp (yi, bhy[i]) || dcomp (yf, bhy[i]))

	    yc = bhy[i];

	  else
	    {

	      while (yi + _dy < bhy[i])

		yi += _dy;

	      while (yf - _dy > bhy[i])

		yf -= _dy;


	      yc = fabs (yi - bhy[i]) < fabs (yf - bhy[i]) ? yi : yf;

	    }



	  if (dcomp (zi, bhz[i]) || dcomp (zf, bhz[i]))

	    zc = bhz[i];

	  else
	    {

	      while (zi + _dz < bhz[i])

		zi += _dz;

	      while (zf - _dz > bhz[i])

		zf -= _dz;


	      zc = fabs (zi - bhz[i]) < fabs (zf - bhz[i]) ? zi : zf;

	    }




	  D.push_back (domain (Nx, Ny, Nz, _dx, _dy, _dz,
			       xc, yc, zc, order, sym_x, sym_y, sym_z));




	}

      num_dom_in_level.push_back (boxes);


    }




  size_t l,
    i,
    j;

  bool restore;

  for (l = mgrid_levels + 1; l < levels; l++)
    {

      restore = false;

      for (i = 0; i < num_dom_in_level[l]; i++)

	for (j = 0; j < num_dom_in_level[l]; j++)

	  if (i != j)
	    {

	      size_t ibox = 0;

	      for (size_t ll = 0; ll < l; ll++)

		ibox += num_dom_in_level[ll];

	      ibox += i;


	      size_t jbox = 0;

	      for (size_t ll = 0; ll < l; ll++)

		jbox += num_dom_in_level[ll];

	      jbox += j;




	      if (D[ibox] && D[jbox])
		{


		  D[ibox] = D[ibox] + D[jbox];

		  D.erase (D.begin () + jbox);

		  num_dom_in_level[l]--;


		  restore = true;





		}


	    }
      if (restore)

	l--;

    }






  size_t ibox = 0;

  for (size_t l = 0; l < levels; l++)

    for (size_t i = 0; i < num_dom_in_level[l]; i++)
      {

	D[ibox].Set_level (l);

	ibox++;

      }






}





















void
elliptic::App_Ell_Operator (size_t ll, size_t ibox)
{



  for (size_t iv = 0; iv < var_names.size (); iv++)
    {



      stringstream stage_var;

      stage_var << iv;

      ffunction *R = vars.Get_ffunction (ll, ibox, "r" + stage_var.str ());

      R->Reset_Norms ();

      vars.Set_Iterator_ibox (total_names, INSIDE, ll, ibox);


      do
	{



	  *R = OPERATOR (iv, ll, ibox);


	  R->Add_To_Norms ();


	}
      while (vars.End (total_names));




      App_Symmetry (R);

      R->Sync ();



    }


}









/*
  double
  elliptic::Smooth (size_t ll, size_t ibox, size_t nu)
  {



  double gLinf = 0;


  for (size_t iv = 0; iv < var_names.size (); iv++)
  {





  ffunction *U = vars.Get_ffunction (ll, ibox, var_names[iv]);


  double d_x = U->Get_domain ()->Get_dx ();

  double d_y = U->Get_domain ()->Get_dy ();

  double d_z = U->Get_domain ()->Get_dz ();

  double d_r2 = 1.5 * (d_x * d_x + d_y * d_y + d_z * d_z);

  stringstream stage_var;

  stage_var << iv;

  ffunction *R = vars.Get_ffunction (ll, ibox, "r" + stage_var.str ());

  ffunction *S = vars.Get_ffunction (ll, ibox, "rho" + stage_var.str ());


  for (size_t kk = 0; kk < nu; kk++)
  {






  vars.Set_Iterator_ibox (total_names, INSIDE, ll, ibox);


  R->Reset_Norms ();





  do
  {


  double x = U->Get_X ();

  double y = U->Get_X ();

  double z = U->Get_X ();

  *R = OPERATOR (iv, ll, ibox) - S->Get_val ();


  if (x * x + y * y + z * z >= d_r2)
  {


  double dr = duOPERATOR (iv, ll, ibox);


  *U -= R->Get_val () / dr;

  }
  else
  {

  *U = (this->*Boundary[iv]) (x, y, z);

  }


  R->Add_To_Norms ();


  }
  while (vars.End ());



  it++;


  U->Sync ();

  App_Symmetry (U);


  R->Sync_Norms ();


  double Linf = R->Get_Norm_LInf ();

  gLinf = Linf > gLinf ? Linf : gLinf;

  }


  }



  return (gLinf);


  }

*/




double elliptic::Smooth( size_t ll, size_t ibox, size_t nu){



  double gLinf = 0;

    
  for(size_t iv = 0; iv < var_names.size(); iv++)
    {
   
  



      ffunction *U = vars.Get_ffunction(ll,ibox,var_names[iv]);

        

      stringstream stage_var;

      stage_var << iv;

      ffunction *R = vars.Get_ffunction(ll,ibox,"r"+stage_var.str());

      ffunction *S = vars.Get_ffunction(ll,ibox,"rho"+stage_var.str());


      for(size_t kk=0; kk < nu; kk++)
        {
      


            

 
	  vars.Set_Iterator_ibox(total_names,INSIDE,ll,ibox);


	  R->Reset_Norms();


        

    
	  do
	    {



	      double tmp1 = OPERATOR(iv,ll,ibox);

	      double tmp2 = S->Get_val(); 
          
            
	      *R = OPERATOR(iv,ll,ibox) - S->Get_val();

                
	      double dr = duOPERATOR(iv,ll,ibox);


      
	      *U -= R->Get_val() / dr ;


                        
	      R->Add_To_Norms();

			      
	      if(!( U->Get_val() == U->Get_val())){

		cout << "NAN check error.. " <<  R->Get_val() 
		     << " iv = " << iv 
		     << " ll = " << ll
		     << " ibox = " << ibox 
		     << " Operator: " << OPERATOR(iv,ll,ibox)
		     << " Source: " << S->Get_val()
		     << " dr: " << dr << endl
		     <<endl;

		exit(0);

	      }
	      


      
	    }while ( vars.End() );


    
	  it++;

    
	  U->Sync();

	  App_Symmetry(U);

        
	  R->Sync_Norms();
                
        
	  double Linf = R->Get_Norm_LInf();


	  gLinf = Linf > gLinf ? Linf : gLinf;

        }

        
    }


    
  return(gLinf);


}








void elliptic::Vcycle ()
{

  double res = 1.0;

  size_t vcycle = 0;

  int signo_prev = 1.0;

  int signo = 1.0;


  double ll =0.0;

  double res_prev = 1.0;
  double res_prev_prev;

  res = 1e12;

  bool reset = false;

  double factor = 1.1;

  if(fit_approx)

    Get_Puncture_App();

  size_t lev = levels-1;

  do
    {

      res_prev_prev = res_prev;

      res_prev = res;

      res = FASN (levels - 1, levels - 1, false, vcycle);

      vcycle++;

      reset = false;
      
      
      if(equation==NSO_ID && vcycle%2==0)
	//if(equation==NSO_ID )
	{

	  cout << "Reset rhoADM ..." << endl;

	  vector<string> var_n;

	  var_n.push_back("psi");

	  var_n.push_back("matter_rho");

	  var_n.push_back("matter_epsl");

	  var_n.push_back("rhoADM");

	  vars.Set_Iterator(var_n,ALL);

	  double twopi=8.0*atan(1.0);

	  do
	    {
	    
	      vars.Set_val (twopi*pow(vars.Get_val("psi"),8)*vars.Get_val("matter_rho")*(1.0+vars.Get_val("matter_epsl")), "rhoADM");
	    

	      
	    }
	  while (vars.End (var_n));


	}

      

    }
  while (res > eps_solv && res < res_prev);
  
  
  for (size_t j = 0; bhmp[j] != 0.0 && j < NumBH; j++)

    cout << "mp [" << j << "] = " << bhmp[j] << endl;

  for (size_t j = 0; bhqp[j] != 0.0 && j < NumBH; j++)

    cout << "qp [" << j << "] = " << bhqp[j] << endl;


}





void
elliptic::Wcycle ()
{

  double res = 1;

  size_t vcycle = 0;



  double res_prev;


  res = 1e12;


  do
    {

      res_prev = res;

      res = FASN (levels - 1, levels - 1, true, vcycle);

      vcycle++;


    }
  while (res > eps_solv && res < res_prev);




}







void
elliptic::Gauss_Seidel ()
{


  double norm;

  double gnorm;




  for (size_t l = 0; l < levels; l++)
    {


      if (l > 0)

	for (size_t iv = 0; iv < var_names.size (); iv++)

	  vars.Transferir (l - 1, var_names[iv], l, var_names[iv], ALL,
			   LAGRANGE, order);


      do
	{

	  gnorm = 0;




	  for (size_t i = 0; i < vars.Get_num_dom (l); i++)
	    {


	      
	      if (l == mgrid_levels)

		App_Bound (l, i);
	      
	      
	      norm = Smooth (l, i);



	      gnorm = norm > gnorm ? norm : gnorm;


	    }


	  it++;

	  Print_Cycle (l, levels - 1, gnorm, it);



	}
      while (norm > eps_solv);


    }



}













double
elliptic::FASN (const size_t k, const size_t k_top, bool Wcycle,
		size_t icycle)
{




  m_interpol method_prolong = LINEAR;

  //    size_t ord_prolong = order-1;    


  //    m_interpol method_restrict = AVERAGE;

  //    size_t ord_restrict = order+1;    


  //    m_interpol method_transfere = HERMITE;

  size_t ord_transfere = order + 1;



  double norm_k;


  double norm_km1;



  size_t cycle = icycle;




  if (k == 0)

    exit (0);


  do
    {







      norm_k = 0;

      for (size_t i = 0; i < vars.Get_num_dom (k); i++)
	{

	  double norm;

	  if (k <= mgrid_levels)

	    App_Bound (k, i);

	  norm = Smooth (k, i, nu1);

	  norm_k = norm > norm_k ? norm : norm_k;

	}



      Print_Pre_Cycle (k, k_top, norm_k, cycle);



      //I^0_1 u^1

      for (size_t iv = 0; iv < var_names.size (); iv++)

	vars.Transferir (k, var_names[iv], k - 1, var_names[iv], INSIDE,
			 LAGRANGE, 3, false);




      /*

      for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	for (size_t iv = 0; iv < source_names.size (); iv++)

	  App_Symmetry (vars.Get_ffunction (k - 1, j, source_names[iv]));



      if (k <= mgrid_levels)

	for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	  App_Bound (k - 1, j);



      for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	for (size_t iv = 0; iv < var_names.size (); iv++)

	  App_Symmetry (vars.Get_ffunction (k - 1, j, var_names[iv]));

      



	  */


      for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)
	{


	  App_Ell_Operator (k - 1, j);

	  if (k <= mgrid_levels)

	    App_Boundary_Operator (k - 1, j);



	}

      






      for (size_t i = 0; i < vars.Get_num_dom (k); i++)
	{

	  App_Ell_Operator (k, i);

	  if (k <= mgrid_levels)

	    App_Boundary_Operator (k, i);



	}










      for (size_t i = 0; i < vars.Get_num_dom (k); i++)

	for (size_t iv = 0; iv < var_names.size (); iv++)
	  {

	    stringstream stage_var;

	    stage_var << iv;

	    ffunction *R = vars.Get_ffunction (k, i, "r" + stage_var.str ());

	    ffunction *S = vars.Get_ffunction (k, i, "rho" + stage_var.str ());

	    ffunction *W = vars.Get_ffunction (k, i, "w" + stage_var.str ());

	    App_Symmetry (W);


	    *W = (*S - *R);


	    App_Symmetry (W);


	  }










      for (size_t iv = 0; iv < var_names.size (); iv++)
	{

	  stringstream stage_var;

	  stage_var << iv;


	  if (k - 1 <= mgrid_levels)

	    vars.Transferir (k, "w" + stage_var.str (), k - 1,
			     "w" + stage_var.str (), ALL, LAGRANGE, 2);

	  else

	    vars.Transferir (k, "w" + stage_var.str (), k - 1,
			     "w" + stage_var.str (), ALL, LAGRANGE, 2, false);



	  for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	    App_Symmetry (vars.
			  Get_ffunction (k - 1, j, "w" + stage_var.str ()));

	}






      for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	for (size_t iv = 0; iv < var_names.size (); iv++)
	  {

	    stringstream stage_var;

	    stage_var << iv;

	    ffunction *R =
	      vars.Get_ffunction (k - 1, j, "r" + stage_var.str ());

	    ffunction *S =
	      vars.Get_ffunction (k - 1, j, "rho" + stage_var.str ());

	    ffunction *W =
	      vars.Get_ffunction (k - 1, j, "w" + stage_var.str ());






	    if (k <= mgrid_levels)
	      {

		*S = (*W + *R);


	      }

	    else
	      {

		*W += *R;


		for (size_t i = 0; i < vars.Get_num_dom (k); i++)
		  {

		    ffunction *Wk =
		      vars.Get_ffunction (k, i, "w" + stage_var.str ());

		    S->Set_In (*W, Wk->Get_domain (), false);





		  }
	      }

	    App_Symmetry (S);


	  }








      for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	for (size_t iv = 0; iv < var_names.size (); iv++)
	  {

	    stringstream stage_var;

	    stage_var << iv;

	    ffunction *U = vars.Get_ffunction (k - 1, j, var_names[iv]);

	    ffunction *ukm1 =
	      vars.Get_ffunction (k - 1, j, "ww" + stage_var.str ());




	    *ukm1 = *U;

	  }






      if (k == 1)
	{



	  for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)
	    {

	      double norm_prev;

	      norm_km1 = 1;

	      do
		{

		  norm_prev = norm_km1;

		  App_Bound (k - 1, j);

		  norm_km1 = Smooth (k - 1, j);


		}
	      while (norm_km1 > eps_coarse && norm_prev != norm_km1);

	    }



	  Print_Pre_Cycle (k - 1, k_top, norm_km1, cycle);




	}
      else

	norm_km1 = FASN (k - 1, k_top, Wcycle, cycle);







      for (size_t j = 0; j < vars.Get_num_dom (k - 1); j++)

	for (size_t iv = 0; iv < var_names.size (); iv++)
	  {

	    stringstream stage_var;

	    stage_var << iv;

	    ffunction *U = vars.Get_ffunction (k - 1, j, var_names[iv]);

	    ffunction *ukm1 =
	      vars.Get_ffunction (k - 1, j, "ww" + stage_var.str ());

	    ffunction *V =
	      vars.Get_ffunction (k - 1, j, "v" + stage_var.str ());

	    *V = (*U - *ukm1);


	    App_Symmetry (V);


	  }






      for (size_t iv = 0; iv < var_names.size (); iv++)
	{

	  stringstream stage_var;

	  stage_var << iv;



	  if (k <= mgrid_levels)

	    vars.Transferir (k - 1, "v" + stage_var.str (), k,
			     "v" + stage_var.str (), ALL, method_prolong, 3,
			     true, true);

	  else

	    vars.Transferir (k - 1, "v" + stage_var.str (), k,
			     "v" + stage_var.str (), INSIDE, method_prolong,
			     ord_transfere);


	  for (size_t i = 0; i < vars.Get_num_dom (k); i++)

	    App_Symmetry (vars.Get_ffunction (k, i, "v" + stage_var.str ()));



	}




      for (size_t iv = 0; iv < var_names.size (); iv++)
	{

	  stringstream stage_var;

	  stage_var << iv;

	  for (size_t i = 0; i < vars.Get_num_dom (k); i++)


	    *vars.Get_ffunction (k, i, var_names[iv]) +=
	      *vars.Get_ffunction (k, i, "v" + stage_var.str ());



	  if (k > mgrid_levels)

	    vars.Transferir (k - 1, var_names[iv], k, var_names[iv], BOUNDARY,
			     LAGRANGE, order + 1);

	}






      norm_k = 0;

      for (size_t i = 0; i < vars.Get_num_dom (k); i++)
	{

	  double norm;

	  if (k <= mgrid_levels)

	    App_Bound (k, i);

	  norm = Smooth (k, i, nu2);

	  norm_k = norm > norm_k ? norm : norm_k;

	}




      Print_Post_Cycle (k, k_top, norm_k, cycle);

      if(k==k_top && (equation==PUNCTURES || equation==PUNCTURES_SCALAR))

	Compute_ADM_Mass();


      if(k==k_top && equation==BOWEN_YORK )

	Compute_Solution();



      cycle++;

      it++;


    }
  while (norm_k > eps_solv && Wcycle);



  return (norm_k);







}








void elliptic::ApproximateID ()
{


  for (size_t lev = 0; lev < levels; lev++)

    for (size_t ibox = 0; ibox < vars.Get_num_dom (lev); ibox++)
      {


	double r,
	  xi,
	  yi,
	  zi;


	double J2,
	  Jx,
	  Jy,
	  Jz;

	double P2,
	  Px,
	  Py,
	  Pz;


	double sbl;

	double pbl;


	double uj,
	  u0j,
	  u2j;

	double up,
	  u0p,
	  u2p;

	double l1,
	  R,
	  muj,
	  mup,
	  PL2j,
	  PL2p;

	double Rx,
	  Ry,
	  Rz;

	double l2,
	  l3,
	  l4,
	  l5;

	double PXJ[3],
	  uc;


	double imass,
	  imass2;



	ffunction *U = vars.Get_ffunction (lev, ibox, "u");

	double dr_at_p = 0.5*sqrt(U->Get_dx()*U->Get_dx()+U->Get_dy()*U->Get_dy()+U->Get_dz()*U->Get_dz());


	U->Set_Iterator (ALL);

	do
	  {

	    *U = 0.0;

	  }
	while (U->End ());



	for (size_t i = 0; bhmp[i] != 0.0 && i < NumBH; i++)
	  {



	    sbl =
	      sqrt (bhsx[i] * bhsx[i] + bhsy[i] * bhsy[i] +
		    bhsz[i] * bhsz[i]);

	    imass = 1.0 / bhmp[i];

	    imass2 = imass * imass;

	    Jx = 4.0 * bhsx[i] * imass2;

	    Jy = 4.0 * bhsy[i] * imass2;

	    Jz = 4.0 * bhsz[i] * imass2;

	    J2 = Jx * Jx + Jy * Jy + Jz * Jz;





	    pbl =
	      sqrt (bhpx[i] * bhpx[i] + bhpy[i] * bhpy[i] +
		    bhpz[i] * bhpz[i]);


	    Px = 2.0 * bhpx[i] * imass;

	    Py = 2.0 * bhpy[i] * imass;

	    Pz = 2.0 * bhpz[i] * imass;


	    P2 = Px * Px + Py * Py + Pz * Pz;


	    U->Set_Iterator (ALL);

	    do
	      {


		xi = U->Get_x () - bhx[i];

		yi = U->Get_y () - bhy[i];

		zi = U->Get_z () - bhz[i];

		r = sqrt (xi * xi + yi * yi + zi * zi);


		if (fabs (r) <= 0.00000000001)

		  r = 0.0000000001;

		Rx = 2.0 * xi * imass;

		Ry = 2.0 * yi * imass;

		Rz = 2.0 * zi * imass;


		R = 2.0 * r * imass;

		l1 = 1.0 / (1.0 + R);




		l2 = l1 * l1;

		l3 = l2 * l1;

		l4 = l3 * l1;

		l5 = l4 * l1;



		u0j = 0.025 * (l1 + l2 + l3 - 4.0 * l4 + 2.0 * l5);


		u2j = -0.05 * l5;

		if (sbl != 0.0 && r != 0.0)

		  muj =
		    (bhsx[i] * xi + bhsy[i] * yi + bhsz[i] * zi) / (r * sbl);

		else

		  muj = 0.0;

		PL2j = 0.5 * (3.0 * muj * muj - 1.0);



		uj = J2 * (u0j + u2j * R * R * PL2j);




		u0p = 0.15625 * (l1 - 2.0 * l2 + 2.0 * l3 - l4 + 0.2 * l5);


		u2p = 0.0125 * (15.0 * l1 + 132.0 * l2 + 53.0 * l3 + 96.0 * l4
				+ 82 * l5 + 84.0 * (l5 +
						    log (l1) / R) / R) / R;


		if (pbl != 0.0 && r != 0.0)

		  mup =
		    (bhpx[i] * xi + bhpy[i] * yi + bhpz[i] * zi) / (r * pbl);

		else

		  mup = 0.0;

		PL2p = 0.5 * (3.0 * mup * mup - 1.0);

		up = P2 * (u0p + u2p * PL2p);



		PXJ[0] = Py * Jz - Pz * Jy;

		PXJ[1] = Pz * Jx - Px * Jz;

		PXJ[2] = Px * Jy - Py * Jx;


		uc =
		  0.0125 * (PXJ[0] * Rx + PXJ[1] * Ry + PXJ[2] * Rz) * (1.0 +
									5.0 *
									R +
									10.0 *
									R *
									R) *
		  l5;


		*U += uj + up + uc;


		  

	      }
	    while (U->End ());
	    
	  }


	cout << " == Approximate ID: level -" << lev << "- box -" << ibox << "-" << endl << endl;



      }



  Get_Puncture_App();

}





void elliptic::Get_Puncture_App()
{



  size_t lev = levels-1;

      for (size_t i = 0; bhmp[i] != 0.0 && i < NumBH; i++)
	u_at_p[i] = 0.0;


  for (size_t ibox = 0; ibox < vars.Get_num_dom (lev); ibox++)
    {


      double r, xi, yi, zi;

      double J2, Jx, Jy, Jz;

      double P2, Px, Py, Pz;


      double sbl;

      double pbl;


      double uj, u0j, u2j;

      double up, u0p, u2p;

      double l1, R, muj, mup, PL2j, PL2p;

      double Rx, Ry, Rz;

      double l2, l3, l4, l5;

      double PXJ[3],
	uc;


      double imass, imass2;



      ffunction *U = vars.Get_ffunction (lev, ibox, "u");


      for (size_t i = 0; bhmp[i] != 0.0 && i < NumBH; i++)
	{



	  sbl = sqrt (bhsx[i] * bhsx[i] + bhsy[i] * bhsy[i] + bhsz[i] * bhsz[i]);

	  imass = 1.0 / bhmp[i];

	  imass2 = imass * imass;

	  Jx = 4.0 * bhsx[i] * imass2;

	  Jy = 4.0 * bhsy[i] * imass2;

	  Jz = 4.0 * bhsz[i] * imass2;

	  J2 = Jx * Jx + Jy * Jy + Jz * Jz;





	  pbl = sqrt (bhpx[i] * bhpx[i] + bhpy[i] * bhpy[i] + bhpz[i] * bhpz[i]);


	  Px = 2.0 * bhpx[i] * imass;

	  Py = 2.0 * bhpy[i] * imass;

	  Pz = 2.0 * bhpz[i] * imass;


	  P2 = Px * Px + Py * Py + Pz * Pz;



	  for (size_t j = 0; bhmp[j] != 0.0 && j < NumBH; j++)
	    if(U->Is_In(bhx[j],bhy[j],bhz[j]))
	      {

		size_t indx_i = U->Get_i(bhx[j]);
		size_t indx_j = U->Get_j(bhy[j]);
		size_t indx_k = U->Get_k(bhz[j]);
	  

		xi = U->Get_x(indx_i,indx_j,indx_k) - bhx[i];
	    
		yi = U->Get_y(indx_i,indx_j,indx_k) - bhy[i];

		zi = U->Get_z(indx_i,indx_j,indx_k) - bhz[i];


		r = sqrt (xi * xi + yi * yi + zi * zi);

		
		if (fabs (r) <= 0.00000000001)
		  
		  r = 0.0000000001;
	    
		Rx = 2.0 * xi * imass;

		Ry = 2.0 * yi * imass;

		Rz = 2.0 * zi * imass;


		R = 2.0 * r * imass;

		l1 = 1.0 / (1.0 + R);




		l2 = l1 * l1;

		l3 = l2 * l1;

		l4 = l3 * l1;

		l5 = l4 * l1;



		u0j = 0.025 * (l1 + l2 + l3 - 4.0 * l4 + 2.0 * l5);


		u2j = -0.05 * l5;

		if (sbl != 0.0 && r != 0.0)

		  muj = (bhsx[i] * xi + bhsy[i] * yi + bhsz[i] * zi) / (r * sbl);

	  else

	    muj = 0.0;

	  PL2j = 0.5 * (3.0 * muj * muj - 1.0);



	  uj = J2 * (u0j + u2j * R * R * PL2j);




	  u0p = 0.15625 * (l1 - 2.0 * l2 + 2.0 * l3 - l4 + 0.2 * l5);


	  u2p = 0.0125 * (15.0 * l1 + 132.0 * l2 + 53.0 * l3 + 96.0 * l4
			  + 82 * l5 + 84.0 * (l5 +
					      log (l1) / R) / R) / R;


	  if (pbl != 0.0 && r != 0.0)

	    mup =
	      (bhpx[i] * xi + bhpy[i] * yi + bhpz[i] * zi) / (r * pbl);

	  else

	    mup = 0.0;

	  PL2p = 0.5 * (3.0 * mup * mup - 1.0);

	  up = P2 * (u0p + u2p * PL2p);



	  PXJ[0] = Py * Jz - Pz * Jy;

	  PXJ[1] = Pz * Jx - Px * Jz;

	  PXJ[2] = Px * Jy - Py * Jx;


	  uc = 0.0125 * (PXJ[0] * Rx + PXJ[1] * Ry + PXJ[2] * Rz) * (1.0 + 5.0 * R + 10.0 * R * R) * l5;


	  u_at_p[j] += uj + up + uc;

	      }
		  

 		    

	}





    }



  for (size_t i = 0;  i < NumBH; i++)
    cout << " == Approximate ID: puncture " << i <<  endl
	 << " x : " << bhx[i] << "\t"
	 << " y : " << bhy[i] << "\t"
	 << " z : " << bhz[i] << "\t"
	 << " u : " << u_at_p[i] << endl;


}


void elliptic::Compute_ADM_Mass(){


  for(int i=0; i<bhmp.size(); i++){

    ffunction *U = vars.Get_ffunction (levels-1, i, "u");

    U->Set_Interpol(LAGRANGE,4);

    double bhmj = 0;    

    for(int j=0; j<bhmp.size(); j++) 
      if(j!=i)
	bhmj += 0.5*bhmp[j]/sqrt( (bhx[i]-bhx[j])*(bhx[i]-bhx[j]) + (bhy[i]-bhy[j])*(bhy[i]-bhy[j]) + (bhz[i]-bhz[j])*(bhz[i]-bhz[j]));

    ADM_mass[i] = bhmp[i]*(1+bhmj+U->Interpol_In(bhx[i],bhy[i],bhz[i])); 
    
  }

  cout << "Puncure method: " << endl;

  for(int i=0; i< ADM_mass.size(); i++)
    cout << "ADM mass[" << i << "] = " << ADM_mass[i] << endl;

  for (size_t lev = 0; lev < levels; lev++)

    for (size_t ibox = 0; ibox < vars.Get_num_dom (lev); ibox++)
      {

	ffunction *U = vars.Get_ffunction (lev, ibox, "u");
	ffunction *C = vars.Get_ffunction (lev, ibox, "c");
	ffunction *dU = vars.Get_ffunction (lev, ibox, "du");

	double cx = U->Get_domain()->Get_Center_x();
	double cy = U->Get_domain()->Get_Center_y();
	double cz = U->Get_domain()->Get_Center_z();

	U->Set_Iterator (ALL);
	dU->Set_Iterator(ALL);
	C->Set_Iterator (ALL);

	do{
	  
	  double x = U->Get_X()-cx;
	  double y = U->Get_Y()-cy;
	  double z = U->Get_Z()-cz;

	  double ir = 1/sqrt(x*x+y*y+z*z);
	  
	  *dU = ir*( x*(C->Dx()+U->Dx()) + y*(C->Dy()+U->Dy()) + z*(C->Dz()+U->Dz()));
	  

	}while(U->End() && dU->End() && C->End() );
	

	double r0 = 0.375*min(U->Get_domain()->Get_Length_x(),min(U->Get_domain()->Get_Length_y(),U->Get_domain()->Get_Length_z()));
	

	int N_theta = 4*25+1;
	int N_phi = N_theta;

	double dtheta = M_PI/(N_theta-1);
	double dphi = 2*M_PI/(N_phi-1);

	double Sph_data[N_theta][N_phi];

	for(int itheta = 0; itheta < N_theta; itheta++){

	  double theta = itheta*dtheta;

	  for(int iphi = 0; iphi < N_phi; iphi++)
	    {
	
	      double phi = iphi*dphi;
	      
	      double sph_x = cx+r0*cos(phi)*sin(theta);
	      double sph_y = cy+r0*sin(phi)*sin(theta);
	      double sph_z = cz+r0*cos(theta);

	      Sph_data[itheta][iphi] = dU->Interpol_In(sph_x,sph_y,sph_z)*sin(theta);


	    }


	}



	double ADMmass = 0;
	double IPhi0, IPhi1, IPhi2, IPhi3, IPhi4;
	/*			
	for(int i = 0; i < N_theta-1; i+=2){

	  IPhi0 = 0;
	  IPhi1 = 0;
	  IPhi2 = 0;

	  for(int j = 0; j < N_phi-1; j+=2){
	    
	    IPhi0 += Sph_data[i][j]+4*Sph_data[i][j+1]+Sph_data[i][j+2];

	    IPhi1 += Sph_data[i+1][j]+4*Sph_data[i+1][j+1]+Sph_data[i+1][j+2];

	    IPhi2 += Sph_data[i+2][j]+4*Sph_data[i+2][j+1]+Sph_data[i+2][j+2];


	  }	

	  IPhi0 *= dphi/3;
	  IPhi1 *= dphi/3;
	  IPhi2 *= dphi/3;

	  ADMmass += (IPhi0+4*IPhi1+IPhi2)*dtheta/3;

	}
	
	ADMmass *= -0.5/M_PI;
	*/


	
	for(int i = 0; i < N_theta-3; i+=4){

	  IPhi0=0;
	  IPhi1=0;
	  IPhi2=0;
	  IPhi3=0;
	  IPhi4=0;

	  for(int j = 0; j < N_phi-3; j+=4){
	    
	    IPhi0 += 7*Sph_data[i][j]+32*Sph_data[i][j+1]+12*Sph_data[i][j+2]+32*Sph_data[i][j+3]+7*Sph_data[i][j+4];

	    IPhi1 += 7*Sph_data[i+1][j]+32*Sph_data[i+1][j+1]+12*Sph_data[i+1][j+2]+32*Sph_data[i+1][j+3]+7*Sph_data[i+1][j+4];

	    IPhi2 += 7*Sph_data[i+2][j]+32*Sph_data[i+2][j+1]+12*Sph_data[i+2][j+2]+32*Sph_data[i+2][j+3]+7*Sph_data[i+2][j+4];

	    IPhi3 += 7*Sph_data[i+3][j]+32*Sph_data[i+3][j+1]+12*Sph_data[i+3][j+2]+32*Sph_data[i+3][j+3]+7*Sph_data[i+3][j+4];

	    IPhi4 += 7*Sph_data[i+4][j]+32*Sph_data[i+4][j+1]+12*Sph_data[i+4][j+2]+32*Sph_data[i+4][j+3]+7*Sph_data[i+4][j+4];


	  }	

	  IPhi0 *= 2*dphi/45;
	  IPhi1 *= 2*dphi/45;
	  IPhi2 *= 2*dphi/45;
	  IPhi3 *= 2*dphi/45;
	  IPhi4 *= 2*dphi/45;


	  ADMmass += 7*IPhi0+32*IPhi1+12*IPhi2+32*IPhi3+7*IPhi4;

	}
	

	ADMmass *= -r0*r0*dtheta/(45*M_PI);
	
	cout << "Level " << lev << " box " << ibox << " r0= "<< r0 <<" AMDmass= " << ADMmass << endl;


      }

}




void elliptic::Compute_Solution(){


  /*
  for (size_t lev = 0; lev < levels; lev++)

    for (size_t ibox = 0; ibox < vars.Get_num_dom (lev); ibox++)
      {

	ffunction *V1 = vars.Get_ffunction (lev, ibox, "V1");
	ffunction *V2 = vars.Get_ffunction (lev, ibox, "V2");
	ffunction *V3 = vars.Get_ffunction (lev, ibox, "V3");
	ffunction *lambda = vars.Get_ffunction (lev, ibox, "lambda");

	ffunction *X1 = vars.Get_ffunction (lev, ibox, "X1");
	ffunction *X2 = vars.Get_ffunction (lev, ibox, "X2");
	ffunction *X3 = vars.Get_ffunction (lev, ibox, "X3");

	vars.Set_Iterator (var_names, ALL);
	vars.Set_Iterator (coeff_names, ALL);


	do{

	  *X1 = V1->Get_val()-0.25*lambda->Dx();	  
	  *X2 = V2->Get_val()-0.25*lambda->Dy();	  
	  *X3 = V3->Get_val()-0.25*lambda->Dz();	  
	  
	    
	}while(vars.End(coeff_names) && vars.End(var_names) );
	

      }

      */
}

