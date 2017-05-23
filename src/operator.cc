//
// operator.cc
//  
// Made by Pablo Galaviz
// Login   <pablo@NerV>
// 
// Started on  Tue Oct 27 11:15:41 2009 Pablo Galaviz
// Started on  Tue Oct 27 11:15:41 2009 Pablo Galaviz
//



#include "elliptic.h"







void
elliptic::Set_Equation ()
{



  total_names.clear ();



  var_names.clear ();

  coeff_names.clear ();

  source_names.clear ();

  aux_names.clear ();


  Operator.clear ();

  duOperator.clear ();


  switch (equation)
    {


    case POISSON:



      var_names.push_back ("u");



      Operator.push_back (&elliptic::poisson_OP);

      duOperator.push_back (&elliptic::poisson_duOP);



      break;

    case BOUNDARY_TEST:


      var_names.push_back ("u");


      Operator.push_back (&elliptic::poisson_OP);

      duOperator.push_back (&elliptic::poisson_duOP);





      break;


    case BRILL_WAVES:


      var_names.push_back ("u");


      coeff_names.push_back ("Lq");


      Operator.push_back (&elliptic::brill_OP);

      duOperator.push_back (&elliptic::brill_duOP);





      break;



    case PUNCTURES:


      var_names.push_back ("u");


      coeff_names.push_back ("b");

      coeff_names.push_back ("c");

      coeff_names.push_back ("du");


      Operator.push_back (&elliptic::puncture_OP);

      duOperator.push_back (&elliptic::puncture_duOP);

      Info_Punctures ();



      break;


    case PUNCTURES_SCALAR:


      var_names.push_back ("u");


      coeff_names.push_back ("b");

      coeff_names.push_back ("c");

      coeff_names.push_back ("d");

      coeff_names.push_back ("e");

      coeff_names.push_back ("du");


      Operator.push_back (&elliptic::puncture_SF_OP);

      duOperator.push_back (&elliptic::puncture_SF_duOP);

      Info_Punctures ();



      break;


    case GAUSS_BONNET:


      var_names.push_back ("psi");


      coeff_names.push_back ("phi");

      coeff_names.push_back ("phi_dx");
      coeff_names.push_back ("phi_dy");
      coeff_names.push_back ("phi_dz");

      coeff_names.push_back ("phi_dxdx");
      coeff_names.push_back ("phi_dydx");
      coeff_names.push_back ("phi_dzdx");

      coeff_names.push_back ("phi_dxdy");
      coeff_names.push_back ("phi_dydy");
      coeff_names.push_back ("phi_dzdy");

      coeff_names.push_back ("phi_dxdz");
      coeff_names.push_back ("phi_dydz");
      coeff_names.push_back ("phi_dzdz");

      
      Operator.push_back (&elliptic::puncture_GB_OP);

      duOperator.push_back (&elliptic::puncture_GB_duOP);

      Info_Punctures ();



      break;


    case PUNCTURES_EM:


      var_names.push_back ("u");

      coeff_names.push_back ("A1");

      coeff_names.push_back ("A2");

      coeff_names.push_back ("A3");

      coeff_names.push_back ("B");

      coeff_names.push_back ("B1x");

      coeff_names.push_back ("B1y");

      coeff_names.push_back ("B1z");

      coeff_names.push_back ("B2x");

      coeff_names.push_back ("B2y");

      coeff_names.push_back ("B2z");

      coeff_names.push_back ("C");

      coeff_names.push_back ("D");

      coeff_names.push_back ("E");



      Operator.push_back (&elliptic::puncture_EM_OP);

      duOperator.push_back (&elliptic::puncture_EM_duOP);

      Info_Punctures ();



      break;




    case NSO_ID:


      var_names.push_back ("psi");

      coeff_names.push_back ("rhoADM");

      coeff_names.push_back ("psi0");


      coeff_names.push_back ("matter_rho");

      coeff_names.push_back ("matter_epsl");

      coeff_names.push_back ("alpha");


      Operator.push_back (&elliptic::NSO_OP);

      duOperator.push_back (&elliptic::NSO_duOP);

      // Info_NSO();



      break;



      /*
    case METRIC_TEST:


      var_names.push_back ("u");

      coeff_names.push_back ("rho");


      Operator.push_back (&elliptic::MT_OP);

      duOperator.push_back (&elliptic::MT_duOP);

      

      break;
      */




    case MULTIBOX_TEST:


      var_names.push_back ("u");


      Operator.push_back (&elliptic::poisson_OP);

      duOperator.push_back (&elliptic::poisson_duOP);





      break;



    case GENERIC:


      var_names.push_back ("u");


      Operator.push_back (&elliptic::poisson_OP);

      duOperator.push_back (&elliptic::poisson_duOP);





      break;



    case SYSTEM_TEST:



      var_names.push_back ("u");

      var_names.push_back ("v");


      coeff_names.push_back ("au");

      coeff_names.push_back ("av");


      Operator.push_back (&elliptic::system_OP_u);

      Operator.push_back (&elliptic::system_OP_v);

      duOperator.push_back (&elliptic::system_duOP_u);

      duOperator.push_back (&elliptic::system_duOP_v);





      break;


    case TRUMPET_1pLOG:



      var_names.push_back ("xi");



      Operator.push_back (&elliptic::trumpet_OP);

      duOperator.push_back (&elliptic::trumpet_duOP);



      break;



    case BOWEN_YORK:



      var_names.push_back ("V1");

      var_names.push_back ("V2");

      var_names.push_back ("V3");

      var_names.push_back ("lambda");


      Operator.push_back (&elliptic::BowenYork_OP_X);

      Operator.push_back (&elliptic::BowenYork_OP_Y);

      Operator.push_back (&elliptic::BowenYork_OP_Z);

      Operator.push_back (&elliptic::BowenYork_OP_L);


      duOperator.push_back (&elliptic::BowenYork_duOP_X);

      duOperator.push_back (&elliptic::BowenYork_duOP_Y);

      duOperator.push_back (&elliptic::BowenYork_duOP_Z);

      duOperator.push_back (&elliptic::BowenYork_duOP_L);




      break;




    case BY_EM:


      var_names.push_back ("u"); //0

      var_names.push_back ("V1"); //1

      var_names.push_back ("V2"); //2

      var_names.push_back ("V3"); //3

      var_names.push_back ("lambda"); //4


      coeff_names.push_back ("A1");//0

      coeff_names.push_back ("A2");//1

      coeff_names.push_back ("A3");//2

      coeff_names.push_back ("B");//3

      coeff_names.push_back ("B1x");//4

      coeff_names.push_back ("B1y");//5

      coeff_names.push_back ("B1z");//6

      coeff_names.push_back ("B2x");//7

      coeff_names.push_back ("B2y");//8

      coeff_names.push_back ("B2z");//9

      coeff_names.push_back ("C");//10

      coeff_names.push_back ("D");//11

      coeff_names.push_back ("E");//12

      coeff_names.push_back ("BiBi");//13 


      coeff_names.push_back ("Jx"); //14

      coeff_names.push_back ("Jy"); //15

      coeff_names.push_back ("Jz"); //16

      coeff_names.push_back ("eta"); //17

      coeff_names.push_back ("phi2"); //18


      coeff_names.push_back ("Axx");//19
      coeff_names.push_back ("Ayy");//20
      coeff_names.push_back ("Azz");//21

      coeff_names.push_back ("Axy");//22
      coeff_names.push_back ("Axz");//23
      coeff_names.push_back ("Ayz");//24


      Operator.push_back (&elliptic::BY_EM_OP_U);

      Operator.push_back (&elliptic::BY_EM_OP_X);

      Operator.push_back (&elliptic::BY_EM_OP_Y);

      Operator.push_back (&elliptic::BY_EM_OP_Z);

      Operator.push_back (&elliptic::BY_EM_OP_L);



      duOperator.push_back (&elliptic::BY_EM_duOP_U);

      duOperator.push_back (&elliptic::BY_EM_duOP_X);

      duOperator.push_back (&elliptic::BY_EM_duOP_Y);

      duOperator.push_back (&elliptic::BY_EM_duOP_Z);

      duOperator.push_back (&elliptic::BY_EM_duOP_L);




      Info_Punctures ();



      break;




    }



  size_t indx = 0;


  for (size_t i = 0; i < var_names.size (); i++)
    {

      total_names.push_back (var_names[i]);

      var_index.push_back (indx);

      indx++;


      stringstream stage_vars;

      stage_vars << i;

      aux_names.push_back ("r" + stage_vars.str ());

      aux_names.push_back ("v" + stage_vars.str ());

      aux_names.push_back ("w" + stage_vars.str ());

      aux_names.push_back ("ww" + stage_vars.str ());


      source_names.push_back ("rho" + stage_vars.str ());


    }




  for (size_t i = 0; i < coeff_names.size (); i++)
    {
      total_names.push_back (coeff_names[i]);

      coeff_index.push_back (indx);

      indx++;

    }



  for (size_t i = 0; i < source_names.size (); i++)

    total_names.push_back (source_names[i]);



  for (size_t i = 0; i < aux_names.size (); i++)

    total_names.push_back (aux_names[i]);



  vars.Make (D, dom, total_names);


  

     cout << "var names: " << endl;

     for(size_t i=0; i < var_names.size(); i++)

     cout << var_names[i] << " index: " << var_index[i] << "\n";

     cout << endl;


     cout << "source names: " << endl;

     for(size_t i=0; i < source_names.size(); i++)

     cout << source_names[i] << "\n";

     cout << endl;



     cout << "coeff names: " << endl;

     for(size_t i=0; i < coeff_names.size(); i++)

     cout << coeff_names[i] << " index: " << coeff_index[i] << "\n";

     cout << endl;


     cout << "aux names: " << endl;

     for(size_t i=0; i < aux_names.size(); i++)

     cout << aux_names[i] << "\n";

     cout << endl;



     cout << "total names: " << endl;

     for(size_t i=0; i < total_names.size(); i++)

     cout << total_names[i] << "\t";

     cout << endl;

   


}






void
elliptic::Set_Coefficients ()
{





  switch (equation)
    {


    case POISSON:

      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (Gaussian_source (x, y, z), source_names);

	}
      while (vars.End (source_names));



      Boundary.push_back (&elliptic::Gaussian);


      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);



      break;

    case BOUNDARY_TEST:


      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (test_boundary_source (x, y, z), "rho0");

	}
      while (vars.End (source_names));




      Boundary.push_back (&elliptic::test_boundary);



      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);




      break;


    case BRILL_WAVES:

      vars.Set_Iterator ("Lq", ALL);

      do
	{

	  double x = vars.Get_x ("Lq");

	  double y = vars.Get_y ("Lq");

	  double z = vars.Get_z ("Lq");

	  vars.Set_val (brill (x, y, z), "Lq");

	}
      while (vars.End ("Lq"));



      vars.Set_Iterator (source_names, ALL);

      do
	{

	  vars.Set_val (0.0, "rho0");

	}
      while (vars.End (source_names));



      Boundary.push_back (&elliptic::one);


      falloff_n.push_back (1);

      boundary_coeff.push_back (1.0);



      break;



    case PUNCTURES:

      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("b");

	  double y = vars.Get_y ("b");

	  double z = vars.Get_z ("b");

	  vars.Set_val (punctures_b (x, y, z), "b");

	  vars.Set_val (punctures_c (x, y, z), "c");

	}
      while (vars.End (coeff_names));


      vars.Set_Iterator (source_names, ALL);

      do
	{


	  vars.Set_val (0.0, "rho0");

	}
      while (vars.End (source_names));


      Boundary.push_back (&elliptic::zero);



      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);




      break;



    case PUNCTURES_SCALAR:

      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("b");

	  double y = vars.Get_y ("b");

	  double z = vars.Get_z ("b");

	  vars.Set_val (punctures_b (x, y, z), "b");

	  vars.Set_val (punctures_c (x, y, z), "c");
	  
	  /*
	  double r2=x*x+y*y+z*z;
	  r2 = r2==0? 0.000001 : r2;
	  double rho2=x*x+y*y;
	  rho2 = rho2==0? 0.000001 : rho2;
	  double r=sqrt(r2);

	  vars.Set_val (r*pow(2*z*sqrt(x*x+y*y)/r2,2)*pow((x*x-y*y)/rho2,2), "c");
	  */
	  vars.Set_val (punctures_scalar_d (x, y, z), "d");

	  vars.Set_val (punctures_scalar_e (x, y, z), "e");

	}
      while (vars.End (coeff_names));

      //Compute_ADM_Mass();
      
      //exit(0);

      vars.Set_Iterator (source_names, ALL);


      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (test_ps(x,y,z), "rho0");



	}
      while (vars.End (source_names));


      Boundary.push_back (&elliptic::zero);



      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);




      break;


    case GAUSS_BONNET:

      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("phi");

	  double y = vars.Get_y ("phi");

	  double z = vars.Get_z ("phi");

	  vars.Set_val (phi (x, y, z), "phi");

	  vars.Set_val (phi_dx (x, y, z), "phi_dx");
	  vars.Set_val (phi_dy (x, y, z), "phi_dy");
	  vars.Set_val (phi_dz (x, y, z), "phi_dz");

	  vars.Set_val (phi_dxdx (x, y, z), "phi_dxdx");
	  vars.Set_val (phi_dydx (x, y, z), "phi_dydx");
	  vars.Set_val (phi_dzdx (x, y, z), "phi_dzdx");

	  vars.Set_val (phi_dxdy (x, y, z), "phi_dxdy");
	  vars.Set_val (phi_dydy (x, y, z), "phi_dydy");
	  vars.Set_val (phi_dzdy (x, y, z), "phi_dzdy");

	  vars.Set_val (phi_dxdz (x, y, z), "phi_dxdz");
	  vars.Set_val (phi_dydz (x, y, z), "phi_dydz");
	  vars.Set_val (phi_dzdz (x, y, z), "phi_dzdz");
	  

	}
      while (vars.End (coeff_names));

      //Compute_ADM_Mass();
      
      //exit(0);

      vars.Set_Iterator (source_names, ALL);


      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (0, "rho0");



	}
      while (vars.End (source_names));


      Boundary.push_back (&elliptic::zero);



      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);




      break;

    case METRIC_TEST:

      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("b");

	  double y = vars.Get_y ("b");

	  double z = vars.Get_z ("b");

	  vars.Set_val (punctures_b (x, y, z), "b");

	  vars.Set_val (punctures_c (x, y, z), "c");
	  
	  /*
	  double r2=x*x+y*y+z*z;
	  r2 = r2==0? 0.000001 : r2;
	  double rho2=x*x+y*y;
	  rho2 = rho2==0? 0.000001 : rho2;
	  double r=sqrt(r2);

	  vars.Set_val (r*pow(2*z*sqrt(x*x+y*y)/r2,2)*pow((x*x-y*y)/rho2,2), "c");
	  */
	  vars.Set_val (punctures_scalar_d (x, y, z), "d");

	  vars.Set_val (punctures_scalar_e (x, y, z), "e");

	}
      while (vars.End (coeff_names));

      //Compute_ADM_Mass();
      
      //exit(0);

      vars.Set_Iterator (source_names, ALL);


      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (test_ps(x,y,z), "rho0");



	}
      while (vars.End (source_names));


      Boundary.push_back (&elliptic::zero);



      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);




      break;

    case PUNCTURES_EM:

      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("A1");

	  double y = vars.Get_y ("A1");

	  double z = vars.Get_z ("A1");

	  //index ---> var 
	  // 0  A1
	  // 1  A2
	  // 2  A3
	  // 3 B
	  // 4 B1x
	  // 5 B1y
	  // 6 B1z 
	  // 7 B2x
	  // 8 B2y
	  // 9 B2z 
	  // 10 C
	  // 11 D 
	  // 12 E

	  vars.Set_val (EM_BH_A1 (x, y, z), "A1");
	  vars.Set_val (EM_BH_A2 (x, y, z), "A2");
	  vars.Set_val (EM_BH_A3 (x, y, z), "A3");

	  vars.Set_val (EM_BH_B (x, y, z), "B");

	  vars.Set_val (EM_BH_B1x (x, y, z), "B1x");
	  vars.Set_val (EM_BH_B1y (x, y, z), "B1y");
	  vars.Set_val (EM_BH_B1z (x, y, z), "B1z");
	  
	  vars.Set_val (EM_BH_B2x (x, y, z), "B2x");
	  vars.Set_val (EM_BH_B2y (x, y, z), "B2y");
	  vars.Set_val (EM_BH_B2z (x, y, z), "B2z");

	  vars.Set_val (EM_BH_C (x, y, z), "C");

	  vars.Set_val (EM_BH_D (x, y, z), "D");

	  vars.Set_val (EM_BH_E (x, y, z), "E");

	}
      while (vars.End (coeff_names));

       //Compute_ADM_Mass();
      
      //exit(0);

      vars.Set_Iterator (source_names, ALL);


      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (0.0, "rho0");



	}
      while (vars.End (source_names));


      Boundary.push_back (&elliptic::one);


      falloff_n.push_back (1);

      boundary_coeff.push_back (1.0);



      break;



    case NSO_ID:

      set_TOV();


      add_perturbation();

      exit(0);

      vars.Set_Iterator (source_names, ALL);

      vars.Set_Iterator("psi0",ALL);

      vars.Set_Iterator("psi",ALL);


      do
	{


	  vars.Set_val (0.0, "rho0");
	  
	  vars.Set_val(vars.Get_val("psi0"),"psi");
	  
	  //vars.Set_val(1.0,"psi");


	}
      while (vars.End (source_names)&&vars.End ("psi0")&&vars.End ("psi") );


      Boundary.push_back (&elliptic::one);



      falloff_n.push_back (1);

      boundary_coeff.push_back (1.0);


      break;



    case MULTIBOX_TEST:

      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (test_mfunct_source (x, y, z), "rho0");

	}
      while (vars.End (source_names));




      Boundary.push_back (&elliptic::test_mfunct);


      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);



      break;



    case GENERIC:

      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (test_mfunct_source (x, y, z), "rho0");

	}
      while (vars.End (source_names));



      Boundary.push_back (&elliptic::test_mfunct);



      falloff_n.push_back (1);

      boundary_coeff.push_back (0.0);



      break;






    case SYSTEM_TEST:


      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (system_u_source (x, y, z), "rho0");

	  vars.Set_val (system_v_source (x, y, z), "rho1");



	}
      while (vars.End (source_names));



      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("au");

	  double y = vars.Get_y ("au");

	  double z = vars.Get_z ("au");


	  vars.Set_val (system_u_a (x, y, z), "au");

	  vars.Set_val (system_v_a (x, y, z), "av");


	}
      while (vars.End (coeff_names));






      Boundary.push_back (&elliptic::system_u);

      Boundary.push_back (&elliptic::system_v);


      falloff_n.push_back (1);

      falloff_n.push_back (2);


      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (1.0);


      break;





    case TRUMPET_1pLOG:

      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  vars.Set_val (Trumpet_source (x, y, z), source_names);

	}
      while (vars.End (source_names));



      Boundary.push_back (&elliptic::Trumpet);


      falloff_n.push_back (1);

      boundary_coeff.push_back (1.0);



      break;



    case BOWEN_YORK:

      /*
      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("sol_X1");

	  double y = vars.Get_y ("sol_X1");

	  double z = vars.Get_z ("sol_X1");

	  vars.Set_val (Bowen_York_X1 (x, y, z), "sol_X1");

	  vars.Set_val (Bowen_York_X2 (x, y, z), "sol_X2");

	  vars.Set_val (Bowen_York_X3 (x, y, z), "sol_X3");


	}
      while (vars.End (coeff_names));
      */

      vars.Set_Iterator (source_names, ALL);

      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  double dx = vars.Get_dx ("rho0");

	  double dy = vars.Get_dy ("rho0");

	  double dz = vars.Get_dz ("rho0");


	  vars.Set_val (Bowen_York_rho0(x,y,z,dx,dy,dz), "rho0");

	  vars.Set_val (Bowen_York_rho1(x,y,z,dx,dy,dz), "rho1");

	  vars.Set_val (Bowen_York_rho2(x,y,z,dx,dy,dz), "rho2");

	  //vars.Set_val (Bowen_York_rho0_test(x,y,z), "rho0");

	  //vars.Set_val (Bowen_York_rho1_test(x,y,z), "rho1");

	  //vars.Set_val (Bowen_York_rho2_test(x,y,z), "rho2");

	  vars.Set_val (0.0, "rho3");


	}
      while (vars.End (source_names));




      Boundary.push_back (&elliptic::zero);

      Boundary.push_back (&elliptic::zero);

      Boundary.push_back (&elliptic::zero);

      Boundary.push_back (&elliptic::zero);


      falloff_n.push_back (2);

      falloff_n.push_back (2);

      falloff_n.push_back (2);

      falloff_n.push_back (2);


      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (0.0);



      break;



    case BY_EM:

      

      vars.Set_Iterator (coeff_names, ALL);

      do
	{

	  double x = vars.Get_x ("A1");

	  double y = vars.Get_y ("A1");

	  double z = vars.Get_z ("A1");

	  //index ---> var 
	  // 0 A1
	  // 1 A2
	  // 2 A3
	  // 3 B
	  // 4 B1x
	  // 5 B1y
	  // 6 B1z 
	  // 7 B2x
	  // 8 B2y
	  // 9 B2z 
	  // 10 C
	  // 11 D 
	  // 12 E
	  // 13 BiBi
	  //14 Jx
	  //15 Jy
	  //16 Jz
	  //17 eta
	  //18 phi2/4
	  //19 Axx
	  //20 Ayy
	  //21 Azz
	  //22 Axy
	  //23 Axz
	  //24 Ayz

	  vars.Set_val (BY_EM_A1 (x, y, z), "A1");
	  vars.Set_val (BY_EM_A2 (x, y, z), "A2");
	  vars.Set_val (BY_EM_A3 (x, y, z), "A3");

	  vars.Set_val (BY_EM_B (x, y, z), "B");

	  vars.Set_val (BY_EM_B1x (x, y, z), "B1x");
	  vars.Set_val (BY_EM_B1y (x, y, z), "B1y");
	  vars.Set_val (BY_EM_B1z (x, y, z), "B1z");
	  
	  vars.Set_val (BY_EM_B2x (x, y, z), "B2x");
	  vars.Set_val (BY_EM_B2y (x, y, z), "B2y");
	  vars.Set_val (BY_EM_B2z (x, y, z), "B2z");

	  vars.Set_val (BY_EM_C (x, y, z), "C");

	  vars.Set_val (BY_EM_D (x, y, z), "D");

	  vars.Set_val (BY_EM_E (x, y, z), "E");


	  vars.Set_val (BY_EM_BiBi (x, y, z), "BiBi");

	  vars.Set_val (BY_EM_Jx (x, y, z), "Jx");

	  vars.Set_val (BY_EM_Jy (x, y, z), "Jy");

	  vars.Set_val (BY_EM_Jz (x, y, z), "Jz");


	  vars.Set_val (BY_EM_eta (x, y, z), "eta");

	  vars.Set_val (0.25*pow(BY_EM_phi (x, y, z),2), "phi2");


	  double Axx=0, Ayy=0, Azz=0;
	  double Axy=0, Axz=0, Ayz=0;

	  for (size_t i = 0; i < NumBH; i++)
	    {

	      double X = x - bhx[i];

	      double Y = y - bhy[i];

	      double Z = z - bhz[i];

	      double x2 = X * X;

	      double y2 = Y * Y;

	      double z2 = Z * Z;

	      double r2 = x2 + y2 + z2;

	      if (r2 < 2 * numeric_limits < double >::epsilon ())

		r2 = 2 * numeric_limits < double >::epsilon ();


	      double r = sqrt (r2);

	      double ir5 = 1.0 / (r2 * r2 * r);

	      Axx += 1.5 * ir5 * (bhpx[i] * X * (x2 + r2) +
				  bhpy[i] * Y * (x2 - r2) +
				  bhpz[i] * Z * (x2 - r2) +
				  4.0 * X * (bhsy[i] * Z - bhsz[i] * Y));

	      Ayy += 1.5 * ir5 * (bhpx[i] * X * (y2 - r2) +
				  bhpy[i] * Y * (y2 + r2) +
				  bhpz[i] * Z * (y2 - r2) +
				  4.0 * Y * (bhsz[i] * X - bhsx[i] * Z));


	      Azz += 1.5 * ir5 * (bhpx[i] * X * (z2 - r2) +
				  bhpy[i] * Y * (z2 - r2) +
				  bhpz[i] * Z * (z2 + r2) +
				  4.0 * Z * (bhsx[i] * Y - bhsy[i] * Z));

	      Axy += 1.5 * ir5 * (bhpx[i] * Y * (r2 + x2) +
				  bhpy[i] * X * (r2 + y2) +
				  bhpz[i] * X * Y * Z -
				  2.0 * bhsx[i] * X * Z +
				  2.0 * bhsy[i] * Y * Z +
				  2.0 * bhsz[i] * (x2 - y2));

	      Axz += 1.5 * ir5 * (bhpx[i] * Z * (r2 + x2) +
				  bhpy[i] * X * Y * Z +
				  bhpz[i] * X * (r2 + z2) +
				  2.0 * bhsx[i] * X * Y +
				  2.0 * bhsy[i] * (z2 - x2) -
				  2.0 * bhsz[i] * Y * Z);

	      Ayz += 1.5 * ir5 * (bhpx[i] * X * Y * Z +
				  bhpy[i] * Z * (r2 + y2) +
				  bhpz[i] * Y * (r2 + z2) +
				  2.0 * bhsx[i] * (y2 - z2) -
				  2.0 * bhsy[i] * X * Y + 2.0 * bhsz[i] * X * Z);


	    }


	  vars.Set_val (Axx, "Axx");
	  vars.Set_val (Ayy, "Ayy");
	  vars.Set_val (Azz, "Azz");

	  vars.Set_val (Axy, "Axy");
	  vars.Set_val (Axz, "Axz");
	  vars.Set_val (Ayz, "Ayz");



	}
      while (vars.End (coeff_names));
      
       //Compute_ADM_Mass();
      
      //exit(0);

      vars.Set_Iterator (source_names, ALL);


      do
	{

	  double x = vars.Get_x ("rho0");

	  double y = vars.Get_y ("rho0");

	  double z = vars.Get_z ("rho0");

	  double dx = vars.Get_dx ("rho0");

	  double dy = vars.Get_dy ("rho0");

	  double dz = vars.Get_dz ("rho0");


	  vars.Set_val (0.0, "rho0");

	  vars.Set_val (0.0, "rho1");

	  vars.Set_val (0.0, "rho2");

	  vars.Set_val (0.0, "rho3");

	  vars.Set_val (0.0, "rho4");



	}
      while (vars.End (source_names));


      // boundary u 
      Boundary.push_back (&elliptic::one);

      falloff_n.push_back (1);

      boundary_coeff.push_back (1.0);


      // boundary V1, V2, V3, lambda

      Boundary.push_back (&elliptic::zero);

      Boundary.push_back (&elliptic::zero);

      Boundary.push_back (&elliptic::zero);

      Boundary.push_back (&elliptic::zero);


      falloff_n.push_back (2);

      falloff_n.push_back (2);

      falloff_n.push_back (2);

      falloff_n.push_back (2);


      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (0.0);

      boundary_coeff.push_back (0.0);



      break;






    }


  Set_Boundary ();





}










//========================== functions definition =========================//




//============================== Brill waves =============================//


double
elliptic::brill (double x, double y, double z)
{



  const double Amp = 1.9;


  double z2 = z * z;

  double rho2 = x * x + y * y;

  double qrr =
    2.0 * Amp * exp (-rho2 - z2) * (1.0 - 5.0 * rho2 + 2.0 * rho2 * rho2);


  double qzz = 2.0 * Amp * exp (-rho2 - z2) * rho2 * (2.0 * z2 - 1.0);


  return (0.25 * (qrr + qzz));



}


//============================= Test function =============================//


double
elliptic::test_funct (double x, double y, double z)
{


  double r = sqrt (x * x + y * y + z * z);

/*
    double r2 = r*r;        

    double a = 0.05;

    double b = 2.0;
            
    
    return( exp( -a*r2 ) * cos(b*r) );
*/

  return (r * r * r * r * r);



}


double
elliptic::test_funct_source (double x, double y, double z)
{


  double r = sqrt (x * x + y * y + z * z);

/*
  double r2 = r*r;        

    double a = 0.05;

    double b = 2.0;
    
    return( exp( -a*r2 ) * ( r * cos( b * r ) * ( 4.0 * a * a * r2 - b * b - 6.0 * a ) +
                             2.0 * b * ( 2*a*r2-1 ) * sin( b * r ) ) / r );
*/


  return (30.0 * r * r * r);


}



//=================== Test function (for Robin boundary) ==================//


double
elliptic::test_boundary (double x, double y, double z)
{




  double r = sqrt (x * x + y * y + z * z);


  if (!dcomp (r, 0.0))

    return (0.004 * tanh (r) / r);

  else

    return (0.004);


}


double
elliptic::test_boundary_source (double x, double y, double z)
{




  double r = sqrt (x * x + y * y + z * z);


  if (!dcomp (r, 0.0))

    return (-2.0 * 0.004 * tanh (r) / (r * cosh (r) * cosh (r)));

  else

    return (-2.0 * 0.004);

}




//===================== Test function (for 2 boxes) =======================//


double
elliptic::test_mfunct (double x, double y, double z)
{



  double r2_1 = (x - bhx[0]) * (x - bhx[0]) + y * y + z * z;

  double r2_2 = (x - bhx[1]) * (x - bhx[1]) + y * y + z * z;


  double a = 0.05;


  return (exp (-a * r2_1) + exp (-a * r2_2));


}


double
elliptic::test_mfunct_source (double x, double y, double z)
{


  double r2_1 = (x - bhx[0]) * (x - bhx[0]) + y * y + z * z;

  double r2_2 = (x - bhx[1]) * (x - bhx[1]) + y * y + z * z;


  double a = 0.05;


  return (2.0 * a * exp (-a * r2_1) * ((2 * a * r2_1 - 3)) +
	  2.0 * a * exp (-a * r2_2) * ((2 * a * r2_2 - 3)));



}



//============================= Gaussian test =============================//


double
elliptic::Gaussian_source (double x, double y, double z)
{


  double r2 = x * x + y * y + z * z;

  double a = 0.004;

  return (a*exp (-r2)*(a+exp(0.5*r2)*(r2-3.0)));



}


double
elliptic::Gaussian (double x, double y, double z)
{

  double r2 = x * x + y * y + z * z;

  double a = 0.004;

  return (a*exp (-0.5*r2 ));



}



//=============================== 1/r^2 test ==============================//


double
elliptic::One_over_r2_source (double x, double y, double z)
{

  double r2 = x * x + y * y + z * z;


  return (-2.0 / (r2 * r2));


}


double
elliptic::One_over_r2 (double x, double y, double z)
{

  double r2 = x * x + y * y + z * z;

  return (2.0 / r2);


}




//===================== Test function (non-linear) ========================//




double
elliptic::Non_linear_source (double x, double y, double z)
{


  double A = 0.005;


  double Sz2 = 0.01;


  double r2 = x * x + y * y + z * z;

  double r = sqrt (r2);


  double b = 2.25 * Sz2 * (x * x + y * y) / pow (r, 8);

  double c = 1 + 0.5 / r;

  b = 5;
  c = 1;

  return (2.0 * A * exp (-r2) * (2.0 * r2 - 3.0) +
	  b / pow (A * exp (-r2) + c, 7));


/*
    return( 2.0 * A * exp(-r2) * (2.0 * r2 - 3.0) +

             pow(  A * exp(-r2) , 5 )  

            );
*/

}


double
elliptic::Non_linear (double x, double y, double z)
{


  double A = 0.005;



  double r2 = x * x + y * y + z * z;

  return (A * exp (-r2));



}

double
elliptic::Non_linear_b (double x, double y, double z)
{


  double Sz2 = 0.01;


  double r = sqrt (x * x + y * y + z * z);


  double b = 2.25 * Sz2 * (x * x + y * y) / pow (r, 8);


  b = 5;


  return (b);



}


double
elliptic::Non_linear_c (double x, double y, double z)
{



  double r = sqrt (x * x + y * y + z * z);



  double c = 1 + 0.5 / r;


  c = 1;

  return (c);



}



double
elliptic::Non_linear_d (double x, double y, double z)
{



  double r = sqrt (x * x + y * y + z * z);


  double pi2 = 8 * atan (1);

  if (r < 3)


    return (-pi2 * 0.1);

  else

    return (0.0);


}






//=============================== system test ==============================//



double
elliptic::system_u_a (double x, double y, double z)
{

  return (sin (sqrt (x * x + y * y + z * z)));

}



double
elliptic::system_v_a (double x, double y, double z)
{

  return (x * x - y * z);

}


double
elliptic::system_u_source (double x, double y, double z)
{

  double a = 0.004;

  double b = 1.0;

  double r2 = x * x + y * y + z * z;

  if (r2 < 2 * numeric_limits < double >::epsilon ())

    r2 = 2 * numeric_limits < double >::epsilon ();

  double r = sqrt (r2);

  double ir4 = 1.0 / (r2 * r2);

  return (a * tanh (r) *
	  (-2.0 * r2 * r * 1 / (cosh (r) * cos (r)) +
	   a * b * sin (r) * tanh (r) * (r2 - tanh (r2))) * ir4);



}



double
elliptic::system_v_source (double x, double y, double z)
{

  double a = 0.8;

  double b = 0.5;

  double r2 = x * x + y * y + z * z;

  if (r2 < 2 * numeric_limits < double >::epsilon ())

    r2 = 2 * numeric_limits < double >::epsilon ();

  double r = sqrt (r2);

  double r4 = r2 * r2;

  double ir6 = 1.0 / (r4 * r2);

  return (b *
	  (a * a * r2 * (x * x - y * z) * (1 / (cosh (r2) * cosh (r2))) *
	   (r2 * cosh (r2) - sinh (r2)) * pow (tanh (r),
					       2) - 2 * r2 * tanh (r2) +
	   2 * r2 * pow ((1 / (cosh (r2) * cosh (r2))),
			 2) * (r2 + 4 * r4 * tanh (r2))) * ir6);


}




double
elliptic::system_u (double x, double y, double z)
{


  double r = sqrt (x * x + y * y + z * z);

  double a = 0.004;

  if (!dcomp (r, 0.0))

    return (a * tanh (r) / r);

  else

    return (a);



}

double
elliptic::system_v (double x, double y, double z)
{

  double r2 = x * x + y * y + z * z;

  double b = 1.0;

  if (!dcomp (r2, 0.0))

    return (b * (1 - tanh (r2) / r2));

  else

    return (b);



}




//=============================== punctures ==============================//


double
elliptic::punctures_c (double X, double Y, double Z)
{


  double cval = 1;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      cval += 0.5 * bhmp[i] / r;

    }


  return (cval);


}



double
elliptic::punctures_b (double X, double Y, double Z)
{



  double A[3][3];


  for (int i = 0; i < 3; i++)

    for (int j = 0; j < 3; j++)

      A[i][j] = 0.0;



  for (size_t i = 0; i < NumBH; i++)
    {

      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];

      double x2 = x * x;

      double y2 = y * y;

      double z2 = z * z;

      double r2 = x2 + y2 + z2;

      if (r2 < 2 * numeric_limits < double >::epsilon ())

	r2 = 2 * numeric_limits < double >::epsilon ();


      double r = sqrt (r2);

      double ir5 = 1.0 / (r2 * r2 * r);

      A[0][0] += 1.5 * ir5 * (bhpx[i] * x * (x2 + r2) +
			      bhpy[i] * y * (x2 - r2) +
			      bhpz[i] * z * (x2 - r2) +
			      4.0 * x * (bhsy[i] * z - bhsz[i] * y));

      A[1][1] += 1.5 * ir5 * (bhpx[i] * x * (y2 - r2) +
			      bhpy[i] * y * (y2 + r2) +
			      bhpz[i] * z * (y2 - r2) +
			      4.0 * y * (bhsz[i] * x - bhsx[i] * z));


      A[2][2] += 1.5 * ir5 * (bhpx[i] * x * (z2 - r2) +
			      bhpy[i] * y * (z2 - r2) +
			      bhpz[i] * z * (z2 + r2) +
			      4.0 * z * (bhsx[i] * y - bhsy[i] * x));

      A[0][1] += 1.5 * ir5 * (bhpx[i] * y * (r2 + x2) +
			      bhpy[i] * x * (r2 + y2) +
			      bhpz[i] * x * y * z -
			      2.0 * bhsx[i] * x * z +
			      2.0 * bhsy[i] * y * z +
			      2.0 * bhsz[i] * (x2 - y2));

      A[0][2] += 1.5 * ir5 * (bhpx[i] * z * (r2 + x2) +
			      bhpy[i] * x * y * z +
			      bhpz[i] * x * (r2 + z2) +
			      2.0 * bhsx[i] * x * y +
			      2.0 * bhsy[i] * (z2 - x2) -
			      2.0 * bhsz[i] * y * z);

      A[1][2] += 1.5 * ir5 * (bhpx[i] * x * y * z +
			      bhpy[i] * z * (r2 + y2) +
			      bhpz[i] * y * (r2 + z2) +
			      2.0 * bhsx[i] * (y2 - z2) -
			      2.0 * bhsy[i] * x * y + 2.0 * bhsz[i] * x * z);


    }


  A[1][0] = A[0][1];

  A[2][0] = A[0][2];

  A[2][1] = A[1][2];


  double bval = 0.0;

  for (int i = 0; i < 3; i++)

    for (int j = 0; j < 3; j++)

      bval += A[i][j] * A[i][j];





  return (0.125 * bval);



}


//======================= punctures+scalar field ============================//


double elliptic::punctures_scalar_d (double X, double Y, double Z)
{

  double r = sqrt (X * X + Y * Y + Z * Z);

  if(ps_case == 1)
    return d_scalar_phi_cI(r);
  else
    if(ps_case == 2)
      return d_scalar_phi_cII(r);
    else
      if(ps_case == 3)
	return d_scalar_phi_ground(r);
      else
	return d_scalar_phi_cIV(r);


}


double elliptic::punctures_scalar_e (double X, double Y, double Z)
{


  if(ps_case == 2 || ps_case ==3 || ps_case ==4){

    double r = sqrt (X * X + Y * Y + Z * Z);

    double ff;

    if(ps_case ==2) ff = ps_e2*scalar_phi_cII(r);
    if(ps_case ==3) ff = ps_e2*scalar_phi_ground(r);
    if(ps_case ==4) ff = ps_e2*scalar_phi_cIV(r);

    return ps_e1*(1-exp(ff))*(1-exp(ff))*exp(-2*ff);
  }
  
  else

    return 0.0;

}


double
elliptic::test_ps (double x, double y, double z)
{

  double r2=x*x+y*y+z*z;
  
  //return -6*exp(-r2)+4*r2*exp(-r2)+punctures_b(x,y,z)/pow(punctures_c(x,y,z)+exp(-r2),7)+punctures_scalar_d(x,y,z)*(punctures_c(x,y,z)+exp(-r2))+punctures_scalar_e(x,y,z)*pow(punctures_c(x,y,z)+exp(-r2),5);
  return 0.0;
}

//===========================================NSO-TOV==============================


/* from bamba  */
/* set TOV star using dookie source code */
int elliptic::chebycoe_Extremes( valarray<double> &u, int n, int inv ) 
{
  
    int k, j, isignum, N=n-1;
    double fac, sum, PioN;
    
    valarray<double> c;

    c.resize( n );  
  
    PioN=4.0*atan(1.0)/N;

    if (inv == 0) {
    
      fac=2.0/N;
      isignum=1;
      
      for (j=0; j<n; j++) {
      
	sum=0.5*(u[0]+u[N]*isignum);
        
	for (k=1; k<N; k++)
	  sum += u[k]*cos(PioN*j*k);
	
	c[j]=fac*sum*isignum; 
	isignum = -isignum;
      }
       
      c[N] = 0.5*c[N];
    }

    else {

        for (j=0; j<n; j++) {

	  sum=-0.5*u[0];
	  isignum=1;
	  
	  for (k=0; k<n; k++){
          
	    sum += u[k]*cos(PioN*j*k)*isignum;
	    isignum = -isignum;
            
	  }
	  
	  c[j]=sum;
        }
    }
    for (j=0;j<n;j++)
      
      u[j]=c[j];

    return (1);

}

double elliptic::chebyeva( double a, double b, valarray<double> &c, int m, double x ) 
{
  
  double d=0.0,dd=0.0,sv,y,y2;
  int j;
  
  if ((x-a)*(x-b) > 0.0){ 

    //    printf( "ERROR in [interp.c] Routine 'chebyeva' :  x not in range [%g,%g].\n",a,b );
    
    error_exit("in [interp.c] Routine 'chebyeva':  x not in range","","",false);

    return (0);
  }

  y2=2.0*(y=(2.0*x-a-b)/(b-a));

  for (j=m-1;j>=1;j--) {

    sv=d;

    d=y2*d-dd+c[j];

    dd=sv;

  }

  return( y*d-dd+0.5*c[0] );

}



int elliptic::set_TOV()
{

    


  fstream File (tov_filename.c_str (),fstream::in);


  if (File.is_open ())
    {

      cout << "File found: " << tov_filename << endl;

      string readstr;

      string param;
      int dim, cols;
      double M,R ;



      File >> readstr; 
      
      param = readstr.substr(0,readstr.find("="));

      dim=atoi( readstr.substr(readstr.find("=")+1).c_str() );

      cout << param <<" \t: " << dim << endl;


      File >> readstr; 
      
      param = readstr.substr(0,readstr.find("="));

      cols=atoi( readstr.substr(readstr.find("=")+1).c_str() );

      cout << param <<" \t: " << cols << endl;



      File >> readstr; 
      
      param = readstr.substr(0,readstr.find("="));

      M=atof( readstr.substr(readstr.find("=")+1).c_str() );

      cout.precision(16);
      cout << param <<" \t: " << M << endl;


      File >> readstr; 
      
      param = readstr.substr(0,readstr.find("="));

      R=atof( readstr.substr(readstr.find("=")+1).c_str() );

      cout.precision(16);
      cout << param <<" \t: " << R << endl;

      File >> readstr; 


      string dummystr;

      getline(File, dummystr); 

      //cout << dummystr << endl;

      valarray <double> rs(dim), lap(dim), Dlap(dim), psi4(dim), Dpsi4(dim);

      valarray <double> rho_(dim), pres(dim), epsl_(dim);

      File.setf (ios_base::scientific, ios_base::floatfield);


      for(int ii=0; ii < dim; ii++)   
	{
      
	  File >> rs[ii]
	       >> lap[ii]
	       >> Dlap[ii]
	       >> psi4[ii]
	       >> Dpsi4[ii]
	       >> rho_[ii]
	       >> pres[ii]
	       >> epsl_[ii];
	   
	}


      /*
      for(int ii=0; ii < dim; ii++)   

	cout << rs[ii]    << "\t"
	     << lap[ii]   << "\t"
	     << Dlap[ii]  << "\t"
	     << psi4[ii]  << "\t"
	     << Dpsi4[ii] << "\t"
	     << rho_[ii]  << "\t"
	     << pres[ii]  << "\t"
	     << epsl_[ii] << endl;
      */


      chebycoe_Extremes( lap,   dim, 0 );  
      chebycoe_Extremes( Dlap,  dim, 0 );
      chebycoe_Extremes( psi4,  dim, 0 );
      chebycoe_Extremes( Dpsi4, dim, 0 );
      chebycoe_Extremes( rho_,  dim, 0 );
      chebycoe_Extremes( pres,  dim, 0 );
      chebycoe_Extremes( epsl_, dim, 0 );
    
      /*
      
      for(int ii=0; ii < dim; ii++)   

	cout << rs[ii]    << "\t"
	     << lap[ii]   << "\t"
	     << Dlap[ii]  << "\t"
	     << psi4[ii]  << "\t"
	     << Dpsi4[ii] << "\t"
	     << rho_[ii]  << "\t"
	     << pres[ii]  << "\t"
	     << epsl_[ii] << endl;
      
      */


      // to test the symmetries => one can test a non staggered grid :)
      double Dx=0.;
      double Dy=0.;
      double Dz=0.;
      //if (Getv("stagger_TOV", "yes")) Dx=Dy=Dz=-0.5*level->dx;

      


      vars.Set_Iterator (coeff_names, ALL);

      double twopi=8.0*atan(1.0);

      do
	{

	  double x = vars.Get_x ("psi0");

	  double y = vars.Get_y ("psi0");

	  double z = vars.Get_z ("psi0");

	  double r = sqrt((x+Dx)*(x+Dx) + (y+Dy)*(y+Dy) + (z+Dz)*(z+Dz));
	  
	  if (r<R) {


	    vars.Set_val (pow(chebyeva( 0.0, 1.0, psi4, dim, r/R ),.25), "psi0");

	    vars.Set_val (chebyeva(0.0, 1.0, lap, dim, r/R ), "alpha");

	    vars.Set_val (chebyeva( 0.0, 1.0, rho_, dim, r/R ), "matter_rho");

	    vars.Set_val (chebyeva( 0.0, 1.0, epsl_, dim, r/R ), "matter_epsl");

	  }
	  else{

	    vars.Set_val (1.0+0.5*M/r, "psi0");

	    vars.Set_val ((1.0-0.5*M/r)/(1.0+0.5*M/r), "alpha");

	    vars.Set_val (0.0, "matter_rho");

	    vars.Set_val (0.0, "matter_epsl");

	  }


	  vars.Set_val (twopi*pow(vars.Get_val("psi0"),8)*vars.Get_val("matter_rho")*(1.0+vars.Get_val("matter_epsl")), "rhoADM");

	  // set psi --> 1
	  //vars.Set_val (twopi*vars.Get_val("matter_rho")*(1.0+vars.Get_val("matter_epsl")), "rhoADM");

	}
      while (vars.End (coeff_names));



      File.close();

    }
  else 

    error_exit("Failed to open file",tov_filename,"TOV equilibrium solution.");



  return (0);
}




/*

double elliptic::plegendre(int l, int m, double theta)
{
  //see numerical recipes
  double x = cos(theta);
  int i,ll;
  double fact,oldfact,pll,pmm,pmmp1,omx2;
    
  pmm = 1.;
  if (m) {
    omx2 = (1.-x)*(1.+x);
    fact = 1.;
    for (i=1; i<=m; i++) {
      pmm *= omx2*fact/(fact+1.);
      fact += 2.;
    }
  }
  pmm = sqrt((2*m+1)*pmm/(4.*PI));
  if (m && 1) pmm = -pmm;
  if (l==m) {
    return pmm;
  } else {
    pmmp1 = x*sqrt(2.*m+3.)*pmm;
    if (l==(m+1)) {
      return pmmp1;
    } else {
      oldfact = sqrt(2.*m+3.);
      for (ll=m+2; ll<=l; ll++) {
	fact = sqrt((4.*ll*ll-1.)/(ll*ll-m*m));
	pll = (x*pmmp1-pmm/oldfact)*fact;
	oldfact = fact;
	pmm = pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }
}
*/


 
double elliptic::plegendre(int l, int m, double x)
{

  double fact,pll,pmm,pmmp1,somx2;
  
  int i,ll;
  

  pmm=1.0; 

  if (m > 0) {

    somx2=sqrt((1.0-x)*(1.0+x));

    fact=1.0;


    for (i=1;i<=m;i++) {

      pmm *= -fact*somx2;

      fact += 2.0;
    }
  }
  
  if (l == m)
  
    return pmm;

  else { 
          
      pmmp1=x*(2*m+1)*pmm;

      if (l == (m+1))

	return pmmp1;

      else { 
    
	for (ll=m+2;ll<=l;ll++) {
	  pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	  pmm=pmmp1;
	  pmmp1=pll;
	}
	return pll;
      }
  }
}
 
 
//fast factorial for small numbers 

double elliptic::ffact(double n){


  double f[] = {1,         1,          2,     
		6,         24,         120,
		720,       5040,       40320, 
		362880,    3628800,    39916800,   
		479001600, 6227020800, 87178291200};

  if(n<=14)
    
    return(f[(int)n]);

  else
    
    return(fact(n));



}


double elliptic::fact(double n)
{

  
  if(n < 0)
    {
      error_exit("Computing a negative factorial");
      
    }
  if (n<=14)

      return(ffact(n));
    
    else
    {

        n = n*fact(n-1);

        return(n);
 
    }

}


void elliptic::computeSphericalHarmonic(double theta, double phi, double *Yr, double *Yi, int l, int m)
{
    if ((l<m) || (l<0) || (m<0)) error_exit("wrong l or m inside SphericalHarmonic");
    
    double Plm = sqrt((2.0*l+1)*ffact(l-m)/(4*PI*ffact(l+m))) * plegendre(l,m, cos(theta));
    
    *Yr = Plm * cos(((double)m)*phi);
    *Yi = Plm * sin(((double)m)*phi);


}


/*
void elliptic::add_perturbation()
{


  ofstream File ("SH.xg", fstream::out);

  File.setf (ios_base::scientific, ios_base::floatfield);

  File.precision (12);

  
  for(int l=0; l <= 3; l++)
    
    for(int m=-l; m <= l; m++)
      
      cout << l << "\t"
	   << m << "\t"
	   << sqrt((2.0*l+1)*ffact(l-m)/(4*PI*ffact(l+m))) << endl;
  

  
  double Yr, Yi;

  for(double theta=0; theta < PI; theta += 0.1) 
    
    for(double phi=0; phi < 2*PI; phi += 0.1) 
      {

   
	computeSphericalHarmonic(theta,phi,&Yr,&Yi,2,0);

	File << theta << "\t"
	     << phi << "\t"
	     << Yr << "\t"
	     << Yi <<  endl;
	//     << plegendre(2,0,cos(theta)) << endl;

      }


}
*/


void elliptic::add_perturbation()
{





}


//===========================================Trumpet==============================


double
elliptic::Trumpet_source (double x, double y, double z)
{

  double r2 = x * x + y * y + z * z;

/*
    if(r2 != 0.0)
        
        return( -(6.0+r2+2*exp(r2)*(2*r2-3) )* exp(-r2) / (exp(r2)-1) );

    else

        return( 1.0 );
    
    return( -( (6.0+r2)*exp(-r2)-(4*r2-6) ) * exp(-r2) / (exp(-r2)+1) );
*/

  if (r2 < 2 * numeric_limits < double >::epsilon ())

    return (0.0);

  else

    return (-5.0 * r2 / (exp (r2) - 1) + exp (-r2) * (r2 + 6));



}



double
elliptic::Trumpet (double x, double y, double z)
{

  double r2 = x * x + y * y + z * z;


  return (1 - exp (-r2));



}






double
a2R (double a, double M)
{


  double ac = sqrt (10.0) - 3.0;

  double Rc = M / (4.0 * ac);

  double Cexpa = 0.5 * (Rc * Rc * Rc) * exp (a - ac) * pow (M, 4);

  double S =
    a >
    ac ? (1 - a) / (Rc * (1 - ac)) : (1.0 / Rc - 1.0 / 1.312) * a / ac +
    1.0 / 1.312;

  double F,
    dF;

  double a2m1 = a * a - 1.0;



  do
    {
      F = a2m1 + 2.0 * S - Cexpa * pow (S, 4);

      dF = 2 - 4 * Cexpa * pow (S, 3);

      S -= F / dF;

    }
  while (abs (F) > 1e-8);


  return (1.0 / S);



}





double
r2a (double r, double M)
{

  double rc = 0.30345204271479997 * M;

  double a = sqrt (10.0) - 3.0;

  double ds = abs (rc - r) / 1000;

  if (r > rc)

    for (double s = rc; s <= r; s += ds)
      {
	a += rhs_arM (a, r, M) * ds;
      }

  else

    for (double s = rc; s >= r; s += ds)
      {
	a += rhs_arM (a, r, M) * ds;
      }


  return (a);

}


double
rhs_arM (double a, double r, double M)
{
  double R = a2R (a, M);

  return (a * (6 * M + 4 * R * (a * a - 1)) /
	  (r * (2 * M + R * (a * a - 2 * a - 1))));

}



//=============================== EM_BH ==============================//


double elliptic::EM_BH_phi (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      cval += bhqp[i] / r;

    }


  return (cval);


}



double elliptic::EM_BH_phi_dx (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r3 = pow (x * x + y * y + z * z,1.5);

      if (r3 < 2 * numeric_limits < double >::epsilon ())
	r3 = 2 * numeric_limits < double >::epsilon ();

      cval -= bhqp[i]*x / r3;

    }


  return (cval);


}


double elliptic::EM_BH_phi_dy (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r3 = pow (x * x + y * y + z * z,1.5);

      if (r3 < 2 * numeric_limits < double >::epsilon ())
	r3 = 2 * numeric_limits < double >::epsilon ();

      cval -= bhqp[i]*y / r3;

    }


  return (cval);


}

double elliptic::EM_BH_phi_dz (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r3 = pow (x * x + y * y + z * z,1.5);

      if (r3 < 2 * numeric_limits < double >::epsilon ())
	r3 = 2 * numeric_limits < double >::epsilon ();

      cval -= bhqp[i]*z / r3;

    }


  return (cval);


}



double elliptic::EM_BH_eta (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      cval += bhmp[i] / r;

    }


  return (0.5*cval);


}




double elliptic::EM_BH_eta_dx (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r3 = pow (x * x + y * y + z * z,1.5);

      if (r3 < 2 * numeric_limits < double >::epsilon ())
	r3 = 2 * numeric_limits < double >::epsilon ();

      cval -= bhmp[i]*x / r3;

    }


  return (0.5*cval);


}


double elliptic::EM_BH_eta_dy (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r3 = pow (x * x + y * y + z * z,1.5);

      if (r3 < 2 * numeric_limits < double >::epsilon ())
	r3 = 2 * numeric_limits < double >::epsilon ();

      cval -= bhmp[i]*y / r3;

    }


  return (0.5*cval);


}

double elliptic::EM_BH_eta_dz (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r3 = pow (x * x + y * y + z * z,1.5);

      if (r3 < 2 * numeric_limits < double >::epsilon ())
	r3 = 2 * numeric_limits < double >::epsilon ();

      cval -= bhmp[i]*z / r3;

    }


  return (0.5*cval);


}



//============================Bowen-York =============================//


double elliptic::Bowen_York_rho0(double X,double Y,double Z, double dx, double dy, double dz){

  double cval = 0;

  double aa = 0.01;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);
      double dr = sqrt (dx * dx + dy * dy + dz * dz);

      if (r <= dr)
	cval +=  3*bhpx[i]/(4*M_PI*dr*dr*dr);


    }


  return (cval);

}



double elliptic::Bowen_York_rho1(double X,double Y,double Z, double dx, double dy, double dz){

  double cval = 0;

  double aa = 0.01;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);
      double dr = sqrt (dx * dx + dy * dy + dz * dz);

      if (r <= dr)
	cval +=  3*bhpy[i]/(4*M_PI*dr*dr*dr);

    }


  return (cval);

}



double elliptic::Bowen_York_rho2(double X,double Y,double Z, double dx, double dy, double dz){

  double cval = 0;

  double aa = 0.01;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);
      double dr = sqrt (dx * dx + dy * dy + dz * dz);

      if (r <= dr)
      	cval +=  3*bhpz[i]/(4*M_PI*dr*dr*dr);

    }


  return (cval);

}




double elliptic::Bowen_York_rho0_test(double X,double Y,double Z){

  double r2 = X*X+Y*Y+Z*Z;

  //return 4*exp(-(BYT1+BYT2+BYT3)*r2)*(-5*BYT1*exp((BYT2+BYT3)*r2)+X*(  BYT2*BYT2*exp((BYT1+BYT3)*r2)*Y+BYT3*BYT3*exp((BYT1+BYT2)*r2)*Z )+BYT1*BYT1*exp((BYT2+BYT3)*r2)*(X*X+3*r2) )/3;


  //  return 2*BYT1*exp(-BYT1*r2)*(2*BYT1*r2-3) - 2*(BYT1*exp(-BYT1*r2)*X + BYT2*exp(-BYT2*r2)*Y + BYT3*exp(-BYT3*r2)*Z)/3;

  return 2*BYT1*exp(-BYT1*r2)*(2*BYT1*r2-3);

}


double elliptic::Bowen_York_rho1_test(double X,double Y,double Z){

  double r2 = X*X+Y*Y+Z*Z;

  //return 4*exp(-(BYT1+BYT2+BYT3)*r2)*(-5*BYT2*exp((BYT1+BYT3)*r2)+Y*(  BYT1*BYT1*exp((BYT2+BYT3)*r2)*X+BYT3*BYT3*exp((BYT1+BYT2)*r2)*Z )+BYT2*BYT2*exp((BYT1+BYT3)*r2)*(Y*Y+3*r2) )/3;

  //return 2*BYT2*exp(-BYT2*r2)*(2*BYT2*r2-3) - 2*(BYT1*exp(-BYT1*r2)*X + BYT2*exp(-BYT2*r2)*Y + BYT3*exp(-BYT3*r2)*Z)/3;

  return 2*BYT2*exp(-BYT2*r2)*(2*BYT2*r2-3);

}



double elliptic::Bowen_York_rho2_test(double X,double Y,double Z){

  double r2 = X*X+Y*Y+Z*Z;
  
  //return 4*exp(-(BYT1+BYT2+BYT3)*r2)*(-5*BYT3*exp((BYT1+BYT2)*r2)+Z*(  BYT1*BYT1*exp((BYT2+BYT3)*r2)*X+BYT2*BYT2*exp((BYT1+BYT3)*r2)*Y )+BYT3*BYT3*exp((BYT1+BYT2)*r2)*(Z*Z+3*r2) )/3;

  //return 2*BYT3*exp(-BYT3*r2)*(2*BYT3*r2-3)- 2*(BYT1*exp(-BYT1*r2)*X + BYT2*exp(-BYT2*r2)*Y + BYT3*exp(-BYT3*r2)*Z)/3;

    return 2*BYT3*exp(-BYT3*r2)*(2*BYT3*r2-3);

}







//=============================== BY_EM ==============================//


double elliptic::BY_EM_phi (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt (x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();

      cval += bhqp[i] / (r - (bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i]  );

    }


  return (cval);


}



double elliptic::BY_EM_phi_dx (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt(x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();
      if(bhmp[i]!=0)
	cval -= bhqp[i]*(x-bhpx[i]*r/bhmp[i]) / (r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2));

    }


  return (cval);


}




double elliptic::BY_EM_phi_dy (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt(x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();
      if(bhmp[i]!=0)
	cval -= bhqp[i]*(y-bhpy[i]*r/bhmp[i]) / (r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2));

    }


  return (cval);


}

double elliptic::BY_EM_phi_dz (double X, double Y, double Z)
{


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt(x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();
      if(bhmp[i]!=0)
	cval -= bhqp[i]*(z-bhpz[i]*r/bhmp[i]) / (r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2));

    }


  return (cval);


}



double elliptic::BY_EM_Bfield_x(double X,double Y,double Z){


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt(x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();
      if(bhmp[i]!=0)
	cval += y*bhqp[i]*(z-bhpz[i]*r/bhmp[i]) / (r*r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2)) - z*bhqp[i]*(y-bhpy[i]*r/bhmp[i]) / (r*r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2));

    }


  return (cval);




}





double elliptic::BY_EM_Bfield_y(double X,double Y,double Z){


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt(x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();
      if(bhmp[i]!=0)
	cval += z*bhqp[i]*(x-bhpx[i]*r/bhmp[i]) / (r*r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2)) - x*bhqp[i]*(z-bhpz[i]*r/bhmp[i]) / (r*r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2));

    }


  return (cval);




}



double elliptic::BY_EM_Bfield_z(double X,double Y,double Z){


  double cval = 0;

  for (size_t i = 0; i < NumBH; i++)
    {


      double x = X - bhx[i];

      double y = Y - bhy[i];

      double z = Z - bhz[i];


      double r = sqrt(x * x + y * y + z * z);

      if (r < 2 * numeric_limits < double >::epsilon ())
	r = 2 * numeric_limits < double >::epsilon ();
      if(bhmp[i]!=0)
	cval += x*bhqp[i]*(y-bhpy[i]*r/bhmp[i]) / (r*r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2)) - y*bhqp[i]*(x-bhpx[i]*r/bhmp[i]) / (r*r*pow(r-(bhpx[i]*x+bhpy[i]*y+bhpz[i]*z)/bhmp[i],2));

    }


  return (cval);




}
