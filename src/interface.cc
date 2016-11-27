//
// programa.cc
//  
// Made by Pablo Galaviz
// Login   <pablo@NerV>
// 
// Started on  Thu Nov 26 17:34:26 2009 Pablo Galaviz
// Started on  Thu Nov 26 17:34:26 2009 Pablo Galaviz
//



//  This file is part of Olliptic
//
//  Olliptic is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  any later version.
//
//  Olliptic is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with Olliptic.  If not, see <http://www.gnu.org/licenses/>.
//



#include "interface.h"



interface::interface (int argc, char *argv[])
{




  //== join parameters ==//

  stringstream param_tmp;

  for (int i = 0; i < argc; i++)

    param_tmp << argv[i] << " ";

  param = param_tmp.str ();



  //== Set default parameters ==//



  by_domain = false;

  by_points = false;

  by_grid_size = false;







  //== Add double parameters ==//



  add_double_param ("dxyz", 0.5,
		    "Grid size [dx=dy=dz]. (default : 0.5)",
		    "-d", "--grid-size");

  add_double_param ("dx", 0.5,
		    "Grid size in x direction. (default : 0.5)",
		    "-dx", "--grid-size-x");

  add_double_param ("dy", 0.5,
		    "Grid size in y direction. (default : 0.5)",
		    "-dy", "--grid-size-y");

  add_double_param ("dz", 0.5,
		    "Grid size in z direction. (default : 0.5)",
		    "-dz", "--grid-size-z");


  add_double_param ("L", 2.0,
		    "Domain length [lx=ly=lz]. (default : 2.0)",
		    "-L", "--length");

  add_double_param ("Lx", 2.0,
		    "Domain length in x direction. (default : 2.0)",
		    "-Lx", "--length-x");


  add_double_param ("Ly", 2.0,
		    "Domain length in y direction. (default : 2.0)",
		    "-Ly", "--length-y");

  add_double_param ("Lz", 2.0,
		    "Domain length in z direction. (default : 2.0)",
		    "-Lz", "--length-z");


  add_double_param ("Dprint", -1.0,
		    "Grid size for interpolated data. (default : L/1000)",
		    "-Dp", "--d-print");


  add_double_param ("tolerance", 1e-8,
		    "Residual tolerance || r^h ||_inf < s (default : 1.0e-8)",
		    "-s", "--tolerance");


  add_double_param ("ps_r0", 120,
		    "Radius of the scalar field perturbation.",
		    "-r0", "--radius_sf_r0");

  add_double_param ("ps_sigma", 8,
		    "Sigma of the scalar field perturbation.",
		    "-sigma", "--sigma_sf");

  add_double_param ("ps_phi0", 0.025,
		    "Amplitude of the scalar field perturbation.",
		    "-phi0", "--phi0_sf");


  add_double_param ("ps_a2", 1,
		    "Scalar field perturbation parameter a2. [1]",
		    "-a2", "--a2_sf");



  //== Add logic parameters ==//

  add_logic_param ("io_bam", false,
		   "Equal to yes if Olliptic should produce ID for bam (default : no).",
		   "-iob", "--io-bam");


  add_logic_param ("io_Zcode", false,
		   "Equal to yes if Olliptic should produce ID for Zcode (default : no).",
		   "-ioZ", "--io-Zcode");

  add_logic_param ("print_coeff", false,
		   "Print coefficients (default : no).",
		   "-pc", "--print-coeff");

  add_logic_param ("print_1d", true,
		   "Print 1d output (default : yes).", "-p1d", "--print-1d");

  add_logic_param ("print_2d", false,
		   "Print 2d output (default : no).", "-p2d", "--print-2d");

  add_logic_param ("print_3d", false,
		   "Print 3d output (default : no).", "-p3d", "--print-3d");

  add_logic_param ("print_1d_it", false,
		   "Print 1d output for each iteration (default : no).",
		   "-p1d_it", "--print-1d-iteration");

  add_logic_param ("print_interpol", true,
		   "Print 1d output with fix mesh (given by -dtp) (default : yes).",
		   "-pint", "--print-interpol");

  add_logic_param ("interpol_output", true,
		   "Interpolate 1d and 2d output (default : yes).",
		   "-int", "--interpol-output");


  add_logic_param ("fit_approx", false,
		   "Fit approximate initial data (default : no).",
		   "-fapp", "--fit-appox");


  //== Add integer parameters =//


  add_int_param ("nxyz", 30,
		 "Number of grid points [nx=ny=nz]. (default : 30)",
		 "-n", "-N");


  add_int_param ("nx", 30,
		 "Number of grid points in x direction. (default : 30)",
		 "-nx", "--grid-points-x");

  add_int_param ("ny", 30,
		 "Number of grid points in y direction. (default : 30)",
		 "-ny", "--grid-points-y");

  add_int_param ("nz", 30,
		 "Number of grid points in z direction. (default : 30)",
		 "-nz", "--grid-points-z");


  add_int_param ("levels", 1,
		 "Number of refinement levels. (default : 1 )",
		 "-l", "--levels");



  add_int_param ("order", 4,
		 "Order of stencil [2,4,6,8,10]. (default : 4)",
		 "-O", "--order");




  add_int_param ("nu1", 3,
		 "nu1 pre-smooth iterations. (default : 3)",
		 "-nu1", "--pre-smooth");




  add_int_param ("nu2", 3,
		 "nu2 post-smooth iterations. (default : 3)",
		 "-nu2", "--post-smooth");



  add_int_param ("NumBH", 1,
		 "Number of punctures. (default : 1)", "-NumBH", "--num-bh");


  add_int_param ("ps_case", 1,
		 "Puncture + scalar field case [1,2,3]. (default : 1)",
		 "-psc", "-ps_case");


  //== Add string parameters =//




  add_string_param ("dir_name", "Output",
		    "Output dirctory. (default : Output)", "-o", "--out-dir");

  add_string_param ("file_bam", "none",
		    "Bam initial data file (ex. punctures_u.xyz).",
		    "-iof", "--io-bam-file");

  add_string_param ("file_TOV", "none",
		    "Name of the equilibrium TOV ID.",
		    "-tovf", "--tov-file");


  add_string_param ("help", "h",
		    "Displays this usage screen.", "-h", "--help");


  add_string_param ("equation", "none",
		    string_equation (), "-eq", "--equation");

  add_string_param ("method", "Vcycle", string_method (), "-m", "--method");



  add_string_param ("output", "yes",
		    "Standar output: yes, partial or no. [default : yes]",
		    "-so", "--std-output");


  add_string_param ("domain_sym", "full",
		    string_domain (), "-ds", "--domain-sym");


  add_string_param ("boundary", "Dirichlet",
		    string_boundary (), "-B", "--boundary");



  equation = POISSON;



  boundary = DIRICHLET;





  bhmp.assign (1, 0);

  bhqp.assign (1, 0);

  bhx.assign (1, 0);

  bhy.assign (1, 0);

  bhz.assign (1, 0);


  bhpx.assign (1, 0);

  bhpy.assign (1, 0);

  bhpz.assign (1, 0);


  bhsx.assign (1, 0);

  bhsy.assign (1, 0);

  bhsz.assign (1, 0);









  exec = true;

  def_eq = false;



  //== read the user command line ==//

  read_console (argc, argv);


  //== check parameters ==//

  proc_parameters ();


  //== Print welcome ==//

  if (!(output == NO))

    print_welcome ();


  //== if something is wrong stop ==//

  if (exec == false || def_eq == false)
    {

      print_help ();

      exit (0);

    }



}









void
interface::read_console (int argc, char *argv[])
{



  //== Read and fill the parameters 

  if (argc > 2)

    for (int i = 1; i < argc; i++)
      {

	char *finalPtr;





//=============== Reading double parameters ====================



	for (size_t iv = 0; iv < double_par_names.size (); iv++)
	  {

	    if (strcmp
		(argv[i], double_option1[double_par_names[iv]].c_str ()) == 0
		|| strcmp (argv[i],
			   double_option2[double_par_names[iv]].c_str ()) ==
		0)
	      {
		if (i + 1 < argc)
		  {


		    double val = strtod (argv[i + 1], &finalPtr);

		    //if (val >= 0){

			double_par[double_par_names[iv]] = val;


			if (double_par_names[iv] == "L" ||
			    double_par_names[iv] == "Lx" ||
			    double_par_names[iv] == "Ly" ||
			    double_par_names[iv] == "Lz")
			  {

			    if (double_par_names[iv] == "L")
			      {
				double_par["Lx"] = val;

				double_par["Ly"] = val;

				double_par["Lz"] = val;


			      }


			    by_domain = true;

			  }



			if (double_par_names[iv] == "dxyz" ||
			    double_par_names[iv] == "dx" ||
			    double_par_names[iv] == "dy" ||
			    double_par_names[iv] == "dz")
			  {

			    if (double_par_names[iv] == "dxyz")
			      {
				double_par["dx"] = val;

				double_par["dy"] = val;

				double_par["dz"] = val;


			      }


			    by_grid_size = true;

			  }


			/*

		      }
		    else
		      {

			cout << endl << endl <<
			  " Error -- Options for float number parameter" <<
			  double_option1[double_par_names[iv]] <<
			  " must be: [ positive double precision float] " <<
			  endl;

			exec = false;

		      }
			*/

		  }

		else
		  {

		    cout << endl << endl <<
		      " Error -- Missing float value for: " <<
		      double_option1[double_par_names[iv]] << endl;

		    exec = false;

		  }

	      }


	  }

















//=============== Reading integer parameters ====================



	for (size_t iv = 0; iv < int_par_names.size (); iv++)
	  {

	    if (strcmp (argv[i], int_option1[int_par_names[iv]].c_str ()) == 0
		|| strcmp (argv[i],
			   int_option2[int_par_names[iv]].c_str ()) == 0)
	      {
		if (i + 1 < argc)
		  {


		    int val = strtoul (argv[i + 1], &finalPtr, 10);

		    if (val >= 0)
		      {

			int_par[int_par_names[iv]] = size_t (val);


			if (int_par_names[iv] == "nxyz" ||
			    int_par_names[iv] == "nx" ||
			    int_par_names[iv] == "ny" ||
			    int_par_names[iv] == "nz")
			  {

			    if (int_par_names[iv] == "nxyz")
			      {
				int_par["nx"] = size_t (val);;

				int_par["ny"] = size_t (val);;

				int_par["nz"] = size_t (val);;


			      }


			    by_points = true;
			  }



		      }
		    else
		      {

			cout << endl << endl <<
			  " Error -- Options for integer parameter" <<
			  int_option1[int_par_names[iv]] <<
			  " must be: [ positive integer ] " << endl;

			exec = false;

		      }


		  }

		else
		  {

		    cout << endl << endl <<
		      " Error -- Missing integer value for: " <<
		      int_option1[int_par_names[iv]] << endl;

		    exec = false;

		  }

	      }


	  }



//=============== Reading logic parameters ====================


	for (size_t iv = 0; iv < logic_par_names.size (); iv++)
	  {

	    if (strcmp (argv[i], logic_option1[logic_par_names[iv]].c_str ())
		== 0
		|| strcmp (argv[i],
			   logic_option2[logic_par_names[iv]].c_str ()) == 0)
	      {
		if (i + 1 < argc)
		  {
		    if (strcmp (argv[i + 1], "yes") == 0 ||
			strcmp (argv[i + 1], "YES") == 0 ||
			strcmp (argv[i + 1], "true") == 0 ||
			strcmp (argv[i + 1], "no") == 0 ||
			strcmp (argv[i + 1], "NO") == 0 ||
			strcmp (argv[i + 1], "false") == 0)
		      {


			if (strcmp (argv[i + 1], "yes") == 0 ||
			    strcmp (argv[i + 1], "YES") == 0 ||
			    strcmp (argv[i + 1], "true") == 0)

			  logic_par[logic_par_names[iv]] = true;

			else

			  logic_par[logic_par_names[iv]] = false;


		      }
		    else
		      {

			cout << endl << endl <<
			  " Error -- Options for logic parameter: " <<
			  logic_option1[logic_par_names[iv]] <<
			  " must be: [yes, YES, true, no, NO, false] " <<
			  endl;

			exec = false;

		      }


		  }

		else
		  {

		    cout << endl << endl <<
		      " Error -- Missing logic value for: " <<
		      logic_option1[logic_par_names[iv]] << endl;

		    exec = false;

		  }

	      }


	  }






//=============== Reading string parameters ====================



	for (size_t iv = 0; iv < string_par_names.size (); iv++)
	  {

	    if (strcmp
		(argv[i], string_option1[string_par_names[iv]].c_str ()) == 0
		|| strcmp (argv[i],
			   string_option2[string_par_names[iv]].c_str ()) ==
		0)
	      {
		if (string_par_names[iv] != "help")
		  {
		    if (i + 1 < argc)
		      {

			string_par[string_par_names[iv]] = argv[i + 1];

		      }

		    else
		      {

			cout << endl << endl <<
			  " Error -- Missing text value for: " <<
			  string_option1[string_par_names[iv]] << endl;

			exec = false;

		      }
		  }

	      }


	  }





	if (strcmp (argv[i], "--help") == 0 || strcmp (argv[i], "-h") == 0)

	  exec = false;











      }


//== for punctures we must fill vectors by default set to zero ==//

  size_t NumBH = Get_integer_par ("NumBH");

  bhmp.assign (NumBH, 0);

  bhqp.assign (NumBH, 0);

  bhx.assign (NumBH, 0);

  bhy.assign (NumBH, 0);

  bhz.assign (NumBH, 0);


  bhpx.assign (NumBH, 0);

  bhpy.assign (NumBH, 0);

  bhpz.assign (NumBH, 0);


  bhsx.assign (NumBH, 0);

  bhsy.assign (NumBH, 0);

  bhsz.assign (NumBH, 0);


//== now we will read the "true" values ==// 

  for (int i = 1; i < argc; i++)
    {

      char *finalPtr;




      for (size_t j = 0; j < NumBH; j++)
	{


	  stringstream num;

	  num << j + 1;

	  string bhparam;

	  bhparam = "-bhmp" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhmp[j] = strtod (argv[i + 1], &finalPtr);


	  bhparam = "-bhqp" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhqp[j] = strtod (argv[i + 1], &finalPtr);




	  bhparam = "-bhx" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhx[j] = strtod (argv[i + 1], &finalPtr);



	  bhparam = "-bhy" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhy[j] = strtod (argv[i + 1], &finalPtr);




	  bhparam = "-bhz" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhz[j] = strtod (argv[i + 1], &finalPtr);





	  bhparam = "-bhpx" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhpx[j] = strtod (argv[i + 1], &finalPtr);



	  bhparam = "-bhpy" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhpy[j] = strtod (argv[i + 1], &finalPtr);




	  bhparam = "-bhpz" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhpz[j] = strtod (argv[i + 1], &finalPtr);




	  bhparam = "-bhsx" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhsx[j] = strtod (argv[i + 1], &finalPtr);



	  bhparam = "-bhsy" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhsy[j] = strtod (argv[i + 1], &finalPtr);




	  bhparam = "-bhsz" + num.str ();

	  if (strcmp (argv[i], bhparam.c_str ()) == 0)

	    bhsz[j] = strtod (argv[i + 1], &finalPtr);







	}

    }



//=========================== End read console =============================//



}







void
Valid_N (size_t & N)
{


  size_t n_V[63] = { 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 17, 19, 21,
    23, 25, 27, 29, 33, 37, 41, 45,
    49, 53, 57, 65, 73, 81, 89, 97,
    105, 113, 129, 145, 161, 177, 193, 209,
    225, 257, 289, 321, 353, 385, 417, 449,
    513, 577, 641, 705, 769, 833, 897, 1153,
    1281, 1409, 1665, 1793, 2305, 2817, 3329
  };				//valid grid size


  bool valid = false;


  size_t suggestion_i = n_V[61],
    suggestion_f = n_V[62],
    i;				// Store valid size sugention 



  for (int j = 0; j < 63; j++)

    if (n_V[j] == N)

      valid = true;



  if (!valid)
    {

      // Error in the axis: search for suggestion sizes

      i = 0;			//Array index

      while (n_V[i] < N && i < 63)
	{

	  suggestion_i = n_V[i];	// lower size

	  suggestion_f = n_V[i + 1];	// upper size

	  i++;


	}



      if (N - suggestion_i < suggestion_f - N)

	N = suggestion_i;

      else

	N = suggestion_f;

    }




}







void
interface::proc_parameters ()
{






  if (Get_integer_par ("levels") < 1)
    {

      cout << "\n\n*** Error: levels must be >= 1" << endl;

      exec = false;



    }



  size_t order = Get_integer_par ("order");


  if (!(order == 2 || order == 4 || order == 6 || order == 8))
    {

      cout << "\n\n*** Error: order must be [2,4,6,8]." << endl;

      exec = false;



    }




  if (by_domain)
    {


      if (by_points && by_grid_size)
	{

	  cout << "\n\n*** Error: Domain definition" << endl;

	  exec = false;

	}
      else

       if (by_grid_size)
	{

	  double Lx = double_par["Lx"];

	  double Ly = double_par["Ly"];

	  double Lz = double_par["Lz"];


	  double dx = double_par["dx"];

	  double dy = double_par["dy"];

	  double dz = double_par["dz"];


	  int_par["nx"] = size_t (fabs (Lx) / dx + 1);

	  int_par["ny"] = size_t (fabs (Ly) / dy + 1);

	  int_par["nz"] = size_t (fabs (Lz) / dz + 1);


	}

      else
	{


	  double Lx = double_par["Lx"];

	  double Ly = double_par["Ly"];

	  double Lz = double_par["Lz"];


	  size_t nx = int_par["nx"];

	  size_t ny = int_par["ny"];

	  size_t nz = int_par["nz"];


	  Valid_N (nx);

	  Valid_N (ny);

	  Valid_N (nz);


	  double_par["dx"] = fabs (Lx) / double (nx - 1.0);

	  double_par["dy"] = fabs (Ly) / double (ny - 1.0);

	  double_par["dz"] = fabs (Lz) / double (nz - 1.0);


	  int_par["nx"] = nx;

	  int_par["ny"] = ny;

	  int_par["nz"] = nz;


	}






    }
  else
    {



      size_t nx = int_par["nx"];

      size_t ny = int_par["ny"];

      size_t nz = int_par["nz"];


      double dx = double_par["dx"];

      double dy = double_par["dy"];

      double dz = double_par["dz"];



      Valid_N (nx);

      Valid_N (ny);

      Valid_N (nz);


      double_par["Lx"] = dx * (nx - 1.0);

      double_par["Ly"] = dy * (ny - 1.0);

      double_par["Lz"] = dz * (nz - 1.0);


      int_par["nx"] = nx;

      int_par["ny"] = ny;

      int_par["nz"] = nz;



    }


  bool mass_inv = false;

  size_t NumBH = Get_integer_par ("NumBH");

  for (size_t i = 0; i < NumBH; i++)

    if (bhmp[i] <= 0)

      mass_inv = true;






  if (mass_inv && equation == PUNCTURES)
    {

      cout << "\n\n*** Error: Wrong punture mass parameters" << endl;

      exec = false;


    }


  if (method == APP && equation != PUNCTURES)
    {

      cout << "\n\n*** Error: App puncture only works with puncture equation"
	<< endl;

      exec = false;


    }




  if (logic_par["io_bam"] && string_par["file_bam"] == "none")
    {


      cout << "\n\n*** Error: Missing file." << endl;

      exec = false;


    }




//== EQUATION PARAMETER ==//

  if (string_par["equation"] == "Poisson")
    {

      equation = POISSON;

      def_eq = true;

    }



  if (string_par["equation"] == "Boundary-test")
    {

      equation = BOUNDARY_TEST;

      def_eq = true;
    }




  if (string_par["equation"] == "Brill-waves")
    {
      equation = BRILL_WAVES;

      def_eq = true;
    }



  if (string_par["equation"] == "Punctures")
    {

      equation = PUNCTURES;

      def_eq = true;
    }

  if (string_par["equation"] == "Punctures_scalar_field")
    {

      equation = PUNCTURES_SCALAR;

      def_eq = true;
    }



  if (string_par["equation"] == "Punctures_EM")
    {

      equation = PUNCTURES_EM;

      def_eq = true;
    }


  if (string_par["equation"] == "Multibox-test")
    {
      equation = MULTIBOX_TEST;


      def_eq = true;

    }


  if (string_par["equation"] == "Generic")
    {
      equation = GENERIC;


      def_eq = true;
    }



  if (string_par["equation"] == "System-test")
    {
      equation = SYSTEM_TEST;

      def_eq = true;

    }



  if (string_par["equation"] == "Trumpet-1+log")
    {
      equation = TRUMPET_1pLOG;

      def_eq = true;

    }


  if (string_par["equation"] == "NSO-ID")
    {
      equation = NSO_ID;

      def_eq = true;

    }


  if (string_par["equation"] == "Metric-Test")
    {
      equation = METRIC_TEST;

      def_eq = true;

    }


  if (string_par["equation"] == "Bowen-York")
    {
      equation = BOWEN_YORK;

      def_eq = true;

    }



  if (string_par["equation"] == "BY-EM")
    {
      equation = BY_EM;

      def_eq = true;

    }


//== END EQUATION PARAMETER ==//





//== METHOD PARAMETER ==//

  method = VCYCLE;


  if (string_par["method"] == "Gauss-Seidel")

    method = GAUSS_SEIDEL;

  if (string_par["method"] == "Full-Multigrid")

    method = FMG;

  if (string_par["method"] == "Vcycle")

    method = VCYCLE;

  if (string_par["method"] == "Wcycle")

    method = WCYCLE;

  if (string_par["method"] == "AppPuncture")

    method = APP;

//== END METHOD PARAMETER ==//




//==OUTPUT PARAMETER ==//

  output = YES;

  if (string_par["output"] == "yes")

    output = YES;

  if (string_par["output"] == "partial")

    output = PARTIAL;

  if (string_par["output"] == "no")

    output = NO;


//== END OUTPUT PARAMETER ==//





//== SYMMETRY PARAMETER ==//

  domain_sym = FULL;

  if (string_par["domain_sym"] == "full")

    domain_sym = FULL;


  if (string_par["domain_sym"] == "octant")

    domain_sym = OCTANT;

  if (string_par["domain_sym"] == "quadrant-xy")

    domain_sym = QUADRANT_XY;

  if (string_par["domain_sym"] == "quadrant-xz")

    domain_sym = QUADRANT_XZ;

  if (string_par["domain_sym"] == "quadrant-yz")

    domain_sym = QUADRANT_YZ;


  if (string_par["domain_sym"] == "bitant-x")

    domain_sym = BITANT_X;

  if (string_par["domain_sym"] == "bitant-")

    domain_sym = BITANT_Y;

  if (string_par["domain_sym"] == "bitant-z")

    domain_sym = BITANT_Z;


  if (string_par["domain_sym"] == "rotate-x")

    domain_sym = ROTATE_X;


  if (string_par["domain_sym"] == "rotate-y")

    domain_sym = ROTATE_Y;

  if (string_par["domain_sym"] == "rotate-z")

    domain_sym = ROTATE_Z;

//== END SYMMETRY PARAMETER ==//






//==BOUNDARY PARAMETER ==//

  boundary = DIRICHLET;

  if (string_par["boundary"] == "Dirichlet")

    boundary = DIRICHLET;

  if (string_par["boundary"] == "Robin")

    boundary = ROBIN;

  if (string_par["boundary"] == "Asymptotic")

    boundary = ASYMPTOTIC;



//== END BOUNDARY PARAMETER ==//















  string dir = string_par["dir_name"] + "/";

  string file_name = string_par["dir_name"];



#ifdef OLLIN_MPI

  const int my_rank = MPI::COMM_WORLD.Get_rank ();


  if (my_rank == 0)
    {

#endif





      string rmdir = "rm -r " + file_name + "_prev";


      int sys_out = system (rmdir.c_str ());



      string mvdir = "mv " + dir + " " + file_name + "_prev";


      sys_out = system (mvdir.c_str ());


      string mkdir = "mkdir " + dir;


      sys_out = system (mkdir.c_str ());







      string par_run = "echo " + param + " > " + dir + file_name + ".run";

      sys_out = system (par_run.c_str ());

      string make_run = "chmod +x " + dir + file_name + ".run";

      sys_out = system (make_run.c_str ());



#ifdef OLLIN_MPI

    }

#endif


  string_par["dir_name"] = dir;




}










void
interface::add_logic_param (string name_par,
			    bool def,
			    string description,
			    string option_key1, string option_key2)
{



  logic_par_names.push_back (name_par);

  logic_par[name_par] = def;

  logic_help[name_par] = description;

  logic_option1[name_par] = option_key1;

  logic_option2[name_par] = option_key2;




}





void
interface::add_int_param (string name_par,
			  size_t def,
			  string description,
			  string option_key1, string option_key2)
{



  int_par_names.push_back (name_par);

  int_par[name_par] = def;

  int_help[name_par] = description;

  int_option1[name_par] = option_key1;

  int_option2[name_par] = option_key2;




}






void
interface::add_double_param (string name_par,
			     double def,
			     string description,
			     string option_key1, string option_key2)
{



  double_par_names.push_back (name_par);

  double_par[name_par] = def;

  double_help[name_par] = description;

  double_option1[name_par] = option_key1;

  double_option2[name_par] = option_key2;




}






void
interface::add_string_param (string name_par,
			     string def,
			     string description,
			     string option_key1, string option_key2)
{



  string_par_names.push_back (name_par);

  string_par[name_par] = def;

  string_help[name_par] = description;

  string_option1[name_par] = option_key1;

  string_option2[name_par] = option_key2;




}





string
interface::string_equation ()
{



  stringstream text;

  text << "Problem to solve:" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Poisson" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Boundary-test" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Brill-waves" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Punctures" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Punctures_scalar_field" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Punctures_EM" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Multibox-test" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Generic" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "System-test" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Trumpet-1+log" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "NSO-ID" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Metric-Test" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Bowen-York" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "BY-EM" << endl;

  return (text.str ());


}





string
interface::string_method ()
{



  stringstream text;

  text << "Method for solve the elliptic equation (default : Vcycle)." <<
    endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Gauss-Seidel" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Full-Multigrid" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Vcycle" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Wcycle" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "AppPuncture" << endl;





  return (text.str ());


}







string
interface::string_domain ()
{



  stringstream text;

  text << "Domain symmetry (default : full)." << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "full" << "- (xi,xf)X(yi,yf)X(zi,zf)" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "octant" << "- (0, xf)X(0, yf)X(0, zf)" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "quadrant-xy" << "- (0, xf)X(0, yf)X(zi,zf)" << endl;


  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "quadrant-xz" << "- (0, xf)X(yi,yf)X(0, zf)" << endl;


  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "quadrant-yz" << "- (xi,xf)X(0, yf)X(0, zf)" << endl;


  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "bitant-x" << "- (0, xf)X(yi,yf)X(zi,zf)" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "bitant-y" << "- (xi,xf)X(0, yf)X(zi,zf)" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "bitant-z" << "- (xi,xf)X(yi,yf)X(0, zf)" << endl;


  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "rotate-x" << "- rotation around X axis." << endl;


  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "rotate-y" << "- rotation around Y axis." << endl;


  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  ";
  text.setf (ios_base::left, ios_base::adjustfield);
  text.width (14);
  text << "rotate-z" << "- rotation around Z axis." << endl;


  return (text.str ());


}





string
interface::string_boundary ()
{



  stringstream text;

  text << "Boundary condition, (default : Dirichlet)" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Dirichlet   - Fixed points at boundary u = A" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " << "Robin       - Solution u = A + B / r^n" << endl;

  text.setf (ios_base::right, ios_base::adjustfield);
  text.width (46);
  text << "  " <<
    "Asymptotic  - The numerical solution is matching to a function at the bondary u <-- A + B(u) / r^n"
    << endl;




  return (text.str ());


}
