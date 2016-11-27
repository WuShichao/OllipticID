//
// output.cc
//  
// Made by Pablo Galaviz
// Login   <pablo@NerV>
// 
// Started on  Tue Oct 27 12:44:49 2009 Pablo Galaviz
// Started on  Tue Oct 27 12:44:49 2009 Pablo Galaviz
//



#include "elliptic.h"

#include "interface.h"




void error_exit(string error_string1,
                string error_string2,
                string error_string3,
                bool warning )
{


    //== Here we just print a message of error or a warning, if we
    //have an error Olliptic finish.
    
#ifdef OLLIN_MPI


    if (MPI::COMM_WORLD.Get_rank () == 0)
    {

#endif


        cerr << endl
             << endl;
        
        if(warning)
        
            cerr << "\t| Warning " << endl;         
        
        else
        
            cerr << "\t| Error " << endl;

        int s1 = error_string1.size(); 

        int s2 = error_string2.size(); 

        int s3 = error_string3.size(); 

        
        int n = max(max(s1,s2),s3);
                
        cerr << "\t|";
        for(int i = 0; i <= n; i++)
            cerr << "-";

        
        cerr << endl;
        
        if(error_string1 != "") cerr << "\t| " << error_string1 << endl;

        if(error_string2 != "") cerr << "\t| " << error_string2 << endl;

        if(error_string3 != "") cerr << "\t| " << error_string3 << endl;



        cerr << "\t|";
        for(int i = 0; i <= n; i++)
            cerr << "-";
        cerr << endl;
        
        
        cerr << "\t| Aborting execution " << endl
             << endl;


#ifdef OLLIN_MPI

    }
    
#endif

    
    if(!warning)
    {


#ifdef OLLIN_MPI

    
        //== Finish MPI ==//    

        MPI::COMM_WORLD.Barrier ();

        MPI::Finalize ();


    
#endif
        
        exit(1);

    }     
}



void
interface::print_welcome ()
{



#ifdef OLLIN_MPI


    if (MPI::COMM_WORLD.Get_rank () == 0)
    {

#endif



        cout << endl
             <<
            "--------------------------------------------------------------------"
             << endl << endl <<
            "   /\\  __`\\ /\\_ \\  /\\_ \\    __        /\\ \\__  __            "
             << endl <<
            "   \\ \\ \\/\\ \\\\//\\ \\ \\//\\ \\  /\\_\\  _____\\ \\ ,_\\/\\_\\    ___    "
             << endl <<
            "    \\ \\ \\ \\ \\ \\ \\ \\  \\ \\ \\ \\/\\ \\/\\ '__`\\ \\ \\/\\/\\ \\  /'___\\  "
             << endl <<
            "     \\ \\ \\_\\ \\ \\_\\ \\_ \\_\\ \\_\\ \\ \\ \\ \\L\\ \\ \\ \\_\\ \\ \\/\\ \\__/  "
             << endl <<
            "      \\ \\_____\\/\\____\\/\\____\\\\ \\_\\ \\ ,__/\\ \\__\\\\ \\_\\ \\____\\ "
             << endl <<
            "       \\/_____/\\/____/\\/____/ \\/_/\\ \\ \\/  \\/__/ \\/_/\\/____/ "
             << endl <<
            "                                   \\ \\_\\                    " <<
            endl <<
            "                                    \\/_/                    " <<
            endl << endl <<
            "--------------------------------------------------------------------"
             << endl << endl << endl;


        cout << "                   ======  Olliptic 9.10  ======" << endl
             << endl
             << "                       Author: Pablo Galaviz    " << endl
             << endl
             << "                     pablo.galaviz@uni-jena.de  " << endl
             << endl
             << "                   =============================" << endl
             << endl;




#ifdef OLLIN_MPI

    }

#endif


}






void
interface::print_help ()
{



#ifdef OLLIN_MPI

    const int my_rank = MPI::COMM_WORLD.Get_rank ();


    if (my_rank == 0)
    {

#endif




        cout << endl << "Usage: Olliptic [Option] {value} " << endl << endl;

//=======================================================================================================//



        cout << "       Float number parameters.             Valid options:" <<
            endl <<
            "                                            [ Positive double ]" <<
            endl <<
            "       --------------------------------------------------------------------"
             << endl << endl;




        for (size_t i = 0; i < double_par_names.size (); i++)
	{

            cout << "       ";
            cout.width (8);
            cout.setf (ios_base::left, ios_base::adjustfield);
            cout << double_option1[double_par_names[i]];
            cout.width (24);
            cout << double_option2[double_par_names[i]] << "  :  " <<
                double_help[double_par_names[i]] << endl;


	}

        cout << endl << endl;

//===================================== end double parameters help ======================================



        cout << "       Integer parameters.                  Valid options:" <<
            endl <<
            "                                            [ Positive integer ]" <<
            endl <<
            "       --------------------------------------------------------------------"
             << endl << endl;




        for (size_t i = 0; i < int_par_names.size (); i++)
	{

            cout << "       ";
            cout.width (8);
            cout.setf (ios_base::left, ios_base::adjustfield);
            cout << int_option1[int_par_names[i]];
            cout.width (24);
            cout << int_option2[int_par_names[i]] << "  :  " <<
                int_help[int_par_names[i]] << endl;


	}

        cout << endl << endl;

//===================================== end integer parameters help ======================================


        cout << "       Logic parameters.                    Valid options:" <<
            endl <<
            "                                            [yes, YES, true, no, NO, false]"
             << endl <<
            "       --------------------------------------------------------------------"
             << endl << endl;



        for (size_t i = 0; i < logic_par_names.size (); i++)
	{

            cout << "       ";
            cout.width (8);
            cout.setf (ios_base::left, ios_base::adjustfield);
            cout << logic_option1[logic_par_names[i]];
            cout.width (24);
            cout << logic_option2[logic_par_names[i]] << "  :  " <<
                logic_help[logic_par_names[i]] << endl;


	}

        cout << endl << endl;


//===================================== end logic parameters help ======================================



        cout << "       Text parameters.                     Valid options:" <<
            endl << "                                            [ text ]" << endl
             <<
            "       --------------------------------------------------------------------"
             << endl << endl;




        for (size_t i = 0; i < string_par_names.size (); i++)
	{

            cout << "       ";
            cout.width (8);
            cout.setf (ios_base::left, ios_base::adjustfield);
            cout << string_option1[string_par_names[i]];
            cout.width (24);
            cout << string_option2[string_par_names[i]] << "  :  " <<
                string_help[string_par_names[i]] << endl;


	}

        cout << endl << endl;

//===================================== end string parameters help ======================================



        cout << "       Puncture parameters.                 Valid options:" <<
            endl << "                                            [ double ]" <<
            endl <<
            "       --------------------------------------------------------------------"
             << endl << endl <<
            "       -bhmp[i]                          :  Black hole mass parameter i "
             << endl <<
            "       -bhqp[i]                          :  Black hole charge parameter i "
             << endl <<
            "       -bhx[i]                           :  Black hole initial position x direction  "
             << endl <<
            "       -bhy[i]                           :  Black hole initial position y direction  "
             << endl <<
            "       -bhz[i]                           :  Black hole initial position z direction  "
             << endl <<
            "       -bhpx[i]                          :  Black hole linear momentum x direction  "
             << endl <<
            "       -bhpy[i]                          :  Black hole linear momentum y direction  "
             << endl <<
            "       -bhpz[i]                          :  Black hole linear momentum z direction  "
             << endl <<
            "       -bhsx[i]                          :  Black hole spin x direction  "
             << endl <<
            "       -bhsy[i]                          :  Black hole spin y direction  "
             << endl <<
            "       -bhsz[i]                          :  Black hole spin z direction  "
             << endl << endl << endl

             << " ======  Olliptic 9.10  ======" << endl << endl
             << "     Author: Pablo Galaviz    " << endl << endl
             << "   pablo.galaviz@uni-jena.de  " << endl << endl
            << " =============================" << endl << endl
             << "   License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html> " << endl
             << "   This is free software: you are free to change and redistribute it." << endl
             << "   There is NO WARRANTY, to the extent permitted by law."
             << endl << endl;

#ifdef OLLIN_MPI

    }

#endif

    exec = false;

}







void elliptic::IO_bam (string file){ 
  

  size_t boxes = 6;

  string ending[] = { "", "a", "b", "c", "d", "e" };


  
  size_t level = 0;

  
  vector < string > name_files;

  size_t num_arch;


  do
    {


      stringstream sl;

      sl << level;
      
      num_arch = 0;

      for (size_t ibox = 0; ibox < boxes; ibox++)
	{


	  string iofile = file + sl.str () + ending[ibox];

	  ifstream File (iofile.c_str ());


	  if (File.is_open ())
	    {

	      num_arch++;

	      if (my_rank == 0 && output == YES)


		cout << "File found: " << iofile << endl;


	      File.close ();

	      name_files.push_back (iofile);


	    }


	  File.close ();


	}

      level++;



    }
  while (num_arch > 0);



  for (size_t ibox = 0; ibox < name_files.size (); ibox++)

    vars.Print_IO_BAM (name_files[ibox], "u", LAGRANGE, order + 1);


}



void elliptic::IO_Zcode(string file){


  size_t boxes = 6;

  size_t level = 0;


  vector < string > name_files;

  size_t num_arch;


  string preL;

  /*
  stringstream ls_str;

  ls_str << "ls -1 " << file << "* > tmp.lst";

  cout << "searching for files: " << file << "*" << endl;

  int sout = system(ls_str.str().c_str());

  cout << "...done" << endl; 

  ifstream File_tmp ("tmp.lst");


  if (File_tmp.is_open ()){

    string namef;

    name_files.resize(0);

    do{

      File_tmp >> namef;   

      name_files.push_back (namef);
    

    }while(!File_tmp.eof());

    File_tmp.close();

    sout = system("rm tmp.lst");

  }
  else
    cout << "non files found..." << endl; 

  name_files.pop_back();

  cout << name_files.size() << " file(s) found: " << endl;

  for(int i=0; i < name_files.size(); i++)
    cout << name_files[i] << endl; 

  */

  
  do
    {


      stringstream sl;

      sl << level;

      if(level < 10)
	preL = "0";
      else
	preL = "";


      num_arch = 0;

      for (size_t ibox = 0; ibox < boxes; ibox++){
	  
	stringstream sb;

	sb << ibox;

	string iofile = file + preL + sl.str () + "-0" + sb.str() + ".mgid";


	ifstream File (iofile.c_str ());


	if (File.is_open ())
	  {
	      
	    num_arch++;

	    if (my_rank == 0 && output == YES)


	      cout << "File found: " << iofile << endl;


	    File.close ();

	    name_files.push_back (iofile);


	  }


      }

      level++;

    } while (num_arch > 0);



  vector < string > shell_files;


  shell_files.push_back( file + "SH-xm.mgid");
  shell_files.push_back( file + "SH-xp.mgid");

  shell_files.push_back( file + "SH-ym.mgid");
  shell_files.push_back( file + "SH-yp.mgid");

  shell_files.push_back( file + "SH-zm.mgid");
  shell_files.push_back( file + "SH-zp.mgid");


  for(int i=0; i<shell_files.size(); i++){

    ifstream File (shell_files[i].c_str ());

    if (File.is_open ())
      {
	      
	num_arch++;

	if (my_rank == 0 && output == YES)


	  cout << "File found: " << shell_files[i] << endl;


	File.close ();

	name_files.push_back (shell_files[i]);


      }

  }



  for (size_t ibox = 0; ibox < name_files.size (); ibox++)

    vars.Print_IO_Zcode (name_files[ibox], "u", LAGRANGE, order + 1);
  


}





void
elliptic::Print_Solution (double Dxyz, double Lx, double Ly, double Lz)
{


    if (my_rank == 0 && output == YES)

        cout << "Printing..." << endl;

    if (print_1d)
    {
        for (size_t l = 0; l < vars.Get_levels (); l++)

            for (size_t i = 0; i < var_names.size (); i++)
            {


                vars.Print_1D (dir_name, var_names[i], l, CUT_X, 0, 0,
                               interpol_output, NULL, LAGRANGE, order + 1);

                vars.Print_1D (dir_name, var_names[i], l, CUT_Y, 0, 0,
                               interpol_output, NULL, LAGRANGE, order + 1);

                vars.Print_1D (dir_name, var_names[i], l, CUT_Z, 0, 0,
                               interpol_output, NULL, LAGRANGE, order + 1);



            }




    }



    if (print_2d)

        for (size_t l = 0; l < vars.Get_levels (); l++)

            for (size_t i = 0; i < var_names.size (); i++)
            {

                vars.Print_2D (dir_name, var_names[i], l, CUT_X, 0,
                               interpol_output);

                vars.Print_2D (dir_name, var_names[i], l, CUT_Y, 0,
                               interpol_output);

                vars.Print_2D (dir_name, var_names[i], l, CUT_Z, 0,
                               interpol_output);

            }



    if (print_3d)

        for (size_t l = 0; l < vars.Get_levels (); l++)

            for (size_t i = 0; i < var_names.size (); i++)

                vars.Print_3D (dir_name, var_names[i], l);







    if (print_interpol)
        for (size_t i = 0; i < var_names.size (); i++)
        {


            vars.Print_Interpol_1D (dir_name, var_names[i], Lx, Dxyz, CUT_X, 0, 0,
                                    NULL, LAGRANGE, order + 1);

            vars.Print_Interpol_1D (dir_name, var_names[i], Ly, Dxyz, CUT_Y, 0, 0,
                                    NULL, LAGRANGE, order + 1);

            vars.Print_Interpol_1D (dir_name, var_names[i], Lz, Dxyz, CUT_Z, 0, 0,
                                    NULL, LAGRANGE, order + 1);


            vars.Print_Interpol_2D (dir_name, var_names[i], Lx, Ly, -0.1, -0.1, CUT_Z, 0,
                                    NULL, LAGRANGE, order + 1);




        }



    if( equation==BOWEN_YORK )
      Print_Coeff ();


    if (my_rank == 0 && output == YES)

        cout << "done" << endl << flush;



}






void
elliptic::Print_Coeff ()
{



    if (print_1d)

        for (size_t l = 0; l < vars.Get_levels (); l++)
        {

            for (size_t i = 0; i < coeff_names.size (); i++)
            {


                vars.Print_1D (dir_name, coeff_names[i], l, CUT_X, 0, 0);

                vars.Print_1D (dir_name, coeff_names[i], l, CUT_Y, 0, 0);

                vars.Print_1D (dir_name, coeff_names[i], l, CUT_Z, 0, 0);



            }


            for (size_t i = 0; i < source_names.size (); i++)
            {


                vars.Print_1D (dir_name, source_names[i], l, CUT_X, 0, 0);

                vars.Print_1D (dir_name, source_names[i], l, CUT_Y, 0, 0);

                vars.Print_1D (dir_name, source_names[i], l, CUT_Z, 0, 0);



            }


        }


    if (print_2d)

        for (size_t l = 0; l < vars.Get_levels (); l++)
        {
            for (size_t i = 0; i < coeff_names.size (); i++)
            {

                vars.Print_2D (dir_name, coeff_names[i], l, CUT_X, 0);

                vars.Print_2D (dir_name, coeff_names[i], l, CUT_Y, 0);

                vars.Print_2D (dir_name, coeff_names[i], l, CUT_Z, 0);

            }


            for (size_t i = 0; i < source_names.size (); i++)
            {


                vars.Print_2D (dir_name, source_names[i], l, CUT_X, 0);

                vars.Print_2D (dir_name, source_names[i], l, CUT_Y, 0);

                vars.Print_2D (dir_name, source_names[i], l, CUT_Z, 0);



            }


        }

    if (print_3d)

        for (size_t l = 0; l < vars.Get_levels (); l++)
        {
            for (size_t i = 0; i < coeff_names.size (); i++)

                vars.Print_3D (dir_name, coeff_names[i], l);


            for (size_t i = 0; i < source_names.size (); i++)

                vars.Print_3D (dir_name, source_names[i], l);

        }



}







void
elliptic::Info_Punctures ()
{




    if (bhx.size () == NumBH &&
        bhy.size () == NumBH &&
        bhz.size () == NumBH &&
        bhpx.size () == NumBH &&
        bhpy.size () == NumBH &&
        bhpz.size () == NumBH && bhsx.size () == NumBH && bhsy.size () == NumBH)
    {


#ifdef OLLIN_MPI



        if (my_rank == 0)
	{

#endif


            cout << "Calculating initial data for " << NumBH <<
                " puncture(s), with parameters: " << endl;

            cout.setf (ios_base::scientific, ios_base::floatfield);
            cout.precision (7);

            for (size_t i = 0; i < NumBH; i++)
	    {

                cout <<
                    "--------------------------------------------------------------------"
                     << endl << "Puncture " << i + 1 << "\t  |\tMass parameter : " << bhmp[i] << "\t |\tCharge parameter : " << bhqp[i] << endl <<
                    "____________________________________________________________________"
                     << endl;


                cout.setf (ios_base::left, ios_base::adjustfield);
                cout.width (18);
                cout << " " << "|";

                cout.width (10);
                cout.setf (ios_base::internal, ios_base::adjustfield);
                cout << "x";

                cout.width (8);
                cout.setf (ios_base::internal, ios_base::adjustfield);
                cout << "|";

                cout.width (8);
                cout.setf (ios_base::internal, ios_base::adjustfield);
                cout << "y";

                cout.width (8);
                cout.setf (ios_base::internal, ios_base::adjustfield);
                cout << "|";


                cout.width (8);
                cout.setf (ios_base::internal, ios_base::adjustfield);
                cout << "z" << endl;



                cout.setf (ios_base::left, ios_base::adjustfield);
                cout.width (18);
                cout << "Position" << ":";

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhx[i];

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhy[i];

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhz[i] << endl;



                cout.setf (ios_base::left, ios_base::adjustfield);
                cout.width (18);
                cout << "Linear momentum" << ":";


                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhpx[i];

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhpy[i];

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhpz[i] << endl;


                cout.setf (ios_base::left, ios_base::adjustfield);
                cout.width (18);
                cout << "Angular momentum" << ":";


                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhsx[i];

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhsy[i];

                cout.width (16);
                cout.setf (ios_base::right, ios_base::adjustfield);
                cout << bhsz[i] << endl
                     <<
                    "--------------------------------------------------------------------"
                     << endl << endl;

	    }



#ifdef OLLIN_MPI

	}

#endif






    }


    else
    {


#ifdef OLLIN_MPI


        if (my_rank == 0)
	{

#endif



            cout << "\n\n*** Error: Punture parameters" << endl;


            exit (1);


#ifdef OLLIN_MPI



	}

#endif









    }





}





void
elliptic::Print_Info ()
{






    string order_str;

    string method_str;

    string equation_str;

    string boundary_str;

    string dom_str;

    switch (equation)
    {


        case POISSON:

            equation_str = "Poisson";

            break;

        case BOUNDARY_TEST:

            equation_str = "Robin boundary test";

            break;

        case BRILL_WAVES:

            equation_str = "Brill waves";

            break;

        case PUNCTURES:

            equation_str = "Punctures ID";

            break;

        case NSO_ID:

            equation_str = "Neutron star oscilations ID";

            break;

        case MULTIBOX_TEST:

            equation_str = "Multibox test";

            break;

        case GENERIC:

            equation_str = "Generic";

            break;

        case SYSTEM_TEST:

            equation_str = "System of equations test";

            break;


        case TRUMPET_1pLOG:

            equation_str = "Trumpet 1+log black-hole initial data";

            break;

        case BOWEN_YORK:

            equation_str = "Bowen-York momentum constrain test";

            break;


    }


    switch (method)
    {



        case GAUSS_SEIDEL:

            method_str = "Gauss - Seidel";

            break;

        case FMG:

            method_str = "Full Multigrid";

            break;

        case VCYCLE:

            method_str = "Multigrid V-Cycles";

            break;

        case WCYCLE:

            method_str = "Multigrid W-Cycle";

            break;

        case APP:

            method_str = "Analitical Approximate";

            break;

    }

    switch (order)
    {

        case 2:

            order_str = "2nd";

            break;


        case 4:

            order_str = "4th";

            break;

        case 6:

            order_str = "6th";

            break;

        case 8:

            order_str = "8th";

            break;


        case 10:

            order_str = "10th";

            break;



    }


    switch (boundary)
    {


        case DIRICHLET:

            boundary_str = "Dirichlet";

            break;

        case ROBIN:

            boundary_str = "Robin";

            break;

        case ASYMPTOTIC:

            boundary_str = "Asymptotic";

            break;

    }



    switch (dom)
    {


        case FULL:

            dom_str = "Full domain";

            break;


        case OCTANT:

            dom_str = "Octant domain: x < 0, y < 0 and z < 0 ";

            break;

        case BITANT_X:

            dom_str = "Bitant domain: x < 0";

            break;

        case BITANT_Y:

            dom_str = "Bitant domain: y < 0";

            break;

        case BITANT_Z:

            dom_str = "Bitant domain: z < 0";

            break;

        case QUADRANT_XY:

            dom_str = "Quadrant domain: x < 0 and y < 0";

            break;

        case QUADRANT_XZ:

            dom_str = "Quadrant domain: x < 0 and z < 0";

            break;


        case QUADRANT_YZ:

            dom_str = "Quadrant domain: y < 0 and z < 0";

            break;

        default:

            dom_str = "Unknow symmetry, set to full domain.";



    }

    stringstream std_output;

    if (output == YES || output == PARTIAL)
    {

        std_output <<
            " ---------------------------- Olliptic -----------------------------"
                   << endl << "Equation\t:\t" << equation_str << endl << "Method\t\t:\t"
                   << method_str << endl << "Tolerance " << endl << " Fine grid\t:\t" <<
            eps_solv << endl << " Coarse grid\t:\t" << eps_coarse << endl <<
            "Smooth sweeps " << endl << " Pre\t\t:\t" << nu1 << endl <<
            " Post\t\t:\t" << nu2 << endl << "Levels\t\t:\t" << levels -
            mgrid_levels << endl << "Boxes\t\t:\t" << boxes << endl <<
            "Order\t\t:\t" << order_str << endl << endl << "Boundary\t:\t" <<
            boundary_str << endl << "Symmetry\t:\t" << dom_str << endl << endl <<
            "Output to\t:\t" << dir_name << "*" << endl;





        for (size_t l = 0; l < levels; l++)

            for (size_t i = 0; i < vars.Get_num_dom (l); i++)
            {

                std_output <<
                    "\n  __________________________________________________________________"
                           << endl << "  level: " << l << " | " << "box: " << i << " | ";

                vars.Get_ffunction (l, i, var_names[0])->Print_Info (std_output);

            }

        std_output <<
            "--------------------------------------------------------------------"
                   << endl << "  Init Multigrid iteration" << endl << flush;

    }




    stringstream file_output;

    file_output << dir_name << "stdout" << my_rank << ".out";

    string outfile = file_output.str ();

    stdout_file.open (outfile.c_str (), fstream::app);


    stdout_file << std_output.str ();


    stdout_file.close ();

    if (my_rank == 0)

        cout << std_output.str ();



}








void
elliptic::Print_Pre_Cycle (size_t k, size_t k_top, double norm, size_t cycle)
{




    if (my_rank == 0 && output == YES)
    {

        cout.setf (ios_base::scientific, ios_base::floatfield);

        cout.precision (8);


        if (k == k_top)
	{

            cout << endl << endl;

            cout << "Cycle: " << cycle << endl;


            for (size_t j = 0; j <= k_top; j++)

                cout << "========================";

            cout << endl;

	}


        for (size_t j = 0; j < k; j++)

            cout << "|\t\t\t";

        cout << "|Loo[" << k << "] : " << norm;


        for (size_t j = k + 1; j <= k_top + 1; j++)

            cout << "|\t\t\t";

        cout << endl;

    }


    if (output == YES)
    {

        stdout_file.open (stdout_file_name.c_str (), fstream::app);


        stdout_file.setf (ios_base::scientific, ios_base::floatfield);

        stdout_file.precision (8);



        if (k == k_top)
	{

            stdout_file << endl << endl << "Cycle: " << cycle << endl;


            for (size_t j = 0; j <= k_top; j++)

                stdout_file << "========================";

            stdout_file << endl;

	}




        for (size_t j = 0; j < k; j++)

            stdout_file << "|\t\t\t";

        stdout_file << "|Loo[" << k << "] : " << norm;


        for (size_t j = k + 1; j <= k_top + 1; j++)

            stdout_file << "|\t\t\t";

        stdout_file << endl << flush;


        stdout_file.close ();


        norm_file.open (norm_file_name.c_str (), fstream::app);

        norm_file.setf (ios_base::scientific, ios_base::floatfield);

        norm_file.precision (18);


        norm_file << cycle << "\t" << norm << endl << flush;

        norm_file.close ();


    }



    if (print_1dt)

        for (size_t iv = 0; iv < var_names.size (); iv++)

            vars.Print_t1D (dir_name, var_names[iv], k, it);




}









void
elliptic::Print_Post_Cycle (size_t k, size_t k_top, double norm, size_t cycle)
{




    if (my_rank == 0 && output == YES)
    {


        for (size_t j = 0; j < k; j++)

            cout << "|\t\t\t";

        cout << "|Loo[" << k << "] : " << norm;


        for (size_t j = k + 1; j <= k_top + 1; j++)

            cout << "|\t\t\t";

        cout << endl;





        if (k == k_top)
	{


            for (size_t j = 0; j <= k_top; j++)

                cout << "========================";


            cout << endl << norm << "  -->  " << eps_solv << endl;


	}


    }




    if (output == YES)
    {


        stdout_file.open (stdout_file_name.c_str (), fstream::app);


        stdout_file.setf (ios_base::scientific, ios_base::floatfield);

        stdout_file.precision (8);




        for (size_t j = 0; j < k; j++)

            stdout_file << "|\t\t\t";

        stdout_file << "|Loo[" << k << "] : " << norm;


        for (size_t j = k + 1; j <= k_top + 1; j++)

            stdout_file << "|\t\t\t";

        stdout_file << endl;




        if (k == k_top)
	{


            for (size_t j = 0; j <= k_top; j++)

                stdout_file << "========================";


            stdout_file << endl << norm << "  -->  " << eps_solv << endl;




	}

        stdout_file.close ();




    }



    if (print_1dt)

        for (size_t iv = 0; iv < var_names.size (); iv++)

            vars.Print_t1D (dir_name, var_names[iv], k, it);



}










void
elliptic::Print_Cycle (size_t k, size_t k_top, double norm, size_t cycle)
{




    if (my_rank == 0 && output == YES)
    {



        cout.setf (ios_base::scientific, ios_base::floatfield);

        cout.precision (8);


        //if( k == 0 ){

        cout << endl << endl;

        cout << "Cycle: " << cycle << endl;


        for (size_t j = 0; j <= k_top; j++)

            cout << "========================";

        cout << endl;

        //}


        for (size_t j = 0; j < k; j++)

            cout << "|\t\t\t";

        cout << "|Loo[" << k << "] : " << norm;


        for (size_t j = k + 1; j <= k_top + 1; j++)

            cout << "|\t\t\t";

        cout << endl;


        //if( k == k_top ){


        for (size_t j = 0; j <= k_top; j++)

            cout << "========================";


        cout << endl << norm << "  -->  " << eps_solv << endl;


        //}


    }




    if (output == YES)
    {




        stdout_file.open (stdout_file_name.c_str (), fstream::app);


        stdout_file.setf (ios_base::scientific, ios_base::floatfield);

        stdout_file.precision (8);



        if (k == 0)
	{

            stdout_file << endl << endl << "Cycle: " << cycle << endl;


            for (size_t j = 0; j <= k_top; j++)

                stdout_file << "========================";

            stdout_file << endl;

	}




        for (size_t j = 0; j < k; j++)

            stdout_file << "|\t\t\t";

        stdout_file << "|Loo[" << k << "] : " << norm;


        for (size_t j = k + 1; j <= k_top + 1; j++)

            stdout_file << "|\t\t\t";

        stdout_file << endl << flush;




        norm_file.open (norm_file_name.c_str (), fstream::app);

        norm_file.setf (ios_base::scientific, ios_base::floatfield);

        norm_file.precision (18);


        norm_file << cycle << "\t" << norm << endl << flush;

        norm_file.close ();





        if (k == k_top)
	{


            for (size_t j = 0; j <= k_top; j++)

                stdout_file << "========================";


            stdout_file << endl << norm << "  -->  " << eps_solv << endl;




	}

        stdout_file.close ();




    }



    if (print_1dt){

        for (size_t iv = 0; iv < var_names.size (); iv++)

            vars.Print_t1D (dir_name, var_names[iv], k, it);

	
        for (size_t iv = 0; iv < aux_names.size (); iv++)

            vars.Print_t1D (dir_name, aux_names[iv], k, it);


    }

}
