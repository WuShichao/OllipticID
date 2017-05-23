//
// interface.h
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


#ifndef INTERFACE_H
#define INTERFACE_H



//============================ Standar libraries =========================//


#include <vector>

#include <string>

#include <sstream>

#include <iostream>

#include <stdlib.h>

#include <math.h>

#include <map>

#include <cstring>


//=============================== MPI library ============================//


#ifdef OLLIN_MPI

#include <mpi.h>

#endif


//============================== End headers =============================//





//== Boundary ==//
enum t_boundary
{ DIRICHLET, ROBIN, ASYMPTOTIC };


//== Method ==//
enum t_method
{ GAUSS_SEIDEL, FMG, VCYCLE, WCYCLE, APP };

//== Elliptic problem ==//
enum t_equation
{ POISSON,
  BOUNDARY_TEST,
  BRILL_WAVES,
  PUNCTURES,
  MULTIBOX_TEST,
  GENERIC,
  SYSTEM_TEST,
  NSO_ID,
  METRIC_TEST,
  TRUMPET_1pLOG,
  PUNCTURES_SCALAR,
  PUNCTURES_EM,
  BOWEN_YORK,
  BY_EM,
  GAUSS_BONNET
};

//== Domain symmetry ==//
enum t_domain
{ FULL,
  BITANT_X, BITANT_Y, BITANT_Z,
  QUADRANT_XY, QUADRANT_XZ, QUADRANT_YZ,
  OCTANT, ROTATE_X, ROTATE_Y, ROTATE_Z
};


//== standard output ==// 
enum t_output
{ YES, PARTIAL, NO };


//== end ==//


using namespace std;




class interface
{



    string param;


//== user parameters ==//


//== Here we define the basic kind of parameters: logic, integers,
//   doubles (floats) and strings. 

//== A vector of string (*_par_names) stores the name of the parameters,
//   a map *_par is used to link the name of the parameter with its value,
//   *_help store a short string with a help comment, *_option1 and
//   *_option2 define strings used as keys in the command line.
    
    

//== Logic parameters

    vector < string > logic_par_names;

    map < string, bool > logic_par;

    map < string, string > logic_help;

    map < string, string > logic_option1;

    map < string, string > logic_option2;




//== Positive integers

    vector < string > int_par_names;

    map < string, size_t > int_par;

    map < string, string > int_help;

    map < string, string > int_option1;

    map < string, string > int_option2;




//== Positive doubles 

    vector < string > double_par_names;

    map < string, double >double_par;

    map < string, string > double_help;

    map < string, string > double_option1;

    map < string, string > double_option2;




// strings

    vector < string > string_par_names;

    map < string, string > string_par;

    map < string, string > string_help;

    map < string, string > string_option1;

    map < string, string > string_option2;



//= Here we define special variables to store parameters with keyword
 

    t_boundary boundary;

    t_method method;

    t_equation equation;

    t_domain domain_sym;

    t_output output;



//== Parameters for puncture initial data


//== mass parameter of each puncture
    vector < double >bhmp;

//== charge parameter of each puncture
    vector < double >bhqp;

//== Position
    vector < double >bhx;

    vector < double >bhy;

    vector < double >bhz;

//== Linear momentum
    vector < double >bhpx;

    vector < double >bhpy;

    vector < double >bhpz;

//== Spin
    vector < double >bhsx;

    vector < double >bhsy;

    vector < double >bhsz;



//== True if the user gives Lx, Ly and Lz
    bool by_domain;

//== True if the user gives nx, ny and nz
    bool by_points;

//== True if the user gives dx, dy and dz
    bool by_grid_size;

//== True if everything is right
    bool exec;

//== True if the user defines an equation to solve.
    bool def_eq;



//== Methods for add parameters, provide parameter's name, default
//   value, a short description of the functionality and two strings
//   defining the keys to set the parameters.  
    
    void add_logic_param (string name_par, bool def, string description,
                          string option_key1, string option_key2);

    void add_int_param (string name_par, size_t def, string description,
                        string option_key1, string option_key2);

    void add_double_param (string name_par, double def, string description,
                           string option_key1, string option_key2);

    void add_string_param (string name_par, string def, string description,
                           string option_key1, string option_key2);

//== Methods to generate the string to describe the functionality of
//   some of the parameter. 
    
    string string_equation ();

    string string_method ();

    string string_domain ();

    string string_boundary ();

//== Method to print a Welcome string to the standard output:
    
    void print_welcome (); //== Definition in output.cc


public:



//== Constructor, use as argument the command line from main.
    interface (int argc, char *argv[]);


//== Destructor, nothing to do.
    ~interface (){};


//== Print a text with the usage.
    void print_help (); //== Definition in output.cc


//== Process the parameters.
    void proc_parameters ();

//== Read the command line text and catch the parameters.
    void read_console (int argc, char *argv[]);



//== Methods to return the parameters. ==//

//== To get the value of the parameter we provide the name of the
//   parameter (defined in add_*_param) and the function returns the
//   value given after the keyword

    inline const bool Get_logic_par (string par)
        {
            return (logic_par[par]);
        }

    inline const size_t Get_integer_par (string par)
        {
            return (int_par[par]);
        }

    inline const double Get_double_par (string par)
        {
            return (double_par[par]);
        }

    inline const string Get_string_par (string par)
        {
            return (string_par[par]);
        }



    inline vector < double >Get_bhmp () const
        {
            return (bhmp);
        };

    inline vector < double >Get_bhqp () const
        {
            return (bhqp);
        };

    inline vector < double >Get_bhx () const
        {
            return (bhx);
        };

    inline vector < double >Get_bhy () const
        {
            return (bhy);
        };

    inline vector < double >Get_bhz () const
        {
            return (bhz);
        };


    inline vector < double >Get_bhpx () const
        {
            return (bhpx);
        };

    inline vector < double >Get_bhpy () const
        {
            return (bhpy);
        };

    inline vector < double >Get_bhpz () const
        {
            return (bhpz);
        };


    inline vector < double >Get_bhsx () const
        {
            return (bhsx);
        };

    inline vector < double >Get_bhsy () const
        {
            return (bhsy);
        };

    inline vector < double >Get_bhsz () const
        {
            return (bhsz);
        };



    inline t_output Get_TOutput () const
        {
            return (output);
        }

    inline t_boundary Get_Boundary () const
        {
            return (boundary);
        }

    inline t_method Get_Method () const
        {
            return (method);
        }

    inline t_equation Get_Equation () const
        {
            return (equation);
        }

    inline t_domain Get_Domain_Sym () const
        {
            return (domain_sym);
        }



    //== End ==//


};


//== Function to check for a valid number of points for multigrid method.
void Valid_N (size_t & N);


//== Function to abort the execution after an error

void error_exit(string error_string1,
                string error_string2="",
                string error_string3="",
                bool warning = false); //== Definition in output.cc
         

#endif // INTERFACE_H
