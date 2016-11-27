//============================= Olliptic main.cc =========================//
//  
//  Olliptic (description... fill)     
//
//
//
//
//
//============================ Pablo Galaviz 2009 ========================//









//============================ Standar libraries =========================//


#include <math.h>

#include <fstream>

#include <iostream>

#include <string>

#include <sstream>

//#include <curses.h>


//============================= Olliptic classes =========================//


#include "interface.h"


#include "elliptic.h"



//=============================== MPI library ============================//


#ifdef OLLIN_MPI

#include <mpi.h>

#endif


//============================== End headers =============================//















using namespace std;


//================================== main ==================================//
//  (description of ... fill)
//
//
//
//================================== ---- ==================================//



int
main (int argc, char *argv[])
{




#ifdef OLLIN_MPI

  //== Initi MPI library ==//

  MPI::Init ();

#endif




  interface Oll (argc, argv);	



  elliptic E (&Oll);		//== Ellitic class with default parameters ==//





  E.Print_Info ();



  //== Solve the equations ==//

  
  E.Solve ();



  //== Print the solution



  E.Print_Solution (Oll.Get_double_par ("Dprint"),
		    Oll.Get_double_par ("Lx"),
		    Oll.Get_double_par ("Ly"), Oll.Get_double_par ("Lz"));


  
  if (Oll.Get_logic_par ("io_bam"))

    E.IO_bam (Oll.Get_string_par ("file_bam"));




  if (Oll.Get_logic_par ("io_Zcode"))

    E.IO_Zcode (Oll.Get_string_par ("file_bam"));








#ifdef MPI_OLLIN

  //== Finish MPI ==//    

  MPI::COMM_WORLD.Barrier ();


  MPI::Finalize ();

#endif



  //== End ==//

  return (0);



}
