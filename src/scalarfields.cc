//
// scalarfields.cc
//  
// Made by Pablo Galaviz
// Login   <pablo@NerV>
// 
// Started on  Thu Oct 22 13:35:47 2009 Pablo Galaviz
// Started on  Thu Oct 22 13:35:47 2009 Pablo Galaviz
//


#include "scalarfields.h"



void
scalarfields::Make (domain * Dom, vector < string > name)
{

  D = Dom;

  NameFields.clear ();

  NumSFs = 0;

  u.clear ();


  if (name.size () == 0)
    {

      New_Field ("u");

    }
  else
    {

      for (size_t i = 0; i < name.size (); i++)

	New_Field (name[i]);
    }

}





void
scalarfields::Make (domain * Dom, string name)
{

  D = Dom;

  NameFields.clear ();

  NumSFs = 0;

  u.clear ();

  New_Field (name);

}
