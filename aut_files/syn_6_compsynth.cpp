/* ==========================================================================
* phonebook.c
* example for illustrating how to manipulate structure and cell array
*
* takes a (MxN) structure matrix and returns a new structure (1x1)
* containing corresponding fields: for string input, it will be (MxN)
* cell array; and for numeric (noncomplex, scalar) input, it will be (MxN)
* vector of numbers with the same classID as input, such as int, double
* etc..
*
* This is a MEX-file for MATLAB.
* Copyright 1984-2011 The MathWorks, Inc.
*==========================================================================*/

#include "mex.h"
#include <stdlib.h> 
#include <iostream>
#include "string.h"
#include <iterator>
#include "matrix.h"
#include <sstream>
#include "D:\IEEE_300_Bus\libFAUDES\Include\libfaudes.h"
// make faudes namespace available

using namespace std;
using std::string;

using namespace faudes;

#define ROWS 1

#define COLUMNS 30000

#define ELEMENTS 30000



#define MAXCHARS 80   /* max length of string contained in each field */

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char *Automaton1;              /* Automaton 1 object file name */
    char *Automaton2;               /*  Automaton 2 object file name */
    size_t ncols;                   /* size of matrix */
    char *outMatrix;              /* output matrix */

	
	/* check proper input and output */
	

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }

   /* get the values of the input  */
     Automaton1 = (char*)mxGetData(prhs[0]);
     Automaton2 = (char*)mxGetData(prhs[1]);
	

    // define Automates
     // Aotomaton1
    char plant[12];
    plant[0] = Automaton1[0]; 
    plant[1] = Automaton1[2]; 
    plant[2] = Automaton1[4];
	plant[3] = Automaton1[6];	
	plant[4] = Automaton1[8];	
	plant[5] = Automaton1[10];
	plant[6] = Automaton1[12];
	plant[7] = Automaton1[14];
	plant[8] = Automaton1[16];
  	plant[9] = Automaton1[18];
   	plant[10] = Automaton1[20];
    plant[11] = '\0';
     // Aotomaton2
    char spec[15];
    spec[0] = Automaton2[0]; 
    spec[1] = Automaton2[2];
    spec[2] = Automaton2[4];
    spec[3] = Automaton2[6];
    spec[4] = Automaton2[8];
    spec[5] = Automaton2[10];
    spec[6] = Automaton2[12];
    spec[7] = Automaton2[14];
    spec[8] = Automaton2[16];
    spec[9] = Automaton2[18];
    spec[10] = Automaton2[20];
    spec[11] = Automaton2[22];
    spec[12] = Automaton2[24];
    spec[13] = Automaton2[26];
    spec[14] = '\0';
    
    // we ay need to create sytamtic way to generate names for automates as results of parallel operation 
    
    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = (char*)mxGetData(plhs[0]);

    
    
    //////////////////////////////
	// DES controller synthesis // 
	//////////////////////////////

	// Compose plant dynamics from two very simple machines 
//   Generator tempgen, machinea, machineb;
  Generator ccplant(plant); 
  System cplant(ccplant); 
 
//   tempgen.Read("libFAUDES/Examples/synthesis/data/verysimplemachine.gen");
//   tempgen.Version("1",machinea);
//   tempgen.Version("2",machineb);
//   Parallel(machinea,machineb,cplant);
  // Declare controllable events
  // maybe define controlability in the defintion of basic automates not here
  EventSet contevents;
// I beleive self loops of tranmission lines are not controllable we only consider self loops of generators and loads
  contevents.Insert("k15"); contevents.Insert("k16"); contevents.Insert("k17"); contevents.Insert("k18"); contevents.Insert("k19"); contevents.Insert("f10"); contevents.Insert("b06");


// the below scentence may need to moved when declaring the result pf parallel operation
   cplant.SetControllable(contevents);
//  Write to file
//   cplant.Write("tmp_cplant12.gen");
  // Report to console
  std::cout << "################################\n";
  std::cout << "# tutorial, plant model \n";
  cplant.DWrite();
  std::cout << "################################\n";
  // Read specification 
  Generator specification(spec);
//    specification.Read("libFAUDES/Examples/synthesis/data/buffer.gen");
 
  InvProject(specification,cplant.Alphabet()); 
  specification.Name("simple machines specification");
  // Write to file
  specification.Write("tmp_specification12.gen");
  // Report to console
  std::cout << "################################\n";
  std::cout << "# tutorial, specification \n";
  specification.DWrite();
  std::cout << "################################\n";
  // Run synthesis algorithm
  System supervisor;
  SupConNB(cplant,specification,supervisor);
  supervisor.Name("simple machines supervisor");
  supervisor.Write("tmp_supervisor12.gen");
  // Report to console
  std::cout << "################################\n";
  std::cout << "# tutorial, supervisor\n";
  supervisor.DWrite();
  std::cout << "################################\n";
  return ;
}
