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
#include "libFAUDES\include\libfaudes.h"
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
    char st1[12];
    st1[0] = Automaton1[0]; 
    st1[1] = Automaton1[2]; 
    st1[2] = Automaton1[4];
	st1[3] = Automaton1[6];	
	st1[4] = Automaton1[8];	
	st1[5] = Automaton1[10];
	st1[6] = Automaton1[12];
	st1[7] = Automaton1[14];
	st1[8] = Automaton1[16];
  	st1[9] = Automaton1[18];
   	st1[10] = Automaton1[20];
    st1[11] = '\0';
     // Aotomaton2
    char ev1[12];
    ev1[0] = Automaton2[0]; 
    ev1[1] = Automaton2[2];
    ev1[2] = Automaton2[4];
    ev1[3] = Automaton2[6];
    ev1[4] = Automaton2[8];
    ev1[5] = Automaton2[10];
    ev1[6] = Automaton2[12];
    ev1[7] = Automaton2[14];
    ev1[8] = Automaton2[16];
    ev1[9] = Automaton2[18];
    ev1[10] = Automaton2[20];
    ev1[11] = '\0';
    
    // we ay need to create sytamtic way to generate names for automates as results of parallel operation 
    
    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = (char*)mxGetData(plhs[0]);

    
    
    ////////////////////////////
	// synchronous composition
	////////////////////////////

	// read generators
	Generator Com1(st1);
	Generator Com2(ev1);
	//Generator line3("line3.gen");
	//Generator line1("line1.gen");



	// perform composition -----> we may need to perform composition through a for loop  
	Generator parallel_Com1Com2; // parallel operation of component 1 and component 2 
	//Generator parallel_g1L3;
	//Generator parallel_g1L1;
	

		
	Parallel(Com1, Com2, parallel_Com1Com2);
	// write a Generator to file 
	//g1.Write("parallel_g1L2.gen");
	// read generator
	//Generator parallel_g1L2("parallel_g1L2.gen");
	
	//Parallel(parallel_g1L2, line3, parallel_g1L3);
	// write a Generator to file 
	//g1.Write("parallel_g1L3.gen");
	// read generator
	//Generator parallel_g1L3("parallel_g1L3.gen");

    //Parallel(parallel_g1L2, line1, parallel_g1L1);
	parallel_Com1Com2.Write("ComCom.gen");
	

	// write result and operands for html docu
	//parallel_g1.Write("tmp_parallel_g1.gen");
	//parallel_g2.Write("tmp_parallel_g2.gen");
	//parallel_g1g2.Write("tmp_parallel_g1g2.gen");

	// inspect result on console

	std::cout << "################################\n";
    parallel_Com1Com2.Write();
    std::cout << "################################\n";
	//parallel_g1L1.Write();
	//std::cout << "################################\n";

	
	
	std::cout << "# tutorial, handcraft generator \n";
	// wasseem: transfer g1.write to variable // 
	std::stringstream ss;
	//change the underlying buffer and save the old buffer
	auto old_buf = std::cout.rdbuf(ss.rdbuf());
	// 
	parallel_Com1Com2.Write();

	// wasseem: transfer g1.write to variable ----> part II // 

	std::cout.rdbuf(old_buf); //reset

	std::cout << "<redirected-output>\n"
		<< ss.str()    // ss is the 
		<< "</redirected-output>" << std::endl;


	std::cout << "################################\n";

	// record test case

	//FAUDES_TEST_DUMP("generator1", generator1);
    
    
   
	
		
	
		/// take vriables to matlab workspace

		double  *pointer;          /* pointer to real data in new array */
		mwSize index;
		//const double data[] = { 2.1, 3.4, 2.3, 2.45 };  /* existing data */
		string data = ss.str().c_str();

		/* Check for proper number of arguments. */
		//if (nrhs != 0) {
		//	mexErrMsgIdAndTxt("MATLAB:arrayFillGetPr:rhs", "This function takes no input arguments.");
		//}

		/* Create an m-by-n mxArray; you will copy existing data into it */
		plhs[0] = mxCreateNumericMatrix(ROWS, COLUMNS, mxDOUBLE_CLASS, mxREAL);
		pointer = mxGetPr(plhs[0]);

		/* Copy data into the mxArray */

		for (index = 0; index < ELEMENTS; index++) {
			pointer[index] = data[index];
		}

		//pointer[] = data; 
		//mxFree(classIDflags);
        
        //FAUDES_TEST_DIFF()
	return;
}
