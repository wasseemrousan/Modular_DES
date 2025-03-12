/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/


/// will try to emply this code for Libfaudes example.////


#include "mex.h"
#include <stdlib.h> 
#include <iostream>
#include <string>
#include <iterator>
#include "matrix.h"
#include <sstream>
#include "libfaudes.h"

// make faudes namespace available

using namespace std;
using std::string;

using namespace faudes;

#define ROWS 1

#define COLUMNS  30000

#define ELEMENTS 30000



#define MAXCHARS 80   /* max length of string contained in each field */



/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    char *states;              /* input states */
    char *events;               /*  input events */
    size_t ncols;                   /* size of matrix */
    char *outMatrix;              /* output matrix */


    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
   
    
    /* get the values of the input  */
     states = (char*)mxGetData(prhs[0]);
     events = (char*)mxGetData(prhs[1]);

     

    // at this point we need to convert char to string  in C
   
    // define strings in C from the inMatrix and multiplier varaibles 
     
    // define states
    char st1[6];
    st1[0] = states[0]; 
    st1[1] = states[2]; 
    st1[2] = states[4];
	st1[3] = states[6];
    st1[4] = states[8];	
    st1[5] = '\0';
    
    char st2[6];
    st2[0] = states[10]; 
    st2[1] = states[12]; 
    st2[2] = states[14]; 
	st2[3] = states[16];
    st2[4] = states[18];
    st2[5] = '\0';
    
//     char st3[5];
//     st3[0] = states[16]; 
//     st3[1] = states[18]; 
//     st3[2] = states[20];
// 	st3[3] = states[22];	
//     st3[4] = '\0';
    
    // from the same row of states - get the line id, which will be used in naming the line for .gen file and parallel operation
    //cid: component id 
    char cid[12];
    cid[0] = states[20]; 
    cid[1] = states[22]; 
    cid[2] = states[24];
    cid[3] = states[26];
    cid[4] = states[28];
	cid[5] = states[30];
    cid[6] = states[32];
    // also define '.gen' extension  in the same variable 
    cid[7] = '.';
    cid[8] = 'g';
    cid[9] = 'e';
    cid[10] = 'n';
    cid[11] = '\0';
    
    
    
    // define events 
    char ev1[5];
    ev1[0] = events[0]; 
    ev1[1] = events[2];
    ev1[2] = events[4];
    ev1[3] = events[6];
    ev1[4] = '\0';

    char ev2[5];
    ev2[0] = events[8]; 
    ev2[1] = events[10];
    ev2[2] = events[12];
    ev2[3] = events[14];
    ev2[4] = '\0';
    
    char ev3[5];
    ev3[0] = events[16];
    ev3[1] = events[18];
    ev3[2] = events[20];
    ev3[3] = events[22];
    ev3[4] = '\0';
    
//     char ev4[4];
//     ev4[0] = events[18]; 
//     ev4[1] = events[20];
//     ev4[2] = events[22];
//     ev4[3] = '\0';
      
    
    
    
    
    

    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = (char*)mxGetData(plhs[0]);


    
    
    ////////////////////////////////////
    //// DES part begins here //////////
    ////////////////////////////////////
    
    
    
    ///////////////////////////////////////////////

	// Constructors (part 1) and filling example///

	//////////////////////////////////////////////

	// at first we create an empty Generator object

	Generator g2;
	


	// do some random "user interaction" stuff with the Generator g1




	// g2 is for line 2 :between bus 1 and bus 3 

	g2.InsState(st1);
	g2.InsState(st2);
// 	g2.InsState(st3);
    

	 g2.InsEvent(ev1);
	 g2.InsEvent(ev2);
	 g2.InsEvent(ev3);
// 	 g2.InsEvent(ev4);
    
          

    
	g2.SetTransition( st1 , ev1 , st2 );
    g2.SetTransition(st1,ev2,st1);
// 	g2.SetTransition(st2,ev2,st3);
// 	g2.SetTransition(st3,ev3,st2);
	g2.SetTransition(st2,ev3,st1);

	g2.SetInitState(st1);
	g2.SetMarkedState(st1);



	
    
 
	// the autoamata will be wriiten into afile and all the Autmatas then will have parallel composition 
	// write a Generator to file 
	g2.Write(cid);
    
    
    // inspect result on console

	std::cout << "################################\n";

	std::cout << "# tutorial, handcraft generator \n";
    
    
    g2.Write(); 
    
	// all the code segemnts after that from MEX are not needed, since our goal is to genearte the automaton in .gen object only. 
    
	return;
}
