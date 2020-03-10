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

#include "mex.h"
#include "LegendreLobattoAndDerivative.c"
#include "math.h"

#define PI 3.14159265358979323846
/* The computational routine */
void LegendreGausLobattoNodesAndWeights(int N, double *xi, double *w)
{
    int i, m, k, nit = 100;
    static double q, qp, LN;
    double delta;
    
    if (N==1){
        xi[0] = -1;
        xi[1] =  1;
        w[0]  =  1;
        w[1]  = w[0];
    }
    else{
        
        /* Setting up the extremities */
        xi[0] = -1;
        xi[N] =  1;
        w[0]  = 2/((double)N*((double)N+1));
        w[N]  = w[0];
        
        /* Setting up the middle */
        LegendreLobattoAndDerivative(N, 0, &q, &qp, &LN);
        
        m = (int)ceil((N+1)/2);
        xi[m] = 0;
        w[m]  = 2/(N*(N+1)*pow(LN,2));
        
        /* Looping through the length of x*/
        for (i=1; i<=(floor((N+1)/2) - 1); i++){
            xi[i] = - cos((i + 1/4)*PI/N - 3/(8*N*PI*(i + 1/4)));
            
            /* Newton iteration */
            for (k=1; k<=nit; k++){
            LegendreLobattoAndDerivative(N, xi[i], &q, &qp, &LN);
            
            delta = -q/qp;
            
            xi[i] = xi[i] + delta;
            if(fabs(delta) < fabs(1e-12*xi[i]))  break;
            }
            
            
            w[i] = 2/(N*(N+1)*pow(LN,2));
            
            xi[N-i] = -xi[i];
            w[N-i]  =  w[i];
        }
    }
    
    
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    int multiplier;              /* input scalar */
    double *outMatrixxi;              /* output matrix */
    double *outMatrixw;              /* output matrix */
    size_t ncols;                   /* size of matrix */

    /* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:LegendreGausLobattoNodesAndWeights:nrhs","One input required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:LegendreGausLobattoNodesAndWeights:nlhs","Two outputs required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:LegendreGausLobattoNodesAndWeights:notScalar","Input multiplier must be a scalar.");
    }
    
    /* get the value of the scalar input  */
    multiplier = mxGetScalar(prhs[0]);
    ncols      = multiplier + 1;


    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize)ncols,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)ncols,1,mxREAL);

    /* get a pointer to the real data in the output matrices */
    outMatrixxi  = mxGetPr(plhs[0]);
    outMatrixw   = mxGetPr(plhs[1]);

    /* call the computational routine */
    LegendreGausLobattoNodesAndWeights(multiplier,outMatrixxi,outMatrixw);
}
