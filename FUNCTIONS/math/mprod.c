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

/* The computational routine */
void LegendrePolynomialm(double *A, int a, double *B, int b, double *C, int colA, int rowA, int colB, int rowB)
{
    int i, j, k, ii, jj; double sum = 0;
    
    // A and B full
    if (a==1 & b==1){
        for (i=0; i<rowA; i++){
            for(j=0; j<colB; j++){
                for (k = 0; k<rowB; k++){
                    sum = sum + *(A + colA*i + k) * *(B + colB*k + j);
                }
                *(C + rowB*i + j) = sum;
                sum = 0;
            }
        }
    }
    // A and B diagonal
    if (a==0 && b==0){
        for (i=0; i<rowA; i++){
            C[rowA*i + i] = A[rowA*i + i] * B[rowA*i + i];
        }
    }
    // A full and B diagonal
    if (a==1 & b==0){
        for (i=0; i<rowA; i++){
            for(j=0; j<colB; j++){
                ii = rowA*i;
                C[ii + j] = A[ii + j] * B[ii + i];
            }
        }
    }
    // A diagonal and B full
    if (a==0 & b==1){
        for (i=0; i<rowA; i++){
            for(j=0; j<colB; j++){
                ii = rowA*i; jj = rowA*j;
                C[ii + j] = A[jj + j] * B[ii + j];
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
//     int multiplier;              /* input scalar */
    double *inMatrixA;           /* 1xN input matrix 1 */
    int     typea;               /* type of a */
    double *inMatrixB;           /* 1xN input matrix 2 */
    int     typeb;               /* type of a */
    size_t cola;                /* size of matrix 1*/
    size_t rowa;                /* size of matrix 1*/
    size_t colb;                /* size of matrix 1*/
    size_t rowb;                /* size of matrix 1*/
    double *outMatrix;            /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:LegendrePolynomialm:nrhs","Four inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:LegendrePolynomialm:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
//     if( !mxIsDouble(prhs[0]) || 
//          mxIsComplex(prhs[0]) ||
//          mxGetNumberOfElements(prhs[0])!=1 ) {
//         mexErrMsgIdAndTxt("MyToolbox:LegendrePolynomialm:notScalar","Input multiplier must be a scalar.");
//     }
    
    /* make sure the second input argument is type double */
//     if( !mxIsDouble(prhs[1]) || 
//          mxIsComplex(prhs[1])) {
//         mexErrMsgIdAndTxt("MyToolbox:LegendrePolynomialm:notDouble","Input matrix must be type double.");
//     }
    
    /* check that number of rows in second input argument is 1 */
//     if(mxGetM(prhs[1])!=1) {
//         mexErrMsgIdAndTxt("MyToolbox:LegendrePolynomialm:notRowVector","Input must be a row vector.");
//     }
    

    /* create a pointer to the real data in the input matrix 1 */
    inMatrixA = mxGetPr(prhs[0]);

    /* get the type of the matrix 1  */
    typea = mxGetScalar(prhs[1]);
    
    /* create a pointer to the real data in the input matrix 2 */
    inMatrixB = mxGetPr(prhs[2]);

    /* get the type of the matrix 1  */
    typeb = mxGetScalar(prhs[3]);
    
    
    /* get dimensions of the input matrix  1*/
    cola = mxGetN(prhs[0]);
    rowa = mxGetM(prhs[0]);
    colb = mxGetN(prhs[2]);
    rowb = mxGetM(prhs[2]);
    

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(cola,rowb,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    LegendrePolynomialm(inMatrixA,typea,inMatrixB,typeb,outMatrix,(mwSize)cola,(mwSize)rowa,(mwSize)colb,(mwSize)rowb);
}
