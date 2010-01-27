#include "mex.h"

void trapezoidUnequal(const double *x, const double *y, const int len,
					 	double *out)
{
	int loop;
	double sum=0;

	*out = 0;
	for (loop = 0;loop < len-1; loop++)
	{
		sum += ((*(y+loop)+*(y+loop+1))/2)*(*(x+loop+1)-*(x+loop));
		*(out+loop+1) = sum;
	}
}

void trapezoidEqual(const double h, const double *y, const int len,
				   	double *out)
{
	int loop;
	double sum=0;

	*out = 0;
	for (loop = 0;loop < len-1; loop++)
	{
		sum += ((*(y+loop)+*(y+loop+1))/2)*h;
		*(out+loop+1) = sum;
	}
}

void simpson(const double h, const double *y, const int len,
			 	double *out)
{
	int loop = 2;
	double sum = 0;
	
	*out = 0;
	
	*(out+1) = (h/12)*(5**(y)+8**(y+1)-*(y+2));

	sum = *(out+1);
	for (loop = 2; loop < len-1; loop++)
	{
		sum += (h/12)*(5**(y+loop-1)+8**(y+loop)-*(y+loop+1));
		*(out+loop) = sum;
	}
	*(out+len-1) = *(out+len-2) + (h/12)*(5**(y+len-1)+8**(y+len-2)-*(y+len-3));
}

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{
	 /*Declare variables*/
	 int len, method;
	 char* methodStr;
	 double *inX, *inY, *integral;

	  /* Check for proper number of arguments. */
     if (nrhs != 3)
	 {
	 	mexErrMsgTxt("DiscreteInt(): Three inputs required");
	 }
	 if (nlhs != 1)
	 {
	 	mexErrMsgTxt("DiscreteInt(): Function only has one ouput");
	 }		 

	 /* Check for proper types*/
	 if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
	 {
		mexErrMsgTxt("DiscreteInt(): x must be type double");
	 }
	 if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
	 {
		mexErrMsgTxt("DiscreteInt(): y must be type double");
	 }
	 if (mxGetClassID(prhs[2]) != mxCHAR_CLASS)
	 {
		mexErrMsgTxt("DiscreteInt(): method must be a string");
	 }
	 
	 /* Parse the method string*/
	 methodStr = (char*)mxCalloc(mxGetN(prhs[2])+1,sizeof(char));
	 mxGetString(prhs[2],methodStr,mxGetN(prhs[2])+1);
	 if (strcmp(methodStr,"trapezoid") == 0)
	 {
	 	method = 0;
	 }
	 else if (strcmp(methodStr,"simpson") == 0)
	 {
	 	method = 1;
	 }
	 else
	 {
	 	mexErrMsgTxt("DiscreteInt(): Invalid method variable");
	 }
	 
	 /* Check for valid dimensions*/
	 if (mxGetM(prhs[0]) != 1)
	 {
	 	mexErrMsgTxt("DiscreteInt(): x must be a row vector");
	 }
	 if (mxGetM(prhs[1]) != 1)
	 {
	 	mexErrMsgTxt("DiscreteInt(): y must be a row vector");
	 }


	 /* Check for array lengths*/
	 if (method == 0 && mxGetN(prhs[0]) != 1 && mxGetN(prhs[0]) != mxGetN(prhs[1]))
	 {
	 	 mexErrMsgTxt("DiscreteInt(): Invalid x vector size");
	 }
	 else if (method == 1 && mxGetN(prhs[0]) != 1)
	 {
	 	 mexErrMsgTxt("DiscreteInt(): Simpson's rule requires a scalar step size");
	 }
			
	 /* Assign scalars*/
	 len = mxGetN(prhs[1]);

	 /* Creates output vectors*/
	 plhs[0] = mxCreateDoubleMatrix(1,len, mxREAL);
       
     /* Assign pointers to each input and output. */
     integral = mxGetPr(plhs[0]);
	 inX = mxGetPr(prhs[0]);
	 inY = mxGetPr(prhs[1]);

     /* Call the appropriate subroutine. */
     switch (method)
     {
     	case 0:		if (mxGetN(prhs[0]) > 1)
     					trapezoidUnequal(inX,inY,len,integral);
     				else
     					trapezoidEqual(*inX,inY,len,integral);
     				break;
     	case 1:		simpson(*inX,inY,len,integral);
     }	
}
