#include "FFT.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
	const int N = 128;
	//declaring two arrays for real and imaginary
	double hw_X_R [N];
	double hw_X_I [N];
    int i;
	//Test inputs with imaginary parts given as zeros and real part as 
	// 0,1,2,........127
	//Inputs are type casted to double to get better precision which is found in matlab
    for (i=0; i < N; i++)
	{
        hw_X_R[i] = (double) i;
        hw_X_I[i] = (double) 0.0;
    }
	//Passing the inputs to FFT fucntion
    FFT(hw_X_R,hw_X_I);
  	
	// Printing the results
  	for (i = 0; i < N; ++i) 
	{
    	printf("(%0.4f, %0.4f)\n",hw_X_R[i], hw_X_I[i]);
  	}
   return 0;
}