#include "FFT.h"
#include <stdio.h>

double W_real[64] = {
		1,0.998795456205172,0.995184726672197,0.989176509964781,0.980785280403230,0.970031253194544,0.956940335732209,0.941544065183021,0.923879532511287,0.903989293123443,0.881921264348355,0.857728610000272,0.831469612302545,0.803207531480645,0.773010453362737,0.740951125354959,0.707106781186548,0.671558954847018,0.634393284163646,0.595699304492434,0.555570233019602,0.514102744193222,0.471396736825998,0.427555093430282,0.382683432365090,0.336889853392220,0.290284677254462,0.242980179903264,0.195090322016128,0.146730474455362,0.0980171403295608,0.0490676743274181,6.12323399573677e-17,-0.0490676743274180,-0.0980171403295607,-0.146730474455362,-0.195090322016128,-0.242980179903264,-0.290284677254462,-0.336889853392220,-0.382683432365090,-0.427555093430282,-0.471396736825998,-0.514102744193222,-0.555570233019602,-0.595699304492433,-0.634393284163645,-0.671558954847018,-0.707106781186548,-0.740951125354959,-0.773010453362737,-0.803207531480645,-0.831469612302545,-0.857728610000272,-0.881921264348355,-0.903989293123443,-0.923879532511287,-0.941544065183021,-0.956940335732209,-0.970031253194544,-0.980785280403230,-0.989176509964781,-0.995184726672197,-0.998795456205172
};

double W_imag[64] = {
    0,-0.0490676743274180,-0.0980171403295606,-0.146730474455362,-0.195090322016128,-0.242980179903264,-0.290284677254462,-0.336889853392220,-0.382683432365090,-0.427555093430282,-0.471396736825998,-0.514102744193222,-0.555570233019602,-0.595699304492433,-0.634393284163646,-0.671558954847018,-0.707106781186548,-0.740951125354959,-0.773010453362737,-0.803207531480645,-0.831469612302545,-0.857728610000272,-0.881921264348355,-0.903989293123443,-0.923879532511287,-0.941544065183021,-0.956940335732209,-0.970031253194544,-0.980785280403230,-0.989176509964781,-0.995184726672197,-0.998795456205172,-1,-0.998795456205172,-0.995184726672197,-0.989176509964781,-0.980785280403230,-0.970031253194544,-0.956940335732209,-0.941544065183021,-0.923879532511287,-0.903989293123443,-0.881921264348355,-0.857728610000272,-0.831469612302545,-0.803207531480645,-0.773010453362737,-0.740951125354959,-0.707106781186548,-0.671558954847019,-0.634393284163646,-0.595699304492434,-0.555570233019602,-0.514102744193222,-0.471396736825998,-0.427555093430282,-0.382683432365090,-0.336889853392220,-0.290284677254462,-0.242980179903264,-0.195090322016129,-0.146730474455362,-0.0980171403295608,-0.0490676743274180
};

void butterfly(double X_real[2], double X_imag[2], double w_r, double w_img)
{
	double E1_real, E1_imag, E2_real, E2_imag;
	double temp1, temp2, temp3, temp4;

    // Storing the real and imaginary values in a temporary variables
	E1_real = X_real[0];
	E2_real = X_real[1];
	E1_imag = X_imag[0];
	E2_imag = X_imag[1];

    // Calculating the Twiddle factors for the butterfly structure
    // and storing it back to X_real and X_imag array
	temp3 = E1_real + E2_real;
	X_real[0] = temp3;
	temp4 = E1_imag + E2_imag;
	X_imag[0] = temp4;

    // Cross multiplication for computing
	temp1 = E1_real - E2_real;
	temp2 = E1_imag - E2_imag;
	X_real[1] = temp1*w_r - temp2*w_img;
    X_imag[1] = temp1*w_img + temp2*w_r;
}


// Main function to calculate the FFT which internatlly calls the butterfly function
void FFT(double X_R[128], double X_I[128])
{
    //N point FFT
	const int N = 128;
	int Head, length, temp1, i, j,k;
    double x_r[2],x_i[2];
    double temp_r[N];
    double temp_i[N];
    double w_r, w_i;
    int index, input_1,input_2;

    // i indicates number of stage minus 1
    // ------------or---------------
    // 2^x = N, then i = x-1
    FFT_label0:for(i = 0; i < 7; i++)
    {
        // temp is used to shift the i as a power of 2
		temp1 =	1<<i;
        //length is used to half the total length after every iteration making it radix structre
        length = N/temp1;
		// j indicates number of sub-operations for that specific stage
        FFT_label2: for(j = 0; j < temp1; j++)
        {
            //Head is to find the starting index of that particular iteration
            Head = length*j;
            {
                 FFT_label3: for(k = 0; k < (length/2); k++)
                 {
                    index = k*(1<<i);
                    //Exploring The symmetries of twiddle factors
	                if (index<64)
	                {
	                    w_r = W_real[index];
	                    w_i = W_imag[index];
	                }
                    else
                    {
                        w_r = -W_real[index-64];
                        w_i = -W_imag[index-64];
                    }
                    input_1 = k + Head;
            		input_2 = input_1+(length/2);
                    //pack two numbers and send for butterfly computation
                    x_r[0] = X_R[input_1];
                    x_r[1] = X_R[input_2];
                    x_i[0] = X_I[input_1];
                    x_i[1] = X_I[input_2];

                    butterfly(x_r,x_i,w_r,w_i);

                    X_R[input_1]= x_r[0] ;
                    X_R[input_2]= x_r[1] ;
                    X_I[input_1]= x_i[0] ;
                    X_I[input_2]= x_i[1] ;

            	}
            }
        }
	}
//Reorder the bits after FFT computation so as to display them in the correct order
   for(i = 0; i < N; i++ )
    {
    	index=(i&0x01) << 6 | (i&0x02) << 4 | (i&0x04)<< 2 | (i&0x08)<< 0| (i&0x10)>>2|(i&0x20)>>4|(i&0x40)>>6;
        temp_r[index] = X_R[i];
        temp_i[index] = X_I[i];
    }
	for(i = 0; i < N; i++)
    {
    	X_R[i] = temp_r[i];
    	X_I[i] = temp_i[i];
    }
}


int main(void)
{
	const int N = 128;
	double hw_X_R [N];
	double hw_X_I [N];
    int i;
    for (i=0; i < N; i++)
	{
        hw_X_R[i] = (double) i;
        hw_X_I[i] = (double) i+1;
    }
    FFT(hw_X_R,hw_X_I);
  	
  	for (i = 0; i < N; ++i) 
	{
    	printf("(%0.4f, %0.4f)\n",hw_X_R[i], hw_X_I[i]);
  	}
   return 0;
}
