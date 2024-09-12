// SCL decoding of polar codes
// Author: Yaoyu Tao
// 05/12/2017
// implementation file

#include <stdlib.h>  
#include <stdio.h>  
#include <malloc.h>  
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "../header/polar_scl.h"

//#define ENABLE_DEBUG

double gaussian(double mean, double std){

	double u1, u2, w, mult, x1;
	
	do{
		u1 = -1 + ((double) rand ()/RAND_MAX)*2;
      		u2 = -1 + ((double) rand ()/RAND_MAX)*2;
      		w = pow (u1, 2) + pow (u2, 2);
	}while(w>=1 || w==0);

	mult = sqrt((-2*log(w))/w);
	x1 = u1*mult;

	return (mean + std*x1);
}

double *addNoise(int *encoded_bits, int block_len, int info_len, double ebno){
	
	double *received_llr = malloc(block_len*sizeof(double));
	double snr_scale = pow(10.0, ebno/20.0)*sqrt(((double) (info_len))/((double) (block_len)));
	double N0 = 1.0;
	double gaussian_num;	
	double received_signal;

#ifdef ENABLE_DEBUG	
	//FILE *fp;
	//fp = fopen("./gaussian.txt", "w");	
#endif

	for (int i = 0; i < block_len; i++){		

		gaussian_num = gaussian(0.0, N0);
#ifdef ENABLE_DEBUG	
		//fprintf(fp, "%.4f ", gaussian_num);
#endif
		received_signal = snr_scale*(2.0*encoded_bits[i]-1.0) + sqrt(N0/2.0)*gaussian_num;
		received_llr[i] = -4*snr_scale*received_signal/N0;
	}
#ifdef ENABLE_DEBUG
	//fclose(fp);
#endif

#ifdef ENABLE_DEBUG
	printf("\nReceived llr:\n");
	for (int i = 0; i < block_len; i++){
		if (i==block_len-1 || (i+1)%32==0){
			printf("%.3f\n", received_llr[i]);
		}else{
			printf("%.3f ", received_llr[i]);
		}
	}
#endif

	return received_llr;
}

