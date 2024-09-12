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

int *polar_encode_core(int *info_bits_padded, int array_size){
	
	int *coded_bits = malloc(array_size*sizeof(int));

	//printf("array_size = %d\n", array_size);

	if (array_size == 1){
		coded_bits[0] = info_bits_padded[0];
	}else{
		int *u1u2 = malloc(array_size/2*sizeof(int));
		int *u2 = malloc(array_size/2*sizeof(int));

		for (int i = 0; i < array_size/2; i++){
			u1u2[i] = (info_bits_padded[2*i]+info_bits_padded[2*i+1])%2;
			u2[i] = info_bits_padded[2*i+1];
		}
		
		int *lower_bits = polar_encode_core(u1u2, array_size/2);
		int *upper_bits = polar_encode_core(u2, array_size/2);

		memcpy(coded_bits, lower_bits, array_size/2*sizeof(int));
		memcpy(coded_bits+array_size/2, upper_bits, array_size/2*sizeof(int));

		free(u1u2);
		free(u2);
		free(lower_bits);
		free(upper_bits);
	}

	return coded_bits;
}

void print_polar_encode(struct polar_code *code, int *info_bits, int *info_bits_padded, int info_bits_count){
	
	printf("\nInfo bits:\n");
	for (int i = 0; i < code->info_len; i++){
		if (i==code->info_len-1 || (i+1)%32==0){
			printf("%d\n", info_bits[i]);
		}else{
			printf("%d ", info_bits[i]);
		}
	}

	printf("\nInfo bits after padding (info bits count = %d):\n", info_bits_count);
	for (int i = 0; i < code->block_len; i++){
		if (i==code->block_len-1 || (i+1)%32==0){
			printf("%d\n", info_bits_padded[i]);
		}else{
			printf("%d ", info_bits_padded[i]);
		}
	}

}

int *polar_encode(int *info_bits, struct polar_code *code, int *channel_order_sorted){
	
	// info bits padding
	int *info_bits_padded = malloc(code->block_len*sizeof(int));
	int info_bits_count = 0;
	for(int i = 0; i < code->block_len; i++){
		if(code->frozen_bits[i]) {
			info_bits_padded[i] = 0;		
		}else{
			info_bits_padded[channel_order_sorted[info_bits_count]] = info_bits[info_bits_count];
			info_bits_count++;
		}
	}

	int *coded_bits = polar_encode_core(info_bits_padded, code->block_len);

#ifdef ENABLE_DEBUG
	print_polar_encode(code, info_bits, info_bits_padded, info_bits_count);
#endif	
	free(info_bits_padded);
   	return coded_bits;

}
