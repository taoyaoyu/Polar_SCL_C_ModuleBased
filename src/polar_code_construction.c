#include <stdlib.h>  
#include <stdio.h>  
#include <malloc.h>  
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "../header/polar_scl.h"

//#define ENABLE_DEBUG

// functions
void bit_reverse_order(struct polar_code *code) {
	//printf("\nBit reversal start...\n");
	for (uint16_t i = 0; i < code->block_len; ++i) {
		uint16_t to_be_reversed = i;
		//printf("to_be_reversed = %d ", to_be_reversed);
		code->bit_reverse_order[i] = (uint16_t) ((to_be_reversed & 1) << (code->n-1));
		for (uint8_t j = (uint8_t) (code->n-1); j>0; --j){
			to_be_reversed >>= 1;
			code->bit_reverse_order[i] += (to_be_reversed & 1) << (j-1);
		}
		//printf(" reversal_result = %d\n", code->bit_reverse_order[i]);
	}
}

int* initialize_frozen_bits (struct polar_code *code) {
	
	// calculate channel polarization
	double *channels = malloc(code->block_len*sizeof(double));
	for (int i = 0; i < code->block_len; i++){
		channels[i] = code->epsilon;
	}
	for (int i = 0; i < code->n; i++){
		double *c1 = malloc(code->block_len/2*sizeof(double));
		double *c2 = malloc(code->block_len/2*sizeof(double));
		for (int j = 0; j < code->block_len; j+=2){
			c1[j/2] = channels[j];
			c2[j/2] = channels[j+1];
		}
		for (int j = 0; j < code->block_len; j++){
			if (j < code->block_len/2){
				channels[j] = c1[j]+c2[j]-c1[j]*c2[j];
			}else{
				channels[j] = c1[j-code->block_len/2]*c2[j-code->block_len/2];
			}
		}
		free(c1);
		free(c2);
	}
	
	bit_reverse_order(code);

	// apply bit reversal
	double *channels_bitrev = malloc(code->block_len*sizeof(double));
	for (int i = 0; i < code->block_len; i++){
		channels_bitrev[i] = channels[code->bit_reverse_order[i]];
	}

	double *channels_bitrev_sorted = malloc(code->block_len*sizeof(double));
	for (int i = 0; i < code->block_len; i++){
		channels_bitrev_sorted[i] = channels_bitrev[i];
	}
	qsort(channels_bitrev_sorted, code->block_len, sizeof(double), compare_double);

	int *channel_order_sorted = malloc(code->block_len*sizeof(int));
	for (int i = 0; i < code->block_len; i++){
		for (int j = 0; j < code->block_len; j++){
			if(channels_bitrev_sorted[i] == channels_bitrev[j]){
				int exist = 0;
				for (int k = 0; k < i; k++){
					if (channel_order_sorted[k] == j){
						exist = 1;
						break;
					}
				}
				if(exist != 1){
					channel_order_sorted[i] = j;
					break;
				}
			}
		}
	}

 	for (int i = 0; i < code->info_len; i++) {
        	code->frozen_bits[channel_order_sorted[i]] = 0;
		code->info_bits_index[i] = channel_order_sorted[i];
    	}
    	for (int i = code->info_len; i < code->block_len; i++) {
       		code->frozen_bits[channel_order_sorted[i]] = 1;
    	}

#ifdef ENABLE_DEBUG
	print_initialization(code, channels, channels_bitrev, channels_bitrev_sorted, channel_order_sorted);
#endif	

	free(channels);
	free(channels_bitrev);
	free(channels_bitrev_sorted);

	return 	channel_order_sorted;

}

