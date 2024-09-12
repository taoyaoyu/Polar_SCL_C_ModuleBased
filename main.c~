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

#include "./header/polar_scl.h"

//#define ENABLE_DEBUG

int main(int argc, char* argv[]){
	
	// configure simulation
	long int num_runs = 1000000000000;
	long int num_max_err = atoi(argv[1]);
	double ebno = atof(argv[2]);
	int list_size = atoi(argv[3]);

	int enable_test_llr = 0;
	int test_info_bits[4] = {0, 0, 1, 0};		
	double test_llr[8] = {-4.365, 9.399, -3.252, 1.060, 12.112, -8.441, -3.234, 4.382};

	int integer_bits = atoi(argv[4]);
	int fractional_bits = atoi(argv[5]);
	int use_log_approx = 1;
	int use_quantization = 1;

	double max_value = pow(2.0, integer_bits-1) - pow(0.5, fractional_bits);
	double step = pow(0.5, fractional_bits);

	// configure polar code struct
	int n = 10; 					
	int block_len = pow(2,n);			
	int info_len = 1<<(n-1);			
	int crc_size = 0;			
	double epsilon = 0.32;			
	int *frozen_bits = malloc(block_len*sizeof(int));
	int *info_bits_index = malloc(info_len*sizeof(int));	
	int *bit_reverse_order = malloc(block_len*sizeof(int));	

#ifdef ENABLE_DEBUG
	printf("\nConfiguring polar code...\n");
	printf("\tn = %d\n\tblock_len = %d\t\n\tinfo_len = %d\n\tcrc_size = %d\n\tlist_size = %d\n", n, block_len, info_len, crc_size, list_size);
#endif

	struct polar_code code = { 
		.n = n,							// number of layers
		.block_len = block_len,			// block length
		.info_len = info_len,			// info bits length
		.crc_size = crc_size,			// crc coding size
		.epsilon = epsilon,			// epsilon for channel polarization
		.frozen_bits = frozen_bits,		// positions of frozen bits
		.info_bits_index = info_bits_index,
		.bit_reverse_order = bit_reverse_order,
		.list_size = list_size			// list size for list decoding
	};
	

	// initialize parameters
	int *info_bits = malloc(info_len*sizeof(int));
	int *channel_order_sorted;

	int *encoded_bits;
	double *received_llr;
	int *decoded_bits;

	if (enable_test_llr == 1){
		encoded_bits = malloc(block_len*sizeof(int));
		received_llr = malloc(block_len*sizeof(double));
		decoded_bits = malloc(info_len*sizeof(int));
	}

	// initialization
	channel_order_sorted = initialize_frozen_bits(&code);

	int frame_err_count = 0;
	int bit_err_count = 0;

	double fer = 0.0;
	double ber = 0.0;

	srand((unsigned) time(NULL));

	// start simulations
	long int sim = 0;
	for (sim = 0; sim < num_runs; sim++){

		// generate infomation bits to be send

		for (int i = 0; i < info_len; i++){
			if (enable_test_llr == 0){		
				info_bits[i] = rand() % 2;   // for simulations
			}else{
				info_bits[i] = test_info_bits[i];	// for verification 
			}
		}

		// encoder
		encoded_bits = polar_encode(info_bits, &code, channel_order_sorted);

		// for simulations
		if (enable_test_llr == 0){
			received_llr = addNoise(encoded_bits, block_len, info_len, ebno);
			//received_llr = malloc(block_len*sizeof(double));
			//for (int i = 0; i < block_len; i++){
			//	received_llr[i] = 0;
			//}
		}else{
			for (int i = 0; i < block_len; i++){
				received_llr[i] = test_llr[i];
			}
		}

		if (use_quantization == 1){
			for (int i = 0; i < block_len; i++){
				received_llr[i] = quantize(received_llr[i], step, max_value);
			}
		}

#ifdef ENABLE_DEBUG
		print_polar_debug(code, encoded_bits, received_llr);
#endif	

		decoded_bits = decode_scl(received_llr, &code, step, max_value, use_log_approx, use_quantization);

		//decoded_bits = malloc(info_len*sizeof(int));
		//for (int i = 0; i < info_len; i++){
		//	decoded_bits[i] = 0;
		//}
	
		if (decoded_bits == NULL){
			printf("\nError: decoded bits array is NULL!\n\n");
			exit(1);
		}else{
#ifdef ENABLE_DEBUG
			printf("\n\nPolar decisions:\n");
			for (int i = 0; i < code.info_len; i++){
				if (i==code.info_len-1 || (i+1)%32==0){
					printf("%d\n", decoded_bits[i]);
				}else{
					printf("%d ", decoded_bits[i]);
				}
			}
			printf("Send:\n");
			for (int i = 0; i < code.info_len; i++){
				if (i==code.info_len-1 || (i+1)%32==0){
					printf("%d\n", info_bits[i]);
				}else{
					printf("%d ", info_bits[i]);
				}
			}
#endif
		}

		int err_count = 0;
		for (int i = 0; i < info_len; i++){
			if(decoded_bits[i] != info_bits[i]){
				err_count++;
			}
		}
		if (err_count != 0){
			frame_err_count++;
			bit_err_count = bit_err_count + err_count;
		}

		if (!enable_test_llr && sim%1000 == 0){
			printf ("\nsimulation %d\tframe_err = %d\tbit_err = %d\n", sim, frame_err_count, bit_err_count);
		}
	
		if (frame_err_count >= num_max_err){
			break;
		}
		
		if (enable_test_llr == 0){
			free(encoded_bits);
			free(received_llr);
			free(decoded_bits);
		}
	}

	free(frozen_bits);
	free(info_bits_index);	
	free(bit_reverse_order);

	free(info_bits);
	free(channel_order_sorted);

	fer = ((double) frame_err_count)/((double) sim);
	ber = ((double) bit_err_count)/((double) sim*info_len);

	printf("\n\nTotal Simulation Count = %d\nFER = %e\nBER = %e\n\n", sim, fer, ber);

	return 0;

}

