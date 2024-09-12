#include <stdlib.h>  
#include <stdio.h>  
#include <malloc.h>  
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "../header/polar_scl.h"

#define ENABLE_DEBUG


void print_initialization(struct polar_code *code, double *channels, double *channels_bitrev, double *channels_bitrev_sorted, int *channel_order_sorted){

	printf("\nChannel vector:\n");
	for (int i = 0; i < code->block_len; i++){
		if (i==code->block_len-1 || (i+1)%16==0){
			printf("%.3f\n", channels[i]);
		}else{
			printf("%.3f ", channels[i]);
		}
	}
	
	printf("\nBit reversal vector:\n");
	for (int i = 0; i < code->block_len; i++){
		if (i==code->block_len-1 || (i+1)%16==0){
			printf("%d\n", code->bit_reverse_order[i]);
		}else{
			printf("%d ", code->bit_reverse_order[i]);
		}
	}

	printf("\nChannel vector after bit reversal:\n");
	for (int i = 0; i < code->block_len; i++){
		if (i==code->block_len-1 || (i+1)%16==0){
			printf("%.3f\n", channels_bitrev[i]);
		}else{
			printf("%.3f ", channels_bitrev[i]);
		}
	}

	printf("\nSorted channel vector after bit reversal:\n");
	for (int i = 0; i < code->block_len; i++){
		if (i==code->block_len-1 || (i+1)%16==0){
			printf("%.3f\n", channels_bitrev_sorted[i]);
		}else{
			printf("%.3f ", channels_bitrev_sorted[i]);
		}
	}

	printf("\nChannel order sorted:\n");
	for (int i = 0; i < code->block_len; i++){
		if (i==code->block_len-1 || (i+1)%16==0){
			printf("%d\n", channel_order_sorted[i]);
		}else{
			printf("%d ", channel_order_sorted[i]);
		}
	}

	printf("\nFrozen bits: (info_len = %d)\n", code->info_len);
	int frozen_bits_count = 0;
	for (int i = 0; i < code->block_len; i++){
		if (code->frozen_bits[i] == 1){
			frozen_bits_count++;
		}
		if (i==code->block_len-1 || (i+1)%32==0){
			printf("%d\n", code->frozen_bits[i]);
		}else{
			printf("%d ", code->frozen_bits[i]);
		}
	}
	printf("Frozen bit count = %d\n", frozen_bits_count);
}


void print_scl_data_structure(int n, int block_len, int list_size, struct scl_data_structure *scl){
	
	printf("\n\n=========================SCL data structures:=========================\n");

	printf("\ninactive_path_indices (stack):\n");
	for (int i = 0; i < scl->inactive_path_indices->size; i++){
		if (i==scl->inactive_path_indices->size-1 || (i+1)%32==0){
			printf("%d\n", scl->inactive_path_indices->array[i]);
		}else{
			printf("%d ", scl->inactive_path_indices->array[i]);
		}
	}

	printf("\nactive_path_array:\n");
	for (int i = 0; i < list_size; i++){
		if (i==list_size-1 || (i+1)%32==0){
			printf("%d\n", scl->active_path_array[i]);
		}else{
			printf("%d ", scl->active_path_array[i]);
		}
	}

	printf("\npathIdx_to_arrayIdx:\n");
	for (int i = 0; i < (n+1); i++){
		printf("Tree level %d: 	", i);
		for (int j = 0; j < list_size; j++){
			printf("%d ", scl->pathIdx_to_arrayIdx[i][j]);
		}
		printf("\n");
	}

	printf("\ninactive_array_indices:\n");
	for (int i = 0; i < (n+1); i++){
		printf("Tree level %d: 	", i);
		for (int j = 0; j < scl->inactive_array_indices[i]->size; j++){
			printf("%d ", scl->inactive_array_indices[i]->array[j]);
		}
		printf("\n");
	}
	
	printf("\narray_ref_count:\n");
	for (int i = 0; i < (n+1); i++){
		printf("Tree level %d: 	", i);
		for (int j = 0; j < list_size; j++){
			printf("%d ", scl->array_ref_count[i][j]);
		}
		printf("\n");
	}
	
	printf("\nllr_scl:\n");
	for (int i = 0; i < list_size*(2*block_len-1); i++){
		if (i==list_size*(2*block_len-1)-1 || (i+1)%32==0){
			printf("%.3f\n", scl->llr_scl[i]);
		}else{
			printf("%.3f ", scl->llr_scl[i]);
		}
	}

	printf("\nllr_path_metric:\n");
	for (int i = 0; i < list_size; i++){
		if (i==list_size-1 || (i+1)%32==0){
			printf("%.3f\n", scl->llr_path_metric[i]);
		}else{
			printf("%.3f ", scl->llr_path_metric[i]);
		}
	}

	printf("\nc_scl:\n");
	for (int i = 0; i < list_size*(2*block_len-1); i++){
		for (int j = 0; j < 2; j++){
			printf("%d ", scl->c_scl[i][j]);
		}
		printf("\n");
	}

	printf("\ni_scl:\n");
	for (int i = 0; i < list_size; i++){
		for (int j = 0; j < block_len; j++){
			printf("%d ", scl->i_scl[i][j]);
		}
		printf("\n");
	}

	printf("\nlambda_offset:\n");
	for (int i = 0; i < (n+1); i++){
		if (i==(n+1)-1 || (i+1)%32==0){
			printf("%d\n", scl->lambda_offset[i]);
		}else{
			printf("%d ", scl->lambda_offset[i]);
		}
	}

	printf("\nlist_offset:\n");
	for (int i = 0; i < (list_size+1); i++){
		if (i==(list_size+1)-1 || (i+1)%32==0){
			printf("%d\n", scl->list_offset[i]);
		}else{
			printf("%d ", scl->list_offset[i]);
		}
	}
}


void print_polar_debug(struct polar_code code, int *encoded_bits, double *received_llr){
	
	printf("\nEncoded bits:\n");
	for (int i = 0; i < code.block_len; i++){
		if (i==code.block_len-1 || (i+1)%32==0){
			printf("%d\n", encoded_bits[i]);
		}else{
			printf("%d ", encoded_bits[i]);
		}
	}

	printf("\nReceived LLRs:\n");
	for (int i = 0; i < code.block_len; i++){
		if (i==code.block_len-1 || (i+1)%32==0){
			printf("%.3f\n", received_llr[i]);
		}else{
			printf("%.3f ", received_llr[i]);
		}
	}

	printf("\nHard decisions:");
	printf("\nRecv:\n");
	for (int i = 0; i < code.block_len; i++){
		if (i==code.block_len-1 || (i+1)%32==0){
			printf("%d\n", received_llr[i]<0);
		}else{
			printf("%d ", received_llr[i]<0);
		}
	}
	printf("Send:\n");
	for (int i = 0; i < code.block_len; i++){
		if (i==code.block_len-1 || (i+1)%32==0){
			printf("%d\n", encoded_bits[i]);
		}else{
			printf("%d ", encoded_bits[i]);
		}
	}
	printf("Miss:\n");
	int missed_bit_count = 0;
	for (int i = 0; i < code.block_len; i++){
		if (i==code.block_len-1 || (i+1)%32==0){
			printf("%d\n", !((received_llr[i]<0 && encoded_bits[i]==1) || (received_llr[i]>=0 && encoded_bits[i]==0)));
		}else{
			printf("%d ", !((received_llr[i]<0 && encoded_bits[i]==1) || (received_llr[i]>=0 && encoded_bits[i]==0)));
		}
		if(!((received_llr[i]<0 && encoded_bits[i]==1) || (received_llr[i]>=0 && encoded_bits[i]==0))){
			missed_bit_count++;
		}
	}
	printf("Hard decision errors: %d\n", missed_bit_count);
}

void print_unfrozen_bits_calculation(int rho, int **cont_forks, double **prob_forks, double *prob, int list_size){
	
	// empty
	printf("\nrho = %d\n", rho);

	printf("\nprob_forks:\n");
	for (int i = 0; i < list_size; i++){
		printf("%.3f %.3f\n", prob_forks[i][0], prob_forks[i][1]);
	}

	printf("\nprob:\n");
	for (int i = 0; i < 2*list_size; i++){
		printf("%.3f ", prob[i]);
	}
	printf("\n");
	printf("\ncont_forks:\n");
	for (int i = 0; i < list_size; i++){
		printf("%d %d\n", cont_forks[i][0], cont_forks[i][1]);
	}
	printf("\n\n");
}

