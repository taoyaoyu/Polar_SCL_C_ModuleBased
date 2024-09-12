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

int *decode_scl(double *llr, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization){

	int list_size = code->list_size;
	int n = code->n;
	int block_len = code->block_len;

#ifdef ENABLE_DEBUG
	printf("\nReceived llr's:\n");
	for (int i = 0; i < block_len; i++){
		if (i==block_len-1 || (i+1)%32==0){
			printf("%.3f\n", llr[i]);
		}else{
			printf("%.3f ", llr[i]);
		}
	}
#endif

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initial scl data structures
        //////////////////////////////////////////////////////////////////////////////////////////////////////

	struct stack inactive_path_indices = {.array = NULL, .size = 0};	// stack
	int *active_path_array = malloc(list_size*sizeof(int));
	int **pathIdx_to_arrayIdx = malloc((n+1)*sizeof(int *));
	struct stack **inactive_array_indices = malloc((n+1)*sizeof(struct stack*));	// each element is stack
	int **array_ref_count = malloc((n+1)*sizeof(int *));
	double *llr_scl = malloc(list_size*(2*block_len-1)*sizeof(double));
	double *llr_path_metric = malloc(list_size*sizeof(double));
	int **c_scl = malloc(list_size*(2*block_len-1)*sizeof(int *));
	int **i_scl = malloc(list_size*sizeof(int *));

	int *lambda_offset = malloc((n+1)*sizeof(int));
	int *list_offset = malloc((list_size+1)*sizeof(int));

	struct scl_data_structure scl = {
		.inactive_path_indices = &inactive_path_indices,	// stack
		.active_path_array = active_path_array,
		.pathIdx_to_arrayIdx = pathIdx_to_arrayIdx,
		.inactive_array_indices = inactive_array_indices,	// each element is stack
		.array_ref_count = array_ref_count,
		.llr_scl = llr_scl,
		.llr_path_metric = llr_path_metric,
		.c_scl = c_scl,
		.i_scl = i_scl,

		.lambda_offset = lambda_offset,
		.list_offset = list_offset
	};

	init_data_structure_scl(n, block_len, list_size, &scl);

#ifdef ENABLE_DEBUG
	printf("\n\n\nAfter Init:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif	

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Assign initial path
        //////////////////////////////////////////////////////////////////////////////////////////////////////
	
	if (scl.inactive_path_indices->size == 0){
		printf("\nAssign initial path inactive_path_indices size == 0\n");
		exit(1);
	}

	int  l_index = top(scl.inactive_path_indices);
    	pop(scl.inactive_path_indices);
    	scl.active_path_array[l_index] = 1;
	
#ifdef ENABLE_DEBUG
	printf("\nl_index = %d\n", l_index);
#endif

	// Associate arrays with path index
    	for (int lambda = 0; lambda < n + 1; lambda++) {

		if (scl.inactive_array_indices[lambda]->size == 0){
			printf("\nAssign initial path inactive_array_indices size == 0\n");
			exit(1);
		}

        	int s = top(scl.inactive_array_indices[lambda]);
        	pop(scl.inactive_array_indices[lambda]);
        	scl.pathIdx_to_arrayIdx[lambda][l_index] = s;
        	scl.array_ref_count[lambda][s] = 1;
    	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Get array pointer P
        //////////////////////////////////////////////////////////////////////////////////////////////////////

    	int s_index = get_array_ptr_p(0, l_index, &scl, code); 

#ifdef ENABLE_DEBUG
	printf("\ns_index = %d\n", s_index);
#endif
	
	int llr_count = 0;
	for (int i = get_i_scl(0, 0, s_index, scl.lambda_offset, scl.list_offset); i < get_i_scl(0, block_len-1, s_index, scl.lambda_offset, scl.list_offset)+1; i++){
		scl.llr_scl[i] = llr[llr_count];
		llr_count++;
	}


#ifdef ENABLE_DEBUG
	printf("\n\n\nAfter assigning init path and get_array_ptr_p and llr_scl:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif	

    	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// Start SCL decoding process
        //////////////////////////////////////////////////////////////////////////////////////////////////////
	
	for (int phi = 0; phi < block_len; phi++){
		

#ifdef ENABLE_DEBUG
	printf("\n\n================= Bit index = %d ==================\n\n", phi);
#endif	


		recursively_calc_p_scl(n, phi, &scl, code, step, max_value, use_log_approx, use_quantization);

#ifdef ENABLE_DEBUG
	printf("\n\nAfter recursively_calc_p_scl:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif	
	
		if (code->frozen_bits[phi] == 1){
			continue_paths_frozen_bit(phi, &scl, code, step, max_value, use_log_approx, use_quantization);
#ifdef ENABLE_DEBUG
	printf("\n\nAfter continue_paths_frozen_bit:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif	

		}else{
			continue_paths_unfrozen_bit(phi, &scl, code, step, max_value, use_log_approx, use_quantization);
#ifdef ENABLE_DEBUG
	printf("\n\nAfter continue_paths_unfrozen_bit:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif	
		}

		if (phi%2 == 1){
			recursively_update_c_scl(n, phi, &scl, code);
		}

#ifdef ENABLE_DEBUG
	printf("\n\nAfter recursively_update_c_scl:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif

		// update llr_metric difference for quantization case
		if (use_quantization){
			double min_llr_path_metric = min_array(scl.llr_path_metric, list_size);
			for (int i = 0; i < list_size; i++){
				scl.llr_path_metric[i] = quantize(scl.llr_path_metric[i]-min_llr_path_metric, step, max_value);
			}
		}

#ifdef ENABLE_DEBUG
	printf("\n\nAfter llr_metric difference for quantization case:");
	print_scl_data_structure(n, block_len, list_size, &scl);
#endif
	}

	l_index = find_most_probable_path(&scl, code);
#ifdef ENABLE_DEBUG
	printf("\nfind most probable path l_index = %d\n", l_index);
#endif
	int c_m = get_array_ptr_p(n, l_index, &scl, code);
	int *info = malloc(block_len*sizeof(int));
	for (int i = 0; i < block_len; i++){
		info[i] = scl.i_scl[c_m][i];
	} 

	int *u = malloc(code->info_len*sizeof(int));
	for (int i = 0; i < code->info_len; i++){
		u[i] = info[code->info_bits_index[i]];
	}  

	free_scl_data_structure(n, block_len, list_size, &scl);
	free(info);
	
	return u;
}

int get_array_ptr_p(int lambda, int l_index, struct scl_data_structure *scl, struct polar_code *code){
	
	//printf("\n\nEntering get_array_ptr_p...\n");

	int n = code->n;
	int list_size = code->list_size;
	int s = scl->pathIdx_to_arrayIdx[lambda][l_index];
    	int s_p;

#ifdef ENABLE_DEBUG
	int *i_s_p;
	int *i_s;
#endif	

    	if (scl->array_ref_count[lambda][s] == 1) {
        	s_p = s;
    	}else {
		
		if (scl->inactive_array_indices[lambda]->size == 0){
			printf("\nget_array_ptr_p inactive_array_indices size == 0\n");
			exit(1);
		}

        	s_p = top(scl->inactive_array_indices[lambda]);
        	pop(scl->inactive_array_indices[lambda]);
		
#ifdef ENABLE_DEBUG	
        	i_s_p = malloc(pow(2, n-lambda)*sizeof(int));
		for (int i = 0; i < pow(2, n-lambda); i++){
			i_s_p[i] = scl->lambda_offset[lambda] + scl->list_offset[s_p] + i;
		}

		printf("\ns = %d\n", s);
		printf("\ni_s_p:\n");
		for (int i = 0; i < pow(2, n-lambda); i++){
			if (i==pow(2, n-lambda)-1 || (i+1)%32==0){
				printf("%d\n", i_s_p[i]);
			}else{
				printf("%d ", i_s_p[i]);
			}
		}	

		i_s = malloc(pow(2, n-lambda)*sizeof(int));
		for (int i = 0; i < pow(2, n-lambda); i++){
			i_s[i] = scl->lambda_offset[lambda] + scl->list_offset[s] + i;
		}


		printf("\ni_s:\n");
		for (int i = 0; i < pow(2, n-lambda); i++){
			if (i==pow(2, n-lambda)-1 || (i+1)%32==0){
				printf("%d\n", i_s[i]);
			}else{
				printf("%d ", i_s[i]);
			}
		}
#endif	
		

		for (int i = 0; i < pow(2, n-lambda); i++){
			scl->c_scl[scl->lambda_offset[lambda] + scl->list_offset[s_p] + i][0] = scl->c_scl[scl->lambda_offset[lambda] + scl->list_offset[s] + i][0];
			scl->c_scl[scl->lambda_offset[lambda] + scl->list_offset[s_p] + i][1] = scl->c_scl[scl->lambda_offset[lambda] + scl->list_offset[s] + i][1];
			scl->llr_scl[scl->lambda_offset[lambda] + scl->list_offset[s_p] + i] = scl->llr_scl[scl->lambda_offset[lambda] + scl->list_offset[s] + i];
		}
		
        	scl->array_ref_count[lambda][s]--;
        	scl->array_ref_count[lambda][s_p] = 1;
        	scl->pathIdx_to_arrayIdx[lambda][l_index] = s_p;

    	}

    	return s_p; 
}


void recursively_calc_p_scl (int lambda, int phi, struct scl_data_structure *scl, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization) {
	
	int n = code->n;
	int list_size = code->list_size;
	
	if (lambda == 0){
		return;
	}

	int psi = phi>>1;

	if (phi%2 == 0){
		recursively_calc_p_scl(lambda-1, psi, scl, code, step, max_value, use_log_approx, use_quantization);
	}

	int *p_index_3_base_list = malloc(list_size*sizeof(int));
	int *l_index_1_list = malloc(list_size*sizeof(int));

	for (int i = 0; i < list_size; i++){
		p_index_3_base_list[i] = 0;
		l_index_1_list[i] = 0;
	}

	for (int l_index = 0; l_index < list_size; l_index++){
		
		if (scl->active_path_array[l_index] == 0){
			continue;
		}

		l_index_1_list[l_index] = get_array_ptr_p(lambda, l_index, scl, code);
                int l_index_2 = get_array_ptr_p(lambda-1, l_index, scl, code);
                int l_index_3 = get_array_ptr_p(lambda, l_index, scl, code);
                
                int p_index_1_base = scl->lambda_offset[lambda-1] + scl->list_offset[l_index_2];
                p_index_3_base_list[l_index] = scl->lambda_offset[lambda] + scl->list_offset[l_index_1_list[l_index]];
                int c_index_3_base = scl->lambda_offset[lambda] + scl->list_offset[l_index_3];

		for (int beta = 0; beta < pow(2, n-lambda); beta++){
			
			int p_index_1 = p_index_1_base + 2 * beta;
                    	int p_index_2 = p_index_1_base + 2 * beta + 1;
                    	int p_index_3 = p_index_3_base_list[l_index] + beta;

			if (phi%2 == 0){
				
				if (max(abs(scl->llr_scl[p_index_1]), abs(scl->llr_scl[p_index_2])) < 40){
					if (use_log_approx == 1){
						scl->llr_scl[p_index_3] = f_func(scl->llr_scl[p_index_1], scl->llr_scl[p_index_2], use_quantization, step, max_value);
					}else{
						scl->llr_scl[p_index_3] = log( (exp(scl->llr_scl[p_index_1] + scl->llr_scl[p_index_2]) + 1) / (exp(scl->llr_scl[p_index_1])  + exp(scl->llr_scl[p_index_2]) ) );
					}
				}else{
					scl->llr_scl[p_index_3] = sign( scl->llr_scl[p_index_1]) * sign(scl->llr_scl[p_index_2]) * min(abs(scl->llr_scl[p_index_2]), abs(scl->llr_scl[p_index_1]));
					if (use_quantization==1){
						scl->llr_scl[p_index_3] = quantize(scl->llr_scl[p_index_3], step, max_value);
					}
				}

			} else {
				int u_p = scl->c_scl[c_index_3_base + beta][0];
				if (use_log_approx == 1){
					scl->llr_scl[p_index_3] = g_func(scl->llr_scl[p_index_1], scl->llr_scl[p_index_2], u_p, use_quantization, step, max_value);
				}else{
					scl->llr_scl[p_index_3] = pow(-1,u_p) * scl->llr_scl[p_index_1] + scl->llr_scl[p_index_2];
					if (use_quantization==1){
						scl->llr_scl[p_index_3] = quantize(scl->llr_scl[p_index_3], step, max_value);
					}
				}
			}
		}
	}

	free(p_index_3_base_list);
	free(l_index_1_list);
}

void recursively_update_c_scl(int lambda, int phi, struct scl_data_structure *scl, struct polar_code *code){

	if (phi%2 == 0){
		printf("Error: phi should always be odd in recursively_update_c_scl function call\n");
		exit(1);
	}
	
	int list_size = code->list_size;
	int n = code->n;
	int psi = floor(phi/2);

	for (int l_index = 0; l_index < list_size; l_index++){

		if (scl->active_path_array[l_index] == 0){	// if not active path, ignore and continue
			continue;
		}
		
		int l_index_1 = get_array_ptr_p(lambda, l_index, scl, code);
		int l_index_2 = get_array_ptr_p(lambda-1, l_index, scl, code);

		int p_index_1 = scl->lambda_offset[lambda] + scl->list_offset[l_index_1];
		int p_index_2 = scl->lambda_offset[lambda-1] + scl->list_offset[l_index_2];

		for (int beta = 0; beta < pow(2, n-lambda); beta++){
			scl->c_scl[p_index_2 + 2*beta][psi%2] = (scl->c_scl[p_index_1+beta][0] + scl->c_scl[p_index_1+beta][1])%2;
			scl->c_scl[p_index_2 + 2*beta + 1][psi%2] = scl->c_scl[p_index_1+beta][1];
		}
	}

	if (psi%2 == 1){
		recursively_update_c_scl(lambda-1, psi, scl, code);
	}

}

void continue_paths_frozen_bit(int phi, struct scl_data_structure *scl, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization ){

	int n = code->n;
	int list_size = code->list_size;

	for (int l_index = 0; l_index < list_size; l_index++){
		if (scl->active_path_array[l_index] == 0){
			continue;
		}

		int l_index_1 = get_array_ptr_p(n, l_index, scl, code);
		scl->c_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]][phi%2] = 0;

		if (use_log_approx){
			if (use_quantization) {
				scl->llr_path_metric[l_index_1] = quantize(scl->llr_path_metric[l_index_1]-scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]], step, max_value);
			}else{
				scl->llr_path_metric[l_index_1] = scl->llr_path_metric[l_index_1]-scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]];
			}
		}else{
			scl->llr_path_metric[l_index_1] = scl->llr_path_metric[l_index_1] + log(1+exp(-scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]]));
		}
	}
}

int clone_path(int l_index, int n, struct scl_data_structure *scl){

	if (scl->inactive_path_indices->size == 0){
		printf("\nclone_path inactive_path_indices size == 0\n");
		exit(1);
	}

	int l_p_index = top(scl->inactive_path_indices);
	pop(scl->inactive_path_indices);
	scl->active_path_array[l_p_index] = 1;

	scl->llr_path_metric[l_p_index] = scl->llr_path_metric[l_index];

	for (int lambda = 0; lambda < n+1; lambda++){
		int s = scl->pathIdx_to_arrayIdx[lambda][l_index];
		scl->pathIdx_to_arrayIdx[lambda][l_p_index] = s;
		scl->array_ref_count[lambda][s]++;
	}

	return l_p_index;
}

void kill_path(int l_index, int n, struct scl_data_structure *scl){

	scl->active_path_array[l_index] = 0;
	push(scl->inactive_path_indices, l_index);

	scl->llr_path_metric[l_index] = 0;

	for (int lambda = 0; lambda < n+1; lambda++){
		
		int s = scl->pathIdx_to_arrayIdx[lambda][l_index];
		scl->array_ref_count[lambda][s]--;
		if (scl->array_ref_count[lambda][s] == 0){
			push(scl->inactive_array_indices[lambda], s);	
		}
	}

}

void continue_paths_unfrozen_bit(int phi, struct scl_data_structure *scl, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization ){

	int n = code->n;
	int list_size = code->list_size;
	int block_len = code->block_len;

	double **prob_forks = malloc(list_size*sizeof(double *));
	int **cont_forks = malloc(list_size*sizeof(int *));
	double *prob = malloc(2*list_size*sizeof(double));

	for (int i = 0; i < list_size; i++){
		prob_forks[i] = malloc(2*sizeof(double));
		prob_forks[i][0] = -REALMAX;
		prob_forks[i][1] = -REALMAX;
	
		cont_forks[i] = malloc(2*sizeof(int));
		cont_forks[i][0] = 0;
		cont_forks[i][1] = 0;
	}
			
	int index = 0;

	for (int l_index = 0; l_index < list_size; l_index++){

		if (scl->active_path_array[l_index] == 1){

			int l_index_1 = get_array_ptr_p(n, l_index, scl, code);

			if (use_log_approx == 1){

				double temp1 = -scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]];	
				double temp2 = scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]];

				if (use_quantization){
					prob_forks[l_index][0] = quantize(-(scl->llr_path_metric[l_index_1]+temp1), step, max_value);
					prob_forks[l_index][1] = quantize(-(scl->llr_path_metric[l_index_1]+temp2), step, max_value);
				}else{
					prob_forks[l_index][0] = -(scl->llr_path_metric[l_index_1]+temp1);
					prob_forks[l_index][1] = -(scl->llr_path_metric[l_index_1]+temp2);
				}
			
			}else{
				
				double temp1 = log(1 + exp(-scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]]));
                        	double temp2 = log(1 + exp(scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]]));
                        	prob_forks[l_index][0] = - (scl->llr_path_metric[l_index_1] + temp1);  
                        	prob_forks[l_index][1] =  - ( scl->llr_path_metric[l_index_1] + temp2); 
			}
			
			index++;
		}
	}

#ifdef ENABLE_DEBUG
	printf("\n\nAfter first half of continue_paths_unfrozen_bit:");
	print_scl_data_structure(n, block_len, list_size, scl);

	printf("\nprob_forks:\n");
	for (int i = 0; i < list_size; i++){
		printf("%.3f %.3f\n", prob_forks[i][0], prob_forks[i][1]);
	}
#endif
	
	// pick L highest prob leaf nodes in forks
	int rho = ((int) min(((double) 2*index), ((double) list_size)));
	
	for (int i = 0; i < 2*list_size; i++){
		prob[i] = prob_forks[i%list_size][((int) i/list_size)];
	}
	qsort(prob, 2*list_size, sizeof(double), compare_double_descend);
		
	double threshold = prob[rho];
	int num_populated = 0;

	for (int l_index = 0; l_index < list_size; l_index++){
		for (int j_index = 0; j_index < 2; j_index++){
			if (num_populated == rho){
				break;
			}

			if (prob_forks[l_index][j_index] > threshold){
				cont_forks[l_index][j_index] = 1;
				num_populated++;
			}		
		}
	}		
	
	if (num_populated < rho){
		for (int l_index = 0; l_index < list_size; l_index++){		
			for (int j_index = 0; j_index < 2; j_index++){
				if (num_populated == rho){
					break;
				}
				if (prob_forks[l_index][j_index] == threshold){
					cont_forks[l_index][j_index] = 1;
					num_populated++;
				}
			}
		}
	}


	for (int l_index = 0; l_index < list_size; l_index++){	

		if (scl->active_path_array[l_index] == 0){
			continue;
		}

		if(cont_forks[l_index][0] == 0 && cont_forks[l_index][1] == 0) {
			kill_path(l_index, n, scl);
		}
	}


#ifdef ENABLE_DEBUG
	printf("\n\nAfter kill path of continue_paths_unfrozen_bit:");
	print_scl_data_structure(n, block_len, list_size, scl);
	print_unfrozen_bits_calculation(rho, cont_forks, prob_forks, prob, list_size);
#endif

	for (int l_index = 0; l_index < list_size; l_index++){		
		
		if (cont_forks[l_index][0] == 0 && cont_forks[l_index][1] == 0){
			continue;
		}

		int l_index_1 = get_array_ptr_p(n, l_index, scl, code);
		if (cont_forks[l_index][0] == 1 && cont_forks[l_index][1] == 1) {
			
			scl->c_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]][phi%2] = 0;
			scl->i_scl[l_index_1][phi] = 0;

			int l_p = clone_path(l_index, n, scl);
#ifdef ENABLE_DEBUG
	printf("\n\nInside loop before update llr path metric: l_index_1 = %d l_p = %d\n\n", l_index_1, l_p);
	print_scl_data_structure(n, block_len, list_size, scl);
#endif

			int l_index_2 = get_array_ptr_p(n, l_p, scl, code);
#ifdef ENABLE_DEBUG
	printf("\n\nInside loop before update llr path metric: l_index_1 = %d\n\n", l_index_1);
	print_scl_data_structure(n, block_len, list_size, scl);
#endif			
			for (int i = 0; i < phi; i++){		
				scl->i_scl[l_index_2][i] = scl->i_scl[l_index_1][i];
			}
			scl->c_scl[scl->lambda_offset[n]+scl->list_offset[l_index_2]][phi%2] = 1;
			scl->i_scl[l_index_2][phi] = 1;

#ifdef ENABLE_DEBUG
	printf("\n\nInside loop before update llr path metric: l_index_1 = %d\n\n", l_index_1);
	print_scl_data_structure(n, block_len, list_size, scl);
#endif


			// update scl->llr_path_metric
			if (use_log_approx==1){
				if (use_quantization == 1){
					scl->llr_path_metric[l_index] = quantize(scl->llr_path_metric[l_index]-scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]], step, max_value);
					scl->llr_path_metric[l_p] = quantize(scl->llr_path_metric[l_p]+scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_2]], step, max_value);
				}else{
					scl->llr_path_metric[l_index] = scl->llr_path_metric[l_index]-scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]];
					scl->llr_path_metric[l_p] = scl->llr_path_metric[l_p]+scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_2]];
				}
			}else{
				scl->llr_path_metric[l_index] = scl->llr_path_metric[l_index] + log(1 + exp(-scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]]));
				scl->llr_path_metric[l_p] = scl->llr_path_metric[l_p] + log(1 + exp(scl->llr_scl[scl->lambda_offset[n]+scl->list_offset[l_index_2]]));
			}

		}else{

			if (cont_forks[l_index][0] == 1){
				
				scl->c_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]][phi%2] = 0;
				scl->i_scl[l_index_1][phi] = 0;

				if (use_log_approx==1){
					if (use_quantization == 1){
						scl->llr_path_metric[l_index] = quantize(scl->llr_path_metric[l_index] - scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]], step, max_value); 
					}else{
						scl->llr_path_metric[l_index] = scl->llr_path_metric[l_index] - scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]];
					}
				}else{
					scl->llr_path_metric[l_index] = scl->llr_path_metric[l_index] + log(1 + exp(-scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]]));
				}

			}else{
				
				scl->c_scl[scl->lambda_offset[n]+scl->list_offset[l_index_1]][phi%2] = 1;
				scl->i_scl[l_index_1][phi] = 1;

				if (use_log_approx==1){
					if (use_quantization == 1){
						scl->llr_path_metric[l_index] = quantize(scl->llr_path_metric[l_index] + scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]], step, max_value); 
					}else{
						scl->llr_path_metric[l_index] = scl->llr_path_metric[l_index] + scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]];
					}
				}else{
					scl->llr_path_metric[l_index] = scl->llr_path_metric[l_index] + log(1 + exp(scl->llr_scl[scl->lambda_offset[n] + scl->list_offset[l_index_1]]));
				}
			}			
		}
	}

	for (int i = 0; i < list_size; i++){
		free(cont_forks[i]);
		free(prob_forks[i]);	
	}
	free(prob_forks);
	free(cont_forks);
	free(prob);
}


int find_most_probable_path( struct scl_data_structure *scl, struct polar_code *code ){
	
	int list_size = code->list_size;
	int n = code->n;

	int l_p_index = 0;
	double p_max = REALMAX;
	int path_with_crc = 0;

	for (int l_index = 0; l_index < list_size; l_index++){	
		
		if (scl->active_path_array[l_index] == 0){
			continue;
		}
		
		int c_index = get_array_ptr_p(n, l_index, scl, code);

		if (p_max > scl->llr_path_metric[l_index]) {
			p_max = scl->llr_path_metric[l_index];
			l_p_index = l_index;
		}
	}
	return l_p_index;
}

void init_data_structure_scl(int n, int block_len, int list_size, struct scl_data_structure *scl){

	for (int i = 0; i < (n+1); i++){
		scl->pathIdx_to_arrayIdx[i] = malloc(list_size*sizeof(int));
	}

	for (int i = 0; i < (n+1); i++){
		struct stack *temp = malloc(sizeof(struct stack));
		scl->inactive_array_indices[i] = temp;
		scl->inactive_array_indices[i]->array = NULL;
		scl->inactive_array_indices[i]->size = 0;
	}	

	for (int i = 0; i < (n+1); i++){
		scl->array_ref_count[i] = malloc(list_size*sizeof(int));
	}

	for (int i = 0; i < list_size*(2*block_len-1); i++){
		scl->c_scl[i] = malloc(2*sizeof(int));
	}

	for (int i = 0; i < list_size; i++){
		scl->i_scl[i] = malloc(block_len*sizeof(int));
	} 
	
	for (int i = 0; i < (n+1); i++){
		scl->lambda_offset[i] = pow(2, n-i)-1;
	}

	for (int i = 0; i < (list_size+1); i++){
		scl->list_offset[i] = i*(2*block_len-1);
	}

	for (int lambda = 0; lambda < n+1; lambda++){
		for (int i_list = 0; i_list < list_size; i_list++){
			push(scl->inactive_array_indices[lambda], i_list);
		}
	}

	for (int i_list = 0; i_list < list_size; i_list++){
		scl->active_path_array[i_list] = 0;
		push(scl->inactive_path_indices, i_list);
	}

}

void free_scl_data_structure(int n, int block_len, int list_size, struct scl_data_structure *scl){

	free(scl->inactive_path_indices->array);
	free(scl->active_path_array);

	for (int i = 0; i < (n+1); i++){
		free(scl->pathIdx_to_arrayIdx[i]);
	}
	free(scl->pathIdx_to_arrayIdx);
	
	for (int i = 0; i < (n+1); i++){
		free(scl->inactive_array_indices[i]->array);
		free(scl->inactive_array_indices[i]);
	}
	free(scl->inactive_array_indices);
	
	for (int i = 0; i < (n+1); i++){
		free(scl->array_ref_count[i]);
	}
	free(scl->array_ref_count);

	free(scl->llr_scl);
	free(scl->llr_path_metric);

	for (int i = 0; i < list_size*(2*block_len-1); i++){
		free(scl->c_scl[i]);
	}
	free(scl->c_scl);
	
	for (int i = 0; i < list_size; i++){
		free(scl->i_scl[i]);
	}
	free(scl->i_scl);

	free(scl->lambda_offset);
	free(scl->list_offset);
}

int get_i_scl(int lambda, int beta, int list_index, int *lambda_offset, int *list_offset){
	return (beta + lambda_offset[lambda] + list_offset[list_index]);
}


