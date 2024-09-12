// SCL decoding of polar codes
// Author: Yaoyu Tao
// 05/12/2017
// header file

#define REALMAX 100000

// struct for polar codes
struct polar_code {
	int n; 				// number of layers
	int block_len;			// block length
	int info_len;			// info bits length
	int crc_size;			// crc coding size
	double epsilon;			// epsilon for channel polarization

	int *frozen_bits;		// positions of frozen bits
	int *info_bits_index;
	int *bit_reverse_order;		

	int list_size;			// list size for list decoding
};

// struct for successive cancellation list decoding
struct scl_data_structure {
	struct stack *inactive_path_indices;	// stack
	int *active_path_array;
	int **pathIdx_to_arrayIdx;
	struct stack **inactive_array_indices;	// each element is stack
	int **array_ref_count;
	double *llr_scl;
	double *llr_path_metric;
	int **c_scl;
	int **i_scl;

	int *lambda_offset;
	int *list_offset;
};

// struct for stack
struct stack {
	int *array;
	int size;
};

// custom helper functions
double min(double a, double b);
double max(double a, double b);
double min_array(double *array, int array_size);
int sign(double in);
int compare_double( const void* a, const void* b );
int compare_double_descend( const void* a, const void* b );
int compare_int( const void* a, const void* b );

// helper functions for stack
void push(struct stack *a, int value);
void pop(struct stack *a);
int top(struct stack *a);

// simulation related helper functions
double quantize(double llr, double step, double max_value);
double f_func(double a, double b, int use_quantize, double step, double max_value);
double g_func(double a, double b, int u_p, int use_quantize, double step, double max_value);
double gaussian(double mean, double std);
double *addNoise(int *encoded_bits, int block_len, int info_len, double ebno);

// polar codes functions
void bit_reverse_order(struct polar_code *code);
int* initialize_frozen_bits (struct polar_code *code);
int *polar_encode_core(int *info_bits_padded, int array_size);
int *polar_encode(int *info_bits, struct polar_code *code, int *channel_order_sorted);
int get_array_ptr_p(int lambda, int l_index, struct scl_data_structure *scl, struct polar_code *code);
void recursively_calc_p_scl (int lambda, int phi, struct scl_data_structure *scl, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization);
void recursively_update_c_scl(int lambda, int phi, struct scl_data_structure *scl, struct polar_code *code);
void continue_paths_frozen_bit(int phi, struct scl_data_structure *scl, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization );
int clone_path(int l_index, int n, struct scl_data_structure *scl);
void kill_path(int l_index, int n, struct scl_data_structure *scl);
void continue_paths_unfrozen_bit(int phi, struct scl_data_structure *scl, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization );
int find_most_probable_path( struct scl_data_structure *scl, struct polar_code *code );
void init_data_structure_scl(int n, int block_len, int list_size, struct scl_data_structure *scl);
int get_i_scl(int lambda, int beta, int list_index, int *lambda_offset, int *list_offset);
int *decode_scl(double *llr, struct polar_code *code, double step, double max_value, int use_log_approx, int use_quantization);
void free_scl_data_structure(int n, int block_len, int list_size, struct scl_data_structure *scl);

// debug printer functions
void print_initialization(struct polar_code *code, double *channels, double *channels_bitrev, double *channels_bitrev_sorted, int *channel_order_sorted);
void print_polar_debug(struct polar_code code, int *encoded_bits, double *received_llr);
void print_unfrozen_bits_calculation(int rho, int **cont_forks, double **prob_forks, double *prob, int list_size);
void print_scl_data_structure(int n, int block_len, int list_size, struct scl_data_structure *scl);
void print_polar_encode(struct polar_code *code, int *info_bits, int *info_bits_padded, int info_bits_count);

