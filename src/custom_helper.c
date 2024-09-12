#include <stdlib.h>  
#include <stdio.h>  
#include <malloc.h>  
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <string.h>

#include "../header/polar_scl.h"

#define ENABLE_DEBUG

double min(double a, double b){
	if(a > b){
		return b;
	}else{
		return a;
	}
}

double max(double a, double b){
	if(a > b){
		return a;
	}else{
		return b;
	}
}

double min_array(double *array, int array_size){
	double min = array[0];
	for (int i = 0; i < array_size; i++){
		if (array[i] < min){
			min = array[i];
		}
	}
	return min;
}

int sign(double in){
	if (in < 0){
		return -1;	
	}else{
		return 1;
	}
}

int compare_double( const void* a, const void* b )
{
    if( *(double*)a == *(double*)b ) return 0;
    return *(double*)a < *(double*)b ? -1 : 1;
}

int compare_double_descend( const void* a, const void* b )
{
    if( *(double*)a == *(double*)b ) return 0;
    return *(double*)a < *(double*)b ? 1 : -1;
}

int compare_int( const void* a, const void* b )
{
    if( *(int*)a == *(int*)b ) return 0;
    return *(int*)a < *(int*)b ? -1 : 1;
}

void push(struct stack *a, int value){
	if(a->size == 0){
		a->array = malloc(sizeof(int));
		a->array[0] = value;
		a->size = 1;
	}else{
		a->array = realloc(a->array, a->size+1);
		a->array[a->size] = value;
		a->size++;
	}
}

void pop(struct stack *a){
	if(a->size == 0){
		printf("\nError: stack is empty in pop\n");
		exit(1);
	}else{
		a->array = realloc(a->array, a->size-1);
		a->size--;
	}
}

int top(struct stack *a){

	if(a->size == 0){
		printf("\nError: stack is empty in top\n");
		exit(1);
	}else{
		return a->array[a->size-1];
	}
}

double quantize(double llr, double step, double max_value){
	
        if(llr >= max_value){
            return max_value;
        }else if (llr <= -max_value){
            return -max_value;
        }else{
            return round(llr/step)*step;
        }
}

// f function in log approximation
double f_func(double a, double b, int use_quantize, double step, double max_value){
	double result = sign(a)*sign(b)*min(abs(a), abs(b));
	if (use_quantize==1){
		result = quantize(result, step, max_value);
	}
	return result;
}

// g function in log approximation
double g_func(double a, double b, int u_p, int use_quantize, double step, double max_value){
	double result = pow(-1, u_p)*a + b;
	if (use_quantize==1){
		result = quantize(result, step, max_value);
	}
	return result; 
}

