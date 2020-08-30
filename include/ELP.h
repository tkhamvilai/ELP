/*
 * ELP.h
 *
 *  Created on: Jul 3, 2020
 *      Author: tkhamvilai
 */

#ifndef ELP_H_
#define ELP_H_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tim.h"

#ifdef USE_RTOS
#define elp_malloc pvPortMalloc
#define elp_free vPortFree
#else
#define elp_malloc malloc
#define elp_free free
#endif /* USE_FREERTOS */

#if __CORTEX_M == 4
typedef float real_t;
#elif __CORTEX_M == 7
typedef double real_t;
#endif

#define VERBOSE 0
#define MATRIX_CHECK 0

#define M 3 /* # of constraints (including bounds) */
#define N 6 /* # of variables (including slacks) */
#define NNZ_A 12

#define ALPHA 10.f /* rate for improving a logarithmic barrier approximation */
#define BETA 0.5f /* rate for decreasing a step size */
#define ELP_MAXITR_INNER 5
#define ELP_MAXITR_OUTER 2

#define ELP_STATUS_SUCCESS 0
#define ELP_STATUS_SIZE_MISMATCH 1
#define ELP_STATUS_ARGUMENT_ERROR 2

typedef uint8_t elp_status;

typedef struct matrix_t
{
	uint8_t is_dense;
	size_t nRows;
	size_t nCols;
	size_t nnz;
	size_t *i;
	size_t *j;
	real_t *val;
} matrix_t;

typedef struct elp_t
{
	matrix_t A;
	matrix_t b;
	matrix_t c;

	matrix_t x;
	matrix_t w;

	matrix_t g;
	matrix_t h;

	matrix_t AHAt;
	matrix_t c_temp;
	matrix_t g_temp;
	matrix_t x_temp;

	size_t *P;
	size_t *etree;
	matrix_t L;
	real_t *D;

	int etree_flag;
} elp_t;

void elp_prob_init(elp_t *prob);
void elp_interior_point_inner(elp_t *prob);
void elp_interior_point_outer(elp_t *prob);

elp_status elp_matrix_init(matrix_t *m, size_t nRows, size_t nCols, size_t nnz, size_t *i, size_t *j, real_t *val, uint8_t is_dense);
elp_status elp_matrix_diag_trans_mult(const matrix_t *m_in, const matrix_t *d_in, matrix_t *m_out);
elp_status elp_matrix_vector_mult_add(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out);
elp_status elp_matrix_vector_mult_sub(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out);
elp_status elp_matrix_trans_vector_mult_add(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out);
elp_status elp_matrix_diag_vector_mult(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out);
elp_status elp_vector_scale(const matrix_t *v_in, real_t scale, matrix_t *v_out);
elp_status elp_vector_scale_add(const matrix_t *v_in, real_t scale, matrix_t *v_out);
elp_status elp_vector_recip(const matrix_t *v_in, real_t scale, matrix_t *v_out);
elp_status elp_vector_subtract(const matrix_t *v_in1, const matrix_t *v_in2, matrix_t *v_out);
elp_status elp_vector_min(const matrix_t *v_in, real_t *val, size_t *ind);

void elp_ldlt_symbolic(const matrix_t *A, matrix_t *L, size_t *etree);
void elp_ldlt_numeric(const matrix_t *A, const size_t *etree, matrix_t *L, real_t *D);
void elp_lsolve(const matrix_t *L, real_t *x);
void elp_dsolve(const real_t *L, real_t *x);
void elp_ltsolve(const matrix_t *L, real_t *x);

#endif /* ELP_H_ */
