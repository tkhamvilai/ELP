/*
 * ELP.c
 *
 *  Created on: Jul 3, 2020
 *      Author: tkhamvilai
 */

#include "ELP.h"

void elp_interior_point_inner(elp_t *prob)
{
	for(size_t itr = ELP_MAXITR_INNER; itr != 0; itr--)
	{
		/* calculate g = c - 1/x + A'w */
		memcpy(prob->g.val, prob->c_temp.val, N * sizeof(real_t));
		assert(elp_vector_recip(&prob->x, -1.f, &prob->g) == ELP_STATUS_SUCCESS);
		assert(elp_matrix_trans_vector_mult_add(&prob->A, &prob->w, &prob->g) == ELP_STATUS_SUCCESS);

		/* calculate h = b - Ax */
		memcpy(prob->h.val, prob->b.val, M * sizeof(real_t));
		assert(elp_matrix_vector_mult_sub(&prob->A, &prob->x, &prob->h) == ELP_STATUS_SUCCESS);

		/* calculate h + A*H_inv*g */
		assert(elp_matrix_diag_vector_mult(&prob->x, &prob->g, &prob->g_temp) == ELP_STATUS_SUCCESS);
		assert(elp_matrix_vector_mult_add(&prob->A, &prob->g_temp, &prob->h) == ELP_STATUS_SUCCESS);

		/* calculate A*H_inv*A' */
		assert(elp_matrix_diag_trans_mult(&prob->A, &prob->x, &prob->AHAt) == ELP_STATUS_SUCCESS);

		/* solve (A*H_inv*A')*dw = (h + A*H_inv*g) */
		if(!prob->etree_flag)
		{
			elp_ldlt_symbolic(&prob->AHAt, &prob->L, prob->etree);
			prob->etree_flag = 1;
		}
		elp_ldlt_numeric(&prob->AHAt, prob->etree, &prob->L, prob->D);
		elp_lsolve(&prob->L, prob->h.val);
		elp_dsolve(prob->D, prob->h.val);
		elp_ltsolve(&prob->L, prob->h.val);

		/* calculate dx = H_inv*(A'*dw - g) */
		memset(prob->g_temp.val, 0, N * sizeof(real_t));
		assert(elp_matrix_trans_vector_mult_add(&prob->A, &prob->h, &prob->g_temp) == ELP_STATUS_SUCCESS);
		assert(elp_vector_subtract(&prob->g_temp, &prob->g, &prob->g_temp) == ELP_STATUS_SUCCESS);
		assert(elp_matrix_diag_vector_mult(&prob->x, &prob->g_temp, &prob->g_temp) == ELP_STATUS_SUCCESS);

		/* choose a step size that does not bring x outside the domain */
		real_t t = 1.f, min_x;
		do
		{
			memcpy(prob->x_temp.val, prob->x.val, N * sizeof(real_t));
			assert(elp_vector_scale_add(&prob->g_temp, t, &prob->x_temp) == ELP_STATUS_SUCCESS); /* t*dx */
			assert(elp_vector_min(&prob->x_temp, &min_x, NULL) == ELP_STATUS_SUCCESS); /* min(x_temp) */
			if(min_x <= 0) t *= BETA;
		}
		while(min_x <= 0);

		memcpy(prob->x.val, prob->x_temp.val, N * sizeof(real_t)); /* update x */
		assert(elp_vector_scale_add(&prob->h, t, &prob->w) == ELP_STATUS_SUCCESS); /* update w */
	}
	return;
}

void elp_interior_point_outer(elp_t *prob)
{
	float t = 1.f;
	for(size_t itr = ELP_MAXITR_OUTER; itr != 0; itr--)
	{
		assert(elp_vector_scale(&prob->c, t, &prob->c_temp) == ELP_STATUS_SUCCESS);
		elp_interior_point_inner(prob);
		t *= ALPHA;
	}
	return;
}

void elp_prob_init(elp_t *prob)
{
	size_t *A_i = (size_t*)elp_malloc(NNZ_A * sizeof(size_t)); assert(A_i != NULL);
	size_t *A_j = (size_t*)elp_malloc(NNZ_A * sizeof(size_t)); assert(A_j != NULL);
	real_t *A_val = (real_t*)elp_malloc(NNZ_A * sizeof(real_t)); assert(A_val != NULL);
	assert(elp_matrix_init(&prob->A, M, N, NNZ_A, A_i, A_j, A_val, 0) == ELP_STATUS_SUCCESS);

	real_t *b_val = (real_t*)elp_malloc(M * sizeof(real_t)); assert(b_val != NULL);
	memset(b_val, 0, M * sizeof(real_t));
	assert(elp_matrix_init(&prob->b, M, 1, M, NULL, NULL, b_val, 1) == ELP_STATUS_SUCCESS);

	real_t *c_val = (real_t*)elp_malloc(N * sizeof(real_t)); assert(c_val != NULL);
	memset(c_val, 0, N * sizeof(real_t));
	assert(elp_matrix_init(&prob->c, N, 1, N, NULL, NULL, c_val, 1) == ELP_STATUS_SUCCESS);

	real_t *x_val = (real_t*)elp_malloc(N * sizeof(real_t)); assert(x_val != NULL);
	for(unsigned int i = N; i != 0; i--) x_val[i-1] = 1.f; /* any non-zero initial point */
	assert(elp_matrix_init(&prob->x, N, 1, N, NULL, NULL, x_val, 1) == ELP_STATUS_SUCCESS);

	real_t *w_val = (real_t*)elp_malloc(M * sizeof(real_t)); assert(w_val != NULL);
	memset(w_val, 0, M * sizeof(real_t));
	assert(elp_matrix_init(&prob->w, M, 1, M, NULL, NULL, w_val, 1) == ELP_STATUS_SUCCESS);

	real_t *g_val = (real_t*)elp_malloc(N * sizeof(real_t)); assert(g_val != NULL);
	assert(elp_matrix_init(&prob->g, N, 1, N, NULL, NULL, g_val, 1) == ELP_STATUS_SUCCESS);

	real_t *h_val = (real_t*)elp_malloc(M * sizeof(real_t)); assert(h_val != NULL);
	memset(h_val, 0, M * sizeof(real_t));
	assert(elp_matrix_init(&prob->h, M, 1, M, NULL, NULL, h_val, 1) == ELP_STATUS_SUCCESS);

	size_t *AHAt_i = (size_t*)elp_malloc( (M*(M+1) >> 1U) * sizeof(size_t)); assert(AHAt_i != NULL);
	size_t *AHAt_j = (size_t*)elp_malloc( (M*(M+1) >> 1U) * sizeof(size_t)); assert(AHAt_j != NULL);
	real_t *AHAt_val = (real_t*)elp_malloc( (M*(M+1) >> 1U) * sizeof(real_t)); assert(AHAt_val != NULL);
	assert(elp_matrix_init(&prob->AHAt, M, M, (M*(M+1) >> 1U), AHAt_i, AHAt_j, AHAt_val, 0) == ELP_STATUS_SUCCESS);

	prob->P = (size_t*)elp_malloc(M * sizeof(size_t)); assert(prob->P != NULL);
	prob->etree = (size_t*)elp_malloc(M * sizeof(size_t)); assert(prob->etree != NULL);
	prob->D = (real_t*)elp_malloc(M * sizeof(real_t)); assert(prob->D != NULL);

	size_t *L_i = (size_t*)elp_malloc( (M*(M-1) >> 1U) * sizeof(size_t)); assert(L_i != NULL);
	memset(L_i, 0, (M*(M-1) >> 1U) * sizeof(size_t));
	size_t *L_j = (size_t*)elp_malloc( (M*(M-1) >> 1U) * sizeof(size_t)); assert(L_j != NULL);
	memset(L_j, 0, (M*(M-1) >> 1U) * sizeof(size_t));
	real_t *L_val = (real_t*)elp_malloc( (M*(M-1) >> 1U) * sizeof(real_t)); assert(L_val != NULL);
	assert(elp_matrix_init(&prob->L, M, M, (M*(M-1) >> 1U), L_i, L_j, L_val, 0) == ELP_STATUS_SUCCESS);

	real_t *c_temp_val = (real_t*)elp_malloc(N * sizeof(real_t)); assert(c_temp_val != NULL);
	assert(elp_matrix_init(&prob->c_temp, N, 1, N, NULL, NULL, c_temp_val, 1) == ELP_STATUS_SUCCESS);

	real_t *g_temp_val = (real_t*)elp_malloc(N * sizeof(real_t)); assert(g_temp_val != NULL);
	assert(elp_matrix_init(&prob->g_temp, N, 1, N, NULL, NULL, g_temp_val, 1) == ELP_STATUS_SUCCESS);

	real_t *x_temp_val = (real_t*)elp_malloc(N * sizeof(real_t)); assert(x_temp_val != NULL);
	assert(elp_matrix_init(&prob->x_temp, N, 1, N, NULL, NULL, x_temp_val, 1) == ELP_STATUS_SUCCESS);

	prob->etree_flag = 0;

	return;
}



