/*
 * ldl.c
 *
 *  Created on: Jul 6, 2020
 *      Author: tkhamvilai
 */

#include "ELP.h"

void elp_ldlt_symbolic(const matrix_t *A, matrix_t *L, size_t *etree)
{
	size_t flag[M], Lnz[M];
	size_t pos_1, pos_2, mk;

	pos_1 = A->nnz;
	for (size_t k = M; k != 0; k--)
	{
		mk = M - k;
		etree[mk] = M;
		flag[mk] = M - k;
		Lnz[mk] = 0;
		while (1)
		{
			pos_2 = A->j[pos_1 - 1];
			if(pos_2 < mk)
			{
				while(flag[pos_2] != mk)
				{
					if(etree[pos_2] == M) etree[pos_2] = mk;
					Lnz[pos_2]++;
					flag[pos_2] = mk;
					pos_2 = etree[pos_2];
				}
			}
			pos_1--;
			if(A->i[pos_1 - 1] != mk || pos_1 == 0U) break;
		}
	}

	for (size_t k = 0; k < M; k++)
	{
		L->j[k + 1] = L->j[k] + Lnz[k];
	}
	return;
}

void elp_ldlt_numeric(const matrix_t *A, const size_t *etree, matrix_t *L, real_t *D)
{
	real_t yi, l_ki;
	real_t y[M] = {0};
	size_t stack[M] = {0}, flag[M], Lnz[M];
	size_t pos_1, pos_2, pos_3, len, top, mk, p;
	pos_1 = A->nnz;
	for (size_t k = M; k != 0; k--)
	{
		mk = M - k;
		top = M;
		flag[mk] = M - k;
		Lnz[mk] = 0;
		while (1)
		{
			pos_2 = A->j[pos_1 - 1];
			len = 0;
			if(pos_2 <= mk)
			{
				y[pos_2] += A->val[pos_1 - 1];
				while(flag[pos_2] != mk)
				{
					stack[len++] = pos_2;
					flag[pos_2] = mk;
					pos_2 = etree[pos_2];
				}
				while(len > 0) stack[--top] = stack[--len];
			}
			pos_1--;
			if(A->i[pos_1 - 1] != mk || pos_1 == 0U) break;
		}

		D[mk] = y[mk];
		y[mk] = 0.0;
		while(top < M)
		{
			pos_2 = stack[top];
			yi = y[pos_2];
			y[pos_2] = 0.0;

			pos_3 = L->j[pos_2] + Lnz[pos_2];
			for(p = L->j[pos_2]; p < pos_3; p++)
			{
				y[L->i[p]] -= L->val[p] * yi;
			}
			l_ki = yi / D [pos_2] ;
			D[mk] -= l_ki * yi ;
			L->i[p] = mk ;
			L->val[p] = l_ki;
			Lnz[pos_2]++;
			top++;
		}
		assert(D[mk] != 0.0);
	}
	return;
}

void elp_lsolve(const matrix_t *L, real_t *x)
{
	size_t j, pos_1, pos_2;
	for(j = 0; j < M; j++)
	{
		pos_2 = L->j[j+1];
		for(pos_1 = L->j[j]; pos_1 < pos_2; pos_1++)
		{
			x[L->i[pos_1]] -= L->val[pos_1] * x[j];
		}
	}
}

void elp_dsolve(const real_t *D, real_t *x)
{
	size_t block_cnt;
	block_cnt = M >> 2U;

	while (block_cnt != 0U)
	{
		*x++ /= *D++;
		*x++ /= *D++;
		*x++ /= *D++;
		*x++ /= *D++;
		block_cnt--;
	}
	block_cnt = M % 0x4U;

	while (block_cnt != 0U)
	{
		*x++ /= *D++;
		block_cnt--;
	}
}

void elp_ltsolve(const matrix_t *L, real_t *x)
{
	size_t j, pos_1, pos_2;
	for(j = M; j != 0; j--)
	{
		pos_2 = L->j[j];
		for(pos_1 = L->j[j-1]; pos_1 < pos_2; pos_1++)
		{
			x[j-1] -= L->val[pos_1] * x[L->i[pos_1]];
		}
	}
}
