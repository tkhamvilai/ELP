/*
 * sparse.c
 *
 *  Created on: Jul 3, 2020
 *      Author: tkhamvilai
 */

#include "ELP.h"

elp_status elp_matrix_init(matrix_t *m, size_t nRows, size_t nCols, size_t nnz, size_t *i, size_t *j, real_t *val, uint8_t is_dense)
{
#if MATRIX_CHECK
	if( nRows < 1 || nCols < 1 || nnz < 1 || val == NULL || (is_dense != 0 && is_dense != 1) )
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
#endif /* MATRIX_CHECK */
	m->nRows = nRows;
	m->nCols = nCols;
	m->nnz = nnz;
	m->val = val;
	m->is_dense = is_dense;
	if(!is_dense)
	{
#if MATRIX_CHECK
		if( (nRows > 1 && i == NULL) || (nCols > 1 && j == NULL) )
		{
			return ELP_STATUS_SIZE_MISMATCH;
		}
#endif /* MATRIX_CHECK */
		m->i = i;
		m->j = j;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_matrix_diag_trans_mult(const matrix_t *m_in, const matrix_t *d_in, matrix_t *m_out)
{
#if MATRIX_CHECK
	if(m_in == NULL || d_in == NULL || m_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (m_out->nCols != m_out->nRows) || (m_in->nCols != d_in->nRows) || (m_in->nRows != m_out->nRows) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */
	size_t *m_out_i = m_out->i;
	size_t *m_out_j = m_out->j;
	real_t *m_out_val = m_out->val;

	size_t pos_1, pos_2, pos_1_temp, pos_2_temp, i, j;
	real_t sum;
	pos_1 = m_in->nnz;
	while (pos_1 != 0U)
	{
		i = m_in->i[pos_1 - 1];
		pos_2 = pos_1;
		while (pos_2 != 0U)
		{
			j = m_in->i[pos_2 - 1];
			pos_1_temp = pos_1;
			pos_2_temp = pos_2;
			sum = 0;
			while (1)
			{
				if(m_in->j[pos_1_temp - 1] > m_in->j[pos_2_temp - 1])
				{
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
				}
				else if(m_in->j[pos_1_temp - 1] < m_in->j[pos_2_temp - 1])
				{
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}
				else
				{
					sum += m_in->val[pos_1_temp - 1] * m_in->val[pos_2_temp - 1] * d_in->val[m_in->j[pos_1_temp - 1]] * d_in->val[m_in->j[pos_1_temp - 1]];
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}

				if(m_in->j[pos_1_temp - 1] > m_in->j[pos_2_temp - 1])
				{
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
				}
				else if(m_in->j[pos_1_temp - 1] < m_in->j[pos_2_temp - 1])
				{
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}
				else
				{
					sum += m_in->val[pos_1_temp - 1] * m_in->val[pos_2_temp - 1] * d_in->val[m_in->j[pos_1_temp - 1]] * d_in->val[m_in->j[pos_1_temp - 1]];
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}

				if(m_in->j[pos_1_temp - 1] > m_in->j[pos_2_temp - 1])
				{
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
				}
				else if(m_in->j[pos_1_temp - 1] < m_in->j[pos_2_temp - 1])
				{
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}
				else
				{
					sum += m_in->val[pos_1_temp - 1] * m_in->val[pos_2_temp - 1] * d_in->val[m_in->j[pos_1_temp - 1]] * d_in->val[m_in->j[pos_1_temp - 1]];
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}

				if(m_in->j[pos_1_temp - 1] > m_in->j[pos_2_temp - 1])
				{
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
				}
				else if(m_in->j[pos_1_temp - 1] < m_in->j[pos_2_temp - 1])
				{
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}
				else
				{
					sum += m_in->val[pos_1_temp - 1] * m_in->val[pos_2_temp - 1] * d_in->val[m_in->j[pos_1_temp - 1]] * d_in->val[m_in->j[pos_1_temp - 1]];
					pos_1_temp--;
					if(m_in->i[pos_1_temp - 1] != i || pos_1_temp == 0U) break;
					pos_2_temp--;
					if(m_in->i[pos_2_temp - 1] != j || pos_2_temp == 0U) break;
				}
			}
			if (sum != 0)
			{
				*m_out_i++ = i;
				*m_out_j++ = j;
				*m_out_val++ = sum;
			}
			while(1)
			{
				pos_2--;
				if(m_in->i[pos_2 - 1] != j || pos_2 == 0U) break;
			}
		}
		while(1)
		{
			pos_1--;
			if(m_in->i[pos_1 - 1] != i || pos_1 == 0U) break;
		}
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_matrix_trans_vector_mult_add(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(m_in == NULL || v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (m_in->nRows != v_in->nRows) || (m_in->nCols != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */
	size_t *m_in_i = m_in->i;
	size_t *m_in_j = m_in->j;
	real_t *m_in_val = m_in->val;
	real_t *v_in_val = v_in->val;
	real_t out_temp[m_in->nnz];

	size_t block_cnt, nz = 0;
	block_cnt = m_in->nnz >> 2U;
	while (block_cnt != 0U)
	{
		out_temp[nz++] = (*m_in_val++) * v_in_val[*m_in_i++];
		out_temp[nz++] = (*m_in_val++) * v_in_val[*m_in_i++];
		out_temp[nz++] = (*m_in_val++) * v_in_val[*m_in_i++];
		out_temp[nz++] = (*m_in_val++) * v_in_val[*m_in_i++];
		block_cnt--;
	}
	block_cnt = m_in->nnz % 0x4U;

	while (block_cnt != 0U)
	{
		out_temp[nz++] = (*m_in_val++) * v_in_val[*m_in_i++];
		block_cnt--;
	}

	nz = 0;
	block_cnt = m_in->nnz >> 2U;
	while (block_cnt != 0U)
	{
		v_out->val[*m_in_j++] += out_temp[nz++];
		v_out->val[*m_in_j++] += out_temp[nz++];
		v_out->val[*m_in_j++] += out_temp[nz++];
		v_out->val[*m_in_j++] += out_temp[nz++];
		block_cnt--;
	}
	block_cnt = m_in->nnz % 0x4U;

	while (block_cnt != 0U)
	{
		v_out->val[*m_in_j++] += out_temp[nz++];
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_matrix_vector_mult_add(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(m_in == NULL || v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (m_in->nCols != v_in->nRows) || (m_in->nRows != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */
	size_t *m_in_i = m_in->i;
	size_t *m_in_j = m_in->j;
	real_t *m_in_val = m_in->val;
	real_t *v_in_val = v_in->val;

	size_t block_cnt;
	block_cnt = m_in->nnz >> 2U;
	while (block_cnt != 0U)
	{
		v_out->val[*m_in_i++] += (*m_in_val++) * v_in_val[*m_in_j++];
		v_out->val[*m_in_i++] += (*m_in_val++) * v_in_val[*m_in_j++];
		v_out->val[*m_in_i++] += (*m_in_val++) * v_in_val[*m_in_j++];
		v_out->val[*m_in_i++] += (*m_in_val++) * v_in_val[*m_in_j++];
		block_cnt--;
	}
	block_cnt = m_in->nnz % 0x4U;

	while (block_cnt != 0U)
	{
		v_out->val[*m_in_i++] += (*m_in_val++) * v_in_val[*m_in_j++];
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_matrix_vector_mult_sub(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(m_in == NULL || v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (m_in->nCols != v_in->nRows) || (m_in->nRows != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */
	size_t *m_in_i = m_in->i;
	size_t *m_in_j = m_in->j;
	real_t *m_in_val = m_in->val;
	real_t *v_in_val = v_in->val;

	size_t block_cnt;
	block_cnt = m_in->nnz >> 2U;
	while (block_cnt != 0U)
	{
		v_out->val[*m_in_i++] -= (*m_in_val++) * v_in_val[*m_in_j++];
		v_out->val[*m_in_i++] -= (*m_in_val++) * v_in_val[*m_in_j++];
		v_out->val[*m_in_i++] -= (*m_in_val++) * v_in_val[*m_in_j++];
		v_out->val[*m_in_i++] -= (*m_in_val++) * v_in_val[*m_in_j++];
		block_cnt--;
	}
	block_cnt = m_in->nnz % 0x4U;

	while (block_cnt != 0U)
	{
		v_out->val[*m_in_i++] -= (*m_in_val++) * v_in_val[*m_in_j++];
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_matrix_diag_vector_mult(const matrix_t *m_in, const matrix_t *v_in, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(m_in == NULL || v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (m_in->nRows != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */
	real_t *m_in_val = m_in->val;
	real_t *v_in_val = v_in->val;
	real_t *v_out_val = v_out->val;

	size_t block_cnt;
	block_cnt = m_in->nRows >> 2U;
	while (block_cnt != 0U)
	{
		(*v_out_val++) = (*m_in_val) * (*m_in_val) * (*v_in_val++); m_in_val++;
		(*v_out_val++) = (*m_in_val) * (*m_in_val) * (*v_in_val++); m_in_val++;
		(*v_out_val++) = (*m_in_val) * (*m_in_val) * (*v_in_val++); m_in_val++;
		(*v_out_val++) = (*m_in_val) * (*m_in_val) * (*v_in_val++); m_in_val++;
		block_cnt--;
	}
	block_cnt = m_in->nRows % 0x4U;

	while (block_cnt != 0U)
	{
		(*v_out_val++) = (*m_in_val) * (*m_in_val) * (*v_in_val++);
		m_in_val++;
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_vector_scale(const matrix_t *v_in, real_t scale, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (v_in->nRows != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */

	real_t *v_in_val = v_in->val;
	real_t *v_out_val = v_out->val;

	size_t block_cnt;
	block_cnt = v_in->nRows >> 2U;

	while (block_cnt != 0U)
	{
		*v_out_val++ = scale * (*v_in_val++);
		*v_out_val++ = scale * (*v_in_val++);
		*v_out_val++ = scale * (*v_in_val++);
		*v_out_val++ = scale * (*v_in_val++);
		block_cnt--;
	}
	block_cnt = v_in->nRows % 0x4U;

	while (block_cnt != 0U)
	{
		*v_out_val++ = scale * (*v_in_val++);
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_vector_scale_add(const matrix_t *v_in, real_t scale, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (v_in->nRows != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */

	real_t *v_in_val = v_in->val;
	real_t *v_out_val = v_out->val;

	size_t block_cnt;
	block_cnt = v_in->nRows >> 2U;

	while (block_cnt != 0U)
	{
		*v_out_val++ += scale * (*v_in_val++);
		*v_out_val++ += scale * (*v_in_val++);
		*v_out_val++ += scale * (*v_in_val++);
		*v_out_val++ += scale * (*v_in_val++);
		block_cnt--;
	}
	block_cnt = v_in->nRows % 0x4U;

	while (block_cnt != 0U)
	{
		*v_out_val++ += scale * (*v_in_val++);
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_vector_recip(const matrix_t *v_in, real_t scale, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(v_in == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (v_in->nRows != v_out->nRows) || (v_in->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */

	real_t *v_in_val = v_in->val;
	real_t *v_out_val = v_out->val;

	size_t block_cnt;
	block_cnt = v_in->nRows >> 2U;

	while (block_cnt != 0U)
	{
		*v_out_val++ += scale / (*v_in_val++);
		*v_out_val++ += scale / (*v_in_val++);
		*v_out_val++ += scale / (*v_in_val++);
		*v_out_val++ += scale / (*v_in_val++);
		block_cnt--;
	}
	block_cnt = v_in->nRows % 0x4U;

	while (block_cnt != 0U)
	{
		*v_out_val++ += scale / (*v_in_val++);
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_vector_subtract(const matrix_t *v_in1, const matrix_t *v_in2, matrix_t *v_out)
{
#if MATRIX_CHECK
	if(v_in1 == NULL || v_in2 == NULL || v_out == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (v_in1->nRows != v_out->nRows) || (v_in2->nRows != v_out->nRows) || (v_in1->nRows != v_in2->nRows) || (v_in1->nCols != 1) || (v_in2->nCols != 1) || (v_out->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */

	real_t *v_in1_val = v_in1->val;
	real_t *v_in2_val = v_in2->val;
	real_t *v_out_val = v_out->val;

	size_t block_cnt;
	block_cnt = v_in1->nRows >> 2U;

	while (block_cnt != 0U)
	{
		*v_out_val++ = (*v_in1_val++) - (*v_in2_val++);
		*v_out_val++ = (*v_in1_val++) - (*v_in2_val++);
		*v_out_val++ = (*v_in1_val++) - (*v_in2_val++);
		*v_out_val++ = (*v_in1_val++) - (*v_in2_val++);
		block_cnt--;
	}
	block_cnt = v_in1->nRows % 0x4U;

	while (block_cnt != 0U)
	{
		*v_out_val++ = (*v_in1_val++) - (*v_in2_val++);
		block_cnt--;
	}
	return ELP_STATUS_SUCCESS;
}

elp_status elp_vector_min(const matrix_t *v_in, real_t *val, size_t *ind)
{
#if MATRIX_CHECK
	if(v_in == NULL)
	{
		return ELP_STATUS_ARGUMENT_ERROR;
	}
	if( (v_in->nRows < 1) || (v_in->nCols != 1) )
	{
		return ELP_STATUS_SIZE_MISMATCH;
	}
#endif /* MATRIX_CHECK */
	real_t *v_in_val = v_in->val;
	real_t min_val1, min_val2, out;
	size_t block_cnt, out_ind, cnt;
	cnt = 0U;
	out_ind = 0U;
	out = *v_in_val++;

	block_cnt = (v_in->nRows - 1) >> 2U;
	while(block_cnt != 0U)
	{
		min_val1 = *v_in_val++;
		min_val2 = *v_in_val++;
		if (out > min_val1)
		{
			out = min_val1;
			out_ind = cnt + 1U;
		}
		if (out > min_val2)
		{
			out = min_val2;
			out_ind = cnt + 2U;
		}

		min_val1 = *v_in_val++;
		min_val2 = *v_in_val++;
		if (out > min_val1)
		{
			out = min_val1;
			out_ind = cnt + 3U;
		}
		if (out > min_val2)
		{
			out = min_val2;
			out_ind = cnt + 4U;
		}

		cnt += 4U;
		block_cnt--;
	}

	block_cnt = (v_in->nRows - 1) % 4U;

	while(block_cnt != 0U)
	{
		min_val1 = *v_in_val++;
		if (out > min_val1)
		{
			out = min_val1;
			out_ind = v_in->nRows - block_cnt;
		}
		block_cnt--;
	}
	*val = out;
	*ind = out_ind;
	return ELP_STATUS_SUCCESS;
}
