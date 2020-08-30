/*
 * test.c
 *
 *  Created on: Jul 3, 2020
 *      Author: tkhamvilai
 */

#include "elp_demo.h"

void elp_demo1()
{
	elp_t *prob = (elp_t*)elp_malloc(sizeof(elp_t));
	elp_prob_init(prob);

	size_t cnt = 0;
	prob->A.i[cnt] = 0; prob->A.j[cnt] = 0; prob->A.val[cnt] = 1; cnt++;
	prob->A.i[cnt] = 0; prob->A.j[cnt] = 1; prob->A.val[cnt] = 1; cnt++;
	prob->A.i[cnt] = 0; prob->A.j[cnt] = 2; prob->A.val[cnt] = 1; cnt++;
	prob->A.i[cnt] = 0; prob->A.j[cnt] = 3; prob->A.val[cnt] = 1; cnt++;
	prob->A.i[cnt] = 1; prob->A.j[cnt] = 0; prob->A.val[cnt] = 10; cnt++;
	prob->A.i[cnt] = 1; prob->A.j[cnt] = 1; prob->A.val[cnt] = 4; cnt++;
	prob->A.i[cnt] = 1; prob->A.j[cnt] = 2; prob->A.val[cnt] = 5; cnt++;
	prob->A.i[cnt] = 1; prob->A.j[cnt] = 4; prob->A.val[cnt] = 1; cnt++;
	prob->A.i[cnt] = 2; prob->A.j[cnt] = 0; prob->A.val[cnt] = 2; cnt++;
	prob->A.i[cnt] = 2; prob->A.j[cnt] = 1; prob->A.val[cnt] = 2; cnt++;
	prob->A.i[cnt] = 2; prob->A.j[cnt] = 2; prob->A.val[cnt] = 6; cnt++;
	prob->A.i[cnt] = 2; prob->A.j[cnt] = 5; prob->A.val[cnt] = 1; cnt++;
	assert(cnt == NNZ_A);

	cnt = 0;
	prob->b.val[cnt++] = 100;
	prob->b.val[cnt++] = 600;
	prob->b.val[cnt++] = 300;

	cnt = 0;
	prob->c.val[cnt++] = -10;
	prob->c.val[cnt++] = -6;
	prob->c.val[cnt++] = -4;

	__HAL_TIM_SET_COUNTER(&htim9,0);
	elp_interior_point_outer(prob);
	printf("\r\nELP demo1: %lu us\r\n", __HAL_TIM_GET_COUNTER(&htim9));

	printf("ELP solution (including slacks): ");
	for(unsigned int i = 0; i < N; i++)
	{
		printf("%d, ", (int)(prob->x.val[i]));
	}
	printf("\r\n");

	return;
}

void elp_demo_main()
{
	elp_demo1();
	return;
}
