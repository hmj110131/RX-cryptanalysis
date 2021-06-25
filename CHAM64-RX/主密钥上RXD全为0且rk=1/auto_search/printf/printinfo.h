/*
 * printinfo.h
 *
 *  Created on: 2017年8月23日
 *      Author: hmj110131
 */

#ifndef PRINTF_PRINTINFO_H_
#define PRINTF_PRINTINFO_H_

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>

#include "ciphers.h"
#include "search.h"
#include "typedef.h"
#include "globalvar.h"

FILE* simon_diff_trail;
FILE* best_Bn;
FILE* speck_best_Bn;
FILE* hight_best_Bn;
FILE* simon_linear_trail;
FILE* sparx_best_Bn;
FILE* sparx_128_best_Bn;
FILE* best_Bn_linear;
FILE* result_print;
clock_t time_start, time_ARX_DDT,time_ARX_cLAT, time_Round, time_finish, time_xor_talbe;
double  run_time;
/**
 * @Print the final result of DC or LC.
 *
 */
//void print_resoult(u16 *xin, u16 *yin,u16 *xout,u16 *yout,u16 rounds, u16 *r_w,u16 w_total);
void print_resoult(u16 search_round);
void test_print(void);

#endif /* PRINTF_PRINTINFO_H_ */
