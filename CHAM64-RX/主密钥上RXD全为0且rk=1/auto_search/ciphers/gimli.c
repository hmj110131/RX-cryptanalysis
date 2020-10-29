/*
 * gimli.c
 *
 *  Created on: 2018年1月18日
 *      Author: hmj110131
 */
#include "gimli.h"
#include "typedef.h"
#include <stdio.h>





u32 GIMLI_state_word[25][12] = {0};   // 24 round with 25 states.



u64 gimli_hwval_space[33][65536] = {0}; // hamming weight space value with corresponding hw.0-32
u64 gimli_hwspace_num[33] = {0};   // the number of corresponding weight 0-32.






/**
 * Compute the tuple or difference of x y z.
 *
 */
u32 gimli_OR_table[8][4] = {
		  //xyz = 0 0 0
		  {0,0,0,0},
		  //xyz = 0 0 1
		  {0,0,0,0},
		  //xyz = 0 1 0
		  {0,0,0,0},
		  //xyz = 0 1 1
		  {0,0,0,0},
		  //xyz = 1 0 0
		  {0,1,2,3},
		  //xyz = 1 0 1
		  {0,0,0,0},
		  //xyz = 1 1 0
		  {2,3,4,5},
		  //xyz = 1 1 1
		  {1,2,4,7}

};



