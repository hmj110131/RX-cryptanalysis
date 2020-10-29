/*
 * gimli.h
 *
 *  Created on: 2018年1月18日
 *      Author: hmj110131
 */
#ifndef CIPHERS_GIMLI_H_
#define CIPHERS_GIMLI_H_

#include <stdio.h>
#include "typedef.h"
// Macros
#define GIMLI_ROUNDS 24 /**< Max. number of rounds */

#define GIMLI_ROUND_CONSTANT  0x9E377900  // Round constant RC.

#define GIMLI_LROT_CONST_24 24  // Left rotation 24 bits.
#define GIMLI_LROT_CONST_9  9   // Left rotation 9 bits.

#define GIMLI_Lshift_1  1  // Left shift 1 bits.
#define GIMLI_Lshift_2  2  // Left shift 2 bits.
#define GIMLI_Lshift_3  3  // Left shift 3 bits.

#define GIMLI_SmallSwap_round(x)  (x == 24) || (x == 20) || (x == 16) || (x == 12) || (x == 8) ||(x == 4)
#define GIMLI_BigSwap_round(x)    (x == 22) || (x == 18) || (x == 14) || (x == 10) || (x == 6) ||(x == 2)


extern u32 GIMLI_state_word[25][12]; // 24 round with 25 states.

extern u64 gimli_hwval_space[33][65536]; // hamming weight space with corresponding hw.
extern u64 gimli_hwspace_num[33];





/**
 * Compute the tuple or difference of x y z.
 *
 */
extern u32 gimli_OR_table[8][4];


/**
 * @Search for the Differential of GIMLI.
 * @ MOde:9
 */
u16 GIMLI_differential_trail_search_entry (u16 search_round);
void gimli_hw_space(void);  //GIMLI search hw_val space.
u16 gimli_round_1(u16 search_round);
u16 gimli_round_r(u16 search_round);
u16 gimli_round_N(u16 search_round);



#endif /* CIPHERS_GIMLI_H_ */
