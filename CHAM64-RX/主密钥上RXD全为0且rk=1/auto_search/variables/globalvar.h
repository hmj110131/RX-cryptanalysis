/*
 * globalvar.h
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 */
#ifndef VARIABLES_GLOBALVAR_H_
#define VARIABLES_GLOBALVAR_H_

#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include "search.h"
#include "typedef.h"



/* The initialization input information of the cipher and mode in the program.*/
extern s8 str_cipher[10],str_blocksize[10];
extern int str_mode;
extern int sc_blocksize;
extern u64 blocksize_len;
extern u64 word_len;
extern int sc_maxBnw;
extern int sc_rounds;

/* The cipher variant with block size is nBytes. */
extern u16 nBytes;
/* The input and output difference*/
extern volatile u64 input_diff_L;  // Left part of input difference.
extern volatile u64 input_diff_R;  // Right part of input difference.
extern volatile u64 output_diff_L;  // Left part of output difference.
extern volatile u64 output_diff_R;  // Right part of output difference.
extern volatile u64 input_mask_L;  // Left part of input mask.
extern volatile u64 input_mask_R;  // Right part of input mask.
extern volatile u64 output_mask_L;  // Left part of output mask.
extern volatile u64 output_mask_R;  // Right part of output mask.
extern volatile u64 x_in[100];  // x of r round.
extern volatile u64 y_in[100];  // y of r round.
extern volatile u64 beta_in[100];  // beta of r round.
extern volatile u64 n_x_in[100];  // x of r round.
extern volatile u64 n_y_in[100];  // y of r round.
extern volatile u64 n_beta_in[100];  // beta of r round.
extern volatile u16 P_w[64]; // the current best weight of current round i.
extern  u16 p_sumof_r_P_w[64]; // p_sumof_r is the sum of r-1 rounds weight.
extern volatile u16 n_P_w[100]; // the current best weight of current round n.
extern volatile u16 n_P_bestofR_w[100]; // the best weight of r round.
extern volatile u16 Bn_w;  // weight of expected n rounds which is dynamicly changed to close to the best weight.
extern volatile u16 DDTAB[256][256]; // DDTAB is the difference distribution table of AND operation within 4 bits.
extern  u64 DDTAB_8bit[65536][257];  // DDTAB is the difference distribution table of AND operation within 8 bits.
extern u64 HMspace[16][65535]; // hamming weight space with corresponding hw.
extern u64 HMspace_num[256];

extern u64 gama_cnt[256];


extern long double prob; // The total probability of differential Characteristic or Linear Hull.//
extern long double prob_w; // The weight of probability of differential Characteristic or Linear Hull.//
extern u64 count;   // The count number of the trails of differential Characteristic or Linear Hull.//
extern u64 trail_nmu;
extern u64 wt_nmu[256];
extern u64 wt_max;
extern u64 wt_ctr_max;
//extern u64 U_base_y[65];


extern u64 flag_stop;



////fro sepck///////////
extern volatile u16 P_sepck[100];  // pro. of each round of speck.
extern volatile u16 P_sepck_bit[100];  // pro. of each round of speck.


extern  u64 Bit_Align;
extern  u64 V_MSB ;
extern  u64 ValueMax_Align;

extern  u64 cDDT_v[8][65536][8][256];
extern  u64 cDDT_n[8][65536][9]; //0-8
extern  u64 cDDT_wt_max[8][65536];
extern  u64 cDDT_wt_min[8][65536];
extern  u64 cDDT_AB_wt_min[65536];
extern  u64 MSB_cDDT_v[8][65536][8][256];
extern  u64 MSB_cDDT_n[8][65536][9]; //0-8
extern  u64 MSB_cDDT_wt_max[8][65536];
extern  u64 MSB_cDDT_wt_min[8][65536];
extern u64 bDDT_v[8][65536][8][256];
extern u64 bDDT_n[8][65536][9];
extern u64 bDDT_wt_max[8][65536];
extern u64 bDDT_wt_min[8][65536];
extern u64 MSB_bDDT_v[8][65536][8][256];
extern u64 MSB_bDDT_n[8][65536][9];
extern u64 MSB_bDDT_wt_max[8][65536] ;
extern u64 MSB_bDDT_wt_min[8][65536] ;


//extern char C0[65];

extern u16 set_A_4[4] ;
extern u16 set_A_3[3] ;
extern u16 set_B_4[4] ;
extern u16 set_B_3[3] ;
extern u16 set_AB_6[6];

extern u16 flag;

///////////以下为线性逼近表的使用参数//////////////////
extern u16 cLAT_W[256][2][9][65536];
extern u16 cLAT_U[256][2][9][65536];
extern u16 cLAT_WU_numb[256][2][9];
extern u16 cLAT_wtcor_min[256][2];
extern u16 cLAT_UVW_bro[256][256][256][2];




#endif /* VARIABLES_GLOBALVAR_H_ */
