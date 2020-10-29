/*
 * gloablvar.c
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 */

#include <stdio.h>
#include "ciphers.h"
#include "search.h"
#include "typedef.h"
#include "globalvar.h"
#include "printinfo.h"



/* The initialization input information of the cipher and mode in the program.*/
s8 str_cipher[10] = "simon";
s8 str_blocksize[10] = "32";
int sc_blocksize = 32;
u64 blocksize_len = 16;
u64 word_len = 16;
int sc_rounds = 0;
int sc_maxBnw= 2;
int str_mode = 0;

/* The cipher variant with block size is nBytes. */
u16 nBytes = 2;
/* The initialization input difference*/
volatile u64 input_diff_L = 0;  // Left part of input difference.
volatile u64 input_diff_R = 0;  // Right part of input difference.
volatile u64 output_diff_L = 0;  // Left part of output difference.
volatile u64 output_diff_R = 0;  // Right part of output difference.
volatile u64 input_mask_L = 0;  // Left part of input mask.
volatile u64 input_mask_R = 0;  // Right part of input mask.
volatile u64 output_mask_L = 0;  // Left part of output mask.
volatile u64 output_mask_R = 0;  // Right part of output mask.
volatile u64 x_in[100] = {0};
volatile u64 n_x_in[100] = {0};
volatile u64 y_in[100] = {0};
volatile u64 n_y_in[100] = {0};
volatile u64 beta_in[100] = {0};  // beta of r round.
volatile u64 n_beta_in[100] = {0};  // beta of r round.
volatile u16 P_w[64] = {0};  // the current best weight of current round i.
 u16 p_sumof_r_P_w[64]={0}; // p_sumof_r is the sum of r-1 rounds weight.

volatile u16 n_P_w[100] = {0};  // the current best weight of current round i.
volatile u16 n_P_bestofR_w[100] = {0};  //{2,0,2,4,6,8,12,14,18,20,25,30,34};// the best weight of n round.
volatile u16 Bn_w = 2;  // The pre set weight of r round which the best weight of r round tend to be not bigger than this.
volatile u16 DDTAB[256][256] = {0};  // DDTAB is the difference distribution table of AND operation within 4 bits.
u64 DDTAB_8bit[65536][257] = {0};  // DDTAB is the difference distribution table of AND operation within 8 bits.
u64 HMspace[16][65535] = {0}; // hamming weight space value with corresponding hw.
u64 HMspace_num[256] = {0};   // the number of corresponding weight.


u64 gama_cnt[256] = {0};


long double prob = 0.0;
long double prob_w = 0.0;
u64 count = 0;
u64 trail_nmu = 0;
u64 wt_nmu[256] = {0};
u64 wt_max = 0;
u64 wt_ctr_max = 0;

//u64 U_base_y[65] = {0};  // non zero vector in U_base construct the bases of Ub.//

u64 flag_stop=0;



////fro sepck///////////
volatile u16 P_sepck[100] = {0};  // pro. of each round of speck.
volatile u16 P_sepck_bit[100] = {0};  // pro. of each round of speck.


u64 Bit_Align = 0xFFFF;
u64 V_MSB = 0x8000;
u64 ValueMax_Align = 0x7FFF;


u64 cDDT_v[8][65536][8][256] = {0};
u64 cDDT_n[8][65536][9] = {0}; //0-8
u64 cDDT_wt_max[8][65536] = {0};
u64 cDDT_wt_min[8][65536] = {0};
u64 cDDT_AB_wt_min[65536] = {0};


u64 MSB_cDDT_v[8][65536][8][256] = {0};
u64 MSB_cDDT_n[8][65536][9] = {0}; //0-8
u64 MSB_cDDT_wt_max[8][65536] = {0};
u64 MSB_cDDT_wt_min[8][65536] = {0};


u64 bDDT_v[8][65536][8][256] = {0};
u64 bDDT_n[8][65536][9] = {0}; //0-8
u64 bDDT_wt_max[8][65536] = {0};
u64 bDDT_wt_min[8][65536] = {0};
u64 MSB_bDDT_v[8][65536][8][256] = {0};
u64 MSB_bDDT_n[8][65536][9] = {0}; //0-8
u64 MSB_bDDT_wt_max[8][65536] = {0};
u64 MSB_bDDT_wt_min[8][65536] = {0};


//char C0[65] = {0};
u16 set_A_4[4] = {0,3,5,6};
u16 set_A_3[3] = {3,5,6};
u16 set_B_4[4] = {1,2,4,7};
u16 set_B_3[3] = {1,2,4};
u16 set_AB_6[6] = {3,5,6,1,2,4};

u16 flag = 0;


///////////以下为线性逼近表的使用参数//////////////////
u16 cLAT_W[256][2][9][65536] = {0};
u16 cLAT_U[256][2][9][65536] = {0};
u16 cLAT_WU_numb[256][2][9]={0};
u16 cLAT_wtcor_min[256][2]={0};
u16 cLAT_UVW_bro[256][256][256][2] = {0};















