/*
 * CHAM.h
 *
 *  Created on: 2019年4月29日
 *      Author: hmj110131
 */

#ifndef CIPHERS_CHAM_H_
#define CIPHERS_CHAM_H_

#include <stdio.h>
#include "typedef.h"

//1从偶数轮开始，0从奇数轮开始,默认CHAM差分特征是从奇数轮开始的
#define CHAM_even_round  0     //宏定义，预编译搜索的是奇数还是偶数轮的CHAM 0/1,
#define CHAM_64_or_128 64    //宏定义，预编译搜索的是CHAM64还是CHAM128 64/128




extern  u32 CHAM_even;  //是否搜索的是CHAM的从偶数轮开始的截断最优路径
volatile u16 CHAM_Pw_even_cor[80]; // 记录CHAM的从偶数轮开始的截断最优路径，的最优线性相关性重量。
volatile u16 CHAM_Pw_even_wt[80]; // 记录CHAM的从偶数轮开始的截断最优差分路径，的最优差分概率重量

volatile u16 CHAM_first_four_wt[4];  //前4轮的相关性重量

extern volatile u32 CHAM_Alpha[80];  //每轮的第模加的输入输出差分
extern volatile u32 CHAM_Beta[80];
extern volatile u32 CHAM_Gamma[80];
extern volatile u32 CHAM_X[4][80]; //每轮的输入差分

extern volatile u16 CHAM_add_cw[80]; //每轮的模加的相关性重量
extern volatile u16 CHAM_64_u[80];  //每轮的第模加的输入输出掩码
extern volatile u16 CHAM_64_v[80];
extern volatile u16 CHAM_64_w[80];
extern volatile u16 CHAM_64_X[4][80]; //每轮的输入掩码

extern volatile u32 CHAM_128_u[80];  //每轮的第模加的输入输出掩码
extern volatile u32 CHAM_128_v[80];
extern volatile u32 CHAM_128_w[80];
extern volatile u32 CHAM_128_X[4][80]; //每轮的输入掩码


extern u16 CHAM_64_IN_Mask[4];
extern u16 CHAM_64_OUT_Mask[4];
extern u32 CHAM_128_IN_Mask[4];
extern u32 CHAM_128_OUT_Mask[4];

extern u8 CHAM_a_beta[256][8][9][32768];  //alpha carry wt num
extern u8 CHAM_a_gama[256][8][9][32768];  //alpha carry wt num
extern u32 CHAM_a_bg_numb[256][8][9];   // alpha carry wt
extern u8 CHAM_a_wt_min[256][8];   // alpha carry
extern u8 CHAM_a_wt_max[256][8];   // alpha carry

extern u8 msb_CHAM_a_beta[256][8][9][32768];  //alpha carry wt num
extern u8 msb_CHAM_a_gama[256][8][9][32768];  //alpha carry wt num
extern u32 msb_CHAM_a_bg_numb[256][8][9];   // alpha carry wt
extern u8 msb_CHAM_a_wt_min[256][8];   // alpha carry
extern u8 msb_CHAM_a_wt_max[256][8];   // alpha carry

void fixed_Alpha_get_betagamma(void);






////////////////////////////////////////////////////////////
u32 CHAM_64_Diff_trail_search_entry(u16 search_round);
u32 CHAM_128_Diff_trail_search_entry(u16 search_round);

u32 CHAM_Diff_trail_round_1(u16 search_round);
u16 CHAM_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 CHAM_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 CHAM_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u32 CHAM_Diff_trail_round_23(u16 search_round, u16 cur_round, u32 *input_x);
u32 CHAM_Diff_trail_round_r(u16 search_round, u16 cur_round, u32 *input_x);
u32 CHAM_Diff_trail_round_N(u16 search_round, u32 *input_x);






////////////////////////////////////////////////////////////
u32 CHAM_64_Linear_trail_search_entry(u16 search_round);
u32 CHAM_R1_tuples(u16 search_round);
u16 CHAM_R1_tuples_ADD_Last(u16 search_round, u16 Cw_1, char *posi );
u16 CHAM_R1_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_R1_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_R2_tuples(u16 search_round);
u16 CHAM_R2_tuples_ADD_Last(u16 search_round,u16 Cw_2, char *posi );
u16 CHAM_R2_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_R2_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_R3_tuples(u16 search_round);
u16 CHAM_R3_tuples_ADD_Last(u16 search_round,u16 Cw_3, char *posi );
u16 CHAM_R3_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_3, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_R3_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_3, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_R4_tuples(u16 search_round);
u16 CHAM_R4_tuples_ADD_Last(u16 search_round,u16 Cw_4, char *posi );
u16 CHAM_R4_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_4, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_R4_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_4, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_First_4_round_masks_comb(u16 search_round);
u32 CHAM_round_r_odd(u16 search_round, u16 cur_round,u16 *Mask_r);
u32 CHAM_round_r_even(u16 search_round, u16 cur_round,u16 *Mask_r);
u32 CHAM_round_N_odd(u16 search_round,u16 *Mask_r);
u32 CHAM_round_N_even(u16 search_round,u16 *Mask_r);

///////////
u16 CHAM64_linear_Hull_search_entry (u16 search_round);
u32 CHAM64_linear_Hull_round_r_odd(u16 search_round, u16 cur_round,u16 *Mask_r);
u32 CHAM64_linear_Hull_round_r_even(u16 search_round, u16 cur_round,u16 *Mask_r);
u32 CHAM64_linear_Hull_round_N_odd(u16 search_round,u16 *Mask_r);
u32 CHAM64_linear_Hull_round_N_even(u16 search_round,u16 *Mask_r);

u16 print_CHAM_64_resoult(u16 search_round);
u16 print_CHAM_64_linear_resoult(u16 search_round);



////////////////////////////////////////////////////////////
u32 CHAM_128_Linear_trail_search_entry(u16 search_round);
u32 CHAM_128_R1_R4_wt(u16 search_round, u16 F4_thrd);
u32 CHAM_128_R1_tuples(u16 search_round);
u16 CHAM_128_R1_tuples_ADD_Last(u16 search_round, u16 Cw_1, char *posi );
u16 CHAM_128_R1_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_128_R1_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_128_R2_tuples(u16 search_round);
u16 CHAM_128_R2_tuples_ADD_Last(u16 search_round,u16 Cw_2, char *posi );
u16 CHAM_128_R2_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_128_R2_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_128_R3_tuples(u16 search_round);
u16 CHAM_128_R3_tuples_ADD_Last(u16 search_round,u16 Cw_3, char *posi );
u16 CHAM_128_R3_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_3, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_128_R3_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_3, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_128_R4_tuples(u16 search_round);
u16 CHAM_128_R4_tuples_ADD_Last(u16 search_round,u16 Cw_4, char *posi );
u16 CHAM_128_R4_tuples_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_4, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 CHAM_128_R4_tuples_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Cw_4, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

u32 CHAM_128_First_4_round_masks_comb(u16 search_round);
u32 CHAM_128_round_r(u16 search_round, u16 cur_round,u32 *Mask_r);
u32 CHAM_128_round_N(u16 search_round,u32 *Mask_r);

u32 CHAM128_Linear_Hull_search_entry(u16 search_round);
u32 CHAM128_Linear_Hull_round_r_odd(u16 search_round, u16 cur_round,u32 *Mask_r);
u32 CHAM128_Linear_Hull_round_r_even(u16 search_round, u16 cur_round,u32 *Mask_r);
u32 CHAM128_Linear_Hull_round_N_odd(u16 search_round,u32 *Mask_r);
u32 CHAM128_Linear_Hull_round_N_even(u16 search_round,u32 *Mask_r);





////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
//////////////////////之前的搜索剪枝条件有问题//////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////







u32 print_CHAM_128_resoult(u16 search_round);
u32 print_CHAM_128_linear_resoult(u16 search_round);

#endif /* CIPHERS_CHAM_H_ */
