/*
 * sparx.h
 *
 *  Only for SPARX-64.
 *  Created on: 2018年8月28日
 *      Author: hmj110131
 */

#ifndef CIPHERS_SPARX_H_
#define CIPHERS_SPARX_H_

//Only for SPARX-64.
#include <stdio.h>
#include "typedef.h"


/////////用于SPARX-64的变量
extern u16 Comb_16bit_cnt[16];
extern u16 Comb_16bit_alpha[16][12000];
extern u16 Comb_16bit_beta[16][12000];
extern u16 Comb_16bit_gamma[16][12000];

extern u16 X_r[4][20];
extern u16 wt_l[20];
extern u16 wt_r[20];

extern u16 X_in_0, X_in_1, X_in_2, X_in_3;
extern u16 X_in_4, X_in_5, X_in_6, X_in_7;
extern u16 Y_out_0, Y_out_1, Y_out_2, Y_out_3;
extern u16 Y_out_4, Y_out_5, Y_out_6, Y_out_7;

extern u16 IN_X_Mask[8];
extern u16 OUT_Y_Mask[8];


extern u16 IN_L_X[5][4];

/////////用于SPARX-128的变量
extern u16 X_128_r[8][20];
extern u16 wt_block[4][20];
extern u16 IN128_L_X[5][8];
extern u16 FIRST;
//用于SPARX128的线性模加的输入输出掩码
extern  u16 R1_mask_UVW[12];
extern  u16 R1_mask_UVW_wt[4];
extern  u16 R2_mask_UVW[12];
extern  u16 R2_mask_UVW_wt[4];
extern u16 R1_wt_part[4];
extern u16 R2_wt_part[4];


u8 L2_Switch_per3Round(u16 *x0, u16 *x1, u16 *x2, u16 *x3 );

u16 SPARX64_round_1_Left(u16 search_round);
u16 SPARX64_input_MSB_Left(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi ,u16 P_1_r);
u16 SPAXR64_input_Middle_Left(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi,u16 P_1_r);
u16 SPARX64_input_Last_Left(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 P_1_r);
u16 SPARX64_round_1_Right(u16 search_round, u16 x0, u16 x1, u16 P_1_l,u16 P_1_r);
u16 SPARX64_input_MSB_Right(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 x0,u16 x1,u16 P_1_l);
u16 SPAXR64_input_Middle_Right(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi,u16 x0,u16 x1,u16 P_1_l);
u16 SPARX64_input_Last_Right(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 x0,u16 x1,u16 P_1_l);
u16 SPARX64_round_1_diff(u16 search_round);
u16 SPARX64_round_2to6(u16 search_round,u16 cur_round,u16 x,u16 y);
u16 SPARX64_round_2toR(u16 search_round,u16 cur_round,u16 x0,u16 x1,u16 x2,u16 x3);
u16 SPARX64_round_2toR_diff(u16 search_round,u16 cur_round,u16 x0,u16 x1,u16 x2,u16 x3);
u16 SPARX64_round_7toR(u16 search_round,u16 cur_round,u16 x0,u16 x1,u16 x2,u16 x3);
u16 SPARX64_round_N(u16 search_round,u16 x0,u16 x1,u16 x2,u16 x3);
u16 SPARX64_round_N_diff(u16 search_round,u16 x0,u16 x1,u16 x2,u16 x3);
u16 SPARX64_differential_trail_search_entry(u16 search_round);
u16 SPARX64_diff_search_entry(u16 search_round);



//用于SPARX-128的函数
u8 L4_Switch_per4Round(u16 *x0, u16 *x1, u16 *x2, u16 *x3, u16 *x4, u16 *x5, u16 *x6, u16 *x7);

u16 SPARX_128_differential_trail_search_entry(u16 search_round);
u16 SPARX_128_diff_search_entry(u16 search_round);

u16 SPARX_128_round_1(u16 search_round);

u16 SPARX_128_round_1_BLOCK_0(u16 search_round,
		u16 w_B_0, u16 w_B_1,u16 w_B_2,u16 w_B_3);
u16 SPARX_128_round_1_BLOCK_0_MSB(u16 search_round,
		u16 w_B_0, u16 w_B_1,u16 w_B_2,u16 w_B_3,
		char *posi);
u16 SPARX_128_round_1_BLOCK_0_Middle(u16 search_round,
		u16 w_B_0, u16 w_B_1,u16 w_B_2,u16 w_B_3,
		u16 alpha, u16 beta, u16 gamma,
		char *posi,u16 cur_posi);
u16 SPARX_128_round_1_BLOCK_0_Last(u16 search_round,
		u16 w_B_0, u16 w_B_1,u16 w_B_2,u16 w_B_3,
		u16 alpha, u16 beta, u16 gamma,
		char *posi);


u16 SPARX_128_round_1_BLOCK_1(u16 search_round,
		u16 w_B_1,u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1);
u16 SPARX_128_round_1_BLOCK_1_MSB(u16 search_round,
		u16 w_B_1,u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,
		char *posi);
u16 SPARX_128_round_1_BLOCK_1_Middle(u16 search_round,
		u16 w_B_1,u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,
		u16 alpha, u16 beta, u16 gamma,
		char *posi,u16 cur_posi);
u16 SPARX_128_round_1_BLOCK_1_Last(u16 search_round,
		u16 w_B_1,u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,
		u16 alpha, u16 beta, u16 gamma,
		char *posi);


u16 SPARX_128_round_1_BLOCK_2(u16 search_round,
		u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3);
u16 SPARX_128_round_1_BLOCK_2_MSB(u16 search_round,
		u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		char *posi);
u16 SPARX_128_round_1_BLOCK_2_Middle(u16 search_round,
		u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		u16 alpha, u16 beta, u16 gamma,
		char *posi,u16 cur_posi);
u16 SPARX_128_round_1_BLOCK_2_Last(u16 search_round,
		u16 w_B_2,u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		u16 alpha, u16 beta, u16 gamma,
		char *posi);


u16 SPARX_128_round_1_BLOCK_3(u16 search_round,
		u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,u16 x_in_4, u16 x_in_5);
u16 SPARX_128_round_1_BLOCK_3_MSB(u16 search_round,
		u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,u16 x_in_4, u16 x_in_5,
		char *posi);
u16 SPARX_128_round_1_BLOCK_3_Middle(u16 search_round,
		u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,u16 x_in_4, u16 x_in_5,
		u16 alpha, u16 beta, u16 gamma,
		char *posi,u16 cur_posi);
u16 SPARX_128_round_1_BLOCK_3_Last(u16 search_round,
		u16 w_B_3,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,u16 x_in_4, u16 x_in_5,
		u16 alpha, u16 beta, u16 gamma,
		char *posi);


u16 SPARX_128_round_1_diff(u16 search_round);

u16 SPARX_128_round_r(u16 search_round, u16 cur_round,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		u16 x_in_4, u16 x_in_5,u16 x_in_6, u16 x_in_7);

u16 SPARX_128_round_r_diff(u16 search_round, u16 cur_round,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		u16 x_in_4, u16 x_in_5,u16 x_in_6, u16 x_in_7);


u16 SPARX_128_round_N(u16 search_round,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		u16 x_in_4, u16 x_in_5,u16 x_in_6, u16 x_in_7);
u16 SPARX_128_round_N_diff(u16 search_round,
		u16 x_in_0, u16 x_in_1,u16 x_in_2, u16 x_in_3,
		u16 x_in_4, u16 x_in_5,u16 x_in_6, u16 x_in_7);

/////////////////////////////////线性特征搜索
u16 SPARX64_linear_trail_search_entry (u16 search_round);
//
u16 SPARX64_linear_round_1_Left(u16 search_round);
u16 SPARX64_R1_mask_Last_left(u16 search_round,u16 Corr_1, char *posi );
u16 SPARX64_R1_mask_Middle_left(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 SPARX64_R1_mask_MSB_left(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
//
u16 SPARX64_linear_round_1_Right(u16 search_round, u16 w1_left,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left );
u16 SPARX64_R1_mask_Last_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u16 Corr_1, char *posi );
u16 SPARX64_R1_mask_Middle_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 SPARX64_R1_mask_MSB_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
//
u16 SPARX64_linear_round_2_Left(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right);
u16 SPARX64_R2_mask_Last_left(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u16 Corr_1, char *posi);
u16 SPARX64_R2_mask_Middle_left(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 SPARX64_R2_mask_MSB_left(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
//
u16 SPARX64_linear_round_2_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u64 tmp_U_left_2, u64 tmp_V_left_2, u64 tmp_W_left_2,
		u16 w2_left);
u16 SPARX64_R2_mask_Last_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u64 tmp_U_left_2, u64 tmp_V_left_2, u64 tmp_W_left_2,
		u16 Corr_1, char *posi );
u16 SPARX64_R2_mask_Middle_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u64 tmp_U_left_2, u64 tmp_V_left_2, u64 tmp_W_left_2,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 SPARX64_R2_mask_MSB_Right(u16 search_round,
		u64 tmp_U_left,u64 tmp_V_left,u64 tmp_W_left,
		u64 tmp_U_right,u64 tmp_V_right,u64 tmp_W_right,
		u64 tmp_U_left_2, u64 tmp_V_left_2, u64 tmp_W_left_2,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);


u16 SPARX64_linear_round_r(u16 search_round, u16 cur_round,
		u64 X0, u64 X1, u64 X2, u64 X3);


u16 SPARX64_linear_round_N(u16 search_round,
		u64 X0, u64 X1, u64 X2, u64 X3);


u16 SPARX64_linear_Hull_search_entry (u16 search_round);
u16 SPARX64_linear_Hull_round_r(u16 search_round, u16 cur_round,
		u64 X0, u64 X1, u64 X2, u64 X3);
u16 SPARX64_linear_Hull_round_N(u16 search_round, u64 X0, u64 X1, u64 X2, u64 X3);





u16 SPARX128_linear_trail_search_entry (u16 search_round);
u16 SPARX128_linear_R1_wt(u16 search_round);
u16 SPARX128_linear_R2_wt(u16 search_round);
///
u16 SPARX128_linear_trail_R1 (u16 search_round, u16 part);
u16 SPARX128_R1_mask_Last(u16 search_round,
		u16 part, u16 part_wt,
		char *posi );
u16 SPARX128_R1_mask_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 part, u16 part_wt, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 SPARX128_R1_mask_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 part, u16 part_wt, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

///
u16 SPARX128_linear_trail_R2 (u16 search_round, u16 part);
u16 SPARX128_R2_mask_Last(u16 search_round,
		u16 part, u16 part_wt,
		char *posi );
u16 SPARX128_R2_mask_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 part, u16 part_wt, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 SPARX128_R2_mask_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 part, u16 part_wt, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);

///
u16 SPARX128_linear_trail_r (u16 search_round, u16 cur_round, u16 *mask_Xr);
u16 SPARX128_linear_trail_N(u16 search_round, u16 *mask_Xr );


u16 SPARX128_linear_Hull_search_entry (u16 search_round);
u16 SPARX128_linear_Hull_r (u16 search_round, u16 cur_round, u16 *mask_Xr );
u16 SPARX128_linear_Hull_N(u16 search_round, u16 *mask_Xr );






u16 print_SPARX64_resoult(u16 search_round);
u16 print_SPARX_128_resoult(u16 search_round);








#endif /* CIPHERS_SPARX_H_ */
