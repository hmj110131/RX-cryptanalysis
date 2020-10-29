/*
 * Alzette.h
 *
 *  Created on: 2020年4月22日
 *      Author: Mingjiang Huang
 */

#ifndef CIPHERS_ALZETTE_H_
#define CIPHERS_ALZETTE_H_

#include <stdio.h>
#include "typedef.h"
#include <time.h>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
#define Alzette_rounds  4     //宏定义，r=4

extern u16  Alzette_rol_a[16];     //定义每轮的循环参数
extern u16  Alzette_rol_b[16];     //定义每轮的循环参数
extern u32  Alzette_Constant[8];     //定义Alzette的实例版本，共8个版本

volatile u16 Alzette_Pw_cor[16]; // 记录Alzette的最优线性相关性重量。
volatile u16 Alzette_Pw_wt[16]; // 记录Alzette的最优差分概率重量

volatile u32 Alzette_X[2][16]; //每轮的输入差分左右两部分
volatile u32 Alzette_Alpha[16];  //每轮的第模加的输入输出差分
volatile u32 Alzette_Beta[16];
volatile u32 Alzette_Gamma[16];
volatile u32 Alzette_u[16];  //每轮的第模加的输入输出掩码 //output u
volatile u32 Alzette_v[16];  //input v
volatile u32 Alzette_w[16];  //input w

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//Alzette每轮的用到的循环移位操作,all are rotated to right
u32 rol_right_32(u32 x, u32 indx);
u32 rol_left_32(u32 x, u32 indx);
u32 Alzette_rol_right_31(u32 input);
u32 Alzette_rol_left_31(u32 input);
u32 Alzette_rol_right_24(u32 input);
u32 Alzette_rol_right_17(u32 input);
//u32 Alzette_rol_right_17(u32 input);
u32 Alzette_rol_right_0(u32 input);
//u32 Alzette_rol_right_31(u32 input);
//u32 Alzette_rol_right_24(u32 input);
u32 Alzette_rol_right_16(u32 input);


////////打印得到的差分和线性最优路径
void print_Alzette_diff_resoult(u16 search_round);
void print_Alzette_linear_resoult(u16 search_round);


//////////////////////////search for differential/////////////////////////////
u32 Alzette_Diff_trail_search_entry(u16 search_round);
u32 Alzette_Diff_trail_round_1(u16 search_round);
u16 Alzette_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 Alzette_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 Alzette_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u32 Alzette_Diff_trail_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);
u32 Alzette_Diff_trail_round_N(u16 search_round, u64 x, u64 y);

//////////////////////////search for linear/////////////////////////////
u32 Alzette_Linear_trail_search_entry(u16 search_round);
u32 Alzette_Linear_trail_round_1(u16 search_round);
u16 Alzette_ADD_mask_MSB
(u16 search_round,u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 Alzette_ADD_mask_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc,char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 Alzette_ADD_mask_Last(u16 search_round,u16 Corr_1, char *posi );
u32 Alzette_Linear_trail_round_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W);
u16 Alzette_ADD_mask_Last_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W,
		u16 Corr_2, char *posi );
u16 Alzette_ADD_mask_Middle_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 Alzette_ADD_mask_MSB_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u32 Alzette_Linear_trail_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);
u32 Alzette_Linear_round_N(u16 search_round, u64 x, u64 y);







#endif /* CIPHERS_ALZETTE_H_ */
