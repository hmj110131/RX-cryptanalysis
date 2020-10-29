/*
 * hight.h
 *
 *  Created on: 2018年8月13日
 *      Author: hmj110131
 */

#ifndef CIPHERS_HIGHT_H_
#define CIPHERS_HIGHT_H_


#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include "search.h"
#include "typedef.h"



extern u8 x_block[33][8];
extern u8 y_block[33][8];

extern u8 F_0_Table[256];
extern u8 F_1_Table[256];
extern u8 wt_F_0_Table[256];
extern u8 wt_F_1_Table[256];

extern u8 wt_x_F1[256][256];

extern u8 Var_F_0_Table[256];
extern u8 Var_F_1_Table[256];

extern u32 Input_Comb_cnt[8];
extern u8 Input_Comb_alpha[8][1400000];
extern u8 Input_Comb_beta[8][1400000];
extern u8 Input_Comb_gamma[8][1400000];

extern u32 fixed_alpha_wt_cnt[256][8];
extern u8 fixed_alpha_wt_beta[256][8][256];
extern u8 fixed_alpha_wt_gamma[256][8][256];
extern u16 fixed_alpha_wt_min[256];
extern u16 fixed_alpha_wt_max[256];

extern u16 Var_XOR_cnt[256];
extern u8 Var_XOR_A[256][256];
extern u8 Var_XOR_B[256][256];

extern u8 P0_wt_HIGHT[32];
extern u8 P1_wt_HIGHT[32];
extern u8 P2_wt_HIGHT[32];
extern u8 P3_wt_HIGHT[32];

extern u8 wt_inc_T[8][128];
extern u8 wt_inc_T_cnt[8];


extern u8 HIGHT_cDDT_v[65536][8][256];
extern u8 HIGHT_cDDT_n[65536][8];
extern u8 HIGHT_cDDT_wt_max[65536];
extern u8 HIGHT_cDDT_wt_min[65536];

extern u8 wt_3_n[32][12000000];
extern u8 wt_2_n[32][12000000] ;
extern u8 wt_1_n[32][12000000];
extern u8 wt_0_n[32][12000000];
extern u64 wt_n_max[32];


extern u8 HIGHT_x[8];
extern u8 HIGHT_y[8];


///////////////////////////
u8 rol_left(u8 x, u8 indx);


u8 hight_F_0(u8 x);
u8 hight_F_1(u8 x);

u8 Construct_F_0(void);
u8 Construct_F_1(void);

u8 Construct_Var_F_0(void);
u8 Construct_Var_F_1(void);

u8 Construct_Input_Comb(void);
u8 Construct_Input_Comb_MSB(u8 alpha, u8 beta, u8 gamma,u8 P_1, char *posi );
u8 Construct_Input_Comb_Middle(u8 alpha,u8 beta,u8 gamma,u8 P_1,char *posi,u8 cur_posi);
u8 Construct_Input_Comb_Last(u8 alpha, u8 beta, u8 gamma,u8 P_1, char *posi );

u8 Construct_fixedAlpha_betagamma_Comb(void);

u8 Construct_XOR_var_Table(void);

u8 Construct_wt_inc_Table(void);

void ARX_carry_DDTm_8bit(void);

void HIGHT_first_wt_Comb(void);

u8 Construct_HIGHT_Tables(void);

u8 print_HIGHT_resoult(u16 search_round);





/////////////////////////////////////////////
u8 HIGHT_round_1(u16 search_round);
u8 HIGHT_round_1_diff(u16 search_round);
u8 HIGHT_round_2(u16 search_round,u8 x7,u8 x6,u8 x5, u8 x4, u8 x3, u8 x2,u8 x1,u8 x0);
u8 HIGHT_round_r(u16 search_round,u16 cur_round,u8 x7, u8 x6, u8 x5, u8 x4,u8 x3,u8 x2,u8 x1,u8 x0);
u8 HIGHT_round_N(u16 search_round,u8 x7, u8 x6, u8 x5, u8 x4,u8 x3,u8 x2,u8 x1,u8 x0);
u8 HIGHT_round_N_diff(u16 search_round,u8 x7, u8 x6, u8 x5, u8 x4,u8 x3,u8 x2,u8 x1,u8 x0);
u8 HIGHT_round_r_diff(u16 search_round,u16 cur_round,u8 x7, u8 x6, u8 x5, u8 x4,u8 x3,u8 x2,u8 x1,u8 x0);
u16 HIGHT_differential_trail_search_entry (u16 search_round);
u16 HIGHT_diff_search_entry(u16 search_round);  // 默认只搜索Hight64的差分












#endif /* CIPHERS_HIGHT_H_ */
