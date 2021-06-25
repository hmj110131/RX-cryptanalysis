/*
 * RX_Siphash.h
 *
 *  Created on: 2020年7月31日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_SIPHASH_H_
#define RX_CRYPTANALYSIS_RX_SIPHASH_H_

#include <stdio.h>
#include "typedef.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include "Alzette.h"
#include <nmmintrin.h>


////////////////////////////////搜索Siphash的最优RX特征的代码//////////////////////////////
volatile float RX_Siphash_P_hr_w[25]; //每半轮中，对应RX差分概率重量
volatile float RX_Siphash_P_hr_w_tmp[25]; //每轮对应RX差分概率重量,临时的
volatile float RX_Siphash_ADD_P_xor_w[2][25]; //每半轮中2个模块加法的对应RX差分概率重量
volatile float RX_Siphash_ADD_P_xor_w_tmp[2][25]; //每轮中，4个模块加法的对应RX差分概率重量
volatile u64 Siphash_R1_ADD_i_RX_abc[2][3]; //2个模块加法的对应a/b/c
volatile u64 Siphash_RX_X[4][25]; //每半轮中，4个分支的输入RX差分
volatile u64 Siphash_RX_X_tmp[4][25]; //每轮中，4个分支的输入RX差分
volatile u64 Siphash_Zeta[2][25]; //每半轮的输zeta
volatile u64 Siphash_zeta_tmp[2][25]; //每半轮的临时zeta
volatile float Siphash_Pr_zeta[2][25]; //每轮RX差分概率
volatile float Siphash_Pr_zeta_tmp[2][25]; //每轮RX差分概率



////////////////////////////////搜索Siphash的最优RX特征的代码//////////////////////////////
u64 Siphash_RX_trail_search_entry(u16 search_round);
u64 Siphash_RX_ADD_wt_inc(u16 search_round);
u64 Siphash_RX_ADD_i_tuples(u16 search_round,u16 ADD_part);
u64 RX_Siphash_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 ADD_part );
u64 RX_Siphash_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi,u16 ADD_part);
u64 RX_Siphash_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 ADD_part);
u64 Siphash_RX_comb_RX(u16 search_round);

u64 Siphash_RX_Middle_Rounds(u16 search_round, u16 current_round, u64 *Xr);


u64 Siphash_RX_N_Rounds(u16 search_round,u32 *Xr);






u64 print_Siphash_RX_resoult(u16 search_round);

#endif /* RX_CRYPTANALYSIS_RX_SIPHASH_H_ */
