/*
 * RX_Chaskey.h
 *
 *  Created on: 2020年7月29日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_CHASKEY_H_
#define RX_CRYPTANALYSIS_RX_CHASKEY_H_


#include <stdlib.h>
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include "printinfo.h"
#include "RX_Chaskey.h"
#include "XOR_offset_table.h"
#include <nmmintrin.h>
#include "chaskey.h"



////////////////////////////////搜索Chaskey的最优RX特征的代码//////////////////////////////
volatile float RX_Chaskey_P_hr_w[25]; //每半轮中，对应RX差分概率重量
volatile float RX_Chaskey_P_hr_w_tmp[25]; //每轮对应RX差分概率重量,临时的
volatile float RX_Chaskey_ADD_P_xor_w[2][25]; //每半轮中2个模块加法的对应RX差分概率重量
volatile float RX_Chaskey_ADD_P_xor_w_tmp[2][25]; //每轮中，4个模块加法的对应RX差分概率重量
volatile u32 Chaskey_R1_ADD_i_RX_abc[2][3]; //2个模块加法的对应a/b/c
volatile u32 RX_X[4][25]; //每半轮中，4个分支的输入RX差分
volatile u32 RX_X_tmp[4][25]; //每轮中，4个分支的输入RX差分
volatile u32 Chaskey_Zeta[2][25]; //每半轮的输zeta
volatile u32 Chaskey_zeta_tmp[2][25]; //每半轮的临时zeta
volatile float Chaskey_Pr_zeta[2][25]; //每轮RX差分概率
volatile float Chaskey_Pr_zeta_tmp[2][25]; //每轮RX差分概率


////////////////////////////////搜索Chaskey的最优RX特征的代码//////////////////////////////
u32 Chaskey_RX_trail_search_entry(u16 search_round);
u32 Chaskey_RX_ADD_wt_inc(u16 search_round);
u32 Chaskey_RX_ADD_i_tuples(u16 search_round,u16 ADD_part);
u16 RX_Chaskey_input_MSB
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi,u16 ADD_part );
u16 RX_Chaskey_input_Middle
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi, u16 cur_posi,u16 ADD_part);
u16 RX_Chaskey_input_Last
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi,u16 ADD_part);
u32 Chaskey_RX_comb_RX(u16 search_round);

u32 Chaskey_RX_Middle_Rounds(u16 search_round, u16 current_round, u32 *Xr);
u32 Chaskey_RX_N_Rounds(u16 search_round,u32 *Xr);

u16 print_Chaskey_RX_resoult(u16 search_round);

////////////////////////////////以下为Chaskey的RX差分特征验证，实验代码///////////
u32 ADD_RX_check(u32 rk);
u32 ADD_RX_delta_check(u32 X0, u32 X1, u32 Y0, u32 Y1, u32 rk);




#endif /* RX_CRYPTANALYSIS_RX_CHASKEY_H_ */
