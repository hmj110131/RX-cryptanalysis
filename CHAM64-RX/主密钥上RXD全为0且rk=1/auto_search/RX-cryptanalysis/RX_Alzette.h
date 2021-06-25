/*
 * RX_Alzette.h
 *
 *  Created on: 2020年7月24日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_ALZETTE_H_
#define RX_CRYPTANALYSIS_RX_ALZETTE_H_


#include <stdio.h>
#include "typedef.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include "Alzette.h"
#include <nmmintrin.h>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
volatile u32 Alzette_RX[2][16]; //每轮的输入RX差分左右两部分
volatile u32 Alzette_Zeta[16]; //每轮的输zeta
volatile u32 Alzette_zeta_tmp[16]; //每轮的临时zeta
volatile float Alzette_Pr_RX[16]; //每轮RX差分概率
volatile float Alzette_Pr_add[16]; //每轮RX差分概率
volatile float Alzette_Pr_zeta[16]; //每轮RX差分概率
volatile float Alzette_Pr_zeta_tmp[16]; //每轮RX差分概率


////////////////////////////////计算Alzette的 8个版本 常数 的RX-差分固定值/////////////////////////
#define Alzette_C_i  0 //共8个版本，旋转0--7
extern u32  Alzette_C_RX[8];     //定义Alzette的实例版本常数RX差分，共8个版本
u32 Alzette_Constant_RX(u16 rx_k);







//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
////////////////////////////////搜索Alzette的最优 RX-差分 特征的代码/////////////////////////
u32 Alzette_RX_Diff_trail_search_entry(u16 search_round);

u32 Alzette_RX_Diff_trail_round_1(u16 search_round);
u16 RX_Alzette_input_MSB
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi );
u16 RX_Alzette_input_Middle
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 RX_Alzette_input_Last
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi );


u32 Alzette_RX_Diff_trail_round_r(u16 search_round, u16 cur_round, u32 x, u32 y);
u32 Alzette_RX_Diff_trail_round_N(u16 search_round, u32 x, u32 y);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//



void print_Alzette_RX_diff_resoult(u16 search_round);





#endif /* RX_CRYPTANALYSIS_RX_ALZETTE_H_ */
