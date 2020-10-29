/*
 * RX_CHAM.h
 *
 *  Created on: 2020年9月8日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_CHAM_H_
#define RX_CRYPTANALYSIS_RX_CHAM_H_


#include <stdio.h>
#include "typedef.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include <nmmintrin.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
volatile u16 CHAM64_const_rxd[80]; //每轮的轮常数的RXD
volatile u32 CHAM128_const_rxd[80]; //每轮的轮常数的RXD

volatile float RX_n_P_bestofR_w[80];  //奇数轮的最优rxdp
volatile float CHAM_RX_Pw_even_wt[80]; //偶数轮的最优rxdp

volatile u32 CHAM_RXD[4][80]; //每轮的输入差分
volatile u32 CHAM_Zeta[80]; //每轮的输zeta
volatile u32 CHAM_zeta_tmp[80]; //每轮的临时zeta
volatile float CHAM_Pr_zeta[80]; //每轮RX差分概率
volatile float CHAM_Pr_zeta_tmp[80]; //每轮RX差分临时概率
volatile float CHAM_Pr_RX[80]; //每轮RX差分概率
volatile float CHAM_Pr_add[80]; //每轮RX差分概率


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Compute_CHAM64_const_rxd(void);
void Compute_CHAM128_const_rxd(void);

u32 CHAM_64_RX_trail_search_entry(u16 search_round);
u32 CHAM_128_RX_trail_search_entry(u16 search_round);
u32 CHAM_RX_trail_round_1(u16 search_round);
u16 CHAM_RX_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 CHAM_RX_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 CHAM_RX_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u32 CHAM_RX_trail_round_23(u16 search_round, u16 cur_round, u32 *input_x);
u32 CHAM_RX_trail_round_r(u16 search_round, u16 cur_round, u32 *input_x);
u32 CHAM_RX_trail_round_N(u16 search_round, u32 *input_x);


u16 print_CHAM_64_RXD_resoult(u16 search_round);
u32 print_CHAM_128_RXD_resoult(u16 search_round);


#endif /* RX_CRYPTANALYSIS_RX_CHAM_H_ */
