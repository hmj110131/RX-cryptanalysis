/*
 * RX_speck_key.h
 *
 *  Created on: 2020年9月9日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_SPECK_KEY_H_
#define RX_CRYPTANALYSIS_RX_SPECK_KEY_H_

#include <stdio.h>
#include "typedef.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include <nmmintrin.h>


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
volatile u64 Speck_key_const_rxd[20]; //每轮的密钥编排部分的轮常数的RXD
volatile u64 speck_key_RXD_tmp[4][20]; //密钥部分，每轮的输入差分RXD，实际从第2轮开始
volatile u64 speck_key_RXD[4][20]; //密钥部分，每轮的输入差分左右两部分RXD，实际从第2轮开始

volatile u64 speck_Zeta[16]; //每轮的输zeta
volatile u64 speck_zeta_tmp[16]; //每轮的临时zeta
volatile float speck_Pr_RX[16]; //每轮RX差分概率
volatile float speck_Pr_XOR[16]; //每轮RX差分概率
volatile float speck_Pr_zeta[16]; //每轮RX差分概率
volatile float speck_Pr_zeta_tmp[16]; //每轮RX差分概率


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//对于Speck的密钥编排部分的RXD路径的搜索，分为3种情况，即对应：
//speck2n/2n    speck2n/3n    speck2n/4n
//密钥部分字节数m=2/3/4的情况，搜索程序的前几轮略有不同
#define speck_m_keyword  3 ////speck2n/2n: 2  speck2n/3n: 3   speck2n/4n: 4
#define speck_m_bits  144 ////speck的密钥比特长度和密钥的字个数需对应
#define speck_keyword_len  speck_m_bits/speck_m_keyword ////speck的密钥比特的字长度
//speck32:  64     m=4
//speck48:  72,96  m=3/4
//speck64:  96,128   m=3/4
//speck96:  96,144     m=2/3
//speck128:  128,192,256   m=2/3/4
void Compute_RX_speck_key_const_rxd(void);  //不用考虑分组长度，轮常数的RXD都不会影响
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
u16 RX_sepck_key_trail_search_entry (u16 search_round);
/////////
#if (speck_m_keyword  == 2)  //speck2n/2n: 2
//%%%%//speck2n/2n %%%%//
u16 RX_sepck_2n2n_key_round_1(u16 search_round);
u16 RX_sepck_2n2n_key_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_2n2n_key_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 RX_sepck_2n2n_key_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_2n2n_key_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 RX_speck_2n2n_key_round_N(u16 search_round, u64 x, u64 y);


#elif (speck_m_keyword  == 3)  //speck2n/3n: 3
//%%%%//speck2n/3n %%%%//
u16 RX_sepck_2n3n_key_round_1(u16 search_round);
u16 RX_sepck_2n3n_key_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_2n3n_key_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 RX_sepck_2n3n_key_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_2n3n_key_round_2(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 RX_sepck_2n3n_key_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 RX_speck_2n3n_key_round_N(u16 search_round, u64 x, u64 y);

#elif (speck_m_keyword  == 4)  //speck2n/4n: 4
//%%%%//speck2n/4n %%%%//
u16 RX_sepck_2n4n_key_round_1(u16 search_round);
u16 RX_sepck_2n4n_key_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_2n4n_key_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 RX_sepck_2n4n_key_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_2n4n_key_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 RX_speck_2n4n_key_round_N(u16 search_round, u64 x, u64 y);

#endif










//%%%%%%%%%%%//
void print_RX_speck_key_resoult(u16 search_round);












#endif /* RX_CRYPTANALYSIS_RX_SPECK_KEY_H_ */
