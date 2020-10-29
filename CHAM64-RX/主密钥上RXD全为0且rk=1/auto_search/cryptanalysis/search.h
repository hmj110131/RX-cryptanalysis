/*
 * search.h
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 */

#ifndef CRYPTANALYSIS_SEARCH_H_
#define CRYPTANALYSIS_SEARCH_H_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "typedef.h"
#include "ciphers.h"
#include "globalvar.h"
/* The search modes be used. */
#define mode0 0;
#define mode1 1;
#define mode2 2;










/**
 * @deal with the preparatory information!
 * @param[in]
 * @param[out]  #ciper #mode
 */
int pre_processing_info(void);
/**
 * @def The varibits
 *
 */
u64 varibits(u64 input);
/**
 * @def The doublebits
 *
 */
u64 doublebits(u64 input);
/**
 * @def Differential probability Hamming weight compute.
 *
 */
u16 HM_weight(u64 input);
/**
 * @def  SIMON Differential probability Hamming weight compute.
 *
 */
u16 SIMON_DP_weight_compute(u64 alpha);
u16 SIMON_SCorr_weight_compute(u64 beta);
u16 SIMON_SCorr_weight_check(u64 beta,u64 *U_base);
//

/**
 * @def Differential probability valid check.
 *
 */
u16 SIMON_DP_weight_check(u64 alpha, u64 gama);
/**
 * @ HW space with the hamming weight increase.
 * @
 */
void HW_space_construct(void);
/**
 * @ Max probability of the corresponding input difference.
 * @ From Liuzhengbing's paper Theory 3.
 */
u16 Pmax_compute(u16 input_alpha_hw);
u16 sCorrelation_max_compute(u16 input_beta_hw);
/**
 * @ Construct the DDTA of simon to find the property gama.
 * @
 */
void DDTA_construct(void);
void ARX_carry_DDTm_construct(void);
void ARX_borrow_DDTm_construct(void);

u64 Speck_XDP_validcheck(u64 alpha, u64 beta, u64 gamma);
u16 Speck_XDP_compute(u64 alpha, u64 beta, u64 gamma);
// 只返回最大概率重量该gamma对应的差分概率为最大值
u16 Speck_XDP_Max(u64 alpha, u64 beta, u64 *gamma_tmp);
u64 Speck_XDP_Max_Gamma(u64 alpha, u64 beta);
u16 Speck_block_wt_compute(u64 alpha, u64 beta, u64 gamma, u64 asign_bits);
u64 Combination_cal(u64 M0, u64 N0, char *T );
u16 Cal_nonezero_subscript(u64 input, u16 *sub );
/**
 * @ The entry of the search program,and set the expected Bn_w
 * @
 */
u64 search_entry(u16 search_round);


/**
 * @Search for the min weight of DP of the first Round of SIMON completely.
 * @
 */
u16 round_1(u16 search_round);
/**
 * @Search for the Differential of SIMON.
 * @
 */
u16 round_1_diff(u16 search_round);
/**
 * @Search for the FIRST Round of SIMON completely.
 * @Round 1.
 */
u16 round_1_j(u16 search_round);
u16 round_1_Lineartrail(u16 search_round);
u16 round_1_LinearHull(u16 search_round);
/**
 * @Search for the min weight of DP of the 2 Round of SIMON completely.
 * @
 */
u16 round_2(u16 search_round);
/**
 * @Search for the FIRST Round of SIMON completely.
 * @Round 2.
 */
u16 round_2_diff(u16 search_round);
/**
 * @Jump to 2rd Round of SIMON completely.
 * @
 */
u16 round_2_j(u16 search_round);
u16 round_2_Lineartrail(u16 search_round);
u16 round_2_LinearHull(u16 search_round);

/**
 * @Search for the min weight of DP of the r>2 Round of SIMON completely.
 * @
 */
u16 round_r(u16 search_round, u16 cur_round);
/**
 * @Search for the FIRST Round of SIMON completely.
 * @Round r.
 */
u16 round_r_diff(u16 search_round, u16 cur_round);
u16 round_r_Lineartrail(u16 search_round, u16 cur_round);
u16 round_r_LinearHull(u16 search_round, u16 cur_round);
/**
 * @Search for the min weight of DP of the last Nth Round of SIMON completely.
 * @
 */
u16 round_N(u16 search_round);
/**
 * @Search for the FIRST Round of SIMON completely.
 * @Round LAST.
 */
u16 round_N_diff(u16 search_round);
u16 round_N_Lineartrail(u16 search_round);
u16 round_N_LinearHull(u16 search_round);




///////////////////////////////////////////
//2019年01月23日15:25:54
//添加对SPECK的线性分析的代码
void ARX_cLAT_construct(void);





/**
 * @Search for the min weight of DP of SIMON of N rounds.
 * @ MOde:0
 */
u16 findMinDPWeightDifferentialTrail (u16 search_round);
/**
 * @Search for the Differential of SIMON.
 * @ MOde:1
 */
u16 find_Differntial_characteristic (u16 search_round);
/**
 * @Search for the min weight of LC of SIMON of N rounds.
 * @ MOde:2
 */
u16 findMinDPWeightLinearTrail(u16 search_round);
/**
 * @Search for the Differential of SIMON.
 * @ MOde:3
 */
u16 find_linearHull_characteristic (u16 search_round);














#endif /* CRYPTANALYSIS_SEARCH_H_ */
