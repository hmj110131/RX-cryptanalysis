/*
 * chaskey.h
 *
 *  Created on: 2018年11月29日
 *      Author: hmj110131
 */

#ifndef CIPHERS_CHASKEY_H_
#define CIPHERS_CHASKEY_H_


#include <stdio.h>
#include "typedef.h"


extern u16 RX_wt_block[4][20];



u32 rol_left32(u32 x, u32 indx);
u32 chaskey_rol_left_5(u32 input);
u32 chaskey_rol_left_8(u32 input);
u32 chaskey_rol_left_16(u32 input);
u32 chaskey_rol_left_7(u32 input);
u32 chaskey_rol_left_13(u32 input);




u32 chaskey_round_1_left(u16 search_round);
u32 chaskey_round_1_right(u16 search_round);
u32 chaskey_round_1_Down_2add(u16 search_round);




////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
//extern u16 Chaskey_R1_ADD_wt[4];
extern volatile u64 Chaskey_R1_ADD_i_UVW[4][3];
extern volatile u16 Chaskey_Round_ADD_wt[5][4];
extern volatile u32 XV[4][20];
extern volatile u32 XV_M[4][20];

////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
u32 Chaskey_Linear_trail_search_entry(u16 search_round);
u32 Chaskey_R1_ADD_wt_inc(u16 search_round);
u32 Chaskey_R1_ADD_i_tuples(u16 search_round,u16 ADD_part);
u16 Chaskey_R1_ADD_Last(u16 search_round,
		u16 part, u16 Corr_1, char *posi );
u16 Chaskey_R1_ADD_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 part, u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 Chaskey_R1_ADD_MSB(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 part, u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u32 Chaskey_R1_comb_MASKs(u16 search_round);

u32 Chaskey_Middle_Rounds(u16 search_round, u16 current_round, u32 *mask_Xr);
u32 Chaskey_Middle_Rounds_UperHalf(u16 search_round, u16 current_round, u32 *mask_Xr);
u32 Chaskey_Middle_Rounds_DownHalf(u16 search_round, u16 current_round,
		u16 wt_cmpr, u16 uper_wt, u32 *mask_Xr);


u32 Chaskey_N_Rounds_UperHalf(u16 search_round, u32 *mask_Xr);
u32 Chaskey_N_Rounds_DownHalf(u16 search_round,
		u16 wt_cmpr, u16 uper_wt, u32 *mask_Xr);




u16 print_Chaskey_resoult(u16 search_round);




u32 Chaskey_1_Rounds(u32 *X, u32 *Y);






#endif /* CIPHERS_CHASKEY_H_ */
