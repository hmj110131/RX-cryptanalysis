/*
 * RX_speck_data.h
 *
 *  Created on: 2020年9月9日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_SPECK_DATA_H_
#define RX_CRYPTANALYSIS_RX_SPECK_DATA_H_

#include <stdio.h>
#include "typedef.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include <nmmintrin.h>
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
u16 RX_sepck_data_trail_search_entry (u16 search_round);
u16 RX_sepck_data_round_1(u16 search_round);
u16 RX_sepck_data_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_data_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi);
u16 RX_sepck_data_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 RX_sepck_data_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);




void print_RX_speck_data_resoult(u16 search_round);











#endif /* RX_CRYPTANALYSIS_RX_SPECK_DATA_H_ */
