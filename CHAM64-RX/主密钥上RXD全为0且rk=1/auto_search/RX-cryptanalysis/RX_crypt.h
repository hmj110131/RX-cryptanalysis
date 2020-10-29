/*
 * RX_crypt.h
 *
 *  Created on: 2020年7月23日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_RX_CRYPT_H_
#define RX_CRYPTANALYSIS_RX_CRYPT_H_

#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "ciphers.h"
#include "globalvar.h"
#include <time.h>
#include "XOR_offset_table.h"
#include "RX_Alzette.h"
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
///////////以下为RX-dif使用搜索算法框架的全局概率重量参数//////////////////
extern volatile float RX_Expected_w;
extern volatile float RX_P_w[100];
extern volatile float RX_P_bestofR_w[100];
extern volatile float RX_P_sumofR_w[100];
//RX分析的旋转参数，需要提前输入
extern u16 RX_k;
//基于概率重量由大到小逼近，选择r-1轮的输出差分为限定条件
extern u64 F_Alpha, F_Beta;
//基于概率重量由大到小逼近，统计动态变化更新期望概率重量的次数
extern u64 Num_Bn_Update;





//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int pre_processing_RX_info(void);
u16 findMinDPWeight_RX_Trail(u16 search_round);
u16 findMinDPWeight_Relatekey_RX_Trail(u16 search_round);
u16 findMinDPWeight_Weakkey_Relatekey_RX_Trail(u16 search_round);
u16 findMinDPWeight_Openkey_Relatekey_RX_Trail(u16 search_round);









//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//



#endif /* RX_CRYPTANALYSIS_RX_CRYPT_H_ */
