/*
 * verify.h
 *
 *  Created on: 2021年6月22日
 *      Author: hmj110131
 */

#ifndef RX_CRYPTANALYSIS_VERIFY_H_
#define RX_CRYPTANALYSIS_VERIFY_H_


#include <stdlib.h>
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include <time.h>
#include "search.h"
#include "globalvar.h"
#include "RX_crypt.h"
#include "printinfo.h"
#include "XOR_offset_table.h"
#include <nmmintrin.h>


u64 Alpha;
u64 Beta;
u64 Gamma;

int pre_processing_TEST_info(void);
long double comput_Pr(u64 a, u64 b, u64 c);
u64 Number_of_xy(u64 alpha, u64 beta, u64 gamma);
u64 Verify_Pr(u64 alpha, u64 beta, u64 gamma);
u64 g_validcheck(u64 alpha, u64 beta, u64 gamma);
u64 w_compute(u64 alpha, u64 beta, u64 gamma);





#endif /* RX_CRYPTANALYSIS_VERIFY_H_ */
