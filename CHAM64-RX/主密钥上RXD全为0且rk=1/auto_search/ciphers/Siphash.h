/*
 * Siphash.h
 *
 *  Created on: 2020年7月31日
 *      Author: hmj110131
 */

#ifndef CIPHERS_SIPHASH_H_
#define CIPHERS_SIPHASH_H_

#include <stdio.h>
#include "typedef.h"
#include <nmmintrin.h>


u64 rol_left_64(u64 x, u64 indx);
u64 siphash_rol_left_13(u64 input);
u64 siphash_rol_left_16(u64 input);
u64 siphash_rol_left_32(u64 input);
u64 siphash_rol_left_17(u64 input);
u64 siphash_rol_left_21(u64 input);






#endif /* CIPHERS_SIPHASH_H_ */
