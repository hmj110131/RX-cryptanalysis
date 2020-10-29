/*
 * Siphash.c
 *
 *  Created on: 2020年7月31日
 *      Author: hmj110131
 */
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include <time.h>
#include "Siphash.h"




u64 rol_left_64(u64 x, u64 indx)
{
	return ((x << indx) | (x >> (64-indx))); // & 0xFFFFFFFFFFFFFFFF;
}
u64 siphash_rol_left_13(u64 input)
{
	return rol_left_64(input,13); // & 0xFFFFFFFF;
}
u64 siphash_rol_left_16(u64 input)
{
	return rol_left_64(input,16); // & 0xFFFFFFFF;
}
u64 siphash_rol_left_32(u64 input)
{
	return rol_left_64(input,32); // & 0xFFFFFFFF;
}
u64 siphash_rol_left_17(u64 input)
{
	return rol_left_64(input,17); // & 0xFFFFFFFF;
}
u64 siphash_rol_left_21(u64 input)
{
	return rol_left_64(input,21); // & 0xFFFFFFFF;
}



















