/*
 * ciphers.c
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 */
#include <stdio.h>
#include <stdlib.h>
#include "ciphers.h"
#include "typedef.h"
#include "globalvar.h"

u8 cipher_name = simon;
u8 cipher_variant = simon32_64;
u8 cipher_wordsize = simon32_wordsize;
u8 cipher_rounds = simon32_64_rounds;


/*the Rotation parameter of simon-LIKE cipher, abc */
u64 rol_a = 8;	/* a */
u64 rol_b = 1;   /* b */
u64 rol_c = 2;	/* c */
u64 rol_delt = 7;   /* a-b */


/**
 * @def The left Rotation in 16 bits
 *
 */
u16 left_rotation_16(u16 input, u8 rotation_abc)
{
	u16 tmp = 0;

	tmp = ((input) << rotation_abc) |  ((input) >> (16 - rotation_abc));

	return tmp;
}
/**
 * @def The left Rotation in 32 bits
 *
 */
u32 left_rotation_32(u32 input, u8 rotation_abc)
{
	u32 tmp = 0;

	tmp = ((input) << rotation_abc) |  ((input) >> (32 - rotation_abc));

	return tmp;
}
/**
 * @def The left Rotation in 64 bits
 *
 */
u64 left_rotation_64(u64 input, u8 rotation_abc)
{
	u64 tmp = 0;

	tmp = ((input) << rotation_abc) |  ((input) >> (64 - rotation_abc));

	return tmp;
}
/**
 * @def The right Rotation in 16 bits.
 *
 */
u16 right_rotation_16(u16 input, u8 rotation_abc)
{
	u16 tmp = 0;

	tmp = ((input) >> rotation_abc) |  ((input) << (16 - rotation_abc));

	return tmp;
}
/**
 * @def The right Rotation in 32bits
 *
 */
u32 right_rotation_32(u32 input, u8 rotation_abc)
{
	u32 tmp = 0;

	tmp = ((input) >> rotation_abc) |  ((input) << (32 - rotation_abc));

	return tmp;
}
/**
 * @def The right Rotation in 64 bits.
 *
 */
u64 right_rotation_64(u64 input, u8 rotation_abc)
{
	u64 tmp = 0;

	tmp = ((input) >> rotation_abc) |  ((input) << (64 - rotation_abc));

	return tmp;
}

u64 LEFT_rotation_64(u64 x, u64 s,u64 n)
{
	u64 tmp = 0;
	u64 tmp_asign = 0;

	tmp = ROTATE_LEFT(x,s,n);

	switch(n)
	{
	case 16:
		tmp_asign = tmp & 0xFFFF;
		break;
	case 24:
		tmp_asign = tmp & 0x00FFFFFF;
		break;
	case 32:
		tmp_asign = tmp & 0xFFFFFFFF;
		break;
	case 48:
		tmp_asign = tmp & 0xFFFFFFFFFFFF;
		break;
	default:  //128
		tmp_asign = tmp & 0xFFFFFFFFFFFFFFFF;
		break;
	}

	return tmp_asign;
}
u64 RIGHT_rotation_64(u64 x, u64 s,u64 n)
{
	u64 tmp = 0;
	u64 tmp_asign = 0;

	tmp = ROTATE_RIGHT(x,s,n);

	switch(n)
	{
	case 16:
		tmp_asign = tmp & 0x0000FFFF;
		break;
	case 24:
		tmp_asign = tmp & 0x00FFFFFF;
		break;
	case 32:
		tmp_asign = tmp & 0xFFFFFFFF;
		break;
	case 48:
		tmp_asign = tmp & 0xFFFFFFFFFFFF;
		break;
	case 64://128
		tmp_asign = tmp & 0xFFFFFFFFFFFFFFFF;
		break;
	default:
		break;
	}

	return tmp_asign;
}


//在32位系统上执行64位数据的循环左移位操作
u64 ROL_left64(u64 input, u64 rol_bit) //rol_bit小于32位
{
	u64 rol_tmp = 0;
	u32 tmp1 = 0;
	u32 tmp2 = 0;
	u32 *p = (u32*)&input;

	printf("P+1: 0x%016x \n", *(p+1) );
	printf("PP: 0x%016x \n", *p);
	// (((x) << s) | ((x) >> (n-s)))
	tmp1 = *(p+1) << rol_bit;  // 高32位左移
	tmp2 = *p >> (32 - rol_bit);
	*(p+1) = tmp1 | tmp2;  //高32位

	tmp1 = *p << rol_bit;  // 低32位左移
	tmp2 = (*p+1) >> (32 - rol_bit);
	*p = tmp1 | tmp2;  //高32位

	printf("P+1: 0x%016x \n", *(p+1) );
	printf("PP: 0x%016x \n", *p);

	rol_tmp = *(p+1)*pow(2,32) + *p;

	return rol_tmp;
}

//在32位系统上执行64位数据的循环右移位操作
u64 ROL_right64(u64 input, u64 rol_bit) //rol_bit小于32位
{
	u64 rol_tmp = 0;
	u32 tmp1 = 0;
	u32 tmp2 = 0;
	u32 *p = (u32*)&input;

	printf("P+1: 0x%016x \n", *(p+1) );
	printf("PP: 0x%016x \n", *p);
	//  ((x) >> (s)) | ((x) << (n - s))
	tmp1 = *(p+1) >> rol_bit;  // 高32位左移
	tmp2 = *p << (32 - rol_bit);
	*(p+1) = tmp1 | tmp2;  //高32位

	tmp1 = *p >> rol_bit;  // 低32位左移
	tmp2 = *(p+1) << (32 - rol_bit);
	*p = tmp1 | tmp2;  //高32位

	printf("P+1: 0x%016x \n", *(p+1) );
	printf("PP: 0x%016x \n", *p);

	printf("input: 0x%llx  \n", input);

	rol_tmp = (*(p+1)) * pow(2,32) + *p;

	return rol_tmp;
}






