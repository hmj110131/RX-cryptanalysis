/*
 * ciphers.h
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 */

#ifndef CIPHERS_CIPHERS_H_
#define CIPHERS_CIPHERS_H_

#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "globalvar.h"

#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include "math.h"





/* The current the specific cipher for searching.*/
 u8 cipher_name;
#define simon 0;
#define simeck 1;
#define speck 2;
#define hight 3;
#define sparx 4;
#define chaskey 5;
#define cham 6;
#define alzette 7;
#define siphash 8;

//......#define ciphername 3-8;
#define gimli 9;
//////////////////



 /*the Rotation parameter of simon-LIKE cipher, abc */
 extern u64 rol_a;	/* a */
 extern u64 rol_b;   /* b */
 extern u64 rol_c;	/* c */
 extern u64 rol_delt;   /* a-b */

 /*the Rotation parameter of simon cipher, abc=182,used in 16 bit cpu */
 #define ROTATE_LEFT(x,s,n)  (((x) << s) | ((x) >> (n-s)))
 #define ROTATE_RIGHT(x, s, n) (((x) >> (s)) | ((x) << (n - s)))

 /*the Rotation parameter of simon cipher, abc=182,used in 32 bit cpu */
 #define rol1_32(x)	( ((x)<<1) | ((x)>>31) )  /* a */
 #define rol8_32(x)	( ((x)<<8) | ((x)>>24) )  /* b */
 #define rol2_32(x)	( ((x)<<2) | ((x)>>30) )  /* c */
 /*the Rotation parameter of simeck cipher, abc=051,used in 32 bit cpu */
 #define rol0_32(x)	( x )                     /* a */
 #define rol5_32(x)	( ((x)<<5) | ((x)>>27) )  /* b */
 #define rol1_32(x)	( ((x)<<1) | ((x)>>31) )  /* c */
 /*the Left Rotation parameter of simon cipher, abc=182,used in 64 bit cpu */
 #define rol1_64(x)	( ((x)<<1) | ((x)>>63) )  /* a */
 #define rol8_64(x)	( ((x)<<8) | ((x)>>56) )  /* b */
 #define rol2_64(x)	( ((x)<<2) | ((x)>>62) )  /* c */
 /*the Left Rotation parameter of simeck cipher, abc=051,used in 64 bit cpu */
 #define rol0_64(x)	( x )                     /* a */
 #define rol5_64(x)	( ((x)<<5) | ((x)>>59) )  /* b */
 #define rol1_64(x)	( ((x)<<1) | ((x)>>63) )  /* c */


/* The current variant of the specific ciphers.*/
extern u8 cipher_variant;
#define simon32_64 0
#define simon48_72 1
#define simon48_96 2
#define simon64_96 3
#define simon64_128 4
#define simon96_96 5
#define simon96_144 6
#define simon128_128 7
#define simon128_192 8
#define simon128_256 9
#define simeck32_64 10
#define simeck48_96 11
#define simeck64_128 12
#define speck32_64 13
#define speck48_72 14
#define speck48_96 15
#define speck64_96 16
#define speck64_128 17
#define speck96_96 18
#define speck96_144 19
#define speck128_128 20
#define speck128_192 21
#define speck128_256 22
/*the block size of specific block cipher variants. */
extern u8 cipher_wordsize;
#define simon32_wordsize 16
#define simon48_wordsize 24
#define simon64_wordsize 32
#define simon96_wordsize 48
#define simon128_wordsize 64
#define simeck32_wordsize 16
#define simeck48_wordsize 24
#define simeck64_wordsize 32
#define speck32_wordsize 16
#define speck48_wordsize 24
#define speck64_wordsize 32simon32/64_rounds
#define speck96_wordsize 48
#define speck128_wordsize 64
/*the total number of Round of specific block cipher variants. */
extern u8 cipher_rounds;
#define simon32_64_rounds 32
#define simon48_72_rounds 36
#define simon48_96_rounds 36
#define simon64_96_rounds 42
#define simon64_128_rounds 44
#define simon96_96_rounds 52
#define simon96_144_rounds 54
#define simon128_128_rounds 68
#define simon128_192_rounds 69
#define simon128_256_rounds 72
#define simeck32_64_rounds 32
#define simeck48_96_rounds 36
#define simeck64_128_rounds 44
#define speck32_64_rounds 22
#define speck48_72_rounds 22
#define speck48_96_rounds 23
#define speck64_96_rounds 26
#define speck64_128_rounds 27
#define speck96_96_rounds 28
#define speck96_144_rounds 29
#define speck128_128_rounds 32
#define speck128_192_rounds 33
#define speck128_256_rounds 34
/*the Rotation parameter of simon cipher, abc=182. */
#define simon_rol_a 8	/* a */
#define simon_rol_b 1   /* b */
#define simon_rol_c 2	/* c */
#define simon_rol_delt 7   /* a-b */
/*the Rotation parameter of simeck cipher, abc=051. */
#define simeck_rol_a 5	/* a */
#define simeck_rol_b 0   /* b */
#define simeck_rol_c 1	/* c */
#define simeck_rol_delt 5   /* a-b */




/**
 * @def The left Rotation in 16 bits
 *
 */
u16 left_rotation_16(u16 input, u8 rotation_abc);
/**
 * @def The left Rotation in 32 bits
 *
 */
u32 left_rotation_32(u32 input, u8 rotation_abc);
/**
 * @def The left Rotation in 64 bits
 *
 */
u64 left_rotation_64(u64 input, u8 rotation_abc);
/**
 * @def The right Rotation in 16 bits.
 *
 */
u16 right_rotation_16(u16 input, u8 rotation_abc);
/**
 * @def The rightt Rotation in 32 bits.
 *
 */
u32 right_rotation_32(u32 input, u8 rotation_abc);
/**
 * @def The right Rotation in 64 bits.
 *
 */
u64 right_rotation_64(u64 input, u8 rotation_abc);


u64 LEFT_rotation_64(u64 x, u64 s,u64 n);
u64 RIGHT_rotation_64(u64 x, u64 s,u64 n);

u64 ROL_left64(u64 input, u64 rol_bit); //rol_bit小于32位
u64 ROL_right64(u64 input, u64 rol_bit); //rol_bit小于32位

#endif /* CIPHERS_CIPHERS_H_ */
