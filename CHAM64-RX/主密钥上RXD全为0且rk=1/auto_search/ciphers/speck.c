/*
 * speck.c
 *
 *  Created on: 2017年9月5日
 *      Author: hmj110131
 */
#include "speck.h"
#include "typedef.h"
#include <stdio.h>
#include <assert.h>
#include "math.h"
#include <time.h>
#include "search.h"
#include "ciphers.h"
#include "search.h"
#include "printinfo.h"


/*  // SPECK96的wt=0的3个差分状态
 0x000000000080    0x000000000000     -0
 0x800000000000    0x800000000000
------------------------------------------
 0x000000000080    0x800000000000     -0
 0x000000000000    0x000000000004
------------------------------------------
 0x000000000000    0x800000000000     -0
 0x800000000000    0x800000000004
*/
/*  //SPECK128的wt=0的差分状态
----------------------------------------------------------
03     0x0000000000000000    0x8000000000000000     -0
04     0x8000000000000000    0x8000000000000004
----------------------------------------------------------
03     0x0000000000000080    0x8000000000000000     -0
04     0x0000000000000000    0x0000000000000004
----------------------------------------------------------
03     0x0000000000000080    0x0000000000000000     -0
04     0x8000000000000000    0x8000000000000000
*/




//新的trail-core扩展法用到的变量
u16 P_up[24] = {0};
u16 P_down[24] = {0};
u16 round_up=0,round_down=0;
u64 up_x[24] = {0};
u64 up_y[24] = {0};
u16 up_wt_min = 100;

u64 down_x[24] = {0};
u64 down_y[24] = {0};
u16 down_wt_min = 100;
u16 Bn_up=0, Bn_down=0;
u16 Bn_TC = 200;


//原始的tc扩展的用到的变量

u64 up_trail_l[12] = {0};
u64 up_trail_r[12] = {0};
u16 up_P_wt[12] = {0};
u64 opt_up_trail_l[12] = {0};
u64 opt_up_trail_r[12] = {0};
u16 opt_up_P_wt[12] = {0};



u64 down_trail_l[12] = {0};
u64 down_trail_r[12] = {0};
u16 down_P_wt[12] = {0};
u64 opt_down_trail_l[12] = {0};
u64 opt_down_trail_r[12] = {0};
u16 opt_down_P_wt[12] = {0};

/*
//For SPECK//
u16 core_round = 0;   // 需要预先设置
u64 in_core_l = 0;    // 需要预先设置
u64 in_core_r = 0;    // 需要预先设置
u64 core_x[16] = {0};   //trail-core 最长15轮
u64 core_y[16] = {0};
u64 out_core_l = 0;   // 需要预先设置
u64 out_core_r = 0;   // 需要预先设置
u16 P_core[16] = {0};
*/

/**
 * Get the size of the key in bits depending on the word size
 *
 * \param word_size word size in bits
 */
u32 speck_get_keysize(u32 word_size)
{
	u32 m = 0;
  switch(word_size) {
  case 16:
	 m = 64;
	 break;
  case 24:
	 m = 96;
	 break;
  case 32:
	 m = 96;
	 //	 m = 128;
	 break;
  case 48:
	 m = 144;
	 break;
  case 64:
	 m = 256;
	 break;
  default:
	 break;
  }
  return m;
}
/**
 * Compute the number of key words depending on the word size
 *
 * \param word_size word size
 * \param key_size key size in bits
 */
u32 speck_compute_nkeywords(u32 word_size, u32 key_size)
{
  if(word_size == 16) {
	 assert((key_size == 64));
  }
  if(word_size == 24) {
	 assert((key_size == 72) || (key_size == 96));
  }
  if(word_size == 32) {
	 assert((key_size == 96) || (key_size == 128));
  }
  if(word_size == 48) {
	 assert((key_size == 96) || (key_size == 144));
  }
  if(word_size == 64) {
	 assert((key_size == 128) || (key_size == 192) || (key_size == 256));
  }
  u32 m = key_size / word_size;
  return m;
}

/**
 * Get the rotation constants.
 */
void speck_get_rot_const(u32 word_size, u32* alpha, u32* beta)
{
  if(word_size == 16) {
	 *alpha = SPECK_RIGHT_ROT_CONST_16BITS;
	 *beta = SPECK_LEFT_ROT_CONST_16BITS;
  } else {
	 *alpha = SPECK_RIGHT_ROT_CONST;
	 *beta = SPECK_LEFT_ROT_CONST;
  }
}

/**
 * Compute the number of rounds for Speck and the index of the z-sequence
 * \param word_size word size
 * \param nkey_words number of key words
 * \return number of rounds
 */
u32 speck_compute_nrounds(u32 word_size, u32 nkey_words)
{
	u32 nrounds = 0;

  switch(word_size) {
  case 16:
	 nrounds = 22;
	 break;
  case 24:
	 if(nkey_words == 3) {
		nrounds = 22;
	 }
	 if(nkey_words == 4) {
		nrounds = 23;
	 }
	 break;
  case 32:
	 if(nkey_words == 3) {
		nrounds = 26;
	 }
	 if(nkey_words == 4) {
		nrounds = 27;
	 }
	 break;
  case 48:
	 if(nkey_words == 2) {
		nrounds = 28;
	 }
	 if(nkey_words == 3) {
		nrounds = 29;
	 }
	 break;
  case 64:
	 if(nkey_words == 2) {
		nrounds = 32;
	 }
	 if(nkey_words == 3) {
		nrounds = 33;
	 }
	 if(nkey_words == 4) {
		nrounds = 34;
	 }
	 break;
  default:
	 break;
  }
  return nrounds;
}
/**
 * Speck key expansion procedure.
 * \param key original key (with enough space for the expanded key)
 * \param nrounds number of rounds
 * \param nkey_words number of key words
 * \param alpha right rotation constant
 * \param beta left rotation constant
 */
void speck_key_expansion(WORD_T key[SPECK_MAX_NROUNDS], u32 nrounds, u32 nkey_words,
		u32 alpha, u32 beta)
{
	u32 T = nrounds;
	u32 m = nkey_words;
  WORD_T L[SPECK_MAX_NROUNDS] = {0};

  for(u32 i = 1; i < m; i++) { // l[m-2], ..., l[0]
	 L[i - 1] = key[i];
  }

  for(u32 i = 0; i < (T - 1); i++) {
	 L[i + m - 1] = ADD(key[i], RROT(L[i], alpha)) ^ i;
	 key[i + 1] = LROT(key[i], beta) ^ L[i + m - 1];
  }
}

/**
 * Speck encryption procedure.
 * \param key expanded key
 * \param nrounds number of rounds
 * \param alpha right rotation constant
 * \param beta left rotation constant
 * \param x_in first plaintext word
 * \param y_in second plaintext word
 */
void speck_encrypt(WORD_T key[SPECK_MAX_NROUNDS], u32 nrounds,
		u32 alpha, u32 beta, WORD_T* x_in, WORD_T* y_in)
{
  WORD_T T = nrounds;
  WORD_T x = *x_in;
  WORD_T y = *y_in;

  for(WORD_T i = 0; i < T; i++) {
#if 0									  // DEBUG
	 printf("[%s:%d] %2d: %8X %8X\n", __FILE__, __LINE__, i, x, y);
#endif
	 x = ADD(RROT(x, alpha), y) ^ key[i];
	 y = LROT(y, beta) ^ x;
  }
#if 0									  // DEBUG
  printf("[%s:%d] %2d: %8X %8X\n", __FILE__, __LINE__, T, x, y);
#endif
  *x_in = x;
  *y_in = y;
}

void speck_decrypt(WORD_T key[SPECK_MAX_NROUNDS], u32 nrounds,
		u32 alpha, u32 beta,WORD_T* x_in, WORD_T* y_in)
{
  WORD_T T = nrounds;
  WORD_T x = *x_in;
  WORD_T y = *y_in;

  for(WORD_T i = 0; i < T; i++) {
#if 0									  // DEBUG
	 printf("[%s:%d] %2d: %8X %8X\n", __FILE__, __LINE__, i, x, y);
#endif
	 y = RROT((y ^ x), beta);
	 x = LROT(SUB((x ^ key[T - i - 1]), y), alpha); // apply keys in reverse order
  }
#if 0									  // DEBUG
  printf("[%s:%d] %2d: %8X %8X\n", __FILE__, __LINE__, T, x, y);
#endif
  *x_in = x;
  *y_in = y;
}



u16 speck_trailCore_Extendsion(u16 search_round)
{
	u16 i = 0;
	u16 extend_rounds = 0;
	u16 up_r_flag = 0;
	u16 up_round = 0;
	u16 down_round = 0;
	u16 diff_wt = 100;
	u16 up_wt = 0, down_wt = 0;
	/*
	u16 core_round = 9;   // 需要预先设置
	u64 core_x[10] = {0}; // 需要预先设置
	u64 core_y[10] = {0}; // 需要预先设置
	u16 P_core[10] = {0}; // 需要预先设置
	u64 in_core_l = 0;    // 需要预先设置
	u64 in_core_r = 0;    // 需要预先设置
	u64 out_core_l = 0;   // 需要预先设置
	u64 out_core_r = 0;   // 需要预先设置
	*/
	u64 in_trail_l[20] = {0};
	u64 in_trail_r[20] = {0};
	u64 out_trail_l[20] = {0};
	u64 out_trail_r[20] = {0};
	u64 in_trail_l_tmp[20] = {0};
	u64 in_trail_r_tmp[20] = {0};
	u64 out_trail_l_tmp[20] = {0};
	u64 out_trail_r_tmp[20] = {0};
	u64 diff_trail_x[30] = {0};
	u64 diff_trail_y[30] = {0};
	u16 P_up[20] = {0};
	u16 P_down[20] = {0};
	u16 P_up_tmp[20] = {0};
	u16 P_down_tmp[20] = {0};
	u16 P_trail[30] = {0};



	//For SPECK96 of 7 rounds.
	u16 core_round = 7;   // 需要预先设置
	u64 in_core_l = 0x000092400040;    // 需要预先设置
	u64 in_core_r = 0x400000104200;    // 需要预先设置
	u64 core_x[16] = {0x000092400040,
					  0x000000820200,
					  0x000000009000,
					  0x000000000080,
					  0x800000000000,
					  0x808000000000,
					  0x800080000004,
					  0x808080800020,
					  0,0,0,0,0,0,0,0
					  };
	u64 core_y[16] = {0x400000104200,
			          0x000000001202,
					  0x000000000010,
					  0x000000000000,
					  0x800000000000,
					  0x808000000004,
					  0x840080000020,
					  0xa08480800124,
					  0,0,0,0,0,0,0,0
					  };
	u64 out_core_l = 0x808080800020;   // 需要预先设置
	u64 out_core_r = 0xa08480800124;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,0,0,0,0,0,0,0,0,0};

/*
	//For SPECK96 of 8 rounds.
	u16 core_round = 8;   // 需要预先设置
	u64 in_core_l = 0x000092400040;    // 需要预先设置
	u64 in_core_r = 0x400000104200;    // 需要预先设置
	u64 core_x[16] = {0x000092400040,
			          0x000000820200,
					  0x000000009000,
					  0x000000000080,
					  0x800000000000,
					  0x808000000000,
					  0x800080000004,
					  0x808080800020,
					  0x800400008124,
					  0,0,0,0,0,0,0
					  };
	u64 core_y[16] = {0x400000104200,
			          0x000000001202,
					  0x000000000010,
					  0x000000000000,
					  0x800000000000,
					  0x808000000004,
					  0x840080000020,
					  0xa08480800124,
					  0x842004008801,
					  0,0,0,0,0,0,0
					  };
	u64 out_core_l = 0x800400008124;   // 需要预先设置
	u64 out_core_r = 0x842004008801;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,9,0,0,0,0,0,0,0,0};


//For SPECK96 of 9 rounds.
	u16 core_round = 9;   // 需要预先设置
	u64 in_core_l = 0x000092400040;    // 需要预先设置
	u64 in_core_r = 0x400000104200;    // 需要预先设置
	u64 core_x[16] = {0x000092400040,
			          0x000000820200,
					  0x000000009000,
					  0x000000000080,
					  0x800000000000,
					  0x808000000000,
					  0x800080000004,
					  0x808080800020,
					  0x800400008124,
					  0x20a000008880,
					  0,0,0,0,0,0
					  };
	u64 core_y[16] = {0x400000104200,
			          0x000000001202,
					  0x000000000010,
					  0x000000000000,
					  0x800000000000,
					  0x808000000004,
					  0x840080000020,
					  0xa08480800124,
					  0x842004008801,
					  0x01a02004c88c,
					  0,0,0,0,0,0
					  };
	u64 out_core_l = 0x20a000008880;   // 需要预先设置
	u64 out_core_r = 0x01a02004c88c;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,9,9,0,0,0,0,0,0,0};

	//For SPECK96 of 10 rounds.
		u16 core_round = 10;   // 需要预先设置
		u64 in_core_l = 0x00800a080800;    // 需要预先设置
		u64 in_core_r = 0x0800124a0848;    // 需要预先设置
		u64 core_x[16] = {0x00800a080800,
						  0x000092400040,
						  0x000000820200,
						  0x000000009000,
						  0x000000000080,
						  0x800000000000,
						  0x808000000000,
						  0x800080000004,
						  0x808080800020,
						  0x800400008124,
						  0xa0a000008880,
						  0,0,0,0,0
						  };
		u64 core_y[16] = {0x0800124a0848,
						  0x400000104200,
				          0x000000001202,
						  0x000000000010,
						  0x000000000000,
						  0x800000000000,
						  0x808000000004,
						  0x840080000020,
						  0xa08480800124,
						  0x842004008801,
						  0x81a02004c88c, //0x01a02004c88c
						  0,0,0,0,0
						  };
		u64 out_core_l = 0x20a000008880;   // 需要预先设置
		u64 out_core_r = 0x01a02004c88c;   // 需要预先设置
		u16 P_core[16] = {10,6,4,2,0,1,3,5,9,9,0,0,0,0,0,0};

	//For SPECK96 of 12 rounds.
		u16 core_round = 12;   // 需要预先设置
		u64 in_core_l = 0x900f00480900;    // 需要预先设置
		u64 in_core_r = 0x011003084009;    // 需要预先设置
		u64 core_x[16] = {0x900f00480900,
						  0x00800a080800,
						  0x000092400040,
						  0x000000820200,
						  0x000000009000,
						  0x000000000080,
						  0x800000000000,
						  0x808000000000,
						  0x800080000004,
						  0x808080800020,
						  0x800400008124,
						  0xa0a000008880,
						  0x80808004c81c,
						  0,0,0
						  };
		u64 core_y[16] = {0x011003084009,
						  0x0800124a0848,
						  0x400000104200,
				          0x000000001202,
						  0x000000000010,
						  0x000000000000,
						  0x800000000000,
						  0x808000000004,
						  0x840080000020,
						  0xa08480800124,
						  0x842004008801,
						  0x81a02004c88c,
						  0x8d8180228c7c,
						  0,0,0
						  };
		u64 out_core_l = 0x80808004c81c;   // 需要预先设置
		u64 out_core_r = 0x8d8180228c7c;   // 需要预先设置
		u16 P_core[16] = {11,10,6,4,2,0,1,3,5,9,9,12,0,0,0,0};
*/


/*  //有问题的
	//For SPECK96 f 8-->12 rounds to extend.,又8轮tc得到的12轮作为新的tc.
	u16 core_round = 12;   // 需要预先设置
	u64 in_core_l = 0x000092400040;    // 需要预先设置
	u64 in_core_r = 0x400000104200;    // 需要预先设置
	u64 core_x[16] = {0x000092400040,
			          0x000000820200,
					  0x000000009000,
					  0x000000000080,
					  0x800000000000,
					  0x808000000000,
					  0x800080000004,
					  0x808080800020,
					  0x800400008124,
					  0x842004038880,
					  0x000000004404,
					  0x000000000020,
					  0x000000000300,
					  0,0,0
					  };
	u64 core_y[16] = {0x400000104200,
			          0x000000001202,
					  0x000000000010,
					  0x000000000000,
					  0x800000000000,
					  0x808000000004,
					  0x840080000020,
					  0xa08480800124,
					  0x842004008801,
					  0x00000000c88c,
					  0x000000000064,
					  0x000000000300,
					  0x000000001b00,
					  0,0,0
					  };
	u64 out_core_l = 0x000000000300;   // 需要预先设置
	u64 out_core_r = 0x000000001b00;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,9,9,9,3,2,0,0,0,0};
*/

/*
	//For SPECK128 of 8 rounds.
	u16 core_round = 8;   // 需要预先设置
	u64 in_core_l = 0x0000000092400040;    // 需要预先设置
	u64 in_core_r = 0x4000000000104200;    // 需要预先设置
	u64 core_x[16] = {0x0000000092400040,
			          0x0000000000820200,
					  0x0000000000009000,
					  0x0000000000000080,
					  0x8000000000000000,
					  0x8080000000000000,
					  0x8000800000000004,
					  0x8080808000000020,
					  0x8004000080000124,
					  0,0,0,0,0,0,0
					  };
	u64 core_y[16] = {0x4000000000104200,
			          0x0000000000001202,
					  0x0000000000000010,
					  0x0000000000000000,
					  0x8000000000000000,
					  0x8080000000000004,
					  0x8400800000000020,
					  0xa084808000000124,
					  0x8420040080000801,
					  0,0,0,0,0,0,0
					  };
	u64 out_core_l = 0x8004000080000124;   // 需要预先设置
	u64 out_core_r = 0x8420040080000801;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,9,0,0,0,0,0,0,0,0};
*/
/*
	//For SPECK128 of 9 rounds.
	u16 core_round = 9;   // 需要预先设置
	u64 in_core_l = 0x0000000092400040;    // 需要预先设置
	u64 in_core_r = 0x4000000000104200;    // 需要预先设置
	u64 core_x[16] = {0x0000000092400040,
			          0x0000000000820200,
					  0x0000000000009000,
					  0x0000000000000080,
					  0x8000000000000000,
					  0x8080000000000000,
					  0x8000800000000004,
					  0x8080808000000020,
					  0x8004000080000124,
					  0x20a0000080800800,
					  0,0,0,0,0,0
					  };
	u64 core_y[16] = {0x4000000000104200,
			          0x0000000000001202,
					  0x0000000000000010,
					  0x0000000000000000,
					  0x8000000000000000,
					  0x8080000000000004,
					  0x8400800000000020,
					  0xa084808000000124,
					  0x8420040080000801,
					  0x01a020048080480c,
					  0,0,0,0,0,0
					  };
	u64 out_core_l = 0x20a0000080800800;   // 需要预先设置
	u64 out_core_r = 0x01a020048080480c;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,9,9,0,0,0,0,0,0,0};

	//For SPECK128 of 10 rounds.
	u16 core_round = 10;   // 需要预先设置
	u64 in_core_l = 0x000000800a080800;    // 需要预先设置
	u64 in_core_r = 0x08000000124a0848;    // 需要预先设置
	u64 core_x[16] = {0x000000800a080800,
			          0x0000000092400040,
			          0x0000000000820200,
					  0x0000000000009000,
					  0x0000000000000080,
					  0x8000000000000000,
					  0x8080000000000000,
					  0x8000800000000004,
					  0x8080808000000020,
					  0x8004000080000124,
					  0xa0a0000080800800,
					  0,0,0,0,0
					  };
	u64 core_y[16] = {0x08000000124a0848,
					  0x4000000000104200,
			          0x0000000000001202,
					  0x0000000000000010,
					  0x0000000000000000,
					  0x8000000000000000,
					  0x8080000000000004,
					  0x8400800000000020,
					  0xa084808000000124,
					  0x8420040080000801,
					  0x81a020048080480c,
					  0,0,0,0,0
					  };
	u64 out_core_l = 0xa0a0000080800800;   // 需要预先设置
	u64 out_core_r = 0x81a020048080480c;   // 需要预先设置
	u16 P_core[16] = {10,6,4,2,0,1,3,5,9,9,0,0,0,0,0,0};

	//For SPECK128 of 12 rounds.
	u16 core_round = 12;   // 需要预先设置
	u64 in_core_l = 0x0000900700480900;    // 需要预先设置
	u64 in_core_r = 0x0100001003084009;    // 需要预先设置
	u64 core_x[16] = {0x0000900700480900,
					  0x000000800a080800,
			          0x0000000092400040,
			          0x0000000000820200,
					  0x0000000000009000,
					  0x0000000000000080,
					  0x8000000000000000,
					  0x8080000000000000,
					  0x8000800000000004,
					  0x8080808000000020,
					  0x8004000080000124,
					  0xa0a0000080800800,
					  0x8000800480004804,
					  0,0,0
					  };
	u64 core_y[16] = {0x0100001003084009,
					  0x08000000124a0848,
					  0x4000000000104200,
			          0x0000000000001202,
					  0x0000000000000010,
					  0x0000000000000000,
					  0x8000000000000000,
					  0x8080000000000004,
					  0x8400800000000020,
					  0xa084808000000124,
					  0x8420040080000801,
					  0x81a020048080480c,
					  0x8d01802084020860,
					  0,0,0
					  };
	u64 out_core_l = 0x8000800480004804;   // 需要预先设置
	u64 out_core_r = 0x8d01802084020860;   // 需要预先设置
	u16 P_core[16] = {11,10,6,4,2,0,1,3,5,9,9,13,0,0,0,0};
*/

/*
////有问题的
	//For SPECK128 of 8-->12 rounds to extend.,又8轮tc得到的12轮作为新的tc.
	u16 core_round = 12;   // 需要预先设置
	u64 in_core_l = 0x0000000092400040;    // 需要预先设置
	u64 in_core_r = 0x4000000000104200;    // 需要预先设置
	u64 core_x[16] = {0x0000000092400040,
			          0x0000000000820200,
					  0x0000000000009000,
					  0x0000000000000080,
					  0x8000000000000000,
					  0x8080000000000000,
					  0x8000800000000004,
					  0x8080808000000020,
					  0x8004000080000124,
					  0x8420040080000800,
					  0x0000000000004804,
					  0x0000000000000824,
					  0x0000000000004904,
					  0,0,0
					  };
	u64 core_y[16] = {0x4000000000104200,
			          0x0000000000001202,
					  0x0000000000000010,
					  0x0000000000000000,
					  0x8000000000000000,
					  0x8080000000000004,
					  0x8400800000000020,
					  0xa084808000000124,
					  0x8420040080000801,
					  0x000000000000480c,
					  0x0000000000000864,
					  0x0000000000004b04,
					  0x0000000000001124,
					  0,0,0
					  };
	u64 out_core_l = 0x0000000000004904;   // 需要预先设置
	u64 out_core_r = 0x0000000000001124;   // 需要预先设置
	u16 P_core[16] = {6,4,2,0,1,3,5,9,6,4,5,6,0,0,0,0};

*/

	/*
 /////////////////////////////////////////////
	//For SPECK96 of wt=0 的第1个差分状态
		u16 core_round = 1;   // 需要预先设置
		u64 in_core_l = 0x000000000080;    // 需要预先设置
		u64 in_core_r = 0x000000000000;    // 需要预先设置
		u64 core_x[16] = {0x000000000080,
						  0x800000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 core_y[16] = {0x000000000000,
						  0x800000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 out_core_l = 0x800000000000;   // 需要预先设置
		u64 out_core_r = 0x800000000000;   // 需要预先设置
		u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		//For SPECK96 of wt=0 的第 2 个差分状态
		u16 core_round = 1;   // 需要预先设置
		u64 in_core_l = 0x000000000080;    // 需要预先设置
		u64 in_core_r = 0x800000000000;    // 需要预先设置
		u64 core_x[16] = {0x000000000080,
						  0x000000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 core_y[16] = {0x800000000000,
						  0x000000000004,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 out_core_l = 0x000000000000;   // 需要预先设置
		u64 out_core_r = 0x000000000004;   // 需要预先设置
		u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		//For SPECK96 of wt=0 的第 3 个差分状态
		u16 core_round = 1;   // 需要预先设置
		u64 in_core_l = 0x000000000000;    // 需要预先设置
		u64 in_core_r = 0x800000000000;    // 需要预先设置
		u64 core_x[16] = {0x000000000000,
						  0x800000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 core_y[16] = {0x800000000000,
						  0x800000000004,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 out_core_l = 0x800000000000;   // 需要预先设置
		u64 out_core_r = 0x800000000004;   // 需要预先设置
		u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
*/

/*
	//For SPECK128 of wt=0的第 1 个差分状态
	u16 core_round = 1;   // 需要预先设置
	u64 in_core_l = 0x0000000000000000;    // 需要预先设置
	u64 in_core_r = 0x8000000000000000;    // 需要预先设置
	u64 core_x[16] = {0x0000000000000000,
					  0x8000000000000000,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	u64 core_y[16] = {0x8000000000000000,
					  0x8000000000000004,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	u64 out_core_l = 0x8000000000000000;   // 需要预先设置
	u64 out_core_r = 0x8000000000000004;   // 需要预先设置
	u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	//For SPECK128 of wt=0的第 2 个差分状态
	u16 core_round = 1;   // 需要预先设置
	u64 in_core_l = 0x0000000000000080;    // 需要预先设置
	u64 in_core_r = 0x8000000000000000;    // 需要预先设置
	u64 core_x[16] = {0x0000000000000080,
					  0x0000000000000000,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	u64 core_y[16] = {0x8000000000000000,
			     	 0x0000000000000004,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	u64 out_core_l = 0x0000000000000000;   // 需要预先设置
	u64 out_core_r = 0x0000000000000004;   // 需要预先设置
	u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	//For SPECK128 of wt=0的第 3 个差分状态
	u16 core_round = 1;   // 需要预先设置
	u64 in_core_l = 0x0000000000000080;    // 需要预先设置
	u64 in_core_r = 0x0000000000000000;    // 需要预先设置
	u64 core_x[16] = {0x0000000000000080,
					  0x8000000000000000,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	u64 core_y[16] = {0x0000000000000000,
					  0x8000000000000000,
					  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	u64 out_core_l = 0x8000000000000000;   // 需要预先设置
	u64 out_core_r = 0x8000000000000000;   // 需要预先设置
	u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
*/



	extend_rounds = search_round - core_round;
	printf("SPECK_%d trail-core round: %d \n",sc_blocksize, core_round);

	printf("Constructing ARX Tables... \n");
	ARX_carry_DDTm_construct();
	ARX_borrow_DDTm_construct();
	time_ARX_DDT = clock();
	run_time =  (time_ARX_DDT - time_start) / CLOCKS_PER_SEC;
	printf("Time of Construct SPECK tables: %.2f seconds.  \n", run_time);


	for(up_round=0; up_round <= extend_rounds;up_round++)  // <= 不搜索向上extend_rounds,向下0轮的情况.
	{
		down_round = extend_rounds - up_round;

		up_wt_min = 100;   // up的轮数概率重量和的最小值IIe{456]2015
		down_wt_min = 100;  // down的轮数概率重量和的最小值
		upper_core_trail(up_round,1,in_core_l,in_core_r,0);
		down_core_trail(down_round,1,out_core_l,out_core_r,0);

		printf("up_round: %d  up_wt_min:%d   down_wt_min:%d \n",up_round,up_wt_min,down_wt_min);

		if((up_wt_min + down_wt_min) < diff_wt)
		{
			up_wt = up_wt_min;
			down_wt = down_wt_min;
			diff_wt = up_wt + down_wt;
			up_r_flag = up_round;

			for(i=0; i < up_r_flag; i++) //暂时取出前r0轮和后r1轮的路径
			{
				in_trail_l[i] = opt_up_trail_l[i+1];
				in_trail_r[i] = opt_up_trail_r[i+1];
				P_up[i] = opt_up_P_wt[i+1];
			}
			for(i=0; i < extend_rounds - up_r_flag; i++)
			{
				out_trail_l[i] = opt_down_trail_l[i+1];
				out_trail_r[i] = opt_down_trail_r[i+1];
				P_down[i] = opt_down_P_wt[i+1];
			}
		}
	}

	for(i=0; i < up_r_flag; i++)
	{
		diff_trail_x[i] = in_trail_l[up_r_flag - 1 -i];
		diff_trail_y[i] = in_trail_r[up_r_flag - 1 -i];
		P_trail[i] = P_up[up_r_flag - 1 -i];
	}
	for(i=up_r_flag; i <= core_round + up_r_flag; i++)
	{
		diff_trail_x[i] = core_x[i - up_r_flag];
		diff_trail_y[i] = core_y[i - up_r_flag];
		P_trail[i] = P_core[i - up_r_flag];
		diff_wt += P_core[i - up_r_flag];
	}
	P_trail[core_round + up_r_flag] = P_down[0];
	for(i=core_round + up_r_flag + 1; i <= extend_rounds + core_round; i++)
	{
		diff_trail_x[i] = out_trail_l[i - core_round - up_r_flag - 1];
		diff_trail_y[i] = out_trail_r[i - core_round - up_r_flag - 1];
		P_trail[i] = P_down[i - core_round - up_r_flag];
	}

	switch(blocksize_len)
	{
	case 48:
		printf("round-------left-------------right--------weight \n");
		for(i=0;i < extend_rounds + core_round;i++)
		{
			printf("%02d     0x%012llx    0x%012llx     -%d \n",i+1,diff_trail_x[i],diff_trail_y[i],P_trail[i]); //,gama_cnt[r_num]);
		}
		printf("%02d     0x%012llx    0x%012llx     NULL \n",(extend_rounds + core_round +1),diff_trail_x[extend_rounds + core_round],diff_trail_y[extend_rounds + core_round]);

		break;
	case 64:
		printf("round---------left---------------right-------------weight \n");
		for(i=0;i < extend_rounds + core_round;i++)
		{
			printf("%02d     0x%016llx    0x%016llx     -%d \n",i+1,diff_trail_x[i],diff_trail_y[i],P_trail[i]); //,gama_cnt[r_num]);
		}
		printf("%02d     0x%016llx    0x%016llx     NULL \n",(extend_rounds + core_round +1),diff_trail_x[extend_rounds + core_round],diff_trail_y[extend_rounds + core_round]);
		break;
	default:  //no
		printf("error ! \n");
		break;
	}
	printf("%d round Total trail weight:  -%d \n",extend_rounds + core_round,diff_wt );
	printf("------------------------------------------------------------ \n");  //----gama_cnt \n");

	time_finish = clock();
	run_time =  (double)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
	printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours. \n", run_time,run_time/60.0,run_time/3600.0 );
	printf("Auto-search END! \n");
	printf("|************************************************************************|\n");
	return 1;
}

u16 upper_core_trail(u16 upper_rounds,u16 cur_rounds,u64 left_x, u64 right_y,u16 p_sumofr)
{
	u16 P_up_wt = 100;
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 in_gamma = 0;
	u64 in_beta  = 0;
	u64 up_gamma = 0;
	u64 up_beta  = 0;
	u64 i = 0, j = 0;
	u64 CB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 alpha_temp = 0;
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;  //w0=0,

	if(upper_rounds == 0)
	{
		up_wt_min = 0;
		return 0;
	}
	in_gamma = left_x & Bit_Align;;
	in_beta =  (ROTATE_RIGHT((right_y ^ left_x), 3, blocksize_len)) & Bit_Align;;


	switch(sc_blocksize)
	{
	case 96:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			gamma_bloc[j] = (in_gamma >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (in_beta >> (8*j))  & 0xFF; //8 bit
			CB_block[j] = ((gamma_bloc[j] << 8) | beta_bloc[j]);
		}
		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((beta_bloc[j-1] >> 7) << 1) + (gamma_bloc[j-1] >> 7);
		}

		for(i0 = bDDT_wt_min[carry[0]][CB_block[0]]; i0 <= bDDT_wt_max[carry[0]][CB_block[0]]; i0++)
		{
			if(i0 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < bDDT_n[carry[0]][CB_block[0]][i0]; j0++)
		{
			alpha_bloc[0] = bDDT_v[carry[0]][CB_block[0]][i0][j0];
			carry[1] = ((alpha_bloc[0] >> 7) << 2) + carry_tmp[1]; //
		for(i1 = bDDT_wt_min[carry[1]][CB_block[1]]; i1 <= bDDT_wt_max[carry[1]][CB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < bDDT_n[carry[1]][CB_block[1]][i1]; j1++)
		{
			alpha_bloc[1] = bDDT_v[carry[1]][CB_block[1]][i1][j1];
			carry[2] = ((alpha_bloc[1] >> 7) << 2) + carry_tmp[2]; // gamma MSB
		for(i2 = bDDT_wt_min[carry[2]][CB_block[2]]; i2 <= bDDT_wt_max[carry[2]][CB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + p_sumofr > up_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < bDDT_n[carry[2]][CB_block[2]][i2]; j2++)
		{
			alpha_bloc[2] = bDDT_v[carry[2]][CB_block[2]][i2][j2];
			carry[3] = ((alpha_bloc[2] >> 7) << 2) + carry_tmp[3]; // gamma MSB
		for(i3 = bDDT_wt_min[carry[3]][CB_block[3]]; i3 <= bDDT_wt_max[carry[3]][CB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < bDDT_n[carry[3]][CB_block[3]][i3]; j3++)
		{
			alpha_bloc[3] = bDDT_v[carry[3]][CB_block[3]][i3][j3];
			carry[4] = ((alpha_bloc[3] >> 7) << 2) + carry_tmp[4]; // gamma MSB
		for(i4 = bDDT_wt_min[carry[4]][CB_block[4]]; i4 <= bDDT_wt_max[carry[4]][CB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + p_sumofr > up_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < bDDT_n[carry[4]][CB_block[4]][i4]; j4++)
		{
			alpha_bloc[4] = bDDT_v[carry[4]][CB_block[4]][i4][j4];
			carry[5] = ((alpha_bloc[4] >> 7) << 2) + carry_tmp[5]; // gamma MSB
		for(i5 = MSB_bDDT_wt_min[carry[5]][CB_block[5]]; i5 <= MSB_bDDT_wt_max[carry[5]][CB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 + p_sumofr > up_wt_min){break;}  //>=则表示不记录同为最小概率重量的路径,不要覆盖,>则表示要覆盖
		for(j5=0; j5 < MSB_bDDT_n[carry[5]][CB_block[5]][i5]; j5++)
		{
			alpha_bloc[5] = MSB_bDDT_v[carry[5]][CB_block[5]][i5][j5];

			alpha_temp = (alpha_bloc[5] << 40) | (alpha_bloc[4] << 32)
					| (alpha_bloc[3] << 24) | (alpha_bloc[2] << 16)
					| (alpha_bloc[1] <<8) | alpha_bloc[0];

			up_gamma = ROTATE_LEFT(alpha_temp,8,blocksize_len) & Bit_Align;
			up_P_wt[cur_rounds] = w5;
			up_trail_l[cur_rounds] = up_gamma;
			up_trail_r[cur_rounds] = in_beta;

			if(cur_rounds < upper_rounds)
			{
				P_up_wt = upper_core_trail(upper_rounds,cur_rounds+1,
						up_gamma,in_beta,up_P_wt[cur_rounds] + p_sumofr);
			}
			else  // 已经是扩展到r0的最上面一轮了
			{
				for(i=1;i<= upper_rounds;i++)  //只记录一次
				{
					opt_up_trail_l[i] = up_trail_l[i];
					opt_up_trail_r[i] = up_trail_r[i];
					opt_up_P_wt[i] = up_P_wt[i];
				}
				up_wt_min = up_P_wt[cur_rounds] + p_sumofr; //覆盖原来的最小值
				P_up_wt = up_wt_min;
				break; //找到1个路径就跳出for,继续判断下一次的for循环,因为概率已经确定
				//return P_up_wt;
			}
		}}}}}}}}}}}}
		break;

	case 128:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			gamma_bloc[j] = (in_gamma >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (in_beta >> (8*j))  & 0xFF; //8 bit
			CB_block[j] = ((gamma_bloc[j] << 8) | beta_bloc[j]);
		}
		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((beta_bloc[j-1] >> 7) << 1) + (gamma_bloc[j-1] >> 7);
		}


		for(i0 = bDDT_wt_min[carry[0]][CB_block[0]]; i0 <= bDDT_wt_max[carry[0]][CB_block[0]]; i0++)
		{
			if(i0 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < bDDT_n[carry[0]][CB_block[0]][i0]; j0++)
		{
			alpha_bloc[0] = bDDT_v[carry[0]][CB_block[0]][i0][j0];
			carry[1] = ((alpha_bloc[0] >> 7) << 2) + carry_tmp[1]; //
		for(i1 = bDDT_wt_min[carry[1]][CB_block[1]]; i1 <= bDDT_wt_max[carry[1]][CB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < bDDT_n[carry[1]][CB_block[1]][i1]; j1++)
		{
			alpha_bloc[1] = bDDT_v[carry[1]][CB_block[1]][i1][j1];
			carry[2] = ((alpha_bloc[1] >> 7) << 2) + carry_tmp[2]; // gamma MSB
		for(i2 = bDDT_wt_min[carry[2]][CB_block[2]]; i2 <= bDDT_wt_max[carry[2]][CB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + p_sumofr > up_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < bDDT_n[carry[2]][CB_block[2]][i2]; j2++)
		{
			alpha_bloc[2] = bDDT_v[carry[2]][CB_block[2]][i2][j2];
			carry[3] = ((alpha_bloc[2] >> 7) << 2) + carry_tmp[3]; // gamma MSB
		for(i3 = bDDT_wt_min[carry[3]][CB_block[3]]; i3 <= bDDT_wt_max[carry[3]][CB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < bDDT_n[carry[3]][CB_block[3]][i3]; j3++)
		{
			alpha_bloc[3] = bDDT_v[carry[3]][CB_block[3]][i3][j3];
			carry[4] = ((alpha_bloc[3] >> 7) << 2) + carry_tmp[4]; // gamma MSB
		for(i4 = bDDT_wt_min[carry[4]][CB_block[4]]; i4 <= bDDT_wt_max[carry[4]][CB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + p_sumofr > up_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < bDDT_n[carry[4]][CB_block[4]][i4]; j4++)
		{
			alpha_bloc[4] = bDDT_v[carry[4]][CB_block[4]][i4][j4];
			carry[5] = ((alpha_bloc[4] >> 7) << 2) + carry_tmp[5]; // gamma MSB
		for(i5 = bDDT_wt_min[carry[5]][CB_block[5]]; i5 <= bDDT_wt_max[carry[5]][CB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 + p_sumofr > up_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < bDDT_n[carry[5]][CB_block[5]][i5]; j5++)
		{
			alpha_bloc[5] = bDDT_v[carry[5]][CB_block[5]][i5][j5];
			carry[6] = ((alpha_bloc[5] >> 7) << 2) + carry_tmp[6]; // gamma MSB
		for(i6 = bDDT_wt_min[carry[6]][CB_block[6]]; i6 <= bDDT_wt_max[carry[6]][CB_block[6]]; i6++)
		{
			w6 = w5 + i6;
			if(w6 + p_sumofr > up_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < bDDT_n[carry[6]][CB_block[6]][i6]; j6++)
		{
			alpha_bloc[6] = bDDT_v[carry[6]][CB_block[6]][i6][j6];
			carry[7] = ((alpha_bloc[6] >> 7) << 2) + carry_tmp[7]; // gamma MSB
		for(i7 = MSB_bDDT_wt_min[carry[7]][CB_block[7]]; i7 <= MSB_bDDT_wt_max[carry[7]][CB_block[7]]; i7++)
		{
			w7 = w6 + i7;
			if(w7 + p_sumofr >= up_wt_min ){break;}  // 一次还是所有
		for(j7=0; j7 < MSB_bDDT_n[carry[7]][CB_block[7]][i7]; j7++)
		{
			alpha_bloc[7] = MSB_bDDT_v[carry[7]][CB_block[7]][i7][j7];
			alpha_temp = (alpha_bloc[7] << 56) | (alpha_bloc[6] << 48)
					| (alpha_bloc[5] << 40) | (alpha_bloc[4] << 32)
					| (alpha_bloc[3] << 24) | (alpha_bloc[2] << 16)
					| (alpha_bloc[1] <<8) | alpha_bloc[0];


			up_gamma = ROTATE_LEFT(alpha_temp,8,blocksize_len) & Bit_Align;
			up_P_wt[cur_rounds] = w7;
			up_trail_l[cur_rounds] = up_gamma;
			up_trail_r[cur_rounds] = in_beta;

			if(cur_rounds < upper_rounds)
			{
				P_up_wt = upper_core_trail(upper_rounds,cur_rounds+1,
						up_gamma,in_beta,up_P_wt[cur_rounds] + p_sumofr);
			}
			else  // 已经是扩展到r0的最上面一轮了
			{
				for(i=1;i<= upper_rounds;i++)  //只记录一次
				{
					opt_up_trail_l[i] = up_trail_l[i];
					opt_up_trail_r[i] = up_trail_r[i];
					opt_up_P_wt[i] = up_P_wt[i];
				}
				up_wt_min = up_P_wt[cur_rounds] + p_sumofr;
				P_up_wt = up_wt_min;
				break;
			}
		}}}}}}}}}}}}}}}}

		break;
	default:
		break;
	}

	return P_up_wt;
}


u16 down_core_trail(u16 down_rounds, u16 cur_rounds,u64 left_x, u64 right_y, u16 p_sumofr)
{
	u16 P_down_wt = 100;
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 in_alpha = 0;
	u64 in_gamma = 0;
	u64 in_beta  = 0;
	u64 down_gamma = 0;
	u64 down_beta  = 0;
	u64 i = 0, j = 0;
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 gamma_temp = 0;
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;  //w0=0,


	if(down_rounds == 0)
	{
		down_wt_min = 0;
		return 0;
	}

	in_alpha = (ROTATE_RIGHT(left_x,8, blocksize_len)) & Bit_Align;
	in_beta = right_y & Bit_Align;;

	//printf("in_alpha:%x   in_alpha:%x \n",in_alpha,in_beta);

	switch(sc_blocksize)
	{
	case 96:
		for(j=0;j<6;j++)  //nBytes
		{
			alpha_bloc[j] = (in_alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (in_beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
		}
		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); //
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + p_sumofr > down_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
			carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]]; i4 <= cDDT_wt_max[carry[4]][AB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + p_sumofr > down_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
		{
			gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
			carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5 <= MSB_cDDT_wt_max[carry[5]][AB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 + p_sumofr > down_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
		{
			gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];

			gamma_temp = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];

			down_beta = (gamma_temp ^ ROTATE_LEFT(in_beta,3,blocksize_len)) & Bit_Align;
			down_P_wt[cur_rounds] = w5;
			down_trail_l[cur_rounds] = gamma_temp;
			down_trail_r[cur_rounds] = down_beta;

			if(cur_rounds < down_rounds)
			{
				P_down_wt = down_core_trail(down_rounds,cur_rounds+1,
						gamma_temp,down_beta,down_P_wt[cur_rounds] + p_sumofr);
			}
			else  // 已经是扩展到r0的最后一轮了
			{
				for(i=1;i<= down_rounds;i++)  //只记录一次
				{
					opt_down_trail_l[i] = down_trail_l[i];
					opt_down_trail_r[i] = down_trail_r[i];
					opt_down_P_wt[i] = down_P_wt[i];
				}
				down_wt_min = down_P_wt[cur_rounds] + p_sumofr;
				P_down_wt = down_wt_min;
				break; //找到1个路径就跳出for,继续判断下一次的for循环
				//return P_down_wt;
			}
		}}}}}}}}}}}}
		break;

	case 128:
		for(j=0;j<8;j++)  //nBytes
		{
			alpha_bloc[j] = (in_alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (in_beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
		}
		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); //
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + p_sumofr > down_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
			carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]]; i4 <= cDDT_wt_max[carry[4]][AB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + p_sumofr > down_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
		{
			gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
			carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = cDDT_wt_min[carry[5]][AB_block[5]]; i5 <= cDDT_wt_max[carry[5]][AB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 + p_sumofr > down_wt_min){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
		{
			gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
			carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7);  // gamma MSB
		for(i6 = cDDT_wt_min[carry[6]][AB_block[6]]; i6 <= cDDT_wt_max[carry[6]][AB_block[6]]; i6++)
		{
			w6 = w5 + i6;
			if(w6 + p_sumofr > down_wt_min ){break;}  //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
		{
			gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];
			carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7);// gamma MSB
		for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7 <= MSB_cDDT_wt_max[carry[7]][AB_block[7]]; i7++)
		{
			w7 = w6 + i7;
			if(w7 + p_sumofr >= down_wt_min ){break;}  // 一次还是所有
		for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
		{
			gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];

			gamma_temp = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
					| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] << 8 ) | gamma_bloc[0];

			down_beta = (gamma_temp ^ ROTATE_LEFT(in_beta,3,blocksize_len)) & Bit_Align;
			down_P_wt[cur_rounds] = w7;
			down_trail_l[cur_rounds] = gamma_temp;
			down_trail_r[cur_rounds] = down_beta;

			if(cur_rounds < down_rounds)
			{
				P_down_wt = down_core_trail(down_rounds,cur_rounds+1,
						gamma_temp,down_beta,down_P_wt[cur_rounds] + p_sumofr);
			}
			else  // 已经是扩展到r1的最后一轮了
			{
				for(i=1;i<= down_rounds;i++)  //只记录一次
				{
					opt_down_trail_l[i] = down_trail_l[i];
					opt_down_trail_r[i] = down_trail_r[i];
					opt_down_P_wt[i] = down_P_wt[i];
				}
				down_wt_min = down_P_wt[cur_rounds] + p_sumofr;
				P_down_wt = down_wt_min;
				break; //找到1个路径就跳出for,继续判断下一次的for循环
				//return P_down_wt;
			}
		}}}}}}}}}}}}}}}}
		break;
	default:
		break;
	}

	return P_down_wt;
}




u16 speck_TC_extend_entry(u16 search_round)
{

	u16 i = 0;
	u16 tc_wt_thrd = 4;   ///该值的大小直接决定搜索的复杂度
	u16 tc_state_wt = 0;


	speck_best_Bn = fopen ("../tmp/speck_best_Bn_wt.xlsx", "a+"); //  "w+"); //

	bit_align_fun();
	ARX_carry_DDTm_construct();

	time_ARX_DDT = clock();
	run_time =  (time_ARX_DDT - time_start) / CLOCKS_PER_SEC;
	printf("Time of Construct ARX_DDT: %.2f seconds.  \n", run_time);

	Bn_w = 0;
	if(search_round > 2)
	{
		for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
		{
			fscanf(speck_best_Bn,"%d",&n_P_bestofR_w[i]);  //read.
			//printf("n_P_bestofR_w[%d]: %d \n",i,n_P_bestofR_w[i]);  // write.
		}
		Bn_w = n_P_bestofR_w[search_round-1] - 1;
	}

	for(tc_state_wt=0; tc_state_wt < tc_wt_thrd; tc_state_wt++)
	{
		sepck_TC_Comb(search_round,tc_state_wt);

	}

	return 1;
}



u16 sepck_TC_Comb(u16 search_round,u16 tc_wt)
{
	u64 M0 = 0;
	u64 N0 = 15;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;


	if(tc_wt < 1)
	{
		sepck_TC_MSB(search_round,tc_wt, C0);
	}
	else
	{
		for ( i=0; i<=(N0-M0); i++) A0[i] = 0;
		for ( i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for ( i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1;
		T0[1] = 0;
		F0[N0] = N0 - M0 + 1;
		I0 = N0 - M0; L0 = N0;

		sepck_TC_MSB(search_round,tc_wt, C0);

		do{
			if (I0 == 0)
			{
				break;
			}
			else
			{
				if (T0[I0] < 0)
				{
					if ((-T0[I0]) != (I0-1))
					{
						T0[I0-1] = T0[I0];
					}
					T0[I0] = I0-1;
				}
				if ( A0[I0]==0 )
				{
					X0 = I0;
					Y0 = F0[L0];
					if (A0[I0-1] == 1)
					{
						F0[I0] = F0[I0 - 1];
					}
					else
					{
						F0[I0] = I0;
					}
					if (F0[L0] == L0)
					{
						L0 = I0; I0 = T0[I0];
						goto CHANGE1;
					}
					if (L0 == N0)
					{
						T0[F0[N0]] = -I0 - 1;
						T0[I0 + 1] = T0[I0];
						I0 = F0[N0];
						F0[N0] = F0[N0] + 1;
						goto CHANGE1;
					}
					T0				[L0] = -I0-1;
					T0[I0+1] = T0[I0];
					F0[L0] = F0[L0] + 1;
					I0 = L0;
					goto CHANGE1;
				}
				Y0 = I0;
				if (I0 != L0)
				{
					F0[L0] = X0 = F0[L0] - 1; F0[I0 - 1] = F0[I0];
					if (L0 == N0)
					{
						if (I0 == (F0[N0] - 1))
						{
							I0 = T0[I0];
							goto CHANGE1;
						}
						T0[F0[N0]-1] = -I0-1;
						T0[I0+1] = T0[I0];
						I0 = F0[N0] - 1;
						goto CHANGE1;
					}
					T0[L0] = -I0 -1; T0[I0 + 1] = T0[I0]; I0 = L0; goto CHANGE1;
				}
				X0 = N0;
				F0[L0 - 1] = F0[L0];
				F0[N0] = N0;
				L0 = N0;
				if (I0 == N0 - 1)
				{
					I0 = T0[N0 - 1];
					goto CHANGE1;
				}
				T0[N0 - 1] = -I0 - 1;
				T0[I0 + 1] = T0[I0];
				I0 = N0 - 1;
	CHANGE1:
			A0[X0] = 1;
			A0[Y0] = 0;
			H0[X0] = Z0 = H0[Y0];
			C0[Z0] = X0;
			}

			sepck_TC_MSB(search_round,tc_wt, C0);


		} while(1);
		///////////////////////////////////////
	}



	return 1;
}




u16 sepck_TC_MSB(u16 search_round,u16 tc_wt,char *posi)
{
	u16 state = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;


	if(tc_wt == 0 )
	{
		for(j_abc=0; j_abc<3; j_abc++)
		{
			if((set_A_3[j_abc] & 0x4) != 0) // alpha bit j.
			{
				Input_alpha = V_MSB;
			}
			else
			{
				Input_alpha = 0;  // & Bit_Align;
			}

			if((set_A_3[j_abc] & 0x2) != 0) // beta bit j.
			{
				Input_beta = V_MSB;
			}
			else
			{
				Input_beta = 0; // & Bit_Align;
			}

			if((set_A_3[j_abc] & 0x1) != 0) // gamma bit j.
			{
				Input_gamma = V_MSB;
			}
			else
			{
				Input_gamma = 0;  // & Bit_Align;
			}

			speck_TC_extend_func(search_round,tc_wt,Input_alpha, Input_beta,Input_gamma);
		}
	}
	else
	{
		for(j_abc = 0; j_abc < 8;j_abc++)   // MSB
		{
			if((j_abc & 0x4) != 0) // alpha bit j.
			{
				Input_alpha = V_MSB;
			}
			else
			{
				Input_alpha = 0;  // & Bit_Align;
			}

			if((j_abc & 0x2) != 0) // beta bit j.
			{
				Input_beta = V_MSB;
			}
			else
			{
				Input_beta = 0; // & Bit_Align;
			}

			if((j_abc & 0x1) != 0) // gamma bit j.
			{
				Input_gamma = V_MSB;
			}
			else
			{
				Input_gamma = 0;  // & Bit_Align;
			}

			if((j_abc==1) || (j_abc==2) || (j_abc==4) || (j_abc==7) )
			{
				Input_alpha |= 0x7F;
				Input_beta  |= 0x7F;
				Input_gamma |= 0x7F;
			}

			sepck_TC_Middle(search_round,tc_wt,Input_alpha, Input_beta,Input_gamma,posi,tc_wt);

		}
	}
	return state;
}



u16 sepck_TC_Middle(u16 search_round,u16 tc_wt,u64 tc_alpha, u64 tc_beta,u64 tc_gamma,char *posi,u16 cur_posi)
{
	u16 state = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	u64 bit_i = 1;

	indx_tmp = posi[cur_posi] - 1;


	// Middle bit positionsof of alpha/beta/gamma
	if( cur_posi > 1)  // Middle bit positionsof of alpha/beta/gamma
	{
		for(j_abc = 1; j_abc < 7; j_abc++)   // MSB //for speck32/48/64
		//for(i = 0; i < 3; i++)   // MSB //for speck96/128
		{
			//j_abc = set_A_3[i];  //for speck96/128

			if((j_abc & 0x4) != 0) // alpha bit j.
			{
				Input_alpha = tc_alpha | (bit_i << indx_tmp);
			}
			else
			{
				//Input_alpha = alpha; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
				Input_alpha = (tc_alpha & ( ~(bit_i << indx_tmp)) ); // & Bit_Align;
			}

			if((j_abc & 0x2) != 0) // alpha bit j.
			{
				Input_beta = tc_beta | (bit_i << indx_tmp); // & Bit_Align;
			}
			else
			{
				//Input_beta = beta; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
				Input_beta = (tc_beta & ( ~(bit_i << indx_tmp)) ); // & Bit_Align;
			}

			if((j_abc & 0x1) != 0) // alpha bit j.
			{
				Input_gamma = tc_gamma | (bit_i << indx_tmp); // & Bit_Align;
			}
			else
			{
				//Input_gamma = gamma; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
				Input_gamma = (tc_gamma & ( ~(bit_i << indx_tmp)) ); // & Bit_Align;
			}


			//set A_3 of Middle
			if( (j_abc==3) || (j_abc==5) || (j_abc==6) )
			{
				// setA对应到下一个不全同比特位置,全设置为0才满足差分有效性条件.
				for(j_last=posi[cur_posi - 1]; j_last < indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
				{
					Input_alpha &= ~(bit_i << j_last);
					Input_beta  &= ~(bit_i << j_last);
					Input_gamma &= ~(bit_i << j_last);
				}
			}     //对于speck96和speck128,不考虑导致后续比特全wei 1情况
			else //set B_3 of Middle //if((j_abc==1) || (j_abc==2) || (j_abc==4) )
			{
				// setB_3对应到下一个不全同比特位置,全设置为 1 才满足差分有效性条件.
				for(j_last=posi[cur_posi - 1]; j_last < indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
				{
					Input_alpha ^= (bit_i << j_last);
					Input_beta  ^= (bit_i << j_last);
					Input_gamma ^= (bit_i << j_last);
				}
			}

			sepck_TC_Middle(search_round,tc_wt,Input_alpha, Input_beta,Input_gamma,posi,cur_posi-1);

		}
	}
	else 	//call last position.
	{
		sepck_TC_Last(search_round,tc_wt, Input_alpha, Input_beta,Input_gamma,posi );

	}

return state;
}



u16 sepck_TC_Last(u16 search_round,u16 tc_wt,u64 tc_alpha, u64 tc_beta,u64 tc_gamma,char *posi )
{
	u16 state = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	u64 bit_i = 1;


	indx_tmp = posi[1] - 1;

	// Last bit position of alpha/beta/gamma // set_A_3
	for(j_abc = 0; j_abc < 3;j_abc++)// Last bit position.
	{
		if((set_A_3[j_abc] & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = tc_alpha | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			Input_alpha = tc_alpha; // & ( ~(1 << indx_tmp))) & Bit_Align;
			//Input_alpha = (alpha & ( ~(1 << indx_tmp))) & Bit_Align;
		}

		if((set_A_3[j_abc] & 0x2) != 0) // alpha bit j.
		{
			Input_beta = tc_beta | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			Input_beta = tc_beta; // & ( ~(1 << indx_tmp))) & Bit_Align;
			//Input_beta = (beta & ( ~(1 << indx_tmp))) & Bit_Align;
		}

		if((set_A_3[j_abc] & 0x1) != 0) // alpha bit j.
		{
			Input_gamma = tc_gamma | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			Input_gamma = tc_gamma; // & ( ~(1 << indx_tmp)) ) & Bit_Align;
			//Input_gamma = (gamma & ( ~(1 << indx_tmp)) ) & Bit_Align;
		}

		speck_TC_extend_func(search_round,tc_wt,Input_alpha, Input_beta,Input_gamma);

	}
	return state;
}


u16 speck_TC_extend_func(u16 search_round,u16 tc_wt,u64 tc_alpha, u64 tc_beta,u64 tc_gamma)
{
	u16 i = 0;
	u16 round_up=0, round_down=0;


	printf(" >>>>>>>>>>>========Trail-Core Extension========<<<<<<<<<<   \n");

	//遍历向下和向上的轮数之和为待搜索的总轮数
	for(round_up=0;round_up < search_round;round_up++)
	{
		round_down = search_round - round_up - 1;
		if(round_up > 0)
		{
			Bn_up = n_P_bestofR_w[round_up] - 1;
		}
		else
		{
			Bn_up = 0;
		}

		if(round_down > 0)
		{
			Bn_down = n_P_bestofR_w[round_down] - 1;
		}
		else
		{
			Bn_down = 0;
		}

		////////////////搜索向下的轮数的最小概率重量/////////////
		if(round_down > 0)
		{
			down_x[1] = tc_gamma;  // output part of the r-th round..
			down_y[1] = tc_gamma ^ (ROTATE_LEFT(tc_beta,rol_b,blocksize_len) & Bit_Align);

			do
			{
				Bn_down = Bn_down + 1;
				printf("Searching round_down =: %d   Bn_down =: %d \n",round_down, Bn_down);

				speck_TC_extend_down(round_down,1,tc_gamma,tc_beta);

				time_Round = clock();
				run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
				printf("Time: %.2f seconds.  \n", run_time);
			}while(Bn_down !=  down_wt_min);   //end conditions.
		}

		/////////////////搜索向上upper的轮数的最小概率重量/////////////
		if(round_up > 0)
		{
			do
			{
				Bn_up = Bn_up + 1;
				printf("Searching round_up =: %d   Bn_up =: %d \n",round_up, Bn_up);

				speck_TC_extend_upper(round_up,1,tc_alpha,tc_beta);

				time_Round = clock();
				run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
				printf("Time: %.2f seconds.  \n", run_time);
			}while(Bn_up !=  up_wt_min);   //end conditions.
		}


		if((Bn_up + Bn_down) < Bn_TC)  // 更新最优路径
		{
			Bn_TC = Bn_up + Bn_down;

			//记录当前最优路径
			for(i=round_up; i>0; i-- )
			{
				n_x_in[round_up - i + 1] = up_x[round_up - i + 1];
				n_y_in[round_up - i + 1] = up_y[round_up - i + 1];
				n_P_w[round_up - i + 1] = P_up[round_up - i + 1];
			}
			n_x_in[round_up + 1] = down_x[1];
			n_y_in[round_up + 1] = down_y[1];
			n_P_w[round_up - i + 1] = tc_wt;
			for(i=2; i<= round_down+1; i++ )
			{
				n_x_in[round_up + i] = down_x[i];
				n_y_in[round_up + i] = down_y[i];
				n_P_w[round_up + i] = P_down[i];
			}

		}
	}
	return 1;
}



u16 speck_TC_extend_upper(u16 upper_rounds,u16 cur_round,u64 in_alpha, u64 in_beta)
{
	u16 best = 0;
	u64 add_alpha = 0;
	u64 add_beta = 0;
	u64 add_gamma = 0;
	u64 up_alpha = 0;
	u64 up_beta = 0;
	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 i = 0, j = 0;
	u16 w_cmp = 0;
	u16 TC_w_cmp = 0;
	u64 xor_beta_gamma = 0;
	u16 w_xor_min = 0;
	u8 w_xor_block[8] = {0};
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u64 CB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;


	up_beta = (in_beta ^ in_alpha) & Bit_Align;
	add_beta = (ROTATE_RIGHT(up_beta, 3, blocksize_len)) & Bit_Align;// beta
	up_beta = add_beta;
	add_gamma = in_alpha; //gamma
	xor_beta_gamma = (add_gamma ^ add_beta) & ValueMax_Align ;
	w_xor_min = HM_weight(xor_beta_gamma );

	for(i=1; i < cur_round; i++)
	{
		p_sumof_r += P_up[i];
	}
	if(upper_rounds == 1 )
	{
		p_sumof_r = 0;
	}

	//向上的轮路需要满足剪枝条件
	w_cmp = Bn_up - p_sumof_r - n_P_bestofR_w[upper_rounds - cur_round];
	if ( w_xor_min > w_cmp)  ////???????????从Bn由小到大逐渐增加
	{
		return 0;
	}
	//向下加上向上的最优概率若大于之前搜索得到的期望轮数的概率，则退出，已经非最优了。
	TC_w_cmp = Bn_TC - Bn_down - p_sumof_r - n_P_bestofR_w[upper_rounds - cur_round];
	if(w_xor_min > TC_w_cmp )
	{
		return 0;
	}

	switch(sc_blocksize)
	{
	case 96:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			gamma_bloc[j] = (add_gamma >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (add_beta >> (8*j))  & 0xFF; //8 bit
			CB_block[j] = ((gamma_bloc[j] << 8) | beta_bloc[j]);
		}
		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((beta_bloc[j-1] >> 7) << 1) + (gamma_bloc[j-1] >> 7);
		}

		for(j=1;j<nBytes;j++)  //nBytes
		{
			w_xor_block[j] = (u8)HM_weight((xor_beta_gamma >> (j*8)));
		}


		for(i0 = bDDT_wt_min[carry[0]][CB_block[0]]; i0 <= bDDT_wt_max[carry[0]][CB_block[0]]; i0++)
		{
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < bDDT_n[carry[0]][CB_block[0]][i0]; j0++)
		{
			alpha_bloc[0] = bDDT_v[carry[0]][CB_block[0]][i0][j0];
			carry[1] = ((alpha_bloc[0] >> 7) << 2) + carry_tmp[1]; //
		for(i1 = bDDT_wt_min[carry[1]][CB_block[1]]; i1 <= bDDT_wt_max[carry[1]][CB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1+ w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < bDDT_n[carry[1]][CB_block[1]][i1]; j1++)
		{
			alpha_bloc[1] = bDDT_v[carry[1]][CB_block[1]][i1][j1];
			carry[2] = ((alpha_bloc[1] >> 7) << 2) + carry_tmp[2]; // gamma MSB
		for(i2 = bDDT_wt_min[carry[2]][CB_block[2]]; i2 <= bDDT_wt_max[carry[2]][CB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < bDDT_n[carry[2]][CB_block[2]][i2]; j2++)
		{
			alpha_bloc[2] = bDDT_v[carry[2]][CB_block[2]][i2][j2];
			carry[3] = ((alpha_bloc[2] >> 7) << 2) + carry_tmp[3]; // gamma MSB
		for(i3 = bDDT_wt_min[carry[3]][CB_block[3]]; i3 <= bDDT_wt_max[carry[3]][CB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < bDDT_n[carry[3]][CB_block[3]][i3]; j3++)
		{
			alpha_bloc[3] = bDDT_v[carry[3]][CB_block[3]][i3][j3];
			carry[4] = ((alpha_bloc[3] >> 7) << 2) + carry_tmp[4]; // gamma MSB
		for(i4 = bDDT_wt_min[carry[4]][CB_block[4]]; i4 <= bDDT_wt_max[carry[4]][CB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < bDDT_n[carry[4]][CB_block[4]][i4]; j4++)
		{
			alpha_bloc[4] = bDDT_v[carry[4]][CB_block[4]][i4][j4];
			carry[5] = ((alpha_bloc[4] >> 7) << 2) + carry_tmp[5]; // gamma MSB
		for(i5 = MSB_bDDT_wt_min[carry[5]][CB_block[5]]; i5 <= MSB_bDDT_wt_max[carry[5]][CB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < MSB_bDDT_n[carry[5]][CB_block[5]][i5]; j5++)
		{
			alpha_bloc[5] = MSB_bDDT_v[carry[5]][CB_block[5]][i5][j5];

			add_alpha = (alpha_bloc[5] << 40) | (alpha_bloc[4] << 32)
					| (alpha_bloc[3] << 24) | (alpha_bloc[2] << 16)
					| (alpha_bloc[1] <<8) | alpha_bloc[0];

			up_alpha = ROTATE_LEFT(add_alpha,8,blocksize_len) & Bit_Align;
			P_up[cur_round] = w5;

			if(upper_rounds == cur_round )
			{
				up_wt_min = p_sumof_r + w5;
				if(up_wt_min ==  Bn_up )
				{
					best = 1;
				}
			}
			else
			{
				best = speck_TC_extend_upper(upper_rounds,cur_round+1, up_alpha, up_beta );
			}

#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			up_x[1] = up_alpha;  // output part of the r-th round..
			up_y[1] = up_beta;
			return 1;
		}
#endif

		}}}}}}}}}}}}
		break;

	case 128:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			gamma_bloc[j] = (in_alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (in_beta >> (8*j))  & 0xFF; //8 bit
			CB_block[j] = ((gamma_bloc[j] << 8) | beta_bloc[j]);
		}
		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((beta_bloc[j-1] >> 7) << 1) + (gamma_bloc[j-1] >> 7);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			w_xor_block[j] = (u8)HM_weight((xor_beta_gamma >> (j*8)));
		}


		for(i0 = bDDT_wt_min[carry[0]][CB_block[0]]; i0 <= bDDT_wt_max[carry[0]][CB_block[0]]; i0++)
		{
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < bDDT_n[carry[0]][CB_block[0]][i0]; j0++)
		{
			alpha_bloc[0] = bDDT_v[carry[0]][CB_block[0]][i0][j0];
			carry[1] = ((alpha_bloc[0] >> 7) << 2) + carry_tmp[1]; //
		for(i1 = bDDT_wt_min[carry[1]][CB_block[1]]; i1 <= bDDT_wt_max[carry[1]][CB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < bDDT_n[carry[1]][CB_block[1]][i1]; j1++)
		{
			alpha_bloc[1] = bDDT_v[carry[1]][CB_block[1]][i1][j1];
			carry[2] = ((alpha_bloc[1] >> 7) << 2) + carry_tmp[2]; // gamma MSB
		for(i2 = bDDT_wt_min[carry[2]][CB_block[2]]; i2 <= bDDT_wt_max[carry[2]][CB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < bDDT_n[carry[2]][CB_block[2]][i2]; j2++)
		{
			alpha_bloc[2] = bDDT_v[carry[2]][CB_block[2]][i2][j2];
			carry[3] = ((alpha_bloc[2] >> 7) << 2) + carry_tmp[3]; // gamma MSB
		for(i3 = bDDT_wt_min[carry[3]][CB_block[3]]; i3 <= bDDT_wt_max[carry[3]][CB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < bDDT_n[carry[3]][CB_block[3]][i3]; j3++)
		{
			alpha_bloc[3] = bDDT_v[carry[3]][CB_block[3]][i3][j3];
			carry[4] = ((alpha_bloc[3] >> 7) << 2) + carry_tmp[4]; // gamma MSB
		for(i4 = bDDT_wt_min[carry[4]][CB_block[4]]; i4 <= bDDT_wt_max[carry[4]][CB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < bDDT_n[carry[4]][CB_block[4]][i4]; j4++)
		{
			alpha_bloc[4] = bDDT_v[carry[4]][CB_block[4]][i4][j4];
			carry[5] = ((alpha_bloc[4] >> 7) << 2) + carry_tmp[5]; // gamma MSB
		for(i5 = bDDT_wt_min[carry[5]][CB_block[5]]; i5 <= bDDT_wt_max[carry[5]][CB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < bDDT_n[carry[5]][CB_block[5]][i5]; j5++)
		{
			alpha_bloc[5] = bDDT_v[carry[5]][CB_block[5]][i5][j5];
			carry[6] = ((alpha_bloc[5] >> 7) << 2) + carry_tmp[6]; // gamma MSB
		for(i6 = bDDT_wt_min[carry[6]][CB_block[6]]; i6 <= bDDT_wt_max[carry[6]][CB_block[6]]; i6++)
		{
			w6 = w5 + i6;
			if(w6 + w_xor_block[7] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < bDDT_n[carry[6]][CB_block[6]][i6]; j6++)
		{
			alpha_bloc[6] = bDDT_v[carry[6]][CB_block[6]][i6][j6];
			carry[7] = ((alpha_bloc[6] >> 7) << 2) + carry_tmp[7]; // gamma MSB
		for(i7 = MSB_bDDT_wt_min[carry[7]][CB_block[7]]; i7 <= MSB_bDDT_wt_max[carry[7]][CB_block[7]]; i7++)
		{
			w7 = w6 + i7;
			if(w7 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j7=0; j7 < MSB_bDDT_n[carry[7]][CB_block[7]][i7]; j7++)
		{
			alpha_bloc[7] = MSB_bDDT_v[carry[7]][CB_block[7]][i7][j7];
			add_alpha = (alpha_bloc[7] << 56) | (alpha_bloc[6] << 48)
					| (alpha_bloc[5] << 40) | (alpha_bloc[4] << 32)
					| (alpha_bloc[3] << 24) | (alpha_bloc[2] << 16)
					| (alpha_bloc[1] <<8) | alpha_bloc[0];

			up_alpha = ROTATE_LEFT(add_alpha,8,blocksize_len) & Bit_Align;
			P_up[cur_round] = w7;

			if(upper_rounds == cur_round )
			{
				up_wt_min = p_sumof_r + w5;
				if(up_wt_min ==  Bn_up )
				{
					best = 1;
				}
			}
			else
			{
				best = speck_TC_extend_upper(upper_rounds,cur_round+1, up_alpha, up_beta );
			}

#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			up_x[1] = up_alpha;  // output part of the r-th round..
			up_y[1] = up_beta;
			return best;
		}
#endif
		}}}}}}}}}}}}}}}}

		break;
	default:
		break;
	}
	return best;
}


u16 speck_TC_extend_down(u16 down_rounds,u16 cur_round,u64 in_gamma, u64 in_beta)
{
	u16 best = 0;
	u64 add_alpha = 0;
	u64 add_beta = 0;
	u64 add_gamma = 0;
//	u64 up_alpha = 0;
//	u64 up_beta = 0;
	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 xor_alpha_beta = 0;
	u16 w_xor_min = 0;
	u16 w_cmp = 0;
	u16 TC_w_cmp = 0;
	u64 i = 0, j = 0;
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u8 w_xor_block[8] = {0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;


/*
	if(down_rounds == cur_round)
	{
		best = speck_round_N(down_rounds, in_gamma, in_beta );
		return best;
	}
*/
	add_beta = (in_gamma ^ ROTATE_LEFT(in_beta,rol_b,blocksize_len) )& Bit_Align;
	add_alpha = (ROTATE_RIGHT(in_gamma,rol_a, blocksize_len)) & Bit_Align;
	xor_alpha_beta = (add_alpha ^ add_beta) & ValueMax_Align ;
	w_xor_min = HM_weight(xor_alpha_beta );



	for(i=1;i<cur_round;i++)
	{
		p_sumof_r += P_down[i]; // p_sumof_r is the sum of r-1 rounds weight.
	}
	if(down_rounds == 1 )
	{
		p_sumof_r = 0;
	}


	//向下轮传播中满足向下的轮靠近最优期望概率
	w_cmp = Bn_down - p_sumof_r - n_P_bestofR_w[down_rounds - cur_round];
	if ( w_xor_min > w_cmp)  ////???????????从Bn由小到大逐渐增加
	{
		return 0;
	}
	//向下加上向上的最优概率若大于之前搜索得到的期望轮数的概率，则退出，已经非最优了。
	TC_w_cmp = Bn_TC - Bn_up - 1 - p_sumof_r - n_P_bestofR_w[down_rounds - cur_round];
	if(w_xor_min > TC_w_cmp )
	{
		return 0;
	}

	switch(sc_blocksize)
	{
	case 96:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (add_alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (add_beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			w_xor_block[j] = (u8)HM_weight((xor_alpha_beta >> (j*8)));
		}

		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
				//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
				carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
				gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
				//carry[4] = ((alpha_bloc[3] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[3] >> 7) << 1)   //beta MSB
				carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]]; i4 <= cDDT_wt_max[carry[4]][AB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5]> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
			{
				gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
				//carry[5] = ((alpha_bloc[4] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[4] >> 7) << 1)   //beta MSB
				carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5 <= MSB_cDDT_wt_max[carry[5]][AB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.

			for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
			{
			gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];
			add_gamma = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];

			P_down[cur_round + 1] = w5;  //

			if(down_rounds == cur_round)
			{
				down_wt_min = p_sumof_r + w5;
				if(down_wt_min ==  Bn_down )
				{
					best = 1;
				}
			}
			else
			{
				best = speck_TC_extend_down(down_rounds,cur_round +1,add_gamma,add_beta);
			}

#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			down_x[cur_round +1] = add_gamma;  // output part of the r-th round..
			down_y[cur_round +1] = add_gamma ^ (ROTATE_LEFT(add_beta, rol_b,blocksize_len) & Bit_Align);
			return 1;
		}
#endif
		}}}}}}}}}}}}
		break;

	case 128:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (add_alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (add_beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			w_xor_block[j] = (u8)HM_weight((xor_alpha_beta >> (j*8)));
		}

		carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//gamma_block_tmp[0] = gamma_bloc[0];
			//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			//gamma_block_tmp[1] = gamma_bloc[1] << 8;
			//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			//gamma_block_tmp[2] = gamma_bloc[2] << 16;
			//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break
		for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
			//gamma_block_tmp[3] = gamma_bloc[3] << 24;
			//carry[4] = ((alpha_bloc[3] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[3] >> 7) << 1)   //beta MSB
			carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]]; i4 <= cDDT_wt_max[carry[4]][AB_block[4]]; i4++)
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
		{
			gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
			//gamma_block_tmp[4] = gamma_bloc[4] << 32;
			//carry[5] = ((alpha_bloc[4] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[4] >> 7) << 1)   //beta MSB
			carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = cDDT_wt_min[carry[5]][AB_block[5]]; i5 <= cDDT_wt_max[carry[5]][AB_block[5]]; i5++)
		{
			w5 = w4 + i5;
			if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
		{
			gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
			//gamma_block_tmp[5] = gamma_bloc[5] << 40;
			//carry[6] = ((alpha_bloc[5] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[5] >> 7) << 1)   //beta MSB
			carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7); // gamma MSB
		for(i6 = cDDT_wt_min[carry[6]][AB_block[6]]; i6 <= cDDT_wt_max[carry[6]][AB_block[6]]; i6++)
		{
			w6 = w5 + i6;
			if(w6 + w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
		{
			gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];
			//gamma_block_tmp[6] = gamma_bloc[6] << 48;
			//carry[7] = ((alpha_bloc[6] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[6] >> 7) << 1)   //beta MSB
			carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7); // gamma MSB
		for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7 <= MSB_cDDT_wt_max[carry[7]][AB_block[7]]; i7++)
		{
			w7 = w6 + i7;
			if(w7> w_cmp ){break;}  //break 要比 continue 高效,直接结束.

		for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
		{
			gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];
			add_gamma = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
					| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];

			P_down[cur_round + 1] = w7;  //

			if(down_rounds == cur_round)
			{
				down_wt_min = p_sumof_r + w7;
				if(down_wt_min ==  Bn_down )
				{
					best = 1;
				}
			}
			else
			{
				best = speck_TC_extend_down(down_rounds,cur_round +1,add_gamma,add_beta);
			}

#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			down_x[cur_round +1] = add_gamma;  // output part of the r-th round..
			down_y[cur_round +1] = add_gamma ^ (ROTATE_LEFT(add_beta, rol_b,blocksize_len) & Bit_Align);
			return best;
		}
#endif

		}}}}}}}}}}}}}}}}
		break;

	default:
		break;

	}

	return best;
}






///////////////////////////////////////////////////////
u16 TC_ext_entry(u16 search_round)
{
	u16 best = 0;
	u16 i = 0;
	u16 extend_rounds = 0;
	u16 opt_wt_MIN = 100;   	u16 diff_wt = 100;
	u16 opt_up_round = 100;
	u64 in_trail_l[20] = {0};
	u64 in_trail_r[20] = {0};
	u64 out_trail_l[20] = {0};
	u64 out_trail_r[20] = {0};

	u64 in_trail_l_tmp[20] = {0};
	u64 in_trail_r_tmp[20] = {0};
	u64 out_trail_l_tmp[20] = {0};
	u64 out_trail_r_tmp[20] = {0};
	u64 diff_trail_x[30] = {0};
	u64 diff_trail_y[30] = {0};
	u16 P_up[20] = {0};
	u16 P_down[20] = {0};
	u16 P_up_tmp[20] = {0};
	u16 P_down_tmp[20] = {0};
	u16 P_trail[30] = {0};


	 /////////////////////////////////////////////
		//For SPECK96 of wt=0 的第1个差分状态
			u16 core_round = 1;   // 需要预先设置
			u64 in_core_l = 0x000000000080;    // 需要预先设置
			u64 in_core_r = 0x000000000000;    // 需要预先设置
			u64 core_x[16] = {0x000000000080,
							  0x800000000000,
							  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			u64 core_y[16] = {0x000000000000,
							  0x800000000000,
							  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			u64 out_core_l = 0x800000000000;   // 需要预先设置
			u64 out_core_r = 0x800000000000;   // 需要预先设置
			u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	/*
			//For SPECK96 of wt=0 的第 2 个差分状态
			u16 core_round = 1;   // 需要预先设置
			u64 in_core_l = 0x000000000080;    // 需要预先设置
			u64 in_core_r = 0x800000000000;    // 需要预先设置
			u64 core_x[16] = {0x000000000080,
							  0x000000000000,
							  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			u64 core_y[16] = {0x800000000000,
							  0x000000000004,
							  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			u64 out_core_l = 0x000000000000;   // 需要预先设置
			u64 out_core_r = 0x000000000004;   // 需要预先设置
			u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

			//For SPECK96 of wt=0 的第 3 个差分状态
			u16 core_round = 1;   // 需要预先设置
			u64 in_core_l = 0x000000000000;    // 需要预先设置
			u64 in_core_r = 0x800000000000;    // 需要预先设置
			u64 core_x[16] = {0x000000000000,
							  0x800000000000,
							  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			u64 core_y[16] = {0x800000000000,
							  0x800000000004,
							  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			u64 out_core_l = 0x800000000000;   // 需要预先设置
			u64 out_core_r = 0x800000000004;   // 需要预先设置
			u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	*/

	/*
		//For SPECK128 of wt=0的第 1 个差分状态
		u16 core_round = 1;   // 需要预先设置
		u64 in_core_l = 0x0000000000000000;    // 需要预先设置
		u64 in_core_r = 0x8000000000000000;    // 需要预先设置
		u64 core_x[16] = {0x0000000000000000,
						  0x8000000000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 core_y[16] = {0x8000000000000000,
						  0x8000000000000004,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 out_core_l = 0x8000000000000000;   // 需要预先设置
		u64 out_core_r = 0x8000000000000004;   // 需要预先设置
		u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		//For SPECK128 of wt=0的第 2 个差分状态
		u16 core_round = 1;   // 需要预先设置
		u64 in_core_l = 0x0000000000000080;    // 需要预先设置
		u64 in_core_r = 0x8000000000000000;    // 需要预先设置
		u64 core_x[16] = {0x0000000000000080,
						  0x0000000000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 core_y[16] = {0x8000000000000000,
				     	 0x0000000000000004,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 out_core_l = 0x0000000000000000;   // 需要预先设置
		u64 out_core_r = 0x0000000000000004;   // 需要预先设置
		u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		//For SPECK128 of wt=0的第 3 个差分状态
		u16 core_round = 1;   // 需要预先设置
		u64 in_core_l = 0x0000000000000080;    // 需要预先设置
		u64 in_core_r = 0x0000000000000000;    // 需要预先设置
		u64 core_x[16] = {0x0000000000000080,
						  0x8000000000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 core_y[16] = {0x0000000000000000,
						  0x8000000000000000,
						  0,0,0,0,0,0,0,0,0,0,0,0,0,0};
		u64 out_core_l = 0x8000000000000000;   // 需要预先设置
		u64 out_core_r = 0x8000000000000000;   // 需要预先设置
		u16 P_core[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	*/


		extend_rounds = search_round - core_round;
		printf("Constructing ARX Tables... \n");
		ARX_carry_DDTm_construct();
		ARX_borrow_DDTm_construct();
		time_ARX_DDT = clock();
		run_time =  (time_ARX_DDT - time_start) / CLOCKS_PER_SEC;
		printf("Time of Construct SPECK tables: %.2f seconds.  \n", run_time);

		up_wt_min = 100;   // up的轮数概率重量和的最小值IIe{456]2015
		down_wt_min = 100;  // down的轮数概率重量和的最小值

		//向上搜索R轮的最小概率重量的路径
		upper_core_trail(search_round,1,in_core_l,in_core_r,0);  //
		//向下搜索R轮的最小概率重量的路径
		down_core_trail(search_round,1,out_core_l,out_core_r,0);
/*
		//opt_wt_MIN = 100;
		for(i=1;i<= search_round-1; i++ )
		{
			if(opt_up_wt_min[i] + opt_down_wt_min[search_round - i] < opt_wt_MIN )
			{
				opt_wt_MIN = opt_up_wt_min[i] + opt_down_wt_min[search_round - i];
				opt_up_round = i;
			}
		}
*/


		// 比较概率重量最小的路径 //并输出路径
		if(up_wt_min < opt_wt_MIN)  // 0+向上R轮为最优
		{
			for(i=0; i < opt_up_round; i++) //暂时取出前r0轮和后r1轮的路径
			{
				in_trail_l[i] = opt_up_trail_l[i+1];
				in_trail_r[i] = opt_up_trail_r[i+1];
				P_up[i] = opt_up_P_wt[i+1];
			}
		}
		else
		if( down_wt_min < opt_wt_MIN)  // 0+向下R轮为最优
		{
			for(i=0; i < search_round - opt_up_round; i++)
			{
				out_trail_l[i] = opt_down_trail_l[i+1];
				out_trail_r[i] = opt_down_trail_r[i+1];
				P_down[i] = opt_down_P_wt[i+1];
			}
		}
		else   // 还是由前面的（1--> R-1 轮的构造的为最优）
		{
			for(i=0; i < opt_up_round; i++) //暂时取出前r0轮和后r1轮的路径
			{
				in_trail_l[i] = opt_up_trail_l[i+1];
				in_trail_r[i] = opt_up_trail_r[i+1];
				P_up[i] = opt_up_P_wt[i+1];
			}
			for(i=0; i < search_round - opt_up_round; i++)
			{
				out_trail_l[i] = opt_down_trail_l[i+1];
				out_trail_r[i] = opt_down_trail_r[i+1];
				P_down[i] = opt_down_P_wt[i+1];
			}
		}


	//组合生成完整的最优差分路径 opt_wt_MIN
	for(i=0; i < opt_up_round; i++)
	{
		diff_trail_x[i] = in_trail_l[opt_up_round - 1 -i];
		diff_trail_y[i] = in_trail_r[opt_up_round - 1 -i];
		P_trail[i] = P_up[opt_up_round - 1 -i];
	}
	for(i=opt_up_round; i <= core_round + opt_up_round; i++)
	{
		diff_trail_x[i] = core_x[i - opt_up_round];
		diff_trail_y[i] = core_y[i - opt_up_round];
		P_trail[i] = P_core[i - opt_up_round];
		diff_wt += P_core[i - opt_up_round];
	}
	P_trail[core_round + opt_up_round] = P_down[0];
	for(i=core_round + opt_up_round + 1; i <= extend_rounds + core_round; i++)
	{
		diff_trail_x[i] = out_trail_l[i - core_round - opt_up_round - 1];
		diff_trail_y[i] = out_trail_r[i - core_round - opt_up_round - 1];
		P_trail[i] = P_down[i - core_round - opt_up_round];
	}






// 打印最优差分路径
switch(blocksize_len)
{
printf("|****************************Trail-Core Extention******************************|\n");
case 48:
	printf("round-------left-------------right--------weight \n");
	for(i=0;i < extend_rounds + core_round;i++)
	{
		printf("%02d     0x%012llx    0x%012llx     -%d \n",i+1,diff_trail_x[i],diff_trail_y[i],P_trail[i]); //,gama_cnt[r_num]);
	}
	printf("%02d     0x%012llx    0x%012llx     NULL \n",(extend_rounds + core_round +1),diff_trail_x[extend_rounds + core_round],diff_trail_y[extend_rounds + core_round]);

	break;
case 64:
	printf("round---------left---------------right-------------weight \n");
	for(i=0;i < extend_rounds + core_round;i++)
	{
		printf("%02d     0x%016llx    0x%016llx     -%d \n",i+1,diff_trail_x[i],diff_trail_y[i],P_trail[i]); //,gama_cnt[r_num]);
	}
	printf("%02d     0x%016llx    0x%016llx     NULL \n",(extend_rounds + core_round +1),diff_trail_x[extend_rounds + core_round],diff_trail_y[extend_rounds + core_round]);
	break;
default:  //no
	printf("error ! \n");
	break;
}
printf("%d round Total trail weight:  -%d \n",extend_rounds + core_round,diff_wt );
printf("------------------------------------------------------------ \n");  //----gama_cnt \n");

time_finish = clock();
run_time =  (double)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours. \n", run_time,run_time/60.0,run_time/3600.0 );
printf("Auto-search END! \n");
printf("|************************************************************************|\n");



	return best;
}

u16 TC_ext_UP(u16 upper_rounds,u16 cur_rounds,u64 left_x, u64 right_y,u16 p_sumofr)
{
	u16 best = 0;



	return best;
}

u16 TC_ext_DOWN(u16 down_rounds,u16 cur_rounds,u64 left_x, u64 right_y,u16 p_sumofr)
{
	u16 best = 0;




	return best;
}

















