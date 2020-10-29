/*
 * speck.h
 *
 *  Created on: 2017年9月5日
 *      Author: hmj110131
 */

#ifndef CIPHERS_SPECK_H_
#define CIPHERS_SPECK_H_

#include "typedef.h"
#include "math.h"

#define SPECK_BYTE 128  //SPECK的分组版本，32/48/64/96/128,比switch块

#define SPECK_MAX_NKEY_WORDS 4
#define SPECK_KEY_LEN_BITS 128
#define SPECK_MAX_NROUNDS 34
#define SPECK_RIGHT_ROT_CONST 8
#define SPECK_LEFT_ROT_CONST 3
#define SPECK_RIGHT_ROT_CONST_16BITS 7
#define SPECK_LEFT_ROT_CONST_16BITS 2


#ifndef WORD_T // abstract word type
#if (WORD_SIZE <= 32)
#define WORD_T u32
#else
#define WORD_T uint64_t
#endif // #if (WORD_SIZE <= 32)
#endif // #ifndef WORD

// Macros
#define NROUNDS_MAX 10 /**< Max. number of rounds */
#ifndef WORD_SIZE
#define WORD_SIZE 3 /**< Word size in  bits. */
#endif

#ifndef MASK
#if(WORD_SIZE <= 32)
#define MASK (0xffffffffUL >> (32 - WORD_SIZE)) /**< A mask for the WORD_SIZE LS bits of a 32-bit word. */
#define MASK_NO_MSB (0xffffffffUL >> (32 - (WORD_SIZE - 1)))
#else // #if(WORD_SIZE > 32)
#define MASK (0xffffffffffffffffULL >> (64 - WORD_SIZE)) /**< A mask for the WORD_SIZE LS bits of a 32-bit word. */
#define MASK_NO_MSB (0xffffffffffffffffULL >> (64 - (WORD_SIZE - 1)))
#endif // #if(WORD_SIZE <= 32)
#endif

#ifndef LROT
#define LROT(x,r) (((x << r) | (x >> (WORD_SIZE - r))) & MASK) /**< Rotate \p x by \p r positions to the left; \p x is of size \ref WORD_SIZE */
#endif
#ifndef RROT
#define RROT(x,r) (((x >> r) | (x << (WORD_SIZE - r))) & MASK) /**< Rotate \p x by \p r positions to the right; \p x is of size \ref WORD_SIZE */
#endif

#ifndef MOD
#define MOD (1ULL << WORD_SIZE) /**< The value 2^{WORD_SIZE}. */
#endif

#ifndef ADD
#define ADD(x,y) ((x + y) & MASK) /**< The ADD operation on words of size \ref WORD_SIZE */
#endif

#ifndef SUB
#if(WORD_SIZE < 64)
#define SUB(x,y) ((WORD_T)(x - y + MOD) & MASK) /**< The modular subtraction (SUB) operation on words of size \ref WORD_SIZE */
#else // #if(WORD_SIZE == 64)
#define SUB(x,y) ((WORD_T)(x - y)) /**< The modular subtraction (SUB) operation on words of size 64-bit */
#endif // #if(WORD_SIZE < 64)
#endif // #ifndef SUB



/**
 * Get the size of the key in bits depending on the word size
 *
 * \param word_size word size in bits
 */
u32 speck_get_keysize(u32 word_size);
/**
 * Compute the number of key words depending on the word size
 *
 * \param word_size word size
 * \param key_size key size in bits
 */
u32 speck_compute_nkeywords(u32 word_size, u32 key_size);
/**
 * Get the rotation constants.
 */
void speck_get_rot_const(u32 word_size, u32* alpha, u32* beta);
/**
 * Compute the number of rounds for Speck and the index of the z-sequence
 * \param word_size word size
 * \param nkey_words number of key words
 * \return number of rounds
 */
u32 speck_compute_nrounds(u32 word_size, u32 nkey_words);


u16 speck_trailCore_Extendsion(u16 search_round );
//u16 speck_trailCore_upper(u16 upper_rounds, u64 left_x, u64 right_y, u16 core_wt);
//u16 speck_trailCore_down(u16 down_rounds, u64 left_x, u64 right_y, u16 core_wt);
//u16 upper_core_trail(u16 upper_rounds,u16 cur_rounds, u16 pre_wt, u64 left_x, u64 right_y,
//		u64 *upper_trail_l, u64 *upper_trail_r, u16 *P_up,u16 p_sumofr);
//u16 down_core_trail(u16 down_rounds, u16 cur_rounds, u16 pre_wt, u64 left_x, u64 right_y,
//		u64 *down_trail_l, u64 *down_trail_r,u16 *P_down,u16 p_sumofr);


u16 upper_core_trail(u16 upper_rounds,u16 cur_rounds,u64 left_x, u64 right_y,u16 p_sumofr);
u16 down_core_trail(u16 down_rounds, u16 cur_rounds,u64 left_x, u64 right_y, u16 p_sumofr);


u16 speck_TC_extend_entry(u16 search_round);
u16 sepck_TC_Comb(u16 search_round,u16 tc_wt);
u16 sepck_TC_MSB(u16 search_round,u16 tc_wt,char *posi );
u16 sepck_TC_Middle(u16 search_round,u16 tc_wt,u64 tc_alpha, u64 tc_beta,u64 tc_gamma,char *posi, u16 cur_posi);
u16 sepck_TC_Last(u16 search_round,u16 tc_wt,u64 tc_alpha, u64 tc_beta,u64 tc_gamma,char *posi );
u16 speck_TC_extend_func(u16 search_round,u16 tc_wt,u64 tc_alpha, u64 tc_beta,u64 tc_gamma);
u16 speck_TC_extend_upper(u16 upper_rounds,u16 cur_round,u64 in_alpha, u64 in_beta);
u16 speck_TC_extend_down(u16 down_rounds,u16 cur_round,u64 in_gamma, u64 in_beta);



//原始的tc扩展的用到的变量
extern u16 P_up[24];
extern u16 P_down[24];
extern u16 round_up,round_down;
extern u64 up_x[24];
extern u64 up_y[24];
extern u16 up_wt_min;
extern u64 down_x[24];
extern u64 down_y[24];
extern u16 down_wt_min;
extern u16 Bn_up, Bn_down;
extern u16 Bn_TC;


//原始的tc扩展的用到的变量
extern u64 up_trail_l[12];
extern u64 up_trail_r[12];
extern u16 up_P_wt[12];
extern u64 opt_up_trail_l[12];
extern u64 opt_up_trail_r[12];
extern u16 opt_up_P_wt[12];



extern u64 down_trail_l[12];
extern u64 down_trail_r[12];
extern u16 down_P_wt[12];
extern u64 opt_down_trail_l[12];
extern u64 opt_down_trail_r[12];
extern u16 opt_down_P_wt[12];

/*
//For SPECK128 of 8 rounds.
extern u16 core_round;  // = 8;   // 需要预先设置
extern u64 in_core_l; // = 0x0000000092400040;    // 需要预先设置
extern u64 in_core_r; // = 0x4000000000104200;    // 需要预先设置
extern u64 core_x[16];
extern u64 core_y[16];
extern u64 out_core_l;  // = 0x8004000080000124;   // 需要预先设置
extern u64 out_core_r; // = 0x8420040080000801;   // 需要预先设置
extern u16 P_core[16];  // = {6,4,2,0,1,3,5,9,0,0};
*/




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


/**
 * @Search for the Differential of speck.
 * @ MOde:9
 */
u16 sepck_differential_trail_search_entry (u16 search_round);
void bit_align_fun(void);
u16 sepck_round_1(u16 search_round  );
u16 sepck_round_1_diff(u16 search_round  );
//u16 xdp_speck(u64 alpha, u64 beta, u64 gama);
u16 sepck_input_MSB(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi);
u16 sepck_input_Middle(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 cur_posi);
u16 sepck_input_Last(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi );
u16 sepck_round_r(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 sepck_round_r_diff(u16 search_round, u16 cur_round, u64 x, u64 y);
/*
u16 speck_gamma_xor_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 curr_wt_xor,u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 wt_or,u16 *or_sub_1,
		u16 wt_cmp);
u16 speck_gamma_and_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 wt_or,u16 *or_sub_1,
		u16 wt_cmp);
u16 speck_gamma_or_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 wt_or,u16 *or_sub_1,
		u16 wt_cmp);
u16 speck_gamma_and_or_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 and_or_sub_1_num,u16 *and_or_sub_1,
		u16 wt_cmp);
u16 speck_gamma(u16 curr_sub,u16 sub_1_num, u16 *sub_1,
		u16 curr_sub_0,u16 sub_0_num, u16 *sub_0,
		u16 search_round, u16 cur_round, u64 Alpha, u64 Beta, u64 Gamma,u16 sum_wt_bf,u16 wt_xor,u16 wt_and);
u16 speck_gamma_combine(u16 curr_sub,u16 sub_1_num, u16 *sub_1,
		u16 curr_sub_0,u16 sub_0_num, u16 *sub_0,
		u16 search_round, u16 cur_round, u64 Alpha, u64 Beta, u64 Gamma,u16 sum_wt_bf,u16 P_r);
*/
u16 speck_round_N(u16 search_round, u64 x, u64 y);
u16 speck_round_N_diff(u16 search_round, u64 x, u64 y);
/*
u16 speck_N_gamma_xor_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 curr_wt_xor,u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 wt_or,u16 *or_sub_1,
		u16 wt_cmp);
u16 speck_N_gamma_and_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 wt_or,u16 *or_sub_1,
		u16 wt_cmp);

u16 speck_N_gamma_or_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 wt_or,u16 *or_sub_1,
		u16 wt_cmp);
u16 speck_N_gamma_and_or_inc(u16 search_round, u16 cur_round,
		u64 Alpha, u64 Beta, u64 Gamma,
		u16 wt_xor, u16 *xor_sub_1,
		u16 wt_and,u16 *and_sub_1,
		u16 and_or_sub_1_num,u16 *and_or_sub_1,
		u16 wt_cmp);
u16 speck_gamma_combine_last(u16 curr_sub,u16 sub_1_num, u16 *sub_1,
		u16 curr_sub_0,u16 sub_0_num, u16 *sub_0,
		u16 search_round, u16 cur_round, u64 Alpha, u64 Beta, u64 Gamma,u16 sum_wt_bf,u16 P_N);
u16 speck_gamma_wt_inc_last(u16 curr_sub_0,u16 sub_0_num, u16 *sub_0,
u16 round_N_sepck(u16 search_round,u16 cur_round, u16 bit_position, u64 input_a, u64 input_b, u64 output_c);
		*/


u16 sepck_linear_trail_search_entry (u16 search_round);

u16 sepck_round_1_linear(u16 search_round);
u16 sepck_ADD_mask_MSB
(u16 search_round,u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 sepck_ADD_mask_Middle(u16 search_round,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_1, u16 Corr_tmp_inc,char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 sepck_ADD_mask_Last(u16 search_round,u16 Corr_1, char *posi );

u16 sepck_round_2_linear(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W);
u16 sepck_ADD_mask_Last_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W,
		u16 Corr_2, char *posi );
u16 sepck_ADD_mask_Middle_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);
u16 sepck_ADD_mask_MSB_2(u16 search_round,
		u64 r1_tmp_U, u64 r1_tmp_V, u64 r1_tmp_W,
		u64 tmp_U, u64 tmp_V, u64 tmp_W,
		u16 Corr_2, u16 Corr_tmp_inc, char *posi,
		u16 current_index, u16 pre_cor_0or1);


u16 sepck_round_r_linear(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 speck_round_N_linear(u16 search_round, u64 x, u64 y);


u16 sepck_linear_Hull_search_entry (u16 search_round);
u16 sepck_round_r_linear_hull(u16 search_round, u16 cur_round, u64 x, u64 y);
u16 speck_round_N_linear_hull(u16 search_round, u64 x, u64 y);




u16 TC_ext_entry(u16 search_round);
u16 TC_ext_UP(u16 upper_rounds,u16 cur_rounds,u64 left_x, u64 right_y,u16 p_sumofr);
u16 TC_ext_DOWN(u16 down_rounds,u16 cur_rounds,u64 left_x, u64 right_y,u16 p_sumofr);







#endif /* CIPHERS_SPECK_H_ */
