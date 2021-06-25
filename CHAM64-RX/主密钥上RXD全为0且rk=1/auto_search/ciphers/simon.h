/*
 * simon.h
 *
 *  Created on: 2017年9月5日
 *      Author: hmj110131
 */
#ifndef CIPHERS_SIMON_H_
#define CIPHERS_SIMON_H_
#include "typedef.h"


// Macros
#define NROUNDS_MAX 10 /**< Max. number of rounds */
#ifndef WORD_SIZE
#define WORD_SIZE 3 /**< Word size in  bits. */
#endif

#ifndef NROUNDS
#define NROUNDS 4 /**< Number of rounds in reduced-round versions of target ciphers. */
#endif

#ifndef ALL_WORDS
#define ALL_WORDS (1ULL << WORD_SIZE) /**< Total number of words of size WORD_SIZE. */
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

#define SIMON_ZSEQ_LEN 62
#define SIMON_MAX_NROUNDS 72


#define SIMON_LROT_CONST_S 1
#define SIMON_LROT_CONST_T 8
#define SIMON_LROT_CONST_U 2
#define SIMON_NPAIRS (1ULL << 20)

#ifndef LROT
#define LROT(x,r) (((x << r) | (x >> (WORD_SIZE - r))) & MASK) /**< Rotate \p x by \p r positions to the left; \p x is of size \ref WORD_SIZE */
#endif
#ifndef RROT
#define RROT(x,r) (((x >> r) | (x << (WORD_SIZE - r))) & MASK) /**< Rotate \p x by \p r positions to the right; \p x is of size \ref WORD_SIZE */
#endif



extern u32 g_simon_zseq[5][62];
/**
 * Compute the number of key words depending on the word size
 *
 * \param word_size word size
 * \param key_size key size in bits
 */
u32 simon_compute_nkeywords(u32 word_size, u32 key_size);
/**
 * Simon key expansion procedure.
 * \param key original key (with enough space for the expanded key)
 * \param Z the z-sequence (\eref g_simon_zseq)
 * \param zseq_j index of the z-seqence
 * \param nrounds number of rounds
 * \param nkey_words number of key words
 */
void simon_key_expansion(u32 key[], u32 Z[5][62], u32 zseq_j,u32 nrounds, u32 nkey_words);
/**
 * Compute the number of rounds for Simon and the index of the z-sequence
 * \param word_size word size
 * \param nkey_words number of key words
 * \param zseq_j index of the z-sequence \ref g_simon_zseq
 * \return number of rounds
 */
u32 simon_compute_nrounds(u32 word_size, u32 nkey_words, u32* zseq_j);
/**
 * Simon encryption procedure.
 * \param key expanded key
 * \param nrounds number of rounds
 * \param x_in first plaintext word
 * \param y_in second plaintext word
 */
void simon_encrypt(u32 key[], u32 nrounds,u32* x_in, u32* y_in);






#endif /* CIPHERS_SIMON_H_ */
