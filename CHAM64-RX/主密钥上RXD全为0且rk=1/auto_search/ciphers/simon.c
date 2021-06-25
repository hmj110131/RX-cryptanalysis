/*
 * simon.c
 *
 *  Created on: 2017年9月5日
 *      Author: hmj110131
 */
#include "simon.h"
#include "typedef.h"
#include <stdio.h>
#include <assert.h>





/**
 * Pre-computed z_j sequences (o <= j < 5) used in the key schedule of Simon.
 */
u32 g_simon_zseq[5][62] =  {
  // z_0
  {1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0,  // 31
	1,1,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,1,1,1,0,0,1,1,0},
  // z_1
  {1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,0,
	1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,1,0,0,0,0,1,0,1,1,0,1,0},
  // z_2
  {1,0,1,0,1,1,1,1,0,1,1,1,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,
	0,1,0,1,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,1,1,0,1,1,0,0,1,1},
  // z_3
  {1,1,0,1,1,0,1,1,1,0,1,0,1,1,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,0,
	0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,1,1,0,0,1,1,0,1,0,0,0,0,1,1,1,1},
  // z_4
  {1,1,0,1,0,0,0,1,1,1,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,
	0,0,1,0,1,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,0,1,1,1,0,1,1,1,1}
};
/**
 * Compute the number of key words depending on the word size
 *
 * \param word_size word size
 * \param key_size key size in bits
 */
u32 simon_compute_nkeywords(u32 word_size, u32 key_size)
{
  if(word_size == 16)
  {
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
 * Simon key expansion procedure.
 * \param key original key (with enough space for the expanded key)
 * \param Z the z-sequence (\eref g_simon_zseq)
 * \param zseq_j index of the z-seqence
 * \param nrounds number of rounds
 * \param nkey_words number of key words
 */
void simon_key_expansion(u32 key[], u32 Z[5][62], u32 zseq_j,u32 nrounds, u32 nkey_words)
{
	u32 T = nrounds;
	u32 m = nkey_words;
	u32 r1 = 3;				  // rot const
	u32 r2 = 1;				  // rot const
	u32 xconst = 3;
	u32 j = zseq_j;

  assert(m <= T);
  assert(key[m] == 0);
  assert(j < 5);

 for(u32 i = m; i < T; i++) {
	 u32 tmp = RROT(key[i - 1], r1);
	 if(m == 4) {
		tmp ^= key[i - 3];		  // !
	 }
	 tmp ^= RROT(tmp, r2);
	 u32 k = (i - m) % SIMON_ZSEQ_LEN;
	 u32 inv_key = (~(key[i - m])) & MASK;
#if 0									  // DEBUG
	 printf("[%s:%d] %8X %8X\n", __FILE__, __LINE__, key[i-m], (~key[i-m]) & MASK);
	 print_binary(key[i-m]);
	 printf("\n");
	 print_binary(inv_key);
	 printf("\n");
#endif
	 key[i] = inv_key ^ tmp ^ Z[j][k] ^ xconst;
  }
}
/**
 * Compute the number of rounds for Simon and the index of the z-sequence
 * \param word_size word size
 * \param nkey_words number of key words
 * \param zseq_j index of the z-sequence \ref g_simon_zseq
 * \return number of rounds
 */
u32 simon_compute_nrounds(u32 word_size, u32 nkey_words, u32* zseq_j)
{
	u32 nrounds = 0;
  *zseq_j = 6;					  // invalid value (for error-check)

  switch(word_size) {
  case 16:
	 nrounds = 32;
	 *zseq_j = 0;
	 break;
  case 24:
	 nrounds = 36;
	 if(nkey_words == 3) {
		*zseq_j = 0;
	 }
	 if(nkey_words == 4) {
		*zseq_j = 1;
	 }
	 break;
  case 32:
	 if(nkey_words == 3) {
		nrounds = 42;
		*zseq_j = 2;
	 }
	 if(nkey_words == 4) {
		nrounds = 44;
		*zseq_j = 3;
	 }
	 break;
  case 48:
	 if(nkey_words == 2) {
		nrounds = 52;
		*zseq_j = 2;
	 }
	 if(nkey_words == 3) {
		nrounds = 54;
		*zseq_j = 3;
	 }
	 break;
  case 64:
	 if(nkey_words == 2) {
		nrounds = 68;
		*zseq_j = 2;
	 }
	 if(nkey_words == 3) {
		nrounds = 69;
		*zseq_j = 3;
	 }
	 if(nkey_words == 4) {
		nrounds = 72;
		*zseq_j = 4;
	 }
	 break;
  default:
	 break;
  }
  return nrounds;
}

/**
 * Simon encryption procedure.
 * \param key expanded key
 * \param nrounds number of rounds
 * \param x_in first plaintext word
 * \param y_in second plaintext word
 */
void simon_encrypt(u32 key[], u32 nrounds,u32* x_in, u32* y_in)
{
	u32 T = nrounds;
	u32 x = *x_in;
	u32 y = *y_in;

  // left rotation constants
	u32 r1 = 1;
	u32 r2 = 8;
	u32 r3 = 2;

  for(u32 i = 0; i < T; i++) {
#if 0									  // DEBUG
	 printf("[%s:%d] %2d: %8X %8X\n", __FILE__, __LINE__, i, x, y);
#endif
	 u32 tmp = x;
	 //	 u32 f = (LROT(x, r1) & LROT(x, r2))) ^ LROT(x, r3);
	 //x = (y ^ (LROT(x, r1) & LROT(x, r2))) ^ LROT(x, r3) ^ key[i];
	 y = tmp;
  }
#if 0									  // DEBUG
  printf("[%s:%d] %2d: %8X %8X\n", __FILE__, __LINE__, T, x, y);
#endif
  *x_in = x;
  *y_in = y;
}


/*
	u64 x = 0;
	u64 cnt = 0;
	u64 a = 0;
	u64 b = 0;
	u64 c = 0;
	u64 d = 0;

	for(x=0;x<=0xFFFF;x++)
	{
		a = ROTATE_LEFT(x,(rol_a - rol_c),blocksize_len) & 0xFFFF;
		b = RIGHT_rotation_64(x,(rol_c - rol_b),blocksize_len)  & 0xFFFF;
		//b = ROTATE_LEFT(x,15,blocksize_len);
		c = ROTATE_LEFT(x,rol_c,blocksize_len) & 0xFFFF;

		d = ROTATE_LEFT(x,rol_a,blocksize_len) & 0xFFFF;
		d &= ROTATE_LEFT(x,rol_b,blocksize_len) & 0xFFFF;
		d ^= ROTATE_LEFT(x,rol_c,blocksize_len) & 0xFFFF;


		if(0x1== d)  //(a&b) )
		{
			printf("cntX: %x  x: %x  \n", cnt,x);
			cnt += 1;
		}
	}
*/




/*
	u64 blocksize_len =24;  //
	u64 cnt = 0;
	u64 ORvalue =0;
	u64 ANDvalue =0xFFFFFF;  //

	//the Rotation parameter of simon-LIKE cipher, abc//
	rol_a = 5;
	rol_b = 0;
	rol_c = 1;
	rol_delt = 5;

	void printf_itoa(int integer)
	{
		int i=0;
		int len = blocksize_len;// sizeof(integer) * 8;

		for(i=len-1;i>=0;i--)
		{
			if((integer & (1<<i)) != 0)
			{
				printf("1");
			}
			else
			{
				printf("0");
			}

			if((i % 4) == 0)
				printf(" ");
		}
		printf("\n");
	}

*/



/*
	// 向上推r轮
	void ROUNDBACK_2(u64 X,u64 Y,int rounds)
	{
		u64 alpha = Y;  //x
		u64 Sa_LSB8,Sb_LSB8;
		u64 Sa,Sb;
		u64 AB_indx[16] = {0};  // Sa and Sb
		u16 gama_num[]={0};
		u64 gama_0 = 0,gama_1 = 0,gama_2 = 0,gama_3 = 0,gama_4 = 0,gama_5 = 0,gama_6 = 0,gama_7 = 0;
		u16 i=0,i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0,j=0,k=0;
		u64 gama = 0;  //gama.andout.

		//////look up the DDTAB table to find the gama.//////
		Sa = ROTATE_LEFT(alpha,rol_a,blocksize_len);
		Sb = ROTATE_LEFT(alpha,rol_b,blocksize_len);

		//x0char  A0[N0 + 1], T0[N0 + 1], F0[N0 + 1], H0[N0 + 1],  C0[N0 + 1], X0, Y0, I0, L0, Z0;
		 *
		Sa_LSB8 = (Sa) & 0xFF; //8 bit
		Sb_LSB8 = (Sb) & 0xFF; //8 bit
		AB_indx[0] = (Sa_LSB8 << 8) | Sb_LSB8;
		gama_num[0] = DDTAB_8bit[AB_indx[0]][256];
		//x1
		Sa_LSB8 = Sa & 0xFF00; //8 bit
		Sb_LSB8 = ( Sb >> 8) & 0xFF; //8 bit
		AB_indx[1] = Sa_LSB8  | Sb_LSB8;
		gama_num[1] = DDTAB_8bit[AB_indx[1]][256];

		for(i1=0;i1<gama_num[1];i1++)
		{
			gama_1 = DDTAB_8bit[AB_indx[1]][i1];  //gama_num 4-7 bits
			for(i0=0;i0 < gama_num[0];i0++)
			{
				gama_0 = DDTAB_8bit[AB_indx[0]][i0];  //gama_num 0-3 bitssearch_round
				gama = (gama_1 <<8) | gama_0;

				if(SIMON_DP_weight_check(alpha,gama) == 1) //valid gama.
				{
				gama = gama ^ ROTATE_LEFT(alpha,rol_c,blocksize_len) ^ X;

				if(rounds > 2)
				{
					ROUNDBACK_2(Y,gama,rounds-1);
				}
				else
				{
					printf_itoa(gama);
					cnt++;
					ORvalue =ORvalue | gama;
					ANDvalue = ANDvalue & gama;
				}
				}
		}
		}
		}

		u64 X=0x2, Y=0x5;
		u64 Sa,Sb;
		u64 alpha = 0;  //x
		u64 Sa_LSB8,Sb_LSB8;
		u64 AB_indx[16] = {0};  // Sa and Sb
		u64 gama_0 = 0,gama_1 = 0,gama_2 = 0,gama_3 = 0,gama_4 = 0,gama_5 = 0,gama_6 = 0,gama_7 = 0;
		u16 gama_num[]={0};
		u16 i=0,i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0,j=0,k=0;
		u64 gama = 0;  //gama.andout.




		alpha = Y;
		//////look up the DDTAB table to find the gama.//////
		Sa = ROTATE_LEFT(alpha,rol_a,blocksize_len);
		Sb = ROTATE_LEFT(alpha,rol_b,blocksize_len);

		//x0
		Sa_LSB8 = (Sa) & 0xFF; //8 bit
		Sb_LSB8 = (Sb) & 0xFF; //8 bit
		AB_indx[0] = (Sa_LSB8 << 8) | Sb_LSB8;
		gama_num[0] = DDTAB_8bit[AB_indx[0]][256];
		//x1
		Sa_LSB8 = Sa & 0xFF00; //8 bit
		Sb_LSB8 = ( Sb >> 8) & 0xFF; //8 bit
		AB_indx[1] = Sa_LSB8  | Sb_LSB8;
		gama_num[1] = DDTAB_8bit[AB_indx[1]][256];

		for(i1=0;i1<gama_num[1];i1++)
		{
			gama_1 = DDTAB_8bit[AB_indx[1]][i1];  //gama_num 4-7 bits
			for(i0=0;i0 < gama_num[0];i0++)
			{
				gama_0 = DDTAB_8bit[AB_indx[0]][i0];  //gama_num 0-3 bitssearch_round
				gama = (gama_1 <<8) | gama_0;
			if(SIMON_DP_weight_check(alpha,gama) == 1) //valid gama.
				{
				gama = gama ^ ROTATE_LEFT(alpha,rol_c,blocksize_len) ^ X;

				printf_itoa(gama);
				//ROUNDBACK_2(Y,gama,2);
				}
			}
		}
*/


/////////////////

	/*

	//向下推r轮
void ROUNDforward_2(u64 X,u64 Y,int rounds)
{
		u64 alpha = X;  //x
		u64 Sa_LSB8,Sb_LSB8;
		u64 Sa,Sb;
		u64 AB_indx[16] = {0};  // Sa and Sb
		u16 gama_num[]={0};
		u64 gama_0 = 0,gama_1 = 0,gama_2 = 0,gama_3 = 0,gama_4 = 0,gama_5 = 0,gama_6 = 0,gama_7 = 0;
		u16 i=0,i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0,j=0,k=0;
		u64 gama = 0;  //gama.andout.

		//////look up the DDTAB table to find the gama.//////
		Sa = ROTATE_LEFT(alpha,rol_a,blocksize_len);
		Sb = ROTATE_LEFT(alpha,rol_b,blocksize_len);

		//x0
		Sa_LSB8 = (Sa) & 0xFF; //8 bit
		Sb_LSB8 = (Sb) & 0xFF; //8 bit
		AB_indx[0] = (Sa_LSB8 << 8) | Sb_LSB8;
		gama_num[0] = DDTAB_8bit[AB_indx[0]][256];
		//x1
		Sa_LSB8 = Sa & 0xFF00; //8 bit
		Sb_LSB8 = ( Sb >> 8) & 0xFF; //8 bit
		AB_indx[1] = Sa_LSB8  | Sb_LSB8;
		gama_num[1] = DDTAB_8bit[AB_indx[1]][256];

		for(i1=0;i1<gama_num[1];i1++)
		{
			gama_1 = DDTAB_8bit[AB_indx[1]][i1];  //gama_num 4-7 bits
			for(i0=0;i0 < gama_num[0];i0++)
			{
				gama_0 = DDTAB_8bit[AB_indx[0]][i0];  //gama_num 0-3 bitssearch_round
				gama = (gama_1 <<8) | gama_0;

				if(SIMON_DP_weight_check(alpha,gama) == 1) //valid gama.
				{
				gama = gama ^ ROTATE_LEFT(alpha,rol_c,blocksize_len) ^ Y;

				if(rounds > 2)
				{
					ROUNDforward_2(gama,X,rounds-1);
				}
				else
				{
					printf_itoa(gama);
					cnt++;
					ORvalue =ORvalue | gama;
					ANDvalue = ANDvalue & gama;
				}
				}

			}
		}
}


		u64 X=0x5, Y=0x2;
		u64 Sa,Sb;
		u64 alpha = 0;  //x
		u64 Sa_LSB8,Sb_LSB8;
		u64 AB_indx[16] = {0};  // Sa and Sb
		u64 gama_0 = 0,gama_1 = 0,gama_2 = 0,gama_3 = 0,gama_4 = 0,gama_5 = 0,gama_6 = 0,gama_7 = 0;
		u16 gama_num[]={0};
		u16 i=0,i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0,j=0,k=0;
		u64 gama = 0;  //gama.andout.

		alpha = X;
		//////look up the DDTAB table to find the gama.//////
		Sa = ROTATE_LEFT(alpha,rol_a,blocksize_len);
		Sb = ROTATE_LEFT(alpha,rol_b,blocksize_len);

		//x0
		Sa_LSB8 = (Sa) & 0xFF; //8 bit
		Sb_LSB8 = (Sb) & 0xFF; //8 bit
		AB_indx[0] = (Sa_LSB8 << 8) | Sb_LSB8;				//printf_itoa(gama);
		gama_num[0] = DDTAB_8bit[AB_indx[0]][256];
		//x1
		Sa_LSB8 = Sa & 0xFF00; //8 bit
		Sb_LSB8 = ( Sb >> 8) & 0xFF; //8 bit
		AB_indx[1] = Sa_LSB8  | Sb_LSB8;
		gama_num[1] = DDTAB_8bit[AB_indx[1]][256];

		for(i1=0;i1<gama_num[1];i1++)
		{
			gama_1 = DDTAB_8bit[AB_indx[1]][i1];  //gama_num 4-7 bits
			for(i0=0;i0 < gama_num[0];i0++)
			{
				gama_0 = DDTAB_8bit[AB_indx[0]][i0];  //gama_num 0-3 bitssearch_round
				gama = (gama_1 <<8) | gama_0;

				if(SIMON_DP_weight_check(alpha,gama) == 1) //valid gama.
				{
					gama = gama ^ ROTATE_LEFT(alpha,rol_c,blocksize_len) ^ Y;

					printf("gama: %d \n",gama);
				//printf_itoa(gama);
				//ROUNDforward_2(gama,X,2);
			}
			}
		}

		//cnt = gama_num[1] * gama_num[0];
		printf("cnt: %d \n",cnt);
		printf("ORvalue: %d \n",ORvalue);
		printf_itoa(ORvalue);
		printf("ANDvalue: %d \n",ANDvalue);
		printf_itoa(ANDvalue);
*/




