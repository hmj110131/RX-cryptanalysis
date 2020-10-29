/*
 * CHAM.c
 *
 *  Created on: 2019年4月29日
 *      Author: hmj110131
 */

#include "CHAM.h"
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"


 u32 CHAM_even = 0;  //是否搜索的是CHAM的从偶数轮开始的截断差分, 默认是0，即模式是从奇数轮开始的截断路径
 volatile u16 CHAM_Pw_even_cor[80] = {0};
 volatile u16 CHAM_Pw_even_wt[80] = {0};

volatile u32 CHAM_Alpha[80]={0};  //每轮的第模加的输入输出差分
volatile u32 CHAM_Beta[80]={0};
volatile u32 CHAM_Gamma[80]={0};
volatile u32 CHAM_X[4][80]={0}; //每轮的输入差分


volatile u16 CHAM_first_four_wt[4]={0};  //前4轮的相关性重量


volatile u16 CHAM_64_u[80]={0};  //每轮的第模加的输入输出掩码
volatile u16 CHAM_64_v[80]={0};
volatile u16 CHAM_64_w[80]={0};
volatile u16 CHAM_64_X[4][80]={0}; //每轮的输入掩码

volatile u32 CHAM_128_u[80]={0};  //每轮的第模加的输入输出掩码
volatile u32 CHAM_128_v[80]={0};
volatile u32 CHAM_128_w[80]={0};
volatile u32 CHAM_128_X[4][80]={0}; //每轮的输入掩码


u16 CHAM_64_IN_Mask[4] = {0};
u16 CHAM_64_OUT_Mask[4] = {0};
u32 CHAM_128_IN_Mask[4] = {0};
u32 CHAM_128_OUT_Mask[4] = {0};

u8 CHAM_a_beta[256][8][9][32768]= {0};  //alpha carry wt num
u8 CHAM_a_gama[256][8][9][32768]= {0};  //alpha carry wt num
u32 CHAM_a_bg_numb[256][8][9]= {0};   // alpha carry wt
u8 CHAM_a_wt_min[256][8]= {0};   // alpha carry
u8 CHAM_a_wt_max[256][8]= {0};   // alpha carry

u8 msb_CHAM_a_beta[256][8][9][32768]= {0};  //alpha carry wt num
u8 msb_CHAM_a_gama[256][8][9][32768]= {0};  //alpha carry wt num
u32 msb_CHAM_a_bg_numb[256][8][9]= {0};   // alpha carry wt
u8 msb_CHAM_a_wt_min[256][8]= {0};   // alpha carry
u8 msb_CHAM_a_wt_max[256][8]= {0};   // alpha carry





u16 print_CHAM_64_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("round----XV0-------XV1--------XV2--------XV3-----Add-wt \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%04llx    0x%04llx     0x%04llx     0x%04llx    -%d \n",
				r_num,
				CHAM_X[0][r_num],CHAM_X[1][r_num],CHAM_X[2][r_num],CHAM_X[3][r_num],
				P_w[r_num]);
	}
	printf("%02d     0x%04llx    0x%04llx     0x%04llx     0x%04llx    -   \n",
			r_num,
			CHAM_X[0][search_round+1],CHAM_X[1][search_round+1],
			CHAM_X[2][search_round+1],CHAM_X[3][search_round+1]);

	printf("------------------ \n");
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
	{
		printf("%d Round CHAM-64 Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
	}
#else
	{
		printf("%d Round CHAM-64 Total Weight: -%d \n",search_round, CHAM_Pw_even_wt[search_round]);
	}
#endif
	return 1;
}




u32 print_CHAM_128_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("round----XV0-----------XV1------------XV2------------XV3---------wt \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx    -%d \n",
				r_num,
				CHAM_X[0][r_num],CHAM_X[1][r_num],
				CHAM_X[2][r_num],CHAM_X[3][r_num],
				P_w[r_num]);
	}
	printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx    -   \n",
			r_num,
			CHAM_X[0][search_round+1],CHAM_X[1][search_round+1],
			CHAM_X[2][search_round+1],CHAM_X[3][search_round+1]);

	printf("------------------ \n");
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
	{
		printf("%d Round CHAM-128 Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
	}
#else
	{
		printf("%d Round CHAM-128 Total Weight: -%d \n",search_round, CHAM_Pw_even_wt[search_round]);
	}
#endif
	return 1;
}




u16 print_CHAM_64_linear_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("round----XV0-------XV1--------XV2--------XV3-----Add-wt \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%04llx    0x%04llx     0x%04llx     0x%04llx    -%d \n",
				r_num,
				CHAM_64_X[0][r_num],CHAM_64_X[1][r_num],CHAM_64_X[2][r_num],CHAM_64_X[3][r_num],
				P_w[r_num]);
	}
	printf("%02d     0x%04llx    0x%04llx     0x%04llx     0x%04llx    -   \n",
			r_num,
			CHAM_64_X[0][search_round+1],CHAM_64_X[1][search_round+1],
			CHAM_64_X[2][search_round+1],CHAM_64_X[3][search_round+1]);

	printf("------------------ \n");
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{
		printf("%d Round CHAM-64 odd-satrt Total Weight: -%d \n",
				search_round, n_P_bestofR_w[search_round]);
			}
#else
			{
				printf("%d Round CHAM-64 even-start Total Weight: -%d \n",
						search_round, CHAM_Pw_even_cor[search_round]);
			}
#endif

	return 1;
}




u32 print_CHAM_128_linear_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("round----XV0-----------XV1------------XV2------------XV3---------wt \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx    -%d \n",
				r_num,
				CHAM_128_X[0][r_num],CHAM_128_X[1][r_num],
				CHAM_128_X[2][r_num],CHAM_128_X[3][r_num],
				P_w[r_num]);
	}
	printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx    -   \n",
			r_num,
			CHAM_128_X[0][search_round+1],CHAM_128_X[1][search_round+1],
			CHAM_128_X[2][search_round+1],CHAM_128_X[3][search_round+1]);

	printf("------------------ \n");
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{
		printf("%d Round CHAM-128 odd-start Total Weight: -%d \n",
				search_round, n_P_bestofR_w[search_round]);
			}
#else
			{
				printf("%d Round CHAM-128 even-start Total Weight: -%d \n",
						search_round, CHAM_Pw_even_cor[search_round]);
			}
#endif
	return 1;
}


void fixed_Alpha_get_betagamma(void)
{
	u16 A,B,C;  //ALPHA,BETA
	u16 carrybits = 0;
	u16 alpha_inv = 0;
	u16 beta_inv = 0;
	u16 gamma_inv = 0;
	u16 alpha_tmp = 0;
	u16 beta_tmp = 0;
	u16 gamma_tmp= 0;
	u16 eq = 0xFFFF;
	//u64 num = 0;
	u16 wt = 0;
	u16 wt_cnt[9] = {0};
	u16 i = 0;
	u16 WT = 0;
	u16 WT_cnt[9] = {0};
	u8 flag1 = 0,flag2 = 0;

	memset(CHAM_a_beta,0,sizeof(CHAM_a_beta)); //内存清零函数
	memset(CHAM_a_gama,0,sizeof(CHAM_a_gama)); //内存清零函数
	memset(CHAM_a_bg_numb,0,sizeof(CHAM_a_bg_numb)); //内存清零函数
	memset(CHAM_a_wt_min,0xFF,sizeof(CHAM_a_wt_min)); //内存清零函数

	memset(msb_CHAM_a_beta,0,sizeof(msb_CHAM_a_beta)); //内存清零函数
	memset(msb_CHAM_a_gama,0,sizeof(msb_CHAM_a_gama)); //内存清零函数
	memset(msb_CHAM_a_bg_numb,0,sizeof(msb_CHAM_a_bg_numb)); //内存清零函数
	memset(msb_CHAM_a_wt_min,0xFF,sizeof(msb_CHAM_a_wt_min)); //内存清零函数


	for(A=0;A<=0xFF;A++)   //8bit DDTA
	{
		alpha_tmp = A<<1;

		for(carrybits=0; carrybits<8; carrybits++)
		{
			for(i=0; i<9;i++ )
			{
				wt_cnt[i] = 0;
				WT_cnt[i] = 0;
			}

			for(B=0;B<=0xFF;B++)
			{
				//AB = ((A << 8) | B);
				beta_tmp = B<<1;

			for(C=0; C<=0xFF;C++ )
			{
				gamma_tmp = C << 1;

				if((carrybits & 0x4) != 0) // alpha  LSB bit.
				{
					alpha_inv = ~(alpha_tmp | 1) ; // & 0xFFFF;
				}
				else
				{
					alpha_inv = ~alpha_tmp ; // & 0xFFFF;
				}

				if((carrybits & 0x2) != 0) // beta LSB bit.
				{
					beta_inv = (beta_tmp | 1) ; // & 0xFFFF;
				}
				else
				{
					beta_inv = beta_tmp; // & 0xFFFF;
				}

				if((carrybits & 0x1) != 0) // gamma LSB bit.
				{
					gamma_inv = (gamma_tmp | 1);
				}
				else
				{
					gamma_inv = gamma_tmp;
				}

				eq = (alpha_inv ^ beta_inv) & (alpha_inv ^ gamma_inv);
				eq = eq & (A ^ B ^ C ^ beta_inv);

				if( (eq & 0xFF) == 0) //
				{
					//compute wt= Speck_block_wt_compute(A,B,C)
					wt = Speck_block_wt_compute(A,B,C,0xFF);
					CHAM_a_beta[A][carrybits][wt][wt_cnt[wt]] = (u8)B;
					CHAM_a_gama[A][carrybits][wt][wt_cnt[wt]] = (u8)C;
					wt_cnt[wt]++ ;

					WT = Speck_block_wt_compute(A,B,C,0x7F); //HM_weight((~(((~A) ^ B) & ((~A) ^ C))) & 0x7F);
					msb_CHAM_a_beta[A][carrybits][WT][WT_cnt[WT]] = (u8)B;
					msb_CHAM_a_gama[A][carrybits][WT][WT_cnt[WT]] = (u8)C;
					WT_cnt[WT]++;
			}}}


			flag1 = 0;
			flag2 = 0;
			for(i=0; i<9;i++ )
			{
				CHAM_a_bg_numb[A][carrybits][i] = wt_cnt[i];
				msb_CHAM_a_bg_numb[A][carrybits][i] = WT_cnt[i];

				if(wt_cnt[i] != 0)
				{
					if(flag1 ==0)
					{
						CHAM_a_wt_min[A][carrybits] = i;
						flag1 = 1;
					}
					CHAM_a_wt_max[A][carrybits] = i;
				}

				if(WT_cnt[i] != 0)
				{
					if(flag2 ==0)
					{
						msb_CHAM_a_wt_min[A][carrybits] = i;
						flag2 = 1;
					}
					msb_CHAM_a_wt_max[A][carrybits] = i;
				}
			}

			//printf("cDDT[carrybits][AB][256]: %d \n",cDDT[carrybits][AB][256]);
		}
	}
}











