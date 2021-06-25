/*
 * verify.c
 *
 *  Created on: 2021年6月22日
 *      Author: hmj110131
 */


#include <stdlib.h>
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include <time.h>
#include "globalvar.h"
#include "search.h"
#include "RX_crypt.h"
#include "printinfo.h"
#include "XOR_offset_table.h"
#include "verify.h"
#include <nmmintrin.h>


u64 Alpha = 0;
u64 Beta = 0;
u64 Gamma = 0;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int pre_processing_TEST_info(void)
{
	printf("|**************************Verification of Theorem 3.2*************************|\n");
	printf( "Enter the parameters : #Blocksize  #rotaion-parameter  #alpha-beta-gamma \n");

	// Get blocksize of addition//
	printf("Blocksize of Addition: ");
	scanf("%d",&blocksize_len); //Input from terminal
//	printf("%d \n",blocksize_len); ///Fixed in the program

	if( blocksize_len > 16)
	{
		printf("Blocksize of Addition is TOO LARGE. Running with very long time. ");
	}

	/* Get rotation-parameter.*/
	printf("rotation-parameter: ");
	scanf("%d",&RX_k); //Input from terminal
//	printf("%d \n", RX_k); ///Fixed in the program
	if(RX_k >=  blocksize_len)
	{
		printf("Error input rotation-parameter! \n");
		return 0;
	}

	/* Get alpha beta gamma.*/
	printf("Alpha: 0x");
	scanf("%x",&Alpha); //Input from terminal
//	printf("%d \n", Alpha); ///Fixed in the program

	printf("Beta:  0x");
	scanf("%x",&Beta); //Input from terminal
//	printf("%d \n", Beta); ///Fixed in the program

	printf("Gamma: 0x");
	scanf("%x",&Gamma); //Input from terminal
//	printf("%x \n", Gamma); ///Fixed in the program


	Verify_Pr(Alpha, Beta, Gamma);

	return 1;
}


u64 Verify_Pr(u64 alpha, u64 beta, u64 gamma)
{
	u64 T = 0;
	u64 Num = 0;
	u64 total_num = 0;
	long double  Pr_T =0;
	long double  Pr_E =0;

	Pr_T = comput_Pr(alpha, beta, gamma);
	total_num = pow(2, 2*blocksize_len);
	Num = Pr_T * total_num;
	printf("------------------------\n");
	printf("Total of (x,y) pairs: %d \n", total_num);
	printf("Number of (x,y) with theoretical Pr_T: %d \n", Num);
	printf("Pr_T = %Lf \n", Pr_T);

	T = Number_of_xy(alpha, beta, gamma);
	Pr_E = (long double)T/total_num;
	printf("Real number of (x,y) with experiment Pr_E: %d  \n", T);
	printf("Pr_E = %Lf \n", Pr_E);

	return 1;
}


long double comput_Pr(u64 a, u64 b, u64 c)
{
	long double  Pr = 0;
	int zeta_num = 0;
	u64 gama = 0;
	u64 F = 0;
	long double  weight = 0;

	RX_zeta_total = 1 + blocksize_len + RX_k*(blocksize_len - RX_k);  //Calculate the total number of zeta
	printf("------------------------\n");
	printf("RX_zeta_total: %d  \n", RX_zeta_total);
	printf("Constructing RX_Ek[4] ... \n"); // Pr(C_L, C_R)
	Compute_RX_Ek(RX_k, (u16)blocksize_len, Ek); //Compute r(C_L, C_R)
	Ek_min = Compute_RX_Ek_min(Ek);

	printf("Constructing RX_XOR_offset_Table ... \n");
	Compute_RX_V_U(RX_k, blocksize_len, RX_v,RX_v_w, RX_u,RX_u_w);
	Sorted_RX_XOR_offset_Table(RX_k, (u16)blocksize_len,
			RX_v,RX_v_w, RX_u,RX_u_w,
			Ek,
			RX_zeta, RX_zeta_w, RX_zeta_i);

	for(zeta_num=0; zeta_num < RX_zeta_total; zeta_num++)
	{
		gama = c ^ RX_zeta[RX_zeta_i[zeta_num]];
		F = g_validcheck(a,b,gama);
		if(F == 0)
		{
			printf("------------------------\n");
			weight = w_compute(a,b,gama);
			//printf("weight_xor: %Lf \n", weight);
			weight += RX_zeta_w[zeta_num];
			//printf("RX_zeta_w: %f \n", RX_zeta_w[zeta_num]);
			//printf("zeta: %x  \n", RX_zeta[RX_zeta_i[zeta_num]]);
			Pr += pow(2,-weight);
			//weight = log(Pr) / log(2); //prob_w = log(prob) / log(2);
			printf("Pr: %Lf \n",Pr);
		}
	}
	return Pr;
}

/////
u64 Number_of_xy(u64 alpha, u64 beta, u64 gamma)
{
	u64 T = 0;
	u64 x=0, y=0, z=0, zt=0;
	u64 xt=0, yt=0, ztt=0;
	u64 mask =0;
	int i = 0;
	u64 N = 0;
	N = (u64)pow(2,blocksize_len);

	for(i=0; i < blocksize_len; i++ )
	{
		mask ^= 1<<i;
		//printf("mask: %x   \n", mask);
	}

	for(x=0; x < N; x++)
	{
		xt = (ROTATE_LEFT(x, RX_k, blocksize_len) & mask) ^ alpha;
		for(y=0; y < N; y++)
		{
			z = (x + y)& mask;

			yt = (ROTATE_LEFT(y, RX_k, blocksize_len) & mask) ^ beta;
			zt = (xt + yt) & mask;

			ztt = (ROTATE_LEFT(z, RX_k, blocksize_len) & mask) ^ gamma;
			if(zt == ztt)
			{
				//printf("zt: %x    ztt: %x \n",zt, ztt);
				T++;
			}
		}
	}
	return T;
}


/////
u64 g_validcheck(u64 alpha, u64 beta, u64 gamma)
{
	u64 eq = 0xF;
	u64 alpha_inv = 0;
	u64 beta_inv = 0;
	u64 mask =0;
	int i = 0;

	for(i=0; i < blocksize_len; i++ )
	{
		mask ^= 1<<i;
		//printf("mask: %x   \n", mask);
	}

	alpha_inv = ~(alpha << 1);
	beta_inv = beta << 1;

	eq = (alpha_inv ^ beta_inv) & (alpha_inv ^ (gamma <<1));
	eq = eq & (alpha ^ beta ^ gamma ^ beta_inv);
	eq = eq & mask;

	return eq;
}

u64 w_compute(u64 alpha, u64 beta, u64 gamma)
{
	u64 alpha_inv = 0;
	u64 eq = 0xFFFFFF;
	int i = 0;
	u64 mask = 0;

	alpha_inv = ~alpha;
	eq = (alpha_inv ^ beta)	& (alpha_inv ^ gamma );

	for(i=0; i < blocksize_len-1; i++ )
	{
		mask ^= 1<<i;
	}
	eq = (~eq) & mask;

	return (u16)_mm_popcnt_u64(eq);
}


