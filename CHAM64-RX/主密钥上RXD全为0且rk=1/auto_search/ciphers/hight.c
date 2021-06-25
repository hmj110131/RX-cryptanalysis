/*
 * hight.c
 *
 *  Created on: 2018年8月13日
 *      Author: hmj110131
 */

#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "math.h"

#include "hight.h"
#include "ciphers.h"
#include <nmmintrin.h>
#include "printinfo.h"
#include <time.h>
#include "search.h"
#include "globalvar.h"




u8 x_block[33][8] = {0};
u8 y_block[33][8] = {0};

u8 F_0_Table[256] = {0};
u8 F_1_Table[256] = {0};
u8 wt_F_0_Table[256] = {0};
u8 wt_F_1_Table[256] = {0};

u8 wt_x_F1[256][256] = {0};

u8 Var_F_0_Table[256] = {0};
u8 Var_F_1_Table[256] = {0};

u32 Input_Comb_cnt[8] = {0};
u8 Input_Comb_alpha[8][1400000] = {0};
u8 Input_Comb_beta[8][1400000] = {0};
u8 Input_Comb_gamma[8][1400000] = {0};


u32 fixed_alpha_wt_cnt[256][8] = {0};
u8 fixed_alpha_wt_beta[256][8][256] = {0};
u8 fixed_alpha_wt_gamma[256][8][256] = {0};
u16 fixed_alpha_wt_min[256] = {0};
u16 fixed_alpha_wt_max[256] = {0};


u16 Var_XOR_cnt[256] = {0};
u8 Var_XOR_A[256][256] = {0};
u8 Var_XOR_B[256][256] = {0};

u8 P0_wt_HIGHT[32] = {0};
u8 P1_wt_HIGHT[32] = {0};
u8 P2_wt_HIGHT[32] = {0};
u8 P3_wt_HIGHT[32] = {0};

u8 wt_inc_T[8][128] = {0};
u8 wt_inc_T_cnt[8] = {0};


u8 HIGHT_cDDT_v[65536][8][256] = {0};
u8 HIGHT_cDDT_n[65536][8] = {0};
u8 HIGHT_cDDT_wt_max[65536] = {0};
u8 HIGHT_cDDT_wt_min[65536] = {0};


u8 wt_3_n[32][12000000] = {0};
u8 wt_2_n[32][12000000] = {0};
u8 wt_1_n[32][12000000] = {0};
u8 wt_0_n[32][12000000] = {0};
u64 wt_n_max[32] = {0};

u8 HIGHT_x[8] = {0};
u8 HIGHT_y[8] = {0};


////////////////////////////////////////////////////////////////////
//#define ROTATE_LEFT(x,s,n)  (((x) << s) | ((x) >> (n-s)))
u8 rol_left(u8 x, u8 indx)
{
	return ((x << indx) | (x >> (8-indx))) & 0xFF;
}

//////////////////////////////
u8 hight_F_0(u8 x)
{
	/*
	u8 y = 0;
	u8 tmp =x;
	u8 a1 = 0;
	u8 a2 = 0;
	u8 a3 = 0;

	a1 = rol_left(tmp,1);
	a2 = rol_left(tmp,2);
	a3 = rol_left(tmp,7);
	y = a1 ^ a2 ^a3;
	*/

/*
	y = ROTATE_LEFT(x,1,8);
	y ^= ROTATE_LEFT(x,2,8);
	y ^= ROTATE_LEFT(x,7,8);
*/
	return 	(rol_left(x,1) ^ rol_left(x,2) ^ rol_left(x,7)) & 0xFF;
}


u8 hight_F_1(u8 x)
{
	/*
	u8 y = 0;
	y = ROTATE_LEFT(x,3,8) ^ ROTATE_LEFT(x,4,8) ^ ROTATE_LEFT(x,6,8);
	return y;
	*/
	return (rol_left(x,3) ^ rol_left(x,4) ^ rol_left(x,6)) & 0xFF;
}

u8 Construct_F_0(void)
{
	u16 i =0;

	for(i=0; i <= 0xFF;i++)
	{
		F_0_Table[i] = hight_F_0((u8)i);
		wt_F_0_Table[i] = (u8)_mm_popcnt_u64(F_0_Table[i] & 0x7F);

		//printf("F_0_Table[%d]: %x   wt_F_0_Table[%d]: %x \n", i, F_0_Table[i], i, wt_F_0_Table[i]);
	}
	return 1;
}

u8 Construct_F_1(void)
{
	u16 i=0;

	for(i=0; i<= 255;i++)
	{
		F_1_Table[i] = hight_F_1((u8)i);
		wt_F_1_Table[i] = (u8)_mm_popcnt_u64(F_1_Table[i] & 0x7F);

		//printf("F_1_Table[%d]: %x   wt_F_1_Table[%d]: %x \n", i, F_1_Table[i], i, wt_F_1_Table[i] );
	}
	return 1;
}


u8 Consturct_X_F1_wt(void)
{
	u16 i = 0;
	u16 j = 0;

	for(i=0; i<256; i++)
	{
		for(j=0; j<256;j++)
		{
			wt_x_F1[i][j] = (u8)_mm_popcnt_u64((i ^ F_1_Table[j]) & 0x7F);  //

			//printf("wt_x_F1[%d][%d]:  %d  \n", i, j, wt_x_F1[i][j]);
		}
	}
	return 0;
}



u8 Construct_Var_F_0(void)
{
	//u8 cnt = 0;
	u16 x = 0;
	u16 y = 0;



	for(y=0; y<= 255;y++)  //y
	{
		//cnt = 0;
		for(x=0; x<256;x++)
		{
			if(y == hight_F_0((u8)x))
			{
				Var_F_0_Table[(u8)y] = (u8)x;
				//printf("Var_F_0_Table[%x]:  %x  \n", y, x);
				//printf("F_0_Table[%x]:  %x  \n", x,y);
				break;
				//cnt++;
			}
		}
		//Var_F_0_Table[y] = cnt;
	}
	return 1;
}

u8 Construct_Var_F_1(void)
{
	//u8 cnt = 0;
	u16 x = 0;
	u16 y = 0;

	for(y=0; y<256;y++)  //y
	{
		//cnt = 0;
		for(x=0; x<= 255;x++)
		{
			if(y == hight_F_1((u8)x))
			{
				Var_F_1_Table[(u8)y] = (u8)x;
				//printf("Var_F_1_Table[%x]:  %x  \n", y, x);
				break;
				//cnt++;
			}
		}
		//Var_F_1_Table[y] = cnt;
		//printf("%d::: %d   \n ",y, cnt );
	}
	return 1;
}

u8 Construct_Input_Comb(void)
{
	u8 wt = 0;
	u8 M0 = 0;
	u8 N0 = 7;
	char A0[9], T0[9], F0[9], H0[9], C0[9];
	char X0, Y0, I0, L0, Z0;
	int i=0;
	u8 j_abc = 0;
	u8 Input_alpha = 0;
	u8 Input_beta = 0;
	u8 Input_gamma = 0;


	for(j_abc = 0; j_abc < 4;j_abc++)
	{
		if((set_A_4[j_abc] & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = 0x80;
		}
		else
		{
			Input_alpha = 0;  // & Bit_Align;
		}

		if((set_A_4[j_abc] & 0x2) != 0) // beta bit j.
		{
			Input_beta = 0x80;
		}
		else
		{
			Input_beta = 0; // & Bit_Align;
		}

		if((set_A_4[j_abc] & 0x1) != 0) // gamma bit j.
		{
			Input_gamma = 0x80;
		}
		else
		{
			Input_gamma = 0;  // & Bit_Align;
		}

		Input_Comb_alpha[0][Input_Comb_cnt[0]] = Input_alpha;
		Input_Comb_beta[0][Input_Comb_cnt[0]]  = Input_beta;
		Input_Comb_gamma[0][Input_Comb_cnt[0]] = Input_gamma;
		Input_Comb_cnt[0]++;
	}

	for(wt = 1; wt<8; wt++)
	{
		M0 = wt;

		for ( i=0; i<=(N0-M0); i++) A0[i] = 0;
		for ( i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for ( i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1;
		T0[1] = 0;
		F0[N0] = N0 - M0 + 1;
		I0 = N0 - M0; L0 = N0;
		///第一种组合模式
		Construct_Input_Comb_MSB(Input_alpha, Input_beta, Input_gamma,wt,C0);

		do
		{
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
						goto CHANGE_comb;
					}
					if (L0 == N0)
					{
						T0[F0[N0]] = -I0 - 1;
						T0[I0 + 1] = T0[I0];
						I0 = F0[N0];
						F0[N0] = F0[N0] + 1;
						goto CHANGE_comb;
					}
					T0				[L0] = -I0-1;
					T0[I0+1] = T0[I0];
					F0[L0] = F0[L0] + 1;
					I0 = L0;
					goto CHANGE_comb;
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
							goto CHANGE_comb;
						}
						T0[F0[N0]-1] = -I0-1;
						T0[I0+1] = T0[I0];
						I0 = F0[N0] - 1;
						goto CHANGE_comb;
					}
					T0[L0] = -I0 -1; T0[I0 + 1] = T0[I0]; I0 = L0;
					goto CHANGE_comb;
				}
				X0 = N0;
				F0[L0 - 1] = F0[L0];
				F0[N0] = N0;
				L0 = N0;
				if (I0 == N0 - 1)
				{
					I0 = T0[N0 - 1];
					goto CHANGE_comb;
				}
				T0[N0 - 1] = -I0 - 1;
				T0[I0 + 1] = T0[I0];
				I0 = N0 - 1;
	CHANGE_comb:
			A0[X0] = 1;
			A0[Y0] = 0;
			H0[X0] = Z0 = H0[Y0];
			C0[Z0] = X0;
			}
			//剩余组合模式输出
			Construct_Input_Comb_MSB(Input_alpha, Input_beta, Input_gamma,wt,C0);

		} while(1);
	}

/*
	for(i=0; i<8; i++)
	{
		printf("Input_Comb_cnt[%d]:  %d  \n", i, Input_Comb_cnt[i]);
	}
*/

	return 1;
}



u8 Construct_Input_Comb_MSB(u8 alpha, u8 beta, u8 gamma,u8 P_1, char *posi )
{
	u8 j_abc = 0;
	u8 Input_alpha = 0;
	u8 Input_beta = 0;
	u8 Input_gamma = 0;

	for(j_abc = 0; j_abc < 4;j_abc++)   // MSB
	{
		if((set_A_4[j_abc] & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = 0x80;
		}
		else
		{
			Input_alpha = 0;  // & Bit_Align;
		}

		if((set_A_4[j_abc] & 0x2) != 0) // beta bit j.
		{
			Input_beta = 0x80;
		}
		else
		{
			Input_beta = 0; // & Bit_Align;
		}

		if((set_A_4[j_abc] & 0x1) != 0) // gamma bit j.
		{
			Input_gamma = 0x80;
		}
		else
		{
			Input_gamma = 0;  // & Bit_Align;
		}

		Construct_Input_Comb_Middle(Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);

	}
	return 1;
}

u8 Construct_Input_Comb_Middle(u8 alpha,u8 beta,u8 gamma,u8 P_1,char *posi,u8 cur_posi)
{
	u8 j_abc = 0;
	u8 Input_alpha = 0;
	u8 Input_beta = 0;
	u8 Input_gamma = 0;
	u8 indx_tmp = 0;
	u8 bit_i = 1;
	u8 i = 0;
	u16 j_last = 0;


	indx_tmp = posi[cur_posi] - 1;

	// Middle bit positionsof of alpha/beta/gamma
	if( cur_posi > 1)  // Middle bit positionsof of alpha/beta/gamma
	{
		//for(j_abc = 1; j_abc < 7; j_abc++)   // MSB //for speck32/48/64
		for(i = 0; i < 3; i++)   // MSB //for speck96/128
		{
			j_abc = set_A_3[i];  //for speck96/128
			//j_abc = set_AB_6[i];  //for all set

			if((j_abc & 0x4) != 0) // alpha bit j.
			{
				Input_alpha = alpha | (bit_i << indx_tmp);
			}
			else
			{
				Input_alpha = alpha; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
				//Input_alpha = (alpha & ( ~(1 << indx_tmp)) ); // & Bit_Align;
			}

			if((j_abc & 0x2) != 0) // alpha bit j.
			{
				Input_beta = beta | (bit_i << indx_tmp); // & Bit_Align;
			}
			else
			{
				Input_beta = beta; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
				//Input_beta = (beta & ( ~(1 << indx_tmp)) ); // & Bit_Align;
			}

			if((j_abc & 0x1) != 0) // alpha bit j.
			{
				Input_gamma = gamma | (bit_i << indx_tmp); // & Bit_Align;
			}
			else
			{
				Input_gamma = gamma; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
				//Input_gamma = (gamma & ( ~(1 << indx_tmp)) ); // & Bit_Align;
			}


/*
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
			}        //对于speck96和speck128,不考虑导致后续比特全wei 1情况
			else //set B_3 of Middle //if((j_abc==1) || (j_abc==2) || (j_abc==4) )
			if((j_abc==1) || (j_abc==2) || (j_abc==4))
			{
				// setB_3对应到下一个不全同比特位置,全设置为 1 才满足差分有效性条件.
				for(j_last=posi[cur_posi - 1]; j_last < indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
				{
					Input_alpha ^= (bit_i << j_last);
					Input_beta  ^= (bit_i << j_last);
					Input_gamma ^= (bit_i << j_last);
				}
			}
*/
			Construct_Input_Comb_Middle(Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
		}
	}
	else 	//call last position.
	{
		Construct_Input_Comb_Last(alpha,beta,gamma,P_1,posi );
	}

	return 1;
}

u8 Construct_Input_Comb_Last(u8 alpha, u8 beta, u8 gamma,u8 P_1, char *posi )
{
	u8 j_abc = 0;
	u8 Input_alpha = 0;
	u8 Input_beta = 0;
	u8 Input_gamma = 0;
	u8 indx_tmp = 0;
	u8 bit_i = 1;
	u16 j_last = 0;


	indx_tmp = posi[1] - 1;

	// Last bit position of alpha/beta/gamma // set_A_3
	for(j_abc = 0; j_abc < 3;j_abc++)// Last bit position.
	{
		if((set_A_3[j_abc] & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = alpha | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			Input_alpha = alpha; // & ( ~(1 << indx_tmp))) & Bit_Align;
		}

		if((set_A_3[j_abc] & 0x2) != 0) // alpha bit j.
		{
			Input_beta = beta | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			Input_beta = beta; // & ( ~(1 << indx_tmp))) & Bit_Align;
		}

		if((set_A_3[j_abc] & 0x1) != 0) // alpha bit j.
		{
			Input_gamma = gamma | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			Input_gamma = gamma; // & ( ~(1 << indx_tmp)) ) & Bit_Align;
		}
/*
	 	 // 若考虑set_B的情况,则需要在此将最低非全0比特位后面的比特全部清零
		if(posi[1] > 1)
			for(j_last=0; j_last< indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
			{
			Input_alpha &= ~(bit_i << j_last);
			Input_beta  &= ~(bit_i << j_last);
			Input_gamma &= ~(bit_i << j_last);
			}
*/

		Input_Comb_alpha[P_1][Input_Comb_cnt[P_1]] = Input_alpha;
		Input_Comb_beta[P_1][Input_Comb_cnt[P_1]]  = Input_beta;
		Input_Comb_gamma[P_1][Input_Comb_cnt[P_1]] = Input_gamma;
		Input_Comb_cnt[P_1]++;
	}

	return 1;
}

u8 Construct_fixedAlpha_betagamma_Comb(void)
{
	u16 A,B,C;  //ALPHA,BETA,GAMMA
	u8 alpha_tmp = 0;
	u8 beta_tmp = 0;
	u8 gamma_tmp= 0;
	u8 eq = 1;
	u8 wt = 0;
	u8 wt_cnt[8] = {0};
	u8 i = 0;
	u8 flag = 0;
	//u8 A_tmp =0;


	memset(fixed_alpha_wt_min,0xFF,sizeof(fixed_alpha_wt_min)); //初始化fixed_alpha_wt_min为全F。
	memset(fixed_alpha_wt_max,0,sizeof(fixed_alpha_wt_max)); //初始化fixed_alpha_wt_max为全F。



	for(A=0;A<=0xFF;A++)   //For all alpha.
	{
		alpha_tmp = ~(A << 1);

		for(i=0; i<8;i++ )
		{
			wt_cnt[i] = 0;
		}

		for(B=0;B<=0xFF; B++)
		{
			beta_tmp = B<<1;
			for(C=0; C<=0xFF; C++ )
			{
				gamma_tmp = C << 1;

				eq = (alpha_tmp ^ beta_tmp) & (alpha_tmp ^ gamma_tmp);
				eq = eq & (A ^ B ^ C ^ beta_tmp);

				if((eq & 0xFF) == 0) // valid alpha, beta, gamma.
				{
					wt = (u8)_mm_popcnt_u64((~((~A ^ B) & (~A ^ C))) & 0x7F); //faster compute hw.

					fixed_alpha_wt_beta[A][wt][wt_cnt[wt]]  = (u8)B;
					fixed_alpha_wt_gamma[A][wt][wt_cnt[wt]] = (u8)C;
					wt_cnt[wt]++;
				}
			}
		}

		flag = 0; //A_tmp = (u8)(A & 0xFF);
		for(i=0; i<8;i++ )
		{
			fixed_alpha_wt_cnt[(u8)A][i] = wt_cnt[i];

			if(wt_cnt[i] != 0)
			{
				if(flag==0)
				{
					fixed_alpha_wt_min[(u8)A] = i;
					flag = 1;
				}
				fixed_alpha_wt_max[(u8)A] = i;
			}
			//printf("fixed_alpha_wt_cnt[%x][%d]:  %d  \n", A,i,fixed_alpha_wt_cnt[(u8)A][i] );
		}
	}
	return 1;
}

u8 Construct_XOR_var_Table(void)
{
	u16 A,B,C;  //ALPHA ^ BETA = GAMMA
	u16 cnt = 0;

	for(C=0; C<=0xFF;C++ )
	{
		cnt = 0;
		for(B=0;B<=0xFF;B++)
		{
			for(A=0;A<=0xFF;A++)   //
			{
				if(C == (A ^ B))
				{
					Var_XOR_A[C][cnt] = (u8)A;
					Var_XOR_B[C][cnt] = (u8)B;
					cnt++;
				}
			}
		}

		Var_XOR_cnt[C] = cnt;
		//printf("%d::: %d   \n ",C, cnt );
	}
	return 1;
}

u8 Construct_wt_inc_Table(void)
{
	u16 A = 0;
	u8 cnt = 0;
	u8 weight_t = 0;

	for(A = 0; A <= 0xFF; A++)
	{
		weight_t = (u8)_mm_popcnt_u64(A & 0x7F);

		wt_inc_T[weight_t][wt_inc_T_cnt[weight_t]] = A;
		wt_inc_T_cnt[weight_t]++;
	}
/*
	for(cnt=0; cnt<8;cnt++ )
	{
		printf("Weight:  %d ", cnt);
		for(weight_t=0; weight_t< wt_inc_T_cnt[cnt];weight_t++)
		{
			printf("%d  ",wt_inc_T[cnt][weight_t]);
		}
		printf("\n");
	}
*/
	return 1;
}

void ARX_carry_DDTm_8bit(void)
{
	u16 A,B,C;  //ALPHA,BETA
	u16 AB = 0;
	u16 alpha_tmp = 0;
	u16 beta_tmp = 0;
	u16 gamma_tmp= 0;
	u16 eq = 0xFF;
	u8 wt = 0;
	u64 wt_cnt[8] = {0};
	u8 flag = 0;
	u64 i = 0;

	memset(HIGHT_cDDT_wt_min,0xFF,sizeof(HIGHT_cDDT_wt_min)); //初始化HIGHT_cDDT_wt_min为全F。
	memset(HIGHT_cDDT_wt_max,0,sizeof(HIGHT_cDDT_wt_max)); //初始化HIGHT_cDDT_wt_max为全F。


	for(A=0;A <= 0xFF;A++)   //8bit DDTA
	{
		alpha_tmp = A<<1;
		for(B=0;B <= 0xFF;B++)
		{
			AB = ((A << 8) | (B & 0xFF));
			beta_tmp = B<<1;
			for(i=0; i<8;i++ )
			{
				wt_cnt[i] = 0;
			}
			for(C=0; C <= 0xFF;C++ )
			{
				gamma_tmp = C << 1;

				eq = ((~alpha_tmp) ^ beta_tmp) & ((~alpha_tmp) ^ gamma_tmp);
				eq = eq & (A ^ B ^ C ^ beta_tmp);

				if((eq & 0xFF) == 0) //
				{
					wt = (u8)_mm_popcnt_u64((~((~A ^ B) & (~A ^ C))) & 0x7F);

					/////组合输入输出,随着wt增加进行排序
					Input_Comb_alpha[wt][Input_Comb_cnt[wt]] = (u8)A;
					Input_Comb_beta[wt][Input_Comb_cnt[wt]]  = (u8)B;
					Input_Comb_gamma[wt][Input_Comb_cnt[wt]] = (u8)C;
					Input_Comb_cnt[wt]++;
					/////组合输入输出

					HIGHT_cDDT_v[AB][wt][wt_cnt[wt]] = (u8)C;
					wt_cnt[wt]++;
				}
			}

			flag = 0;
			for(i=0; i<8;i++)
			{
				HIGHT_cDDT_n[AB][i] = wt_cnt[i];

				if(wt_cnt[i] != 0)
				{
					if(flag==0)
					{
						HIGHT_cDDT_wt_min[AB] = i;
						flag = 1;
					}
					HIGHT_cDDT_wt_max[AB] = i;
				}
			}
		}}

/*
	for(i=0;i<8;i++)
	{
		printf("Input_Comb_cnt[%d]::: %d \n",i, Input_Comb_cnt[i]);
	}
*/


}

void HIGHT_first_wt_Comb(void)
{
	u8 Wt_0 = 0;
	u8 Wt_1 = 0;
	u8 Wt_2 = 0;
	u8 Wt_3 = 0;
	u8 Wt_n = 0;


	for(Wt_n = 0; Wt_n <= 28; Wt_n++){   //first round max = 8;
		for(Wt_3 = 0; Wt_3 < 8; Wt_3++)
			for(Wt_2 = 0; Wt_2 < 8; Wt_2++)
				for(Wt_1=0;Wt_1 < 8;Wt_1++)
					for(Wt_0 = 0; Wt_0 < 8; Wt_0++)
					{
						if((Wt_3 + Wt_2 + Wt_1 + Wt_0) == Wt_n )
						{
							wt_3_n[Wt_n][wt_n_max[Wt_n]] = Wt_3;
							wt_2_n[Wt_n][wt_n_max[Wt_n]] = Wt_2;
							wt_1_n[Wt_n][wt_n_max[Wt_n]] = Wt_1;
							wt_0_n[Wt_n][wt_n_max[Wt_n]] = Wt_0;
							wt_n_max[Wt_n]++;
						}
			}
		//printf("wt_n_max[%d]::: %d   \n ",Wt_n, wt_n_max[Wt_n] );
	}
}


//构造HIGHT的查找表非常耗时
u8 Construct_HIGHT_Tables(void)
{
	Construct_F_0();
	Construct_F_1();
	Construct_Var_F_0();
	Construct_Var_F_1();
//	Construct_Input_Comb();
	Construct_fixedAlpha_betagamma_Comb();
	//Construct_XOR_var_Table();
	//Construct_wt_inc_Table();
	Consturct_X_F1_wt();
	ARX_carry_DDTm_8bit();
	HIGHT_first_wt_Comb();

/*
	u8 i =0 ;u16 j = 0;
	for(i=0;i<8;i++)
	{
		printf("--------- Input_Comb_cnt[%d]: %d--------- \n" ,i,Input_Comb_cnt[i]);

		for(j=0;j<256;j++)
		printf("fixed_alpha_wt_cnt[%d][%d]: %d \n" ,j, i,fixed_alpha_wt_cnt[j][i]);
	}
*/

	return 1;
}


u8 print_HIGHT_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("---------------------------HIGHT Optima differentail trails------------------------------- \n");  //----gama_cnt \n");
	fprintf(result_print,"---------------------------HIGHT Optima differentail trails------------------------------- \n");
	printf("Round---x7-----x6-----x5-----x4-----x3-----x2-----x1-----x0----wt3---wt2---wt1---wt0---P_wt \n");  //----gama_cnt \n");
	fprintf(result_print,"Round---x7-----x6-----x5-----x4-----x3-----x2-----x1-----x0----wt3---wt2---wt1---wt0---P_wt \n");


	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d    0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x    -%d    -%d    -%d    -%d    -%d \n",
				r_num,x_block[r_num][7],x_block[r_num][6],x_block[r_num][5],x_block[r_num][4],
				x_block[r_num][3],x_block[r_num][2],x_block[r_num][1],x_block[r_num][0],
				P3_wt_HIGHT[r_num],P2_wt_HIGHT[r_num],P1_wt_HIGHT[r_num],P0_wt_HIGHT[r_num],
				P_w[r_num]);
		fprintf(result_print,"%02d    0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x    -%d    -%d    -%d    -%d    -%d \n",
				r_num,x_block[r_num][7],x_block[r_num][6],x_block[r_num][5],x_block[r_num][4],
				x_block[r_num][3],x_block[r_num][2],x_block[r_num][1],x_block[r_num][0],
				P3_wt_HIGHT[r_num],P2_wt_HIGHT[r_num],P1_wt_HIGHT[r_num],P0_wt_HIGHT[r_num],
				P_w[r_num]);
	}
	printf("%02d    0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   null  null  null  null   NULL \n",
			search_round+1,x_block[search_round+1][7],x_block[search_round+1][6],
			x_block[search_round+1][5],x_block[search_round+1][4],x_block[search_round+1][3],
			x_block[search_round+1][2],x_block[search_round+1][1],x_block[search_round+1][0]
			);
	fprintf(result_print,"%02d    0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   0x%02x   null  null  null  null   NULL \n",
			search_round+1,x_block[search_round+1][7],x_block[search_round+1][6],
			x_block[search_round+1][5],x_block[search_round+1][4],x_block[search_round+1][3],
			x_block[search_round+1][2],x_block[search_round+1][1],x_block[search_round+1][0]
			);

	printf("------------------------------------------------------------------------------------------- \n");
	fprintf(result_print,"------------------------------------------------------------------------------------------- \n");
	printf("***************\n");
	fprintf(result_print,"***************\n");
	printf("%d Round HIGHT Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
	fprintf(result_print,"%d Round HIGHT Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
	return 0;
}










































