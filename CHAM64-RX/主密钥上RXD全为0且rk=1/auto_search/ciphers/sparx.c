/*
 * sparx.c
 *
 *  Created on: 2018年8月28日
 *      Author: hmj110131
 */


#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "math.h"

#include "ciphers.h"
#include <nmmintrin.h>
#include "printinfo.h"
#include <time.h>
#include "sparx.h"
#include "globalvar.h"


/////////用于SPARX-64的变量
u16 Comb_16bit_cnt[16] = {0};
u16 Comb_16bit_alpha[16][12000] = {0};
u16 Comb_16bit_beta[16][12000] = {0};
u16 Comb_16bit_gamma[16][12000] = {0};

u16 X_r[4][20] = {0};
u16 wt_l[20] = {0};
u16 wt_r[20] = {0};

u16 X_in_0=0,X_in_1=0,X_in_2=0,X_in_3=0;
u16 X_in_4=0,X_in_5=0,X_in_6=0,X_in_7=0;
u16 Y_out_0=0,Y_out_1=0,Y_out_2=0,Y_out_3=0;
u16 Y_out_4=0,Y_out_5=0,Y_out_6=0,Y_out_7=0;



u16 IN_X_Mask[8] = {0};
u16 OUT_Y_Mask[8] = {0};



u16 IN_L_X[5][4] = {0};

///////////用于SPARX-128的变量
u16 X_128_r[8][20] = {0};
u16 wt_block[4][20] = {0};
u16 IN128_L_X[5][8] = {0};
u16 FIRST = 0;

u16 R1_mask_UVW[12] = {0};
u16 R1_mask_UVW_wt[4] = {0};
u16 R2_mask_UVW[12] = {0};
u16 R2_mask_UVW_wt[4] ={0};
u16 R1_wt_part[4] = {0};
u16 R2_wt_part[4] = {0};


u8 L2_Switch_per3Round(u16 *x0, u16 *x1, u16 *x2, u16 *x3)
{
	u16 y0 = 0;
	u16 y1 = 0;
	u16 y2 = 0;
	u16 y3 = 0;
	u16 tmp = 0;
	u16 tmp_xor = 0;
	u16 x0_tmp = 0;
	u16 x1_tmp = 0;

	y0 = *x0;
	y1 = *x1;
	y2 = *x2;
	y3 = *x3;

	tmp = y0 ^ y1;
	tmp_xor = ((tmp << 8) | (tmp >> 8));

	x0_tmp = tmp_xor ^ y0 ^ y2;
	x1_tmp = tmp_xor ^ y1 ^ y3;

	*x2 = y0;
	*x3 = y1;
	*x0 = x0_tmp;
	*x1 = x1_tmp;


	return 0;
}

u8 L4_Switch_per4Round(u16 *x0, u16 *x1, u16 *x2, u16 *x3, u16 *x4, u16 *x5, u16 *x6, u16 *x7)
{
	u16 y0 = 0;
	u16 y1 = 0;
	u16 y2 = 0;
	u16 y3 = 0;
	u16 y4 = 0;
	u16 y5 = 0;
	u16 y6 = 0;
	u16 y7 = 0;
	u16 tmp = 0;
	u16 tmp_xor = 0;
	u16 x0_tmp = 0;
	u16 x1_tmp = 0;
	u16 x2_tmp = 0;
	u16 x3_tmp = 0;


	y0 = *x0;  // left
	y1 = *x1;
	y2 = *x2;
	y3 = *x3;
	y4 = *x4;  // right
	y5 = *x5;
	y6 = *x6;
	y7 = *x7;

	tmp = y0 ^ y1 ^ y2 ^ y3;
	tmp_xor =  ((tmp << 8) | (tmp >> 8));

	x0_tmp = tmp_xor ^ y2 ^ y4;
	x1_tmp = tmp_xor ^ y1 ^ y5;
	x2_tmp = tmp_xor ^ y0 ^ y6;
	x3_tmp = tmp_xor ^ y3 ^ y7;

	*x4 = y0;
	*x5 = y1;
	*x6 = y2;
	*x7 = y3;
	*x0 = x0_tmp;
	*x1 = x1_tmp;
	*x2 = x2_tmp;
	*x3 = x3_tmp;

	return 0;
}

u8 Construct_16bit_Comb(void)
{
	u8 wt = 0;
	u8 M0 = 0;
	u8 N0 = 7;
	char A0[9], T0[9], F0[9], H0[9], C0[9];
	char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;
	u16 Input_alpha = 0;
	u16 Input_beta = 0;
	u16 Input_gamma = 0;


	for(j_abc = 0; j_abc < 4;j_abc++)
	{
		if((set_A_4[j_abc] & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = 0x8000;
		}
		else
		{
			Input_alpha = 0;  // & Bit_Align;
		}

		if((set_A_4[j_abc] & 0x2) != 0) // beta bit j.
		{
			Input_beta = 0x8000;
		}
		else
		{
			Input_beta = 0; // & Bit_Align;
		}

		if((set_A_4[j_abc] & 0x1) != 0) // gamma bit j.
		{
			Input_gamma = 0x8000;
		}
		else
		{
			Input_gamma = 0;  // & Bit_Align;
		}

		Comb_16bit_alpha[M0][Comb_16bit_cnt[M0]] = Input_alpha;
		Comb_16bit_beta[M0][Comb_16bit_cnt[M0]]  = Input_beta;
		Comb_16bit_gamma[M0][Comb_16bit_cnt[M0]] = Input_gamma;
		Comb_16bit_cnt[M0]++;
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
	return 1;
}

u8 Construct_16bit_Comb_MSB(u16 alpha, u16 beta, u16 gamma,u8 P_1, char *posi )
{
	u16 j_abc = 0;
	u16 Input_alpha = 0;
	u16 Input_beta = 0;
	u16 Input_gamma = 0;

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

u8 Construct_16bit_Comb_Middle(u16 alpha,u16 beta,u16 gamma,u8 P_1,char *posi,u8 cur_posi)
{
	u16 j_abc = 0;
	u16 Input_alpha = 0;
	u16 Input_beta = 0;
	u16 Input_gamma = 0;
	u16 indx_tmp = 0;
	u16 bit_i = 1;
	u16 i = 0;
	u16 j_last = 0;


	indx_tmp = posi[cur_posi] - 1;

	// Middle bit positionsof of alpha/beta/gamma
	if( cur_posi > 1)  // Middle bit positionsof of alpha/beta/gamma
	{
		//for(j_abc = 1; j_abc < 7; j_abc++)   // MSB //for speck32/48/64
		for(i = 0; i < 6; i++)   // MSB //for speck96/128
		{
			//j_abc = set_A_3[i];  //for speck96/128
			j_abc = set_AB_6[i];  //for all

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
					Input_alpha &= ~(1 << j_last);
					Input_beta  &= ~(1 << j_last);
					Input_gamma &= ~(1 << j_last);
				}
			}     //对于speck96和speck128,不考虑导致后续比特全wei 1情况
			else //set B_3 of Middle //if((j_abc==1) || (j_abc==2) || (j_abc==4) )
			*/
			if((j_abc==1) || (j_abc==2) || (j_abc==4))
			{
				// setB_3对应到下一个不全同比特位置,全设置为 1 才满足差分有效性条件.
				for(j_last=posi[cur_posi - 1]; j_last < indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
				{
					Input_alpha ^= (1 << j_last);
					Input_beta  ^= (1 << j_last);
					Input_gamma ^= (1 << j_last);
				}
			}

			Construct_Input_Comb_Middle(Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
		}
	}
	else 	//call last position.
	{
		Construct_Input_Comb_Last(alpha,beta,gamma,P_1,posi );
	}

	return 1;
}

u8 Construct_16bit_Comb_Last(u16 alpha, u16 beta, u16 gamma,u8 P_1, char *posi )
{
	u16 j_abc = 0;
	u16 Input_alpha = 0;
	u16 Input_beta = 0;
	u16 Input_gamma = 0;
	u16 indx_tmp = 0;
	u16 bit_i = 1;
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

	 	 // 若考虑set_B的情况,则需要在此将最低非全0比特位后面的比特全部清零
		if(posi[1] > 1)
			for(j_last=0; j_last< indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
			{
			Input_alpha &= ~(bit_i << j_last);
			Input_beta  &= ~(bit_i << j_last);
			Input_gamma &= ~(bit_i << j_last);
			}

		Comb_16bit_alpha[P_1][Comb_16bit_cnt[P_1]] = Input_alpha;
		Comb_16bit_beta[P_1][Comb_16bit_cnt[P_1]]  = Input_beta;
		Comb_16bit_gamma[P_1][Comb_16bit_cnt[P_1]] = Input_gamma;
		Comb_16bit_cnt[P_1]++;


	}

	return 1;
}






u16 print_SPARX64_resoult(u16 search_round)
{
	u16 r_num = 0;
	u16 i_L = 0;

	printf("round----Xr0-------Xr1--------Xr2-------Xr3----w_l----w_r----Pw \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%04x    0x%04x    0x%04x    0x%04x    -%d    -%d     -%d \n",
				r_num,X_r[0][r_num],X_r[1][r_num],X_r[2][r_num],X_r[3][r_num],wt_l[r_num],wt_r[r_num],P_w[r_num]); //,gama_cnt[r_num]);
	}

	printf("%02d     0x%04x    0x%04x    0x%04x    0x%04x    -     -    NULL \n",
			search_round+1,X_r[0][search_round+1],X_r[1][search_round+1],X_r[2][search_round+1],X_r[3][search_round+1]); //,gama_cnt[r_num]);
	printf("------------------ \n");
	for(i_L = 1; i_L <= search_round/3;i_L++ )
	{
		printf("%02d-L   0x%04x    0x%04x    0x%04x    0x%04x    -     -      - \n",
				3*i_L,IN_L_X[i_L][0],IN_L_X[i_L][1],IN_L_X[i_L][2],IN_L_X[i_L][3]); //,gama_cnt[r_num]);
	}

	printf("------------------ \n");
	printf("%d Round SPARX64 Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);

	return 1;
}



u16 print_SPARX_128_resoult(u16 search_round)
{
	u16 r_num = 0;
	u16 i_L = 0;

	printf("round---Xr0----Xr1-----Xr2-----Xr3-----Xr4-----Xr5-----Xr6-----Xr7----w_0---w_1--w_2--w_3----Pw \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d   0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x    -%d   -%d   -%d   -%d     -%d \n",
				r_num,
				X_128_r[0][r_num],X_128_r[1][r_num],X_128_r[2][r_num],X_128_r[3][r_num],
				X_128_r[4][r_num],X_128_r[5][r_num],X_128_r[6][r_num],X_128_r[7][r_num],
				wt_block[0][r_num],wt_block[1][r_num],wt_block[2][r_num],wt_block[3][r_num],
				P_w[r_num]); //,gama_cnt[r_num]);
	}

	printf("%02d   0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x    -    -    -    -     NULL \n",
			search_round+1,
			X_128_r[0][search_round+1],X_128_r[1][search_round+1],X_128_r[2][search_round+1],X_128_r[3][search_round+1],
			X_128_r[4][search_round+1],X_128_r[5][search_round+1],X_128_r[6][search_round+1],X_128_r[7][search_round+1]
			);

	printf("------------------ \n");
	for(i_L = 1; i_L <= search_round/4;i_L++ )
	{
		printf("%02d-L 0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x  0x%04x    -    -    -    -      - \n",
				4*i_L,
				IN128_L_X[i_L][0],IN128_L_X[i_L][1],IN128_L_X[i_L][2],IN128_L_X[i_L][3],
				IN128_L_X[i_L][4],IN128_L_X[i_L][5],IN128_L_X[i_L][6],IN128_L_X[i_L][7]
				);
	}


	printf("------------------ \n");
	printf("%d Round SPARX-128 Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);

	return 1;
}












