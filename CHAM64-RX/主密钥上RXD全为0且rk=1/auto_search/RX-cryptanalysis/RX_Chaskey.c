/*
 * RX_Chaskey.c
 *
 *  Created on: 2020年7月29日
 *      Author: hmj110131
 */
#include <stdlib.h>
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include <time.h>
#include "search.h"
#include "RX_crypt.h"
#include "printinfo.h"
#include "RX_Chaskey.h"
#include "XOR_offset_table.h"
#include <nmmintrin.h>
#include "chaskey.h"


u16 RX_wt_block[4][20] = {0};

////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
volatile float RX_Chaskey_P_hr_w[25] = {0}; //每半轮对应RX差分概率重量
volatile float RX_Chaskey_P_hr_w_tmp[25] = {0}; //每轮对应RX差分概率重量
volatile float RX_Chaskey_ADD_P_xor_w[2][25] = {0}; //每轮中，4个模块加法的对应RX差分概率重量
volatile float RX_Chaskey_ADD_P_xor_w_tmp[2][25] = {0}; //每轮中，4个模块加法的对应RX差分概率重量
volatile u32 Chaskey_R1_ADD_i_RX_abc[2][3] = {0}; //4个模块加法的对应abc
volatile u32 RX_X[4][25] = {0};; //每轮中，4个分支的输入RX差分
volatile u32 RX_X_tmp[4][25] = {0}; //每轮中，4个分支的输入RX差分
volatile u32 Chaskey_Zeta[2][25]={0}; //每半轮的输zeta
volatile u32 Chaskey_zeta_tmp[2][25]={0}; //每半轮的临时zeta
volatile float Chaskey_Pr_zeta[2][25]={0}; //每半轮RX差分概率
volatile float Chaskey_Pr_zeta_tmp[2][25]={0}; //每轮RX差分概率



////////////////////////////////搜索Chaskey的最优RX特征的代码//////////////////////////////
u32 Chaskey_RX_trail_search_entry(u16 search_round)
{
	u32 best = 0;
	u16 i = 0;
	FILE* CHASKEY_RX_Bn_wt;


	printf("--Note: For Chaskey, the 'round' here means half-round.--\n");
	if(sc_blocksize != 128)
	{
		printf("The block size of CHASKEY should be 128 bits. \n");
		return 0;
	}

	// For chaskey modular additon is 32 bits.
	blocksize_len = 32;
	Bit_Align = 0xFFFFFFFF;
	ValueMax_Align = 0x7FFFFFFF;
	V_MSB = 0x80000000;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	RX_zeta_total = 1 + blocksize_len + RX_k*(blocksize_len -RX_k);  //计算Zeta的总个数
	printf("RX_zeta_total: %d  \n",RX_zeta_total);
	printf("Constructing RX_Ek[4] ... \n");
	Compute_RX_Ek(RX_k, (u16)blocksize_len,  Ek); //计算ek的4个值
	printf("k=%d  RX_Ek[0]: %f RX_Ek[1]: %f RX_Ek[2]: %f RX_Ek[3]: %f \n",
			RX_k,Ek[0],Ek[1],Ek[2],Ek[3]);
	Ek_min = Compute_RX_Ek_min(Ek);
	printf("Constructing RX_XOR_offset_Table ... \n");
	Compute_RX_V_U(RX_k, blocksize_len, RX_v,RX_v_w, RX_u,RX_u_w);
	Sorted_RX_XOR_offset_Table(RX_k, (u16)blocksize_len,
			RX_v,RX_v_w, RX_u,RX_u_w,
			Ek,
			RX_zeta, RX_zeta_w, RX_zeta_i);
	//-------------取前r-1轮最优概率重量，搜索期望概率重量的最大值----------------
	CHASKEY_RX_Bn_wt = fopen ("../tmp/Chaskey_RX_Bn_wt.xlsx", "a+"); //  "w+"); //
	for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
	{
		fscanf(CHASKEY_RX_Bn_wt,"%f", &RX_P_bestofR_w[i]);  //read.
	}
	//-------记录时间-----
		time_xor_talbe = clock();
		run_time =  (time_xor_talbe - time_start) / CLOCKS_PER_SEC;
		printf("Time of Preprocessing info: %.2f seconds.  \n", run_time);
		printf("|--------------------------------------------------|\n");
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		printf("Constructing cDDT Lookup Tables ... \n");
		//构造cDDT
		ARX_carry_DDTm_construct();  // FOR cDDT
		//fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
		time_ARX_DDT = clock();
		run_time =  (time_ARX_DDT - time_xor_talbe) / CLOCKS_PER_SEC;
		printf("Time of Construct Chaskey cDDT: %.2f seconds.  \n", run_time);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		int RX_Expected_w_int = 0;
			RX_Expected_w_int = (int)RX_P_bestofR_w[search_round - 1];
			RX_Expected_w = RX_Expected_w_int;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		printf("*****************************Searching processing**************************|\n");
		do{
			RX_Expected_w = RX_Expected_w + 1; // 从r-1轮的概率重量整数部分+1开始
			printf("Update: %d   RX_Bn_w: %f  --Max expected weight-- \n",
					Num_Bn_Update, RX_Expected_w);

			// Search Entry.
			best = Chaskey_RX_ADD_wt_inc(search_round); //效率高

			time_Round = clock();
			run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
			printf("Time: %.2f seconds.\n", run_time);

		}while( Num_Bn_Update == 0 );   //end conditions. // 概率重量逼近次数大于1，且返回值表示搜索到最优路径

		fprintf(CHASKEY_RX_Bn_wt,"%f \r\n",RX_Expected_w);  // write.RX_P_bestofR_w[search_round]
		printf("**************************************************************************|\n");
		print_Chaskey_RX_resoult(search_round);  //打印

		time_finish = clock();
		run_time = (float)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
		printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours == %.2f days. \n", run_time,run_time/60.0,run_time/3600.0,run_time/86400.0 );
		printf("Auto-search %d-Round Chaskey optimal RX-differential trails END! \n", search_round);
		printf("|************************************************************************|\n");
		fclose(CHASKEY_RX_Bn_wt);

	return 0;
}


u32 Chaskey_RX_ADD_wt_inc(u16 search_round)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 wt_round_1 = 0;
	u16 wt_Block_0 = 0,wt_Block_1 = 1; //CHASKEy的RX差分特征，只需要构造首轮前半轮的两个模架
	u16 i = 0;

	// 限定SPARX-128的首轮的概率重量 //blocksize_len-1
	for(wt_round_1 = 1; wt_round_1 < blocksize_len-1; wt_round_1++)  //0::n-1  //各个分块最多15比特的概率重量:61
	{
		// thr_d : the bits number of alpha/beta/gamma that the position not equal at the same time.
		if((wt_round_1 + 2*Ek_min + RX_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
		{
			return 0;  //Return the upper procedure::
		}
		P_w[1] = wt_round_1;  //P_w为每半轮中，两个模加的XOR差分概率重量的和

		for(wt_Block_0=0; wt_Block_0 <= wt_round_1; wt_Block_0++) //15
		{
			for(wt_Block_1=0; wt_Block_1 <= wt_round_1; wt_Block_1++) //15
			{
				if((wt_Block_0 + wt_Block_1) == wt_round_1)
				{
					RX_Chaskey_ADD_P_xor_w_tmp[0][1] = wt_Block_0;
					RX_Chaskey_ADD_P_xor_w_tmp[1][1] = wt_Block_1;

					//根据第1轮的前半轮中两个模加的，每个模加的RX差分概率重量构造对应的输入输出XOR差分
					best = Chaskey_RX_ADD_i_tuples(search_round,0);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
{
	////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	if(best == 2)
	{
		updated = 2;
		return updated;
	}
	else
	{
		updated = 1;  //准备通知已经更新
	}
}
#endif
	}}}}
	return best;
}


//ADD_part = 0 /1 .
u32 Chaskey_RX_ADD_i_tuples(u16 search_round,u16 ADD_part)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 thr_d = 0;
	u32 Input_alpha = 0;
	u32 Input_beta = 0;
	u32 Input_gamma = 0;

	u32 zeta_num = 0;
	u32 updata_cur_num = 0;

	u32 M0 = 0;
	u32 N0 = blocksize_len-1;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;

	updata_cur_num = Num_Bn_Update;

	///////////////////round 1 entry/////////////////////
	///////////2个模加的输入输出XOR差分构造完毕，则调到第1轮的下半轮//////////////
	if(ADD_part >= 2)
	{
		best = Chaskey_RX_comb_RX(search_round);
		return best;
	}


	M0 = RX_Chaskey_ADD_P_xor_w_tmp[ADD_part][1];
	thr_d = M0;
	//P_w[1] = M0; //相关性重量的记录数组

	if(M0 == 0)  /// abd [0 -> n-2] 没有完全不相同的比特,且MSB之后只能为0
	{
		for(j_abc = 0; j_abc < 4;j_abc++)
		{
			if((set_A_4[j_abc] & 0x4) != 0) // alpha bit j.
			{
				Input_alpha = V_MSB;
			}
			else
			{
				Input_alpha = 0;  // & Bit_Align;
			}

			if((set_A_4[j_abc] & 0x2) != 0) // beta bit j.
			{
				Input_beta = V_MSB;
			}
			else
			{
				Input_beta = 0; // & Bit_Align;
			}

			if((set_A_4[j_abc] & 0x1) != 0) // gamma bit j.
			{
				Input_gamma = V_MSB;
			}
			else
			{
				Input_gamma = 0;  // & Bit_Align;
			}

			Chaskey_R1_ADD_i_RX_abc[ADD_part][0] = Input_alpha;
			Chaskey_R1_ADD_i_RX_abc[ADD_part][1] = Input_beta;
			Chaskey_R1_ADD_i_RX_abc[ADD_part][2] = Input_gamma;

			//直接跳转到第1轮各个部分加法的内部构造
			best = Chaskey_RX_ADD_i_tuples(search_round,ADD_part+1);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(best == 2)
		{
			updated = 2;
			return updated;
		}
		else
		{
			updated = 1;  //准备通知已经更新
		}
	}
#endif
	}
}
	else  //M0 > 0.
	{
		for ( i=0; i<=(N0-M0); i++) A0[i] = 0;
		for ( i=N0-M0+1; i<=N0; i++) A0[i] = 1;
		for ( i = 1; i<=M0; i++) { C0[i] = N0 - M0 + i; H0[N0-M0+i] = i; }
		T0[N0-M0] = -1;
		T0[1] = 0;
		F0[N0] = N0 - M0 + 1;
		I0 = N0 - M0; L0 = N0;


		best = RX_Chaskey_input_MSB(search_round,
				Input_alpha,
				Input_beta,
				Input_gamma,thr_d,C0,ADD_part);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(best == 2)
		{
			updated = 2;
			return updated;
		}
		else
		{
			updated = 1;  //准备通知已经更新
		}
	}
#endif

	////////////////////
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

	/////////////////////////////////
						best = RX_Chaskey_input_MSB(search_round,
								Input_alpha,
								Input_beta,
								Input_gamma,thr_d,C0,ADD_part);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(best == 2)
		{
			updated = 2;
			return updated;
		}
		else
		{
			updated = 1;  //准备通知已经更新
		}
	}
#endif
	///////////////////////////////////////
						} while(1);
					}
	return updated;
}


u16 RX_Chaskey_input_MSB
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi,u16 ADD_part )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 Input_alpha = 0;
	u32 Input_beta = 0;
	u32 Input_gamma = 0;
	u16 j_abc = 0;


	// MSB of of alpha/beta/gamma
	for(j_abc = 0; j_abc < 8;j_abc++)   // MSB
	{
		//set A of MSB.
		if((j_abc==0) || (j_abc==3) || (j_abc==5) || (j_abc==6) )
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

			state = RX_Chaskey_input_Middle(search_round,
					Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1,ADD_part);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(state == 2)
		{
			updated = 2;
			return updated;
		}
		else
		{
			updated = 1;  //准备通知已经更新
		}
	}
#endif
		}
		else //set B of MSB //if((j_abc==1) || (j_abc==2) || (j_abc==4) || (j_abc==7) )
		{
			//对于speck96和speck128,不考虑导致后续比特全确定情况
			//若为{0,3,5,6}则除MSB外其它后面比特全为0
			//若为{1、2、4,7}则除MSB外其它后面比特全为1
		if((j_abc & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = Bit_Align;
		}
		else
		{
			Input_alpha = ValueMax_Align;  // & Bit_Align;
		}

		if((j_abc & 0x2) != 0) // beta bit j.
		{
			Input_beta = Bit_Align;
		}
		else
		{
			Input_beta = ValueMax_Align; // & Bit_Align;
		}

		if((j_abc & 0x1) != 0) // gamma bit j.
		{
			Input_gamma = ValueMax_Align;
		}
		else
		{
			Input_gamma = Bit_Align;  // & Bit_Align;
		}

		state = RX_Chaskey_input_Middle(search_round,
				Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1,ADD_part);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(state == 2)
		{
			updated = 2;
			return updated;
		}
		else
		{
			updated = 1;  //准备通知已经更新
		}
	}
#endif
		}
}
	return updated;
}



u16 RX_Chaskey_input_Middle
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi, u16 cur_posi,u16 ADD_part)
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 Input_alpha = 0;
	u32 Input_beta = 0;
	u32 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u32 indx_tmp = 0;
	//u16 i = 0;
	u32 bit_i = 1;


	indx_tmp = posi[cur_posi] - 1;

// 考虑到产生概率重量的位置为全部可能的状态{3/5/6,1/2/4}
// Middle bit positionsof of alpha/beta/gamma
if( cur_posi > 1)  // Middle bit positionsof of alpha/beta/gamma
{
	//只有{3/5/6,1/2/4}会导致该位置产生概率重量
	for(j_abc = 1; j_abc < 7; j_abc++)   //概率重量位置的可能取值组合
	{
		//判断该位置对应的alpha、beta、gamma的比特值为0还是1；
		if((j_abc & 0x4) != 0) // alpha bit j.
		{
			Input_alpha = alpha | (bit_i << indx_tmp);
		}
		else
		{
			//Input_alpha = alpha; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
			Input_alpha = (alpha & ( ~(bit_i << indx_tmp)) ); // & Bit_Align;
		}

		if((j_abc & 0x2) != 0) // alpha bit j.
		{
			Input_beta = beta | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			//Input_beta = beta; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
			Input_beta = (beta & ( ~(bit_i << indx_tmp)) ); // & Bit_Align;
		}

		if((j_abc & 0x1) != 0) // alpha bit j.
		{
			Input_gamma = gamma | (bit_i << indx_tmp); // & Bit_Align;
		}
		else
		{
			//Input_gamma = gamma; // & ( ~(1 << indx_tmp)) ); // & Bit_Align;
			Input_gamma = (gamma & ( ~(bit_i << indx_tmp)) ); // & Bit_Align;
		}


		//判断该位置是属于{3/5/6}，后续比特位置全清0
		if( (j_abc==3) || (j_abc==5) || (j_abc==6) )
		{
			// setA对应到下一个不全同比特位置,全设置为0才满足差分有效性条件.
			for(j_last=posi[cur_posi - 1]; j_last < indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
			{
				Input_alpha &= ~(bit_i << j_last);
				Input_beta  &= ~(bit_i << j_last);
				Input_gamma &= ~(bit_i << j_last);
			}
		}
		else //判断该位置是属于{1/2/4}，后续比特位置全部置1；
		{
			// setB_3对应到下一个不全同比特位置,全设置为 1 才满足差分有效性条件.
			for(j_last=posi[cur_posi - 1]; j_last < indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
			{
				Input_alpha ^= (bit_i << j_last);
				Input_beta  ^= (bit_i << j_last);
				Input_gamma ^= (bit_i << j_last);
			}
		}


/*  //只考虑了产生概率重量的比特位置为{3/5/6}的情况，可能会忽略部分组合
//indx_tmp = 1 << (posi[cur_posi] - 1);
indx_tmp = posi[cur_posi] - 1;

// Middle bit positionsof of alpha/beta/gamma
if( cur_posi > 1)  // Middle bit positionsof of alpha/beta/gamma
{

	//for(j_abc = 1; j_abc < 7; j_abc++)   // MSB //for speck32/48/64
	for(i = 0; i < 3; i++)   // MSB //for speck96/128
	{
		j_abc = set_A_3[i];  //for speck96/128

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
*/
		state = RX_Chaskey_input_Middle(search_round,
				Input_alpha,Input_beta, Input_gamma,
				P_1,posi,cur_posi-1,ADD_part);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
{
	////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	if(state == 2)
	{
		updated = 2;
		return updated;
	}
	else
	{
		updated = 1;  //准备通知已经更新
	}
}
#endif
	}
}
else 	//call last position.
{
	state = RX_Chaskey_input_Last(search_round,
			alpha,beta,	gamma,P_1,posi,ADD_part );
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
{
	////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	if(state == 2)
	{
		updated = 2;
		return updated;
	}
	else
	{
		updated = 1;  //准备通知已经更新
	}
}
#endif
}
return updated;
}

u16 RX_Chaskey_input_Last
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi,u16 ADD_part )
{
u16 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
u32 Input_alpha = 0;
u32 Input_beta = 0;
u32 Input_gamma = 0;
u16 j_abc = 0;
u16 j_last = 0;
u32 indx_tmp = 0;
u32 bit_i = 1;

u32 zeta_num = 0;
u32 updata_cur_num = 0;

indx_tmp = posi[1] - 1;


// 最后的一个产生概率重量的比特位置只能为{3/5/6}，在该位置之后的更低的比特位置都全部清0；
for(j_abc = 0; j_abc < 3;j_abc++)// Last bit position.
{
	if((set_A_3[j_abc] & 0x4) != 0) // alpha bit j.
	{
		Input_alpha = alpha | (bit_i << indx_tmp); // & Bit_Align;
	}
	else
	{
		//Input_alpha = alpha; // & ( ~(1 << indx_tmp))) & Bit_Align;
		Input_alpha = (alpha & ( ~(bit_i << indx_tmp))) & Bit_Align;
	}

	if((set_A_3[j_abc] & 0x2) != 0) // alpha bit j.
	{
		Input_beta = beta | (bit_i << indx_tmp); // & Bit_Align;
	}
	else
	{
		//Input_beta = beta; // & ( ~(1 << indx_tmp))) & Bit_Align;
		Input_beta = (beta & ( ~(bit_i << indx_tmp))) & Bit_Align;
	}

	if((set_A_3[j_abc] & 0x1) != 0) // alpha bit j.
	{
		Input_gamma = gamma | (bit_i << indx_tmp); // & Bit_Align;
	}
	else
	{
		//Input_gamma = gamma; // & ( ~(1 << indx_tmp)) ) & Bit_Align;
		Input_gamma = (gamma & ( ~(bit_i << indx_tmp)) ) & Bit_Align;
	}

	 	 // 需要在此将产生概率重量位置更低的比特位置全部清零
	if(posi[1] > 1)  //若产生概率重量的最后一个位置是LSB，则不用清0
		for(j_last=0; j_last< indx_tmp;j_last++ )  // Last position of alpha/beta/gamma set to all 0s.
		{
			Input_alpha &= ~(bit_i << j_last);
			Input_beta  &= ~(bit_i << j_last);
			Input_gamma &= ~(bit_i << j_last);
		}

	Chaskey_R1_ADD_i_RX_abc[ADD_part][0] = Input_alpha;
	Chaskey_R1_ADD_i_RX_abc[ADD_part][1] = Input_beta;
	Chaskey_R1_ADD_i_RX_abc[ADD_part][2] = Input_gamma;

	//直接跳转到第1轮各个部分加法的内部构造
	state = Chaskey_RX_ADD_i_tuples(search_round,ADD_part+1);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
if(state != 0 )  //判断本轮是否要更新RX_Expected_w
{
	if(state ==2)
	{
		updated = 2;  //准备通知上一轮最新的RX_Expected_w
	return updated;
	}
	else
	{
		updated = 1;
	}
}
#endif
}
return updated;
}


////%%%%%%%%%%%%%%%%%%%%%%%
u32 Chaskey_RX_comb_RX(u16 search_round)
{
	u32 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新

	u32 out_X[4];
	u32 out_RX[4];

	u32 zeta_num_0 = 0,zeta_num_1 = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;

	w_cmp_rx = RX_Expected_w - P_w[1] -RX_P_bestofR_w[search_round -1];


	//首轮的输入RX差分
	RX_X_tmp[0][1] = Chaskey_R1_ADD_i_RX_abc[0][0]; //x0 =alpha
	RX_X_tmp[1][1] = Chaskey_R1_ADD_i_RX_abc[0][1]; //x1= beta
	RX_X_tmp[2][1] = Chaskey_R1_ADD_i_RX_abc[1][0]; //x2 =alpha
	RX_X_tmp[3][1] = Chaskey_R1_ADD_i_RX_abc[1][1]; //x3= beta

	//输出XOR差分，注意Chaskey的x的排序，中间轮最左边为X1
	out_X[1] = Chaskey_R1_ADD_i_RX_abc[0][2] ^ chaskey_rol_left_5(RX_X_tmp[1][1]);
	out_X[2] = chaskey_rol_left_16(Chaskey_R1_ADD_i_RX_abc[0][2]);
	out_X[0] = Chaskey_R1_ADD_i_RX_abc[1][2]; // OUT_X0 = 第二个莫加的gamma
	out_X[3] = out_X[0] ^ chaskey_rol_left_8(RX_X_tmp[3][1]);

	//用zeta补偿XOR差分，得到RX差分
	for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
	{
		for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
		{
			if((RX_zeta_w[zeta_num_0] + RX_zeta_w[zeta_num_1]) >= w_cmp_rx)
			{
				return 0;
			}

			Chaskey_Pr_zeta_tmp[0][1] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
			Chaskey_Pr_zeta_tmp[1][1] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
			Chaskey_zeta_tmp[0][1] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
			Chaskey_zeta_tmp[1][1] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1
			RX_Chaskey_P_hr_w_tmp[1] = P_w[1] + Chaskey_Pr_zeta_tmp[0][1]
											+ Chaskey_Pr_zeta_tmp[1][1];
			RX_P_sumofR_w[1] = RX_Chaskey_P_hr_w_tmp[1];


			out_RX[1] = out_X[1] ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[2] = out_X[2] ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[0] = out_X[0] ^ RX_zeta[RX_zeta_i[zeta_num_1]];
			out_RX[3] = out_X[3] ^ RX_zeta[RX_zeta_i[zeta_num_1]];

			RX_X_tmp[0][2] = out_RX[0];
			RX_X_tmp[1][2] = out_RX[1];
			RX_X_tmp[2][2] = out_RX[2];
			RX_X_tmp[3][2] = out_RX[3];
			//进入第1轮的下半轮,即第2个半轮
			state = Chaskey_RX_Middle_Rounds(search_round, 2, out_RX);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
{
	////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	if(state == 2)
	{
		updated = 2;
		return updated;
	}
	updated = 1;  //准备通知已经更新
}
#endif
		}
	}
	return updated;
}




u32 Chaskey_RX_Middle_Rounds(u16 search_round, u16 current_round, u32 *Xr)
{
	u32 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 xor_0 = 0, xor_1 = 0;
	u16 w_xor_0 = 0,w_xor_1 = 0;
	u16 w_xor = 0;
	u32 i = 0,j = 0;
	u32 A0_alpha_bloc[4] = {0};
	u32 A0_beta_bloc[4] = {0};
	u32 A0_gamma_bloc[4] = {0};
	u32 A0_gamma_temp = 0;
	u32 A0_AB_block[4] = {0};
	u32 A0_carry[4] ={0};
	u32 A0_carry_tmp[4] ={0};
	u32 A0_i0=0,A0_i1=0,A0_i2=0,A0_i3=0;
	u32 A0_j0=0,A0_j1=0,A0_j2=0,A0_j3=0;
	u32 A0_w1=0,A0_w2=0,A0_w3=0;
	u8 A0_w_xor_block[8] = {0};
	u32 A1_alpha_bloc[4] = {0};
	u32 A1_beta_bloc[4] = {0};
	u32 A1_gamma_bloc[4] = {0};
	u32 A1_gamma_temp = 0;
	u32 A1_AB_block[4] = {0};
	u32 A1_carry[4] ={0};
	u32 A1_carry_tmp[4] ={0};
	u32 A1_i0=0,A1_i1=0,A1_i2=0,A1_i3=0;
	u32 A1_j0=0,A1_j1=0,A1_j2=0,A1_j3=0;
	u32 A1_w1=0,A1_w2=0,A1_w3=0;
	u8 A1_w_xor_block[8] = {0};

	u32 out_X1_tmp=0, out_X3_tmp=0;
	u32 out_RX[4];

	u32 zeta_num_0 = 0,zeta_num_1 = 0;
	float w_cmp = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;
	if(search_round == current_round)
	{
		updated = Chaskey_RX_N_Rounds(search_round, Xr);
		return updated;
	}


	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[current_round -1]
				- RX_P_bestofR_w[search_round - current_round];
	w_cmp = w_cmp_rx - 2*Ek_min;

	xor_0 = (Xr[0] ^ Xr[1]) & ValueMax_Align ;
	xor_1 = (Xr[2] ^ Xr[3]) & ValueMax_Align ;
	w_xor_0 = HM_weight(xor_0);
	w_xor_1 = HM_weight(xor_1);
	w_xor = w_xor_0 + w_xor_1;
	if (w_xor >= w_cmp)
	{
		return 0;
	}

	if(current_round % 2 == 0)  //偶半轮，即某轮的下半轮，旋转参数为7/13/16
	{
		//rc=16;
		out_X1_tmp = chaskey_rol_left_7(Xr[1]);
		out_X3_tmp = chaskey_rol_left_13(Xr[3]);
	}
	else //奇半轮，即某轮的上半轮，旋转参数为5/8/16
	{
		//rc=16;
		out_X1_tmp = chaskey_rol_left_5(Xr[1]);
		out_X3_tmp = chaskey_rol_left_8(Xr[3]);
	}


	////////////////////////////对应32bit字长版本的切块预处理////////////////////////////
		for(i=0;i<nBytes;i++)  //nBytes
		{
			A0_alpha_bloc[i] = (Xr[0] >> (8*i)) & 0xFF; //8 bit
			A0_beta_bloc[i]  = (Xr[1] >> (8*i))  & 0xFF; //8 bit
			A0_AB_block[i] = ((A0_alpha_bloc[i] << 8) | A0_beta_bloc[i]);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A0_w_xor_block[j] = (u8)_mm_popcnt_u64((xor_0 >> (j*8)));
		}
		A0_carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A0_carry_tmp[j] = ((A0_alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((A0_beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}
	////////////////////////////对应32bit字长版本的切块预处理////////////////////////////
		for(i=0;i<nBytes;i++)  //nBytes
		{
			A1_alpha_bloc[i] = (Xr[2] >> (8*i)) & 0xFF; //8 bit
			A1_beta_bloc[i]  = (Xr[3] >> (8*i))  & 0xFF; //8 bit
			A1_AB_block[i] = ((A1_alpha_bloc[i] << 8) | A1_beta_bloc[i]);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A1_w_xor_block[j] = (u8)_mm_popcnt_u64((xor_1 >> (j*8)));
		}
		A1_carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A1_carry_tmp[j] = ((A1_alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((A1_beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}


/////////////////////第一个模加
		for(A0_i0 = cDDT_wt_min[A0_carry[0]][A0_AB_block[0]];
				A0_i0 <= cDDT_wt_max[A0_carry[0]][A0_AB_block[0]]; A0_i0++)
		{
			if(A0_i0 + A0_w_xor_block[1] +w_xor_1>= w_cmp ){break;}
		for(A0_j0=0; A0_j0 < cDDT_n[A0_carry[0]][A0_AB_block[0]][A0_i0]; A0_j0++)
		{
			A0_gamma_bloc[0] = cDDT_v[A0_carry[0]][A0_AB_block[0]][A0_i0][A0_j0];
			A0_carry[1] = A0_carry_tmp[1] + (A0_gamma_bloc[0] >> 7); // gamma MSB

			for(A0_i1 = cDDT_wt_min[A0_carry[1]][A0_AB_block[1]];
					A0_i1 <= cDDT_wt_max[A0_carry[1]][A0_AB_block[1]]; A0_i1++)
			{
				A0_w1 = A0_i0 + A0_i1;
				if(A0_w1 + A0_w_xor_block[2]  +w_xor_1>= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A0_j1=0; A0_j1 < cDDT_n[A0_carry[1]][A0_AB_block[1]][A0_i1]; A0_j1++)
			{
				A0_gamma_bloc[1] = cDDT_v[A0_carry[1]][A0_AB_block[1]][A0_i1][A0_j1];
				A0_carry[2] = A0_carry_tmp[2] + (A0_gamma_bloc[1] >> 7); // gamma MSB

			for(A0_i2 = cDDT_wt_min[A0_carry[2]][A0_AB_block[2]];
					A0_i2 <= cDDT_wt_max[A0_carry[2]][A0_AB_block[2]]; A0_i2++)
			{
				A0_w2 = A0_w1 + A0_i2;
				if(A0_w2 + A0_w_xor_block[3] +w_xor_1 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A0_j2=0; A0_j2 < cDDT_n[A0_carry[2]][A0_AB_block[2]][A0_i2]; A0_j2++)
			{
				A0_gamma_bloc[2] = cDDT_v[A0_carry[2]][A0_AB_block[2]][A0_i2][A0_j2];
				//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
				A0_carry[3] = A0_carry_tmp[3] + (A0_gamma_bloc[2] >> 7); // gamma MSB

			for(A0_i3 = cDDT_wt_min[A0_carry[3]][A0_AB_block[3]];
					A0_i3 <= MSB_cDDT_wt_max[A0_carry[3]][A0_AB_block[3]]; A0_i3++)
			{
				A0_w3 = A0_w2 + A0_i3;
				if(A0_w3  +w_xor_1>= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(A0_j3=0; A0_j3 < MSB_cDDT_n[A0_carry[3]][A0_AB_block[3]][A0_i3]; A0_j3++)
				{
					A0_gamma_bloc[3] = MSB_cDDT_v[A0_carry[3]][A0_AB_block[3]][A0_i3][A0_j3];
					A0_gamma_temp = (A0_gamma_bloc[3] << 24) | (A0_gamma_bloc[2] << 16)
							   | (A0_gamma_bloc[1] <<8) | A0_gamma_bloc[0];


/////////////////////第二个模加
		for(A1_i0 = cDDT_wt_min[A1_carry[0]][A1_AB_block[0]];
				A1_i0 <= cDDT_wt_max[A1_carry[0]][A1_AB_block[0]]; A1_i0++)
		{
			if(A0_w3 + A1_i0 + A1_w_xor_block[1] >= w_cmp ){break;}
		for(A1_j0=0; A1_j0 < cDDT_n[A1_carry[0]][A1_AB_block[0]][A1_i0]; A1_j0++)
		{
			A1_gamma_bloc[0] = cDDT_v[A1_carry[0]][A1_AB_block[0]][A1_i0][A1_j0];
			A1_carry[1] = A1_carry_tmp[1] + (A1_gamma_bloc[0] >> 7); // gamma MSB

			for(A1_i1 = cDDT_wt_min[A1_carry[1]][A1_AB_block[1]];
					A1_i1 <= cDDT_wt_max[A1_carry[1]][A1_AB_block[1]]; A1_i1++)
			{
				A1_w1 = A1_i0 + A1_i1;
				if(A0_w3 + A1_w1 + A1_w_xor_block[2] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j1=0; A1_j1 < cDDT_n[A1_carry[1]][A1_AB_block[1]][A1_i1]; A1_j1++)
			{
				A1_gamma_bloc[1] = cDDT_v[A1_carry[1]][A1_AB_block[1]][A1_i1][A1_j1];
				A1_carry[2] = A1_carry_tmp[2] + (A1_gamma_bloc[1] >> 7); // gamma MSB

			for(A1_i2 = cDDT_wt_min[A1_carry[2]][A1_AB_block[2]];
					A1_i2 <= cDDT_wt_max[A1_carry[2]][A1_AB_block[2]]; A1_i2++)
			{
				A1_w2 = A1_w1 + A1_i2;
				if(A0_w3 + A1_w2 + A1_w_xor_block[3] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j2=0; A1_j2 < cDDT_n[A1_carry[2]][A1_AB_block[2]][A1_i2]; A1_j2++)
			{
				A1_gamma_bloc[2] = cDDT_v[A1_carry[2]][A1_AB_block[2]][A1_i2][A1_j2];
				A1_carry[3] = A1_carry_tmp[3] + (A1_gamma_bloc[2] >> 7); // gamma MSB
				for(A1_i3 = cDDT_wt_min[A1_carry[3]][A1_AB_block[3]];
						A1_i3 <= MSB_cDDT_wt_max[A1_carry[3]][A1_AB_block[3]];A1_i3++)
				{
					A1_w3 = A1_w2 + A1_i3;
					if(A0_w3 + A1_w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.

//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
			for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
				for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
				{
					w_cmp = w_cmp_rx - RX_zeta_w[zeta_num_0] - RX_zeta_w[zeta_num_1];
					if(A0_w3 + A1_w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.

					for(A1_j3=0; A1_j3 < MSB_cDDT_n[A1_carry[3]][A1_AB_block[3]][A1_i3]; A1_j3++)
					{
						A1_gamma_bloc[3] = MSB_cDDT_v[A1_carry[3]][A1_AB_block[3]][A1_i3][A1_j3];
						A1_gamma_temp = (A1_gamma_bloc[3] << 24) | (A1_gamma_bloc[2] << 16)
								   | (A1_gamma_bloc[1] <<8) | A1_gamma_bloc[0];


			RX_Chaskey_ADD_P_xor_w_tmp[0][current_round] = A0_w3;
			RX_Chaskey_ADD_P_xor_w_tmp[1][current_round] = A1_w3;
			P_w[current_round] = A0_w3 + A1_w3;
			Chaskey_Pr_zeta_tmp[0][current_round] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
			Chaskey_Pr_zeta_tmp[1][current_round] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
			Chaskey_zeta_tmp[0][current_round] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
			Chaskey_zeta_tmp[1][current_round] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1
			RX_Chaskey_P_hr_w_tmp[current_round] = P_w[current_round]
											    + Chaskey_Pr_zeta_tmp[0][current_round]
												+ Chaskey_Pr_zeta_tmp[1][current_round];
			RX_P_sumofR_w[current_round] = RX_Chaskey_P_hr_w_tmp[current_round]
											+ RX_P_sumofR_w[current_round -1];


			out_RX[2] =  A0_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[1] = out_X1_tmp ^ out_RX[2];
			out_RX[2] = chaskey_rol_left_16(out_RX[2]);
			out_RX[0] = A1_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_1]];
			out_RX[3] = out_RX[0] ^ out_X3_tmp;

			RX_X_tmp[0][current_round +1] = out_RX[0];
			RX_X_tmp[1][current_round +1] = out_RX[1];
			RX_X_tmp[2][current_round +1] = out_RX[2];
			RX_X_tmp[3][current_round +1] = out_RX[3];

			//进入下一个半轮,
			state = Chaskey_RX_Middle_Rounds(search_round, current_round+1, out_RX);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
{
	////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	if(state == 2)
	{
		updated = 2;
		return updated;
	}
	updated = 1;  //准备通知已经更新
}
#endif
		}
	}}}}}}}}}}}}}}}}

	return updated;
}


u32 Chaskey_RX_N_Rounds(u16 search_round,u32 *Xr)
{
	u32 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 xor_0 = 0, xor_1 = 0;
	u16 w_xor_0 = 0,w_xor_1 = 0;
	u16 w_xor = 0;
	u32 i = 0,j = 0;
	u32 A0_alpha_bloc[4] = {0};
	u32 A0_beta_bloc[4] = {0};
	u32 A0_gamma_bloc[4] = {0};
	u32 A0_gamma_temp = 0;
	u32 A0_AB_block[4] = {0};
	u32 A0_carry[4] ={0};
	u32 A0_carry_tmp[4] ={0};
	u32 A0_i0=0,A0_i1=0,A0_i2=0,A0_i3=0;
	u32 A0_j0=0,A0_j1=0,A0_j2=0,A0_j3=0;
	u32 A0_w1=0,A0_w2=0,A0_w3=0;
	u8 A0_w_xor_block[8] = {0};
	u32 A1_alpha_bloc[4] = {0};
	u32 A1_beta_bloc[4] = {0};
	u32 A1_gamma_bloc[4] = {0};
	u32 A1_gamma_temp = 0;
	u32 A1_AB_block[4] = {0};
	u32 A1_carry[4] ={0};
	u32 A1_carry_tmp[4] ={0};
	u32 A1_i0=0,A1_i1=0,A1_i2=0,A1_i3=0;
	u32 A1_j0=0,A1_j1=0,A1_j2=0,A1_j3=0;
	u32 A1_w1=0,A1_w2=0,A1_w3=0;
	u8 A1_w_xor_block[8] = {0};

	u32 out_X1_tmp=0, out_X3_tmp=0;
	u32 out_RX[4];

	u32 zeta_num_0 = 0,zeta_num_1 = 0;
	float w_cmp = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;
	float Bn_cmp = 0.0;


	updata_cur_num = Num_Bn_Update;

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[search_round -1];
	w_cmp = w_cmp_rx - 2*Ek_min;
	xor_0 = (Xr[0] ^ Xr[1]) & ValueMax_Align ;
	xor_1 = (Xr[2] ^ Xr[3]) & ValueMax_Align ;
	w_xor_0 = HM_weight(xor_0);
	w_xor_1 = HM_weight(xor_1);
	w_xor = w_xor_0 + w_xor_1;
	if (w_xor >= w_cmp)
	{
		return 0;
	}

	if(search_round % 2 == 0)  //偶半轮，即某轮的下半轮，旋转参数为7/13/16
	{
		//rc=16;
		out_X1_tmp = chaskey_rol_left_7(Xr[1]);
		out_X3_tmp = chaskey_rol_left_13(Xr[3]);
	}
	else //奇半轮，即某轮的上半轮，旋转参数为5/8/16
	{
		//rc=16;
		out_X1_tmp = chaskey_rol_left_5(Xr[1]);
		out_X3_tmp = chaskey_rol_left_8(Xr[3]);
	}
///////////

	////////////////////////////对应32bit字长版本的切块预处理////////////////////////////
		for(i=0;i<nBytes;i++)  //nBytes
		{
			A0_alpha_bloc[i] = (Xr[0] >> (8*i)) & 0xFF; //8 bit
			A0_beta_bloc[i]  = (Xr[1] >> (8*i))  & 0xFF; //8 bit
			A0_AB_block[i] = ((A0_alpha_bloc[i] << 8) | A0_beta_bloc[i]);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A0_w_xor_block[j] = (u8)_mm_popcnt_u64((xor_0 >> (j*8)));
		}
		A0_carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A0_carry_tmp[j] = ((A0_alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((A0_beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}
	////////////////////////////对应32bit字长版本的切块预处理////////////////////////////
		for(i=0;i<nBytes;i++)  //nBytes
		{
			A1_alpha_bloc[i] = (Xr[2] >> (8*i)) & 0xFF; //8 bit
			A1_beta_bloc[i]  = (Xr[3] >> (8*i))  & 0xFF; //8 bit
			A1_AB_block[i] = ((A1_alpha_bloc[i] << 8) | A1_beta_bloc[i]);
		}
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A1_w_xor_block[j] = (u8)_mm_popcnt_u64((xor_1 >> (j*8)));
		}
		A1_carry[0] = 0;
		for(j=1;j<nBytes;j++)  //nBytes
		{
			A1_carry_tmp[j] = ((A1_alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((A1_beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}
///////////


/////////////////////第一个模加
		for(A0_i0 = cDDT_wt_min[A0_carry[0]][A0_AB_block[0]];
				A0_i0 <= cDDT_wt_max[A0_carry[0]][A0_AB_block[0]]; A0_i0++)
		{
			if(A0_i0 + A0_w_xor_block[1] + w_xor_1 >= w_cmp ){break;}
		for(A0_j0=0; A0_j0 < cDDT_n[A0_carry[0]][A0_AB_block[0]][A0_i0]; A0_j0++)
		{
			A0_gamma_bloc[0] = cDDT_v[A0_carry[0]][A0_AB_block[0]][A0_i0][A0_j0];
			A0_carry[1] = A0_carry_tmp[1] + (A0_gamma_bloc[0] >> 7); // gamma MSB

			for(A0_i1 = cDDT_wt_min[A0_carry[1]][A0_AB_block[1]];
					A0_i1 <= cDDT_wt_max[A0_carry[1]][A0_AB_block[1]]; A0_i1++)
			{
				A0_w1 = A0_i0 + A0_i1;
				if(A0_w1 + A0_w_xor_block[2] +w_xor_1 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A0_j1=0; A0_j1 < cDDT_n[A0_carry[1]][A0_AB_block[1]][A0_i1]; A0_j1++)
			{
				A0_gamma_bloc[1] = cDDT_v[A0_carry[1]][A0_AB_block[1]][A0_i1][A0_j1];
				A0_carry[2] = A0_carry_tmp[2] + (A0_gamma_bloc[1] >> 7); // gamma MSB
			for(A0_i2 = cDDT_wt_min[A0_carry[2]][A0_AB_block[2]];
					A0_i2 <= cDDT_wt_max[A0_carry[2]][A0_AB_block[2]]; A0_i2++)
			{
				A0_w2 = A0_w1 + A0_i2;
				if(A0_w2 + A0_w_xor_block[3] +w_xor_1 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A0_j2=0; A0_j2 < cDDT_n[A0_carry[2]][A0_AB_block[2]][A0_i2]; A0_j2++)
			{
				A0_gamma_bloc[2] = cDDT_v[A0_carry[2]][A0_AB_block[2]][A0_i2][A0_j2];
				//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
				A0_carry[3] = A0_carry_tmp[3] + (A0_gamma_bloc[2] >> 7); // gamma MSB

				for(A0_i3 = cDDT_wt_min[A0_carry[3]][A0_AB_block[3]];
						A0_i3 <= MSB_cDDT_wt_max[A0_carry[3]][A0_AB_block[3]]; A0_i3++)
				{
					A0_w3 = A0_w2 + A0_i3;
					if(A0_w3 +w_xor_1 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A0_j3=0; A0_j3 < MSB_cDDT_n[A0_carry[3]][A0_AB_block[3]][A0_i3]; A0_j3++)
					{
						A0_gamma_bloc[3] = MSB_cDDT_v[A0_carry[3]][A0_AB_block[3]][A0_i3][A0_j3];
						A0_gamma_temp = (A0_gamma_bloc[3] << 24) | (A0_gamma_bloc[2] << 16)
								   | (A0_gamma_bloc[1] <<8) | A0_gamma_bloc[0];

/////////////////////第二个模加
		for(A1_i0 = cDDT_wt_min[A1_carry[0]][A1_AB_block[0]];
				A1_i0 <= cDDT_wt_max[A1_carry[0]][A1_AB_block[0]]; A1_i0++)
		{
			if(A0_w3 + A1_i0 + A1_w_xor_block[1] >= w_cmp ){break;}
		for(A1_j0=0; A1_j0 < cDDT_n[A1_carry[0]][A1_AB_block[0]][A1_i0]; A1_j0++)
		{
			A1_gamma_bloc[0] = cDDT_v[A1_carry[0]][A1_AB_block[0]][A1_i0][A1_j0];
			A1_carry[1] = A1_carry_tmp[1] + (A1_gamma_bloc[0] >> 7); // gamma MSB

		for(A1_i1 = cDDT_wt_min[A1_carry[1]][A1_AB_block[1]];
					A1_i1 <= cDDT_wt_max[A1_carry[1]][A1_AB_block[1]]; A1_i1++)
		{
			A1_w1 = A1_i0 + A1_i1;
			if(A0_w3 + A1_w1 + A1_w_xor_block[2] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(A1_j1=0; A1_j1 < cDDT_n[A1_carry[1]][A1_AB_block[1]][A1_i1]; A1_j1++)
		{
			A1_gamma_bloc[1] = cDDT_v[A1_carry[1]][A1_AB_block[1]][A1_i1][A1_j1];
			A1_carry[2] = A1_carry_tmp[2] + (A1_gamma_bloc[1] >> 7); // gamma MSB

			for(A1_i2 = cDDT_wt_min[A1_carry[2]][A1_AB_block[2]];
					A1_i2 <= cDDT_wt_max[A1_carry[2]][A1_AB_block[2]]; A1_i2++)
			{
				A1_w2 = A1_w1 + A1_i2;
				if(A0_w3 + A1_w2 + A1_w_xor_block[3] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j2=0; A1_j2 < cDDT_n[A1_carry[2]][A1_AB_block[2]][A1_i2]; A1_j2++)
			{
				A1_gamma_bloc[2] = cDDT_v[A1_carry[2]][A1_AB_block[2]][A1_i2][A1_j2];
				A1_carry[3] = A1_carry_tmp[3] + (A1_gamma_bloc[2] >> 7); // gamma MSB
				for(A1_i3 = cDDT_wt_min[A1_carry[3]][A1_AB_block[3]];
						A1_i3 <= MSB_cDDT_wt_max[A1_carry[3]][A1_AB_block[3]];A1_i3++)
				{
					A1_w3 = A1_w2 + A1_i3;
					if(A0_w3 + A1_w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A1_j3=0; A1_j3 < MSB_cDDT_n[A1_carry[3]][A1_AB_block[3]][A1_i3]; A1_j3++)
					{
						A1_gamma_bloc[3] = MSB_cDDT_v[A1_carry[3]][A1_AB_block[3]][A1_i3][A1_j3];
						A1_gamma_temp = (A1_gamma_bloc[3] << 24) | (A1_gamma_bloc[2] << 16)
								   | (A1_gamma_bloc[1] <<8) | A1_gamma_bloc[0];

//%%%%%%%%基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减%%%%%%%%%%%%%%%%%%%%%%
//	for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
//		for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
		{
	//		zeta_num_0=0; zeta_num_1=0;  //启发式，只取EK_min
	//		w_cmp = w_cmp_rx - RX_zeta_w[zeta_num_0] - RX_zeta_w[zeta_num_1];
	//		if(A0_w3 + A1_w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.

			RX_Chaskey_ADD_P_xor_w_tmp[0][search_round] = A0_w3;
			RX_Chaskey_ADD_P_xor_w_tmp[1][search_round] = A1_w3;
			P_w[search_round] = A0_w3 + A1_w3;
			Chaskey_Pr_zeta_tmp[0][search_round] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
			Chaskey_Pr_zeta_tmp[1][search_round] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
			Chaskey_zeta_tmp[0][search_round] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
			Chaskey_zeta_tmp[1][search_round] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1


			RX_Chaskey_P_hr_w_tmp[search_round] = P_w[search_round]
											    + Chaskey_Pr_zeta_tmp[0][search_round]
												+ Chaskey_Pr_zeta_tmp[1][search_round];
//			RX_P_sumofR_w[search_round] = RX_Chaskey_P_hr_w_tmp[search_round]
//											+ RX_P_sumofR_w[search_round -1];


			out_RX[2] =  A0_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[1] = out_X1_tmp ^ out_RX[2];
			out_RX[2] = chaskey_rol_left_16(out_RX[2]);
			out_RX[0] = A1_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_1]];
			out_RX[3] = out_RX[0] ^ out_X3_tmp;
			RX_X[0][search_round +1] = out_RX[0];
			RX_X[1][search_round +1] = out_RX[1];
			RX_X[2][search_round +1] = out_RX[2];
			RX_X[3][search_round +1] = out_RX[3];

		//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
		for(i=1; i <= search_round; i++ )
		{
			//记录更新后的每个半轮的RX差分
			RX_X[0][i] = RX_X_tmp[0][i];
			RX_X[1][i] = RX_X_tmp[1][i];
			RX_X[2][i] = RX_X_tmp[2][i];
			RX_X[3][i] = RX_X_tmp[3][i];
			//记录更新后的每个半轮的zeta
			Chaskey_Zeta[0][i] = Chaskey_zeta_tmp[0][i];
			Chaskey_Zeta[1][i] = Chaskey_zeta_tmp[1][i];
			//记录更新后的每个半轮的zeta的概率
			Chaskey_Pr_zeta[0][i] = Chaskey_Pr_zeta_tmp[0][i];
			Chaskey_Pr_zeta[1][i] = Chaskey_Pr_zeta_tmp[1][i];
			//记录更新后的每个半轮的模加的XOR差分概率
			RX_Chaskey_ADD_P_xor_w[0][i] = RX_Chaskey_ADD_P_xor_w_tmp[0][i];
			RX_Chaskey_ADD_P_xor_w[1][i] = RX_Chaskey_ADD_P_xor_w_tmp[1][i];
			//记录更新后的每个半轮的RX总的概率重量，包括XOR和zeta
			RX_Chaskey_P_hr_w[i] = RX_Chaskey_P_hr_w_tmp[i];
		}
//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
					Num_Bn_Update++;  //更新次数记录
					Bn_cmp = 0;
					for(i=1; i <= search_round; i++ )
					{
						Bn_cmp += RX_Chaskey_P_hr_w[i];
					}
					RX_Expected_w = Bn_cmp;
					w_cmp = RX_Expected_w - 2*Ek_min - RX_P_sumofR_w[search_round -1];
					state = 1;
					updated = 1;


		//%%%%%打印更新次数等信息
		time_Round = clock();
		run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
		printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
		printf("Time: %.2f seconds. ==>>\n", run_time);

		}
			}}}}}}}}}}}}}}}}

	return updated;
}





u16 print_Chaskey_RX_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("-r----X0--------X1-------X2-------X3------zeta_0---zeta_1----P_zeta_0<->P_zeta_3---P_xor_0---P_xor_1---Pw_hr------Pw_r\n");
	for(r_num=1;r_num <= search_round; )
	{
		printf("%02d  %08llx %08llx %08llx %08llx  %08llx %08llx   %f   %f   %f  %f  %f     -\n",
				r_num,
				RX_X[0][r_num],RX_X[1][r_num],RX_X[2][r_num],RX_X[3][r_num],
				Chaskey_Zeta[0][r_num],Chaskey_Zeta[1][r_num],
				Chaskey_Pr_zeta[0][r_num],Chaskey_Pr_zeta[1][r_num],
				RX_Chaskey_ADD_P_xor_w[0][r_num],RX_Chaskey_ADD_P_xor_w[1][r_num],
				RX_Chaskey_P_hr_w[r_num]); //,gama_cnt[r_num]);

//printf("hr----XV0--------XV1---------XV2---------XV3------zeta_0------zeta_1------P_zeta_0---P_zeta_0----P_xor_0---P_xor_1----Pw_hr \n");
		printf("%02d  %08llx %08llx %08llx %08llx  %08llx %08llx   %f   %f   %f  %f  %f  %f \n",
				r_num+1,
				RX_X[0][r_num+1],RX_X[1][r_num+1],RX_X[2][r_num+1],RX_X[3][r_num+1],
				Chaskey_Zeta[0][r_num+1],Chaskey_Zeta[1][r_num+1],
				Chaskey_Pr_zeta[0][r_num+1],Chaskey_Pr_zeta[1][r_num+1],
				RX_Chaskey_ADD_P_xor_w[0][r_num+1],RX_Chaskey_ADD_P_xor_w[1][r_num+1],
				RX_Chaskey_P_hr_w[r_num+1],
			    RX_Chaskey_P_hr_w[r_num] + RX_Chaskey_P_hr_w[r_num+1] ); //

		r_num = r_num + 2;
	}
	printf("out %08llx %08llx %08llx %08llx   \n",
			RX_X[0][search_round+1],RX_X[1][search_round+1],RX_X[2][search_round+1],RX_X[3][search_round+1]);

	printf("------------------ \n");
	printf("%d-half-Round Chaskey Total Weight: -%f \n",search_round, RX_Expected_w);

	return 1;
}




//	ADD_RX_check(1);
////////////////////////////////以下为Chaskey的RX差分特征验证，实验代码///////////
u32 ADD_RX_check(u32 rk)
{
	u32 X_V = 0, X_U=0;
	u32 X[4] = {0}, Y[4]={0};
	u32 x0=0, x1 =0;
	u32 y0=0, y1=0;
	u64 Nx = 0xFF;
	u64 Ny = 0xFF;
	u64 N =0;  //计数
	u64 t=0;
	float pr = 0;
	u32 tx0=0, tx1 =0;




	t = (Nx+1) * (Ny+1);

	for(X_V=0xFFFF; X_V< 0xFFFF+ Nx; X_V++)  //0xFFF1FFFF  =  Nx
	{
		x0 = X_V;  //(X_V  << 1) ^ 0x1;
		y0 = rol_left32(x0,1) ^ 0x00000001;
		for(X_U=0xFFF1FFFF; X_U<= 0xFFF1FFFF + Ny; X_U++)  //0xFFFFFFFF
		{
			x1 = X_U;  //(X_U  << 1) ^ 0x1;
			y1 = rol_left32(x1,1);

			//printf("---x0: %x   x1: %x -----", x0,x1);
			//printf("---y0: %x   y1: %x -----", y0,y1);
			//N += ADD_RX_delta_check(x0,x1,y0,y1,  1);

			for(tx0 = 0x0A0A; tx0 < 0x0A0F; tx0++)
				for(tx1 = 0x0A0A; tx1 < 0x0A0F; tx1++)
				{
					X[0] = tx0; Y[0] = rol_left32(tx0,1);
					X[1] = tx1; Y[1] = rol_left32(tx1,1);

			//X[0]=0x10;X[1]=0x0F;
			X[2]=x0;X[3]=x1;
			//Y[0]=0x20;Y[1]=0x1E;
			Y[2]=y0;Y[3]=y1;
			Chaskey_1_Rounds(X,Y);
			if(Y[0] == rol_left32(X[0],1))
				if(Y[1] == rol_left32(X[1],1))
					if(Y[2] == rol_left32(X[2],1))
						if(Y[3] == rol_left32(X[3],1))
			{
				N++;
			}
		}
		}
	}

	pr = (float)N/t;
	printf("t = %d------ N =  %d ----- Pr = %f------\n", t, N, pr);
	printf("X[0] = %x--X[1] = %x--X[2] = %x--X[3] = %x\n", X[0],X[1],X[2],X[3]);
	printf("Y[0] = %x--Y[1] = %x--Y[2] = %x--Y[3] = %x\n", Y[0],Y[1],Y[2],Y[3]);
return N;
}


u32 ADD_RX_delta_check(u32 X0, u32 X1, u32 Y0, u32 Y1, u32 rk)
{
	u32 cnt=0;
	u32 tmp_X = 0, tmp_Y = 0;
	u32 delta = 0xFF;

	tmp_X = (X0 + X1) & 0xFFFFFFFF;
	tmp_Y = (Y0 + Y1) & 0xFFFFFFFF;

	delta = rol_left32(tmp_X,rk) ^ tmp_Y;
	printf("tmp_X: %x---tmp_Y: %x----delta: %x------\n",tmp_X, tmp_Y,delta);

	//zeta =v||u; v=0, u=1,则 高n-1比特相加无进位，低1比特加法有进位
	if((delta ^ 0x00000000) == 0)  //判断RX差分是否为0
	{
		cnt =1;
	}

	return cnt;
}



