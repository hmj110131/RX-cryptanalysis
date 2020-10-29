/*
 * RX_Siphash.c
 *
 *  Created on: 2020年7月31日
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
#include "XOR_offset_table.h"
#include <nmmintrin.h>
#include "Siphash.h"
#include "RX_Siphash.h"





////////////////////////////////搜索Siphash的最优RX特征的代码//////////////////////////////
volatile float RX_Siphash_P_hr_w[25] = {0};; //每半轮中，对应RX差分概率重量
volatile float RX_Siphash_P_hr_w_tmp[25] = {0};; //每轮对应RX差分概率重量,临时的
volatile float RX_Siphash_ADD_P_xor_w[2][25] = {0};; //每半轮中2个模块加法的对应RX差分概率重量
volatile float RX_Siphash_ADD_P_xor_w_tmp[2][25] = {0};; //每轮中，4个模块加法的对应RX差分概率重量
volatile u64 Siphash_R1_ADD_i_RX_abc[2][3] = {0};; //2个模块加法的对应a/b/c
volatile u64 Siphash_RX_X[4][25] = {0};; //每半轮中，4个分支的输入RX差分
volatile u64 Siphash_RX_X_tmp[4][25] = {0};; //每轮中，4个分支的输入RX差分
volatile u64 Siphash_Zeta[2][25] = {0};; //每半轮的输zeta
volatile u64 Siphash_zeta_tmp[2][25] = {0};; //每半轮的临时zeta
volatile float Siphash_Pr_zeta[2][25] = {0};; //每轮RX差分概率
volatile float Siphash_Pr_zeta_tmp[2][25] = {0};; //每轮RX差分概率




////////////////////////////////搜索Siphash的最优RX特征的代码//////////////////////////////
u64 Siphash_RX_trail_search_entry(u16 search_round)
{
	u64 best = 0;
	u16 i = 0;
	FILE* Siphash_RX_Bn_wt;


	printf("--Note: For Siphash, the 'round' here means half-round of SipRound.--\n");
	if(sc_blocksize != 256)
	{
		printf("The block size of Siphash should be 256-bits. \n");
		return 0;
	}

	// For Siphash, word length of  modular additon is: 64 bits.
	blocksize_len = 64;
	nBytes = 8;
//	word_len = 64;  //模加的字长
	Bit_Align = 0xFFFFFFFFFFFFFFFF;
	ValueMax_Align = 0x7FFFFFFFFFFFFFFF;
	V_MSB = 0x8000000000000000;
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
	Siphash_RX_Bn_wt = fopen ("../tmp/Siphash_RX_Bn_wt.xlsx", "a+"); //  "w+"); //
	for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
	{
		fscanf(Siphash_RX_Bn_wt,"%f", &RX_P_bestofR_w[i]);  //read.
	}
	//-------记录时间-----
		time_xor_talbe = clock();
		run_time =  (time_xor_talbe - time_start) / CLOCKS_PER_SEC;
		printf("Time of Preprocessing Siphash info: %.2f seconds.  \n", run_time);
		printf("|--------------------------------------------------|\n");
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
		printf("Constructing cDDT Lookup Tables ... \n");
		//构造cDDT
		ARX_carry_DDTm_construct();  // FOR cDDT
		//fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
		time_ARX_DDT = clock();
		run_time =  (time_ARX_DDT - time_xor_talbe) / CLOCKS_PER_SEC;
		printf("Time of Construct Siphash cDDT: %.2f seconds.  \n", run_time);
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
			best = Siphash_RX_ADD_wt_inc(search_round); //效率高

			time_Round = clock();
			run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
			printf("Time: %.2f seconds.\n", run_time);

		}while( Num_Bn_Update == 0 );   //end conditions. // 概率重量逼近次数大于1，且返回值表示搜索到最优路径

		fprintf(Siphash_RX_Bn_wt,"%f \r\n",RX_Expected_w);  // write.RX_P_bestofR_w[search_round]
		printf("**************************************************************************|\n");
		print_Siphash_RX_resoult(search_round);  //打印

		time_finish = clock();
		run_time = (float)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
		printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours == %.2f days. \n", run_time,run_time/60.0,run_time/3600.0,run_time/86400.0 );
		printf("Auto-search %d-Round Siphash optimal RX-differential trails END! \n", search_round);
		printf("|************************************************************************|\n");
		fclose(Siphash_RX_Bn_wt);

	return 0;
}


u64 Siphash_RX_ADD_wt_inc(u16 search_round)
{
	u64 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 wt_round_1 = 0;
	u16 wt_Block_0 = 0,wt_Block_1 = 1; //Siphash的RX差分特征，只需要构造首轮前半轮的两个模架
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
						RX_Siphash_ADD_P_xor_w_tmp[0][1] = wt_Block_0;
						RX_Siphash_ADD_P_xor_w_tmp[1][1] = wt_Block_1;

						//根据第1轮的前半轮中两个模加的，每个模加的RX差分概率重量构造对应的输入输出XOR差分
						best = Siphash_RX_ADD_i_tuples(search_round,0);

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
u64 Siphash_RX_ADD_i_tuples(u16 search_round,u16 ADD_part)
{
	u64 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 thr_d = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;

	u64 zeta_num = 0;
	u64 updata_cur_num = 0;

	u64 M0 = 0;
	u64 N0 = blocksize_len-1;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;

	updata_cur_num = Num_Bn_Update;

	///////////////////round 1 entry/////////////////////
	///////////2个模加的输入输出XOR差分构造完毕，则调到第1轮的下半轮//////////////
	if(ADD_part >= 2)
	{
		best = Siphash_RX_comb_RX(search_round);
		return best;
	}


	M0 = RX_Siphash_ADD_P_xor_w_tmp[ADD_part][1];
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

			Siphash_R1_ADD_i_RX_abc[ADD_part][0] = Input_alpha;
			Siphash_R1_ADD_i_RX_abc[ADD_part][1] = Input_beta;
			Siphash_R1_ADD_i_RX_abc[ADD_part][2] = Input_gamma;

			//直接跳转到第1轮各个部分加法的内部构造
			best = Siphash_RX_ADD_i_tuples(search_round,ADD_part+1);

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


		best = RX_Siphash_input_MSB(search_round,
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
						best = RX_Siphash_input_MSB(search_round,
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

u64 RX_Siphash_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 ADD_part )
{
	u64 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
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

			state = RX_Siphash_input_Middle(search_round,
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

		state = RX_Siphash_input_Middle(search_round,
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


u64 RX_Siphash_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi,u16 ADD_part)
{
	u64 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	//u16 i = 0;
	u64 bit_i = 1;


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
		state = RX_Siphash_input_Middle(search_round,
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
	state = RX_Siphash_input_Last(search_round,
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

u64 RX_Siphash_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi,u16 ADD_part)
{
	u64 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	u64 bit_i = 1;

	u64 zeta_num = 0;
	u64 updata_cur_num = 0;

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

	Siphash_R1_ADD_i_RX_abc[ADD_part][0] = Input_alpha;
	Siphash_R1_ADD_i_RX_abc[ADD_part][1] = Input_beta;
	Siphash_R1_ADD_i_RX_abc[ADD_part][2] = Input_gamma;

	//直接跳转到第1轮各个部分加法的内部构造
	state = Siphash_RX_ADD_i_tuples(search_round,ADD_part+1);

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
u64 Siphash_RX_comb_RX(u16 search_round)
{
	u64 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新

	u64 out_X[4];
	u64 out_RX[4];

	u64 zeta_num_0 = 0,zeta_num_1 = 0;
	float w_cmp_rx = 0;
	u64 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;

	w_cmp_rx = RX_Expected_w - P_w[1] -RX_P_bestofR_w[search_round -1];


	//首轮的输入RX差分
	Siphash_RX_X_tmp[0][1] = Siphash_R1_ADD_i_RX_abc[0][0]; //x0 =alpha
	Siphash_RX_X_tmp[1][1] = Siphash_R1_ADD_i_RX_abc[0][1]; //x1= beta
	Siphash_RX_X_tmp[2][1] = Siphash_R1_ADD_i_RX_abc[1][0]; //x2 =alpha
	Siphash_RX_X_tmp[3][1] = Siphash_R1_ADD_i_RX_abc[1][1]; //x3= beta

	//输出XOR差分，注意Siphash的x的排序，中间轮最左边为X1
	out_X[1] = Siphash_R1_ADD_i_RX_abc[0][2] ^siphash_rol_left_13(Siphash_RX_X_tmp[1][1]);
	out_X[2] = siphash_rol_left_32(Siphash_R1_ADD_i_RX_abc[0][2]);
	out_X[0] = Siphash_R1_ADD_i_RX_abc[1][2]; // OUT_X0 = 第二个莫加的gamma
	out_X[3] = out_X[0] ^ siphash_rol_left_16(Siphash_RX_X_tmp[3][1]);


	//用zeta补偿XOR差分，得到RX差分
	for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
	{
		for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
		{
			if((RX_zeta_w[zeta_num_0] + RX_zeta_w[zeta_num_1]) >= w_cmp_rx)
			{
				return 0;
			}

			Siphash_Pr_zeta_tmp[0][1] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
			Siphash_Pr_zeta_tmp[1][1] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
			Siphash_zeta_tmp[0][1] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
			Siphash_zeta_tmp[1][1] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1
			RX_Siphash_P_hr_w_tmp[1] = P_w[1] + Siphash_Pr_zeta_tmp[0][1]
											+ Siphash_Pr_zeta_tmp[1][1];
			RX_P_sumofR_w[1] = RX_Siphash_P_hr_w_tmp[1];


			out_RX[1] = out_X[1] ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[2] = out_X[2] ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[0] = out_X[0] ^ RX_zeta[RX_zeta_i[zeta_num_1]];
			out_RX[3] = out_X[3] ^ RX_zeta[RX_zeta_i[zeta_num_1]];

			Siphash_RX_X_tmp[0][2] = out_RX[0];
			Siphash_RX_X_tmp[1][2] = out_RX[1];
			Siphash_RX_X_tmp[2][2] = out_RX[2];
			Siphash_RX_X_tmp[3][2] = out_RX[3];
			//进入第1轮的下半轮,即第2个半轮
			state = Siphash_RX_Middle_Rounds(search_round, 2, out_RX);
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

u64 Siphash_RX_Middle_Rounds(u16 search_round, u16 current_round, u64 *Xr)
{
	u64 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 xor_0 = 0, xor_1 = 0;
	u64 w_xor_0 = 0,w_xor_1 = 0;
	u64 w_xor = 0;
	u64 i = 0,j = 0;
	u64 A0_alpha_bloc[8] = {0};
	u64 A0_beta_bloc[8] = {0};
	u64 A0_gamma_bloc[8] = {0};
	u64 A0_gamma_temp = 0;
	u64 A0_AB_block[8] = {0};
	u64 A0_carry[8] ={0};
	u64 A0_carry_tmp[8] ={0};
	u64 A0_i0=0,A0_i1=0,A0_i2=0,A0_i3=0,A0_i4=0,A0_i5=0,A0_i6=0,A0_i7=0;
	u64 A0_j0=0,A0_j1=0,A0_j2=0,A0_j3=0,A0_j4=0,A0_j5=0,A0_j6=0,A0_j7=0;
	u32 A0_w1=0,A0_w2=0,A0_w3=0,A0_w4=0,A0_w5=0,A0_w6=0,A0_w7=0;
	u8 A0_w_xor_block[8] = {0};
	u64 A1_alpha_bloc[8] = {0};
	u64 A1_beta_bloc[8] = {0};
	u64 A1_gamma_bloc[8] = {0};
	u64 A1_gamma_temp = 0;
	u64 A1_AB_block[8] = {0};
	u64 A1_carry[8] ={0};
	u64 A1_carry_tmp[8] ={0};
	u64 A1_i0=0,A1_i1=0,A1_i2=0,A1_i3=0,A1_i4=0,A1_i5=0,A1_i6=0,A1_i7=0;
	u64 A1_j0=0,A1_j1=0,A1_j2=0,A1_j3=0,A1_j4=0,A1_j5=0,A1_j6=0,A1_j7=0;
	u64 A1_w1=0,A1_w2=0,A1_w3=0,A1_w4=0,A1_w5=0,A1_w6=0,A1_w7=0;
	u8 A1_w_xor_block[8] = {0};

	u64 out_X1_tmp=0, out_X3_tmp=0;
	u64 out_RX[4];

	u64 zeta_num_0 = 0,zeta_num_1 = 0;
	float w_cmp = 0;
	float w_cmp_rx = 0;
	u64 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;
	if(search_round == current_round)
	{
		updated = Siphash_RX_N_Rounds(search_round, Xr);
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
		out_X1_tmp = siphash_rol_left_17(Xr[1]);
		out_X3_tmp = siphash_rol_left_21(Xr[3]);
	}
	else //奇半轮，即某轮的上半轮，旋转参数为5/8/16
	{
		//rc=16;
		out_X1_tmp = siphash_rol_left_13(Xr[1]);
		out_X3_tmp = siphash_rol_left_16(Xr[3]);
	}


	////////////////////////////对应64bit字长版本的切块预处理////////////////////////////
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
	////////////////////////////对应64bit字长版本的切块预处理////////////////////////////
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
				if(A0_w3 + A0_w_xor_block[4] +w_xor_1>= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A0_j3=0; A0_j3 < MSB_cDDT_n[A0_carry[3]][A0_AB_block[3]][A0_i3]; A0_j3++)
			{
				A0_gamma_bloc[3] = cDDT_v[A0_carry[3]][A0_AB_block[3]][A0_i3][A0_j3];
				A0_carry[4] = A0_carry_tmp[4] + (A0_gamma_bloc[3] >> 7); // gamma MSB
				for(A0_i4 = cDDT_wt_min[A0_carry[4]][A0_AB_block[4]];
						A0_i4 <= cDDT_wt_min[A0_carry[4]][A0_AB_block[4]]; A0_i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
				{
					A0_w4 = A0_w3 + A0_i4;
					if(A0_w4 + A0_w_xor_block[5] +w_xor_1> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(A0_j4=0; A0_j4 < cDDT_n[A0_carry[4]][A0_AB_block[4]][A0_i4]; A0_j4++)
				{
					A0_gamma_bloc[4] = cDDT_v[A0_carry[4]][A0_AB_block[4]][A0_i4][A0_j4];
					A0_carry[5] = A0_carry_tmp[5] + (A0_gamma_bloc[4] >> 7); // gamma MSB
				for(A0_i5 = cDDT_wt_min[A0_carry[5]][A0_AB_block[5]];
						A0_i5 <= cDDT_wt_min[A0_carry[5]][A0_AB_block[5]]; A0_i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
				{
					A0_w5 = A0_w4 + A0_i5;
					if(A0_w5 + A0_w_xor_block[6] +w_xor_1 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(A0_j5=0; A0_j5 < cDDT_n[A0_carry[5]][A0_AB_block[5]][A0_i5]; A0_j5++)
				{
					A0_gamma_bloc[5] = cDDT_v[A0_carry[5]][A0_AB_block[5]][A0_i5][A0_j5];
					A0_carry[6] = A0_carry_tmp[6] + (A0_gamma_bloc[5] >> 7); // gamma MSB
				for(A0_i6 = cDDT_wt_min[A0_carry[6]][A0_AB_block[6]];
						A0_i6 <= cDDT_wt_min[A0_carry[6]][A0_AB_block[6]]; A0_i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
				{
					A0_w6 = A0_w5 + A0_i6;
					if(A0_w6 + A0_w_xor_block[7] +w_xor_1 > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
				for(A0_j6=0; A0_j6 < cDDT_n[A0_carry[6]][A0_AB_block[6]][A0_i6]; A0_j6++)
				{
					A0_gamma_bloc[6] = cDDT_v[A0_carry[6]][A0_AB_block[6]][A0_i6][A0_j6];
					A0_carry[7] = A0_carry_tmp[7] + (A0_gamma_bloc[6] >> 7); // gamma MSB
				for(A0_i7 = MSB_cDDT_wt_min[A0_carry[7]][A0_AB_block[7]];
						A0_i7 <= MSB_cDDT_wt_min[A0_carry[7]][A0_AB_block[7]]; A0_i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
				{
					A0_w7 = A0_w6 + A0_i7;
					if(A0_w7 +w_xor_1> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
///////




/////////////////////第二个模加
		for(A1_i0 = cDDT_wt_min[A1_carry[0]][A1_AB_block[0]];
				A1_i0 <= cDDT_wt_max[A1_carry[0]][A1_AB_block[0]]; A1_i0++)
		{
			if(A0_w7 + A1_i0 + A1_w_xor_block[1] >= w_cmp ){break;}
		for(A1_j0=0; A1_j0 < cDDT_n[A1_carry[0]][A1_AB_block[0]][A1_i0]; A1_j0++)
		{
			A1_gamma_bloc[0] = cDDT_v[A1_carry[0]][A1_AB_block[0]][A1_i0][A1_j0];
			A1_carry[1] = A1_carry_tmp[1] + (A1_gamma_bloc[0] >> 7); // gamma MSB

			for(A1_i1 = cDDT_wt_min[A1_carry[1]][A1_AB_block[1]];
					A1_i1 <= cDDT_wt_max[A1_carry[1]][A1_AB_block[1]]; A1_i1++)
			{
				A1_w1 = A1_i0 + A1_i1;
				if(A0_w7 + A1_w1 + A1_w_xor_block[2] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j1=0; A1_j1 < cDDT_n[A1_carry[1]][A1_AB_block[1]][A1_i1]; A1_j1++)
			{
				A1_gamma_bloc[1] = cDDT_v[A1_carry[1]][A1_AB_block[1]][A1_i1][A1_j1];
				A1_carry[2] = A1_carry_tmp[2] + (A1_gamma_bloc[1] >> 7); // gamma MSB

			for(A1_i2 = cDDT_wt_min[A1_carry[2]][A1_AB_block[2]];
					A1_i2 <= cDDT_wt_max[A1_carry[2]][A1_AB_block[2]]; A1_i2++)
			{
				A1_w2 = A1_w1 + A1_i2;
				if(A0_w7 + A1_w2 + A1_w_xor_block[3] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j2=0; A1_j2 < cDDT_n[A1_carry[2]][A1_AB_block[2]][A1_i2]; A1_j2++)
			{
				A1_gamma_bloc[2] = cDDT_v[A1_carry[2]][A1_AB_block[2]][A1_i2][A1_j2];
				A1_carry[3] = A1_carry_tmp[3] + (A1_gamma_bloc[2] >> 7); // gamma MSB
			for(A1_i3 = cDDT_wt_min[A1_carry[3]][A1_AB_block[3]];
						A1_i3 <= MSB_cDDT_wt_max[A1_carry[3]][A1_AB_block[3]];A1_i3++)
			{
					A1_w3 = A1_w2 + A1_i3;
					if(A0_w7 + A1_w3 +A1_w_xor_block[4] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.

			for(A1_j3=0; A1_j3 < MSB_cDDT_n[A1_carry[3]][A1_AB_block[3]][A1_i3]; A1_j3++)
			{
				A1_gamma_bloc[3] = cDDT_v[A1_carry[3]][A1_AB_block[3]][A1_i3][A1_j3];
				A1_carry[4] = A1_carry_tmp[4] + (A1_gamma_bloc[3] >> 7); // gamma MSB
			for(A1_i4 = cDDT_wt_min[A1_carry[4]][A1_AB_block[4]];
							A1_i4 <= cDDT_wt_min[A1_carry[4]][A1_AB_block[4]];A1_i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
			{
				A1_w4 = A1_w3 + A1_i4;
				if(A0_w7 +A1_w4 + A1_w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j4=0; A1_j4 < cDDT_n[A1_carry[4]][A1_AB_block[4]][A1_i4]; A1_j4++)
			{
				A1_gamma_bloc[4] = cDDT_v[A1_carry[4]][A1_AB_block[4]][A1_i4][A1_j4];
				A1_carry[5] = A1_carry_tmp[5] + (A1_gamma_bloc[4] >> 7); // gamma MSB
			for(A1_i5 = cDDT_wt_min[A1_carry[5]][A1_AB_block[5]];
					A1_i5 <= cDDT_wt_min[A1_carry[5]][A1_AB_block[5]]; A1_i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
			{
				A1_w5 = A1_w4 + A1_i5;
				if(A0_w7 + A1_w5 + A1_w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(A1_j5=0; A1_j5 < cDDT_n[A1_carry[5]][A1_AB_block[5]][A1_i5]; A1_j5++)
			{
				A1_gamma_bloc[5] = cDDT_v[A1_carry[5]][A1_AB_block[5]][A1_i5][A1_j5];
				A1_carry[6] = A1_carry_tmp[6] + (A1_gamma_bloc[5] >> 7); // gamma MSB
			for(A1_i6 = cDDT_wt_min[A1_carry[6]][A1_AB_block[6]];
					A1_i6 <= cDDT_wt_min[A1_carry[6]][A1_AB_block[6]]; A1_i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
			{
				A1_w6 = A1_w5 + A1_i6;
				if(A0_w7 + A1_w6 + A1_w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
			for(A1_j6=0; A1_j6 < cDDT_n[A1_carry[6]][A1_AB_block[6]][A1_i6]; A1_j6++)
			{
				A1_gamma_bloc[6] = cDDT_v[A1_carry[6]][A1_AB_block[6]][A1_i6][A1_j6];
				A1_carry[7] = A1_carry_tmp[7] + (A1_gamma_bloc[6] >> 7); // gamma MSB
			for(A1_i7 = MSB_cDDT_wt_min[A1_carry[7]][A1_AB_block[7]];
					A1_i7 <= MSB_cDDT_wt_min[A1_carry[7]][A1_AB_block[7]]; A1_i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
			{
				A1_w7 = A1_w6 + A1_i7;

////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
		for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
			for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
			{
				w_cmp = w_cmp_rx - RX_zeta_w[zeta_num_0] - RX_zeta_w[zeta_num_1];
				if(A0_w7 + A1_w7 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				RX_Siphash_ADD_P_xor_w_tmp[0][current_round] = A0_w7;
				RX_Siphash_ADD_P_xor_w_tmp[1][current_round] = A1_w7;
				P_w[current_round] = A0_w7 + A1_w7;

				Siphash_Pr_zeta_tmp[0][current_round] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
				Siphash_Pr_zeta_tmp[1][current_round] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
				Siphash_zeta_tmp[0][current_round] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
				Siphash_zeta_tmp[1][current_round] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1
				RX_Siphash_P_hr_w_tmp[current_round] = P_w[current_round]
												    + Siphash_Pr_zeta_tmp[0][current_round]
													+ Siphash_Pr_zeta_tmp[1][current_round];
				RX_P_sumofR_w[current_round] = RX_Siphash_P_hr_w_tmp[current_round]
												+ RX_P_sumofR_w[current_round -1];

			for(A0_j7=0; A0_j7 < MSB_cDDT_n[A0_carry[7]][A0_AB_block[7]][A0_i7]; A0_j7++)
			{
				A0_gamma_bloc[7] = MSB_cDDT_v[A0_carry[7]][A0_AB_block[7]][A0_i7][A0_j7];

				//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
				A0_gamma_temp = (A0_gamma_bloc[7] << 56) | (A0_gamma_bloc[6] << 48)
						| (A0_gamma_bloc[5] << 40) | (A0_gamma_bloc[4] << 32)
						| (A0_gamma_bloc[3] << 24) | (A0_gamma_bloc[2] << 16)
						| (A0_gamma_bloc[1] <<8) | A0_gamma_bloc[0];
			for(A1_j7=0; A1_j7 < MSB_cDDT_n[A1_carry[7]][A1_AB_block[7]][A1_i7]; A1_j7++)
			{
				A1_gamma_bloc[7] = MSB_cDDT_v[A1_carry[7]][A1_AB_block[7]][A1_i7][A1_j7];

				//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
				A1_gamma_temp = (A1_gamma_bloc[7] << 56) | (A1_gamma_bloc[6] << 48)
						| (A1_gamma_bloc[5] << 40) | (A1_gamma_bloc[4] << 32)
						| (A1_gamma_bloc[3] << 24) | (A1_gamma_bloc[2] << 16)
						| (A1_gamma_bloc[1] <<8) | A1_gamma_bloc[0];

			out_RX[2] =  A0_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_0]];
			out_RX[1] = out_X1_tmp ^ out_RX[2];
			out_RX[2] = siphash_rol_left_32(out_RX[2]);
			out_RX[0] = A1_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_1]];
			out_RX[3] = out_RX[0] ^ out_X3_tmp;


			Siphash_RX_X_tmp[0][current_round +1] = out_RX[0];
			Siphash_RX_X_tmp[1][current_round +1] = out_RX[1];
			Siphash_RX_X_tmp[2][current_round +1] = out_RX[2];
			Siphash_RX_X_tmp[3][current_round +1] = out_RX[3];

			//进入下一个半轮,
			state = Siphash_RX_Middle_Rounds(search_round, current_round+1, out_RX);
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
	}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
	return updated;
}


u64 Siphash_RX_N_Rounds(u16 search_round,u32 *Xr)
{
	u64 state = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 xor_0 = 0, xor_1 = 0;
	u64 w_xor_0 = 0,w_xor_1 = 0;
	u64 w_xor = 0;
	u64 i = 0,j = 0;
	u64 A0_alpha_bloc[8] = {0};
	u64 A0_beta_bloc[8] = {0};
	u64 A0_gamma_bloc[8] = {0};
	u64 A0_gamma_temp = 0;
	u64 A0_AB_block[8] = {0};
	u64 A0_carry[8] ={0};
	u64 A0_carry_tmp[8] ={0};
	u64 A0_i0=0,A0_i1=0,A0_i2=0,A0_i3=0,A0_i4=0,A0_i5=0,A0_i6=0,A0_i7=0;
	u64 A0_j0=0,A0_j1=0,A0_j2=0,A0_j3=0,A0_j4=0,A0_j5=0,A0_j6=0,A0_j7=0;
	u32 A0_w1=0,A0_w2=0,A0_w3=0,A0_w4=0,A0_w5=0,A0_w6=0,A0_w7=0;
	u8 A0_w_xor_block[8] = {0};
	u64 A1_alpha_bloc[8] = {0};
	u64 A1_beta_bloc[8] = {0};
	u64 A1_gamma_bloc[8] = {0};
	u64 A1_gamma_temp = 0;
	u64 A1_AB_block[8] = {0};
	u64 A1_carry[8] ={0};
	u64 A1_carry_tmp[8] ={0};
	u64 A1_i0=0,A1_i1=0,A1_i2=0,A1_i3=0,A1_i4=0,A1_i5=0,A1_i6=0,A1_i7=0;
	u64 A1_j0=0,A1_j1=0,A1_j2=0,A1_j3=0,A1_j4=0,A1_j5=0,A1_j6=0,A1_j7=0;
	u64 A1_w1=0,A1_w2=0,A1_w3=0,A1_w4=0,A1_w5=0,A1_w6=0,A1_w7=0;
	u8 A1_w_xor_block[8] = {0};

	u64 out_X1_tmp=0, out_X3_tmp=0;
	u64 out_RX[4];

	u64 zeta_num_0 = 0,zeta_num_1 = 0;
	float w_cmp = 0;
	float w_cmp_rx = 0;
	u64 updata_cur_num = 0;
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
		out_X1_tmp = siphash_rol_left_17(Xr[1]);
		out_X3_tmp = siphash_rol_left_21(Xr[3]);
	}
	else //奇半轮，即某轮的上半轮，旋转参数为5/8/16
	{
		//rc=16;
		out_X1_tmp = siphash_rol_left_13(Xr[1]);
		out_X3_tmp = siphash_rol_left_16(Xr[3]);
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
						if(A0_w3 + A0_w_xor_block[4] +w_xor_1>= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A0_j3=0; A0_j3 < MSB_cDDT_n[A0_carry[3]][A0_AB_block[3]][A0_i3]; A0_j3++)
					{
						A0_gamma_bloc[3] = cDDT_v[A0_carry[3]][A0_AB_block[3]][A0_i3][A0_j3];
						A0_carry[4] = A0_carry_tmp[4] + (A0_gamma_bloc[3] >> 7); // gamma MSB
						for(A0_i4 = cDDT_wt_min[A0_carry[4]][A0_AB_block[4]];
								A0_i4 <= cDDT_wt_min[A0_carry[4]][A0_AB_block[4]]; A0_i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
						{
							A0_w4 = A0_w3 + A0_i4;
							if(A0_w4 + A0_w_xor_block[5] +w_xor_1> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
						for(A0_j4=0; A0_j4 < cDDT_n[A0_carry[4]][A0_AB_block[4]][A0_i4]; A0_j4++)
						{
							A0_gamma_bloc[4] = cDDT_v[A0_carry[4]][A0_AB_block[4]][A0_i4][A0_j4];
							A0_carry[5] = A0_carry_tmp[5] + (A0_gamma_bloc[4] >> 7); // gamma MSB
						for(A0_i5 = cDDT_wt_min[A0_carry[5]][A0_AB_block[5]];
								A0_i5 <= cDDT_wt_min[A0_carry[5]][A0_AB_block[5]]; A0_i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
						{
							A0_w5 = A0_w4 + A0_i5;
							if(A0_w5 + A0_w_xor_block[6] +w_xor_1 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
						for(A0_j5=0; A0_j5 < cDDT_n[A0_carry[5]][A0_AB_block[5]][A0_i5]; A0_j5++)
						{
							A0_gamma_bloc[5] = cDDT_v[A0_carry[5]][A0_AB_block[5]][A0_i5][A0_j5];
							A0_carry[6] = A0_carry_tmp[6] + (A0_gamma_bloc[5] >> 7); // gamma MSB
						for(A0_i6 = cDDT_wt_min[A0_carry[6]][A0_AB_block[6]];
								A0_i6 <= cDDT_wt_min[A0_carry[6]][A0_AB_block[6]]; A0_i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
						{
							A0_w6 = A0_w5 + A0_i6;
							if(A0_w6 + A0_w_xor_block[7] +w_xor_1 > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
						for(A0_j6=0; A0_j6 < cDDT_n[A0_carry[6]][A0_AB_block[6]][A0_i6]; A0_j6++)
						{
							A0_gamma_bloc[6] = cDDT_v[A0_carry[6]][A0_AB_block[6]][A0_i6][A0_j6];
							A0_carry[7] = A0_carry_tmp[7] + (A0_gamma_bloc[6] >> 7); // gamma MSB
						for(A0_i7 = MSB_cDDT_wt_min[A0_carry[7]][A0_AB_block[7]];
								A0_i7 <= MSB_cDDT_wt_min[A0_carry[7]][A0_AB_block[7]]; A0_i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
						{
							A0_w7 = A0_w6 + A0_i7;
							if(A0_w7 +w_xor_1> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		///////




		/////////////////////第二个模加
				for(A1_i0 = cDDT_wt_min[A1_carry[0]][A1_AB_block[0]];
						A1_i0 <= cDDT_wt_max[A1_carry[0]][A1_AB_block[0]]; A1_i0++)
				{
					if(A0_w7 + A1_i0 + A1_w_xor_block[1] >= w_cmp ){break;}
				for(A1_j0=0; A1_j0 < cDDT_n[A1_carry[0]][A1_AB_block[0]][A1_i0]; A1_j0++)
				{
					A1_gamma_bloc[0] = cDDT_v[A1_carry[0]][A1_AB_block[0]][A1_i0][A1_j0];
					A1_carry[1] = A1_carry_tmp[1] + (A1_gamma_bloc[0] >> 7); // gamma MSB

					for(A1_i1 = cDDT_wt_min[A1_carry[1]][A1_AB_block[1]];
							A1_i1 <= cDDT_wt_max[A1_carry[1]][A1_AB_block[1]]; A1_i1++)
					{
						A1_w1 = A1_i0 + A1_i1;
						if(A0_w7 + A1_w1 + A1_w_xor_block[2] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A1_j1=0; A1_j1 < cDDT_n[A1_carry[1]][A1_AB_block[1]][A1_i1]; A1_j1++)
					{
						A1_gamma_bloc[1] = cDDT_v[A1_carry[1]][A1_AB_block[1]][A1_i1][A1_j1];
						A1_carry[2] = A1_carry_tmp[2] + (A1_gamma_bloc[1] >> 7); // gamma MSB

					for(A1_i2 = cDDT_wt_min[A1_carry[2]][A1_AB_block[2]];
							A1_i2 <= cDDT_wt_max[A1_carry[2]][A1_AB_block[2]]; A1_i2++)
					{
						A1_w2 = A1_w1 + A1_i2;
						if(A0_w7 + A1_w2 + A1_w_xor_block[3] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A1_j2=0; A1_j2 < cDDT_n[A1_carry[2]][A1_AB_block[2]][A1_i2]; A1_j2++)
					{
						A1_gamma_bloc[2] = cDDT_v[A1_carry[2]][A1_AB_block[2]][A1_i2][A1_j2];
						A1_carry[3] = A1_carry_tmp[3] + (A1_gamma_bloc[2] >> 7); // gamma MSB
					for(A1_i3 = cDDT_wt_min[A1_carry[3]][A1_AB_block[3]];
								A1_i3 <= MSB_cDDT_wt_max[A1_carry[3]][A1_AB_block[3]];A1_i3++)
					{
							A1_w3 = A1_w2 + A1_i3;
							if(A0_w7 + A1_w3 +A1_w_xor_block[4] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.

					for(A1_j3=0; A1_j3 < MSB_cDDT_n[A1_carry[3]][A1_AB_block[3]][A1_i3]; A1_j3++)
					{
						A1_gamma_bloc[3] = cDDT_v[A1_carry[3]][A1_AB_block[3]][A1_i3][A1_j3];
						A1_carry[4] = A1_carry_tmp[4] + (A1_gamma_bloc[3] >> 7); // gamma MSB
					for(A1_i4 = cDDT_wt_min[A1_carry[4]][A1_AB_block[4]];
									A1_i4 <= cDDT_wt_min[A1_carry[4]][A1_AB_block[4]];A1_i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
					{
						A1_w4 = A1_w3 + A1_i4;
						if(A0_w7 +A1_w4 + A1_w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A1_j4=0; A1_j4 < cDDT_n[A1_carry[4]][A1_AB_block[4]][A1_i4]; A1_j4++)
					{
						A1_gamma_bloc[4] = cDDT_v[A1_carry[4]][A1_AB_block[4]][A1_i4][A1_j4];
						A1_carry[5] = A1_carry_tmp[5] + (A1_gamma_bloc[4] >> 7); // gamma MSB
					for(A1_i5 = cDDT_wt_min[A1_carry[5]][A1_AB_block[5]];
							A1_i5 <= cDDT_wt_min[A1_carry[5]][A1_AB_block[5]]; A1_i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
					{
						A1_w5 = A1_w4 + A1_i5;
						if(A0_w7 + A1_w5 + A1_w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(A1_j5=0; A1_j5 < cDDT_n[A1_carry[5]][A1_AB_block[5]][A1_i5]; A1_j5++)
					{
						A1_gamma_bloc[5] = cDDT_v[A1_carry[5]][A1_AB_block[5]][A1_i5][A1_j5];
						A1_carry[6] = A1_carry_tmp[6] + (A1_gamma_bloc[5] >> 7); // gamma MSB
					for(A1_i6 = cDDT_wt_min[A1_carry[6]][A1_AB_block[6]];
							A1_i6 <= cDDT_wt_min[A1_carry[6]][A1_AB_block[6]]; A1_i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
					{
						A1_w6 = A1_w5 + A1_i6;
						if(A0_w7 + A1_w6 + A1_w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
					for(A1_j6=0; A1_j6 < cDDT_n[A1_carry[6]][A1_AB_block[6]][A1_i6]; A1_j6++)
					{
						A1_gamma_bloc[6] = cDDT_v[A1_carry[6]][A1_AB_block[6]][A1_i6][A1_j6];
						A1_carry[7] = A1_carry_tmp[7] + (A1_gamma_bloc[6] >> 7); // gamma MSB
					for(A1_i7 = MSB_cDDT_wt_min[A1_carry[7]][A1_AB_block[7]];
							A1_i7 <= MSB_cDDT_wt_min[A1_carry[7]][A1_AB_block[7]]; A1_i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
					{
						A1_w7 = A1_w6 + A1_i7;

		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
				//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
//				for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
//					for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
					{
						zeta_num_0=0; zeta_num_1=0;
						w_cmp = w_cmp_rx - RX_zeta_w[zeta_num_0] - RX_zeta_w[zeta_num_1];
						if(A0_w7 + A1_w7 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
						w_cmp = A0_w7 + A1_w7;
						RX_Siphash_ADD_P_xor_w_tmp[0][search_round] = A0_w7;
						RX_Siphash_ADD_P_xor_w_tmp[1][search_round] = A1_w7;
						P_w[search_round] = A0_w7 + A1_w7;

						Siphash_Pr_zeta_tmp[0][search_round] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
						Siphash_Pr_zeta_tmp[1][search_round] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
						Siphash_zeta_tmp[0][search_round] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
						Siphash_zeta_tmp[1][search_round] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1
						RX_Siphash_P_hr_w_tmp[search_round] = P_w[search_round]
														    + Siphash_Pr_zeta_tmp[0][search_round]
															+ Siphash_Pr_zeta_tmp[1][search_round];
						RX_P_sumofR_w[search_round] = RX_Siphash_P_hr_w_tmp[search_round]
														+ RX_P_sumofR_w[search_round -1];

					//for(A0_j7=0; A0_j7 < MSB_cDDT_n[A0_carry[7]][A0_AB_block[7]][A0_i7]; A0_j7++)
					{
						A0_j7=0; //取1个即可
						A0_gamma_bloc[7] = MSB_cDDT_v[A0_carry[7]][A0_AB_block[7]][A0_i7][A0_j7];

						//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
						A0_gamma_temp = (A0_gamma_bloc[7] << 56) | (A0_gamma_bloc[6] << 48)
								| (A0_gamma_bloc[5] << 40) | (A0_gamma_bloc[4] << 32)
								| (A0_gamma_bloc[3] << 24) | (A0_gamma_bloc[2] << 16)
								| (A0_gamma_bloc[1] <<8) | A0_gamma_bloc[0];
					//for(A1_j7=0; A1_j7 < MSB_cDDT_n[A1_carry[7]][A1_AB_block[7]][A1_i7]; A1_j7++)
					{
						A1_j7=0;  //取1个即可
						A1_gamma_bloc[7] = MSB_cDDT_v[A1_carry[7]][A1_AB_block[7]][A1_i7][A1_j7];

						//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
						A1_gamma_temp = (A1_gamma_bloc[7] << 56) | (A1_gamma_bloc[6] << 48)
								| (A1_gamma_bloc[5] << 40) | (A1_gamma_bloc[4] << 32)
								| (A1_gamma_bloc[3] << 24) | (A1_gamma_bloc[2] << 16)
								| (A1_gamma_bloc[1] <<8) | A1_gamma_bloc[0];


//%%%%%%%%基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减%%%%%%%%%%%%%%%%%%%%%%
	//	for(zeta_num_0=0; zeta_num_0 < RX_zeta_total ; zeta_num_0++)
	//	for(zeta_num_1=0; zeta_num_1 < RX_zeta_total ; zeta_num_1++)
		{
	//		zeta_num_0=0; zeta_num_1=0;  //启发式，只取EK_min
	//		w_cmp = w_cmp_rx - RX_zeta_w[zeta_num_0] - RX_zeta_w[zeta_num_1];
	//		if(A0_w3 + A1_w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.

		RX_Siphash_ADD_P_xor_w_tmp[0][search_round] = A0_w3;
		RX_Siphash_ADD_P_xor_w_tmp[1][search_round] = A1_w3;
		P_w[search_round] = A0_w3 + A1_w3;
		Siphash_Pr_zeta_tmp[0][search_round] = RX_zeta_w[zeta_num_0];  //Pr_zeta0
		Siphash_Pr_zeta_tmp[1][search_round] = RX_zeta_w[zeta_num_1];  //Pr_zeta1
		Siphash_zeta_tmp[0][search_round] = RX_zeta[RX_zeta_i[zeta_num_0]]; //zeta0
		Siphash_zeta_tmp[1][search_round] = RX_zeta[RX_zeta_i[zeta_num_1]]; //zeta1


		RX_Siphash_P_hr_w_tmp[search_round] = P_w[search_round]
										    + Siphash_Pr_zeta_tmp[0][search_round]
											+ Siphash_Pr_zeta_tmp[1][search_round];

		out_RX[2] =  A0_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_0]];
		out_RX[1] = out_X1_tmp ^ out_RX[2];
		out_RX[2] = siphash_rol_left_32(out_RX[2]);
		out_RX[0] = A1_gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num_1]];
		out_RX[3] = out_RX[0] ^ out_X3_tmp;
		Siphash_RX_X[0][search_round +1] = out_RX[0];
		Siphash_RX_X[1][search_round +1] = out_RX[1];
		Siphash_RX_X[2][search_round +1] = out_RX[2];
		Siphash_RX_X[3][search_round +1] = out_RX[3];

		//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
		for(i=1; i <= search_round; i++ )
		{
			//记录更新后的每个半轮的RX差分
			Siphash_RX_X[0][i] = Siphash_RX_X_tmp[0][i];
			Siphash_RX_X[1][i] = Siphash_RX_X_tmp[1][i];
			Siphash_RX_X[2][i] = Siphash_RX_X_tmp[2][i];
			Siphash_RX_X[3][i] = Siphash_RX_X_tmp[3][i];
			//记录更新后的每个半轮的zeta
			Siphash_Zeta[0][i] = Siphash_zeta_tmp[0][i];
			Siphash_Zeta[1][i] = Siphash_zeta_tmp[1][i];
			//记录更新后的每个半轮的zeta的概率
			Siphash_Pr_zeta[0][i] = Siphash_Pr_zeta_tmp[0][i];
			Siphash_Pr_zeta[1][i] = Siphash_Pr_zeta_tmp[1][i];
			//记录更新后的每个半轮的模加的XOR差分概率
			RX_Siphash_ADD_P_xor_w[0][i] = RX_Siphash_ADD_P_xor_w_tmp[0][i];
			RX_Siphash_ADD_P_xor_w[1][i] = RX_Siphash_ADD_P_xor_w_tmp[1][i];
			//记录更新后的每个半轮的RX总的概率重量，包括XOR和zeta
			RX_Siphash_P_hr_w[i] = RX_Siphash_P_hr_w_tmp[i];
		}
		//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
				Num_Bn_Update++;  //更新次数记录
				Bn_cmp = 0;
				for(i=1; i <= search_round; i++ )
				{
					Bn_cmp += RX_Siphash_P_hr_w[i];
				}
				RX_Expected_w = Bn_cmp;
//				w_cmp = RX_Expected_w - Siphash_Pr_zeta_tmp[0][search_round]
//							- Siphash_Pr_zeta_tmp[1][search_round]
//						    - RX_P_sumofR_w[search_round -1];
//if(P_w[search_round] >= w_cmp)

				state = 1;
				updated = 1;

				//%%%%%打印更新次数等信息
				time_Round = clock();
				run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
				printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
				printf("Time: %.2f seconds. ==>>\n", run_time);

				}
				}
}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
	return updated;
}





u64 print_Siphash_RX_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("-r----X0--------X1-------X2-------X3------zeta_0---zeta_1----P_zeta_0<->P_zeta_3---P_xor_0---P_xor_1---Pw_hr------Pw_r\n");
	printf("-------------------\n");
	for(r_num=1;r_num <= search_round; )
	{
		printf("%02d  %16.16llx %16.16llx %16.16llx %16.16llx   %16.16llx %16.16llx \n",
				r_num,
				Siphash_RX_X[0][r_num],Siphash_RX_X[1][r_num],Siphash_RX_X[2][r_num],Siphash_RX_X[3][r_num],
				Siphash_Zeta[0][r_num],Siphash_Zeta[1][r_num]);
		printf("P_zeta_0 = %f   P_zeta_1 = %f   P_xor_0 = %f   P_xor_1 = %f   Pw_hr = %f     ----\n",
				Siphash_Pr_zeta[0][r_num],Siphash_Pr_zeta[1][r_num],
				RX_Siphash_ADD_P_xor_w[0][r_num],RX_Siphash_ADD_P_xor_w[1][r_num],
				RX_Siphash_P_hr_w[r_num]);

//printf("hr----XV0--------XV1---------XV2---------XV3------zeta_0------zeta_1------P_zeta_0---P_zeta_0----P_xor_0---P_xor_1----Pw_hr \n");
		printf("%02d  %16.16llx %16.16llx %16.16llx %16.16llx   %16.16llx %16.16llx \n",
				r_num+1,
				Siphash_RX_X[0][r_num+1],Siphash_RX_X[1][r_num+1],Siphash_RX_X[2][r_num+1],Siphash_RX_X[3][r_num+1],
				Siphash_Zeta[0][r_num+1],Siphash_Zeta[1][r_num+1]);

		printf("P_zeta_2 = %f   P_zeta_3 = %f   P_xor_2 = %f   P_xor_3 = %f   Pw_hr = %f   Pw_r = %f \n",
				Siphash_Pr_zeta[0][r_num+1],Siphash_Pr_zeta[1][r_num+1],
				RX_Siphash_ADD_P_xor_w[0][r_num+1],RX_Siphash_ADD_P_xor_w[1][r_num+1],
				RX_Siphash_P_hr_w[r_num+1],
				RX_Siphash_P_hr_w[r_num] + RX_Siphash_P_hr_w[r_num+1] ); //

		r_num = r_num + 2;
		printf("-------------------\n" );
	}
	printf("out %16.16llx %16.16llx %16.16llx %16.16llx   \n",
			Siphash_RX_X[0][search_round+1],Siphash_RX_X[1][search_round+1],
			Siphash_RX_X[2][search_round+1],Siphash_RX_X[3][search_round+1]);

	printf("------------------\n");
	printf("%d-half-Round Siphash Total Weight: -%f \n",search_round, RX_Expected_w);

	return 1;


}











