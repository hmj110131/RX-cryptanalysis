/*
 * RX_Alzette.c
 *
 *  Created on: 2020年7月24日
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
#include "Alzette.h"
#include "RX_Alzette.h"
#include "XOR_offset_table.h"
#include <nmmintrin.h>


volatile u32 Alzette_RX[2][16]={0}; //每轮的输入差分左右两部分
volatile u32 Alzette_Zeta[16]={0}; //每轮的输zeta
volatile u32 Alzette_zeta_tmp[16]={0}; //每轮的临时zeta
volatile float Alzette_Pr_RX[16]={0}; //每轮RX差分概率
volatile float Alzette_Pr_add[16]={0}; //每轮RX差分概率
volatile float Alzette_Pr_zeta[16]={0}; //每轮RX差分概率

u32  Alzette_C_RX[8]={0};     //定义Alzette的实例版本常数RX差分，共8个版本
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
u32 Alzette_Constant_RX(u16 rx_k)
{
	u16 i=0;

	for(i=0; i<8; i++)
	{
		Alzette_C_RX[i] = Alzette_Constant[i] ^
				(rol_left_32(Alzette_Constant[i], rx_k));
	}
	return 0;
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
////////////////////////////////搜索Alzette的最优 差分 特征的代码/////////////////////////
u32 Alzette_RX_Diff_trail_search_entry(u16 search_round)
{
	u32 status = 0;
	u16 i = 0, j=0;
	u32 tmp=0,tmp_i=0;
	u16 t_index = 0;
	FILE* Alzette_RX_diff_Bn;
	float Bn_wt_inc_max = 0;
	float zeta_tmp = 0;
//	char strValue[64];


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	// For Alzette modular additon is 32 bits.
	blocksize_len = 32;
	nBytes = 4;  //分为几个8比特的块
	Bit_Align = 0xFFFFFFFF;
	ValueMax_Align = 0x7FFFFFFF;
	V_MSB = 0x80000000;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	printf("Constructing Alzette_Constant_RX ... \n");
	Alzette_Constant_RX (RX_k);  //计算Alzette常数的RX差分
//	printf("|------------------------- k=%d -----------------------|\n", RX_k);
//	for(i=0; i<8; i++)printf("Alzette_C_RX[%d]: %x  \n",i,Alzette_C_RX[i]);
//	printf("k=%d   Alzette_C[%d]: %x   Alzette_C_RX[%d]: %x \n",
//			RX_k, Alzette_C_i, Alzette_Constant[Alzette_C_i], Alzette_C_i,Alzette_C_RX[Alzette_C_i]);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	RX_zeta_total = 1 + blocksize_len + RX_k*(blocksize_len -RX_k);  //计算Zeta的总个数
	printf("RX_zeta_total: %d  \n",RX_zeta_total);
	printf("Constructing RX_Ek[4] ... \n");
	Compute_RX_Ek(RX_k, (u16)blocksize_len,  Ek); //计算ek的4个值
	Ek_min = Compute_RX_Ek_min(Ek);
//	printf("k=%d  RX_Ek[0]: %f RX_Ek[1]: %f RX_Ek[2]: %f RX_Ek[3]: %f \n",
//			RX_k,Ek[0],Ek[1],Ek[2],Ek[3]);

	printf("Constructing RX_XOR_offset_Table ... \n");
	Compute_RX_V_U(RX_k, blocksize_len, RX_v,RX_v_w, RX_u,RX_u_w);
	Sorted_RX_XOR_offset_Table(RX_k, (u16)blocksize_len,
			RX_v,RX_v_w, RX_u,RX_u_w,
			Ek,
			RX_zeta, RX_zeta_w, RX_zeta_i);
//-------------取前r-1轮最优概率重量，搜索期望概率重量的最大值----------------
	Alzette_RX_diff_Bn = fopen ("../tmp/Alzette_RX_diff_Bn.xlsx", "a+"); //  "w+"); //
	for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
	{
		fscanf(Alzette_RX_diff_Bn,"%f", &RX_P_bestofR_w[i]);  //read.
//		fscanf(Alzette_RX_diff_Bn,"%s", strValue);  //read.
//		printf("strValue[%d]: %s \n",i,strValue);  // write.
//		RX_P_bestofR_w[i] = atof(strValue);
//		printf("RX_P_bestofR_w[%d]: %f \n",i,RX_P_bestofR_w[i]);  // write.
	}
/*
	Bn_wt_inc_max = wt_inc_app(F_Alpha,
			rol_right_32(F_Beta,Alzette_rol_a[search_round-1]),
			ValueMax_Align); //计算最大概率重量增量
	printf("Bn_wt_inc_max: %f   \n",Bn_wt_inc_max);
//	printf("Bn_wt_inc_max: %f  RX_P_bestofR_w[search_round - 1]:%f \n",Bn_wt_inc_max, RX_P_bestofR_w[search_round - 1]);
	//以前一轮的最优RX差分概率重量，加上以其输出RX差分扩展1轮的概率重量为最大的概率重量，即概率下界逼近最优
	RX_Expected_w = Bn_wt_inc_max + RX_P_bestofR_w[search_round - 1];
//	RX_Expected_w = 8;  //!!!!!启发式，可以通过设定直接的开始期望概率重量!!!!!//
	printf("RX_Expected_w: %f  \n", RX_Expected_w);
	*/
//-------记录时间-----
	time_xor_talbe = clock();
	run_time =  (time_xor_talbe - time_start) / CLOCKS_PER_SEC;
	printf("Time of Preprocessing Alzette info: %.2f seconds.  \n", run_time);
	printf("|--------------------------------------------------|\n");
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	printf("Constructing cDDT Lookup Tables ... \n");
	//构造cDDT
	ARX_carry_DDTm_construct();  // FOR cDDT
	//fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
	time_ARX_DDT = clock();
	run_time =  (time_ARX_DDT - time_xor_talbe) / CLOCKS_PER_SEC;
	printf("Time of Construct Alzette cDDT: %.2f seconds.  \n", run_time);
	//printf("|--------------------------------------------------|\n");
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int RX_Expected_w_int = 0;
	RX_Expected_w_int = (int)RX_P_bestofR_w[search_round - 1];
	RX_Expected_w = RX_Expected_w_int;
	RX_Expected_w = 42;  ////!!!!!启发式，可以通过设定直接的开始期望概率重量!!!!!//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
printf("*****************************Searching processing**************************|\n");
	do{
		RX_Expected_w = RX_Expected_w + 1; // 从r-1轮的概率重量整数部分+1开始
		printf("Update: %d   RX_Bn_w: %f  --Max expected weight-- \n",
				Num_Bn_Update, RX_Expected_w);

		// Search Entry.
		status = Alzette_RX_Diff_trail_round_1(search_round); //效率高
//		if((status != 2)|| (Num_Bn_Update == 0))
		{
			//printf("******* Not updated anymore. status:%d******* \n",status);
//		printf("***** Not updated anymore. status:%d  Num_Bn_Update: %d ***** \n",
//					status, Num_Bn_Update);
		}

		time_Round = clock();
		run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
		printf("Time: %.2f seconds.\n", run_time);

	}while( Num_Bn_Update == 0 );   //end conditions. // 概率重量逼近次数大于1，且返回值表示搜索到最优路径

	fprintf(Alzette_RX_diff_Bn,"%f \r\n",RX_Expected_w);  // write.RX_P_bestofR_w[search_round]
	printf("**************************************************************************|\n");
	print_Alzette_RX_diff_resoult(search_round);  //打印

	time_finish = clock();
	run_time = (float)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
	printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours == %.2f days. \n", run_time,run_time/60.0,run_time/3600.0,run_time/86400.0 );
	printf("Auto-search %d-Round Alzette optimal RX-differential trails END! \n", search_round);
	printf("|************************************************************************|\n");
	fclose(Alzette_RX_diff_Bn);

	return 0;  //status
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
u32 Alzette_RX_Diff_trail_round_1(u16 search_round)
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
 	printf("**Noted**: for(thr_d = 0; thr_d < blocksize_len-1; thr_d++) \n"); //打印启发信息
// 可以选择是否限定第1轮的概率重量? thr_d的搜索范围直接影响复杂度，可考虑限定。
		for(thr_d = 0; thr_d < blocksize_len-1; thr_d++)  //0::n-1 //blocksize_len-1
		{
			// thr_d : the bits number of alpha/beta/gamma that the position not equal at the same time.
			if ((thr_d +Ek_min + RX_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{	//统一退出位置
				return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
			}

			M0 = thr_d;  //三个alpha,beta和gamma同时考虑
			P_w[1] = M0;  //pro;
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

					Alzette_X[0][0] = Input_alpha;
					Alzette_X[1][0] = Alzette_rol_left_31(Input_beta);

					for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
					{
						RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
						if ((RX_P_w[1] + RX_P_bestofR_w[search_round - 1])
								>= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
						{
							break;
						}
						RX_P_sumofR_w[1] = RX_P_w[1];
						Alzette_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
						Alzette_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

					Alzette_X[0][1] = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]];
					Alzette_X[1][1] = Alzette_X[1][0] ^ (Alzette_rol_right_24(Alzette_X[0][1]));
					Alzette_X[0][1] ^= Alzette_C_RX[Alzette_C_i];
				//	Alzette_X[0][1] ^= 3;  //测试C=1

				best = Alzette_RX_Diff_trail_round_r(search_round, 2,
								Alzette_X[0][1],
								Alzette_X[1][1]);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		if(updata_cur_num < Num_Bn_Update)
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知已经更新
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if((RX_P_w[1] +	RX_P_bestofR_w[search_round - 1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
		{
			updated = 2;
			break;
		}
	}}
#endif
					}
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

/*
					if(Bn_w == 5 ){
				//第一种模式输出 //These bits position are SET. LSB-->MSB.
					for( i=1; i<=M0; i++)
					{
						printf("%d\t",C0[i]);
					}
				printf("\n");
					}
*/
					best = RX_Alzette_input_MSB(search_round,
							Input_alpha,
							Input_beta,
							Input_gamma,thr_d,C0);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(best == 2)
		{
			updated = 2;
			return updated;
		}
		updated = 1;  //准备通知已经更新
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
			best = RX_Alzette_input_MSB(search_round,
					Input_alpha,
					Input_beta,
					Input_gamma,thr_d,C0);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		if(best == 2)
		{
			updated = 2;
			return updated;
		}
		updated = 1;  //准备通知已经更新
	}
#endif
///////////////////////////////////////
					} while(1);
				}
		}
	return updated;
}



u16 RX_Alzette_input_MSB
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 Input_alpha = 0;
	u32 Input_beta = 0;
	u32 Input_gamma = 0;
	u16 j_abc = 0;

/*  //只考虑了set-A
	for(j_abc = 0; j_abc < 4;j_abc++)   // MSB
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
*/


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

			state = RX_Alzette_input_Middle(search_round,
					Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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

		state = RX_Alzette_input_Middle(search_round,
				Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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

//由高位到低位逐步判断
u16 RX_Alzette_input_Middle
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi, u16 cur_posi)
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
			state = RX_Alzette_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
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
	else 	//call last position.
	{
		state = RX_Alzette_input_Last(search_round,
				alpha,beta,	gamma,P_1,posi );
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
	return updated;
}

u16 RX_Alzette_input_Last
(u16 search_round,u32 alpha, u32 beta, u32 gamma,u16 P_1, char *posi )
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

		//基于构造的第1轮的输入输出差分组合，调用下一轮
		Alzette_X[0][0] = Input_alpha;
		Alzette_X[1][0] = Alzette_rol_left_31(Input_beta);

		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
			if((RX_P_w[1] +	RX_P_bestofR_w[search_round - 1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				break;
			}
			RX_P_sumofR_w[1] = RX_P_w[1];
			Alzette_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
			Alzette_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

		Alzette_X[0][1] = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]];
		Alzette_X[1][1] = Alzette_X[1][0] ^	(Alzette_rol_right_24(Alzette_X[0][1]));
		Alzette_X[0][1] ^= Alzette_C_RX[Alzette_C_i];
	//	Alzette_X[0][1] ^= 3;  //测试C=1

		state = Alzette_RX_Diff_trail_round_r(search_round, 2,
					Alzette_X[0][1],
					Alzette_X[1][1]);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(state == 1 )  //判断本轮是否要更新RX_Expected_w
	{
		if(updata_cur_num < Num_Bn_Update)
		{
			updata_cur_num = Num_Bn_Update;
			updated = 1;  //准备通知上一轮最新的RX_Expected_w

		//更新了期望概率重量RX_Expected_w//先考虑是否直接退出本次的PW
		if((P_w[1] + RX_zeta_w[0] +	RX_P_bestofR_w[search_round - 1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
		{
			updated = 2;  //准备通知上一轮，准备结束搜索
			return updated;  //这里是跳出for(j_abc = 0;
		}
		////更新了期望概率重量RX_Expected_w，考虑是否搜剩下的组合%%%%%%%%%%%%%%%%%%%%%%%
		if((RX_P_w[1] +	RX_P_bestofR_w[search_round - 1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
		{
			break;   //这里是退出for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		}
	}}
#endif
		}
	}
	return updated;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

u32 Alzette_RX_Diff_trail_round_r(u16 search_round, u16 cur_round, u32 x, u32 y)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 Alpha = 0;
	u32 Beta = 0;
	u32 gamma_temp = 0;
	u32 xor_alpha_beta = 0;
	u16 w_xor = 0;
	float w_cmp = 0;
	u32 i = 0,j = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u32 alpha_bloc[8] = {0};
	u32 beta_bloc[8] = {0};
	u32 gamma_bloc[8] = {0};
	u32 AB_block[8] = {0};
	u32 carry[8] ={0};
	u32 carry_tmp[8] ={0};
	u8 w_xor_block[8] = {0};
	u32 i0=0,i1=0,i2=0,i3=0;
	u32 j0=0,j1=0,j2=0,j3=0;
	u32 w1=0,w2=0,w3=0;

	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制

	updata_cur_num = Num_Bn_Update;

	if(search_round == cur_round)
	{
		updated = Alzette_RX_Diff_trail_round_N(search_round, x, y );
		return updated;
	}

	Alpha = x;
	Beta =rol_right_32(y,Alzette_rol_a[cur_round-1]);
	xor_alpha_beta = (Beta ^ Alpha) & ValueMax_Align ;
	w_xor = HM_weight(xor_alpha_beta );

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
				- RX_P_bestofR_w[search_round - cur_round];
	w_cmp = w_cmp_rx - Ek_min;
	if ( w_xor >= w_cmp)
	{
		return 0;
	}


////////////////////////////选择对应64bit分组版本////////////////////////////
	for(i=0;i<nBytes;i++)  //nBytes
	{
		alpha_bloc[i] = (Alpha >> (8*i)) & 0xFF; //8 bit
		beta_bloc[i]  = (Beta >> (8*i))  & 0xFF; //8 bit
		AB_block[i] = ((alpha_bloc[i] << 8) | beta_bloc[i]);
	}
	for(j=1;j<nBytes;j++)  //nBytes
	{
		w_xor_block[j] = (u8)_mm_popcnt_u64((xor_alpha_beta >> (j*8)));
	}
	carry[0] = 0;
	for(j=1;j<nBytes;j++)  //nBytes
	{
		carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
			+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
	}


/////////////////////
	for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
	{
		if(i0 + w_xor_block[1] >= w_cmp ){break;}
	for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
	{
		gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
		//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
		//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
		carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

	for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
	{
		w1 = i0 + i1;
		if(w1 + w_xor_block[2] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
	for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
	{
		gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
		//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
		//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
		carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

	for(i2 = cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
	{
		w2 = w1 + i2;
		if(w2 + w_xor_block[3] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
	for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
	{
		gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
		//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
		//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
		carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB

	for(i3 = MSB_cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= MSB_cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
	{
		w3 = w2 + i3;
		if(w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		/////////////////////
		//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			w_cmp = w_cmp_rx - RX_zeta_w[zeta_num];
			if(w3 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			Alzette_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
			Alzette_zeta_tmp[cur_round] =RX_zeta[RX_zeta_i[zeta_num]];
			P_w[cur_round] = w3; //i0 + i1 +i2 +i3;
			//p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];
			RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
			RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

			//若增加zeta变量后，概率仍然满足条件，再取遍历最后一个block并组合可能输出RX差分
		for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
			gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					   | (gamma_bloc[1] <<8) | gamma_bloc[0];

		Alzette_X[0][cur_round] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]];
		Alzette_X[1][cur_round] = rol_right_32(Alzette_X[0][cur_round],
										Alzette_rol_b[cur_round-1])^y;
		Alzette_X[0][cur_round] ^= Alzette_C_RX[Alzette_C_i];
	//	Alzette_X[0][cur_round] ^= 3;  //测试C=1

		best = Alzette_RX_Diff_trail_round_r(search_round,
				cur_round +1,
				Alzette_X[0][cur_round],
				Alzette_X[1][cur_round]);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best == 1 )  //判断本轮是否要更新RX_Expected_w
	{
		if(updata_cur_num < Num_Bn_Update)   //快速退出机制
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知上一轮最新的RX_Expected_w
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%//快速退出机制//%%%%%%%%
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_P_bestofR_w[search_round - cur_round];
		w_cmp = w_cmp_rx - RX_zeta_w[zeta_num];
		if(w3 >= w_cmp ){
			breakfor = 1;
			break;
		}  //这里有问题，只应该退出i3的嵌套
	}}
#endif
		}
		if(breakfor ==1 ){
			break;} //退出for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		}
		if(breakfor ==1 ){
			break;} //退出for(i3 = cDDT_wt_min
	}
	if(breakfor ==1 ){
		if(w2 + w_xor_block[3] >= w_cmp ){break;}  //退出for(j2=0;
		else
		{ breakfor = 0; } //清楚breakfor，继续执行j2中的gamma
		} //退出for(j2
	}
/*	if(breakfor ==1 ){
		if(w2 + w_xor_block[3] >= w_cmp ){break;} //退出w2
		else
		{	breakfor = 0; } //清楚breakfor，继续增长w2
		} //退出for(w2=
	*/
	}
/*	if(breakfor ==1 ){
		if(w1 + w_xor_block[2] >= w_cmp ){break;}  //退出j1
		else
		{breakfor = 0; } //清楚breakfor
		} //退出for(i1=
	*/
	}
/*	if(breakfor ==1 ){
		if(w1 + w_xor_block[2] >= w_cmp ){break;}  //退出w1
		else
		{breakfor = 0;  }//清楚breakfor
		} //退出for(w1=
*/
	}
/*	if(breakfor ==1 ){
		if(i0 + w_xor_block[1] >= w_cmp ){break;} //退出j0
		else
		{breakfor = 0; } //清楚breakfor
		} //退出for(j=
	*/
	}
/*	if(breakfor ==1 ){
		if(i0 + w_xor_block[1] >= w_cmp ){break;} //退出i0
		else
		{breakfor = 0;  } //清楚breakfor
		} //退出for(i0=
	*/
	}
	return updated;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
u32 Alzette_RX_Diff_trail_round_N(u16 search_round, u32 x, u32 y)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u32 Alpha = 0;
	u32 Beta = 0;
	u32 gamma_temp = 0;
	u32 xor_alpha_beta = 0;
	u16 w_xor = 0;
	float w_cmp = 0;
	u32 i = 0,j = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u32 alpha_bloc[8] = {0};
	u32 beta_bloc[8] = {0};
	u32 gamma_bloc[8] = {0};
	u32 AB_block[8] = {0};
	u32 carry[8] ={0};
	u32 carry_tmp[8] ={0};
	u8 w_xor_block[8] = {0};
	u32 i0=0,i1=0,i2=0,i3=0;
	u32 j0=0,j1=0,j2=0,j3=0;
	u32 w1=0,w2=0,w3=0;

	u32 zeta_num = 0;
	float w_cmp_rx = 0;

//	u32 breakfor = 0; //快速退出机制

	Alpha = x;
	Beta = rol_right_32(y,Alzette_rol_a[search_round-1]);
	xor_alpha_beta = (Beta ^ Alpha) & ValueMax_Align ;
	w_xor = HM_weight(xor_alpha_beta );

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[search_round -1];
	w_cmp = w_cmp_rx - Ek_min;
	if (w_xor >= w_cmp)
	{
		return 0;
	}

	/*
//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	P_w[search_round] = w_xor; //异或的xdp差分为概率上届,这里可能有问题*******************
	RX_P_w[search_round] = P_w[search_round] + Ek_min;
	Alzette_Pr_zeta_tmp[search_round] = Ek_min;
	Alzette_zeta_tmp[search_round] = 0;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 概率重量更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%//
	if(RX_P_w[search_round] >= w_cmp_rx) // 不比原来的小，则不更新
	{
		return 0;
	}
	//RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1]
	//											+ RX_P_w[search_round];
	Alzette_X[0][search_round] = Beta ^ Alpha; //这里有大问题*******************
	Alzette_X[1][search_round] = rol_right_32(Alzette_X[0][search_round],
									Alzette_rol_b[search_round-1])^y;
	Alzette_X[0][search_round] ^= Alzette_C_RX[Alzette_C_i];

//%%%%%记录更新次数
	Num_Bn_Update++;  //更新次数记录
	RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
	best = 1;
//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
	for(i=0; i <= search_round; i++ )
	{
		Alzette_RX[0][i] = Alzette_X[0][i];
		Alzette_RX[1][i] = Alzette_X[1][i];

		Alzette_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
		Alzette_Pr_add[i] = P_w[i];   //每轮的模加的XOR差分概率重量
		Alzette_Pr_zeta[i] =  Alzette_Pr_zeta_tmp[i];
		Alzette_Zeta[i] = Alzette_zeta_tmp[i];
	}

//	printf("==>> RX_Bn_w: %f   Num_Bn_Update: %d ", RX_Expected_w, Num_Bn_Update);
	time_Round = clock();
	run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
	printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
	printf("Time: %.2f seconds. ==>>\n", run_time);
//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
	*/


////////////////////////////选择对应64bit分组版本////////////////////////////
	for(i=0;i<nBytes;i++)  //nBytes
	{
		alpha_bloc[i] = (Alpha >> (8*i)) & 0xFF; //8 bit
		beta_bloc[i]  = (Beta >> (8*i))  & 0xFF; //8 bit
		AB_block[i] = ((alpha_bloc[i] << 8) | beta_bloc[i]);
	}
	for(j=1;j<nBytes;j++)  //nBytes
	{
		w_xor_block[j] = (u8)_mm_popcnt_u64((xor_alpha_beta >> (j*8)));
	}
	carry[0] = 0;
	for(j=1;j<nBytes;j++)  //nBytes
	{
		carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
			+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
	}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
/////////////////////
	for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
	{
//		i0 = cDDT_wt_min[carry[0]][AB_block[0]]; ///启发式
		if(i0 + w_xor_block[1] >= w_cmp ){break;}
	for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
	{
//		j0=0; ///启发式
		gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
		carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

	for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
	{
//		i1 = cDDT_wt_min[carry[1]][AB_block[1]]; ///启发式
		w1 = i0 + i1;
		if(w1 + w_xor_block[2] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
	for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
	{
//		j1=0;  ///启发式
		gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
		carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

	for(i2 = cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
	{
//		i2 = cDDT_wt_min[carry[2]][AB_block[2]]; ///启发式
		w2 = w1 + i2;
		if(w2 + w_xor_block[3] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
	for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
	{
//		j2=0; ///启发式
		gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
		carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB

//	for(i3 = MSB_cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= MSB_cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
	{
		i3 = MSB_cDDT_wt_min[carry[3]][AB_block[3]];  ///启发式
		w3 = w2 + i3;
		 //概率满足Matsui条件,找到期望相关性重量的线性路径
		if(w3 >= w_cmp ){continue;}  //break 要比 continue 高效,直接结束.
//		if(w3 >= w_cmp ) {return 0;}///启发式
//		for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
			 	j3=0;
				P_w[search_round] = w3; //i0 + i1 +i2 +i3;

				gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
				gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
						   | (gamma_bloc[1] <<8) | gamma_bloc[0];
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%组合zeta%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//基于构造的第i轮的输入输出XOR差分组合，结合zeta构造RX输出差分
//	for(zeta_num=0; zeta_num < RX_zeta_total; zeta_num++)
	{
		zeta_num=0;  //只考虑最后一轮的CL，CR=（0,0)的情况  ///启发式,最后一轮只考虑Ek_min

		RX_P_w[search_round] = P_w[search_round] + RX_zeta_w[zeta_num];
//		if(RX_P_w[search_round] >= w_cmp_rx) // Pmax_compute(xi) is the ri round min weight.
		{
			//return 0; ///启发式
//			break;//只更新比原来RX_Expected_w小的
		}

//		w_cmp = w_cmp_rx - RX_zeta_w[zeta_num];
		Alzette_Pr_zeta_tmp[search_round] = RX_zeta_w[zeta_num];
		Alzette_zeta_tmp[search_round] = RX_zeta[RX_zeta_i[zeta_num]];

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 概率重量更新 %%%%%%%%%%%%%%%%%%%%%%%%%%%//
		//RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1]
		//											+ RX_P_w[search_round];
		Alzette_X[0][search_round] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]];
		Alzette_X[1][search_round] = rol_right_32(Alzette_X[0][search_round],
										Alzette_rol_b[search_round-1])^y;
		Alzette_X[0][search_round] ^= Alzette_C_RX[Alzette_C_i];
//		Alzette_X[0][search_round] ^= 3;

//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		Num_Bn_Update++;  //更新次数记录
		RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
		w_cmp = w3;
		best = 1;
//		breakfor =1; //快速退出机制
		updated = 1;
//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
	for(i=0; i <= search_round; i++ )
	{
		Alzette_RX[0][i] = Alzette_X[0][i];
		Alzette_RX[1][i] = Alzette_X[1][i];

		Alzette_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
		Alzette_Pr_add[i] = P_w[i];   //每轮的模加的XOR差分概率重量
		Alzette_Pr_zeta[i] =  Alzette_Pr_zeta_tmp[i];
		Alzette_Zeta[i] =  Alzette_zeta_tmp[i];
	}
//%%%%%打印更新次数等信息

//	printf("==>> RX_Bn_w: %f   Num_Bn_Update: %d ", RX_Expected_w, Num_Bn_Update);
	time_Round = clock();
	run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
	printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
	printf("Time: %.2f seconds. ==>>\n", run_time);

	print_Alzette_RX_diff_resoult(search_round);  //打印
		}}}
	}
/*	if(breakfor ==1 ){
		if(w2 + w_xor_block[3] >= w_cmp ){break;}  //退出for(i2 =
	else
	{
		breakfor = 0;  //清楚breakfor，继续执行j2中的gamma
		}
	}
	*/
	}
/*	if(breakfor ==1 ){
		if(w1 + w_xor_block[2] >= w_cmp ){break;} //退出for(j1=0;
	else
	{
		breakfor = 0;  //清楚breakfor
		}
	}
	*/
	}
/*	if(breakfor ==1 ){
		if(w1 + w_xor_block[2] >= w_cmp ){break;} //退出for(i1=0;
	else
	{
		breakfor = 0;  //清楚breakfor，
		}
	}
	*/
	}
/*	if(breakfor ==1 ){
		if(i0 + w_xor_block[1] >= w_cmp ){break;}//退出for(j0=0;
	else
	{
		breakfor = 0;  //清楚breakfor，
		}
	}
	*/
	}
/*	if(breakfor ==1 ){
		if(i0 + w_xor_block[1] >= w_cmp ){break;}//退出for(i0=0;
	else
	{
		breakfor = 0;  //清楚breakfor，
		}
	}
	*/
	}
////////////////////////////
return  updated;
}






void print_Alzette_RX_diff_resoult(u16 search_round)
{
	u16 r_num = 0;

printf("|-----------------------------------------------------------------------|\n");
printf("round------x----------y----------zeta-----Pr_zeta-------Pr_xor-----Pr_RX \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     %08llx    %08llx   %08llx   %f    %f    %f \n",
				r_num-1,
				Alzette_RX[0][r_num-1], Alzette_RX[1][r_num-1],Alzette_Zeta[r_num],
				Alzette_Pr_zeta[r_num], Alzette_Pr_add[r_num], Alzette_Pr_RX[r_num]);
	}
	printf("%02d     %08llx    %08llx      NULL        NULL         NULL      NULL\n",
			search_round,
			Alzette_RX[0][search_round], Alzette_RX[1][search_round]);

	printf("%d Round Total Weight: -%f \n",search_round, RX_Expected_w);
	printf("|-----------------------------------------------------------------------|\n");
}



