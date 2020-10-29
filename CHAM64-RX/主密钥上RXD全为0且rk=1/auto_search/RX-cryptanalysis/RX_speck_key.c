/*
 * RX_speck_key.c
 *
 *  Created on: 2020年9月9日
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
#include "RX_speck_key.h"
#include "speck.h"
#include "CHAM.h"  //用到固定一个输入差分，查询另外一个输入差分和输出差分

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
volatile u64 Speck_key_const_rxd[20]={0}; //每轮的密钥编排部分的轮常数的RXD
volatile u64 speck_key_RXD_tmp[4][20]={0}; //密钥部分,左右，每轮的输入差分左右两部分RXD，实际从第2轮开始
volatile u64 speck_key_RXD[4][20]={0}; //密钥部分,左右，每轮的输入差分左右两部分RXD，实际从第2轮开始

volatile u64 speck_Zeta[16]={0}; //每轮的输zeta
volatile u64 speck_zeta_tmp[16]={0}; //每轮的临时zeta
volatile float speck_Pr_RX[16]={0}; //每轮RX差分概率
volatile float speck_Pr_XOR[16]={0}; //每轮RX差分概率
volatile float speck_Pr_zeta[16]={0}; //每轮RX差分概率
volatile float speck_Pr_zeta_tmp[16]={0}; //每轮RX差分概率


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Compute_RX_speck_key_const_rxd(void)  //不用考虑分组长度，轮常数的RXD都不会影响
{
	u16 i = 0;
	for(i=0; i<20; i++)
	{
		Speck_key_const_rxd[i+1] = i ^ (ROTATE_LEFT(i, RX_k, 64));
	}
/*
#if ((speck_m_bits  == 64))  ////speck32:  64     m=4

#elif (speck_m_bits  == 72)  //speck48:  72,96  m=3/4
#if ((speck_m_keyword  == 3))//speck48/72  m=3

#elif ((speck_m_keyword  == 4))//speck48/96  m=4

#endif

#elif (speck_m_bits  == 4)  //speck2n/4n: 4

#endif
*/
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
u16 RX_sepck_key_trail_search_entry (u16 search_round)
{
	FILE* RX_speck_key_Bn;
	u16 i = 0;


	Compute_RX_speck_key_const_rxd();
	//%%%%%%%%%%%%%%%%%%%%%%%预计算XOR-offset表相关信息%%%%%%%%%%%%%%%%%%%%%%%%%//
	RX_zeta_total = 1 + blocksize_len + RX_k*(blocksize_len -RX_k);  //计算Zeta的总个数
	printf("RX_zeta_total: %d  \n",RX_zeta_total);
	printf("Constructing RX_Ek[4] ... \n");
	Compute_RX_Ek(RX_k, (u16)blocksize_len,  Ek); //计算ek的4个值
	Ek_min = Compute_RX_Ek_min(Ek);

	printf("Constructing RX_XOR_offset_Table ... \n");
	Compute_RX_V_U(RX_k, blocksize_len, RX_v,RX_v_w, RX_u,RX_u_w);
	Sorted_RX_XOR_offset_Table(RX_k, (u16)blocksize_len,
			RX_v,RX_v_w, RX_u,RX_u_w,
			Ek,
			RX_zeta, RX_zeta_w, RX_zeta_i);
	//-------------取前r-1轮最优概率重量，搜索期望概率重量的最大值----------------//
	RX_speck_key_Bn = fopen ("../tmp/RX_speck_key_Bn_wt.xlsx", "a+"); //  "w+"); //
	if(search_round > 1)
	{
		for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
		{
			fscanf(RX_speck_key_Bn,"%f",&RX_P_bestofR_w[i]);  //read.
		}
	}
	fclose(speck_best_Bn);
		//-------记录时间-----
		time_xor_talbe = clock();
		run_time =  (time_xor_talbe - time_start) / CLOCKS_PER_SEC;
		printf("Time of Preprocessing info: %.2f seconds.  \n", run_time);
		printf("|--------------------------------------------------|\n");
		//%%%%%%%%%%%%%%%%%%%%%%%%%%构造cDDT、XOR-offset表%%%%%%%%%%%%%%%%%%//
		printf("Constructing cDDT Lookup Tables ... \n");
		//构造cDDT
#if (speck_m_keyword  == 2)  //speck2n/2n: 2
	bit_align_fun();
	ARX_carry_DDTm_construct();
	//fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
#elif (speck_m_keyword  == 3)  //speck2n/3n: 3
	bit_align_fun();
	ARX_carry_DDTm_construct();
	fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
#elif (speck_m_keyword  == 4)  //speck2n/4n: 4
	bit_align_fun();
	ARX_carry_DDTm_construct();
	fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
#endif
	time_ARX_DDT = clock();
	run_time =  (time_ARX_DDT - time_start) / CLOCKS_PER_SEC;
	printf("Time of Construct ARX_DDT: %.2f seconds.  \n", run_time);

printf("*****************************Searching processing**************************|\n");
//%%%%%%%%%%%%%%%%%%%%%%%%%%搜索程序入口%%%%%%%%%%%%%%%%%%//
//	Bn_w = 47;  //直接设定开始的概率重量,快速达到期望的概率重量
	do
	{
		//Bn_w = Bn_w + 1;
		RX_Expected_w = speck_m_bits;  //密钥编排部分最大概率重量为密钥长度
		printf("Update: %d   RX_Bn_w: %f  --Max expected weight-- \n",
				Num_Bn_Update, RX_Expected_w);

#if (speck_m_keyword  == 2)  //speck2n/2n: 2
		// Search Entry of speck2n/2n.
		RX_sepck_2n2n_key_round_1(search_round );
#elif (speck_m_keyword  == 3)  //speck2n/3n: 3
		// Search Entry of speck2n/2n.
		RX_sepck_2n3n_key_round_1(search_round );
#elif (speck_m_keyword  == 4)  //speck2n/4n: 4
		// Search Entry of speck2n/2n.
		RX_sepck_2n3n_key_round_1(search_round );
#endif

		time_Round = clock();
		run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
		printf("Time: %.2f seconds.  \n", run_time);
	}while( Num_Bn_Update == 0 );   //end conditions. // 概率重量逼近次数大于1，且返回值表示搜索到最优路径
	printf("**************************************************************************|\n");

	fprintf(speck_best_Bn,"%d \n",n_P_bestofR_w[search_round]);  // write.

	for(i=0; i <= search_round; i++)
	{
		n_P_w[i] = P_w[i];
	}
	print_resoult(search_round);

	time_finish = clock();
	run_time =  (double)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
	printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours. \n", run_time,run_time/60.0,run_time/3600.0 );
	printf("Auto-search_speck END! \n");
	printf("|************************************************************************|\n");
	return 1;
}

/////////
#if (speck_m_keyword  == 2)  //speck2n/2n: 2
/////////
u16 RX_sepck_2n2n_key_round_1(u16 search_round)
{
	u16 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 thr_d = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;

	u64 M0 = 0;
	u64 N0 = blocksize_len-1;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;

	u32 zeta_num = 0;
	u32 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;

///////////////////round 1 entry/////////////////////
	// 可以选择是否限定第1轮的概率重量
	//printf("**Noted**: for(thr_d = 0; thr_d < ?; thr_d++) \n"); //打印启发信息
	for(thr_d = 0; thr_d < blocksize_len-1 ;thr_d++)  //0::n-1 //blocksize_len-1
	{
		// thr_d : the bits number of alpha/beta/gamma that the position not equal at the same time.
		if ((thr_d +Ek_min+ RX_P_bestofR_w[search_round - 1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
		{  //统一退出位置
			return updated;//Return the upper procedure
		}

		M0 = thr_d;  //三个alpha,beta和gamma同时考虑
		P_w[1] = M0;  //pro; //Speck_XDP_compute(Input_Alpha,Input_Beta,Input_Gamma );
		if(M0 == 0)  /// abd [0 -> n-2] 没有完全不相同的比特,且MSB之后只能为0
		{
			for(j_abc = 1; j_abc < 4;j_abc++)
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

				//主密钥,左侧
				speck_key_RXD_tmp[0][0] = Input_beta; //主密钥,右侧k0
				speck_key_RXD_tmp[1][0] = ROTATE_LEFT(Input_alpha,rol_a,blocksize_len) & Bit_Align;//主密钥,左侧k1

				speck_key_RXD_tmp[0][1] = Input_beta;  //第1轮输出，实际第2轮的key的RXD
				speck_key_RXD_tmp[1][1] = speck_key_RXD_tmp[1][0]; //第1轮输出，实际第2轮的key的RXD


				for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
					if ((RX_P_w[1] + RX_P_bestofR_w[search_round - 1])
							>= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
					{
						break;
					}
					RX_P_sumofR_w[1] = RX_P_w[1];
					speck_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
					speck_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

					speck_key_RXD_tmp[1][2] = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD
					speck_key_RXD_tmp[0][2] = speck_key_RXD_tmp[1][2] ^
						(ROTATE_LEFT(Input_beta,rol_b,blocksize_len) & Bit_Align);  //第2轮输出RXD

				best = RX_sepck_2n2n_key_round_r(search_round, 2,
						speck_key_RXD_tmp[1][2],speck_key_RXD_tmp[0][2]);
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
	}}
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
				best = RX_sepck_2n2n_key_input_MSB(search_round,
						Input_alpha, Input_beta, Input_gamma,thr_d,C0);
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
		best = RX_sepck_2n2n_key_input_MSB(search_round,Input_alpha, Input_beta, Input_gamma,thr_d,C0);
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



u16 RX_sepck_2n2n_key_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
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

			state = RX_sepck_2n2n_key_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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

		state = RX_sepck_2n2n_key_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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
	return state;
}

//由高位到低位逐步判断
u16 RX_sepck_2n2n_key_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi)
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	u16 i = 0;
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


			state = RX_sepck_2n2n_key_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
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
		state = RX_sepck_2n2n_key_input_Last(search_round,
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
	return state;
}

u16 RX_sepck_2n2n_key_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;

	u64 bit_i = 1;

	u32 zeta_num = 0;
	u32 updata_cur_num = 0;


	indx_tmp = posi[1] - 1;

	// 最后的一个产生概率重量的比特位置只能为{3、5/6}，在该位置之后的更低的比特位置都全部清0；
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


		//主密钥,左侧
		speck_key_RXD_tmp[0][0] = Input_beta; //主密钥,右侧
		speck_key_RXD_tmp[1][0] = ROTATE_LEFT(Input_alpha,rol_a,blocksize_len) & Bit_Align;//主密钥,左侧

		speck_key_RXD_tmp[0][1] = Input_beta;  //第1轮输出，实际第2轮的key的RXD
		speck_key_RXD_tmp[1][1] = speck_key_RXD_tmp[1][0]; //第1轮输出，实际第2轮的key的RXD

		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
			if ((RX_P_w[1] + RX_P_bestofR_w[search_round - 1])
					>= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				break;
			}
			RX_P_sumofR_w[1] = RX_P_w[1];
			speck_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
			speck_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

			speck_key_RXD_tmp[1][2] = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD
			speck_key_RXD_tmp[0][2] = speck_key_RXD_tmp[1][2] ^
				(ROTATE_LEFT(Input_beta,rol_b,blocksize_len) & Bit_Align);  //第2轮输出RXD

		state = RX_sepck_2n2n_key_round_r(search_round, 2,
				speck_key_RXD_tmp[1][2],speck_key_RXD_tmp[0][2]);
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


u16 RX_sepck_2n2n_key_round_r(u16 search_round, u16 cur_round, u64 x, u64 y)
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 xor_alpha_beta = 0;
//	u64 and_alpha_beta = 0;
	u16 w_xor = 0;
	float w_cmp = 0;
	u16 i = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
//	u64 rol_left = 0;
//	u64 rol_right = 0;
	u64 alpha_bloc[8] = {0};
//	u64 gamma_block_tmp[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 j = 0;
	u64 beta_temp = 0;
	u64 gamma_temp = 0;
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};
//	u16 wt_purning = 63;
//	u16 wt_purning_chk = 63;


	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制
	updata_cur_num = Num_Bn_Update;

	if(search_round == cur_round)
	{
		state = RX_speck_2n2n_key_round_N(search_round, x, y );
		return state;
	}

	Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align;
	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;

//	rol_right = ((x) >> (rol_a)) | ((x) << (blocksize_len - rol_a)); //right
//	rol_left = (((y) << rol_b) | ((y) >> (blocksize_len - rol_b)));
//	Beta = (x ^ rol_left) & Bit_Align;
//	Alpha = rol_right & Bit_Align; // alpha = gama_uper <<< rol_a;
	//Beta = (x ^ ROTATE_LEFT(y,rol_b,blocksize_len) )& Bit_Align;
	//Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align; // alpha = gama_uper <<< rol_a;
	xor_alpha_beta = (Beta ^ Alpha) & ValueMax_Align ;
	w_xor = HM_weight(xor_alpha_beta );
	//w_xor = Speck_XDP_Max(Alpha, Beta);

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
				- RX_P_bestofR_w[search_round - cur_round];
	w_cmp = w_cmp_rx - Ek_min;
	if ( w_xor >= w_cmp)
	{
		return 0;
	}

////////////////////////////选择SPECK的对应密钥长度版本////////////////////////////
#if (speck_m_bits  == 96)  //SPECK-96/96
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
				i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]]
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
				i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
				carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
				i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  // cDDT_wt_max[carry[2]][AB_block[2]]
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
				i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
				gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
				carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
				i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++)  //cDDT_wt_max[carry[4]][AB_block[4]]
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5]> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
			{
				gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
				carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]];
				i5 <= MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //MSB_cDDT_wt_max[carry[5]][AB_block[5]]
		{
			w5 = w4 + i5;
			if(w5 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			{
				P_w[cur_round] = w5;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
				RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
				if (RX_P_w[cur_round] >= w_cmp_rx){break;} // Pmax_compute(xi) is the ri round min weight.
				speck_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
				speck_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
				RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

			for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
			{
			gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];
			gamma_temp = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];

			speck_key_RXD_tmp[1][cur_round+1] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]]; //下轮输出RXD
			speck_key_RXD_tmp[1][cur_round+1] ^= Speck_key_const_rxd[cur_round]; //轮常数
			speck_key_RXD_tmp[0][cur_round+1] = speck_key_RXD_tmp[1][cur_round+1] ^next_Beta_tmp;  //下一轮输出RXD

			state = RX_sepck_2n2n_key_round_r(search_round, cur_round+1,
					speck_key_RXD_tmp[1][cur_round],speck_key_RXD_tmp[0][cur_round]);
	#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
		if(state == 1 )  //判断本轮是否要更新RX_Expected_w
		{
			if(updata_cur_num < Num_Bn_Update)
			{
				updata_cur_num = Num_Bn_Update;
				updated = 1;  //准备通知上一轮最新的RX_Expected_w

			//更新了期望概率重量RX_Expected_w//先考虑是否直接退出本次的PW
			if((P_w[cur_round] + RX_zeta_w[0] +	RX_P_bestofR_w[search_round - cur_round]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				updated = 2;  //准备通知上一轮，准备结束搜索
				return updated;  //这里是跳出for(j_abc = 0;
			}
			////更新了期望概率重量RX_Expected_w，考虑是否搜剩下的组合%%%%%%%%%%%%%%%%%%%%%%%
			if((RX_P_w[cur_round] +	RX_P_bestofR_w[search_round - cur_round]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				break;   //这里是退出for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			}
		}}
	#endif
			}
		}
		}}}}}}}}}}}
////////////////////////////
#elif (speck_m_bits  == 128)  //SPECK-48

		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
			//w_xor_block[j] = cDDT_AB_wt_min[AB_block[j]];
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

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
				i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]];
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//gamma_block_tmp[0] = gamma_bloc[0];
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
				i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
				i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  //cDDT_wt_max[carry[2]][AB_block[2]]
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
				i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break
		for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
			carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
				i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
		{
			gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
			carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = cDDT_wt_min[carry[5]][AB_block[5]];
				i5 <= cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
		{
			w5 = w4 + i5;
			if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
		{
			gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
			carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7); // gamma MSB
		for(i6 = cDDT_wt_min[carry[6]][AB_block[6]];
				i6 <= cDDT_wt_min[carry[6]][AB_block[6]]; i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
		{
			w6 = w5 + i6;
			if(w6 + w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
		{
			gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];
			carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7); // gamma MSB
		for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]];
				i7 <= MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
		{
			w7 = w6 + i7;
			if(w7> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			{
				P_w[cur_round] = w7;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
				RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
				if (RX_P_w[cur_round] >= w_cmp_rx){break;} // Pmax_compute(xi) is the ri round min weight.
				speck_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
				speck_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
				RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

				for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
				{
					gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];
					//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
					gamma_temp = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
							| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
							| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
							| (gamma_bloc[1] <<8) | gamma_bloc[0];

					speck_key_RXD_tmp[1][cur_round+1] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]]; //下轮输出RXD
					speck_key_RXD_tmp[1][cur_round+1] ^= Speck_key_const_rxd[cur_round]; //轮常数
					speck_key_RXD_tmp[0][cur_round+1] = speck_key_RXD_tmp[1][cur_round+1] ^next_Beta_tmp;  //下一轮输出RXD

			state = RX_sepck_2n2n_key_round_r(search_round, cur_round+1,
					speck_key_RXD_tmp[1][cur_round],speck_key_RXD_tmp[0][cur_round]);
	#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
		if(state == 1 )  //判断本轮是否要更新RX_Expected_w
		{
			if(updata_cur_num < Num_Bn_Update)
			{
				updata_cur_num = Num_Bn_Update;
				updated = 1;  //准备通知上一轮最新的RX_Expected_w

			//更新了期望概率重量RX_Expected_w//先考虑是否直接退出本次的PW
			if((P_w[cur_round] + RX_zeta_w[0] +
					RX_P_bestofR_w[search_round - cur_round]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				updated = 2;  //准备通知上一轮，准备结束搜索
				return updated;  //这里是跳出for(j_abc = 0;
			}
			////更新了期望概率重量RX_Expected_w，考虑是否搜剩下的组合%%%%%%%%%%%%%%%%%%%%%%%
			if((RX_P_w[cur_round] +	RX_P_bestofR_w[search_round - cur_round])
					>= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				break;   //这里是退出for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			}
		}}
	#endif
			}
		}
		}}}}}}}}}}}}}}}
#endif
	return state;
}



u16 RX_speck_2n2n_key_round_N(u16 search_round, u64 x, u64 y)
{
	u16 best = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 xor_alpha_beta = 0;
	u16 w_xor = 0;
	float w_cmp = 0;
	u16 i = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 j = 0;
	u64 gamma_temp = 0;
//	u64 gamma_block_tmp[8] = {0};
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};

	u32 zeta_num = 0;
	float w_cmp_rx = 0;



	Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align;
	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;


//	rol_right = ((x) >> (rol_a)) | ((x) << (blocksize_len - rol_a)); //right
//	rol_left = (((y) << rol_b) | ((y) >> (blocksize_len-rol_b)));
//	Beta = (x ^ rol_left) & Bit_Align;
//	Alpha = rol_right & Bit_Align; // alpha = gama_uper <<< rol_a;
	//Beta = x ^ (ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align);
	//Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align; // alpha = gama_uper <<< rol_a;
	xor_alpha_beta = (Beta ^ Alpha)  & ValueMax_Align;
	w_xor = HM_weight(xor_alpha_beta);

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[search_round -1];
	w_cmp = w_cmp_rx - Ek_min;
	if (w_xor >= w_cmp)
	{
		return 0;
	}

	////////////////////////////选择SPECK的对应密钥长度版本////////////////////////////
#if (speck_m_bits  == 96)  //SPECK-96/96
	for(j=0;j<nBytes;j++)  //nBytes
	{
		alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
		beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
		AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

			for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
					i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]]
			{
				//w0 = i0;
				if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
				{
					gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
					carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
			for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
					i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
			{
				w1 = i0 + i1;
				if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
				{
					gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
					carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
			for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
					i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  // cDDT_wt_max[carry[2]][AB_block[2]]
			{
				w2 = w1 + i2;
				if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
				{
				gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
				carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
			for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
					i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
			{
				w3 = w2 + i3;
				if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
				{
					gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
					carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
			for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
					i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++)  //cDDT_wt_max[carry[4]][AB_block[4]]
			{
				w4 = w3 + i4;
				if(w4 + w_xor_block[5]> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
				{
					gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
					carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
			for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]];
					i5 <= MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //MSB_cDDT_wt_max[carry[5]][AB_block[5]]
			{
				w5 = w4 + i5;
				if(w5 >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				//for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					zeta_num=0;
					P_w[search_round] = w5;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
					RX_P_w[search_round] = P_w[search_round] + RX_zeta_w[zeta_num];
					if (RX_P_w[search_round] >= w_cmp_rx){break;} // Pmax_compute(xi) is the ri round min weight.
					speck_Pr_zeta_tmp[search_round] = RX_zeta_w[zeta_num];
					speck_zeta_tmp[search_round] = RX_zeta[RX_zeta_i[zeta_num]];
					RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1] + RX_P_w[search_round];

				//for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
				{
						j5=0;
				gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];
				gamma_temp = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
						| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
						| (gamma_bloc[1] <<8) | gamma_bloc[0];

				speck_key_RXD[1][search_round+1] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]]; //下轮输出RXD
				speck_key_RXD[1][search_round+1] ^= Speck_key_const_rxd[search_round]; //轮常数
				speck_key_RXD[0][search_round+1] = speck_key_RXD[1][search_round+1] ^next_Beta_tmp;  //下一轮输出RXD

				//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
						Num_Bn_Update++;  //更新次数记录
						RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
						w_cmp = w5;
						best = 1;
				//		breakfor =1; //快速退出机制
						updated = 1;
				//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
					for(i=0; i <= search_round; i++ )
					{
						speck_key_RXD[0][i] = speck_key_RXD_tmp[0][i];
						speck_key_RXD[1][i] = speck_key_RXD_tmp[1][i];

						speck_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
						speck_Pr_XOR[i] = P_w[i];   //每轮的模加的XOR差分概率重量
						speck_Pr_zeta[i] =  speck_Pr_zeta_tmp[i];
						speck_Zeta[i] =  speck_zeta_tmp[i];
					}

					time_Round = clock();
					run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
					printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
					printf("Time: %.2f seconds. ==>>\n", run_time);

					print_RX_speck_key_resoult(search_round);  //打印
				}
			}
			}}}}}}}}}}}
	////////////////////////////
#elif (speck_m_bits  == 128)  //SPECK-128/128
			for(j=0;j<nBytes;j++)  //nBytes
			{
				alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
				beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
				AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
				//w_xor_block[j] = cDDT_AB_wt_min[AB_block[j]];
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

			for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
					i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]];
			{
				//w0 = i0;
				if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				//gamma_block_tmp[0] = gamma_bloc[0];
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
			for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
					i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
			{
				w1 = i0 + i1;
				if(w1 + w_xor_block[2] > w_cmp ){break;}  //break
			for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
				carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
			for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
					i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  //cDDT_wt_max[carry[2]][AB_block[2]]
			{
				w2 = w1 + i2;
				if(w2 + w_xor_block[3] > w_cmp ){break;}  //break
			for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
				gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
				carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
			for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
					i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
			{
				w3 = w2 + i3;
				if(w3 + w_xor_block[4] > w_cmp ){break;}  //break
			for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
				gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
				carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
			for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
					i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
			{
				w4 = w3 + i4;
				if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
			{
				gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
				carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
			for(i5 = cDDT_wt_min[carry[5]][AB_block[5]];
					i5 <= cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
			{
				w5 = w4 + i5;
				if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
			{
				gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
				carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7); // gamma MSB
			for(i6 = cDDT_wt_min[carry[6]][AB_block[6]];
					i6 <= cDDT_wt_min[carry[6]][AB_block[6]]; i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
			{
				w6 = w5 + i6;
				if(w6 + w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
			for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
			{
				gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];
				carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7); // gamma MSB
			for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]];
					i7 <= MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
			{
				w7 = w6 + i7;
				if(w7>= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				//for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					zeta_num=0;
					P_w[search_round] = w7;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
					RX_P_w[search_round] = P_w[search_round] + RX_zeta_w[zeta_num];
					if (RX_P_w[search_round] >= w_cmp_rx){break;} // Pmax_compute(xi) is the ri round min weight.
					speck_Pr_zeta_tmp[search_round] = RX_zeta_w[zeta_num];
					speck_zeta_tmp[search_round] = RX_zeta[RX_zeta_i[zeta_num]];
					RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1] + RX_P_w[search_round];

				//for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
				{
						j7=0;
						gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];
						//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
						gamma_temp = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
								| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
								| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
								| (gamma_bloc[1] <<8) | gamma_bloc[0];

				speck_key_RXD[1][search_round+1] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]]; //下轮输出RXD
				speck_key_RXD[1][search_round+1] ^= Speck_key_const_rxd[search_round]; //轮常数
				speck_key_RXD[0][search_round+1] = speck_key_RXD[1][search_round+1] ^next_Beta_tmp;  //下一轮输出RXD

				//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
						Num_Bn_Update++;  //更新次数记录
						RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
						w_cmp = w7;
						best = 1;
				//		breakfor =1; //快速退出机制
						updated = 1;
				//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
					for(i=0; i <= search_round; i++ )
					{
						speck_key_RXD[0][i] = speck_key_RXD_tmp[0][i];
						speck_key_RXD[1][i] = speck_key_RXD_tmp[1][i];

						speck_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
						speck_Pr_XOR[i] = P_w[i];   //每轮的模加的XOR差分概率重量
						speck_Pr_zeta[i] =  speck_Pr_zeta_tmp[i];
						speck_Zeta[i] =  speck_zeta_tmp[i];
					}

					time_Round = clock();
					run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
					printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
					printf("Time: %.2f seconds. ==>>\n", run_time);

					print_RX_speck_key_resoult(search_round);  //打印
				}
			}
			}}}}}}}}}}}}}}}
	#endif
return  updated;
}
/////////
#elif (speck_m_keyword  == 3)  //speck2n/3n: 3
/////////
u16 RX_sepck_2n3n_key_round_1(u16 search_round)
{
	u16 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 thr_d = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;

	u64 M0 = 0;
	u64 N0 = blocksize_len-1;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;

	u32 zeta_num = 0;
	u32 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;

///////////////////round 1 entry/////////////////////
	// 可以选择是否限定第1轮的概率重量
	//printf("**Noted**: for(thr_d = 0; thr_d < ?; thr_d++) \n"); //打印启发信息
	for(thr_d = 0; thr_d < blocksize_len-1 ;thr_d++)  //0::n-1 //blocksize_len-1
	{
		// thr_d : the bits number of alpha/beta/gamma that the position not equal at the same time.
		if ((thr_d +Ek_min+ RX_P_bestofR_w[search_round - 1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
		{  //统一退出位置
			return updated;//Return the upper procedure
		}

		M0 = thr_d;  //三个alpha,beta和gamma同时考虑
		P_w[1] = M0;  //pro; //Speck_XDP_compute(Input_Alpha,Input_Beta,Input_Gamma );
		if(M0 == 0)  /// abd [0 -> n-2] 没有完全不相同的比特,且MSB之后只能为0
		{
			for(j_abc = 1; j_abc < 4;j_abc++)
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

				speck_key_RXD_tmp[0][0] = Input_beta; //主密钥,右侧k0
				speck_key_RXD_tmp[1][0] = ROTATE_LEFT(Input_alpha,rol_a,blocksize_len) & Bit_Align;//主密钥,k1
				//speck_key_RXD_tmp[2][0] = ?  //k2

				speck_key_RXD_tmp[0][1] = speck_key_RXD_tmp[0][0]; //第1轮输出，实际第2轮的key的RXD
				speck_key_RXD_tmp[1][1] = speck_key_RXD_tmp[1][0];  //第1轮输出，实际第2轮的key的RXD
				//speck_key_RXD_tmp[2][1] = ?  //k2

				for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
					if ((RX_P_w[1] + RX_P_bestofR_w[search_round - 1])
							>= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
					{
						break;
					}
					RX_P_sumofR_w[1] = RX_P_w[1];
					speck_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
					speck_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

					speck_key_RXD_tmp[2][2] = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD
					speck_key_RXD_tmp[0][2] = speck_key_RXD_tmp[2][2] ^
						(ROTATE_LEFT(Input_beta,rol_b,blocksize_len) & Bit_Align);  //第2轮输出RXD

				best = RX_sepck_2n3n_key_round_2(search_round, 2,
						speck_key_RXD_tmp[2][2],speck_key_RXD_tmp[0][2]);
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
	}}
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
				best = RX_sepck_2n3n_key_input_MSB(search_round,
						Input_alpha, Input_beta, Input_gamma,thr_d,C0);
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
		best = RX_sepck_2n3n_key_input_MSB(search_round,Input_alpha, Input_beta, Input_gamma,thr_d,C0);
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




u16 RX_sepck_2n3n_key_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
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

			state = RX_sepck_2n3n_key_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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

		state = RX_sepck_2n3n_key_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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
	return state;
}

//由高位到低位逐步判断
u16 RX_sepck_2n3n_key_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi)
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	u16 i = 0;
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


			state = RX_sepck_2n3n_key_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
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
		state = RX_sepck_2n3n_key_input_Last(search_round,
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
	return state;
}

u16 RX_sepck_2n3n_key_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;

	u64 bit_i = 1;

	u32 zeta_num = 0;
	u32 updata_cur_num = 0;


	indx_tmp = posi[1] - 1;

	// 最后的一个产生概率重量的比特位置只能为{3、5/6}，在该位置之后的更低的比特位置都全部清0；
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

		speck_key_RXD_tmp[0][0] = Input_beta; //主密钥,右侧k0
		speck_key_RXD_tmp[1][0] = ROTATE_LEFT(Input_alpha,rol_a,blocksize_len) & Bit_Align;//主密钥,k1
		//speck_key_RXD_tmp[2][0] = ?  //k2

		speck_key_RXD_tmp[0][1] = speck_key_RXD_tmp[0][0]; //第1轮输出，实际第2轮的key的RXD
		speck_key_RXD_tmp[1][1] = speck_key_RXD_tmp[1][0];  //第1轮输出，实际第2轮的key的RXD
		//speck_key_RXD_tmp[2][1] = ?  //k2

		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
			if ((RX_P_w[1] + RX_P_bestofR_w[search_round - 1])
					>= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
			{
				break;
			}
			RX_P_sumofR_w[1] = RX_P_w[1];
			speck_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
			speck_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

			speck_key_RXD_tmp[2][2] = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD
			speck_key_RXD_tmp[0][2] = speck_key_RXD_tmp[2][2] ^
				(ROTATE_LEFT(Input_beta,rol_b,blocksize_len) & Bit_Align);  //第2轮输出RXD

		state = RX_sepck_2n3n_key_round_2(search_round, 2,
				speck_key_RXD_tmp[2][2],speck_key_RXD_tmp[0][2]);
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




u16 RX_sepck_2n3n_key_round_2(u16 search_round, u16 cur_round, u64 x, u64 y)
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 K0_bloc[8] = {0};
	u8  K0_bloc_wt_min[8] = {0};
	u64 K0_bloc_msb[8] = {0};
	u64 i = 0, j = 0;
	u64 alpha_bloc[8] = {0};
	u64 gama_block[8] = {0};
	u8 carry_block[4] = {0};
	u64 alpha_tmp = 0;
	u64 gamma_tmp = 0;

	u64 gamma_block_tmp[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;


	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};


	float w_cmp = 0;
	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制

	updata_cur_num = Num_Bn_Update;

	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;
	speck_key_RXD_tmp[1][3] = speck_key_RXD_tmp[2][2]; //第2轮输出RXD,k1

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
				- RX_P_bestofR_w[search_round - cur_round];
	w_cmp = w_cmp_rx - Ek_min;


	////////////////////////////
#if (speck_m_bits  == 72)  //SPECK-48/72
	K0_bloc[0] = Beta & 0xFF;   //取出每个8比特
	K0_bloc[1] = (Beta >> 8) & 0xFF;   //取出每个8比特
	K0_bloc[2] = (Beta >> 16) & 0xFF;   //取出每个8比特
	//确定每个输入差分beta的每8比特的子块所对应的可能的最小差分概率重量
	K0_bloc_wt_min[0] = CHAM_a_wt_min[K0_bloc[0]][0];
	K0_bloc_wt_min[1] = 0xFF;
	K0_bloc_wt_min[2] = 0xFF;
	for(j=0;j<2;j++)
	{
		if( (K0_bloc[j] & 0x80) == 0)
		{
			for(i=0;i<4;i++) //carry 为 0/1/2/3
			{
				if(K0_bloc_wt_min[j+1] > CHAM_a_wt_min[K0_bloc[j+1]][i])
				{
					K0_bloc_wt_min[j+1] = CHAM_a_wt_min[K0_bloc[j+1]][i];
				}
			}
		}
		else
		{
			for(i=4;i<8;i++)  //carry 为 4/5/6/7
			{
				if(K0_bloc_wt_min[j+1] > CHAM_a_wt_min[K0_bloc[j+1]][i])
				{
					K0_bloc_wt_min[j+1] = CHAM_a_wt_min[K0_bloc[j+1]][i];
				}
			}
		}
	}

	if(K0_bloc_wt_min[0] + K0_bloc_wt_min[1]+ K0_bloc_wt_min[2] > w_cmp)
	{return 0;}

	//注意： alpha和Beta是可以互换的
	K0_bloc_msb[0] = ((K0_bloc[0] >> 7) << 2);  //alpha 最低8比特的 MSB
	K0_bloc_msb[1] = ((K0_bloc[1] >> 7) << 2);  //alpha 8-16比特的 MSB
//	K0_bloc_msb[2] = ((K0_bloc[2] >> 7) << 2);  //alpha 16-24比特的 MSB

	//block 0  //可以控制每个block的概率重量范围, 例如[min, min +1]
	for(i0=CHAM_a_wt_min[K0_bloc[0]][0];
			i0<=CHAM_a_wt_min[K0_bloc[0]][0] +1; i0++)  //X0_block_wt_min[0] // <=8 //CHAM_a_wt_max[X0_block[0]][0]
	{
		if(i0 + K0_bloc_wt_min[1] + K0_bloc_wt_min[2] > w_cmp ){break;}
		for(j0=0; j0<CHAM_a_bg_numb[K0_bloc[0]][0][i0]; j0++)
		{
			alpha_bloc[0] = CHAM_a_beta[K0_bloc[0]][0][i0][j0];
			gama_block[0] = CHAM_a_gama[K0_bloc[0]][0][i0][j0];
			carry_block[0] = K0_bloc_msb[0] +
					((alpha_bloc[0] >> 7) << 1) + (gama_block[0] >> 7);

	//block 1
	for(i1=CHAM_a_wt_min[K0_bloc[1]][carry_block[0]];
					i1<=CHAM_a_wt_min[K0_bloc[1]][carry_block[0]] +1; i1++) // CHAM_a_wt_max[X0_block[1]][carry_block[0]]
	{
		if(i0 + i1 + K0_bloc_wt_min[2] > w_cmp ){break;}
		for(j1=0; j1<CHAM_a_bg_numb[K0_bloc[1]][carry_block[0]][i1]; j1++)
		{
			alpha_bloc[1] = CHAM_a_beta[K0_bloc[1]][carry_block[0]][i1][j1];
			gama_block[1] = CHAM_a_gama[K0_bloc[1]][carry_block[0]][i1][j1];
			carry_block[1] = K0_bloc_msb[1] +
					((alpha_bloc[1] >> 7) << 1) + (gama_block[1] >> 7);

	//block 2
	for(i2=msb_CHAM_a_wt_min[K0_bloc[2]][carry_block[1]];
							i2<=msb_CHAM_a_wt_min[K0_bloc[2]][carry_block[1]] +1; i2++)  //CHAM_a_wt_max[X0_block[2]][carry_block[1]]
	{
		P_w[cur_round] = i0 + i1 + i2;  //// 当前轮的概率重量
		if(P_w[cur_round] > w_cmp ){break;}
		//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
			speck_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
			speck_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
			RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
			RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

			for(j2=0; j2<msb_CHAM_a_bg_numb[K0_bloc[2]][carry_block[1]][i2]; j2++)
			{
				alpha_bloc[2] = CHAM_a_beta[K0_bloc[2]][carry_block[1]][i2][j2];
				gama_block[2] = CHAM_a_gama[K0_bloc[2]][carry_block[1]][i2][j2];


				alpha_tmp = ((alpha_bloc[2] << 16) |
							(alpha_bloc[1] << 8) | alpha_bloc[0]) & Bit_Align;
				gamma_tmp = ((gama_block[2] << 16) |
							(gama_block[1] << 8) | gama_block[0]) & Bit_Align;

				//整理前两轮的临时输入输出RXD
				speck_key_RXD_tmp[2][0] = (ROTATE_LEFT(alpha_tmp,rol_a,blocksize_len) & Bit_Align);  //第2轮输出RXD

				speck_key_RXD_tmp[2][1] = speck_key_RXD_tmp[2][0]; //第0/1轮输入RXD,k2

				speck_key_RXD_tmp[1][2] = speck_key_RXD_tmp[2][0]; //第1轮输出RXD,k1

				//speck_key_RXD_tmp[1][3] = speck_key_RXD_tmp[2][2]; //第2轮输出RXD,k1
				speck_key_RXD_tmp[2][3] = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD,k1
				speck_key_RXD_tmp[2][3] ^= Speck_key_const_rxd[cur_round]; //轮常数; //第2轮输出RXD,k1
				speck_key_RXD_tmp[0][3] = speck_key_RXD_tmp[2][3] ^ next_Beta_tmp; //第2轮输出RXD,k1

				state = RX_sepck_2n3n_key_round_r(search_round, 3,
						speck_key_RXD_tmp[1][3],speck_key_RXD_tmp[0][3]);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(state == 1 )  //判断本轮是否要更新RX_Expected_w
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
			//breakfor = 1;
			break;
		}  //这里有问题，只应该退出i3的嵌套
	}}
#endif
			}
		}}}}}}
	////////////////////////////
#elif (speck_m_bits  == 96)  //SPECK-64/96
	K0_bloc[0] = Beta & 0xFF;   //取出每个8比特  4X8=32
	K0_bloc[1] = (Beta >> 8) & 0xFF;   //取出每个8比特
	K0_bloc[2] = (Beta >> 16) & 0xFF;   //取出每个8比特
	K0_bloc[2] = (Beta >> 24) & 0xFF;   //取出每个8比特
	//确定每个输入差分beta的每8比特的子块所对应的可能的最小差分概率重量
	K0_bloc_wt_min[0] = CHAM_a_wt_min[K0_bloc[0]][0];
	K0_bloc_wt_min[1] = 0xFF;
	K0_bloc_wt_min[2] = 0xFF;
	K0_bloc_wt_min[3] = 0xFF;
	for(j=0;j<3;j++)
	{
		if( (K0_bloc[j] & 0x80) == 0)
		{
			for(i=0;i<4;i++) //carry 为 0/1/2/3
			{
				if(K0_bloc_wt_min[j+1] > CHAM_a_wt_min[K0_bloc[j+1]][i])
				{
					K0_bloc_wt_min[j+1] = CHAM_a_wt_min[K0_bloc[j+1]][i];
				}
			}
		}
		else
		{
			for(i=4;i<8;i++)  //carry 为 4/5/6/7
			{
				if(K0_bloc_wt_min[j+1] > CHAM_a_wt_min[K0_bloc[j+1]][i])
				{
					K0_bloc_wt_min[j+1] = CHAM_a_wt_min[K0_bloc[j+1]][i];
				}
			}
		}
	}

	if(K0_bloc_wt_min[0] + K0_bloc_wt_min[1]
		+ K0_bloc_wt_min[2]+ K0_bloc_wt_min[3]
		> w_cmp){return 0;}


	//注意： alpha和Beta是可以互换的
	K0_bloc_msb[0] = ((K0_bloc[0] >> 7) << 2);  //alpha 最低8比特的 MSB
	K0_bloc_msb[1] = ((K0_bloc[1] >> 7) << 2);  //alpha 8-16比特的 MSB
	K0_bloc_msb[2] = ((K0_bloc[2] >> 7) << 2);  //alpha 16-24比特的 MSB
//	K0_bloc_msb[3] = ((K0_bloc[3] >> 7) << 2);  //alpha 16-24比特的 MSB

	//block 1
	for(i1=CHAM_a_wt_min[K0_bloc[1]][carry_block[0]];
					i1<=CHAM_a_wt_min[K0_bloc[1]][carry_block[0]] +1; i1++) // CHAM_a_wt_max[X0_block[1]][carry_block[0]]
	{
		if(i0 + i1 + K0_bloc_wt_min[2] > w_cmp ){break;}
		for(j1=0; j1<CHAM_a_bg_numb[K0_bloc[1]][carry_block[0]][i1]; j1++)
		{
			alpha_bloc[1] = CHAM_a_beta[K0_bloc[1]][carry_block[0]][i1][j1];
			gama_block[1] = CHAM_a_gama[K0_bloc[1]][carry_block[0]][i1][j1];
			carry_block[1] = K0_bloc_msb[1] +
					((alpha_bloc[1] >> 7) << 1) + (gama_block[1] >> 7);

	//block 2
	for(i2=CHAM_a_wt_min[K0_bloc[2]][carry_block[1]];
							i2<=CHAM_a_wt_min[K0_bloc[2]][carry_block[1]] +1; i2++)  //CHAM_a_wt_max[X0_block[2]][carry_block[1]]
	{
		if(i0 + i1 + i2 +
				K0_bloc_wt_min[3] > w_cmp ){break;}
		for(j2=0; j2<CHAM_a_bg_numb[K0_bloc[2]][carry_block[1]][i2]; j2++)
		{
			alpha_bloc[2] = CHAM_a_beta[K0_bloc[2]][carry_block[1]][i2][j2];
			gama_block[2] = CHAM_a_gama[K0_bloc[2]][carry_block[1]][i2][j2];
			carry_block[2] = K0_bloc_msb[2] +
					((alpha_bloc[2] >> 7) << 1) + (gama_block[2] >> 7);

	//block 3
	for(i3=msb_CHAM_a_wt_min[K0_bloc[3]][carry_block[2]];
					i3<=msb_CHAM_a_wt_min[K0_bloc[3]][carry_block[2]] +1; i3++)  //<8 //msb_CHAM_a_wt_max[X0_block[3]][carry_block[2]]
	{
		P_w[cur_round] = i0 + i1 + i2 + i3;  //// 当前轮的概率重量
		if(P_w[cur_round] > w_cmp ){break;}
		//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
			speck_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
			speck_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
			RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
			RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

			for(j3=0; j3< msb_CHAM_a_bg_numb[K0_bloc[3]][carry_block[2]][i3]; j3++)
			{
				alpha_bloc[3] = msb_CHAM_a_beta[K0_bloc[3]][carry_block[2]][i3][j3];
				gama_block[3] = msb_CHAM_a_gama[K0_bloc[3]][carry_block[2]][i3][j3];

				alpha_tmp = ((alpha_bloc[3] << 24) | (alpha_bloc[2] << 16) |
						(alpha_bloc[1] << 8) | alpha_bloc[0]) & Bit_Align;
				gamma_tmp = ((gama_block[3] << 24) | (gama_block[2] << 16) |
						(gama_block[1] << 8) | gama_block[0]) & Bit_Align;

				//整理前两轮的临时输入输出RXD
				speck_key_RXD_tmp[2][0] = (ROTATE_LEFT(alpha_tmp,rol_a,blocksize_len) & Bit_Align);  //第2轮输出RXD

				speck_key_RXD_tmp[2][1] = speck_key_RXD_tmp[2][0]; //第0/1轮输入RXD,k2

				speck_key_RXD_tmp[1][2] = speck_key_RXD_tmp[2][0]; //第1轮输出RXD,k1

				//speck_key_RXD_tmp[1][3] = speck_key_RXD_tmp[2][2]; //第2轮输出RXD,k1
				speck_key_RXD_tmp[2][3] = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD,k1
				speck_key_RXD_tmp[2][3] ^= Speck_key_const_rxd[cur_round]; //轮常数; //第2轮输出RXD,k1
				speck_key_RXD_tmp[0][3] = speck_key_RXD_tmp[2][3] ^ next_Beta_tmp; //第2轮输出RXD,k1

				state = RX_sepck_2n3n_key_round_r(search_round, 3,
						speck_key_RXD_tmp[1][3],speck_key_RXD_tmp[0][3]);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(state == 1 )  //判断本轮是否要更新RX_Expected_w
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
			//breakfor = 1;
			break;
		}  //这里有问题，只应该退出i3的嵌套
	}}
#endif
			}
}}}}}}
////////////////////////////
#elif (speck_m_bits  == 144)  //SPECK-96/144





	////////////////////////////
#elif (speck_m_bits  == 192)  //SPECK-128/192


#endif
	return updated;
}

u16 RX_sepck_2n3n_key_round_r(u16 search_round, u16 cur_round, u64 x, u64 y)
{
	u16 best = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 xor_alpha_beta = 0;
//	u64 and_alpha_beta = 0;
	u16 w_xor = 0;
	float	w_cmp = 0;
	u16 i = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
//	u64 rol_left = 0;
//	u64 rol_right = 0;
	u64 alpha_bloc[8] = {0};
//	u64 gamma_block_tmp[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 j = 0;
	u64 beta_temp = 0;
	u64 gamma_temp = 0;
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};
//	u16 wt_purning = 63;
//	u16 wt_purning_chk = 63;


	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制
	updata_cur_num = Num_Bn_Update;


	if(search_round == cur_round)
	{
		best = RX_speck_2n3n_key_round_N(search_round, x, y );
		return best;
	}

	Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align;
	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;
	speck_key_RXD_tmp[1][cur_round+1] = speck_key_RXD_tmp[2][cur_round]; //输出RXD,k1=输入k2

//	rol_right = ((x) >> (rol_a)) | ((x) << (blocksize_len - rol_a)); //right
//	rol_left = (((y) << rol_b) | ((y) >> (blocksize_len - rol_b)));
//	Beta = (x ^ rol_left) & Bit_Align;
//	Alpha = rol_right & Bit_Align; // alpha = gama_uper <<< rol_a;
	//Beta = (x ^ ROTATE_LEFT(y,rol_b,blocksize_len) )& Bit_Align;
	//Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align; // alpha = gama_uper <<< rol_a;
	xor_alpha_beta = (Beta ^ Alpha) & ValueMax_Align ;
	w_xor = HM_weight(xor_alpha_beta );
	//w_xor = Speck_XDP_Max(Alpha, Beta);

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
				- RX_P_bestofR_w[search_round - cur_round];
	w_cmp = w_cmp_rx - Ek_min;
	if ( w_xor >= w_cmp)
	{
		return 0;
	}


////////////////////////////选择SPECK的对应分组版本////////////////////////////
#if (speck_m_bits  == 72)  //SPECK-48/72
		for(i=0;i<speck_m_keyword;i++)  //speck_m_keyword
		{
			alpha_bloc[i] = (Alpha >> (8*i)) & 0xFF; //8 bit
			beta_bloc[i]  = (Beta >> (8*i))  & 0xFF; //8 bit
			AB_block[i] = ((alpha_bloc[i] << 8) | beta_bloc[i]);
		}
		for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
		{
			w_xor_block[j] = (u8)_mm_popcnt_u64((xor_alpha_beta >> (j*8)));
		}
		carry[0] = 0;
		for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
		{
			carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}

////////Recombining////////////
		for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 + w_xor_block[1] > w_cmp ){break;}
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

		for(i1 = 0; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

		for(i2 = 0; i2 <= MSB_cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			P_w[cur_round] = w1 + i2;
			if(P_w[cur_round] > w_cmp ){break;}
			//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
			for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			{
				if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
				speck_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
				speck_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
				RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
				RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

				for(j2=0; j2 < MSB_cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
				{
					gamma_bloc[2] = MSB_cDDT_v[carry[2]][AB_block[2]][i2][j2];
					gamma_temp = (gamma_bloc[2] << 16) | (gamma_bloc[1] <<8)
							     | gamma_bloc[0];

					//整理临时输出RXD
					//speck_key_RXD_tmp[1][cur_round+1] = speck_key_RXD_tmp[2][cur_round]; //输出RXD,k1=输入k2
					speck_key_RXD_tmp[2][cur_round+1] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD,k1
					speck_key_RXD_tmp[2][cur_round+1] ^= Speck_key_const_rxd[cur_round]; //轮常数; //第2轮输出RXD,k1
					speck_key_RXD_tmp[0][cur_round+1] = speck_key_RXD_tmp[2][3] ^ next_Beta_tmp; //第2轮输出RXD,k1

					best = RX_sepck_2n3n_key_round_r(search_round, cur_round+1,
							speck_key_RXD_tmp[1][cur_round+1],
							speck_key_RXD_tmp[0][cur_round+1]);
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
			//breakfor = 1;
			break;
		}  //这里有问题，只应该退出i3的嵌套
	}}
#endif
		}
	}
}}}}}
////////////////////////////
#elif (speck_m_bits  == 96)  //SPECK-64/96
		for(i=0;i<speck_m_keyword;i++)  //speck_m_keyword
		{
			alpha_bloc[i] = (Alpha >> (8*i)) & 0xFF; //8 bit
			beta_bloc[i]  = (Beta >> (8*i))  & 0xFF; //8 bit
			AB_block[i] = ((alpha_bloc[i] << 8) | beta_bloc[i]);
		}
		for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
		{
			w_xor_block[j] = (u8)_mm_popcnt_u64((xor_alpha_beta >> (j*8)));
		}
		carry[0] = 0;
		for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
		{
			carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
				+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
		}


		////////Recombining////////////
				for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
				{
					if(i0 + w_xor_block[1] > w_cmp ){break;}
				for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
				{
					gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
					//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
					//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
					carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

				for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
				{
					w1 = i0 + i1;
					if(w1 + w_xor_block[2] > w_cmp ){break;}
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
					P_w[cur_round]= w2 + i3;
					if(P_w[cur_round] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
					for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
					{
						if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
						speck_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
						speck_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
						RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
						RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

						for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
						{
							gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
							gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
									   | (gamma_bloc[1] <<8) | gamma_bloc[0];

							//整理临时输出RXD
				//speck_key_RXD_tmp[1][cur_round+1] = speck_key_RXD_tmp[2][cur_round]; //输出RXD,k1=输入k2
				speck_key_RXD_tmp[2][cur_round+1] = gamma_temp ^ RX_zeta[RX_zeta_i[zeta_num]]; //第2轮输出RXD,k1
				speck_key_RXD_tmp[2][cur_round+1] ^= Speck_key_const_rxd[cur_round]; //轮常数; //第2轮输出RXD,k1
				speck_key_RXD_tmp[0][cur_round+1] = speck_key_RXD_tmp[2][3] ^ next_Beta_tmp; //第2轮输出RXD,k1

				best = RX_sepck_2n3n_key_round_r(search_round, cur_round+1,
									speck_key_RXD_tmp[1][cur_round+1],
									speck_key_RXD_tmp[0][cur_round+1]);
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
			//breakfor = 1;
			break;
		}  //这里有问题，只应该退出i3的嵌套
	}}
#endif
	}
	}
}}}}}}}
////////////////////////////
#elif (speck_m_bits  == 144)  //SPECK-96/144
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

		for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 + w_xor_block[1] > w_cmp ){break;}
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

		for(i1 = 0; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

		for(i2 = 0; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB

		for(i3 = 0; i3 <= MSB_cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
			gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					   | (gamma_bloc[1] <<8) | gamma_bloc[0];

			beta_temp = gamma_temp ^ next_Beta_tmp;

			P_w[cur_round] = w3; //i0 + i1 +i2 +i3;
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

			best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}}}}}}}}
////////////////////////////
#elif (speck_m_bits  == 192)  //SPECK-128/192
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
				i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]]
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
				i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
				carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
				i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  // cDDT_wt_max[carry[2]][AB_block[2]]
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
				i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
				gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
				carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
				i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++)  //cDDT_wt_max[carry[4]][AB_block[4]]
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5]> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
			{
				gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
				carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]];
				i5 <= MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //MSB_cDDT_wt_max[carry[5]][AB_block[5]]
		{
			w5 = w4 + i5;
			if(w5 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			P_w[cur_round] = w5;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

			for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
			{
			gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];
			gamma_temp = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];

			beta_temp = gamma_temp ^ next_Beta_tmp;

/*
 	 	 	 //限定每一个生产的beta的汉明重量
			wt_purning_chk = HM_weight(beta_temp);
			if(wt_purning_chk <= wt_purning )
			{
				wt_purning = wt_purning_chk;
				best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
			}
*/
			//if(w_xor + 9 > w5 ) //对应的DDT的概率在一个范围内
			{
				best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
			}

			//best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}
		}}}}}}}}}}}

/*
//////////////贪婪算法////////////每一轮的都只是概率最大的分支,该相同概率下的分支会有很多,但是这里只搜了其中一条,因此,遍历
		///的总的分支数就是第1轮的输入组合分支数

		//xor_alpha_beta = Beta ^ Alpha ;
		//w_xor = HM_weight(xor_alpha_beta & ValueMax_Align);
		and_alpha_beta = Beta & Alpha;
		w_and = HM_weight(and_alpha_beta & ValueMax_Align);

		P_w[cur_round] = w_xor + w_and;  // 00-->0;01,10-->1; 11-->0,概率的和
		if ( P_w[cur_round] > w_cmp)
		{
			return 0;
		}

//		gamma_temp |= V_MSB;
		gamma_temp = xor_alpha_beta; // 00-->1; 01,10--->0;11-->1;的情况先不考虑,因此只贪婪的考虑SetA的情况.

		best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
#if 1   //是否找到第一条最优路径就返回？
	if( best == 1 )
	{
		n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
		n_y_in[cur_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
		return 1;
	}
#endif
//////////////贪婪算法////////////
	*/
		////////////////////////////
#elif (SPECK_BYTE  == 128)  //SPECK-128
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
			//w_xor_block[j] = cDDT_AB_wt_min[AB_block[j]];
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

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
				i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]];
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//gamma_block_tmp[0] = gamma_bloc[0];
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
				i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
				i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  //cDDT_wt_max[carry[2]][AB_block[2]]
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
				i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break
		for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
			carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
				i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
		{
			gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
			carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = cDDT_wt_min[carry[5]][AB_block[5]];
				i5 <= cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
		{
			w5 = w4 + i5;
			if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
		{
			gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
			carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7); // gamma MSB
		for(i6 = cDDT_wt_min[carry[6]][AB_block[6]];
				i6 <= cDDT_wt_min[carry[6]][AB_block[6]]; i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
		{
			w6 = w5 + i6;
			if(w6 + w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
		{
			gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];

/*
///////////
			gamma_block_msb[6] = gamma_bloc[6] >> 7;
			if(gamma_block_msb_flag[7] == 0)  //本轮此处第一次进入下一个循环嵌套
			{
				carry[7] = carry_tmp[7] + gamma_block_msb[6]; // gamma MSB

			}
			else //本轮此处不是第一次进入下一个循环嵌套，
			{   //则需要判断gamma_bloc的最高位的值gamma_block_msb[6]以前是否出现过并break了

			}
*/

////////////

			carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7); // gamma MSB
		for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]];
				i7 <= MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
		{
			w7 = w6 + i7;
			if(w7> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			P_w[cur_round] = w7;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

		for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
		{
			gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];
//			gamma_temp = (gamma_bloc[7] << 56) | gamma_block_tmp[6]
//					| gamma_block_tmp[5] | gamma_block_tmp[4]
//					| gamma_block_tmp[3] | gamma_block_tmp[2]
//					| gamma_block_tmp[1] | gamma_bloc[0];
//
			//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
			gamma_temp = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
					| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];
			beta_temp = gamma_temp ^ next_Beta_tmp;

			best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}
		}}}}}}}}}}}}}}}


/*
//////////////贪婪算法////////////每一轮的都只是概率最大的分支,该相同概率下的分支会有很多,但是这里只搜了其中一条,因此,遍历
				///的总的分支数就是第1轮的输入组合分支数

				//xor_alpha_beta = Beta ^ Alpha ;
				//w_xor = HM_weight(xor_alpha_beta & ValueMax_Align);
				and_alpha_beta = Beta & Alpha;
				w_and = HM_weight(and_alpha_beta & ValueMax_Align);

				P_w[cur_round] = w_xor + w_and;  // 00-->0;01,10-->1; 11-->0,概率的和
				if ( P_w[cur_round] > w_cmp)
				{
					return 0;
				}

		//		gamma_temp |= V_MSB;
				gamma_temp = xor_alpha_beta; //对应最小概率了，00-->1; 01,10--->0;11-->1;的情况先不考虑,因此只贪婪的考虑SetA的情况.

				best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
		#if 1   //是否找到第一条最优路径就返回？
			if( best == 1 )
			{
				n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
				n_y_in[cur_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
				return 1;
			}
		#endif
//////////////贪婪算法////////////
*/
#endif
	return updated;
}



u16 RX_speck_2n3n_key_round_N(u16 search_round, u64 x, u64 y)
{
	u16 best = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 xor_alpha_beta = 0; u64 and_alpha_beta = 0;
	u16 w_xor = 0;   u16 w_and = 0;
	float w_cmp = 0;
	u16 i = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u64 rol_left = 0;
	u64 rol_right = 0;
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 j = 0;
	u64 gamma_temp = 0;
//	u64 gamma_block_tmp[8] = {0};
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};

	u32 zeta_num = 0;
	float w_cmp_rx = 0;



	Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align;
	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;

//	rol_right = ((x) >> (rol_a)) | ((x) << (blocksize_len - rol_a)); //right
//	rol_left = (((y) << rol_b) | ((y) >> (blocksize_len-rol_b)));
//	Beta = (x ^ rol_left) & Bit_Align;
//	Alpha = rol_right & Bit_Align; // alpha = gama_uper <<< rol_a;
	//Beta = x ^ (ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align);
	//Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align; // alpha = gama_uper <<< rol_a;
	xor_alpha_beta = (Beta ^ Alpha)  & ValueMax_Align;
	w_xor = HM_weight(xor_alpha_beta);

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[search_round -1];
	w_cmp = w_cmp_rx - Ek_min;
	if (w_xor >= w_cmp)
	{
		return 0;
	}

////////////////////////////选择SPECK的对应分组版本////////////////////////////
#if (speck_m_bits  == 72)  //SPECK-48/72
	for(i=0;i<speck_m_keyword;i++)  //speck_m_keyword
	{
		alpha_bloc[i] = (Alpha >> (8*i)) & 0xFF; //8 bit
		beta_bloc[i]  = (Beta >> (8*i))  & 0xFF; //8 bit
		AB_block[i] = ((alpha_bloc[i] << 8) | beta_bloc[i]);
	}
	for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
	{
		w_xor_block[j] = (u8)_mm_popcnt_u64((xor_alpha_beta >> (j*8)));
	}
	carry[0] = 0;
	for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
	{
		carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
			+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
	}

////////Recombining////////////
	for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
	{
		if(i0 + w_xor_block[1] > w_cmp ){break;}
	for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
	{
		gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
		//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
		//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
		carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

	for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
	{
		w1 = i0 + i1;
		if(w1 + w_xor_block[2] > w_cmp ){break;}
	for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
	{
		gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
		//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
		//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
		carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

//	for(i2 = MSB_cDDT_wt_min[carry[2]][AB_block[2]]; i2 <= MSB_cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
	{
		i2 = MSB_cDDT_wt_min[carry[2]][AB_block[2]];
		P_w[search_round] = w1 + i2;
		if(P_w[search_round] > w_cmp ){break;}
		//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
//		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			zeta_num=0;
			if(P_w[search_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
			speck_Pr_zeta_tmp[search_round] = RX_zeta_w[zeta_num];
			speck_zeta_tmp[search_round] = RX_zeta[RX_zeta_i[zeta_num]];
			RX_P_w[search_round] = P_w[search_round] + RX_zeta_w[zeta_num];
			RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1] + RX_P_w[search_round];

//			for(j2=0; j2 < MSB_cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
				j2=0;
				gamma_bloc[2] = MSB_cDDT_v[carry[2]][AB_block[2]][i2][j2];
				gamma_temp = (gamma_bloc[2] << 16) | (gamma_bloc[1] <<8)
						     | gamma_bloc[0];

				//整理输出RXD
				speck_key_RXD[1][search_round+1] = speck_key_RXD_tmp[2][search_round]; //输出RXD,k1=输入k2
				speck_key_RXD[2][search_round+1] = gamma_temp
						^ RX_zeta[RX_zeta_i[zeta_num]]; //下轮输出RXD
				speck_key_RXD[2][search_round+1] ^= Speck_key_const_rxd[search_round]; //轮常数
				speck_key_RXD[0][search_round+1] = speck_key_RXD[2][search_round+1]
																	^next_Beta_tmp;  //下一轮输出RXD

				//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
						Num_Bn_Update++;  //更新次数记录
						RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
						w_cmp = P_w[search_round];
						best = 1;
				//		breakfor =1; //快速退出机制
						updated = 1;
				//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
					for(i=0; i <= search_round; i++ )
					{
						speck_key_RXD[0][i] = speck_key_RXD_tmp[0][i];
						speck_key_RXD[1][i] = speck_key_RXD_tmp[1][i];
						speck_key_RXD[2][i] = speck_key_RXD_tmp[2][i];

						speck_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
						speck_Pr_XOR[i] = P_w[i];   //每轮的模加的XOR差分概率重量
						speck_Pr_zeta[i] =  speck_Pr_zeta_tmp[i];
						speck_Zeta[i] =  speck_zeta_tmp[i];
					}

					time_Round = clock();
					run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
					printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
					printf("Time: %.2f seconds. ==>>\n", run_time);

					print_RX_speck_key_resoult(search_round);  //打印
	}
}}}}}}
////////////////////////////
#elif (speck_m_bits  == 96)  //SPECK-64/96
	for(i=0;i<speck_m_keyword;i++)  //speck_m_keyword
	{
		alpha_bloc[i] = (Alpha >> (8*i)) & 0xFF; //8 bit
		beta_bloc[i]  = (Beta >> (8*i))  & 0xFF; //8 bit
		AB_block[i] = ((alpha_bloc[i] << 8) | beta_bloc[i]);
	}
	for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
	{
		w_xor_block[j] = (u8)_mm_popcnt_u64((xor_alpha_beta >> (j*8)));
	}
	carry[0] = 0;
	for(j=1;j<speck_m_keyword;j++)  //speck_m_keyword
	{
		carry_tmp[j] = ((alpha_bloc[j-1] >> 7) << 2) //alpha MSB
			+ ((beta_bloc[j-1] >> 7) << 1);   //beta MSB
	}


	////////Recombining////////////
			for(i0 = cDDT_wt_min[carry[0]][AB_block[0]]; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
			{
				if(i0 + w_xor_block[1] > w_cmp ){break;}
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

			for(i1 = cDDT_wt_min[carry[1]][AB_block[1]]; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
			{
				w1 = i0 + i1;
				if(w1 + w_xor_block[2] > w_cmp ){break;}
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

//			for(i3 = MSB_cDDT_wt_min[carry[3]][AB_block[3]]; i3 <= MSB_cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
			{
				i3 = MSB_cDDT_wt_min[carry[3]][AB_block[3]];
				P_w[search_round]= w2 + i3;
				if(P_w[search_round] >= w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
				for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					if(P_w[search_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
					speck_Pr_zeta_tmp[search_round] = RX_zeta_w[zeta_num];
					speck_zeta_tmp[search_round] = RX_zeta[RX_zeta_i[zeta_num]];
					RX_P_w[search_round] = P_w[search_round] + RX_zeta_w[zeta_num];
					RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1] + RX_P_w[search_round];

//					for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
					{
						j3=0;
						gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
						gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
								   | (gamma_bloc[1] <<8) | gamma_bloc[0];

						//整理输出RXD
						speck_key_RXD[1][search_round+1] = speck_key_RXD_tmp[2][search_round]; //输出RXD,k1=输入k2
						speck_key_RXD[2][search_round+1] = gamma_temp
								^ RX_zeta[RX_zeta_i[zeta_num]]; //下轮输出RXD
						speck_key_RXD[2][search_round+1] ^= Speck_key_const_rxd[search_round]; //轮常数
						speck_key_RXD[0][search_round+1] = speck_key_RXD[2][search_round+1]
																			^next_Beta_tmp;  //下一轮输出RXD

						//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
								Num_Bn_Update++;  //更新次数记录
								RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
								w_cmp = P_w[search_round];
								best = 1;
						//		breakfor =1; //快速退出机制
								updated = 1;
						//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
							for(i=0; i <= search_round; i++ )
							{
								speck_key_RXD[0][i] = speck_key_RXD_tmp[0][i];
								speck_key_RXD[1][i] = speck_key_RXD_tmp[1][i];
								speck_key_RXD[2][i] = speck_key_RXD_tmp[2][i];

								speck_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
								speck_Pr_XOR[i] = P_w[i];   //每轮的模加的XOR差分概率重量
								speck_Pr_zeta[i] =  speck_Pr_zeta_tmp[i];
								speck_Zeta[i] =  speck_zeta_tmp[i];
							}

							time_Round = clock();
							run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
							printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
							printf("Time: %.2f seconds. ==>>\n", run_time);

							print_RX_speck_key_resoult(search_round);  //打印
}
}
}}}}}}}
///////////////////////////
#elif (speck_m_bits  == 144)  //SPECK-96/144
////////////////////////////
#elif (speck_m_bits  == 192)  //SPECK-128/192

#endif
return  updated;
}
/////////
#elif (speck_m_keyword  == 4)  //speck2n/4n: 4
/////////
u16 RX_sepck_2n4n_key_round_1(u16 search_round)
{
	u16 best = 0;
	u16 thr_d = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;


	u64 M0 = 0;
	u64 N0 = blocksize_len-1;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;


///////////////////round 1 entry/////////////////////
	if( search_round == 1)
	{
		best = round_1_j(1);
		return best;
	}
	else
	{	// 可以选择是否限定第1轮的概率重量
		for(thr_d = 0; thr_d < blocksize_len-1 ;thr_d++)  //0::n-1 //blocksize_len-1
		{
			// thr_d : the bits number of alpha/beta/gamma that the position not equal at the same time.
			if ((thr_d + n_P_bestofR_w[search_round - 1]) > Bn_w) // Pmax_compute(xi) is the ri round min weight.
			{
				return 0;  //Return the upper procedure::sepck_differential_trail_search_entry
			}

			M0 = thr_d;  //三个alpha,beta和gamma同时考虑
			P_w[1] = M0;  //pro; //Speck_XDP_compute(Input_Alpha,Input_Beta,Input_Gamma );
			p_sumof_r_P_w[1] = P_w[1];
			if(M0 == 0)  /// abd [0 -> n-2] 没有完全不相同的比特,且MSB之后只能为0
			{
				for(j_abc = 1; j_abc < 4;j_abc++)
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

					best = sepck_round_r(search_round, 2,Input_gamma,Input_gamma ^ (ROTATE_LEFT(Input_beta, rol_b,blocksize_len) & Bit_Align));
#if 1   //是否找到第一条最优路径就返回？
		if(best == 1 )
		{
			n_x_in[1] = (ROTATE_LEFT(Input_alpha, rol_a,blocksize_len)) & Bit_Align; // left part of Input of the first round.
			n_y_in[1] = Input_beta;  // right part of Input of the first round.

			n_x_in[2 ] = Input_gamma;  // left part output of fisrt round.
			n_y_in[2 ] = Input_gamma ^ (ROTATE_LEFT(Input_beta, rol_b,blocksize_len) & Bit_Align);
			return 1;
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
					best = sepck_input_MSB(search_round,Input_alpha, Input_beta, Input_gamma,thr_d,C0);
#if 1   //是否找到第一条最优路径就返回？
	if(best == 1 )
	{
		return 1;
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
			best = sepck_input_MSB(search_round,Input_alpha, Input_beta, Input_gamma,thr_d,C0);

#if 1   //是否找到第一条最优路径就返回？
			if(best == 1 )
			{
				return 1;
			}
#endif
///////////////////////////////////////
					} while(1);
				}
		}
	}
	return best;
}



u16 RX_sepck_2n4n_key_input_MSB
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
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

			state = sepck_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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

		state = sepck_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
		}

		//state = sepck_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
#if 1   //是否找到第一条最优路径就返回？
	if(state == 1 )
	{
		return 1;
	}
#endif
}
	return state;
}

//由高位到低位逐步判断
u16 RX_sepck_2n4n_key_input_Middle
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi, u16 cur_posi)
{
	u16 state = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;
	u16 i = 0;
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


			state = sepck_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
#if 1   //是否找到第一条最优路径就返回？
	if(state == 1 )
	{
		return 1;
	}
#endif
		}
	}
	else 	//call last position.
	{
		state = sepck_input_Last(search_round,
				alpha,beta,	gamma,P_1,posi );
#if 1   //是否找到第一条最优路径就返回？
	if(state == 1 )
	{
		return 1;
	}
#endif
	}
	return state;
}

u16 RX_sepck_2n4n_key_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0;
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u16 j_abc = 0;
	u16 j_last = 0;
	u64 indx_tmp = 0;

	u64 bit_i = 1;


	indx_tmp = posi[1] - 1;

	// 最后的一个产生概率重量的比特位置只能为{3、5/6}，在该位置之后的更低的比特位置都全部清0；
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
		//P_w[1] = P_1;
		state = sepck_round_r(search_round, 2,(Input_gamma & Bit_Align),Input_gamma ^ (ROTATE_LEFT(Input_beta, rol_b,blocksize_len) & Bit_Align));
#if 1   //是否找到第一条最优路径就返回？
	if(state == 1 )
	{
		n_x_in[1] = (ROTATE_LEFT(Input_alpha, rol_a, blocksize_len)) & Bit_Align; // left part of Input of the first round.
		n_y_in[1] = Input_beta;  // right part of Input of the first round.

		n_x_in[2 ] = Input_gamma;  // left part output of fisrt round.
		n_y_in[2 ] = Input_gamma ^ (ROTATE_LEFT(Input_beta, rol_b,blocksize_len) & Bit_Align);
		return 1;
	}
#endif
	}
	return state;
}


u16 RX_sepck_2n4n_key_round_r(u16 search_round, u16 cur_round, u64 x, u64 y)
{
	u16 best = 0;
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 xor_alpha_beta = 0;
//	u64 and_alpha_beta = 0;
	u16 w_xor = 0;
	u16	w_cmp = 0;
	u16 i = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
//	u64 rol_left = 0;
//	u64 rol_right = 0;
	u64 alpha_bloc[8] = {0};
//	u64 gamma_block_tmp[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 j = 0;
	u64 beta_temp = 0;
	u64 gamma_temp = 0;
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};
//	u16 wt_purning = 63;
//	u16 wt_purning_chk = 63;

	if(search_round == cur_round)
	{
		best = speck_round_N(search_round, x, y );
		return best;
	}

	Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align;
	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;

//	rol_right = ((x) >> (rol_a)) | ((x) << (blocksize_len - rol_a)); //right
//	rol_left = (((y) << rol_b) | ((y) >> (blocksize_len - rol_b)));
//	Beta = (x ^ rol_left) & Bit_Align;
//	Alpha = rol_right & Bit_Align; // alpha = gama_uper <<< rol_a;
	//Beta = (x ^ ROTATE_LEFT(y,rol_b,blocksize_len) )& Bit_Align;
	//Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align; // alpha = gama_uper <<< rol_a;
	xor_alpha_beta = (Beta ^ Alpha) & ValueMax_Align ;
	w_xor = HM_weight(xor_alpha_beta );
	//w_xor = Speck_XDP_Max(Alpha, Beta);

/*
	for(i=1;i<cur_round;i++)
	{
		p_sumof_r += P_w[i]; // p_sumof_r is the sum of r-1 rounds weight.
	}
*/

	w_cmp = Bn_w - p_sumof_r_P_w[cur_round -1] - n_P_bestofR_w[search_round - cur_round];
	if ( w_xor > w_cmp)
	{
		return 0;
	}

//双目剪枝方法，另外一种判断变量逻辑关系，理论支撑？

////////////////////////////选择SPECK的对应分组版本////////////////////////////
#if (SPECK_BYTE  == 32)  //SPECK-32
		//m0
		alpha_bloc[0] = (Alpha) & 0xFF; //8 bit
		beta_bloc[0]  = (Beta) & 0xFF; //8 bit
		AB_block[0] = (alpha_bloc[0] << 8) | beta_bloc[0];
		carry[0] = 0;
		//gamma_num[0] = cDDT[carry[0]][AB_block[0]][256];
		//m1
		alpha_bloc[1]= (Alpha >> 8) & 0xFF; //8 bit
		beta_bloc[1] = ( Beta >> 8) & 0xFF; //8 bit
		AB_block[1] = (alpha_bloc[1] << 8) | beta_bloc[1];


		for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
					+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
					+ (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = 0; i1 <= MSB_cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			if(i0 + i1 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < MSB_cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = MSB_cDDT_v[carry[1]][AB_block[1]][i1][j1];
			gamma_temp = (gamma_bloc[1] <<8) | gamma_bloc[0];
			beta_temp = gamma_temp ^ next_Beta_tmp;
			P_w[cur_round] = i0 + i1;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];
			best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}}}}
////////////////////////////
#elif (SPECK_BYTE  == 48)  //SPECK-48
	case 48:
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

		for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 + w_xor_block[1] > w_cmp ){break;}
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

		for(i1 = 0; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

		for(i2 = 0; i2 <= MSB_cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 > w_cmp ){break;}
		for(j2=0; j2 < MSB_cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = MSB_cDDT_v[carry[2]][AB_block[2]][i2][j2];
			gamma_temp = (gamma_bloc[2] << 16) | (gamma_bloc[1] <<8) | gamma_bloc[0];
			beta_temp = gamma_temp ^ next_Beta_tmp;

			P_w[cur_round] = w2;  //i0 + i1 +i2;
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

			best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}}}}}}
////////////////////////////
#elif (SPECK_BYTE  == 64)  //SPECK-64
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

		for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
		{
			if(i0 + w_xor_block[1] > w_cmp ){break;}
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

		for(i1 = 0; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

		for(i2 = 0; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
			//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB

		for(i3 = 0; i3 <= MSB_cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
		{
			w3 = w2 + i3;
			if(w3 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
			gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					   | (gamma_bloc[1] <<8) | gamma_bloc[0];

			beta_temp = gamma_temp ^ next_Beta_tmp;

			P_w[cur_round] = w3; //i0 + i1 +i2 +i3;
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

			best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}}}}}}}}
////////////////////////////
#elif (SPECK_BYTE  == 96)  //SPECK-96
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
				i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]]
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
				i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
				carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
				i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  // cDDT_wt_max[carry[2]][AB_block[2]]
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
				i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
				gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
				carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
				i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++)  //cDDT_wt_max[carry[4]][AB_block[4]]
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5]> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
			{
				gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
				carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]];
				i5 <= MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //MSB_cDDT_wt_max[carry[5]][AB_block[5]]
		{
			w5 = w4 + i5;
			if(w5 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			P_w[cur_round] = w5;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

			for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
			{
			gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];
			gamma_temp = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];

			beta_temp = gamma_temp ^ next_Beta_tmp;

/*
 	 	 	 //限定每一个生产的beta的汉明重量
			wt_purning_chk = HM_weight(beta_temp);
			if(wt_purning_chk <= wt_purning )
			{
				wt_purning = wt_purning_chk;
				best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
			}
*/
			//if(w_xor + 9 > w5 ) //对应的DDT的概率在一个范围内
			{
				best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
			}

			//best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}
		}}}}}}}}}}}

/*
//////////////贪婪算法////////////每一轮的都只是概率最大的分支,该相同概率下的分支会有很多,但是这里只搜了其中一条,因此,遍历
		///的总的分支数就是第1轮的输入组合分支数

		//xor_alpha_beta = Beta ^ Alpha ;
		//w_xor = HM_weight(xor_alpha_beta & ValueMax_Align);
		and_alpha_beta = Beta & Alpha;
		w_and = HM_weight(and_alpha_beta & ValueMax_Align);

		P_w[cur_round] = w_xor + w_and;  // 00-->0;01,10-->1; 11-->0,概率的和
		if ( P_w[cur_round] > w_cmp)
		{
			return 0;
		}

//		gamma_temp |= V_MSB;
		gamma_temp = xor_alpha_beta; // 00-->1; 01,10--->0;11-->1;的情况先不考虑,因此只贪婪的考虑SetA的情况.

		best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
#if 1   //是否找到第一条最优路径就返回？
	if( best == 1 )
	{
		n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
		n_y_in[cur_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
		return 1;
	}
#endif
//////////////贪婪算法////////////
	*/
		////////////////////////////
#elif (SPECK_BYTE  == 128)  //SPECK-128
		for(j=0;j<nBytes;j++)  //nBytes
		{
			alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
			beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
			AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
			//w_xor_block[j] = cDDT_AB_wt_min[AB_block[j]];
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

		for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
				i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]];
		{
			//w0 = i0;
			if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
		{
			gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
			//gamma_block_tmp[0] = gamma_bloc[0];
			carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
		for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
				i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
		{
			w1 = i0 + i1;
			if(w1 + w_xor_block[2] > w_cmp ){break;}  //break
		for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
		{
			gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
			carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
		for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
				i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  //cDDT_wt_max[carry[2]][AB_block[2]]
		{
			w2 = w1 + i2;
			if(w2 + w_xor_block[3] > w_cmp ){break;}  //break
		for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
		{
			gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
			carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
		for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
				i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
		{
			w3 = w2 + i3;
			if(w3 + w_xor_block[4] > w_cmp ){break;}  //break
		for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
		{
			gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
			carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
		for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
				i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
		{
			w4 = w3 + i4;
			if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
		{
			gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
			carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
		for(i5 = cDDT_wt_min[carry[5]][AB_block[5]];
				i5 <= cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
		{
			w5 = w4 + i5;
			if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
		for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
		{
			gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
			carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7); // gamma MSB
		for(i6 = cDDT_wt_min[carry[6]][AB_block[6]];
				i6 <= cDDT_wt_min[carry[6]][AB_block[6]]; i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
		{
			w6 = w5 + i6;
			if(w6 + w_xor_block[7] > w_cmp ){break;}  //break //break 要比 continue 高效,直接结束.
		for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
		{
			gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];

/*
///////////
			gamma_block_msb[6] = gamma_bloc[6] >> 7;
			if(gamma_block_msb_flag[7] == 0)  //本轮此处第一次进入下一个循环嵌套
			{
				carry[7] = carry_tmp[7] + gamma_block_msb[6]; // gamma MSB

			}
			else //本轮此处不是第一次进入下一个循环嵌套，
			{   //则需要判断gamma_bloc的最高位的值gamma_block_msb[6]以前是否出现过并break了

			}
*/

////////////

			carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7); // gamma MSB
		for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]];
				i7 <= MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
		{
			w7 = w6 + i7;
			if(w7> w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			P_w[cur_round] = w7;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);
			p_sumof_r_P_w[cur_round] = p_sumof_r_P_w[cur_round -1] + P_w[cur_round];

		for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
		{
			gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];
//			gamma_temp = (gamma_bloc[7] << 56) | gamma_block_tmp[6]
//					| gamma_block_tmp[5] | gamma_block_tmp[4]
//					| gamma_block_tmp[3] | gamma_block_tmp[2]
//					| gamma_block_tmp[1] | gamma_bloc[0];
//
			//这个gamma太特么多了，导致分支数太多，搜索复杂度太高
			gamma_temp = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
					| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
					| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
					| (gamma_bloc[1] <<8) | gamma_bloc[0];
			beta_temp = gamma_temp ^ next_Beta_tmp;

			best = sepck_round_r(search_round,cur_round +1,gamma_temp, beta_temp);
#if 1   //是否找到第一条最优路径就返回？
		if( best == 1 )
		{
			n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[cur_round +1] = beta_temp; // Y
			return 1;
		}
#endif
		}
		}}}}}}}}}}}}}}}


/*
//////////////贪婪算法////////////每一轮的都只是概率最大的分支,该相同概率下的分支会有很多,但是这里只搜了其中一条,因此,遍历
				///的总的分支数就是第1轮的输入组合分支数

				//xor_alpha_beta = Beta ^ Alpha ;
				//w_xor = HM_weight(xor_alpha_beta & ValueMax_Align);
				and_alpha_beta = Beta & Alpha;
				w_and = HM_weight(and_alpha_beta & ValueMax_Align);

				P_w[cur_round] = w_xor + w_and;  // 00-->0;01,10-->1; 11-->0,概率的和
				if ( P_w[cur_round] > w_cmp)
				{
					return 0;
				}

		//		gamma_temp |= V_MSB;
				gamma_temp = xor_alpha_beta; //对应最小概率了，00-->1; 01,10--->0;11-->1;的情况先不考虑,因此只贪婪的考虑SetA的情况.

				best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
		#if 1   //是否找到第一条最优路径就返回？
			if( best == 1 )
			{
				n_x_in[cur_round +1] = gamma_temp;  // output part of the r-th round..
				n_y_in[cur_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
				return 1;
			}
		#endif
//////////////贪婪算法////////////
*/
#endif
	return best;
}



u16 RX_speck_2n4n_key_round_N(u16 search_round, u64 x, u64 y)
{
	u16 best = 0;
	u64 Alpha = 0;
	u64 Beta = 0;
	u64 next_Beta_tmp = 0;
	u64 xor_alpha_beta = 0; u64 and_alpha_beta = 0;
	u16 w_xor = 0;   u16 w_and = 0;
	u16 w_cmp = 0;
	u16 i = 0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	u64 rol_left = 0;
	u64 rol_right = 0;
	u64 alpha_bloc[8] = {0};
	u64 beta_bloc[8] = {0};
	u64 gamma_bloc[8] = {0};
	u64 AB_block[8] = {0};
	u64 carry[8] ={0};
	u64 carry_tmp[8] ={0};
	u64 i0=0,i1=0,i2=0,i3=0,i4=0,i5=0,i6=0,i7=0;
	u64 j0=0,j1=0,j2=0,j3=0,j4=0,j5=0,j6=0,j7=0;
	u64 j = 0;
	u64 gamma_temp = 0;
//	u64 gamma_block_tmp[8] = {0};
	u64 w1=0,w2=0,w3=0,w4=0,w5=0,w6=0,w7=0;
	u8 w_xor_block[8] = {0};


	Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align;
	Beta = y & Bit_Align;
	next_Beta_tmp = ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align;


//	rol_right = ((x) >> (rol_a)) | ((x) << (blocksize_len - rol_a)); //right
//	rol_left = (((y) << rol_b) | ((y) >> (blocksize_len-rol_b)));
//	Beta = (x ^ rol_left) & Bit_Align;
//	Alpha = rol_right & Bit_Align; // alpha = gama_uper <<< rol_a;
	//Beta = x ^ (ROTATE_LEFT(y,rol_b,blocksize_len) & Bit_Align);
	//Alpha = (ROTATE_RIGHT(x,rol_a, blocksize_len)) & Bit_Align; // alpha = gama_uper <<< rol_a;
	xor_alpha_beta = (Beta ^ Alpha)  & ValueMax_Align;
	w_xor = HM_weight(xor_alpha_beta);

	//w_xor = Speck_XDP_Max(Alpha, Beta, &gamma_temp);
/*
	for(i=1; i < search_round;i++)
	{
		p_sumof_r += P_w[i]; // p_sumof_r is the sum of r-1 rounds weight.
	}
*/
	w_cmp = Bn_w - p_sumof_r_P_w[search_round -1];
	if (w_xor > w_cmp)
	{
		return 0;
	}

	/*
	if (w_xor == w_cmp)
	{
		P_w[search_round] = w_xor;
		//gamma_temp = Speck_XDP_Max_Gamma(Alpha, Beta);
		best  = 1;
		n_P_bestofR_w[search_round] = Bn_w;
		n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
		n_y_in[search_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
		return 1;
	}
	else
	{
		return 0;
	}
	*/

	////////////////////////////选择SPECK的对应分组版本////////////////////////////
	#if (speck_m_bits  == 72)  //SPECK-48/72




	////////////////////////////
	#elif (speck_m_bits  == 96)  //SPECK-64/96

	///////////////////////////
	#elif (speck_m_bits  == 144)  //SPECK-96/144
	////////////////////////////
	#elif (speck_m_bits  == 192)  //SPECK-128/192

	#endif


////////////////////////////选择SPECK的对应分组版本////////////////////////////
#if (SPECK_BYTE  == 32)  //SPECK-32
			//m0
			alpha_bloc[0] = (Alpha) & 0xFF; //8 bit
			beta_bloc[0]  = (Beta) & 0xFF; //8 bit
			AB_block[0] = (alpha_bloc[0] << 8) | beta_bloc[0];
			carry[0] = 0;
			//m1
			alpha_bloc[1]= (Alpha >> 8) & 0xFF; //8 bit
			beta_bloc[1] = ( Beta >> 8) & 0xFF; //8 bit
			AB_block[1] = (alpha_bloc[1] << 8) | beta_bloc[1];


			for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
			{
				if(i0 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
						+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
						+ (gamma_bloc[0] >> 7); // gamma MSB
			for(i1 = 0; i1 <= MSB_cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
			{
				if(i0 + i1 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j1=0; j1 < MSB_cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				P_w[search_round] = i0 + i1;  //Speck_XDP_compute(Alpha, Beta, gamma_temp);

				if (P_w[search_round] == w_cmp) //概率满足Matsui条件
				{
					gamma_bloc[1] = MSB_cDDT_v[carry[1]][AB_block[1]][i1][j1];
					gamma_temp = (gamma_bloc[1] <<8) | gamma_bloc[0];

					//printf("gamma[%d]: %x   \n", j2, gamma_temp);

					//best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
#if 1   //是否找到第一条最优路径就返回？
			//if( best == 1 )
			{
				best  = 1;
				n_P_bestofR_w[search_round] = Bn_w;
				n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
				n_y_in[search_round +1] = gamma_temp ^ next_Beta_tmp;
				return 1;
			}
#endif
		}}}}}
////////////////////////////
#elif (SPECK_BYTE  == 48)  //SPECK-48
			for(j=0;j<nBytes;j++)  //nBytes
				{
					alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
					beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
					AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

					for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
					{
						if(i0 + w_xor_block[1] > w_cmp ){break;}
					for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
					{
						gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
						//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
						carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
					for(i1 = 0; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
					{
						w1 = i0 + i1;
						if(w1 + w_xor_block[2] > w_cmp ){break;}
					for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
					{
						gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
						//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
						carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

					for(i2 = 0; i2 <= MSB_cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
					{
						w2 = w1 + i2;
						if(w2 > w_cmp ){break;}
					for(j2=0; j2 < MSB_cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
					{
						gamma_bloc[2] = MSB_cDDT_v[carry[2]][AB_block[2]][i2][j2];
						gamma_temp = (gamma_bloc[2] << 16) | (gamma_bloc[1] <<8) | gamma_bloc[0];

						P_w[search_round] = w2; //i0 + i1 +i2;
						if (P_w[search_round] == w_cmp) //概率满足Matsui条件
						{
				#if 1   //是否找到第一条最优路径就返回？
						//if( best == 1 )
						{
							best  = 1;
							n_P_bestofR_w[search_round] = Bn_w;
							n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
							n_y_in[search_round +1] = gamma_temp ^ next_Beta_tmp;
							return 1;
						}
				#endif
						}
					}
					}
					}
					}
					}
					}
////////////////////////////
#elif (SPECK_BYTE  == 64)  //SPECK-64
			for(j=0;j<nBytes;j++)  //nBytes
			{
				alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
				beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
				AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

			for(i0 = 0; i0 <= cDDT_wt_max[carry[0]][AB_block[0]]; i0++)
			{
				if(i0 + w_xor_block[1] > w_cmp ){break;}
					for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
					{
						gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
						//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
						carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

					for(i1 = 0; i1 <= cDDT_wt_max[carry[1]][AB_block[1]]; i1++)
					{
						w1 = i0 + i1;
						if(w1 + w_xor_block[2] > w_cmp ){break;}
					for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
					{
						gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
						//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
						carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

					for(i2 = 0; i2 <= cDDT_wt_max[carry[2]][AB_block[2]]; i2++)
					{
						w2 = w1 + i2;
						if(w2 + w_xor_block[3] > w_cmp ){break;}
					for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
					{
						gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
						//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
						carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB

					for(i3 = 0; i3 <= MSB_cDDT_wt_max[carry[3]][AB_block[3]]; i3++)
					{
						w3 = w2 + i3;
						if(w3 > w_cmp ){break;}
					for(j3=0; j3 < MSB_cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
					{
						gamma_bloc[3] = MSB_cDDT_v[carry[3]][AB_block[3]][i3][j3];
						gamma_temp = (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
								   | (gamma_bloc[1] <<8) | gamma_bloc[0];
						P_w[search_round] = w3; //i0 + i1 +i2 +i3;
						if (P_w[search_round] == w_cmp) //概率满足Matsui条件
						{
			#if 1   //是否找到第一条最优路径就返回？
					//if( best == 1 )
					{
						best  = 1;
						n_P_bestofR_w[search_round] = Bn_w;
						n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
						n_y_in[search_round +1] = gamma_temp ^ next_Beta_tmp;
						return 1;
/*
						if(flag == 0)
						{
							flag =1;
							return 0;
						}
						else
						{
							return 1;
						}
*/
					}
			#endif
					}}}}}}}}}
////////////////////////////
#elif (SPECK_BYTE  == 96)  //SPECK-96
				for(j=0;j<nBytes;j++)  //nBytes
				{
					alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
					beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
					AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

					for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
							i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++) //cDDT_wt_max[carry[0]][AB_block[0]]
					{
						if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
					{
						gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
						//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
						carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB
					for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
							i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++) //cDDT_wt_max[carry[1]][AB_block[1]]
					{
						w1 = i0 + i1;
						if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
					{
						gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
						//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
						carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB
					for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
							i2 <= cDDT_wt_min[carry[2]][AB_block[2]];i2++) // cDDT_wt_max[carry[2]][AB_block[2]];
					{
						w2 = w1 + i2;
						if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
					{
						gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
						//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
						carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
					for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
							i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++) //cDDT_wt_max[carry[3]][AB_block[3]]
					{
						w3 = w2 + i3;
						if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.

					for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
					{
						gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
						//carry[4] = ((alpha_bloc[3] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[3] >> 7) << 1)   //beta MSB
						carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB
					for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
							i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++) //cDDT_wt_max[carry[4]][AB_block[4]]
					{
						w4 = w3 + i4;
						if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
					for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
					{
						gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
						//carry[5] = ((alpha_bloc[4] >> 7) << 2) //alpha MSB
						//		+ ((beta_bloc[4] >> 7) << 1)   //beta MSB
						carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB
					for(i5 = MSB_cDDT_wt_min[carry[5]][AB_block[5]];
							i5 <= MSB_cDDT_wt_min[carry[5]][AB_block[5]]; i5++) //MSB_cDDT_wt_max[carry[5]][AB_block[5]];
					{
						w5 = w4 + i5;
						if(w5 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
						P_w[search_round] = w5; //Speck_XDP_compute(Alpha, Beta, gamma_temp);
					for(j5=0; j5 < MSB_cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
					{
						gamma_bloc[5] = MSB_cDDT_v[carry[5]][AB_block[5]][i5][j5];
						gamma_temp = (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
								| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
								| (gamma_bloc[1] <<8) | gamma_bloc[0];

						if (P_w[search_round] == w_cmp) //概率满足Matsui条件
						//if ((P_w[search_round] <= w_cmp) && (P_w[search_round] != 0))
						{
#if 1   //是否找到第一条最优路径就返回？
		//if( best == 1 )
		{
			best  = 1;
			n_P_bestofR_w[search_round] = Bn_w;
			n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
			n_y_in[search_round +1] = gamma_temp ^ next_Beta_tmp;
			return 1;
		}
#endif
	    }}}}}}}}}}}}}
/*
//////////////贪婪算法////////////每一轮的都只是概率最大的分支,该相同概率下的分支会有很多,但是这里只搜了其中一条,因此,遍历
							///的总的分支数就是第1轮的输入组合分支数

							//xor_alpha_beta = Beta ^ Alpha ;
							//w_xor = HM_weight(xor_alpha_beta & ValueMax_Align);
							and_alpha_beta = Beta & Alpha;
							w_and = HM_weight(and_alpha_beta & ValueMax_Align);

							P_w[search_round] = w_xor + w_and;  // 00-->0;01,10-->1; 11-->0,概率的和
							if ( P_w[search_round] == w_cmp)
							{
							gamma_temp = xor_alpha_beta; // 00-->1; 01,10--->0;11-->1;的情况先不考虑,因此只贪婪的考虑SetA的情况.

							//best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
					#if 1   //是否找到第一条最优路径就返回？
						//if( best == 1 )
						{
							best  = 1;
							n_P_bestofR_w[search_round] = Bn_w;
							n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
							n_y_in[search_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
							return 1;
						}
					#endif
							}
//////////////贪婪算法////////////
*/
////////////////////////////
#elif (SPECK_BYTE  == 128)  //SPECK-128
			for(j=0;j<nBytes;j++)  //nBytes
			{
				alpha_bloc[j] = (Alpha >> (8*j)) & 0xFF; //8 bit
				beta_bloc[j]  = (Beta >> (8*j))  & 0xFF; //8 bit
				AB_block[j] = ((alpha_bloc[j] << 8) | beta_bloc[j]);
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

			for(i0 = cDDT_wt_min[carry[0]][AB_block[0]];
					i0 <= cDDT_wt_min[carry[0]][AB_block[0]]; i0++)  //cDDT_wt_max[carry[0]][AB_block[0]]
			{
				if(i0 + w_xor_block[1] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j0=0; j0 < cDDT_n[carry[0]][AB_block[0]][i0]; j0++)
			{
				gamma_bloc[0] = cDDT_v[carry[0]][AB_block[0]][i0][j0];
				//gamma_block_tmp[0] = gamma_bloc[0];
				//carry[1] = ((alpha_bloc[0] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[0] >> 7) << 1)   //beta MSB
				carry[1] = carry_tmp[1] + (gamma_bloc[0] >> 7); // gamma MSB

			for(i1 = cDDT_wt_min[carry[1]][AB_block[1]];
					i1 <= cDDT_wt_min[carry[1]][AB_block[1]]; i1++)  //cDDT_wt_max[carry[1]][AB_block[1]]
			{
				w1 = i0 + i1;
				if(w1 + w_xor_block[2] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j1=0; j1 < cDDT_n[carry[1]][AB_block[1]][i1]; j1++)
			{
				gamma_bloc[1] = cDDT_v[carry[1]][AB_block[1]][i1][j1];
				//gamma_block_tmp[1] = gamma_bloc[1] << 8;
				//carry[2] = ((alpha_bloc[1] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[1] >> 7) << 1)   //beta MSB
				carry[2] = carry_tmp[2] + (gamma_bloc[1] >> 7); // gamma MSB

			for(i2 = cDDT_wt_min[carry[2]][AB_block[2]];
					i2 <= cDDT_wt_min[carry[2]][AB_block[2]]; i2++)  //cDDT_wt_max[carry[2]][AB_block[2]]
			{
				w2 = w1 + i2;
				if(w2 + w_xor_block[3] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j2=0; j2 < cDDT_n[carry[2]][AB_block[2]][i2]; j2++)
			{
				gamma_bloc[2] = cDDT_v[carry[2]][AB_block[2]][i2][j2];
				//gamma_block_tmp[2] = gamma_bloc[2] << 16;
				//carry[3] = ((alpha_bloc[2] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[2] >> 7) << 1)   //beta MSB
				carry[3] = carry_tmp[3] + (gamma_bloc[2] >> 7); // gamma MSB
			for(i3 = cDDT_wt_min[carry[3]][AB_block[3]];
					i3 <= cDDT_wt_min[carry[3]][AB_block[3]]; i3++)  //cDDT_wt_max[carry[3]][AB_block[3]]
			{
				w3 = w2 + i3;
				if(w3 + w_xor_block[4] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j3=0; j3 < cDDT_n[carry[3]][AB_block[3]][i3]; j3++)
			{
				gamma_bloc[3] = cDDT_v[carry[3]][AB_block[3]][i3][j3];
				//gamma_block_tmp[3] = gamma_bloc[3] << 24;
				//carry[4] = ((alpha_bloc[3] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[3] >> 7) << 1)   //beta MSB
				carry[4] = carry_tmp[4] + (gamma_bloc[3] >> 7); // gamma MSB

			for(i4 = cDDT_wt_min[carry[4]][AB_block[4]];
					i4 <= cDDT_wt_min[carry[4]][AB_block[4]]; i4++)  //cDDT_wt_max[carry[4]][AB_block[4]]
			{
				w4 = w3 + i4;
				if(w4 + w_xor_block[5] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j4=0; j4 < cDDT_n[carry[4]][AB_block[4]][i4]; j4++)
			{
				gamma_bloc[4] = cDDT_v[carry[4]][AB_block[4]][i4][j4];
				//gamma_block_tmp[4] = gamma_bloc[4] << 32;
				//carry[5] = ((alpha_bloc[4] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[4] >> 7) << 1)   //beta MSB
				carry[5] = carry_tmp[5] + (gamma_bloc[4] >> 7); // gamma MSB

			for(i5 = cDDT_wt_min[carry[5]][AB_block[5]];
					i5 <= cDDT_wt_min[carry[5]][AB_block[5]]; i5++)  //cDDT_wt_max[carry[5]][AB_block[5]]
			{
				w5 = w4 + i5;
				if(w5 + w_xor_block[6] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j5=0; j5 < cDDT_n[carry[5]][AB_block[5]][i5]; j5++)
			{
				gamma_bloc[5] = cDDT_v[carry[5]][AB_block[5]][i5][j5];
				//gamma_block_tmp[5] = gamma_bloc[5] << 40;
				//carry[6] = ((alpha_bloc[5] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[5] >> 7) << 1)   //beta MSB
				carry[6] = carry_tmp[6] + (gamma_bloc[5] >> 7); // gamma MSB

			for(i6 = cDDT_wt_min[carry[6]][AB_block[6]];
					i6 <= cDDT_wt_min[carry[6]][AB_block[6]]; i6++)  //cDDT_wt_max[carry[6]][AB_block[6]]
			{
				w6 = w5 + i6;
				if(w6 + w_xor_block[7] > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
			for(j6=0; j6 < cDDT_n[carry[6]][AB_block[6]][i6]; j6++)
			{
				gamma_bloc[6] = cDDT_v[carry[6]][AB_block[6]][i6][j6];
				//gamma_block_tmp[6] = gamma_bloc[6] << 48;
				//carry[7] = ((alpha_bloc[6] >> 7) << 2) //alpha MSB
				//		+ ((beta_bloc[6] >> 7) << 1)   //beta MSB
				carry[7] = carry_tmp[7] + (gamma_bloc[6] >> 7); // gamma MSB

			for(i7 = MSB_cDDT_wt_min[carry[7]][AB_block[7]];
					i7 <= MSB_cDDT_wt_min[carry[7]][AB_block[7]]; i7++)  //MSB_cDDT_wt_max[carry[7]][AB_block[7]]
			{
				w7 = w6 + i7;
				if(w7 > w_cmp ){break;}  //break 要比 continue 高效,直接结束.
				P_w[search_round] = w7; //Speck_XDP_compute(Alpha, Beta, gamma_temp);
			for(j7=0; j7 < MSB_cDDT_n[carry[7]][AB_block[7]][i7]; j7++)
			{
				gamma_bloc[7] = MSB_cDDT_v[carry[7]][AB_block[7]][i7][j7];
//								gamma_temp = (gamma_bloc[7] << 56) | gamma_block_tmp[6]
//								| gamma_block_tmp[5] | gamma_block_tmp[4]
//								| gamma_block_tmp[3] | gamma_block_tmp[2]
//								| gamma_block_tmp[1] | gamma_bloc[0];
//
				gamma_temp = (gamma_bloc[7] << 56) | (gamma_bloc[6] << 48)
						| (gamma_bloc[5] << 40) | (gamma_bloc[4] << 32)
						| (gamma_bloc[3] << 24) | (gamma_bloc[2] << 16)
						| (gamma_bloc[1] <<8) | gamma_bloc[0];

				if (P_w[search_round] == w_cmp) //概率满足Matsui条件
				{
	#if 1   //是否找到第一条最优路径就返回？
			//if( best == 1 )
			{
				best  = 1;
				n_P_bestofR_w[search_round] = Bn_w;
				n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
				n_y_in[search_round +1] = gamma_temp ^ next_Beta_tmp;
				return 1;
			}
	#endif
				}
			}}}}}}}}}}}}}}}}

/*
//////////////贪婪算法////////////每一轮的都只是概率最大的分支,该相同概率下的分支会有很多,但是这里只搜了其中一条,因此,遍历
							///的总的分支数就是第1轮的输入组合分支数

							//xor_alpha_beta = Beta ^ Alpha ;
							//w_xor = HM_weight(xor_alpha_beta & ValueMax_Align);
							and_alpha_beta = Beta & Alpha;
							w_and = HM_weight(and_alpha_beta & ValueMax_Align);

							P_w[search_round] = w_xor + w_and;  // 00-->0;01,10-->1; 11-->0,概率的和
							if ( P_w[search_round] == w_cmp)
							{

							gamma_temp = xor_alpha_beta; // 00-->1; 01,10--->0;11-->1;的情况先不考虑,因此只贪婪的考虑SetA的情况.

							//best = sepck_round_r(search_round,cur_round +1,gamma_temp,Beta);
					#if 1   //是否找到第一条最优路径就返回？
						//if( best == 1 )
						{
							best  = 1;
							n_P_bestofR_w[search_round] = Bn_w;
							n_x_in[search_round +1] = gamma_temp;  // output part of the r-th round..
							n_y_in[search_round +1] = gamma_temp ^ (ROTATE_LEFT(Beta, rol_b,blocksize_len) & Bit_Align);
							return 1;
						}
					#endif
//////////////贪婪算法////////////
*/
#endif
return  best;
}
/////////
#endif






/**
 * @Print the final result of RXD in speck key schedule part.
 *
 */
void print_RX_speck_key_resoult(u16 search_round)
{
	u16 r_num = 0;

	switch(blocksize_len)
	{
	case 16:
		printf("round----left-----right----weight \n");  //----gama_cnt \n");
		for(r_num=1;r_num <= search_round;r_num++)
		{
			printf("%02d     0x%04llx    0x%04llx     -%d \n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num]); //,gama_cnt[r_num]);
			fprintf(simon_diff_trail,"%02d  \t 0x%04llx \t 0x%04llx \t -%d \t -%d \n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num],n_P_bestofR_w[r_num]);//,gama_cnt[r_num]); // write.
		}
		printf("%02d     0x%04llx    0x%04llx     NULL \n",(search_round+1),n_x_in[search_round+1],n_y_in[search_round+1]);
		fprintf(simon_diff_trail,"%02d  \t 0x%04llx \t 0x%04llx \t NULL \t NULL\n",search_round+1,n_x_in[search_round+1],n_y_in[search_round+1]); // write.
		break;
	case 24:
		printf("round----left--------right----weight \n");
		for(r_num=1;r_num <= search_round;r_num++)
		{
			printf("%02d     0x%06llx    0x%06llx     -%d \n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num]);
			fprintf(simon_diff_trail,"%02d  \t 0x%06llx \t 0x%06llx \t -%d \t -%d\n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num],n_P_bestofR_w[r_num]); // write.
		}
		printf("%02d     0x%06llx    0x%06llx     NULL \n",(search_round+1),n_x_in[search_round+1],n_y_in[search_round+1]);
		fprintf(simon_diff_trail,"%02d  \t 0x%06llx \t 0x%06llx \t NULL \t NULL\n",search_round+1,n_x_in[search_round+1],n_y_in[search_round+1]); // write.
		break;
	case 32:
		printf("round-----left---------right------weight \n");
		for(r_num=1;r_num <= search_round;r_num++)
		{
			printf("%02d     0x%08llx    0x%08llx     -%d \n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num]);
			fprintf(simon_diff_trail,"%02d  \t 0x%08llx \t 0x%08llx \t -%d \t -%d\n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num],n_P_bestofR_w[r_num]); // write.
		}
		printf("%02d     0x%08llx    0x%08llx     NULL \n",(search_round+1),n_x_in[search_round+1],n_y_in[search_round+1]);
		fprintf(simon_diff_trail,"%02d  \t 0x%08llx \t 0x%08llx \t NULL \t NULL\n",search_round+1,n_x_in[search_round+1],n_y_in[search_round+1]); // write.
		break;
	case 48:
		printf("round-------left-------------right--------weight \n");
		for(r_num=1;r_num <= search_round;r_num++)
		{
			printf("%02d     0x%012llx    0x%012llx     -%d \n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num]);
			fprintf(simon_diff_trail,"%02d  \t 0x%012llx \t 0x%012llx \t -%d \t -%d\n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num],n_P_bestofR_w[r_num]); // write.
		}
		printf("%02d     0x%012llx    0x%012llx     NULL \n",(search_round+1),n_x_in[search_round+1],n_y_in[search_round+1]);
		fprintf(simon_diff_trail,"%02d  \t 0x%012llx \t 0x%012llx \t NULL \t NULL\n",search_round+1,n_x_in[search_round+1],n_y_in[search_round+1]); // write.
		break;
	default:  //64
		printf("round---------left---------------right-------------weight \n");
		for(r_num=1;r_num <= search_round;r_num++)
		{
			printf("%02d     0x%016llx    0x%016llx     -%d \n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num]);
			fprintf(simon_diff_trail,"%02d  \t 0x%016llx \t 0x%016llx \t -%d \t -%d\n",r_num,n_x_in[r_num],n_y_in[r_num],n_P_w[r_num],n_P_bestofR_w[r_num]); // write.
		}
		printf("%02d     0x%016llx    0x%016llx     NULL \n",(search_round+1),n_x_in[search_round+1],n_y_in[search_round+1]);
		fprintf(simon_diff_trail,"%02d  \t 0x%016llx \t 0x%016llx \t NULL \t NULL\n",search_round+1,n_x_in[search_round+1],n_y_in[search_round+1]); // write.
		break;
	}

	printf("%d Round Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
}






