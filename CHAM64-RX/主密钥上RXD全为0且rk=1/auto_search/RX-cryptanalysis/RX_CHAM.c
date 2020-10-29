/*
 * RX_CHAM.c
 *
 *  Created on: 2020年9月8日
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
#include "RX_CHAM.h"
#include "CHAM.h"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
volatile u16 CHAM64_const_rxd[80]={0}; //每轮的轮常数的RXD
volatile u32 CHAM128_const_rxd[80]={0}; //每轮的轮常数的RXD

volatile float RX_n_P_bestofR_w[80] = {0};
volatile float CHAM_RX_Pw_even_wt[80] = {0};

volatile u32 CHAM_RXD[4][80]={0}; //每轮的输入差分
volatile u32 CHAM_Zeta[80]={0}; //每轮的输zeta
volatile u32 CHAM_zeta_tmp[80]={0}; //每轮的临时zeta
volatile float CHAM_Pr_zeta[80]={0}; //每轮RX差分概率
volatile float CHAM_Pr_zeta_tmp[80]={0}; //每轮RX差分临时概率
volatile float CHAM_Pr_RX[80]={0}; //每轮RX差分概率
volatile float CHAM_Pr_add[80]={0}; //每轮RX差分概率




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
void Compute_CHAM64_const_rxd(void)
{
	u16 i = 0;

	for(i=0; i<80; i++)
	{
		CHAM64_const_rxd[i+1] = i ^ (ROTATE_LEFT(i, RX_k, 16));
	}
}
void Compute_CHAM128_const_rxd(void)
{
	u16 i = 0;

	for(i=0; i<80; i++)
	{
		CHAM128_const_rxd[i+1] = i ^ (ROTATE_LEFT(i, RX_k, 32));
	}
}


/////////////搜索CHAM的固定密钥中RXD固定为0的 RX 特征的代码//////////////////////////////
u32 CHAM_64_RX_trail_search_entry(u16 search_round)
{
	u32 best = 0;
	u16 i = 0;
	FILE* CHAM_64_RX_Bn;
	FILE* CHAM_64_RX_Bn_even;
	u16 Bn_wt_expected = 0;

	CHAM_64_RX_Bn = fopen ("../tmp/CHAM_64_RX_Bn_wt.xlsx", "a+"); //  "w+"); //
	CHAM_64_RX_Bn_even = fopen ("../tmp/CHAM_64_RX_Bn_wt_even.xlsx", "a+"); //  "w+"); //

	if(sc_blocksize != 64)
	{
		printf("The block size of CHAM64 should be 64 bits. \n");
		return 0;
	}
	// For CHAM-64 modular additon is 16 bits.
	blocksize_len = 16;
	nBytes = 2;
	Bit_Align = 0xFFFF;
	ValueMax_Align = 0x7FFF;
	V_MSB = 0x8000;

	Compute_CHAM64_const_rxd(); //提前计算各个轮常数对应的RXD

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
		if(search_round > 1)
		{
			for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
			{
				fscanf(CHAM_64_RX_Bn,"%f",&RX_n_P_bestofR_w[i]);  //read.
				fscanf(CHAM_64_RX_Bn_even,"%f",&CHAM_RX_Pw_even_wt[i]);  //read.
			}
		}
	//-------记录时间-----
		time_xor_talbe = clock();
		run_time =  (time_xor_talbe - time_start) / CLOCKS_PER_SEC;
		printf("Time of Preprocessing info: %.2f seconds.  \n", run_time);
		printf("|--------------------------------------------------|\n");
	//%%%%%%%%%%%%%%%%%%%%%%%%%%构造cDDT%%%%%%%%%%%%%%%%%%%%%%%%//
		printf("Constructing cDDT Lookup Tables ... \n");
		//构造cDDT
		ARX_carry_DDTm_construct();  // FOR cDDT
		fixed_Alpha_get_betagamma();   //for 234 round
		//fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
		time_ARX_DDT = clock();
		run_time =  (time_ARX_DDT - time_xor_talbe) / CLOCKS_PER_SEC;
		printf("Time of Construct CHAM cDDTs: %.2f seconds.  \n", run_time);

//%%%%%%%%%%%%%%%%%判断开始的奇偶轮%%%%%%%%%%%%%%%%%%%%%%%//
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径,=0为奇数轮开始的正常搜索
	{
		printf("===============Search odd-round============= \n");
		//Bn_w = n_P_bestofR_w[search_round -1] - 1;
		//	Bn_w = 0;
			p_sumof_r_P_w[search_round] = 0xFFFF;

			int RX_Expected_w_int = 0;
				RX_Expected_w_int = (int)RX_n_P_bestofR_w[search_round - 1];
				RX_Expected_w = RX_Expected_w_int;
				RX_Expected_w = 63;  ////!!!!!启发式，可以通过设定直接的开始期望概率重量!!!!!//
	}
#else
	{
		printf("===============Search even-round============= \n");
		//Bn_w = CHAM_Pw_even_wt[search_round -1] - 1;
		//	Bn_w = 0;
			p_sumof_r_P_w[search_round] = 0xFFFF;

			int RX_Expected_w_int = 0;
				RX_Expected_w_int = (int)CHAM_RX_Pw_even_wt[search_round - 1];
				RX_Expected_w = RX_Expected_w_int;
				RX_Expected_w = 63;  ////!!!!!启发式，可以通过设定直接的开始期望概率重量!!!!!//
	}
#endif

printf("*****************************Searching processing**************************|\n");
//	Bn_w = 63;  //直接设定开始的概率重量,快速达到期望的概率重量
			do
			{
				RX_Expected_w = RX_Expected_w + 1; // 从r-1轮的概率重量整数部分+1开始
				printf("Update: %d   RX_Bn_w: %f  --Max expected weight-- \n",
						Num_Bn_Update, RX_Expected_w);
				// Search Entry.
				best = CHAM_RX_trail_round_1(search_round); //效率高
/*
				u32 Output_X[4] = {0};
				Output_X[0] = 0x0010;
				Output_X[1] = 0x1020;
				Output_X[2] = 0x2800;
				Output_X[3] = 0x0000;
				P_w[1] = 1;  // P_w[2] = 2; P_w[3] = 3;   //pro;
				p_sumof_r_P_w[1] = P_w[1]; //  p_sumof_r_P_w[2] = 3; p_sumof_r_P_w[3] = 6;
				best = CHAM_Diff_trail_round_23(search_round, 2, Output_X);  //第1轮的输出
	//			best = CHAM_Diff_trail_round_r(search_round, 4, Output_X); //r lun //第3轮的输出
*/

				time_Round = clock();
				run_time =  (time_Round - time_ARX_cLAT) / CLOCKS_PER_SEC;
				printf("Time: %.2f seconds.  \n", run_time);
/*
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
				{
					//Bn_wt_expected =  p_sumof_r_P_w[search_round];
					Bn_wt_expected =  n_P_bestofR_w[search_round];
				}
#else
				{
					//Bn_wt_expected =  p_sumof_r_P_w[search_round];
					Bn_wt_expected =  CHAM_Pw_even_wt[search_round];
				}
#endif
				*/
			}while( Num_Bn_Update == 0 );   //end conditions. // 概率重量逼近次数大于1，且返回值表示搜索到最优路径
			printf("**************************************************************************|\n");
			print_CHAM_64_RXD_resoult(search_round);  //打印CHAM的RXDT

			time_finish = clock();
			run_time = (double)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
			printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours == %.2f days. \n", run_time,run_time/60.0,run_time/3600.0,run_time/86400.0 );
			printf("Auto-search CHAM-64 RX-differential trails END! \n");
			printf("|************************************************************************|\n");
			fclose(CHAM_64_RX_Bn);
			fclose(CHAM_64_RX_Bn_even);
	return best;
}

////////////////////////////////
u32 CHAM_128_RX_trail_search_entry(u16 search_round)
{
	u32 best = 0;
	u16 i = 0;
	FILE* CHAM_128_RX_Bn;
	FILE* CHAM_128_RX_Bn_even;
	u16 Bn_wt_expected = 0;


	CHAM_128_RX_Bn = fopen ("../tmp/CHAM_128_RX_Bn_wt.xlsx", "a+"); //  "w+"); //
	CHAM_128_RX_Bn_even = fopen ("../tmp/CHAM_128_RX_Bn_wt_even.xlsx", "a+"); //  "w+"); //
	if(sc_blocksize != 128)
	{
		printf("The block size of CHAM128 should be 128 bits. \n");
		return 0;
	}
	// For CHAM-128 modular additon is 32 bits.
	blocksize_len = 32;
	nBytes = 4;
	Bit_Align = 0xFFFFFFFF;
	ValueMax_Align = 0x7FFFFFFF;
	V_MSB = 0x80000000;

	Compute_CHAM128_const_rxd();



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
		if(search_round > 1)
		{
			for(i=1; i < search_round; i++)  //搜索n轮，需要先搜索前n-1轮得到P_bestofR_w。
			{
				fscanf(CHAM_128_RX_Bn,"%f",&n_P_bestofR_w[i]);  //read.
				fscanf(CHAM_128_RX_Bn_even,"%f",&CHAM_Pw_even_wt[i]);  //read.
			}
		}
		//-------记录时间-----
			time_xor_talbe = clock();
			run_time =  (time_xor_talbe - time_start) / CLOCKS_PER_SEC;
			printf("Time of Preprocessing info: %.2f seconds.  \n", run_time);
			printf("|--------------------------------------------------|\n");
		//%%%%%%%%%%%%%%%%%%%%%%%%%%构造cDDT%%%%%%%%%%%%%%%%%%%%%%%%//
			printf("Constructing cDDT Lookup Tables ... \n");
			//构造cDDT
			ARX_carry_DDTm_construct();  // FOR cDDT
			fixed_Alpha_get_betagamma();   //for 234 round
			//fixed_Alpha_get_betagamma();   //对于模加一个输入差分已知的情况，构造这类cDDT
			time_ARX_DDT = clock();
			run_time =  (time_ARX_DDT - time_xor_talbe) / CLOCKS_PER_SEC;
			printf("Time of Construct CHAM cDDTs: %.2f seconds.  \n", run_time);

	//%%%%%%%%%%%%%%%%%判断开始的奇偶轮%%%%%%%%%%%%%%%%%%%%%%%//
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径,=0为奇数轮开始的正常搜索
	{
		printf("===============Search odd-round============= \n");
		//Bn_w = n_P_bestofR_w[search_round -1] - 1;
		//	Bn_w = 0;
			p_sumof_r_P_w[search_round] = 0xFFFF;

			int RX_Expected_w_int = 0;
				RX_Expected_w_int = (int)RX_n_P_bestofR_w[search_round - 1];
				RX_Expected_w = RX_Expected_w_int;
				RX_Expected_w = 42;  ////!!!!!启发式，可以通过设定直接的开始期望概率重量!!!!!//
	}
#else
	{
		printf("===============Search even-round============= \n");
		//Bn_w = CHAM_Pw_even_wt[search_round -1] - 1;
		//	Bn_w = 0;
			p_sumof_r_P_w[search_round] = 0xFFFF;

			int RX_Expected_w_int = 0;
				RX_Expected_w_int = (int)CHAM_RX_Pw_even_wt[search_round - 1];
				RX_Expected_w = RX_Expected_w_int;
				RX_Expected_w = 42;  ////!!!!!启发式，可以通过设定直接的开始期望概率重量!!!!!//
	}
#endif

	printf("*****************************Searching processing**************************|\n");
//	  	Bn_w = 1;  //直接设定开始的概率重量,快速达到期望的概率重量
			do
			{
				RX_Expected_w = RX_Expected_w + 1; // 从r-1轮的概率重量整数部分+1开始
				printf("Update: %d   RX_Bn_w: %f  --Max expected weight-- \n",
						Num_Bn_Update, RX_Expected_w);
				// Search Entry.
				best = CHAM_RX_trail_round_1(search_round); //效率高
/*
				u32 Output_X[4] = {0};
				Output_X[0] = 0x10200000;
				Output_X[1] = 0x08002000;
				Output_X[2] = 0x00000000;
				Output_X[3] = 0x40000000;
				P_w[1] = 3;   P_w[2] = 1;  P_w[3] = 2;   //pro;
				p_sumof_r_P_w[1] = P_w[1];   p_sumof_r_P_w[2] = 4;  p_sumof_r_P_w[3] = 6;
	//			best = CHAM_Diff_trail_round_23(search_round, 2, Output_X);  //第1轮的输出
				best = CHAM_Diff_trail_round_r(search_round, 4, Output_X); //r lun //第3轮的输出
*/

				time_Round = clock();
				run_time =  (time_Round - time_ARX_cLAT) / CLOCKS_PER_SEC;
				printf("Time: %.2f seconds.  \n", run_time);
/*
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
				{
					Bn_wt_expected = n_P_bestofR_w[search_round];
				}
#else
				{
					Bn_wt_expected = CHAM_Pw_even_wt[search_round];
				}
#endif
*/
			}while( Num_Bn_Update == 0 );   //end conditions. // 概率重量逼近次数大于1，且返回值表示搜索到最优路径
			printf("**************************************************************************|\n");
			print_CHAM_128_RXD_resoult(search_round);  //打印CHAM的RXDT

			time_finish = clock();
			run_time = (double)(time_finish - time_ARX_DDT) / CLOCKS_PER_SEC;
			printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours == %.2f days. \n", run_time,run_time/60.0,run_time/3600.0,run_time/86400.0 );
			printf("Auto-search CHAM-64 RX-differential trails END! \n");
			printf("|************************************************************************|\n");
			fclose(CHAM_128_RX_Bn);
			fclose(CHAM_128_RX_Bn_even);
	return best;
}




u32 CHAM_RX_trail_round_1(u16 search_round)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 thr_d = 0;
	u32 Input_alpha = 0;
	u32 Input_beta = 0;
	u32 Input_gamma = 0;
	u32 Out_RXD_tmp = 0;
	u16 ROL_a =1, ROL_b = 8;  //当搜索过程中间轮为 奇数 轮，默认旋转参数为1,8
	u32 Output_X[4] = {0};

	u32 M0 = 0;
	u64 N0 = blocksize_len-1;
	char A0[65], T0[65], F0[65], H0[65], C0[65];
	volatile char X0, Y0, I0, L0, Z0;
	int i=0;
	u16 j_abc = 0;

	u32 zeta_num = 0;
	u32 updata_cur_num = 0;

	updata_cur_num = Num_Bn_Update;

#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
	{
		ROL_a =1; ROL_b = 8;
	}
#else
	{
		ROL_a =8, ROL_b = 1;
	}
#endif
///////////////////round 1 entry/////////////////////
	if( search_round == 1)
	{
		return best;
	}
	else
	{
		// 可以选择是否限定第1轮的概率重量
		//printf("**Noted**: for(thr_d = 0; thr_d < ?; thr_d++) \n"); //打印启发信息
		for(thr_d = 0; thr_d < blocksize_len-1; thr_d++)  //0::n-1 // blocksize_len-1
		{
			// thr_d : the bits number of alpha/beta/gamma that the position not equal to each other at the same time.
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{	//奇数开始轮，则链接偶数轮
				if ((thr_d +Ek_min + CHAM_RX_Pw_even_wt[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{	//统一退出位置
					return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
				}
			}
#else
			{	//偶数轮开始，则链接奇数轮
				if ((thr_d +Ek_min + RX_n_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{	//统一退出位置
					return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
				}
			}
#endif


			M0 = thr_d;  //三个alpha,beta和gamma同时考虑
			P_w[1] = M0;  //pro;
			//p_sumof_r_P_w[1] = P_w[1];
			if(M0 == 0)  /// abd [0 -> n-2] 没有完全不相同的比特,且MSB之后只能为0
			{
				for(j_abc = 0; j_abc < 4;j_abc++)  //*******注意在CHAM中首轮不能为全0的输入输出差分******
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


					for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
					{
						RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{	//奇数开始轮，则链接偶数轮
				if ((RX_P_w[1]+ CHAM_RX_Pw_even_wt[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{	//统一退出位置
					return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
				}
			}
#else
			{	//偶数轮开始，则链接奇数轮
				if ((RX_P_w[1]+ RX_n_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{	//统一退出位置
					return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
				}
			}
#endif

			RX_P_sumofR_w[1] = RX_P_w[1];
			CHAM_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
			CHAM_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

			Output_X[0] = ROTATE_RIGHT(Input_beta, ROL_a, blocksize_len) & Bit_Align;
			Out_RXD_tmp = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]];
			Output_X[3] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;

			CHAM_X[0][1] = Input_alpha;
			CHAM_X[1][1] = Output_X[0];

			best = CHAM_RX_trail_round_23(search_round, 2, Output_X);

#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		if(updata_cur_num < Num_Bn_Update)
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知已经更新
		////更新期望概率重量RX_Expected_w，逐渐逼近,搜首轮其它概率重量%%%%%%%%%%%%%%%%%
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{	//奇数开始轮，则链接偶数轮
				if ((RX_P_w[1]+ CHAM_RX_Pw_even_wt[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{
					updated = 2;
					break;
				}
			}
#else
			{	//偶数轮开始，则链接奇数轮
				if ((RX_P_w[1]+ RX_n_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{
					updated = 2;
					break;
				}
			}
#endif
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

		best = CHAM_RX_input_MSB(search_round,Input_alpha, Input_beta, Input_gamma,thr_d,C0);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一轮有更新
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
								goto cham_r1;
							}
							if (L0 == N0)
							{
								T0[F0[N0]] = -I0 - 1;
								T0[I0 + 1] = T0[I0];
								I0 = F0[N0];
								F0[N0] = F0[N0] + 1;
								goto cham_r1;
							}
							T0				[L0] = -I0-1;
							T0[I0+1] = T0[I0];
							F0[L0] = F0[L0] + 1;
							I0 = L0;
							goto cham_r1;
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
									goto cham_r1;
								}
								T0[F0[N0]-1] = -I0-1;
								T0[I0+1] = T0[I0];
								I0 = F0[N0] - 1;
								goto cham_r1;
							}
							T0[L0] = -I0 -1; T0[I0 + 1] = T0[I0]; I0 = L0; goto cham_r1;
						}
						X0 = N0;
						F0[L0 - 1] = F0[L0];
						F0[N0] = N0;
						L0 = N0;
						if (I0 == N0 - 1)
						{
							I0 = T0[N0 - 1];
							goto cham_r1;
						}
						T0[N0 - 1] = -I0 - 1;
						T0[I0 + 1] = T0[I0];
						I0 = N0 - 1;
			cham_r1:
					A0[X0] = 1;
					A0[Y0] = 0;
					H0[X0] = Z0 = H0[Y0];
					C0[Z0] = X0;
					}
/////////////////////////////////
			best = CHAM_RX_input_MSB(search_round,Input_alpha, Input_beta, Input_gamma,thr_d,C0);
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
	}
	return best;
}



u16 CHAM_RX_input_MSB
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

			state = CHAM_RX_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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

		state = CHAM_RX_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,P_1);
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
u16 CHAM_RX_input_Middle
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

			state = CHAM_RX_input_Middle(search_round,Input_alpha,Input_beta, Input_gamma,P_1,posi,cur_posi-1);
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
		state = CHAM_RX_input_Last(search_round,
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

u16 CHAM_RX_input_Last
(u16 search_round,u64 alpha, u64 beta, u64 gamma,u16 P_1, char *posi )
{
	u16 state = 0,updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u64 Input_alpha = 0;
	u64 Input_beta = 0;
	u64 Input_gamma = 0;
	u32 Out_RXD_tmp = 0;
	u32 Output_X[4] = {0};
	u16 ROL_a =1, ROL_b = 8;  //当搜索过程中间轮为 奇数 轮，默认旋转参数为1,8

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

		//基于构造的第1轮的输入输出差分组合，调用下一轮
		for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
		{
			RX_P_w[1] = P_w[1] + RX_zeta_w[zeta_num];
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{	//奇数开始轮，则链接偶数轮
				if ((RX_P_w[1]+ CHAM_RX_Pw_even_wt[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{	//统一退出位置
					return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
				}
				ROL_a =1; ROL_b = 8;
			}
#else
			{	//偶数轮开始，则链接奇数轮
				if ((RX_P_w[1]+ RX_n_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{	//统一退出位置
					return updated;  //Return the upper procedure::sepck_differential_trail_search_entry
				}
				ROL_a =8, ROL_b = 1;
			}
#endif
			RX_P_sumofR_w[1] = RX_P_w[1];
			CHAM_Pr_zeta_tmp[1] = RX_zeta_w[zeta_num];
			CHAM_zeta_tmp[1] = RX_zeta[RX_zeta_i[zeta_num]];

			Output_X[0] = ROTATE_RIGHT(Input_beta, ROL_a, blocksize_len) & Bit_Align;
			Out_RXD_tmp = Input_gamma ^ RX_zeta[RX_zeta_i[zeta_num]];
			Output_X[3] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;

			CHAM_X[0][1] = Input_alpha;
			CHAM_X[1][1] = Output_X[0];

			state = CHAM_RX_trail_round_23(search_round, 2, Output_X);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(state != 0 )  //判断是结束搜索，还是进行遍历其它组合输入输出XOR差分
	{
		if(updata_cur_num < Num_Bn_Update)
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知已经更新
		////更新期望概率重量RX_Expected_w，逐渐逼近,搜首轮其它概率重量%%%%%%%%%%%%%%%%%
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
			{	//奇数开始轮，则链接偶数轮
				if ((RX_P_w[1]+ CHAM_RX_Pw_even_wt[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{
					updated = 2;
					break;
				}
			}
#else
			{	//偶数轮开始，则链接奇数轮
				if ((RX_P_w[1]+ RX_n_P_bestofR_w[search_round -1]) >= RX_Expected_w) // Pmax_compute(xi) is the ri round min weight.
				{
					updated = 2;
					break;
				}
			}
#endif
	}}
#endif
		}
	}
	return state;
}


u32 CHAM_RX_trail_round_23(u16 search_round, u16 cur_round, u32 *input_x)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 ROL_a =1, ROL_b = 8;  //当搜索过程中间轮为 奇数 轮，默认旋转参数为1,8
	u32 beta_tmp=0, gamma_tmp=0;
	u32 Input_RX0 = 0;
	u32 Out_RXD_tmp = 0;
	u32 Output_X[4] = {0};
	u32 beta_block[4] = {0};
	u32 gama_block[4] = {0};
	u8 carry_block[4] = {0};
	u16 X0_block[4] = {0};
	u16 X0_block_msb[4] = {0};
	u8 X0_block_wt_min[4] = {0};
	u16 i=0, j=0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	float w_cmp = 0;
	u16 i0=0,i1=0,i2=0,i3=0;
	u16 j0=0,j1=0,j2=0,j3=0;
	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制

	updata_cur_num = Num_Bn_Update;



#if (CHAM_even_round == 0 ) //当前搜索参数（1，8）开始的奇数轮的最优路径
	if((cur_round & 0x01) == 0)  //当前偶数轮%2
	{
		ROL_a = 8; ROL_b = 1;

		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
	else  // 奇数轮
	{
		//ROL_a =1; ROL_b = 8; //默认的值就是1,8
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
#else  //CHAM_even_round=1, 当前搜索参数（8,1）开始的偶数轮的最优路径
	if((cur_round & 0x01) == 0)  //偶数轮
	{
		//ROL_a = 1; ROL_b = 8;  //默认的值就是1,8
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
	else  // 奇数轮
	{
		ROL_a = 8; ROL_b = 1;

		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
#endif



#if (CHAM_64_or_128 == 64)
		Input_RX0 = input_x[0]^ CHAM64_const_rxd[cur_round];
		X0_block[0] = Input_RX0 & 0xFF;   //取出每个8比特
		X0_block[1] = (Input_RX0 >> 8) & 0xFF;   //取出每个8比特

		//确定每个输入差分alpha的每8比特的子块所对应的可能的最小差分概率重量
		X0_block_wt_min[0] = CHAM_a_wt_min[X0_block[0]][0];
		X0_block_wt_min[1] = 0xFF;
		if( (X0_block[0] & 0x80) == 0)
		{
			for(i=0;i<4;i++)  //carry 为 0/1/2/3  carry=alpha||beta||gamma,即alpha的这位为0
			{
				if(X0_block_wt_min[1] > msb_CHAM_a_wt_min[X0_block[1]][i])
				{
					X0_block_wt_min[1] = msb_CHAM_a_wt_min[X0_block[1]][i];
				}
			}
		}
		else
		{
			for(i=4;i<8;i++) //carry 为 4/5/6/7  carry=alpha||beta||gamma,即alpha的这位为1
			{
				if(X0_block_wt_min[1] > msb_CHAM_a_wt_min[X0_block[1]][i])
				{
					X0_block_wt_min[1] = msb_CHAM_a_wt_min[X0_block[1]][i];
				}
			}
		}

		if( X0_block_wt_min[0] + X0_block_wt_min[1] > w_cmp)
		{return 0;}
/*
	if( _mm_popcnt_u64(input_x[0] & 0x7FFF) > w_cmp)
	{return 0;}
	X0_block_wt_min[1] = _mm_popcnt_u64(input_x[0] & 0x7F);
*/

		X0_block_msb[0] = ((X0_block[0] >> 7) << 2);  //alpha 最低8比特的 MSB
		//整理给下一轮的输入差分
		Output_X[1] = input_x[2];
		Output_X[2] = input_x[3];

		//block 0  //可以控制每个block的概率重量范围, 例如[min, min +1]
		for(i0=CHAM_a_wt_min[X0_block[0]][0];
				i0<=CHAM_a_wt_min[X0_block[0]][0] +1; i0++)  //X0_block_wt_min[0] // <=8 //CHAM_a_wt_max[X0_block[0]][0]
		{
			if(i0 + X0_block_wt_min[1] > w_cmp ){break;}
			for(j0=0; j0<CHAM_a_bg_numb[X0_block[0]][0][i0]; j0++)
			{
				beta_block[0] = CHAM_a_beta[X0_block[0]][0][i0][j0];
				gama_block[0] = CHAM_a_gama[X0_block[0]][0][i0][j0];
				carry_block[0] = X0_block_msb[0] +
						((beta_block[0] >> 7) << 1) + (gama_block[0] >> 7);

		//block 1
		for(i1=msb_CHAM_a_wt_min[X0_block[1]][carry_block[0]];
				i1<=msb_CHAM_a_wt_min[X0_block[1]][carry_block[0]] +1; i1++) //X0_block_wt_min[1]  //<8//msb_CHAM_a_wt_max[X0_block[1]][carry_block[0]]
			{
				P_w[cur_round] = i0 + i1 ;  //// 当前轮的概率重量
				if(P_w[cur_round] > w_cmp ){break;}
				//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
				for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
					CHAM_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
					CHAM_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
					RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
					RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

				for(j1=0; j1<msb_CHAM_a_bg_numb[X0_block[1]][carry_block[0]][i1]; j1++)
				{
					beta_block[1] = msb_CHAM_a_beta[X0_block[1]][carry_block[0]][i1][j1];
					gama_block[1] = msb_CHAM_a_gama[X0_block[1]][carry_block[0]][i1][j1];

					beta_tmp = ((beta_block[1] << 8) | beta_block[0]) & Bit_Align;
					gamma_tmp = ((gama_block[1] << 8) | gama_block[0]) & Bit_Align;

					//整理给下一轮的输入差分
					Output_X[0] = ROTATE_RIGHT(beta_tmp, ROL_a, blocksize_len) & Bit_Align;
					Out_RXD_tmp = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]];
					Output_X[3] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;

					if(cur_round == 2 ) //记录第2轮的输入
					{
						CHAM_X[2][1] = Output_X[0];

						CHAM_X[0][cur_round] = input_x[0];
						CHAM_X[1][cur_round] = Output_X[0];
						//CHAM_X[2][cur_round] = CHAM_X[1][3];
						CHAM_X[3][cur_round] = input_x[3];

					//	Output_X[2] = input_x[3];
					}
					if(cur_round == 3 ) //记录第3轮的输出
					{
						CHAM_X[3][1] = Output_X[0];

						CHAM_X[2][2] = Output_X[0];

						CHAM_X[0][cur_round] = input_x[0];
						CHAM_X[1][cur_round] = Output_X[0];
						CHAM_X[2][cur_round] = input_x[2];
						CHAM_X[3][cur_round] = input_x[3];

					//	Output_X[1] = input_x[2];
					//	Output_X[2] = input_x[3];
					}

					if(cur_round >= 3 ) //第4轮之后跳转的CHAM_Diff_trail_round_r进行两个输入差分的查询输出差分的过程
					{
						best = CHAM_RX_trail_round_r(search_round, cur_round+1, Output_X);
					}
					else
					{
						best = CHAM_RX_trail_round_23(search_round, cur_round+1, Output_X);
					}
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best == 1 )  //判断本轮是否要更新RX_Expected_w
	{
		if(updata_cur_num < Num_Bn_Update)   //快速退出机制
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知上一轮最新的RX_Expected_w
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%//快速退出机制//%%%%%%%%
#if (CHAM_even_round == 0 ) //当前搜索参数（1，8）开始的奇数轮的最优路径
	if((cur_round & 0x01) == 0)  //当前偶数轮%2
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
#else  //CHAM_even_round=1, 当前搜索参数（8,1）开始的偶数轮的最优路径
	if((cur_round & 0x01) == 0)  //偶数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
#endif
		if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx ){
			breakfor = 1;
			break;
		}  //退出j1的循环
	}}
#endif
	}
		}}}}
///CHAM_64_or_128 == 128
#else
		Input_RX0 = input_x[0]^ CHAM128_const_rxd[cur_round];
		for(i=0; i< nBytes; i++)
		{
			X0_block[i] = (Input_RX0 >> (8*i)) & 0xFF;   //取出每个8比特
		}
		//确定每个输入差分alpha的每8比特的子块所对应的可能的最小差分概率重量
		X0_block_wt_min[0] = CHAM_a_wt_min[X0_block[0]][0];
		X0_block_wt_min[1] = 0xFF;
		X0_block_wt_min[2] = 0xFF;
		X0_block_wt_min[3] = 0xFF;
		for(j=0;j<3;j++)
		{
			if( (X0_block[j] & 0x80) == 0)
			{
				for(i=0;i<4;i++) //carry 为 0/1/2/3
				{
					if(X0_block_wt_min[j+1] > CHAM_a_wt_min[X0_block[j+1]][i])
					{
						X0_block_wt_min[j+1] = CHAM_a_wt_min[X0_block[j+1]][i];
					}
				}
			}
			else
			{
				for(i=4;i<8;i++)  //carry 为 4/5/6/7
				{
					if(X0_block_wt_min[j+1] > CHAM_a_wt_min[X0_block[j+1]][i])
					{
						X0_block_wt_min[j+1] = CHAM_a_wt_min[X0_block[j+1]][i];
					}
				}
			}
		}

		if(X0_block_wt_min[0] + X0_block_wt_min[1]+
				X0_block_wt_min[2] + X0_block_wt_min[3] > w_cmp)
		{return 0;}

		X0_block_msb[0] = ((X0_block[0] >> 7) << 2);  //alpha 最低8比特的 MSB
		X0_block_msb[1] = ((X0_block[1] >> 7) << 2);  //alpha 8-16比特的 MSB
		X0_block_msb[2] = ((X0_block[2] >> 7) << 2);  //alpha 16-24比特的 MSB

		//整理给下一轮的输入差分
		Output_X[1] = input_x[2];
		Output_X[2] = input_x[3];

		//block 0
		for(i0=CHAM_a_wt_min[X0_block[0]][0];
				i0<= CHAM_a_wt_min[X0_block[0]][0] +1; i0++)   //CHAM_a_wt_max[X0_block[0]][0]
		{
			if(i0 + X0_block_wt_min[1] +
					X0_block_wt_min[2] + X0_block_wt_min[3] > w_cmp ){break;}
			for(j0=0; j0<CHAM_a_bg_numb[X0_block[0]][0][i0]; j0++)
			{
				beta_block[0] = CHAM_a_beta[X0_block[0]][0][i0][j0];
				gama_block[0] = CHAM_a_gama[X0_block[0]][0][i0][j0];
				carry_block[0] = X0_block_msb[0] +
						((beta_block[0] >> 7) << 1) + (gama_block[0] >> 7);

		//block 1
		for(i1=CHAM_a_wt_min[X0_block[1]][carry_block[0]];
				i1<=CHAM_a_wt_min[X0_block[1]][carry_block[0]] +1; i1++) // CHAM_a_wt_max[X0_block[1]][carry_block[0]]
			{
				if(i0 + i1 +
						X0_block_wt_min[2] + X0_block_wt_min[3] > w_cmp ){break;}
				for(j1=0; j1<CHAM_a_bg_numb[X0_block[1]][carry_block[0]][i1]; j1++)
				{
					beta_block[1] = CHAM_a_beta[X0_block[1]][carry_block[0]][i1][j1];
					gama_block[1] = CHAM_a_gama[X0_block[1]][carry_block[0]][i1][j1];
					carry_block[1] = X0_block_msb[1] +
							((beta_block[1] >> 7) << 1) + (gama_block[1] >> 7);

		//block 2
		for(i2=CHAM_a_wt_min[X0_block[2]][carry_block[1]];
				i2<=CHAM_a_wt_min[X0_block[2]][carry_block[1]] +1; i2++)  //CHAM_a_wt_max[X0_block[2]][carry_block[1]]
		{
			if(i0 + i1 + i2 +
					X0_block_wt_min[3] > w_cmp ){break;}
			for(j2=0; j2<CHAM_a_bg_numb[X0_block[2]][carry_block[1]][i2]; j2++)
			{
				beta_block[2] = CHAM_a_beta[X0_block[2]][carry_block[1]][i2][j2];
				gama_block[2] = CHAM_a_gama[X0_block[2]][carry_block[1]][i2][j2];
				carry_block[2] = X0_block_msb[2] +
						((beta_block[2] >> 7) << 1) + (gama_block[2] >> 7);

		//block 3
		for(i3=msb_CHAM_a_wt_min[X0_block[3]][carry_block[2]];
				i3<=msb_CHAM_a_wt_min[X0_block[3]][carry_block[2]] +1; i3++)  //<8 //msb_CHAM_a_wt_max[X0_block[3]][carry_block[2]]
		{
			P_w[cur_round] = i0 + i1 + i2 + i3;  //// 当前轮的概率重量
			if(P_w[cur_round] > w_cmp ){break;}
			//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
			for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			{
				if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
				CHAM_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
				CHAM_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
				RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
				RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

			for(j3=0; j3< msb_CHAM_a_bg_numb[X0_block[3]][carry_block[2]][i3]; j3++)
			{
				beta_block[3] = msb_CHAM_a_beta[X0_block[3]][carry_block[2]][i3][j3];
				gama_block[3] = msb_CHAM_a_gama[X0_block[3]][carry_block[2]][i3][j3];

				beta_tmp = ((beta_block[3] << 24) | (beta_block[2] << 16) |
						(beta_block[1] << 8) | beta_block[0]) & Bit_Align;
				gamma_tmp = ((gama_block[3] << 24) | (gama_block[2] << 16) |
						(gama_block[1] << 8) | gama_block[0]) & Bit_Align;

				//整理给下一轮的输入差分
				Output_X[0] = ROTATE_RIGHT(beta_tmp, ROL_a, blocksize_len) & Bit_Align;
				Out_RXD_tmp = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]];
				Output_X[3] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;

				if(cur_round == 2 ) //记录第2轮的输入
				{
					CHAM_X[2][1] = Output_X[0];

					CHAM_X[0][cur_round] = input_x[0];
					CHAM_X[1][cur_round] = Output_X[0];
					//CHAM_X[2][cur_round] = CHAM_X[1][3];
					CHAM_X[3][cur_round] = input_x[3];

				//	Output_X[2] = input_x[3];
				}
				if(cur_round == 3 ) //记录第3轮的输出
				{
					CHAM_X[3][1] = Output_X[0];

					CHAM_X[2][2] = Output_X[0];

					CHAM_X[0][cur_round] = input_x[0];
					CHAM_X[1][cur_round] = Output_X[0];
					CHAM_X[2][cur_round] = input_x[2];
					CHAM_X[3][cur_round] = input_x[3];

				//	Output_X[1] = input_x[2];
				//	Output_X[2] = input_x[3];
				}

				if(cur_round >= 3 ) //第4轮之后跳转的CHAM_Diff_trail_round_r进行两个输入差分的查询输出差分的过程
				{
					best = CHAM_RX_trail_round_r(search_round, cur_round+1, Output_X);
				}
				else
				{
					best = CHAM_RX_trail_round_23(search_round, cur_round+1, Output_X);
				}
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best == 1 )  //判断本轮是否要更新RX_Expected_w
	{
		if(updata_cur_num < Num_Bn_Update)   //快速退出机制
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知上一轮最新的RX_Expected_w
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%//快速退出机制//%%%%%%%%
#if (CHAM_even_round == 0 ) //当前搜索参数（1，8）开始的奇数轮的最优路径
	if((cur_round & 0x01) == 0)  //当前偶数轮%2
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
#else  //CHAM_even_round=1, 当前搜索参数（8,1）开始的偶数轮的最优路径
	if((cur_round & 0x01) == 0)  //偶数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
#endif
		if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx ){
			breakfor = 1;
			break;
		}  //退出j1的循环
	}}
#endif
		}
			}
		}}}}}}}
#endif
		return updated;
}


u32 CHAM_RX_trail_round_r(u16 search_round, u16 cur_round, u32 *input_x)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 ROL_a = 1, ROL_b = 8;
	u32 gamma_tmp=0;
	u32 Output_X[4] = {0};
	u32 Input_RX0 = 0;
	u32 Out_RXD_tmp = 0;
	u32 gama_block[4] = {0};
	u8 carry_block[4] = {0};
	u8 X0_block[4] = {0};
	u32 X1_tmp = 0;
	u8 X1_block[4] = {0};
	u8 X0X1_block_msb[4] = {0};
	u32 AB_block[4] = {0};
	u8 AB_block_wt_min[4] = {0};
	u32 xor_alpha_beta = 0;
//	u16 i=0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	float w_cmp = 0;
	u16 i0=0,i1=0,i2=0,i3=0;
	u16 j0=0,j1=0,j2=0,j3=0;
//	u8 Block_wt_max[4] = {8,8,8,8};

	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制

	updata_cur_num = Num_Bn_Update;

	if( cur_round >= search_round)  //跳转到最后一轮
	{
		best = CHAM_RX_trail_round_N(search_round, input_x);
		return best;
	}
	CHAM_X[0][cur_round] = input_x[0];
	CHAM_X[1][cur_round] = input_x[1];
	CHAM_X[2][cur_round] = input_x[2];
	CHAM_X[3][cur_round] = input_x[3];


#if (CHAM_even_round == 0 ) //当前搜索参数（1，8）开始的奇数轮的最优路径
	if((cur_round & 0x01) == 0)  //当前偶数轮%2
	{
		ROL_a = 8; ROL_b = 1;

		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
	else  // 奇数轮
	{
		//ROL_a =1; ROL_b = 8; //默认的值就是1,8
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
#else  //CHAM_even_round=1, 当前搜索参数（8,1）开始的偶数轮的最优路径
	if((cur_round & 0x01) == 0)  //偶数轮
	{
		//ROL_a = 1; ROL_b = 8;  //默认的值就是1,8
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
	else  // 奇数轮
	{
		ROL_a = 8; ROL_b = 1;

		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
		w_cmp = w_cmp_rx - Ek_min;
	}
#endif


	//处理CHAM-64 或者 CHAM-128
#if (CHAM_64_or_128 == 64)
		X1_tmp = ROTATE_LEFT(input_x[1], ROL_a, blocksize_len); // & Bit_Align;
		//确定每个输入差分alpha||beta的每8比特的子块所对应的可能的最小差分概率重量
		Input_RX0 = input_x[0]^ CHAM64_const_rxd[cur_round];
		xor_alpha_beta = (X1_tmp ^ Input_RX0); // & 0x7FFF;
		if( _mm_popcnt_u64(xor_alpha_beta & 0x7FFF) > w_cmp)
		{return 0;}

/*		for(i=0; i< nBytes; i++)
		{
			X0_block[i] = (input_x[0] >> (8*i)) & 0xFF;   //取出每个8比特
			X1_block[i] = (X1_tmp >> (8*i)) & 0xFF;   //取出每个8比特
			AB_block[i] = X1_block[i]*256 + X0_block[i];
			AB_block_wt_min[i] = (u8)_mm_popcnt_u64( (xor_alpha_beta >> (i*8)) & 0xFF );
		}
*/
		X0_block[0] = Input_RX0 & 0xFF;   //取出每个8比特
		X1_block[0] = X1_tmp & 0xFF;   //取出每个8比特
		AB_block[0] = X1_block[0] *256 + X0_block[0];
		AB_block_wt_min[0] = (u8)_mm_popcnt_u64(xor_alpha_beta & 0xFF);
		X0_block[1] = (Input_RX0 >> 8) & 0xFF;   //取出每个8比特
		X1_block[1] = (X1_tmp >> 8) & 0xFF;   //取出每个8比特
		AB_block[1] = X1_block[1] *256 + X0_block[1];
		AB_block_wt_min[1] = (u8)_mm_popcnt_u64(xor_alpha_beta & 0x7F00);

		X0X1_block_msb[0] = ((X0_block[0] >> 7) << 2) + ((X1_block[0] >> 7) << 1);  //alpha 最低8比特的 MSB
		//整理给下一轮的输入差分
		Output_X[0] = input_x[1];
		Output_X[1] = input_x[2];
		Output_X[2] = input_x[3];

		//block 0
		for(i0=cDDT_wt_min[0][AB_block[0]];
				i0<= cDDT_wt_min[0][AB_block[0]] ; i0++)  //<=8  // //cDDT_wt_max[0][AB_block[0]]
		{
			if(i0 + AB_block_wt_min[1] > w_cmp){break;}
			for(j0=0; j0 < cDDT_n[0][AB_block[0]][i0]; j0++)
			{
				gama_block[0] = cDDT_v[0][AB_block[0]][i0][j0];
				carry_block[0] = X0X1_block_msb[0] + (gama_block[0] >> 7);

		//block 1  MSB_cDDT_v
		for(i1=MSB_cDDT_wt_min[carry_block[0]][AB_block[1]];
				i1<=MSB_cDDT_wt_min[carry_block[0]][AB_block[1]] ; i1++)  //<8 //Block_wt_max[1] //MSB_cDDT_wt_max[carry_block[0]][AB_block[1]]
			{
				P_w[cur_round] = i0 + i1 ;  //// 当前轮的概率重量
				if(P_w[cur_round] > w_cmp ){break;}
				//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
				for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
				{
					if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
					CHAM_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
					CHAM_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
					RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
					RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

					for(j1=0; j1 < MSB_cDDT_n[carry_block[0]][AB_block[1]][i1]; j1++)
					{
						gama_block[1] = MSB_cDDT_v[carry_block[0]][AB_block[1]][i1][j1];
						gamma_tmp = ((gama_block[1] << 8) | gama_block[0]) & Bit_Align;
						//整理给下一轮的输入差分
						Out_RXD_tmp = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]];
						Output_X[3] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;

						best = CHAM_RX_trail_round_r(search_round, cur_round+1, Output_X);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best == 1 )  //判断本轮是否要更新RX_Expected_w
	{
		if(updata_cur_num < Num_Bn_Update)   //快速退出机制
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知上一轮最新的RX_Expected_w
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%//快速退出机制//%%%%%%%%
#if (CHAM_even_round == 0 ) //当前搜索参数（1，8）开始的奇数轮的最优路径
	if((cur_round & 0x01) == 0)  //当前偶数轮%2
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
#else  //CHAM_even_round=1, 当前搜索参数（8,1）开始的偶数轮的最优路径
	if((cur_round & 0x01) == 0)  //偶数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
#endif
		if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx ){
			breakfor = 1;
			break;
		}  //退出j1的循环
	}}
#endif
	}
		}}}}
#else //CHAM-128
		X1_tmp = ROTATE_LEFT(input_x[1], ROL_a, blocksize_len); // & Bit_Align;
		//确定每个输入差分alpha||beta的每8比特的子块所对应的可能的最小差分概率重量
		Input_RX0 = input_x[0]^ CHAM128_const_rxd[cur_round];
		//确定每个输入差分alpha||beta的每8比特的子块所对应的可能的最小差分概率重量
		xor_alpha_beta = (X1_tmp ^ Input_RX0) & Bit_Align;
/*		for(i=0; i< nBytes; i++)
		{
			X0_block[i] = (input_x[0] >> (8*i)) & 0xFF;   //取出每个8比特
			X1_block[i] = (X1_tmp >> (8*i)) & 0xFF;   //取出每个8比特
			AB_block[i] = X1_block[i]*256 + X0_block[i];
			AB_block_wt_min[i] = (u8)_mm_popcnt_u64( (xor_alpha_beta >> (i*8)) & 0xFF );
		}
*/
		X0_block[0] = Input_RX0 & 0xFF;   //取出每个8比特
		X1_block[0] = X1_tmp & 0xFF;   //取出每个8比特
		AB_block[0] = X1_block[0]*256 + X0_block[0];
		AB_block_wt_min[0] = (u8)_mm_popcnt_u64(xor_alpha_beta & 0xFF );
		X0_block[1] = (Input_RX0 >> 8) & 0xFF;   //取出每个8比特
		X1_block[1] = (X1_tmp >> 8) & 0xFF;   //取出每个8比特
		AB_block[1] = X1_block[1]*256 + X0_block[1];
		AB_block_wt_min[1] = (u8)_mm_popcnt_u64( xor_alpha_beta & 0x0000FF00 );
		X0_block[2] = (Input_RX0 >> 16) & 0xFF;   //取出每个8比特
		X1_block[2] = (X1_tmp >> 16) & 0xFF;   //取出每个8比特
		AB_block[2] = X1_block[2]*256 + X0_block[2];
		AB_block_wt_min[2] = (u8)_mm_popcnt_u64( xor_alpha_beta & 0x00FF0000 );
		X0_block[3] = (Input_RX0 >> 24) & 0xFF;   //取出每个8比特
		X1_block[3] = (X1_tmp >> 24) & 0xFF;   //取出每个8比特
		AB_block[3] = X1_block[3]*256 + X0_block[3];
		AB_block_wt_min[3] = (u8)_mm_popcnt_u64( xor_alpha_beta & 0x7F000000 );


		//粗略剪枝条件
		if( AB_block_wt_min[0] + AB_block_wt_min[1] +
				AB_block_wt_min[2] + AB_block_wt_min[3] > w_cmp)
		{return 0;}

		X0X1_block_msb[0] = ((X0_block[0] >> 7) << 2) + ((X1_block[0] >> 7) << 1);  //alpha 最低8比特的 MSB
		X0X1_block_msb[1] = ((X0_block[1] >> 7) << 2) + ((X1_block[1] >> 7) << 1);  //alpha 最低8比特的 MSB
		X0X1_block_msb[2] = ((X0_block[1] >> 7) << 2) + ((X1_block[2] >> 7) << 1);  //alpha 最低8比特的 MSB

		//整理给下一轮的输入差分
		Output_X[0] = input_x[1];
		Output_X[1] = input_x[2];
		Output_X[2] = input_x[3];

		//block 0
		for(i0=cDDT_wt_min[0][AB_block[0]];
				i0<=cDDT_wt_min[0][AB_block[0]] ; i0++)
		{
			if(i0 + AB_block_wt_min[1] +
					AB_block_wt_min[2] + AB_block_wt_min[3] > w_cmp ){break;}
			for(j0=0; j0<cDDT_n[0][AB_block[0]][i0]; j0++)
			{
				gama_block[0] = cDDT_v[0][AB_block[0]][i0][j0];
				carry_block[0] = X0X1_block_msb[0] + (gama_block[0] >> 7);

		//block 1
		for(i1=cDDT_wt_min[carry_block[0]][AB_block[1]];
				i1<=cDDT_wt_min[carry_block[0]][AB_block[1]] ; i1++) ///cDDT_wt_max[carry_block[0]][AB_block[1]]
			{
				if(i0 + i1 +
						AB_block_wt_min[2] + AB_block_wt_min[3] > w_cmp ){break;}
				for(j1=0; j1<cDDT_n[carry_block[0]][AB_block[1]][i1]; j1++)
				{
					gama_block[1] = cDDT_v[carry_block[0]][AB_block[1]][i1][j1];
					carry_block[1] = X0X1_block_msb[1] + (gama_block[1] >> 7);

		//block 2
		for(i2=cDDT_wt_min[carry_block[1]][AB_block[2]];
				i2<=cDDT_wt_min[carry_block[1]][AB_block[2]] ; i2++)  ////
		{
			if(i0 + i1 + i2 +
					AB_block_wt_min[3] > w_cmp ){break;}
			for(j2=0; j2<cDDT_n[carry_block[1]][AB_block[2]][i2]; j2++)
			{
				gama_block[2] = cDDT_v[carry_block[1]][AB_block[2]][i2][j2];
				carry_block[2] = X0X1_block_msb[2] + (gama_block[2] >> 7);

		//block 3
		for(i3=MSB_cDDT_wt_min[carry_block[2]][AB_block[3]];
				i3<=MSB_cDDT_wt_min[carry_block[2]][AB_block[3]] ; i3++)  //
		{
			P_w[cur_round] = i0 + i1 + i2 + i3;  //// 当前轮的概率重量
			if(P_w[cur_round] > w_cmp ){break;}
			//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
			for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			{
				if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
				CHAM_Pr_zeta_tmp[cur_round] = RX_zeta_w[zeta_num];
				CHAM_zeta_tmp[cur_round] = RX_zeta[RX_zeta_i[zeta_num]];
				RX_P_w[cur_round] = P_w[cur_round] + RX_zeta_w[zeta_num];
				RX_P_sumofR_w[cur_round] = RX_P_sumofR_w[cur_round -1] + RX_P_w[cur_round];

				for(j3=0; j3<MSB_cDDT_n[carry_block[2]][AB_block[3]][i3]; j3++)
				{
					gama_block[3] = MSB_cDDT_v[carry_block[2]][AB_block[3]][i3][j3];
					gamma_tmp = ((gama_block[3] << 24) | (gama_block[2] << 16) |
							(gama_block[1] << 8) | gama_block[0]) & Bit_Align;

					//整理给下一轮的输入差分
					Out_RXD_tmp = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]];
					Output_X[3] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;

					best = CHAM_RX_trail_round_r(search_round, cur_round+1, Output_X);
#if 1   //更新期望概率重量RX_Expected_w,best为1表示来自下一路有更新
	if(best == 1 )  //判断本轮是否要更新RX_Expected_w
	{
		if(updata_cur_num < Num_Bn_Update)   //快速退出机制
		{
			updata_cur_num = Num_Bn_Update;
		updated = 1;  //准备通知上一轮最新的RX_Expected_w
		////更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%//快速退出机制//%%%%%%%%
#if (CHAM_even_round == 0 ) //当前搜索参数（1，8）开始的奇数轮的最优路径
	if((cur_round & 0x01) == 0)  //当前偶数轮%2
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
#else  //CHAM_even_round=1, 当前搜索参数（8,1）开始的偶数轮的最优路径
	if((cur_round & 0x01) == 0)  //偶数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- CHAM_RX_Pw_even_wt[search_round - cur_round];
	}
	else  // 奇数轮
	{
		w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[cur_round -1]
					- RX_n_P_bestofR_w[search_round - cur_round];
	}
#endif
		if(P_w[cur_round] + RX_zeta_w[zeta_num] >= w_cmp_rx ){
			breakfor = 1;
			break;
		}  //退出j1的循环
	}}
#endif
		}
		}}}}}}}}
#endif
		return updated;
}



u32 CHAM_RX_trail_round_N(u16 search_round, u32 *input_x)
{
	u32 best = 0, updated = 0;  //updated通知上一轮更新了，best是来自下一轮的更新
	u16 ROL_a = 8, ROL_b = 1;
	u32 gamma_tmp=0;
//	u32 Output_X[4] = {0};
	u32 Input_RX0 = 0;
	u32 Out_RXD_tmp = 0;
	u32 gama_block[4] = {0};
	u8 carry_block[4] = {0};
	u8 X0_block[4] = {0};
	u32 X1_tmp = 0;
	u8 X1_block[4] = {0};
	u8 X0X1_block_msb[4] = {0};
	u32 AB_block[4] = {0};
	u8 AB_block_wt_min[4] = {0};
	u32 xor_alpha_beta = 0;
	u16 i=0;
//	u16 p_sumof_r = 0; // p_sumof_r is the sum of r-1 rounds weight.
	float w_cmp = 0;
	u16 i0=0,i1=0,i2=0,i3=0;
	u16 j0=0,j1=0,j2=0,j3=0;

	u32 zeta_num = 0;
	float w_cmp_rx = 0;
	u32 updata_cur_num = 0;

	u32 breakfor = 0;  //快速退出机制

	updata_cur_num = Num_Bn_Update;



#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
	{
		if((search_round &0x01) == 0)  //偶数轮
		{
			ROL_a = 8; ROL_b = 1;
		}
		else  //奇数轮
		{
			ROL_a =1; ROL_b = 8;
		}
	}
#else // 搜索的路径为从奇数轮开始的
	{
		if((search_round &0x01) == 0)  //偶数轮
		{
			ROL_a = 1; ROL_b = 8;
		}
		else  //奇数轮
		{
			ROL_a =8; ROL_b = 1;
		}
	}
#endif

	w_cmp_rx = RX_Expected_w - RX_P_sumofR_w[search_round -1];
	w_cmp = w_cmp_rx - Ek_min;


//处理CHAM-64 或者 CHAM-128
#if (CHAM_64_or_128 == 64)
	X1_tmp = ROTATE_LEFT(input_x[1], ROL_a, blocksize_len); // & Bit_Align;
	//确定每个输入差分alpha||beta的每8比特的子块所对应的可能的最小差分概率重量
	Input_RX0 = input_x[0]^ CHAM64_const_rxd[search_round];
	xor_alpha_beta = (X1_tmp ^ Input_RX0); // & 0x7FFF;
	if( _mm_popcnt_u64(xor_alpha_beta & 0x7FFF) > w_cmp)
	{return 0;}

	CHAM_X[0][search_round] = input_x[0];
	CHAM_X[1][search_round] = input_x[1];
	CHAM_X[2][search_round] = input_x[2];
	CHAM_X[3][search_round] = input_x[3];

/*		for(i=0; i< nBytes; i++)
		{
			X0_block[i] = (input_x[0] >> (8*i)) & 0xFF;   //取出每个8比特
			X1_block[i] = (X1_tmp >> (8*i)) & 0xFF;   //取出每个8比特
			AB_block[i] = X1_block[i]*256 + X0_block[i];
			AB_block_wt_min[i] = (u8)_mm_popcnt_u64( (xor_alpha_beta >> (i*8)) & 0xFF );
		}
*/
		X0_block[0] = Input_RX0 & 0xFF;   //取出每个8比特
		X1_block[0] = X1_tmp & 0xFF;   //取出每个8比特
		AB_block[0] = X1_block[0] *256 + X0_block[0];
		AB_block_wt_min[0] = (u8)_mm_popcnt_u64(xor_alpha_beta & 0xFF);
		X0_block[1] = (Input_RX0 >> 8) & 0xFF;   //取出每个8比特
		X1_block[1] = (X1_tmp >> 8) & 0xFF;   //取出每个8比特
		AB_block[1] = X1_block[1] *256 + X0_block[1];
		AB_block_wt_min[1] = (u8)_mm_popcnt_u64(xor_alpha_beta & 0x7F00);

		X0X1_block_msb[0] = ((X0_block[0] >> 7) << 2) + ((X1_block[0] >> 7) << 1);  //alpha 最低8比特的 MSB

		//最后一轮，只用考虑到wt取最小值的情况即可
		for(i0=cDDT_wt_min[0][AB_block[0]];
				i0<= cDDT_wt_min[0][AB_block[0]]; i0++) //8 //Block_wt_max[0]
		{
			if(i0 + AB_block_wt_min[1] > w_cmp ){break;}
			for(j0=0; j0 < cDDT_n[0][AB_block[0]][i0]; j0++)
			{
				gama_block[0] = cDDT_v[0][AB_block[0]][i0][j0];
				carry_block[0] = (X0X1_block_msb[0] + (gama_block[0] >> 7));

//		for(i1=MSB_cDDT_wt_min[carry_block[0]][AB_block[1]];
//				i1<=MSB_cDDT_wt_min[carry_block[0]][AB_block[1]]; i1++) //8 //Block_wt_max[1]
			{
				i1=MSB_cDDT_wt_min[carry_block[0]][AB_block[1]];
			P_w[search_round] = i0 + i1 ;  //// 当前轮的概率重量
			if(P_w[search_round] >= w_cmp ){continue;}
			//基于RX_XOR_offset_Table表，遍历可能的zeta，P_zeta单调递减
		//	for(zeta_num=0; zeta_num < RX_zeta_total ; zeta_num++)
			{
				zeta_num=0;
				if(P_w[search_round] + RX_zeta_w[zeta_num] >= w_cmp_rx){break;}
				CHAM_Pr_zeta_tmp[search_round] = RX_zeta_w[zeta_num];
				CHAM_zeta_tmp[search_round] = RX_zeta[RX_zeta_i[zeta_num]];
				RX_P_w[search_round] = P_w[search_round] + RX_zeta_w[zeta_num];
				RX_P_sumofR_w[search_round] = RX_P_sumofR_w[search_round -1] + RX_P_w[search_round];

			//	for(j1=0; j1 < MSB_cDDT_n[carry_block[0]][AB_block[1]][i1]; j1++)
				{
					gama_block[1] = MSB_cDDT_v[carry_block[0]][AB_block[1]][i1][j1];
					gamma_tmp = ((gama_block[1] << 8) | gama_block[0]) & Bit_Align;

					//整理给最后一轮的输输出差分
					CHAM_RXD[0][search_round+1] = input_x[1];
					CHAM_RXD[1][search_round+1] = input_x[2];
					CHAM_RXD[2][search_round+1] = input_x[3];
					Out_RXD_tmp = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]];
					CHAM_RXD[3][search_round+1] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;
//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
					Num_Bn_Update++;  //更新次数记录
					RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
					best = 1;
					//		breakfor =1; //快速退出机制
					updated = 1;
				}
				//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
					for(i=1; i <= search_round; i++ )
					{
						CHAM_RXD[0][i] = CHAM_X[0][i];
						CHAM_RXD[1][i] = CHAM_X[1][i];
						CHAM_RXD[2][i] = CHAM_X[2][i];
						CHAM_RXD[3][i] = CHAM_X[3][i];

						CHAM_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
						CHAM_Pr_add[i] = P_w[i];   //每轮的模加的XOR差分概率重量
						CHAM_Pr_zeta[i] =  CHAM_Pr_zeta_tmp[i];
						CHAM_Zeta[i] =  CHAM_zeta_tmp[i];
					}
					//CHAM_RXD[0][search_round+1] = CHAM_X[0][search_round+1];
					//CHAM_RXD[1][search_round+1] = CHAM_X[1][search_round+1];
					//CHAM_RXD[2][search_round+1] = CHAM_X[2][search_round+1];
					//CHAM_RXD[3][search_round+1] = CHAM_X[3][search_round+1];
					//%%%%%打印更新次数等信息
					//	printf("==>> RX_Bn_w: %f   Num_Bn_Update: %d ", RX_Expected_w, Num_Bn_Update);
					time_Round = clock();
					run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
					printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
					printf("Time: %.2f seconds. ==>>\n", run_time);

					print_CHAM_64_RXD_resoult(search_round);  //打印
		}}}}

#else  //搜索CHAM128
		X1_tmp = ROTATE_LEFT(input_x[1], ROL_a, blocksize_len); // & Bit_Align;
		//确定每个输入差分alpha||beta的每8比特的子块所对应的可能的最小差分概率重量
		Input_RX0 = input_x[0]^ CHAM128_const_rxd[search_round];
		//确定每个输入差分alpha||beta的每8比特的子块所对应的可能的最小差分概率重量
		xor_alpha_beta = (X1_tmp ^ Input_RX0) & Bit_Align;

		for(i=0; i< nBytes; i++)
		{
			X0_block[i] = (Input_RX0 >> (8*i)) & 0xFF;   //取出每个8比特
			X1_block[i] = (X1_tmp >> (8*i)) & 0xFF;   //取出每个8比特
			AB_block[i] = X1_block[i]*256 + X0_block[i];
			AB_block_wt_min[i] = (u8)_mm_popcnt_u64( (xor_alpha_beta >> (i*8)) & 0xFF );
		}

		//粗略剪枝条件
		if( AB_block_wt_min[0] + AB_block_wt_min[1] +
				AB_block_wt_min[2] + AB_block_wt_min[3] > w_cmp)
		{return 0;}

		CHAM_X[0][search_round] = input_x[0];
		CHAM_X[1][search_round] = input_x[1];
		CHAM_X[2][search_round] = input_x[2];
		CHAM_X[3][search_round] = input_x[3];

		X0X1_block_msb[0] = ((X0_block[0] >> 7) << 2) + ((X1_block[0] >> 7) << 1);  //alpha 最低8比特的 MSB
		X0X1_block_msb[1] = ((X0_block[1] >> 7) << 2) + ((X1_block[1] >> 7) << 1);  //alpha 最低8比特的 MSB
		X0X1_block_msb[2] = ((X0_block[1] >> 7) << 2) + ((X1_block[2] >> 7) << 1);  //alpha 最低8比特的 MSB


		//block 0
		for(i0=cDDT_wt_min[0][AB_block[0]];
				i0<=cDDT_wt_min[0][AB_block[0]]; i0++)
		{
			if(i0 + AB_block_wt_min[1] +
					AB_block_wt_min[2] + AB_block_wt_min[3] > w_cmp ){break;}
			for(j0=0; j0<cDDT_n[0][AB_block[0]][i0]; j0++)
			{
				gama_block[0] = cDDT_v[carry_block[0]][AB_block[0]][i0][j0];
				carry_block[0] = X0X1_block_msb[0] + (gama_block[0] >> 7);

		//block 1
		for(i1=cDDT_wt_min[carry_block[0]][AB_block[1]];
				i1<=cDDT_wt_min[carry_block[0]][AB_block[1]]; i1++)
			{
				if(i0 + i1 +
						AB_block_wt_min[2] + AB_block_wt_min[3]  > w_cmp ){break;}
				for(j1=0; j1<cDDT_n[carry_block[0]][AB_block[1]][i1]; j1++)
				{
					gama_block[1] = cDDT_v[carry_block[0]][AB_block[1]][i1][j1];
					carry_block[1] = X0X1_block_msb[1] + (gama_block[1] >> 7);

		//block 2
		for(i2=cDDT_wt_min[carry_block[1]][AB_block[2]];
				i2<=cDDT_wt_min[carry_block[1]][AB_block[2]]; i2++)
		{
			if(i0 + i1 + i2 +
					AB_block_wt_min[3] > w_cmp ){break;}
			for(j2=0; j2<cDDT_n[carry_block[1]][AB_block[2]][i2]; j2++)
			{
				gama_block[2] = cDDT_v[carry_block[1]][AB_block[2]][i2][j2];
				carry_block[2] = X0X1_block_msb[2] + (gama_block[2] >> 7);

		//block 3
//		for(i3=MSB_cDDT_wt_min[carry_block[2]][AB_block[3]];
//				i3<=MSB_cDDT_wt_min[carry_block[2]][AB_block[3]]; i3++)
		{
			i3=MSB_cDDT_wt_min[carry_block[2]][AB_block[3]];
		P_w[search_round] = i0 + i1 + i2 + i3;  //// 当前轮的概率重量
		if(P_w[search_round] >= w_cmp ){continue;}
//		for(j3=0; j3<MSB_cDDT_n[carry_block[2]][AB_block[3]][i3]; j3++)
		gama_block[3] = MSB_cDDT_v[carry_block[2]][AB_block[3]][i3][j3];
		gamma_tmp = ((gama_block[3] << 24) | (gama_block[2] << 16) |
				(gama_block[1] << 8) | gama_block[0]) & Bit_Align;

		//整理给最后一轮的输输出差分
		CHAM_RXD[0][search_round+1] = input_x[1];
		CHAM_RXD[1][search_round+1] = input_x[2];
		CHAM_RXD[2][search_round+1] = input_x[3];
		Out_RXD_tmp = gamma_tmp ^ RX_zeta[RX_zeta_i[zeta_num]];
		CHAM_RXD[3][search_round+1] = ROTATE_LEFT(Out_RXD_tmp, ROL_b, blocksize_len) & Bit_Align;
//%%%%%%%%%更新期望概率重量RX_Expected_w，逐渐逼近%%%%%%%%%%%%%%%%%%%%%%%
		Num_Bn_Update++;  //更新次数记录
		RX_Expected_w = RX_P_w[search_round] + RX_P_sumofR_w[search_round -1];
		best = 1;
		//	breakfor =1; //快速退出机制
		updated = 1;

		//%%%%%记录每次更新后的路径的详情，RX差分和每轮的概率重量等
		for(i=1; i <= search_round; i++ )
		{
			CHAM_RXD[0][i] = CHAM_X[0][i];
			CHAM_RXD[1][i] = CHAM_X[1][i];
			CHAM_RXD[2][i] = CHAM_X[2][i];
			CHAM_RXD[3][i] = CHAM_X[3][i];

			CHAM_Pr_RX[i] = RX_P_w[i];  //每轮的总的RX概率重量
			CHAM_Pr_add[i] = P_w[i];   //每轮的模加的XOR差分概率重量
			CHAM_Pr_zeta[i] =  CHAM_Pr_zeta_tmp[i];
			CHAM_Zeta[i] =  CHAM_zeta_tmp[i];
		}
		//CHAM_RXD[0][search_round+1] = CHAM_X[0][search_round+1];
		//CHAM_RXD[1][search_round+1] = CHAM_X[1][search_round+1];
		//CHAM_RXD[2][search_round+1] = CHAM_X[2][search_round+1];
		//CHAM_RXD[3][search_round+1] = CHAM_X[3][search_round+1];
		//%%%%%打印更新次数等信息
		//	printf("==>> RX_Bn_w: %f   Num_Bn_Update: %d ", RX_Expected_w, Num_Bn_Update);
		time_Round = clock();
		run_time =  (time_Round - time_ARX_DDT) / CLOCKS_PER_SEC;
		printf("Update: %d   RX_Bn_w: %f  --", Num_Bn_Update, RX_Expected_w);
		printf("Time: %.2f seconds. ==>>\n", run_time);

		print_CHAM_128_RXD_resoult(search_round);  //打印
			}}}}}}}
#endif
		return  updated;
}




u16 print_CHAM_64_RXD_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("round----X0------X1-------X2-------X3-----zeta-----Pr_zeta----Pr_XOR-----Pr-RX \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%04llx  0x%04llx   0x%04llx   0x%04llx   0x%04llx   %f  %f  %f\n",
				r_num,
				CHAM_RXD[0][r_num],CHAM_RXD[1][r_num],CHAM_RXD[2][r_num],CHAM_RXD[3][r_num],
				CHAM_Zeta[r_num],CHAM_Pr_zeta[r_num],
				CHAM_Pr_add[r_num],CHAM_Pr_RX[r_num]);
	}
	printf("%02d     0x%04llx  0x%04llx   0x%04llx   0x%04llx    -   \n",
			r_num,
			CHAM_RXD[0][search_round+1],CHAM_RXD[1][search_round+1],
			CHAM_RXD[2][search_round+1],CHAM_RXD[3][search_round+1]);

	printf("------------------ \n");
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
	{
		printf("%d Round CHAM-64 Total Weight: -%f \n",search_round, RX_Expected_w);
	}
#else
	{
		printf("%d Round CHAM-64 Total Weight: -%f \n",search_round,RX_Expected_w);
	}
#endif
	return 1;
}





u32 print_CHAM_128_RXD_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("round----X0------X1------X2------X3------zeta-----Pr_zeta----Pr_XOR-----Pr-RX\n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d  0x%08llx  0x%08llx  0x%08llx  0x%08llx   0x%08llx   %f   %f   %f \n",
				r_num,
				CHAM_RXD[0][r_num],CHAM_RXD[1][r_num],CHAM_RXD[2][r_num],CHAM_RXD[3][r_num],
				CHAM_Zeta[r_num],CHAM_Pr_zeta[r_num],
				CHAM_Pr_add[r_num],CHAM_Pr_RX[r_num]);
	}
	printf("%02d     0x%08llx  0x%08llx  0x%08llx  0x%08llx       \n",
			r_num,
			CHAM_RXD[0][search_round+1],CHAM_RXD[1][search_round+1],
			CHAM_RXD[2][search_round+1],CHAM_RXD[3][search_round+1]);

	printf("------------------ \n");
#if(CHAM_even_round == 0) //结束条件判断是搜索的奇数还是偶数轮开始的路径
	{
		printf("%d Round CHAM-128 Total Weight: -%f \n",search_round, RX_Expected_w);
	}
#else
	{
		printf("%d Round CHAM-128 Total Weight: -%f \n",search_round, RX_Expected_w);
	}
#endif
	return 1;
}





