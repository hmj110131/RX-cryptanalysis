/*
 * RX_crypt.c
 *
 *  Created on: 2020年7月23日
 *      Author: hmj110131
 */

#include <stdio.h>
#include <stdlib.h>
#include <table.h>
#include "math.h"
#include "typedef.h"
#include "ciphers.h"
#include "search.h"
#include "globalvar.h"
#include "printinfo.h"
#include <time.h>
#include <nmmintrin.h>
#include "RX_crypt.h"
#include "XOR_offset_table.h"
#include "RX_Alzette.h"
#include "RX_Alzette.h"
#include "RX_Chaskey.h"
#include "RX_Siphash.h"
#include "RX_CHAM.h"

#if 1
#include "Alzette.h"
#endif

#if 1
#include "speck.h"
#endif

#if 1
#include "sparx.h"
#endif
#if 1
#include "Siphash.h"
#endif


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
///////////以下为RX-dif使用的全局概率重量参数//////////////////
volatile float RX_Expected_w = 0;  //最优路径的概率重量，该值会动态更新
volatile float RX_P_w[100]={0};   //记录RX路径的每轮的概率值
volatile float RX_P_bestofR_w[100] = {0}; //用于记录每轮中的概率重量
volatile float RX_P_sumofR_w[100]={0};;
//RX分析的旋转参数，需要提前输入
u16 RX_k = 1;
//基于概率重量由大到小逼近，选择r-1轮的输出差分为限定条件
u64 F_Alpha=0, F_Beta=0;
//基于概率重量由大到小逼近，统计动态变化更新期望概率重量的次数
u64 Num_Bn_Update = 0;






//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
int pre_processing_RX_info(void)
{
	//ADD_RX_check(1); //return 0;
	printf("|*******************************RX-Cryptanalysis****************************|\n");
	printf("Auto-search Begin... \n");
	printf( "Enter the parameters : #ciper #blocksize #rounds  #mode \n");

	////////用nohup后台执行时，需要先写死参数，编译后直接执行/////////
		strcpy(str_cipher, "cham");    //设定待搜索的密码算法
		strcpy(str_blocksize, "64");   //设定分组长度  //注意需要到对应.h文件中改宏定义的分组长度选择编译相应代码
		sc_rounds = 16;   //设定要搜索的轮数

		/*|***************RX差分特征的模式*************|*/
		RX_k = 1;  //旋转参数
		str_mode = 0;  //搜索模式
		//基于概率重量由大到小逼近，选择r-1轮的输出差分为限定条件,概率重量动态逼近的搜索方法
//		F_Alpha = 0xd823f3a6;  //1轮的最优输出RX差分
//		F_Beta  = 0x01;

//		F_Alpha = 0xd823f3a7; //2轮的最优输出RX差分
//		F_Beta  = 0x67cf8038;

///////需要在下面代码中开启scanf函数，从终端获取输入参数//////////
		/* Get the input ciper name.*/
		printf("ciper: ");
	//	scanf("%s",str_cipher); //从终端输入
		printf("%s \n",str_cipher); ///在程序里面写死


	if( strcmp(str_cipher,"speck") == 0)
	{
		cipher_name = speck;
		/*the Rotation parameter of simon-LIKE cipher, abc */
		rol_a = 8;	/* a */
		rol_b = 3;   /* b */
		rol_delt = 5;   /* a-b */
	}
	else
	if( strcmp(str_cipher,"sparx") == 0)
	{
		cipher_name = sparx;
		rol_a = 7;	/* a */
		rol_b = 2;   /* b */
		//printf("Note: blocksize Now only fixed for 64 or 128 bits.  \n");
	}
	else
	if( strcmp(str_cipher,"alzette") == 0)
	{
		cipher_name = alzette;
		rol_a = 31;	/* a */ //first round rol_constants
		rol_b = 24;   /* b */
		//printf("Note: blocksize only fixed with 64-bits for Alzette.  \n");
	}
	else
	if( strcmp(str_cipher,"chaskey") == 0)
	{
		cipher_name = chaskey;
		//rol_a = 31;	/* a */ //first round rol_constants
		//rol_b = 24;   /* b */
		//printf("Note: blocksize only fixed with 64-bits for Alzette.  \n");
	}
	else
	if( strcmp(str_cipher,"siphash") == 0)
	{
		cipher_name = siphash;
		//rol_a = 31;	/* a */ //first round rol_constants
		//rol_b = 24;   /* b */
		//printf("Note: blocksize only fixed with 64-bits for Alzette.  \n");
	}
	else
	if( strcmp(str_cipher,"cham") == 0)
	{
		cipher_name = cham;
		//rol_a = 31;	/* a */ //first round rol_constants
		//rol_b = 24;   /* b */
		//printf("Note: blocksize only fixed with 64-bits for Alzette.  \n");
	}
	else
	{
		printf("Errors. Input with wrong cipher name! \n");
		return 0;
	}

	/* Get the input ciper variant with corresponding blocksize.*/
	printf("blocksize: ");
//	scanf("%s",str_blocksize); //从终端输入
	printf("%s \n",str_blocksize); ///在程序里面写死

	if( strcmp(str_blocksize,"32") == 0)
	{
		sc_blocksize = 32;
		blocksize_len = 16;
		nBytes = 2;

		if( strcmp(str_cipher,"speck") == 0)
		{
			rol_a = 7;
			rol_b = 2;
			rol_delt = 5;
		}
	}
	else
	if( strcmp(str_blocksize,"48") == 0)
	{
		sc_blocksize = 48;
		blocksize_len = 24;
		nBytes = 3;
		if( strcmp(str_cipher,"speck") == 0)
		{
			rol_a = 8;
			rol_b = 3;
		}
	}
	else
	if( strcmp(str_blocksize,"64") == 0)
	{
		sc_blocksize = 64;
		blocksize_len = 32;
		nBytes = 4;
		if( strcmp(str_cipher,"speck") == 0)
		{
			rol_a = 8;
			rol_b = 3;
		}
	}
	else
	if( strcmp(str_blocksize,"96") == 0)
	{
		sc_blocksize = 96;
		blocksize_len = 48;
		nBytes = 6;
		if( strcmp(str_cipher,"speck") == 0)
		{
			rol_a = 8;
			rol_b = 3;
		}
	}
	else
	if( strcmp(str_blocksize,"128") == 0)
	{
		sc_blocksize = 128;
		blocksize_len = 64;
		nBytes = 8;
		if( strcmp(str_cipher,"speck") == 0)
		{
			rol_a = 8;
			rol_b = 3;
		}

		if( strcmp(str_cipher,"chaskey") == 0)
		{
			sc_blocksize = 128;
			blocksize_len = 32;
			nBytes = 4;
		}
	}
	else
	if( strcmp(str_blocksize,"256") == 0)
	{
		sc_blocksize = 128;
		blocksize_len = 64;
		nBytes = 8;
		if( strcmp(str_cipher,"siphash") == 0)
		{
			sc_blocksize = 256;
			blocksize_len = 64;
			nBytes = 4;
		}
	}
	else
	{
		printf("Error input with blocksize of the cipher! \n");
		return 0;
	}


	/* Get the rotational parameter RX_k. //获取旋转参数 */
	printf("rotational parameter: ");
//	scanf("%s",str_blocksize); //从终端输入
	printf("%d \n",RX_k); ///在程序里面写死



	/* Get the input expected r round max weight of DP.*/
	printf("rounds: ");
//	scanf("%d",&sc_rounds); //从终端输入
	printf("%d \n",sc_rounds); ///在程序里面写死了

	if((sc_rounds < 1) || (sc_rounds > 100))
	{
		printf("Error input search Rounds! \n");
		return 0;
	}


	/* Get the input search mode.*/
	printf("mode: %d \n", str_mode); //在程序里面写死  //0, //5 for SPARX64-diff ....
//	scanf("%d",&str_mode); //从终端输入


	switch(str_mode)
	{
	case 0: // Optimal RX-diff trail. 秘钥为常数,纯RX差分路径
		findMinDPWeight_RX_Trail((u16)sc_rounds);
		break;

	case 1: // Related Key RX-diff trail. 相关秘钥情形下的，RX差分路径
		findMinDPWeight_Relatekey_RX_Trail((u16)sc_rounds);
		break;

	case 2: // weak-key class Related Key RX-diff trail. 弱秘钥情形下，相关秘钥情形下的，RX差分路径
		findMinDPWeight_Weakkey_Relatekey_RX_Trail((u16)sc_rounds);
		break;

	case 3: // Open-key class Related Key RX-diff trail. 开放秘钥情形下，相关秘钥情形下的，RX差分路径
		findMinDPWeight_Openkey_Relatekey_RX_Trail((u16)sc_rounds);
		break;


	default:
		printf("Error search mode,does not have this search mode! \n");
		break;
	}

	return 1;
}



u16 findMinDPWeight_RX_Trail(u16 search_round)
{

	printf("----------------Search Optimal pure RX-Differential Trail-------------------\n");

	if(search_round <= 1)  // 0/1 rounds.
	{
		printf("Error! -----Please set search round > 1.-----\n");
		return 1;
	}


	if( strcmp(str_cipher,"speck") == 0)
	{
		// SPECK differential trail search.
		//sepck_differential_trail_search_entry(search_round);
	}
	else
	if( strcmp(str_cipher,"sparx") == 0)
	{
		//str_cipher[10] = "SPARX64";

		if( sc_blocksize == 64 )
		{
			// SPARX-64 optimal differential trail search.
		//	SPARX64_differential_trail_search_entry(search_round);
		}
		else
		{
			// SPARX-128 optimal differential trail search.
		//	SPARX_128_differential_trail_search_entry(search_round);
		}
	}
	else
	if( strcmp(str_cipher,"cham") == 0)
		{
			if( sc_blocksize == 64 )
			{
			// CHAM-64 optimal differential trail search.
				CHAM_64_RX_trail_search_entry( search_round);
			}
			else
			{
				// CHAM-128 optimal differential trail search.
				CHAM_128_RX_trail_search_entry( search_round);
			}
		}
	else
	if( strcmp(str_cipher,"alzette") == 0)
	{
		if( sc_blocksize == 64 )
		{
		// Alzette optimal differential trail search.
			Alzette_RX_Diff_trail_search_entry( search_round);
		}
		else
		{
			printf("Error! The block size of Alzette(ARX-Box) should be 64 bits. \n");
		}
	}
	else
	if( strcmp(str_cipher,"chaskey") == 0)
	{
		if( sc_blocksize == 128 )
		{
		// Chaskey optimal differential trail search.
			Chaskey_RX_trail_search_entry( search_round);
		}
		else
		{
			printf("Error! The block size of Chaskey(ARX-Box) should be 128 bits. \n");
		}
	}
	else
	if( strcmp(str_cipher,"siphash") == 0)
	{
		if( sc_blocksize == 256 )
		{
		// Chaskey optimal differential trail search.
			Siphash_RX_trail_search_entry( search_round);
		}
		else
		{
			printf("Error! The block size of Chaskey(ARX-Box) should be 128 bits. \n");
		}
	}

	time_finish = clock();
	run_time =  (float)(time_finish - time_start) / CLOCKS_PER_SEC;
	printf( "Time cost: %.2f seconds ==  %.2f minutes == %.2f hours. \n", run_time,run_time/60.0,run_time/3600.0 );
	printf("Auto-RX-search END! \n");
	printf("|************************************************************************|\n");
	return 1;
}


u16 findMinDPWeight_Relatekey_RX_Trail(u16 search_round)
{

	return 1;
}


u16 findMinDPWeight_Weakkey_Relatekey_RX_Trail(u16 search_round)
{



	return 1;
}


u16 findMinDPWeight_Openkey_Relatekey_RX_Trail(u16 search_round)
{


	return 1;
}









