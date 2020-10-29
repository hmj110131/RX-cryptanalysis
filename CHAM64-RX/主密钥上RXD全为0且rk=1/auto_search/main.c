/***
 *
/										    \\\\\|////
/											\\	- -	//
/											 (	@ @	)
+------------------------------------------oOOo-(_)-oOOo-----------------------+
|						@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
|								 Name: auto_search
|								 Function:  Search for the DC&LC of Simon
|								 Version: V1.0
|								 Created On: 2017-7-28
|								 Author:	Mingjiang Huang
|								 Organization:	Chinese Academy of Sciences
|						 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
|														 													|
+-------------------------------------------oooO-----Oooo----------------------+
/											(	)	 (   )
/											 \  )	 (  /
/											  \_)	 (_/
 *Revised:
 * * main.c
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 *
 *
 *
 ***/
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <table.h>
#include <time.h>
#include "math.h"
#include "ciphers.h"
#include "search.h"
#include "typedef.h"
#include "globalvar.h"
#include "printinfo.h"
#include "ciphers.h"
#include "hight.h"
#include "RX_crypt.h"



//这里设置本程序的搜索为XOR差分或线性特征，还是RX差分特征
#define DC_LC 0 //本程序为搜索XOR差分特征或者线性特征
#define RX_DC 1 //本程序为搜索RX差分特征

#define USEMPI 0     //是否开启并行MPI计算,1-yes.
#define USEOPENMP 0  //是否开启并行OpenMP计算,1-yes.
#define Thread_N 4   //定义使用的进程最大数量

#if (USEMPI==1)
#include </opt/mpich/include/mpi.h>
#endif

#if (USEOPENMP==1)
#include <omp.h>
#endif


/** Main Function which is the main entry of the auto-search program.
 ** Calls function pre_processing_info() in the main after initialization,
 ** Iterating on @sc_rounds R (the number of round).
 ** Find the best differential of linear trail of the specific cipher.
 ** */
int main(int argc, char *argv[])
{


#if (USEMPI==1)   /* Parallel Computation. */
	int comm_rank; //进程号
	int comm_size;  //进程数目

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&comm_rank);  //进程号
	MPI_Comm_size(MPI_COMM_WORLD,&comm_size);  //进程数

	//MPI_Barrier(MPI_COMM_WORLD);  //进程间同步//进程锁
	//char processor_name[MPI_MAX_PROCESSOR_NAME];
	//MPI_Get_processor_name(processor_name,&comm_size); //处理器名称

	if (comm_size != Thread_N)
	{
		printf( "Number of processor is not same as expected. comm_rank:%d comm_size: %d \n",comm_rank, comm_size);
		MPI_Finalize();  //并行计算结束

		return 0;
	}
#endif


#if (USEMPI==1)
	int thread_id;
	thread_id = comm_rank;  //进程号
//#else
//	int thread_id;
//	printf("Please input the index of the thread(0~%d).\n", Thread_N-1);
//	scanf("%d", &thread_id);
#endif

#if (USEOPENMP==1)
	omp_set_num_threads(4);
#endif

/*
	int i=0;
	wt_nmu[17] = 1;
	wt_nmu[18] = 0;
	wt_nmu[19] = 3;
	wt_nmu[20] = 1;
	wt_nmu[21] = 9;
	wt_nmu[22] = 14;
	wt_nmu[23] = 68;
	wt_nmu[24] = 171;

	for(i=0;i<=24;i++)  // 超过wt_max的也统计
	{
		prob = prob + (wt_nmu[i] * pow(2,-2*i)); //注意相关性和线性平方相关的区别
		prob_w = log(prob) / log(2); //prob_w = log(prob) / log(2);
		printf("wt:-%d   wt_num:%d   prob_sum:%Lf \n",i,wt_nmu[i],prob_w);
	}
*/



/*
	Bit_Align = 0xFFFFFFFFFFFF;
	//ValueMax = 0xFFFFFFFFFFFF;
	ValueMax_Align = 0x7FFFFFFFFFFF;
	printf("PPPPPP: %d   \n", Speck_XDP_compute(0x080800FE0808,0x0800EE4A0848,0x000772400040));
	printf("PPXXX: %d   \n", HM_weight(0x080800FE0808 ^ 0x0800EE4A0848));
	printf("PPXXX: 0x%016x   \n", 0x080800FE0808 ^ 0x0800EE4A0848);
*/
	//Bit_Align = 0xFFFFFFFFFFFF;
	//ValueMax_Align = 0x7FFFFFFFFFFF;
	//printf("Pwt: %d   \n", Speck_XDP_validcheck(0x400000924000,0x400000104200,0x820200));
	//printf("PPwwwt: %x \n", Speck_XDP_compute(0x008000000000,0x800000000000,0x808000000000) );
	//printf("PPXXX: %x \n", Speck_XDP_validcheck(0x008000000000,0x800000000000,0x808000000000) );
	//Input_alpha = 0x400000924000;
	//Input_gamma = 0x400000924000;
	//Input_beta  = 0x400000104200;


#if (DC_LC == 1)  //注意在本文件开头的宏命令中定义#define DC_LC 1才能开启
///////////////搜索ARX、AND-RX密码的差分和线性特征/////////////////////
	if( pre_processing_info() == 0)
		return 0;
//#else
#elif (RX_DC == 1)   //注意在本文件开头的宏命令中定义#define DC_LC 1才能开启
///////////////搜索ARX密码的RX差分特征/////////////////////
	if( pre_processing_RX_info() == 0)
		return 0;

#endif




#if (USEMPI==1)
	MPI_Finalize();
#endif

	return 1; /* All functions are normally returned 1, and return 0 with unexpected error. */
}
