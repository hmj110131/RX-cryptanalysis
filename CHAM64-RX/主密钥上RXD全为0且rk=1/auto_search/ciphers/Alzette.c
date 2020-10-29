/*
 * Alzette.c
 *
 *  Created on: 2020年4月22日
 *      Author: Mingjiang Huang
 */


#include "Alzette.h"
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include <time.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//重复4次Alzette
u16  Alzette_rol_a[16]={31,17,0,24,31,17,0,24,31,17,0,24,31,17,0,24};     //定义每轮的循环参数
u16  Alzette_rol_b[16]={24,17,31,16,24,17,31,16,24,17,31,16,24,17,31,16};     //定义每轮的循环参数
//定义Alzette的实例版本，共8个版本
u32  Alzette_Constant[8]={0xB7E15162, 0xBF715880, 0x38B4DA56, 0x324E7738,
                          0xBB1185EB, 0x4F7C7B57, 0xCFBFA1C8, 0xC2B3293D};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
volatile u16 Alzette_Pw_cor[16] = {0}; // 记录Alzette的最优线性相关性重量。
volatile u16 Alzette_Pw_wt[16] = {0}; // 记录Alzette的最优差分概率重量

volatile u32 Alzette_X[2][16]={0}; //每轮的输入差分左右两部分
volatile u32 Alzette_Alpha[16]={0};  //每轮的第模加的输入输出差分
volatile u32 Alzette_Beta[16]={0};
volatile u32 Alzette_Gamma[16]={0};
volatile u32 Alzette_u[16]={0};  //每轮的第模加的输入输出掩码 //output u
volatile u32 Alzette_v[16]={0};  //input v
volatile u32 Alzette_w[16]={0};  //input w




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//Alzette每轮的用到的循环移位操作,all are rotated to right
u32 rol_right_32(u32 x, u32 indx)
{
	return ((x >> indx) | (x << (32-indx))); // & 0xFFFFFFFF;
}
u32 rol_left_32(u32 x, u32 indx)
{
	return ((x << indx) | (x >> (32-indx))); // & 0xFFFFFFFF;
}
u32 Alzette_rol_right_31(u32 input)
{
	return rol_right_32(input,31); // & 0xFFFFFFFF;
}
u32 Alzette_rol_left_31(u32 input)
{
	return rol_left_32(input,31); // & 0xFFFFFFFF;
}
u32 Alzette_rol_right_24(u32 input)
{
	return rol_right_32(input,24);
}
u32 Alzette_rol_right_17(u32 input)
{
	return rol_right_32(input,17);
}
u32 Alzette_rol_right_0(u32 input)
{
	return input;
}
u32 Alzette_rol_right_16(u32 input)
{
	return rol_right_32(input,16);
}
























void print_Alzette_diff_resoult(u16 search_round)
{
	u16 r_num = 0;

	printf("|-------------------------------------------|\n");
	printf("round-----left---------right------weight \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%08llx    0x%08llx     -%d \n",r_num-1,Alzette_X[0][r_num-1],Alzette_X[1][r_num-1],P_w[r_num]);
	}
	printf("%02d     0x%08llx    0x%08llx     NULL \n",search_round,Alzette_X[0][search_round],Alzette_X[1][search_round]);

	printf("%d Round Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
	printf("|-------------------------------------------|\n");
}

void print_Alzette_linear_resoult(u16 search_round)
{
	u16 r_num = 0;
	printf("|-------------------------------------------|\n");
	printf("round-----left---------right------weight \n");
	for(r_num=0;r_num < search_round;r_num++)
	{
		printf("%02d     0x%08llx    0x%08llx     -%d \n",r_num,Alzette_X[0][r_num],Alzette_X[1][r_num],P_w[r_num+1]);
	}
	printf("%02d     0x%08llx    0x%08llx     NULL \n",search_round,Alzette_X[0][search_round],Alzette_X[1][search_round]);

	printf("%d Round Total Correlation Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);
	printf("|-------------------------------------------|\n");
}










