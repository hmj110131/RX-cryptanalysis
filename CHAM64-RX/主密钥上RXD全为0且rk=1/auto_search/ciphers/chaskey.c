/*
 * chaskey.c
 *
 *  Created on: 2018年11月29日
 *      Author: hmj110131
 */
#include <stdio.h>
#include "typedef.h"
#include "ciphers.h"
#include "chaskey.h"




u32 rol_left32(u32 x, u32 indx)
{
	return ((x << indx) | (x >> (32-indx))); // & 0xFFFFFFFF;
}
u32 chaskey_rol_left_5(u32 input)
{
	return rol_left32(input,5); // & 0xFFFFFFFF;
}
u32 chaskey_rol_left_8(u32 input)
{
	return rol_left32(input,8); // & 0xFFFFFFFF;
}
u32 chaskey_rol_left_16(u32 input)
{
	return rol_left32(input,16); // & 0xFFFFFFFF;
}
u32 chaskey_rol_left_7(u32 input)
{
	return rol_left32(input,7); // & 0xFFFFFFFF;
}
u32 chaskey_rol_left_13(u32 input)
{
	return rol_left32(input,13); // & 0xFFFFFFFF;
}



u32 chaskey_round_1_left(u16 search_round)
{
	u32 best = 0;


	return best;
}

u32 chaskey_round_1_right(u16 search_round)
{
	u32 best = 0;


	return best;
}

u32 chaskey_round_1_Down_2add(u16 search_round)
{
	u32 best = 0;


	return best;

}




////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////
//u16 Chaskey_R1_ADD_wt[4] = {0};
volatile u64 Chaskey_R1_ADD_i_UVW[4][3] = {0}; //4个模块加法的对应UVW
volatile u16 Chaskey_Round_ADD_wt[5][4] = {0};
volatile u32 XV[4][20] = {0};
volatile u32 XV_M[4][20] = {0};


////////////////////////////////搜索Chaskey的最优线性特征的代码//////////////////////////////





u16 print_Chaskey_resoult(u16 search_round)
{
	u16 r_num = 0;
	u16 i_L = 0;

	printf("round----XV0-----------XV1------------XV2------------XV3---------w_0--w_1--w_2--w_3---Pw \n");  //----gama_cnt \n");
	for(r_num=1;r_num <= search_round;r_num++)
	{
		printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx    -%d   -%d   -%d   -%d    -%d \n",
				r_num,XV[0][r_num],XV[1][r_num],XV[2][r_num],XV[3][r_num],
				Chaskey_Round_ADD_wt[r_num][0],Chaskey_Round_ADD_wt[r_num][1],
				Chaskey_Round_ADD_wt[r_num][2],Chaskey_Round_ADD_wt[r_num][3],
				P_w[r_num]); //,gama_cnt[r_num]);
	}

	printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx     -    -    -    -   NULL \n",
			search_round+1,
			XV[0][search_round+1],XV[1][search_round+1],
			XV[2][search_round+1],XV[3][search_round+1]); //,gama_cnt[r_num]);
	printf("Tmp-----MV0-----------MV1------------MV2------------MV3--------- \n");
	for(i_L = 1; i_L <= search_round;i_L++ )
	{
		printf("%02d     0x%08llx    0x%08llx     0x%08llx     0x%08llx   \n",
				i_L,
				XV_M[0][i_L],XV_M[1][i_L],XV_M[2][i_L],XV_M[3][i_L]); //,gama_cnt[r_num]);
	}


	printf("------------------ \n");
	printf("%d Round Chaskey Total Weight: -%d \n",search_round, n_P_bestofR_w[search_round]);

	return 1;
}







u32 Chaskey_1_Rounds(u32 *X, u32 *Y)
{
	u32 cnt=0;
	u32 v0 =0, v1=0, v2=0, v3=0;
	u32 v0_t =0, v1_t=0, v2_t=0, v3_t=0;
	u32 u0 =0, u1=0, u2=0, u3=0;
	u32 u0_t =0, u1_t=0, u2_t=0, u3_t=0;
	u32 tmp_left = 0, tmp_right = 0;

	u32 out_X[4] = {0};
	u32 out_Y[4] = {0};


///////////X
	v0 = X[0];
	v1 = X[1];
	v2 = X[2];
	v3 = X[0];

	tmp_left = (v0 + v1); // & 0xFFFFFFFF;
	tmp_right = (v2 + v3); // & 0xFFFFFFFF;

	v1_t = chaskey_rol_left_5(v1) ^ tmp_left;
	v0_t = tmp_right;
	v2_t = chaskey_rol_left_16(tmp_left);
	v3_t = v0_t ^ chaskey_rol_left_8(v3);

	tmp_left = (v0_t + v1_t); // & 0xFFFFFFFF;
	tmp_right = (v2_t + v3_t); // & 0xFFFFFFFF;

	out_X[1] = chaskey_rol_left_7(v1_t) ^ tmp_left;
	out_X[0] = tmp_right;
	out_X[2] = chaskey_rol_left_16(tmp_left);
	out_X[3] = out_X[0] ^ chaskey_rol_left_13(v3_t);

	X[0]= out_X[0];
	X[1]= out_X[1];
	X[2]= out_X[2];
	X[3]= out_X[3];
///////////Y
	u0 = Y[0];
	u1 = Y[1];
	u2 = Y[2];
	u3 = Y[0];

	tmp_left = (u0 + u1); // & 0xFFFFFFFF;
	tmp_right = (u2 + u3); // & 0xFFFFFFFF;

	u1_t = chaskey_rol_left_5(u1) ^ tmp_left;
	u0_t = tmp_right;
	u2_t = chaskey_rol_left_16(tmp_left);
	u3_t = u0_t ^ chaskey_rol_left_8(u3);


	tmp_left = (u0_t + u1_t); // & 0xFFFFFFFF;
	tmp_right = (u2_t + u3_t); // & 0xFFFFFFFF;

	out_Y[1] = chaskey_rol_left_7(u1_t) ^ tmp_left;
	out_Y[0] = tmp_right;
	out_Y[2] = chaskey_rol_left_16(tmp_left);
	out_Y[3] = out_Y[0] ^ chaskey_rol_left_13(u3_t);

	Y[0] = out_Y[0];
	Y[1] = out_Y[1];
	Y[2] = out_Y[2];
	Y[3] = out_Y[3];
///////////////////////////比较out_和out_Y///////////////////

	return cnt;
}














