/***
|				@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     
|				 Name: auto_search                                   
|				 Function:  Searching DC&LC&RX&RK of ARX ciphers      
|				 Version: V1.0                                         
|				 Created since: 2017-7-28                             
|			         Author: Mingjiang Huang                         
|			         Email:	huangmingjiang@iie.ac.cn                    
|				 Organization:	Chinese Academy of Sciences        
|				 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       
|														 	
***/
/*
 * RX_offset_table.c
 *
 *  Created on: 2020年7月7日
 *      Author: hmj110131
 */
#include <stdio.h>
#include "RX_offset_table.h"
#include "typedef.h"
#include "math.h"
#include "ciphers.h"
#include <nmmintrin.h>
#include "RX_crypt.h"



//Ek里面的值，对应k和n都需要在这里提前手工计算得到并提前存储。
//循环左移位参数k=1，且假设n很大的情况
extern float  Ek[4]={1.415, 1.415, 3, 3};     //指定k的，对应CL，CR取值的4种情况的概率重量的小数值，当k和n确定，Ek为固定值,这里以k=1为例子。

//RX分析中用到的修正量和偏移量，并记录对应的概率重量，提前计算好
extern u64     RX_u[64]={0};     //RX_offset的右边部分对应的dummy-offset的值
extern float  RX_u_w[64]={0};     //RX_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u64     RX_v[64]={0};     //RX_offset的左边部分对应的dummy-offset的值
extern float  RX_v_w[64]={0};     //RX_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u64     RX_zeta[3600]={0};  //RX_offset偏移量，zeta
extern float  RX_zeta_w[3600]={0};     //RX_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u16     RX_zeta_i[3600]={0};      //RX_offset的右边部分对应的dummy-offset的值各自对应的概率重量经过排序后的下标
extern u64     RX_zeta_total=0;   //rx-xor表中的zeta的个数
float Ek_min = 0;


//计算指定旋转参数k，所对应的Ek的4个值，通过指针返回
u16 Compute_RX_Ek(u16 rx_k, u16 rx_n,  float *ek)
{
	float pp0 = 0,pp1 = 0,pp2 = 0;
	float pr = 0;

	pp0 = pow(2,rx_k-rx_n);
	pp1 = pow(2,-rx_k);
	pp2 = pow(2,-rx_n);

	pr = (1 + pp0 + pp1 + pp2)/4;
	ek[0]= -log(pr) / log(2);

	pr = (1 - pp0 + pp1 - pp2)/4;
	ek[1]= -log(pr) / log(2);

	pr = (1 + pp0 - pp1 - pp2)/4;
	ek[2]= -log(pr) / log(2);

	pr = (1 - pp0 - pp1 + pp2)/4;
	ek[3]= -log(pr) / log(2);

	return 1;
}

//计算指定旋转参数k，所对应的Ek的4个值中的最小值
float Compute_RX_Ek_min(float *ek)
{
	float ek_min = 0;
	float ek_tmp = 0;
	u16 i = 0;

	ek_min = ek[0];
	ek_tmp = ek[0];
	for(i=1; i<4; i++)
	{
		if(ek_tmp > ek[i])
		{
			ek_tmp = ek[i];
		}
	}
	ek_min = ek_tmp;
	return ek_min;
}


//计算指定旋转参数k，构造RX_offset对应的u和v，通过指针返回
u16 Compute_RX_V_U(u16 rx_k, u16 rx_n,
		u64 *rx_v, float *rx_v_w, u64 *rx_u,float *rx_u_w)
{
	u16 i=0, j=0;
	u64 tmp=0;
	u64 tmp_rx=0;
	u64 f =1;

	if(rx_k == 1) //rx_k >=1
	{
		rx_u[0]= 1;
		rx_u_w[0]=0;  //k=1时，默认概率重量为0

		rx_v[0] = ROTATE_LEFT(f, rx_k, rx_n);  //n-k=1时，默认概率重量为0
		rx_v_w[0]=1;  //n-k=1时，默认概率重量为1
		//f= ROTATE_LEFT(f, 1, rx_n);
		for(i=1; i< rx_n-1; i++)
		{
			rx_v[i] = ROTATE_LEFT(rx_v[0], i, rx_n) ^ rx_v[i-1];
		}
		for(i=0; i< rx_n-rx_k-1; i++)
		{
			rx_v_w[i]= i+1;
		}
		rx_v_w[rx_n-rx_k-1] = rx_n-2;
	}
	else
	{
		rx_u[0]= 1;
//		rx_u_w[0]=1;  //k=1时，默认概率重量为0
		for(i=1; i < rx_k; i++)
		{
			rx_u[i] = ROTATE_LEFT(rx_u[0], i, rx_n) ^ rx_u[i-1];
		}
		for(i=0; i < rx_k-1; i++)
		{
			rx_u_w[i] = i+1;
		}
		rx_u_w[rx_k-1] = rx_k-1;

		if((rx_n - rx_k) == 1)  //对于V，若n-k=1的特殊情况
		{
			rx_v[0] = ROTATE_LEFT(f, rx_k, rx_n);  //n-k=1时，默认概率重量为0
			rx_v_w[0]= 0;
		}
		else //对于V，若n-k !=1的一般情况
		{
			rx_v[0] = ROTATE_LEFT(f, rx_k, rx_n);  //n-k !=1时
			for(i=1; i< rx_n-rx_k; i++)
			{
				rx_v[i] = ROTATE_LEFT(rx_v[0], i, rx_n) ^ rx_v[i-1];
			}

			for(i=0; i< rx_n-rx_k-1; i++)
			{
				rx_v_w[i]= i+1;
			}
			rx_v_w[rx_n-rx_k-1] = rx_n-rx_k-1;
		}
	}

/*
	rx_u[0]= 1;
	for(i=1; i < rx_k; i++) //构造U的值
	{
		rx_u[i]= ROTATE_LEFT(1, i, rx_n);
		rx_u[i] ^= rx_u[i-1];
	}
	rx_u_w[0]=0;  //k=1时，默认概率重量为0
	for(i=0; i < rx_k-1; i++) //构造U的值对应的概率重量
	{
		rx_u_w[i] = i+1;
	}
	if(rx_k -1 > 0)
	{
		rx_u_w[rx_k-1] = rx_k-1; //1^k对应的概率重量为k-1
	}
//%%%%%%%%%%%%%%%%%%%%%%%以下为V和V对应概率重量的构造%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	tmp_rx = ROTATE_LEFT(1, rx_k, rx_n); //V的值，为n比特，低k位全为0,后面可以直接用
	rx_v[i] = 1; //n-k=1时，默认概率重量为1
	for(i=1; i < rx_n - rx_k; i++) //构造V的值的数组,需要注意，V的值，为n比特，低k位全为0
	{
		rx_v[i]= ROTATE_LEFT(tmp_rx, i, rx_n);
		rx_v[i] ^= rx_v[i-1];
	}
	rx_v_w[0]=0;  //n-k=1时，默认概率重量为0
	for(i=0; i < rx_n - rx_k -1; i++) //构造V的值对应的概率重量
	{
		rx_v_w[i]=i+1;
	}
	if( rx_n - rx_k -1 > 0)
	{
		rx_v_w[rx_n - rx_k-1] = rx_n - rx_k -1; //1^{n-k}对应的概率重量为n-k-1
	}
*/
	return 1;
}



//基于冒泡排序:从小到大排列，最小下标冒泡排序
//预计算RX_offset对应的按照zeta的概率重量递增的排序表
//返回u64 *zeta, float *zeta_w, u16 *zeta_index
u16 Sorted_RX_offset_Table(u16 rx_k, u16 rx_n,
		u64 *rx_v, float *rx_v_w, u64 *rx_u, float *rx_u_w,
		float *ek,
		u64 *zeta, float *zeta_w, u16 *zeta_index)
{
	u16 i=0,j=0;
	u16 ua=0, vb=0;
	u64 tmp=0,tmp_i=0;
	float zeta_tmp = 0;
	u16 t_index = 0;
	u16 l=0,m=0,flag=0;



//预处理，将所有zeta的可能取值，全部存入zeta和zeta_w中
	//(CL,CR)=(0,0)
	*zeta++   = 0;
	*zeta_w++ = ek[0]; // printf("Ek[0][%d]: %f  \n",0,Ek[0]);

	//(CL,CR)=(0,1)
	if(rx_k==1)
	{
		*zeta++   = rx_u[0];
		*zeta_w++ = ek[1] + rx_u_w[0];
	//printf("zeta_w[%d]: %f  RX_u_w[0]:%f ?? \n", 1,Ek[1],RX_u_w[0]);
	}
	else
	{
		for(i=0; i< rx_k; i++)
		{
			*zeta++   = rx_u[i];
			*zeta_w++ = ek[1] + rx_u_w[i];
	}}

	//(CL,CR)=(1,0)
	if((rx_n-rx_k)==1)
	{
		*zeta++   = rx_v[0];
		*zeta_w++ = ek[2] + rx_v_w[0];
	//printf("zeta_w[%d]: %f  RX_u_w[0]:%f ?? \n", 1,Ek[1],RX_u_w[0]);
	}
	else
	{
		for(i=0; i< rx_n-rx_k; i++)
		{
			*zeta++   = rx_v[i];
			*zeta_w++ = ek[2] + rx_v_w[i];
		}
	}
   //(CL,CR)=(1,1)
//	i=0;
//	for(i=rx_n; i< rx_n + rx_k*(rx_n -rx_k); i++)
//	{
		for(ua=0; ua < rx_k; ua++)
			for(vb=0; vb < rx_n-rx_k; vb++)
			{
				tmp = rx_v[vb] ^ rx_u[ua];
				*zeta++   = tmp;
				*zeta_w++ = ek[3] + rx_v_w[vb] + rx_u_w[ua];
//				i++;
			}
//	}

//		for(i=0; i<1+ rx_n + rx_k*(rx_n -rx_k); i++)printf("zeta_w[%d]: %f  \n",i, zeta_w[i]);


//%%%%%%%%%%%%%%%%%%%%%以下开始为基于冒泡排序法，按照zeta_w由小到大的顺序进行排序%%%%%%%%%%%%//
		for(i=0; i < RX_zeta_total; i++)
		{
			RX_zeta_i[i] = i;
		}

		for(i=0; i < RX_zeta_total-1; i++)
		{
			t_index = i;
			for(j=i+1; j < RX_zeta_total; j++)
			{
				if(RX_zeta_w[t_index] > RX_zeta_w[j])
				{
					t_index = j;
				}
			}
			zeta_tmp = RX_zeta_w[i];
			RX_zeta_w[i] = RX_zeta_w[t_index];
			RX_zeta_w[t_index] = zeta_tmp;

	        tmp_i = RX_zeta_i[i];
	        RX_zeta_i[i] = RX_zeta_i[t_index];
	        RX_zeta_i[t_index] = tmp_i;
		}
/*
			for(i=0; i<63; i++)
			{
				printf("RX_zeta[%d]: %x RX_zeta_w[%d]: %f  RX_zeta_i[%d]: %d \n",
						i,RX_zeta[i], i,RX_zeta_w[i], i,  RX_zeta_i[i] );
			}
			printf("|--------------------------------------------------|\n");
			for(i=0; i<63; i++)
			{
			printf("RX_zeta[%d]: %x RX_zeta_w[%d]: %f  RX_zeta_i[%d]: %d \n",
					RX_zeta_i[i],RX_zeta[RX_zeta_i[i]], i,RX_zeta_w[i], i,  RX_zeta_i[i] );
			}
*/
	//zeta_idex中记录的值，表示概率重量递增的下标顺序
	return 1;
}





















