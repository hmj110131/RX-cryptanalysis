/*
 * RX_offset_table.c
 *
 *  Created on: 2020年7月7日
 *      Author: hmj110131
 */
#include <stdio.h>
#include "XOR_offset_table.h"
#include "typedef.h"
#include "math.h"
#include "ciphers.h"
#include <nmmintrin.h>
#include "RX_crypt.h"



//Ek里面的值，对应k和n都需要在这里提前手工计算得到并提前存储。
//循环左移位参数k=1，且假设n很大的情况
extern float  Ek[4]={1.415, 1.415, 3, 3};     //指定k的，对应CL，CR取值的4种情况的概率重量的小数值，当k和n确定，Ek为固定值,这里以k=1为例子。

//RX分析中用到的修正量和偏移量，并记录对应的概率重量，提前计算好
extern u64     RX_u[64]={0};     //XOR_offset的右边部分对应的dummy-offset的值
extern float  RX_u_w[64]={0};     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u64     RX_v[64]={0};     //XOR_offset的左边部分对应的dummy-offset的值
extern float  RX_v_w[64]={0};     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u64     RX_zeta[3600]={0};  //XOR_offset偏移量，zeta
extern float  RX_zeta_w[3600]={0};     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u16     RX_zeta_i[3600]={0};      //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量经过排序后的下标
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


//计算指定旋转参数k，构造XOR_offset对应的u和v，通过指针返回
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
//预计算XOR_offset对应的按照zeta的概率重量递增的排序表
//返回u64 *zeta, float *zeta_w, u16 *zeta_index
u16 Sorted_RX_XOR_offset_Table(u16 rx_k, u16 rx_n,
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

//计算概率重量的增量逼近上界
float wt_inc_app(u64 Alpha, u64 Beta, u64 len_mask)
{
	float Bn_inc = 0;
	u64 xor_alpha_beta = 0;

	xor_alpha_beta = (Beta ^ Alpha) & len_mask ;

	Bn_inc = (float)_mm_popcnt_u64(xor_alpha_beta);

	Bn_inc += Compute_RX_Ek_min(Ek);  // 扩展一轮的xdp重量最小值和Ek的最小值为RX差分特征的概率增量上界

	return Bn_inc;
}





//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

/* 冒泡排序:从小到大排列 */
/* 1. 从当前元素起，向后依次比较每一对相邻元素，若逆序则交换 */
/* 2. 对所有元素均重复以上步骤，直至最后一个元素 */
/* 3. 需要记录arr[]排序完成后，原始arr[]中的元素的下标变换，并存入数组Seq[]中 */
/* elemType arr[]: 排序目标数组; int len: 元素个数 */
u64 bubbleSort (elemType arr[], int len) {
    elemType temp;
    elemType *Seq;

    int i, j;
    for (i=0; i<len-1; i++) /* 外循环为排序趟数，len个数进行len-1趟 */
        for (j=0; j<len-1-i; j++) { /* 内循环为每趟比较的次数，第i趟比较len-i次 */
            if (arr[j] > arr[j+1]) { /* 相邻元素比较，若逆序则交换（升序为左大于右，降序反之） */
                temp = arr[j];
                arr[j] = arr[j+1];
                arr[j+1] = temp;
            }
            Seq[len-1-i] =i;
        }
    return *Seq;
}

/*
int[] a={9,5,12,3,8,18,7,1,2,35,6};
for (int i = 0; i < a.length-1; i++) {
	//假定index时最小下标
	int index=i;
	for(int j=i+1;j<a.length;j++){
		if(a[index]>a[j]){
			index=j;//始终让index下标的值最小，此时index=j
		}
	}
	int item=a[i];//将最大值给item，然后互换位置,循环遍历
	a[i]=a[index];
	a[index]=item;
//System.out.println(Arrays.toString(a));
//		最小下标法就是每次拿第一个和后面的去进行比较，如果比后面的大，就互换位置，然后继续之
}
*/
u64 bubbleSort_index (elemType arr[], int len) {
	//假定index时最小下标
    elemType temp;
    u64 *Seq;
    int i=0,j=0;
	int index=i;

	for (int i = 0; i < len-1; i++) {
	for(int j=i+1;j<len;j++){
		if(arr[index]>arr[j]){
			index=j;//始终让index下标的值最小，此时index=j
		}
	}
	int item=arr[i];//将最大值给item，然后互换位置,循环遍历
	arr[i]=arr[index];
	arr[index]=item;

	Seq[i] = index;  //记录下标，由小到大，这里只能记录到0到len-2,最后一个默认是最大的P(cl,cr）+Pv+Pu
	}
    return *Seq;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//






















