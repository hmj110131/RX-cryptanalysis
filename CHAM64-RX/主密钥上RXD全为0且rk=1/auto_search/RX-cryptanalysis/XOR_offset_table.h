/*
 * XOR_offset_table.h
 *
 *  Created on: 2020年7月7日
 *      Author: hmj110131
 */
#ifndef RX_CRYPTANALYSIS_XOR_OFFSET_TABLE_H_
#define RX_CRYPTANALYSIS_XOR_OFFSET_TABLE_H_

#include <stdio.h>
#include "typedef.h"
#include "math.h"
#include "ciphers.h"


////////////////////////////////////////////////////////////////
extern float  Ek[4];     //指定k的，对应CL，CR取值的4种情况的概率重量的小数值，当k和n确定，Ek为固定值。
extern u64  RX_u[64];     //XOR_offset的右边部分对应的dummy-offset的值
extern float  RX_u_w[64];     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u64  RX_v[64];     //XOR_offset的左边部分对应的dummy-offset的值
extern float  RX_v_w[64];     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u64  RX_zeta[3600];  //XOR_offset偏移量，zeta
extern float  RX_zeta_w[3600];     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量
extern u16     RX_zeta_i[3600];     //XOR_offset的右边部分对应的dummy-offset的值各自对应的概率重量经过排序后的下标
extern u64  RX_zeta_total;   //rx-xor表中的zeta的个数
//////////////////////////////////////////////////////////////////
//计算指定旋转参数k，所对应的Ek的4个值，通过指针返回
u16 Compute_RX_Ek(u16 rx_k, u16 rx_n,  float *ek);
//计算指定旋转参数k，所对应的Ek的4个值中的最小值
float Compute_RX_Ek_min(float *ek);
extern float Ek_min;
//计算概率重量的增量逼近上界
float wt_inc_app(u64 Alpha, u64 Beta, u64 len_mask);

//计算指定旋转参数k，构造XOR_offset对应的u和v，通过指针返回
u16 Compute_RX_V_U(u16 rx_k, u16 rx_n,
		u64 *rx_v, float *rx_v_w, u64 *rx_u, float *rx_u_w);

//基于冒泡排序:从小到大排列，最小下标冒泡排序
//预计算XOR_offset对应的按照zeta的概率重量递增的排序表
u16 Sorted_RX_XOR_offset_Table(u16 rx_k, u16 rx_n,
		u64 *rx_v, float *rx_v_w, u64 *rx_u, float *rx_u_w,
		float *ek,
		u64 *zeta, float *zeta_w, u16 *zeta_index);













//////////////////////////冗余的函数
#define ARR_LEN 255 /*数组长度上限*/
#define elemType int /*元素类型*/
u64 bubbleSort (elemType arr[], int len);
u64 bubbleSort_index (elemType arr[], int len);

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//









#endif /* RX_CRYPTANALYSIS_XOR_OFFSET_TABLE_H_ */
