/*
 * printinfo.c
 *
 *  Created on: 2017年8月23日
 *      Author: hmj110131
 */


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include "ciphers.h"
#include "search.h"
#include "typedef.h"
#include "globalvar.h"
#include "printinfo.h"



/**
 * @Print the final result of DC or LC.
 *
 */
void print_resoult(u16 search_round)
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



void test_print(void)
{

  int i=0,n,a[32];
  u16 l=0, r=0;
  printf("请输入一个十进制整数.\n");
  scanf("%d",&n);
  while (n>0)
  {
      a[i]=n%2;
      i=i+1;
      n=n/2;
  }
  printf("十进制整数转换为二进制数是:\n");
  for(i--;i>=0;i--)
      printf("%d",a[i]);
  printf("\n");

  r = ((l)<<1) | (l)>>15;

  printf(" %x   %x \n ", l,r);
 }






