//Copyright@ 
|				@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
|				 Name: auto_search
|				 Function:  Searching DC&LC&RX&RK of ARX ciphers 
|				 Version: V1.0 
|				 Created since: 2017-7-28 
|			         Author: Mingjiang Huang 
|			         Email:	huangmingjiang@iie.ac.cn 
|				 Organization:	Chinese Academy of Sciences 
|				 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
/*
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 *      ///cDDT with 8-bit, lookup \gamma, based on \alpha and \beta
 */
void ARX_cDDT_construct(void)
{
	u16 A,B,C;  //ALPHA,BETA
	u16 AB = 0;
	u16 carrybits = 0;
	u16 alpha_inv = 0;
	u16 beta_inv = 0;
	u16 gamma_inv = 0;
	u16 alpha_tmp = 0;
	u16 beta_tmp = 0;
	u16 gamma_tmp= 0;
	u16 eq = 0xFFFF;
	//u64 num = 0;
	u16 wt = 0;
	u16 wt_cnt[9] = {0};
	u16 i = 0;
	u16 WT = 0;
	u16 WT_cnt[9] = {0};
	u8 flag1 = 0,flag2 = 0;


	memset(cDDT_v,0,sizeof(MSB_cDDT_v)); //内存清零函数
	memset(cDDT_n,0,sizeof(MSB_cDDT_n)); //内存清零函数
	memset(MSB_cDDT_v,0,sizeof(MSB_cDDT_v)); //内存清零函数
	memset(MSB_cDDT_n,0,sizeof(MSB_cDDT_n)); //内存清零函数

	memset(cDDT_wt_min,0xFF,sizeof(cDDT_wt_min)); //初始化cDDT_wt_min为全F。
	memset(cDDT_wt_max,0,sizeof(cDDT_wt_max)); //初始化MSB_cDDT_wt_min为全F。
	memset(MSB_cDDT_wt_min,0xFF,sizeof(MSB_cDDT_wt_min)); //初始化MSB_cDDT_wt_min为全F。
	memset(MSB_cDDT_wt_max,0,sizeof(MSB_cDDT_wt_max)); //内存清零函数

	for(A=0;A<=0xFF;A++)   //8bit DDTA
	{
		alpha_tmp = A<<1;
		for(B=0;B<=0xFF;B++)
		{
			AB = ((A << 8) | B);
			beta_tmp = B<<1;
			for(carrybits=0; carrybits<8; carrybits++)
			{
				for(i=0; i<9;i++ )
				{
					wt_cnt[i] = 0;
					WT_cnt[i] = 0;
				}
			for(C=0; C<=0xFF;C++ )
			{
				gamma_tmp = C << 1;

				if((carrybits & 0x4) != 0) // alpha  LSB bit.
				{
					alpha_inv = ~(alpha_tmp | 1) ; // & 0xFFFF;
				}
				else
				{
					alpha_inv = ~alpha_tmp ; // & 0xFFFF;
				}

				if((carrybits & 0x2) != 0) // beta LSB bit.
				{
					beta_inv = (beta_tmp | 1) ; // & 0xFFFF;
				}
				else
				{
					beta_inv = beta_tmp; // & 0xFFFF;
				}

				if((carrybits & 0x1) != 0) // gamma LSB bit.
				{
					gamma_inv = (gamma_tmp | 1);
				}
				else
				{
					gamma_inv = gamma_tmp;
				}

				eq = (alpha_inv ^ beta_inv) & (alpha_inv ^ gamma_inv);
				eq = eq & (A ^ B ^ C ^ beta_inv);

				if( (eq & 0xFF) == 0) //
				{
					//compute wt= Speck_block_wt_compute(A,B,C)
					wt = Speck_block_wt_compute(A,B,C,0xFF);
					cDDT_v[carrybits][AB][wt][wt_cnt[wt]] = (u16)C;  //for gamma numb.
					wt_cnt[wt]++ ;

					WT = Speck_block_wt_compute(A,B,C,0x7F); //HM_weight((~((~A ^ B) & (~A ^ C))) & 0x7F);
					MSB_cDDT_v[carrybits][AB][WT][WT_cnt[WT]] = (u16)C;  //for gamma numb.
					WT_cnt[WT]++;
				}
			}

			flag1 = 0;
			flag2 = 0;
			for(i=0; i<9;i++ )
			{
				cDDT_n[carrybits][AB][i] = wt_cnt[i];
				MSB_cDDT_n[carrybits][AB][i] = WT_cnt[i];

				if(wt_cnt[i] != 0)
				{
					if(flag1 ==0)
					{
						cDDT_wt_min[carrybits][AB] = i;
						flag1 = 1;
					}
					cDDT_wt_max[carrybits][AB] = i;
				}

				if(WT_cnt[i] != 0)
				{
					if(flag2 ==0)
					{
						MSB_cDDT_wt_min[carrybits][AB] = i;
						flag2 = 1;
					}
					MSB_cDDT_wt_max[carrybits][AB] = i;
				}

			}

			if(cDDT_AB_wt_min[AB] > cDDT_wt_min[carrybits][AB])
			{
				cDDT_AB_wt_min[AB] = cDDT_wt_min[carrybits][AB];
			}

			//printf("cDDT[carrybits][AB][256]: %d \n",cDDT[carrybits][AB][256]);
			}
		}
	}
}


/*
 *
 *  Created on: 2017年8月22日
 *      Author: hmj110131
 *      ///Another kind of cDDT with 8-bit, lookup \gamma and \beta, based on \alpha.
 */
void fixed_Alpha_get_betagamma(void)
{
	u16 A,B,C;  //ALPHA,BETA
	u16 carrybits = 0;
	u16 alpha_inv = 0;
	u16 beta_inv = 0;
	u16 gamma_inv = 0;
	u16 alpha_tmp = 0;
	u16 beta_tmp = 0;
	u16 gamma_tmp= 0;
	u16 eq = 0xFFFF;
	//u64 num = 0;
	u16 wt = 0;
	u16 wt_cnt[9] = {0};
	u16 i = 0;
	u16 WT = 0;
	u16 WT_cnt[9] = {0};
	u8 flag1 = 0,flag2 = 0;

	memset(CHAM_a_beta,0,sizeof(CHAM_a_beta)); //内存清零函数
	memset(CHAM_a_gama,0,sizeof(CHAM_a_gama)); //内存清零函数
	memset(CHAM_a_bg_numb,0,sizeof(CHAM_a_bg_numb)); //内存清零函数
	memset(CHAM_a_wt_min,0xFF,sizeof(CHAM_a_wt_min)); //内存清零函数

	memset(msb_CHAM_a_beta,0,sizeof(msb_CHAM_a_beta)); //内存清零函数
	memset(msb_CHAM_a_gama,0,sizeof(msb_CHAM_a_gama)); //内存清零函数
	memset(msb_CHAM_a_bg_numb,0,sizeof(msb_CHAM_a_bg_numb)); //内存清零函数
	memset(msb_CHAM_a_wt_min,0xFF,sizeof(msb_CHAM_a_wt_min)); //内存清零函数


	for(A=0;A<=0xFF;A++)   //8bit DDTA
	{
		alpha_tmp = A<<1;

		for(carrybits=0; carrybits<8; carrybits++)
		{
			for(i=0; i<9;i++ )
			{
				wt_cnt[i] = 0;
				WT_cnt[i] = 0;
			}

			for(B=0;B<=0xFF;B++)
			{
				//AB = ((A << 8) | B);
				beta_tmp = B<<1;

			for(C=0; C<=0xFF;C++ )
			{
				gamma_tmp = C << 1;

				if((carrybits & 0x4) != 0) // alpha  LSB bit.
				{
					alpha_inv = ~(alpha_tmp | 1) ; // & 0xFFFF;
				}
				else
				{
					alpha_inv = ~alpha_tmp ; // & 0xFFFF;
				}

				if((carrybits & 0x2) != 0) // beta LSB bit.
				{
					beta_inv = (beta_tmp | 1) ; // & 0xFFFF;
				}
				else
				{
					beta_inv = beta_tmp; // & 0xFFFF;
				}

				if((carrybits & 0x1) != 0) // gamma LSB bit.
				{
					gamma_inv = (gamma_tmp | 1);
				}
				else
				{
					gamma_inv = gamma_tmp;
				}

				eq = (alpha_inv ^ beta_inv) & (alpha_inv ^ gamma_inv);
				eq = eq & (A ^ B ^ C ^ beta_inv);

				if( (eq & 0xFF) == 0) //
				{
					//compute wt= Speck_block_wt_compute(A,B,C)
					wt = Speck_block_wt_compute(A,B,C,0xFF);
					CHAM_a_beta[A][carrybits][wt][wt_cnt[wt]] = (u8)B;
					CHAM_a_gama[A][carrybits][wt][wt_cnt[wt]] = (u8)C;
					wt_cnt[wt]++ ;

					WT = Speck_block_wt_compute(A,B,C,0x7F); //HM_weight((~(((~A) ^ B) & ((~A) ^ C))) & 0x7F);
					msb_CHAM_a_beta[A][carrybits][WT][WT_cnt[WT]] = (u8)B;
					msb_CHAM_a_gama[A][carrybits][WT][WT_cnt[WT]] = (u8)C;
					WT_cnt[WT]++;
			}}}


			flag1 = 0;
			flag2 = 0;
			for(i=0; i<9;i++ )
			{
				CHAM_a_bg_numb[A][carrybits][i] = wt_cnt[i];
				msb_CHAM_a_bg_numb[A][carrybits][i] = WT_cnt[i];

				if(wt_cnt[i] != 0)
				{
					if(flag1 ==0)
					{
						CHAM_a_wt_min[A][carrybits] = i;
						flag1 = 1;
					}
					CHAM_a_wt_max[A][carrybits] = i;
				}

				if(WT_cnt[i] != 0)
				{
					if(flag2 ==0)
					{
						msb_CHAM_a_wt_min[A][carrybits] = i;
						flag2 = 1;
					}
					msb_CHAM_a_wt_max[A][carrybits] = i;
				}
			}

			//printf("cDDT[carrybits][AB][256]: %d \n",cDDT[carrybits][AB][256]);
		}
	}
}



