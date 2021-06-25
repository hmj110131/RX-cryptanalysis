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
 *      Created on: 2018年7月23日
 *      Author: hmj110131
 *      ///cLAT with 8-bit, lookup \u and \w, based on \v
 *      //构造分块的cLAT
 */
void ARX_cLAT_construct(void)
{
	//u16 e_cor[2];
	//e_cor[0] = 0x0, e_cor[1] =0x01;
	volatile u16 ie = 0, je = 0, k =0;
	u16 U,V,W;
	u16 Z=0, A=0, B =0, C =0;
	u16 i_index = 0;
	u16 Cj[8], MT[8];
	u16 wt_cor = 0;
	u16 flag_A=0xFF, flag_B =0xFF;
	u16 bit_v = 0x01;

	u16 i=0;

/*
	memset(cLAT_W,0,sizeof(cLAT_W)); //内存清零函数
	memset(cLAT_U,0,sizeof(cLAT_U)); //内存清零函数
	memset(cLAT_WU_numb,0,sizeof(cLAT_WU_numb)); //内存清零函数
	memset(cLAT_wtcor_min,0xFF,sizeof(cLAT_wtcor_min)); //内存清零函数
	memset(cLAT_UVW_bro,0,sizeof(cLAT_UVW_bro)); //内存清零函数
*/


	for(k=0; k<8; k++)
	{
		MT[k] = 0;
	}

	for(ie=0; ie<2; ie++) //作为索引值，进行查询,级联上一个block
	{
	for(V=0; V<256; V++) //模加的输入掩码，作为索引值，进行查询 
	{
		cLAT_wtcor_min[V][ie] = 0xFF;
		wt_cor = 0; // 记录V对应的U和W下的相关性重量

		for(W=0;W<256;W++) //模加的另外一个输入掩码，通过查询得到
		{
			for(U=0;U<256;U++) //模加的输出掩码，通过查询得到
			{
				A = U ^ V;
				B = U ^ W;
				C = U ^ V ^ W;
				for(je=0; je<8; je++)
				{
					Cj[je] = (C >> (7-je)) & 0x1;
				}

				wt_cor = 0; // 记录V对应的U和W下的相关性重量
				if(ie == 1) //block的MSB
				{
					wt_cor++;
					Z = 0x80;  // Z | (bit_v << 7);
					MT[0] = 1;  //0+ie, 0bit,对于完整分组MSB的最高位必须为0
				}
				else
				{
					Z = 00;  // Z | (bit_v << 7);
					MT[0] = 0;  //0+ie, 0bit,对于完整分组MSB的最高位必须为0
				}

				for(i_index=1; i_index<8; i_index++) //1-->7 bit
				{ //我们的目标Mt向量数组，向量中由0-->n-1对应MT分块向量的MSB-->LSB的比特值 异或上 ie
					MT[i_index] = (Cj[i_index -1] + MT[i_index -1]) & 0x1; //% 2;
					if(MT[i_index] == 1 ) //产生相关性重量的位置
					{
						wt_cor++;  //wt_cor为 0 -->8
						Z = Z | (bit_v << (7-i_index)); //构造Z=Mnt向量
					}
				}
				Z = Z & 0xFF; //8bit
				flag_A = A & (~(A&Z));
				flag_B = B & (~(B&Z));

				if( (flag_A==0) && (flag_B ==0) ) // 指标函数，有效性判断 
				{
					cLAT_W[V][ie][wt_cor][cLAT_WU_numb[V][ie][wt_cor]] = W;
					cLAT_U[V][ie][wt_cor][cLAT_WU_numb[V][ie][wt_cor]] = U;
					cLAT_WU_numb[V][ie][wt_cor]++;

					//MT[8]
					cLAT_UVW_bro[U][V][W][ie] = (MT[7] + Cj[7]) & 0x01; //% 2; //向右侧下一个分块的递进比特

					if(cLAT_wtcor_min[V][ie] > wt_cor )  //记录V输入掩码下，对应所有UW的最小相关性重量
					{
						cLAT_wtcor_min[V][ie] = wt_cor;
					}
				}
				else
				{
					cLAT_UVW_bro[U][V][W][ie] = 0xFF; //生成对应的UW为无效的输入和输出,判断没通过
				}
	}}
	}}

	/*
		//printf("cLAT_wtcor_min[V][ie]: %d   \n ", cLAT_wtcor_min[0x08][0] );
		printf("cLAT_WU_numb: %d   \n ", cLAT_WU_numb[0x08][0][1] );
		for(i=0; i< cLAT_WU_numb[0x08][0][1]; i++)  //	cLAT_WU_numb: 191
		{
			printf("%d  W: %x \t ",i, cLAT_W[0x08][0][1][i] );
			printf("U: %x \n ", cLAT_U[0x08][0][1][i] );
		}

		//if(V == 0x08)
		//if(cLAT_WU_numb[0x08][0][1] >= 2 )
		{
		//	printf("Z: %x  W: %x  U: %x ie: %d wt_cor: %d \n ", Z, W, U, ie, wt_cor );
		}
		*/

}






