//#include "StdAfx.h"
#include <math.h>
#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include "tinex.h"
#define  PI 3.1415926

//Retinex Ä£ÐÍÈ¥Îí
// ŸßÓÐRetinexºÍ¶à³ß¶ÈRetinexÁœÖÖ·œ·š


/***************************************************************/
/*ŽËº¯ÊýÊµÏÖÍŒÏñµÄ²¹ÁãµÄ¹ŠÄÜ£¬Ê¹µÃÍŒÏñµÄ³€¿íÊÇ2µÄÕûÊýŽÎÃÝ*/
unsigned long sizeCanculator(unsigned long height)
{

	unsigned long numh;
	numh=1;
	while(numh < height)
	{
		numh *= 2;
	}
	return numh;
}

///////////////////////////////////////////////////////////////////////////
//ŽËº¯ÊýÀŽÊµÏÖ¶þžµÀïÒ¶±ä»»
//Ô­ÊŒÊýŸÝdataÒÔžŽÊýÐÎÊœ±íÊŸ£¬ÓÃÒ»Î¬µÄdoubleÏòÁ¿±£Žæ£¬ÆäÖÐ2rÎªÊµ²¿£¬2r+1ÎªÐé²¿¡£
//×îÖÕ±äÁ¿µÄœá¹û±£ŽæÔÚdataÖÐ
//height,width ·Ö±ð±íÊŸÍŒÏñµÄžß¶ÈºÍ¿í¶È
//isign±íÊŸÎªÕý·Ž±ä»»£¬=1Ê±ÊÇÕý±ä»»£¬=-1Ê±Îª·Ž±ä»»
////////////////////////////////////////////////////////////////////////////////////
void Fourier2(double *data,unsigned long height, unsigned long width, int isign)
{
	unsigned long x_index, i_x, i_x_rev, i_y, i_y_rev, jx1, jy1;
	unsigned long k1, k2, k_inblock, k_rev, ntot, block_length, block_totalnum, block_index, i_Lf_block;
	double tempr, tempi, w_temp, theta, wr, wi, wpr, wpi;

	ntot = height * width;//ÍŒÏñÏñËØµã×ÜÊý
	//µûÐÎËã·š²œÖè1£¬±äÖ·
	for(i_x = 0; i_x < width; i_x++)
	{
		x_index = i_x;
		//ŒÆËãi_x(ÏñËØÁÐÏÂ±ê)µÄµ¹ÐòÊýi_x_rev
		i_x_rev=0;

		jx1=width/2;

		//theRew(jx1, width, x_index);

		while(jx1>=1)
		{
			if (x_index >= jx1)
			{
				i_x_rev += (width/2)/jx1;
				x_index -= jx1;
			}
			jx1 /= 2;
		}

		if (i_x<i_x_rev)
		{
			//ÖðÐÐœ»»»Data(k_inblock)ÓëData(k_rev)
			for (i_y=0; i_y<height; i_y++)
			{
				k_inblock = (i_y*width + i_x)*2;
				k_rev = (i_y*width + i_x_rev)*2;

				tempr = data[k_inblock];
				tempi = data[k_inblock+1];
				data[k_inblock] = data[k_rev];
				data [k_inblock+1] = data[k_rev+1];
				data[k_rev] =tempr;
				data[k_rev+1] = tempi;
			}
		}

	}

 //µûÐÎËã·š²œÖè2£¬œ»²æŒÆËã
	//¿é³€¶Èblock_length³õÖµÊÇ2£¬ÒÔºó²»¶Ï³Ë2
	for(block_length=2;block_length <= height;block_length *=2)
	{
		block_totalnum=ntot/block_length;//¿é×ÜžöÊý
		theta = -isign*PI*2/block_length;
		wr = 1.0-2.0*sin(0.5*theta)*sin(0.5*theta);
		wi = sin(theta);
		wpr =1.0;
		wpi =0.0;
		//Öð¿éœøÐÐµÚL²ãµÄµûÐÎŒÆËã£º¿éÄÚÏÂ±êŽÓ0µœ¿é³€¶ÈµÄÒ»°ë
		for(block_index=0; block_index<block_length/2;block_index++)
		{
			//ÏÈ¿éºóÐÐµÄµûÐÎŒÆËã£ši_Lf_block¿é±àºÅË÷Òý£©
				for(i_Lf_block=0;i_Lf_block<block_totalnum;i_Lf_block++)
			{
				k1 = (block_length * i_Lf_block + block_index)*2;
				k2 = k1 +(block_length/2)*2;

				tempr = data[k2]*wpr-data[k2+1]*wpi;
				tempi = data[k2]*wpi+data[k2+1]*wpr;
				data[k2]= data[k1]-tempr;
				data[k2+1]=data[k1+1]-tempi;
				data[k1] += tempr;
				data[k1+1] += tempi;
			}
			//Ðý×ªÒò×ÓÔÚÉÏÒ»ŽÎµÄ»ùŽ¡ÉÏ³ËÒÔµÚL²ãµÄµ¥Î»±ŸÔ­žù
			w_temp = wpr;
			wpr = w_temp*wr-wpi*wi;
			wpi = wpi*wr +w_temp*wi;
		}
	}


	//×ÝÏò±ä»»
	//µûÐÎËã·š²œÖè1£º±äÖ·
	for(i_y=0;i_y<height;i_y++)
	{
		x_index= i_y;
		i_y_rev=0;
		jy1=height/2;
		while(jy1>=1)
		{
			if (x_index>=jy1)
			{
				i_y_rev += (height/2)/jy1;
				x_index -= jy1;
			}
			jy1/=2;
		}
		//=======================================//
		if(i_y < i_y_rev)
		{
			for(i_x=0;i_x<width;i_x++)
			{
				k_inblock=(i_y*width+i_x)*2;
				k_rev=(i_y_rev*width+i_x)*2;

				tempr = data[k_inblock];
				tempi = data[k_inblock+1];
				data[k_inblock]=data[k_rev];
				data[k_inblock+1]=data[k_rev+1];
				data[k_rev] = tempr;
				data[k_rev+1]= tempi;
 
			}
		}
	}
/////////////////////////////////////////////////////
	//////////////µûÐÎËã·š²œÖè2£ºœ»²æŒÆËã
	for(block_length=2; block_length <= height; block_length *= 2)
	{
		theta = -isign*PI*2/block_length;
		wr = 1.0-2.0*sin(0.5*theta)*sin(0.5*theta);
		wi = sin(theta);
		wpr = 1.0;
		wpi = 0.0;
		//Öð¿éœøÐÐµÚL²ãµÄµûÖµŒÆËã
		for(block_index=0;block_index<block_length/2;block_index++)
		{
			
			for (i_y=0;i_y<height;i_y +=block_length)
			{
				for (i_x=0;i_x<width;i_x++)
				{
					k1 = (i_x+(block_index+i_y)*width)*2;
					k2 = k1+(block_length/2)*width*2;
					tempr = data[k2]*wpr-data[k2+1]*wpi;
					tempi = data[k2]*wpi+data[k2+1]*wpr;
					data[k2]= data[k1]-tempr;
					data[k2+1]=data[k1+1]-tempi;
					data[k1] += tempr;
					data[k1+1] += tempi;

				}
			}
			//Ðý×ªÒò×ÓÔÚÉÏÒ»ŽÎµÄ»ùŽ¡ÉÏ³ËÒÔµÚL²ãµÄµ¥Î»±ŸÔ­žù
			w_temp =wpr;
			wpr = w_temp*wr-wpi*wi;
			wpi = wpi*wr+w_temp*wi;
		}
	}
////////////////////////////////////////////////////////////////////////////////
	if (isign == 1)
	{
		//Õýœ»±ä»»œá¹ûµÄËõ·ÅŽŠÀí£ºdata=data/M
		for(i_x=0; i_x<width;i_x++)
		{
			for(i_y=0;i_y<height;i_y++)
			{
				data[(i_y*width+i_x)*2] /=height*width;
				data[(i_y*width+i_x)*2+1] /=height*width;
			}
				
		}
	}
}
/***************************************************************/


/////////////////////////////////////////////////////
/*Guassian function calculation for whole image
  input: height:image height
         width:image width
  output: Guassian value for each point*/
// ŒÆËãžßË¹º¯Êý
// height£ºÍŒÏñžß
// width: ÍŒÏñ¿í
// tao: ÉèÖÃžßË¹º¯ÊýµÄtao²ÎÊý
//////////////////////////////////////////////////////
//高斯卷集核扩展并归一化
float *CalGuassian(int height,int width,int tao)
{
	float *guass=(float*)malloc(sizeof(float)*height*width);
	//CPU function code
	float center_x = (float)width/2;
	float center_y = (float)height/2;
	double sum = 0.0;
	for (int i =0;i<height;i++)
	{
		for (int j = 0;j<width;j++)
		{
			float ty = i - center_y;
			float tx = j - center_x;
//			float tmp = (float)exp(-(double)(ty*ty+tx*tx)/(2*tao*tao))/(2*PI*tao*tao);
			float tmp = (float)exp(-(double)(ty*ty+tx*tx)/(tao*tao));//ŒÆËãžßË¹º¯Êý²ÉÑùÖµ
			guass[i*width+j] = tmp;
			sum +=tmp;
		}
	}

	for(int k =0;k<height*width;k++)
	{
		guass[k]=guass[k]/sum;
	}
	return guass;
}


/////////////////////////////////////////////////////
/*SSR:Single Scale Retinex implementation for defogging 
 input:  pic: original image data
		 height: image height
		 width:image width
		 bandcount:image bands count
 output: float* result*/
//µ¥Ò»RetinexÄ£ÐÍÈ¥Îí
//pic: Ô­ÊŒÍŒÏñ
//height£º ÍŒÏñžß
//width£º  ÍŒÏñ¿í
//bandcount: ÍŒÏñ²š¶ÎÊý
//////////////////////////////////////////////////  
//main
float* OnRetinex(float* pic, int height,int width,int bandcount)
{
	float *result = (float*)malloc(sizeof(float)*height*width*bandcount);
	double log_data;
	double low_I;
	int tao = 200;
	float* guass = CalGuassian(height,width,tao);
//	Cbutter filter;
	unsigned long numh,numw;//²ÎÊý numh,numw·Ö±ð±íÊŸ²¹ÁãºóµÄÍŒÏñµÄžß¶ÈºÍ¿í¶È
	unsigned long i,j;
	double *data,*zero_pic,*d_guass,*zero_guass;//žŽÊý»¯µÃœá¹ûŽæŽ¢ÓÚdataÖÐ£¬2rÎªÊµ²¿£¬2r+1ÎªÐé²¿
	numh = sizeCanculator(height);
	numw = sizeCanculator(width);//ŒÆËã²¹ÁãºóµÄ³ßŽç
	double max=0.0;
	double min=100000;
	zero_pic = (double *)malloc(numh*numw*sizeof(double));
	zero_guass = (double *)malloc(numh*numw*sizeof(double));
	data = (double *)malloc(numh*numw*2*sizeof(double));
	d_guass = (double *)malloc(numh*numw*2*sizeof(double));
/**/
		//²¹Áã
	for (j=0;j<numh;j++)
	{
		for (i=0;i<numw;i++)
		{
			zero_pic[j*numw+i]=0;
			zero_guass[j*numw+i]=0;
			if ((j<height)&&(i<width))
			{
				zero_pic[j*numw+i]=(double)pic[j*width+i];
				zero_guass[j*numw+i]=(double)guass[j*width+i];
			}
		}
	}
	//žŽÊý»¯²¢È¡¶ÔÊý
	for(i=0;i<numh*numw;i++)
	{
		data[2*i]=zero_pic[i];
		data[2*i+1]=0;
		d_guass[2*i]=zero_guass[i];
		d_guass[2*i+1]=0;
	}

	Fourier2(data,numh,numw,1);//žµÀïÒ¶Õý±ä»»
	Fourier2(d_guass,numh,numw,1);
/**/
	//CPU function code
	for (i=0;i<numh;i++)
	{
		for (j=0;j<numw*2;j++)
		{
			data[i*numw*2+j]=data[i*numw*2+j]*d_guass[i*numw*2+j];
			if(data[i*numw*2+j]<0)
			{
				data[i*numw*2+j]=0;
			}
		}
	}

	Fourier2(data,numh,numw,-1);//žµÀïÒ¶·Ž±ä»»£¬µÃµœÈëÉä·ÖÁ¿data
	float *tmp = (float*)malloc(sizeof(float)*height*width);
/**/
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			log_data=log10((double)pic[i*width+j]);
			low_I = log10(data[i*2*numw+j*2]+1);
//			tmp[i*width+j]=(float)exp(log_data-low_I);//»ñÈ¡·ŽÉä·ÖÁ¿tmp
			tmp[i*width+j]=(float)(log_data-low_I);
			if(max<tmp[i*width+j]) 
			{
				max = tmp[i*width+j];
			}
			if (min>tmp[i*width+j])
			{
				min = tmp[i*width+j];
			}
		}
	}
	printf("host :max = %4.4f\tmin = %4.4f", max, min);
/**/
	for (i=0;i<height;i++)
	{
		for (j = 0;j<width;j++)
		{
			result[i*width+j]=(float)(tmp[i*width+j]-min)*255/(max-min);//À­ÉìÏÔÊŸ
		}
	}
/**/
	free(data);
	free(zero_pic);
	free(d_guass);
	free(zero_guass);

	free(tmp);
	return result;
}
////////////////////////////////////////////////
/*HIS2RGB transformation
  type = 1: RGB TO HIS, return I V1 V2
       = -1:HIS TO RGB, return R  G  B*/
/////////////////////////////////////////////////
BYTE* OnHIS(BYTE* data,int height,int width,int type)
{
	BYTE *ISH = (BYTE*)malloc(sizeof(BYTE)*height*width*3);
	switch(type)
	{
		case 1:		
			float I ;
			float v1;
			float v2;
			for (int i=0;i<height*width;i++)
			{
				I = (float)(data[i]+data[height*width+i]+data[height*width*2+i])/sqrt(3.0);
				v1 = (float)(data[i]+data[height*width+i]-2*data[height*width*2+i])/sqrt(6.0);
				v2 = (float)(data[i]-data[height*width+i])/sqrt(2.0);
				ISH[i] = (BYTE)I;
				ISH[height*width+i]= (BYTE)v1;
				ISH[height*width*2+i]=(BYTE)v2;
				/*ISH[height*width+i]= (BYTE)sqrt(v1*v1+v2*v2);
				ISH[height*width*2+i]=(BYTE)atan(v2/v1);*/
			}
			return ISH;
			break;
		case -1:		
			float R;
			float G;
			float B;
			for (int j = 0;j<height*width*3;j++)
			{
				R = (float)data[j]/sqrt(3.0)+data[height*width+j]/sqrt(6.0)+data[height*width*2+j]/sqrt(2.0);
				G = (float)data[j]/sqrt(3.0)+data[height*width+j]/sqrt(6.0)-data[height*width*2+j]/sqrt(2.0);
				B = (float)data[j]/sqrt(3.0)-2*data[height*width+j]/sqrt(6.0);
			
                ISH[j]=(BYTE)R;
				ISH[height*width+j]=(BYTE)G;
				ISH[height*width*2+j]=(BYTE)B;
			   
			}
			return ISH;
			break;
	}
}

///////////////////////////////////////////////////////////
/* Multiple Scale Retinex implementation for defogging
  use three tao value: 20,40,60 to gain different Gaussian function
  three weights: 0.3,0.4,0.3
  input:  data:original image
          height:image height
		  width:image width
  output:MSR result*/
//¶à³ß¶ÈRetinexº¯ÊýÈ¥Îí
//pic: Ô­ÊŒÍŒÏñ
//height: ÍŒÏñžß
//width: ÍŒÏñ¿í
//////////////////////////////////////////////////////////////////
float* OnMSR(float*pic,int height,int width)
{
	float *result = (float*)malloc(sizeof(float)*height*width);
	double log_data;
	double low_I; 
	double mid_I;
	double high_I;

	///////////////////////////////////////////////////////////////
	//1,set guassian parameters and calculate guassian function
	int tao_low = 20;//·Ö±ðÉèÖÃŒÆËãÈýžöžßË¹º¯ÊýµÄtao²ÎÊý
	int tao_mid = 40;
	int tao_high = 60;
	float weight_low = 0.3f;//ÉèÖÃÈýžöretinexº¯ÊýµÄÈšÖµ
	float weight_mid = 0.4f;
	float weight_high = 0.3f;

	float* guass_low = CalGuassian(height,width,tao_low);//ŒÆËãžßË¹º¯ÊýºË
	float* guass_mid = CalGuassian(height,width,tao_mid);
	float* guass_high = CalGuassian(height,width,tao_high);
//	Cbutter filter;
	unsigned long numh,numw;//²ÎÊý numh,numw·Ö±ð±íÊŸ²¹ÁãºóµÄÍŒÏñµÄžß¶ÈºÍ¿í¶È
	unsigned long i,j;
	float *zero_pic,*zero_guass;//žŽÊý»¯µÃœá¹ûŽæŽ¢ÓÚdataÖÐ£¬2rÎªÊµ²¿£¬2r+1ÎªÐé²¿
	double  *d_guass,*d_guass_mid,*d_guass_high;
	numh = sizeCanculator(height);
	numw = sizeCanculator(width);//ŒÆËã²¹ÁãºóµÄ³ßŽç
	float max=0.0;
	float min=100000;
	zero_pic = (float *)malloc(numh*numw*sizeof(float));
	zero_guass = (float *)malloc(numh*numw*sizeof(float));
    float *zero_guass_mid = (float *)malloc(numh*numw*sizeof(float));
	float *zero_guass_high = (float *)malloc(numh*numw*sizeof(float));
	double *data = (double *)malloc(numh*numw*2*sizeof(double));
	d_guass = (double *)malloc(numh*numw*2*sizeof(double)); 
	d_guass_mid = (double *)malloc(numh*numw*2*sizeof(double)); 
	d_guass_high = (double *)malloc(numh*numw*2*sizeof(double)); 
	//²¹Áã
	for (j=0;j<numh;j++)
	{
		for (i=0;i<numw;i++)
		{
			zero_pic[j*numw+i]=0;
			zero_guass[j*numw+i]=0;
			zero_guass_mid[j*numw+i]=0;
			zero_guass_high[j*numw+i]=0;
			if ((j<height)&&(i<width))
			{
				zero_pic[j*numw+i]=pic[j*width+i];
				zero_guass[j*numw+i]=guass_low[j*width+i];
				zero_guass_mid[j*numw+i]=guass_mid[j*width+i];
				zero_guass_high[j*numw+i]=guass_high[j*width+i];
			}
		}
	}
	//žŽÊý»¯²¢È¡¶ÔÊý
	for(i=0;i<numh*numw;i++)
	{
		data[2*i]=(double)zero_pic[i];
		data[2*i+1]=0;
		d_guass[2*i]=(double)zero_guass[i];
		d_guass[2*i+1]=0;
		d_guass_mid[2*i]=(double)zero_guass_mid[i];
		d_guass_mid[2*i+1]=0;
		d_guass_high[2*i]=(double)zero_guass_high[i];
		d_guass_high[2*i+1]=0;
	}
	Fourier2(data,numh,numw,1);//žµÀïÒ¶Õý±ä»»
	Fourier2(d_guass,numh,numw,1);
	Fourier2(d_guass_mid,numh,numw,1);
	Fourier2(d_guass_high,numh,numw,1);
	free(zero_pic);
	free(zero_guass);
	free(zero_guass_mid);
	free(zero_guass_high);
  
	for (i=0;i<numh;i++)
	{
		for (j=0;j<numw*2;j++)
		{
			d_guass[i*numw*2+j]=data[i*numw*2+j]*d_guass[i*numw*2+j];
			d_guass_mid[i*numw*2+j]=data[i*numw*2+j]*d_guass_mid[i*numw*2+j];
			d_guass_high[i*numw*2+j]=data[i*numw*2+j]*d_guass_high[i*numw*2+j];
			if(d_guass[i*numw*2+j]<0)
			{
				d_guass[i*numw*2+j]=0;
			}
			if(d_guass_mid[i*numw*2+j]<0)
			{
				d_guass_mid[i*numw*2+j]=0;
			}
			if(d_guass_high[i*numw*2+j]<0)
			{
				d_guass_high[i*numw*2+j]=0;
			}

		}
	}
	Fourier2(d_guass,numh,numw,-1);//žµÀïÒ¶·Ž±ä»»£¬µÃµœÈëÉä·ÖÁ¿data
	Fourier2(d_guass_mid,numh,numw,-1);
	Fourier2(d_guass_high,numh,numw,-1);
	float *tmp = (float*)malloc(sizeof(float)*height*width);
	for (i=0;i<height;i++)
	{
		for (j=0;j<width;j++)
		{
			log_data=log10((double)pic[i*width+j]);
			low_I = log10(d_guass[i*2*numw+j*2]+1);
			mid_I = log10(d_guass_mid[i*2*numw+j*2]+1);
			high_I = log10(d_guass_high[i*2*numw+j*2]+1);
			tmp[i*width+j]=exp(weight_low*(log_data-low_I)+weight_mid*(log_data-mid_I)+weight_high*(log_data-high_I));//»ñÈ¡·ŽÉä·ÖÁ¿tmp
			if(max<tmp[i*width+j]) 
			{
				max=tmp[i*width+j];
			}
			if (min>tmp[i*width+j])
			{
				min =tmp[i*width+j];
			}
		}
	}
	for (i=0;i<height;i++)
	{
		for (j = 0;j<width;j++)
		{
			result[i*width+j]=(tmp[i*width+j]-min)*250/(max-min);//À­ÉìÏÔÊŸ
		}
	}
	free(d_guass);
	free(d_guass_mid);
	free(d_guass_high);
	free(tmp);
	free(data);
	return result;
}
