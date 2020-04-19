#include "StdAfx.h"
#include "DarkChannel.h"
#include "Matrix.h"
#include <math.h>
#include <string.h>
#include "Guide_Fliter.h"


void DirectSelectSort(float *data, ssize_t *index, ssize_t size)
{
	int k;
	float temp;
	ssize_t tempId;
	for(int i = 0; i <= (size / 1000) + 1; i++)
	{
		k = i;
		for(int j = i + 1; j < size; j++)
		{
			if(*(data + k) < *(data + j))
				k = j;
		}
		if(k != i)
		{
			temp = *(data + k);
			*(data + k) = *(data + i);
			*(data + i) = temp;
			tempId = *(index + k);
			*(index + k) = *(index + i);
			*(index + i) = tempId;
		}
	}
}


CDarkChannel::CDarkChannel(void)
{
}

CDarkChannel::~CDarkChannel(void)
{
}

//////////////////////////////////////////////////////////
/* calculate dark channel*/
//////////////////////////////////////////////////////////
float* CDarkChannel::CalDarkChannel(float*data,int height,int width,int nbands,int size)
{
	int i,j,m,n;
	float*result = (float*)malloc(sizeof(float)*height*width);

	int h_num = height/size;
	int w_num = width/size;

	for (i=0;i<(h_num+1);i++)
	{
		for (j=0;j<(w_num+1);j++)
		{
			int min_i = i*size;
			int max_i = i*size+size;
			int min_j = j*size;
			int max_j = j*size+size;
			float min_tmp = 1.0;
			float tmp = 1.0;
			for(m = min_i;m<max_i;m++)
			{
				for (n = min_j;n<max_j;n++)
				{
					if (m>=0&&m<height&&n>=0&&n<width)
					{
						if (nbands==3)
						{
							float R = *(data + m*width+n);
							float G = *(data + width * height + m * width + n);
							float B = *(data + 2 * width * height + m * width + n);
							tmp=min(B,min(R,G));
						}
						else if(nbands==1)
						{
							tmp = *(data+m*width+n);
						}
						if (tmp<min_tmp)
						{
							min_tmp = tmp;
						}
					}					
				}
			}
			for(m = min_i;m<max_i;m++)
			{
				for (n = min_j;n<max_j;n++)
				{
					if (m>=0&&m<height&&n>=0&&n<width)
					{
						*(result+m*width+n)=min_tmp;
					}
				}
			}
		}
	}

	return result;
}

//////////////////////////////////////////////////////
/*estimate atmosphere light*/
//////////////////////////////////////////////////////////
float* CDarkChannel::GetLight(float* channel,float*data,int height,int width,int nbands)
{

    float max_channel = 0;
	int max_i;
		float* light=(float*)malloc(sizeof(float) * 3);
	ssize_t *index = (ssize_t *)malloc(height * width *sizeof(ssize_t));


	for(ssize_t i = 0; i < height * width; i++)
		*(index + i) = i;

	DirectSelectSort(channel, index, height * width);

	float maxVal = 0;
	int idMax = height * width / 1000;
	if (nbands==3)
	{
		for(int i = 0; i < idMax; i++)
		{
			if(	*(data + *(index + i)) + 
				*(data + *(index + i) + height * width) + 
				*(data + *(index + i) + height * width * 2) > maxVal)
			{
				maxVal = *(data + *(index + i)) + 
					*(data + *(index + i) + height * width) + 
					*(data + *(index + i) + height * width * 2);
				max_i = i;
			}
		}
		light[0] = *(data + *(index + max_i));
		light[1] = *(data + height * width + *(index + max_i));
		light[2] = *(data + 2 * height * width + *(index + max_i));
	}
	else if (nbands==1)
	{
		for(int i = 0; i < idMax; i++)
		{
			if( *(data + *(index + i)) > maxVal )
			{
				maxVal = *(data + *(index + i));
				max_i = i;
			}
		}
		light[0]=*(data + *(index + max_i));
	}
	free(index);
	index = NULL;

	return light;
}

//////////////////////////////////////////////////////////////
/* estimate the transmission*/
//////////////////////////////////////////////////////////////
float* CDarkChannel::CalTransmission(float* channel, float* light,int height,int width,int nbands)
{
	float weight = 0.95f;
	float* Trans = (float*)malloc(sizeof(float)*height*width);
	float tmp_R,tmp_G,tmp_B ;
	if (nbands==3)
	{
		for (int i = 0;i<height;i++)
		{
			for (int j = 0;j<width;j++)
			{
				tmp_R = 1-weight*(*(channel+i*width+j))/light[0];
				tmp_G = 1-weight*(*(channel+i*width+j))/light[1];
				tmp_B = 1-weight*(*(channel+i*width+j))/light[2];
				Trans[i*width+j] = min(tmp_R,min(tmp_G,tmp_B));
			}
		}
	}
	if (nbands==1)
	{
		for (int i = 0;i<height;i++)
		{
			for (int j = 0;j<width;j++)
			{
				tmp_R = 1-weight*(*(channel+i*width+j))/light[0];
				
				Trans[i*width+j] =tmp_R;
			}
		}
	}
	return Trans;
}

/////////////////////////////////////////////////////////////////////
/*recovering the scene radiance*/
////////////////////////////////////////////////////////////////////
BYTE * CDarkChannel::CalRadiance(float* Trans,float *data,float *light,int height,int width,int nbands)
{
	float t0 = 0.1f;
	//float*result = (float*)malloc(sizeof(float)*height*width*3);
	BYTE*result = (BYTE*)malloc(sizeof(BYTE)*height*width*3);
	float result_R,result_G,result_B;
	if (nbands==3)
	{
		for (int i =0;i<height;i++)
		{
			for (int j = 0;j<width;j++)
			{
				result_R=((data[i*width+j])-light[0])/max(Trans[i*width+j],t0)+light[0];
				result_G=((data[width*height+i*width+j])-light[1])
					/max(Trans[i*width+j],t0)+light[1];
				result_B=((data[2*width*height+i*width+j])-light[2])
					/max(Trans[i*width+j],t0)+light[2];
				if (result_R<0)
				{
					result_R=0;
				}
				if (result_G<0)
				{
					result_G=0;
				}
				if (result_B<0)
				{
					result_B=0;
				}
				if (result_R>1)
				{
					result_R=1;
				}
				if (result_G>1)
				{
					result_G=1;
				}
				if (result_B>1)
				{
					result_B=1;
				}
				result[i*width+j] = result_R*255;
				result[width*height+i*width+j] = result_G*255;
				result[2*width*height+i*width+j]=result_B*255;

			}
		}
	}
    if (nbands==1)
    {
		for (int i =0;i<height;i++)
		{
			for (int j = 0;j<width;j++)
			{
				result_R=((data[i*width+j])-light[0])/max(Trans[i*width+j],t0)+light[0];
				if (result_R<0)
				{
					result_R=0;
				}
		
				if (result_R>1)
				{
					result_R=1;
				}
				result[i*width+j] = result_R*255;
			}
		}
    }
	return result;
}

////////////////////////////////////////////////////////
/* dark channel defogging implementation*/
//////////////////////////////////////////////////////
BYTE* CDarkChannel::OnDarkChannel(BYTE*data,int height,int width,int nbands)
{
	int size = 5;
	float* darkchannel = (float*)malloc(sizeof(float)*height*width);
	float* fData = (float*)malloc(sizeof(float)*height*width*nbands);//
    int maxI=0;
	int new_height = height;
	int new_width = width;
    for (int i=0;i<height*width*nbands;i++)
    {
		fData[i]=(float)data[i]/255;
    }

 

	darkchannel = CalDarkChannel(fData,new_height,new_width,nbands,size);
    float* light = GetLight(darkchannel,fData,new_height,new_width,nbands);
    float *Trans = CalTransmission(darkchannel,light,new_height,new_width,nbands);
	float*Trans2=(float*)malloc(sizeof(float)*height*width);
	free(darkchannel);

 	//	Trans = GuidedFilter(Trans,fData,new_height,new_width,nbands);

	Guide_Fliter filter;
	Trans2=filter.piq_GuideFliter(Trans,fData,Trans2,width,height,nbands);
    
	BYTE*result=(BYTE*)malloc(sizeof(BYTE)*new_height*new_width*nbands);



	result = CalRadiance(Trans2,fData,light,new_height,new_width,nbands);


	free(fData);

     return result;

}


/////////////////////////////////////////////////////////////////
/* Using guided filter to do fast matting for transmission*/
///////////////////////////////////////////////////////////////////
float* CDarkChannel::GuidedFilter(float*Trans,float*data,int height,int width,int nbands)
{
	int i,j,m,n;
	int size =20;//
	int allsize = (size*2+1)*(size*2+1);
	float epsilon = 0.0001f;
	float* A = (float*)malloc(sizeof(float)*height*width*nbands);//
	float* B = (float*)malloc(sizeof(float)*height*width*nbands);
	float*result = (float*)malloc(sizeof(float)*height*width*nbands);
	
	float *subpixels_trans = (float*)malloc(sizeof(float)*allsize);//
	
	float *subpixels_R = (float*)malloc(sizeof(float)*allsize);//store sub window pixels
	float *subpixels_G = (float*)malloc(sizeof(float)*allsize);
	float *subpixels_B = (float*)malloc(sizeof(float)*allsize);//store sub window transmission for each band
	float mean_T;
	////////////////////////////////////////////////////////////
	if (nbands==3)
	{
		float mean_R,mean_G,mean_B,var_R,var_G,var_B;
		float sum_R, sum_G,sum_B;
		int index;
		for (i=0;i<height;i++)
		{
			for (j=0;j<width;j++)
		 {
			 int min_y = i-size;
			 int max_y = i+size;
			 int min_x = j-size;
			 int max_x = j+size;
		
			 index = 0;
			 sum_R = 0.0f;
			 sum_G = 0.0f;
			 sum_B = 0.0f;
			 mean_R = 0.0f;
			 mean_G = 0.0f;
			 mean_B = 0.0f;
			 mean_T = 0.0f;
			 var_R = 0.0f;
			 var_G = 0.0f;
			 var_B = 0.0f;
			 for (m = min_y;m<max_y;m++)
			 {
				 for (n = min_x;n<max_x;n++)
				 {
					
					 if (m>0&&m<height&&n>0&&n<width)
					 { 
						 subpixels_R[index]=data[m*width+n];
						 subpixels_G[index]=data[height*width+m*width+n];
						 subpixels_B[index]=data[2*height*width+m*width+n];
						 subpixels_trans[index]=Trans[m*width+n];

						 sum_R += subpixels_R[index]*subpixels_trans[index];
						 sum_G += subpixels_G[index]*subpixels_trans[index];
						 sum_B += subpixels_B[index]*subpixels_trans[index];

						 mean_R += subpixels_R[index];
						 mean_G += subpixels_G[index];
						 mean_B += subpixels_B[index];
						 mean_T += subpixels_trans[index];
						 index++;

					 }
					
				 }
			 }
			 mean_R = mean_R/index;
			 mean_G = mean_G/index;
			 mean_B = mean_B/index;
			
			 float sum = 0.0f;
			 for(int ii =0;ii<index;ii++)
			 {
				 var_R +=(subpixels_R[ii]-mean_R)*(subpixels_R[ii]-mean_R);
				 var_G +=(subpixels_G[ii]-mean_G)*(subpixels_G[ii]-mean_G);
				 var_B +=(subpixels_B[ii]-mean_B)*(subpixels_B[ii]-mean_B);

			 }
			 var_R = var_R/index;
			 var_G = var_G/index;
			 var_B = var_R/index;
			  mean_T = mean_T/index;

			  A[i*width+j] = (sum_R/index-mean_R*mean_T)/(var_R+epsilon);
			  A[height*width+i*width+j] = (sum_G/index-mean_G*mean_T)/(var_G+epsilon);
			  A[2*height*width+i*width+j] = (sum_B/index-mean_B*mean_T)/(var_B+epsilon);
              B[i*width+j] = mean_T - A[i*width+j]*mean_R;
			  B[height*width+i*width+j] = mean_T - A[height*width+i*width+j]*mean_G;
			  B[2*height*width+i*width+j] = mean_T - A[2*height*width+i*width+j]*mean_B;
			  /////////////////////////////////////////////////////////////////////////
             

		 }// end of j
		}//end of i

          /////////////////////////////////////////////q = a*I + b. gain matting result
		  for (i = 0;i<height;i++)
		  {
			  for (j = 0;j<width;j++)
			  {
				  int sum = 0;
				  float sumA_R = 0.0f;
				  float sumB_R = 0.0f;
				  float sumA_G = 0.0f;
				  float sumB_G = 0.0f;
				  float sumA_B = 0.0f;
				  float sumB_B = 0.0f;
				  for(m = (i-size);m<(i+size);m++)
				  {
					  for (n = (j-size);n<(j+size);n++)
					  {
                          if (m>=0&&m<height&&n>=0&&n<width)
                          {
							  sumA_R += A[m*width+n];
							  sumB_R += B[m*width+n];
							  sumA_G += A[height*width+m*width+n];
							  sumB_G += B[height*width+m*width+n];
							  sumA_B += A[2*height*width+m*width+n];
							  sumB_B += B[2*height*width+m*width+n];
							  sum++;
                          }
					  }
				  }

				  result[i*width+j]=sumA_R*data[i*width+j]/sum+sumB_R/sum;
				  result[height*width+i*width+j]=sumA_G*data[height*width+i*width+j]/sum+sumB_G/sum;
				  result[2*height*width+i*width+j]=sumA_B*data[2*height*width+i*width+j]/sum+sumB_B/sum;
			
			  }
		  }

        free(A);
		free(B);
		free(subpixels_B);
		free(subpixels_G);
		free(subpixels_R);
		free(subpixels_trans);

		return result;
	
	}

}
