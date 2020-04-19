// 注意：rgb色彩信息按BSQ方式排列，即r,g,b三个矩阵分开排列

#include "StdAfx.h"
#include "Guide_Fliter.h"
//#include <Windows.h>

Guide_Fliter::Guide_Fliter(void)
{
	r_win = 20;
	yipu = 0.001;
}

Guide_Fliter::~Guide_Fliter(void)
{
}
//函数名称：integralImage
//功    能：计算累积图
//调用函数：无
//被调函数：main(), windowAverage, windowVar
//参数说明：
//输    入：输入图像pic
//输    出：累积图p_sum
void Guide_Fliter::integralImage(float *pic,double *p_sum)
{
	unsigned long i1,j1;

	p_sum[0] = pic[0];

	// 计算第一列
	for (j1=1;j1<height;j1++)
	{
		p_sum[j1*width] = p_sum[(j1-1)*width] + pic[j1*width];
	}

	// 计算第一行
	for (j1=1;j1<width;j1++)
	{
		p_sum[j1] = p_sum[j1-1] + pic[j1];
	}

	for (i1=1;i1<height;i1++)
	{
		for (j1=1;j1<width;j1++)
		{
			p_sum[i1*width+j1] = pic[i1*width+j1]+p_sum[(i1-1)*width+j1]
			+p_sum[i1*width+(j1-1)]-p_sum[(i1-1)*width+(j1-1)];
		}
	}

}

//函数名称：dim1wAvg
//功    能：计算窗口内像素的均值
//调用函数：integralImage
//调用函数：windowAverage
//参数说明：
//输    入：图像累积图ps，图像的宽和高分别为w1、h1，窗口半径r
//输    出：均值图ave
//方法说明：
//1.在计算平均值时，未延拓图像，直接在图像的内部求和平均
//2.分别按9种情况处理：中间、四个角、四条边
//3.中心为(i1,j1)、半径为r时，window的覆盖范围是[i1-r,i1+r]×[j1-r,j1+r]，尺寸为(2r+1)×(2r+1)

void Guide_Fliter::dim1winAvg(float *pic,float *avg)
{
	long i1,j1,w1,h1,r;
	long x1,x2,y1,y2;
	double temp =0.0;
	double *ps =NULL;

	w1 = width;
	h1 = height;
	r = r_win;

	if ( r>=(w1/2) || r>=(h1/2) )
	{
		printf("窗口半径太大！\n\n");
		exit(1);
	}
	ps = (double*)malloc(w1*h1*sizeof(double));
	integralImage(pic,ps);  // 累积图

	for (i1=0;i1<h1;i1++)
	{
		y1 = i1-r-1;
		y2 = i1+r;

		for (j1=0;j1<w1;j1++)
		{
			x1 = j1-r-1;
			x2 = j1+r;

			if(y1<0)
			{
				if(x1<0) temp = ps[y2*w1+x2]/((x2+1)*(y2+1));// 左上角
				else if(x2 >=w1) temp = (ps[y2*w1+w1-1]-ps[y2*w1+x1])/((w1-x1-1)*(y2+1));// 右上角
				else temp = (ps[y2*w1+x2]-ps[y2*w1+x1])/((x2-x1)*(y2+1));// 上边
			}
			else if(y2 >= h1)
			{
				if(x2>=w1) temp = (ps[h1*w1-1]+ps[y1*w1+x1]-ps[y1*w1+w1-1]-ps[(h1-1)*w1+x1])/((w1-1-x1)*(h1-1-y1));// 右下角
				else if(x1<0) temp = (ps[(h1-1)*w1+x2] -ps[y1*w1+x2])/((x2+1)*(h1-1-y1));// 左下角
				else temp = (ps[(h1-1)*w1+x2]+ps[y1*w1+x1]-ps[y1*w1+x2]-ps[(h1-1)*w1+x1])/((x2-x1)*(h1-1-y1));//下边
			}
			else
			{
				if(x1<0)  temp = (ps[y2*w1+x2] -ps[y1*w1+x2])/((x2+1)*(y2-y1));// 左边
				else if(x2>=w1) temp = (ps[y2*w1+w1-1]+ps[y1*w1+x1]-ps[y1*w1+w1-1]-ps[y2*w1+x1])/((w1-1-x1)*(y2-y1));// 右边
				else temp = (ps[y2*w1+x2]+ps[y1*w1+x1]-ps[y1*w1+x2]-ps[y2*w1+x1])/((x2-x1)*(y2-y1));// 中间
			}
			avg[i1*w1+j1] = (float)temp;
		}
	} // end of for_loop

	free(ps);
	ps = NULL;
}

void Guide_Fliter::windowAverage(float *pic,float *p_avg)
{
	dim1winAvg(&pic[0],&p_avg[0]);                             // r分量
	dim1winAvg(&pic[width*height],&p_avg[width*height]);       // g分量
	dim1winAvg(&pic[width*height*2],&p_avg[width*height*2]);   // b分量
}


//函数名称：dim1winCov
//功    能：计算窗口内像素的方差
//调用函数：integralImage，dim1winAvg
//被调函数：windowCov
//参数说明：略
//方法说明：为方便使用累积图方法，需要将方差变形

void Guide_Fliter::dim1winCov(float *pic1,float *avg1,float *pic2,float *avg2,float *cov)
{
	unsigned long i1,k1;
	float *quad_avg= NULL,*pic_quad=NULL;

	quad_avg =(float*)malloc(width*height*sizeof(float));
	pic_quad =(float*)malloc(width*height*sizeof(float));

	k1 = width*height;
	
	for (i1=0;i1<k1;i1++)
	{
		pic_quad[i1] = pic1[i1]*pic2[i1];// 计算像素平方图
	}

	dim1winAvg(pic_quad,quad_avg);

	for (i1=0;i1<k1;i1++)
	{
		cov[i1] = quad_avg[i1]- avg1[i1]*avg2[i1];
	}

	free(quad_avg);	quad_avg = NULL;
	free(pic_quad);	pic_quad = NULL;
}


//函数名称：matrix_Cov
//功    能：指导图像的协方差矩阵
//调用函数：dim1winCov
//被调函数：piq_Fliter3
//参数说明：略
//注意事项：本函数用于彩色图像！rgb色彩信息按BSQ方式排列！

void Guide_Fliter::matrix_Cov(float *img,float *i_avg,float *cov)
{
	unsigned long k1,k2,k3,k4,k5,k6,k7,k8,hw;

	k1 = width*height;
	k2 = k1*2;
	k3 = k1*3;
	k4 = k1*4;
	k5 = k1*5;
	k6 = k1*6;
	k7 = k1*7;
	k8 = k1*8;

	dim1winCov(&img[0],&i_avg[0],&img[0],&i_avg[0],&cov[0]);
	dim1winCov(&img[k1],&i_avg[k1],&img[k1],&i_avg[k1],&cov[k4]);
	dim1winCov(&img[k2],&i_avg[k2],&img[k2],&i_avg[k2],&cov[k8]);
	dim1winCov(&img[0],&i_avg[0],&img[k1],&i_avg[k1],&cov[k1]);
	dim1winCov(&img[0],&i_avg[0],&img[k2],&i_avg[k2],&cov[k2]);
	dim1winCov(&img[k1],&i_avg[k1],&img[k2],&i_avg[k2],&cov[k5]);

	for (hw=0;hw<k1;hw++)
	{
		cov[k3+hw] = cov[k1+hw];
		cov[k6+hw] = cov[k2+hw];
		cov[k7+hw] = cov[k5+hw];
	}
}

//函数名称：winCov_Ip
//功    能：计算mat(img协方差矩阵+yiputh*U)的逆矩阵：inv_mat
//调用函数：无
//被调函数：piq_Fliter3
//参数说明：
//输    入: cov_img、yipu；
//输    出: inv_mat
//方法说明：针对3维对称实数矩阵，采用代数余子式方法计算

void Guide_Fliter::inv_r3Mat(float *cov_img,float *inv_mat)
{

	float m[9],mval;
	unsigned long i1,k1,k2,k3,k4,k5,k6,k7,k8;

	k1 = width*height;
	k2 = k1*2;
	k3 = k1*3;
	k4 = k1*4;
	k5 = k1*5;
	k6 = k1*6;
	k7 = k1*7;
	k8 = k1*8;

	for (i1=0;i1<k1;i1++)
	{
		// mat=cov_img(img协方差矩阵+yiputh*U)
		m[0]=cov_img[i1] +yipu;
		m[1]=cov_img[i1+k1];
		m[2]=cov_img[i1+k2];
		m[3]=cov_img[i1+k3];
		m[4]=cov_img[i1+k4]+yipu;
		m[5]=cov_img[i1+k5];
		m[6]=cov_img[i1+k6];
		m[7]=cov_img[i1+k7];
		m[8]=cov_img[i1+k8]+yipu;
		
		// ===================矩阵求逆=======================
		mval =  m[0]*m[4]*m[8]+m[1]*m[5]*m[6]+m[2]*m[3]*m[7]
		       -m[2]*m[4]*m[6]-m[0]*m[5]*m[7]-m[1]*m[3]*m[8];

		if (mval ==0.0)
		{
			printf("矩阵不可逆\n");
			exit(1);
		}
		else
		{ 
			inv_mat[i1]    =(m[4]*m[8]-m[5]*m[7])/mval;
			inv_mat[i1+k1] =(m[5]*m[6]-m[3]*m[8])/mval;
			inv_mat[i1+k2] =(m[3]*m[7]-m[4]*m[6])/mval;

			inv_mat[i1+k3] =inv_mat[i1+k1];
			inv_mat[i1+k4] =(m[0]*m[8]-m[2]*m[6])/mval;
			inv_mat[i1+k5] =(m[1]*m[6]-m[0]*m[7])/mval;

			inv_mat[i1+k6] =inv_mat[i1+k2];
			inv_mat[i1+k7] =inv_mat[i1+k5];
			inv_mat[i1+k8] =(m[0]*m[4]-m[1]*m[3])/mval;
		}// ================================================
	}
}


//函数名称：winCov_Ip
//功    能：计算指导图像与输入图像的窗内协方差
//调用函数：dim1winCov
//被调函数：piq_Fliter1
//参数说明：输入->略；输出->无返回
//方法说明：ak = covIp/(cov_img+yipu)
//注意事项：本函数用于彩色图像！rgb色彩信息按BSQ方式排列！
void Guide_Fliter::winCov_Ip(float *pic,float *p_avg,float *img,float *i_avg,float *covIp)
{
	unsigned long k2,k1;

	k1 = width*height;
	k2 = 2*k1;
	
	dim1winCov(img,i_avg,pic,p_avg,covIp);
	dim1winCov(&img[k1],&i_avg[k1],pic,p_avg,&covIp[k1]);
	dim1winCov(&img[k2],&i_avg[k2],pic,p_avg,&covIp[k2]);
}


//函数名称：qout_Caculater1
//功    能：指导滤波计算
//调用函数：dim1winAvg
//被调函数：piq_Fliter1
//参数说明：输入->略；输出->无返回
//方法说明：ak = covIp/(cov_img+yipu)
//注意事项：rgb色彩信息按BSQ方式排列！
void Guide_Fliter::qout_Caculater1(float *img,float *i_avg,float *p_avg,float *cov_img,float *covIp,float *q_out)
{
	unsigned long i1,k1= width*height;
	float *ak,*bk,*ak_avg,*bk_avg;

	ak     = (float*)malloc(width*height*sizeof(float));
	bk     = (float*)malloc(width*height*sizeof(float));
	ak_avg = (float*)malloc(width*height*sizeof(float));
	bk_avg = (float*)malloc(width*height*sizeof(float));

	for (i1=0;i1<k1;i1++)
	{
		ak[i1] = covIp[i1]/(cov_img[i1]+yipu);
		bk[i1] = p_avg[i1]-ak[i1]*i_avg[i1];
	}

	dim1winAvg(ak,ak_avg);                // ak均值
	dim1winAvg(bk,bk_avg);                // bk均值

	for (i1=0;i1<k1;i1++)
	{
		q_out[i1] = ak_avg[i1]*img[i1]+bk_avg[i1];
	}

	free(ak);ak = NULL;
	free(bk);bk = NULL;
	free(ak_avg);ak_avg = NULL;
	free(bk_avg);bk_avg = NULL;

}


//函数名称：qout_Caculater3
//功    能：指导滤波计算
//调用函数：windowAverage
//被调函数：piq_Fliter3
//参数说明：输入->略；输出->无返回
//方法说明：ak = inv_mat*covIp
//注意事项：rgb色彩信息按BSQ方式排列！

void Guide_Fliter::qout_Caculater3(float *img,float *i_avg,float *p_avg,float *inv_mat,float *covIp,float *q_out)
{
	unsigned long i1,k1,k2,k3,k4,k5,k6,k7,k8;
	float *ak,*bk,*ak_avg,*bk_avg;

	k1 = width*height;
	k2 = k1*2;
	k3 = k1*3;
	k4 = k1*4;
	k5 = k1*5;
	k6 = k1*6;
	k7 = k1*7;
	k8 = k1*8;

	ak     = (float*)malloc(width*height*3*sizeof(float));
	ak_avg = (float*)malloc(width*height*3*sizeof(float));
	bk     = (float*)malloc(width*height*sizeof(float));
	bk_avg = (float*)malloc(width*height*sizeof(float));

	for (i1=0;i1<k1;i1++)
	{
		ak[i1]    = inv_mat[i1]*covIp[i1]+inv_mat[i1+k1]*covIp[i1+k1]+inv_mat[i1+k2]*covIp[i1+k2];
		ak[i1+k1] = inv_mat[i1+k3]*covIp[i1]+inv_mat[i1+k4]*covIp[i1+k1]+inv_mat[i1+k5]*covIp[i1+k2];
		ak[i1+k2] = inv_mat[i1+k6]*covIp[i1]+inv_mat[i1+k7]*covIp[i1+k1]+inv_mat[i1+k8]*covIp[i1+k2];

		bk[i1] = p_avg[i1]-(ak[i1]*i_avg[i1]+ak[i1+k1]*i_avg[i1+k1]+ak[i1+k2]*i_avg[i1+k2]);
	}

	windowAverage(ak,ak_avg);             // ak均值
	dim1winAvg(bk,bk_avg);                // bk均值

	for (i1=0;i1<k1;i1++)
	{
		q_out[i1] = ak_avg[i1]*img[i1]+ak_avg[i1+k1]*img[i1+k1]+ak_avg[i1+k2]*img[i1+k2]+bk_avg[i1];
	}

	free(ak);ak = NULL;
	free(bk);bk = NULL;
	free(ak_avg);ak_avg = NULL;
	free(bk_avg);bk_avg = NULL;

}


//函数名称：piq_Fliter1
//功    能：实现灰度图像的指导滤波
//调用函数：dim1winAvg,dim1winCov,qout_Caculater1
//被调函数：piq_GuideFliter
//参数说明：
//输    入：输入图像*pic,指导图像 *img,指针 *q_out
//输    出：滤波结果，float型指针q_out

float* Guide_Fliter::piq_Fliter1(float *pic,float *img,float *q_out)
{
	float *p_avg,*i_avg,*cov_img,*covIp;

	p_avg       = (float*)malloc(width*height*sizeof(float));
	i_avg       = (float*)malloc(width*height*sizeof(float));

	cov_img     = (float*)malloc(width*height*sizeof(float));
	covIp       = (float*)malloc(width*height*sizeof(float));

	dim1winAvg(pic,p_avg);
	dim1winAvg(img,i_avg);

	dim1winCov(img,i_avg,img,i_avg,cov_img);
	dim1winCov(pic,p_avg,img,i_avg,covIp);

	qout_Caculater1(img,i_avg,p_avg,cov_img,covIp,q_out);//计算输出图像q

	free(p_avg);p_avg = NULL;
	free(covIp);covIp = NULL;
	free(i_avg);i_avg = NULL;
	free(cov_img);cov_img = NULL;

	return q_out;
}


//函数名称：piq_Fliter3
//功    能：实现3通道彩色图像的指导滤波
//调用函数：windowAverage,matrix_Cov,winCov_Ip,inv_r3Mat,qout_Caculater3
//被调函数：piq_GuideFliter
//参数说明：
//输    入：输入图像*pic,指导图像 *img,指针 *q_out
//输    出：滤波结果，float型指针q_out

float* Guide_Fliter::piq_Fliter3(float *pic,float *img,float *q_out)
{
	float *p_avg,*i_avg,*cov_img,*inv_mat,*covIp;

	p_avg       = (float*)malloc(width*height*sizeof(float));
	i_avg       = (float*)malloc(width*height*3*sizeof(float));

	cov_img     = (float*)malloc(width*height*9*sizeof(float));
	inv_mat     = (float*)malloc(width*height*9*sizeof(float));
	covIp       = (float*)malloc(width*height*3*sizeof(float));

	dim1winAvg(pic,p_avg);
	windowAverage(img,i_avg);

	matrix_Cov(img,i_avg,cov_img);
	winCov_Ip(pic,p_avg,img,i_avg,covIp);
	inv_r3Mat(cov_img,inv_mat);// 计算mat(img协方差矩阵+yiputh*U)的逆矩阵：inv_mat

	qout_Caculater3(img,i_avg,p_avg,inv_mat,covIp,q_out);//计算输出图像q

	free(p_avg);p_avg = NULL;
	free(covIp);covIp = NULL;
	free(i_avg);i_avg = NULL;
	free(cov_img);cov_img = NULL;
	free(inv_mat);inv_mat = NULL;

	return q_out;
}


//函数名称：piq_GuideFliter
//功    能：指导滤波计算
//调用函数：piq_Fliter1,piq_Fliter3
//被调函数：main()
//参数说明：
//输    入：按次序输入，指导图像I、输入图像p、输出图像q的文件路径
//          图像宽度w、图像高度h、色彩通道数d
//输    出：float型指针q_out，q_out指向指导滤波图像的数据
//注意事项：rgb色彩信息按BSQ方式排列！

float* Guide_Fliter::piq_GuideFliter(float *pic,float *img,float *q_out,unsigned long w,unsigned long h,int d)
{
	float *ptemp = NULL;

	width = w;
	height = h;

	if (d==1)
	{
		ptemp = piq_Fliter1(pic,img,q_out);
	}
	else if (d==3)
	{
		ptemp =piq_Fliter3(pic,img,q_out);
	} 
	else
	{	
		printf("色彩通道数有误，程序无法处理\n");
	}
	return ptemp;
}
