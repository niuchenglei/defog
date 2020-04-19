#pragma once


class Guide_Fliter
{
public:
	Guide_Fliter(void);
	~Guide_Fliter(void);
private:
	unsigned long width, height,r_win;
	float yipu;
private:
	// ====================滤波计算函数case1：灰度=============================
	void integralImage(float *pic,double *p_sum);
	void dim1winAvg(float *pic,float *avg);
	void dim1winCov(float *pic1,float *avg1,float *pic2,float *avg2,float *cov);
	void qout_Caculater1(float *img, float *i_avg, float *p_avg, float *cov_img, float *covIp, float *q_out);
	float* piq_Fliter1(float *pic,float *img,float *q_out);

	// ====================滤波计算函数case2：24位彩色========================
	void windowAverage(float *pic,float *avg);
	void matrix_Cov(float *pic1,float *avg1,float *cov);
	void inv_r3Mat(float *cov_img,float *inv_mat);
	void winCov_Ip(float *pic1,float *avg1,float *pic2,float *avg2,float *cov);
	void qout_Caculater3(float *img, float *i_avg, float *p_avg, float *inv_mat, float *covIp, float *q_out);
	float* piq_Fliter3(float *pic,float *img,float *q_out);
	//==========================================================================
public:
	float* piq_GuideFliter(float *pic,float *img,float *q_out,unsigned long w,unsigned long h,int d);
};

