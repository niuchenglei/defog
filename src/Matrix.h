#pragma once
///////////////////////////////
/*矩阵类*/
///////////////////////////////
class CMatrix
{
public:
	CMatrix(void);
	~CMatrix(void);
	 float* SubMatrix(float *A1, float *A2,int a1,int a2,int b1,int b2);//矩阵相减
	 double* SubMatrix(double *A1, double *A2,int a1,int a2,int b1,int b2);

	 double* PlusMatrix(double *A1,double *A2,int rows,int cols);//矩阵相加
	 float* PlusMatrix(float *A1,float *A2,int rows,int cols);

	 float* MultiplyMatrix(float *firstMatrix, float *secondMatrix,int h1,int w1,int h2,int w2 );//矩阵相乘
	 double* MultiplyMatrix(double *firstMatrix, double *secondMatrix,int h1,int w1,int h2,int w2 );

	 float* TransposeMaxtrix(float* dMatrix,int row,int col);//矩阵转置
	 double* TransposeMaxtrix(double* dMatrix,int row,int col);

	 double MatrixValue( float* MatrixList, int Level );//行列式

	 float* ReverseMatrix(float* dMatrix, int Level);//矩阵逆

	 float* CovMatrix(float *dMatrix,int bands,int num,float *meanMatrix);//协方差矩阵
	 float* rinv(float *a, int n);//矩阵逆
	 float* Eigenvec(float *data,float *eigenvector, float *eigenvalue, int dimension);//矩阵特征向量
private:
	float* sstq(int dimension, float *b,float *c,float *q);
	float* strq(float *data, int dimension,float *q,float *b, float *c);
};
