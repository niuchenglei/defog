#include "StdAfx.h"
#include "Matrix.h"
#include <math.h>

CMatrix::CMatrix(void)
{
}

CMatrix::~CMatrix(void)
{
}

/////////////////////////////////////////////////
/*矩阵加法
//rows:行
//cols:列*/
////////////////////////////////////////////////
 float* CMatrix::PlusMatrix(float *A1,float *A2,int rows,int cols)
 {
	 //矩阵相减
	 float *B = new float[rows*cols];
	 for (int i = 0; i < rows; i++)
	 {
		 for (int j = 0; j < cols;j++ )
			{
				B[i*cols+j] = A1[i*cols+j] + A2[i*cols+j];
			}

	 }
	 return B;
 }

 /////////////////////////////////////////////////
 /*矩阵加法
 //rows:行
 //cols:列*/
 ////////////////////////////////////////////////
 double* CMatrix::PlusMatrix(double *A1,double *A2,int rows,int cols)
 {
	 double *B = new double[rows*cols];
	 for (int i = 0; i < rows; i++)
	 {
		 for (int j = 0; j < cols;j++ )
			{
				B[i*cols+j] = A1[i*cols+j] + A2[i*cols+j];
			}

	 }
	 return B;
 }
/////////////////////////////
/* 矩阵减法
//a1:A1矩阵行；
//a2： A2矩阵行
//b1: A1矩阵列
//b2: B2矩阵列*/
///////////////////////////////
 float* CMatrix::SubMatrix(float *A1, float *A2,int a1,int a2,int b1,int b2)
{
	//判断矩阵的长短是否一致
	/*if ((a1 != a2))
	{
		return NULL;
	}*/
	if (a2 ==1)
	{
		//矩阵相减
		float *B = new float[a1*b1];
		for (int i = 0; i < a1; i++)
		{
			for (int j = 0; j < b1;j++ )
			{
				B[i*b1+j] = A1[i*b1+j] - A2[j];
			}

		}
		return B;

	}
	else if (b2 ==1)/////////////////////////////////减去列向量也可以。
	{
		//矩阵相减
		float *B = new float[a1*b1];
		for (int i = 0; i < a1; i++)
		{
			for (int j = 0; j < b1;j++ )
			{
				B[i*b1+j] = A1[i*b1+j] - A2[j];
			}

		}
		return B;
	}

	else
	{
		//矩阵相减
		float *B = new float[a1*b1];
		for (int i = 0; i < a1; i++)
		{
			for (int j = 0; j < b1;j++ )
			{
				B[i*b1+j] = A1[i*b1+j] - A2[i*b1+j];
			}

		}
		return B;
	}
	
	
}


 double* CMatrix::SubMatrix(double *A1, double *A2,int a1,int a2,int b1,int b2)
 {
	 //判断矩阵的长短是否一致
	 /*if ((a1 != a2))
	 {
	 return NULL;
	 }*/
	 if (a2 ==1)
	 {
		 //矩阵相减
		 double *B = new double[a1*b1];
		 for (int i = 0; i < a1; i++)
		 {
			 for (int j = 0; j < b1;j++ )
			 {
				 B[i*b1+j] = A1[i*b1+j] - A2[j];
			 }

		 }
		 return B;

	 }
	 else if (b2 ==1)/////////////////////////////////减去列向量也可以。
	 {
		 //矩阵相减
		 double *B = new double[a1*b1];
		 for (int i = 0; i < a1; i++)
		 {
			 for (int j = 0; j < b1;j++ )
			 {
				 B[i*b1+j] = A1[i*b1+j] - A2[j];
			 }

		 }
		 return B;
	 }

	 else
	 {
		 //矩阵相减
		 double *B = new double[a1*b1];
		 for (int i = 0; i < a1; i++)
		 {
			 for (int j = 0; j < b1;j++ )
			 {
				 B[i*b1+j] = A1[i*b1+j] - A2[i*b1+j];
			 }

		 }
		 return B;
	 }


 }

 //////////////////////////////////////
 /// 矩阵乘法
 /// firstMatrix,secondMatrix:两个输入矩阵
 /// h1:第一矩阵行
 /// w1:第一矩阵列
 /// h2:第二矩阵行
 /// w2:第二矩阵列
 /////////////////////////////////////
  float* CMatrix::MultiplyMatrix(float *firstMatrix, float *secondMatrix,int h1,int w1,int h2,int w2 )
 {
	 
	 float *resultMatrix = new float[h1*w2];
	 //判断相乘矩阵是否合法，即第一个矩阵的列要等于第二个矩阵的行
	 if (w1 != h2)
	 {
		 return NULL;
	 }
	 //求结果矩阵
	 for (int rowIndex = 0; rowIndex < h1; rowIndex++)
	 {
		 for (int colIndex = 0; colIndex < w2; colIndex++)
		 {
			 //初始化结果矩阵的元素
			 resultMatrix[rowIndex*w2+colIndex] = 0;

			 for (int i = 0; i < w1; i++)
			 {
				 //求结果矩阵的元素值
				 resultMatrix[rowIndex*w2+colIndex] += firstMatrix[rowIndex*w1+ i] * secondMatrix[i*w2 +colIndex];
			 }
		 }
	 }
	 return resultMatrix;
 }


  double* CMatrix::MultiplyMatrix(double *firstMatrix, double *secondMatrix,int h1,int w1,int h2,int w2 )
  {

	  double *resultMatrix = new double[h1*w2];
	  //判断相乘矩阵是否合法，即第一个矩阵的列要等于第二个矩阵的行
	  if (w1 != h2)
	  {
		  return NULL;
	  }
	  //求结果矩阵
	  for (int rowIndex = 0; rowIndex < h1; rowIndex++)
	  {
		  for (int colIndex = 0; colIndex < w2; colIndex++)
		  {
			  //初始化结果矩阵的元素
			  resultMatrix[rowIndex*w2+colIndex] = 0;

			  for (int i = 0; i < w1; i++)
			  {
				  //求结果矩阵的元素值
				  resultMatrix[rowIndex*w2+colIndex] += firstMatrix[rowIndex*w1+ i] * secondMatrix[i*w2 +colIndex];
			  }
		  }
	  }
	  return resultMatrix;
  }


  ///////////////////////////////////////
  /// 求转置矩阵
  /// row:矩阵行
  /// col:矩阵列
  //////////////////////////////////////
   float*  CMatrix::TransposeMaxtrix(float* dMatrix,int row,int col)
  {
	
	  float *resultMaxtrix = new float[col*row];

	  for (int i = 0; i < row;i++ )
	  {
		  for (int j = 0; j < col;j++ )
		  {
			  resultMaxtrix[j*row+i] = dMatrix[i*col+j];
		  }
	  }
	  return resultMaxtrix;
  }


   ///////////////////////////////////////
   /// 求转置矩阵
   /// row:矩阵行
   /// col:矩阵列
   //////////////////////////////////////
   double*  CMatrix::TransposeMaxtrix(double* dMatrix,int row,int col)
   {

	   double * resultMaxtrix = new double[col*row];

	   for (int i = 0; i < row;i++ )
	   {
		   for (int j = 0; j < col;j++ )
		   {
			   resultMaxtrix[j*row+i] = dMatrix[i*col+j];
		   }
	   }
	   return resultMaxtrix;
   }

   /////////////////////////////////////////////
   ////求逆矩阵
   /////////////////////////////////////////////
    float*  CMatrix::ReverseMatrix(float* dMatrix, int Level)
   {
	   double dMatrixValue = MatrixValue( dMatrix, Level );
	   if( dMatrixValue == 0 ) return NULL;
	   float* dReverseMatrix = new float[Level*2*Level];
	   float x, c;
	   // Init Reverse matrix 
	   for( int i = 0; i < Level; i++ )
	   {
		   for( int j = 0; j < 2 * Level; j++ )
		   {
			   if( j < Level )
				   dReverseMatrix[i*2 * Level+j] = dMatrix[i,j];
			   else
				   dReverseMatrix[i*2 * Level+j] = 0;
		   }
		   dReverseMatrix[i*2 * Level+Level + i ] = 1;
	   }
	   for( int i = 0, j = 0; i < Level && j < Level; i++, j++ )
	   {
		   if( dReverseMatrix[i*2 * Level+j] == 0 )
		   {
			   int m = i;
			   for( ; dMatrix[m,j] == 0; m++ );
			   if( m == Level )
				   return NULL;
			   else
			   {

				   // Add i-row with m-row

				   for( int n = j; n < 2 * Level; n++ )

					   dReverseMatrix[i*2 * Level+n] += dReverseMatrix[m*2 * Level+n];
			   }
		   }
		   // Format the i-row with "1" start
		   x = dReverseMatrix[i*2 * Level+j];
		   if( x != 1 )
		   {
			   for( int n = j; n < 2 * Level; n++ )
				   if( dReverseMatrix[i*2 * Level+n] != 0 )
					   dReverseMatrix[i*2 * Level+n] /= x;
		   }
		   // Set 0 to the current column in the rows after current row
		   for( int s = Level - 1; s > i;s-- )
		   {
			   x = dReverseMatrix[s*2 * Level+j];
			   for( int t = j; t < 2 * Level; t++ )
				   dReverseMatrix[s*2 * Level+t] -= ( dReverseMatrix[i*2 * Level+t]* x );
		   }
	   }
	   // Format the first matrix into unit-matrix

	   for( int i = Level - 2; i >= 0; i-- )
	   {
		   for( int j = i + 1 ; j < Level; j++ )
			   if( dReverseMatrix[i*2 * Level+j] != 0 )
			   {
				   c = dReverseMatrix[i*2 * Level+j];
				   int n;
				   for(  n = j; n < 2*Level; n++ )
				   {
					   dReverseMatrix[i*2 * Level+n] -= ( c * dReverseMatrix[j*2 * Level+n] );
				   }
			   }
	   }

	   float* dReturn = new float[Level, Level];
	   for( int i = 0; i < Level; i++ )
	   {
		   for( int j = 0; j < Level; j++ )
		   {
			   dReturn[i*Level+j] = dReverseMatrix[i*2 * Level+j+Level];
		   }
	   }
	   return dReturn;

   }



   /////////////////////////////////////////
   ///矩阵行列式
   ///level：矩阵列
   /////////////////////////////////////////
    double CMatrix::MatrixValue( float* MatrixList, int Level )
   {
	   double* dMatrix = new double[Level*Level];
	   for( int i = 0; i < Level; i++ )
	   {
		   for (int j = 0; j < Level; j++)
		   {
			   dMatrix[i*Level+j] = MatrixList[i*Level+j];
		   }
	   }
	   float c, x;
	   double k = 1;
	   for( int i = 0, j = 0; i < Level && j < Level; i++, j++ )
	   {
		   if( dMatrix[i*Level+j] == 0 )
		   {
			   int m = i;
			   for( ; dMatrix[m*Level+j] == 0; m++ );
			   if( m == Level )
				   return 0;
			   else
			   {
				   // Row change between i-row and m-row
				   for( int n = j; n < Level; n++ )
				   {
					   c = dMatrix[i*Level+n];
					   dMatrix[i*Level+n] = dMatrix[m*Level+n];
					   dMatrix[m*Level+n] = c;
				   }

				   // Change value pre-value
				   k *= (-1);
			   }
		   }
		   // Set 0 to the current column in the rows after current row
		   for( int s = Level - 1; s > i;s-- )
		   {
			   x = dMatrix[s*Level+j];
			   for( int t = j; t < Level; t++ )
				   dMatrix[s*Level+t] -= dMatrix[i*Level+t]* ( x/dMatrix[i*Level+j] );
		   }
	   }
	   double sn = 1;
	   for( int i = 0; i < Level; i++ )
	   {
		   if( dMatrix[i*Level+i] != 0 )
			   sn *= dMatrix[i*Level+i];
		   else
			   return 0;
	   }
	   return k*sn;
   }

	 //////////////////////////////////////
	////*计算协方差矩阵
	////dMatrix:
	////bands:
	////num:
	////mean[]:
	////////////////////////////////////
	float *  CMatrix::CovMatrix(float *dMatrix,int bands,int num,float* meanMatrix)
	{
		float *dtMatrix = TransposeMaxtrix(dMatrix,bands,num);//transpose data matrix
		
		int i,j,k;
	
		float * subMatrix = SubMatrix(dtMatrix,meanMatrix,num,1,bands,bands);//sub the mean value
	//	float * subMatrix = SubMatrix(dtMatrix,meanMatrix,num,num,bands,bands);
		
		float * result = new float[bands*bands];//store final result
		for (i = 0;i<bands;i++)
		{
			for(j = 0;j<bands;j++)
			{
				int firstNum = i;
				int secondNum= j;
				float firstCol ;
				float secondCol;
				double sum = 0 ;
				for (k = 0;k<num;k++)
				{
					firstCol=subMatrix[k*bands + i];
					secondCol=subMatrix[k*bands+ j];
					sum +=double(firstCol*secondCol);
				}
				result[i*bands+j]=sum/(num-1);
                
			}
		}

         return result;
	}


	/************************************************************************/
	/* a为输入矩阵，n为输入矩阵维数   ,矩阵求逆                                                                  */
	/************************************************************************/
    float* CMatrix::rinv(float *a, int n)
	{ 
		int *is,*js,i,j,k,l,u,v;
		double d,p;
		is=(int*)malloc(n*sizeof(int));
		js=(int*)malloc(n*sizeof(int));
		for (k=0; k<=n-1; k++)
		{ 
			d=0.0;
			for (i=k; i<=n-1; i++)
				for (j=k; j<=n-1; j++)
				{
					l=i*n+j; p=fabs(a[l]);
					if (p>d)
					{ 
						d=p; is[k]=i; js[k]=j;
					}
				}
				if (d+1.0==1.0)
				{ 
					free(is); free(js); printf("err**not inv\n");
					return(0);
				}
				if (is[k]!=k)
					for (j=0; j<=n-1; j++)
					{ 
						u=k*n+j; v=is[k]*n+j;
						p=a[u]; a[u]=a[v]; a[v]=p;
					}
				if (js[k]!=k)
					for (i=0; i<=n-1; i++)
					{
						u=i*n+k; v=i*n+js[k];
						p=a[u]; a[u]=a[v]; a[v]=p;
					}
					l=k*n+k;
					a[l]=1.0/a[l];
					for (j=0; j<=n-1; j++)
					if (j!=k)
					{ 
						u=k*n+j; a[u]=a[u]*a[l];
					}
					for (i=0; i<=n-1; i++)
					if (i!=k)
					for (j=0; j<=n-1; j++)
						if (j!=k)
						{
							u=i*n+j;
							a[u]=a[u]-a[i*n+k]*a[k*n+j];
						}
					for (i=0; i<=n-1; i++)
						if (i!=k)
						{
							u=i*n+k; a[u]=-a[u]*a[l];
						}
		}
		for (k=n-1; k>=0; k--)
		{
			if (js[k]!=k)
			for (j=0; j<=n-1; j++)
			{
				u=k*n+j; v=js[k]*n+j;
				p=a[u]; a[u]=a[v]; a[v]=p;
			}
			if (is[k]!=k)
				for (i=0; i<=n-1; i++)
				{
					u=i*n+k; v=i*n+is[k];
					p=a[u]; a[u]=a[v]; a[v]=p;
				}
		}
		free(is); free(js);
		return a;
	}



	/************************************************************************/
	/*计算对称矩阵的特征值和特征向量。输入矩阵为data[dimension][dimension];
	返回值存储在矩阵eigenvector中，eigenvector为特征向量矩阵，b为特征值矩阵*/
	/************************************************************************/
	float *CMatrix::Eigenvec(float *data,float *eigenvector, float *eigenvalue, int dimension)
	{

		float *c = (float*)malloc(sizeof(float) * dimension);//三对角元素的副对角线元素
		eigenvector = strq(data, dimension, eigenvector, eigenvalue, c);
		eigenvector = sstq(dimension, eigenvalue, c, eigenvector);
		free(c);
		return eigenvector;
	}
	/************************************************************************/
	/* 计算对称阵特征值特征向量first step：将对称阵约化为三对角矩阵（私有） ，利用household方法
	data[dimension][dimension]:输入的对称阵
	q：householder 变换矩阵
	b：存储主对角线矩阵
	c：存储副对角线矩阵
	输出：q*/
	/************************************************************************/
	float *CMatrix::strq(float *data, int dimension,float *q,float *b, float *c)
	{ 
		int i, j, k, u;
		float h, f, g, h2;
		for (i = 0; i<= dimension-1; i++ )
			for (j = 0; j<= dimension-1; j++ )
			{ 
				u = i *  dimension + j;
				q[u] = data[u];
			}
			for (i = dimension-1; i>= 1; i--)
			{ 
				h = 0.0;
				if (i>1)
					for (k = 0; k<= i-1; k++ )
					{
						u = i *  dimension + k;
						h = h + q[u] * q[u];
					}
					if (h + 1.0 == 1.0)
					{ 
						c[i] = 0.0;
						if (i == 1) 
							c[i] = q[i * dimension + i-1];
						b[i] = 0.0;
					}
					else
					{
						c[i] = sqrt(h);
						u = i * dimension + i-1;
						if (q[u] > 0.0) c[i] =- c[i];
						h = h-q[u] * c[i];
						q[u] = q[u] - c[i];
						f = 0.0;
						for (j = 0; j<= i-1; j++ )
						{
							q[j * dimension + i] = q[i * dimension + j] / h;
							g = 0.0;
							for (k = 0; k<= j; k++ )
								g = g + q[j * dimension + k] *  q[i * dimension + k];
							if (j + 1<= i - 1)
								for (k = j + 1; k<= i - 1; k++ )
									g = g + q[k * dimension + j] *  q[i * dimension + k];
							c[j] = g / h;
							f = f + g * q[j * dimension + i];
						}
						h2 = f / (h + h);
						for (j = 0; j<= i-1; j++ )
						{
							f = q[i *  dimension + j];
							g = c[j]-h2 * f;
							c[j] = g;
							for (k = 0; k<= j; k++ )
							{
								u = j * dimension + k;
								q[u] = q[u] - f * c[k]-g * q[i * dimension + k];
							}
						}
						b[i] = h;
					}
			}
			for (i = 0; i<= dimension-2; i++ ) 
				c[i] = c[i + 1];
			c[dimension-1] = 0.0;
			b[0] = 0.0;
			for (i = 0; i<= dimension-1; i++ )
			{ 
				if ((b[i]!= 0.0)&&(i-1 >= 0))
					for (j = 0; j<= i-1; j++ )
					{ 
						g = 0.0;
						for (k = 0; k<= i-1; k++ )
							g = g + q[i * dimension + k] *  q[k * dimension + j];
						for (k = 0; k<= i-1; k++ )
						{ 
							u = k * dimension + j;
							q[u] = q[u]-g * q[k * dimension + i];
						}
					}
					u = i * dimension + i;
					b[i] = q[u]; q[u] = 1.0;
					if (i-1>= 0)
						for (j = 0; j<= i-1; j++ )
						{
							q[i * dimension + j] = 0.0; q[j * dimension + i] = 0.0;
						}
			}
			return q;
	}
	/************************************************************************/
	/*  j计算三对角矩阵的特征值和特征向量second step。
	dimension为q的维数。
	参数如上个函数*/
	/************************************************************************/
	float *CMatrix::sstq(int dimension, float *b,float *c,float *q)
	{
		int i, j, k, m, it, u, v;
		float d, f, h, g, p, r, e, s;
		c[dimension-1] = 0.0; d = 0.0; f = 0.0;
		for (j = 0; j <= dimension-1; j++)
		{ 
			it = 0;
			h =  0.000001 * (fabs(b[j]) + fabs(c[j]));
			if (h > d)
				d = h;
			m = j;
			while ((m <= dimension-1) && (fabs(c[m]) > d)) 
				m = m + 1;
			if (m != j)
			{ 
				do
				{
					/*if (it==l)
					{ 
					printf("fail\n");
					return(NULL);
					}*/
					it = it + 1;
					g = b[j];
					p = (b[j + 1]- g) / (2.0 * c[j]);
					r = sqrt(p *  p + 1.0);
					if (p>= 0.0) 
						b[j] = c[j] / (p + r);
					else
						b[j] = c[j] / (p-r);
					h = g - b[j];
					for (i = j + 1; i<= dimension-1; i ++ )
						b[i] = b[i]-h;
					f = f + h; p = b[m]; e = 1.0; s = 0.0;
					for (i = m-1; i>= j; i--)
					{
						g = e *  c[i]; h = e *  p;
						if (fabs(p)>= fabs(c[i]))
						{
							e = c[i] / p; r = sqrt(e *  e + 1.0);
							c[i + 1] = s *  p *  r; s = e / r; e = 1.0 / r;
						}
						else
						{ 
							e = p / c[i]; r = sqrt(e *  e + 1.0);
							c[i + 1] = s *  c[i] *  r;
							s = 1.0 / r; e = e / r;
						}
						p = e *  b[i]-s *  g;
						b[i + 1] = h + s *  (e *  g + s *  b[i]);
						for (k = 0; k<= dimension-1; k++ )
						{ 
							u = k *  dimension + i + 1; v = u-1;
							h = q[u]; q[u] = s *  q[v] + e *  h;
							q[v] = e *  q[v]-s *  h;
						}
					}
					c[j] = s * p; b[j] = e * p;
				}
				while (fabs(c[j]) > d);
			}
			b[j] = b[j] + f;
		}
		for (i = 0; i<= dimension-1; i++ )
		{ 
			k = i; p = b[i];
			if (i + 1<= dimension-1)
			{ 
				j = i + 1;
				while ((j<= dimension-1)&&(b[j]<= p))
				{ k = j; p = b[j]; j = j + 1;}
			}
			if (k!= i)
			{
				b[k] = b[i]; b[i] = p;
				for (j = 0; j<= dimension-1; j++ )
				{ 
					u = j *  dimension + i; v = j *  dimension + k;
					p = q[u]; q[u] = q[v]; q[v] = p;
				}
			}
		}

		return q;
	}

