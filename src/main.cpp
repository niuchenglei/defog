#include "StdAfx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gdal_priv.h"
#include "DarkChannel.h"
#include "tinex.h"

/** /
void sort(float *data, size_t *index, size_t left, size_t right)
{
	float temp = *(data + left);
	float tempId = *(index + left);
	size_t p = left;
	size_t i = left, j = right;
	while(i < j)
	{
		while(j >= p && *(data + j) <= temp)
			j--;
		if(j >= p)
		{
			*(data + p) = *(data + j);
			*(index + p) = *(index + j);
			p = j;
		}
		while(i <= p && *(data + i) >= temp)
			i++;
		if(i <= p)
		{
			*(data + p) = *(data + i);
			*(index + p) = *(index + i);
			p = i;
		}
	}
	*(data + p) = temp;
	*(index + p) = tempId;
	if(p - left > 1)
		sort(data, index, left, p - 1);
	if(right - p > 1)
		sort(data, index, p + 1, right);
}
/**/

int dcmodel(const char* src, const char* dst){
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
/*
	float *data = (float *)malloc(2000 * 1500 * sizeof(float));
	size_t *index = (size_t *)malloc(2000 * 1500 * sizeof(size_t));
	for(size_t i = 0; i < 3000000; i++)
		*(index + i) = i;
	for(size_t i = 0; i < 3000000; i++)
		*(data + i) = (random() / 10000000) / 100.0;

//	sort(data, index, 0, 200 * 150);
	/***************/
	ssize_t nImgSizeX, nImgSizeY, nBandCount;
	double adfGeoTransform[6];	//地理坐标信息
	//char *strFile = "test1.bmp";
	GDALDataset *poDataset; //GDAL数据集
	GDALAllRegister();
	poDataset = (GDALDataset *) GDALOpen(src, GA_ReadOnly );

	if(poDataset == NULL)
	{
		printf("\n文件打开失败\n");
		exit(0);
	}
	nBandCount = poDataset->GetRasterCount();
	nImgSizeX = poDataset->GetRasterXSize();
	nImgSizeY = poDataset->GetRasterYSize();
	poDataset->GetGeoTransform( adfGeoTransform );
	int bandmap[3] = {1,2,3};
	ssize_t height = nImgSizeX, width = nImgSizeY;
	unsigned char *pic = (unsigned char *)malloc(sizeof(char) * nImgSizeX * nImgSizeY * 3);

	poDataset->RasterIO( GF_Read, 0, 0, nImgSizeX, nImgSizeY, pic, nImgSizeX, nImgSizeY, GDT_Byte, 3, bandmap, 0, 0, 0 );

	/***************/
clock_t a = clock();
	CDarkChannel *test = new CDarkChannel();
	pic = test->OnDarkChannel(pic, nImgSizeX, nImgSizeY, nBandCount);
clock_t b = clock();
printf("time = %lf", (double)(b - a)/1000000.0);
	/***************/
	char *strFilePath = "test1.tif", *fomat="GTiff";
	GDALDriver *poDriver;

	poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
	if(poDriver == NULL)
	{
		printf("\n格式不匹配\n");
		exit(0);
	}
	char **papszMetadata = poDriver->GetMetadata();
//	poDataset = poDriver->CreateCopy(strFilePath, poDataset, 1, papszMetadata, NULL, NULL);
	poDataset = poDriver->Create(dst, nImgSizeX, nImgSizeY, nBandCount, GDT_Byte, papszMetadata);	//GDALRasterBand
	poDataset->RasterIO(GF_Write, 0, 0, nImgSizeX, nImgSizeY, pic, nImgSizeX, nImgSizeY, GDT_Byte, nBandCount, 0,0,0,0);

	GDALClose(poDataset);
	/***************/
//	printf("min(4.0, 40.0) : %.2f\t max(4.0, 40.0) : %.2f\t", min(4.0,40.0), max(4.0, 40.0));
/** /
	for(int i = 0; i < 10; i++)
		printf("%4.2f\t", *(data + i));
	sort(data, index, 0, 9);
	printf("\n");
	for(int i = 0; i < 10; i++)
		printf("%.2f\t", *(data + i));
	printf("\n");
	for(int i = 0; i < 10; i++)
		printf("%4d\t", *(index + i));
/**/
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	return 0;
}
int retinex(const char* src, const char* dst){
	
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	int nImgSizeX, nImgSizeY, nBandCount;
	double adfGeoTransform[6];	//地理坐标信息
	char *strFile = "test1.bmp";
	GDALDataset *poDataset; //GDAL数据集
	GDALAllRegister();
	poDataset = (GDALDataset *) GDALOpen(src, GA_ReadOnly );

	if(poDataset == NULL)
	{
		printf("\n文件打开失败\n");
		exit(0);
	}
	nBandCount = poDataset->GetRasterCount();
	nImgSizeX = poDataset->GetRasterXSize();
	nImgSizeY = poDataset->GetRasterYSize();
	poDataset->GetGeoTransform( adfGeoTransform );
	int bandmap[1] = {1};
	int height = nImgSizeX, width = nImgSizeY;
	float *pic = (float *)malloc(sizeof(float) * nImgSizeX * nImgSizeY );
	float *temp = pic;
	unsigned char *test = (unsigned char *)malloc(sizeof(char) * nImgSizeX * nImgSizeY );
	unsigned char *tempTest = test;

	poDataset->RasterIO( GF_Read, 0, 0, nImgSizeX, nImgSizeY, test, nImgSizeX, nImgSizeY, GDT_Byte, 1, bandmap, 0, 0, 0 );

	for(int i = 0; i < height * width; i++)
	{
		*temp++ = (float)*tempTest++;
	}

	//SYS_TIME_USE("cpu time\n")
	printf("cpu time\n");
	{
		temp = OnRetinex(pic, height, width, 1);
	}
	char *strFilePath = "test1.tif", *fomat="GTiff";
	GDALDataset *poDataset2;
	GDALDriver *poDriver;
	tempTest = test;
	for(int i = 0; i < height * width; i++)
	{
		*tempTest++ = (unsigned char)*temp++;
	}

	poDriver = GetGDALDriverManager()->GetDriverByName(fomat);
	if(poDriver == NULL)
	{
		printf("\n格式不匹配\n");
		exit(0);
	}
	char **papszMetadata = poDriver->GetMetadata();
//	poDataset2 = poDriver->CreateCopy(strFilePath, poDataset, 1, papszMetadata, NULL, NULL);
	poDataset2 = poDriver->Create(dst, nImgSizeX, nImgSizeY, nBandCount, GDT_Byte, papszMetadata);	//GDALRasterBand
	poDataset2->RasterIO(GF_Write, 0, 0, nImgSizeX, nImgSizeY, test, nImgSizeX, nImgSizeY, GDT_Byte, nBandCount, 0,0,0,0);

/**/
	GDALClose(poDataset);
	GDALClose(poDataset2);
	printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
	return 0;
}

int main(int argc, char **argv)
{
	if(argc != 4)
		return 0;
	if(!strcmp(argv[1], "0"))
		dcmodel(argv[2], argv[3]);
	else
		retinex(argv[2], argv[3]);
}
