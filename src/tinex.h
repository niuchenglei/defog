#pragma once
//////////////////////////
/*RetinexÄ£ÐÍÈ¥Îí */
////////////////////////////
typedef unsigned char BYTE;

float *CalGuassian(int height,int width,int tao);//ŒÆËãžßË¹º¯Êý
BYTE* OnHIS(BYTE* data,int height,int width,int type);

float* OnRetinex(float* data, int height,int width,int bandcount);//º¯ÊýœÓ¿Ú
float* OnMSR(float*pic,int height,int width);//¶à³ß¶ÈRetinexº¯ÊýœÓ¿Ú