#pragma once

typedef unsigned char BYTE;

class CDarkChannel
{
public:
	CDarkChannel(void);
	~CDarkChannel(void);
	BYTE* OnDarkChannel(BYTE*data,int height,int width,int nbands);
private:
	float* CalDarkChannel(float*data,int height,int width,int nbands,int size);
	float* GetLight(float* channel,float* data,int height,int width,int nbands);
	float* CalTransmission(float* channel, float* light,int height,int width,int nbands);
	BYTE* CalRadiance(float* Trans,float *data,float *light,int height,int width,int nbands);

	float* GuidedFilter(float*Trans,float*data,int height,int width,int nbands);
};
