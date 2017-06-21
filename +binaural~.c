#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include "m_pd.h"
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif
#include "impulsesc.h"
#include "impulsesd.h"
#include "impulsekemar.h"
#define max(a,b) (a>b?a:b)
#define BINAURAL_BLOCK_SIZE 128
#define	ELEVATED	0
#define	EAR_LEVEL	1
#define	LOWERED		2

#define TRUE 1.0f
#define FALSE 0.0f


/* ------------------------ binaural~ ----------------------------- */

static t_class *binaural_class;

// some convenient constants
enum
{
	kLevelIncrements = 96,
	kSizeFFT = 1024,
	kHalfSizeFFT = 512
};


typedef struct 
{
	float azimuth;
	float elevation;
	float *nearEar;
	float *farEar;
}hrtfPosition;

typedef struct 
{
	long numPositions;
	float gain;
	hrtfPosition *position;
}hrtf;



typedef struct _binaural
{
	
	t_object x_obj;
	t_float x_f; //current sample
	t_float azimuth;
	t_float gain;

	t_float sampleRate;
	t_int   bufferPosition;
	
	long filterSet;
	float currentAzimuth;
	int lastFilterFilled;
	long height;
	long filterSwitch;
	float mix, mixIncrement;
	float oneOverSizeImpulse;
	float azimuthOne, azimuthTwo, lastAzimuth;
	float *inBuffer, *outBufferL, *outBufferR;
	float *inputDouble;
	float *outputRightA, *outputLeftA;
	float *outputRightB, *outputLeftB;
	float *impulseLeftA, *impulseRightA;
	float *impulseLeftB, *impulseRightB;
	float *overlapLeftA, *overlapRightA;
	float *overlapLeftB, *overlapRightB;
	hrtf filterOne, filterTwo;
	float impulseScale;
	long halfSizeFFT, sizeFFT;
	long sizeImpulse;
	long sampleNumber;

	float *inSpectra, *impulseLeftASpectra, *impulseRightASpectra, *impulseLeftBSpectra, *impulseRightBSpectra,
			*outLeftSpectra, *outRightSpectra;
	hrtf *currentFilter;
} t_binaural;


//method declarations
void binaural_initHRTF(t_binaural *x);
void binaural_setFilterSet(t_binaural *x);
void binaural_ProcessMotion(t_binaural *x);
void binaural_Process(t_binaural *x);
short binaural_newImpulse(t_binaural *x);
void binaural_findImpulse(t_binaural *x, float *pImpulseL, float *pImpulseR, float pAzimuth);
void binaural_copyImpulse(float source[], float target[], long delay);
void binaural_mixImpulse(float sourceA[], float sourceB[], float target[], float percentB);


static void binaural_azimuth(t_binaural *x, t_float value)
{
	if(value<-180.0||value>180.0)
	{
		error("angle value must be between -180 and 180");
		return;
	}
	else
	{
		x->azimuth = value;
	}
}

static void binaural_gain(t_binaural *x, t_float value)
{	
	
	x->gain = powf(10.f, (value * 0.05f));
}

static void binaural_filterSet(t_binaural *x, t_float value)
{
	if (value < .5){
		x->filterSet = 0;
	}else{
		x->filterSet = 1;
	}
	binaural_setFilterSet(x);
}

void binaural_setFilterSet(t_binaural *x)
{
	x->filterSwitch = TRUE;
	switch(x->filterSet)
	{
		case 0:
			break;
		case 1:
			break;
	}
}

void binaural_ProcessMotion(t_binaural *x)
{
	long	i;
	long	n;
// this version of binauralProcess filters 2 positions and cross fades between them.
// we will cheat here - if we are going to cross a position, we would have to cross
// fade between 3 filters. too much processing though. instead, i will cross fade up
// to the position, then next block will start at the position and go a little faster
// hopefully this cheat will be clean.
		
	// a - get ready for a new position	
	binaural_newImpulse(x);

	// b - get new samples and FFT
	memset(x->inputDouble, 0, x->sizeFFT * sizeof(float));
	for(i = 0; i < x->sizeImpulse; i++)
		*(x->inputDouble + i) = *(x->inBuffer + i) * x->impulseScale * x->currentFilter->gain;

	mayer_realfft(x->sizeFFT, x->inputDouble);
	memcpy(x->inSpectra, x->inputDouble, x->sizeFFT * sizeof(float));

	// c - perform complex multiplication
	// left low angle
    x->outLeftSpectra[0] = x->impulseLeftASpectra[0] * x->inSpectra[0];	// DC Component
    x->outLeftSpectra[x->halfSizeFFT] = x->impulseLeftASpectra[x->halfSizeFFT] * x->inSpectra[x->halfSizeFFT];	// Nyquist Frequency
    for(i = 1; i < x->halfSizeFFT; ++i)
    {
    	x->outLeftSpectra[i] = x->impulseLeftASpectra[i] * x->inSpectra[i] - x->impulseLeftASpectra[x->sizeFFT - i] * x->inSpectra[x->sizeFFT - i];
    	x->outLeftSpectra[x->sizeFFT - i] = x->impulseLeftASpectra[i] * x->inSpectra[x->sizeFFT - i] + x->impulseLeftASpectra[x->sizeFFT - i] * x->inSpectra[i];
    }

	mayer_realifft(x->sizeFFT, x->outLeftSpectra);
	memcpy(x->outputLeftA, x->outLeftSpectra, x->sizeFFT * sizeof(float));
				
	// left high angle
    x->outLeftSpectra[0] = x->impulseLeftBSpectra[0] * x->inSpectra[0];	// DC Component
    x->outLeftSpectra[x->halfSizeFFT] = x->impulseLeftBSpectra[x->halfSizeFFT] * x->inSpectra[x->halfSizeFFT];	// Nyquist Frequency
    for(i = 1; i < x->halfSizeFFT; ++i)
    {
    	x->outLeftSpectra[i] = x->impulseLeftBSpectra[i] * x->inSpectra[i] - x->impulseLeftBSpectra[x->sizeFFT - i] * x->inSpectra[x->sizeFFT - i];
    	x->outLeftSpectra[x->sizeFFT - i] = x->impulseLeftBSpectra[i] * x->inSpectra[x->sizeFFT - i] + x->impulseLeftBSpectra[x->sizeFFT - i] * x->inSpectra[i];
    }

	mayer_realifft(x->sizeFFT, x->outLeftSpectra);
	memcpy(x->outputLeftB, x->outLeftSpectra, x->sizeFFT * sizeof(float));
	
	// right low angle
    x->outRightSpectra[0] = x->impulseRightASpectra[0] * x->inSpectra[0];	// DC Component
    x->outRightSpectra[x->halfSizeFFT] = x->impulseRightASpectra[x->halfSizeFFT] * x->inSpectra[x->halfSizeFFT];	// Nyquist Frequency
    for(i = 1; i < x->halfSizeFFT; ++i)
    {
    	x->outRightSpectra[i] = x->impulseRightASpectra[i] * x->inSpectra[i] - x->impulseRightASpectra[x->sizeFFT - i] * x->inSpectra[x->sizeFFT - i];
    	x->outRightSpectra[x->sizeFFT - i] = x->impulseRightASpectra[i] * x->inSpectra[x->sizeFFT - i] + x->impulseRightASpectra[x->sizeFFT - i] * x->inSpectra[i];
    }

	mayer_realifft(x->sizeFFT, x->outRightSpectra);
	memcpy(x->outputRightA, x->outRightSpectra, x->sizeFFT * sizeof(float));
	
	// right high angle
    x->outRightSpectra[0] = x->impulseRightBSpectra[0] * x->inSpectra[0];	// DC Component
    x->outRightSpectra[x->halfSizeFFT] = x->impulseRightBSpectra[x->halfSizeFFT] * x->inSpectra[x->halfSizeFFT];	// Nyquist Frequency
    for(i = 1; i < x->halfSizeFFT; ++i)
    {
    	x->outRightSpectra[i] = x->impulseRightBSpectra[i] * x->inSpectra[i] - x->impulseRightBSpectra[x->sizeFFT - i] * x->inSpectra[x->sizeFFT - i];
    	x->outRightSpectra[x->sizeFFT - i] = x->impulseRightBSpectra[i] * x->inSpectra[x->sizeFFT - i] + x->impulseRightBSpectra[x->sizeFFT - i] * x->inSpectra[i];
    }

	mayer_realifft(x->sizeFFT, x->outRightSpectra);
	memcpy(x->outputRightB, x->outRightSpectra, x->sizeFFT * sizeof(float));

	for(n = 0; n < (x->sizeImpulse - 1); n++)
	{
		x->outputLeftA[n] = x->outputLeftA[n] + x->overlapLeftA[n];
		x->overlapLeftA[n] = x->outputLeftA[x->sizeImpulse + n];
		x->outputRightA[n] = x->outputRightA[n] + x->overlapRightA[n];
		x->overlapRightA[n] = x->outputRightA[x->sizeImpulse + n];
		x->outputLeftB[n] = x->outputLeftB[n] + x->overlapLeftB[n];
		x->overlapLeftB[n] = x->outputLeftB[x->sizeImpulse + n];
		x->outputRightB[n] = x->outputRightB[n] + x->overlapRightB[n];
		x->overlapRightB[n] = x->outputRightB[x->sizeImpulse + n];
	}
	// e - write stuff out	
	for(i = 0; i < x->sizeImpulse; i++)
	{
		*(x->outBufferL+i) = (*(x->outputLeftA+i) + ((*(x->outputLeftB+i) - *(x->outputLeftA+i)) * x->mix)) * x->gain;
		*(x->outBufferR+i) = (*(x->outputRightA+i) + ((*(x->outputRightB+i) - *(x->outputRightA+i)) * x->mix)) * x->gain;
		x->mix += x->mixIncrement;
	}
}


short binaural_newImpulse(t_binaural *x)
{
	long n;
	float rotatedAzimuth;
	
	switch(x->filterSet)
	{
		case 0:
			x->currentFilter = &x->filterOne;
			break;
		case 1:
			x->currentFilter = &x->filterTwo;
			break;
	}
		
	// go from -180 - 0 - 180 to 0 - 360
	if(x->azimuth < 0)
		rotatedAzimuth = x->azimuth + 360.0f;
	else
		rotatedAzimuth = x->azimuth;
		
//	if(filterSwitch == false)
//		return(false);
//	filterSwitch = false;

	// limit the speed of change
	if(x->currentAzimuth > rotatedAzimuth)
	{
		if((x->currentAzimuth - rotatedAzimuth) >= 180.0f)
			x->currentAzimuth += 1.0f;
		else
			x->currentAzimuth -= 1.0f;
		if(x->currentAzimuth <= rotatedAzimuth)
			x->currentAzimuth = rotatedAzimuth;
	}
	else
	{
		if((rotatedAzimuth - x->currentAzimuth) >= 180.0f)
			x->currentAzimuth -= 1.0f;
		else
			x->currentAzimuth += 1.0f;
		if(x->currentAzimuth >= rotatedAzimuth)
			x->currentAzimuth = rotatedAzimuth;
	}
	// without limits
//	currentAzimuth = rotatedAzimuth;
	
	while(x->currentAzimuth >= 360.0f)
		x->currentAzimuth -= 360.0f;
	while(x->currentAzimuth < 0.0f)
		x->currentAzimuth += 360.0f;

	if(x->lastFilterFilled == 0)
	{
		binaural_findImpulse(x, x->impulseLeftB, x->impulseRightB, x->currentAzimuth);
		for(n = x->sizeImpulse; n<x->sizeFFT; n++)
		{
			x->impulseLeftB[n] = 0.0;
			x->impulseRightB[n] = 0.0;
		}

		mayer_realfft(x->sizeFFT, x->impulseLeftB);
		memcpy(x->impulseLeftBSpectra, x->impulseLeftB, x->sizeFFT * sizeof(float));
		
		mayer_realfft(x->sizeFFT, x->impulseRightB);
		memcpy(x->impulseRightBSpectra, x->impulseRightB, x->sizeFFT * sizeof(float));

		x->mix = 0.0f;
		x->mixIncrement = x->oneOverSizeImpulse;
		x->lastFilterFilled = 1;
	}
	else
	{
		binaural_findImpulse(x, x->impulseLeftA, x->impulseRightA, x->currentAzimuth);
		for(n = x->sizeImpulse; n<x->sizeFFT; n++)
		{
			x->impulseLeftA[n] = 0.0;
			x->impulseRightA[n] = 0.0;
		}
		mayer_realfft(x->sizeFFT, x->impulseLeftA);
		memcpy(x->impulseLeftASpectra, x->impulseLeftA, x->sizeFFT * sizeof(float));
		
		mayer_realfft(x->sizeFFT, x->impulseRightA);
		memcpy(x->impulseRightASpectra, x->impulseRightA, x->sizeFFT * sizeof(float));

		x->mix = 1.0f;
		x->mixIncrement = -x->oneOverSizeImpulse;
		x->lastFilterFilled = 0;
	}
	return(TRUE);
}



void binaural_findImpulse(t_binaural *x, float *pImpulseL, float *pImpulseR, float pAzimuth)
{
	int i, j;
	float ratio;
	float *filterOneL, *filterTwoL, *filterOneR, *filterTwoR;

	if(pAzimuth < 180.0f)
	{
		filterOneR = filterTwoR = x->currentFilter->position[0].nearEar;
		filterOneL = filterTwoL = x->currentFilter->position[0].farEar;
		for(i = 0; i < x->currentFilter->numPositions; i++)
		{
			j = i+1;
			// this will be true if j is 37
			// that is, if currentazimuth is not between i = 35 and j = 36 (180)
			// on the last test
			// actually - this should never occur, since we are testing for CA under 180
			if(j >= x->currentFilter->numPositions)
			{
				// position zero for j is always 180
				j = i;
				if(x->currentFilter->position[i].azimuth <= pAzimuth)
				{
					filterOneR = x->currentFilter->position[i].nearEar;
					filterOneL = x->currentFilter->position[i].farEar;
					x->azimuthOne = x->currentFilter->position[i].azimuth;
					filterTwoR = x->currentFilter->position[j].nearEar;
					filterTwoL = x->currentFilter->position[j].farEar;
					x->azimuthTwo = x->currentFilter->position[j].azimuth;
					ratio = (pAzimuth - x->azimuthOne) / (180.f - x->azimuthOne);
					// this causes it to break out
					i = x->currentFilter->numPositions;
				}
			}
			else
			{
				if((x->currentFilter->position[i].azimuth <= pAzimuth) && (x->currentFilter->position[j].azimuth >= pAzimuth))
				{
					filterOneR = x->currentFilter->position[i].nearEar;
					filterOneL = x->currentFilter->position[i].farEar;
					x->azimuthOne = x->currentFilter->position[i].azimuth;
					filterTwoR = x->currentFilter->position[j].nearEar;
					filterTwoL = x->currentFilter->position[j].farEar;
					x->azimuthTwo = x->currentFilter->position[j].azimuth;
					ratio = (pAzimuth - x->azimuthOne) / (x->azimuthTwo - x->azimuthOne);
					// this causes it to break out
					i = x->currentFilter->numPositions;
				}
			}
		}
	}
	// we are mirroring zimuth above 180 by using 360 - azimut, and L = near, R = far
	else
	{
		filterOneL = filterTwoL = x->currentFilter->position[0].nearEar;
		filterOneR = filterTwoR = x->currentFilter->position[0].farEar;
		for(i = 0; i < x->currentFilter->numPositions; i++)
		{
			j = i+ 1;
			if(j >= x->currentFilter->numPositions)
			{
				// position zero for j is always 180
				j = 0;
				if(x->currentFilter->position[i].azimuth <= (360.0f - pAzimuth))
				{
					filterOneL = x->currentFilter->position[i].nearEar;
					filterOneR = x->currentFilter->position[i].farEar;
					x->azimuthOne = x->currentFilter->position[i].azimuth;
					filterTwoL = x->currentFilter->position[j].nearEar;
					filterTwoR = x->currentFilter->position[j].farEar;
					x->azimuthTwo = x->currentFilter->position[j].azimuth;
					ratio = ((360.0f - pAzimuth) - x->azimuthOne) / (180.f - x->azimuthOne);
					// this causes it to break out
					i = x->currentFilter->numPositions;
				}
			}
			else
			{
				if((x->currentFilter->position[i].azimuth <= (360.0f - pAzimuth)) 
					&& (x->currentFilter->position[j].azimuth >= (360.0f - pAzimuth)))
				{
					// we have a match for filter one and two
					filterOneL = x->currentFilter->position[i].nearEar;
					filterOneR = x->currentFilter->position[i].farEar;
					x->azimuthOne = x->currentFilter->position[i].azimuth;
					filterTwoL = x->currentFilter->position[j].nearEar;
					filterTwoR = x->currentFilter->position[j].farEar;
					x->azimuthTwo = x->currentFilter->position[j].azimuth;
					ratio = ((360.0f - pAzimuth) - x->azimuthOne) / (x->azimuthTwo - x->azimuthOne);
					// la 1 ca 359 - la 361
					i = x->currentFilter->numPositions;
				}
			}
		}
	}
	
	binaural_mixImpulse(filterOneL, filterTwoL, pImpulseL, ratio);
	binaural_mixImpulse(filterOneR, filterTwoR, pImpulseR, ratio);
}

void binaural_copyImpulse(float source[], float target[], long delay)
{
	long m, n;
	
	for(n = 0; n < delay; n++)
		target[n] = 0.0;
	for(m = 0; n < BINAURAL_BLOCK_SIZE; n++, m++)
		target[n] = source[m];
}

void binaural_mixImpulse(float sourceA[], float sourceB[], float target[], float percentB)
{
	long ia, ib, n;

	for(n = 0, ia = 0, ib = 0; n < 256; ia++, ib++, n++)
		target[n] = ((sourceB[ib] - sourceA[ia]) * percentB) + sourceA[ia];
}


void binaural_initHRTF(t_binaural *x)
{
	 
	// read in HRTF KEMAR
	x->filterOne.numPositions = 37;
	x->filterOne.gain = 0.2630891791f;
	//x->filterOne.position = (hrtfPosition *)malloc(sizeof(hrtfPosition*));
	x->filterOne.position = malloc (sizeof (hrtfPosition)*x->filterOne.numPositions);
	// 0 - 1
	x->filterOne.position[0].azimuth = 0.0f;
	x->filterOne.position[0].elevation = 0.0f;
	x->filterOne.position[0].nearEar = H0e000aright;
	x->filterOne.position[0].farEar = H0e000aleft;
	// 5 - 2
	x->filterOne.position[1].azimuth = 5.0f;
	x->filterOne.position[1].elevation = 0.0f;
	x->filterOne.position[1].nearEar = H0e005aright;
	x->filterOne.position[1].farEar = H0e005aleft;
	// 10 - 3
	x->filterOne.position[2].azimuth = 10.0f;
	x->filterOne.position[2].elevation = 0.0f;
	x->filterOne.position[2].nearEar = H0e010aright;
	x->filterOne.position[2].farEar = H0e010aleft;
	// 15 - 4
	x->filterOne.position[3].azimuth = 15.0f;
	x->filterOne.position[3].elevation = 0.0f;
	x->filterOne.position[3].nearEar = H0e015aright;
	x->filterOne.position[3].farEar = H0e015aleft;
	// 20 - 5
	x->filterOne.position[4].azimuth = 20.0f;
	x->filterOne.position[4].elevation = 0.0f;
	x->filterOne.position[4].nearEar = H0e020aright;
	x->filterOne.position[4].farEar = H0e020aleft;
	// 25 - 6
	x->filterOne.position[5].azimuth = 25.0f;
	x->filterOne.position[5].elevation = 0.0f;
	x->filterOne.position[5].nearEar = H0e025aright;
	x->filterOne.position[5].farEar = H0e025aleft;
	// 30 - 7
	x->filterOne.position[6].azimuth = 30.0f;
	x->filterOne.position[6].elevation = 0.0f;
	x->filterOne.position[6].nearEar = H0e030aright;
	x->filterOne.position[6].farEar = H0e030aleft;
	// 35 - 8
	x->filterOne.position[7].azimuth = 35.0f;
	x->filterOne.position[7].elevation = 0.0f;
	x->filterOne.position[7].nearEar = H0e035aright;
	x->filterOne.position[7].farEar = H0e035aleft;
	// 40 - 9
	x->filterOne.position[8].azimuth = 40.0f;
	x->filterOne.position[8].elevation = 0.0f;
	x->filterOne.position[8].nearEar = H0e040aright;
	x->filterOne.position[8].farEar = H0e040aleft;
	// 45 - 10
	x->filterOne.position[9].azimuth = 45.0f;
	x->filterOne.position[9].elevation = 0.0f;
	x->filterOne.position[9].nearEar = H0e045aright;
	x->filterOne.position[9].farEar = H0e045aleft;
	// 50 - 11
	x->filterOne.position[10].azimuth = 50.0f;
	x->filterOne.position[10].elevation = 0.0f;
	x->filterOne.position[10].nearEar = H0e050aright;
	x->filterOne.position[10].farEar = H0e050aleft;
	// 55 - 12
	x->filterOne.position[11].azimuth = 55.0f;
	x->filterOne.position[11].elevation = 0.0f;
	x->filterOne.position[11].nearEar = H0e055aright;
	x->filterOne.position[11].farEar = H0e055aleft;
	// 60 - 13
	x->filterOne.position[12].azimuth = 60.0f;
	x->filterOne.position[12].elevation = 0.0f;
	x->filterOne.position[12].nearEar = H0e060aright;
	x->filterOne.position[12].farEar = H0e060aleft;
	// 65 - 14
	x->filterOne.position[13].azimuth = 65.0f;
	x->filterOne.position[13].elevation = 0.0f;
	x->filterOne.position[13].nearEar = H0e065aright;
	x->filterOne.position[13].farEar = H0e065aleft;
	// 70 - 15
	x->filterOne.position[14].azimuth = 70.0f;
	x->filterOne.position[14].elevation = 0.0f;
	x->filterOne.position[14].nearEar = H0e070aright;
	x->filterOne.position[14].farEar = H0e070aleft;
	// 75 - 16
	x->filterOne.position[15].azimuth = 75.0f;
	x->filterOne.position[15].elevation = 0.0f;
	x->filterOne.position[15].nearEar = H0e075aright;
	x->filterOne.position[15].farEar = H0e075aleft;
	// 80 - 17
	x->filterOne.position[16].azimuth = 80.0f;
	x->filterOne.position[16].elevation = 0.0f;
	x->filterOne.position[16].nearEar = H0e080aright;
	x->filterOne.position[16].farEar = H0e080aleft;
	// 85 - 18
	x->filterOne.position[17].azimuth = 85.0f;
	x->filterOne.position[17].elevation = 0.0f;
	x->filterOne.position[17].nearEar = H0e085aright;
	x->filterOne.position[17].farEar = H0e085aleft;
	// 90 - 19
	x->filterOne.position[18].azimuth = 90.0f;
	x->filterOne.position[18].elevation = 0.0f;
	x->filterOne.position[18].nearEar = H0e090aright;
	x->filterOne.position[18].farEar = H0e090aleft;
	// 95 - 20
	x->filterOne.position[19].azimuth = 95.0f;
	x->filterOne.position[19].elevation = 0.0f;
	x->filterOne.position[19].nearEar = H0e095aright;
	x->filterOne.position[19].farEar = H0e095aleft;
	// 100 - 21
	x->filterOne.position[20].azimuth = 100.0f;
	x->filterOne.position[20].elevation = 0.0f;
	x->filterOne.position[20].nearEar = H0e100aright;
	x->filterOne.position[20].farEar = H0e100aleft;
	// 105 - 22
	x->filterOne.position[21].azimuth = 105.0f;
	x->filterOne.position[21].elevation = 0.0f;
	x->filterOne.position[21].nearEar = H0e105aright;
	x->filterOne.position[21].farEar = H0e105aleft;
	// 110 - 23
	x->filterOne.position[22].azimuth = 110.0f;
	x->filterOne.position[22].elevation = 0.0f;
	x->filterOne.position[22].nearEar = H0e110aright;
	x->filterOne.position[22].farEar = H0e110aleft;
	// 115 - 24
	x->filterOne.position[23].azimuth = 115.0f;
	x->filterOne.position[23].elevation = 0.0f;
	x->filterOne.position[23].nearEar = H0e115aright;
	x->filterOne.position[23].farEar = H0e115aleft;
	// 120 - 25
	x->filterOne.position[24].azimuth = 120.0f;
	x->filterOne.position[24].elevation = 0.0f;
	x->filterOne.position[24].nearEar = H0e120aright;
	x->filterOne.position[24].farEar = H0e120aleft;
	// 125 - 26
	x->filterOne.position[25].azimuth = 125.0f;
	x->filterOne.position[25].elevation = 0.0f;
	x->filterOne.position[25].nearEar = H0e125aright;
	x->filterOne.position[25].farEar = H0e125aleft;
	// 130 - 27
	x->filterOne.position[26].azimuth = 130.0f;
	x->filterOne.position[26].elevation = 0.0f;
	x->filterOne.position[26].nearEar = H0e130aright;
	x->filterOne.position[26].farEar = H0e130aleft;
	// 135 - 28
	x->filterOne.position[27].azimuth = 135.0f;
	x->filterOne.position[27].elevation = 0.0f;
	x->filterOne.position[27].nearEar = H0e135aright;
	x->filterOne.position[27].farEar = H0e135aleft;
	// 140 - 29
	x->filterOne.position[28].azimuth = 140.0f;
	x->filterOne.position[28].elevation = 0.0f;
	x->filterOne.position[28].nearEar = H0e140aright;
	x->filterOne.position[28].farEar = H0e140aleft;
	// 145 - 30
	x->filterOne.position[29].azimuth = 145.0f;
	x->filterOne.position[29].elevation = 0.0f;
	x->filterOne.position[29].nearEar = H0e145aright;
	x->filterOne.position[29].farEar = H0e145aleft;
	// 150 - 31
	x->filterOne.position[30].azimuth = 150.0f;
	x->filterOne.position[30].elevation = 0.0f;
	x->filterOne.position[30].nearEar = H0e150aright;
	x->filterOne.position[30].farEar = H0e150aleft;
	// 155 - 32
	x->filterOne.position[31].azimuth = 155.0f;
	x->filterOne.position[31].elevation = 0.0f;
	x->filterOne.position[31].nearEar = H0e155aright;
	x->filterOne.position[31].farEar = H0e155aleft;
	// 160 - 33
	x->filterOne.position[32].azimuth = 160.0f;
	x->filterOne.position[32].elevation = 0.0f;
	x->filterOne.position[32].nearEar = H0e160aright;
	x->filterOne.position[32].farEar = H0e160aleft;
	// 165 - 34
	x->filterOne.position[33].azimuth = 165.0f;
	x->filterOne.position[33].elevation = 0.0f;
	x->filterOne.position[33].nearEar = H0e165aright;
	x->filterOne.position[33].farEar = H0e165aleft;
	// 170 - 35
	x->filterOne.position[34].azimuth = 170.0f;
	x->filterOne.position[34].elevation = 0.0f;
	x->filterOne.position[34].nearEar = H0e170aright;
	x->filterOne.position[34].farEar = H0e170aleft;
	// 175 - 36
	x->filterOne.position[35].azimuth = 175.0f;
	x->filterOne.position[35].elevation = 0.0f;
	x->filterOne.position[35].nearEar = H0e175aright;
	x->filterOne.position[35].farEar = H0e175aleft;
	// 180 - 37
	x->filterOne.position[36].azimuth = 180.0f;
	x->filterOne.position[36].elevation = 0.0f;
	x->filterOne.position[36].nearEar = H0e180aright;
	x->filterOne.position[36].farEar = H0e180aleft;
/*	
	// read in HRTF D
	filterOne.numPositions = 36;
	filterOne.gain = 0.4510151663f;
	filterOne.position = new hrtfPosition[filterOne.numPositions];
	// 0 - 1
	filterOne.position[0].azimuth = 0.0f;
	filterOne.position[0].elevation = 0.0f;
	filterOne.position[0].nearEar = d000l;
	filterOne.position[0].farEar = d000r;
	// 5 - 2
	filterOne.position[1].azimuth = 5.0f;
	filterOne.position[1].elevation = 0.0f;
	filterOne.position[1].nearEar = d005l;
	filterOne.position[1].farEar = d005r;
	// 10 - 3
	filterOne.position[2].azimuth = 10.0f;
	filterOne.position[2].elevation = 0.0f;
	filterOne.position[2].nearEar = d010l;
	filterOne.position[2].farEar = d010r;
	// 15 - 4
	filterOne.position[3].azimuth = 15.0f;
	filterOne.position[3].elevation = 0.0f;
	filterOne.position[3].nearEar = d015l;
	filterOne.position[3].farEar = d015r;
	// 20 - 5
	filterOne.position[4].azimuth = 20.0f;
	filterOne.position[4].elevation = 0.0f;
	filterOne.position[4].nearEar = d020l;
	filterOne.position[4].farEar = d020r;
	// 25 - 6
	filterOne.position[5].azimuth = 25.0f;
	filterOne.position[5].elevation = 0.0f;
	filterOne.position[5].nearEar = d025l;
	filterOne.position[5].farEar = d025r;
	// 30 - 7
	filterOne.position[6].azimuth = 30.0f;
	filterOne.position[6].elevation = 0.0f;
	filterOne.position[6].nearEar = d030l;
	filterOne.position[6].farEar = d030r;
	// 35 - 8
	filterOne.position[7].azimuth = 35.0f;
	filterOne.position[7].elevation = 0.0f;
	filterOne.position[7].nearEar = d035l;
	filterOne.position[7].farEar = d035r;
	// 40 - 9
	filterOne.position[8].azimuth = 40.0f;
	filterOne.position[8].elevation = 0.0f;
	filterOne.position[8].nearEar = d040l;
	filterOne.position[8].farEar = d040r;
	// 45 - 10
	filterOne.position[9].azimuth = 45.0f;
	filterOne.position[9].elevation = 0.0f;
	filterOne.position[9].nearEar = d045l;
	filterOne.position[9].farEar = d045r;
	// 50 - 11
	filterOne.position[10].azimuth = 50.0f;
	filterOne.position[10].elevation = 0.0f;
	filterOne.position[10].nearEar = d050l;
	filterOne.position[10].farEar = d050r;
	// 55 - 12
	filterOne.position[11].azimuth = 55.0f;
	filterOne.position[11].elevation = 0.0f;
	filterOne.position[11].nearEar = d055l;
	filterOne.position[11].farEar = d055r;
	// 60 - 13
	filterOne.position[12].azimuth = 60.0f;
	filterOne.position[12].elevation = 0.0f;
	filterOne.position[12].nearEar = d060l;
	filterOne.position[12].farEar = d060r;
	// 65 - 14
	filterOne.position[13].azimuth = 65.0f;
	filterOne.position[13].elevation = 0.0f;
	filterOne.position[13].nearEar = d065l;
	filterOne.position[13].farEar = d065r;
	// 70 - 15
	filterOne.position[14].azimuth = 70.0f;
	filterOne.position[14].elevation = 0.0f;
	filterOne.position[14].nearEar = d070l;
	filterOne.position[14].farEar = d070r;
	// 75 - 16
	filterOne.position[15].azimuth = 75.0f;
	filterOne.position[15].elevation = 0.0f;
	filterOne.position[15].nearEar = d075l;
	filterOne.position[15].farEar = d075r;
	// 80 - 17
	filterOne.position[16].azimuth = 80.0f;
	filterOne.position[16].elevation = 0.0f;
	filterOne.position[16].nearEar = d080l;
	filterOne.position[16].farEar = d080r;
	// 85 - 18
	filterOne.position[17].azimuth = 85.0f;
	filterOne.position[17].elevation = 0.0f;
	filterOne.position[17].nearEar = d085l;
	filterOne.position[17].farEar = d085r;
	// 90 - 19
	filterOne.position[18].azimuth = 90.0f;
	filterOne.position[18].elevation = 0.0f;
	filterOne.position[18].nearEar = d090l;
	filterOne.position[18].farEar = d090r;
	// 95 - 20
	filterOne.position[19].azimuth = 95.0f;
	filterOne.position[19].elevation = 0.0f;
	filterOne.position[19].nearEar = d095l;
	filterOne.position[19].farEar = d095r;
	// 100 - 21
	filterOne.position[20].azimuth = 100.0f;
	filterOne.position[20].elevation = 0.0f;
	filterOne.position[20].nearEar = d100l;
	filterOne.position[20].farEar = d100r;
	// 105 - 22
	filterOne.position[21].azimuth = 105.0f;
	filterOne.position[21].elevation = 0.0f;
	filterOne.position[21].nearEar = d105l;
	filterOne.position[21].farEar = d105r;
	// 110 - 23
	filterOne.position[22].azimuth = 110.0f;
	filterOne.position[22].elevation = 0.0f;
	filterOne.position[22].nearEar = d110l;
	filterOne.position[22].farEar = d110r;
	// 115 - 24
	filterOne.position[23].azimuth = 115.0f;
	filterOne.position[23].elevation = 0.0f;
	filterOne.position[23].nearEar = d115l;
	filterOne.position[23].farEar = d115r;
	// 120 - 25
	filterOne.position[24].azimuth = 120.0f;
	filterOne.position[24].elevation = 0.0f;
	filterOne.position[24].nearEar = d120l;
	filterOne.position[24].farEar = d120r;
	// 125 - 26
	filterOne.position[25].azimuth = 125.0f;
	filterOne.position[25].elevation = 0.0f;
	filterOne.position[25].nearEar = d125l;
	filterOne.position[25].farEar = d125r;
	// 130 - 27
	filterOne.position[26].azimuth = 130.0f;
	filterOne.position[26].elevation = 0.0f;
	filterOne.position[26].nearEar = d130l;
	filterOne.position[26].farEar = d130r;
	// 135 - 28
	filterOne.position[27].azimuth = 135.0f;
	filterOne.position[27].elevation = 0.0f;
	filterOne.position[27].nearEar = d135l;
	filterOne.position[27].farEar = d135r;
	// 140 - 29
	filterOne.position[28].azimuth = 140.0f;
	filterOne.position[28].elevation = 0.0f;
	filterOne.position[28].nearEar = d140l;
	filterOne.position[28].farEar = d140r;
	// 145 - 30
	filterOne.position[29].azimuth = 145.0f;
	filterOne.position[29].elevation = 0.0f;
	filterOne.position[29].nearEar = d145l;
	filterOne.position[29].farEar = d145r;
	// 155 - 31
	filterOne.position[30].azimuth = 155.0f;
	filterOne.position[30].elevation = 0.0f;
	filterOne.position[30].nearEar = d155l;
	filterOne.position[30].farEar = d155r;
	// 160 - 32
	filterOne.position[31].azimuth = 160.0f;
	filterOne.position[31].elevation = 0.0f;
	filterOne.position[31].nearEar = d160l;
	filterOne.position[31].farEar = d160r;
	// 165 - 33
	filterOne.position[32].azimuth = 165.0f;
	filterOne.position[32].elevation = 0.0f;
	filterOne.position[32].nearEar = d165l;
	filterOne.position[32].farEar = d165r;
	// 170 - 34
	filterOne.position[33].azimuth = 170.0f;
	filterOne.position[33].elevation = 0.0f;
	filterOne.position[33].nearEar = d170l;
	filterOne.position[33].farEar = d170r;
	// 175 - 35
	filterOne.position[34].azimuth = 175.0f;
	filterOne.position[34].elevation = 0.0f;
	filterOne.position[34].nearEar = d175l;
	filterOne.position[34].farEar = d175r;
	// 180 - 36
	filterOne.position[35].azimuth = 180.0f;
	filterOne.position[35].elevation = 0.0f;
	filterOne.position[35].nearEar = d180l;
	filterOne.position[35].farEar = d180r;
*/	
	// read in HRTF C
	x->filterTwo.numPositions = 7;
	x->filterTwo.gain = 0.5405306362f;
	x->filterTwo.position = malloc (sizeof (hrtfPosition)*x->filterTwo.numPositions);
	// 0 - 1
	x->filterTwo.position[0].azimuth = 0.0f;
	x->filterTwo.position[0].elevation = 0.0f;
	x->filterTwo.position[0].nearEar = a000l;
	x->filterTwo.position[0].farEar = a360l;
	// 30 - 2
	x->filterTwo.position[1].azimuth = 30.0f;
	x->filterTwo.position[1].elevation = 0.0f;
	x->filterTwo.position[1].nearEar = a030l;
	x->filterTwo.position[1].farEar = a330l;
	// 60 - 3
	x->filterTwo.position[2].azimuth = 60.0f;
	x->filterTwo.position[2].elevation = 0.0f;
	x->filterTwo.position[2].nearEar = a060l;
	x->filterTwo.position[2].farEar = a300l;
	// 90 - 4
	x->filterTwo.position[3].azimuth = 90.0f;
	x->filterTwo.position[3].elevation = 0.0f;
	x->filterTwo.position[3].nearEar = a090l;
	x->filterTwo.position[3].farEar = a270l;
	// 120 - 5
	x->filterTwo.position[4].azimuth = 120.0f;
	x->filterTwo.position[4].elevation = 0.0f;
	x->filterTwo.position[4].nearEar = a120l;
	x->filterTwo.position[4].farEar = a240l;
	// 150 - 6
	x->filterTwo.position[5].azimuth = 150.0f;
	x->filterTwo.position[5].elevation = 0.0f;
	x->filterTwo.position[5].nearEar = a150l;
	x->filterTwo.position[5].farEar = a210l;
	// 180 - 7
	x->filterTwo.position[6].azimuth = 180.0f;
	x->filterTwo.position[6].elevation = 0.0f;
	x->filterTwo.position[6].nearEar = a180l;
	x->filterTwo.position[6].farEar = a181l;
}


static void *binaural_new(t_floatarg f)
{
    long sizeConvolution, i;
	float Pi, twoPi;
	
	t_binaural *x = (t_binaural *)pd_new(binaural_class);

    outlet_new(&x->x_obj, gensym("signal"));
	outlet_new(&x->x_obj, gensym("signal"));

	x->mix = 0;
	x->mixIncrement = 0;
	x->outBufferL = x->outBufferR = x->inBuffer = 0;
	x->impulseLeftA = x->impulseRightA = x->impulseLeftB = x->impulseRightB = 0;
	x->outputRightA = x->outputLeftA = x->outputRightB = x->outputLeftB = x->inputDouble = 0;
	x->overlapLeftA = x->overlapRightA = x->overlapLeftB = x->overlapRightB = 0;
	x->lastFilterFilled = 1;

	Pi = 4.0f * atanf(1.0f);
	twoPi = 8.0f * atanf(1.0f);
	x->sizeImpulse = BINAURAL_BLOCK_SIZE;
	x->oneOverSizeImpulse = 1.0f/(float)x->sizeImpulse;
	sizeConvolution = 2 * x->sizeImpulse - 1;
	for(x->sizeFFT = 1; x->sizeFFT < sizeConvolution;)
		x->sizeFFT <<= 1;
	x->halfSizeFFT = x->sizeFFT >> 1;
	x->impulseScale = (2.0f)/(x->sizeFFT);
	x->gain = 1.0f;
	x->bufferPosition = 0;
	x->currentAzimuth = 0.0f;

	x->inBuffer = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outBufferL = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outBufferR = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseLeftA = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseRightA = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseLeftB = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseRightB = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->inputDouble = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outputLeftA = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outputRightA = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outputLeftB = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outputRightB = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->overlapLeftA = (float *) malloc(sizeof(float) * (x->sizeFFT + 2));
	x->overlapRightA = (float *) malloc(sizeof(float) * (x->sizeFFT + 2));
	x->overlapLeftB = (float *) malloc(sizeof(float) * (x->sizeFFT + 2));
	x->overlapRightB = (float *) malloc(sizeof(float) * (x->sizeFFT + 2));
	x->inSpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseLeftASpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseRightASpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseLeftBSpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->impulseRightBSpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outLeftSpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outRightSpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
		
	x->height = EAR_LEVEL;
	/*if(allocMem() == false)
		return(false);*/
	binaural_initHRTF(x);
	for(i = 0; i <(x->sizeFFT+2); i++)
		x->overlapRightA[i] = x->overlapLeftA[i] = x->overlapRightB[i] = x->overlapLeftB[i] = 0.0;
	for(i = 0; i<x->sizeFFT; i++)
	{
		x->inputDouble[i] = x->inBuffer[i] = x->outBufferL[i] = x->outBufferR[i] = 0.0;
	}
	x->currentFilter = &x->filterOne;

	binaural_azimuth(x, f);
	
    return (x);
}


static t_int *binaural_perform(t_int *w)
{

	t_binaural *x = (t_binaural *)(w[1]);
    t_float *in = (t_float *)(w[2]);
    t_float *outL = (t_float *)(w[3]);
	t_float *outR = (t_float *)(w[4]);
    int sampleframes = (int)(w[5]);


	int i, framesLeft, processframes;

	
	framesLeft = sampleframes;
	while (framesLeft > 0)
	{
		// how many frames can we process now
		// with this we insure that we stop on the 
		// BINAURAL_BLOCK_SIZE boundary
		if(framesLeft+x->bufferPosition < BINAURAL_BLOCK_SIZE)
			processframes = framesLeft;
		else
			processframes = BINAURAL_BLOCK_SIZE - x->bufferPosition;
		// copy in the new input..., flush out the previous output
		memcpy(x->inBuffer + x->bufferPosition, in, processframes*sizeof(float));
		for(i=0; i<processframes; i++)
		{
			outL[i] = x->outBufferL[i+x->bufferPosition];
			outR[i] = x->outBufferR[i+x->bufferPosition];
		}		
		x->bufferPosition += processframes;
		// if over the midpoint boundry, we process a new block
		if(x->bufferPosition == BINAURAL_BLOCK_SIZE)
		{
			x->bufferPosition = 0;
			binaural_ProcessMotion(x);
		}
		in += processframes;
		outL += processframes;
		outR += processframes;
		framesLeft -= processframes;
	}
	x->sampleNumber += sampleframes;




  //  while (n-- > 0)
  //  {
  //      *outL++ = *in++; 
		//*outR++ = *in++;
  //  }
    return (w+6);
}

static void binaural_dsp(t_binaural *x, t_signal **sp)
{
	x->sampleRate = sp[0]->s_sr;
	dsp_add(binaural_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
}


static void binaural_free(t_binaural *x) 
{
	if(x->inBuffer != 0) free(x->inBuffer);
	if(x->outBufferL != 0) free(x->outBufferL);
	if(x->outBufferR != 0) free(x->outBufferR);

	if(x->impulseLeftA != 0) free(x->impulseLeftA);
	if(x->impulseRightA != 0) free(x->impulseRightA);
	if(x->impulseLeftB != 0) free(x->impulseLeftB);
	if(x->impulseRightB != 0) free(x->impulseRightB);
	
	if(x->inputDouble != 0) free(x->inputDouble);
	
	if(x->outputLeftA != 0) free(x->outputLeftA);
	if(x->outputRightA != 0) free(x->outputRightA);
	if(x->outputLeftB != 0) free(x->outputLeftB);
	if(x->outputRightB != 0) free(x->outputRightB);
	
	if(x->overlapLeftA != 0) free(x->overlapLeftA);
	if(x->overlapRightA != 0) free(x->overlapRightA);
	if(x->overlapLeftB != 0) free(x->overlapLeftB);
	if(x->overlapRightB != 0) free(x->overlapRightB);
	
	if(x->inSpectra != 0) free(x->inSpectra);
	
	if(x->impulseLeftASpectra != 0) free(x->impulseLeftASpectra);
	if(x->impulseRightASpectra != 0) free(x->impulseRightASpectra);
	if(x->impulseLeftBSpectra != 0) free(x->impulseLeftBSpectra);
	if(x->impulseRightBSpectra != 0) free(x->impulseRightBSpectra);
	
	if(x->outLeftSpectra != 0) free(x->outLeftSpectra);
	if(x->outRightSpectra != 0) free(x->outRightSpectra);

	if(x->filterOne.position != 0) free(x->filterOne.position);
	if(x->filterTwo.position != 0) free(x->filterTwo.position);
	
	x->outBufferL = x->outBufferR = x->inBuffer = 0;
	x->impulseRightA = x->impulseLeftA = x->impulseRightB = x->impulseLeftB = 0;
	x->outputRightA = x->outputLeftA = x->outputRightB = x->outputLeftB = x->inputDouble = 0;
	x->overlapLeftA = x->overlapRightA = x->overlapLeftB = x->overlapRightB = 0;
	x->inSpectra = x->impulseLeftASpectra = x->impulseRightASpectra = x->outLeftSpectra = x->outRightSpectra = 0;
	x->impulseLeftBSpectra = x->impulseRightBSpectra = 0;
}


void setup_0x2bbinaural_tilde(void)
{
    binaural_class = class_new(gensym("+binaural~"), (t_newmethod)binaural_new, (t_method)binaural_free,
    	sizeof(t_binaural), 0, A_DEFFLOAT, 0);
	CLASS_MAINSIGNALIN(binaural_class, t_binaural, x_f);
	class_addmethod(binaural_class, (t_method)binaural_dsp, gensym("dsp"), 0);
	class_addmethod(binaural_class, (t_method)binaural_azimuth, gensym("angle"), A_DEFFLOAT, 0);
	class_addmethod(binaural_class, (t_method)binaural_gain, gensym("gain"), A_DEFFLOAT, 0);
	class_addmethod(binaural_class, (t_method)binaural_filterSet, gensym("filterSet"), A_DEFFLOAT, 0);

}

