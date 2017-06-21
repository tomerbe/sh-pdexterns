#include "m_pd.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

#define TRUE 1.0f
#define FALSE 0.0f

/* ------------------------ spectralcompand~ ----------------------------- */


static t_class *spectralcompand_tilde_class;

// some convenient constants
enum
{
	kLevelIncrements = 96,
	kSizeFFT = 1024,
	kHalfSizeFFT = 512
};


typedef struct _spectralcompand_tilde
	{
		t_object x_obj;         /* obligatory header */
		
		t_float sampleRate;
		t_int   bufferPosition;
		double denormalValue;
		long blockSize;
		// FFT stuff
		float *inBuffer, *outBuffer, *inWindowed;
		float *analysisWindow, *synthesisWindow;
		float *inShift, *outShift;
		long halfSizeFFT, sizeFFT;
		long inputTime, outputTime;
		float pi, twoPi;
		
		// these use the external table
		int x_array_points;
		t_symbol *x_arrayname;					
		t_word *gainTweak;

		// external parameters
		t_float averageLearnValue;
		t_float peakLearnValue;
		t_float invertValue;
		t_float thresholdValue;
		t_float attack;
		t_float release;
		t_float tilt;
		
		// internal parameters
		t_float invert;
		t_float learning;
		t_float learnFrames;
		t_float threshold;
		t_float thresholdSquare;
		t_float ratio;
		t_float attackFactor;
		t_float releaseFactor;
		t_float makeupGain;
		
		// internal tables
		float *threshLearn;

		float *gainTable;
		float *tiltTable;
		float *squareLevel;
		float *innerTweak;
		float *linearTweak;
		
		float *gainIncrements;
		float *levelIncrements;
		
		
	} t_spectralcompand_tilde;


/*these functions clip incoming parameter values*/
void spectralcompand_tilde_set(t_spectralcompand_tilde *x, t_symbol *s) 
{
	t_garray *a;
	int old_array_points;
	int i;
	
	old_array_points = x->x_array_points;
	
	if(!(a = (t_garray *)pd_findbyclass(s, garray_class)))
	{
		pd_error(x, "%s: no such array", s->s_name);
		x->gainTweak = 0;
	}
	else if(!garray_getfloatwords(a, &x->x_array_points, &x->gainTweak))
	{
		pd_error(x, "%s: bad template for spectralcompand", s->s_name);
		x->gainTweak = 0;
	}
	else 
	{
		x->x_arrayname = s;
		garray_resize((t_garray *)pd_findbyclass(x->x_arrayname,garray_class), 513.0);
		garray_redraw((t_garray *)pd_findbyclass(x->x_arrayname,garray_class));
		garray_getfloatwords(a, &x->x_array_points, &x->gainTweak);
		garray_usedindsp(a);
		for(i = 0; i <= x->halfSizeFFT; i++)
			x->innerTweak[i] = x->gainTweak[i].w_float = 1.0f;
	}	
}

void spectralcompand_tilde_learnpeak(t_spectralcompand_tilde *x, t_float value)
{
	int i;
	float threshPeak;
	
	x->peakLearnValue = value;
	// start threshold learning
	if(x->peakLearnValue > 0.5f && x->learning == FALSE)
	{
		for(i = 0; i <= x->halfSizeFFT; i++)
		{
			x->innerTweak[i] = 0.0f;
			x->linearTweak[i] = 1.0f;
			x->threshLearn[i] = 0.0f;
		}
		x->learning = TRUE;
	}
	else if(x->peakLearnValue < 0.5f)
	{
		if(x->learnFrames > 0)
		{
			threshPeak = 0.0f;
			for(i = 0; i < x->halfSizeFFT; i++)
				if(x->threshLearn[i] > threshPeak)
					threshPeak = x->threshLearn[i];
			for(i = 0; i < x->halfSizeFFT; i++)
			{
				x->linearTweak[i] = x->threshLearn[i]/threshPeak;
				x->innerTweak[i] = 20.0f * log10f(x->linearTweak[i]);
				if(x->innerTweak[i] < -96.f)
					x->innerTweak[i] = -96.0f;
				if(x->gainTweak != 0)
					x->gainTweak[i].w_float = x->innerTweak[i];
			}
		}
		x->learning = FALSE;
		x->learnFrames = 0;
	}
	if(x->gainTweak != 0)
		garray_redraw((t_garray *)pd_findbyclass(x->x_arrayname,garray_class));
}

void spectralcompand_tilde_learnavg(t_spectralcompand_tilde *x, t_float value) 
{
	int i;
		
	float threshPeak;

	x->averageLearnValue = value;
	// start threshold learning
	if(x->averageLearnValue > 0.5f && x->learning == FALSE)
	{
		for(i = 0; i <= x->halfSizeFFT; i++)
		{
			x->innerTweak[i] = 0.0f;
			x->linearTweak[i] = 1.0f;
			x->threshLearn[i] = 0.0f;
		}
		x->learning = TRUE;
	}
	else if(x->averageLearnValue < 0.5f)
	{
		if(x->learnFrames > 0)
		{
			threshPeak = 0.0f;
			for(i = 0; i < x->halfSizeFFT; i++)
				if(x->threshLearn[i] > threshPeak)
					threshPeak = x->threshLearn[i];
			for(i = 0; i < x->halfSizeFFT; i++)
			{
				x->linearTweak[i] = x->threshLearn[i]/threshPeak;
				x->innerTweak[i] = 20.0f * log10f(x->linearTweak[i]);
				if(x->innerTweak[i] < -96.f)
					x->innerTweak[i] = -96.0f;
				if(x->gainTweak != 0)
					x->gainTweak[i].w_float = x->innerTweak[i];
			}
		}
		x->learnFrames = 0;
		x->learning = FALSE;
	}	
	if(x->gainTweak != 0)
		garray_redraw((t_garray *)pd_findbyclass(x->x_arrayname,garray_class));
}

// 1 equals TRUE, 0 equals FALSE
void spectralcompand_tilde_invert(t_spectralcompand_tilde *x, t_float value) 
{
	if(value <= 1 && value >= 0) 
		x->invertValue = value;
	else
		error("invert value must be 0 or 1");
	
	if(x->invertValue >= 0.5f)
		x->invert = TRUE;
	else
		x->invert = FALSE;
	
}

void spectralcompand_tilde_thresh(t_spectralcompand_tilde *x, t_float f)
{
	float dbLevel, level;
	long i;

	if(f>0.0 || f<-96.0)
	{
		error("threshold value must be >= -96 and <= 0.");
		return;
	}
	else
		x->thresholdValue = f;		
	x->threshold = powf(10.f, x->thresholdValue * 0.05f);
	x->thresholdSquare = x->threshold * x->threshold;
	if(x->ratio >= 1.0f)
	{
		for(i = 0; i < kLevelIncrements; i++)
		{
			dbLevel = (float)(-96 + i);
			level =  powf(10.f, dbLevel * 0.05f);
			if(level < x->threshold)
				x->gainIncrements[i] = powf(10.0f, (dbLevel - x->thresholdValue) * (x->ratio - 1.0f) * 0.05f);
			else
				x->gainIncrements[i] = 1.0f;
		}
	}
	else
	{
		for(i = 0; i < kLevelIncrements; i++)
		{
			dbLevel = (float)(-96 + i);
			level =  powf(10.f, dbLevel * 0.05f);
			if(level > x->threshold)
				x->gainIncrements[i] = powf(10.0f, (x->thresholdValue - dbLevel) * (1.0f/x->ratio - 1.0f) * 0.05f);
			else
				x->gainIncrements[i] = 1.0f;
		}
	}
}

void spectralcompand_tilde_ratio(t_spectralcompand_tilde *x, t_float value) 
{
	long i;
	float dbLevel, level;
	if(value > 5.0)
	{
		// expansion
		error("ratio value must be >= 0.2 and <= 5.0");
		value = 5.0;
	}
	else if(value < 0.01)
	{
		// compression
		error("ratio value must be >= 0.2 and <= 5.0");
		value = 0.2;
	}
	x->ratio = value;				
	if(x->ratio >= 1.0)
	{
		for(i = 0; i < kLevelIncrements; i++)
		{
			dbLevel = (float)(-96 + i);
			level =  powf(10.f, dbLevel * 0.05f);
			if(level < x->threshold)
				x->gainIncrements[i] = powf(10.0f, (dbLevel - x->thresholdValue) * (x->ratio - 1.0f) * 0.05f);
			else
				x->gainIncrements[i] = 1.0f;
		}
	}
	else
	{
		for(i = 0; i < kLevelIncrements; i++)
		{
			dbLevel = (float)(-96 + i);
			level =  powf(10.f, dbLevel * 0.05f);
			if(level > x->threshold)
				x->gainIncrements[i] = powf(10.0f, (x->thresholdValue - dbLevel) * ((1.0f/x->ratio) - 1.0f) * 0.05f);
			else
				x->gainIncrements[i] = 1.0f;
		}
	}
}

void spectralcompand_tilde_attack(t_spectralcompand_tilde *x, t_float f)
{	
	if(f<0.0)
	{
		error("attack above 0 seconds");
		return;
	}
	else
		x->attack = f;	
	if(x->attack == 0.0f)
		x->attackFactor = 1.0f/0.001f;
	else
		x->attackFactor = powf(1.0f/0.1f, ((float)(x->blockSize)/(x->attack*x->sampleRate)));
}

void spectralcompand_tilde_release(t_spectralcompand_tilde *x, t_float f)
{	
	if(f<0.0)
	{
		error("attack above 0 seconds");
		return;
	}
	else
		x->release = f;	
	if(x->release == 0.0f)
		x->releaseFactor = 0.001f;
	else
		x->releaseFactor = powf(0.1f, ((float)(x->blockSize)/(x->release * x->sampleRate)));
}

void spectralcompand_tilde_tilt(t_spectralcompand_tilde *x , t_float value)
{
    int i;
    float tiltBasis;
    
    if ((value <= 6) && (value >= -6))  
	{ 
		x->tilt = value;	
		tiltBasis = log10f((float)x->halfSizeFFT) * 20.0f * 0.5f;
		
		for(i=0; i <= x->halfSizeFFT; i++)
			x->tiltTable[i] = powf(10.0f, ((log10f((float)i+1) * 20.0f) -  tiltBasis) * (x->tilt/tiltBasis));
		
	} 
	else
		error("tilt between -6 and 6 db/oct.");
}

void spectralcompand_tilde_gain(t_spectralcompand_tilde *x , t_float value)
{
	if ((value <= 24) && (value >= -24))  
		x->makeupGain = powf(10.f, value * 0.05f);
	else
		error("gain between -24 and 24 db");
}




void spectralcompand_tilde_initHammingWindows ( t_spectralcompand_tilde *x ) 
{
	long N, k;
	N = x->sizeFFT;
	for ( k = 0; k < N; k++ )
		x->analysisWindow[k] = x->synthesisWindow[k] = (float) (.54f - (.46f * cosf(x->twoPi * k / (N - 1)) ) );

}


void spectralcompand_tilde_scaleWindows (t_spectralcompand_tilde *x) {
	long i;
	float a, b, sum, analFactor, synthFactor;
	
	a = 0.54f;
	b = 0.46f;
	sum = 0.0f;
	for (i = 0; i < x->sizeFFT; i++)
		sum += x->analysisWindow[i];
	
	synthFactor = analFactor = 2.0f/sum;
	
	for (i = 0; i < x->sizeFFT; i++) {
		x->analysisWindow[i] *= analFactor;
		x->synthesisWindow[i] *= synthFactor;
	}

	sum = 0.0;
	for (i = 0; i < x->sizeFFT; i += x->blockSize)
		sum += x->synthesisWindow[i] * x->synthesisWindow[i];

	sum = 0.5f/(sum * x->sizeFFT);

	for (i = 0; i < x->sizeFFT; i++)
		x->synthesisWindow[i] *= sum;
}




void *spectralcompand_tilde_new(t_symbol *table) 
{
	long i;
	t_spectralcompand_tilde *x = (t_spectralcompand_tilde *)pd_new(spectralcompand_tilde_class);

	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"),gensym("interptVal")); // make this float inlet?
	outlet_new(&x->x_obj, gensym("signal"));
		
	// add method to set table names
		
	//	expects 1 table names 

	x->gainTweak = 0;

	x->makeupGain = 1.0f;							/* might be neccesary to do some gain compensation */
	x->bufferPosition = 0;
	x->inputTime = x->outputTime = 0;
	x->sampleRate = 44100.0f;		/* just to make sure */
	x->sizeFFT = kSizeFFT;
	x->blockSize = x->sizeFFT >> 2;
	x->halfSizeFFT = x->sizeFFT >> 1;
	x->pi = 4.0f * atanf(1.0f);
	x->twoPi = 8.0f * atanf(1.0f);

	
	// lots of memory to set up
	//  first - preset all pointers to zero
	x->outBuffer = x->inBuffer = x->inWindowed = x->inShift = x->outShift = 0;
	// second - allocate a huge number of pointers
		
	x->inBuffer = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->inWindowed = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->inShift = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outShift = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->outBuffer = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->synthesisWindow = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->analysisWindow = (float *) malloc(sizeof(float) * x->sizeFFT);
	
	// third - zero out all the memory
	for(i = 0; i<x->sizeFFT; i++)
		x->inBuffer[i] = x->outBuffer[i] = x->inShift[i] = x->outShift[i] = 0.0;
	
	// spectralcompand specific stuff
	x->gainTable = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->squareLevel = (float *) malloc(sizeof(float) * x->sizeFFT);		
	x->tiltTable = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->threshLearn = (float *) malloc(sizeof(float) * x->sizeFFT);

	for(i = 0; i < x->sizeFFT; i++)
	{
		x->threshLearn[i] = 0.0f;
		x->squareLevel[i] = 0.0f;
		x->gainTable[i] = 1.0f;
	}
	x->levelIncrements = (float *) malloc(sizeof(float) * kLevelIncrements);
	x->gainIncrements = (float *) malloc(sizeof(float) * kLevelIncrements);
	for(i = 0; i < kLevelIncrements; i++)
		x->levelIncrements[i] =  powf(10.f, (float)(-96 + i) * 0.1f);

	// spectralcompand external specific stuff
	x->innerTweak = (float *) malloc(sizeof(float) * 513);
	x->linearTweak = (float *) malloc(sizeof(float) * 513);
	for(i = 0; i < 513; i++)
	{
		x->innerTweak[i] = 0.0f;
		x->linearTweak[i] = 1.0f;
	}
	
	//        fourth - set up the FFT and windows
	spectralcompand_tilde_initHammingWindows(x);		// 
	spectralcompand_tilde_scaleWindows(x);				// these we keep
	
	x->averageLearnValue = 0.0f;
	x->peakLearnValue = 0.0;
	x->learning = FALSE;
	x->learnFrames = 0;
	spectralcompand_tilde_invert(x, 0.0);
	spectralcompand_tilde_thresh(x, -30.0);
	spectralcompand_tilde_attack(x, 0.1);
	spectralcompand_tilde_release(x, 0.5);
	spectralcompand_tilde_tilt(x, 0.0);
	spectralcompand_tilde_gain(x, 0.0);	
	if(table) 
	{  
		x->x_array_points = 0;	
		spectralcompand_tilde_set(x,table); 
	}
	
	return (x);
}


void spectralcompand_tilde_processSpect(t_spectralcompand_tilde *x) 
{	
	long i, j;
	float gain;
	
	// first - get the levels for all bands and identify the peaks
	// left channel
	x->squareLevel[0] = x->inWindowed[0] * x->inWindowed[0];
	x->squareLevel[x->halfSizeFFT] = x->inWindowed[x->halfSizeFFT] * x->inWindowed[x->halfSizeFFT];
	for(i = 1; i < x->halfSizeFFT; i++)
	{
		x->squareLevel[i] = x->inWindowed[i] * x->inWindowed[i] + x->inWindowed[x->sizeFFT - i] * x->inWindowed[x->sizeFFT - i];
	}
	if(x->learning == TRUE)
	{
		for(i = 0; i<x->halfSizeFFT; i++)
		{
			if(x->peakLearnValue > 0.5f)
			{
				if(x->squareLevel[i] > x->threshLearn[i])
					x->threshLearn[i] = x->squareLevel[i];
			}
			else
				x->threshLearn[i] = x->threshLearn[i] + (x->squareLevel[i] * 0.05f);
		}
		x->learnFrames++;
		if(x->learnFrames > 40)
			x->learning = FALSE;
	}
	
	// now we reset the gain table, multiplying by the attackFactor or releaseFactor depending on whether 
	// above or below the threshold
	for(i = 0; i <= x->halfSizeFFT; i++)
	{
			
		if((((x->squareLevel[i] * x->linearTweak[i] * x->tiltTable[i] * x->tiltTable[i]) < (x->thresholdSquare)) && (x->ratio >= 1.0f)) ||
		   (((x->squareLevel[i] * x->linearTweak[i] * x->tiltTable[i] * x->tiltTable[i]) > (x->thresholdSquare )) && (x->ratio < 1.0f)))
		{
			x->gainTable[i] *= x->releaseFactor;
			for(j = 0; j < kLevelIncrements; j++)
				if((x->squareLevel[i] * x->linearTweak[i]) < x->levelIncrements[j])
				{
					gain = x->gainIncrements[j];
					break;
				}
			if(x->gainTable[i] < gain)
				x->gainTable[i] = gain;
		}
		else
		{
			x->gainTable[i] *= x->attackFactor;
			if(x->gainTable[i] > 1.0f)
				x->gainTable[i] = 1.0f;
		}
	}
	if(x->invert == TRUE)
	{
	    x->outShift[0] = x->inWindowed[0] * (1.0f - x->gainTable[0]);	// DC Component
	    x->outShift[x->halfSizeFFT] = x->inWindowed[x->halfSizeFFT] * (1.0f - x->gainTable[x->halfSizeFFT]);	// Nyquist Frequency
	    for(i = 1; i < x->halfSizeFFT; ++i)
	    {
	    	x->outShift[i] = x->inWindowed[i] * (1.0f - x->gainTable[i]);
	    	x->outShift[x->sizeFFT - i] = x->inWindowed[x->sizeFFT - i] * (1.0f - x->gainTable[i]);
	    }
    }
	else
	{
		// c - gain the spectra
	    x->outShift[0] = x->inWindowed[0] * x->gainTable[0];	// DC Component
	    x->outShift[x->halfSizeFFT] = x->inWindowed[x->halfSizeFFT] * x->gainTable[x->halfSizeFFT];	// Nyquist Frequency
	    for(i = 1; i < x->halfSizeFFT; ++i)
	    {
	    	x->outShift[i] = x->inWindowed[i] * x->gainTable[i];
	    	x->outShift[x->sizeFFT - i] = x->inWindowed[x->sizeFFT - i] * x->gainTable[i];
	    }
	}
}


void spectralcompand_tilde_block(t_spectralcompand_tilde *x) {
	long i;
	//	long j;
	long maskFFT = x->sizeFFT - 1;
	float tweakSum;
	
	if(x->gainTweak != 0)
	{
		tweakSum = 0.0f;
		for(i = 0; i < 513; i++)
			tweakSum += (x->innerTweak[i] - x->gainTweak[i].w_float);
		if(tweakSum != 0.0f)
		{
			for(i = 0; i < 513; i++)
			{
				x->innerTweak[i] = x->gainTweak[i].w_float;
				x->linearTweak[i] = powf(10.f, x->innerTweak[i] * 0.05f);
			}
		}
	}

			
		
	// shift data in the outBuffer toward the beginning of the buffer
	memcpy(x->outBuffer, x->outBuffer+x->blockSize, (x->sizeFFT - x->blockSize) * sizeof(float));
	// zero out the end of the outBuffer
	memset(x->outBuffer+(x->sizeFFT - x->blockSize), 0, x->blockSize *  sizeof(float));
		
	// shift data in the inShift buffer toward the beginning of the buffer
	memcpy(x->inShift, x->inShift+x->blockSize, (x->sizeFFT - x->blockSize) * sizeof(float));
	// put new samples in the end of the buffer
	memcpy(x->inShift + (x->sizeFFT - x->blockSize), x->inBuffer, x->blockSize * sizeof(float));

	//window our input samples in preparation for FFT
	for(i = 0; i < x->sizeFFT; i++) {
		*(x->inWindowed + x->inputTime) = *(x->inShift + i) * *(x->analysisWindow + i);
		++(x->inputTime);
		x->inputTime = x->inputTime & maskFFT;
	}
	
	mayer_realfft(x->sizeFFT, x->inWindowed);
		
	// use averaged table later, for now, this is just straight table 1		
	spectralcompand_tilde_processSpect(x);
	
	mayer_realifft(x->sizeFFT, x->outShift);
		
	// now copy the output into the output buffer, multiplying by the window
	for(i = 0; i < x->sizeFFT; i++) {
		*(x->outBuffer + i) += *(x->outShift + x->outputTime) * *(x->synthesisWindow + i);
		++x->outputTime;
		// this uses a mask to do the modulo addressing
		x->outputTime = x->outputTime & maskFFT;
	}
}
/* called to start DSP.  Here we call Pd back to add our perform
 routine to a linear callback list which Pd in turn calls to grind
 out the samples. */


/* this is the actual performance routine which acts on the  samples.
   It's called with a single pointer "w" which is our location in the
   DSP call list.  We return a new "w" which will point to the next item
   after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
   (no use to us), w[1] and w[2] are the input and output vector  
   locations, and w[3] is the number of points to calculate. */
static t_int *spectralcompand_tilde_perform(t_int *w) {
		
	t_spectralcompand_tilde *x = (t_spectralcompand_tilde *)(w[1]);
	t_sample *in = (t_float *)(w[2]);
	//	t_float *freqResponse = (t_float *)(w[3]);
	t_sample *out = (t_float *)(w[3]);
	int n = (int)(w[4]);
		
	// i like counting from zero, so i use sample to count the offset from
	// the start of the in and out blocks
	int i;
	//	long frames;
	long framesLeft, processframes;
	long bandsPerIndex;
		
	// gotta copy things from the filter input first off
	// since the block size doesn't always match the number of frequency  bands
	// lets adjust in our copy
	bandsPerIndex = kHalfSizeFFT/n;
		
		
	// framesLeft is the number of samples that have to be copied from *in
	// this loop continues to copy from in until no samples are left

	framesLeft = n;
		
	while ( framesLeft > 0 ) {
		if (framesLeft + x->bufferPosition < x->blockSize)
			processframes = framesLeft;
		else
			processframes = x->blockSize - x->bufferPosition;
		// flush out previous output, copy in new input...
		// x->bufferPosition is used as a way to keep track of position in both
		// x->inBuffer and x->outBuffer
			
		memcpy(x->inBuffer+(x->bufferPosition), in, processframes *  sizeof(float));

		for (i=0; i<processframes; i++)	
		{
			out[i] = x->outBuffer[i+x->bufferPosition] * x->makeupGain;
			if(out[i] > 1000.0f)
				out[i] = 1000.0f;
			else if(out[i] < -1000.0f)
				out[i] = -1000.0f;
		}
		// increment in and out pointers        
		out += processframes;
		in += processframes;
		// increment bufferPostion, if the bufferPosition hits the blockSize  (1/4 of FFT size)
		// perform another FFT.

		x->bufferPosition += processframes;
		if (x->bufferPosition >= x->blockSize){

			x->bufferPosition = 0;
			spectralcompand_tilde_block(x);
		}
		// decrement framesLeft by the number of frames (samples) processed
		framesLeft -= processframes;
	}
	return (w + 5);
}


void spectralcompand_tilde_dsp(t_spectralcompand_tilde *x, t_signal **sp) {
	x->sampleRate = sp[0]->s_sr;
		
	dsp_add(spectralcompand_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}
   	
	
// since we allocated some memory, we need a delete function
static void spectralcompand_tilde_free(t_spectralcompand_tilde *x) 
{
	if(x->inBuffer != 0) free(x->inBuffer);
	if(x->inWindowed != 0) free(x->inWindowed);
	if(x->outBuffer != 0) free(x->outBuffer);
	if(x->inShift != 0) free(x->inShift);
	if(x->outShift != 0) free(x->outShift);
	
	if(x->analysisWindow != 0) free(x->analysisWindow);
	if(x->synthesisWindow != 0) free(x->synthesisWindow);
	
	if(x->gainTable != 0) free(x->gainTable);
	if(x->tiltTable != 0) free(x->tiltTable);
	
	if(x->threshLearn != 0) free(x->threshLearn);
	if(x->squareLevel != 0) free(x->squareLevel);		
	if(x->levelIncrements != 0) free(x->levelIncrements);
	if(x->gainIncrements != 0) free(x->gainIncrements);		
	// spectralcompand external specific stuff
	if(x->innerTweak != 0) free(x->innerTweak);		
	if(x->linearTweak != 0) free(x->linearTweak);			
}

/* this routine, which must have exactly this name (with the "~"  replaced 
 * by "_tilde) is called when the code is first loaded, and tells Pd  how to build the "class". */
      

void setup_0x2bspectralcompand_tilde(void) {
	spectralcompand_tilde_class = class_new(gensym("+spectralcompand~"),  
				                 (t_newmethod)spectralcompand_tilde_new, 
								 (t_method)spectralcompand_tilde_free,
								 sizeof(t_spectralcompand_tilde), 
								 CLASS_DEFAULT, 
								 A_DEFSYMBOL,
								 0);
	/* this is magic to declare that the leftmost, "main" inlet
	 * takes signals; other signal inlets are done differently... */

	/* also installs delay_time as the leftmost inlet float */	

	CLASS_MAINSIGNALIN(spectralcompand_tilde_class, t_spectralcompand_tilde, makeupGain);
	/* here we tell Pd about the "dsp" method, which is called back when DSP is turned on. */

	class_addmethod(spectralcompand_tilde_class, (t_method) spectralcompand_tilde_dsp, gensym("dsp"), (t_atomtype)0);
	
	/* a set function to set a threshold table */
	class_addmethod(spectralcompand_tilde_class, (t_method) spectralcompand_tilde_set, gensym("set"), A_SYMBOL, 0);

	class_addmethod(spectralcompand_tilde_class, (t_method) spectralcompand_tilde_learnpeak, gensym("learnpeak"), A_DEFFLOAT,0);
	class_addmethod(spectralcompand_tilde_class, (t_method) spectralcompand_tilde_learnavg, gensym("learnavg"), A_DEFFLOAT,0);
	class_addmethod(spectralcompand_tilde_class, (t_method) spectralcompand_tilde_invert, gensym("invert"), A_DEFFLOAT,0);
		
	class_addmethod(spectralcompand_tilde_class, (t_method)spectralcompand_tilde_thresh, gensym("threshold"), A_DEFFLOAT, 0);
	class_addmethod(spectralcompand_tilde_class, (t_method)spectralcompand_tilde_ratio, gensym("ratio"), A_DEFFLOAT,0);
	class_addmethod(spectralcompand_tilde_class, (t_method)spectralcompand_tilde_attack, gensym("attack"), A_DEFFLOAT,0);
	class_addmethod(spectralcompand_tilde_class, (t_method)spectralcompand_tilde_release, gensym("release"), A_DEFFLOAT,0);
	class_addmethod(spectralcompand_tilde_class, (t_method)spectralcompand_tilde_tilt, gensym("tilt"), A_DEFFLOAT, 0); 		
	class_addmethod(spectralcompand_tilde_class, (t_method)spectralcompand_tilde_gain, gensym("gain"), A_DEFFLOAT, 0); 		
}