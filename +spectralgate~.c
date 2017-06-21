/*
 *  spectralgate~.c
 *
 *
 */


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

static t_class *spectralgate_class;

enum
{
	kGate,
	kDuck
};


enum
{
	kSizeFFT = 1024,
	kHalfSizeFFT = 512
};


typedef struct _spectralgate
{
	t_object x_obj;         /* obligatory header */	
	t_float sampleRate;
	long blockSize;
	
	t_int   bufferPosition;

	t_float gain;
	
	// FFT stuff
	float *inBuffer, *outBuffer, *inWindowed;
	float *analysisWindow, *synthesisWindow;
	float *inShift, *outShift;
	long halfSizeFFT, sizeFFT;
	long inputTime, outputTime;
	float pi, twoPi;

	/* external table */

	t_word *gainTweak;
	int x_array_points;
	t_symbol *x_arrayname;
	
	 
	/* Parameter variables  from setParameter */
		
	t_float typeValue;
	t_float learnValue;
	t_float learnFrames;
	
	t_float * threshTable;
	t_float * threshLearn;

	t_float learning;
	t_float resetValue;
	t_float threshAverage;
	
	t_float peakTrackValue;
	t_float peakTrack;
	t_float thresholdValue;
	t_float threshold;
	t_float thresholdSquare;
	t_float attackValue;
	t_float attack;
	t_float releaseValue;
	t_float release;
	t_float gateGain;
	t_float highGain;
	t_float lowGain;
	
	t_float makeupGainValue;
	t_float makeupGain;
	
	t_float tiltValue;
	t_float tilt;
	t_float tiltBasis;

	t_float releaseFactor;
	t_float attackFactor;
	
	t_float type;
	
	/* internal tables */
	t_float * tiltTable;
	t_float * squareLevel;
	t_float * gainTable;	
	t_float *innerTweak;	
	
} t_spectralgate;

void spectralgate_tilde_set(t_spectralgate *x, t_symbol *s) {
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
		pd_error(x, "%s: bad template for spectralgate", s->s_name);
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
			x->innerTweak[i] = x->gainTweak[i].w_float = -48.0f;
	}	
}   	


void spectralgate_tilde_type(t_spectralgate *x, t_float value) 
{
	x->typeValue = value;
	
	if(x->typeValue < 0.5f)
		x->type = kGate;
	else
		x->type = kDuck;
}

void spectralgate_tilde_learn (t_spectralgate *x, t_float value) 
{
	int i;
	float threshPeak;
	
	x->learnValue = value;
	// start threshold learning
	if(x->learnValue > 0.5f && x->learning == FALSE)
	{
		for(i = 0; i <= x->halfSizeFFT; i++)
		{
			x->innerTweak[i] = 0.0f;
			x->threshTable[i] = 1.0f;
			x->threshLearn[i] = 0.0f;
		}
		x->learning = TRUE;
	}
	else if(x->learnValue < 0.5f)
	{
		if(x->learnFrames > 0)
		{
			threshPeak = 0.0f;
			for(i = 0; i < x->halfSizeFFT; i++)
				if(x->threshLearn[i] > threshPeak)
					threshPeak = x->threshLearn[i];
			for(i = 0; i < x->halfSizeFFT; i++)
			{
				x->threshTable[i] = x->threshLearn[i]/threshPeak;
				x->innerTweak[i] = 20.0f * log10f(x->threshTable[i]);
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

void spectralgate_tilde_reset(t_spectralgate *x, t_float value)
{
	int i;
	x->resetValue = value;
	for(i = 0; i <= x->halfSizeFFT; i++) 
	{
		x->innerTweak[i] = 0.0f;
		x->threshTable[i] = 1.0f;
		x->threshLearn[i] = 0.0f;
		if(x->gainTweak != 0)
			x->gainTweak[i].w_float = x->innerTweak[i];
	}
	if(x->gainTweak != 0)
		garray_redraw((t_garray *)pd_findbyclass(x->x_arrayname,garray_class));
}

void spectralgate_tilde_peaktrack(t_spectralgate *x, t_float value) {
	x->peakTrackValue = value;
	
	if(x->peakTrackValue < 0.5f)
		x->peakTrack = FALSE;
	else
		x->peakTrack = TRUE;
}


void spectralgate_tilde_threshold(t_spectralgate *x, t_float value) {
	
	if(value>0.0 || value<-96.0)
	{
		error("threshold value must be >= -96 and <= 0.");
		return;
	}
	
	x->thresholdValue = value;
	x->threshold = powf(10.f, x->thresholdValue * 0.05f);
	x->thresholdSquare = x->threshold * x->threshold;
}

void spectralgate_tilde_attack(t_spectralgate *x, t_float value) {
	
	if(value <0.0)
	{
		error("attack above 0 seconds");
		return;
	}
	else 
		x->attack = value;
	if(x->attack == 0.0f)
		x->attackFactor = 1.0f/0.001f;
	else
		x->attackFactor = powf(1.0f/0.1f, ((float)(x->blockSize)/(x->attack*x->sampleRate)));
	x->attackFactor *= x->attackFactor;
}

void spectralgate_tilde_release(t_spectralgate *x, t_float value) {	
	if(value<0.0)
	{
		error("release above 0 seconds");
		return;
	}
	
	else
		x->release = value;	
	if(x->release == 0.0f)
		x->releaseFactor = 0.001f;
	else
		x->releaseFactor = powf(0.1f, ((float)(x->blockSize)/(x->release * x->sampleRate)));
	x->releaseFactor = x->releaseFactor * x->releaseFactor;
}

void spectralgate_tilde_gain(t_spectralgate *x, t_float value)
{
	
	if ((value <= 60) && (value >= -60))
	{
		x->gateGain = powf(10.f, value * 0.05f);
		if(x->gateGain >= 1.0f) 
		{
			x->highGain = x->gateGain;
			x->lowGain = 1.0f;
		}
		else 
		{
			x->highGain = 1.0f;
			x->lowGain = x->gateGain;
		}
	}
	else
		error("gain between -60 and 60 db");
}
void spectralgate_tilde_makeupgain(t_spectralgate *x, t_float value) 
{
	
	if ((value <= 24) && (value >= -24))  
		x->makeupGain = powf(10.f, value * 0.05f);
	else
		error("makeupgain between -24 and 24 db");
}

void spectralgate_tilde_tilt(t_spectralgate *x, t_float value) 
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

void spectralgate_tilde_processSpect(t_spectralgate *x)
{
	long i;
	
	float	*triggerSpectra;
	float	*localThresh;
	
	triggerSpectra = x->inWindowed;
	localThresh = x->threshTable;
	
	
	// if peak detection is on, the threshold is based on the peak value plus the threshold value
	// the range of threshold values changes in the GUI from 0 to -96 (no peak) to +/- 48 (peak),
	// so we can just combine the threshold with the peak
	
	// first - get the levels for all bands and identify the peaks
	// left channel
	
	float peakSquared = x->squareLevel[0] = triggerSpectra[0] * triggerSpectra[0];
	x->squareLevel[x->halfSizeFFT] = triggerSpectra[x->halfSizeFFT] * triggerSpectra[x->halfSizeFFT];
	if(x->squareLevel[x->halfSizeFFT] > peakSquared)
		peakSquared = x->squareLevel[x->halfSizeFFT];
	for(i = 1; i < x->halfSizeFFT; i++)
	{
		x->squareLevel[i] = triggerSpectra[i] * triggerSpectra[i] + triggerSpectra[x->sizeFFT - i] * triggerSpectra[x->sizeFFT - i];
		if(x->squareLevel[i] > peakSquared)
			peakSquared = x->squareLevel[i];
	}
	
	if(x->learning == TRUE)
	{
		for(i = 0; i < x->halfSizeFFT; i++)
		{
			if(x->squareLevel[i] > x->threshLearn[i])
				x->threshLearn[i] = x->squareLevel[i];
			
		}
		x->learnFrames++;
		if(x->learnFrames > 40)
			x->learnValue = 0.0f;
	}
	if(x->peakTrack == FALSE)
	{
		peakSquared = peakSquared = 1.0f;
	}
	
	// now we reset the gain table, multiplying by the attackFactor or releaseFactor depending on whether 
	// above or below the threshold
	// left channel
	for(i = 0; i <= x->halfSizeFFT; i++)
	{
		if(((x->squareLevel[i] < (x->thresholdSquare * localThresh[i] * peakSquared * x->tiltTable[i] * x->tiltTable[i])) && (x->type == kGate))
		   || ((x->squareLevel[i] > (x->thresholdSquare * localThresh[i] * peakSquared * x->tiltTable[i] * x->tiltTable[i])) && (x->type == kDuck)) )
			x->gainTable[i] *= x->releaseFactor;
		else
			x->gainTable[i] *= x->attackFactor;
		if(x->gainTable[i] > x->highGain)
			x->gainTable[i] = x->highGain;
		else if(x->gainTable[i] < x->lowGain)
			x->gainTable[i] = x->lowGain;
	}
	
	
	// c - gain the spectra
    x->outShift[0] = x->inWindowed[0]* x->gainTable[0];	// DC Component
    x->outShift[x->halfSizeFFT] = x->inWindowed[x->halfSizeFFT] * x->gainTable[x->halfSizeFFT];	// Nyquist Frequency
    for(i = 1; i < x->halfSizeFFT; ++i)
    {
    	x->outShift[i] = x->inWindowed[i] * x->gainTable[i];
    	x->outShift[x->sizeFFT - i] = x->inWindowed[x->sizeFFT - i] * x->gainTable[i];
    }
	
}

void spectralgate_tilde_block(t_spectralgate *x) {
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
				x->threshTable[i] = powf(10.f, x->innerTweak[i] * 0.05f);
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
	
	
	spectralgate_tilde_processSpect(x);
	
	mayer_realifft(x->sizeFFT, x->outShift);
	
	// now copy the output into the output buffer, multiplying by the window
	for(i = 0; i < x->sizeFFT; i++) {
		*(x->outBuffer + i) += *(x->outShift + x->outputTime) * *(x->synthesisWindow + i);
		++x->outputTime;
		// this uses a mask to do the modulo addressing
		x->outputTime = x->outputTime & maskFFT;
	}
}



t_int *spectralgate_perform(t_int *w) {
	
	t_spectralgate *x = (t_spectralgate *)(w[1]);
	t_float *in = (t_float *)(w[2]);
	t_float *out = (t_float *)(w[3]);
	int n = (int)(w[4]);
	
	// i like counting from zero, so i use sample to count the offset from
	// the start of the in and out blocks
	int i;
//	int  j;
	//	long frames;
	long framesLeft, processframes;
	

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
		
		for (i=0; i<processframes; i++)	{
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
			spectralgate_tilde_block(x);
		}
		// decrement framesLeft by the number of frames (samples) processed
		framesLeft -= processframes;
	}
	

	return (w + 5);
}	



void spectralgate_tilde_initHammingWindows ( t_spectralgate *x ) {
	long N,k;

	N = x->sizeFFT;
	for (k = 0; k < x->sizeFFT; k++ )
		x->analysisWindow[k] = x->synthesisWindow[k] = (float) (.54f - (.46f * cosf(x->twoPi * k / (x->sizeFFT - 1)) ) );
	
}


void spectralgate_tilde_scaleWindows (t_spectralgate *x) {
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




void spectralgate_dsp(t_spectralgate *x, t_signal **sp) {
	x->sampleRate = sp[0]->s_sr;
	dsp_add(spectralgate_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}


void *spectralgate_new(t_symbol *table1) {
	long i;
	t_spectralgate *x = (t_spectralgate *)pd_new(spectralgate_class);

	// an inlet to bring in filter coefficients

	//	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("signal"),0); // make this float inlet?

	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"),gensym("interptVal")); // make this float inlet?

	outlet_new(&x->x_obj, gensym("signal"));

	x->sizeFFT = kSizeFFT;
	x->blockSize = x->sizeFFT >> 2;
	x->halfSizeFFT = x->sizeFFT >> 1;
	x->sampleRate = 44100.0f;		/* just to make sure */
	x->makeupGain = 1.0f;
	x->pi = 4.0f * atanf(1.0f);
	x->twoPi = 8.0f * atanf(1.0f);


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

	x->threshTable = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->innerTweak = (float * ) malloc(sizeof(float) * 513);

	for(i = 0; i < 513; i++)
	{
		x->innerTweak[i] = 0.0f;
		x->threshTable[i] = 1.0f;
	}
	
	x->learnValue = 0.0f;
	x->learning = FALSE;
	x->learnFrames = 0;
	

	spectralgate_tilde_initHammingWindows( x );		// 
	spectralgate_tilde_scaleWindows(x);				// these we keep
	
	spectralgate_tilde_peaktrack(x, 0.0);
	spectralgate_tilde_threshold(x, -6.0);
	spectralgate_tilde_attack(x, 1.0);
	spectralgate_tilde_release(x, 1.0);
	spectralgate_tilde_gain(x, 0.0);
	spectralgate_tilde_makeupgain(x, 0.0);
	spectralgate_tilde_tilt(x, 0.0);
	
	if(table1) {  
		x->x_array_points = 0;	
		spectralgate_tilde_set(x,table1); 
	}
	
	
	return (x);
}

void spectralgate_free(t_spectralgate *x) {
	if(x->inBuffer) free(x->inBuffer);
	if(x->inShift) free(x->inShift);
	if(x->outShift) free(x->outShift);
	if(x->outBuffer) free(x->outBuffer);

	if(x->tiltTable) free(x->tiltTable);
	if(x->analysisWindow) free(x->analysisWindow);
	if(x->synthesisWindow) free(x->synthesisWindow);
	if(x->threshLearn) free(x->threshLearn);	
	if(x->squareLevel) free(x->squareLevel) ;
	if(x->gainTable) free(x->gainTable);
	if(x->innerTweak) free(x->innerTweak);
	if(x->threshTable) free(x->threshTable);
	
}

void setup_0x2bspectralgate_tilde(void) {
	spectralgate_class = class_new(gensym("+spectralgate~"),  
								  (t_newmethod)spectralgate_new, 
								  (t_method)spectralgate_free,
								  sizeof(t_spectralgate), 
								  CLASS_DEFAULT, 
								  A_DEFSYMBOL,
								  A_DEFSYMBOL,
								  0);
	/* this is magic to declare that the leftmost, "main" inlet
	 * takes signals; other signal inlets are done differently... */
	
	/* also installs delay_time as the leftmost inlet float */	
	
	CLASS_MAINSIGNALIN(spectralgate_class, t_spectralgate, makeupGain);

	/* here we tell Pd about the "dsp" method, which is called back when DSP is turned on. */
	
	class_addmethod(spectralgate_class, (t_method) spectralgate_dsp,  
					gensym("dsp"), (t_atomtype)0); 
	
	
	
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_set, gensym("set"), A_SYMBOL, 0);

	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_type, gensym("type"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_learn, gensym("learn"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_reset, gensym("reset"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_peaktrack, gensym("peaktrack"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_threshold, gensym("threshold"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_attack, gensym("attack"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_release, gensym("release"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_gain, gensym("gain"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_makeupgain, gensym("makeupgain"), A_DEFFLOAT, 0);
	class_addmethod(spectralgate_class, (t_method) spectralgate_tilde_tilt, gensym("tilt"), A_DEFFLOAT, 0);
}	

