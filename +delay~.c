/*

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 1:45PM - removed unused variables.

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "m_pd.h"
#define LFOTABLESIZE 1048576
#define ENVARRAYSIZE 16384
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *delay_tilde_class;

typedef struct _delay_tilde
{
	t_object x_obj;
	float x_twoPi;
	float x_sr;
	float x_n;
	float x_oneOverFrames;
	
	float x_envelopeModulate;
	float x_fEnvelopeModulate;
	long x_envelopePosition;
	long x_envelopeArraySize;
    float x_envelopeArrayTotal;
	float  x_envelopeArray[ENVARRAYSIZE];
	
	long x_readPos;
	long x_readPos2;
	long x_writePos;
	long x_delaySize;
	long x_delayMask;
	
	float x_time;
	float x_fTime;
	float x_timeInc;
	
	float x_feedback;
	float x_fFeedback;
	
	float x_lfoOut;
	float x_lfoOut2;
	float x_lfoRandInc;
	float x_lfoPhase;
	float x_lfoPhase2;
	float x_lfoIncrement;
	float x_lfoTableSize;

	float x_modSpeed;
	float x_fModSpeed;
	float x_modDepth;
	float x_fModDepth;	
	float x_depthInc;
	int x_depthIncTicks;
	float x_modPhase;
	float x_fModPhase;
	float x_modShape;
	float x_fModShape;
	
	float x_delay[LFOTABLESIZE];
	float x_lfoWaveform[LFOTABLESIZE];
	float x_sineWaveform[LFOTABLESIZE];
	float x_squareWaveform[LFOTABLESIZE];
	float x_triWaveform[LFOTABLESIZE];
	float x_rampUpWaveform[LFOTABLESIZE];
	float x_rampDownWaveform[LFOTABLESIZE];
		
	float x_filterFreq;
	float x_fFilterFreq;
	float x_resonance;
	float x_fResonance;

	int x_feedbackSwitch;
	int x_freeze;
	int x_twoHead;

	float X1;
	float X2;
	float Y1;
	float Y2;	
	float dcOut1;
	float dcIn1;
	
    float x_f;
} t_delay_tilde;

/***************** UTILITY FUNCTIONS *****************/
static float delay_tilde_envelopeFollower(t_delay_tilde *x, float in)
{
    if(x->x_envelopePosition >= x->x_envelopeArraySize)
            x->x_envelopePosition = 0;
	if(x->x_envelopePosition < 0)
            x->x_envelopePosition = 0;
// the in is the sum of two channels
// FIRST: rectify the input
	if(in < 0.0)
		in = -in;
		
// SECOND: add to array to calculate mean
    x->x_envelopeArrayTotal = x->x_envelopeArrayTotal - x->x_envelopeArray[x->x_envelopePosition] + in;
    
    x->x_envelopeArray[x->x_envelopePosition] = in;
    x->x_envelopePosition++;
    
// THIRD: mean is total/arraysize
    return(x->x_envelopeArrayTotal/(float)x->x_envelopeArraySize);
};

// from http://www.musicdsp.org/archive.php?classid=5#93
static float delay_tilde_hermite(float xx, float yy0, float yy1, float yy2, float yy3)
{
        // 4-point, 3rd-order Hermite (x-form)
        float c0 = yy1;
        float c1 = 0.5f * (yy2 - yy0);
        float y0my1 = yy0 - yy1;
        float c3 = (yy1 - yy2) + 0.5f * (yy3 - y0my1 - yy2);
        float c2 = y0my1 + c1 - c3;

        return ((c3 * xx + c2) * xx + c1) * xx + c0;
}
/***************** UTILITY FUNCTIONS *****************/

static void delay_tilde_time(t_delay_tilde *x, t_floatarg t)
{	
	if(t>5000.0 || t<1.0)
		error("time value must be >= 1.0 ms and <= 5000.0 ms.");
	else
	{	// must fix this - it's scaled twice below
		x->x_time = t/1000.0; // convert from ms to sec
		x->x_timeInc = (x->x_time - x->x_fTime)/(x->x_sr * 0.1f);
	}
}

static void delay_tilde_feedback(t_delay_tilde *x, t_floatarg f)
{	
	if(f>200.0 || f<0.0)
		error("feedback value must be >= 0%% and <= 200%%.");
	else
		x->x_feedback = f/200.0;
}

static void delay_tilde_modSpeed(t_delay_tilde *x, t_floatarg f)
{	
	if(f>10.0 || f<0.01)
		error("lfoFreq value must be >= 0.01 Hz and <= 10.0 Hz.");
	else
		x->x_modSpeed = f;
}

static void delay_tilde_modDepth(t_delay_tilde *x, t_floatarg d)
{	
	float diff;
	
	diff = (x->x_modDepth - x->x_fModDepth);
	
	if(d>100.0 || d<0.0)
		error("lfoDepth value must be >= 0.0%% and <= 100.0%%.");
	else
	{
		x->x_modDepth = d/100.0;
		x->x_depthInc = diff/x->x_sr;
		x->x_depthIncTicks = (int)x->x_sr/x->x_n;
	}
}

static void delay_tilde_modShape(t_delay_tilde *x, t_floatarg s)
{
	int i;
	
	if(s>5.0 || s<0.0)
		error("lfoShape value must be either 0, 1, 2, 3, 4, or 5.");
	else
	{
		x->x_modShape = s/5.0;

		if(x->x_modShape < 0.1666666)
			for(i=0; i<x->x_lfoTableSize; i++)
				x->x_lfoWaveform[i] = x->x_sineWaveform[i];
		else if(x->x_modShape < 0.33333333)
			for(i=0; i<x->x_lfoTableSize; i++)
				x->x_lfoWaveform[i] = x->x_triWaveform[i];
		else if(x->x_modShape < 0.5)
			for(i=0; i<x->x_lfoTableSize; i++)
				x->x_lfoWaveform[i] = x->x_squareWaveform[i];
		else if(x->x_modShape < 0.66666666)
			for(i=0; i<x->x_lfoTableSize; i++)
				x->x_lfoWaveform[i] = x->x_rampUpWaveform[i];
		else if(x->x_modShape < 0.83333333)
			for(i=0; i<x->x_lfoTableSize; i++)
				x->x_lfoWaveform[i] = x->x_rampDownWaveform[i];
	};
				
}

static void delay_tilde_lfoPhase(t_delay_tilde *x, t_floatarg p)
{	
	if(p>180.0 || p < -180.0)
		error("lfoPhase value must be >= -180 deg and <= 180 deg.");
	else
		x->x_lfoPhase = (180.0+p)/360.0;
}

static void delay_tilde_filterFreq(t_delay_tilde *x, t_floatarg ff)
{	
	if(ff>20000.0 || ff<20.0)
		error("filterFreq value must be >= 20 Hz and <= 20000 Hz.");
	else
		x->x_filterFreq = ff;
}

static void delay_tilde_resonance(t_delay_tilde *x, t_floatarg r)
{	
	if(r>100 || r<0)
		error("resonance value must be >= 0%% and <= 100%%.");
	else
		x->x_resonance = r/100.0;
}

static void delay_tilde_feedbackSwitch(t_delay_tilde *x, t_floatarg f)
{
	if(f>1 || f<0)
		error("feedbackSwitch value must be 0 or 1.");
	else
		x->x_feedbackSwitch = (int)f;
}

static void delay_tilde_freeze(t_delay_tilde *x, t_floatarg f)
{
	if(f>1 || f<0)
		error("freeze value must be 0 or 1.");
	else
	x->x_freeze = (int)f;
}

static void delay_tilde_heads(t_delay_tilde *x, t_floatarg h)
{	
	if(h>2.0 || h<1.0)
		error("head value must be 1 or 2.");
	else
		x->x_twoHead = (int)h-1;
}

static void *delay_tilde_new()
{
	t_delay_tilde *x = (t_delay_tilde *)pd_new(delay_tilde_class);
	int i;

	outlet_new(&x->x_obj, &s_signal);
	
	x->x_twoPi = 2.0*M_PI;
	x->x_sr = 44100.0;
	x->x_n = 64.0;
	x->x_oneOverFrames = 1.0f/x->x_n;
	
	x->x_envelopeModulate = x->x_fEnvelopeModulate = 0.0;
	x->x_envelopePosition = 0;
	x->x_envelopeArraySize = ENVARRAYSIZE;
	x->x_envelopeArrayTotal = 0.0;
	
	for(i=0; i<ENVARRAYSIZE; i++)
		x->x_envelopeArray[i] = 0.0;
		
	x->x_readPos = x->x_readPos2 = x->x_writePos = 0.0;
	x->x_delaySize = 1048576;
	x->x_delayMask = 1048576-1;
	

	x->x_time = x->x_fTime = 0.1; // 100 ms
	x->x_timeInc = 0.0;

	x->x_feedback = x->x_fFeedback = 0.5;
	
	x->x_lfoOut = x->x_lfoOut2 = 0.0;
	x->x_lfoRandInc = 0.0;
	x->x_lfoPhase = x->x_lfoPhase2 = 0.75;
	x->x_lfoIncrement = 0.0;
	x->x_lfoTableSize = 1048576;	

	x->x_modSpeed = x->x_fModSpeed = 0.2;
	x->x_modDepth = x->x_fModDepth = 0.0;
	x->x_depthInc = 0.0;
	x->x_depthIncTicks = 0;
	x->x_modPhase = x->x_fModPhase = 0.75;
 	x->x_modShape = x->x_fModShape = 0.0;

	// init delay
	for(i=0; i<x->x_lfoTableSize; i++)
		x->x_delay[i] = 0.0;

	// init lfoWaveform
	for(i=0; i<x->x_lfoTableSize; i++)
		x->x_lfoWaveform[i] = 0.0;

	for(i=0; i<x->x_lfoTableSize; i++)
		x->x_sineWaveform[i] = sin(x->x_twoPi * (float)i/(float)x->x_lfoTableSize);

	for(i=0; i<(x->x_lfoTableSize/2.0); i++)
		x->x_squareWaveform[i] = 1.0;

	for(; i<x->x_lfoTableSize; i++)
		x->x_squareWaveform[i] = -1.0;

	for(i=0; i<(x->x_lfoTableSize/2.0); i++)
		x->x_triWaveform[i] = ((float)i/(float)x->x_lfoTableSize * 4.0) - 1.0f;

	for(; i<x->x_lfoTableSize; i++)
		x->x_triWaveform[i] = 3.0f - ((float)i/(float)x->x_lfoTableSize * 4.0);
		
	for(i=0; i<x->x_lfoTableSize; i++)
		x->x_rampUpWaveform[i] = ((float)i/(float)x->x_lfoTableSize * 2.0) - 1.0f;
	
	for(i=0; i<x->x_lfoTableSize; i++)
		x->x_rampDownWaveform[i] = 1.0 - ((float)i/(float)x->x_lfoTableSize * 2.0);

	for(i=0; i<x->x_lfoTableSize; i++)
		x->x_lfoWaveform[i] = x->x_sineWaveform[i];

	x->x_filterFreq = x->x_fFilterFreq = 20000.0;
	x->x_resonance = x->x_fResonance = 0.5;
	x->x_feedbackSwitch = 1.0; // 1 is on, 0 is off
	
	x->x_freeze = 0;
	x->x_twoHead = 0; // 1 or 2 heads

	// biquad filter stuff
	x->X1 = x->X2 = x->Y1 = x->Y2 = 0.0;
	x->dcOut1 = x->dcIn1 = 0.0f;
	
	return(void *)x;
};

static t_int *delay_tilde_perform(t_int *w)
{
	t_delay_tilde *x = (t_delay_tilde *)(w[1]);
	
	t_sample *in1 = (t_sample *)(w[2]);
	t_sample *out = (t_sample *)(w[3]);
	int n = (int)(w[4]);
	
	long frame;
	float delayTime, delayTimeL;
	float output = 0.0;
	float output2 = 0.0;
	float delayTime2L = 0.0;
	float preClip, preFilter, randTarget;
	float lfoPeriodPerSample;
	float lfoFreq;
	// feedback filter variables
	float filterFreq, filterReson, f0, C;
	float A1, A2, A3, B1, B2;
	float dcOut0, dcIn0, dcR;

	float feedbackInc, speedInc, tempTime, envModInc, lfoPhaseInc;
	float clipFactor = 1.0f;
	float envelopeModulate;
	float feedback;
	float delayFracL = 0.0;
	float delayFrac2L = 0.0;
	long delayLong;
	
	// set up for parameter interpolation
	// delay time interpolation set in setParameter
	feedbackInc = (x->x_feedback - x->x_fFeedback) * x->x_oneOverFrames;
	speedInc = (x->x_modSpeed - x->x_fModSpeed) * x->x_oneOverFrames;
	envModInc = (x->x_envelopeModulate - x->x_fEnvelopeModulate) * x->x_oneOverFrames;
	lfoPhaseInc = (x->x_modPhase - x->x_fModPhase) * x->x_oneOverFrames;
	
	x->x_fModShape = x->x_modShape;
	x->x_fResonance = x->x_resonance;
	x->x_fFilterFreq = x->x_filterFreq;

	dcR = 1.0f - (126.0f/x->x_sr);
	// filter coefficients get calculated only once per block
	filterFreq = x->x_fFilterFreq;
	// calculate new frequency coefficient and damping value
	filterReson = (1.0 - x->x_fResonance) * (sqrt(2.0) - 0.1) + 0.1;

	if(filterFreq > x->x_sr * 0.5)
		filterFreq = x->x_sr * 0.5;
	f0 = filterFreq/x->x_sr;

	if(f0 < 0.1) 
		C = 1.0 / (f0 * M_PI);
	else 
		C = tan((0.5 - f0) * M_PI);
	A1 = 1.0 / ( 1.0 + filterReson * C + C * C);
	A2 = 2.0 * A1;
	A3 = A1;
	B1 = 2.0 * ( 1.0 - C * C) * A1;
	B2 = ( 1.0 - filterReson * C + C * C) * A1;


	// delay processing loop
	for(frame=0; frame<n; frame++)
	{
		// 0 - run the envelope follower
		envelopeModulate = delay_tilde_envelopeFollower(x, *(in1+frame)) * ((x->x_fEnvelopeModulate*6.0));
		if(envelopeModulate > 1.0f)
			envelopeModulate = 1.0f;
		
		// 1 - calculate lfo output
		lfoFreq = x->x_fModSpeed;

		// derive the phase increment from the frequency

		lfoPeriodPerSample = lfoFreq/x->x_sr;
		x->x_lfoIncrement = lfoPeriodPerSample * x->x_lfoTableSize;
		
		// we need a separate bit of code for random ramps
		// we use random when fModShape is between 0.75 and 1.0
		
		if(x->x_fModShape > 0.83333333333f)
		{
			x->x_lfoOut += x->x_lfoRandInc;
		}
		else
		{
			x->x_lfoOut = *(x->x_lfoWaveform + (long)x->x_lfoPhase) * x->x_fModDepth;
		}

		// increment the phase
		x->x_lfoPhase += x->x_lfoIncrement;

		if(x->x_lfoPhase >= (float)x->x_lfoTableSize)
		{
			x->x_lfoPhase -= (float)x->x_lfoTableSize;
			// pick a random number, divide it by the range of randomness
			randTarget = (float)rand()/(float)RAND_MAX; // this is between 0.0 and 1.0
			// change the range to -1.0 to 1.0
			randTarget = ((randTarget * 2.0) - 1.0)  * x->x_fModDepth;
			// find a persample increment (target - current)/number of samples in period
			x->x_lfoRandInc = (randTarget - x->x_lfoOut) * lfoPeriodPerSample;
		}
		// check to see if phase has wrapped
		// that is, the new phase is less than the old

		if(x->x_twoHead == 1)
		{
			if(x->x_fModShape > 0.83333333333f)
			{
				x->x_lfoOut2 += x->x_fModDepth;
				if(x->x_lfoOut2 > x->x_fModDepth)
					x->x_lfoOut2 -= x->x_fModDepth * 2.0;
				else if(x->x_lfoOut2 < -x->x_fModDepth)
					x->x_lfoOut2 += x->x_fModDepth * 2.0;
			}
			else
			{
				x->x_lfoPhase2 = x->x_lfoPhase + (float)x->x_lfoTableSize/2.0;
				if(x->x_lfoPhase2 >= (float)x->x_lfoTableSize)
					x->x_lfoPhase2 -= (float)x->x_lfoTableSize;
				x->x_lfoOut2 = *(x->x_lfoWaveform + (long)x->x_lfoPhase2) * x->x_fModDepth;
			}
		}

		// 2 - determine the delay time and the readPos
		delayTime = x->x_fTime * x->x_sr;

		// a two sample buffer for the hermite interpolation
		delayTimeL = delayTime + (delayTime * x->x_lfoOut);
		delayTimeL = delayTimeL + 2.0 - (delayTimeL * envelopeModulate);
		
		if(x->x_twoHead == 1)
		{
			delayTime2L = delayTime + (delayTime * x->x_lfoOut2);
			delayTime2L = delayTime2L + 2.0 - (delayTime2L * envelopeModulate);
		}

		delayLong = (long)delayTimeL;
		delayFracL = delayTimeL - (float)delayLong;
		delayFracL = 1.0 - delayFracL;
		x->x_readPos = x->x_writePos - delayLong;
		x->x_readPos &= x->x_delayMask;

		if(x->x_twoHead == 1)
		{
			delayLong = (long)delayTime2L;
			delayFrac2L = delayTime2L - (float)delayLong;
			delayFrac2L = 1.0 - delayFrac2L;
			x->x_readPos2 = x->x_writePos - delayLong;
			x->x_readPos2 &= x->x_delayMask;
		}
		
		if(x->x_feedbackSwitch == 1)
			feedback = ((x->x_fFeedback * 2.0f));
		else
			feedback = -((x->x_fFeedback * 2.0f));
			
		// 3 - calculate left channel
		if(x->x_freeze == 0)
		{
			// our hermite spline function
			// we grab the 2 points before and 2 points after the
			// actual delay point
			// in brackets just for local variable scope
			{
				float inm1 = x->x_delay[(x->x_readPos - 1) & x->x_delayMask];
				float inm0   = x->x_delay[(x->x_readPos + 0) & x->x_delayMask];
				float inp1 = x->x_delay[(x->x_readPos + 1) & x->x_delayMask];
				float inp2 = x->x_delay[(x->x_readPos + 2) & x->x_delayMask];
				output = delay_tilde_hermite(delayFracL, inm1, inm0, inp1, inp2);
			}
			if(x->x_twoHead == 1)
			{
				float inm1 = x->x_delay[(x->x_readPos2 - 1) & x->x_delayMask];
				float inm0   = x->x_delay[(x->x_readPos2 + 0) & x->x_delayMask];
				float inp1 = x->x_delay[(x->x_readPos2 + 1) & x->x_delayMask];
				float inp2 = x->x_delay[(x->x_readPos2 + 2) & x->x_delayMask];
				output2 = delay_tilde_hermite(delayFrac2L, inm1, inm0, inp1, inp2);
				output = (output * (*(x->x_triWaveform + (long)x->x_lfoPhase)+1.0f))
					+ (output2 * (*(x->x_triWaveform + (long)x->x_lfoPhase2)+1.0f));
				output *= 0.5f;
			}

//			output = *(x->x_delay + x->x_readPos);
			preFilter = ((output * feedback) + *(in1+frame));
// here's where we filter
			preClip = A1 * preFilter + A2 * x->X1 + A3 * x->X2 - B1 * x->Y1 - B2 * x->Y2;
			x->X2 = x->X1;
			x->X1 = preFilter;
			x->Y2 = x->Y1;
			x->Y1 = preClip;
// bypass the filter if the frequency is high
			if(x->x_fFilterFreq > 19000)
				preClip = preFilter;

			
			if(feedback > 1.0f)
			{
			dcIn0 = preClip;
			dcOut0 = dcIn0 - x->dcIn1 + dcR * x->dcOut1;
			x->dcOut1 = dcOut0;
			x->dcIn1 = dcIn0;
			preClip = dcOut0;
			}

			// the new arctan softclip			
			*(x->x_delay+x->x_writePos) = atan( preClip * clipFactor )/clipFactor + 0.000001f; 
			*(out+frame) = output;
		}
		else
		{
			// freeze just copies delay out to delay in with no
			// input signal added
			// feedback doesn't matter
			output = *(x->x_delay+x->x_writePos) = *(x->x_delay + x->x_readPos);
			*(out+frame) = output;
		}		
		
		x->x_writePos++;
		x->x_writePos &= x->x_delayMask;

		// interpolate
		x->x_fFeedback += feedbackInc;
		x->x_fModSpeed += speedInc;
		x->x_fModDepth += x->x_depthInc;
		x->x_fEnvelopeModulate += envModInc;
		x->x_fModPhase += lfoPhaseInc;

		tempTime = x->x_fTime + x->x_timeInc;

		// a test to see if tempTime has gone past x->x_time
		if((tempTime - x->x_time) * (x->x_fTime - x->x_time) < 0.0)
		{
			x->x_fTime = x->x_time;
			x->x_timeInc = 0.0f;
		}
		else
			x->x_fTime = tempTime;
	}
		
	x->x_fEnvelopeModulate = x->x_envelopeModulate;
	x->x_fFeedback = x->x_feedback;
	x->x_fModSpeed = x->x_modSpeed;
	
	if( x->x_depthIncTicks==0 )
	{
		x->x_fModDepth = x->x_modDepth;
		x->x_depthInc = 0.0;
	}
	else
		x->x_depthIncTicks--;
		
	x->x_fModPhase = x->x_modPhase;
	
	return(w + 5);
};

static void delay_tilde_dsp(t_delay_tilde *x, t_signal **sp)
{
	dsp_add(
		delay_tilde_perform,
		4,
		x,
		sp[0]->s_vec,
		sp[1]->s_vec,
		sp[0]->s_n
	);

	if(sp[0]->s_n != x->x_n || sp[0]->s_sr != x->x_sr)
	{
		x->x_n = sp[0]->s_n;
		x->x_sr = sp[0]->s_sr;

		x->x_oneOverFrames = 1.0f/x->x_n;
	};
};

void setup_0x2bdelay_tilde(void)
{
	delay_tilde_class = 
	class_new(
		gensym("+delay~"),
		(t_newmethod)delay_tilde_new,
		0,
		sizeof(t_delay_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(delay_tilde_class, t_delay_tilde, x_f);

	class_addmethod(
		delay_tilde_class,
		(t_method)delay_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_time,
		gensym("time"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_feedback,
		gensym("feedback"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_modSpeed,
		gensym("lfoFreq"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_modDepth,
		gensym("lfoDepth"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_modShape,
		gensym("lfoShape"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_lfoPhase,
		gensym("lfoPhase"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_filterFreq,
		gensym("filterFreq"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_resonance,
		gensym("resonance"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_feedbackSwitch,
		gensym("feedbackSwitch"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_freeze,
		gensym("freeze"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		delay_tilde_class, 
        (t_method)delay_tilde_heads,
		gensym("heads"),
		A_DEFFLOAT,
		0
	);
}