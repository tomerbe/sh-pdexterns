/*

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 2:00PM - removed unused variables.

*/

#include <math.h>
#include "m_pd.h"
#define DELAYSIZE 1048576

static t_class *pitchdelay_tilde_class;

typedef struct _pitchdelay_tilde
{
	t_object x_obj;
	float x_sr;
	float x_n;
	float x_delay[DELAYSIZE];
	float x_triWaveform[DELAYSIZE];
	float x_time;
	float x_fTime;
	float x_timeInc;
	float x_pitchFactor[35];
	int x_pitchShift;
	float x_feedback;
	float x_fFeedback;
	float x_octave;
	float x_loopDepth;
	float x_fLoopDepth;
	float x_depthInc;
	int x_depthIncTicks;
	float x_rampOutA;
	float x_rampOutB;
	float x_rampPhaseA;
	float x_rampPhaseB;
	float x_rampIncrement;
	long x_rampTableSize;
	long x_readPosA;
	long x_readPosB;
	long x_writePos;
	long x_delayMask;
	float x_dcIn1;
	float x_dcOut1;
    float x_f;
} t_pitchdelay_tilde;

/***************** UTILITY FUNCTIONS *****************/
// from http://www.musicdsp.org/archive.php?classid=5#93
static float pitchdelay_tilde_hermite(float xx, float yy0, float yy1, float yy2, float yy3)
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


static void pitchdelay_tilde_time(t_pitchdelay_tilde *x, t_floatarg t)
{	
	if(t>5000.0 || t<1.0)
		error("time value must be >= 1.0 ms and <= 5000.0 ms.");
	else
	{
		x->x_time = t/1000.0;
		x->x_timeInc = (x->x_time - x->x_fTime) /(x->x_sr * 0.5f);
	}
}

// need discrete just/eq temp steps
static void pitchdelay_tilde_pitchFactor(t_pitchdelay_tilde *x, t_floatarg pf)
{	
	if( pf>34 || pf<0 )
		error("pitchFactor values must between 0 and 1200 cents.");
	else
		x->x_pitchShift = pf; // 0-34 for indexing into pitchFactor
}

static void pitchdelay_tilde_feedback(t_pitchdelay_tilde *x, t_floatarg f)
{	
	if(f>200.0 || f<0.0)
		error("feedback value must be >= 0%% and <= 200%%.");
	else
		x->x_feedback = f/200.0;
}

static void pitchdelay_tilde_octave(t_pitchdelay_tilde *x, t_floatarg o)
{	
	if(o > 3.0 || o < -3.0)
		error("octave value must be >= -3.0 and <= 3.0.");
	else
		x->x_octave = o;
}

static void pitchdelay_tilde_loopDepth(t_pitchdelay_tilde *x, t_floatarg d)
{	
	float diff;
	
	diff = (x->x_loopDepth - x->x_fLoopDepth);
	
	if(d>100.0 || d<10.0)
		error("lfoDepth value must be >= 10.0%% and <= 100.0%%.");
	else
	{
		x->x_loopDepth = d/100.0;
		x->x_depthInc = diff/x->x_sr;
		x->x_depthIncTicks = (int)x->x_sr/x->x_n;
	}
}

static void *pitchdelay_tilde_new()
{	
	t_pitchdelay_tilde *x = (t_pitchdelay_tilde *)pd_new(pitchdelay_tilde_class);
	int i;

	outlet_new(&x->x_obj, &s_signal);

	x->x_sr = 44100.0;
	x->x_n = 64.0;
	
	for(i=0; i<DELAYSIZE; i++)
		x->x_delay[i] = 0.0;

	for(i = 0; i < DELAYSIZE/2; i++)
		x->x_triWaveform[i] = ((float)i/(float)DELAYSIZE * 4.0) - 1.0f;
	for(; i < DELAYSIZE; i++)
		x->x_triWaveform[i] = 3.0f - ((float)i/(float)DELAYSIZE * 4.0);

	x->x_delayMask = DELAYSIZE-1;

	x->x_time = x->x_fTime = 0.3;
	x->x_timeInc = 0.0;
	
	x->x_pitchFactor[0] = 1.0f;
	x->x_pitchFactor[1] = pow(2.0, 1.0/12.0);	//1.05946
	x->x_pitchFactor[2] = 16.0f/15.0f;			//1.06666
	x->x_pitchFactor[3] = 10.0f/9.0f;			//1.11111
	x->x_pitchFactor[4] = pow(2.0, 2.0/12.0);	//1.12246
	x->x_pitchFactor[5] = 9.0f/8.0f;				//1.125
	x->x_pitchFactor[6] = 8.0f/7.0f;				//1.14285
	x->x_pitchFactor[7] = 7.0f/6.0f;				//1.16666
	x->x_pitchFactor[8] = pow(2.0, 3.0/12.0);	//1.18920
	x->x_pitchFactor[9] = 6.0f/5.0f;				//1.2
	x->x_pitchFactor[10] = 5.0f/4.0f;			//1.25
	x->x_pitchFactor[11] = pow(2.0, 4.0/12.0);	//1.25992
	x->x_pitchFactor[12] = 9.0f/7.0f;			//1.28571
	x->x_pitchFactor[13] = 4.0f/3.0f;			//1.33333
	x->x_pitchFactor[14] = pow(2.0, 5.0/12.0);	//1.33483
	x->x_pitchFactor[15] = 7.0f/5.0f;			//1.4
	x->x_pitchFactor[16] = pow(2.0, 6.0/12.0);	//1.41421
	x->x_pitchFactor[17] = 10.0f/7.0f;			//1.42857
	x->x_pitchFactor[18] = pow(2.0, 7.0/12.0);	//1.49830
	x->x_pitchFactor[19] = 3.0f/2.0f;			//1.5
	x->x_pitchFactor[20] = 14.0f/9.0f;			//1.55555
	x->x_pitchFactor[21] = 11.0f/7.0f;			//1.57142
	x->x_pitchFactor[22] = pow(2.0, 8.0/12.0);	//1.58740
	x->x_pitchFactor[23] = 8.0f/5.0f;			//1.6
	x->x_pitchFactor[24] = 5.0f/3.0f;			//1.66666
	x->x_pitchFactor[25] = pow(2.0, 9.0/12.0);	//1.68179
	x->x_pitchFactor[26] = 12.0f/7.0f;			//1.71428
	x->x_pitchFactor[27] = 7.0f/4.0f;			//1.75
	x->x_pitchFactor[28] = 16.0f/9.0f;			//1.77777
	x->x_pitchFactor[29] = pow(2.0, 10.0/12.0);	//1.78179
	x->x_pitchFactor[30] = 9.0f/5.0f;			//1.8
	x->x_pitchFactor[31] = 13.0f/7.0f;			//1.85714
	x->x_pitchFactor[32] = 15.0f/8.0f;			//1.875
	x->x_pitchFactor[33] = pow(2.0, 11.0/12.0);	//1.88774
	x->x_pitchFactor[34] = 2.0f/1.0f;			//2.0
	
	x->x_pitchShift = 0;
	x->x_feedback = x->x_fFeedback = 0.5;
	x->x_octave = 0;
	x->x_loopDepth = x->x_fLoopDepth = 0.25;
	x->x_depthInc = 0.0;
	x->x_depthIncTicks = 0;
	
	x->x_rampTableSize = DELAYSIZE;
	x->x_rampIncrement = 0.0f;
	x->x_rampPhaseA = 0.0f;
	x->x_rampPhaseB = 0.0f;
	x->x_rampOutA = 0.0f;
	x->x_rampOutB = 0.0f;
	x->x_readPosA = x->x_readPosB = x->x_writePos = 0.0;

	x->x_dcIn1 = x->x_dcOut1 = 0.0;
	
	return(void *)x;
};

static t_int *pitchdelay_tilde_perform(t_int *w)
{
	t_pitchdelay_tilde *x = (t_pitchdelay_tilde *)(w[1]);
	
	t_sample *in1 = (t_sample *)(w[2]);
	t_sample *out = (t_sample *)(w[3]);
	int n = (int)(w[4]);
	int index;

	long frame;
	float outputA, delayTime, delayTimeA, delayTimeB;
	float outputB;
	float preClip;
	
	float dcOut0, dcIn0, dcR;
	float feedbackInc, tempTime;
	float oneOverFrames = 1.0/(float)n;
	float clipFactor = 1.0f;
	float feedback;
	float minDelayMod, maxDelayMod;

	float delayFracA;
	float delayFracB;
	float octave;
	long delayLong;
	
	dcR = 1.0f - (126.0f/x->x_sr);
	// set up for parameter interpolation
	// delay time interpolation set in setParameter
	feedbackInc = (x->x_feedback - x->x_fFeedback) * oneOverFrames;
	
	index = x->x_pitchShift;

	octave = x->x_octave;
	x->x_rampIncrement = (x->x_pitchFactor[index] * pow(2.0, octave)) - 1.0f;
	

	// delay processing loop
	for(frame = 0; frame < n; frame++)
	{
		// 2 - determine the delay time and the readPos
		delayTime = x->x_fTime * x->x_sr;
		
		if(x->x_fLoopDepth < 0.01f)
			x->x_fLoopDepth = 0.01f;
		
		x->x_rampOutA += x->x_rampIncrement;
		x->x_rampOutB = x->x_rampOutA + (delayTime * x->x_fLoopDepth);
		
		// calculate lowest and highest delay
		maxDelayMod = delayTime * x->x_fLoopDepth;
		minDelayMod = -maxDelayMod;
		while(x->x_rampOutA > maxDelayMod)
			x->x_rampOutA -= (2.0f * maxDelayMod);
		while(x->x_rampOutB > maxDelayMod)
			x->x_rampOutB -= (2.0f * maxDelayMod);
		while(x->x_rampOutA <= minDelayMod)
			x->x_rampOutA += (2.0f * maxDelayMod);
		while(x->x_rampOutB <= minDelayMod)
			x->x_rampOutB += (2.0f * maxDelayMod);
		
		x->x_rampPhaseA = x->x_rampTableSize*(x->x_rampOutA + maxDelayMod)/(2.0f*maxDelayMod);
		x->x_rampPhaseB = x->x_rampTableSize*(x->x_rampOutB + maxDelayMod)/(2.0f*maxDelayMod);
		while(x->x_rampPhaseA>=x->x_rampTableSize)
			x->x_rampPhaseA -= x->x_rampTableSize;
		while(x->x_rampPhaseA<0.0)
			x->x_rampPhaseA += x->x_rampTableSize;
		while(x->x_rampPhaseB>=x->x_rampTableSize)
			x->x_rampPhaseB -= x->x_rampTableSize;
		while(x->x_rampPhaseB<0.0)
			x->x_rampPhaseB += x->x_rampTableSize;
		// a two sample buffer for the hermite interpolation
		delayTimeA = delayTime - x->x_rampOutA + 2.0f;
		delayTimeB = delayTime - x->x_rampOutB + 2.0f;
		

		delayLong = (long)delayTimeA;
		delayFracA = delayTimeA - (float)delayLong;
		delayFracA = 1.0 - delayFracA;
		x->x_readPosA = x->x_writePos - delayLong;
		x->x_readPosA &= x->x_delayMask;


		// MIDI clock locking code for the delay time

		delayLong = (long)delayTimeB;
		delayFracB = delayTimeB - (float)delayLong;
		delayFracB = 1.0 - delayFracB;
		x->x_readPosB = x->x_writePos - delayLong;
		x->x_readPosB &= x->x_delayMask;
			
		feedback = ((x->x_fFeedback * 2.0f));
			
			// 3 calculate output
			// our hermite spline function
			// we grab the 2 points before and 2 points after the
			// actual delay point
			// in brackets just for local variable scope
			{
				float inm1 = x->x_delay[(x->x_readPosA - 1) & x->x_delayMask];
				float in   = x->x_delay[(x->x_readPosA + 0) & x->x_delayMask];
				float inp1 = x->x_delay[(x->x_readPosA + 1) & x->x_delayMask];
				float inp2 = x->x_delay[(x->x_readPosA + 2) & x->x_delayMask];
				outputA = pitchdelay_tilde_hermite(delayFracA, inm1, in, inp1, inp2);
			}
			{
				float inm1 = x->x_delay[(x->x_readPosB - 1) & x->x_delayMask];
				float in   = x->x_delay[(x->x_readPosB + 0) & x->x_delayMask];
				float inp1 = x->x_delay[(x->x_readPosB + 1) & x->x_delayMask];
				float inp2 = x->x_delay[(x->x_readPosB + 2) & x->x_delayMask];
				outputB = pitchdelay_tilde_hermite(delayFracB, inm1, in, inp1, inp2);
				outputA = (outputA * (*(x->x_triWaveform + (long)x->x_rampPhaseA)+1.0f))
					+ (outputB * (*(x->x_triWaveform + (long)x->x_rampPhaseB)+1.0f));
				outputA *= 0.5f;
			}
			preClip = ((outputA * feedback) + *(in1+frame));

		dcIn0 = preClip;
		dcOut0 = dcIn0 - x->x_dcIn1 + dcR * x->x_dcOut1;
		x->x_dcOut1 = dcOut0;
		x->x_dcIn1 = dcIn0;
		preClip = dcOut0;

		// the new arctan softclip			
		*(x->x_delay+x->x_writePos) = atan( preClip * clipFactor )/clipFactor + 0.000001f; 
		*(out+frame) = outputA;
				
		x->x_writePos++;
		x->x_writePos &= x->x_delayMask;
		
		// interpolate
		x->x_fFeedback += feedbackInc;
		x->x_fLoopDepth += x->x_depthInc;
		
		tempTime = x->x_fTime + x->x_timeInc;

		// a test to see if tempTime has gone past pTime
		if((tempTime - x->x_time) * (x->x_fTime - x->x_time) < 0.0)
		{
			x->x_fTime = x->x_time;
			x->x_timeInc = 0.0f;
		}
		else
			x->x_fTime = tempTime;
	
	};

	if( x->x_depthIncTicks==0 )
	{
		x->x_fLoopDepth = x->x_loopDepth;
		x->x_depthInc = 0.0;
	}
	else
		x->x_depthIncTicks--;
		
	x->x_fFeedback = x->x_feedback;
	
	return(w + 5);
};

static void pitchdelay_tilde_dsp(t_pitchdelay_tilde *x, t_signal **sp)
{
	dsp_add(
		pitchdelay_tilde_perform,
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
	};
};

void setup_0x2bpitchdelay_tilde(void)
{
	pitchdelay_tilde_class = 
	class_new(
		gensym("+pitchdelay~"),
		(t_newmethod)pitchdelay_tilde_new,
		0,
		sizeof(t_pitchdelay_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(pitchdelay_tilde_class, t_pitchdelay_tilde, x_f);

	class_addmethod(
		pitchdelay_tilde_class,
		(t_method)pitchdelay_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		pitchdelay_tilde_class, 
        (t_method)pitchdelay_tilde_time,
		gensym("time"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		pitchdelay_tilde_class, 
        (t_method)pitchdelay_tilde_pitchFactor,
		gensym("pitchFactor"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		pitchdelay_tilde_class, 
        (t_method)pitchdelay_tilde_feedback,
		gensym("feedback"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		pitchdelay_tilde_class, 
        (t_method)pitchdelay_tilde_octave,
		gensym("octave"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		pitchdelay_tilde_class, 
        (t_method)pitchdelay_tilde_loopDepth,
		gensym("loopDepth"),
		A_DEFFLOAT,
		0
	);
}