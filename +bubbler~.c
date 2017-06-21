/*

June 20, 2010: added a second inlet for MIDI note on pairs, and internal logic for which notes are active.

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 2:10PM - changing only the longs that are necessary: x_randSeed, and randSeed in the randomUlong function. Removed unused variables and the unnecessary if statements that were at line 541.

*/

#include <stdint.h> // this won't work on windows...
#include <math.h>
#include "m_pd.h"
#define DELAYSIZE 2097152
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

static t_class *bubbler_tilde_class;

typedef struct singleGrain
{
	float start;
	float length;
	float end;
	float attackend;
	float decaystart;
	float divisor;
	float delay;
	float readPos;
	float increment;
	int active;
} t_singleGrain;

typedef struct _bubbler_tilde
{
	t_object x_obj;
	
	uint32_t x_randSeed;
	float x_sr;
	float x_n;	
	
	t_singleGrain x_grain[8];
	
	float x_delay[DELAYSIZE];
	long x_delayMask;
	long x_delaySize;
	
	long x_readPos;
	long x_writePos;
	long x_currentSample;
	long x_lastGrainSample;

	float x_time;
	float x_fTime;
	float x_timeInc;
	float x_timeVariation;
	float x_fTimeVariation;
	
	float x_feedback;
	float x_fFeedback;
	float x_filterFreq;
	float x_fFilterFreq;
	float x_resonance;
	float x_fResonance;
	
	float x_density;
	float x_fDensity;
	
	float x_grainSpacing;
	float x_grainSize;
	float x_fGrain_size;
	float x_grainSizeVariation;
	float x_fGrainSizeVariation;
	float x_grainStartVariation;
	float x_fGrainStartVariation;
	float x_grainReversal;
	float x_fGrainReversal;
	
	float x_octave;
	float x_fOctave;
	float x_octaveVariation;
	float x_fOctaveVariation;
	
	int x_just;
	
	long x_numNotes;
	long x_pitchOn[12];
	float x_pitch[12];
	float x_scale[2][13];

	float X1;
	float X2;
	float Y1;
	float Y2;
	
	float dcOut1;
	float dcIn1;
    
    float x_f;
    
} t_bubbler_tilde;


/***************** UTILITY FUNCTIONS *****************/
static void bubbler_tilde_randomULong(uint32_t *randSeed)
{
   /* Change this for different random sequences. */
   *randSeed = (*randSeed * 196314165) + 907633515;
};

// from http://www.musicdsp.org/archive.php?classid=5#93
static float bubbler_tilde_hermite(float xx, float yy0, float yy1, float yy2, float yy3)
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


static void bubbler_tilde_time(t_bubbler_tilde *x, t_floatarg t)
{	
	if(t>10000.0 || t<1.0)
		error("time value must be >= 1.0 ms and <= 10000.0 ms.");
	else
	{
		x->x_time = t/1000.0; //convert to seconds
		x->x_timeInc = (x->x_time - x->x_fTime)/(x->x_sr * 0.1f);
	}
}

static void bubbler_tilde_time_vari(t_bubbler_tilde *x, t_floatarg tv)
{	
	if(tv>100.0 || tv<0.0)
		error("timeVariation value must be >= 0%% and <= 100%%.");
	else
		x->x_timeVariation = tv/100;
}

static void bubbler_tilde_feedback(t_bubbler_tilde *x, t_floatarg f)
{	
	if(f>200.0 || f<0.0)
		error("feedback value must be >= 0%% and <= 200%%.");
	else
		x->x_feedback = f/200;
}

static void bubbler_tilde_filter_freq(t_bubbler_tilde *x, t_floatarg ff)
{	
	if(ff>20000.0 || ff<20.0)
		error("filterFreq value must be >= 20 Hz and <= 20000 Hz.");
	else
		x->x_filterFreq = ff;
}

static void bubbler_tilde_resonance(t_bubbler_tilde *x, t_floatarg r)
{	
	if(r>100 || r<0)
		error("resonance value must be >= 0%% and <= 100%%.");
	else
		x->x_resonance = r/100.0;
}

static void bubbler_tilde_density(t_bubbler_tilde *x, t_floatarg d)
{	
	if(d>200 || d<0)
		error("density value must be >= 0%% and <= 200%%.");
	else
		x->x_density = d/200.0;
}

static void bubbler_tilde_grain_size(t_bubbler_tilde *x, t_floatarg g)
{	
	if(g>50 || g<0)
		error("grainSize value must be >= 0%% and <= 50%%.");
	else
		x->x_grainSize = g/50.0;
}

static void bubbler_tilde_grain_start_vari(t_bubbler_tilde *x, t_floatarg gs)
{	
	if(gs>100 || gs<0)
		error("grainStartVariation value must be >= 0%% and <= 100%%.");
	else
		x->x_grainStartVariation = gs/100.0;
}

static void bubbler_tilde_octave(t_bubbler_tilde *x, t_floatarg o)
{	
	if(o > 4.0 || o < -4.0)
		error("octave value must be >= -4.0 and <= 4.0.");
	else
		x->x_octave = o;
}

static void bubbler_tilde_oct_vari(t_bubbler_tilde *x, t_floatarg ov)
{	
	if(ov>2.0 || ov<0.0)
		error("octaveVariation value must be >= 0.0 and <= 2.0.");
	else
		x->x_octaveVariation = ov;
}

static void bubbler_tilde_grain_reversal(t_bubbler_tilde *x, t_floatarg gr)
{	
	if(gr>100 || gr<0)
		error("grainReversal value must be >= 0%% and <= 100%%.");
	else
		x->x_grainReversal = gr/100.0;
}

static void bubbler_tilde_12tet(t_bubbler_tilde *x)
{	
	x->x_just = 0;
	post("12tet intonation.");
}

static void bubbler_tilde_just(t_bubbler_tilde *x)
{	
	x->x_just = 1;
	post("just intonation.");
}

static void bubbler_tilde_pitchOn(t_bubbler_tilde *x, t_float note, t_float vel)
{	
	int i, mod_note;
	t_float octave;

	if(note<0)
	{
		error("note values must be positive");
		return;
	}
	else
		mod_note = (int)note%12;
	
	if(vel<0)
		vel=0.0;

	if(x->x_pitch[mod_note] == 0.0f && vel>0)
	{	
		octave = (int)(note/12) - 5;
		x->x_octave = octave;
		
		//add a note
		x->x_pitch[mod_note] = 1.0f;
		
		// safety: clip numNotes to 0-11 for pitchOn
		if(x->x_numNotes < 0)
			x->x_numNotes = 0;
		else if(x->x_numNotes > 11)
			x->x_numNotes = 11;
			
		x->x_pitchOn[x->x_numNotes] = mod_note;
		x->x_numNotes++;
	}
	else if(x->x_pitch[mod_note] == 1.0f && vel==0)
	{		
		// delete a note
		x->x_pitch[mod_note] = 0.0f;
		for(i=0; i<12; i++)
		{	
			if(x->x_pitchOn[i] != mod_note)
				;
			else
			{
				for(; i<11; i++)
					x->x_pitchOn[i] = x->x_pitchOn[i+1];
				
				x->x_pitchOn[i] = 0.0;
				break;
			}
		}
		
		x->x_numNotes--;
	}
}

static void *bubbler_tilde_new()
{
	t_bubbler_tilde *x = (t_bubbler_tilde *)pd_new(bubbler_tilde_class);
	unsigned long i;

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("list"), gensym("pitchOn"));
	outlet_new(&x->x_obj, &s_signal);

	bubbler_tilde_randomULong(&x->x_randSeed);
	x->x_sr = 44100.0;
	x->x_n = 64.0;
	
	for(i=0; i<8; i++)
	{
		x->x_grain[i].start = 0;
		x->x_grain[i].length = 0;
		x->x_grain[i].end = 0;
		x->x_grain[i].attackend = 0;
		x->x_grain[i].decaystart = 0;
		x->x_grain[i].divisor = 0;
		x->x_grain[i].delay = 0;
		x->x_grain[i].readPos = 0;
		x->x_grain[i].increment = 0;
		x->x_grain[i].active = 0;
	};

	for(i=0; i<DELAYSIZE; i++)
		x->x_delay[i] = 0.0;
		
	x->x_delaySize = DELAYSIZE;
	x->x_delayMask = x->x_delaySize - 1;
	
	x->x_readPos = x->x_writePos = 0.0;
	x->x_currentSample = x->x_lastGrainSample = 0.0;

	x->x_time = x->x_fTime = 8.32;
	x->x_timeInc = (x->x_time - x->x_fTime) /(x->x_sr * 0.1f);
	x->x_timeVariation = x->x_fTimeVariation = 0.0;
	
	x->x_feedback = x->x_fFeedback = 0.0;
	x->x_filterFreq = x->x_fFilterFreq = 20000.0;
	x->x_resonance = x->x_fResonance = 0.5;
	
	x->x_density = x->x_fDensity = 0.8;
	x->x_grainSpacing = 0.0;
	x->x_grainSize = x->x_fGrain_size = 0.2;
	x->x_grainSizeVariation = x->x_fGrainSizeVariation = 0.0;
	x->x_grainStartVariation = x->x_fGrainStartVariation = 0.0;
	x->x_grainReversal = x->x_fGrainReversal = 0.0;
	
	x->x_octave = x->x_fOctave = 0.0;
	x->x_octaveVariation = x->x_fOctaveVariation = 0.0;
	x->x_just = 0; // 0 is 12tet, 1 is just
	
	for(i=0; i<12; i++)
		x->x_pitch[i]=0.0;

	for(i=0; i<12; i++)
		x->x_pitchOn[i]=0.0;
		
		
	for(i=0; i<13; i++)
		x->x_scale[0][i] = pow(2.0, (float)i/12.0);
		
	x->x_scale[1][0] = 1.0;
	x->x_scale[1][1] = 16.0/15.0;
	x->x_scale[1][2] = 9.0/8.0;
	x->x_scale[1][3] = 6.0/5.0;
	x->x_scale[1][4] = 5.0/4.0;
	x->x_scale[1][5] = 4.0/3.0;
	x->x_scale[1][6] = 7.0/5.0;
	x->x_scale[1][7] = 3.0/2.0;
	x->x_scale[1][8] = 8.0/5.0;
	x->x_scale[1][9] = 5.0/3.0;
	x->x_scale[1][10] = 7.0/4.0;
	x->x_scale[1][11] = 15.0/8.0;
	x->x_scale[1][12] = 2.0;
	
	// biquad filter stuff
	x->X1 = x->X2 = x->Y1 = x->Y2 = 0.0;
	x->dcOut1 = x->dcIn1 = 0.0f;
	
	return(void *)x;
};

static t_int *bubbler_tilde_perform(t_int *w)
{
	t_bubbler_tilde *x = (t_bubbler_tilde *)(w[1]);
	
	t_sample *mainInput = (t_sample *)(w[2]);
	t_sample *out = (t_sample *)(w[3]);
	int n = (int)(w[4]);

	long frame, frames;
	long i;
	float output, delayTime;
	float preClip, preFilter;
	float grainSizeVariation, timeVariation, grainSamples;
	float octaveVariation, octave, reversal;
	long pitchOffset;
	float twoPi;

	// feedback filter variables
	float filterFreq, filterReson, f0, C;
	float A1, A2, A3, B1, B2;
	
	float feedbackInc, grainSizeInc, grainSizeVariInc, grainStartVariInc, tempTime, timeVariInc, densityInc, octaveInc, octaveVariationInc, reverseGrainInc;
	float oneOverFrames;
	float clipFactor = 1.0f;
	float feedback;

	float grainOut;
	float floatNote;
	
	float dcOut0, dcIn0, dcR;
	
	float delayFrac;
	long delayLong;

    frames = (long)n;
	oneOverFrames = 1.0/(float)frames;
	
	// set up for parameter interpolation
	// graindelay time interpolation set in setParameter
	feedbackInc = (x->x_feedback - x->x_fFeedback) * oneOverFrames;
	timeVariInc = (x->x_timeVariation - x->x_fTimeVariation) * oneOverFrames;
	grainSizeInc = (x->x_grainSize - x->x_fGrain_size) * oneOverFrames;
	grainSizeVariInc = (x->x_grainSizeVariation - x->x_fGrainSizeVariation) * oneOverFrames;
	grainStartVariInc = (x->x_grainStartVariation - x->x_fGrainStartVariation) * oneOverFrames;
	densityInc = (x->x_density - x->x_fDensity) * oneOverFrames;
	octaveInc = (x->x_octave - x->x_fOctave) * oneOverFrames;
	octaveVariationInc = (x->x_octaveVariation - x->x_fOctaveVariation) * oneOverFrames;
	reverseGrainInc = (x->x_grainReversal - x->x_fGrainReversal) * oneOverFrames;
	
	twoPi = 8.0 * atan(1.0);
	x->x_fResonance = x->x_resonance;
	x->x_fFilterFreq = x->x_filterFreq;
	
	dcR = 1.0f - (126.0f/x->x_sr);
	// filter coefficients get calculated only once per block
	filterFreq = x->x_fFilterFreq;
// 	filterFreq = (9.965784 * x->x_fFilterFreq) + 4.321928;
// 	filterFreq = pow(2.0, filterFreq);
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

	// graindelay processing loop
	for(frame = 0; frame < frames; frame++, x->x_currentSample++)
	{		
		// 2 - determine the graindelay time and the readPos
		delayTime = (x->x_fTime * x->x_sr)  + 2.0;
		
		delayLong = (long)delayTime;
		delayFrac = delayTime - (float)delayLong;
		delayFrac = 1.0 - delayFrac;
		
		// 0 - do we create a new grain
		// chances of grain creation
		// 1 grain every duration(samples)/density samples
		// grainsSample can vary from 0.0 to 1.0
		if (x->x_fGrain_size <= 0.0)
			x->x_fGrain_size = 0.001;
		if (x->x_fDensity <= 0.0)
			x->x_fDensity = 0.00;
		grainSamples = x->x_fGrain_size * 0.5f * delayTime;
		if(x->x_fDensity > 0.0f)
			x->x_grainSpacing = grainSamples/(8.0f * x->x_fDensity);
		else 
			x->x_grainSpacing = 1000000000000000000000000000000000.0f;

		if(x->x_lastGrainSample + x->x_grainSpacing < x->x_currentSample)
		{
			for(i = 0; i < 8; i++)
			{
				if(x->x_grain[i].active == 0)
				{
					long pitch;
					float freqRatio;
					x->x_lastGrainSample = x->x_currentSample;
					x->x_grain[i].active = 1;
					bubbler_tilde_randomULong(&x->x_randSeed);
					x->x_grain[i].start = x->x_currentSample + ((float)x->x_randSeed/4294967296.0) * x->x_grainSpacing * x->x_fGrainStartVariation;
					bubbler_tilde_randomULong(&x->x_randSeed);
					grainSizeVariation = ((((float)x->x_randSeed/4294967296.0) * 2.0) - 1.0) * x->x_fGrainSizeVariation;
					x->x_grain[i].length = grainSamples + (grainSamples * grainSizeVariation);
					x->x_grain[i].end = x->x_grain[i].start + x->x_grain[i].length;
					if(x->x_grain[i].length > 0)
					{
						if(x->x_grain[i].length/x->x_sr > 0.1)
						{
							x->x_grain[i].attackend = x->x_grain[i].start + (0.05 * x->x_sr);
							x->x_grain[i].decaystart = x->x_grain[i].end - (0.05 * x->x_sr);
							x->x_grain[i].divisor = 1.0f/(0.1 * x->x_sr);
						}
						else
						{
							x->x_grain[i].divisor = 1.0f/x->x_grain[i].length;
							x->x_grain[i].attackend = x->x_grain[i].decaystart = 0.0;
						}
					}
					else
						x->x_grain[i].divisor = 1.0f;

					// how much pitch shift
					bubbler_tilde_randomULong(&x->x_randSeed);
					floatNote = ((float)x->x_randSeed/4294967296.0) * (float)x->x_numNotes;
					if(floatNote >= (float)x->x_numNotes)
						floatNote = (float)x->x_numNotes - 1.0f;
					bubbler_tilde_randomULong(&x->x_randSeed);
					octaveVariation = ((((float)x->x_randSeed/4294967296.0) * 2.0) - 1.0) * x->x_fOctaveVariation;
					octave = x->x_octave;
					pitchOffset = (long)(octave + octaveVariation);
					pitch = x->x_pitchOn[(long)floatNote];
					
					// safety: clip pitch to 0-12 range for x_scale
					if(pitch < 0)
						pitch = 0;
					else if(pitch > 12)
						pitch = 12;
						
					freqRatio = x->x_scale[x->x_just][pitch];
					if(x->x_numNotes == 0)
						x->x_grain[i].increment = powf(2.0, (float)pitchOffset);
					else
						x->x_grain[i].increment = freqRatio * powf(2.0, (float)pitchOffset);

					// how much delay for this grain
					bubbler_tilde_randomULong(&x->x_randSeed);
					timeVariation = ((((float)x->x_randSeed/4294967296.0) * 2.0) - 1.0) * x->x_fTimeVariation;
					x->x_grain[i].delay = (delayTime * x->x_grain[i].increment) + (delayTime * timeVariation);
					x->x_grain[i].readPos = x->x_writePos - x->x_grain[i].delay;
					// are we reversed or not
					bubbler_tilde_randomULong(&x->x_randSeed);
					reversal = (((float)x->x_randSeed/4294967296.0));
					if(reversal < x->x_fGrainReversal)
						x->x_grain[i].increment = x->x_grain[i].increment * -1.0;
					
					break;
				}
			}
		}

		output = 0.0f;
		for(i=0; i<8; i++)
		{
			if((x->x_grain[i].active == 1) && (x->x_currentSample >= x->x_grain[i].start))
			{
				x->x_readPos = (long)x->x_grain[i].readPos;
				delayFrac = x->x_grain[i].readPos - (float)x->x_readPos;
				delayFrac = 1.0 - delayFrac;
				x->x_readPos &= x->x_delayMask;
				{
					float inm1, in, inp1, inp2;
					unsigned long inm1Ind, inInd, inp1Ind, inp2Ind;

					inm1Ind = (x->x_readPos - 1) & x->x_delayMask;
					inInd   = (x->x_readPos + 0) & x->x_delayMask;
					inp1Ind = (x->x_readPos + 1) & x->x_delayMask;
					inp2Ind = (x->x_readPos + 2) & x->x_delayMask;
					
					if(inm1Ind >= DELAYSIZE)
					{
						inm1Ind = DELAYSIZE-1;
						post("over");
					}
					else if(inInd >= DELAYSIZE)
					{
						inInd = DELAYSIZE-1;
						post("over");
					}
					else if(inp1Ind >= DELAYSIZE)
					{
						inp1Ind = DELAYSIZE-1;
						post("over");
					}
					else if(inp2Ind >= DELAYSIZE)
					{
						inp2Ind = DELAYSIZE-1;
						post("over");
					};
					
					inm1 = x->x_delay[inm1Ind];
					in   = x->x_delay[inInd];
					inp1 = x->x_delay[inp1Ind];
					inp2 = x->x_delay[inp2Ind];
					
					grainOut = bubbler_tilde_hermite(delayFrac, inm1, in, inp1, inp2);
				}
				if(x->x_grain[i].attackend == x->x_grain[i].decaystart)
					grainOut *= (0.5 - 0.5*cos(twoPi*(x->x_currentSample - x->x_grain[i].start)*x->x_grain[i].divisor));
				else if(x->x_currentSample <= x->x_grain[i].attackend)
					grainOut *= (0.5 - 0.5*cos(twoPi*(x->x_currentSample - x->x_grain[i].start)*x->x_grain[i].divisor));
				else if(x->x_currentSample >= x->x_grain[i].decaystart)
					grainOut *= (0.5 - 0.5*cos(twoPi*(x->x_currentSample - (x->x_grain[i].decaystart+(0.05*x->x_sr)))*x->x_grain[i].divisor));
				output += grainOut * 0.5;

				// in vst this is only if(kMonoMode)
				x->x_grain[i].readPos += x->x_grain[i].increment;					
				while(x->x_grain[i].readPos > (float)x->x_delaySize)
					x->x_grain[i].readPos -= (float)x->x_delaySize;
				while(x->x_grain[i].readPos < 0.0)
					x->x_grain[i].readPos += (float)x->x_delaySize;
			}
		}
		
		feedback = x->x_fFeedback * 2.0f;
		preFilter = ((output * feedback) + *(mainInput+frame));
// here's where we filter
		preClip = A1 * preFilter + A2 * x->X1 + A3 * x->X2 - B1 * x->Y1 - B2 * x->Y2;
		x->X2 = x->X1;
		x->X1 = preFilter;
		x->Y2 = x->Y1;
		x->Y1 = preClip;
// bypass the filter if the frequency is high
		if(x->x_fFilterFreq > 19000) // used to be 0.95 for 0-1 range
			preClip = preFilter;
		
		dcIn0 = preClip;
		dcOut0 = dcIn0 - x->dcIn1 + dcR * x->dcOut1;
		x->dcOut1 = dcOut0;
		x->dcIn1 = dcIn0;
		preClip = dcOut0;
		
// the new arctan softclip			
		*(x->x_delay+x->x_writePos) = atan( preClip * clipFactor )/clipFactor + 0.000001f; 
		*(out+frame) = output;
			
		x->x_writePos++;
		x->x_writePos &= x->x_delayMask;
		
		// interpolate
		x->x_fFeedback += feedbackInc;

		tempTime = x->x_fTime + x->x_timeInc;
		x->x_fTimeVariation += timeVariInc;
		x->x_fGrain_size += grainSizeInc;
		x->x_fGrainSizeVariation += grainSizeVariInc;
		x->x_fGrainStartVariation += grainStartVariInc;
		x->x_fDensity += densityInc;
		x->x_fOctave += octaveInc;
		x->x_fOctaveVariation += octaveVariationInc;
		x->x_fGrainReversal += reverseGrainInc;
		for(i = 0; i<8; i++)
		{
			if(x->x_currentSample > x->x_grain[i].end)
				x->x_grain[i].active = 0;
		}
	
		// a test to see if tempTime has gone past pTime
		if((tempTime - x->x_time) * (x->x_fTime - x->x_time) < 0.0)
		{
			x->x_fTime = x->x_time;
			x->x_timeInc = 0.0f;
		}
		else
			x->x_fTime = tempTime;
	};
	
	x->x_fFeedback = x->x_feedback;
	x->x_fTimeVariation = x->x_timeVariation;
	x->x_fGrain_size = x->x_grainSize;
	x->x_fGrainSizeVariation = x->x_grainSizeVariation;
	x->x_fGrainStartVariation = x->x_grainStartVariation;
	x->x_fDensity = x->x_density;
	x->x_fOctave = x->x_octave;
	x->x_fOctaveVariation = x->x_octaveVariation;
	x->x_fGrainReversal = x->x_grainReversal;
	
	
	return(w + 5);
};

static void bubbler_tilde_dsp(t_bubbler_tilde *x, t_signal **sp)
{
	dsp_add(
		bubbler_tilde_perform,
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

void setup_0x2bbubbler_tilde(void)
{
	bubbler_tilde_class = 
	class_new(
		gensym("+bubbler~"),
		(t_newmethod)bubbler_tilde_new,
		0,
		sizeof(t_bubbler_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(bubbler_tilde_class, t_bubbler_tilde, x_f);

	class_addmethod(
		bubbler_tilde_class,
		(t_method)bubbler_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_time,
		gensym("time"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_time_vari,
		gensym("timeVariation"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_feedback,
		gensym("feedback"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_filter_freq,
		gensym("filterFreq"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_resonance,
		gensym("resonance"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_density,
		gensym("density"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_grain_size,
		gensym("grainSize"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_grain_start_vari,
		gensym("grainStartVariation"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_octave,
		gensym("octave"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_oct_vari,
		gensym("octaveVariation"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_grain_reversal,
		gensym("grainReversal"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_12tet,
		gensym("12tet"),
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_just,
		gensym("just"),
		0
	);
	
	class_addmethod(
		bubbler_tilde_class, 
        (t_method)bubbler_tilde_pitchOn,
		gensym("pitchOn"),
		A_DEFFLOAT,
		A_DEFFLOAT,
		0
	);
}