/*

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 12:50PM - changed to int32_t instead of long for meansquare and levelindex, even though it didn't seem to crash with regular longs. Also included stdint.h to enable. Removed unused variables.

*/

#include <stdint.h>
#include <math.h>
#include "m_pd.h"

static t_class *compand_tilde_class;

typedef struct _compand_tilde
{
	t_object x_obj;
	float x_sr;
	float x_n;
	float x_compGainTable[8192];
	float x_compSquareArray[64];
	float x_compRoot[65536];
	long x_compRMSArraySize;
	float x_compRMSScale;
	float x_compRMSArrayTotal;
	long x_compRMSPosition;
	float x_compSoftLowThreshold;
	float x_compSoftHighThreshold;
	float x_compDBSoftLowThreshold;
	float x_compReleaseLevel;
	float x_compAttackLevel;
	float x_compAttackInc;
	float x_compAttackMult;
	float x_compReleaseMult;
	float x_gain;
	float x_ratio;
	int x_RMSLevelDetect;
	float x_threshold;
	float x_attack;
	float x_release;
	float x_softKnee;
	float x_makeupGain;
    float x_f;
} t_compand_tilde;

/*********************** UTILITY FUNCTIONS ***********************/
static void compand_tilde_compInitGainTable(t_compand_tilde *x)
{
    long i;
    float insignal, indecibel, gaindecibel;

    if(x->x_ratio <= 1.0f)
        // 0 is -96dB or 20 * log10((0 + 1)/65536.0)
        for(i = 0; i < 8192; i++)
        {
            insignal = (i + 1) *.0001220703125f; // i+1/8192
            indecibel = 20.0f * log10(insignal);
            
            if(insignal <= x->x_compSoftLowThreshold)
                x->x_compGainTable[i] = 1.0f;
            else if(insignal > x->x_compSoftLowThreshold && insignal <= x->x_compSoftHighThreshold)
            {
                gaindecibel = ((x->x_threshold - indecibel) 
                    * (1.0f - x->x_ratio)) * ((indecibel - x->x_compDBSoftLowThreshold)/(x->x_softKnee * 12.0f));
                x->x_compGainTable[i] = pow(10.0f, gaindecibel * 0.05f);
            }
            else
            {
                gaindecibel = ((x->x_threshold - indecibel) * (1.0f - x->x_ratio));
                x->x_compGainTable[i] = pow(10.0f, gaindecibel * 0.05f);
            }
        }
    else
        // 0 is -96dB or 20 * log10((0 + 1)/65536.0)
        for(i = 0; i < 8192; i++)
        {
            insignal = (i + 1)/8192.0f;
            indecibel = 20.0f * log10(insignal);
            if(insignal > x->x_compSoftHighThreshold)
                x->x_compGainTable[i] = 1.0f;
            else if(insignal > x->x_compSoftLowThreshold && insignal <= x->x_compSoftHighThreshold)
            {
                gaindecibel = ((indecibel - x->x_threshold) 
                    * (x->x_ratio - 1.0f)) * ((x->x_compSoftHighThreshold - indecibel)/(x->x_softKnee * 12.0f));
                x->x_compGainTable[i] = pow(10.0f, gaindecibel * 0.05f);
            }
            else
            {
                gaindecibel = ((indecibel - x->x_threshold) * (x->x_ratio - 1.0f));
                x->x_compGainTable[i] = pow(10.0f, gaindecibel * 0.05f);
            }
        }
};


static float compand_tilde_compPeakDetect(float in)
{
	if(in < 0)
		return(-in);
	else
		return(in);
}


// we will try 64 samples for rms length
// too short and we will have low frequency distortion
// too long and we will miss peaks
static float compand_tilde_compRMSDetect(t_compand_tilde *x, float in)
{
    int32_t meansquare;
    float square;
    if(x->x_compRMSPosition >= x->x_compRMSArraySize)
            x->x_compRMSPosition = 0;
	if(x->x_compRMSPosition < 0)
            x->x_compRMSPosition = 0;
    
    // FIRST: square the new sample
    // in range -1.0 to 1.0
    // square range 0.0 to 1.0
    square = in * in;
    
    // SECOND: add to array to calculate mean
    // arraytotal range 0.0 to comprmsarrarysize
    x->x_compRMSArrayTotal = x->x_compRMSArrayTotal - x->x_compSquareArray[x->x_compRMSPosition] + square;
    x->x_compSquareArray[x->x_compRMSPosition] = square;
    x->x_compRMSPosition++;
    
    // THIRD: meansquare needs to be scaled up to 65535 for root table lookup
    meansquare = (int32_t)(x->x_compRMSArrayTotal * x->x_compRMSScale);
	if(meansquare > 65535)
		meansquare = 65535;
    // FOURTH: return the root of the mean of squares
    return(x->x_compRoot[meansquare]);
}
/********************* UTILITY FUNCTIONS END *********************/


static void compand_tilde_ratio(t_compand_tilde *x, t_floatarg value)
{	
	if(value > 10.0)
	{
	// expansion
		error("ratio value must be >= 0.01 and <= 10.0");
		value = 10.0;
	}
	else if(value < 0.01)
	{
	// compression
		error("ratio value must be >= 0.01 and <= 10.0");
		value = 0.01;
	};
		
	x->x_ratio = value;				
	compand_tilde_compInitGainTable(x);
}

static void compand_tilde_threshold(t_compand_tilde *x, t_floatarg threshold)
{
	float compLinearThreshold;
	
	
	if(threshold > 0.0 || threshold < -96.0)
		error("threshold value must be >= -96.0dB and <= 0.0dB.");
	else
		x->x_threshold = threshold;


    compLinearThreshold = pow(10.0f, (x->x_threshold * 0.05f));
    
    if(x->x_softKnee == 0.0f)
    {
        x->x_compSoftLowThreshold = x->x_compSoftHighThreshold = compLinearThreshold;
    }
    else
    {
        // softness can go from -6 to +6 around threshold
        x->x_compSoftHighThreshold = pow(10.0f, ((x->x_threshold + (x->x_softKnee * 6.0f)) * 0.05f));
        x->x_compSoftLowThreshold = pow(10.0f, (((20.0f * log10(x->x_compSoftHighThreshold)) - (x->x_softKnee * 12.0f)) * 0.05f));
    }
    x->x_threshold = 20.0f * log10(compLinearThreshold);

    x->x_compDBSoftLowThreshold = 20.0f * log10(x->x_compSoftLowThreshold);

    compand_tilde_compInitGainTable(x);
}

static void compand_tilde_attack(t_compand_tilde *x, t_floatarg attack)
{	
	if(attack > 100.0 || attack < 1.0)
		error("attack value must be >= 1.0ms and <= 100.0ms.");
	else
		x->x_attack = attack/1000;

	x->x_compAttackInc = 1.0f/(x->x_attack * x->x_sr);
	x->x_compAttackMult = pow(10.0f, (1.0f/(x->x_attack * x->x_sr)));
}

static void compand_tilde_release(t_compand_tilde *x, t_floatarg release)
{	
	if(release > 2000.0 || release < 100.0)
		error("release value must be >= 100ms and <= 2000.0ms.");
	else
		x->x_release = release/1000;

	// release is also in ms, NOT seconds.
	x->x_compReleaseMult = pow(0.1f, (1.0f/(x->x_release * x->x_sr)));
}

static void compand_tilde_softKnee(t_compand_tilde *x, t_floatarg softKnee)
{	
	if(softKnee > 1.0 || softKnee < 0.0)
		error("softKnee value must be >= 0.0 and <= 1.0.");
	else
		x->x_softKnee = softKnee;

	compand_tilde_threshold(x, x->x_threshold);
}

static void compand_tilde_makeupGain(t_compand_tilde *x, t_floatarg makeupGain)
{	
	if(makeupGain > 60.0 || makeupGain < 0.0)
		error("makeupGain value must be >= 0.0dB and <= 60.0dB.");
	else
		x->x_makeupGain = pow(10.0f, makeupGain*0.05f);
}

static void compand_tilde_rms_detect(t_compand_tilde *x, t_floatarg rms)
{	
	if(rms > 1.0)
		rms=1.0;
	else if(rms < 0.0)
		rms=0.0;
	else
		x->x_RMSLevelDetect = rms;
	
	if(x->x_RMSLevelDetect==1)
		post("RMS Level Detect Mode.");
	else
		post("Peak Level Detect Mode.");
}

static void *compand_tilde_new()
{	
	t_compand_tilde *x = (t_compand_tilde *)pd_new(compand_tilde_class);
	int i;

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("ratio"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("threshold"));
	outlet_new(&x->x_obj, &s_signal);
	
	x->x_sr = 44100.0;
	x->x_n = 64.0;
	x->x_compReleaseLevel = 0.0;
	x->x_compAttackLevel = 0.0;
    x->x_compRMSArraySize = 64;
    x->x_compRMSScale = 65535.0f/(float)x->x_compRMSArraySize;
	x->x_compRMSArrayTotal = 0.0f;
	x->x_compRMSPosition = 0;

	// init tables
    for(i=0; i<(int)x->x_compRMSArraySize; i++)
    	x->x_compSquareArray[i] = 0.0;

    // range of rms table: from 0 to sqrt(65535) = 255.99/256.0
    // this is 17 bit precision at the top of the range
    // 8 bit precision at the bottom
    for(i=0; i<65536; i++)
	    x->x_compRoot[i] = pow((float)i, 0.5f) / 256.0f;

    for(i=0; i<8192; i++)
	    x->x_compGainTable[i] = 0.0;
          
	x->x_ratio = 39.0f/42.0f;
	x->x_gain = 1.0;
	x->x_RMSLevelDetect = 0;
	x->x_threshold = -40.0;
	x->x_attack = 0.02;
	x->x_compAttackInc = 1.0f/(x->x_attack * 0.001f * x->x_sr);
	x->x_compAttackMult = pow(10.0f, (1.0f/(x->x_attack * 0.001f * x->x_sr)));
	x->x_release = 0.5;
	x->x_compReleaseMult = pow(0.1f, (1.0f/(x->x_release*x->x_sr)));
	x->x_softKnee = 0.0;
	x->x_makeupGain = 1.0;
	x->x_compSoftHighThreshold = pow(10.0f, ((x->x_threshold + (x->x_softKnee * 6.0f)) * 0.05f));
	x->x_compSoftLowThreshold = pow(10.0f, (((20.0f * log10(x->x_compSoftHighThreshold)) - (x->x_softKnee * 12.0f)) * 0.05f));
    x->x_compDBSoftLowThreshold = 20.0f * log10(x->x_compSoftLowThreshold);

	// init compGainTable
	compand_tilde_compInitGainTable(x);
	
	return(void *)x;
};

static t_int *compand_tilde_perform(t_int *w)
{
	t_compand_tilde *x = (t_compand_tilde *)(w[1]);
	
	t_sample *in1 = (t_sample *)(w[2]);
	t_sample *out = (t_sample *)(w[3]);
	int n = (int)(w[4]);

    float level, gainmix, flevelindex;
    int32_t levelindex, sample;
	
	for(sample=0; sample<n; sample++)
	{
		x->x_compReleaseLevel *= x->x_compReleaseMult;
		x->x_compAttackLevel *= x->x_compAttackMult;
		if(x->x_RMSLevelDetect == 1)
			level = compand_tilde_compRMSDetect(x, *(in1+sample));
		else
			level = compand_tilde_compPeakDetect(*(in1+sample));
		// release algorithm before attack algorithm
				
		if(level > 1.0)
			level = 1.0;
		if(level < x->x_compReleaseLevel)
			x->x_compAttackLevel = level = x->x_compReleaseLevel;
		else if (level > x->x_compAttackLevel)
			level = x->x_compReleaseLevel = x->x_compAttackLevel;
		else if (level < 0.00001f)
			x->x_compReleaseLevel = x->x_compAttackLevel = 0.00001f;
		else
			x->x_compReleaseLevel = x->x_compAttackLevel = level;
		
		flevelindex = level * 8191.0f;

		levelindex = (int32_t)flevelindex;
		gainmix = flevelindex - levelindex;
		if(levelindex == 8191)
			x->x_gain = x->x_gain * 0.5f + 0.5f * x->x_compGainTable[levelindex];
		else
			x->x_gain = x->x_gain * 0.5f + 0.5f * (x->x_compGainTable[levelindex] 
				+ (x->x_compGainTable[levelindex+1] - x->x_compGainTable[levelindex]) * gainmix);
		*(out + sample) = *(in1 + sample) * x->x_gain * x->x_makeupGain;
	};

	return(w + 5);
};

static void compand_tilde_dsp(t_compand_tilde *x, t_signal **sp)
{
	dsp_add(
		compand_tilde_perform,
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

// unusual name for setup routine to allow for "+".
void setup_0x2bcompand_tilde(void)
{
	compand_tilde_class = 
	class_new(
		gensym("+compand~"),
		(t_newmethod)compand_tilde_new,
		0,
		sizeof(t_compand_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(compand_tilde_class, t_compand_tilde, x_f);

	class_addmethod(
		compand_tilde_class,
		(t_method)compand_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_ratio,
		gensym("ratio"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_threshold,
		gensym("threshold"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_attack,
		gensym("attack"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_release,
		gensym("release"),
		A_DEFFLOAT,
		0
	);
	
	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_softKnee,
		gensym("softKnee"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_makeupGain,
		gensym("makeupGain"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		compand_tilde_class, 
        (t_method)compand_tilde_rms_detect,
		gensym("rmsDetect"),
		A_DEFFLOAT,
		0
	);
}