/*

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 12:17PM - changed to int32_t instead of long for longTmp lMask and lMask2 in the perform routine. Also included stdint.h to enable this.

*/

#include <stdint.h>
#include <math.h>
#include "m_pd.h"
#define STDBLOCK 256

static t_class *decimate_tilde_class;

typedef struct _decimate_tilde
{
	t_object x_obj;
	float x_sr;
	float x_n;
	int x_block;
	t_sample *x_inBuffer;
	t_sample *x_outBuffer;
	long x_lMask;
	long x_lMask2;
	float x_maskMix;
	float x_averageScale;
	float x_avg;
    float x_f;
} t_decimate_tilde;


static void decimate_tilde_depth(t_decimate_tilde *x, t_floatarg d)
{	
	unsigned long mask, bit;
	long i, bitDepth;
	
	if(d>32.0 || d<1.0)
		error("depth value must be >= 1.0 and <= 24.0.");
	else
	{
		bitDepth = (long)d;
		x->x_maskMix = d - bitDepth;		
		mask = 0;
		for(i = 0, bit = 0x80000000UL; i < bitDepth; i++)
		{
			mask = mask | bit;
			bit = bit>>1;
		}
		x->x_lMask = mask;
		mask = 0;
		for(i = 0, bit = 0x80000000UL; i < bitDepth+1; i++)
		{
			mask = mask | bit;
			bit = bit>>1;
		}
		x->x_lMask2 = mask;
	};	
}

static void decimate_tilde_folds(t_decimate_tilde *x, t_floatarg f)
{	
	if(f>8.0 || f<0.0)
		error("folds value must be between 0 and 8.");
	else
	{
		x->x_avg = powf(2.0, f);
		x->x_averageScale = 1.0/x->x_avg;
	}
}

static void *decimate_tilde_new()
{	
	t_decimate_tilde *x = (t_decimate_tilde *)pd_new(decimate_tilde_class);
	int i;

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("depth"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("folds"));
	outlet_new(&x->x_obj, &s_signal);
	x->x_inBuffer = (t_sample *)getbytes(0);

	x->x_inBuffer = (t_sample *)t_resizebytes(x->x_inBuffer, 0, STDBLOCK * sizeof(t_sample));

	// init signal buffer
	for(i=0; i<STDBLOCK; i++)
		x->x_inBuffer[i] = 0.0;

	x->x_outBuffer = (t_sample *)t_resizebytes(x->x_outBuffer, 0, STDBLOCK * sizeof(t_sample));

	// init out buffer
	for(i=0; i<STDBLOCK; i++)
		x->x_outBuffer[i] = 0.0;
		
	x->x_n = 64;
	x->x_block = STDBLOCK;
	x->x_lMask = 0xffffffff;
	x->x_lMask2 = 0xffffffff;
	x->x_maskMix = 0.0;
	x->x_avg = 1.0;
	x->x_averageScale = 1.0/x->x_avg;

	return(void *)x;
};

static t_int *decimate_tilde_perform(t_int *w)
{
	int i, j, block, avg;
	float tmp, maskMix, avgScale;
	int32_t longTmp, lMask, lMask2;
	t_decimate_tilde *x = (t_decimate_tilde *)(w[1]);
	
	t_sample *in1 = (t_sample *)(w[2]);
	t_sample *out = (t_sample *)(w[3]);
	int n = (int)(w[4]);
	t_sample *pOutBuffer;

	// grab a pointer to outBuffer for the final while(n--)
	pOutBuffer = x->x_outBuffer;
	
	// rename dataspace variables locally
	block = x->x_block;
	avg = (int)x->x_avg;
	avgScale = x->x_averageScale;
	lMask = x->x_lMask;
	lMask2 = x->x_lMask2;
	maskMix = x->x_maskMix;

	// init out buffer
	for(i=0; i<n; i++)
		x->x_outBuffer[i] = 0.0;

	// shift signal buffer contents back.
	for(i=0; i<(block-n); i++)
		x->x_inBuffer[i] = x->x_inBuffer[i+n];
	
	// write new block to end of signal buffer.
	for(i=0; i<n; i++)
		x->x_inBuffer[block-n+i] = in1[i];
	
	for(i=0; i<n; i+=avg)
	{
		for(j=0, tmp=0; j<avg; j++)
			tmp += x->x_inBuffer[i+j];

		tmp *= avgScale; // avoid divides using precomputed avgScale
		
		for(j=0, longTmp=0; j<avg; j++)
		{	
			longTmp = (int32_t)(tmp * 2147483647.0f);
			longTmp = (int32_t)((longTmp & lMask) + ((longTmp & lMask2) - (longTmp & lMask)) * maskMix);
			x->x_outBuffer[i+j] = (float)longTmp * 0.0000000004656612875245797f;
		};
	};	

	
	while(n--)
		*out++ = *pOutBuffer++;
	
	return(w+5);
};

static void decimate_tilde_dsp(t_decimate_tilde *x, t_signal **sp)
{
	int i;
	
	dsp_add(
		decimate_tilde_perform,
		4,
		x,
		sp[0]->s_vec,
		sp[1]->s_vec,
		sp[0]->s_n
	);
	
	// update blocksize if it has changed
	if(x->x_n != sp[0]->s_n)
	{
		if(sp[0]->s_n > STDBLOCK)
		{
			// resize signal buffer
			x->x_inBuffer = (t_sample *)t_resizebytes(x->x_inBuffer, x->x_block * sizeof(t_sample), sp[0]->s_n * sizeof(t_sample));
				
			// init signal buffer
			for(i=0; i<sp[0]->s_n; i++)
				x->x_inBuffer[i] = 0.0;

			// resize out buffer
			x->x_outBuffer = (t_sample *)t_resizebytes(x->x_outBuffer, x->x_block * sizeof(t_sample), sp[0]->s_n * sizeof(t_sample));
				
			// init out buffer
			for(i=0; i<sp[0]->s_n; i++)
				x->x_outBuffer[i] = 0.0;
				
			x->x_block = sp[0]->s_n;
		}
		else if( (sp[0]->s_n <= STDBLOCK) && (x->x_n > STDBLOCK) )
		{
			// resize signal buffer
			x->x_inBuffer = (t_sample *)t_resizebytes(x->x_inBuffer, x->x_block * sizeof(t_sample), STDBLOCK * sizeof(t_sample));
				
			// init signal buffer
			for(i=0; i<STDBLOCK; i++)
				x->x_inBuffer[i] = 0.0;

			// resize out buffer
			x->x_outBuffer = (t_sample *)t_resizebytes(x->x_outBuffer, x->x_block * sizeof(t_sample), STDBLOCK * sizeof(t_sample));
				
			// init out buffer
			for(i=0; i<STDBLOCK; i++)
				x->x_outBuffer[i] = 0.0;
				
			x->x_block = STDBLOCK;
		};

		x->x_n = sp[0]->s_n;
	}

	if(sp[0]->s_sr != x->x_sr)
		x->x_sr = sp[0]->s_sr;
};

static void decimate_tilde_free(t_decimate_tilde *x)
{
	// free the input buffer memory
    t_freebytes(x->x_inBuffer, x->x_block*sizeof(t_sample));

	// free the input buffer memory
    t_freebytes(x->x_outBuffer, x->x_block*sizeof(t_sample));
};

// unusual name for setup routine to allow for "+".
void setup_0x2bdecimate_tilde(void)
{
	decimate_tilde_class = 
	class_new(
		gensym("+decimate~"),
		(t_newmethod)decimate_tilde_new,
		(t_method)decimate_tilde_free,
		sizeof(t_decimate_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(decimate_tilde_class, t_decimate_tilde, x_f);

	class_addmethod(
		decimate_tilde_class,
		(t_method)decimate_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		decimate_tilde_class, 
        (t_method)decimate_tilde_depth,
		gensym("depth"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		decimate_tilde_class, 
        (t_method)decimate_tilde_folds,
		gensym("folds"),
		A_DEFFLOAT,
		0
	);
}