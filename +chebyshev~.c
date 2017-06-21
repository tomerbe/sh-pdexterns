/*

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 12:28PM - removing unused variables and envelope follower stuff entirely.

*/

#include <math.h>
#include "m_pd.h"
#define ENVARRAYSIZE 16384

static t_class *chebyshev_tilde_class;

typedef struct _chebyshev_tilde
{
	t_object x_obj;
	float x_sr;
	float x_n;
	float x_dbOne;
	float x_dbTwo;
	float x_dbThree;
	float x_dbFour;
	float x_dbFive;
	float x_dbSeven;
	float x_gain;
    float x_f; 
} t_chebyshev_tilde;


static void chebyshev_tilde_p1(t_chebyshev_tilde *x, t_floatarg p1)
{	
	if(p1 > 0.0 || p1 < -60.0)
		error("value must be >= -60.0dB and <= 0.0dB.");
	else
		x->x_dbOne = pow(10.0f, p1 * 0.05);
}

static void chebyshev_tilde_p2(t_chebyshev_tilde *x, t_floatarg p2)
{	
	if(p2 > 00.0 || p2 < -60.0)
		error("value must be >= -60.0dB and <= 0.0dB.");
	else
		x->x_dbTwo = pow(10.0f, p2 * 0.05);
}

static void chebyshev_tilde_p3(t_chebyshev_tilde *x, t_floatarg p3)
{	
	if(p3 > 00.0 || p3 < -60.0)
		error("value must be >= -60.0dB and <= 0.0dB.");
	else
		x->x_dbThree = pow(10.0f, p3 * 0.05);
}

static void chebyshev_tilde_p4(t_chebyshev_tilde *x, t_floatarg p4)
{	
	if(p4 > 00.0 || p4 < -60.0)
		error("value must be >= -60.0dB and <= 0.0dB.");
	else
		x->x_dbFour = pow(10.0f, p4 * 0.05);
}

static void chebyshev_tilde_p5(t_chebyshev_tilde *x, t_floatarg p5)
{	
	if(p5 > 00.0 || p5 < -60.0)
		error("value must be >= -60.0dB and <= 0.0dB.");
	else
		x->x_dbFive = pow(10.0f, p5 * 0.05);
}

static void chebyshev_tilde_p7(t_chebyshev_tilde *x, t_floatarg p7)
{	
	if(p7 > 00.0 || p7 < -60.0)
		error("value must be >= -60.0dB and <= 0.0dB.");
	else
		x->x_dbSeven = pow(10.0f, p7 * 0.05);
}

static void *chebyshev_tilde_new()
{	
	t_chebyshev_tilde *x = (t_chebyshev_tilde *)pd_new(chebyshev_tilde_class);
	outlet_new(&x->x_obj, &s_signal);

	x->x_dbOne = 1.0;
	x->x_dbTwo = 0.0;
	x->x_dbThree = 0.0;
	x->x_dbFour = 0.0;
	x->x_dbFive = 0.0;
	x->x_dbSeven = 0.0;
	x->x_gain = 0.0;
	
	return(void *)x;
};

static t_int *chebyshev_tilde_perform(t_int *w)
{
	t_chebyshev_tilde *x = (t_chebyshev_tilde *)(w[1]);
	
	t_sample *in = (t_sample *)(w[2]);
	t_sample *out = (t_sample *)(w[3]);
	int n = (int)(w[4]);

	float gainOne, gainTwo, gainThree, gainFour, gainFive, gainSeven, gainSum;
	float poly1, poly2, poly3, poly4, poly5, poly7, square;
	float dbOne, dbTwo, dbThree, dbFour, dbFive, dbSeven;
	float gain;
	
	//local copies of dataspace variables:
	dbOne = x->x_dbOne;
	dbTwo = x->x_dbTwo;
	dbThree = x->x_dbThree;
	dbFour = x->x_dbFour;
	dbFive = x->x_dbFive;
	dbSeven = x->x_dbSeven;
	gain = x->x_gain;
	
	gainOne = pow(10.0f, ((dbOne * 3.f) - 3.f));
	gainTwo = pow(10.0f, ((dbTwo * 3.f) - 3.f));
	gainThree = pow(10.0f, ((dbThree * 3.f) - 3.f));
	gainFour = pow(10.0f, ((dbFour * 3.f) - 3.f));
	gainFive = pow(10.0f, ((dbFive * 3.f) - 3.f));
	gainSeven = pow(10.0f, ((dbSeven * 3.f) - 3.f));
	
	gainSum = gainOne + gainTwo + gainThree + gainFour + gainFive + gainSeven;

	if(gainSum == 0.0f)
		x->x_gain = gain = 0.0f;
	else
		x->x_gain = gain = 1.0f/gainSum;
	
	while(n--) 
	{
		float inputSample = *in++;
        float outputSample;

		// here's where you do your DSP work
		square = inputSample * inputSample;
		poly1 = inputSample;
		poly2 = (2.f * square) - 1.f;
		poly3 = ((4.f * square) - 3.f) * inputSample;
		poly4 = (((8.f * square) - 8.f) * square) + 1.f;
		poly5 = ((((16.f * square) - 20.f) * square) + 5.f) * inputSample;
//		poly6 = (((((32.f * square) - 48.f) * square) + 18.f) * square) - 1.f;
		poly7 = ((((((64.f * square) - 112.f) * square) + 56.f) * square) - 7.f) * inputSample;
		
		outputSample = (gainOne * poly1) + (gainTwo * poly2) + (gainThree * poly3) + (gainFour * poly4) 
			+ (gainFive * poly5) + (gainSeven * poly7);
		outputSample *= gain; 
						
		*out++ = outputSample;
	};
	
	return(w + 5);
};

static void chebyshev_tilde_dsp(t_chebyshev_tilde *x, t_signal **sp)
{
	dsp_add(
		chebyshev_tilde_perform,
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

void setup_0x2bchebyshev_tilde(void)
{
	chebyshev_tilde_class = 
	class_new(
		gensym("+chebyshev~"),
		(t_newmethod)chebyshev_tilde_new,
		0,
		sizeof(t_chebyshev_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(chebyshev_tilde_class, t_chebyshev_tilde, x_f);

	class_addmethod(
		chebyshev_tilde_class,
		(t_method)chebyshev_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		chebyshev_tilde_class, 
        (t_method)chebyshev_tilde_p1,
		gensym("p1"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		chebyshev_tilde_class, 
        (t_method)chebyshev_tilde_p2,
		gensym("p2"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		chebyshev_tilde_class, 
        (t_method)chebyshev_tilde_p3,
		gensym("p3"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		chebyshev_tilde_class, 
        (t_method)chebyshev_tilde_p4,
		gensym("p4"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		chebyshev_tilde_class, 
        (t_method)chebyshev_tilde_p5,
		gensym("p5"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		chebyshev_tilde_class, 
        (t_method)chebyshev_tilde_p7,
		gensym("p7"),
		A_DEFFLOAT,
		0
	);
}