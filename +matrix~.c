/*

May 20, 2010: using a single makefile for all sources that should work on both mac and windows. removing the sethelpsymbol() and going with the standard whatever~-help.pd convention. also adding static to all function names (except setup).

Nov 27, 2009: 12:35PM - removed unused variables.

*/

#include <math.h>
#include "m_pd.h"
#define ENVARRAYSIZE 16384

static t_class *matrix_tilde_class;

typedef struct _matrix_tilde
{
	t_object x_obj;
	float x_sr;
	float x_n;
	float x_midGain;
	float x_sideGain;
	int x_direction;
    float x_f; 
} t_matrix_tilde;

static void matrix_tilde_midGain(t_matrix_tilde *x, t_floatarg mg)
{	
	if(mg > 12.0 || mg < -12.0)
		error("value must be >= -12.0dB and <= 12.0dB.");
	else
		x->x_midGain = pow(10.0, mg*0.05);
}

static void matrix_tilde_sideGain(t_matrix_tilde *x, t_floatarg sg)
{	
	if(sg > 12.0 || sg < -12.0)
		error("value must be >= -12.0dB and <= 12.0dB.");
	else
		x->x_sideGain = pow(10.0, sg*0.05);
}

static void matrix_tilde_direction(t_matrix_tilde *x, t_floatarg d)
{	
	if(d > 1.0)
		x->x_direction = 1;
	else if(d < 0.0)
		x->x_direction = 0;
	else
		x->x_direction = (int)d;
}

static void *matrix_tilde_new()
{
	t_matrix_tilde *x = (t_matrix_tilde *)pd_new(matrix_tilde_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("midGain"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("sideGain"));
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("direction"));
	outlet_new(&x->x_obj, &s_signal);
	outlet_new(&x->x_obj, &s_signal);

	x->x_sr = 44100.0;
	x->x_n = 64.0;
	x->x_midGain = x->x_sideGain = 1.0;
	x->x_direction = 0;
	
	return(void *)x;
};

static t_int *matrix_tilde_perform(t_int *w)
{
	t_matrix_tilde *x = (t_matrix_tilde *)(w[1]);
	
	t_sample *in1 = (t_sample *)(w[2]);
	t_sample *in2 = (t_sample *)(w[3]);
	t_sample *out1 = (t_sample *)(w[4]);
	t_sample *out2 = (t_sample *)(w[5]);
	int n = (int)(w[6]);

	long framesLeft;
	long frame;
	float mid, side, left, right;

	framesLeft = n;
	
	// LR to MS - false - 0
	// MS to LR - true - 1
	if(x->x_direction == 0)
	{
		// direction 0 is MS encode lr to ms
		for(frame = 0; frame < framesLeft; frame++)
		{
			left = in1[frame];
			right = in2[frame];
			out1[frame] = (left + right) * x->x_midGain * 0.5f; // m = l + r
			out2[frame] = (left - right) * x->x_sideGain * 0.5f; // s = l - r
		}
	}
	else // initial direction - ms to lr
	{
		for(frame = 0; frame < framesLeft; frame++)
		{
			mid = in1[frame];
			side = in2[frame];
			out1[frame] = (mid * x->x_midGain) + (side * x->x_sideGain); // l = m + s
			out2[frame] = (mid * x->x_midGain) - (side * x->x_sideGain); // r = m - s
		}
	}
	
	return(w + 7);
};

static void matrix_tilde_dsp(t_matrix_tilde *x, t_signal **sp)
{
	dsp_add(
		matrix_tilde_perform,
		6,
		x,
		sp[0]->s_vec,
		sp[1]->s_vec,
		sp[2]->s_vec,
		sp[3]->s_vec,
		sp[0]->s_n
	);

	if(sp[0]->s_n != x->x_n || sp[0]->s_sr != x->x_sr)
	{
		x->x_n = sp[0]->s_n;
		x->x_sr = sp[0]->s_sr;
	};
};

void setup_0x2bmatrix_tilde(void)
{
	matrix_tilde_class = 
	class_new(
		gensym("+matrix~"),
		(t_newmethod)matrix_tilde_new,
		0,
		sizeof(t_matrix_tilde),
		CLASS_DEFAULT,
		0
	);

    CLASS_MAINSIGNALIN(matrix_tilde_class, t_matrix_tilde, x_f);

	class_addmethod(
		matrix_tilde_class,
		(t_method)matrix_tilde_dsp,
		gensym("dsp"),
		0
	);

	class_addmethod(
		matrix_tilde_class, 
        (t_method)matrix_tilde_midGain,
		gensym("midGain"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		matrix_tilde_class, 
        (t_method)matrix_tilde_sideGain,
		gensym("sideGain"),
		A_DEFFLOAT,
		0
	);

	class_addmethod(
		matrix_tilde_class, 
        (t_method)matrix_tilde_direction,
		gensym("direction"),
		A_DEFFLOAT,
		0
	);
}