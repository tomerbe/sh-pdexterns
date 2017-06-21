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

/* ------------------------ morphfilter~ ----------------------------- */


static t_class *morphfilter_tilde_class;

// some convenient constants
enum
{
	kSizeFFT = 1024,
	kHalfSizeFFT = 512
};


struct mfilter {

	t_word *x_vec;
};

typedef struct _morphfilter_tilde
	{
		t_object x_obj;         /* obligatory header */
		t_float gain;            /* might be neccesary to do some gain  compensation */
		
		t_float sampleRate;
		t_int   bufferPosition;
		t_int   impulseSize;
		double denormalValue;
		long blockSize;
		// FFT stuff
		float *inBuffer, *outBuffer, *inWindowed;
		float *analysisWindow, *synthesisWindow;
		float *inShift, *outShift;
		long halfSizeFFT, sizeFFT;
		long inputTime, outputTime;
		float pi, twoPi;
		float *filterSpectra;
		
		float *dbTable1;
		float *dbTable2;
		float *filterLearn;
		
		float *gainTable;
		float *tiltTable;
		float *squareLevel;
		
		// These are the user defined tables 
		int x_array_points1;
		int x_array_points2;
		t_symbol *x_arrayname1;
		t_symbol *x_arrayname2;
					
		
		t_float depthFilter;

		t_float tilt;
		
		t_float learning;
		t_float learnValue;
		t_float learnFrames;
		
		int indexFilter;
		t_float numberFilter;
		
		struct mfilter filter[2];
		
	} t_morphfilter_tilde;


/*these functions clip incoming parameter values*/
void morphfilter_tilde_setDepth(t_morphfilter_tilde *x, t_float f)
{
	if(f>2.0 || f<-2.0)
		error("depth value must be >= -2 and <= 2.");
	else
		x->depthFilter = f;		
}


void morphfilter_tilde_morph(t_morphfilter_tilde *x, t_float value) {

	 if(value <= 1 && value >= 0) 
		{ x->numberFilter = value; }
	else
		error("morph value must be >= 0 and <= 1.");
	
	if(x->numberFilter >= 0.5f)
			x->indexFilter = 1;
	else
		x->indexFilter = 0;
	
//	filterTableL = filter[indexFilter].linTableL;
//	filterTableR = filter[indexFilter].linTableR;
//	mixFilter();

}


void morphfilter_tilde_tilt(t_morphfilter_tilde *x , t_float value)
{
    int i;
    float tiltBasis;
    
    if ((value <= 3) && (value >= -3))  
	{ 
		x->tilt = value;	
		tiltBasis = log10f((float)x->halfSizeFFT) * 20.0f * 0.5f;

		for(i=0; i <= x->halfSizeFFT; i++)
			x->tiltTable[i] = powf(10.0f, ((log10f((float)i+1) * 20.0f) -  tiltBasis) * (x->tilt/tiltBasis));
	
	} 
	else
		error("tilt value must be >= -3 and <= 3.");
	
}

	
void morphfilter_tilde_set(t_morphfilter_tilde *x, t_symbol *s1, t_symbol *s2) {
	t_garray *a1;
	t_garray *a2;
	int old_array_points1;
	int old_array_points2;
	
	x->x_arrayname1 = 0;
	x->x_arrayname2 = 0;
	
	old_array_points1 = x->x_array_points1;
	old_array_points2 = x->x_array_points2;	
	
	if(!(a1 = (t_garray *)pd_findbyclass(s1, garray_class)))
	{
		x->filter[0].x_vec = 0;
		pd_error(x, "%s: no such array", s1->s_name);
	}
	else if(!garray_getfloatwords(a1, &x->x_array_points1, &x->filter[0].x_vec))
	{
		x->filter[0].x_vec = 0;
		pd_error(x, "%s: bad template for morphfilter", s1->s_name);
	}
	else {
		garray_usedindsp(a1);
		x->x_arrayname1 = s1;
		garray_resize( (t_garray *)pd_findbyclass(x->x_arrayname1,garray_class), 513.0);
		garray_getfloatwords(a1, &x->x_array_points1, &x->filter[0].x_vec);
	}
		
	if(!(a2 = (t_garray *)pd_findbyclass(s2, garray_class)))
	{
		x->filter[1].x_vec = 0;
		pd_error(x, "%s: no such array", s2->s_name);
	}
	else if(!garray_getfloatwords(a2, &x->x_array_points2, &x->filter[1].x_vec))
	{
		x->filter[1].x_vec = 0;
		pd_error(x, "%s: bad template for morphfilter", s2->s_name);
	}
	else {
		x->x_arrayname2 = s2;
		garray_usedindsp(a2);
		garray_resize( (t_garray *)pd_findbyclass(x->x_arrayname1,garray_class), 513.0);
		garray_getfloatwords(a2, &x->x_array_points2, &x->filter[1].x_vec);
	}
		
}


void morphfilter_tilde_mixFilter(t_morphfilter_tilde *x) {
	// depth needs to be applied to the filter to make a new filtercurve
	// one could mix between two filter powers (f * f) and (f * f * f)
	// for the range between 2 and 3, or just use pow(f, 2.4). the previous
	// should be faster
	// negative depths should invert the filter shape first f = (1.0f - f)
	// zero depth should result in a straight line of 1.0f
	
	int i;
	float dbValue;
	if(x->filter[0].x_vec == 0 || x->filter[1].x_vec == 0)
		return;
	//	if(modified)    //find a way to check if table has been modified - modified is checked in mixFilter, probably doesn't need to be re-checked
	//	{
	/*
	 for(i = 0; i <= x->halfSizeFFT; i++)   //convert to dbValues
	 {
	 //this may cause problems??

	 x->dbTable1[i] = 20.0f *log10f(  x->x_vec1[i].w_float);
	 x->dbTable2[i] = 20.0f *log10f(  x->x_vec2[i].w_float);
	 } */
	//		modified = false;
	//	}
	
	
	for(i = 0; i <= x->halfSizeFFT; i++)	   // mix the two values
	{
		dbValue = x->depthFilter * (x->filter[0].x_vec[i].w_float * (1.0f - x->numberFilter) + x->filter[1].x_vec[i].w_float * x->numberFilter);		
		x->filterSpectra[i] = powf(10.0f, dbValue * 0.05f);
	}
}


void morphfilter_tilde_initHammingWindows ( t_morphfilter_tilde *x ) {
	long N, k;
	N = x->sizeFFT;
	for ( k = 0; k < N; k++ )
		x->analysisWindow[k] = x->synthesisWindow[k] = (float) (.54f - (.46f * cosf(x->twoPi * k / (N - 1)) ) );

}


void morphfilter_tilde_scaleWindows (t_morphfilter_tilde *x) {
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


void morphfilter_tilde_setLearn (t_morphfilter_tilde *x, t_float valIn) 
{
	int i;

	x->learnValue = valIn;
	if(x->filter[0].x_vec == 0 || x->filter[1].x_vec == 0)
		return;
 // start filter learning

	if(x->learnValue > 0.5f && x->learning == FALSE) {
		for(i=0; i <= x->halfSizeFFT; i++) 
		{
			x->filter[ x->indexFilter].x_vec[i].w_float = 0.0f;
			x->filterLearn[i] = 0.0f;
		}

		
		x->learning = TRUE;
	}

	else if(x->learnValue < 0.5f && x->learning == TRUE) {
		x->learning = 0;
		if(x->learnFrames > 0) {

		// we are gonna normalize this shape between 24 and -24
		// lets look in the range 0 to -96
		// without doing the divide by 20, that's from -4.8 to 0
		// to -1.2 to 1.2

			float filterHigh = -4.8f;        // lowest possible
			float filterLow = 0.0f;        // highest possible

			for(i=0; i <= x->halfSizeFFT; i++) {
				x->filterLearn[i] = x->filterLearn[i]/ x->learnFrames;
				x->filterLearn[i] = log10f(x->filterLearn[i] + 0.000001f);
				if(x->filterLearn[i] > filterHigh)
					filterHigh = x->filterLearn[i];
				if(x->filterLearn[i] < filterLow)
					filterLow = x->filterLearn[i];
			}
			
			float filterRangeAdjust = 2.4f/(filterHigh-filterLow);

			for(i=0; i <= x->halfSizeFFT; i++) 
			{
				float temp = (x->filterLearn[i] - filterHigh) * filterRangeAdjust;
				x->filter[x->indexFilter].x_vec[i].w_float = (temp * 20.f) + 24.0f;

			}
		}
		x->learnFrames = 0;
		morphfilter_tilde_mixFilter(x);
		garray_redraw((t_garray *)pd_findbyclass(x->x_arrayname1,garray_class));
		garray_redraw((t_garray *)pd_findbyclass(x->x_arrayname2,garray_class));

	}
}


void *morphfilter_tilde_new(t_symbol *table1, t_symbol *table2) {
	long i;
	t_morphfilter_tilde *x = (t_morphfilter_tilde *)pd_new(morphfilter_tilde_class);
	// an inlet to bring in filter coefficients

	//	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("signal"),0); // make this float inlet?

	inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"),gensym("interptVal")); // make this float inlet?

	outlet_new(&x->x_obj, gensym("signal"));
		
	// add method to set table names
		
	//	expects 2 table names 

	x->x_arrayname1 = 0;
	x->x_arrayname2 = 0;
	x->filter[0].x_vec = x->filter[1].x_vec = 0;
	
	if(table1 && table2) {  
		x->x_array_points1 = 0;	
		x->x_array_points2 = 0;	
		morphfilter_tilde_set(x,table1, table2); 
	}
 
	x->gain = 1.0f;							/* might be neccesary to do some gain compensation */
	x->bufferPosition = 0;
	x->inputTime = x->outputTime = 0;
		
	x->sizeFFT = kSizeFFT;
	x->blockSize = x->sizeFFT >> 2;
	x->halfSizeFFT = x->sizeFFT >> 1;
	x->pi = 4.0f * atanf(1.0f);
	x->twoPi = 8.0f * atanf(1.0f);

	x->learning = FALSE;
	x->learnValue = 0;
	x->learnFrames = 0;
	
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
	x->filterSpectra = (float *) malloc(sizeof(float) * x->sizeFFT);
	
	
	x->gainTable = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->tiltTable = (float *) malloc(sizeof(float) * x->sizeFFT);
	
	x->filterLearn = (float *) malloc(sizeof(float) * x->sizeFFT);
	x->squareLevel = (float *) malloc(sizeof(float) * x->sizeFFT);		

	// third - zero out all the memory
	for(i = 0; i<x->sizeFFT; i++) {
		x->inBuffer[i] = x->outBuffer[i] = x->inShift[i] = x->outShift[i] = 0.0;
		x->filterSpectra[i] = 0.0;
	}
	//        fourth - set up the FFT and windows
	morphfilter_tilde_initHammingWindows( x );		// 
	morphfilter_tilde_scaleWindows(x);				// these we keep
	
	x->depthFilter = 0.0f;
	x->indexFilter = 0;
	x->numberFilter = 0.0f;
	
	x->tilt = 0.0f;
	float tiltBasis = log10f((float) x->halfSizeFFT) * 20.0f * 0.5f;
	for(i = 0; i <= x->halfSizeFFT; i++)
		x->tiltTable[i] = powf(10.0f, ((log10f((float)i+1) * 20.0f) -  tiltBasis) * (x->tilt/tiltBasis));
		
		
	// this is how you access values in table: x->x_vec[i].w_float
		
	return (x);
}


void morphfilter_tilde_processSpect(t_morphfilter_tilde *x) {
	long i;
	
	
	// if peak detection is on, the filter is based on the peak value  plus the filter value
	// the range of filter values changes in the GUI from 0 to -96 (no peak) to +/- 48 (peak),
	// so we can just combine the filter with the peak
	
	// first - get the levels for all bands and identify the peaks
	// left channel

		 
	if(x->learning == TRUE && x->learnValue != 0.0f) {

		
		x->squareLevel[0] = x->inWindowed[0] * x->inWindowed[0];
		x->squareLevel[x->halfSizeFFT] = x->inWindowed[x->halfSizeFFT] *  x->inWindowed[x->halfSizeFFT];

		for(i=1; i < x->halfSizeFFT; i++)
			x->squareLevel[i] = x->inWindowed[i] * x->inWindowed[i] + x->inWindowed[x->sizeFFT - i] * x->inWindowed[x->sizeFFT - i];

		for(i = 0; i< x->halfSizeFFT; i++)
			x->filterLearn[i] += x->squareLevel[i];

		x->learnFrames++;

		if(x->learnFrames > 3000)
			x->learnValue = 0.0f;
		
	}
	
	
	//				x->filter[indexFilter].x_vec[i].w_float =  x->filter[indexFilter].dbTable[i];
	//	if(modified) 
	//	{
	
	
	morphfilter_tilde_mixFilter(x);
	
	
	//	}
	
	
	// now we reset the gain table, multiplying by the attackFactor or  releaseFactor depending on whether
	// above or below the filter
	
	
	for(i = 0; i <= x->halfSizeFFT; i++)
		x->gainTable[i] = (x->gainTable[i] * 0.6f) + (x->filterSpectra[i]  *  x->tiltTable[i] * 0.4f);
	
	
	// c - gain the spectra
	x->outShift[0] = x->inWindowed[0] * x->gainTable[0];        // DC Component
	x->outShift[x->halfSizeFFT] = x->inWindowed[x->halfSizeFFT] *  x->gainTable[x->halfSizeFFT];        // Nyquist Frequency
	
	for(i = 1; i < x->halfSizeFFT; ++i)
	{
		x->outShift[i] = x->inWindowed[i] * x->gainTable[i];
		x->outShift[x->sizeFFT - i] = x->inWindowed[x->sizeFFT - i] * x->gainTable[i];
	}
}


void morphfilter_tilde_block(t_morphfilter_tilde *x) {
	long i;
	//	long j;
	long maskFFT = x->sizeFFT - 1;
		
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
	morphfilter_tilde_processSpect(x);
	
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
static t_int *morphfilter_tilde_perform(t_int *w) {
		
	t_morphfilter_tilde *x = (t_morphfilter_tilde *)(w[1]);
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
		
	//	if(bandsPerIndex <= 1)
	//		memcpy(x->filterSpectra, freqResponse, n);
	//	else
	//        // if there are more bands in the fft than the number of numbers in the freqResponse
	//        // we will spread out the freqResponse to fit the fft
	//	{
	//		for(i=0; i<n; i++)
	//			for(j=0; j<bandsPerIndex; j++)
	//				*(x->filterSpectra+(i*bandsPerIndex)+j) = *(freqResponse+i);
	//	}
		
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
			out[i] = x->outBuffer[i+x->bufferPosition] * x->gain;
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
			morphfilter_tilde_block(x);
		}
		// decrement framesLeft by the number of frames (samples) processed
		framesLeft -= processframes;
	}
	return (w + 5);
}


void morphfilter_tilde_dsp(t_morphfilter_tilde *x, t_signal **sp) {
	x->sampleRate = sp[0]->s_sr;
		
	dsp_add(morphfilter_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}
   	
	
// since we allocated some memory, we need a delete function
static void morphfilter_tilde_free(t_morphfilter_tilde *x) {
	if(x->inBuffer != 0) free(x->inBuffer);
	if(x->inWindowed != 0) free(x->inWindowed);
	if(x->outBuffer != 0) free(x->outBuffer);
	if(x->inShift != 0) free(x->inShift);
	if(x->outShift != 0) free(x->outShift);
	if(x->filterSpectra != 0) free(x->filterSpectra);
	
	if(x->analysisWindow != 0) free(x->analysisWindow);
	if(x->synthesisWindow != 0) free(x->synthesisWindow);
	
	if(x->dbTable1 != 0) free(x->dbTable1);
	if(x->dbTable2 != 0) free(x->dbTable2);
	if(x->gainTable != 0) free(x->gainTable);
	if(x->tiltTable != 0) free(x->tiltTable);
	
	if(x->filterLearn != 0) free(x->filterLearn);
	if(x->squareLevel != 0) free(x->squareLevel);		
/*
	if(x->filter[0].dbTable != 0) free(x->filter[0].dbTable);
	if(x->filter[1].dbTable != 0) free(x->filter[1].dbTable);

	if(x->filter[0].linTable != 0) free(x->filter[0].linTable);
	if(x->filter[1].linTable != 0) free(x->filter[1].linTable);

 */	
	/*
	if(x->filter[0].x_vec != 0) free(x->filter[0].x_vec);
	if(x->filter[1].x_vec != 0) free(x->filter[1].x_vec); */

}

/* this routine, which must have exactly this name (with the "~"  replaced 
 * by "_tilde) is called when the code is first loaded, and tells Pd  how to build the "class". */
      

void setup_0x2bmorphfilter_tilde(void) {
	morphfilter_tilde_class = class_new(gensym("+morphfilter~"),  
				                 (t_newmethod)morphfilter_tilde_new, 
								 (t_method)morphfilter_tilde_free,
								 sizeof(t_morphfilter_tilde), 
								 CLASS_DEFAULT, 
								 A_DEFSYMBOL,
								 A_DEFSYMBOL,
								 0);
	/* this is magic to declare that the leftmost, "main" inlet
	 * takes signals; other signal inlets are done differently... */

	/* also installs delay_time as the leftmost inlet float */	

	CLASS_MAINSIGNALIN(morphfilter_tilde_class, t_morphfilter_tilde, gain);
	/* here we tell Pd about the "dsp" method, which is called back when DSP is turned on. */

	class_addmethod(morphfilter_tilde_class, (t_method) morphfilter_tilde_dsp,  
						gensym("dsp"), (t_atomtype)0);
		
	class_addmethod(morphfilter_tilde_class, (t_method) morphfilter_tilde_set, gensym("set"), A_SYMBOL, A_SYMBOL, 0);

	class_addmethod(morphfilter_tilde_class, (t_method) morphfilter_tilde_setLearn, gensym("learn"), A_DEFFLOAT,0);
	//		class_addmethod(morphfilter_tilde_class, (t_method) morphfilter_tilde_interptVal, gensym("interptVal"), A_DEFFLOAT, 0);
		
	class_addmethod(morphfilter_tilde_class, (t_method)morphfilter_tilde_setDepth, gensym("depth"), A_DEFFLOAT,0);
	class_addmethod(morphfilter_tilde_class, (t_method)morphfilter_tilde_morph, gensym("morph"), A_DEFFLOAT, 0);
	class_addmethod(morphfilter_tilde_class, (t_method)morphfilter_tilde_tilt, gensym("tilt"), A_DEFFLOAT, 0); 		
}