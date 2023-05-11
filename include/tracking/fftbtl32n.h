#pragma once
#ifndef _FFTBTL32N_H_
#define _FFTBTL32N_H_
/*	fftbtl32n.h
 *
 *	5 fev 1992	JMF
 *
 *	header for fftbtl32n.c
 *		set of routines for fast Fourier transforms
 *		in 32 bits with the DOS EXTENDER
 */

#ifdef _WIN32
#define M_PI 3.141592653589793238462643383279502884197
#endif

#include <math.h>
#include <stdlib.h>


typedef struct _fft_plan
{
  int nx;                  // the total number of frames requires to perform calibration
  int n_2, n_4;            // nx/2 and nx/4
  float    *fftbtsin;      // the sinus array
  int 	   *mixbt;         // the mixing array
  int pow_2[31];
} fft_plan;


#ifdef __cplusplus
extern "C" {
#endif

fft_plan* fftbt_init(fft_plan *initial, int npts);


#ifdef __cplusplus
}
#endif


int free_fft_plan_bt(fft_plan *fp);
int fftbtmixing(fft_plan *fp, float *x);
int fftbt(fft_plan *fp, float *x, int df);
void realtr1bt(fft_plan *fp, float *x);
void realtr2bt(fft_plan *fp, float *x, int df);
int fftbtwindow(fft_plan *fp, float *x);
int defftbtwindow(fft_plan *fp, float *x);
int fftbtwindow1(fft_plan *fp, float *x, float smp);
int defftbtwindow1(fft_plan *fp, float *x, float smp);
int fftbtwc(fft_plan *fp, float *x);
int fftbtwc1(fft_plan *fp, float *x, float smp);
int defftbtwc(fft_plan *fp, float *x);
int defftbtwc1(fft_plan *fp, float *x, float smp);
void spec_realbt(fft_plan *fp, float *x, float *spe);
void spec_compbt(fft_plan *fp, float *x, float *spe);
int demodulatebt(fft_plan *fp, float *x, int freq);
void derive_realbt(fft_plan *fp, float *x, float *y);
void ampbt(fft_plan *fp, float *x, float *y);

# endif
