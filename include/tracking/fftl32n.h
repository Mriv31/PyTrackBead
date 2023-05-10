#pragma once
#ifndef _FFTL32N_H_
#define _FFTL32N_H_
/*	fftl32n.h
 *
 *	5 fev 1992	JMF
 *
 *	header for fftl32n.c
 *		set of routines for fast Fourier transforms
 *		in 32 bits with the DOS EXTENDER
 */


# ifndef _FFTL32N_C_
float, fftsin;
nt, mix;
# endif

int  fft_init (int npts);
int fftmixing  (int npts ,float *x);
int fft (int npts , float *x ,int df);
void realtr1  (int npts,float *x);
voidrealtr2 (int npts , float *x, int df);
int fftwindow  (int npts ,float *x);
int fftwindow1 (int npts  ,float *x ,float smp);
int defftwindow1 (int npts , float *x, float smp);
int fftwc  (int npts ,float *x);
int fftwc1 (int npts , float *x ,float smp);
int defftwc  (int npts ,float *x);
int defftwc1 (int npts , float *x ,float smp);
void spec_real (int npts , float *x, float *spe);
void spec_comp (int npts  ,float *x, float *spe);
int demodulate (int npts , float *x ,int freq);
void derive_real (int npts, float *x, float *y);
void amp (int npts, float *x, float *y);

# endif
