#ifndef _FFTBTL32N_C_
#define _FFTBTL32N_C_

/*   fftbtl32n.c
 *
 *	5 fev. 92	VC
 *	5 fev. 92	revision JMF
 *	20 mars 92	rev JMF function defftbtwindow1
 *
 *	A collection of routines for performing Fast Fourier Transforms on
 *	real and complex data. In all routines, npts is the number of
 *	elements in the data set, so for complex input npts is twice the
 *	number of complex data pairs. x is an array of floats which is the
 *	time series on input and is replaced by the fftbt output. The length
 *	of the array is npts. For complex input data, the real and imaginary
 *	parts are alternated in the array, i.e.
 *	x[0] = re{z[0]}, x[1] = im{z[0]}, x[2] = re{z[1]}, x[3] = im{z[1]} ...
 *	The same is true on output (excepted for z[1] in the case of real
 *	data).
 */
# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include "fftbtl32n.h"

/* # include "fftbtl32n.i"*/





# ifdef FORTIFY
# include "fortify.h"
# endif

//static int oldnpts = 0, nptsmax = 0;

/*	fftbt_init()
 *	DESCRIPTION	Calculates the array of sines that will be needed
 *			for the fftbt. This routine is called by fftbt(). The
 *			user may call this routine explicitly before calling
 *			fftbt() for the sake of using the array fftbtsin[] for
 *			windowing prior to calculating the fftbt.
 *
 *	RETURNS		0 on success, 1 on npts which is not a power of two,
 *			2 on array memory could not be allocated.
 */
fft_plan *fftbt_init(fft_plan *initial, int npts)
{
  register int i, j, k;
  int n=1, nbits;
  fft_plan *fp = NULL;

  if (initial != NULL && initial->nx == npts)  return initial;

  for (i = 0; i < 30; i++)
    {/* Check npts is a power of 2 */
      if (npts == n) break;
      n *= 2;
    }
  if (i == 30)    return NULL;

  if (initial == NULL)
    {
      fp = (fft_plan*)calloc(1,sizeof(fft_plan));
      if (fp == NULL)  return NULL;
      fp->mixbt = (int*)calloc(npts,sizeof(int));
      fp->fftbtsin = (float*)calloc(npts,sizeof(float));
    }
  else
    {
      initial->mixbt = (int*)realloc(initial->mixbt, npts*sizeof(int));
      initial->fftbtsin = (float*)realloc(initial->fftbtsin, npts*sizeof(float));
      fp = initial;
    }
  if (fp->mixbt == NULL || fp->fftbtsin == NULL)  return NULL;


  fp->nx = npts;
  fp->n_2 = npts/2;
  fp->n_4 = npts/4;

  for (i = 0 ; i <= fp->n_4 ; i++)
      fp->fftbtsin[i] =  sin ((M_PI * i) / fp->n_2);
  fp->pow_2[0] = fp->n_4;
  i=0;
  do
    {
      i++;
      fp->pow_2[i] = fp->pow_2[i - 1] >> 1;
    } while (fp->pow_2[i] != 1 );
  nbits = i;

  for (j=0; j < fp->n_2; j++)
    {
      k=0;
      for (i=0; i <= nbits; i++)
	{
	  if ((j & fp->pow_2[i] ) > 0)	      k += fp->pow_2[nbits-i];
	}
      fp->mixbt[j] = ( k > j)  ? 2*k : 2*j;
    }
  return fp;
}


int free_fft_plan_bt(fft_plan *fp)
{
  if (fp == NULL) return 1;
  if (fp->mixbt != NULL) free(fp->mixbt);
  if (fp->fftbtsin != NULL) free(fp->fftbtsin);
  free(fp);
  return 0;
}


int fftbtmixing(fft_plan *fp, float *x)
{
	register int i = 0, j, k;
	float temp;

	for (j=0, i=0; i < fp->nx; j++)
	{

		if ( (k = fp->mixbt[j]) > i)
		{
			temp = x[i];
			x[i++] = x[k];
			x[k++] = temp;
			temp = x[i];
			x[i++] = x[k];
			x[k] = temp;
		}
		else
		{
			i++;
			i++;
		}
	}
	return 0;
}

/*	fftbt()
 *	DESCRIPTION	Caluculates a Fast Fourier Transform on the array x.
 *			A forward transform is calculated for df=1, backward
 *			for df=-1. If the array x is complex data, the data
 *			should be placed in x alternating real and imaginary
 *			terms.
 *
 *	RETURNS		0 on success, non-0 on failure (see fftbt_init() above).
 */
int fftbt (fft_plan *fp, float *x, int df)
{
	register int j, k;
	int  nk, kk, k1, k2, kN;
	float temp0, temp1, temp2, temp3, *x1, *x2, *si1, *si2;


	if (fp == NULL)  return 1;
	fftbtmixing(fp, x);
	nk = 4;
	for ( j = 0; j < fp->nx; j += nk)
	{
		x1 = x + j;
		temp0 = x1[0];
		temp1 = x1[1];
		x1[0] += x1[2];
		x1[1] += x1[3];
		x1[2] = temp0 - x1[2];
		x1[3] = temp1 - x1[3];
	}
	nk = 8;
	for (j = 0; j < fp->nx; j += nk)
	{
		x1 = x + j;
		temp0 = x1[0];
		temp1 = x1[1];
		x1[0] += x1[4];
		x1[1] += x1[5];
		x1[4] = temp0 - x1[4];
		x1[5] = temp1 - x1[5];
		temp0 = x1[2];
		temp1 = x1[3];
		x1[6] *= df;
		x1[7] *= df;
		x1[2] += x1[7];
		x1[3] -= x1[6];
		temp2 = x1[6];
		x1[6] = temp0 - x1[7];
		x1[7] = temp1 + temp2;
	}

	for ( nk = 16; nk <= fp->nx; nk=2*nk )
	{
		k1=nk/2;
		k2=nk/4;
		for ( k = 0; k < nk / 8 ;k++)
		{
			kk = fp->nx / k1;
			kk = k * kk;
			kN = fp->n_4 - kk;
			si1 = fp->fftbtsin + kN;
			si2 = fp->fftbtsin + kk;
			for ( j = 0 ; j < fp->nx ; j += nk)
			{
				x1 = x + j + k + k;
				x2 = x1 + k1;
				temp2 = df*si2[0];
				temp3 = - temp2 * x2[0];
				temp2 *= x2[1];
				temp2 += x2[0]*si1[0];
				temp3 += x2[1]*si1[0];
				temp0 = x1[0];
				temp1 = x1[1];
				x1[0] += temp2;
				x1[1] += temp3;
				x2[0] = temp0 - temp2;
				x2[1] = temp1 - temp3;
				x1 += k2;
				x2 += k2;
				temp2 = df*si1[0];
				temp3 = - temp2 * x2[0];
				temp2 *= x2[1];
				temp2 -= x2[0]*si2[0];
				temp3 -= x2[1]*si2[0];
				temp0 = x1[0];
				temp1 = x1[1];
				x1[0] += temp2;
				x1[1] += temp3;
				x2[0] = temp0 - temp2;
				x2[1] = temp1 - temp3;
			}
		}
	}

	if (df ==  1)
	{
		for ( j = 0; j < fp->nx; j++)
			x[j] /= fp->n_2;
	}
	return 0;
}

/*	realtr()
 *	DESCRIPTION	Rearranges the real input data for use by fftbt()
 *			which expects complex data. This routine is not	used
 *			for complex input data. For real input data, it should
 *			be called just before fftbt() for forward	transforms or
 *			just after fftbt() for backward transforms.
 *
 *	RETURNS		Nothing.
 */
void realtr1bt(fft_plan *fp, float *x)
{
	register int i, j;
	float temp;

	if (fp == NULL)  return;
	for (i=1, j=fp->nx-1; i <= fp->n_2; i+=2, j-=2)
	{
		temp = x[i];
		x[i] = x[j];
		x[j] = temp;
	}
}
/*	realtr2()
 *	DESCRIPTION	Converts the complex output of fftbt() to real output.
 *			This routine is not used for complex input data. For
 *			real input data, it should be called just after fftbt()
 *			for forward transforms or just before fftbt() for
 *			backward transforms.
 *
 *	RETURNS		Nothing.
 */
void realtr2bt(fft_plan *fp, float *x, int df)
{
	register int k, nk;
	float temp0, temp1, temp2, temp3, *x1, *x2, *si1, *si2;


	if (fp == NULL)  return;
	if ( df == 1)
	{
		temp0 = (x[0] + x[1]) / 2;
		x[1] = (x[0] - x[1]) / 2;
		x[0] = temp0;
	}
	if ( df ==  - 1)
	{
		temp0 = (x[0] + x[1]) ;
		x[1] = (x[0] - x[1]) ;
		x[0] = temp0;
	}

	for ( k = 1; k <= fp->n_4; k++)
	{
		nk = fp->n_4 - k;
		x1 = x + k + k;
		x2 = x + fp->nx - k - k;
		si1 = fp->fftbtsin + k;
		si2 = fp->fftbtsin + nk;
		temp0 = (x2[0] - x1[0]) / 2;
		temp1 = (x1[1] + x2[1]) / 2;
		temp2 = si1[0] * temp0 + si2[0] * temp1;
		temp3 =  - si2[0] * temp0 + si1[0] * temp1;
		temp0 = (x2[0] + x1[0]) / 2;
		temp1 = (x1[1] - x2[1]) / 2;
		x1[0] = temp0 + temp2;
		x1[1] = temp1 + temp3;
		x2[0] = temp0 - temp2;
		x2[1] =  - temp1 + temp3;
	}
}

/*	fftbtwindow()
 *	DESCRIPTION	Performs a standard cosine window on an entire real
 *			data set. The window does not work for complex data.
 *			fftbt_init() need not be called before this routine
 *			and the sine table is filled by it.
 *
 *	RETURNS		0 on success, non-0 on failure (see fftbt_init()).
 */
int fftbtwindow(fft_plan *fp, float *x)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	for (i=0, j=fp->n_4; i<fp->n_4; i++, j--)
		x[i] *= (1.0 - fp->fftbtsin[j]);
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++)
		x[i] *= (1.0 + fp->fftbtsin[j]);
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, j--)
		x[i] *= (1.0 + fp->fftbtsin[j]);
	for (i=nk, j=0; i < fp->nx ; i++, j++)
		x[i] *= (1.0 - fp->fftbtsin[j]);
	return 0;
}
/*	defftbtwindow()
 *	DESCRIPTION	Performs a standard cosine dewindow on an real
 *			data set excepted point 0. The window does not work for complex data.
 *			fftbt_init() need not be called before this routine
 *			and the sine table is filled by it.
 *
 *	RETURNS		0 on success, non-0 on failure (see fftbt_init()).
 */
int defftbtwindow(fft_plan *fp, float *x)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	for (i=1, j=fp->n_4-1; i<fp->n_4; i++, j--)
		x[i] /= (1.0 - fp->fftbtsin[j]);
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++)
		x[i] /= (1.0 + fp->fftbtsin[j]);
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, j--)
		x[i] /= (1.0 + fp->fftbtsin[j]);
	for (i=nk, j=0; i < fp->nx ; i++, j++)
		x[i] /= (1.0 - fp->fftbtsin[j]);
	return 0;
}


int fftbtwindow1(fft_plan *fp, float *x, float smp)
{
	register int i, j;
	int nk;

	if (fp == NULL)  return 1;
	for (i=0, j=fp->n_4; i<fp->n_4; i++, j--)
		x[i] *= (1.0+smp - fp->fftbtsin[j]);
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++)
		x[i] *= (1.0+smp + fp->fftbtsin[j]);
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, j--)
		x[i] *= (1.0+smp + fp->fftbtsin[j]);
	for (i=nk, j=0; i < fp->nx ; i++, j++)
		x[i] *= (1.0+smp - fp->fftbtsin[j]);
	return 0;
}

int defftbtwindow1(fft_plan *fp, float *x, float smp)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	for (i=0, j=fp->n_4; i<fp->n_4; i++, j--)
		x[i] /= (1.0+smp - fp->fftbtsin[j]);
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++)
		x[i] /= (1.0+smp + fp->fftbtsin[j]);
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, j--)
		x[i] /= (1.0+smp + fp->fftbtsin[j]);
	for (i=nk, j=0; i < fp->nx ; i++, j++)
		x[i] /= (1.0+smp - fp->fftbtsin[j]);
	return 0;
}

int fftbtwc(fft_plan *fp, float *x)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	for (i=0, j=fp->n_4; i<fp->n_4; i++,i++, j--,j--)
	{
		x[i] *= (1.0 - fp->fftbtsin[j]);
		x[i+1] *= (1.0 - fp->fftbtsin[j]);
	}
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++, i++, j++)
	{
		x[i] *= (1.0 + fp->fftbtsin[j]);
		x[i+1] *= (1.0 + fp->fftbtsin[j]);
	}
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] *= (1.0 + fp->fftbtsin[j]);
		x[i+1] *= (1.0 + fp->fftbtsin[j]);
	}
	for (i=nk, j=0; i < fp->nx ; i++, i++, j++, j++)
	{
		x[i] *= (1.0 - fp->fftbtsin[j]);
		x[i+1] *= (1.0 - fp->fftbtsin[j]);
	}
	return 0;
}
int fftbtwc1(fft_plan *fp, float *x, float smp)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	for (i=0, j=fp->n_4; i<fp->n_4; i++,i++, j--,j--)
	{
		x[i] *= (1.0 + smp - fp->fftbtsin[j]);
		x[i+1] *= (1.0 + smp - fp->fftbtsin[j]);
	}
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++, i++, j++)
	{
		x[i] *= (1.0 + smp + fp->fftbtsin[j]);
		x[i+1] *= (1.0 + smp + fp->fftbtsin[j]);
	}
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] *= (1.0 + smp + fp->fftbtsin[j]);
		x[i+1] *= (1.0 + smp + fp->fftbtsin[j]);
	}
	for (i=nk, j=0; i < fp->nx ; i++, i++, j++, j++)
	{
		x[i] *= (1.0 + smp - fp->fftbtsin[j]);
		x[i+1] *= (1.0 + smp - fp->fftbtsin[j]);
	}
	return 0;
}

int defftbtwc(fft_plan *fp, float *x)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	x[0]=0;
	x[1]=1;
	for (i=2, j=fp->n_4 - 2; i<fp->n_4; i++,i++, j--,j--)
	{
		x[i] /= (1.0 - fp->fftbtsin[j]);
		x[i+1] /= (1.0 - fp->fftbtsin[j]);
	}
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++, i++, j++)
	{
		x[i] /= (1.0 + fp->fftbtsin[j]);
		x[i+1] /= (1.0 + fp->fftbtsin[j]);
	}
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] /= (1.0 + fp->fftbtsin[j]);
		x[i+1] /= (1.0 + fp->fftbtsin[j]);
	}
	for (i=nk, j=0; i < fp->nx ; i++, i++, j++, j++)
	{
		x[i] /= (1.0 - fp->fftbtsin[j]);
		x[i+1] /= (1.0 - fp->fftbtsin[j]);
	}
	return 0;
}

int defftbtwc1(fft_plan *fp, float *x, float smp)
{
	register int i, j;
	int  nk;

	if (fp == NULL)  return 1;
	for (i=0, j=fp->n_4 ; i<fp->n_4; i++,i++, j--,j--)
	{
		x[i] /= (1.0+smp - fp->fftbtsin[j]);
		x[i+1] /= (1.0+smp - fp->fftbtsin[j]);
	}
	for (i=fp->n_4, j=0; i < fp->n_2 ; i++, j++, i++, j++)
	{
		x[i] /= (1.0+smp + fp->fftbtsin[j]);
		x[i+1] /= (1.0+smp + fp->fftbtsin[j]);
	}
	nk = fp->n_2 + fp->n_4;
	for (i=fp->n_2, j=fp->n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] /= (1.0+smp + fp->fftbtsin[j]);
		x[i+1] /= (1.0+smp + fp->fftbtsin[j]);
	}
	for (i=nk, j=0; i < fp->nx ; i++, i++, j++, j++)
	{
		x[i] /= (1.0+smp - fp->fftbtsin[j]);
		x[i+1] /= (1.0+smp - fp->fftbtsin[j]);
	}
	return 0;
}
/*	spec_real()
 *
 *	Compute the spectrum for the Fourier transform of real data.
 *	If data was npts points then the spectrum has npts/2 + 1 points
 *
 *	RETURNS		nothing
 */
void spec_realbt (fft_plan *fp, float *x, float *spe)
{
	register int i, k;
	float tmp;

	if (fp == NULL)  return;
	spe [0] =  x[0]*x[0];
	tmp = x[1]*x[1];
	for (i = 1; i < fp->n_2 ; i++)
	{
		k = 2*i;
		spe[i] = x[k] *x[k] + x[k+1]*x[k+1];
	}
	spe [fp->n_2] = tmp;
}
/*	spec_comp()
 *
 *	Compute the spectrum for the Fourier transform of complex data.
 *	If data was npts/2 complex points then the spectrum has npts/2
 *	real points
 *
 *	RETURNS		nothing
 */
void spec_compbt (fft_plan *fp, float *x, float *spe)
{
	register int i, k;
	float *xt,*spet;

	if (fp == NULL)  return;
	xt = x + fp->n_2;
	spet = spe + fp->n_2/2;

	for (i = 0; i < fp->n_2/2 ; i++)
	{
		k = 2*i;
		spet[i] = x[k] *x[k];
		spe[i] = xt[k] *xt[k];
		k++;
		spet[i] += x[k] *x[k];
		spe[i] += xt[k] *xt[k];
	}
}
/*	demodulate a spectrum
*	the initial array is the fftbt of real data, the output is the fftbt
*	of complex data, npts is the number of memory boxes neaded
*/
int demodulatebt (fft_plan *fp, float *x, int freq)
{
	register int i;
	float *x1, *x2, tmp;


	if (fp == NULL)  return 1;
	if ( 2*freq > fp->n_2)
		return 3;

	x1 = x + 2*freq;
	x2 = x + fp->nx -2*freq;


	for (i = 2*freq -1; i >= 0 ; i--)
	{
		tmp=x[i];
		x[i]=x1[i];
		x2[i]=-tmp;
		i--;
		tmp=x[i];
		x[i]=x1[i];
		x2[i]=tmp;
	}
	for (i = 2*freq; i < fp->nx -2*freq; i++)
		x[i]=0;
	return 0;
}
/*	derive a spectrum
*	the initial array is the fftbt of real data, the output is the fftbt
*	of real data, npts is the number of memory boxes neaded
*/
void derive_realbt (fft_plan *fp, float *x, float *y)
{
	register int i,k;
	float tmp;

	if (fp == NULL)  return;
	y[0]=0;
	y[1]=0;

	for (i = 2; i < fp->nx ; i++, i++)
	{
		k = i/2;
		tmp=x[i];
		y[i]=k*x[i+1];
		y[i+1]=-k*tmp;
	}
}

void ampbt(fft_plan *fp, float *x, float *y)
{
	register int i,j;
	float tmp;

	if (fp == NULL)  return;

	for (i = 0; i < fp->n_2; i++)
	{
		j=2*i;
		tmp = x[j]*x[j]+x[j+1]*x[j+1];
		y[i] = sqrt (tmp);
	}
}
# endif
