#ifndef _FFTL32N_C_
#define _FFTL32N_C_

/*   fftl32n.c
 *
 *	5 fev. 92	VC
 *	5 fev. 92	revision JMF
 *	20 mars 92	rev JMF function defftwindow1
 *
 *	A collection of routines for performing Fast Fourier Transforms on
 *	real and complex data. In all routines, npts is the number of
 *	elements in the data set, so for complex input npts is twice the
 *	number of complex data pairs. x is an array of floats which is the
 *	time series on input and is replaced by the fft output. The length
 *	of the array is npts. For complex input data, the real and imaginary
 *	parts are alternated in the array, i.e.
 *	x[0] = re{z[0]}, x[1] = im{z[0]}, x[2] = re{z[1]}, x[3] = im{z[1]} ...
 *	The same is true on output (excepted for z[1] in the case of real
 *	data).
 */
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include "fftl32n.h"

/* # include "fftl32n.i"*/

float 	*fftsin = NULL;
int 	*mix = NULL;

# ifdef FORTIFY
# include "fortify.h"
# endif

static int oldnpts = 0, nptsmax = 0;

/*	fft_init()
 *	DESCRIPTION	Calculates the array of sines that will be needed
 *			for the fft. This routine is called by fft(). The
 *			user may call this routine explicitly before calling
 *			fft() for the sake of using the array fftsin[] for
 *			windowing prior to calculating the fft.
 *
 *	RETURNS		0 on success, 1 on npts which is not a power of two,
 *			2 on array memory could not be allocated.
 */
int fft_init(int npts)
{
	register int i, k;
	int n_2, n_4, n=1, nbits, j;
	static int pow_2[31] = {0};

#ifdef DEBUG
	fprintf(stderr,"Initialisation: fftsin[i]\n");
#endif

	/* Check npts is a power of 2 */
	for (i=0; i<31; i++)
	{
		if (npts == n)
			break;
		n *= 2;
	}
	if (i == 31)
		return 1;

	n_2 = npts / 2;
	n_4 = n_2 / 2;
	if ( oldnpts != 0 && npts > nptsmax)
	{
		if (fftsin != NULL)	free (fftsin);
		if (mix != NULL)	free (mix);
		fftsin = NULL; mix = NULL;
	}
	if ( npts > nptsmax)
	{
		fftsin=(float *)malloc((n_4+1)*sizeof(float));
		mix = (int *)malloc(n_2*sizeof(int));
		if (fftsin == NULL || mix == NULL)
			return 2;
		nptsmax = npts;
	}

	for (i = 0 ; i <= n_4 ; i ++)
	{
		fftsin[i] =  sin ((M_PI * i) / n_2);
	}
	oldnpts = npts;
	pow_2[0] = n_4;
	i=0;
	do
	{
		i++;
		pow_2[i] = pow_2[i - 1] >> 1;
	} while (pow_2[i] != 1 );
	nbits = i;

	for (j=0; j < n_2; j++)
	{
		k=0;
		for (i=0; i <= nbits; i++)
		{
			if ((j & pow_2[i] ) > 0)
			{
				k += pow_2[nbits-i];
			}
		}
		if ( k > j)
			mix[j] = 2*k;
		else
			mix[j] = 2*j;
	}
	return 0;
}

int fftmixing(int npts, float *x)
{
	register int i=0, k;
	int j;
	float temp;

#ifdef DEBUG
	fprintf (stderr,"Mixing()...\n");
	fprintf(stderr,"Initialization: pow_2[i]\n");
#endif


	for (j=0, i=0; i < npts; j++)
	{

		if ( (k = mix[j]) > i)
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

/*	fft()
 *	DESCRIPTION	Caluculates a Fast Fourier Transform on the array x.
 *			A forward transform is calculated for df=1, backward
 *			for df=-1. If the array x is complex data, the data
 *			should be placed in x alternating real and imaginary
 *			terms.
 *
 *	RETURNS		0 on success, non-0 on failure (see fft_init() above).
 */
int fft (int npts, float *x, int df)
{
	register int j, k;
	int n_2, n_4, nk, kk, k1, k2, kN;
	float temp0, temp1, temp2, temp3, *x1, *x2, *si1, *si2;

#ifdef DEBUG
	fprintf (stderr, "FFT()...\n");
#endif

	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}
	fftmixing(npts, x);
	n_2 = npts/2;
	n_4 = n_2/2;
	nk = 4;
	for ( j = 0; j < npts; j += nk)
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
	for (j = 0; j < npts; j += nk)
	{
		x1 = x + j;
		temp0 = x1[0];
		temp1 = x1[1];
		x1[0] += x1[4];         //x1[0] = x1[0] + x1[4]
		x1[1] += x1[5];         //x1[1] = x1[1] + x1[5]
		x1[4] = temp0 - x1[4];  //x1[4] = x1[0] - x1[4]
		x1[5] = temp1 - x1[5];  //x1[5] = x1[1] - x1[5]
		temp0 = x1[2];
		temp1 = x1[3];
		x1[6] *= df;
		x1[7] *= df;
		x1[2] += x1[7];         //x1[2] = x1[2] + df * x1[7]
		x1[3] -= x1[6];         //x1[3] = x1[3] - df * x1[6]
		temp2 = x1[6];
		x1[6] = temp0 - x1[7];  //x1[6] = x1[2] - df * x1[7]
		x1[7] = temp1 + temp2;  //x1[7] = x1[3] + df * x1[6]
	}

	for ( nk = 16; nk <= npts; nk=2*nk )
	{
		k1=nk/2;
		k2=nk/4;
		for ( k = 0; k < nk / 8 ;k++)
		{
			kk = npts / k1;
			kk = k * kk;
			kN = n_4 - kk;
			si1 = fftsin + kN;
			si2 = fftsin + kk;
			for ( j = 0 ; j < npts ; j += nk)
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
		    //x1[0] = x1[0] + x2[0] * si1[0] + x2[1] * df * si2[0];
		    //x1[1] = x1[1] + x2[1] * si1[0] - x2[0] * df * si2[0];
		    //x2[0] = x1[0] - x2[0] * si1[0] - x2[1] * df * si2[0];
		    //x2[1] = x1[1] - x2[1] * si1[0] + x2[0] * df * si2[0];
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
		    //x1[0] = x1[0] - x2[0] * si2[0] + x2[1] * df * si1[0];
		    //x1[1] = x1[1] - x2[1] * si2[0] - x2[0] * df * si1[0];
		    //x2[0] = x1[0] + x2[0] * si1[0] - x2[1] * df * si1[0];
		    //x2[1] = x1[1] - x2[1] * si1[0] + x2[0] * df * si1[0];
			}
		}
	}

	if (df ==  1 )
	{
		for ( j = 0; j < npts; j++)
			x[j] /= n_2;
	}
	return 0;
}

/*	realtr()
 *	DESCRIPTION	Rearranges the real input data for use by fft()
 *			which expects complex data. This routine is not	used
 *			for complex input data. For real input data, it should
 *			be called just before fft() for forward	transforms or
 *			just after fft() for backward transforms.
 *
 *	RETURNS		Nothing.
 */
void realtr1(int npts, float *x)
{
	register int i, j;
	int n_2;
	float temp;

#ifdef DEBUG
	fprintf (stderr,"realtr1()...\n");
#endif

	n_2 = npts/2;
	for (i=1, j=npts-1; i <= n_2; i+=2, j-=2)
	{
		temp = x[i];
		x[i] = x[j];
		x[j] = temp;
	}
}
/*	realtr2()
 *	DESCRIPTION	Converts the complex output of fft() to real output.
 *			This routine is not used for complex input data. For
 *			real input data, it should be called just after fft()
 *			for forward transforms or just before fft() for
 *			backward transforms.
 *
 *	RETURNS		Nothing.
 */
void realtr2(int npts, float *x, int df)
{
	register int k, n_4;
	int nk;
	float temp0, temp1, temp2, temp3, *x1, *x2, *si1, *si2;

#ifdef DEBUG
	fprintf (stderr,"realtr2()\n");
#endif

	n_4 = npts/4;
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

	for ( k = 1; k <= n_4; k++)
	{
		nk = n_4 - k;
		x1 = x + k + k;
		x2 = x + npts - k - k;
		si1 = fftsin + k;
		si2 = fftsin + nk;
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

/*	fftwindow()
 *	DESCRIPTION	Performs a standard cosine window on an entire real
 *			data set. The window does not work for complex data.
 *			fft_init() need not be called before this routine
 *			and the sine table is filled by it.
 *
 *	RETURNS		0 on success, non-0 on failure (see fft_init()).
 */
int fftwindow(int npts, float *x)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Windowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;

	for (i=0, j=n_4; i<n_4; i++, j--)
	{
		x[i] *= (1.0 - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++)
	{
		x[i] *= (1.0 + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, j--)
	{
		x[i] *= (1.0 + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, j++)
	{
		x[i] *= (1.0 - fftsin[j]);
	}
	return 0;
}


int fftwindow1(int npts, float *x, float smp)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Windowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;

	for (i=0, j=n_4; i<n_4; i++, j--)
	{
		x[i] *= (1.0+smp - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++)
	{
		x[i] *= (1.0+smp + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, j--)
	{
		x[i] *= (1.0+smp + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, j++)
	{
		x[i] *= (1.0+smp - fftsin[j]);
	}
	return 0;
}

int defftwindow1(int npts, float *x, float smp)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Dewindowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;

	for (i=0, j=n_4; i<n_4; i++, j--)
	{
		x[i] /= (1.0+smp - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++)
	{
		x[i] /= (1.0+smp + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, j--)
	{
		x[i] /= (1.0+smp + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, j++)
	{
		x[i] /= (1.0+smp - fftsin[j]);
	}
	return 0;
}

int fftwc(int npts, float *x)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Complex Windowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;

	for (i=0, j=n_4; i<n_4; i++,i++, j--,j--)
	{
		x[i] *= (1.0 - fftsin[j]);
		x[i+1] *= (1.0 - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++, i++, j++)
	{
		x[i] *= (1.0 + fftsin[j]);
		x[i+1] *= (1.0 + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] *= (1.0 + fftsin[j]);
		x[i+1] *= (1.0 + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, i++, j++, j++)
	{
		x[i] *= (1.0 - fftsin[j]);
		x[i+1] *= (1.0 - fftsin[j]);
	}
	return 0;
}
int fftwc1(int npts, float *x, float smp)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Complex Windowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;

	for (i=0, j=n_4; i<n_4; i++,i++, j--,j--)
	{
		x[i] *= (1.0 + smp - fftsin[j]);
		x[i+1] *= (1.0 + smp - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++, i++, j++)
	{
		x[i] *= (1.0 + smp + fftsin[j]);
		x[i+1] *= (1.0 + smp + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] *= (1.0 + smp + fftsin[j]);
		x[i+1] *= (1.0 + smp + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, i++, j++, j++)
	{
		x[i] *= (1.0 + smp - fftsin[j]);
		x[i+1] *= (1.0 + smp - fftsin[j]);
	}
	return 0;
}

int defftwc(int npts, float *x)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Complex De_Windowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;
	x[0]=0;
	x[1]=1;
	for (i=2, j=n_4 - 2; i<n_4; i++,i++, j--,j--)
	{
		x[i] /= (1.0 - fftsin[j]);
		x[i+1] /= (1.0 - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++, i++, j++)
	{
		x[i] /= (1.0 + fftsin[j]);
		x[i+1] /= (1.0 + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] /= (1.0 + fftsin[j]);
		x[i+1] /= (1.0 + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, i++, j++, j++)
	{
		x[i] /= (1.0 - fftsin[j]);
		x[i+1] /= (1.0 - fftsin[j]);
	}
	return 0;
}

int defftwc1(int npts, float *x, float smp)
{
	register int i, j;
	int n_2, n_4, nk;

#ifdef DEBUG
	fprintf(stderr,"Complex De_Windowing...\n");
#endif
	if (npts != oldnpts)
	{
		if ((j = fft_init(npts)) != 0)
			return j;
		oldnpts = npts;
	}

	n_2 = npts/2;
	n_4 = n_2/2;

	for (i=0, j=n_4 ; i<n_4; i++,i++, j--,j--)
	{
		x[i] /= (1.0+smp - fftsin[j]);
		x[i+1] /= (1.0+smp - fftsin[j]);
	}
	for (i=n_4, j=0; i < n_2 ; i++, j++, i++, j++)
	{
		x[i] /= (1.0+smp + fftsin[j]);
		x[i+1] /= (1.0+smp + fftsin[j]);
	}
	nk = n_2 + n_4;
	for (i=n_2, j=n_4; i < nk ; i++, i++, j--, j--)
	{
		x[i] /= (1.0+smp + fftsin[j]);
		x[i+1] /= (1.0+smp + fftsin[j]);
	}
	for (i=nk, j=0; i < npts ; i++, i++, j++, j++)
	{
		x[i] /= (1.0+smp - fftsin[j]);
		x[i+1] /= (1.0+smp - fftsin[j]);
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
void spec_real (int npts, float *x, float *spe)
{
	register int i, k;
	int n_2;
	float tmp;

	n_2 = npts/2;
	spe [0] =  x[0]*x[0];
	tmp = x[1]*x[1];
	for (i = 1; i < n_2 ; i++)
	{
		k = 2*i;
		spe[i] = x[k] *x[k] + x[k+1]*x[k+1];
	}
	spe [n_2] = tmp;
}
/*	spec_comp()
 *
 *	Compute the spectrum for the Fourier transform of complex data.
 *	If data was npts/2 complex points then the spectrum has npts/2
 *	real points
 *
 *	RETURNS		nothing
 */
void spec_comp (int npts, float *x, float *spe)
{
	register int i, k;
	int n_2;
	float *xt = NULL,*spet = NULL;

	n_2 = npts/2;
	xt = x + n_2;
	spet = spe + n_2/2;

	for (i = 0; i < n_2/2 ; i++)
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
*	the initial array is the fft of real data, the output is the fft
*	of complex data, npts is the number of memory boxes neaded
*/
int demodulate (int npts, float *x, int freq)
{
	register int i;
	int n_2;
	float *x1 = NULL, *x2 = NULL, tmp;


	n_2 = npts/2;
	if ( 2*freq > n_2)
		return 3;

	x1 = x + 2*freq;
	x2 = x +npts -2*freq;


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
	for (i = 2*freq; i < npts -2*freq; i++)
		x[i]=0;
	return 0;
}
/*	derive a spectrum
*	the initial array is the fft of real data, the output is the fft
*	of real data, npts is the number of memory boxes neaded
*/
void derive_real (int npts, float *x, float *y)
{
	register int i,k;
	float tmp;

	y[0]=0;
	y[1]=0;


	for (i = 2; i <npts ; i++, i++)
	{
		k = i/2;
		tmp=x[i];
		y[i]=k*x[i+1];
		y[i+1]=-k*tmp;
	}
}

void amp(int npts, float *x, float *y)
{
	register int i,j;
	int n_2;
	float tmp;

	n_2 = npts /2;

	for (i = 0; i < n_2; i++)
	{
		j=2*i;
		tmp = x[j]*x[j]+x[j+1]*x[j+1];
		y[i] = sqrt (tmp);
	}
}
# endif
