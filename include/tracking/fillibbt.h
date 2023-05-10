/*	fillibbt.h
 *
 *	5 fev 92	VC
 *			revision JMF
 *
 *	header file for fillib.c
 *		Collection of functions performing lowpass, highpass
 *			and bandpass filtering on the fourier transform of
 *			either real or complex data. Modified for multi threads
 */
#pragma once

typedef struct _filter
{
	int n;		/*size of the filter array in use*/
	int m;		/*max. size of the array, n<m */
	int hp;		/*high-pass freq.*/
	int wh;		/*high-pass width*/
	int lp;		/*low-pass freq.*/
	int wl;		/*low-pass width*/
	int fl;		/*flag 0(1) indicating real(complex) data */
	float *f;  	/*pointer to the array of filter_bt */
} filter_bt;


#ifdef __cplusplus
extern "C" {
#endif

filter_bt* filter_init_bt(filter_bt *initial, int npts);


#ifdef __cplusplus
}
#endif

int free_filter_bt(filter_bt *fil);
int lowpass_smooth_half_bt(filter_bt *ft, int npts, float *x, int cutoff);
int hipass_smooth_half_bt(filter_bt *ft, int npts, float *x, int cutoff);
int bandpass_smooth_half_bt(filter_bt *ft, int npts, float *x, int center, int width);
int lowpass_and_highpass_smooth_half_bt(filter_bt *ft, int npts, float *x, int lp, int lw, int hp, int hw);
int lowpass_smooth_sym_bt(filter_bt *ft, int npts, float *x, int cutoff);
int hipass_smooth_sym_bt(filter_bt *ft, int npts, float *x, int cutoff);
int bandpass_smooth_sym_bt(filter_bt *ft, int npts, float *x, int center, int width);
int lowpass_smooth_dissym_bt(filter_bt *ft, int npts, float *x, int cutoff);
int hipass_smooth_dissym_bt(filter_bt *ft, int npts, float *x, int cutoff);
int bandpass_smooth_dissym_bt(filter_bt *ft, int npts, float *x, int center, int width);
char* get_filer_error_bt(filter_bt *ft, int type);

# define	MALLOC_ERROR			1
# define	HP_TOO_HIGH			2
# define	HP_W_TOO_HIGH			4
# define	HP_W_TOO_LOW			8
# define	LP_TOO_HIGH			16
# define	LP_W_TOO_HIGH			32
# define	LP_W_TOO_LOW			64
# define	NEGATIVE_FREQ_IN_REAL_DATA	128
