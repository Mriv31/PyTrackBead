/*  fillib.c
 *
 *  5 fev 92    VC
 *          revision JMF
 *
 *  DESCRIPTION Collection of function performing lowpass, highpass
 *          and bandpass filtering on the fourier transform of
 *          either real or complex data.
 *
 *
 *  filter_init(npts)
 *
 *  lowpass_smooth_half (npts, x, cutoff)
 *  hipass_smooth_half (npts, x, cutoff)
 *  bandpass_smooth_half (npts, x, center, width)
 *  (for filtering Fourier transform data of real serie)
 *
 *  lowpass_smooth_sym (npts, x, cutoff)
 *  hipass_smooth_sym (npts, x, cutoff)
 *  bandpass_smooth_sym (npts, x, center, width)
 *  (for filtering Fourier transform data of complex serie)
 *
 *  npts is the number of memories required to store the Fourier transform
 *  corresponding to npts points of real data, or npts/2 points of
 *  complex data. npts should be a multiple of 4. npts may be
 *  changed from one call to the other, asking for a smaller
 *  array has no major consequences, increasing npts to a value
 *  larger than anyone previously invoked, reallocates memory and
 *  thus is not not speed efficient. It is wise to call
 *  filter-init with the largest n used later on.
 *
 *
 *  cutoff, center, width characterize the filter shape, their units
 *  are in FFT mode number
 *
 *  When ever one filter is call, the filter parameter are
 *  compared with those stored in memory, if they differ, the
 *  filter array is recomputed, otherwise the Fourier array is
 *  just multiplied by the filter array.
 *
 *
 */
# include <stdio.h>
# include <math.h>
# include <malloc.h>
# include <stdlib.h>
# include "fillibbt.h"
# ifdef FORTIFY
# include "fortify.h"
# endif

//static struct filter fil;

//static int fillib_init=0;

/*  filter_init(npts)
 *  DESCRIPTION Allocate memory to the float array fil.f[]
 *          used for filtering the data Fourier transformed,
 *          this array is static.
 *
 *  RETURNS     0 on success, 1 on  malloc error
 *
 */

filter_bt *filter_init_bt(filter_bt *initial, int npts)
{
    int n_2 = npts / 2;
    filter_bt *ft = NULL;

    if (initial == NULL)    /*first initialisation*/
    {
        ft = (filter_bt *)calloc(1, sizeof(filter_bt));
        ft->n = 0;
        ft->m = 0;
        ft->lp = 0;
        ft->wl = 0;
        ft->hp = 0;
        ft->wh = 0;
    }
    else
    {
        ft = initial;
    }
    if ( ft->n != n_2 )
    {
        if (ft->n != 0 && ft->m < n_2)
        {
            if (ft->f)
            {
                free (ft->f);
            }
            ft->f = NULL;
        }
        if (ft->m < n_2)
        {
            ft->f = (float *)malloc((n_2 + 1) * sizeof(float));
            if (ft->f == NULL)
            {
                fprintf(stderr, "malloc error n");
                return NULL;
            }
            else
            {
                ft->n = n_2;
                ft->m = n_2;
            }
        }
        else
        {
            ft->n = n_2;
        }
    }
    return ft;
}

int free_filter_bt(filter_bt *fil)
{
    if (fil == NULL)
    {
        return 1;
    }
    if (fil->f != NULL)
    {
        free(fil->f);
    }
    free(fil);
    return 0;
}

/*  lowpass_smooth_half (npts, x, cutoff)
 *  DESCRIPTION Subroutine performing a low-pass filtering on a
 *          Fourier transform of a npts real data points stored
 *          in the x array. The filter is built with cosine,
 *          starts at 1 at zero frequency, 1/2 at cutoff,
 *          and reaches 0 at 2*cutoff.
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */

int lowpass_smooth_half_bt (filter_bt *ft, int npts, float *x, int cutoff)
{
    int i, j;
    int n_4, n_2, c;

    c = cutoff;
    n_2 = npts / 2;
    n_4 = npts / 4;

    if ( cutoff > n_4 || cutoff < 0 || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != 0 || ft->wh != 0 || ft->lp != c || ft->wl != c || ft->fl != 0)
    {
        ft->f[0] = 1;
        for ( i = 1; i < 2 * cutoff ; i++)
        {
            ft->f[i] = (1 + cos ( i * M_PI_2 / cutoff)) / 2;
        }
        for ( i = 2 * cutoff; i < npts / 2 ; i++)
        {
            ft->f[i] = 0;
        }
        ft->f[npts / 2] = 0;
        ft->hp = 0;
        ft->wh = 0;
        ft->lp = cutoff;
        ft->wl = cutoff;
        ft->fl = 0;
    }
    if (x != NULL)
    {
        x[0] *= ft->f[0];
        x[1] *= ft->f[n_2];
        for ( i = 1, j = 2; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}

/*  hipass_smooth_half (npts, x, cutoff)
 *  DESCRIPTION Subroutine performing a high-pass filtering on a
 *          Fourier transform of a npts real data points stored
 *          in the x array. The filter is built with cosine,
 *          starts at 0 at zero frequency, 1/2 at cutoff,
 *          and reaches 1 at 2*cutoff.
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */

int hipass_smooth_half_bt (filter_bt *ft, int npts, float *x, int cutoff)
{
    int i, j;
    int n_4, n_2, c;

    c = cutoff;
    n_2 = npts / 2;
    n_4 = npts / 4;

    if ( cutoff > n_4 || cutoff < 0  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->lp != 0 || ft->wl != 0 || ft->hp != c || ft->wh != c || ft->fl != 0)
    {
        ft->f[0] = 0;
        for ( i = 1; i < 2 * cutoff ; i++)
        {
            ft->f[i] = (1 - cos ( i * M_PI_2 / cutoff)) / 2;
        }
        for ( i = 2 * cutoff; i < n_2 ; i++)
        {
            ft->f[i] = 1;
        }
        ft->f[n_2] = 1;
        ft->lp = 0;
        ft->wl = 0;
        ft->hp = cutoff;
        ft->wh = cutoff;
        ft->fl = 0;
    }
    if (x != NULL)
    {
        x[0] *= ft->f[0];
        x[1] *= ft->f[n_2];
        for ( i = 1, j = 2; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}
/*  bandpass_smooth_half (npts, x, center, width)
 *  DESCRIPTION Subroutine performing a band-pass filtering on a
 *          Fourier transform of a npts real data points stored
 *          in the x array. The filter is built with cosine,
 *          has a gain of 1 at center frequency, falls to 1/2
 *          at center +/- width/2 and reaches 0 at
 *          center +/- width.
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int bandpass_smooth_half_bt (filter_bt *ft, int npts, float *x, int center, int width)
{
    int i, j;
    int n_4, n_2, c1, c2, w;

    n_2 = npts / 2;
    n_4 = npts / 4;
    c1 = center - width;
    c2 = center + width;
    w = width;

    if ( center > n_4 || width >= center || center < 0 || width < 0  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != c1 || ft->wh != w || ft->lp != c2 || ft->wl != w || ft->fl != 0)
    {
        for ( i = 0; i < c1 ; i++)
        {
            ft->f[i] = 0;
        }
        for ( i = 0 ; i < 2 * width ; i++)
        {
            ft->f[i + c1] = (1 - cos ( i * M_PI / width)) / 2;
        }
        for ( i = c2; i < n_2 ; i++)
        {
            ft->f[i] = 0;
        }
        ft->f[n_2] = 0;
        ft->hp = c1;
        ft->wh = width;
        ft->lp = c2;
        ft->wl = width;
        ft->fl = 0;
    }
    if (x != NULL)
    {
        x[0] *= ft->f[0];
        x[1] *= ft->f[n_2];
        for ( i = 1, j = 2; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}



/*  lowpass_and_highpass_smooth_half (npts, x, lp, hp)
 *  DESCRIPTION Subroutine performing a band-pass filtering on a
 *          Fourier transform of a npts real data points stored
 *          in the x array. The filter is built with cosine,
 *          has a gain of 1 at center frequency, falls to 1/2
 *          at center +/- width/2 and reaches 0 at
 *          center +/- width.
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int lowpass_and_highpass_smooth_half_bt (filter_bt *ft, int npts, float *x, int lp, int lw, int hp, int hw)
{
    int i, j;
    int  n_2, c0, c1, c2, c3, w;

    n_2 = npts / 2;
    c0 = lp - lw;
    c1 = lp + lw;
    c2 = hp - hw;
    c3 = hp + hw;



    if ( c2 > n_2 || lw <= 0 || lp <= 0 || hw <= 0 || hp < lp|| ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != hp || ft->wh != hw || ft->lp != lp || ft->wl != lw || ft->fl != 0)
    {
        for (i = 0; i < c0 ; i++)
	  {
	    ft->f[i] = 0;
	  }
        for (i = 0, w = 2*lw; i < w ; i++)
	  {
            ft->f[c0+i] = (1 - cos ( i * M_PI / w)) / 2;
	  }
        for (i = c1; i < c2 ; i++)
	  {
	    ft->f[i] = 1;
	  }
        for (i = 0, w = 2*hw ; i < w ; i++)
	  {
	    ft->f[i + c2] = (1 + cos ( i * M_PI / w)) / 2;
	  }
        for (i = c3; i < n_2 ; i++)
	  {
            ft->f[i] = 0;
	  }
        ft->f[n_2] = 0;
        ft->hp = hp;
        ft->wh = hw;
        ft->lp = lp;
        ft->wl = lw;
        ft->fl = 0;
    }
    if (x != NULL)
    {
        x[0] *= ft->f[0];
        x[1] *= ft->f[n_2];
        for ( i = 1, j = 2; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}




/*  lowpass_smooth_sym (npts, x, cutoff)
 *  DESCRIPTION Subroutine performing a low-pass filtering on a
 *          Fourier transform of a npts complex data points stored
 *          in the x array. The filter is built with cosine,
 *          starts at 1 at zero frequency, 1/2 at +/-cutoff,
 *          and reaches 0 at +/-   2*cutoff.
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int lowpass_smooth_sym_bt (filter_bt *ft, int npts, float *x, int cutoff)
{
    int i, j;
    int n_4, n_2;

    n_2 = npts / 2;
    n_4 = npts / 4;

    if ( 2 * cutoff > n_4  || cutoff < 0  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != 0 || ft->wh != 0 || ft->lp != cutoff || ft->wl != cutoff || ft->fl != 1)
    {
        ft->f[0] = 1;
        for ( i = 1; i < 2 * cutoff ; i++)
        {
            ft->f[i] = (1 + cos ( i * M_PI_2 / cutoff)) / 2;
            ft->f[n_2 - i] = ft->f[i];
        }
        for ( i = 2 * cutoff; i < n_4 ; i++)
        {
            ft->f[i] = 0;
            ft->f[n_2 - i] = 0;
        }
        ft->f[n_4] = 0;
        ft->hp = 0;
        ft->wh = 0;
        ft->lp = cutoff;
        ft->wl = cutoff;
        ft->fl = 1;
    }
    if (x != NULL)
    {
        for ( i = 0, j = 0; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}
/*  hipass_smooth_sym (npts, x, cutoff)
 *  DESCRIPTION Subroutine performing a high-pass filtering on a
 *          Fourier transform of a npts complex data points stored
 *          in the x array. The filter is built with cosine,
 *          starts at 0 at zero frequency, 1/2 at +/-cutoff,
 *          and reaches 1 at +/-2*cutoff.
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int hipass_smooth_sym_bt (filter_bt *ft, int npts, float *x, int cutoff)
{
    int i, j;
    int n_4, n_2;

    n_2 = npts / 2;
    n_4 = npts / 4;

    if ( 2 * cutoff > n_4 || cutoff < 0  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->lp != 0 || ft->wl != 0 || ft->hp != cutoff || ft->wh != cutoff || ft->fl != 1)
    {
        ft->f[0] = 0;
        for ( i = 1; i < 2 * cutoff ; i++)
        {
            ft->f[i] = (1 - cos ( i * M_PI_2 / cutoff)) / 2;
            ft->f[n_2 - i] = (1 - cos ( i * M_PI_2 / cutoff)) / 2;
        }
        for ( i = 2 * cutoff; i < n_4 ; i++)
        {
            ft->f[i] = 1;
            ft->f[n_2 - i] = 1;
        }
        ft->f[n_4] = 1;
        ft->lp = 0;
        ft->wl = 0;
        ft->hp = cutoff;
        ft->wh = cutoff;
        ft->fl = 1;
    }
    if (x != NULL)
    {
        for ( i = 0, j = 0; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}
/*  bandpass_smooth_sym (npts, x, center, width)
 *  DESCRIPTION Subroutine performing a band-pass filtering on a
 *          Fourier transform of a npts complex data points stored
 *          in the x array. The filter is built with cosine,
 *          has a gain of 1 at +/- center frequency, falls to 1/2
 *          at +/-(center +/- width/2) and reaches 0 at
 *          +/-(center +/- width).
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int bandpass_smooth_sym_bt (filter_bt *ft, int npts, float *x, int center, int width)
{
    int i, j;
    int n_4, n_2, c1, c2, w;

    n_2 = npts / 2;
    n_4 = npts / 4;
    c1 = center - width;
    c2 = center + width;
    w = width;

    if ( 2 * center > n_4 || width >= center || center < 0 || width < 0  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != c1 || ft->wh != w || ft->lp != c2 || ft->wl != w || ft->fl != 1)
    {
        ft->f[0] = 0;
        for ( i = 1; i < c1 ; i++)
        {
            ft->f[i] = 0;
            ft->f[n_2 - i] = 0;
        }
        for ( i = 0 ; i < 2 * width ; i++)
        {
            ft->f[i + c1] = (1 - cos ( i * M_PI / width)) / 2;
            ft->f[n_2 - i - c1] = ft->f[i + c1];
        }
        for ( i = c2; i < n_4 ; i++)
        {
            ft->f[i] = 0;
            ft->f[n_2 - i] = 0;
        }
        ft->f[n_4] = 0;
        ft->hp = c1;
        ft->wh = width;
        ft->lp = c2;
        ft->wl = width;
        ft->fl = 1;
    }
    if (x != NULL)
    {
        for ( i = 0, j = 0; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}

/*  lowpass_smooth_dissym (npts, x, cutoff)
 *  DESCRIPTION Subroutine performing a low-pass filtering on a
 *          Fourier transform of a npts complex data points stored
 *          in the x array. The filter is built with cosine,
 *          starts at 1 at zero frequency, 1/2 at cutoff,
 *          and reaches 0 at 2*cutoff. All modes which
 *          sign differs from that of cutoff are set to 0
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int lowpass_smooth_dissym_bt (filter_bt *ft, int npts, float *x, int cutoff)
{
    int i, j;
    int n_4, n_2;

    n_2 = npts / 2;
    n_4 = npts / 4;
    j = abs (cutoff);

    if ( 2 * j > n_4  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != 0 || ft->wh != 0 || ft->lp != cutoff || ft->wl != cutoff || ft->fl != 2)
    {
        ft->f[0] = 1;
        for ( i = 1; i < 2 * j ; i++)
        {
            ft->f[i] = (1 + cos ( i * M_PI_2 / j)) / 2;
            ft->f[n_2 - i] = ft->f[i];
        }
        if (cutoff >= 0 )
        {
            for ( i = 1; i < 2 * j ; i++)
            {
                ft->f[n_2 - i] = 0;
            }
        }
        else
        {
            for ( i = 1; i < 2 * j ; i++)
            {
                ft->f[i] = 0;
            }
        }
        for ( i = 2 * j; i < n_4 ; i++)
        {
            ft->f[i] = 0;
            ft->f[n_2 - i] = 0;
        }
        ft->f[n_4] = 0;
        ft->hp = 0;
        ft->wh = 0;
        ft->lp = cutoff;
        ft->wl = cutoff;
        ft->fl = 2;
    }
    if (x != NULL)
    {
        for ( i = 0, j = 0; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}
/*  hipass_smooth_dissym (npts, x, cutoff)
 *  DESCRIPTION Subroutine performing a high-pass filtering on a
 *          Fourier transform of a npts complex data points stored
 *          in the x array. The filter is built with cosine,
 *          starts at 0 at zero frequency, 1/2 at cutoff,
 *          and reaches 1 at 2*cutoff.  All modes which
 *          sign differs from that of cutoff are set to 0
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int hipass_smooth_dissym_bt (filter_bt *ft, int npts, float *x, int cutoff)
{
    int i, j;
    int n_4, n_2;

    n_2 = npts / 2;
    n_4 = npts / 4;
    j = abs (cutoff);

    if ( 2 * j > n_4  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->lp != 0 || ft->wl != 0 || ft->hp != cutoff || ft->wh != cutoff || ft->fl != 2)
    {
        ft->f[0] = 0;
        for ( i = 1; i < 2 * j ; i++)
        {
            ft->f[i] = (1 - cos ( i * M_PI_2 / j)) / 2;
            ft->f[n_2 - i] = ft->f[i];
        }
        for ( i = 2 * j; i < n_4 ; i++)
        {
            ft->f[i] = 1;
            ft->f[n_2 - i] = 1;
        }
        if (cutoff >= 0 )
        {
            for ( i = 1; i < n_4 ; i++)
            {
                ft->f[n_2 - i] = 0;
            }
        }
        else
        {
            for ( i = 1; i < n_4 ; i++)
            {
                ft->f[i] = 0;
            }
        }
        ft->f[n_4] = 1;
        ft->lp = 0;
        ft->wl = 0;
        ft->hp = cutoff;
        ft->wh = cutoff;
        ft->fl = 2;
    }
    if (x != NULL)
    {
        for ( i = 0, j = 0; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}
/*  bandpass_smooth_dissym (npts, x, center, width)
 *  DESCRIPTION Subroutine performing a band-pass filtering on a
 *          Fourier transform of a npts complex data points stored
 *          in the x array. The filter is built with cosine,
 *          has a gain of 1 at center frequency, falls to 1/2
 *          at (center +/- width/2) and reaches 0 at
 *          (center +/- width).
 *
 *  RETURNS     0 on success, 3 on incompatible cutoff frequency,
 *          1 on  malloc error (in filter_init)
 */
int bandpass_smooth_dissym_bt (filter_bt *ft, int npts, float *x, int center, int width)
{
    int i, j;
    int n_4, n_2, c1, c2, w;

    n_2 = npts / 2;
    n_4 = npts / 4;
    c1 = center - width;
    c2 = center + width;
    w = width;

    if (abs(c1) > n_4 || abs(c2) > n_4 || (abs(center) - abs(width)) < 0 || width < 0  || ft == NULL)
    {
        return 3;
    }

    if (ft->n != n_2)
    {
        if ( filter_init_bt (ft, n_2) == NULL)
        {
            return 1;
        }
        ft->n = n_2;
    }

    if (ft->hp != c1 || ft->wh != w || ft->lp != c2 || ft->wl != w || ft->fl != 2)
    {
        ft->f[0] = 0;
        ft->f[1] = 0;
        c1 = abs(center) - abs(w);
        for ( i = 1; i < c1 ; i++)
        {
            ft->f[i] = 0;
            ft->f[n_2 - i] = 0;
        }
        for ( i = 0 ; i < 2 * width ; i++)
        {
            ft->f[i + c1] = (1 - cos ( i * M_PI / width)) / 2;
            ft->f[n_2 - i - c1] = ft->f[i + c1];
        }
        for ( i = abs(center) + abs(w); i < n_4 ; i++)
        {
            ft->f[i] = 0;
            ft->f[n_2 - i] = 0;
        }
        if (center >= 0) for (i = n_4; i < n_2; i++)
            {
                ft->f[i] = 0;
            }
        else for (i = 2; i < n_4; i++)
            {
                ft->f[i] = 0;
            }
        c1 = center - width;
        ft->f[n_4] = 0;
        ft->hp = c1;
        ft->wh = width;
        ft->lp = c2;
        ft->wl = width;
        ft->fl = 2;
    }
    if (x != NULL)
    {
        for ( i = 0, j = 0; i < n_2; i++)
        {
            x[j++] *= ft->f[i];
            x[j++] *= ft->f[i];
        }
    }
    return 0;
}
