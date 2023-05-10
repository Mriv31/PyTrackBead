#include "track_util.h"


int rm_low_modes_in_auto_convolution = 0;
int noise_from_auto_conv = 2;
int auto_conv_norm = 0;
int filter_after_conv = 1;

int deplace(float *x, int nx)
{
    int i, j;
    float tmp;

    for (i = 0, j = nx / 2; j < nx ; i++, j++)
    {
        tmp = x[j];
        x[j] = x[i];
        x[i] = tmp;
    }

    return 0;
}

int check_bead_not_lost_diff(int *x, int *xd, float *x1, int *y, int *yd, float *y1, int cl, float bd_mul, int do_diff)
{
    int i, j;
    int minx, maxx, miny, maxy;
    int minxd, maxxd, minyd, maxyd;

    minxd = maxxd = x[0] - xd[0];
    minyd = maxyd = y[0] - yd[0];

    for (i = 0, minx = maxx = x[0], miny = maxy = y[0]; i < cl; i++)
    {
        j = x[i];
        minx = (j < minx) ? j : minx;
        maxx = (j > maxx) ? j : maxx;

        if (do_diff == 1)
        {
            j -= xd[i];
            minxd = (j < minxd) ? j : minxd;
            maxxd = (j > maxxd) ? j : maxxd;
        }

        x1[i] = (float)j;
        j = y[i];
        miny = (j < miny) ? j : miny;
        maxy = (j > maxy) ? j : maxy;

        if (do_diff == 1)
        {
            j -= yd[i];
            minyd = (j < minyd) ? j : minyd;
            maxyd = (j > maxyd) ? j : maxyd;
        }

        y1[i] = (float)j;
    }

    //test image contrast, if it is too low, the bead is lost
    if (do_diff == 0)
        return !!(((maxx - minx) * bd_mul) < (maxx + minx) && ((maxy - miny) * bd_mul) < (maxy + miny));
    else
        return !!(((maxxd - minxd) * bd_mul) < (maxx + minx) && ((maxyd - minyd) * bd_mul) < (maxy + miny));
}

int fftwindow_flat_top1(fft_plan *fp, int npts, float *x, float smp)
{
    int i, j;
    //int n_2, n_4;
    float tmp; //, *offtsin;

    (void)npts;

    if (fp == NULL)  return 1;

    //offtsin = get_fftsin();

    //if ((j = fft_init(npts)) != 0)  return j;

    //n_2 = npts/2;
    //n_4 = n_2/2;
    x[0] *= smp;
    smp += 1.0;

    for (i = 1, j = fp->n_4 - 2; j >= 0; i++, j -= 2)
    {
        tmp = (smp - fp->fftbtsin[j]);
        x[i] *= tmp;
        x[fp->nx - i] *= tmp;
    }

    for (i = fp->n_4 / 2, j = 0; i < fp->n_4 ; i++, j += 2)
    {
        tmp = (smp + fp->fftbtsin[j]);
        x[i] *= tmp;
        x[fp->nx - i] *= tmp;
    }

    return 0;
}

int find_max1(float *x, int nx, float *Max_pos, float *Max_val, float *Max_deriv, int max_limit)
{
    register int i;
    int  nmax, ret = 0, delta, imin, imax;
    double a, b, c, d, f_2, f_1, f_0, f1, f2, xm = 0;
    int maxs_pos[PROFILE_BUFFER_SIZE];
    int nmaxs = 0;

    if (max_limit < 2) max_limit = nx / 2;

    imin = nx / 2 - max_limit;
    imin = (imin < 0) ? 0 : imin;
    imax = nx / 2 + max_limit;
    imax = (imax >= nx) ? nx : imax;

    maxs_pos[0] = 0;
    // find maximum, handle several points with the same maxima value.
    // TODO : handle correctly the case where maximums are not neighbours
    // find the first maximum
    for (i = imin, nmaxs = 1; i < imax; i++)
    {
      maxs_pos[0] = (i == imin || x[i] > x[maxs_pos[0]]) ? i : maxs_pos[0];
    }
    // then find points having the same value (saturation)
    for (i = imin; i < imax; i++)
    {
        if (i == maxs_pos[0]) continue;
        if (x[i] == x[maxs_pos[0]] && nmaxs < PROFILE_BUFFER_SIZE)
        {
            maxs_pos[nmaxs++] = i;
        }
    }
    // take the average position of the maximum
    for (i = 0, nmax = 0; i < nmaxs; ++i)
    {
        nmax += maxs_pos[i];
    }
    nmax /= nmaxs;

    // get neighbourhood
    f_2 = x[(nx + nmax - 2) % nx];
    f_1 = x[(nx + nmax - 1) % nx];
    f_0 = x[nmax];
    f1 = x[(nmax + 1) % nx];
    f2 = x[(nmax + 2) % nx];

    // polinomiale interpolation
    if (f_1 < f1)
    {
        a = -f_1 / 6 + f_0 / 2 - f1 / 2 + f2 / 6;
        b = f_1 / 2 - f_0 + f1 / 2;
        c = -f_1 / 3 - f_0 / 2 + f1 - f2 / 6;
        d = f_0;
        delta = 0;
    }
    else
    {
        a = b = c = 0;
        a = -f_2 / 6 + f_1 / 2 - f_0 / 2 + f1 / 6;
        b = f_2 / 2 - f_1 + f_0 / 2;
        c = -f_2 / 3 - f_1 / 2 + f_0 - f1 / 6;
        d = f_1;
        delta = -1;
    }


    if (fabs(a) < 1e-8)
    {
        if (b != 0)
        {
            xm =  - c / (2 * b) - (3 * a * c * c) / (4 * b * b * b);
        }
        else
        {
            xm = 0;
            ret |= PB_WITH_EXTREMUM;
        }
    }
    else if ((b * b - 3 * a * c) < 0)
    {
        ret |= PB_WITH_EXTREMUM;
        xm = 0;
    }
    else
        xm = (-b - sqrt(b * b - 3 * a * c)) / (3 * a);

    /*
       for (p = -2; p < 1; p++)
       {
       if (6*a*p+b > 0) ret |= TOO_MUCH_NOISE;
       }
       */
    *Max_pos = (float)(xm + nmax + delta);
    *Max_val = (float)((a * xm * xm * xm) + (b * xm * xm) + (c * xm) + d);

    if (Max_deriv)    *Max_deriv = (float)((6 * a * xm) + (2 * b));

    return ret;
}

int correlate_1d_sig_and_invert(float *x1, int nx, int window, float *fxl1, int lfilter, int remove_dc, fft_plan *fp,
                                filter_bt *fil, float *noise)
{
    int i, j;
    int n_2;
    float moy, smp = 0.05, tmpf;
    float m1, m2, ms, md, re1, re2, im1, im2;

    if (x1 == NULL)             return WRONG_ARGUMENT;

    if (fp == NULL || fil == NULL)          return FFT_NOT_SUPPORTED;

    if (fftbt_init(fp, nx) == NULL)         return FFT_NOT_SUPPORTED;

    //if (fft_init(nx) || filter_init(nx))  return FFT_NOT_SUPPORTED;
    if (fxl1 == NULL)   fxl1 = (float *)calloc(nx, sizeof(float));

    if (fxl1 == NULL)           return OUT_MEMORY;

    n_2 = nx / 2;

    /*  compute fft of x1 */
    for (i = 0; i < nx; i++)
        fxl1[i] = x1[i];

    for (i = 0, j = (3 * nx) / 4, moy = 0; i < nx / 4; i++, j++)
        moy += fxl1[i] + fxl1[j];

    moy = (2 * moy) / nx;

    for (i = 0 ; remove_dc && i < nx; i++)
        fxl1[i] = x1[i] - moy;

    if (window == 1)    fftbtwindow1(fp, fxl1, smp);
    else if (window == 2)   fftwindow_flat_top1(fp, nx, fxl1, smp);

    //for(i = 1, dfxl1[0] = fxl1[0] - fxl1[nx-1]; i < nx; i++)      dfxl1[i] = fxl1[i] - fxl1[i-1];


    realtr1bt(fp, fxl1);
    //for(i = 0; i < nx; i++)       x1[i] = fxl1[i];
    //return 0;
    //win_printf("bef fft fxl1[3] = %f moy %f nx %d filter %d",fxl1[3],moy,nx,filter);
    fftbt(fp, fxl1, 1);
    realtr2bt(fp, fxl1, 1);


    if (noise && (noise_from_auto_conv == 0))
    {
        // we estimate noise on convolution by using the high frequency modes
        for (j = 0, i = 3 * nx / 4, m1 = 0; i < nx; i++, j++)
            m1 += fxl1[i] * fxl1[i];

        m1 = (j) ? m1 / j : m1;
        m1 *= n_2;   // we extend the noise to all modes
        *noise = sqrt(m1);
    }

    //we filter low modes
    for (i = 1; i < rm_low_modes_in_auto_convolution; i++)
        fxl1[i] = 0;

    if (lfilter > 0 && filter_after_conv == 0)
        lowpass_smooth_half_bt(fil, nx, fxl1, lfilter);

    /*  compute normalization without DC and Nyquist*/
    for (i = 0, ms = md = 0; i < nx; i += 2)
    {
        tmpf = fxl1[i] * fxl1[i] + fxl1[i + 1] * fxl1[i + 1];
        ms += tmpf;
        md += tmpf * i * i;  // we compute derivative amplitude
    }

    md *= (M_PI * M_PI) / (nx * nx);
    //m1 /= 2;
    m1 = (ms == 0) ? 1 : ms;

    if (auto_conv_norm == 0)      m1 = 1;
    else if (auto_conv_norm == 1) m1 = (float)1 / m1;
    else if (auto_conv_norm == 2) m1 = (float)1 / sqrt(m1);

    /* if m1 = 0, the normalization still needs to multiply by nx/2 to reach \Sigma_i=1^N Si * Si'
       if we assume shot noise in the profile with variance \sigma_n, the error
       on the convolution error is \sigma_c^2 = 4 * \Sigma_i=1^N Si * Si' * \sigma_n^2

       To determine in the error on the position of the maximum, one needs to compute
       the derivative of the maximum and the noise of this derirative, the error on this
       derivatives equals \sigma_d^2 = 4 * \Sigma_i=1^N \Delta Si * \Delta Si' * \sigma_n^2
       without filtering.
       */

    /*  compute correlation in Fourier space */
    for (i = 2; i < nx; i += 2)
    {
        re1 = fxl1[i];
        re2 = fxl1[i];
        im1 = fxl1[i + 1];
        im2 = -fxl1[i + 1];
        fxl1[i]         = (re1 * re2 + im1 * im2) * m1;
        fxl1[i + 1]   = (im1 * re2 - re1 * im2) * m1;
    }

    fxl1[0] = 2 * (fxl1[0] * fxl1[0]) * m1; /* these too mode are special */
    fxl1[1]     = 0;//2*(fxl1[1] * fxl1[1])/m1;   may be to correct with remove_dc

    if (noise && noise_from_auto_conv ==  1)
    {
        // we estimate noise on convolution by using the high frequency modes
        for (j = 0, i = 3 * nx / 4, m1 = 0; i < nx; i++, j++)
            m1 += fxl1[i] * fxl1[i];

        m1 = (j) ? m1 / j : m1;
        m1 *= n_2;   // we extend the noise to all modes
    }

    if (noise && noise_from_auto_conv == 2)
    {
        // we estimate noise using the derivative of the autoconvolution
        if (auto_conv_norm == 0)      m1 = md;
        else if (auto_conv_norm == 1) m1 = md / ms;
        else if (auto_conv_norm == 2) m1 = md / sqrt(ms);

        m1 = 4 * m1 / n_2;
    }

    if (lfilter > 0 && filter_after_conv)
        lowpass_smooth_half_bt(fil, nx, fxl1, lfilter);

    if (noise && noise_from_auto_conv)
    {
        if (lfilter > 0)
        {
            for (i = 0, m2 = 0; i < n_2; i++)
                m2 += fil->f[i];// * fil->f[i];

            m2 /= n_2;
            m1 *= m2;           // we take in account the filter
        }

        *noise = sqrt(m1);
    }

    /*  get back to real world */
    realtr2bt(fp, fxl1, -1);
    //win_printf("bef 2 fft fxl1[2] = %f",fxl1[2]);
    fftbt(fp, fxl1, -1);
    realtr1bt(fp, fxl1);
    deplace(fxl1, nx);
    return 0;
}


float   find_distance_from_center_2(float *x1, float *lfx1, int cl, int flag, int lfilter, int black_circle,
                                    float *corr, float *deriv, fft_plan *fp, filter_bt *fil, float *noise, int max_limit)
{
    float dx = 0;
    (void)black_circle;
    correlate_1d_sig_and_invert(x1, cl, flag, lfx1, lfilter, 1, fp, fil, noise);//(black_circle)?0:1);
    find_max1(lfx1, cl, &dx, corr, deriv, max_limit);
    dx -= cl / 2;
    dx /= 2;
    return dx;
}


int fill_avg_profiles_from_im_diff(O_i *oi, int xc, int yc, int cl, int cw, int *x, int *y, int *xd, int *yd,
                                   int *x_var, int *y_var, int *xd_var, int *yd_var)
{
    int j, i, k;
    unsigned char *choi = NULL;
    short int *inoi = NULL;
    unsigned short int *inui = NULL;
    int ytt, xll, onx, ony, xco, yco, cl2, cw2, cl22;
    unsigned char **pd = NULL;

    if (oi == NULL || x == NULL || y == NULL)   return 1;

    onx = oi->im.nx;
    ony = oi->im.ny;
    pd = oi->im.pixel;
    cl2 = cl >> 1;
    cw2 = cw >> 1;
    cl22 = cl2 * cl2;
    yco = (yc - cl2 < 0) ? cl2 : yc;
    yco = (yco + cl2 <= ony) ? yco : ony - cl2;
    xco = (xc - cl2 < 0) ? cl2 : xc;
    xco = (xco + cl2 <= onx) ? xco : onx - cl2;

    for (i = 0; i < cl; i++)     x[i] = y[i] = xd[i] = yd[i] = 0;

    if (oi->im.data_type == IS_CHAR_IMAGE)
    {
        for (j = 0, ytt = yco - cw2, xll = xco - cl2; j < cw; j++)
        {
            // we average the central band
            choi = pd[ytt + j] + xll;

            for (i = 0; i < cl; i++)        x[i] += (int)choi[i];
        }


        if (x_var)
        {
            for (i = 0; i < cl; i++)       x_var[i] = 0;

            for (j = 0, ytt = yco - cw2, xll = xco - cl2; j < cw; j++)
            {
                for (i = 0, choi = pd[ytt + j] + xll; i < cl; i++)
                    x_var[i] += ((j & 0x01) == (i & 0x01)) ?  - choi[i] : choi[i];
            }
        }

# ifdef PARA_DIFF

        for (j = 0, ytt = yco, xll = xco - cl2; j < cw2; j++)
        {
            // we remove lateral bands on a parabola
            for (i = 0; i < cl; i++)
            {
                k = ((cl2 - cw) * (i - cl2) * (i - cl2)) / cl22;
                k = cl2 - k;
                xd[i] += (int)pd[ytt - k + j][i + xll];
                xd[i] += (int)pd[ytt + k - j - 1][i + xll];
            }
        }

# else

        for (j = 0, ytt = yco - cw, xll = xco - cl2; j < cw2; j++)
        {
            // we remove lateral bands
            choi = pd[ytt + j] + xll;

            for (i = 0; i < cl; i++)        xd[i] += (int)choi[i];
        }

        for (j = 0, ytt = yco + cw2, xll = xco - cl2; j < cw2; j++)
        {
            // we remove lateral bands
            choi = pd[ytt + j] + xll;

            for (i = 0; i < cl; i++)        xd[i] += (int)choi[i];
        }

# endif

        if (xd_var)
        {
            for (i = 0; i < cl; i++)      xd_var[i] = 0;

            for (j = 0, ytt = yco - cw, xll = xco - cl2; j < cw2; j++)
            {
                // we remove lateral bands
                for (i = 0, choi = pd[ytt + j] + xll; i < cl; i++)
                    xd_var[i] += ((j & 0x01) == (i & 0x01)) ? - choi[i] : choi[i];
            }

            for (j = 0, ytt = yco + cw2, xll = xco - cl2; j < cw2; j++)
            {
                // we remove lateral bands
                for (i = 0, choi = pd[ytt + j] + xll; i < cl; i++)
                    xd_var[i] += ((j & 0x01) == (i & 0x01)) ?  - choi[i] : choi[i];
            }
        }


        for (j = 0, ytt = yco - cl2, xll = xco - cw2; j < cl; j++)
        {
            choi = pd[ytt + j] + xll;

            for (i = 0; i < cw; i++)        y[j] += (int)choi[i];
        }


        if (y_var)
        {
            for (i = 0; i < cl; i++)      y_var[i] = 0;

            for (j = 0, ytt = yco - cl2, xll = xco - cw2; j < cl; j++)
            {
                for (i = 0, choi = pd[ytt + j] + xll; i < cw; i++)
                    y_var[j] += ((j & 0x01) == (i & 0x01)) ? -choi[i] : choi[i];
            }
        }

# ifdef PARA_DIFF

        for (j = 0, ytt = yco - cl2, xll = xco; j < cl; j++)
        {
            // we remove lateral bands on a parabola
            k = ((cl2 - cw) * (j - cl2) * (j - cl2)) / cl22;
            k = cl2 - k;

            for (i = 0; i < cw2; i++)
            {
                yd[j] += (int)pd[ytt + j][xll - k + i];
                yd[j] += (int)pd[ytt + j][xll + k - i - 1];
            }
        }

# else

        for (j = 0, ytt = yco - cl2, xll = xco - cw; j < cl; j++)
        {
            choi = pd[ytt + j] + xll;

            for (i = 0; i < cw2; i++)       yd[j] += (int)choi[i];
        }

        for (j = 0, ytt = yco - cl2, xll = xco + cw2; j < cl; j++)
        {
            choi = pd[ytt + j] + xll;

            for (i = 0; i < cw2; i++)       yd[j] += (int)choi[i];
        }

# endif

        if (yd_var)
        {
            for (i = 0; i < cl; i++)           yd_var[i] = 0;

            for (j = 0, ytt = yco - cl2, xll = xco - cw; j < cl; j++)
            {
                for (i = 0, choi = pd[ytt + j] + xll; i < cw2; i++)
                    yd_var[j] += ((j & 0x01) == (i & 0x01)) ? - choi[i] : choi[i];
            }

            for (j = 0, ytt = yco - cl2, xll = xco + cw2; j < cl; j++)
            {
                for (i = 0, choi = pd[ytt + j] + xll; i < cw2; i++)
                    yd_var[j] += ((j & 0x01) == (i & 0x01)) ? - choi[i] : choi[i];
            }
        }

        //my_set_window_title("       short int profile x0 %d",x[0]);
    }

    else
    {
        //my_set_window_title("       not the right profile type x0 %d",x[0]);
        return 1;
    }

    return 0;
}

int fill_X_avg_profiles_from_im_diff(O_i *oi, int xc, int yc, int cl, int cw, int *x, int *y, int *xd, int *yd,
                                     int *x_var, int *y_var, int *xd_var, int *yd_var)

{
    int j, i;
    int ytt, xll, onx, ony, xco, yco, cl2, cw2, d, cw_2, cw_3, ytij, xlij;//, d2
    unsigned char **pd = NULL;

    if (oi == NULL || x == NULL || y == NULL) return 1;

    onx = oi->im.nx;
    ony = oi->im.ny;
    pd = oi->im.pixel;
    cl2 = cl >> 1;
    cw2 = cw >> 1;
    cw_2 = 2 * cw;
    cw_3 = cw + cw2;
    //d2 = cl2 + cw2 - 1;
    d = cl2 + cw - 1;
    yco = (yc - d < 0) ? d : yc;
    yco = (yco + d < ony) ? yco : ony - d - 1;
    xco = (xc - d < 0) ? d : xc;
    xco = (xco + d < onx) ? xco : onx - d - 1;

    for (i = 0; i < cl; i++)       x[i] = y[i] = xd[i] = yd[i] = 0;

    if (oi->im.data_type == IS_CHAR_IMAGE)
    {
        ytt = yco - d;
        xll = xco - d;
        xll += cw_2 - 1;

        for (i = 0; i < cl; i++)
        {
            for (j = 0; j < cw2; j++)
            {
                ytij = ytt + i + j;
                xlij = xll + i - j;
                ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                xd[i] += (int)pd[ytij][xlij];
            }

            for (j = cw2; j < cw_3; j++)
            {
                ytij = ytt + i + j;
                xlij = xll + i - j;
                ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                x[i] += (int)pd[ytij][xlij];
            }

            for (j = cw_3; j < cw_2; j++)
            {
                ytij = ytt + i + j;
                xlij = xll + i - j;
                ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                xd[i] += (int)pd[ytij][xlij];
            }
        }


        if (x_var)
        {
            for (i = 0; i < cl; i++)
            {
                for (j = cw2, x_var[i] = 0; j < cw_3; j++)
                {
                    ytij = ytt + i + j;
                    xlij = xll + i - j;
                    ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                    xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                    x_var[i] += (j % 2) ? - pd[ytij][xlij] : pd[ytij][xlij];
                }
            }
        }

        if (xd_var)
        {
            for (i = 0; i < cl; i++)
            {
                xd_var[i] = 0;

                for (j = 0; j < cw2; j++)
                {
                    ytij = ytt + i + j;
                    xlij = xll + i - j;
                    ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                    xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                    xd_var[i] += (j % 2) ? - pd[ytij][xlij] : pd[ytij][xlij];
                }

                for (j = cw_3; j < cw_2; j++)
                {
                    ytij = ytt + i + j;
                    xlij = xll + i - j;
                    ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                    xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                    xd_var[i] += (j % 2) ? - pd[ytij][xlij] : pd[ytij][xlij];
                }
            }
        }

        ytt = yco - d;
        ytt += cl + cw_2 - 2;
        xll = xco - d;
        xll += cw_2 - 1;

        for (i = 0; i < cl; i++)
        {
            for (j = 0; j < cw2; j++)
            {
                ytij = ytt - i - j;
                xlij = xll + i - j;
                ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                yd[i] += (int)pd[ytij][xlij];
            }

            for (j = cw2; j < cw_3; j++)
            {
                ytij = ytt - i - j;
                xlij = xll + i - j;
                ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                y[i] += (int)pd[ytij][xlij];
            }

            for (j = cw_3; j < cw_2; j++)
            {
                ytij = ytt - i - j;
                xlij = xll + i - j;
                ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                yd[i] += (int)pd[ytij][xlij];
            }
        }


        if (y_var)
        {
            for (i = 0; i < cl; i++)
            {
                for (j = cw2, y_var[i] = 0; j < cw_3; j++)
                {
                    ytij = ytt - i - j;
                    xlij = xll + i - j;
                    ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                    xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                    y_var[i] += (j % 2) ? - pd[ytij][xlij] : pd[ytij][xlij];
                }
            }
        }

        if (yd_var)
        {
            for (i = 0; i < cl; i++)
            {
                yd_var[i] = 0;

                for (j = 0; j < cw2; j++)
                {
                    ytij = ytt - i - j;
                    xlij = xll + i - j;
                    ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                    xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                    yd_var[i] += (j % 2) ? - pd[ytij][xlij] : pd[ytij][xlij];
                }

                for (j = cw_3; j < cw_2; j++)
                {
                    ytij = ytt - i - j;
                    xlij = xll + i - j;
                    ytij = (ytij < 0) ? 0 : ((ytij >= oi->im.ny) ? oi->im.ny - 1 : ytij);
                    xlij = (xlij < 0) ? 0 : ((xlij >= oi->im.nx) ? oi->im.nx - 1 : xlij);
                    yd_var[i] += (j % 2) ? - pd[ytij][xlij] : pd[ytij][xlij];
                }
            }
        }


    }
    else return 1;

    return 0;
}
