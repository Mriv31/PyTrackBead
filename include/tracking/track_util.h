#pragma once

#include "fillibbt.h"
#include "fftbtl32n.h"

#define FFT_NOT_SUPPORTED 1
#define WRONG_ARGUMENT 2
#define PROFILE_BUFFER_SIZE 3
#define PB_WITH_EXTREMUM 4
#define OUT_MEMORY 5
#define IS_CHAR_IMAGE 19


typedef struct
{
 int nx;
 int ny;
 unsigned char ** pixel;
 int data_type;
} image ;


typedef struct
{
image im;
} O_i;

#ifdef __cplusplus
extern "C" {
#endif

float   find_distance_from_center_2(float *x1, float *lfx1, int cl, int flag, int lfilter, int black_circle,
                                    float *corr, float *deriv, fft_plan *fp, filter_bt *fil, float *noise, int max_limit);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
int fill_avg_profiles_from_im_diff(O_i *oi, int xc, int yc, int cl, int cw, int *x, int *y, int *xd, int *yd,
                                   int *x_var, int *y_var, int *xd_var, int *yd_var);
#ifdef __cplusplus
}

#endif


#ifdef __cplusplus
extern "C" {
#endif
int fill_X_avg_profiles_from_im_diff(O_i *oi, int xc, int yc, int cl, int cw, int *x, int *y, int *xd, int *yd,
                                     int *x_var, int *y_var, int *xd_var, int *yd_var);
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif
int check_bead_not_lost_diff(int *x, int *xd, float *x1, int *y, int *yd, float *y1, int cl, float bd_mul, int do_diff);

#ifdef __cplusplus
}
#endif



int find_max1(float *x, int nx, float *Max_pos, float *Max_val, float *Max_deriv, int max_limit);


int correlate_1d_sig_and_invert(float *x1, int nx, int window, float *fxl1, int lfilter, int remove_dc, fft_plan *fp,
                                filter_bt *fil, float *noise);
