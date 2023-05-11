#pragma once
#include "track_util.h"
#include <vector>
#define TRACKING_BUFFER_SIZE 128
#define NB_FRAMES_BEF_LOST 128




class Bead
{
public:
  Bead(int a,int b, int c, int d);
  ~Bead();
  void Track(O_i *oi);
  float get_x(){return x[ci-1];};
  float get_y(){return y[ci-1];};
  std::vector<float> get_copy_x_array();
  std::vector<float> get_copy_y_array();




private:
  int xc;
  int yc;
  int cl;
  int cw;
  int x0;
  int y0;
  int cross_45;
  int ci;

  float x[TRACKING_BUFFER_SIZE];     // bead position
  float y[TRACKING_BUFFER_SIZE];

  int *x_var, *y_var;
  int *xs_var, *ys_var;

  int ncl, ncw;

  int do_diff_track;
  int dis_filter;
  int max_limit;
  int bd_mul;
  fft_plan *xy_trp;                // the fft plan for x, y tracking
  filter_bt *xy_fil;                  // the filter paln for x, y tracking
  float xnoise[TRACKING_BUFFER_SIZE];
  float ynoise[TRACKING_BUFFER_SIZE];
  float x_v_var, xs_v_var, x_v_var_m, xs_v_var_m;
  float y_v_var, ys_v_var, y_v_var_m, ys_v_var_m;
  float dxcor_dx, dxcor_dx_m, dycor_dy, dycor_dy_m;
  int x_v[2][256], y_v[2][256];     // tmp float buffers to variance in avg profile
  int xs_v[2][256], ys_v[2][256];   // tmp float buffers to v;ariance in avg side profile
  float saved_x,saved_y,last_x_not_lost,last_y_not_lost;
  int not_lost;
  int start_im;
  int type;
  int black_circle;
  int fence_x;
  int fence_y;

  float xf[256], yf[256];           // tmp float buffers to do fft
  float xcor[256], ycor[256];
  float fx1[512], fy1[512];


  int xi[256], yi[256];             // tmp int buffers to avg profile
  int xis[256], yis[256];           // tmp int buffers to avg side profile
        // tmp float buffers to do fft
  float dx;                         // the x shift in profile
  float dy;

  float dcor_dx[TRACKING_BUFFER_SIZE];
  float dcor_dy[TRACKING_BUFFER_SIZE];     // derivative of the convolution
  float xmax, ymax;






};
