#include "Bead.hpp"
#include <cstdio>


int mouse_draged;

Bead::Bead(int clh,int cwh, int xch, int ych)
{
  ncl = cl = clh;
  ncw =cw = cwh;

  black_circle = 0;
  xy_trp = fftbt_init(NULL, cl);
  xy_fil = filter_init_bt(NULL, cl);

  type = 0;

  dis_filter = cl >> 2;
  max_limit = cl >> 1;

  bd_mul = 1;
  do_diff_track = 1;

  fence_x = -1;
  fence_y = -1;

  if (xy_trp == NULL || xy_fil == NULL )
  {
      return; // MR Comebck
  }

  not_lost = -1;
  ci = 0;

  saved_x = last_x_not_lost = (float)xc;
  saved_y = last_y_not_lost = (float)yc;
  xc = x0 = xch;
  yc = y0 = ych;
  mouse_draged = 0;
  not_lost = NB_FRAMES_BEF_LOST;
  start_im = 0;
  cross_45 = 0;

}




void Bead::Track(O_i *oi)
{
cl = ncl;              // we update cross changes
cw = ncw;
int cl2 = cl / 2;
float deriv,corr;
int i;
float xcf, ycf; // new position of the bead



if (ci & 0x01)
{
    x_var = x_v[1];
    y_var = y_v[1];
    xs_var = xs_v[1];
    ys_var = ys_v[1];

}
else
{
    x_var = x_v[0];
    y_var = y_v[0];
    xs_var = xs_v[0];
    ys_var = ys_v[0];

}

if (cross_45)
        {
            fill_X_avg_profiles_from_im_diff(oi, xc, yc, cl, cw, xi, yi, xis, yis, x_var, y_var, xs_var, ys_var);
        }
        else
        {
            fill_avg_profiles_from_im_diff(oi, xc, yc, cl, cw, xi, yi, xis, yis, x_var, y_var, xs_var, ys_var);
        }


int lost_this_time = check_bead_not_lost_diff(xi, xis, xf, yi, yis, yf, cl, bd_mul, do_diff_track);

dx = find_distance_from_center_2(xf, fx1, cl, 0, dis_filter, !black_circle, &corr, &deriv,
                                 xy_trp, xy_fil, xnoise + ci, max_limit);
xnoise[ci] *= (do_diff_track) ? sqrt((x_v_var_m + xs_v_var_m) / 2) : sqrt(x_v_var_m / 2);

for (i = 0; i < cl;  i++)
{
    xcor[i] = fx1[i];
}

dxcor_dx = dcor_dx[ci] = deriv;
xmax = corr;
dy = find_distance_from_center_2(yf, fy1, cl, 0, dis_filter, !black_circle, &corr, &deriv,
                                 xy_trp, xy_fil, ynoise + ci, max_limit);
ynoise[ci] *= (do_diff_track) ? sqrt((y_v_var_m + ys_v_var_m) / 2) : sqrt(y_v_var_m / 2);

for (i = 0; i < cl;  i++)
{
    ycor[i] = fy1[i];
}

dycor_dy = dcor_dy[ci] = deriv;
ymax = corr;

if (cross_45)
{

        x[ci] = xcf = xc + (dx + dy);
        y[ci] = ycf = yc + (dx - dy) - 1;

    dx = (dx + dy);
    dy = (dx - dy) - 1;
}
else
{

        x[ci] = xcf = dx + xc;               // we save its new x, y position
        y[ci] = ycf = yc + dy;

    dx = dx;
    dy = dy;
}
//printf("x = %f\n",x[ci]);
//printf("y = %f\n",y[ci]);

//xc = (int) (xcf + 0.5); Should be update xc ?
//yc = (int) (xcf + 0.5);


ci++;
if (ci == TRACKING_BUFFER_SIZE) ci =0;

}

std::vector<float> Bead::get_copy_x_array()
{
  std::vector<float> a(x, x+TRACKING_BUFFER_SIZE);
  return a;
}
std::vector<float> Bead::get_copy_y_array()
{
  std::vector<float> a(y, y+TRACKING_BUFFER_SIZE);
  return a;
}
