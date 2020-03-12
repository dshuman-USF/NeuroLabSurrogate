#define _GNU_SOURCE
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <float.h>
#include<X11/Xlib.h>
#include <X11/Xutil.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include "unpin.h"
#if !HAVE_ASPRINTF && HAVE_CONFIG_H
int asprintf(char **buffer, char *fmt, ...);
#endif

static int width = 1024;
static int height = 768;
static Display *dpy;
static Window rootwin;
static Window win;
static int scr;
static GC gc;
static Colormap cm;
static XColor yellow, blue, green, red;

static int stim[] = {
  110,183,247,296,363,439,561,671,758,879,937,1011,1268,1332,1445,1612,1953,2193,
  2270,2497,2593,2787,2921,3008,3305,3465,3604,3702,4026,4217,4399,4600,4804,5027,
  5176,5547,5608,5716,5844,6279,6442,6594,6655,6777,7056,7241,7354,7485,7548,7656,
  7877,8013,8064,8151,8402,8494,8524,8666,9104,9400,9808,10133,10238,10378,10561,
  10804,10963,11065,11112,11173,11206,11249,11316,11454,11654,11754,11870,11950,
  12029,12230,12278,12600,12637,12704,12943,13145,13210,13270,13494,13590,13708,
  13820,14017,14644,14755,14846,14943,15007,15266,15883
};
#define SPIKE_COUNT (sizeof stim / sizeof stim[0])
static int spike_count = SPIKE_COUNT;

static int time_spikes[SPIKE_COUNT];
static int min_isi;
static int *lodif;
static int t_est_count;
static int t_est_spikes[SPIKE_COUNT];
static double t_est_ispikes[SPIKE_COUNT];
static int t_est[SPIKE_COUNT];
static double tgen[2 * SPIKE_COUNT];
static double ygen[2 * SPIKE_COUNT];
static int gen_count;
static double isi[SPIKE_COUNT];
static double g, scale;
double max_d;
double min_d = DBL_MAX;
double max_d_loc;
double max_d_top;
double max_d_bot;

int xl0, xl1;
int x0;
int x1;
int x2, x3;
int x4, x5;
int xd01, xd23, xd45;

int y0_;
int y1_;
int y2;
int y3;
int yd23;
int y4;
int y5;
int yd45;
int y6, y7;

double td;

static int
compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

static gsl_error_handler_t *default_gsl_err_handler;
static int cdf_gamma_P_error_count;

static void
gsl_cdf_gamma_P_err (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_errno == GSL_EMAXITER) {
    if (cdf_gamma_P_error_count == 0)
      fprintf (stderr, "gsl_cdf_gamma_P: %s\n", gsl_strerror (gsl_errno));
    cdf_gamma_P_error_count++;
    return;
  }
  default_gsl_err_handler (reason, file, line, gsl_errno);
}
                          
static void
get_t_est (void)
{
  char drop[spike_count];
  int n, from, to;
  int drop_count = spike_count - t_est_count;
  static int last_t_est_count = -1;

  if (t_est_count == last_t_est_count)
    return;
  last_t_est_count = t_est_count;
  memset (drop, 0, sizeof drop);
  int dropped = 0;
  for (n = 0; dropped < drop_count; n++)
    if (lodif[n] != spike_count + 1) {
      drop[lodif[n]-2] = 1;
      dropped++;
    }
  from = to = 0;
  for (from = 0; from < spike_count; from++)
    if (!drop[from]) {
      t_est[to] = stim[from];
      t_est_spikes[to] = time_spikes[from];
      to++;
    }
  from = to = 0;
  int last_t_est = 0;
  int last_t_est_spikes = 0;
  for (from = 0; from < spike_count; from++) {
    while (t_est[to] < stim[from]) {
      last_t_est = t_est[to];
      last_t_est_spikes = t_est_spikes[to];
      to++;
    }
    t_est_ispikes[from] = (last_t_est_spikes + (double)(stim[from] - last_t_est) / (t_est[to] - last_t_est)
                           * (t_est_spikes[to] - last_t_est_spikes));
  }
  isi[0] = t_est_ispikes[0];
  double sum_log = log (isi[0]);
  double sum = isi[0];
  for (n = 1; n < spike_count; n++) {
    isi[n] = t_est_ispikes[n] - t_est_ispikes[n - 1];
    sum_log += log (isi[n]);
    sum += isi[n];
  }
  double s = log (sum / spike_count) - sum_log / spike_count;
  double t = s - 3;
  g = s == 0 ? 1e6 : (3 - s + sqrt (t*t + 24 * s)) / (12 * s);
  scale = 1 / g;
  
  qsort (isi, spike_count, sizeof (double), compare_doubles);
  
  max_d = 0;
  cdf_gamma_P_error_count = 0;
  default_gsl_err_handler = gsl_set_error_handler (gsl_cdf_gamma_P_err);
  for (n = 0; n < spike_count; n++) {
    double actual_lo = (double) n / spike_count;
    double actual_hi = (n + 1.) / spike_count;
    double cdf = gsl_cdf_gamma_P (isi[n], g, scale);
    if (actual_hi - cdf > max_d) {
      max_d = actual_hi - cdf;
      max_d_loc = isi[n];
      max_d_top = actual_hi;
      max_d_bot = cdf;
    }
    if (cdf - actual_lo > max_d) {
      max_d = cdf - actual_lo;
      max_d_loc = isi[n];
      max_d_top = cdf;
      max_d_bot = actual_lo;
    }
  }
  if (max_d < min_d)
    min_d = max_d;
  gsl_set_error_handler (default_gsl_err_handler);

  static gsl_rng *rng;
  if (rng == NULL)
    rng = gsl_rng_alloc (gsl_rng_mt19937); 
  double ysum = 0;
  gen_count = 0;
  while (1) {
    double ranval = gsl_ran_gamma (rng, g, scale);
    ysum += ranval;
    if (ysum > spike_count)
      break;
    if (gen_count == spike_count * 2) {
      fprintf (stderr, "bug at %s line %d\n", __FILE__, __LINE__);
      exit (1);
    }
    ygen[gen_count++] = ysum;
  }
  last_t_est = 0;
  last_t_est_spikes = 0;
  int ne = 0;
  int ng = 0;
  for (ng = 0; ng < gen_count; ng++) {
    while (t_est_spikes[ne] < ygen[ng]) {
      last_t_est = t_est[ne];
      last_t_est_spikes = t_est_spikes[ne];
      ne++;
    }
    tgen[ng] = (last_t_est + (ygen[ng] - last_t_est_spikes) / (t_est_spikes[ne] - last_t_est_spikes)
                * (t_est[ne] - last_t_est));
  }
}

static void
plot (int *stim, int *time_spikes, int local_spike_count)
{
  int n;
  int lastx = x0;
  int lastt = 0;
  int lasty = y3;
  int lastcy = y5;
  double lastcnt = 0;
  for (n = 0; n < local_spike_count; n++) {
    int x = x0 + stim[n] / td * xd01;
    int y = y3 - (time_spikes[n] - lastcnt) / (stim[n] - lastt) * min_isi * yd23;
    int cy = y5 - (double)time_spikes[n] / spike_count * yd45;
    XDrawLine (dpy, win, gc, lastx, y, x, y);
    XDrawLine (dpy, win, gc, lastx, lasty, lastx, y);
    XDrawLine (dpy, win, gc, lastx, lastcy, x, cy);
    lastt = stim[n];
    lasty = y;
    lastx = x;
    lastcy = cy;
    lastcnt = time_spikes[n];
  }
}

static void
plot_df (void)
{
  int n;
  int pxl[x5+1];
  memset (pxl, 0, sizeof pxl);

  static double max_isi = 2;

  if (isi[spike_count - 1] > max_isi)
    max_isi = isi[spike_count - 1];
  int count_max = 1;
  int x;
  for (n = 0; n < spike_count; n++) {
    x = x4 + isi[n] / max_isi * xd45;
    pxl[x]++;
    if (pxl[x] > count_max)
      count_max = pxl[x];
  }
  for (x = x4; x <= x5; x++) {
    if (pxl[x] == 0)
      continue;
    int y = y3 - (double)pxl[x] / count_max * yd23;
    XDrawLine (dpy, win, gc, x, y3, x, y);
  }
  int lastx = x4;
  int lasty = y5;
  for (n = 0; n < spike_count; n++) {
    x = x4 + isi[n] / max_isi * xd45;
    int y = y5 - (n + 1.) / spike_count * yd45;
    XDrawLine (dpy, win, gc, lastx, lasty, x, lasty);
    XDrawLine (dpy, win, gc, x, lasty, x, y);
    lastx = x;
    lasty = y;
  }

  lastx = x4;
  lasty = y5;
  double pdf[x5 + 1];
  double max_pdf = 0;
  for (x = x4 + 1; x <= x5; x++) {
    double isi = (double)(x - x4) / xd45 * max_isi;
    pdf[x] = gsl_ran_gamma_pdf (isi, g, scale);
    if (pdf[x] > max_pdf)
      max_pdf = pdf[x];
  }
  
  int lastyp = y3 - gsl_ran_gamma_pdf (0, g, scale) * yd23;
  cdf_gamma_P_error_count = 0;
  default_gsl_err_handler = gsl_set_error_handler (gsl_cdf_gamma_P_err);
  for (x = x4 + 1; x <= x5; x++) {
    double isi = (double)(x - x4) / xd45 * max_isi;
    double cdf = gsl_cdf_gamma_P (isi, g, scale);
    int y = y5 - cdf * yd45;
    int yp = y3 - pdf[x] / max_pdf * yd23;
    XDrawLine (dpy, win, gc, lastx, lasty, x, lasty);
    XDrawLine (dpy, win, gc, x, lasty, x, y);
    XDrawLine (dpy, win, gc, lastx, lastyp, x, yp);
    lastx = x;
    lasty = y;
    lastyp = yp;
  }
  gsl_set_error_handler (default_gsl_err_handler);
  x = x4 + max_d_loc / max_isi * xd45;
  int yt = y5 - max_d_top * yd45;
  int yb = y5 - max_d_bot * yd45;
  XSetForeground (dpy, gc, red.pixel);
  XDrawLine (dpy, win, gc, x, yt, x, yb);
  XSetForeground (dpy, gc, BlackPixel (dpy, scr));
  
}

static void
paint (void)
{
  XWindowAttributes xwa;

  XGetWindowAttributes (dpy, win, &xwa);
  width = xwa.width;
  height = xwa.height;
  XClearWindow(dpy, win);

  int n;
  xl0 = .01 * width;
  xl1 = .07 * width;
  x0 = .08 * width;
  x1 = .47 * width;
  x2 = .50 * width;
  x3 = .56 * width;
  x4 = .60 * width;
  x5 = .99 * width;

  xd01 = x1 - x0;
  xd23 = x3 - x2;
  xd45 = x5 - x4;

  y0_ = height * .01;
  y1_ = height * .08;
  y2 = height * .09;
  y3 = height * .47;
  yd23 = y3 - y2;
  y4 = height * .48;
  y5 = height * .9;
  yd45 = y5 - y4;
  y6 = height * .92;
  y7 = height * .99;

  td = stim[spike_count - 1];
  get_t_est ();
  for (n = 0; n < spike_count; n++) {
    int x = x0 + stim[n] / td * xd01;
    int cy = y5 - (double)t_est_ispikes[n] / spike_count * yd45;
    XDrawLine (dpy, win, gc, x, y0_, x, y1_);
    XSetForeground (dpy, gc, yellow.pixel);
    XDrawLine (dpy, win, gc, x, y1_, x, cy);
    XDrawLine (dpy, win, gc, x, cy, x2, cy);
    XSetForeground (dpy, gc, BlackPixel (dpy, scr));
    XDrawLine (dpy, win, gc, x2, cy, x3, cy);
  }
  for (n = 0; n < gen_count; n++) {
    int x = x0 + tgen[n] / td * xd01;
    int y = y5 - ygen[n] / spike_count * yd45;
    XDrawLine (dpy, win, gc, xl0, y, xl1, y);
    XSetForeground (dpy, gc, green.pixel);
    XDrawLine (dpy, win, gc, xl1, y, x, y);
    XDrawLine (dpy, win, gc, x, y, x, y6);
    XSetForeground (dpy, gc, BlackPixel (dpy, scr));
    XDrawLine (dpy, win, gc, x, y6, x, y7);
  }
  XSetForeground (dpy, gc, red.pixel);
  plot (t_est, t_est_spikes, t_est_count);
  XSetForeground (dpy, gc, BlackPixel (dpy, scr));
  plot (stim, time_spikes, spike_count);
  
  XDrawLine (dpy, win, gc, x0, y3, x1, y3);
  XDrawLine (dpy, win, gc, x0, y2, x0, y3);
  XDrawLine (dpy, win, gc, x0, y5, x1, y5);
  XDrawLine (dpy, win, gc, x0, y4, x0, y5);
  XDrawLine (dpy, win, gc, x4, y3, x5, y3);
  XDrawLine (dpy, win, gc, x4, y2, x4, y3);
  XDrawLine (dpy, win, gc, x4, y5, x5, y5);
  XDrawLine (dpy, win, gc, x4, y4, x4, y5);
  plot_df ();
  char *s;
  if (min_d == max_d)
    XSetForeground (dpy, gc, red.pixel);
  if (asprintf (&s, "    d: %g", max_d) == -1) exit (1);
  XDrawString (dpy, win, gc, (x4 + x5) / 2 + 5, (y4 + y5) / 2, s, strlen(s));
  free (s);
  if (asprintf (&s, "min d: %g", min_d) == -1) exit (1);
  XDrawString (dpy, win, gc, (x4 + x5) / 2 + 5, (y4 + y5) / 2 + 15, s, strlen(s));
  free (s);
  XSetForeground (dpy, gc, BlackPixel (dpy, scr));
  if (asprintf (&s, "shape: %g", g) == -1) exit (1);
  XDrawString (dpy, win, gc, (x4 + x5) / 2 + 5, (y4 + y5) / 2 + 30, s, strlen(s));
  free (s);
  if (asprintf (&s, "count: %d", t_est_count) == -1) exit (1);
  XDrawString (dpy, win, gc, (x4 + x5) / 2 + 5, (y4 + y5) / 2 + 45, s, strlen(s));
  free (s);

  s = " <- original spike train";
  XDrawString (dpy, win, gc, x1, (y0_ + y1_) / 2, s, strlen(s));
  s = " <- surrogate spike train";
  XDrawString (dpy, win, gc, x1, (y6 + y7) / 2, s, strlen(s));
  s = "  firing rate estimates (black is \"reciprocal ISI\")";
  XDrawString (dpy, win, gc, x0, y3 - 8, s, strlen(s));
  s = "  integral of firing rate estimates";
  XDrawString (dpy, win, gc, x0, y4 + .2 * (y5 - y4), s, strlen(s));
  s = "    ISI distribution and fitted gamma distribution PDF";
  XDrawString (dpy, win, gc, x4, y2 - 8, s, strlen(s));
  s = "empirical CDF and gamma CDF";
  XDrawString (dpy, win, gc, x4 + 8, y4 + 12, s, strlen(s));
  s = "red is max difference";
  XDrawString (dpy, win, gc, x4 + 8, y4 + 30, s, strlen(s));
  s = "| random";
  XDrawString (dpy, win, gc, xl0, y4-7, s, strlen(s));
  s = "v";
  XDrawString (dpy, win, gc, xl0, y4-2, s, strlen(s));
  s = "| ISI's";
  XDrawString (dpy, win, gc, x2, y4-7, s, strlen(s));
  s = "v";

  XDrawString (dpy, win, gc, x2, y4-2, s, strlen(s));
}

static void
init (void)
{
  int n;
  min_isi = stim[1] - stim[0];
  for (n = 2; n < spike_count; n++)
    if (stim[n] - stim[n - 1] < min_isi)
      min_isi = stim[n] - stim[n - 1];
  cm = DefaultColormap(dpy, scr);
  XAllocNamedColor (dpy, cm, "yellow", &yellow, &yellow);
  XAllocNamedColor (dpy, cm, "blue", &blue, &blue);
  XAllocNamedColor (dpy, cm, "green", &green, &green);
  XAllocNamedColor (dpy, cm, "red", &red, &red);
  lodif = unpin (stim, spike_count);
  for (n = 0; n < spike_count; n++)
    time_spikes[n] = n + 1;
}

int
main (int argc, char **argv)
{
  t_est_count = spike_count;
  XEvent e;

  if (argc > 1)
    spike_count = atoi (argv[1]);

  if (!(dpy = XOpenDisplay (NULL))) {
    fprintf (stderr, "ERROR: Could not open display\n");
    exit (1);
  }
  init ();

  scr = DefaultScreen (dpy);
  rootwin = RootWindow (dpy, scr);

  win = XCreateSimpleWindow (dpy, rootwin, 1, 1, width, height, 0, 
                          BlackPixel (dpy, scr), WhitePixel (dpy, scr));
  gc = XCreateGC (dpy, win, 0, NULL);

  XSelectInput (dpy, win, ExposureMask|KeyPressMask);
  XMapWindow (dpy, win);

  while (1) {
    XNextEvent (dpy, &e);
    if (e.type == Expose && e.xexpose.count < 1) {
      paint ();
    }
    else if (e.type == KeyPress) {
#     define BUFSZ 8
      static char buf[BUFSZ];
      XLookupString (&e.xkey, buf, BUFSZ, 0, 0);
      if (buf[0] == 'q')
        break;
      if (buf[0] == '-' && t_est_count > 1)
        t_est_count--;
      if ((buf[0] == '+' || buf[0] == '=') && t_est_count < spike_count)
        t_est_count++;
      paint ();
    }
  }

  XCloseDisplay (dpy);
  return 0;
}
