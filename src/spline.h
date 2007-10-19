///////////////////////////////////////////////////////////////////////////////
//
// spline.h
//
///////////////////////////////////////////////////////////////////////////////
// $Id: spline.h,v 1.3 2007-10-19 16:24:06 fgygi Exp $
void spline(int n, double *x, double *y, double yp_left, double yp_right,
            int bcnat_left, int bcnat_right, double *y2);
void splint (int n, double *xa, double *ya, double *y2a, double x, double *y);
void splintd (int n, double *xa, double *ya, double *y2a,
              double x, double *y, double *dy);
