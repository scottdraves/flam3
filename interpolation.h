/*
    FLAM3 - cosmic recursive fractal flames
    Copyright (C) 1992-2009 Spotworks LLC

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef interpolation_included
#define interpolation_included

#include "private.h"

#define INTERP(x)  do { result->x = 0.0; \
   for (k = 0; k < ncp; k++) result->x += c[k] * cpi[k].x; } while(0)

#define INTERI(x)  do { double tt = 0.0; \
   for (k = 0; k < ncp; k++) tt += c[k] * cpi[k].x; \
   result->x = (int)rint(tt); } while(0)

double adjust_percentage(double in);
double motion_funcs(int funcnum, double timeval);

double smoother(double t);
double get_stagger_coef(double t, double stagger_prc, int num_xforms, int this_xform);

double det_matrix(double s[2][2]);
int id_matrix(double s[3][2]);
int zero_matrix(double s[3][2]);
void copy_matrix(double to[3][2], double from[3][2]);
void clear_matrix(double m[3][2]);
void sum_matrix(double s, double m1[3][2], double m2[3][2]);
void mult_matrix(double s1[2][2], double s2[2][2], double d[2][2]);

int compare_xforms(const void *av, const void *bv);

void interpolate_cmap(flam3_palette cmap, double blend,
              int index0, double hue0, int index1, double hue1);
void interp_and_convert_back(double *c, int ncps, int xfi, double cxang[4][2], 
                             double cxmag[4][2], double cxtrn[4][2],double store_array[3][2]);                             
void convert_linear_to_polar(flam3_genome *cp, int ncps, int xfi, int cflag, 
                             double cxang[4][2], double cxmag[4][2], double cxtrn[4][2]);

void interpolate_catmull_rom(flam3_genome cps[], double t, flam3_genome *result);
void flam3_interpolate_n(flam3_genome *result, int ncp, flam3_genome *cpi, double *c, double stagger);
void establish_asymmetric_refangles(flam3_genome *cp, int ncps);
void flam3_align(flam3_genome *dst, flam3_genome *src, int nsrc);
#endif
