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

#ifndef variations_included
#define variations_included

#include "private.h"



/* Variation functions */
void var0_linear(flam3_iter_helper *, double);
void var1_sinusoidal(flam3_iter_helper *, double);
void var2_spherical(flam3_iter_helper *, double);
void var3_swirl(flam3_iter_helper *, double);
void var4_horseshoe(flam3_iter_helper *, double);
void var5_polar(flam3_iter_helper *, double);
void var6_handkerchief(flam3_iter_helper *, double);
void var7_heart(flam3_iter_helper *, double);
void var8_disc(flam3_iter_helper *, double);
void var9_spiral(flam3_iter_helper *, double);
void var10_hyperbolic(flam3_iter_helper *, double);
void var11_diamond(flam3_iter_helper *, double);
void var12_ex(flam3_iter_helper *, double);
void var13_julia(flam3_iter_helper *, double);
void var14_bent(flam3_iter_helper *, double);
void var15_waves(flam3_iter_helper *, double);
void var16_fisheye(flam3_iter_helper *, double);
void var17_popcorn(flam3_iter_helper *, double);
void var18_exponential(flam3_iter_helper *, double);
void var19_power(flam3_iter_helper *, double);
void var20_cosine(flam3_iter_helper *, double);
void var21_rings(flam3_iter_helper *, double);
void var22_fan(flam3_iter_helper *, double);
void var23_blob(flam3_iter_helper *, double);
void var24_pdj(flam3_iter_helper *, double);
void var25_fan2(flam3_iter_helper *, double);
void var26_rings2(flam3_iter_helper *, double);
void var27_eyefish(flam3_iter_helper *, double);
void var28_bubble(flam3_iter_helper *, double);
void var29_cylinder(flam3_iter_helper *, double);
void var30_perspective(flam3_iter_helper *, double);
void var31_noise(flam3_iter_helper *, double);
void var32_juliaN_generic(flam3_iter_helper *, double);
void var33_juliaScope_generic(flam3_iter_helper *, double);
void var34_blur(flam3_iter_helper *, double);
void var35_gaussian(flam3_iter_helper *, double);
void var36_radial_blur(flam3_iter_helper *, double);
void var37_pie(flam3_iter_helper *, double);
void var38_ngon(flam3_iter_helper *, double);
void var39_curl(flam3_iter_helper *, double);
void var40_rectangles(flam3_iter_helper *, double);
void var41_arch(flam3_iter_helper *, double);
void var42_tangent(flam3_iter_helper *, double);
void var43_square(flam3_iter_helper *, double);
void var44_rays(flam3_iter_helper *, double);
void var45_blade(flam3_iter_helper *, double);
void var46_secant2(flam3_iter_helper *, double);
void var47_twintrian(flam3_iter_helper *, double);
void var48_cross(flam3_iter_helper *, double);
void var49_disc2(flam3_iter_helper *, double);
void var50_supershape(flam3_iter_helper *, double);
void var51_flower(flam3_iter_helper *, double);
void var52_conic(flam3_iter_helper *, double);
void var53_parabola(flam3_iter_helper *, double);
void var54_bent2(flam3_iter_helper *, double);
void var55_bipolar(flam3_iter_helper *, double);
void var56_boarders(flam3_iter_helper *, double);
void var57_butterfly(flam3_iter_helper *, double);
void var58_cell(flam3_iter_helper *, double);
void var59_cpow(flam3_iter_helper *, double);
void var60_curve(flam3_iter_helper *, double);
void var61_edisc(flam3_iter_helper *, double);
void var62_elliptic(flam3_iter_helper *, double);
void var63_escher(flam3_iter_helper *, double);
void var64_foci(flam3_iter_helper *, double);
void var65_lazysusan(flam3_iter_helper *, double);
void var66_loonie(flam3_iter_helper *, double);
void var67_pre_blur(flam3_iter_helper *, double);
void var68_modulus(flam3_iter_helper *, double);
void var69_oscope(flam3_iter_helper *, double);
void var70_polar2(flam3_iter_helper *, double);
void var71_popcorn2(flam3_iter_helper *, double);
void var72_scry(flam3_iter_helper *, double);
void var73_separation(flam3_iter_helper *, double);
void var74_split(flam3_iter_helper *, double);
void var75_splits(flam3_iter_helper *, double);
void var76_stripes(flam3_iter_helper *, double);
void var77_wedge(flam3_iter_helper *, double);
void var78_wedge_julia(flam3_iter_helper *, double);
void var79_wedge_sph(flam3_iter_helper *, double);
void var80_whorl(flam3_iter_helper *, double);
void var81_waves2(flam3_iter_helper *, double);
void var82_exp (flam3_iter_helper *f, double weight);
void var83_log (flam3_iter_helper *f, double weight);
void var84_sin (flam3_iter_helper *f, double weight);
void var85_cos (flam3_iter_helper *f, double weight);
void var86_tan (flam3_iter_helper *f, double weight);
void var87_sec (flam3_iter_helper *f, double weight);
void var88_csc (flam3_iter_helper *f, double weight);
void var89_cot (flam3_iter_helper *f, double weight);
void var90_sinh (flam3_iter_helper *f, double weight);
void var91_cosh (flam3_iter_helper *f, double weight);
void var92_tanh (flam3_iter_helper *f, double weight);
void var93_sech (flam3_iter_helper *f, double weight);
void var94_csch (flam3_iter_helper *f, double weight);
void var95_coth (flam3_iter_helper *f, double weight);
void var96_auger (flam3_iter_helper *f, double weight);
void var97_flux (flam3_iter_helper *f, double weight);
void var98_mobius (flam3_iter_helper *f, double weight);

/* Precalculation functions */
void perspective_precalc(flam3_xform *xf);
void juliaN_precalc(flam3_xform *xf);
void juliaScope_precalc(flam3_xform *xf);
void radial_blur_precalc(flam3_xform *xf);
void waves_precalc(flam3_xform *xf);
void disc2_precalc(flam3_xform *xf);
void supershape_precalc(flam3_xform *xf);
void wedgeJulia_precalc(flam3_xform *xf);

void xform_precalc(flam3_genome *cp, int xi);
int prepare_precalc_flags(flam3_genome *);

int apply_xform(flam3_genome *cp, int fn, double *p, double *q, randctx *rc);
void initialize_xforms(flam3_genome *thiscp, int start_here);
#endif
