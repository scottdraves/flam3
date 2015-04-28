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

#include "interpolation.h"
#include "palettes.h"

double adjust_percentage(double in) {

   if (in==0.0)
      return(0.0);
   else
      return(pow(10.0, -log(1.0/in)/log(2)));

}

double motion_funcs(int funcnum, double timeval) {

   /* motion funcs should be cyclic, and equal to 0 at integral time values */
   /* abs peak values should be not be greater than 1                       */
   if (funcnum==MOTION_SIN) {
      return (sin(2.0*M_PI*timeval));
   } else if (funcnum==MOTION_TRIANGLE) {
      double fr = fmod(timeval,1.0);
      
      if (fr<0) fr+= 1.0;
      
      if (fr<=.25)
         fr = 4.0 * fr;
      else if (fr<=.75)
         fr = -4.0 * fr + 2.0;
      else
         fr = 4.0 * fr - 4.0;
     
      return(fr);
   } else { //if (funcnum==MOTION_HILL) {
      return( (1.0-cos(2.0*M_PI*timeval)) * 0.5);
   }
   
}

double det_matrix(double s[2][2]) {
   return s[0][0] * s[1][1] - s[0][1] * s[1][0];
}

int id_matrix(double s[3][2]) {
   return
      (s[0][0] == 1.0) &&
      (s[0][1] == 0.0) &&
      (s[1][0] == 0.0) &&
      (s[1][1] == 1.0) &&
      (s[2][0] == 0.0) &&
      (s[2][1] == 0.0);
}

int zero_matrix(double s[3][2]) {
   return
      (s[0][0] == 0.0) &&
      (s[0][1] == 0.0) &&
      (s[1][0] == 0.0) &&
      (s[1][1] == 0.0) &&
      (s[2][0] == 0.0) &&
      (s[2][1] == 0.0);
}

void copy_matrix(double to[3][2], double from[3][2]) {

   to[0][0] = from[0][0];
   to[0][1] = from[0][1];
   to[1][0] = from[1][0];
   to[1][1] = from[1][1];
   to[2][0] = from[2][0];
   to[2][1] = from[2][1];
}


void clear_matrix(double m[3][2]) {
   m[0][0] = 0.0;
   m[0][1] = 0.0;
   m[1][0] = 0.0;
   m[1][1] = 0.0;
   m[2][0] = 0.0;
   m[2][1] = 0.0;
}

void sum_matrix(double s, double m1[3][2], double m2[3][2]) {

   m2[0][0] += s * m1[0][0];
   m2[0][1] += s * m1[0][1];
   m2[1][0] += s * m1[1][0];
   m2[1][1] += s * m1[1][1];
   m2[2][0] += s * m1[2][0];
   m2[2][1] += s * m1[2][1];
}

void mult_matrix(double s1[2][2], double s2[2][2], double d[2][2]) {
   d[0][0] = s1[0][0] * s2[0][0] + s1[1][0] * s2[0][1];
   d[1][0] = s1[0][0] * s2[1][0] + s1[1][0] * s2[1][1];
   d[0][1] = s1[0][1] * s2[0][0] + s1[1][1] * s2[0][1];
   d[1][1] = s1[0][1] * s2[1][0] + s1[1][1] * s2[1][1];
}

int compare_xforms(const void *av, const void *bv) {
   flam3_xform *a = (flam3_xform *) av;
   flam3_xform *b = (flam3_xform *) bv;
   double aa[2][2];
   double bb[2][2];
   double ad, bd;

   aa[0][0] = a->c[0][0];
   aa[0][1] = a->c[0][1];
   aa[1][0] = a->c[1][0];
   aa[1][1] = a->c[1][1];
   bb[0][0] = b->c[0][0];
   bb[0][1] = b->c[0][1];
   bb[1][0] = b->c[1][0];
   bb[1][1] = b->c[1][1];
   ad = det_matrix(aa);
   bd = det_matrix(bb);

   if (a->color_speed > b->color_speed) return 1;
   if (a->color_speed < b->color_speed) return -1;
   if (a->color_speed) {
      if (ad < 0) return -1;
      if (bd < 0) return 1;
      ad = atan2(a->c[0][0], a->c[0][1]);
      bd = atan2(b->c[0][0], b->c[0][1]);
   }

   if (ad < bd) return -1;
   if (ad > bd) return 1;
   return 0;
}

void interpolate_cmap(flam3_palette cmap, double blend,
                      int index0, double hue0, int index1, double hue1) {
                 
   flam3_palette p0,p1;
   int i, j, rcode;

   rcode = flam3_get_palette(index0, p0, hue0);
   if (rcode<0)
      fprintf(stderr,"unable to retrieve palette %d, setting to white\n", index0);
   rcode = flam3_get_palette(index1, p1, hue1);
   if (rcode<0)
      fprintf(stderr,"unable to retrieve palette %d, setting to white\n", index1);

   for (i = 0; i < 256; i++) {
      double t[5], s[5];
    
      rgb2hsv(p0[i].color, s);
      rgb2hsv(p1[i].color, t);
      
      s[3] = p0[i].color[3];
      t[3] = p1[i].color[3];
      
      s[4] = p0[i].index;
      t[4] = p1[i].index;

      /* take the shorter way around the hue circle */
//      if (0 == i) {
//	fprintf(stderr, "xxx interpolating between hues, %g %g\n", s[0], t[0]);
 //     }
      
      /* Correct the first hue to go the short way around */
      if ((s[0] - t[0]) > 3.0)  /* first hue much bigger than second hue */
        s[0] -= 6.0;
      if ((s[0] - t[0]) < -3.0)  /* first hue much smaller than second hue */
	s[0] += 6.0;
    
      for (j = 0; j < 5; j++)
         t[j] = ((1.0-blend) * s[j]) + (blend * t[j]);
         
      hsv2rgb(t, cmap[i].color);
      cmap[i].color[3] = t[3];
      cmap[i].index = t[4];      
   }
}

void interp_and_convert_back(double *c, int ncps, int xfi, double cxang[4][2], 
                             double cxmag[4][2], double cxtrn[4][2],double store_array[3][2]) {

   int i,col;
   
   double accang[2],accmag[2];
   double expmag;
   int accmode[2];
   
   accang[0] = 0.0;
   accang[1] = 0.0;
   accmag[0] = 0.0;
   accmag[1] = 0.0;

   accmode[0]=accmode[1]=0;
   
   /* accumulation mode defaults to logarithmic, but in special */
   /* cases we want to switch to linear accumulation            */
   for (col=0; col<2; col++) {
      for (i=0; i<ncps; i++) {
         if (log(cxmag[i][col])<-10)
            accmode[col]=1; // Mode set to linear interp
      }
   }
   
   for (i=0; i<ncps; i++) {
      for (col=0; col<2; col++) {
      
         accang[col] += c[i] * cxang[i][col];
         
         if (accmode[col]==0)
            accmag[col] += c[i] * log(cxmag[i][col]);
         else 
            accmag[col] += c[i] * (cxmag[i][col]);
            
         /* translation is ready to go */
         store_array[2][col] += c[i] * cxtrn[i][col];
      }
   }
   
   /* Convert the angle back to rectangular */
   for (col=0;col<2;col++) {
      if (accmode[col]==0)
         expmag = exp(accmag[col]);
      else
         expmag = accmag[col];
      
      store_array[col][0] = expmag * cos(accang[col]);
      store_array[col][1] = expmag * sin(accang[col]);
   }
   
}

void convert_linear_to_polar(flam3_genome *cp, int ncps, int xfi, int cflag, 
                             double cxang[4][2], double cxmag[4][2], double cxtrn[4][2]) {

   double c1[2],d,t,refang;
   int col,k;
   int zlm[2];

   for (k=0; k<ncps;k++) {

      /* Establish the angles and magnitudes for each component */
      /* Keep translation linear */
      zlm[0]=zlm[1]=0;
      for (col=0;col<2;col++) {
      
         if (cflag==0) {
            c1[0] = cp[k].xform[xfi].c[col][0];
            c1[1] = cp[k].xform[xfi].c[col][1];
            t = cp[k].xform[xfi].c[2][col];            
         } else {
            c1[0] = cp[k].xform[xfi].post[col][0];
            c1[1] = cp[k].xform[xfi].post[col][1];
            t = cp[k].xform[xfi].post[2][col];
         }
         
         cxang[k][col] = atan2(c1[1],c1[0]);
         cxmag[k][col] = sqrt(c1[0]*c1[0] + c1[1]*c1[1]);
         
         if (cxmag[k][col]== 0.0)
            zlm[col]=1;
         
         cxtrn[k][col] = t;
      }
      
      if (zlm[0]==1 && zlm[1]==0)
         cxang[k][0] = cxang[k][1];
      else if (zlm[0]==0 && zlm[1]==1)
         cxang[k][1] = cxang[k][0];
      
   }
   
   /* Make sure the rotation is the shorter direction around the circle */
   /* by adjusting each angle in succession, and rotate clockwise if 180 degrees */
   for (col=0; col<2; col++) {
      for (k=1;k<ncps;k++) {

         /* Adjust angles differently if we have an asymmetric case */   
         if (cp[k].xform[xfi].wind[col]>0 && cflag==0) {

            /* Adjust the angles to make sure that it's within wind:wind+2pi */
            refang = cp[k].xform[xfi].wind[col] - 2*M_PI;

            /* Make sure both angles are within [refang refang+2*pi] */
            while(cxang[k-1][col] < refang)
                  cxang[k-1][col] += 2*M_PI;
            
            while(cxang[k-1][col] > refang + 2*M_PI)
                  cxang[k-1][col] -= 2*M_PI;
                  
            while(cxang[k][col] < refang)
                  cxang[k][col] += 2*M_PI;
            
            while(cxang[k][col] > refang + 2*M_PI)
                  cxang[k][col] -= 2*M_PI;

         } else {

            /* Normal way of adjusting angles */
            d = cxang[k][col]-cxang[k-1][col];
      
            /* Adjust to avoid the -pi/pi discontinuity */
            if (d > M_PI+EPS)
               cxang[k][col] -= 2*M_PI;
            else if (d < -(M_PI-EPS) ) /* Forces clockwise rotation at 180 */
               cxang[k][col] += 2*M_PI;
         }
      }
   }
}

void interpolate_catmull_rom(flam3_genome cps[], double t, flam3_genome *result) {
   double t2 = t * t;
   double t3 = t2 * t;
   double cmc[4];

   cmc[0] = (2*t2 - t - t3) / 2;
   cmc[1] = (3*t3 - 5*t2 + 2) / 2;
   cmc[2] = (4*t2 - 3*t3 + t) / 2;
   cmc[3] = (t3 - t2) / 2;

   flam3_interpolate_n(result, 4, cps, cmc, 0);
}

double smoother(double t) {
  return 3*t*t - 2*t*t*t;
}

double get_stagger_coef(double t, double stagger_prc, int num_xforms, int this_xform) {

   /* max_stag is the spacing between xform start times if stagger_prc = 1.0 */
   double max_stag = (double)(num_xforms-1)/num_xforms;
   
   /* scale the spacing by stagger_prc */
   double stag_scaled = stagger_prc * max_stag;

   /* t ranges from 1 to 0 (the contribution of cp[0] to the blend) */
   /* the first line below makes the first xform interpolate first */
   /* the second line makes the last xform interpolate first */
   double st = stag_scaled * (num_xforms - 1 - this_xform) / (num_xforms-1);
//   double st = stag_scaled * (this_xform) / (num_xforms-1);
   double et = st + (1-stag_scaled);
   
//   printf("t=%f xf:%d st=%f et=%f : : %f\n",t,this_xform,st,et,smoother((t-st)/(1-stag_scaled)));
   
   if (t <= st)
      return (0);
   else if (t >= et)
      return (1);
   else
      return ( smoother((t-st)/(1-stag_scaled)) );

}
   


/* all cpi and result must be aligned (have the same number of xforms,
   and have final xform in the same slot) */
void flam3_interpolate_n(flam3_genome *result, int ncp,
          flam3_genome *cpi, double *c, double stagger) {
   int i, j, k, l, numstd;

   if (flam3_palette_interpolation_sweep != cpi[0].palette_interpolation) {

      /* rgb, hsv or hsv_circular modes. */
      double rgb_fraction = 0.0;  /* Assume that we are in plain hsv mode */
      if (flam3_palette_interpolation_rgb == cpi[0].palette_interpolation)
         rgb_fraction = 1.0;  /* All RGB output */
      else if (flam3_palette_interpolation_hsv_circular == cpi[0].palette_interpolation)
         rgb_fraction = cpi[0].hsv_rgb_palette_blend;

      for (i = 0; i < 256; i++) {
         double col_rgb[3], col_hsv[3];
         double new_rgb[3] = {0, 0, 0};
         double new_hsv[3] = {0, 0, 0};
         double new_count = 0, new_index = 0;
         int alpha1 = 1;
         
         /* Loop over each control point's color at this index */
         for (k = 0; k < ncp; k++) {

            /* Convert to hsv */
            rgb2hsv(cpi[k].palette[i].color, col_hsv);
	    
	    /* Store the rgb */
            for (l = 0; l < 3; l++)
               col_rgb[l] = cpi[k].palette[i].color[l];

	    if (2 == ncp && k == 0 && cpi[0].palette_interpolation == flam3_palette_interpolation_hsv_circular) {
               /* only adjust the first coordinate based on the other control point's hue */
               double second_color[3];
               rgb2hsv(cpi[1].palette[i].color, second_color);

               /* Adjust the hue so that we go the shorter direction around the circle */
               if ((second_color[0] - col_hsv[0]) > 3.0) {
                  col_hsv[0] += 6.0;
               } else if ((second_color[0] - col_hsv[0]) < -3.0) {
                  col_hsv[0] -= 6.0;	          
               }
            }

            for (j = 0; j < 3; j++) {
               new_rgb[j] += c[k] * col_rgb[j];
               new_hsv[j] += c[k] * col_hsv[j];
            }
               
            /* Compute the other two components of the color (count and index) */            
            new_count += c[k] * cpi[k].palette[i].color[3];
            if (cpi[k].palette[i].color[3] != 1.0)
               alpha1 = 0;
            new_index += c[k] * cpi[k].palette[i].index;
            
         }

         if (alpha1 == 1)
            new_count = 1.0;

         /* Convert the new hsv coord to back rgb */
         double new_hsv_rgb[3];
         hsv2rgb(new_hsv, new_hsv_rgb);

         /* Store the interpolated color in the new palette */
         for (l = 0; l < 3; l++)
            result->palette[i].color[l] = rgb_fraction * new_rgb[l] + (1.0-rgb_fraction) * new_hsv_rgb[l];
            
         result->palette[i].color[3] = new_count;
         result->palette[i].index = new_index;

         /* Clip the new color appropriately */
         for (j = 0; j < 4; j++) {
            if (result->palette[i].color[j] < 0.0)
               result->palette[i].color[j] = 0.0;
            if (result->palette[i].color[j] > 1.0)
               result->palette[i].color[j] = 1.0;
         }
         
         if (result->palette[i].index < 0.0)
            result->palette[i].index = 0.0;
         if (result->palette[i].index > 255.0)
            result->palette[i].index = 255.0;
      }
   } else {
      /* Sweep - not the best option for float indices */
      for (i = 0; i < 256; i++) {
         j = (i < (256 * c[0])) ? 0 : 1;
         result->palette[i] = cpi[j].palette[i];
      }
   }

   result->palette_index = flam3_palette_random;
   result->symmetry = 0;
   result->spatial_filter_select = cpi[0].spatial_filter_select;
   result->temporal_filter_type = cpi[0].temporal_filter_type;
   result->palette_mode = cpi[0].palette_mode;

   result->interpolation_type = cpi[0].interpolation_type;
   result->palette_interpolation = cpi[0].palette_interpolation;
   result->hsv_rgb_palette_blend = cpi[0].hsv_rgb_palette_blend;
   INTERP(brightness);
   INTERP(contrast);
   INTERP(highlight_power);
   INTERP(gamma);
   INTERP(vibrancy);
   INTERP(hue_rotation);
   INTERI(width);
   INTERI(height);
   INTERI(spatial_oversample);
   INTERP(center[0]);
   INTERP(center[1]);
   INTERP(rot_center[0]);
   INTERP(rot_center[1]);
   INTERP(background[0]);
   INTERP(background[1]);
   INTERP(background[2]);
   INTERP(pixels_per_unit);
   INTERP(spatial_filter_radius);
   INTERP(temporal_filter_exp);
   INTERP(temporal_filter_width);
   INTERP(sample_density);
   INTERP(zoom);
   INTERP(rotate);
   INTERI(nbatches);
   INTERI(ntemporal_samples);
   INTERP(estimator);
   INTERP(estimator_minimum);
   INTERP(estimator_curve);
   INTERP(gam_lin_thresh);
   
   /* Interpolate the chaos array */
   numstd = cpi[0].num_xforms - (cpi[0].final_xform_index >= 0);
   for (i=0;i<numstd;i++) {
      for (j=0;j<numstd;j++) {
         INTERP(chaos[i][j]);
         if (result->chaos[i][j]<0) result->chaos[i][j]=0;
         //chaos can be > 1
         //if (result->chaos[i][j]>1) result->chaos[i][j]=1.0;
      }
   }

   /* Interpolate each xform */
   for (i = 0; i < cpi[0].num_xforms; i++) {
   
      double csave[2] = {0, 0};     
      double td;
      int all_id;
      int nx = cpi[0].num_xforms-(cpi[0].final_xform_index>=0);
      
      if (ncp==2 && stagger>0 && i!=cpi[0].final_xform_index) {
         csave[0] = c[0];
         csave[1] = c[1];
         c[0] = get_stagger_coef(csave[0],stagger,nx,i);
         c[1] = 1.0-c[0];
      }
      
      
      INTERP(xform[i].density);
      td = result->xform[i].density;
      result->xform[i].density = (td < 0.0) ? 0.0 : td;
      INTERP(xform[i].color);
      if (result->xform[i].color<0) result->xform[i].color=0;
      if (result->xform[i].color>1) result->xform[i].color=1;

      INTERP(xform[i].color_speed);
      if (result->xform[i].color_speed<0) result->xform[i].color_speed=0;
      if (result->xform[i].color_speed>1) result->xform[i].color_speed=1;
      
      INTERP(xform[i].opacity);      
      INTERP(xform[i].animate);
      INTERP(xform[i].blob_low);
      INTERP(xform[i].blob_high);
      INTERP(xform[i].blob_waves);
      INTERP(xform[i].pdj_a);
      INTERP(xform[i].pdj_b);
      INTERP(xform[i].pdj_c);
      INTERP(xform[i].pdj_d);
      INTERP(xform[i].fan2_x);
      INTERP(xform[i].fan2_y);
      INTERP(xform[i].rings2_val);
      INTERP(xform[i].perspective_angle);
      INTERP(xform[i].perspective_dist);
      INTERP(xform[i].julian_power);
      INTERP(xform[i].julian_dist);
      INTERP(xform[i].juliascope_power);
      INTERP(xform[i].juliascope_dist);
      INTERP(xform[i].radial_blur_angle);
      INTERP(xform[i].pie_slices);
      INTERP(xform[i].pie_rotation);
      INTERP(xform[i].pie_thickness);
      INTERP(xform[i].ngon_sides);
      INTERP(xform[i].ngon_power);
      INTERP(xform[i].ngon_circle);
      INTERP(xform[i].ngon_corners);
      INTERP(xform[i].curl_c1);
      INTERP(xform[i].curl_c2);
      INTERP(xform[i].rectangles_x);
      INTERP(xform[i].rectangles_y);
      INTERP(xform[i].amw_amp);
      INTERP(xform[i].disc2_rot);
      INTERP(xform[i].disc2_twist);
      INTERP(xform[i].super_shape_rnd);
      INTERP(xform[i].super_shape_m);
      INTERP(xform[i].super_shape_n1);
      INTERP(xform[i].super_shape_n2);
      INTERP(xform[i].super_shape_n3);
      INTERP(xform[i].super_shape_holes);
      INTERP(xform[i].flower_petals);
      INTERP(xform[i].flower_holes);
      INTERP(xform[i].conic_eccentricity);
      INTERP(xform[i].conic_holes);
      INTERP(xform[i].parabola_height);
      INTERP(xform[i].parabola_width);
      INTERP(xform[i].bent2_x);
      INTERP(xform[i].bent2_y);
      INTERP(xform[i].bipolar_shift);
      INTERP(xform[i].cell_size);
      INTERP(xform[i].cpow_r);
      INTERP(xform[i].cpow_i);
      INTERP(xform[i].cpow_power);
      INTERP(xform[i].curve_xamp);
      INTERP(xform[i].curve_yamp);
      INTERP(xform[i].curve_xlength);
      INTERP(xform[i].curve_ylength);
      INTERP(xform[i].escher_beta);
      INTERP(xform[i].lazysusan_x);
      INTERP(xform[i].lazysusan_y);
      INTERP(xform[i].lazysusan_twist);
      INTERP(xform[i].lazysusan_space);
      INTERP(xform[i].lazysusan_spin);
      INTERP(xform[i].modulus_x);
      INTERP(xform[i].modulus_y);
      INTERP(xform[i].oscope_separation);
      INTERP(xform[i].oscope_frequency);
      INTERP(xform[i].oscope_amplitude);
      INTERP(xform[i].oscope_damping);
      INTERP(xform[i].popcorn2_x);
      INTERP(xform[i].popcorn2_y);
      INTERP(xform[i].popcorn2_c);
      INTERP(xform[i].separation_x);
      INTERP(xform[i].separation_xinside);
      INTERP(xform[i].separation_y);
      INTERP(xform[i].separation_yinside);
      INTERP(xform[i].split_xsize);
      INTERP(xform[i].split_ysize);
      INTERP(xform[i].splits_x);
      INTERP(xform[i].splits_y);
      INTERP(xform[i].stripes_space);
      INTERP(xform[i].stripes_warp);
      INTERP(xform[i].wedge_angle);
      INTERP(xform[i].wedge_hole);
      INTERP(xform[i].wedge_count);
      INTERP(xform[i].wedge_swirl);
      INTERP(xform[i].wedge_julia_angle);
      INTERP(xform[i].wedge_julia_count);
      INTERP(xform[i].wedge_julia_power);
      INTERP(xform[i].wedge_julia_dist);
      INTERP(xform[i].wedge_sph_angle);
      INTERP(xform[i].wedge_sph_hole);
      INTERP(xform[i].wedge_sph_count);
      INTERP(xform[i].wedge_sph_swirl);
      INTERP(xform[i].whorl_inside);
      INTERP(xform[i].whorl_outside);
      INTERP(xform[i].waves2_scalex);
      INTERP(xform[i].waves2_scaley);
      INTERP(xform[i].waves2_freqx);
      INTERP(xform[i].waves2_freqy);
      INTERP(xform[i].auger_sym);
      INTERP(xform[i].auger_weight);
      INTERP(xform[i].auger_freq);
      INTERP(xform[i].auger_scale);
      INTERP(xform[i].flux_spread);
      INTERP(xform[i].mobius_re_a);
      INTERP(xform[i].mobius_im_a);
      INTERP(xform[i].mobius_re_b);
      INTERP(xform[i].mobius_im_b);
      INTERP(xform[i].mobius_re_c);
      INTERP(xform[i].mobius_im_c);
      INTERP(xform[i].mobius_re_d);
      INTERP(xform[i].mobius_im_d);

      for (j = 0; j < flam3_nvariations; j++)
         INTERP(xform[i].var[j]);

      if (flam3_inttype_log == cpi[0].interpolation_type) {
         double cxmag[4][2];  // XXX why only 4? should be ncp
         double cxang[4][2];
         double cxtrn[4][2];

         /* affine part */
         clear_matrix(result->xform[i].c);
         convert_linear_to_polar(cpi,ncp,i,0,cxang,cxmag,cxtrn);
         interp_and_convert_back(c, ncp, i, cxang, cxmag, cxtrn,result->xform[i].c);

         /* post part */
         all_id = 1;
         for (k=0; k<ncp; k++)
            all_id &= id_matrix(cpi[k].xform[i].post);
         
         clear_matrix(result->xform[i].post);
         if (all_id) {
            result->xform[i].post[0][0] = 1.0;
            result->xform[i].post[1][1] = 1.0;
         } else {
            convert_linear_to_polar(cpi,ncp,i,1,cxang,cxmag,cxtrn);
            interp_and_convert_back(c, ncp, i, cxang, cxmag, cxtrn,result->xform[i].post);
         }         
         
      } else {

         /* Interpolate c matrix & post */
         clear_matrix(result->xform[i].c);
         clear_matrix(result->xform[i].post);
         all_id = 1;
         for (k = 0; k < ncp; k++) {
            sum_matrix(c[k], cpi[k].xform[i].c, result->xform[i].c);
            sum_matrix(c[k], cpi[k].xform[i].post, result->xform[i].post);

            all_id &= id_matrix(cpi[k].xform[i].post);

         }
         if (all_id) {
            clear_matrix(result->xform[i].post);
            result->xform[i].post[0][0] = 1.0;
            result->xform[i].post[1][1] = 1.0;
         }
      }
      
      if (ncp==2 && stagger>0 && i!=cpi[0].final_xform_index) {
         c[0] = csave[0];
         c[1] = csave[1];
      }
      
   }
   
}

void establish_asymmetric_refangles(flam3_genome *cp, int ncps) {

   int k, xfi, col;
   
   double cxang[4][2],d,c1[2];

   for (xfi=0; xfi<cp[0].num_xforms; xfi++) {
   
      /* Final xforms don't rotate regardless of their symmetry */
      if (cp[0].final_xform_enable==1 && xfi==cp[0].final_xform_index)
         continue;

      for (k=0; k<ncps;k++) {

         /* Establish the angle for each component */
         /* Should potentially functionalize */
         for (col=0;col<2;col++) {
         
            c1[0] = cp[k].xform[xfi].c[col][0];
            c1[1] = cp[k].xform[xfi].c[col][1];
            
            cxang[k][col] = atan2(c1[1],c1[0]);
         }
      }
      
      for (k=1; k<ncps; k++) {
     
         for (col=0;col<2;col++) {

            int sym0,sym1;
            int padsymflag;

            d = cxang[k][col]-cxang[k-1][col];

            /* Adjust to avoid the -pi/pi discontinuity */
            if (d > M_PI+EPS)
            cxang[k][col] -= 2*M_PI;
            else if (d < -(M_PI-EPS) )
            cxang[k][col] += 2*M_PI;

            /* If this is an asymmetric case, store the NON-symmetric angle    */
            /* Check them pairwise and store the reference angle in the second */
            /* to avoid overwriting if asymmetric on both sides                */
            padsymflag = 0;
         
            sym0 = (cp[k-1].xform[xfi].animate==0 || (cp[k-1].xform[xfi].padding==1 && padsymflag));
            sym1 = (cp[k].xform[xfi].animate==0 || (cp[k].xform[xfi].padding==1 && padsymflag));

            if ( sym1 && !sym0 )
               cp[k].xform[xfi].wind[col] = cxang[k-1][col] + 2*M_PI;
            else if ( sym0 && !sym1 )
               cp[k].xform[xfi].wind[col] = cxang[k][col] + 2*M_PI;

         }
      }
   }
}

void flam3_align(flam3_genome *dst, flam3_genome *src, int nsrc) {
   int i, tfx, tnx, max_nx = 0, max_fx = 0;
   int already_aligned=1;
   int xf,j;
   int ii,fnd;
   double normed;
   int usethisone;
   
   usethisone = (nsrc/2) - 1;
   
   max_nx = src[0].num_xforms - (src[0].final_xform_index >= 0);
   max_fx = src[0].final_xform_enable;
   
   for (i = 1; i < nsrc; i++) {
      tnx = src[i].num_xforms - (src[i].final_xform_index >= 0);
      if (max_nx != tnx) {
         already_aligned = 0;
         if (tnx > max_nx) max_nx = tnx;
      }
      
      tfx = src[i].final_xform_enable;
      if (max_fx != tfx) {
         already_aligned = 0;
         max_fx |= tfx;
      }
   }

   /* Pad the cps to equal xforms */
   for (i = 0; i < nsrc; i++) {
      flam3_copyx(&dst[i], &src[i], max_nx, max_fx);
   }

   /* Skip if this genome is compatibility mode */
   if (dst[usethisone].interpolation_type == flam3_inttype_compat ||
         dst[usethisone].interpolation_type == flam3_inttype_older)
      return;

      
   /* Check to see if there's a parametric variation present in one xform   */
   /* but not in an aligned xform.  If this is the case, use the parameters */
   /* from the xform with the variation as the defaults for the blank one.  */
   
   /* All genomes will have the same number of xforms at this point */
   /* num = max_nx + max_fx */
   for (i = 0; i<nsrc; i++) {


      for (xf = 0; xf<max_nx+max_fx; xf++) {
                  
         /* Loop over the variations to see which of them are set to 0 */
         /* Note that there are no parametric variations < 23 */
         for (j = 23; j < flam3_nvariations; j++) {
         
              if (dst[i].xform[xf].var[j]==0) {
            
                 if (i>0) {
                              
                    /* Check to see if the prior genome's xform is populated */
                    if (dst[i-1].xform[xf].var[j] != 0) {
                  
                       /* Copy the prior genome's parameters and continue */
                       flam3_copy_params(&(dst[i].xform[xf]), &(dst[i-1].xform[xf]), j);
                       continue;
                    }

                 } else if (i<nsrc-1) {

                    /* Check to see if the next genome's xform is populated */
                    if (dst[i+1].xform[xf].var[j] != 0) {
                  
                       /* Copy the next genome's parameters and continue */
                       flam3_copy_params(&(dst[i].xform[xf]), &(dst[i+1].xform[xf]), j);
                       continue;
                    }
                 }
              }
          } /* variations */

          if (dst[i].xform[xf].padding == 1 && !already_aligned) {
         
             /* This is a new xform.  Let's see if we can choose a better 'identity' xform. */
             /* Check the neighbors to see if any of these variations are used: */
             /* rings2, fan2, blob, perspective, julian, juliascope, ngon, curl, super_shape, split */
             /* If so, we can use a better starting point for these */
            
             /* Remove linear from the list */
             dst[i].xform[xf].var[0] = 0.0;
            
             /* Look through all of the 'companion' xforms to see if we get a match on any of these */
             fnd=0;

             /* Only do the next substitution for log interpolation */
             if ( (i==0 && dst[i].interpolation_type == flam3_inttype_log)
                  || (i>0 && dst[i-1].interpolation_type==flam3_inttype_log) ) {

             for (ii=-1; ii<=1; ii+=2) {

                /* Skip if out of bounds */
                if (i+ii<0 || i+ii>=nsrc)
                   continue;
                  
                /* Skip if this is also padding */
                if (dst[i+ii].xform[xf].padding==1)
                   continue;

                /* Spherical / Ngon (trumps all others due to holes)       */
                /* Interpolate these against a 180 degree rotated identity */
                /* with weight -1.                                         */
                /* Added JULIAN/JULIASCOPE to get rid of black wedges      */
                if (dst[i+ii].xform[xf].var[VAR_SPHERICAL]>0 ||
                      dst[i+ii].xform[xf].var[VAR_NGON]>0 || 
                      dst[i+ii].xform[xf].var[VAR_JULIAN]>0 || 
                      dst[i+ii].xform[xf].var[VAR_JULIASCOPE]>0 ||
                      dst[i+ii].xform[xf].var[VAR_POLAR]>0 ||
                      dst[i+ii].xform[xf].var[VAR_WEDGE_SPH]>0 ||
                      dst[i+ii].xform[xf].var[VAR_WEDGE_JULIA]>0) {
                 
                   dst[i].xform[xf].var[VAR_LINEAR] = -1.0;
                   /* Set the coefs appropriately */
                   dst[i].xform[xf].c[0][0] = -1.0;
                   dst[i].xform[xf].c[0][1] = 0.0;
                   dst[i].xform[xf].c[1][0] = 0.0;
                   dst[i].xform[xf].c[1][1] = -1.0;
                   dst[i].xform[xf].c[2][0] = 0.0;
                   dst[i].xform[xf].c[2][1] = 0.0;               
                   fnd=-1;
                }
             }

             }

             if (fnd==0) {

                for (ii=-1; ii<=1; ii+=2) {

                   /* Skip if out of bounds */
                   if (i+ii<0 || i+ii>=nsrc)
                      continue;
                     
                   /* Skip if also padding */
                   if (dst[i+ii].xform[xf].padding==1)
                      continue;

                   /* Rectangles */
                   if (dst[i+ii].xform[xf].var[VAR_RECTANGLES]>0) {
                      dst[i].xform[xf].var[VAR_RECTANGLES] = 1.0;
                      dst[i].xform[xf].rectangles_x = 0.0;
                      dst[i].xform[xf].rectangles_y = 0.0;
                      fnd++;
                   }

                   /* Rings 2 */
                   if (dst[i+ii].xform[xf].var[VAR_RINGS2]>0) {
                      dst[i].xform[xf].var[VAR_RINGS2] = 1.0;
                      dst[i].xform[xf].rings2_val = 0.0;
                      fnd++;
                   }
                  
                   /* Fan 2 */
                   if (dst[i+ii].xform[xf].var[VAR_FAN2]>0) {
                      dst[i].xform[xf].var[VAR_FAN2] = 1.0;
                      dst[i].xform[xf].fan2_x = 0.0;
                      dst[i].xform[xf].fan2_y = 0.0;
                      fnd++;
                   }
               
                   /* Blob */
                   if (dst[i+ii].xform[xf].var[VAR_BLOB]>0) {
                      dst[i].xform[xf].var[VAR_BLOB] = 1.0;
                      dst[i].xform[xf].blob_low = 1.0;
                      dst[i].xform[xf].blob_high = 1.0;
                      dst[i].xform[xf].blob_waves = 1.0;
                      fnd++;
                   }
               
                   /* Perspective */
                   if (dst[i+ii].xform[xf].var[VAR_PERSPECTIVE]>0) {
                      dst[i].xform[xf].var[VAR_PERSPECTIVE] = 1.0;
                      dst[i].xform[xf].perspective_angle = 0.0;
                      /* Keep the perspective distance as-is */
                      fnd++;
                   }
               
                   /* Curl */
                   if (dst[i+ii].xform[xf].var[VAR_CURL]>0) {
                      dst[i].xform[xf].var[VAR_CURL] = 1.0;
                      dst[i].xform[xf].curl_c1 = 0.0;
                      dst[i].xform[xf].curl_c2 = 0.0;
                      fnd++;
                   }

                   /* Super-Shape */
                   if (dst[i+ii].xform[xf].var[VAR_SUPER_SHAPE]>0) {
                      dst[i].xform[xf].var[VAR_SUPER_SHAPE] = 1.0;
                      /* Keep supershape_m the same */
                      dst[i].xform[xf].super_shape_n1 = 2.0;
                      dst[i].xform[xf].super_shape_n2 = 2.0;
                      dst[i].xform[xf].super_shape_n3 = 2.0;
                      dst[i].xform[xf].super_shape_rnd = 0.0;
                      dst[i].xform[xf].super_shape_holes = 0.0;
                      fnd++;
                   }
                }
             }

             /* If we didn't have any matches with those, */
             /* try the affine ones, fan and rings        */
             if (fnd==0) {
            
                for (ii=-1; ii<=1; ii+=2) {

                   /* Skip if out of bounds */
                   if (i+ii<0 || i+ii>=nsrc)
                      continue;                  

                   /* Skip if also a padding xform */
                   if (dst[i+ii].xform[xf].padding==1)
                      continue;
                     
                   /* Fan */
                   if (dst[i+ii].xform[xf].var[VAR_FAN]>0) {
                      dst[i].xform[xf].var[VAR_FAN] = 1.0;
                      fnd++;
                   }

                   /* Rings */
                   if (dst[i+ii].xform[xf].var[VAR_RINGS]>0) {
                      dst[i].xform[xf].var[VAR_RINGS] = 1.0;
                      fnd++;
                   }

                }
               
                if (fnd>0) {
                   /* Set the coefs appropriately */
                   dst[i].xform[xf].c[0][0] = 0.0;
                   dst[i].xform[xf].c[0][1] = 1.0;
                   dst[i].xform[xf].c[1][0] = 1.0;
                   dst[i].xform[xf].c[1][1] = 0.0;
                   dst[i].xform[xf].c[2][0] = 0.0;
                   dst[i].xform[xf].c[2][1] = 0.0;               
                }
             }
                                          
             /* If we still have no matches, switch back to linear */
             if (fnd==0)

                dst[i].xform[xf].var[VAR_LINEAR] = 1.0;

             else if (fnd>0) {

                /* Otherwise, go through and normalize the weights. */
                normed = 0.0;
                for (j = 0; j < flam3_nvariations; j++)
                   normed += dst[i].xform[xf].var[j];
                  
                for (j = 0; j < flam3_nvariations; j++)
                   dst[i].xform[xf].var[j] /= normed;

             }         
          }
       } /* xforms */
   } /* genomes */
                              
}

