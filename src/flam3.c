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

#ifndef _MSC_VER /* VC++ */
#define _GNU_SOURCE
#endif

#include "private.h"
#include "img.h"
#include "config.h"
#include "variations.h"
#include "interpolation.h"
#include "parser.h"
#include "filters.h"
#include "palettes.h"
#include <limits.h>
#include <locale.h>
#include <math.h>
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <errno.h>

#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif

#ifdef __APPLE__
#include <mach/mach.h>
#include <mach/mach_error.h>
#define flam3_os "OSX"
#else
#ifdef _WIN32
#define WINVER 0x0500
#include <windows.h>
#define flam3_os "WIN"
#else
#define flam3_os "LNX"
#endif
#endif



char *flam3_version() {

  if (strcmp(SVN_REV, "exported"))
    return flam3_os "-" VERSION "." SVN_REV;
  return flam3_os "-" VERSION;
}


#define CHOOSE_XFORM_GRAIN 16384
#define CHOOSE_XFORM_GRAIN_M1 16383

#define random_distrib(v) ((v)[random()%vlen(v)])
                                   
                                    
unsigned short *flam3_create_xform_distrib(flam3_genome *cp) {

   /* Xform distrib is created in this function             */   
   int numrows;
   int dist_row,i;
   unsigned short *xform_distrib;
   
   numrows = cp->num_xforms - (cp->final_xform_index>=0) + 1;
   xform_distrib = calloc(numrows*CHOOSE_XFORM_GRAIN,sizeof(unsigned short));
   
   /* First, set up the first row of the xform_distrib (raw weights) */
   flam3_create_chaos_distrib(cp, -1, xform_distrib);
   
   /* Check for non-unity chaos */
   cp->chaos_enable = 1 - flam3_check_unity_chaos(cp);
   
   if (cp->chaos_enable) {
   
      /* Now set up a row for each of the xforms */
      dist_row = 0;
      for (i=0;i<cp->num_xforms;i++) {
      
         if (cp->final_xform_index == i)
            continue;
         else
            dist_row++;
         
         if (flam3_create_chaos_distrib(cp, i, &(xform_distrib[CHOOSE_XFORM_GRAIN*(dist_row)]))) {
            free(xform_distrib);
            return(NULL);
         }
      }
   }
   
   return(xform_distrib);
}

/* Run this on a copy of a genome to get a strip of the output */
int flam3_make_strip(flam3_genome *cp, int nstrips, int stripnum) {

    double old_center[2];
    int j;

    /* Strip out motion elements */         
    for (j=0;j<cp->num_xforms;j++)
    flam3_delete_motion_elements(&cp->xform[j]);

    old_center[0] = cp->center[0];
    old_center[1] = cp->center[1];
    cp->height /= nstrips;
    cp->center[1] = cp->center[1] - ((nstrips - 1) * cp->height) / (2 * cp->pixels_per_unit * pow(2.0, cp->zoom));
    cp->center[1] += cp->height * stripnum / ( cp->pixels_per_unit * pow(2.0,cp->zoom) );
    rotate_by(cp->center, old_center, cp->rotate);
    
    return(0);
}


void rotate_by(double *p, double *center, double by) {
    double r[2];
    double th = by * 2 * M_PI / 360.0;
    double c = cos(th);
    double s = -sin(th);
    p[0] -= center[0];
    p[1] -= center[1];
    r[0] = c * p[0] - s * p[1];
    r[1] = s * p[0] + c * p[1];
    p[0] = r[0] + center[0];
    p[1] = r[1] + center[1];
}


int flam3_check_unity_chaos(flam3_genome *cp) {

   int i,j;
   int num_std;
   int unity=1;
   num_std = cp->num_xforms - (cp->final_xform_index >= 0);
   
   for (i=0;i<num_std;i++) {
      for (j=0;j<num_std;j++) {
         if ( fabs(cp->chaos[i][j]-1.0) > EPS)
            unity=0;
      }
   }
   
   return(unity);
}

int flam3_create_chaos_distrib(flam3_genome *cp, int xi, unsigned short *xform_distrib) {

   /* Xform distrib is a preallocated array of CHOOSE_XFORM_GRAIN chars */
   /* address of array is passed in, contents are modified              */
   double t,r,dr;
   int i,j;
   int num_std;
   
   //fprintf(stdout,"storing at %ld\n",xform_distrib);
   
   num_std = cp->num_xforms - (cp->final_xform_index >= 0);
   
   dr = 0.0;
   for (i = 0; i < num_std; i++) {
      double d = cp->xform[i].density;
                     
      if (xi>=0)
         d *= cp->chaos[xi][i];
         
      //fprintf(stdout,"%f ",d);
      if (d < 0.0) {
         fprintf(stderr, "xform weight must be non-negative, not %g.\n", d);
         return(1);
      }         
      
      dr += d;
   }
   
   //fprintf(stdout,"dr=%f\n",dr);
   
   if (dr == 0.0) {
      fprintf(stderr, "cannot iterate empty flame.\n");
      return(1);
   }
   
   dr = dr / CHOOSE_XFORM_GRAIN;

   j = 0;
   t = cp->xform[0].density;
   if (xi>=0)
     t *= cp->chaos[xi][0];
   r = 0.0;
   for (i = 0; i < CHOOSE_XFORM_GRAIN; i++) {
      while (r >= t) {
         j++;

         if (xi>=0)
            t += cp->xform[j].density*cp->chaos[xi][j];
         else
            t += cp->xform[j].density;
         
      }
      //fprintf(stdout,"%d ",j);
      xform_distrib[i] = j;
      r += dr;
   }
   //fprintf(stdout,"\n---\n");
   
   return(0);
}
/*
 * run the function system described by CP forward N generations.  store
 * the N resulting 4-vectors in SAMPLES.  the initial point is passed in
 * SAMPLES[0..3].  ignore the first FUSE iterations.
 */


int flam3_iterate(flam3_genome *cp, int n, int fuse,  double *samples, unsigned short *xform_distrib, randctx *rc) {
   int i;
   double p[4], q[4];
   int consec = 0;
   int badvals = 0;
   int lastxf=0;
   int fn;
   
   p[0] = samples[0];
   p[1] = samples[1];
   p[2] = samples[2];
   p[3] = samples[3];

   /* Perform precalculations */   
   for (i=0;i<cp->num_xforms;i++)
      xform_precalc(cp,i);

   for (i = -4*fuse; i < 4*n; i+=4) {
   
//         fn = xform_distrib[ lastxf*CHOOSE_XFORM_GRAIN + (((unsigned)irand(rc)) % CHOOSE_XFORM_GRAIN)];
      if (cp->chaos_enable)
         fn = xform_distrib[ lastxf*CHOOSE_XFORM_GRAIN + (((unsigned)irand(rc)) & CHOOSE_XFORM_GRAIN_M1)];
      else
         fn = xform_distrib[ ((unsigned)irand(rc)) & CHOOSE_XFORM_GRAIN_M1 ];
      
      if (apply_xform(cp, fn, p, q, rc)>0) {
         consec ++;
         badvals ++;
         if (consec<5) {
            p[0] = q[0];
            p[1] = q[1];
            p[2] = q[2];
            p[3] = q[3];
            i -= 4;
            continue;
         } else
            consec = 0;
      } else
         consec = 0;

      /* Store the last used transform */
      lastxf = fn+1;

      p[0] = q[0];
      p[1] = q[1];
      p[2] = q[2];
      p[3] = q[3];

      if (cp->final_xform_enable == 1) {
         if (cp->xform[cp->final_xform_index].opacity==1 || 
                flam3_random_isaac_01(rc)<cp->xform[cp->final_xform_index].opacity) {
             apply_xform(cp, cp->final_xform_index, p, q, rc);
             /* Keep the opacity from the original xform */
             q[3] = p[3];
         }
      }

      /* if fuse over, store it */
      if (i >= 0) {
         samples[i] = q[0];
         samples[i+1] = q[1];
         samples[i+2] = q[2];
         samples[i+3] = q[3];
      }
   }
   
   return(badvals);
}

int flam3_xform_preview(flam3_genome *cp, int xi, double range, int numvals, int depth, double *result, randctx *rc) {

   /* We will evaluate the 'xi'th xform 'depth' times, over the following values:           */
   /* x in [-range : range], y in [-range : range], with 2* (2*numvals+1)^2 values returned */ 
   double p[4];
   double incr;
   int outi;
   int xx,yy,dd;
   double oldweight;
   
   outi=0;
   
   oldweight = cp->xform[xi].density;
   cp->xform[xi].density = 1.0;
   
   /* Prepare the function pointers */
   if (prepare_precalc_flags(cp)) {
      cp->xform[xi].density = oldweight;
      return(1);
   }
   
   /* Calculate increment */
   incr = range / (double)numvals;
   
   /* Perform precalculations */
   xform_precalc(cp,xi);
   
   /* Loop over the grid */
   for (xx=-numvals;xx<=numvals;xx++) {
      for (yy=-numvals;yy<=numvals;yy++) {
      
         /* Calculate the input coordinates */
         p[0] = (double)xx * incr;
         p[1] = (double)yy * incr;
         
         /* Loop over the depth */
         for (dd=0;dd<depth;dd++)
            apply_xform(cp, xi, p, p, rc);
         
         result[outi] = p[0];
         result[outi+1] = p[1];
         
         outi += 2;
      }
   }
   cp->xform[xi].density = oldweight;

   return(0);
}         

int flam3_colorhist(flam3_genome *cp, int num_batches, randctx *rc, double *hist) {

  int lp,plp;
  int mycolor;
  unsigned short *xform_distrib;
  int sbs = 10000;
  double sub_batch[4*10000];

  memset(hist,0,256*sizeof(double));
  
  // get into the attractor
  if (prepare_precalc_flags(cp))
     return(1);
      
  xform_distrib = flam3_create_xform_distrib(cp);

  for (lp=0;lp<num_batches;lp++) {
  
    sub_batch[0] = flam3_random_isaac_11(rc);
    sub_batch[1] = flam3_random_isaac_11(rc);
    sub_batch[2] = 0;
    sub_batch[3] = 0;

    
    if (xform_distrib==NULL)
       return(1);
    flam3_iterate(cp, sbs, 20, sub_batch, xform_distrib, rc);
    
    // histogram the colors in the sub_batch array
    for (plp=0;plp<4*sbs;plp+=4) {
      mycolor = (int)(sub_batch[plp+2]*CMAP_SIZE);
      if (mycolor<0) mycolor=0;
      if (mycolor>CMAP_SIZE_M1) mycolor=CMAP_SIZE_M1;
      
      hist[mycolor] += 1;
    }
  }
  
  free(xform_distrib);
  for (plp=0;plp<256;plp++)
    hist[plp] /= (float)(num_batches*sbs);
    
  return(0);
} 
  
flam3_genome *sheep_loop(flam3_genome *cp, double blend) {
   
   flam3_genome *result;
   int i;

   /* Allocate the genome - this must be freed by calling function */
   result = calloc(1,sizeof(flam3_genome));
   
   /* Clear it */
   clear_cp(result,flam3_defaults_on);
   
   /* Copy the original */
   flam3_copy(result,cp);

   /*
    * Insert motion magic here :
    * if there are motion elements, we will modify the contents of
    * the result genome before flam3_rotate is called.
    */
   for (i=0;i<cp->num_xforms;i++) {
      if (cp->xform[i].num_motion>0) {
         /* Apply motion parameters to result.xform[i] using blend parameter */
         apply_motion_parameters(&cp->xform[i], &result->xform[i], blend);
      }
      
      /* Delete the motion parameters from the result */
      flam3_delete_motion_elements(&result->xform[i]);
   }

   /* Rotate the affines */
   flam3_rotate(result, blend*360.0,result->interpolation_type);
   
   return(result);
}




flam3_genome *sheep_edge(flam3_genome *cp, double blend, int seqflag, double stagger) {

   flam3_genome spun[2];
   flam3_genome prealign[2];
   flam3_genome *result;
   int i,si;
   char *ai;

   memset(spun, 0, 2*sizeof(flam3_genome));
   memset(prealign, 0, 2*sizeof(flam3_genome));

   /* Allocate the memory for the result */
   result = calloc(1,sizeof(flam3_genome));

   /*
    * Insert motion magic here :
    * if there are motion elements, we will modify the contents of
    * the prealign genomes before we rotate and interpolate.
    */
  
   for (si=0;si<2;si++) {
      flam3_copy(&prealign[si], &cp[si]);
      for (i=0;i<cp[si].num_xforms;i++) {
         if (cp[si].xform[i].num_motion>0) {
            /* Apply motion parameters to result.xform[i] using blend parameter */
            apply_motion_parameters(&cp[si].xform[i], &prealign[si].xform[i], blend);
         }
      }
   }

   /* Use the un-padded original for blend=0 when creating a sequence */
   /* This keeps the original interpolation type intact               */
   if (seqflag && 0.0 == blend) {
      flam3_copy(result, &prealign[0]);
   } else {

      /* Align what we're going to interpolate */
      flam3_align(spun, prealign, 2);

      spun[0].time = 0.0;
      spun[1].time = 1.0;

      /* Call this first to establish the asymmetric reference angles */
      establish_asymmetric_refangles(spun,2);  

      /* Rotate the aligned xforms */
      flam3_rotate(&spun[0], blend*360.0, spun[0].interpolation_type);
      flam3_rotate(&spun[1], blend*360.0, spun[0].interpolation_type);

      /* Now call the interpolation */
      if (argi("unsmoother",0) == 0)
         flam3_interpolate(spun, 2, smoother(blend), stagger, result);
      else
         flam3_interpolate(spun, 2, blend, stagger, result);
      

     /* Interpolation type no longer needs to be forced to linear mode */
//     if (!seqflag)
//        result.interpolation_type = flam3_inttype_linear;
   }
   
   /* Clear the genomes we used */
   clear_cp(&spun[0],flam3_defaults_on);
   clear_cp(&spun[1],flam3_defaults_on);
   clear_cp(&prealign[0],flam3_defaults_on);
   clear_cp(&prealign[1],flam3_defaults_on);

   /* Make sure there are no motion elements in the result */
   for (i=0;i<result->num_xforms;i++)
      flam3_delete_motion_elements(&result->xform[i]);

   return(result);
}


/* BY is angle in degrees */
void flam3_rotate(flam3_genome *cp, double by, int interpolation_type) {
   int i;
   for (i = 0; i < cp->num_xforms; i++) {
      double r[2][2];
      double T[2][2];
      double U[2][2];
      double dtheta = by * 2.0 * M_PI / 360.0;

      /* Don't rotate xforms with > 0 animate values */
      if (cp->xform[i].animate == 0.0)
         continue;

      if (cp->xform[i].padding == 1) {
         if (interpolation_type == flam3_inttype_compat) {
            /* gen 202 era flam3 did not rotate padded xforms */
            continue;
         } else if (interpolation_type == flam3_inttype_older) {
            /* not sure if 198 era flam3 rotated padded xforms */
            continue;
         } else if (interpolation_type == flam3_inttype_linear) {
            /* don't rotate for prettier symsings */
            continue;
         } else if (interpolation_type == flam3_inttype_log) {
            /* Current flam3: what do we prefer? */
            //continue;
         }
      }

      /* Do NOT rotate final xforms */
      if (cp->final_xform_enable==1 && cp->final_xform_index==i)
         continue;

      r[1][1] = r[0][0] = cos(dtheta);
      r[0][1] = sin(dtheta);
      r[1][0] = -r[0][1];
      T[0][0] = cp->xform[i].c[0][0];
      T[1][0] = cp->xform[i].c[1][0];
      T[0][1] = cp->xform[i].c[0][1];
      T[1][1] = cp->xform[i].c[1][1];
      mult_matrix(r, T, U);
      cp->xform[i].c[0][0] = U[0][0];
      cp->xform[i].c[1][0] = U[1][0];
      cp->xform[i].c[0][1] = U[0][1];
      cp->xform[i].c[1][1] = U[1][1];
   }
}

#define APPMOT(x)  do { addto->x += mot[i].x * motion_funcs(func,freq*blend); } while (0);

void apply_motion_parameters(flam3_xform *xf, flam3_xform *addto, double blend) {

   int i,j,k;
   int freq;
   int func;
   flam3_xform* mot;
   
   mot = xf->motion;

   /* Loop over the motion elements and add their contribution to the original vals */
   for (i=0; i<xf->num_motion; i++) {   
   
      freq = mot->motion_freq;
      func = mot->motion_func;
      
      APPMOT(density); /* Must ensure > 0 after all is applied */
      APPMOT(color); /* Must ensure [0,1] after all is applied */
      
      APPMOT(opacity);      
      APPMOT(color_speed);
      APPMOT(animate);
      APPMOT(blob_low);
      APPMOT(blob_high);
      APPMOT(blob_waves);
      APPMOT(pdj_a);
      APPMOT(pdj_b);
      APPMOT(pdj_c);
      APPMOT(pdj_d);
      APPMOT(fan2_x);
      APPMOT(fan2_y);
      APPMOT(rings2_val);
      APPMOT(perspective_angle);
      APPMOT(perspective_dist);
      APPMOT(julian_power);
      APPMOT(julian_dist);
      APPMOT(juliascope_power);
      APPMOT(juliascope_dist);
      APPMOT(radial_blur_angle);
      APPMOT(pie_slices);
      APPMOT(pie_rotation);
      APPMOT(pie_thickness);
      APPMOT(ngon_sides);
      APPMOT(ngon_power);
      APPMOT(ngon_circle);
      APPMOT(ngon_corners);
      APPMOT(curl_c1);
      APPMOT(curl_c2);
      APPMOT(rectangles_x);
      APPMOT(rectangles_y);
      APPMOT(amw_amp);
      APPMOT(disc2_rot);
      APPMOT(disc2_twist);
      APPMOT(super_shape_rnd);
      APPMOT(super_shape_m);
      APPMOT(super_shape_n1);
      APPMOT(super_shape_n2);
      APPMOT(super_shape_n3);
      APPMOT(super_shape_holes);
      APPMOT(flower_petals);
      APPMOT(flower_holes);
      APPMOT(conic_eccentricity);
      APPMOT(conic_holes);
      APPMOT(parabola_height);
      APPMOT(parabola_width);
      APPMOT(bent2_x);
      APPMOT(bent2_y);
      APPMOT(bipolar_shift);
      APPMOT(cell_size);
      APPMOT(cpow_r);
      APPMOT(cpow_i);
      APPMOT(cpow_power);
      APPMOT(curve_xamp);
      APPMOT(curve_yamp);
      APPMOT(curve_xlength);
      APPMOT(curve_ylength);
      APPMOT(escher_beta);
      APPMOT(lazysusan_x);
      APPMOT(lazysusan_y);
      APPMOT(lazysusan_twist);
      APPMOT(lazysusan_space);
      APPMOT(lazysusan_spin);
      APPMOT(modulus_x);
      APPMOT(modulus_y);
      APPMOT(oscope_separation);
      APPMOT(oscope_frequency);
      APPMOT(oscope_amplitude);
      APPMOT(oscope_damping);
      APPMOT(popcorn2_x);
      APPMOT(popcorn2_y);
      APPMOT(popcorn2_c);
      APPMOT(separation_x);
      APPMOT(separation_xinside);
      APPMOT(separation_y);
      APPMOT(separation_yinside);
      APPMOT(split_xsize);
      APPMOT(split_ysize);
      APPMOT(splits_x);
      APPMOT(splits_y);
      APPMOT(stripes_space);
      APPMOT(stripes_warp);
      APPMOT(wedge_angle);
      APPMOT(wedge_hole);
      APPMOT(wedge_count);
      APPMOT(wedge_swirl);
      APPMOT(wedge_julia_angle);
      APPMOT(wedge_julia_count);
      APPMOT(wedge_julia_power);
      APPMOT(wedge_julia_dist);
      APPMOT(wedge_sph_angle);
      APPMOT(wedge_sph_hole);
      APPMOT(wedge_sph_count);
      APPMOT(wedge_sph_swirl);
      APPMOT(whorl_inside);
      APPMOT(whorl_outside);
      APPMOT(waves2_scalex);
      APPMOT(waves2_scaley);
      APPMOT(waves2_freqx);
      APPMOT(waves2_freqy);
      APPMOT(auger_sym);
      APPMOT(auger_weight);
      APPMOT(auger_freq);
      APPMOT(auger_scale);
      APPMOT(flux_spread);
      APPMOT(mobius_re_a);
      APPMOT(mobius_re_b);
      APPMOT(mobius_re_c);
      APPMOT(mobius_re_d);
      APPMOT(mobius_im_a);
      APPMOT(mobius_im_b);
      APPMOT(mobius_im_c);
      APPMOT(mobius_im_d);

      for (j = 0; j < flam3_nvariations; j++)
         APPMOT(var[j]);
         
      for (j=0; j<3; j++) {
         for (k=0; k<2; k++) {
            APPMOT(c[j][k]);
            APPMOT(post[j][k]);
         }
      }
         
   }
   
   /* Make sure certain params are within reasonable bounds */
   if (addto->color<0) addto->color=0;
   if (addto->color>1) addto->color=1;
   if (addto->density<0) addto->density=0;
   
}
      

/*
 * create a control point that interpolates between the control points
 * passed in CPS.  CPS must be sorted by time.
 */
void flam3_interpolate(flam3_genome cps[], int ncps,
             double time, double stagger, flam3_genome *result) {
   int i1, i2;
   double c[2];
   flam3_genome cpi[4];
   int smoothflag = 0;

   if (1 == ncps) {
      flam3_copy(result, &(cps[0]));
      return;
   }
      
   if (cps[0].time >= time) {
      i1 = 0;
      i2 = 1;
   } else if (cps[ncps - 1].time <= time) {
      i1 = ncps - 2;
      i2 = ncps - 1;
   } else {
      i1 = 0;
      while (cps[i1].time < time)
         i1++;

      i1--;
      i2 = i1 + 1;

   }

   c[0] = (cps[i2].time - time) / (cps[i2].time - cps[i1].time);
   c[1] = 1.0 - c[0];

   memset(cpi, 0, 4*sizeof(flam3_genome));
   
   /* To interpolate the xforms, we will make copies of the source cps  */
   /* and ensure that they both have the same number before progressing */
   if (flam3_interpolation_linear == cps[i1].interpolation) {
       flam3_align(&cpi[0], &cps[i1], 2);
       smoothflag = 0;
     
   } else {
       if (0 == i1) {
          fprintf(stderr, "error: cannot use smooth interpolation on first segment.\n");
          fprintf(stderr, "reverting to linear interpolation.\n");
          flam3_align(&cpi[0], &cps[i1], 2);
          smoothflag = 0;
       }

       if (ncps-1 == i2) {
          fprintf(stderr, "error: cannot use smooth interpolation on last segment.\n");
          fprintf(stderr, "reverting to linear interpolation.\n");
          flam3_align(&cpi[0], &cps[i1], 2);
          smoothflag = 0;
       }

       flam3_align(&cpi[0], &cps[i1-1], 4);
       smoothflag = 1;
   }
   
   /* Clear the destination cp */
   clear_cp(result, 1);
      
   if (cpi[0].final_xform_index >= 0) {
      flam3_add_xforms(result, cpi[0].num_xforms-1, 0, 0);
      flam3_add_xforms(result, 1, 0, 1);
   } else
      flam3_add_xforms(result, cpi[0].num_xforms, 0, 0);
   

   result->time = time;
   result->interpolation = flam3_interpolation_linear;
   result->interpolation_type = cpi[0].interpolation_type;
   result->palette_interpolation = flam3_palette_interpolation_hsv;

   if (!smoothflag) {
       flam3_interpolate_n(result, 2, cpi, c, stagger);
   } else {
       interpolate_catmull_rom(cpi, c[1], result);
       clear_cp(&(cpi[2]),0);
       clear_cp(&(cpi[3]),0);
   }
   
   clear_cp(&(cpi[0]),0);
   clear_cp(&(cpi[1]),0);

}

void flam3_copy_params(flam3_xform *dest, flam3_xform *src, int varn) {

   /* We only want to copy param var coefs for this one */
   if (varn==VAR_BLOB) {
      /* Blob */
      dest->blob_low = src->blob_low;
      dest->blob_high = src->blob_high;
      dest->blob_waves = src->blob_waves;
   } else if (varn==VAR_PDJ) {
      /* PDJ */
      dest->pdj_a = src->pdj_a;
      dest->pdj_b = src->pdj_b;
      dest->pdj_c = src->pdj_c;
      dest->pdj_d = src->pdj_d;
   } else if (varn==VAR_FAN2) {
      /* Fan2 */
      dest->fan2_x = src->fan2_x;
      dest->fan2_y = src->fan2_y;
   } else if (varn==VAR_RINGS2) {
      /* Rings2 */
      dest->rings2_val = src->rings2_val;
   } else if (varn==VAR_PERSPECTIVE) {
      /* Perspective */
      dest->perspective_angle = src->perspective_angle;
      dest->perspective_dist = src->perspective_dist;
      dest->persp_vsin = src->persp_vsin;
      dest->persp_vfcos = src->persp_vfcos;
   } else if (varn==VAR_JULIAN) {
      /* Julia_N */
      dest->julian_power = src->julian_power;
      dest->julian_dist = src->julian_dist;
      dest->julian_rN = src->julian_rN;
      dest->julian_cn = src->julian_cn;
   } else if (varn==VAR_JULIASCOPE) {
      /* Julia_Scope */
      dest->juliascope_power = src->juliascope_power;
      dest->juliascope_dist = src->juliascope_dist;
      dest->juliascope_rN = src->juliascope_rN;
      dest->juliascope_cn = src->juliascope_cn;
   } else if (varn==VAR_RADIAL_BLUR) {
      /* Radial Blur */
      dest->radial_blur_angle = src->radial_blur_angle;
   } else if (varn==VAR_PIE) {
      /* Pie */
      dest->pie_slices = src->pie_slices;
      dest->pie_rotation = src->pie_rotation;
      dest->pie_thickness = src->pie_thickness;
   } else if (varn==VAR_NGON) {
      /* Ngon */
      dest->ngon_sides = src->ngon_sides;
      dest->ngon_power = src->ngon_power;
      dest->ngon_corners = src->ngon_corners;
      dest->ngon_circle = src->ngon_circle;
   } else if (varn==VAR_CURL) {
      /* Curl */
      dest->curl_c1 = src->curl_c1;
      dest->curl_c2 = src->curl_c2;
   } else if (varn==VAR_RECTANGLES) {
      /* Rect */
      dest->rectangles_x = src->rectangles_x;
      dest->rectangles_y = src->rectangles_y;
   } else if (varn==VAR_DISC2) {
      /* Disc2 */
      dest->disc2_rot = src->disc2_rot;
      dest->disc2_twist = src->disc2_twist;
   } else if (varn==VAR_SUPER_SHAPE) {
      /* Supershape */
      dest->super_shape_rnd = src->super_shape_rnd;
      dest->super_shape_m = src->super_shape_m;
      dest->super_shape_n1 = src->super_shape_n1;
      dest->super_shape_n2 = src->super_shape_n2;
      dest->super_shape_n3 = src->super_shape_n3;
      dest->super_shape_holes = src->super_shape_holes;
   } else if (varn==VAR_FLOWER) {
      /* Flower */
      dest->flower_petals = src->flower_petals;
      dest->flower_petals = src->flower_petals;
   } else if (varn==VAR_CONIC) {
      /* Conic */
      dest->conic_eccentricity = src->conic_eccentricity;
      dest->conic_holes = src->conic_holes;
   } else if (varn==VAR_PARABOLA) {
      /* Parabola */
      dest->parabola_height = src->parabola_height;
      dest->parabola_width = src->parabola_width;
   } else if (varn==VAR_BENT2) {
      /* Bent2 */
      dest->bent2_x = src->bent2_x;
      dest->bent2_y = src->bent2_y;
   } else if (varn==VAR_BIPOLAR) {
      /* Bipolar */
      dest->bipolar_shift = src->bipolar_shift;
   } else if (varn==VAR_CELL) {
      /* Cell */
      dest->cell_size = src->cell_size;
   } else if (varn==VAR_CPOW) {
      /* Cpow */
      dest->cpow_i = src->cpow_i;
      dest->cpow_r = src->cpow_r;
      dest->cpow_power = src->cpow_power;
   } else if (varn==VAR_CURVE) {
      /* Curve */
      dest->curve_xamp = src->curve_xamp;
      dest->curve_yamp = src->curve_yamp;
      dest->curve_xlength = src->curve_xlength;
      dest->curve_ylength = src->curve_ylength;
   } else if (varn==VAR_ESCHER) {
      /* Escher */
      dest->escher_beta = src->escher_beta;
   } else if (varn==VAR_LAZYSUSAN) {
      /* Lazysusan */
      dest->lazysusan_x = src->lazysusan_x;
      dest->lazysusan_y = src->lazysusan_y;
      dest->lazysusan_spin = src->lazysusan_spin;
      dest->lazysusan_space = src->lazysusan_space;
      dest->lazysusan_twist = src->lazysusan_twist;
   } else if (varn==VAR_MODULUS) {
      /* Modulus */
      dest->modulus_x = src->modulus_x;
      dest->modulus_y = src->modulus_y;
   } else if (varn==VAR_OSCILLOSCOPE) {
      /* Oscope */
      dest->oscope_separation = src->oscope_separation;
      dest->oscope_frequency = src->oscope_frequency;
      dest->oscope_amplitude = src->oscope_amplitude;
      dest->oscope_damping = src->oscope_damping;
   } else if (varn==VAR_POPCORN2) {
      /* Popcorn2 */
      dest->popcorn2_x = src->popcorn2_x;
      dest->popcorn2_y = src->popcorn2_y;
      dest->popcorn2_c = src->popcorn2_c;
   } else if (varn==VAR_SEPARATION) {
      /* Separation */
      dest->separation_x = src->separation_x;
      dest->separation_y = src->separation_y;
      dest->separation_xinside = src->separation_xinside;
      dest->separation_yinside = src->separation_yinside;
   } else if (varn==VAR_SPLIT) {
      /* Split */
      dest->split_xsize = src->split_xsize;
      dest->split_ysize = src->split_ysize;
   } else if (varn==VAR_SPLITS) {
      /* Splits */
      dest->splits_x = src->splits_x;
      dest->splits_y = src->splits_y;
   } else if (varn==VAR_STRIPES) {
      /* Stripes */
      dest->stripes_space = src->stripes_space;
      dest->stripes_warp = src->stripes_warp;
   } else if (varn==VAR_WEDGE) {
      /* Wedge */
      dest->wedge_angle = src->wedge_angle;
      dest->wedge_hole = src->wedge_hole;
      dest->wedge_count = src->wedge_count;
      dest->wedge_swirl = src->wedge_swirl;
   } else if (varn==VAR_WEDGE_JULIA) {
      /* Wedge_Julia */
      dest->wedge_julia_angle = src->wedge_julia_angle;
      dest->wedge_julia_count = src->wedge_julia_count;
      dest->wedge_julia_power = src->wedge_julia_power;
      dest->wedge_julia_dist = src->wedge_julia_dist;
      dest->wedgeJulia_cf = src->wedgeJulia_cf;
      dest->wedgeJulia_cn = src->wedgeJulia_cn;
      dest->wedgeJulia_rN = src->wedgeJulia_rN;
   } else if (varn==VAR_WEDGE_SPH) {
      /* Wedge_sph */
      dest->wedge_sph_angle = src->wedge_sph_angle;
      dest->wedge_sph_hole = src->wedge_sph_hole;
      dest->wedge_sph_count = src->wedge_sph_count;
      dest->wedge_sph_swirl = src->wedge_sph_swirl;      
   } else if (varn==VAR_WHORL) {
      /* whorl */
      dest->whorl_inside = src->whorl_inside;
      dest->whorl_outside = src->whorl_outside;
   } else if (varn==VAR_WAVES2) {
      /* waves2 */
      dest->waves2_scalex = src->waves2_scalex;
      dest->waves2_scaley = src->waves2_scaley;
      dest->waves2_freqx = src->waves2_freqx;
      dest->waves2_freqy = src->waves2_freqy;
   } else if (varn==VAR_AUGER) {
      /* auger */
      dest->auger_sym = src->auger_sym;
      dest->auger_weight = src->auger_weight;
      dest->auger_freq = src->auger_freq;
      dest->auger_scale = src->auger_scale;
   } else if (varn==VAR_FLUX) {
      /* flux */
      dest->flux_spread = src->flux_spread;
   } else if (varn==VAR_MOBIUS) {
      /* mobius */
      dest->mobius_re_a = src->mobius_re_a;
      dest->mobius_re_b = src->mobius_re_b;
      dest->mobius_re_c = src->mobius_re_c;
      dest->mobius_re_d = src->mobius_re_d;
      dest->mobius_im_a = src->mobius_im_a;
      dest->mobius_im_b = src->mobius_im_b;
      dest->mobius_im_c = src->mobius_im_c;
      dest->mobius_im_d = src->mobius_im_d;
   }
}

/* Motion support functions */
void flam3_add_motion_element(flam3_xform *xf) {

   /* Add one to the xform's count of motion elements */
   xf->num_motion++;
   
   /* Reallocate the motion storage to include the empty space */
   xf->motion = (struct xform *)realloc(xf->motion, xf->num_motion * sizeof(struct xform));
   
   /* Initialize the motion element */
   /* In this case, all elements should be set to 0 */
   memset( &(xf->motion[xf->num_motion-1]), 0, sizeof(struct xform));
   
}   
   
/* Motion support functions */
void flam3_delete_motion_elements(flam3_xform *xf) {

   /* Free the motion elements */
   if (xf->num_motion>0) {
      free(xf->motion);
      xf->num_motion = 0;
   }
   
}

/* Xform support functions */
void flam3_add_xforms(flam3_genome *thiscp, int num_to_add, int interp_padding, int final_flag) {

   int i,j;
   int old_num = thiscp->num_xforms;
   int oldstd,numstd;
   flam3_xform tmp;
   
   oldstd = thiscp->num_xforms - (thiscp->final_xform_index >= 0);
   
   /* !!! must make sure that if final_flag is specified, we don't already have a final xform! !!! */

//   if (thiscp->num_xforms > 0)
      thiscp->xform = (flam3_xform *)realloc(thiscp->xform, (thiscp->num_xforms + num_to_add) * sizeof(flam3_xform));
//   else
//      thiscp->xform = (flam3_xform *)malloc(num_to_add * sizeof(flam3_xform));

   thiscp->num_xforms += num_to_add;

   /* Initialize all the new xforms */
   initialize_xforms(thiscp, old_num);

   /* Set the padding flag for the new xforms */
   if (interp_padding) {
      for (i = old_num ; i < thiscp->num_xforms ; i++)
         thiscp->xform[i].padding=1;
   }
   
   /* If the final xform is not the last xform in the list, make it so */
   if (thiscp->final_xform_index >= 0 && thiscp->final_xform_index != thiscp->num_xforms-1) {
      tmp = thiscp->xform[thiscp->final_xform_index];
      for (i=thiscp->final_xform_index; i < thiscp->num_xforms-1; i++)
         thiscp->xform[i] = thiscp->xform[i+1];
      
      thiscp->final_xform_index = thiscp->num_xforms-1;
      thiscp->xform[thiscp->final_xform_index] = tmp;
   }
   
   if (final_flag) {
      /* Set the final xform index */
      thiscp->final_xform_enable = 1;
      thiscp->final_xform_index = thiscp->num_xforms-1;
   } else {
      /* Handle the chaos array */
      numstd = thiscp->num_xforms - (thiscp->final_xform_index>=0);
      
      /* Pad existing rows */
      for (i=0;i<oldstd;i++) {
         thiscp->chaos[i] = realloc(thiscp->chaos[i], numstd * sizeof(double));
         for (j=oldstd; j<numstd; j++)
            thiscp->chaos[i][j] = 1.0;
      }
      
      /* Add new rows */
      thiscp->chaos = realloc(thiscp->chaos,numstd * sizeof(double *));
      for (i=oldstd; i<numstd; i++) {
         thiscp->chaos[i] = malloc(numstd * sizeof(double));
         for (j=0;j<numstd;j++)
            thiscp->chaos[i][j] = 1.0;
      }
   }
}

void flam3_delete_xform(flam3_genome *thiscp, int idx_to_delete) {

   int i,j;
   int num_std = thiscp->num_xforms - (thiscp->final_xform_index >= 0);

   if (thiscp->final_xform_index != idx_to_delete) {
      /* We're going to delete the nth std xform. */

      /* Delete the nth_std row of the chaos array */
      free(thiscp->chaos[idx_to_delete]);

      /* Shift the pointers down one */
      for (i=idx_to_delete+1;i<num_std;i++)
         thiscp->chaos[i-1] = thiscp->chaos[i];
      
      /* Realloc the pointer array */
      thiscp->chaos = realloc(thiscp->chaos,(num_std-1)*sizeof(double *));
      num_std--;
      
      /* Loop over all of the rows and remove the nth_std element from them */
      for (i=0;i<num_std;i++) {
         for (j=idx_to_delete+1;j<num_std+1;j++) {
            thiscp->chaos[i][j-1] = thiscp->chaos[i][j];
         }
         /* Realloc the vector to have one less element */
         thiscp->chaos[i] = realloc(thiscp->chaos[i],num_std*sizeof(double));
         
      }
   }      
   
   /* Handle the final xform index */
   if (thiscp->final_xform_index == idx_to_delete) {
      thiscp->final_xform_index = -1;
      thiscp->final_xform_enable = 0;
   } else if (thiscp->final_xform_index > idx_to_delete) {
      thiscp->final_xform_index--;
   }

   /* Delete the motion elements of the banished xform */
   flam3_delete_motion_elements(&(thiscp->xform[idx_to_delete]));

   /* Move all of the xforms down one - this does not require manual motion xform adjustment */
   for (i=idx_to_delete; i<thiscp->num_xforms-1; i++)
      thiscp->xform[i] = thiscp->xform[i+1];

   thiscp->num_xforms--;

   /* Reduce the memory storage by one xform */
   thiscp->xform = (flam3_xform *)realloc(thiscp->xform, sizeof(flam3_xform) * thiscp->num_xforms);
   
}

void flam3_copy_xform(flam3_xform *dest, flam3_xform *src) {

   int j;
   
   /* Make sure the dest doesn't have motion already */
   if (dest->num_motion>0)
      flam3_delete_motion_elements(dest);

   /* Copy everything */
   *dest = *src;
   
   /* Reset motion in dest and copy it */
   dest->num_motion=0;
   dest->motion=NULL;

   if (src->num_motion>0) {
      for (j=0;j<src->num_motion;j++)
         flam3_add_motion_element(dest);

      memcpy(dest->motion,src->motion,src->num_motion*sizeof(flam3_xform));
   }
}

/* Copy one control point to another */
void flam3_copy(flam3_genome *dest, flam3_genome *src) {
   
   int i,ii;
   int numstd;

   /* If there are any xforms in dest before the copy, clean them up */
   clear_cp(dest, 1);

   /* Copy main contents of genome */
   memcpy(dest, src, sizeof(flam3_genome));

   /* Only the pointer to the xform was copied, not the actual xforms. */
   /* We need to create new xform memory storage for this new cp       */
   /* This goes for chaos, too.                                        */
   dest->num_xforms = 0;
   dest->final_xform_index = -1;
   dest->xform = NULL;
   dest->chaos = NULL;

   /* Add the standard xforms first */
   numstd = src->num_xforms-(src->final_xform_index>=0);
   flam3_add_xforms(dest, numstd, 0, 0);
   for (i=0;i<numstd;i++)
      flam3_copy_xform(&dest->xform[i], &src->xform[i]);
      
   /* Add the final x if it's present */
   if (src->final_xform_index>=0) {
      i = src->final_xform_index;
      flam3_add_xforms(dest, 1, 0, 1);
      ii = dest->final_xform_index;
      flam3_copy_xform(&dest->xform[ii],&src->xform[i]);      
   }
   
   /* Also, only the pointer to the chaos array was copied.
    * We have to take care of that as well.                 */   
   for (i=0;i<numstd;i++)
      memcpy(dest->chaos[i],src->chaos[i], numstd * sizeof(double));
      
}

void flam3_copyx(flam3_genome *dest, flam3_genome *src, int dest_std_xforms, int dest_final_xform) {

   int i,numsrcstd;
   
   /* If there are any xforms in dest before the copy, clean them up */
   clear_cp(dest, 1);

   /* Copy main contents of genome */
   memcpy(dest, src, sizeof(flam3_genome));

   /* Only the pointer to the xform was copied, not the actual xforms. */
   /* We need to create new xform memory storage for this new cp       */
   /* This goes for chaos, too.                                        */
   dest->num_xforms = 0;
   dest->xform = NULL;
   dest->chaos = NULL;
   dest->final_xform_index = -1;

   /* Add the padded standard xform list */
   /* Set the pad to 1 for these */
   flam3_add_xforms(dest, dest_std_xforms, 1, 0);

   numsrcstd = src->num_xforms - (src->final_xform_index >= 0);

   for(i=0;i<numsrcstd;i++) {

      /* When we copy the old xform, the pad is set to 0 */
      flam3_copy_xform(&dest->xform[i],&src->xform[i]);

      /* Copy the initial chaos from the src - the rest are already 1 */
      memcpy(dest->chaos[i], src->chaos[i], numsrcstd*sizeof(double));
      
   }   
   
   /* Add the final xform if necessary */
   if (dest_final_xform > 0) {
      flam3_add_xforms(dest, dest_final_xform, 1, 1);

      if (src->final_xform_enable > 0) {
      
         i = src->final_xform_index;
         
         flam3_copy_xform(&dest->xform[dest->num_xforms-1],&src->xform[i]);
         
      } else {
         /* Interpolated-against final xforms need animate & color_speed set to 0.0 */
         dest->xform[dest->num_xforms-1].num_motion = 0;
         dest->xform[dest->num_xforms-1].motion=NULL;
         dest->xform[dest->num_xforms-1].animate=0.0;
         dest->xform[dest->num_xforms-1].color_speed=0.0;
      }

   } else {
      dest->final_xform_index = -1;
      dest->final_xform_enable = 0;
   }

}

void clear_cp(flam3_genome *cp, int default_flag) {
    cp->palette_index = flam3_palette_random;
    cp->center[0] = 0.0;
    cp->center[1] = 0.0;
    cp->rot_center[0] = 0.0;
    cp->rot_center[1] = 0.0;
    cp->gamma = 4.0;
    cp->vibrancy = 1.0;
    cp->contrast = 1.0;
    cp->brightness = 4.0;
    cp->symmetry = 0;
    cp->hue_rotation = 0.0;
    cp->rotate = 0.0;
    cp->pixels_per_unit = 50;
    cp->interpolation = flam3_interpolation_linear;
    cp->palette_interpolation = flam3_palette_interpolation_hsv;

    cp->genome_index = 0;
    memset(cp->parent_fname,0,flam3_parent_fn_len);

    if (default_flag==flam3_defaults_on) {
       /* If defaults are on, set to reasonable values */
       cp->highlight_power = -1.0;
       cp->background[0] = 0.0;
       cp->background[1] = 0.0;
       cp->background[2] = 0.0;
       cp->width = 100;
       cp->height = 100;
       cp->spatial_oversample = 1;
       cp->spatial_filter_radius = 0.5;
       cp->zoom = 0.0;
       cp->sample_density = 1;
       /* Density estimation stuff defaulting to ON */
       cp->estimator = 9.0;
       cp->estimator_minimum = 0.0;
       cp->estimator_curve = 0.4;
       cp->gam_lin_thresh = 0.01;
//       cp->motion_exp = 0.0;
       cp->nbatches = 1;
       cp->ntemporal_samples = 1000;
       cp->spatial_filter_select = flam3_gaussian_kernel;
       cp->interpolation_type = flam3_inttype_log;
       cp->temporal_filter_type = flam3_temporal_box;
       cp->temporal_filter_width = 1.0;
       cp->temporal_filter_exp = 0.0;
       cp->palette_mode = flam3_palette_mode_step;

    } else {
       /* Defaults are off, so set to UN-reasonable values. */
       cp->highlight_power = -1.0;
       cp->background[0] = -1.0;
       cp->background[1] = -1.0;
       cp->background[2] = -1.0;
       cp->zoom = 999999999;
       cp->spatial_oversample = -1;
       cp->spatial_filter_radius = -1;
       cp->nbatches = -1;
       cp->ntemporal_samples = -1;
       cp->width = -1;
       cp->height = -1;
       cp->sample_density = -1;
       cp->estimator = -1;
       cp->estimator_minimum = -1;
       cp->estimator_curve = -1;
       cp->gam_lin_thresh = -1;
//       cp->motion_exp = -999;
       cp->nbatches = 0;
       cp->ntemporal_samples = 0;
       cp->spatial_filter_select = -1;
       cp->interpolation_type = -1;
       cp->temporal_filter_type = -1;
       cp->temporal_filter_width = -1;
       cp->temporal_filter_exp = -999;
       cp->palette_mode = -1;
    }

    if (cp->xform != NULL && cp->num_xforms > 0) {
        int i;
        int ns = cp->num_xforms - (cp->final_xform_index>=0);
        
       for (i=0;i<ns;i++) {
          free(cp->chaos[i]);
       }
       free(cp->chaos);
       cp->chaos=NULL;

      for (i=0;i<cp->num_xforms;i++)
         flam3_delete_motion_elements(&cp->xform[i]);
                       
       free(cp->xform);
       cp->xform=NULL;
       
       cp->num_xforms = 0;
    }
    
    cp->final_xform_enable = 0;
    cp->final_xform_index = -1;

}


int flam3_count_nthreads(void) {
   int nthreads;
#ifndef HAVE_LIBPTHREAD
   return(1);
#endif

#ifdef _WIN32
   SYSTEM_INFO sysInfo;
   GetSystemInfo(&sysInfo);
   nthreads = sysInfo.dwNumberOfProcessors;
#else
#ifdef __APPLE__
   kern_return_t    kr;
   host_name_port_t   host;
   unsigned int     size;
   struct host_basic_info   hi;

   host = mach_host_self();
   size = sizeof(hi)/sizeof(int);
   kr = host_info(host, HOST_BASIC_INFO, (host_info_t)&hi, &size);
   if (kr != KERN_SUCCESS) {
       mach_error("host_info():", kr);
       /* set threads to 1 on error */
       nthreads = 1;
   } else
       nthreads = hi.avail_cpus;
#else 
#ifndef _SC_NPROCESSORS_ONLN
   char line[MAXBUF];
   FILE *f = fopen("/proc/cpuinfo", "r");
   if (NULL == f) goto def;
   nthreads = 0;
   while (fgets(line, MAXBUF, f)) {
      if (!strncmp("processor\t:", line, 11))
         nthreads++;
   }
   fclose(f);
   if (nthreads < 1) goto def;
   return (nthreads);
def:
   fprintf(stderr, "could not read /proc/cpuinfo, using one render thread.\n");
   nthreads = 1;
#else
   nthreads = sysconf(_SC_NPROCESSORS_ONLN);
   if (nthreads < 1) nthreads = 1;
#endif
#endif
#endif
   return (nthreads);
}

flam3_genome *flam3_parse_xml2(char *xmldata, char *xmlfilename, int default_flag, int *ncps) {

   xmlDocPtr doc; /* Parsed XML document tree */
   xmlNode *rootnode;
   char *bn;
   int i;
   int loc_all_ncps=0;
   flam3_genome *loc_all_cp=NULL;
   char* locale = NULL;
   char* lorig  = setlocale(LC_NUMERIC, NULL);

   /* Parse XML string into internal document */
   /* Forbid network access during read       */
   doc = xmlReadMemory(xmldata, (int)strlen(xmldata), xmlfilename, NULL, XML_PARSE_NONET);

   /* Check for errors */
   if (doc==NULL) {
      fprintf(stderr, "Failed to parse %s\n", xmlfilename);
      return NULL;
   }

   /* What is the root node of the document? */
   rootnode = xmlDocGetRootElement(doc);
   
   // force use of "C" locale when writing reals.
   // first save away the current settings.
   if (lorig == NULL)
      fprintf(stderr, "error: couldn't get current locale\n");
   else {
      int slen = strlen(lorig) + 1;
      locale = (char*)malloc(slen);
      if (locale != NULL)
         memcpy(locale, lorig, slen);
   }
   if (setlocale(LC_NUMERIC, "C") == NULL)
      fprintf(stderr, "error: couldn't set C locale\n");

   /* Scan for <flame> nodes, starting with this node */
   bn = basename(xmlfilename);

   /* Have to use &loc_all_cp since the memory gets allocated in scan_for_flame_nodes */
   scan_for_flame_nodes(rootnode, bn, default_flag,&loc_all_cp,&loc_all_ncps);

   // restore locale
   if (locale != NULL) {
      if (setlocale(LC_NUMERIC, locale) == NULL)
         fprintf(stderr, "error: couldn't replace locale settings\n");
      free(locale);
   }
   
   xmlFreeDoc(doc);

   *ncps = loc_all_ncps;
   
   /* Check to see if the first control point or the second-to-last */
   /* control point has interpolation="smooth".  This is invalid    */
   /* and should be reset to linear (with a warning).               */
   if (loc_all_ncps>=1) {
      if (loc_all_cp[0].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for first segment.\n"
                        "         switching to linear.\n");
         loc_all_cp[0].interpolation = flam3_interpolation_linear;
      }
   }
   
   if (loc_all_ncps>=2) {
      if (loc_all_cp[(loc_all_ncps)-2].interpolation == flam3_interpolation_smooth) {
         fprintf(stderr,"Warning: smooth interpolation cannot be used for last segment.\n"
                        "         switching to linear.\n");
         loc_all_cp[loc_all_ncps-2].interpolation = flam3_interpolation_linear;
      }
   }
   
   /* Finally, ensure that consecutive 'rotate' parameters never exceed */
   /* a difference of more than 180 degrees (+/-) for interpolation.    */
   /* An adjustment of +/- 360 degrees is made until this is true.      */
   if (*ncps>1) {
   
      for (i=1;i<*ncps;i++) {

         /* Only do this adjustment if we're not in compat mode */
         if (flam3_inttype_compat != loc_all_cp[i-1].interpolation_type
        && flam3_inttype_older != loc_all_cp[i-1].interpolation_type) {
      
            while (loc_all_cp[i].rotate < loc_all_cp[i-1].rotate-180)
               loc_all_cp[i].rotate += 360;
            
            while (loc_all_cp[i].rotate > loc_all_cp[i-1].rotate+180)
               loc_all_cp[i].rotate -= 360;
         }
      }
   }
   
   //Note that concurrent calls to flam3, if in parallel, potentially segfault
   //if this function is called.  technically it's required but it doesn't
   //leak memory continuously.
   //xmlCleanupParser();
   
   return loc_all_cp;
}

flam3_genome * flam3_parse_from_file(FILE *f, char *fname, int default_flag, int *ncps) {
   int i, c, slen = 5000;
   char *s, *snew;
   flam3_genome *ret;

   /* Incrementally read XML file into a string */
   s = malloc(slen);
   i = 0;
   do {
      c = getc(f);
      if (EOF == c)
         break;
      s[i++] = c;
      if (i == slen-1) {
         slen *= 2;
         snew = realloc(s, slen);
         if (snew==NULL) {
            fprintf(stderr,"XML file too large to be read. continuing with partial file.\n");
            break;
         } else
            s = snew;
      }
   } while (1);

   /* Null-terminate the read XML data */
   s[i] = 0;

   /* Parse the XML string */
   if (fname)
      ret = flam3_parse_xml2(s, fname, default_flag, ncps);
   else
      ret = flam3_parse_xml2(s, "stdin", default_flag, ncps);
      
   free(s);
   
   return(ret);

}


void flam3_apply_template(flam3_genome *cp, flam3_genome *templ) {

   /* Check for invalid values - only replace those with valid ones */
   if (templ->background[0] >= 0)
      cp->background[0] = templ->background[0];
   if (templ->background[1] >= 0)
      cp->background[1] = templ->background[1];
   if (templ->background[1] >= 0)
      cp->background[2] = templ->background[2];
   if (templ->zoom < 999999998)
      cp->zoom = templ->zoom;
   if (templ->spatial_oversample > 0)
      cp->spatial_oversample = templ->spatial_oversample;
   if (templ->spatial_filter_radius >= 0)
      cp->spatial_filter_radius = templ->spatial_filter_radius;
   if (templ->sample_density > 0)
      cp->sample_density = templ->sample_density;
   if (templ->nbatches > 0)
      cp->nbatches = templ->nbatches;
   if (templ->ntemporal_samples > 0)
      cp->ntemporal_samples = templ->ntemporal_samples;
   if (templ->width > 0) {
      /* preserving scale should be an option */
      cp->pixels_per_unit = cp->pixels_per_unit * templ->width / cp->width;
      cp->width = templ->width;
   }
   if (templ->height > 0)
      cp->height = templ->height;
   if (templ->estimator >= 0)
      cp->estimator = templ->estimator;
   if (templ->estimator_minimum >= 0)
      cp->estimator_minimum = templ->estimator_minimum;
   if (templ->estimator_curve >= 0)
      cp->estimator_curve = templ->estimator_curve;
   if (templ->gam_lin_thresh >= 0)
      cp->gam_lin_thresh = templ->gam_lin_thresh;
   if (templ->nbatches>0)
      cp->nbatches = templ->nbatches;
   if (templ->ntemporal_samples>0)
      cp->ntemporal_samples = templ->ntemporal_samples;
   if (templ->spatial_filter_select>0)
      cp->spatial_filter_select = templ->spatial_filter_select;
   if (templ->interpolation >= 0)
      cp->interpolation = templ->interpolation;
   if (templ->interpolation_type >= 0)
      cp->interpolation_type = templ->interpolation_type;
   if (templ->temporal_filter_type >= 0)
      cp->temporal_filter_type = templ->temporal_filter_type;
   if (templ->temporal_filter_width > 0)
      cp->temporal_filter_width = templ->temporal_filter_width;
   if (templ->temporal_filter_exp > -900)
      cp->temporal_filter_exp = templ->temporal_filter_exp;
   if (templ->highlight_power >=0)
      cp->highlight_power = templ->highlight_power;
   if (templ->palette_mode >= 0)
      cp->palette_mode = templ->palette_mode;

}

char *flam3_print_to_string(flam3_genome *cp) {

   FILE *tmpflame;
   long stringbytes;
   char *genome_string;
   
   int using_tmpdir = 0;
   char *tmp_path;
   char tmpnam[256];
   
   tmpflame = tmpfile();
   if (NULL==tmpflame) {
#ifdef _WIN32       
       // This might be a permissions problem, so let's try to open a
       // tempfile in the env var TEMP's area instead
       tmp_path = getenv("TEMP");

       if (tmp_path != NULL) {
          strcpy(tmpnam, tmp_path);
          strcat(tmpnam, "\\fr0st.tmp");
          tmpflame = fopen(tmpnam, "w+");
          if (tmpflame != NULL) {
             using_tmpdir = 1;
          }
       }
#endif
       if (using_tmpdir == 0) {
          perror("FLAM3: opening temporary file");
          return (NULL);
       }
   }
   flam3_print(tmpflame,cp,NULL,flam3_dont_print_edits);
   stringbytes = ftell(tmpflame);
   fseek(tmpflame,0L, SEEK_SET);
   genome_string = (char *)calloc(stringbytes+1,1);
   if (stringbytes != fread(genome_string, 1, stringbytes, tmpflame)) {
       perror("FLAM3: reading string from temp file");
   }
   fclose(tmpflame);

   if (using_tmpdir)
      unlink(tmpnam);
   
   return(genome_string);
}
   

void flam3_print(FILE *f, flam3_genome *cp, char *extra_attributes, int print_edits) {
   int i,numstd;
   int flam27_flag;
   char *ai;

   // force use of "C" locale when writing reals.
   // first save away the current settings.
   char* locale = NULL;
   char* lorig  = setlocale(LC_NUMERIC, NULL);
   if (lorig == NULL)
      fprintf(stderr, "error: couldn't get current locale\n");
   else {
      int slen = strlen(lorig) + 1;
      locale = (char*)malloc(slen);
      if (locale != NULL)
      memcpy(locale, lorig, slen);
   }
   if (setlocale(LC_NUMERIC, "C") == NULL)
      fprintf(stderr, "error: couldn't set C locale\n");

   
   flam27_flag = argi("flam27",0);

   fprintf(f, "<flame version=\"FLAM3-%s\" time=\"%g\"", flam3_version(),cp->time);
   
   if (cp->flame_name[0]!=0)
      fprintf(f, " name=\"%s\"",cp->flame_name);
   
   fprintf(f, " size=\"%d %d\"", cp->width, cp->height);
   fprintf(f, " center=\"%g %g\"", cp->center[0], cp->center[1]);
   fprintf(f, " scale=\"%g\"", cp->pixels_per_unit);

   if (cp->zoom != 0.0)
      fprintf(f, " zoom=\"%g\"", cp->zoom);

   fprintf(f, " rotate=\"%g\"", cp->rotate);
   fprintf(f, " supersample=\"%d\"", cp->spatial_oversample);
   fprintf(f, " filter=\"%g\"", cp->spatial_filter_radius);

   /* Need to print the correct kernel to use */
   if (cp->spatial_filter_select == flam3_gaussian_kernel)
      fprintf(f, " filter_shape=\"gaussian\"");
   else if (cp->spatial_filter_select == flam3_hermite_kernel)
      fprintf(f, " filter_shape=\"hermite\"");
   else if (cp->spatial_filter_select == flam3_box_kernel)
      fprintf(f, " filter_shape=\"box\"");
   else if (cp->spatial_filter_select == flam3_triangle_kernel)
      fprintf(f, " filter_shape=\"triangle\"");
   else if (cp->spatial_filter_select == flam3_bell_kernel)
      fprintf(f, " filter_shape=\"bell\"");
   else if (cp->spatial_filter_select == flam3_b_spline_kernel)
      fprintf(f, " filter_shape=\"bspline\"");
   else if (cp->spatial_filter_select == flam3_mitchell_kernel)
      fprintf(f, " filter_shape=\"mitchell\"");
   else if (cp->spatial_filter_select == flam3_blackman_kernel)
      fprintf(f, " filter_shape=\"blackman\"");
   else if (cp->spatial_filter_select == flam3_catrom_kernel)
      fprintf(f, " filter_shape=\"catrom\"");
   else if (cp->spatial_filter_select == flam3_hanning_kernel)
      fprintf(f, " filter_shape=\"hanning\"");
   else if (cp->spatial_filter_select == flam3_hamming_kernel)
      fprintf(f, " filter_shape=\"hamming\"");
   else if (cp->spatial_filter_select == flam3_lanczos3_kernel)
      fprintf(f, " filter_shape=\"lanczos3\"");
   else if (cp->spatial_filter_select == flam3_lanczos2_kernel)
      fprintf(f, " filter_shape=\"lanczos2\"");
   else if (cp->spatial_filter_select == flam3_quadratic_kernel)
      fprintf(f, " filter_shape=\"quadratic\"");

   if (cp->temporal_filter_type == flam3_temporal_box)
      fprintf(f, " temporal_filter_type=\"box\"");
   else if (cp->temporal_filter_type == flam3_temporal_gaussian)
      fprintf(f, " temporal_filter_type=\"gaussian\"");
   else if (cp->temporal_filter_type == flam3_temporal_exp)
      fprintf(f, " temporal_filter_type=\"exp\" temporal_filter_exp=\"%g\"",cp->temporal_filter_exp);

   fprintf(f, " temporal_filter_width=\"%g\"",cp->temporal_filter_width);
   


   fprintf(f, " quality=\"%g\"", cp->sample_density);
   fprintf(f, " passes=\"%d\"", cp->nbatches);
   fprintf(f, " temporal_samples=\"%d\"", cp->ntemporal_samples);
   fprintf(f, " background=\"%g %g %g\"",
      cp->background[0], cp->background[1], cp->background[2]);
   fprintf(f, " brightness=\"%g\"", cp->brightness);
   fprintf(f, " gamma=\"%g\"", cp->gamma);
   
   if (!flam27_flag)
      fprintf(f, " highlight_power=\"%g\"", cp->highlight_power);
      
   fprintf(f, " vibrancy=\"%g\"", cp->vibrancy);
   fprintf(f, " estimator_radius=\"%g\" estimator_minimum=\"%g\" estimator_curve=\"%g\"",
      cp->estimator, cp->estimator_minimum, cp->estimator_curve);
   fprintf(f, " gamma_threshold=\"%g\"", cp->gam_lin_thresh);
   
   if (flam3_palette_mode_step == cp->palette_mode)
      fprintf(f, " palette_mode=\"step\"");
   else if (flam3_palette_mode_linear == cp->palette_mode)
      fprintf(f, " palette_mode=\"linear\"");

   if (flam3_interpolation_linear != cp->interpolation)
       fprintf(f, " interpolation=\"smooth\"");
       
   if (flam3_inttype_linear == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"linear\"");
   else if (flam3_inttype_log == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"log\"");
   else if (flam3_inttype_compat == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"old\"");
   else if (flam3_inttype_older == cp->interpolation_type)
       fprintf(f, " interpolation_type=\"older\"");


   if (flam3_palette_interpolation_hsv != cp->palette_interpolation)
       fprintf(f, " palette_interpolation=\"sweep\"");

   if (extra_attributes)
      fprintf(f, " %s", extra_attributes);

   fprintf(f, ">\n");

   if (cp->symmetry)
      fprintf(f, "   <symmetry kind=\"%d\"/>\n", cp->symmetry);
   
   numstd = cp->num_xforms - (cp->final_xform_index>=0);
   
   for (i = 0; i < cp->num_xforms; i++) {
   
      if (i==cp->final_xform_index)   
         flam3_print_xform(f, &cp->xform[i], 1, numstd, NULL, 0);
      else 
         flam3_print_xform(f, &cp->xform[i], 0, numstd, cp->chaos[i], 0);

   }

   int hexpalette = argi("hexpalette",0);

   if (hexpalette) {

      fprintf(f,"   <palette count=\"256\" format=\"RGB\">");

      for (i=0; i < 256; i++) {

         int r, g, b;
         r = rint(cp->palette[i].color[0] * 255.0);
         g = rint(cp->palette[i].color[1] * 255.0);
         b = rint(cp->palette[i].color[2] * 255.0);

         if (i % 8 == 0) {
            fprintf(f,"\n");
            fprintf(f,"      ");
         }

         fprintf(f,"%2x%2x%2x",r,g,b);

      }

      fprintf(f,"\n");
      fprintf(f,"   </palette>\n");

   } else {
   for (i = 0; i < 256; i++) {
      double r, g, b, a;
      r = (cp->palette[i].color[0] * 255.0);
      g = (cp->palette[i].color[1] * 255.0);
      b = (cp->palette[i].color[2] * 255.0);
      a = (cp->palette[i].color[3] * 255.0);
      
      fprintf(f, "   ");
      
      if (flam27_flag || a==255.0) {
      
         if (flam27_flag && a!=255.0)
            fprintf(stderr,"alpha channel in palette cannot be stored in 2.7 compatibility mode; truncating\n");
      
         if (getenv("intpalette"))
            fprintf(f, "<color index=\"%d\" rgb=\"%d %d %d\"/>", i, (int)rint(r), (int)rint(g), (int)rint(b));
         else {
#ifdef USE_FLOAT_INDICES
            fprintf(f, "<color index=\"%.10g\" rgb=\"%.6g %.6g %.6g\"/>", cp->palette[i].index, r, g, b);
#else
            fprintf(f, "<color index=\"%d\" rgb=\"%.6g %.6g %.6g\"/>", i, r, g, b);
#endif
         }
      } else {
         if (getenv("intpalette"))
            fprintf(f, "   <color index=\"%d\" rgba=\"%d %d %d %d\"/>", i, (int)rint(r), (int)rint(g), (int)rint(b), (int)rint(a));
         else
            fprintf(f, "   <color index=\"%d\" rgba=\"%.6g %.6g %.6g %.6g\"/>", i, r, g, b, a);
      }
//      if (i%4 == 3)
         fprintf(f, "\n");
         
   }
   }

   if (cp->edits != NULL && print_edits==flam3_print_edits) {

      /* We need a custom script for printing these */
      /* and it needs to be recursive               */
      xmlNodePtr elem_node = xmlDocGetRootElement(cp->edits);
      flam3_edit_print(f,elem_node, 1, 1);
   }
   fprintf(f, "</flame>\n");

   if (locale != NULL) {
      if (setlocale(LC_NUMERIC, locale) == NULL)
         fprintf(stderr, "error: couldn't restore locale settings\n");
      free(locale);
   }
}

#define PRINTNON(p) do { if (x->p != 0.0) fprintf(f, #p "=\"%f\" ",x->p); } while(0)
   

void flam3_print_xform(FILE *f, flam3_xform *x, int final_flag, int numstd, double *chaos_row, int motion_flag) {

   int blob_var=0,pdj_var=0,fan2_var=0,rings2_var=0,perspective_var=0;
   int juliaN_var=0,juliaScope_var=0,radialBlur_var=0,pie_var=0,disc2_var=0;
   int ngon_var=0,curl_var=0,rectangles_var=0,supershape_var=0;
   int flower_var=0,conic_var=0,parabola_var=0,bent2_var=0,bipolar_var=0;
   int cell_var=0,cpow_var=0,curve_var=0,escher_var=0,lazys_var=0;
   int modulus_var=0,oscope_var=0,popcorn2_var=0,separation_var=0;
   int split_var=0,splits_var=0,stripes_var=0,wedge_var=0,wedgeJ_var=0;
   int wedgeS_var=0,whorl_var=0,waves2_var=0,auger_var=0,flux_var=0;
   int mobius_var=0;

   int j;
   int lnv;

   int flam27_flag;
   char *ai;
   
   flam27_flag = argi("flam27",0);

   /* Motion flag will not be set if flam27_flag is set */
   if (motion_flag) {
      fprintf(f, "      <motion motion_frequency=\"%d\" ",x->motion_freq);
      if (x->motion_func == MOTION_SIN)
         fprintf(f, "motion_function=\"sin\" ");
      else if (x->motion_func == MOTION_TRIANGLE)
         fprintf(f, "motion_function=\"triangle\" ");
      else if (x->motion_func == MOTION_HILL)
         fprintf(f, "motion_function=\"hill\" ");
   } else {
      if (final_flag)
         fprintf(f, "   <finalxform ");
      else
         fprintf(f, "   <xform weight=\"%g\" ", x->density);
   }
   
   if (!motion_flag || x->color != 0.0)
      fprintf(f, "color=\"%g\" ", x->color);
   
   if (flam27_flag)
      fprintf(f, "symmetry=\"%g\" ", 1.0-2.0*x->color_speed);
   else if (!motion_flag)   
      fprintf(f, "color_speed=\"%g\" ", x->color_speed);
   
   if (!final_flag && !motion_flag && !flam27_flag)
      fprintf(f, "animate=\"%g\" ", x->animate);

   lnv = flam27_flag ? 54:flam3_nvariations;

   for (j = 0; j < lnv; j++) {
      double v = x->var[j];
      if (0.0 != v) {
         fprintf(f, "%s=\"%g\" ", flam3_variation_names[j], v);
         if (j==VAR_BLOB)
            blob_var=1;
         else if (j==VAR_PDJ)
            pdj_var=1;
         else if (j==VAR_FAN2)
            fan2_var=1;
         else if (j==VAR_RINGS2)
            rings2_var=1;
         else if (j==VAR_PERSPECTIVE)
            perspective_var=1;
         else if (j==VAR_JULIAN)
            juliaN_var=1;
         else if (j==VAR_JULIASCOPE)
            juliaScope_var=1;
         else if (j==VAR_RADIAL_BLUR)
            radialBlur_var=1;
         else if (j==VAR_PIE)
            pie_var=1;
         else if (j==VAR_NGON)
            ngon_var=1;
         else if (j==VAR_CURL)
            curl_var=1;
         else if (j==VAR_RECTANGLES)
            rectangles_var=1;
         else if (j==VAR_DISC2)
            disc2_var=1;
         else if (j==VAR_SUPER_SHAPE)
            supershape_var=1;
         else if (j==VAR_FLOWER)
            flower_var=1;
         else if (j==VAR_CONIC)
            conic_var=1;
         else if (j==VAR_PARABOLA)
            parabola_var=1;
         else if (j==VAR_BENT2)
            bent2_var=1;
         else if (j==VAR_BIPOLAR)
            bipolar_var=1;
         else if (j==VAR_CELL)
            cell_var=1;
         else if (j==VAR_CPOW)
            cpow_var=1;
         else if (j==VAR_CURVE)
            curve_var=1;
         else if (j==VAR_ESCHER)
            escher_var=1;
         else if (j==VAR_LAZYSUSAN)
            lazys_var=1;
         else if (j==VAR_MODULUS)
            modulus_var=1;
         else if (j==VAR_OSCILLOSCOPE)
            oscope_var=1;
         else if (j==VAR_POPCORN2)
            popcorn2_var=1;
         else if (j==VAR_SPLIT)
            split_var=1;
         else if (j==VAR_SPLITS)
            splits_var=1;
         else if (j==VAR_STRIPES)
            stripes_var=1;
         else if (j==VAR_WEDGE)
            wedge_var=1;
         else if (j==VAR_WEDGE_JULIA)
            wedgeJ_var=1;
         else if (j==VAR_WEDGE_SPH)
            wedgeS_var=1;
         else if (j==VAR_WHORL)
            whorl_var=1;
         else if (j==VAR_WAVES2)
            waves2_var=1;
         else if (j==VAR_AUGER)
            auger_var=1;
         else if (j==VAR_FLUX)
            flux_var=1;
         else if (j==VAR_MOBIUS)
            mobius_var=1;
      }
   }

   if (!motion_flag) {
      if (blob_var==1) {
         fprintf(f, "blob_low=\"%g\" ", x->blob_low);
         fprintf(f, "blob_high=\"%g\" ", x->blob_high);
         fprintf(f, "blob_waves=\"%g\" ", x->blob_waves);
      }

      if (pdj_var==1) {
         fprintf(f, "pdj_a=\"%g\" ", x->pdj_a);
         fprintf(f, "pdj_b=\"%g\" ", x->pdj_b);
         fprintf(f, "pdj_c=\"%g\" ", x->pdj_c);
         fprintf(f, "pdj_d=\"%g\" ", x->pdj_d);
      }

      if (fan2_var==1) {
         fprintf(f, "fan2_x=\"%g\" ", x->fan2_x);
         fprintf(f, "fan2_y=\"%g\" ", x->fan2_y);
      }

      if (rings2_var==1) {
         fprintf(f, "rings2_val=\"%g\" ", x->rings2_val);
      }

      if (perspective_var==1) {
         fprintf(f, "perspective_angle=\"%g\" ", x->perspective_angle);
         fprintf(f, "perspective_dist=\"%g\" ", x->perspective_dist);
      }

      if (juliaN_var==1) {
         fprintf(f, "julian_power=\"%g\" ", x->julian_power);
         fprintf(f, "julian_dist=\"%g\" ", x->julian_dist);
      }

      if (juliaScope_var==1) {
         fprintf(f, "juliascope_power=\"%g\" ", x->juliascope_power);
         fprintf(f, "juliascope_dist=\"%g\" ", x->juliascope_dist);
      }

      if (radialBlur_var==1) {
         fprintf(f, "radial_blur_angle=\"%g\" ", x->radial_blur_angle);
      }

      if (pie_var==1) {
         fprintf(f, "pie_slices=\"%g\" ", x->pie_slices);
         fprintf(f, "pie_rotation=\"%g\" ", x->pie_rotation);
         fprintf(f, "pie_thickness=\"%g\" ", x->pie_thickness);
      }

      if (ngon_var==1) {
         fprintf(f, "ngon_sides=\"%g\" ", x->ngon_sides);
         fprintf(f, "ngon_power=\"%g\" ", x->ngon_power);
         fprintf(f, "ngon_corners=\"%g\" ", x->ngon_corners);
         fprintf(f, "ngon_circle=\"%g\" ", x->ngon_circle);
      }

      if (curl_var==1) {
         fprintf(f, "curl_c1=\"%g\" ", x->curl_c1);
         fprintf(f, "curl_c2=\"%g\" ", x->curl_c2);
      }

      if (rectangles_var==1) {
         fprintf(f, "rectangles_x=\"%g\" ", x->rectangles_x);
         fprintf(f, "rectangles_y=\"%g\" ", x->rectangles_y);
      }

      if (disc2_var==1) {
         fprintf(f, "disc2_rot=\"%g\" ", x->disc2_rot);
         fprintf(f, "disc2_twist=\"%g\" ", x->disc2_twist);
      }

      if (supershape_var==1) {
         fprintf(f, "super_shape_rnd=\"%g\" ", x->super_shape_rnd);
         fprintf(f, "super_shape_m=\"%g\" ", x->super_shape_m);
         fprintf(f, "super_shape_n1=\"%g\" ", x->super_shape_n1);
         fprintf(f, "super_shape_n2=\"%g\" ", x->super_shape_n2);
         fprintf(f, "super_shape_n3=\"%g\" ", x->super_shape_n3);
         fprintf(f, "super_shape_holes=\"%g\" ", x->super_shape_holes);
      }

      if (flower_var==1) {
         fprintf(f, "flower_petals=\"%g\" ", x->flower_petals);
         fprintf(f, "flower_holes=\"%g\" ", x->flower_holes);
      }

      if (conic_var==1) {
         fprintf(f, "conic_eccentricity=\"%g\" ", x->conic_eccentricity);
         fprintf(f, "conic_holes=\"%g\" ", x->conic_holes);
      }

      if (parabola_var==1) {
         fprintf(f, "parabola_height=\"%g\" ", x->parabola_height);
         fprintf(f, "parabola_width=\"%g\" ", x->parabola_width);
      }
      
      if (bent2_var==1) {
         fprintf(f, "bent2_x=\"%g\" ", x->bent2_x);
         fprintf(f, "bent2_y=\"%g\" ", x->bent2_y);
      }
      
      if (bipolar_var==1) {
         fprintf(f, "bipolar_shift=\"%g\" ", x->bipolar_shift);
      }

      if (cell_var==1) {
         fprintf(f, "cell_size=\"%g\" ", x->cell_size);
      }

      if (cpow_var==1) {
         fprintf(f, "cpow_i=\"%g\" ", x->cpow_i);
         fprintf(f, "cpow_r=\"%g\" ", x->cpow_r);
         fprintf(f, "cpow_power=\"%g\" ", x->cpow_power);
      }

      if (curve_var==1) {
         fprintf(f, "curve_xamp=\"%g\" ", x->curve_xamp);
         fprintf(f, "curve_yamp=\"%g\" ", x->curve_yamp);
         fprintf(f, "curve_xlength=\"%g\" ", x->curve_xlength);
         fprintf(f, "curve_ylength=\"%g\" ", x->curve_ylength);
      }

      if (escher_var==1) {
         fprintf(f, "escher_beta=\"%g\" ", x->escher_beta);
      }

      if (lazys_var==1) {
         fprintf(f, "lazysusan_x=\"%g\" ", x->lazysusan_x);
         fprintf(f, "lazysusan_y=\"%g\" ", x->lazysusan_y);
         fprintf(f, "lazysusan_spin=\"%g\" ", x->lazysusan_spin);
         fprintf(f, "lazysusan_space=\"%g\" ", x->lazysusan_space);
         fprintf(f, "lazysusan_twist=\"%g\" ", x->lazysusan_twist);
      }

      if (modulus_var==1) {
         fprintf(f, "modulus_x=\"%g\" ", x->modulus_x);
         fprintf(f, "modulus_y=\"%g\" ", x->modulus_y);
      }

      if (oscope_var==1) {
         fprintf(f, "oscilloscope_separation=\"%g\" ", x->oscope_separation);
         fprintf(f, "oscilloscope_frequency=\"%g\" ", x->oscope_frequency);
         fprintf(f, "oscilloscope_amplitude=\"%g\" ", x->oscope_amplitude);
         fprintf(f, "oscilloscope_damping=\"%g\" ", x->oscope_damping);
      }

      if (popcorn2_var==1) {
         fprintf(f, "popcorn2_x=\"%g\" ", x->popcorn2_x);
         fprintf(f, "popcorn2_y=\"%g\" ", x->popcorn2_y);
         fprintf(f, "popcorn2_c=\"%g\" ", x->popcorn2_c);
      }

      if (separation_var==1) {
         fprintf(f, "separation_x=\"%g\" ", x->separation_x);
         fprintf(f, "separation_y=\"%g\" ", x->separation_y);
         fprintf(f, "separation_xinside=\"%g\" ", x->separation_xinside);
         fprintf(f, "separation_yinside=\"%g\" ", x->separation_yinside);
      }

      if (split_var==1) {
         fprintf(f, "split_xsize=\"%g\" ", x->split_xsize);
         fprintf(f, "split_ysize=\"%g\" ", x->split_ysize);
      }

      if (splits_var==1) {
         fprintf(f, "splits_x=\"%g\" ", x->splits_x);
         fprintf(f, "splits_y=\"%g\" ", x->splits_y);
      }

      if (stripes_var==1) {
         fprintf(f, "stripes_space=\"%g\" ", x->stripes_space);
         fprintf(f, "stripes_warp=\"%g\" ", x->stripes_warp);
      }

      if (wedge_var==1) {
         fprintf(f, "wedge_angle=\"%g\" ", x->wedge_angle);
         fprintf(f, "wedge_hole=\"%g\" ", x->wedge_hole);
         fprintf(f, "wedge_count=\"%g\" ", x->wedge_count);
         fprintf(f, "wedge_swirl=\"%g\" ", x->wedge_swirl);            
      }

      if (wedgeJ_var==1) {
         fprintf(f, "wedge_julia_angle=\"%g\" ", x->wedge_julia_angle);
         fprintf(f, "wedge_julia_count=\"%g\" ", x->wedge_julia_count);
         fprintf(f, "wedge_julia_power=\"%g\" ", x->wedge_julia_power);
         fprintf(f, "wedge_julia_dist=\"%g\" ", x->wedge_julia_dist);            
      }

      if (wedgeS_var==1) {
         fprintf(f, "wedge_sph_angle=\"%g\" ", x->wedge_sph_angle);
         fprintf(f, "wedge_sph_hole=\"%g\" ", x->wedge_sph_hole);
         fprintf(f, "wedge_sph_count=\"%g\" ", x->wedge_sph_count);
         fprintf(f, "wedge_sph_swirl=\"%g\" ", x->wedge_sph_swirl);            
      }

      if (whorl_var==1) {
         fprintf(f, "whorl_inside=\"%g\" ", x->whorl_inside);
         fprintf(f, "whorl_outside=\"%g\" ", x->whorl_outside);
      }

      if (waves2_var==1) {
         fprintf(f, "waves2_scalex=\"%g\" ", x->waves2_scalex);
         fprintf(f, "waves2_scaley=\"%g\" ", x->waves2_scaley);
         fprintf(f, "waves2_freqx=\"%g\" ", x->waves2_freqx);
         fprintf(f, "waves2_freqy=\"%g\" ", x->waves2_freqy);
      }

      if (auger_var==1) {
         fprintf(f, "auger_sym=\"%g\" ", x->auger_sym);
         fprintf(f, "auger_weight=\"%g\" ", x->auger_weight);
         fprintf(f, "auger_freq=\"%g\" ", x->auger_freq);
         fprintf(f, "auger_scale=\"%g\" ", x->auger_scale);
      }

      if (flux_var==1)
         fprintf(f, "flux_spread=\"%g\" ", x->flux_spread);

      if (mobius_var==1) {
         fprintf(f, "mobius_re_a=\"%g\" ", x->mobius_re_a);
         fprintf(f, "mobius_im_a=\"%g\" ", x->mobius_im_a);
         fprintf(f, "mobius_re_b=\"%g\" ", x->mobius_re_b);
         fprintf(f, "mobius_im_b=\"%g\" ", x->mobius_im_b);
         fprintf(f, "mobius_re_c=\"%g\" ", x->mobius_re_c);
         fprintf(f, "mobius_im_c=\"%g\" ", x->mobius_im_c);
         fprintf(f, "mobius_re_d=\"%g\" ", x->mobius_re_d);
         fprintf(f, "mobius_im_d=\"%g\" ", x->mobius_im_d);
      }

      fprintf(f, "coefs=\"");
      for (j = 0; j < 3; j++) {
         if (j) fprintf(f, " ");
         fprintf(f, "%g %g", x->c[j][0], x->c[j][1]);
      }
      fprintf(f, "\"");

      if (!id_matrix(x->post)) {
         fprintf(f, " post=\"");
         for (j = 0; j < 3; j++) {
            if (j) fprintf(f, " ");
            fprintf(f, "%g %g", x->post[j][0], x->post[j][1]);
         }
         fprintf(f, "\"");
      }


   } else {
      /* For motion, print any parameter if it's nonzero */
      PRINTNON(blob_low);
      PRINTNON(blob_high);
      PRINTNON(blob_waves);

      PRINTNON(pdj_a);
      PRINTNON(pdj_b);
      PRINTNON(pdj_c);
      PRINTNON(pdj_d);

      PRINTNON(fan2_x);
      PRINTNON(fan2_y);

      PRINTNON(rings2_val);

      PRINTNON(perspective_angle);
      PRINTNON(perspective_dist);
      
      PRINTNON(julian_power);
      PRINTNON(julian_dist);      

      PRINTNON(juliascope_power);
      PRINTNON(juliascope_dist);      

      PRINTNON(radial_blur_angle);

      PRINTNON(pie_slices);
      PRINTNON(pie_rotation);
      PRINTNON(pie_thickness);      

      PRINTNON(ngon_sides);
      PRINTNON(ngon_power);
      PRINTNON(ngon_corners);
      PRINTNON(ngon_circle);

      PRINTNON(curl_c1);
      PRINTNON(curl_c2);

      PRINTNON(rectangles_x);
      PRINTNON(rectangles_y);

      PRINTNON(disc2_rot);
      PRINTNON(disc2_twist);

      PRINTNON(super_shape_rnd);
      PRINTNON(super_shape_m);
      PRINTNON(super_shape_n1);
      PRINTNON(super_shape_n2);
      PRINTNON(super_shape_n3);
      PRINTNON(super_shape_holes);

      PRINTNON(flower_petals);
      PRINTNON(flower_holes);

      PRINTNON(conic_eccentricity);
      PRINTNON(conic_holes);

      PRINTNON(parabola_height);
      PRINTNON(parabola_width);

      PRINTNON(bent2_x);
      PRINTNON(bent2_y);

      PRINTNON(bipolar_shift);      

      PRINTNON(cell_size);

      PRINTNON(cpow_i);
      PRINTNON(cpow_r);
      PRINTNON(cpow_power);

      PRINTNON(curve_xamp);
      PRINTNON(curve_yamp);
      PRINTNON(curve_xlength);
      PRINTNON(curve_ylength);

      PRINTNON(escher_beta);

      PRINTNON(lazysusan_x);
      PRINTNON(lazysusan_y);
      PRINTNON(lazysusan_spin);
      PRINTNON(lazysusan_space);
      PRINTNON(lazysusan_twist);

      PRINTNON(modulus_x);
      PRINTNON(modulus_y);

      PRINTNON(oscope_separation);
      PRINTNON(oscope_frequency);
      PRINTNON(oscope_amplitude);
      PRINTNON(oscope_damping);

      PRINTNON(popcorn2_x);
      PRINTNON(popcorn2_y);
      PRINTNON(popcorn2_c);

      PRINTNON(separation_x);
      PRINTNON(separation_y);
      PRINTNON(separation_xinside);
      PRINTNON(separation_yinside);

      PRINTNON(split_xsize);
      PRINTNON(split_ysize);

      PRINTNON(splits_x);
      PRINTNON(splits_y);

      PRINTNON(stripes_space);
      PRINTNON(stripes_warp);

      PRINTNON(wedge_angle);
      PRINTNON(wedge_hole);
      PRINTNON(wedge_count);
      PRINTNON(wedge_swirl);

      PRINTNON(wedge_julia_angle);
      PRINTNON(wedge_julia_count);
      PRINTNON(wedge_julia_power);
      PRINTNON(wedge_julia_dist);

      PRINTNON(wedge_sph_angle);
      PRINTNON(wedge_sph_hole);
      PRINTNON(wedge_sph_count);
      PRINTNON(wedge_sph_swirl);
      
      PRINTNON(whorl_inside);
      PRINTNON(whorl_outside);
      
      PRINTNON(waves2_scalex);
      PRINTNON(waves2_scaley);
      PRINTNON(waves2_freqx);
      PRINTNON(waves2_freqy);

      PRINTNON(auger_sym);
      PRINTNON(auger_weight);
      PRINTNON(auger_freq);
      PRINTNON(auger_scale);

      PRINTNON(flux_spread);

      PRINTNON(mobius_re_a);
      PRINTNON(mobius_im_a);
      PRINTNON(mobius_re_b);
      PRINTNON(mobius_im_b);
      PRINTNON(mobius_re_c);
      PRINTNON(mobius_im_c);
      PRINTNON(mobius_re_d);
      PRINTNON(mobius_im_d);
      
      if (!zero_matrix(x->c)) {
         fprintf(f, "coefs=\"");
         for (j = 0; j < 3; j++) {
            if (j) fprintf(f, " ");
            fprintf(f, "%g %g", x->c[j][0], x->c[j][1]);
         }
         fprintf(f, "\"");
      }

      if (!zero_matrix(x->post)) {
         fprintf(f, " post=\"");
         for (j = 0; j < 3; j++) {
            if (j) fprintf(f, " ");
            fprintf(f, "%g %g", x->post[j][0], x->post[j][1]);
         }
         fprintf(f, "\"");
      }
   }   
      
   if (!final_flag && !motion_flag && !flam27_flag) {
      
      /* Print out the chaos row for this xform */
      int numcols = numstd;

      while (numcols > 0 && chaos_row[numcols-1]==1.0)
         numcols--;
         
      if (numcols>0) {   
         fprintf(f, " chaos=\"");
         for (j=0;j<numcols;j++)
            fprintf(f, "%g ",chaos_row[j]);
         fprintf(f, "\"");
      }
      
               
   }

   if (!flam27_flag && !motion_flag) {
      fprintf(f, " opacity=\"%g\"",x->opacity);
   }

   if (!motion_flag && x->num_motion>0 && !flam27_flag) {
      int nm;
      
      fprintf(f,">\n");
      
      for (nm=0; nm<x->num_motion; nm++)
         flam3_print_xform(f, &(x->motion[nm]), 0, 0, NULL, 1);
         
      if (final_flag)
         fprintf(f,"   </finalxform>\n");
      else
         fprintf(f,"   </xform>\n");

   } else
      fprintf(f, "/>\n");
}


/* returns a uniform variable from 0 to 1 */
double flam3_random01() {
   return (random() & 0xfffffff) / (double) 0xfffffff;
}

double flam3_random11() {
   return ((random() & 0xfffffff) - 0x7ffffff) / (double) 0x7ffffff;
}

/* This function must be called prior to rendering a frame */
void flam3_init_frame(flam3_frame *f) {

   char *ai;
   char *isaac_seed = args("isaac_seed",NULL);
   long int default_isaac_seed = (long int)time(0);  

   /* Clear out the isaac state */
   memset(f->rc.randrsl, 0, RANDSIZ*sizeof(ub4));

   /* Set the isaac seed */
   if (NULL == isaac_seed) {
      int lp;
      /* No isaac seed specified.  Use the system time to initialize. */
      for (lp = 0; lp < RANDSIZ; lp++)
         f->rc.randrsl[lp] = default_isaac_seed;
   } else {
      /* Use the specified string */
      strncpy((char *)&f->rc.randrsl,(const char *)isaac_seed, RANDSIZ*sizeof(ub4));
   }

   /* Initialize the random number generator */
   irandinit(&f->rc,1);
}

/* returns uniform variable from ISAAC rng */
double flam3_random_isaac_01(randctx *ct) {
   return ((int)irand(ct) & 0xfffffff) / (double) 0xfffffff;
}

double flam3_random_isaac_11(randctx *ct) {
   return (((int)irand(ct) & 0xfffffff) - 0x7ffffff) / (double) 0x7ffffff;
}

int flam3_random_bit() {
  /* might not be threadsafe */
  static int n = 0;
  static int l;
  if (0 == n) {
    l = random();
    n = 20;
  } else {
    l = l >> 1;
    n--;
  }
  return l & 1;
}

int flam3_random_isaac_bit(randctx *ct) {
   int tmp = irand(ct);
   return tmp & 1;
}

static double round6(double x) {
  x *= 1e6;
  if (x < 0) x -= 1.0;
  return 1e-6*(int)(x+0.5);
}

/* sym=2 or more means rotational
   sym=1 means identity, ie no symmetry
   sym=0 means pick a random symmetry (maybe none)
   sym=-1 means bilateral (reflection)
   sym=-2 or less means rotational and reflective
*/
void flam3_add_symmetry(flam3_genome *cp, int sym) {
   int i, j, k;
   double a;
   int result = 0;

   if (0 == sym) {
      static int sym_distrib[] = {
         -4, -3,
         -2, -2, -2,
         -1, -1, -1,
         2, 2, 2,
         3, 3,
         4, 4,
      };
      if (random()&1) {
         sym = random_distrib(sym_distrib);
      } else if (random()&31) {
         sym = (random()%13)-6;
      } else {
         sym = (random()%51)-25;
      }
   }

   if (1 == sym || 0 == sym) return;

   cp->symmetry = sym;

   if (sym < 0) {

      i = cp->num_xforms;
      if (cp->final_xform_enable)
         i -= 1;

      flam3_add_xforms(cp,1,0,0);

      cp->xform[i].density = 1.0;
      cp->xform[i].color_speed = 0.0;
      cp->xform[i].animate = 0.0;
      cp->xform[i].var[0] = 1.0;
      for (j = 1; j < flam3_nvariations; j++)
         cp->xform[i].var[j] = 0;
      cp->xform[i].color = 1.0;
      cp->xform[i].c[0][0] = -1.0;
      cp->xform[i].c[0][1] = 0.0;
      cp->xform[i].c[1][0] = 0.0;
      cp->xform[i].c[1][1] = 1.0;
      cp->xform[i].c[2][0] = 0.0;
      cp->xform[i].c[2][1] = 0.0;

      result++;
      sym = -sym;
   }

   a = 2*M_PI/sym;

   for (k = 1; k < sym; k++) {

      i = cp->num_xforms;
      if (cp->final_xform_enable)
         i -= 1;

      flam3_add_xforms(cp, 1, 0,0);

      cp->xform[i].density = 1.0;
      cp->xform[i].color_speed = 0.0;
      cp->xform[i].animate = 0.0;
      cp->xform[i].var[0] = 1.0;
      for (j = 1; j < flam3_nvariations; j++)
         cp->xform[i].var[j] = 0;
      cp->xform[i].color = (sym<3) ? 0.0 : ((k-1.0)/(sym-2.0));
      cp->xform[i].c[0][0] = round6(cos(k*a));
      cp->xform[i].c[0][1] = round6(sin(k*a));
      cp->xform[i].c[1][0] = round6(-cp->xform[i].c[0][1]);
      cp->xform[i].c[1][1] = cp->xform[i].c[0][0];
      cp->xform[i].c[2][0] = 0.0;
      cp->xform[i].c[2][1] = 0.0;

      result++;
   }   
   
   qsort((char *) &cp->xform[cp->num_xforms-result], result,
      sizeof(flam3_xform), compare_xforms);
      
}

void add_to_action(char *action, char *addtoaction) {

   /* action is a flam3_max_action_length array */
   if (action != NULL) {

      int alen = strlen(action);
      int addlen = strlen(addtoaction);

      if (alen+addlen < flam3_max_action_length)
         strcat(action,addtoaction);
      else
         fprintf(stderr,"action string too long, truncating...\n");
   }
}


void flam3_cross(flam3_genome *cp0, flam3_genome *cp1, flam3_genome *out, int cross_mode, randctx *rc, char *action) {

   int i0,i1, i,j, rb;
   char ministr[10];   

   if (cross_mode == CROSS_NOT_SPECIFIED) {
   
      double s = flam3_random_isaac_01(rc);
      
      if (s < 0.1)
         cross_mode = CROSS_UNION;
      else if (s < 0.2)
         cross_mode = CROSS_INTERPOLATE;
      else
         cross_mode = CROSS_ALTERNATE;

   }
   
   if (cross_mode == CROSS_UNION) {
   
      flam3_xform mycopy;
   
      /* Make a copy of cp0 */
      flam3_copy(out, cp0);
      
      for (j=0;j<cp1->num_xforms;j++) {
         /* Skip over the final xform, if it's present.    */
         /* Default behavior keeps the final from parent0. */
         if (cp1->final_xform_index == j)		     
            continue;

         i = out->num_xforms;
         if (out->final_xform_enable)
            i -= 1;

         flam3_add_xforms(out, 1, 0, 0);
         flam3_copy_xform(&out->xform[i],&cp1->xform[j]);
      }
      
      /* Put the final xform last (if there is one) */
      /* We do not need to do complicated xform copies here since we're just moving them around */
      if (out->final_xform_index >= 0) {
         mycopy = out->xform[out->final_xform_index];
         out->xform[out->final_xform_index] = out->xform[out->num_xforms-1];
         out->xform[out->num_xforms-1] = mycopy;
         out->final_xform_index = out->num_xforms-1;
      }
      
      add_to_action(action,"cross union");
      
   } else if (cross_mode == CROSS_INTERPOLATE) {
   
      /* linearly interpolate somewhere between the two */
      flam3_genome parents[2];
      double t = flam3_random_isaac_01(rc);

      memset(parents, 0, 2*sizeof(flam3_genome));

      flam3_copy(&(parents[0]), cp0);
      flam3_copy(&(parents[1]), cp1);

      parents[0].time = 0.0;
      parents[1].time = 1.0;
      flam3_interpolate(parents, 2, t, 0, out);
      
      for (i=0;i<out->num_xforms;i++)
         flam3_delete_motion_elements(&out->xform[i]);
         
      clear_cp(&parents[0],flam3_defaults_on);
      clear_cp(&parents[1],flam3_defaults_on);
      
      sprintf(ministr,"%7.5g",t);
   
      add_to_action(action,"cross interpolate ");
      add_to_action(action,ministr);
      
   } else {
   
      /* alternate mode */
      int got0, got1, used_parent;
      char *trystr;

      trystr = calloc(4 * (cp0->num_xforms + cp1->num_xforms), sizeof(char));

      /* each xform comes from a random parent, possible for an entire parent to be excluded */
      do {

         trystr[0] = 0;
         got0 = got1 = 0;
         rb = flam3_random_isaac_bit(rc);
         sprintf(ministr,"%d:",rb);
         strcat(trystr,ministr);

         /* Copy the parent, sorting the final xform to the end if it's present. */
         if (rb)
            flam3_copyx(out, cp1, cp1->num_xforms - (cp1->final_xform_index > 0), cp1->final_xform_enable);
         else
            flam3_copyx(out, cp0, cp0->num_xforms - (cp0->final_xform_index > 0), cp0->final_xform_enable);

         used_parent = rb;

         /* Only replace non-final xforms */

         for (i = 0; i < out->num_xforms - out->final_xform_enable; i++) {
            rb = flam3_random_isaac_bit(rc);

            /* Replace xform if bit is 1 */
            if (rb==1) {
               if (used_parent==0) {
                  if (i < cp1->num_xforms && cp1->xform[i].density > 0) {
                     flam3_copy_xform(&out->xform[i],&cp1->xform[i]);
                     sprintf(ministr," 1");
                     got1 = 1;
                  } else {
                     sprintf(ministr," 0");
                     got0 = 1;
                  }
               } else {
                  if (i < cp0->num_xforms && cp0->xform[i].density > 0) {
                     flam3_copy_xform(&out->xform[i],&cp0->xform[i]);
                     sprintf(ministr," 0");
                     got0 = 1;
                  } else {
                     sprintf(ministr," 1");
                     got1 = 1;
                  }
               }
            } else {
               sprintf(ministr," %d",used_parent);
               if (used_parent)
                  got1 = 1;
               else
                  got0 = 1;
            }

            strcat(trystr,ministr);
         }
             
         if (used_parent==0 && cp0->final_xform_enable)
            got0 = 1;
         else if (used_parent==1 && cp1->final_xform_enable)
            got1 = 1;
               
      } while ((i > 1) && !(got0 && got1));

      add_to_action(action,"cross alternate ");
      add_to_action(action,trystr);
      
      free(trystr);
   }
   
   /* reset color coords */
   for (i = 0; i < out->num_xforms; i++) {
      out->xform[i].color = i&1;
   }
               
   /* Potentially genetically cross the two colormaps together */
   if (flam3_random_isaac_01(rc) < 0.4) {
                              
      /* Select the starting parent */
      int startParent=flam3_random_isaac_bit(rc);
      int ci;
                        
      add_to_action(action," cmap_cross");
      sprintf(ministr," %d:",startParent);
      add_to_action(action,ministr);
                  
      /* Loop over the entries, switching to the other parent 1% of the time */
      for (ci=0;ci<256;ci++) {
         if (flam3_random_isaac_01(rc)<.01) {
            startParent = 1-startParent;
            sprintf(ministr," %d",ci);
            add_to_action(action,ministr);
         }
                     
         out->palette[ci] = startParent ? cp1->palette[ci]: cp0->palette[ci];
      }
   }

}

void flam3_mutate(flam3_genome *cp, int mutate_mode, int *ivars, int ivars_n, int sym, double speed, randctx *rc, char *action) {

   double randselect;
   flam3_genome mutation;
   int i,j,done;
   char ministr[30];
   
   /* If mutate_mode = -1, choose a random mutation mode */
   if (mutate_mode == MUTATE_NOT_SPECIFIED) {
   
      randselect = flam3_random_isaac_01(rc);
      
      if (randselect < 0.1)
         mutate_mode = MUTATE_ALL_VARIATIONS;
      else if (randselect < 0.3)
         mutate_mode = MUTATE_ONE_XFORM_COEFS;
      else if (randselect < 0.5)
         mutate_mode = MUTATE_ADD_SYMMETRY;
      else if (randselect < 0.6)
         mutate_mode = MUTATE_POST_XFORMS;
      else if (randselect < 0.7)
         mutate_mode = MUTATE_COLOR_PALETTE;
      else if (randselect < 0.8)
         mutate_mode = MUTATE_DELETE_XFORM;
      else
         mutate_mode = MUTATE_ALL_COEFS;

   }
   
   memset(&mutation, 0, sizeof(flam3_genome));
      
   if (mutate_mode == MUTATE_ALL_VARIATIONS) {
   
      add_to_action(action,"mutate all variations");

      do {
         /* Create a random flame, and use the variations */
         /* to replace those in the original              */
         flam3_random(&mutation, ivars, ivars_n, sym, cp->num_xforms);
         for (i = 0; i < cp->num_xforms; i++) {
            for (j = 0; j < flam3_nvariations; j++) {
               if (cp->xform[i].var[j] != mutation.xform[i].var[j]) {
               
                  /* Copy the new var weights */
                  cp->xform[i].var[j] = mutation.xform[i].var[j];

                  /* Copy parameters for this variation only */
                  flam3_copy_params(&(cp->xform[i]),&(mutation.xform[i]),j);

                  done = 1;
               }
            }
         }
      } while (!done);
      
   } else if (mutate_mode == MUTATE_ONE_XFORM_COEFS) {
   
      int modxf;
   
      /* Generate a 2-xform random */
      flam3_random(&mutation, ivars, ivars_n, sym, 2);
      
      /* Which xform do we mutate? */
      modxf = ((unsigned)irand(rc)) % cp->num_xforms;
      
      add_to_action(action,"mutate xform ");
      sprintf(ministr,"%d coefs",modxf);
      add_to_action(action,ministr);
      
      /* if less than 3 xforms, then change only the translation part */
      if (2 >= cp->num_xforms) {
         for (j = 0; j < 2; j++)
            cp->xform[modxf].c[2][j] = mutation.xform[0].c[2][j];
      } else {
         for (i = 0; i < 3; i++)
            for (j = 0; j < 2; j++)
               cp->xform[modxf].c[i][j] = mutation.xform[0].c[i][j];
      }
      
   } else if (mutate_mode == MUTATE_ADD_SYMMETRY) {
   
      add_to_action(action,"mutate symmetry");
      flam3_add_symmetry(cp, 0);
      
   } else if (mutate_mode == MUTATE_POST_XFORMS) {
   
      int b = 1 + ((unsigned)irand(rc))%6;
      int same = ((unsigned)irand(rc))&3; /* 25% chance of using the same post for all of them */
      
      sprintf(ministr,"(%d%s)",b,(same>0) ? " same" : "");
      add_to_action(action,"mutate post xforms ");
      add_to_action(action,ministr);
      for (i = 0; i < cp->num_xforms; i++) {
         int copy = (i > 0) && same;

         if (copy) { /* Copy the post from the first xform to the rest of them */
            for (j = 0; j < 3; j++) {
               cp->xform[i].post[j][0] = cp->xform[0].post[j][0];
               cp->xform[i].post[j][1] = cp->xform[0].post[j][1];
            }

         } else {

            if (b&1) { /* 50% chance */
            
               double f = M_PI * flam3_random_isaac_11(rc);
               double t[2][2];

               t[0][0] = (cp->xform[i].c[0][0] * cos(f) + cp->xform[i].c[0][1] * -sin(f));
               t[0][1] = (cp->xform[i].c[0][0] * sin(f) + cp->xform[i].c[0][1] * cos(f));
               t[1][0] = (cp->xform[i].c[1][0] * cos(f) + cp->xform[i].c[1][1] * -sin(f));
               t[1][1] = (cp->xform[i].c[1][0] * sin(f) + cp->xform[i].c[1][1] * cos(f));

               cp->xform[i].c[0][0] = t[0][0];
               cp->xform[i].c[0][1] = t[0][1];
               cp->xform[i].c[1][0] = t[1][0];
               cp->xform[i].c[1][1] = t[1][1];

               f *= -1.0;

               t[0][0] = (cp->xform[i].post[0][0] * cos(f) + cp->xform[i].post[0][1] * -sin(f));
               t[0][1] = (cp->xform[i].post[0][0] * sin(f) + cp->xform[i].post[0][1] * cos(f));
               t[1][0] = (cp->xform[i].post[1][0] * cos(f) + cp->xform[i].post[1][1] * -sin(f));
               t[1][1] = (cp->xform[i].post[1][0] * sin(f) + cp->xform[i].post[1][1] * cos(f));

               cp->xform[i].post[0][0] = t[0][0];
               cp->xform[i].post[0][1] = t[0][1];
               cp->xform[i].post[1][0] = t[1][0];
               cp->xform[i].post[1][1] = t[1][1];

            }

            if (b&2) { /* 33% chance */
            
               double f = 0.2 + flam3_random_isaac_01(rc);
               double g = 0.2 + flam3_random_isaac_01(rc);

               if (flam3_random_isaac_bit(rc))
                  f = 1.0 / f;
               
               if (flam3_random_isaac_bit(rc))
                  g = f;
               else {               
                  if (flam3_random_isaac_bit(rc))
                     g = 1.0 / g;
               }

               cp->xform[i].c[0][0] /= f;
               cp->xform[i].c[0][1] /= f;
               cp->xform[i].c[1][1] /= g;
               cp->xform[i].c[1][0] /= g;
               cp->xform[i].post[0][0] *= f;
               cp->xform[i].post[1][0] *= f;
               cp->xform[i].post[0][1] *= g;
               cp->xform[i].post[1][1] *= g;
            }

            if (b&4) { /* 16% chance */

               double f = flam3_random_isaac_11(rc);
               double g = flam3_random_isaac_11(rc);

               cp->xform[i].c[2][0] -= f;
               cp->xform[i].c[2][1] -= g;
               cp->xform[i].post[2][0] += f;
               cp->xform[i].post[2][1] += g;
            }
         }
      }
   } else if (mutate_mode == MUTATE_COLOR_PALETTE) {
   
      double s = flam3_random_isaac_01(rc);

      if (s < 0.4) { /* randomize xform color coords */
      
         flam3_improve_colors(cp, 100, 0, 10);
         add_to_action(action,"mutate color coords");
         
      } else if (s < 0.8) { /* randomize xform color coords and palette */
      
         flam3_improve_colors(cp, 25, 1, 10);
         add_to_action(action,"mutate color all");
         
      } else { /* randomize palette only */

         cp->palette_index = flam3_get_palette(flam3_palette_random, cp->palette, cp->hue_rotation);
         /* if our palette retrieval fails, skip the mutation */
         if (cp->palette_index >= 0)
            add_to_action(action,"mutate color palette");
         else
            fprintf(stderr,"failure getting random palette, palette set to white\n");

      }
   } else if (mutate_mode == MUTATE_DELETE_XFORM) {
   
      int nx = ((unsigned)irand(rc))%cp->num_xforms;
      sprintf(ministr,"%d",nx);
      add_to_action(action,"mutate delete xform ");
      add_to_action(action,ministr);

      if (cp->num_xforms > 1)
         flam3_delete_xform(cp,nx);

   } else { /* MUTATE_ALL_COEFS */ 
   
      int x;
      add_to_action(action,"mutate all coefs");
      flam3_random(&mutation, ivars, ivars_n, sym, cp->num_xforms);

      /* change all the coefs by a fraction of the random */
      for (x = 0; x < cp->num_xforms; x++) {
         for (i = 0; i < 3; i++) {
            for (j = 0; j < 2; j++) {
               cp->xform[x].c[i][j] += speed * mutation.xform[x].c[i][j];

            }
         }
         /* Eventually, we can mutate the parametric variation coefs here. */
      }
   }
   
   clear_cp(&mutation,flam3_defaults_on);

}

static int random_var() {
  return random() % flam3_nvariations;
}

static int random_varn(int n) {
   return random() % n;
}

void flam3_random(flam3_genome *cp, int *ivars, int ivars_n, int sym, int spec_xforms) {

   int i, j, nxforms, var, samed, multid, samepost, postid, addfinal=0;
   int finum = -1;
   int n;
   char *ai;
   int f27 = argi("flam27",0);
   int mvar = f27 ? 54 : flam3_nvariations;
   double sum;

   static int xform_distrib[] = {
     2, 2, 2, 2,
     3, 3, 3, 3,
     4, 4, 4,
     5, 5,
     6
   };

   clear_cp(cp,flam3_defaults_on);

   cp->hue_rotation = (random()&7) ? 0.0 : flam3_random01();
   cp->palette_index = flam3_get_palette(flam3_palette_random, cp->palette, cp->hue_rotation);
   if (cp->palette_index < 0)
      fprintf(stderr,"error getting palette from xml file, setting to all white\n");
   cp->time = 0.0;
   cp->interpolation = flam3_interpolation_linear;
   cp->palette_interpolation = flam3_palette_interpolation_hsv;

   /* Choose the number of xforms */
   if (spec_xforms>0) {
      nxforms = spec_xforms;
      flam3_add_xforms(cp,nxforms,0,0);
   } else {
      nxforms = random_distrib(xform_distrib);
      flam3_add_xforms(cp,nxforms,0,0);
      /* Add a final xform 15% of the time */
      addfinal = flam3_random01() < 0.15;
      if (addfinal) {
         flam3_add_xforms(cp,1,0,1);
         nxforms = nxforms + addfinal;
         finum = nxforms-1;
      }
   }   

   /* If first input variation is 'flam3_variation_random' */
   /* choose one to use or decide to use multiple    */
   if (flam3_variation_random == ivars[0]) {
      if (flam3_random_bit()) {
         var = random_varn(mvar);
      } else {
         var = flam3_variation_random;
      }
   } else {
      var = flam3_variation_random_fromspecified;
   }

   samed = flam3_random_bit();
   multid = flam3_random_bit();
   postid = flam3_random01() < 0.6;
   samepost = flam3_random_bit();

   /* Loop over xforms */
   for (i = 0; i < nxforms; i++) {
      int j, k;
      cp->xform[i].density = 1.0 / nxforms;
      cp->xform[i].color = i&1;
      cp->xform[i].color_speed = 0.5;
      cp->xform[i].animate = 1.0;
      for (j = 0; j < 3; j++) {
         for (k = 0; k < 2; k++) {
            cp->xform[i].c[j][k] = flam3_random11();
            cp->xform[i].post[j][k] = (double)(k==j);
         }
      }
      
      if ( i != finum ) {
      
         if (!postid) {

            for (j = 0; j < 3; j++)
            for (k = 0; k < 2; k++) {
               if (samepost || (i==0))
                  cp->xform[i].post[j][k] = flam3_random11();
               else
                  cp->xform[i].post[j][k] = cp->xform[0].post[j][k];
            }
         }

         /* Clear all variation coefs */
         for (j = 0; j < flam3_nvariations; j++)
            cp->xform[i].var[j] = 0.0;

         if (flam3_variation_random != var && 
               flam3_variation_random_fromspecified != var) {

            /* Use only one variation specified for all xforms */
            cp->xform[i].var[var] = 1.0;

         } else if (multid && flam3_variation_random == var) {

           /* Choose a random var for this xform */
             cp->xform[i].var[random_varn(mvar)] = 1.0;

         } else {

            if (samed && i > 0) {

               /* Copy the same variations from the previous xform */
               for (j = 0; j < flam3_nvariations; j++) {
                  cp->xform[i].var[j] = cp->xform[i-1].var[j];
                  flam3_copy_params(&(cp->xform[i]),&(cp->xform[i-1]),j);
               }

            } else {

               /* Choose a random number of vars to use, at least 2 */
               /* but less than flam3_nvariations.Probability leans */
               /* towards fewer variations.                         */
               n = 2;
               while ((flam3_random_bit()) && (n<mvar))
                  n++;
               
               /* Randomly choose n variations, and change their weights. */
               /* A var can be selected more than once, further reducing  */
               /* the probability that multiple vars are used.            */
               for (j = 0; j < n; j++) {
                  if (flam3_variation_random_fromspecified != var)
                     cp->xform[i].var[random_varn(mvar)] = flam3_random01();
                  else
                     cp->xform[i].var[ivars[random_varn(ivars_n)]] = flam3_random01();
               }

               /* Normalize weights to 1.0 total. */
               sum = 0.0;
               for (j = 0; j < flam3_nvariations; j++)
                  sum += cp->xform[i].var[j];
               if (sum == 0.0)
                  cp->xform[i].var[random_var()] = 1.0;
               else {
                  for (j = 0; j < flam3_nvariations; j++)
                     cp->xform[i].var[j] /= sum;
               }
            }
         }
      } else {
         /* Handle final xform randomness. */
         n = 1;
         if (flam3_random_bit()) n++;
         
         /* Randomly choose n variations, and change their weights. */
         /* A var can be selected more than once, further reducing  */
         /* the probability that multiple vars are used.            */
         for (j = 0; j < n; j++) {
            if (flam3_variation_random_fromspecified != var)
               cp->xform[i].var[random_varn(mvar)] = flam3_random01();
            else
               cp->xform[i].var[ivars[random_varn(ivars_n)]] = flam3_random01();
         }

         /* Normalize weights to 1.0 total. */
         sum = 0.0;
         for (j = 0; j < flam3_nvariations; j++)
            sum += cp->xform[i].var[j];
         if (sum == 0.0)
            cp->xform[i].var[random_var()] = 1.0;
         else {
            for (j = 0; j < flam3_nvariations; j++)
               cp->xform[i].var[j] /= sum;
         }
      }  

      /* Generate random params for parametric variations, if selected. */
      if (cp->xform[i].var[VAR_BLOB] > 0) {
         /* Create random params for blob */
         cp->xform[i].blob_low = 0.2 + 0.5 * flam3_random01();
         cp->xform[i].blob_high = 0.8 + 0.4 * flam3_random01();
         cp->xform[i].blob_waves = (int)(2 + 5 * flam3_random01());
      }

      if (cp->xform[i].var[VAR_PDJ] > 0) {
         /* Create random params for PDJ */
         cp->xform[i].pdj_a = 3.0 * flam3_random11();
         cp->xform[i].pdj_b = 3.0 * flam3_random11();
         cp->xform[i].pdj_c = 3.0 * flam3_random11();
         cp->xform[i].pdj_d = 3.0 * flam3_random11();
      }

      if (cp->xform[i].var[VAR_FAN2] > 0) {
         /* Create random params for fan2 */
         cp->xform[i].fan2_x = flam3_random11();
         cp->xform[i].fan2_y = flam3_random11();
      }

      if (cp->xform[i].var[VAR_RINGS2] > 0) {
         /* Create random params for rings2 */
         cp->xform[i].rings2_val = 2*flam3_random01();
      }

      if (cp->xform[i].var[VAR_PERSPECTIVE] > 0) {

         /* Create random params for perspective */
         cp->xform[i].perspective_angle = flam3_random01();
         cp->xform[i].perspective_dist = 2*flam3_random01() + 1.0;

      }

      if (cp->xform[i].var[VAR_JULIAN] > 0) {

         /* Create random params for julian */
         cp->xform[i].julian_power = (int)(5*flam3_random01() + 2);
         cp->xform[i].julian_dist = 1.0;

      }

      if (cp->xform[i].var[VAR_JULIASCOPE] > 0) {

         /* Create random params for juliaScope */
         cp->xform[i].juliascope_power = (int)(5*flam3_random01() + 2);
         cp->xform[i].juliascope_dist = 1.0;

      }

      if (cp->xform[i].var[VAR_RADIAL_BLUR] > 0) {

         /* Create random params for radialBlur */
         cp->xform[i].radial_blur_angle = (2 * flam3_random01() - 1);

      }

      if (cp->xform[i].var[VAR_PIE] > 0) {
         /* Create random params for pie */
         cp->xform[i].pie_slices = (int) 10.0*flam3_random01();
         cp->xform[i].pie_thickness = flam3_random01();
         cp->xform[i].pie_rotation = 2.0 * M_PI * flam3_random11();
      }

      if (cp->xform[i].var[VAR_NGON] > 0) {
         /* Create random params for ngon */
         cp->xform[i].ngon_sides = (int) flam3_random01()* 10 + 3;
         cp->xform[i].ngon_power = 3*flam3_random01() + 1;
         cp->xform[i].ngon_circle = 3*flam3_random01();
         cp->xform[i].ngon_corners = 2*flam3_random01()*cp->xform[i].ngon_circle;
      }

      if (cp->xform[i].var[VAR_CURL] > 0) {
         /* Create random params for curl */
         cp->xform[i].curl_c1 = flam3_random01();
         cp->xform[i].curl_c2 = flam3_random01();
      }

      if (cp->xform[i].var[VAR_RECTANGLES] > 0) {
         /* Create random params for rectangles */
         cp->xform[i].rectangles_x = flam3_random01();
         cp->xform[i].rectangles_y = flam3_random01();
      }

      if (cp->xform[i].var[VAR_DISC2] > 0) {
      /* Create random params for disc2 */
      cp->xform[i].disc2_rot = 0.5 * flam3_random01();
      cp->xform[i].disc2_twist = 0.5 * flam3_random01();

      }

      if (cp->xform[i].var[VAR_SUPER_SHAPE] > 0) {
         /* Create random params for supershape */
         cp->xform[i].super_shape_rnd = flam3_random01();
         cp->xform[i].super_shape_m = (int) flam3_random01()*6;
         cp->xform[i].super_shape_n1 = flam3_random01()*40;
         cp->xform[i].super_shape_n2 = flam3_random01()*20;
         cp->xform[i].super_shape_n3 = cp->xform[i].super_shape_n2;
         cp->xform[i].super_shape_holes = 0.0;
      }

      if (cp->xform[i].var[VAR_FLOWER] > 0) {
         /* Create random params for flower */
         cp->xform[i].flower_petals = 4 * flam3_random01();
         cp->xform[i].flower_holes = flam3_random01();
      }

      if (cp->xform[i].var[VAR_CONIC] > 0) {
         /* Create random params for conic */
         cp->xform[i].conic_eccentricity = flam3_random01();
         cp->xform[i].conic_holes = flam3_random01();
      }

      if (cp->xform[i].var[VAR_PARABOLA] > 0) {
         /* Create random params for parabola */
         cp->xform[i].parabola_height = 0.5 + flam3_random01();
         cp->xform[i].parabola_width = 0.5 + flam3_random01();
      }

      if (cp->xform[i].var[VAR_BENT2] > 0) {
         /* Create random params for bent2 */
         cp->xform[i].bent2_x = 3*(-0.5 + flam3_random01());
         cp->xform[i].bent2_y = 3*(-0.5 + flam3_random01());
      }

      if (cp->xform[i].var[VAR_BIPOLAR] > 0) {
         /* Create random params for bipolar */
         cp->xform[i].bipolar_shift = 2.0 * flam3_random01() - 1;
      }

      if (cp->xform[i].var[VAR_CELL] > 0) {
         /* Create random params for cell */
         cp->xform[i].cell_size = 2.0 * flam3_random01() + 0.5;
      }

      if (cp->xform[i].var[VAR_CPOW] > 0) {
         /* Create random params for cpow */
         cp->xform[i].cpow_r = 3.0 * flam3_random01();
         cp->xform[i].cpow_i = flam3_random01() - 0.5;
         cp->xform[i].cpow_power = (int)(5.0 * flam3_random01());
      }

      if (cp->xform[i].var[VAR_CURVE] > 0) {
         /* Create random params for curve */
         cp->xform[i].curve_xamp = 5 * (flam3_random01()-.5);
         cp->xform[i].curve_yamp = 4 * (flam3_random01()-.5);
         cp->xform[i].curve_xlength = 2 * (flam3_random01()+.5);
         cp->xform[i].curve_ylength = 2 * (flam3_random01()+.5);
      }

      if (cp->xform[i].var[VAR_ESCHER] > 0) {
         /* Create random params for escher */
         cp->xform[i].escher_beta = M_PI * flam3_random11();
      }

      if (cp->xform[i].var[VAR_LAZYSUSAN] > 0) {
         /* Create random params for lazysusan */
         cp->xform[i].lazysusan_x = 2.0*flam3_random11();
         cp->xform[i].lazysusan_y = 2.0*flam3_random11();
         cp->xform[i].lazysusan_spin = M_PI*flam3_random11();
         cp->xform[i].lazysusan_space = 2.0*flam3_random11();
         cp->xform[i].lazysusan_twist = 2.0*flam3_random11();
      }

      if (cp->xform[i].var[VAR_MODULUS] > 0) {
         /* Create random params for modulus */
         cp->xform[i].modulus_x = flam3_random11();
         cp->xform[i].modulus_y = flam3_random11();
      }

      if (cp->xform[i].var[VAR_OSCILLOSCOPE] > 0) {
         /* Create random params for oscope */
         cp->xform[i].oscope_separation = 1.0 + flam3_random11();
         cp->xform[i].oscope_frequency = M_PI * flam3_random11();
         cp->xform[i].oscope_amplitude = 1.0 + 2 * flam3_random01();
         cp->xform[i].oscope_damping = flam3_random01();
      }

      if (cp->xform[i].var[VAR_POPCORN2] > 0) {
         /* Create random params for popcorn2 */
         cp->xform[i].popcorn2_x = 0.2 * flam3_random01();
         cp->xform[i].popcorn2_y = 0.2 * flam3_random01();
         cp->xform[i].popcorn2_c = 5 * flam3_random01();
      }

      if (cp->xform[i].var[VAR_SEPARATION] > 0) {
         /* Create random params for separation */
         cp->xform[i].separation_x = 1 + flam3_random11();
         cp->xform[i].separation_y = 1 + flam3_random11();
         cp->xform[i].separation_xinside = flam3_random11();
         cp->xform[i].separation_yinside = flam3_random11();
      }

      if (cp->xform[i].var[VAR_SPLIT] > 0) {
         /* Create random params for split */
         cp->xform[i].split_xsize = flam3_random11();
         cp->xform[i].split_ysize = flam3_random11();
      }

      if (cp->xform[i].var[VAR_SPLITS] > 0) {
         /* Create random params for splits */
         cp->xform[i].splits_x = flam3_random11();
         cp->xform[i].splits_y = flam3_random11();
      }

      if (cp->xform[i].var[VAR_STRIPES] > 0) {
         /* Create random params for stripes */
         cp->xform[i].stripes_space = flam3_random01();
         cp->xform[i].stripes_warp = 5*flam3_random01();
      }

      if (cp->xform[i].var[VAR_WEDGE] > 0) {
         /* Create random params for wedge */
         cp->xform[i].wedge_angle = M_PI*flam3_random01();
         cp->xform[i].wedge_hole = 0.5*flam3_random11();
         cp->xform[i].wedge_count = floor(5*flam3_random01())+1;
         cp->xform[i].wedge_swirl = flam3_random01();
      }

      if (cp->xform[i].var[VAR_WEDGE_JULIA] > 0) {

         /* Create random params for wedge_julia */
         cp->xform[i].wedge_julia_power = (int)(5*flam3_random01() + 2);
         cp->xform[i].wedge_julia_dist = 1.0;
         cp->xform[i].wedge_julia_count = (int)(3*flam3_random01() + 1);
         cp->xform[i].wedge_julia_angle = M_PI * flam3_random01();

      }

      if (cp->xform[i].var[VAR_WEDGE_SPH] > 0) {
         /* Create random params for wedge_sph */
         cp->xform[i].wedge_sph_angle = M_PI*flam3_random01();
         cp->xform[i].wedge_sph_hole = 0.5*flam3_random11();
         cp->xform[i].wedge_sph_count = floor(5*flam3_random01())+1;
         cp->xform[i].wedge_sph_swirl = flam3_random01();
      }

      if (cp->xform[i].var[VAR_WHORL] > 0) {
         /* Create random params for whorl */
         cp->xform[i].whorl_inside = flam3_random01();
         cp->xform[i].whorl_outside = flam3_random01();
      }

      if (cp->xform[i].var[VAR_WAVES2] > 0) {
         /* Create random params for waves2 */
         cp->xform[i].waves2_scalex = 0.5 + flam3_random01();
         cp->xform[i].waves2_scaley = 0.5 + flam3_random01();
         cp->xform[i].waves2_freqx = 4 * flam3_random01();
         cp->xform[i].waves2_freqy = 4 * flam3_random01();
      }         

      if (cp->xform[i].var[VAR_AUGER] > 0) {
         /* Create random params for auger */
         cp->xform[i].auger_sym = 0;
         cp->xform[i].auger_weight = 0.5 + flam3_random01()/2.0;
         cp->xform[i].auger_freq = floor(5*flam3_random01())+1;
         cp->xform[i].auger_scale = flam3_random01();
      }         

      if (cp->xform[i].var[VAR_FLUX] > 0) {
         /* Create random params for flux */
         cp->xform[i].flux_spread = 0.5 + flam3_random01()/2.0;
      }

      if (cp->xform[i].var[VAR_MOBIUS] > 0) {
         /* Create random params for mobius */
         cp->xform[i].mobius_re_a = flam3_random11();
         cp->xform[i].mobius_im_a = flam3_random11();
         cp->xform[i].mobius_re_b = flam3_random11();
         cp->xform[i].mobius_im_b = flam3_random11();
         cp->xform[i].mobius_re_c = flam3_random11();
         cp->xform[i].mobius_im_c = flam3_random11();
         cp->xform[i].mobius_re_d = flam3_random11();
         cp->xform[i].mobius_im_d = flam3_random11();
      }

   }

   /* Randomly add symmetry (but not if we've already added a final xform) */
   if (sym || (!(random()%4) && !addfinal))
      flam3_add_symmetry(cp, sym);
   else
      cp->symmetry = 0;

   //qsort((char *) cp->xform, (cp->num_xforms-addfinal), sizeof(flam3_xform), compare_xforms);


}


static int sort_by_x(const void *av, const void *bv) {
    double *a = (double *) av;
    double *b = (double *) bv;
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    return 0;
}

static int sort_by_y(const void *av, const void *bv) {
    double *a = (double *) av;
    double *b = (double *) bv;
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    return 0;
}


/* Memory helper functions because 

    Python on Windows uses the MSVCR71.dll version of the C Runtime and 
    mingw uses the MSVCRT.dll version. */

void *flam3_malloc(size_t size) {

   return (malloc(size));
   
}

void flam3_free(void *ptr) {

   free(ptr);
   
}

/*
 * find a 2d bounding box that does not enclose eps of the fractal density
 * in each compass direction.
 */
int flam3_estimate_bounding_box(flam3_genome *cp, double eps, int nsamples,
             double *bmin, double *bmax, randctx *rc) {
   int i;
   int low_target, high_target;
   double min[2], max[2];
   double *points;
   int bv;
   unsigned short *xform_distrib;

   if (nsamples <= 0) nsamples = 10000;

   points = (double *) malloc(sizeof(double) * 4 * nsamples);
   points[0] = flam3_random_isaac_11(rc);
   points[1] = flam3_random_isaac_11(rc);
   points[2] = 0.0;
   points[3] = 0.0;

   if (prepare_precalc_flags(cp))
      return(-1);
   xform_distrib = flam3_create_xform_distrib(cp);
   if (xform_distrib==NULL)
      return(-1);
   bv=flam3_iterate(cp, nsamples, 20, points, xform_distrib, rc);
   free(xform_distrib);
      
   if ( bv/(double)nsamples > eps )
      eps = 3*bv/(double)nsamples;
   
   if ( eps > 0.3 )
      eps = 0.3;
      
   low_target = (int)(nsamples * eps);
   high_target = nsamples - low_target;

   
   min[0] = min[1] =  1e10;
   max[0] = max[1] = -1e10;

   for (i = 0; i < nsamples; i++) {
      double *p = &points[4*i];
      if (p[0] < min[0]) min[0] = p[0];
      if (p[1] < min[1]) min[1] = p[1];
      if (p[0] > max[0]) max[0] = p[0];
      if (p[1] > max[1]) max[1] = p[1];
   }

   if (low_target == 0) {
      bmin[0] = min[0];
      bmin[1] = min[1];
      bmax[0] = max[0];
      bmax[1] = max[1];
      free(points);
      return(bv);
   }

   qsort(points, nsamples, sizeof(double) * 4, sort_by_x);
   bmin[0] = points[4 * low_target];
   bmax[0] = points[4 * high_target];

   qsort(points, nsamples, sizeof(double) * 4, sort_by_y);
   bmin[1] = points[4 * low_target + 1];
   bmax[1] = points[4 * high_target + 1];
   free(points);
   
   return(bv);
}


typedef double bucket_double[5];
typedef double abucket_double[4];
typedef unsigned int bucket_int[5];
typedef unsigned int abucket_int[4];
typedef float bucket_float[5];
typedef float abucket_float[4];

#ifdef HAVE_GCC_64BIT_ATOMIC_OPS
static inline void
double_atomic_add(double *dest, double delta)
{
   uint64_t *int_ptr = (uint64_t *)dest;
   union {
      double dblval;
      uint64_t intval;
   } old_val, new_val;
   int success;

   do {
      old_val.dblval = *dest;
      new_val.dblval = old_val.dblval + delta;
      success = __sync_bool_compare_and_swap(
         int_ptr, old_val.intval, new_val.intval);
   } while (!success);
}
#endif /* HAVE_GCC_64BIT_ATOMIC_OPS */

#ifdef HAVE_GCC_ATOMIC_OPS
static inline void
float_atomic_add(float *dest, float delta)
{
   uint32_t *int_ptr = (uint32_t *)dest;
   union {
      float fltval;
      uint32_t intval;
   } old_val, new_val;
   int success;

   do {
      old_val.fltval = *dest;
      new_val.fltval = old_val.fltval + delta;
      success = __sync_bool_compare_and_swap(
         int_ptr, old_val.intval, new_val.intval);
   } while (!success);
}

static inline void
uint_atomic_add(unsigned int *dest, unsigned int delta)
{
   unsigned int old_val, new_val;
   int success;

   do {
      old_val = *dest;
      if (UINT_MAX - old_val > delta)
         new_val = old_val + delta;
      else
         new_val = UINT_MAX;
      success = __sync_bool_compare_and_swap(
         dest, old_val, new_val);
   } while (!success);
}

static inline void
ushort_atomic_add(unsigned short *dest, unsigned short delta)
{
   unsigned short old_val, new_val;
   int success;

   do {
      old_val = *dest;
      if (USHRT_MAX - old_val > delta)
         new_val = old_val + delta;
      else
         new_val = USHRT_MAX;
      success = __sync_bool_compare_and_swap(
         dest, old_val, new_val);
   } while (!success);
}
#endif /* HAVE_GCC_ATOMIC_OPS */

/* 64-bit datatypes */
#define bucket bucket_double
#define abucket abucket_double
#define abump_no_overflow(dest, delta) do {dest += delta;} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta)  do {dest += delta;} while (0)
#define render_rectangle render_rectangle_double
#define iter_thread iter_thread_double
#define de_thread_helper de_thread_helper_64
#define de_thread de_thread_64
#include "rect.c"
#ifdef HAVE_GCC_64BIT_ATOMIC_OPS
   /* multi-threaded */
   #undef USE_LOCKS
   #undef bump_no_overflow
   #undef render_rectangle
   #undef iter_thread
   #undef de_thread_helper
   #undef de_thread
   #define bump_no_overflow(dest, delta)  double_atomic_add(&dest, delta)
   #define render_rectangle render_rectangle_double_mt
   #define iter_thread iter_thread_double_mt
   #define de_thread_helper de_thread_helper_64_mt
   #define de_thread de_thread_64_mt
   #include "rect.c"
#else /* !HAVE_GCC_64BIT_ATOMIC_OPS */
   #define render_rectangle_double_mt render_rectangle_double
#endif /* HAVE_GCC_64BIT_ATOMIC_OPS */
#undef render_rectangle
#undef iter_thread
#undef add_c_to_accum
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow
#undef de_thread_helper
#undef de_thread

/* 32-bit datatypes */
#define bucket bucket_int
#define abucket abucket_int
#define abump_no_overflow(dest, delta) do { \
   if (UINT_MAX - dest > delta) dest += delta; else dest = UINT_MAX; \
} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta) do { \
   if (UINT_MAX - dest > delta) dest += delta; else dest = UINT_MAX; \
} while (0)
#define render_rectangle render_rectangle_int
#define iter_thread iter_thread_int
#define de_thread_helper de_thread_helper_32
#define de_thread de_thread_32
#include "rect.c"
#ifdef HAVE_GCC_ATOMIC_OPS
   /* multi-threaded */
   #undef USE_LOCKS
   #undef bump_no_overflow
   #undef render_rectangle
   #undef iter_thread
   #undef de_thread_helper
   #undef de_thread
   #define bump_no_overflow(dest, delta)  uint_atomic_add(&dest, delta)
   #define render_rectangle render_rectangle_int_mt
   #define iter_thread iter_thread_int_mt
   #define de_thread_helper de_thread_helper_32_mt
   #define de_thread de_thread_32_mt
   #include "rect.c"
#else /* !HAVE_GCC_ATOMIC_OPS */
   #define render_rectangle_int_mt render_rectangle_int
#endif /* HAVE_GCC_ATOMIC_OPS */
#undef iter_thread
#undef render_rectangle
#undef add_c_to_accum
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow
#undef de_thread_helper
#undef de_thread

/* experimental 32-bit datatypes (called 33) */
#define bucket bucket_int
#define abucket abucket_float
#define abump_no_overflow(dest, delta) do {dest += delta;} while (0)
#define add_c_to_accum(acc,i,ii,j,jj,wid,hgt,c) do { \
   if ( (j) + (jj) >=0 && (j) + (jj) < (hgt) && (i) + (ii) >=0 && (i) + (ii) < (wid)) { \
   abucket *a = (acc) + ( (i) + (ii) ) + ( (j) + (jj) ) * (wid); \
   abump_no_overflow(a[0][0],(c)[0]); \
   abump_no_overflow(a[0][1],(c)[1]); \
   abump_no_overflow(a[0][2],(c)[2]); \
   abump_no_overflow(a[0][3],(c)[3]); \
   } \
} while (0)
/* single-threaded */
#define USE_LOCKS
#define bump_no_overflow(dest, delta) do { \
   if (UINT_MAX - dest > delta) dest += delta; else dest = UINT_MAX; \
} while (0)
#define render_rectangle render_rectangle_float
#define iter_thread iter_thread_float
#define de_thread_helper de_thread_helper_33
#define de_thread de_thread_33
#include "rect.c"
#ifdef HAVE_GCC_ATOMIC_OPS
   /* multi-threaded */
   #undef USE_LOCKS
   #undef bump_no_overflow
   #undef render_rectangle
   #undef iter_thread
   #undef de_thread_helper
   #undef de_thread
   #define bump_no_overflow(dest, delta)  uint_atomic_add(&dest, delta)
   #define render_rectangle render_rectangle_float_mt
   #define iter_thread iter_thread_float_mt
   #define de_thread_helper de_thread_helper_33_mt
   #define de_thread de_thread_33_mt
   #include "rect.c"
#else /* !HAVE_GCC_ATOMIC_OPS */
   #define render_rectangle_float_mt render_rectangle_float
#endif /* HAVE_GCC_ATOMIC_OPS */
#undef iter_thread
#undef render_rectangle
#undef add_c_to_accum
#undef bucket
#undef abucket
#undef bump_no_overflow
#undef abump_no_overflow
#undef de_thread_helper
#undef de_thread


double flam3_render_memory_required(flam3_frame *spec)
{
  flam3_genome *cps = spec->genomes;
  int real_bits = spec->bits;
  int real_bytes;

  if (33 == real_bits) real_bits = 32;

  real_bytes = real_bits / 8;

  return
    (double) cps[0].spatial_oversample * cps[0].spatial_oversample *
    (double) cps[0].width * cps[0].height * real_bytes * 9.0;
}

void bits_error(flam3_frame *spec) {
      fprintf(stderr, "flam3: bits must be 32, 33, or 64 not %d.\n",
         spec->bits);
}

int flam3_render(flam3_frame *spec, void *out,
        int field, int nchan, int trans, stat_struct *stats) {
         
  int retval;
  
  if (spec->nthreads <= 2) {
    /* single-threaded or 2 threads without atomic operations */
    switch (spec->bits) {
    case 32:
      retval = render_rectangle_int(spec, out, field, nchan, trans, stats);
      return(retval);
    case 33:
      retval = render_rectangle_float(spec, out, field, nchan, trans, stats);
      return(retval);
    case 64:
      retval = render_rectangle_double(spec, out, field, nchan, trans, stats);
      return(retval);
    default:
      bits_error(spec);
      return(1);
    }
  } else {
    /* 3+ threads using atomic ops if available */
    switch (spec->bits) {
    case 32:
      retval = render_rectangle_int_mt(spec, out, field, nchan, trans, stats);
      return(retval);
    case 33:
      retval = render_rectangle_float_mt(spec, out, field, nchan, trans, stats);
      return(retval);
    case 64:
      retval = render_rectangle_double_mt(spec, out, field, nchan, trans, stats);
      return(retval);
    default:
      bits_error(spec);
      return(1);
    }
  }
}


void flam3_srandom() {
   unsigned int seed;
   char *s = getenv("seed");

   if (s)
      seed = atoi(s);
   else
      seed = time(0) + getpid();

   srandom(seed);
}


/* correlation dimension, after clint sprott.
   computes slope of the correlation sum at a size scale
   the order of 2% the size of the attractor or the camera. */
double flam3_dimension(flam3_genome *cp, int ntries, int clip_to_camera) {
  double fd;
  double *hist;
  double bmin[2];
  double bmax[2];
  double d2max;
  int lp;
  long int default_isaac_seed = (long int)time(0);
  randctx rc;
  int SBS = 10000;
  int i, n1=0, n2=0, got, nclipped;

  /* Set up the isaac rng */
  for (lp = 0; lp < RANDSIZ; lp++)
     rc.randrsl[lp] = default_isaac_seed;

  irandinit(&rc,1);

  if (ntries < 2) ntries = 3000*1000;

  if (clip_to_camera) {
    double scale, ppux, corner0, corner1;
    scale = pow(2.0, cp->zoom);
    ppux = cp->pixels_per_unit * scale;
    corner0 = cp->center[0] - cp->width / ppux / 2.0;
    corner1 = cp->center[1] - cp->height / ppux / 2.0;
    bmin[0] = corner0;
    bmin[1] = corner1;
    bmax[0] = corner0 + cp->width  / ppux;
    bmax[1] = corner1 + cp->height / ppux;
  } else {
    if (flam3_estimate_bounding_box(cp, 0.0, 0, bmin, bmax, &rc)<0)
       return(-1.0);
       
  }

  d2max =
    (bmax[0] - bmin[0]) * (bmax[0] - bmin[0]) +
    (bmax[1] - bmin[1]) * (bmax[1] - bmin[1]);

  //  fprintf(stderr, "d2max=%g %g %g %g %g\n", d2max,
  //  bmin[0], bmin[1], bmax[0], bmax[1]);

  hist = malloc(2 * ntries * sizeof(double));

  got = 0;
  nclipped = 0;
  while (got < 2*ntries) {
    double subb[40000];
    int i4, clipped;
    unsigned short *xform_distrib;
    subb[0] = flam3_random_isaac_11(&rc);
    subb[1] = flam3_random_isaac_11(&rc);
    subb[2] = 0.0;
    subb[3] = 0.0;
    if (prepare_precalc_flags(cp))
      return(-1.0);
    xform_distrib = flam3_create_xform_distrib(cp);
    if (xform_distrib==NULL)
      return(-1.0);
    flam3_iterate(cp, SBS, 20, subb, xform_distrib, &rc);
    free(xform_distrib);
    i4 = 0;
    for (i = 0; i < SBS; i++) {
      if (got == 2*ntries) break;
      clipped = clip_to_camera &&
   ((subb[i4] < bmin[0]) ||
    (subb[i4+1] < bmin[1]) ||
    (subb[i4] > bmax[0]) ||
    (subb[i4+1] > bmax[1]));
      if (!clipped) {
   hist[got] = subb[i4];
   hist[got+1] = subb[i4+1];
   got += 2;
      } else {
   nclipped++;
   if (nclipped > 10 * ntries) {
       fprintf(stderr, "warning: too much clipping, "
          "flam3_dimension giving up.\n");
       return sqrt(-1.0);
   }
      }
      i4 += 4;
    }
  }
  if (0)
    fprintf(stderr, "cliprate=%g\n", nclipped/(ntries+(double)nclipped));

  for (i = 0; i < ntries; i++) {
    int ri;
    double dx, dy, d2;
    double tx, ty;

    tx = hist[2*i];
    ty = hist[2*i+1];

    do {
      ri = 2 * (random() % ntries);
    } while (ri == i);

    dx = hist[ri] - tx;
    dy = hist[ri+1] - ty;
    d2 = dx*dx + dy*dy;
    if (d2 < 0.004 * d2max) n2++;
    if (d2 < 0.00004 * d2max) n1++;
  }

  fd = 0.434294 * log(n2 / (n1 - 0.5));

  if (0)
    fprintf(stderr, "n1=%d n2=%d\n", n1, n2);

  free(hist);
  return fd;
}

double flam3_lyapunov(flam3_genome *cp, int ntries) {
  double p[4];
  double x, y;
  double xn, yn;
  double xn2, yn2;
  double dx, dy, r;
  double eps = 1e-5;
  int i;
  double sum = 0.0;
  unsigned short *xform_distrib;

  int lp;
  long int default_isaac_seed = (long int)time(0);
  randctx rc;

  /* Set up the isaac rng */
  for (lp = 0; lp < RANDSIZ; lp++)
     rc.randrsl[lp] = default_isaac_seed;

  irandinit(&rc,1);


  if (ntries < 1) ntries = 10000;

  for (i = 0; i < ntries; i++) {
    x = flam3_random_isaac_11(&rc);
    y = flam3_random_isaac_11(&rc);

    p[0] = x;
    p[1] = y;
    p[2] = 0.0;
    p[3] = 0.0;

    // get into the attractor
    if (prepare_precalc_flags(cp))
      return(-1.0);
    xform_distrib = flam3_create_xform_distrib(cp);
    if (xform_distrib==NULL)
      return(-1.0);
      
    flam3_iterate(cp, 1, 20+(random()%10), p, xform_distrib, &rc);
    free(xform_distrib);

    x = p[0];
    y = p[1];

    // take one deterministic step
    srandom(i);

    if (prepare_precalc_flags(cp))
      return(-1.0);
    xform_distrib = flam3_create_xform_distrib(cp);
    if (xform_distrib==NULL)
      return(-1.0);

    flam3_iterate(cp, 1, 0, p, xform_distrib, &rc);
    free(xform_distrib);

    xn = p[0];
    yn = p[1];

    do {
      dx = flam3_random_isaac_11(&rc);
      dy = flam3_random_isaac_11(&rc);
      r = sqrt(dx * dx + dy * dy);
    } while (r == 0.0);
    dx /= r;
    dy /= r;

    dx *= eps;
    dy *= eps;

    p[0] = x + dx;
    p[1] = y + dy;
    p[2] = 0.0;

    // take the same step but with eps
    srandom(i);
    if (prepare_precalc_flags(cp))
      return(-1.0);
    xform_distrib = flam3_create_xform_distrib(cp);
    if (xform_distrib==NULL)
      return(-1.0);

    flam3_iterate(cp, 1, 0, p, xform_distrib, &rc);
    free(xform_distrib);

    xn2 = p[0];
    yn2 = p[1];

    r = sqrt((xn-xn2)*(xn-xn2) + (yn-yn2)*(yn-yn2));

    sum += log(r/eps);
  }
  return sum/(log(2.0)*ntries);
}

