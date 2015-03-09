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

#include "filters.h"


/*
 * filter function definitions
 * from Graphics Gems III code
 * and ImageMagick resize.c
 */


double flam3_spatial_support[flam3_num_spatialfilters] = {

   1.5, /* gaussian */
   1.0, /* hermite */
   0.5, /* box */
   1.0, /* triangle */
   1.5, /* bell */
   2.0, /* b spline */
   2.0, /* mitchell */
   1.0, /* blackman */
   2.0, /* catrom */
   1.0, /* hanning */
   1.0, /* hamming */
   3.0, /* lanczos3 */
   2.0, /* lanczos2 */
   1.5  /* quadratic */
};

double flam3_hermite_filter(double t) {
   /* f(t) = 2|t|^3 - 3|t|^2 + 1, -1 <= t <= 1 */
   if(t < 0.0) t = -t;
   if(t < 1.0) return((2.0 * t - 3.0) * t * t + 1.0);
   return(0.0);
}

double flam3_box_filter(double t) {
   if((t > -0.5) && (t <= 0.5)) return(1.0);
   return(0.0);
}

double flam3_triangle_filter(double t) {
   if(t < 0.0) t = -t;
   if(t < 1.0) return(1.0 - t);
   return(0.0);
}

double flam3_bell_filter(double t) {
   /* box (*) box (*) box */
   if(t < 0) t = -t;
   if(t < .5) return(.75 - (t * t));
   if(t < 1.5) {
      t = (t - 1.5);
      return(.5 * (t * t));
   }
   return(0.0);
}

double flam3_b_spline_filter(double t) {

   /* box (*) box (*) box (*) box */
   double tt;

   if(t < 0) t = -t;
   if(t < 1) {
      tt = t * t;
      return((.5 * tt * t) - tt + (2.0 / 3.0));
   } else if(t < 2) {
      t = 2 - t;
      return((1.0 / 6.0) * (t * t * t));
   }
   return(0.0);
}

double flam3_sinc(double x) {
   x *= M_PI;
   if(x != 0) return(sin(x) / x);
   return(1.0);
}

double flam3_blackman_filter(double x) {
  return(0.42+0.5*cos(M_PI*x)+0.08*cos(2*M_PI*x));
}

double flam3_catrom_filter(double x) {
  if (x < -2.0)
    return(0.0);
  if (x < -1.0)
    return(0.5*(4.0+x*(8.0+x*(5.0+x))));
  if (x < 0.0)
    return(0.5*(2.0+x*x*(-5.0-3.0*x)));
  if (x < 1.0)
    return(0.5*(2.0+x*x*(-5.0+3.0*x)));
  if (x < 2.0)
    return(0.5*(4.0+x*(-8.0+x*(5.0-x))));
  return(0.0);
}

double flam3_mitchell_filter(double t) {
   double tt;

   tt = t * t;
   if(t < 0) t = -t;
   if(t < 1.0) {
      t = (((12.0 - 9.0 * flam3_mitchell_b - 6.0 * flam3_mitchell_c) * (t * tt))
         + ((-18.0 + 12.0 * flam3_mitchell_b + 6.0 * flam3_mitchell_c) * tt)
         + (6.0 - 2 * flam3_mitchell_b));
      return(t / 6.0);
   } else if(t < 2.0) {
      t = (((-1.0 * flam3_mitchell_b - 6.0 * flam3_mitchell_c) * (t * tt))
         + ((6.0 * flam3_mitchell_b + 30.0 * flam3_mitchell_c) * tt)
         + ((-12.0 * flam3_mitchell_b - 48.0 * flam3_mitchell_c) * t)
         + (8.0 * flam3_mitchell_b + 24 * flam3_mitchell_c));
      return(t / 6.0);
   }
   return(0.0);
}

double flam3_hanning_filter(double x) {
  return(0.5+0.5*cos(M_PI*x));
}

double flam3_hamming_filter(double x) {
  return(0.54+0.46*cos(M_PI*x));
}

double flam3_lanczos3_filter(double t) {
   if(t < 0) t = -t;
   if(t < 3.0) return(flam3_sinc(t) * flam3_sinc(t/3.0));
   return(0.0);
}

double flam3_lanczos2_filter(double t) {
   if(t < 0) t = -t;
   if(t < 2.0) return(flam3_sinc(t) * flam3_sinc(t/2.0));
   return(0.0);
}

double flam3_gaussian_filter(double x) {
  return(exp((-2.0*x*x))*sqrt(2.0/M_PI));
}

double flam3_quadratic_filter(double x) {
  if (x < -1.5)
    return(0.0);
  if (x < -0.5)
    return(0.5*(x+1.5)*(x+1.5));
  if (x < 0.5)
    return(0.75-x*x);
  if (x < 1.5)
    return(0.5*(x-1.5)*(x-1.5));
  return(0.0);
}

double flam3_spatial_filter(int knum, double x) {

   if (knum==0)
      return flam3_gaussian_filter(x);
   else if (knum==1)
      return flam3_hermite_filter(x);
   else if (knum==2)
      return flam3_box_filter(x);
   else if (knum==3)
      return flam3_triangle_filter(x);
   else if (knum==4)
      return flam3_bell_filter(x);
   else if (knum==5)
      return flam3_b_spline_filter(x);
   else if (knum==6)
      return flam3_mitchell_filter(x);
   else if (knum==7)
      return flam3_sinc(x)*flam3_blackman_filter(x);
   else if (knum==8)
      return flam3_catrom_filter(x);
   else if (knum==9)
      return flam3_sinc(x)*flam3_hanning_filter(x);
   else if (knum==10)
      return flam3_sinc(x)*flam3_hamming_filter(x);
   else if (knum==11)
      return flam3_lanczos3_filter(x)*flam3_sinc(x/3.0);   
   else if (knum==12)
      return flam3_lanczos2_filter(x)*flam3_sinc(x/2.0);
   else  // if (knum==13)
      return flam3_quadratic_filter(x);
}

int normalize_vector(double *v, int n) {
   double t = 0.0;
   int i;
   for (i = 0; i < n; i++)
      t += v[i];
   if (0.0 == t) return 1;
   t = 1.0 / t;
   for (i = 0; i < n; i++)
      v[i] *= t;
   return 0;
}


int flam3_create_spatial_filter(flam3_frame *spec, int field, double **filter) {

   int sf_kernel = spec->genomes[0].spatial_filter_select;
   int supersample = spec->genomes[0].spatial_oversample;
   double sf_radius = spec->genomes[0].spatial_filter_radius;
   double aspect_ratio = spec->pixel_aspect_ratio;   
   double sf_supp = flam3_spatial_support[sf_kernel];
   
   double fw = 2.0 * sf_supp * supersample * sf_radius / aspect_ratio;
   double adjust, ii, jj;
   
   int fwidth = ((int) fw) + 1;
   int i,j;
   
   
   /* Make sure the filter kernel has same parity as oversample */
   if ((fwidth ^ supersample) & 1)
      fwidth++;

   /* Calculate the coordinate scaling factor for the kernel values */
   if (fw > 0.0)
      adjust = sf_supp * fwidth / fw;
   else
      adjust = 1.0;

   /* Calling function MUST FREE THE RETURNED KERNEL, lest ye leak memory */
   (*filter) = (double *)calloc(fwidth * fwidth,sizeof(double));

   /* fill in the coefs */
   for (i = 0; i < fwidth; i++)
      for (j = 0; j < fwidth; j++) {
      
         /* Calculate the function inputs for the kernel function */
         ii = ((2.0 * i + 1.0) / (double)fwidth - 1.0)*adjust;
         jj = ((2.0 * j + 1.0) / (double)fwidth - 1.0)*adjust;

         /* Scale for scanlines */
         if (field) jj *= 2.0;

         /* Adjust for aspect ratio */
         jj /= aspect_ratio;

         (*filter)[i + j * fwidth] = 
               flam3_spatial_filter(sf_kernel,ii) * flam3_spatial_filter(sf_kernel,jj);
      }


   if (normalize_vector((*filter), fwidth * fwidth)) {
      fprintf(stderr, "Spatial filter value is too small: %g.  Terminating.\n",sf_radius);
      return(-1);
   }   
   
   return (fwidth);
}

flam3_de_helper flam3_create_de_filters(double max_rad, double min_rad, double curve, int ss) {

   flam3_de_helper de;
   double comp_max_radius, comp_min_radius;
   double num_de_filters_d;
   int num_de_filters,de_max_ind;
   int de_row_size, de_half_size;
   int filtloop;
   int keep_thresh=100;

   de.kernel_size=-1;

   if (curve <= 0.0) {
      fprintf(stderr,"estimator curve must be > 0\n");
      return(de);
   }

   if (max_rad < min_rad) {
      fprintf(stderr,"estimator must be larger than estimator_minimum.\n");
      fprintf(stderr,"(%f > %f) ? \n",max_rad,min_rad);
      return(de);
   }

   /* We should scale the filter width by the oversample          */
   /* The '+1' comes from the assumed distance to the first pixel */
   comp_max_radius = max_rad*ss + 1;
   comp_min_radius = min_rad*ss + 1;

   /* Calculate how many filter kernels we need based on the decay function */
   /*                                                                       */
   /*    num filters = (de_max_width / de_min_width)^(1/estimator_curve)    */
   /*                                                                       */
   num_de_filters_d = pow( comp_max_radius/comp_min_radius, 1.0/curve );
   if (num_de_filters_d>1e7) {
      fprintf(stderr,"too many filters required in this configuration (%g)\n",num_de_filters_d);
      return(de);
   }
   num_de_filters = (int)ceil(num_de_filters_d);
         
   /* Condense the smaller kernels to save space */
   if (num_de_filters>keep_thresh) { 
      de_max_ind = (int)ceil(DE_THRESH + pow(num_de_filters-DE_THRESH,curve))+1;
      de.max_filtered_counts = (int)pow( (double)(de_max_ind-DE_THRESH), 1.0/curve) + DE_THRESH;
   } else {
      de_max_ind = num_de_filters;
      de.max_filtered_counts = de_max_ind;
   }

   /* Allocate the memory for these filters */
   /* and the hit/width lookup vector       */
   de_row_size = (int)(2*ceil(comp_max_radius)-1);
   de_half_size = (de_row_size-1)/2;
   de.kernel_size = (de_half_size+1)*(2+de_half_size)/2;

   de.filter_coefs = (double *) calloc (de_max_ind * de.kernel_size,sizeof(double));
   de.filter_widths = (double *) calloc (de_max_ind,sizeof(double));

   /* Generate the filter coefficients */
   de.max_filter_index = 0;
   for (filtloop=0;filtloop<de_max_ind;filtloop++) {

      double de_filt_sum=0.0, de_filt_d;
      double de_filt_h;
      int dej,dek;
      double adjloop;
      int filter_coef_idx;

      /* Calculate the filter width for this number of hits in a bin */
      if (filtloop<keep_thresh)
         de_filt_h = (comp_max_radius / pow(filtloop+1,curve));
      else {
         adjloop = pow(filtloop-keep_thresh,(1.0/curve)) + keep_thresh;
         de_filt_h = (comp_max_radius / pow(adjloop+1,curve));
      }

      /* Once we've reached the min radius, don't populate any more */
      if (de_filt_h <= comp_min_radius) {
         de_filt_h = comp_min_radius;
         de.max_filter_index = filtloop;
      }

      de.filter_widths[filtloop] = de_filt_h;

      /* Calculate norm of kernel separately (easier) */
      for (dej=-de_half_size; dej<=de_half_size; dej++) {
         for (dek=-de_half_size; dek<=de_half_size; dek++) {
            
            de_filt_d = sqrt( (double)(dej*dej+dek*dek) ) / de_filt_h;

            /* Only populate the coefs within this radius */
            if (de_filt_d <= 1.0) {

               /* Gaussian */
               de_filt_sum += flam3_spatial_filter(flam3_gaussian_kernel,
                        flam3_spatial_support[flam3_gaussian_kernel]*de_filt_d);

               /* Epanichnikov */
//             de_filt_sum += (1.0 - (de_filt_d * de_filt_d));
            }
         }
      }

      filter_coef_idx = filtloop*de.kernel_size;

      /* Calculate the unique entries of the kernel */
      for (dej=0; dej<=de_half_size; dej++) {
         for (dek=0; dek<=dej; dek++) {
            de_filt_d = sqrt( (double)(dej*dej+dek*dek) ) / de_filt_h;

            /* Only populate the coefs within this radius */
            if (de_filt_d>1.0)
               de.filter_coefs[filter_coef_idx] = 0.0;
            else {

               /* Gaussian */
               de.filter_coefs[filter_coef_idx] = flam3_spatial_filter(flam3_gaussian_kernel,
                        flam3_spatial_support[flam3_gaussian_kernel]*de_filt_d)/de_filt_sum; 
                            
               /* Epanichnikov */
//             de_filter_coefs[filter_coef_idx] = (1.0 - (de_filt_d * de_filt_d))/de_filt_sum;
            }
                  
            filter_coef_idx ++;
         }
      }

      if (de.max_filter_index>0)
         break;
   }

   if (de.max_filter_index==0)
      de.max_filter_index = de_max_ind-1;

   
   return(de);
}                       

double flam3_create_temporal_filter(int numsteps, int filter_type, double filter_exp, double filter_width,
                                    double **temporal_filter, double **temporal_deltas) {

   double maxfilt = 0.0;
   double sumfilt = 0.0;
   double slpx,halfsteps;
   double *deltas, *filter;
   
   int i;

   /* Allocate memory - this must be freed in the calling routine! */   
   deltas = (double *)malloc(numsteps*sizeof(double));
   filter = (double *)malloc(numsteps*sizeof(double));
   
   /* Deal with only one step */
   if (numsteps==1) {
      deltas[0] = 0;
      filter[0] = 1.0;
      *temporal_deltas = deltas;
      *temporal_filter = filter;
      return(1.0);
   }
      
   /* Define the temporal deltas */   
   for (i = 0; i < numsteps; i++)
      deltas[i] = ((double)i /(double)(numsteps - 1) - 0.5)*filter_width;
      
   /* Define the filter coefs */
   if (flam3_temporal_exp == filter_type) {

      for (i=0; i < numsteps; i++) {

         if (filter_exp>=0)
            slpx = ((double)i+1.0)/numsteps;
         else
            slpx = (double)(numsteps - i)/numsteps;

         /* Scale the color based on these values */
         filter[i] = pow(slpx,fabs(filter_exp));
         
         /* Keep the max */
         if (filter[i]>maxfilt)
            maxfilt = filter[i];
      }

   } else if (flam3_temporal_gaussian == filter_type) {

      halfsteps = numsteps/2.0;
      for (i=0; i < numsteps; i++) {
      
         /* Gaussian */
         filter[i] = flam3_spatial_filter(flam3_gaussian_kernel,
                           flam3_spatial_support[flam3_gaussian_kernel]*fabs(i - halfsteps)/halfsteps);
         /* Keep the max */
         if (filter[i]>maxfilt)
            maxfilt = filter[i];
      }
      
   } else { // (flam3_temporal_box)

      for (i=0; i < numsteps; i++)
         filter[i] = 1.0;
         
	   maxfilt = 1.0;
	   
   }

   /* Adjust the filter so that the max is 1.0, and */
   /* calculate the K2 scaling factor  */
   for (i=0;i<numsteps;i++) {
      filter[i] /= maxfilt;
      sumfilt += filter[i];
   }
         
   sumfilt /= numsteps;
   
   *temporal_deltas = deltas;
   *temporal_filter = filter;
   
   return(sumfilt);
}                                     
 
