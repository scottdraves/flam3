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

#ifdef WIN32
#define WINVER 0x0500
#include <windows.h>
#endif

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#include <limits.h>


#include "private.h"
#include "img.h"
#include "isaacs.h"


int calc_nstrips(flam3_frame *spec) {
  double mem_required;
  double mem_available;
  int nstrips,ninc;
  char *testmalloc;
#ifdef WIN32
  MEMORYSTATUS stat;
  stat.dwLength = sizeof(stat);
  GlobalMemoryStatus(&stat); // may want to go to GlobalMemoryStatusEx eventually
  mem_available = (double)stat.dwTotalPhys;
//  fprintf(stderr,"%lu bytes free memory...\n",(size_t)stat.dwAvailPhys);
//  if (mem_available > 1e9) mem_available = 1e9;
#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
  mem_available =
      (double)sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGESIZE);
#elif defined __APPLE__
#ifdef __LP64__
long physmem;
size_t len = sizeof(physmem);
static int mib[2] = { CTL_HW, HW_MEMSIZE };
#else
unsigned int physmem;
size_t len = sizeof(physmem);
static int mib[2] = { CTL_HW, HW_PHYSMEM };
#endif
if (sysctl(mib, 2, &physmem, &len, NULL, 0) == 0 && len == sizeof(physmem)) {
   mem_available = (double )physmem;
} else {
   fprintf(stderr, "warning: unable to determine physical memory.n");
   mem_available = 2e9;
}
#else
  fprintf(stderr, "warning: unable to determine physical memory.\n");
  mem_available = 2e9;
#endif
#if 0
  fprintf(stderr,"available phyical memory is %lu\n",
	  (unsigned long)mem_available);
#endif
  mem_available *= 0.8;
  if (getenv("use_mem")) {
      mem_available = atof(getenv("use_mem"));
  }
  mem_required = flam3_render_memory_required(spec);
  if (mem_available >= mem_required) return 1;
  nstrips = (int) ceil(mem_required / mem_available);

if (0) {
  /* Attempt to malloc a strip, and if it fails, try adding additional strips */
  ninc=-1;
  testmalloc = NULL;
  while(NULL==testmalloc && ninc<3) {
     ninc++;
	  testmalloc = (char *)malloc((int)ceil(mem_required / (nstrips+ninc)));
  }
  if (NULL==testmalloc) {
     fprintf(stderr,"Unable to allocate memory for render.  Please close some running programs and try to render again.\n");
     exit(1);
  } else {
     free(testmalloc);
     nstrips = nstrips + ninc;
  }
}
  
  return nstrips;
}

int print_progress(void *foo, double fraction, int stage, double eta) {
    fprintf(stderr, "stage=%s progress=%g eta=%g\n", stage?"filtering":"chaos", fraction, eta);
  return 0;
}

int main(int argc, char **argv) {
   flam3_frame f;
   char *ai;
   flam3_genome *cps;
   int ncps;
   int i;
   void *image=NULL;
   FILE *fp;
   char fname[256];
   size_t this_size, last_size = -1;
   double imgmem;
   unsigned int strip;
   double center_y, center_base;
   unsigned int nstrips;
   randctx savectx;
   char *prefix = args("prefix", "");
   char *out = args("out", NULL);
   char *format = getenv("format");
   int verbose = argi("verbose", 1);
   int bits = argi("bits", 33);
   int bpc = argi("bpc",8);
   int seed = argi("seed", 0);
   int transparency = argi("transparency", 0);
   char *inf = getenv("in");
   double qs = argf("qs", 1.0);
   double ss = argf("ss", 1.0);
   double pixel_aspect = argf("pixel_aspect", 1.0);
   int sub_batch_size = argi("sub_batch_size",10000);
   int name_enable = argi("name_enable",0);
   int num_threads = argi("nthreads",0);
   int earlyclip = argi("earlyclip",0);
   FILE *in;
   double zoom_scale;
   unsigned int channels;
   long start_time = (long)time(0);
   flam3_img_comments fpc;
   stat_struct stats;
   char numiter_string[64];
   char badval_string[64];
   char rtime_string[64];

#ifdef WIN32
   
   char *slashloc;
   char exepath[256];
   char palpath[256];
   memset(exepath,0,256);
   memset(palpath,0,256); 

    slashloc = strrchr(argv[0],'\\');
	if (NULL==slashloc) {
	   sprintf(palpath,"flam3_palettes=flam3-palettes.xml");
	} else {
       strncpy(exepath,argv[0],slashloc-argv[0]+1);
	   sprintf(palpath,"flam3_palettes=%sflam3-palettes.xml",exepath);
	}
	putenv(palpath);

#endif         
   
   
   if (1 != argc) {
     docstring();
     exit(0);
   }

   /* Init random number generators */
   flam3_init_frame(&f);
   flam3_srandom();

   /* Set the number of threads */
   if (num_threads==0) {
      num_threads = flam3_count_nthreads();
      if (verbose > 1)
         fprintf(stderr,"Automatically detected %d core(s)...\n",num_threads);
   } else{
      if (verbose)
         fprintf(stderr,"Manually specified %d thread(s)...\n",num_threads);
   }


   if (NULL == format) format = "png";
   if (strcmp(format, "jpg") &&
       strcmp(format, "ppm") &&
       strcmp(format, "png")) {
       fprintf(stderr,
          "format must be either jpg, ppm, or png, not %s.\n",
          format);
       exit(1);
   }

   channels = strcmp(format, "png") ? 3 : 4;

   /* Check for 16-bit-per-channel processing */
   if ( (16 == bpc) && (strcmp(format,"png") != 0)) {
	fprintf(stderr,"Support for 16 bpc images is only present for the png format.\n");
	exit(1);
   } else if (bpc != 8 && bpc != 16) {
	fprintf(stderr,"Unexpected bpc specified (%d)\n",bpc);
	exit(1);
   }
   
   if (pixel_aspect <= 0.0) {
     fprintf(stderr, "pixel aspect ratio must be positive, not %g.\n",
        pixel_aspect);
     exit(1);
   }

   if (inf)
     in = fopen(inf, "rb");
   else
     in = stdin;
   if (NULL == in) {
     perror(inf);
     exit(1);
   }

   cps = flam3_parse_from_file(in, inf, flam3_defaults_on, &ncps);
   if (NULL == cps) {
     fprintf(stderr,"error reading genomes from file\n");
     exit(1);
   }
   
   if (inf)
      fclose(in);

   for (i = 0; i < ncps; i++) {
      /* Force ntemporal_samples to 1 for -render */
      cps[i].ntemporal_samples = 1;
      cps[i].sample_density *= qs;
      cps[i].height = (int)(cps[i].height * ss);
      cps[i].width = (int)(cps[i].width * ss);
      cps[i].pixels_per_unit *= ss;
      if (cps[i].height<=0 || cps[i].width<=0) {
         fprintf(stderr,"output image has dimension <=0, aborting.\n");
         exit(1);
      }
   }

   if (out && (ncps > 1)) {
      fprintf(stderr, "hqi-flame: warning: writing multiple images "
      "to one file.  all but last will be lost.\n");
   }


   for (i = 0; i < ncps; i++) {
      int real_height;

      if (verbose && ncps > 1) {
         fprintf(stderr, "flame = %d/%d ", i+1, ncps);
      }

//      f.temporal_filter_radius = 0.0;
      f.genomes = &cps[i];
      f.ngenomes = 1;
      f.verbose = verbose;
      f.bits = bits;
      f.time = 0.0;
      f.pixel_aspect_ratio = pixel_aspect;
      f.progress = 0;//print_progress;
      f.nthreads = num_threads;
      f.earlyclip = earlyclip;
      f.sub_batch_size = sub_batch_size;
      
      if (16==bpc)
         f.bytes_per_channel = 2;
      else
         f.bytes_per_channel = 1;
         

      if (getenv("nstrips")) {
         nstrips = atoi(getenv("nstrips"));
      } else {
         nstrips = calc_nstrips(&f);
      }

      if (nstrips > cps[i].height) {
         fprintf(stderr, "cannot have more strips than rows but %d>%d.\n",
         nstrips, cps[i].height);
         exit(1);
      }
      
      imgmem = (double)channels * (double)cps[i].width 
               * (double)cps[i].height * f.bytes_per_channel;
      
      if (imgmem > ULONG_MAX) {
         fprintf(stderr,"Image size > ULONG_MAX.  Aborting.\n");
         exit(1);
      }

      this_size = (size_t)channels * (size_t)cps[i].width 
                  * (size_t)cps[i].height * f.bytes_per_channel;
      if (this_size != last_size) {
         if (last_size != -1)
            free(image);
         last_size = this_size;
         image = (void *) calloc(this_size, sizeof(char));
         if (NULL==image) {
            fprintf(stderr,"Error allocating memory for image.  Aborting\n");
            exit(1);
         }
      } else {
         memset(image, 0, this_size);
      }

      cps[i].sample_density *= nstrips;
      real_height = cps[i].height;
      cps[i].height = (int) ceil(cps[i].height / (double) nstrips);
      center_y = cps[i].center[1];
      zoom_scale = pow(2.0, cps[i].zoom);
      center_base = center_y - ((nstrips - 1) * cps[i].height) /
      (2 * cps[i].pixels_per_unit * zoom_scale);
      
      /* Copy off random context to use for each strip */
      memcpy(&savectx, &f.rc, sizeof(randctx));

      for (strip = 0; strip < nstrips; strip++) {
         size_t ssoff = (size_t)cps[i].height * strip * cps[i].width * channels * f.bytes_per_channel;
         void *strip_start = image + ssoff;
         cps[i].center[1] = center_base + cps[i].height * (double) strip / (cps[i].pixels_per_unit * zoom_scale);
         
         if ((cps[i].height * (strip + 1)) > real_height) {
            int oh = cps[i].height;
            cps[i].height = real_height - oh * strip;
            cps[i].center[1] -=
            (oh - cps[i].height) * 0.5 /
            (cps[i].pixels_per_unit * zoom_scale);
         }
         
         /* Use the same random context for each strip */
         memcpy(&f.rc, &savectx, sizeof(randctx));

         if (verbose && nstrips > 1) {
            fprintf(stderr, "strip = %d/%d\n", strip+1, nstrips);
         }
         if (verbose && (1 == nstrips) && (ncps > 1)) {
            fprintf(stderr, "\n");
         }
         cps[i].ntemporal_samples = 1;
         if (flam3_render(&f, strip_start, flam3_field_both, channels, transparency, &stats)) {
            fprintf(stderr,"error rendering image: aborting.\n");
            exit(1);
         }

         if (NULL != out) {
            strcpy(fname,out);
         } else if (name_enable && cps[i].flame_name[0]>0) {
            sprintf(fname, "%s.%s",cps[i].flame_name,format);
         } else {
            sprintf(fname, "%s%05d.%s", prefix, i, format);
         }
         if (verbose) {
            fprintf(stderr, "writing %s...", fname);
         }
         fp = fopen(fname, "wb");
         if (NULL == fp) {
            perror(fname);
            exit(1);
         }

         /* Generate temp file with genome */
         fpc.genome = flam3_print_to_string(f.genomes);
         
         sprintf(badval_string,"%g",stats.badvals/(double)stats.num_iters);
         fpc.badvals = badval_string;
         sprintf(numiter_string,"%g",(double)stats.num_iters);
         fpc.numiters = numiter_string;
         sprintf(rtime_string,"%d",stats.render_seconds);
         fpc.rtime = rtime_string;

         if (!strcmp(format, "png")) {

             write_png(fp, image, cps[i].width, real_height, &fpc, f.bytes_per_channel);            
            
         } else if (!strcmp(format, "jpg")) {
                                      
             write_jpeg(fp, (unsigned char *)image, cps[i].width, real_height, &fpc);
            
         } else {
            fprintf(fp, "P6\n");
            fprintf(fp, "%d %d\n255\n", cps[i].width, real_height);
            if (this_size != fwrite((unsigned char *)image, 1, this_size, fp)) {
		       perror(fname);
	        }
         }
         /* Free string */
         free(fpc.genome);

         fclose(fp);
      }

      /* restore the cps values to their original values */
      cps[i].sample_density /= nstrips;
      cps[i].height = real_height;
      cps[i].center[1] = center_y;

      if (verbose) {
         fprintf(stderr, "done.\n");
      }
   }
   if (verbose && ((ncps > 1) || (nstrips > 1))) {
      long total_time = (long)time(0) - start_time;

      if (total_time > 100)
         fprintf(stderr, "total time = %.1f minutes\n", total_time / 60.0);
      else
         fprintf(stderr, "total time = %ld seconds\n", total_time);
   }
   
   for (i=0;i<ncps;i++) {
   
      xmlFreeDoc(cps[i].edits);
      clear_cp(&cps[i],0);
   
   }
   free(cps);
   
   free(image);
   return 0;
}
