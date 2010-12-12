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

#include "private.h"
#include "isaacs.h"
#include "config.h"

int verbose;

char *get_extras() {
  char *e = getenv("extras");
  return e;
}

void gprint(flam3_genome *cp, int extras) {
    if (getenv("noedits"))
        flam3_print(stdout, cp, extras ? get_extras() : NULL, flam3_dont_print_edits);
    else
        flam3_print(stdout, cp, extras ? get_extras() : NULL, flam3_print_edits);
}


void test_cp(flam3_genome *cp) {
   cp->time = 0.0;
   cp->interpolation = flam3_interpolation_linear;
   cp->palette_interpolation = flam3_palette_interpolation_hsv;
   cp->background[0] = 0.0;
   cp->background[1] = 0.0;
   cp->background[2] = 0.0;
   cp->center[0] = 0.0;
   cp->center[1] = 0.0;
   cp->rotate = 0.0;
   cp->pixels_per_unit = 64;
   cp->width = 128;
   cp->height = 128;
   cp->spatial_oversample = 1;
   cp->spatial_filter_radius = 0.5;
   cp->spatial_filter_select = 0;
   cp->zoom = 0.0;
   cp->sample_density = 1;
   cp->nbatches = 1;
   cp->ntemporal_samples = 1;
   cp->estimator = 0.0;
   cp->estimator_minimum = 0.0;
   cp->estimator_curve = 0.6;
}

flam3_genome *string_to_cp(char *s, int *n) {
  flam3_genome *cp;
  FILE *fp;

  fp = fopen(s, "rb");
  if (NULL == fp) {
    perror(s);
    exit(1);
  }
  cp = flam3_parse_from_file(fp, s, flam3_defaults_on, n);
  if (NULL == cp) {
      fprintf(stderr, "could not read genome from %s.\n", s);
      exit(1);
  }
  return cp;
}

xmlDocPtr create_new_editdoc(char *action, flam3_genome *parent0, flam3_genome *parent1) {

   xmlDocPtr doc = NULL, comment_doc = NULL;
   xmlNodePtr root_node = NULL, node = NULL, nodecopy = NULL;
   xmlNodePtr root_comment = NULL;
   struct tm *localt;
   time_t mytime;
   char *ai;
   char timestring[100];
   char *nick = getenv("nick");
   char *url = getenv("url");
   char *id = getenv("id");
   char *comment = getenv("comment");
   int sheep_gen = argi("sheep_gen",-1);
   int sheep_id = argi("sheep_id",-1);
   char buffer[100];
   char comment_string[100];

   doc = xmlNewDoc( (const xmlChar *)"1.0");

   /* Create the root node, called "edit" */
   root_node = xmlNewNode(NULL, (const xmlChar *)"edit");
   xmlDocSetRootElement(doc,root_node);
   /* Add the edit attributes */

   /* date */
   mytime = time(NULL);
   localt = localtime(&mytime);
   /* XXX use standard time format including timezone */
   strftime(timestring, 100, "%a %b %e %H:%M:%S %z %Y", localt);
   xmlNewProp(root_node, (const xmlChar *)"date", (const xmlChar *)timestring);

   /* nick */
   if (nick) {
      xmlNewProp(root_node, (const xmlChar *)"nick", (const xmlChar *)nick);
   }

   /* url */
   if (url) {
      xmlNewProp(root_node, (const xmlChar *)"url", (const xmlChar *)url);
   }

   if (id) {
      xmlNewProp(root_node, (const xmlChar *)"id", (const xmlChar *)id);
   }

   /* action */
   xmlNewProp(root_node, (const xmlChar *)"action", (const xmlChar *)action);

   /* sheep info */
   if (sheep_gen > 0 && sheep_id > 0) {
      /* Create a child node of the root node called sheep */
      node = xmlNewChild(root_node, NULL, (const xmlChar *)"sheep", NULL);

      /* Create the sheep attributes */
      sprintf(buffer, "%d", sheep_gen);
      xmlNewProp(node, (const xmlChar *)"generation", (const xmlChar *)buffer);

      sprintf(buffer, "%d", sheep_id);
      xmlNewProp(node, (const xmlChar *)"id", (const xmlChar *)buffer);
   }

   /* Check for the parents */
   /* If Parent 0 not specified, this is a randomly generated genome. */
   if (parent0) {
      if (parent0->edits) {
         /* Copy the node from the parent */
         node = xmlDocGetRootElement(parent0->edits);
         nodecopy = xmlCopyNode(node, 1);
         xmlNewProp(nodecopy,(const xmlChar *)"filename", (const xmlChar *)parent0->parent_fname);
         sprintf(buffer,"%d",parent0->genome_index);
         xmlNewProp(nodecopy,(const xmlChar *)"index", (const xmlChar *)buffer);
         xmlAddChild(root_node, nodecopy);
      } else {
         /* Insert a (parent has no edit) message */
         nodecopy = xmlNewChild(root_node, NULL, (const xmlChar *)"edit",NULL);
         xmlNewProp(nodecopy,(const xmlChar *)"filename", (const xmlChar *)parent0->parent_fname);
         sprintf(buffer,"%d",parent0->genome_index);
         xmlNewProp(nodecopy,(const xmlChar *)"index", (const xmlChar *)buffer);

      }
   }

   if (parent1) {

      if (parent1->edits) {
         /* Copy the node from the parent */
         node = xmlDocGetRootElement(parent1->edits);
         nodecopy = xmlCopyNode(node, 1);
         xmlNewProp(nodecopy,(const xmlChar *)"filename", (const xmlChar *)parent1->parent_fname);
         sprintf(buffer,"%d",parent1->genome_index);
         xmlNewProp(nodecopy,(const xmlChar *)"index", (const xmlChar *)buffer);
         xmlAddChild(root_node, nodecopy);
      } else {
         /* Insert a (parent has no edit) message */
         nodecopy = xmlNewChild(root_node, NULL, (const xmlChar *)"edit",NULL);
         xmlNewProp(nodecopy,(const xmlChar *)"filename", (const xmlChar *)parent1->parent_fname);
         sprintf(buffer,"%d",parent1->genome_index);
         xmlNewProp(nodecopy,(const xmlChar *)"index", (const xmlChar *)buffer);
      }
   }

   /* Comment string */
   /* This one's hard, since we have to treat the comment string as   */
   /* a valid XML document.  Create a new document using the comment  */
   /* string as the in-memory document, and then copy all children of */
   /* the root node into the edit structure                           */
   /* Parsing the comment string should be done once and then copied  */
   /* for each call to create_new_editdoc, but that's for later.      */
   if (comment) {

      sprintf(comment_string,"<comm>%s</comm>",comment);

      comment_doc = xmlReadMemory(comment_string, strlen(comment_string), "comment.env", NULL, XML_PARSE_NONET);

      /* Check for errors */
      if (comment_doc==NULL) {
         fprintf(stderr, "Failed to parse comment into XML!\n");
         exit(1);
      }

      /* Loop through the children of the new document and copy */
      /* them into the root_node */
      root_comment = xmlDocGetRootElement(comment_doc);

      for (node=root_comment->children; node; node = node->next) {

         nodecopy = xmlCopyNode(node,1);
         xmlAddChild(root_node, nodecopy);
      }

      /* Free the created document */
      xmlFreeDoc(comment_doc);
   }


   /* return the xml doc */
   return(doc);
}

void offset(flam3_genome *g) {
    char *os = getenv("offset");
    double ox, oy;
    if (NULL == os) return;
    sscanf(os, "%lf:%lf", &ox, &oy);
    g->center[0] += ox / (g->pixels_per_unit * g->spatial_oversample);
    g->center[1] += oy / (g->pixels_per_unit * g->spatial_oversample);
}

void spin(int frame, double blend, flam3_genome *parent, flam3_genome *templ)
{
   flam3_genome *result;
   char action[50];
   xmlDocPtr doc;
   
   /* Spin the parent blend*360 degrees */
   result = sheep_loop(parent,blend);

   /* Apply the template if necessary */  
   if (templ)
      flam3_apply_template(result, templ);

   /* Set genome parameters accordingly */
   result->time = (double)frame;
   result->interpolation = flam3_interpolation_linear;
   result->palette_interpolation = flam3_palette_interpolation_hsv;

   /* Force linear interpolation - unsure if this is still necessary     */
   /* I believe we put this in so that older clients could render frames */
//   result->interpolation_type = flam3_inttype_linear;

   /* Create the edit doc xml */
   sprintf(action,"rotate %g",blend*360.0);
   doc = create_new_editdoc(action, parent, (flam3_genome *)NULL);
   result->edits = doc;

   /* Subpixel jitter */
   offset(result);

   /* Make the name of the flame the time */
   sprintf(result->flame_name,"%f",result->time);

   /* Print the resulting xml */
   gprint(result, 1);

   /* Clear out the xml doc */
   xmlFreeDoc(result->edits);

   /* Clear the result cp */
   clear_cp(result,flam3_defaults_on);
 
   /* Free the cp allocated in flam3_sheep_loop */ 
   free(result);
}

void spin_inter(int frame, double blend, int seqflag, flam3_genome *parents, flam3_genome *templ) {

   flam3_genome *result;
   char action[50];
   xmlDocPtr doc;
   char *ai;
   double stagger = argf("stagger", 0.0);

   /* Interpolate between spun parents */
   result = sheep_edge(parents, blend, seqflag, stagger);

   /* Unsure why we check for random palettes on both ends... */  
   if ((parents[0].palette_index != flam3_palette_random) &&
      (parents[1].palette_index != flam3_palette_random)) {

         result->palette_index = flam3_palette_interpolated;
         result->palette_index0 = parents[0].palette_index;
         result->hue_rotation0 = parents[0].hue_rotation;
         result->palette_index1 = parents[1].palette_index;
         result->hue_rotation1 = parents[1].hue_rotation;
         result->palette_blend = blend;
   }

   /* Apply template if necessary */
   if (templ)
      flam3_apply_template(result, templ);

   /* Set genome attributes */
   result->time = (double)frame;
//   result->interpolation_type = flam3_inttype_linear;

   /* Create the edit doc xml */
   sprintf(action,"interpolate %g",blend*360.0);
   doc = create_new_editdoc(action, &parents[0], &parents[1]);
   result->edits = doc;

   /* Subpixel jitter */
   offset(result);

   /* Make the name of the flame the time */
   sprintf(result->flame_name,"%f",result->time);

   /* Print the genome */  
   gprint(result, 1);

   /* Clean up */
   xmlFreeDoc(result->edits);

   /* Free genome storage */
   clear_cp(result,flam3_defaults_on);
   free(result);
}

void truncate_variations(flam3_genome *g, int max_vars, char *action) {
   int i, j, nvars, smallest;
   double sv=0;
   char trunc_note[30];

   for (i = 0; i < g->num_xforms; i++) {
      double d = g->xform[i].density;

/*      if (0.0 < d && d < 0.001) */

      if (d < 0.001 && (g->final_xform_index != i)) {
         sprintf(trunc_note," trunc_density %d",i);
         //strcat(action,trunc_note);
         add_to_action(action,trunc_note);
         flam3_delete_xform(g, i);

/*         g->xform[i].density = 0.0;
      } else if (d > 0.0) {
*/
      } else {
         do {
            nvars = 0;
            smallest = -1;
            for (j = 0; j < flam3_nvariations; j++) {
               double v = g->xform[i].var[j];
               if (v != 0.0) {
                  nvars++;
                  if (-1 == smallest || fabs(v) < sv) {
                     smallest = j;
                     sv = fabs(v);
                  }
               }
            }
            if (nvars > max_vars) {
               sprintf(trunc_note," trunc %d %d",i,smallest);
               //strcat(action,trunc_note);
               add_to_action(action,trunc_note);
               g->xform[i].var[smallest] = 0.0;
            }
         } while (nvars > max_vars);
      }
   }
}

static double golden_bit(randctx *rc) {
  return flam3_random_isaac_bit(rc)?0.38196:0.61804;
}

int
main(argc, argv)
   int argc;
   char **argv;
{
   int debug = 0;
   int count;
   char *ai;
   unsigned char *image;
   flam3_genome *templ = NULL;
   flam3_genome cp_orig, cp_save;
   int i, j;
   double avg_pix, fraction_black, fraction_white;
   double avg_thresh = argf("avg", 20.0);
   double black_thresh = argf("black", 0.01);
   double white_limit = argf("white", 0.05);
   int nframes = argi("nframes", 100);
   int sym = argi("symmetry", 0);
   int enclosed = argi("enclosed", 1);
   char *clone = getenv("clone");
   char *clone_all = getenv("clone_all");
   char *animate = getenv("animate");
   char *mutate = getenv("mutate");
   char *cross0 = getenv("cross0");
   char *cross1 = getenv("cross1");
   char *method = getenv("method");
   char *inter = getenv("inter");
   char *rotate = getenv("rotate");
   char *strip = getenv("strip");
   char *sequence = getenv("sequence");
   int loops = argi("loops", 1);
   int frame = argi("frame", 0);
   int rep, repeat = argi("repeat", 1);
   double speed = argf("speed", 0.1);
   int bits = argi("bits", 33);
   int ntries = argi("tries", 10);
   char *use_vars = getenv("use_vars");
   char *dont_use_vars = getenv("dont_use_vars");
   flam3_genome *parent0=NULL, *parent1=NULL;
   flam3_genome selp0, selp1;
   flam3_genome *aselp0, *aselp1;
   int parent0_n, parent1_n;
   int num_threads = 1;
   flam3_genome *cp;
   int ncp;

   int ivars[max_specified_vars];
   int novars[max_specified_vars];
   int num_ivars = 0;
   int num_novars = 0;
   char *var_tok;

   flam3_frame f;
   char action[flam3_max_action_length];  
   
   stat_struct stats;



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

   if (argc>1) {
      if (strcmp("--version",argv[1])==0) {
         printf("FLAM3-%s\n",flam3_version());
         exit(0);
      } else {
         printf("unrecognized option %s, aborting.\n",argv[1]);
         exit(-1);
      }
   }

   verbose = argi("verbose", 0);

   memset(&cp_orig, 0, sizeof(flam3_genome));
   memset(&cp_save, 0, sizeof(flam3_genome));
   memset(&selp0, 0, sizeof(flam3_genome));
   memset(&selp1, 0, sizeof(flam3_genome));

   if (1 != argc) {
      docstring();
      exit(0);
   }

   /* Init random number generators */
   flam3_init_frame(&f);
   flam3_srandom();

   //f.temporal_filter_radius = 0.0;
   f.bits = bits;
   f.bytes_per_channel = 1;
   f.earlyclip = 1;
   f.verbose = 0;
   f.genomes = &cp_orig;
   f.ngenomes = 1;
   f.pixel_aspect_ratio = 1.0;
   f.progress = 0;
   f.nthreads = num_threads;
   f.sub_batch_size = 10000;
   test_cp(&cp_orig);  // just for the width & height

   /* Are the variations to be used specified? */
   if (use_vars && dont_use_vars) {
      fprintf(stderr,"use_vars and dont_use_vars cannot both be specified.  Terminating.\n");
      exit(-1);
   }
   
   /* Specify reasonable defaults if nothing is specified */
   if (!use_vars && !dont_use_vars) {
      novars[num_novars++] = VAR_NOISE;
      novars[num_novars++] = VAR_BLUR;
      novars[num_novars++] = VAR_GAUSSIAN_BLUR;
      novars[num_novars++] = VAR_RADIAL_BLUR;
      novars[num_novars++] = VAR_NGON;
      novars[num_novars++] = VAR_SQUARE;
      novars[num_novars++] = VAR_RAYS;
      novars[num_novars++] = VAR_CROSS;
      novars[num_novars++] = VAR_PRE_BLUR;
      novars[num_novars++] = VAR_SEPARATION;
      novars[num_novars++] = VAR_SPLIT;
      novars[num_novars++] = VAR_SPLITS;
      
      

      /* Loop over the novars and set ivars to the complement */
      for (i=0;i<flam3_nvariations;i++) {
         for (j=0;j<num_novars;j++) {
            if (novars[j] == i)
               break;
         }
         if (j==num_novars)
            ivars[num_ivars++] = i;
      }

   } else {
   
      if (use_vars) {
         /* Parse comma-separated list of variations to use */
         var_tok = strtok(use_vars,",");
         ivars[num_ivars++] = atoi(var_tok);
         while(1) {
            var_tok = strtok(NULL,",");

            if (var_tok==NULL)
               break;

            ivars[num_ivars++] = atoi(var_tok);

            if (num_ivars==max_specified_vars) {
               fprintf(stderr,"Maximum number of user-specified variations exceeded.  Truncating.\n");
               break;
            }
         }

         /* Error checking */
         for (i=0;i<num_ivars;i++) {
            if (ivars[i]<0 || ivars[i]>=flam3_nvariations) {
               fprintf(stderr,"specified variation list includes bad value. (%d)\n",ivars[i]);
               exit(1);
            }
         }
      } else if (dont_use_vars) {
         /* Parse comma-separated list of variations NOT to use */
         var_tok = strtok(dont_use_vars,",");
         novars[num_novars++] = atoi(var_tok);
         while(1) {
            var_tok = strtok(NULL,",");

            if (var_tok==NULL)
               break;

            novars[num_novars++] = atoi(var_tok);

            if (num_novars==max_specified_vars) {
               fprintf(stderr,"Maximum number of user-specified variations exceeded.  Truncating.\n");
               break;
            }
         }
         
         /* Loop over the novars and set ivars to the complement */
         for (i=0;i<flam3_nvariations;i++) {
            for (j=0;j<num_novars;j++) {
               if (novars[j] == i)
                  break;
            }
            if (j==num_novars)
               ivars[num_ivars++] = i;
         }
      }
   }
   
   if (1 < (!!mutate + !!(cross0 || cross1) +
       !!inter + !!rotate + !!clone + !!strip )) {
      fprintf(stderr,
      "can only specify one of mutate, clone, cross, rotate, strip, or inter.\n");
      exit(1);
   }
   
   if ( (!cross0) ^ (!cross1) ) {
      fprintf(stderr, "must specify both crossover arguments.\n");
      exit(1);
   }

   if (method && (!cross0 && !mutate)) {
      fprintf(stderr, "cannot specify method unless doing crossover or mutate.\n");
      exit(1);
   }

   if (getenv("template")) {
      char *tf = getenv("template");

      templ = string_to_cp(tf, &ncp);
      if (1 < ncp) {
         fprintf(stderr, "more than one control point in template, "
            "ignoring all but first.\n");
      } else if (0 == ncp) {
         fprintf(stderr, "no control points in template.\n");
         exit(1);
      }
      
   }

   /* Methods for genetic manipulation begin here */

   if (clone_all) {

      cp = string_to_cp(clone_all, &ncp);

      printf("<clone_all version=\"FLAM3-%s\">\n", flam3_version());
      for (i = 0; i < ncp; i++) {
         if (templ) flam3_apply_template(&cp[i], templ);
         offset(&cp[i]);
         gprint(&cp[i], 1);
      }
      printf("</clone_all>\n");
      
      exit(0);
   }
   
   if (animate) {
      flam3_genome interpolated;
      int first_frame,last_frame;
      int ftime,iscp;
      double stagger = argf("stagger", 0.0);
      cp = string_to_cp(animate, &ncp);
      
      
      for (i = 0; i < ncp; i++) {
         if (i > 0 && cp[i].time <= cp[i-1].time) {
            fprintf(stderr, "error: control points must be sorted by time, but %g <= %g, index %d.\n",
            cp[i].time, cp[i-1].time, i);
            exit(1);
         }
         /* Strip out all motion elements here */
         for (j=0;j<cp[i].num_xforms;j++)
            flam3_delete_motion_elements(&cp[i].xform[j]);
                     
      }

      if (!getenv("begin"))
         first_frame = (int) cp[0].time;
      else
         first_frame = argi("begin",0);
            
      if (!getenv("end"))
         last_frame = (int) cp[ncp-1].time;
      else
         last_frame = argi("end",0);
      
      if (last_frame < first_frame) last_frame = first_frame;

      printf("<animate version=\"FLAM3-%s\">\n", flam3_version());

      for (ftime = first_frame; ftime <= last_frame; ftime += 1) {
         iscp=0;
         for (i=0;i<ncp;i++) {
            if ( ftime==cp[i].time ) {
               flam3_copy(&interpolated, &(cp[i]) );
               iscp=1;
            }
         }
         if (iscp==0) {
            flam3_interpolate(cp, ncp, (double)ftime, stagger, &interpolated);
            for (i=0;i<ncp;i++) {
               if ( ftime==cp[i].time-1 ) {
                  iscp=1;
               }
            }
            if (iscp==0)
               interpolated.interpolation_type = flam3_inttype_linear;
         }
         
         if (templ) flam3_apply_template(&interpolated, templ);
         gprint(&interpolated, 1);
      }
      printf("</animate>\n");
      exit(0);
   }


   if (sequence) {
      double blend, spread;
      int seqflag;
      int framecount;

      if (nframes <= 0) {
         fprintf(stderr, "nframes must be positive, not %d.\n", nframes);
         exit(1);
      }

      cp = string_to_cp(sequence, &ncp);

      if (enclosed) printf("<sequence version=\"FLAM3-%s\">\n", flam3_version());
      spread = 1.0/nframes;
      framecount = 0;
#if 1
      for (i = 0; i < ncp; i++) {
         if (loops) {
            for (frame = 0; frame < nframes; frame++) {
               blend = frame/(double)nframes;
               spin(framecount++, blend, &cp[i], templ);
            }
         }
         if (i < ncp-1)
        for (frame = 0; frame < nframes; frame++) {
           if (0==frame || nframes-1==frame)
              seqflag=1;
           else
              seqflag=0;
       blend = frame/(double)nframes;
       spin_inter(framecount++, blend, seqflag, &cp[i], templ);
        }
      }
      spin(framecount, 0.0, &cp[ncp-1], templ);
#else
      if (1) {
     flam3_genome res;
     memset(&res, 0, sizeof(flam3_genome));
     res.final_xform_index = -1;
     flam3_add_xforms(&res, cp[0].num_xforms, 0);

     if (ncp < 4) {
         fprintf(stderr, "catmull-rom requires 4 or more control points.\n");
         exit(1);
     }
     for (i = 0; i < ncp - 3; i++) {
         for (frame = 0; frame < nframes; frame++) {
        blend = frame/(double)nframes;
        interpolate_catmull_rom(cp+i, blend, &res);
        res.time = frame + i * nframes;
        gprint(&res, 0);
        fflush(stdout);
         }
     }
      }
#endif
      if (enclosed) printf("</sequence>\n");
      exit(0);
   }

   if (inter || rotate) {

      double blend, spread;
      char *fname = inter?inter:rotate;
      int ni;

      if (nframes <= 0) {
         fprintf(stderr, "nframes must be positive, not %d.\n", nframes);
         exit(1);
      }

      blend = frame/(double)nframes;
      spread = 1.0/nframes;

      cp = string_to_cp(fname, &ncp);

      if (enclosed) printf("<pick version=\"FLAM3-%s\">\n", flam3_version());
      if (rotate) {
         if (1 != ncp) {
            fprintf(stderr, "rotation requires one control point, not %d.\n", ncp);
            exit(1);
         }
         spin(frame - 1, blend - spread, cp, templ);
         spin(frame    , blend         , cp, templ);
         spin(frame + 1, blend + spread, cp, templ);
      } else {
         if (2 != ncp) {
            fprintf(stderr, "interpolation requires two control points, not %d.\n", ncp);
            exit(1);
         }
         spin_inter(frame - 1, blend - spread, 0, cp, templ);
         spin_inter(frame    , blend         , 0, cp, templ);
         spin_inter(frame + 1, blend + spread, 0, cp, templ);
      }
      if (enclosed) printf("</pick>\n");
      
      for (ni=0;ni<ncp;ni++) {
         xmlFreeDoc(cp[ni].edits);
         clear_cp(&cp[ni],flam3_defaults_on);
      }
      free(cp);
      
      exit(0);
   }

   if (strip) {

      cp = string_to_cp(strip, &ncp);

      if (enclosed) printf("<pick version=\"FLAM3-%s\">\n", flam3_version());

      for (i = 0; i < ncp; i++) {
         double old_center[2];

         /* Strip out motion elements */         
         for (j=0;j<cp[i].num_xforms;j++)
            flam3_delete_motion_elements(&cp[i].xform[j]);
         
         old_center[0] = cp[i].center[0];
         old_center[1] = cp[i].center[1];
         cp[i].height /= nframes;
         cp[i].center[1] = cp[i].center[1] - ((nframes - 1) * cp[i].height) /
            (2 * cp[i].pixels_per_unit * pow(2.0, cp[i].zoom));
         cp[i].center[1] += cp[i].height * frame / ( cp[i].pixels_per_unit * pow(2.0,cp[i].zoom) );
         rotate_by(cp[i].center, old_center, cp[i].rotate);

         if (templ) flam3_apply_template(&cp[i], templ);
         offset(&cp[i]);
         gprint(&cp[i], 1);
      }

      if (enclosed) printf("</pick>\n");
      exit(0);
   }

   /* pick a control point until it looks good enough */
   if (repeat <= 0) {
     fprintf(stderr, "repeat must be positive, not %d.\n", repeat);
     exit(1);
   }

   if (enclosed) printf("<pick version=\"FLAM3-%s\">\n", flam3_version());
   image = (unsigned char *) malloc(3 * cp_orig.width * cp_orig.height);

   for (rep = 0; rep < repeat; rep++) {
   
      if (verbose)
         fprintf(stderr, "flame = %d/%d..", rep+1, repeat);

      count = 0;

      if (clone) {

         parent0 = string_to_cp(clone, &parent0_n);
         /* Action is 'clone' with trunc_vars concat */
         sprintf(action,"clone");
         
         if (getenv("clone_action"))
            sprintf(action,"clone %s", getenv("clone_action"));

         flam3_copy(&selp0, &(parent0[random()%parent0_n]));
         flam3_copy(&cp_save, &selp0);
         aselp0 = &selp0;
         aselp1 = NULL;
         truncate_variations(&cp_save, 5, action);

         cp_save.edits = create_new_editdoc(action, aselp0, aselp1);

      } else {
         int did_color;

         do {
         
            int random_mode=0;
         
            if (verbose) fprintf(stderr, ".");
            did_color = 0;
            f.time = (double) 0.0;
            action[0] = 0;

            if (mutate) {
               int mutmeth;

               parent0 = string_to_cp(mutate, &parent0_n);
               flam3_copy(&selp0, &(parent0[((unsigned)irand(&f.rc))%parent0_n]));
               flam3_copy(&cp_orig, &selp0);
               aselp0 = &selp0;
               aselp1 = NULL;

               if (NULL == getenv("method"))
                  mutmeth = MUTATE_NOT_SPECIFIED;
               else if (!strcmp(method,"all_vars"))
                  mutmeth = MUTATE_ALL_VARIATIONS;
               else if (!strcmp(method,"one_xform"))
                  mutmeth = MUTATE_ONE_XFORM_COEFS;
               else if (!strcmp(method,"add_symmetry"))
                  mutmeth = MUTATE_ADD_SYMMETRY;
               else if (!strcmp(method,"post_xforms"))
                  mutmeth = MUTATE_POST_XFORMS;
               else if (!strcmp(method,"color_palette"))
                  mutmeth = MUTATE_COLOR_PALETTE;
               else if (!strcmp(method,"delete_xform"))
                  mutmeth = MUTATE_DELETE_XFORM;
               else if (!strcmp(method,"all_coefs"))
                  mutmeth = MUTATE_ALL_COEFS;
               else {
                  fprintf(stderr,"method '%s' not defined for mutate.  defaulting to random.\n",method);
                  mutmeth = MUTATE_NOT_SPECIFIED;
               }

               flam3_mutate(&cp_orig, mutmeth, ivars, num_ivars, sym, speed, &f.rc, action);
               
               /* Scan string returned for 'mutate color' */
               if ( strstr(action,"mutate color") )
                  did_color = 1;
                  
               if (cp_orig.flame_name[0]) {
                  char tm[flam3_name_len+1];
                  strncpy(tm, cp_orig.flame_name, flam3_name_len);
                  snprintf(cp_orig.flame_name, flam3_name_len, "mutation %d of %s", rep, tm);
               }

            } else if (cross0) {
               int i0, i1;
               int crossmeth;

               parent0 = string_to_cp(cross0, &parent0_n);
               parent1 = string_to_cp(cross1, &parent1_n);

               i0 = ((unsigned)irand(&f.rc))%parent0_n;
               i1 = ((unsigned)irand(&f.rc))%parent1_n;

               flam3_copy(&selp0, &(parent0[i0]));
               flam3_copy(&selp1, &(parent1[i1]));

               aselp0 = &selp0;
               aselp1 = &selp1;

               if (NULL == getenv("method"))
                  crossmeth = CROSS_NOT_SPECIFIED;
               else if (!strcmp(method,"union"))
                  crossmeth = CROSS_UNION;
               else if (!strcmp(method,"interpolate"))
                  crossmeth = CROSS_INTERPOLATE;
               else if (!strcmp(method,"alternate"))
                  crossmeth = CROSS_ALTERNATE;
               else {
                  fprintf(stderr,"method '%s' not defined for cross.  defaulting to random.\n",method);
                  crossmeth = CROSS_NOT_SPECIFIED;
               }

               flam3_cross(&parent0[i0], &parent1[i1], &cp_orig, crossmeth, &f.rc, action);

               if (parent0[i0].flame_name[0] || parent1[i1].flame_name[0]) {
                  snprintf(cp_orig.flame_name, flam3_name_len, "%d of %s x %s", 
                           rep, parent0[i0].flame_name, parent1[i1].flame_name);
               }

            } else {
               sprintf(action,"random");
               random_mode=1;
               flam3_random(&cp_orig, ivars, num_ivars, sym, 0);


               aselp0 = NULL;
               aselp1 = NULL;
            }

            /* Adjust bounding box half the time */
            if (flam3_random_bit(&f.rc) || random_mode) {
               double bmin[2], bmax[2];
               flam3_estimate_bounding_box(&cp_orig, 0.01, 100000, bmin, bmax, &f.rc);
               if (flam3_random_isaac_01(&f.rc) < 0.3) {
                  cp_orig.center[0] = (bmin[0] + bmax[0]) / 2.0;
                  cp_orig.center[1] = (bmin[1] + bmax[1]) / 2.0;
                  add_to_action(action," recentered");
               } else {
                  double mix0, mix1;
                  if (flam3_random_isaac_bit(&f.rc)) {
                     mix0 = golden_bit(&f.rc) + flam3_random_isaac_11(&f.rc)/5;
                     mix1 = golden_bit(&f.rc);
                     add_to_action(action," reframed0");
                  } else if (flam3_random_isaac_bit(&f.rc)) {
                     mix0 = golden_bit(&f.rc);
                     mix1 = golden_bit(&f.rc) + flam3_random_isaac_11(&f.rc)/5;
                     add_to_action(action," reframed1");
                  } else {
                     mix0 = golden_bit(&f.rc) + flam3_random_isaac_11(&f.rc)/5;
                     mix1 = golden_bit(&f.rc) + flam3_random_isaac_11(&f.rc)/5;
                     add_to_action(action," reframed2");
                  }
                  cp_orig.center[0] = mix0 * bmin[0] + (1-mix0)*bmax[0];
                  cp_orig.center[1] = mix1 * bmin[1] + (1-mix1)*bmax[1];
               }
               cp_orig.rot_center[0] = cp_orig.center[0];
               cp_orig.rot_center[1] = cp_orig.center[1];
               cp_orig.pixels_per_unit = cp_orig.width / (bmax[0] - bmin[0]);
            }

            truncate_variations(&cp_orig, 5, action);

            if (!did_color && random()&1) {
               if (debug)
                  fprintf(stderr,"improving colors...\n");

               flam3_improve_colors(&cp_orig, 100, 0, 10);
               //strcat(action," improved colors");
               add_to_action(action," improved colors");
            }

            cp_orig.edits = create_new_editdoc(action, aselp0, aselp1);
            flam3_copy(&cp_save, &cp_orig);
            test_cp(&cp_orig);
            if (flam3_render(&f, image, flam3_field_both, 3, 0, &stats)) { 
               fprintf(stderr,"error rendering test image: aborting.\n");
               exit(1);
            }

            if (1) {
               int n, tot, totb, totw;
               n = cp_orig.width * cp_orig.height;
               tot = 0;
               totb = 0;
               totw = 0;
               for (i = 0; i < 3*n; i+=3) {
               
                  tot += (image[i]+image[i+1]+image[i+2]);
                  if (0 == image[i] && 0 == image[i+1] && 0 == image[i+2]) totb++;
                  if (255 == image[i] && 255 == image[i+1] && 255 == image[i+2]) totw++;

                  // printf("%d ", image[i]);
               }

               avg_pix = (tot / (double)(3*n));
               fraction_black = totb / (double)n;
               fraction_white = totw / (double)n;

               if (debug)
                  fprintf(stderr,
                     "avg_pix=%g fraction_black=%g fraction_white=%g n=%g\n",
                     avg_pix, fraction_black, fraction_white, (double)n);

            } else {
               avg_pix = avg_thresh + 1.0;
               fraction_black = black_thresh + 1.0;
               fraction_white = white_limit - 1.0;
            }
            
            clear_cp(&cp_orig,flam3_defaults_on);

            count++;
         } while ((avg_pix < avg_thresh ||
                   fraction_black < black_thresh ||
                   fraction_white > white_limit) &&
                   count < ntries);

         if (ntries == count) {
            fprintf(stderr, "warning: reached maximum attempts, giving up.\n");
         }

      }

      if (templ)
         flam3_apply_template(&cp_save,templ);

      cp_save.time = (double)rep;

      if (1) {
         char *maxforms = getenv("maxforms");
         if (maxforms && atoi(maxforms)) {
            cp_save.symmetry = 0;
            while (cp_save.num_xforms > atoi(maxforms))
               flam3_delete_xform(&cp_save, cp_save.num_xforms - 1);
         }
      }

      gprint(&cp_save, 1);
      fflush(stdout);

      /* Free created documents */
      /* (Only free once, since the copy is a ptr to the original) */
      xmlFreeDoc(cp_save.edits);
      clear_cp(&cp_save,0);

      if (verbose) {
         fprintf(stderr, "\ndone.  action = %s\n", action);
      }

   }
   if (enclosed) printf("</pick>\n");
   free(image);
      
   return 0;
}
