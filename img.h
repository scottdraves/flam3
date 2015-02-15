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


#include <stdio.h>
#include "flam3.h"

#define FLAM3_PNG_COM 8

#ifdef WIN32
   #define snprintf _snprintf
#endif

typedef struct {

   char *genome;
   char *badvals;
   char *numiters;
   char *rtime;

} flam3_img_comments;


void write_jpeg(FILE *file, unsigned char *image, int width, int height, flam3_img_comments *fpc);
void write_png(FILE *file, void *image, int width, int height, flam3_img_comments *fpc, int bpc);

/* returns RGBA pixel array or NULL on failure */
unsigned char *read_png(FILE *file, int *width, int *height);
unsigned char *read_jpeg(FILE *file, int *width, int *height);
