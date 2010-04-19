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

#ifndef parser_included
#define parser_included

#include "private.h"

int flam3_atoi(char *nstr);
double flam3_atof(char *nstr);
int var2n(const char *s);
int flam3_parse_hexformat_colors(char *colstr, flam3_genome *cp, int numcolors, int chan);

void scan_for_flame_nodes(xmlNode *cur_node, char *parent_file, int default_flag, flam3_genome **all_cp, int *all_ncps);
int parse_flame_element(xmlNode *flame_node, flam3_genome *loc_current_cp);
int parse_xform_xml(xmlNode *chld_node,flam3_xform *this_xform, int *num_xaos, 
                    flam3_chaos_entry **xaos, int numstd, int motionxf);
void flam3_edit_print(FILE *f, xmlNodePtr editNode, int tabs, int formatting);
int flam3_interp_missing_colors(flam3_genome *cp);
#endif
