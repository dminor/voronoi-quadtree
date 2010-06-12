/*
Copyright (c) 2010 Daniel Minor 

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <shapefil.h>

#include "voronoi_quadtree.h"

void render_point(FILE *f, int site, double x, double y)
{
    fprintf(f, "colour-site-%d\n", site); 
    fprintf(f, "%.0f %.0f 2 0 360 arc\n", x, y);
    fprintf(f, "fill\n"); 
}

void render_line(FILE *f, int site, SHPObject *obj)
{ 
    fprintf(f, "colour-site-%d\n", site); 
    fprintf(f, "newpath\n");
    fprintf(f, "%.0f %.0f moveto\n", obj->padfX[0], obj->padfY[0]);

    for (int j = 0; j < obj->nVertices; ++j) { 
        fprintf(f, "%.0f %.0f lineto\n", obj->padfX[j], obj->padfY[j]);
    }

    fprintf(f, "stroke\n"); 
}

void render_poly(FILE *f, int site, SHPObject *obj)
{ 
    fprintf(f, "colour-site-%d\n", site); 
    fprintf(f, "newpath\n");
    fprintf(f, "%.0f %.0f moveto\n", obj->padfX[0], obj->padfY[0]);

    for (int j = 0; j < obj->nVertices; ++j) { 
        fprintf(f, "%.0f %.0f lineto\n", obj->padfX[j], obj->padfY[j]);
    }

    fprintf(f, "closepath\n");
    fprintf(f, "fill\n"); 
}

void render_voronoi_quadtree(FILE *f, struct VoronoiQuadtreeNode *node)
{
    if (node->ne == 0 && node->nw == 0 && node->sw == 0 && node->se == 0) {
        fprintf(f, "colour-site-%d\n", node->site_id);
        fprintf(f, "%.1f %.1f %.1f %.1f box\n", node->x1, node->x2, node->y1, node->y2);
    } else {
        render_voronoi_quadtree(f, node->ne);
        render_voronoi_quadtree(f, node->nw);
        render_voronoi_quadtree(f, node->sw);
        render_voronoi_quadtree(f, node->se);
    }
}

int main(int argc, char **argv)
{ 
    //handle command line
    if (argc < 2) {
        printf("usage: render-quadtree <filename> [maximum depth]\n");
        exit(1);
    }

    int max_depth = 8;
    if (argc == 3) max_depth = atoi(argv[2]);

    //read data from shapefile
    SHPHandle shp = SHPOpen(argv[1], "rb");
    if (!shp) {
        printf("could not open shapefile: %s\n", argv[1]);
        exit(1);
    }

    int entcount;
    int type;
    SHPGetInfo(shp, &entcount, &type, 0, 0); 

    if (entcount <= 0) {
        printf("error: shapefile is empty\n");
        exit(1); 
    }

    //setup output file
    FILE *f = fopen("voronoi-diagram.ps", "w");
    fprintf(f, "%\n");

    //define box function for later
    fprintf(f, "/box {\n");
    fprintf(f, "    /y2 exch def\n");
    fprintf(f, "    /y1 exch def\n");
    fprintf(f, "    /x2 exch def\n");
    fprintf(f, "    /x1 exch def\n");
    fprintf(f, "    gsave\n");
    fprintf(f, "    newpath\n");
    fprintf(f, "    x1 y1 moveto\n");
    fprintf(f, "    x1 y2 lineto\n");
    fprintf(f, "    x2 y2 lineto\n");
    fprintf(f, "    x2 y1 lineto\n");
    fprintf(f, "    closepath\n");
    fprintf(f, "    stroke\n");
    fprintf(f, "    grestore\n");
    fprintf(f, "} def\n");

    //setup colours
    for (int i = 0; i < entcount; ++i) {
        fprintf(f, "/colour-site-%d {%.1f %.1f %.1f setrgbcolor } def\n", i, (double)rand()/(double)RAND_MAX, (double)rand()/(double)RAND_MAX,(double)rand()/(double)RAND_MAX);
    }

    struct VoronoiQuadtreeNode *node = 0;
    if (type == SHPT_POINT) {

        struct Site *sites = malloc(entcount * sizeof(struct Site));

        for (int i = 0; i < entcount; ++i) {
            SHPObject *obj = SHPReadObject(shp, i);
            render_point(f, i, obj->padfX[0], obj->padfY[0]);

            sites[i].n = 1;
            sites[i].xs = malloc(sizeof(double));
            sites[i].xs[0] = obj->padfX[0];
            sites[i].ys = malloc(sizeof(double)); 
            sites[i].ys[0] = obj->padfY[0];
            sites[i].metric = pt_euclidean_metric;

            SHPDestroyObject(obj);
        }

        node = build_voronoi_quadtree(entcount, sites, max_depth);

        for (int i = 0; i < entcount; ++i) {
            free(sites[i].xs);
            free(sites[i].ys);
        }

        free(sites);

    } else if (type == SHPT_ARC || type == SHPT_POLYGON) {

        struct Site *sites = malloc(entcount * sizeof(struct Site));

        for (int i = 0; i < entcount; ++i) {
            SHPObject *obj = SHPReadObject(shp, i);
            if (type == SHPT_ARC) render_line(f, i, obj);
            else render_poly(f, i, obj);

            sites[i].n = obj->nVertices;
            sites[i].xs = malloc(sites[i].n * sizeof(double));
            sites[i].ys = malloc(sites[i].n * sizeof(double));

            for (int j = 0; j < obj->nVertices; ++j) { 
                sites[i].xs[j] = obj->padfX[j];
                sites[i].ys[j] = obj->padfY[j];
            }

            if (type == SHPT_ARC) sites[i].metric = line_euclidean_metric;
            else sites[i].metric = poly_euclidean_metric;

            SHPDestroyObject(obj); 
        }

        node = build_voronoi_quadtree(entcount, sites, max_depth);

        for (int i = 0; i < entcount; ++i) {
            free(sites[i].xs);
            free(sites[i].ys);
        }

        free(sites);


    } else { 
        printf("error: shapefile does not contain point, line or areal data\n");
        exit(1); 
    }

    SHPClose(shp);

    if (node) render_voronoi_quadtree(f, node);
    fclose(f);

    return 0;
}
