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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <shapefil.h>

#include "voronoi_quadtree.h"

struct Site {

    SHPObject *obj;
    size_t id;

    Site() : obj(0), id(0)
    {
    }

    virtual ~Site()
    {
        if (obj) SHPDestroyObject(obj);
    } 

    void set_obj(SHPObject *obj)
    {
        this->obj = obj;
    }
};

double pt_distance_to_line(double x, double y, double a_x, double a_y, double b_x, double b_y)
{
    double vba_x = b_x - a_x;
    double vba_y = b_y - a_y;

    double vpta_x = x - a_x;
    double vpta_y = y - a_y;

    double denom = (vba_x*vba_x + vba_y*vba_y);
    if (denom == 0.0) {
        return (a_x-x)*(a_x-x) + (a_y-y)*(a_y-y);
    }

    double t = (vba_x*vpta_x + vba_y*vpta_y) / denom;
    if (t > 1.0) t = 1.0;
    if (t < 0.0) t = 0.0;

    double line_x = a_x + t * vba_x;
    double line_y = a_y + t * vba_y;

    return (line_x-x)*(line_x-x) + (line_y-y)*(line_y-y); 
} 

double metric(Site *site, double *pt)
{
    SHPObject *obj = site->obj;
    double distance = std::numeric_limits<double>::max(); 

    double &x = pt[0];
    double &y = pt[1];

    if (obj) {
        switch(obj->nSHPType) {

            case SHPT_POINT:
                distance = (x-obj->padfX[0])*(x-obj->padfX[0])
                    + (y-obj->padfY[0])*(y-obj->padfY[0]);
                break;

            case SHPT_ARC:
                for (int i = 0; i < obj->nVertices; ++i) {
                    double d = pt_distance_to_line(x, y,
                        obj->padfX[i], obj->padfY[i],
                        obj->padfX[i + 1], obj->padfY[i + 1]); 

                    if (d < distance) {
                        distance = d; 
                    }
                } 
                break;

            case SHPT_POLYGON:
                for (int i = 0; i < obj->nVertices; ++i) {

                    double d; 
                    if (i + 1 == obj->nVertices) {
                        d = pt_distance_to_line(x, y,
                            obj->padfX[i], obj->padfY[i],
                            obj->padfX[0], obj->padfY[0]); 
                    } else { 
                        d = pt_distance_to_line(x, y,
                            obj->padfX[i], obj->padfY[i],
                            obj->padfX[i + 1], obj->padfY[i + 1]); 
                    }

                    if (d < distance) { 
                        distance = d; 
                    }
                } 
                break; 
        }
    }

    return distance;
}

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

void render_voronoi_quadtree(FILE *f, VoronoiQuadtree<Site>::Node *node)
{
    if (node->nodes == 0) {
        fprintf(f, "colour-site-%d\n", node->site->id);
        fprintf(f, "%.1f %.1f %.1f %.1f box\n",
            node->mid[0] - node->radius, node->mid[0] + node->radius,
            node->mid[1] - node->radius, node->mid[1] + node->radius);
    } else {
        render_voronoi_quadtree(f, node->nodes[0]);
        render_voronoi_quadtree(f, node->nodes[1]);
        render_voronoi_quadtree(f, node->nodes[2]);
        render_voronoi_quadtree(f, node->nodes[3]);
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
    double min[4];
    double max[4];
    SHPGetInfo(shp, &entcount, &type, min, max); 

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
        fprintf(f, "/colour-site-%d {%.1f %.1f %.1f setrgbcolor } def\n", i,
            (double)rand()/(double)RAND_MAX,
            (double)rand()/(double)RAND_MAX,
            (double)rand()/(double)RAND_MAX);
    }

    Site *sites = new Site[entcount];

    for (int i = 0; i < entcount; ++i) {
        sites[i].id = i;
        sites[i].obj = SHPReadObject(shp, i);

        switch(type) {

            case SHPT_POINT:
                render_point(f, i,
                    sites[i].obj->padfX[0], sites[i].obj->padfY[0]);
                break;

            case SHPT_ARC:
                render_line(f, i, sites[i].obj);
                break;

            case SHPT_POLYGON:
                render_poly(f, i, sites[i].obj);
                break;

        }
    }

    double mid[2];
    mid[0] = (min[0] + max[0])*0.5;
    mid[1] = (min[1] + max[1])*0.5;

    double radius = std::max((max[0] - min[0])*0.5, (max[1] - min[1])*0.5); 

    VoronoiQuadtree<Site> *qt = new VoronoiQuadtree<Site>(2, mid, radius,
        sites, entcount, max_depth, metric);
    render_voronoi_quadtree(f, qt->root);
    delete qt;

    delete[] sites;
    SHPClose(shp);

    fclose(f);
 
    return 0;
}
