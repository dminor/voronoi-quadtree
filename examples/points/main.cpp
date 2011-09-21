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

#include "voronoi_quadtree.h"

struct Point {

    double v[2]; 
    int id;

    double &operator[](const int &index) 
    { 
        return v[index]; 
    } 

    const double &operator[](const int &index) const 
    { 
        return v[index]; 
    } 
};

double metric(Point *pt1, double *pt2)
{ 
    return (pt2[0]-(*pt1)[0])*(pt2[0]-(*pt1)[0]) +
        (pt2[1]-(*pt1)[1])*(pt2[1]-(*pt1)[1]);
}

void render_point(FILE *f, int site, double x, double y)
{
    fprintf(f, "colour-site-%d\n", site); 
    fprintf(f, "%.0f %.0f 2 0 360 arc\n", x, y);
    fprintf(f, "fill\n"); 
}

void render_voronoi_quadtree(FILE *f, VoronoiQuadtree<Point>::Node *node)
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

Point *read_points(FILE *f, int *pt_count)
{
    int err;
    char buf[80];
    err = fscanf(f, "%d", pt_count);
    fgets(buf, 80, f);

    if (err != 1 || pt_count < 0) {
        fprintf(stderr, "error: invalid point count %d\n", *pt_count);
        return 0;
    }

    Point *pts = new Point[*pt_count];

    double x, y;
    for (int i = 0; i < *pt_count; ++i) {
        err = fscanf(f, "%lf, %lf", &x, &y);
        if (err != 2) {
            fprintf(stderr, "error: invalid points in file\n"); 
        }
        pts[i][0] = x;
        pts[i][1] = y; 
        pts[i].id = i;
    }

    return pts; 
}


int main(int argc, char **argv)
{ 
    //handle command line
    if (argc < 2) {
        fprintf(stderr, "usage: render-quadtree <filename> [maximum depth]\n");
        exit(1);
    }

    int max_depth = 8;
    if (argc == 3) max_depth = atoi(argv[2]);

    //read data
    int ptcount;
    FILE *f = fopen(argv[1], "r");
    if (!f) {
        fprintf(stderr, "error: could not open points file: %s\n", argv[1]);
        exit(1); 
    }

    Point *pts = read_points(f, &ptcount); 
    fclose(f);

    if (!pts) {
        exit(0);
    }

    //setup output
    fprintf(stdout, "%%\n");

    //define box function for later
    fprintf(stdout, "/box {\n");
    fprintf(stdout, "    /y2 exch def\n");
    fprintf(stdout, "    /y1 exch def\n");
    fprintf(stdout, "    /x2 exch def\n");
    fprintf(stdout, "    /x1 exch def\n");
    fprintf(stdout, "    gsave\n");
    fprintf(stdout, "    newpath\n");
    fprintf(stdout, "    x1 y1 moveto\n");
    fprintf(stdout, "    x1 y2 lineto\n");
    fprintf(stdout, "    x2 y2 lineto\n");
    fprintf(stdout, "    x2 y1 lineto\n");
    fprintf(stdout, "    closepath\n");
    fprintf(stdout, "    stroke\n");
    fprintf(stdout, "    grestore\n");
    fprintf(stdout, "} def\n");

    //setup colours
    for (int i = 0; i < ptcount; ++i) {
        fprintf(stdout, "/colour-site-%d {%.1f %.1f %.1f setrgbcolor } def\n", i,
            (double)rand()/(double)RAND_MAX,
            (double)rand()/(double)RAND_MAX,
            (double)rand()/(double)RAND_MAX);
    }

    //get bounds and render points
    double min[2], max[2];
    min[0] = pts[0][0]; min[1] = pts[0][1];
    max[0] = pts[0][0]; max[1] = pts[0][1];

    for (int i = 1; i < ptcount; ++i) {
        render_point(stdout, i, pts[i][0], pts[i][1]);
        if (pts[i][0] < min[0]) min[0] = pts[i][0];
        if (pts[i][1] < min[1]) min[1] = pts[i][1];
        if (pts[i][0] > max[0]) max[0] = pts[i][0];
        if (pts[i][1] > max[1]) max[1] = pts[i][1];
    }

    double mid[2];
    mid[0] = (min[0] + max[0])*0.5;
    mid[1] = (min[1] + max[1])*0.5;

    double radius = std::max((max[0] - min[0])*0.5, (max[1] - min[1])*0.5); 

    VoronoiQuadtree<Point> *qt = new VoronoiQuadtree<Point>(2, mid, radius,
        pts, ptcount, max_depth, metric);

    render_voronoi_quadtree(stdout, qt->root);
    delete qt;

    delete[] pts;
 
    return 0;
}
