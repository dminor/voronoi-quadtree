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

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <sys/stat.h>
#include <unistd.h>

#include "voronoi_quadtree.h"

double pt_euclidean_metric(double x, double y, struct Site *site)
{ 
    return (x-site->xs[0])*(x-site->xs[0]) + (y-site->ys[0])*(y-site->ys[0]);
}

static double pt_distance_to_line(double x, double y, double a_x, double a_y, double b_x, double b_y)
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

double line_euclidean_metric(double x, double y, struct Site *site)
{ 
    double distance = DBL_MAX;
    for (int i = 0; i < site->n; ++i) {
        double d = pt_distance_to_line(x, y, site->xs[i], site->ys[i], site->xs[i + 1], site->ys[i + 1]); 

        if (d < distance) {
            distance = d; 
        }
    }

    return distance;
}

double poly_euclidean_metric(double x, double y, struct Site *site)
{ 
    double distance = DBL_MAX;
    for (int i = 0; i < site->n; ++i) {

        double d; 
        if (i + 1 == site->n) {
            d = pt_distance_to_line(x, y, site->xs[i], site->ys[i], site->xs[0], site->ys[0]); 
        } else { 
            d = pt_distance_to_line(x, y, site->xs[i], site->ys[i], site->xs[i + 1], site->ys[i + 1]); 
        }

        if (d < distance) { 
            distance = d; 
        }
    }

    return distance;
}

static int closest_site(double x, double y, int n, struct Site *sites)
{
    int site = -1;

    double distance = DBL_MAX;
    for (int i = 0; i < n; ++i) {
        double d = sites[i].metric(x, y, &sites[i]); 
        if (d < distance) {
            site = i;
            distance = d; 
        }
    }

    return site;
}

static struct VoronoiQuadtreeNode *voronoi_quadtree_worker(double x1, double x2, double y1, double y2, int n, struct Site *sites, int depth, int max_depth)
{
    struct VoronoiQuadtreeNode *node = malloc(sizeof(struct VoronoiQuadtreeNode)); 

    int cp1 = closest_site(x1, y1, n, sites);
    int cp2 = closest_site(x1, y2, n, sites);
    int cp3 = closest_site(x2, y1, n, sites);
    int cp4 = closest_site(x2, y2, n, sites);

    if (cp1 == cp2 && cp2 == cp3 && cp3 == cp4 || depth == max_depth) {

        node->ne = node->nw = node->sw = node->se = 0;
        node->x1 = x1;
        node->x2 = x2;
        node->y1 = y1;
        node->y2 = y2;

        node->site_id = cp1;

    } else {

        double xmid = x1 + 0.5 * (x2 - x1);
        double ymid = y1 + 0.5 * (y2 - y1);

        node->ne = voronoi_quadtree_worker(xmid, x2, y1, ymid, n, sites, depth+1, max_depth);
        node->nw = voronoi_quadtree_worker(x1, xmid, y1, ymid, n, sites, depth+1, max_depth);
        node->sw = voronoi_quadtree_worker(x1, xmid, ymid, y2, n, sites, depth+1, max_depth);
        node->se = voronoi_quadtree_worker(xmid, x2, ymid, y2, n, sites, depth+1, max_depth);
    }

    return node;

}

struct VoronoiQuadtreeNode *build_voronoi_quadtree(int n, struct Site *sites, int max_depth)
{
    if (n == 0) return 0;

    //determine bounds
    double x1 = sites[0].xs[0], x2 = sites[0].xs[0], y1 = sites[0].ys[0], y2 = sites[0].ys[0];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < sites[i].n; ++j) {
            if (sites[i].xs[j] < x1) x1 = sites[i].xs[j];
            if (sites[i].xs[j] > x2) x2 = sites[i].xs[j];
            if (sites[i].ys[j] < y1) y1 = sites[i].ys[j];
            if (sites[i].ys[j] > y2) y2 = sites[i].ys[j];
        }
    }

    return voronoi_quadtree_worker(x1, x2, y1, y2, n, sites, 0, max_depth);

}

