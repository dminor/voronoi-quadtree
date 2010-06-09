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

typedef int (*closest_pt_fn)(unsigned int n, double x, double y, double *xs, double *ys, int *sites);

static int closest_pt(unsigned int n, double x, double y, double *xs, double *ys, int *sites)
{
    int site = -1;

    double distance = DBL_MAX;
    for (int i = 0; i < n; ++i) {
        double d = (x-xs[i])*(x-xs[i]) + (y-ys[i])*(y-ys[i]);
        if (d < distance) {
            site = i;
            distance = d; 
        }
    }

    return site;
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

static int closest_pt_lines(unsigned int n, double x, double y, double *xs, double *ys, int *sites)
{
    int site = -1;

    int current_site = sites[0];

    double distance = DBL_MAX;
    for (int i = 0; i < n; ++i) {
        //if end of object, skip to next
        if (i + 1 != n && sites[i + 1] != current_site) {
            current_site = sites[i + 1];
        } else { 
            double d = pt_distance_to_line(x, y, xs[i], ys[i], xs[i + 1], ys[i + 1]); 

            if (d < distance) {
                site = sites[i];
                distance = d; 
            }
        }
    }

    return site;
}

static int closest_pt_polys(unsigned int n, double x, double y, double *xs, double *ys, int *sites)
{
    int site = -1;

    int start_index = 0;
    int current_site = sites[0];

    double distance = DBL_MAX;
    for (int i = 0; i < n; ++i) {

        double d;

        //if next site is different, we need to test line that closes polygon
        if (i + 1 == n || sites[i + 1] != current_site) {
            d = pt_distance_to_line(x, y, xs[i], ys[i], xs[start_index], ys[start_index]); 
            if (i + 1 != n) {
                current_site = sites[i + 1];
                start_index = i + 1;
            }
        } else { 
            d = pt_distance_to_line(x, y, xs[i], ys[i], xs[i + 1], ys[i + 1]); 
        }

        if (d < distance) {
            site = sites[i];
            distance = d; 
        }
    }

    return site;
}

static struct VoronoiQuadtreeNode *voronoi_quadtree_worker(closest_pt_fn fn, double x1, double x2, double y1, double y2, int n, int depth, int max_depth, double *xs, double *ys, int *sites)
{
    struct VoronoiQuadtreeNode *node = malloc(sizeof(struct VoronoiQuadtreeNode)); 

    int cp1 = fn(n, x1, y1, xs, ys, sites);
    int cp2 = fn(n, x1, y2, xs, ys, sites);
    int cp3 = fn(n, x2, y1, xs, ys, sites);
    int cp4 = fn(n, x2, y2, xs, ys, sites);

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

        node->ne = voronoi_quadtree_worker(fn, xmid, x2, y1, ymid, n, depth+1, max_depth, xs, ys, sites); 
        node->nw = voronoi_quadtree_worker(fn, x1, xmid, y1, ymid, n, depth+1, max_depth, xs, ys, sites); 
        node->sw = voronoi_quadtree_worker(fn, x1, xmid, ymid, y2, n, depth+1, max_depth, xs, ys, sites); 
        node->se = voronoi_quadtree_worker(fn, xmid, x2, ymid, y2, n, depth+1, max_depth, xs, ys, sites); 
    }

    return node;

}

struct VoronoiQuadtreeNode *build_voronoi_quadtree_pts(int n, double *xs, double *ys, int max_depth)
{ 
    //determine bounds
    double x1 = xs[0], x2 = xs[0], y1 = ys[0], y2 = ys[0];
    for (int i = 0; i < n; ++i) {
        if (xs[i] < x1) x1 = xs[i];
        if (xs[i] > x2) x2 = xs[i];
        if (ys[i] < y1) y1 = ys[i];
        if (ys[i] > y2) y2 = ys[i];
    }

    return voronoi_quadtree_worker(closest_pt, x1, x2, y1, y2, n, 0, max_depth, xs, ys, 0);
}

struct VoronoiQuadtreeNode *build_voronoi_quadtree_lines(int n, double *xs, double *ys, int *sites, int max_depth)
{
    //determine bounds
    double x1 = xs[0], x2 = xs[0], y1 = ys[0], y2 = ys[0];
    for (int i = 0; i < n; ++i) {
        if (xs[i] < x1) x1 = xs[i];
        if (xs[i] > x2) x2 = xs[i];
        if (ys[i] < y1) y1 = ys[i];
        if (ys[i] > y2) y2 = ys[i];
    }

    return voronoi_quadtree_worker(closest_pt_lines, x1, x2, y1, y2, n, 0, max_depth, xs, ys, sites);

}

struct VoronoiQuadtreeNode *build_voronoi_quadtree_polys(int n, double *xs, double *ys, int *sites, int max_depth)
{
    //determine bounds
    double x1 = xs[0], x2 = xs[0], y1 = ys[0], y2 = ys[0];
    for (int i = 0; i < n; ++i) {
        if (xs[i] < x1) x1 = xs[i];
        if (xs[i] > x2) x2 = xs[i];
        if (ys[i] < y1) y1 = ys[i];
        if (ys[i] > y2) y2 = ys[i];
    }

    return voronoi_quadtree_worker(closest_pt_polys, x1, x2, y1, y2, n, 0, max_depth, xs, ys, sites);

}
