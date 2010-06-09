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

#ifndef VORONOI_QUADTREE_H_ 
#define VORONOI_QUADTREE_H_ 

struct VoronoiQuadtreeNode {
    struct VoronoiQuadtreeNode *ne, *nw, *sw, *se;
    double x1, x2, y1, y2;
    int site_id;
};

struct VoronoiQuadtreeNode *build_voronoi_quadtree_pts(int n, double *xs, double *ys, int max_depth);
struct VoronoiQuadtreeNode *build_voronoi_quadtree_lines(int n, double *xs, double *ys, int *sites, int max_depth);
struct VoronoiQuadtreeNode *build_voronoi_quadtree_polys(int n, double *xs, double *ys, int *sites, int max_depth);

#endif
