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

#include <limits>

template<class Site> class VoronoiQuadtree {

public:

    struct Node {
        Node *ne, *nw, *sw, *se;
        double x1, x2, y1, y2;
        Site *site;
    };

    VoronoiQuadtree(double x1, double x2, double y1, double y2, Site *sites, size_t n, size_t max_depth, double (*metric)(Site *site, double x, double y)) : sites(sites), n(n), max_depth(max_depth), metric(metric)
    {        
        root = worker(x1, x2, y1, y2, 0); 
    }

    virtual ~VoronoiQuadtree()
    { 
    }

    Node *root;
    double (*metric)(Site *site, double x, double y);

private:

    Site *closest_site(double x, double y)
    {
        Site *site = 0;

        double distance = std::numeric_limits<double>::max();
        for (int i = 0; i < n; ++i) {
            double d = metric(&sites[i], x, y);
            if (d < distance) {
                site = &sites[i];
                distance = d; 
            }
        }

        return site; 
    }

    Node *worker(double x1, double x2, double y1, double y2, int depth)
    {
        Node *node = new Node; 

        Site *site1 = closest_site(x1, y1);
        Site *site2 = closest_site(x1, y2);
        Site *site3 = closest_site(x2, y1);
        Site *site4 = closest_site(x2, y2);

        if (site1 == site2 && site2 == site3 && site3 == site4 || depth == max_depth) {

            node->ne = node->nw = node->sw = node->se = 0;
            node->x1 = x1;
            node->x2 = x2;
            node->y1 = y1;
            node->y2 = y2;

            node->site = site1;

        } else {

            double xmid = x1 + 0.5 * (x2 - x1);
            double ymid = y1 + 0.5 * (y2 - y1);

            node->ne = worker(xmid, x2, y1, ymid, depth+1);
            node->nw = worker(x1, xmid, y1, ymid, depth+1);
            node->sw = worker(x1, xmid, ymid, y2, depth+1);
            node->se = worker(xmid, x2, ymid, y2, depth+1);
        }

        return node; 

    }

    Site *sites;
    size_t n;
    size_t max_depth;
};

#endif
