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
        Node **nodes;
        double *mid;
        double radius;
        Site *site;

        Node() : nodes(0), mid(0), site(0)
        {
        }

        virtual ~Node()
        {
            if (mid) delete[] mid;
            if (nodes) delete[] nodes;
        }
    };

    VoronoiQuadtree(int dim, double *mid, double radius,
        Site *sites, int n, int max_depth,
        double (*metric)(Site *site, double *pt)) :
        dim(dim), nnodes(2*dim), sites(sites), n(n),
        max_depth(max_depth), metric(metric)
    {        
        root = worker(mid, radius, 0); 
    }

    virtual ~VoronoiQuadtree()
    { 
        delete_worker(root);
    }

    Node *root;
    double (*metric)(Site *site, double *pt);

private:

    Site *closest_site(double *pt)
    {
        Site *site = 0;

        double distance = std::numeric_limits<double>::max();
        for (int i = 0; i < n; ++i) {
            double d = metric(&sites[i], pt);
            if (d < distance) {
                site = &sites[i];
                distance = d; 
            }
        }

        return site; 
    }

    Node *worker(double *mid, double radius, int depth)
    { 
        Node *node = new Node; 

        Site *closest = 0, *last_closest = 0;

        double *pt = new double[dim];

        bool all_same = true;

        if (depth != max_depth) {
            for (int i = 0; i < nnodes && all_same; ++i) { 
                for (int d = 0; d < dim; ++d) {
                    if (i & (1 << d)) {
                        pt[d] = mid[d] + radius; 
                    } else { 
                        pt[d] = mid[d] - radius; 
                    } 
                }

                closest = closest_site(pt);

                if (last_closest != 0 && closest != last_closest) {
                    all_same = false;
                    break;
                }

                last_closest = closest;
            }
        } else {
            closest = closest_site(mid);
        }

        if (all_same) {

            node->mid = new double[dim];

            for (int d = 0; d < dim; ++d) {
                node->mid[d] = mid[d];
                node->radius = radius;
            }

            node->site = closest;

        } else {

            node->nodes = new Node *[nnodes];

            double new_radius = 0.5*radius;
            for (int i = 0; i < nnodes; ++i) { 
                for (int d = 0; d < dim; ++d) {
                    if (i & (1 << d)) {
                        pt[d] = mid[d] + new_radius; 
                    } else { 
                        pt[d] = mid[d] - new_radius; 
                    } 
                }

                node->nodes[i] = worker(pt, new_radius, depth + 1); 
            } 
        }

        delete[] pt;

        return node; 

    }

    void delete_worker(Node *n)
    { 
        if (n->nodes) {
            for (int i = 0; i < nnodes; ++i) {
                if (n->nodes[i]) delete_worker(n->nodes[i]);
            }
        } 

        delete n;
    }

    Site *sites;
    int dim;
    int nnodes;
    int n;
    int max_depth;
};

#endif
