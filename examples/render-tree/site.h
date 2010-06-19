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

#ifndef SITE_H_
#define SITE_H_

#include <limits>

struct Site {

    double *xs;
    double *ys;
    size_t n;

    int id;

    Site() : n(0) {};

    virtual ~Site()
    {
        if (n) {
            delete[] xs;
            delete[] ys;
        }
    }; 

    virtual double distance_to(double x, double y) = 0; 

};

struct PointSite : public Site {

    double distance_to(double x, double y)
    { 
        return (x-this->xs[0])*(x-this->xs[0]) + (y-this->ys[0])*(y-this->ys[0]);
    }

};

struct LineSite : public Site {

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

    double distance_to(double x, double y)
    { 
        double distance = std::numeric_limits<double>::max(); 
        for (int i = 0; i < n; ++i) {
            double d = pt_distance_to_line(x, y, xs[i], ys[i], xs[i + 1], ys[i + 1]); 

            if (d < distance) {
                distance = d; 
            }
        }

        return distance;
    }

};

struct PolygonSite : public LineSite {

    double distance_to(double x, double y)
    {
        double distance = std::numeric_limits<double>::max(); 
        for (int i = 0; i < n; ++i) {

            double d; 
            if (i + 1 == n) {
                d = pt_distance_to_line(x, y, xs[i], ys[i], xs[0], ys[0]); 
            } else { 
                d = pt_distance_to_line(x, y, xs[i], ys[i], xs[i + 1], ys[i + 1]); 
            }

            if (d < distance) { 
                distance = d; 
            }
        }

        return distance;
    }
}; 

#endif
