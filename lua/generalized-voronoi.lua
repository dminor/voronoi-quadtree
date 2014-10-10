--[[
Copyright (c) 2011 Daniel Minor 

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
--]]

gd = require('gd')
shapelib = require('shapelib')

MAX_DEPTH = 6 

function readpoints(filename)
    --read points from a text file

    local pts = {}
    local file = assert(io.open(filename, 'r'))

    local line
    for line in file:lines() do 
        local pt = {} 
        _, _, pt.x, pt.y = line:find('([.%w]+), ([.%w]+)')

        if pt.x and pt.y then 
            pt.x = tonumber(pt.x)
            pt.y = tonumber(pt.y)
            table.insert(pts, pt)
        end
    end

    file:close()

    return pts
end

function manhattan_metric(pt1, pt2)
    -- L1 metric
    return math.abs(pt1.x-pt2.x) + math.abs(pt1.y-pt2.y)
end

function euclidean_metric(pt1, pt2)
    -- L2 metric
    return (pt2.x-pt1.x)*(pt2.x-pt1.x) + (pt2.y-pt1.y)*(pt2.y-pt1.y) 
end

function closest_pt(metric)
    --closure for closest pt function for provided metric

    return function (pt, pts)
        local closest = 1 
        local distance = metric(pt, pts[1])

        for k,v in ipairs(pts) do
            local d = metric(pt, v)
            if d < distance then
                closest = k 
                distance = d
            end 
        end 

        return closest, distance 
    end 
end 

function closest_pt_on_line(pt, a, b)
    -- return closest point on line segment (a,b) to pt

    local vba_x = b.x - a.x
    local vba_y = b.y - a.y

    local vpta_x = pt.x - a.x
    local vpta_y = pt.y - a.y

    local denom = (vba_x*vba_x + vba_y*vba_y)
    if denom == 0.0 then
        return a
    end

    local t = (vba_x*vpta_x + vba_y*vpta_y) / denom
    if t > 1.0 then t = 1.0 end
    if t < 0.0 then t = 0.0 end

    return {x=a.x + t * vba_x,y=a.y + t * vba_y} 

end

function closest_pt_poly(pt, polys) 
    --return closest point on polygon to pt

    local closest = 1 
    local distance = euclidean_metric(pt, polys[1][1])
    local closest_pt = closest_pt(euclidean_metric)

    for k,v in ipairs(polys) do

        pts = {} 
        for i=1,#v-1 do
            table.insert(pts, closest_pt_on_line(pt, v[i], v[i+1])) 
        end
        table.insert(pts, closest_pt_on_line(pt, v[1], v[#v])) 

        local d
        _, d = closest_pt(pt, pts)

        if d < distance then
            closest = k
            distance = d
        end
    end

    return closest, distance 
end

function build_quadtree(x1, x2, y1, y2, closest_fn, objs, depth) 
    --recursively build quadtree approximating Voronoi diagram

    depth = depth or 0

    local qt = {}

    local cp1 = closest_fn({x=x1,y=y1}, objs)
    local cp2 = closest_fn({x=x1,y=y2}, objs)
    local cp3 = closest_fn({x=x2,y=y1}, objs)
    local cp4 = closest_fn({x=x2,y=y2}, objs)

    -- leaf node, store square and value
    if (cp1 == cp2 and cp2 == cp3 and cp3 == cp4) or depth > MAX_DEPTH then
        qt.square = {}
        qt.square.x1 = x1
        qt.square.x2 = x2
        qt.square.y1 = y1
        qt.square.y2 = y2

        qt.closest_site = cp1 
    else

        local xmid = x1 + 0.5 * (x2 - x1)
        local ymid = y1 + 0.5 * (y2 - y1)

        qt.ne = build_quadtree(xmid, x2, y1, ymid, closest_fn, objs, depth + 1)
        qt.nw = build_quadtree(x1, xmid, y1, ymid, closest_fn, objs, depth + 1)
        qt.sw = build_quadtree(x1, xmid, ymid, y2, closest_fn, objs, depth + 1)
        qt.se = build_quadtree(xmid, x2, ymid, y2, closest_fn, objs, depth + 1)

    end

    return qt

end

function render_voronoi_diagram(im, qt, colours, transform)
    --render voronoi diagram to image

    if qt.square then
        if transform then
            local x1; local y1
            x1, y1 = transform(qt.square.x1, qt.square.y1)

            local x2; local y2
            x2, y2 = transform(qt.square.x2, qt.square.y2)

            im:rectangle(x1, y1, x2, y2, colours[qt.closest_site])
        else
            im:rectangle(qt.square.x1, qt.square.y1, qt.square.x2, qt.square.y2, colours[qt.closest_site]) 
        end 
    else
        render_voronoi_diagram(im, qt.ne, colours, transform)
        render_voronoi_diagram(im, qt.nw, colours, transform)
        render_voronoi_diagram(im, qt.sw, colours, transform)
        render_voronoi_diagram(im, qt.se, colours, transform)
    end 
end

function transform_fn(x_offset, y_offset, x_scale, y_scale, height) 
    -- transform point from (x,y) to image (u,v)
    return function(x, y)
        return (x-x_offset)/x_scale, height-(y-y_offset)/y_scale 
    end
end

function voronoi_pts(filename)
    -- generate voronoi diagram based upon points in filename

    local pts = readpoints(filename)

    --find bounds of region
    local x1 = pts[1].x 
    local x2 = pts[1].x 
    local y1 = pts[1].y 
    local y2 = pts[1].y 

    for k,v in ipairs(pts) do
        if v.x < x1 then
            x1 = v.x
        end
        if v.x > x2 then
            x2 = v.x
        end
        if v.y < y1 then
            y1 = v.y
        end
        if v.y > y2 then
            y2 = v.y
        end
    end

    local qt = build_quadtree(x1, x2, y1, y2, closest_pt(euclidean_metric), pts, 0) 

    local im = gd.createTrueColor(x2-x1,y2-y1) 
    local colours = {}

    for k, v in ipairs(pts) do
        colours[k] = im:colorAllocate(math.random(255), math.random(255), math.random(255)) 
    end

    render_voronoi_diagram(im, qt, colours)
    im:png('voronoi.png') 
end

function render_polys(im, polys, colours, transform) 
    --render polygons from shapefiles 

    for k, v in ipairs(polys) do
        local pts = {}
        for _, pt in ipairs(v) do
            local x; local y
            x, y = transform(pt.x, pt.y)
            table.insert(pts, {x,y})
        end
        im:filledPolygon(pts, colours[k]) 
    end 
end

function voronoi_polys(filename)
    -- generate voronoi diagram based upon shapefile areal features

    local shp = shapelib.open(filename, 'rb')

    local entcount; local min_x; local min_y; local max_x; local max_y
    entcount, _, x1, y1, _, _, x2, y2 = shapelib.getinfo(shp)

    local polys = {}

    for i=1,entcount do
        local poly = {}
        local o = shapelib.readobject(shp, i)
        for v=1,o.Vertices do
            local pt = {x=o.X[v], y=o.Y[v]}
            table.insert(poly, pt)
        end

        table.insert(polys, poly)
    end 

    shapelib.close(shp)

    local qt = build_quadtree(x1, x2, y1, y2, closest_pt_poly, polys, 0) 

    local im = gd.createTrueColor(1000, 1000)
    local colours = {}

    for k, v in ipairs(polys) do
        colours[k] = im:colorAllocate(math.random(255), math.random(255), math.random(255)) 
    end

    render_voronoi_diagram(im, qt, colours, transform_fn(x1, y1, 5, 5, 1000))
    render_polys(im, polys, colours, transform_fn(x1, y1, 5, 5, 1000))
    im:png('voronoi.png') 

end

if #arg ~= 2 then
    print('usage: generalized-voronoi.lua (-pts|-poly) filename')
    return
end

if arg[1] == '-pts' then
    voronoi_pts(arg[2])
else
    voronoi_polys(arg[2]) 
end
 
