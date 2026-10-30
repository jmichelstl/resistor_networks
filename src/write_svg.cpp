/*
This program reads data describing a planar line graph by specifying pairs
of end points for edges, and creates an svg file intended to describe a tensile
test specimen. The outer boundary will be detected and moified to provide tabs
that can be grasped by penumatic grips. Discreet polygons will be detected and
grouped to produce efficient directives for a two-axis cutting device.
*/

#include <iostream>
#include <unordered_map>
#include "network_utils.hpp"
#include <vector>
#include <cmath>
#include <list>
#include <sstream>
#include <algorithm>

using namespace std;

#define PAD .5
#define CCW 0
#define CW 1

inline bool in_range(int idx, int max){
    return idx >= 0 && idx < max;
}

bool plist_comp(vector<Point> list1, vector<Point> list2){
    double meanx1 = 0, meanx2 = 0, meany1 = 0, meany2 = 0;

    for(Point p : list1){
        meanx1 += p.x;
        meany1 += p.y;
    }
    meanx1 /= list1.size();
    meany1 /= list1.size();

    for(Point p : list2){
        meanx2 += p.x;
        meany2 += p.y;
    }
    meanx2 /= list2.size();
    meany2 /= list2.size();

    if(meany1 < meany2) return true;
    if(meany2 < meany1) return false;
    if(meanx1 < meanx2) return true;
    return false;
}

/*
*Produce a more compact list of points for the planar line graph that eliminates
groups of three or more colinear points, and generate the corresponding edges.
*/
void get_svg_file_info(NetworkComplex nc, vector<vector<int>> cycle_list, vector<vector<Point>> &polygons, vector<Point> &border){
    int start_index, curr_index, next_index, cycle_index = 0;
    int min_cycle = 0;
    Point prev_point, curr_point, next_point;
    double slope1, slope2, dx1, dy1, dx2, dy2, angle1, angle2;
    double min = FLT_MAX, cycle_min;
    size_t length;

    for(vector<int> next_cycle : cycle_list){
        //List of points describing the next polygon in a set
        vector<Point> next_poly;
        length = next_cycle.size();

        //Locate a starting index at a corner of a polygon, and not along
        //an edge
        start_index = pick_start(next_cycle, nc);

        //Iterate over indices in the current cycle, beginning with the
        //strategically chosen starting point, and coallesce edges, then add
        //end points to the minimal map and update edge list
        curr_index = start_index;
        prev_point = nc.points[next_cycle[(start_index + length - 1) % length]];
        cycle_min = prev_point.y;

        do{
            next_index = (curr_index + 1) % length;
            curr_point = nc.points[next_cycle[curr_index]];
            next_point = nc.points[next_cycle[next_index]];
            if(curr_point.y < cycle_min) cycle_min = curr_point.y;

            slope1 = (curr_point.y-prev_point.y)/(curr_point.x-prev_point.x);
            slope2 = (next_point.y-curr_point.y)/(next_point.x-curr_point.x);
            if(abs(slope2 - slope1) > FLOAT_TOL){
                next_poly.push_back(curr_point);
            }

            prev_point = curr_point;
            curr_index = next_index;
        }while(curr_index != start_index);

        //Determine whether this cycle has the lowest y coordinate of any cycle
        //so far
        if(cycle_min < min){
            min = cycle_min;
            min_cycle = cycle_index;
        }

        polygons.push_back(next_poly);
        cycle_index ++;
    }

    //The border will contain the lowest point in the planar line graph. This
    //polygon is set aside so that grips may be added and so that it can be
    //color coded
    border = polygons[min_cycle];
    polygons.erase(polygons.begin() + min_cycle, polygons.begin() + min_cycle + 1);

    sort(polygons.begin(), polygons.end(), plist_comp);
}

string poly_path(vector<Point> poly, double minx, double maxy){
    double xpos, ypos;
    ostringstream pstream;
    int index, len;

    len = poly.size();

    pstream << "<path d=\"";

    xpos = poly[0].x - minx;
    ypos = maxy - poly[0].y;
    pstream << "M " << xpos << " " << ypos;

    for(index = 1; index < len; index ++){
        xpos = poly[index].x - minx;
        ypos = maxy - poly[index].y;
        pstream << " L " << xpos << " " << ypos;
    }

    pstream << " z\" fill=\"none\" stroke=\"red\" stroke-width=\".005\"/>";

    return pstream.str();
}

void get_min_max_xy(vector<Point> border, double &minx, double &miny, double &maxx, double &maxy){
    minx = miny = FLT_MAX;
    maxx = maxy = FLT_MIN;

    for(Point next_point : border){
        if(next_point.x < minx) minx = next_point.x;
        if(next_point.y < miny) miny = next_point.y;
        if(next_point.x > maxx) maxx = next_point.x;
        if(next_point.y > maxy) maxy = next_point.y;
    }
}

string arc_string(double xr, double yr, int sense, double dx, double dy){
    ostringstream astream;
    astream << " a" << xr << "," << yr << ",0,0," << sense << "," << dx << "," << dy;
    return astream.str();
}

string get_border_path(vector<Point> bpoints, double minx, double miny, double maxx, double maxy, double crad, double tap_w, double tap_h, double tab_h){

    //Indices of border corners
    int bottom_l, bottom_r, top_l, top_r, index, incr, len, diff;

    //Information about corner positions
    double bottom_minx = FLT_MAX, bottom_maxx = FLT_MIN;
    double top_minx = FLT_MAX, top_maxx = FLT_MIN;

    //Current point along the border path
    double xpos, ypos, bwidth, twidth;
    ostringstream bstream;

    bstream << "<path d=\"";

    len = bpoints.size();

    //Find positions of corners in border
    index = 0;
    for(Point bp : bpoints){
        if(abs(bp.y - miny) < FLOAT_TOL){
            if(bp.x < bottom_minx){
                bottom_minx = bp.x;
                bottom_l = index;
            }
            if(bp.x > bottom_maxx){
                bottom_maxx = bp.x;
                bottom_r = index;
            }
        }

        if(abs(bp.y - maxy) < FLOAT_TOL){
            if(bp.x < top_minx){
                top_minx = bp.x;
                top_l = index;
            }
            if(bp.x > top_maxx){
                top_maxx = bp.x;
                top_r = index;
            }
        }
        index ++;
    }

    bwidth = bottom_maxx - bottom_minx + 2*(crad + tap_w);
    twidth = top_maxx - top_minx + 2*(crad + tap_w);
    minx = minx - PAD - crad - tap_w;
    maxy = maxy + crad + tap_h + tab_h + PAD;

    //Traverse the border in a clockwise manner, adding the left edge, the
    //top tab, the right edge and finally the bottom tab.
    
    //First determine wheter corners are in clockwise or counter-clockwise
    //order
    diff = (bottom_r + len - bottom_l)%len - (top_l + len - bottom_l)%len;
    incr = diff > 0 ? 1 : -1; 

    //Start at the bottom left corner
    xpos = bpoints[bottom_l].x - minx;
    ypos = maxy - bpoints[bottom_l].y;
    bstream << "M " << xpos << " " << ypos;

    //Add left side
    index = bottom_l;
    do{
        index = (index + incr + len) % len;
        xpos = bpoints[index].x - minx;
        ypos = maxy - bpoints[index].y;
        bstream << " L " << xpos << " " << ypos;
    }while(index != top_l);

    //Add the top tab
    if(tap_w > 0){
        bstream << arc_string(tap_w, tap_h, CCW, -tap_w, -tap_h);
    }
    if(crad > 0){
        bstream << arc_string(crad, crad, CW, -crad, -crad);
    }
    bstream << " l 0 " << -tab_h;
    bstream << " l " << twidth << " 0"; 
    bstream << " l 0 " << tab_h;
    if(crad > 0){
        bstream << arc_string(crad, crad, CW, -crad, crad);
    }
    if(tap_w > 0){
        bstream << arc_string(tap_w, tap_h, CCW, -tap_w, tap_h);
    }

    //Add the right side
    index = top_r;
    do{
        index = (index + incr + len) % len;
        xpos = bpoints[index].x - minx;
        ypos = maxy - bpoints[index].y;
        bstream << " L " << xpos << " " << ypos;
    }while(index != bottom_r);

    //Add the bottom tab
    if(tap_w > 0){
        bstream << arc_string(tap_w, tap_h, CCW, tap_w, tap_h);
    }
    if(crad > 0){
        bstream << arc_string(crad, crad, CW, crad, crad);
    }
    bstream << " l 0 " << tab_h;
    bstream << " l -" << bwidth << " 0"; 
    bstream << " l 0 " << -tab_h;
    if(crad > 0){
        bstream << arc_string(crad, crad, CW, crad, -crad);
    }
    if(tap_w > 0){
        bstream << arc_string(tap_w, tap_h, CCW, tap_w, -tap_h);
    }

    bstream << "\" fill=\"none\" stroke=\"blue\" stroke-width=\".005\"/>";

    return bstream.str();
}

void write_svg(vector<Point> border, vector<vector<Point>> polygons, double tap_w, double tap_h, double crad, double tab_h, string file_name){

    double width, height, minx, miny, maxx, maxy;
    string svg_open, path;
    ofstream svg_file;

    svg_file.open(file_name, ofstream::out);
    if(! svg_file.is_open()) return;

    //Determine minimum and maximum x and y values for the document
    get_min_max_xy(border, minx, miny, maxx, maxy);

    width = maxx - minx + 2 * (crad + tap_w + PAD);
    height = maxy - miny + 2 * (tap_h + crad + tab_h + PAD);

    //Write a generic svg file header
    svg_file << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
    svg_file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n"; 
    svg_file << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";

    //Open top-level svg tag and declare document properties
    svg_file << "<svg width=\"" << width << "cm\" height=\""<< height <<"cm\"";
    svg_file << " viewBox=\"0 0 " << width << " " << height << "\"";
    svg_file << " version = \"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";

    //Open group tag to group all paths in the design together
    svg_file << "<g>\n";

    //Create the path for the border and add it as the first path
    svg_file << "\t" << get_border_path(border, minx, miny, maxx, maxy, crad, tap_w, tap_h, tab_h) << "\n";

    //Now add a path for each interior polygon
    minx = minx - PAD - crad - tap_w;
    maxy = maxy + crad + tap_h + tab_h + PAD;

    for(vector<Point> nextPoly : polygons){
        svg_file << "\t" << poly_path(nextPoly, minx, maxy) << "\n";
    }

    //Close the group
    svg_file << "</g>\n";

    //Close the SVG document
    svg_file << "</svg>\n";

    svg_file.close();
}

/*
*Prompt the user for the name of a data file and an output SVG file. Create
an open input stream, obtain a map of points to indices, and a map of indices
to lists of other indices. Use the two maps to build a table of cycles, then
convert the table of cycles to polygon data.
*/
bool export_svg(NetworkComplex nc, double tap_w, double tap_h, double crad, double tab_h, string file_name){
    ofstream svg_file;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    vector<vector<Point>> polygons;
    vector<Point> border;

    //Obtain a map from unique vertices to neighboring vertices
    build_plc_map(nc, emap);

    //Use the map to build cycles described by integer indices in cyclical order
    cycle_list = build_cycles(emap);

    //Get final data to construct an SVG file
    get_svg_file_info(nc, cycle_list, polygons, border); 

    write_svg(border, polygons, tap_w, tap_h, crad, tab_h, file_name);

    return 0;
}
