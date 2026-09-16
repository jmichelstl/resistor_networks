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
If a point is not currently in an index map, add it and map it to the integer
index one greater than the current size of the map.
*/
void update_maps(Point p, map<Point, int> &pmap, map<int, list<int>> &emap){
    int index;

    if(pmap.find(p) == pmap.end()){
        index = pmap.size() + 1;
        pmap.insert(make_pair(p, index));
        emap.insert(make_pair(index, list<int>()));
    }
}

inline bool in_range(int idx, int max){
    return idx >= 0 && idx < max;
}

/*
This function processes a file describing a planar line graph by specifying
pairs of endpoints of line segments.

Arguments:
datfile - An open stream to a file describing the planar line graph
imap - A map from vertices in the planar line graph to integer indices
emap - A map in which each pair is a mapping from a point's integer index to
the indices of all other points to which that point is connected
*/
bool read_edges(ifstream &datfile, map<int, Point> &imap, map<int, list<int>> &emap){

    string nextline;
    Point p1, p2;
    int num_pts, num_edges, edge_count = 0, idx1, idx2;
    double x, y;

    //Attempt to read a header
    if(! datfile.eof()){
        getline(datfile, nextline);
        if(sscanf(nextline.c_str(), "%d %d", &num_pts, &num_edges) < 2){
            cerr << "No header could be read.\n";
	        datfile.close();
	        return false;
        }

        if(num_pts <= 0 || num_edges <= 0){
            cerr << "An invalid header was found.\n";
            datfile.close();
	        return false;
        }
    }

    //Attempt to read the points in the planar line graph
    while(! datfile.eof() && imap.size() < num_pts){
        getline(datfile, nextline);

        if(sscanf(nextline.c_str(), "%lf %lf", &x, &y) < 2){
            cout << "A valid point could not be read.\n";
            datfile.close();
            return false;
        }

        imap.insert(make_pair(imap.size(), Point(x, y)));
    }

    if(imap.size() < num_pts){
        cerr << "Too few points were read.\n";
        datfile.close();
        return false;
    }

    for(int idx = 0; idx < num_pts; idx++){
        emap.insert(make_pair(idx, list<int>()));
    }

    while(! datfile.eof() && edge_count < num_edges){
        getline(datfile, nextline);

        if(sscanf(nextline.c_str(), "%d %d", &idx1, &idx2) < 2){
            cerr << "A valid edge could not be read.\n";
            datfile.close();
            return false;
        }

        if(! in_range(idx1, num_pts) || ! in_range(idx2, num_pts)){
            cerr << "An invalid index was encountered.\n";
            datfile.close();
	        return false;
        }

        emap[idx1].push_back(idx2);
        emap[idx2].push_back(idx1);

        edge_count ++;
    }

    datfile.close();

    if(edge_count < num_edges){
        cerr << "Too few edges were read.\n";
        return false;
    }

    return true;
}

vector<vector<int>> build_cycles(map<int, list<int>> emap){
    int curr_index, next_index;
    map<int, list<int>>::iterator emap_pos;
    vector<int> curr_cycle;
    vector<vector<int>> cycle_list;
    bool cycle_complete;

    //Build cycles by iterating over entries in the map from indices to
    //lists of neighbor indices
    while(! emap.empty()){
        //Start a new cycle by finding the first element of the map and
        //making the key the first element in the cycle
        vector<int> curr_cycle;
        curr_index = emap.begin()->first;
        cycle_complete = false;

        //Add successive elements to the current cycle by finding the first
        //neighbor listed for the current vertex in the cycle
        do{
            next_index = emap[curr_index].front();
            curr_cycle.push_back(next_index);

            //Remove the current index and the next index from each other's
            //lists of neighbors
            emap[curr_index].remove(next_index);
            if(emap[curr_index].empty()){
                emap_pos = emap.find(curr_index);
                emap.erase(emap_pos);
            }

            emap[next_index].remove(curr_index);
            if(emap[next_index].empty()){
                emap_pos = emap.find(next_index);
                emap.erase(emap_pos);
                //When the list to which next_index is mapped has no elements,
                //the current cycle is complete
                cycle_complete = true;
            }
            curr_index = next_index;
        }while(! cycle_complete);
        cycle_list.push_back(curr_cycle);
    }

    return cycle_list;
}

/*
*Locate a starting point in a cycle that is not part of a colinear group of
*three or more points
*/
int pick_start(vector<int> cycle, map<int, Point> imap){
    int start_index, curr_index, next_index_1, next_index_2;
    Point curr_p, next_p_1, next_p_2;
    double slope1, slope2;
    size_t length = cycle.size();

    for(curr_index = 0; curr_index < length; curr_index ++){
        next_index_1 = (curr_index + 1) % length;
        next_index_2 = (curr_index + 2) % length;
        curr_p = imap[cycle[curr_index]];
        next_p_1 = imap[cycle[next_index_1]];
        next_p_2 = imap[cycle[next_index_2]];

        //If the slopes of the next two edges are different, choose the
        //next index as the starting index
        slope1 = (next_p_1.y - curr_p.y) / (next_p_1.x - curr_p.x);
        slope2 = (next_p_2.y - next_p_1.y) / (next_p_2.x - next_p_1.x);

        if(abs(slope2 - slope1) > FLOAT_TOL){
            start_index = next_index_1;
            break;
        }
    }

    return start_index;
}

/*
*Produce a more compact list of points for the planar line graph that eliminates
groups of three or more colinear points, and generate the corresponding edges.
*/
void get_svg_file_info(map<int, Point> imap, vector<vector<int>> cycle_list, vector<vector<Point>> &polygons, vector<Point> &border){
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
        start_index = pick_start(next_cycle, imap);

        //Iterate over indices in the current cycle, beginning with the
        //strategically chosen starting point, and coallesce edges, then add
        //end points to the minimal map and update edge list
        curr_index = start_index;
        prev_point = imap[next_cycle[(start_index + length - 1) % length]];
        cycle_min = prev_point.y;

        do{
            next_index = (curr_index + 1) % length;
            curr_point = imap[next_cycle[curr_index]];
            next_point = imap[next_cycle[next_index]];
            if(curr_point.y < cycle_min) cycle_min = curr_point.y;

            slope1 = (curr_point.y-prev_point.y)/(curr_point.x-prev_point.x);
            slope2 = (next_point.y-curr_point.y)/(next_point.x-curr_point.x);
            if(abs(slope2 - slope1) > FLOAT_TOL){
                next_poly.push_back(curr_point);
            }

            prev_point = curr_point;
            curr_index = next_index;
        }while(curr_index != start_index);

        //cout << "Survived.\n";

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
            //cout << "bp.x: " << bp.x << "\n";
            if(bp.x < bottom_minx){
                bottom_minx = bp.x;
                bottom_l = index;
            }
            if(bp.x > bottom_maxx){
                bottom_maxx = bp.x;
                //cout << "Now it's " << bottom_maxx << "\n";
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

    //cout << "Maxx: " << bottom_maxx << " Minx: " << bottom_minx << "\n";
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

void write_svg(vector<Point> border, vector<vector<Point>> polygons){

    double width, height, minx, miny, maxx, maxy;
    double tap_w, tap_h, crad, tab_h;
    vector<double> response;
    string svg_open, path;
    bool valid = false;
    ofstream svg_file;

    open_output_file("Enter the name for the SVG file: ", svg_file);
    if(! svg_file.is_open()) return;

    //Obtain information about the grip portions
    do{
        response = getdoubles("Enter taper width and height, corner radius and tab height: ");
        if(response.size() >= 4 && response[0] >= 0 && response[1] >= 0 && response[2] >= 0 && response[3] >= 0) valid = true;
        else cerr << "Invalid response. Enter four non-negative numbers.\n";
    }while(! valid);

    tap_w = response[0];
    tap_h = response[1];
    crad = response[2];
    tab_h = response[3];

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
int main(int argc, char **argv){
    ifstream datfile;
    ofstream svg_file;
    map<int, Point> imap;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    vector<vector<Point>> polygons;
    vector<Point> border;

    //Obtain network data and obtain maps from points to integers and from
    //integers to lists of neighbor indices
    open_dat_file("Enter the data file describing the network: ", datfile);
    if(! datfile.is_open()){
        cerr << "The specified file could not be opened for reading.\n";
        exit(1);
    }

    if(! read_edges(datfile, imap, emap)){
        cerr << "A valid planar line graph could not be read.\n";
        exit(2);
    }

    //Use maps to build cycles described by integer indices in cyclical order
    cycle_list = build_cycles(emap);

    //Get final data to construct an SVG file
    get_svg_file_info(imap, cycle_list, polygons, border); 

    write_svg(border, polygons);

    return 0;
}
