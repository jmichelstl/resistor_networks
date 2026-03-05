/*Author: Jonathan MichenyThis program prompts a user for boundaries and rules for producing lattice
sites and nearest neighbor connections. From here, a set of line segments
joining nearest neighbor locations is produced. Connections can be removed
from the network in two ways; one method is completely random, while the other
uses randomness but also guarantees that the network remains fully connected.
*/

#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <fstream>
#include <stack>
#include "network_utils.hpp"
#include "hnet_bvh.hpp"
#include <unistd.h>
#include <ctype.h>
#include <float.h>
#include <random>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <tuple>
#include <time.h>
#include <memory>

using namespace std;

//Data type for storing ordered pairs of indices into lists of geometrical
//objects
struct IdxPair{

    IdxPair(int arg1, int arg2){
        idx1 = min(arg1, arg2);
	idx2 = max(arg1, arg2);
    }

    int idx1, idx2;

    //Utility functions establishing equality and ordering relations
    bool operator == (const IdxPair& pair) const {
        return pair.idx1 == idx1 && pair.idx2 == idx2;
    }

    friend bool operator < (const IdxPair &p1, const IdxPair &p2){
        if(p2.idx2 > p1.idx2) return true;
        else if(p2.idx2 < p1.idx2) return false;
        else if(p2.idx1 > p1.idx1) return true;
        else return false;
    }
};

//Structure for stitching together grains from a disordered large-scale network
struct GrainFusion{

    GrainFusion(){}

    //Endpoints of the edge shared by the two grains
    Point p1, p2;
    //Lists of points on either side of the grain boundary
    set<Point> list1, list2;
};

struct NetData{

    NetData(vector<vector<double>> rls, map<int,vector<vector<double>>> nn, double w, double tk, bool prtct) : rules(rls), nns(nn), width(w), to_keep(tk), protect(prtct) {
        displace = false;
        sdev = 0;
    }
    
    NetData(vector<vector<double>> rls, map<int,vector<vector<double>>> nn, double w, double tk, double dev, bool prtct, bool disp) : rules(rls), nns(nn), width(w), to_keep(tk), sdev(dev), protect(prtct), displace(disp) {}

    vector<vector<double>> rules;
    map<int,vector<vector<double>>> nns;
    double width;
    double to_keep;
    double sdev;
    bool protect;
    bool displace;
};

/*
 * I introduced this struct with the idea that I might use it to develop a
 * parallel version of the code for stitching together different large-scale
 * grains when making a diluted, hierarchical network.
 * TODO: Either implement parallel hierarhical network creation, or replace
 * this superfluous struct.
 */
struct MST_JOB{

    MST_JOB(vector<vector<Edge>> *arg1, int arg2, int arg3, int *arg4, vector<Edge> *arg5, vector<Edge> *arg6, vector<int> *arg7, unordered_map<Point, int> *arg8, int *arg9) : collection(arg1), begin(arg2), end(arg3), edge_count(arg4), pool(arg5), retain(arg6), mst_table(arg7), point_map(arg8), canonical_count(arg9){}

    vector<vector<Edge>> *collection; 
    int begin;
    int end;
    int *edge_count;
    vector<Edge> *pool;
    vector<Edge> *retain; 
    vector<int> *mst_table;
    unordered_map<Point, int> *point_map;
    int *canonical_count;
};

//Utility function to create a random seed
unsigned get_random_seed(){

    unsigned seed;
    FILE *ranfile;
    int num_read = 0;
    time_t the_time;

    ranfile = fopen("/dev/urandom", "r");
    if(ranfile != NULL){
        num_read = fread(&seed, sizeof(unsigned), 1, ranfile);
        fclose(ranfile);
    }

    if(num_read > 0){
        return seed;
    }

    time(&the_time);
    return (unsigned int) the_time;
}

void print_network(NetworkComplex nc, string message){    
    string nextline, filename;
    FILE *out = NULL;
    vector<string> tokens;

    while(out == NULL){
        cout << message;
        getline(cin, nextline);
        tokens = split(nextline, ' ');

        if(tokens.size() > 0){
            filename = tokens[0];
            out = fopen(filename.c_str(), "w");
            if(out == NULL){
                if(! yesno("The file could not be opened. Try again?")){
		    return;
		}
            }
        }

	else{
	    if(! yesno("No file name was read. Try again?")){
                return;
            }
	}
    }

    for(Edge e : nc.edges){
        if(e.idx1 >= nc.points.size() || e.idx2 >= nc.points.size()){
            cerr << "Out of bounds edge: (" << e.idx1 << ", " << e.idx2 << ")\n";
        }
        Point p1 = nc.points[e.idx1];
	Point p2 = nc.points[e.idx2];
        fprintf(out, "%10.8lf %10.8lf \n%10.8lf %10.8lf \n\n", p1.x, p1.y, p2.x, p2.y);
    }

    fclose(out);
}

void print_network_compact(NetworkComplex nc, string message){ 
    string nextline, filename;
    FILE *out = NULL;
    vector<string> tokens;

    while(out == NULL){
        cout << message;
        getline(cin, nextline);
        tokens = split(nextline, ' ');
	if(tokens.size() > 0){
            filename = tokens[0];
            out = fopen(filename.c_str(), "w");
            if(out == NULL){
                if(! yesno("The file could not be opened. Try again?")){
		    return;
		}
            }
	}
	else{
	    if(! yesno("No file name was read. Try again?")){
                return;
            }
        }
    }

    fprintf(out, "%ld %ld\n", nc.points.size(), nc.edges.size());

    for(Point p : nc.points){
        fprintf(out, "%1.10lf\t%1.10lf\n", p.x, p.y);
    }

    for(Edge e : nc.edges){
        fprintf(out, "%u\t%u\n", e.idx1, e.idx2);
    }

    fclose(out);
}

void print_polygons_verbose(FILE *poly_report, PolygonComplex pc){

    for(NetPolygon np : pc.polygons){
        for(int iter = 0; iter < np.vertices.size() - 1; iter++){
            Point p1 = pc.points[np.vertices[iter]];
            Point p2 = pc.points[np.vertices[(iter+1)%np.vertices.size()]];
            fprintf(poly_report, "%lf\t%lf\n%lf\t%lf\n\n", p1.x, p1.y, p2.x, p2.y);
        }
    }

    fclose(poly_report);
}

void print_polygons_concise(FILE *poly_report, PolygonComplex pc){

    fprintf(poly_report, "%d\t%d\n", pc.points.size(), pc.polygons.size());
    for(Point p : pc.points){
        fprintf(poly_report, "%lf\t%lf\n", p.x, p.y);
    }

    for(NetPolygon poly : pc.polygons){
        for(int iter = 0; iter < poly.vertices.size() - 1; iter++){
            fprintf(poly_report, "%d ", poly.vertices[iter]);
        }
        fprintf(poly_report, "%d\n", *(poly.vertices.rbegin()));
    }

    fclose(poly_report);
}

void makeedges(vector<vector<double>> rules, map<int,vector<vector<double>>> nns, vector<double> bnds, NetworkComplex &nc){
    double x, y, x2, y2;
    int ruleiter = 0, rulesize, xiter, base = 0, index, numrules = 0;
    vector<double> rule;
    Point p1, p2;
    map<Point, unsigned int> pmap;
    map<Point, unsigned int>::iterator piter;

    nc.clear();

    for(vector<double> vec : rules){
        numrules += vec.size() - 2;
    }

    y = bnds[1];
    while(y <= bnds[3]){
        xiter = 0;
        rule = rules[ruleiter%rules.size()];
        rulesize = rule.size() - 2;
        x = rule[0] + bnds[0];

        while(x <= bnds[2]){
           p1 = Point(x,y);
           piter = pmap.find(p1);
           if(piter != pmap.end()){
               p1.x = piter->first.x;
               p1.y = piter->first.y;
           }
           else{
               pmap.insert(make_pair(p1, pmap.size()));
               nc.points.push_back(p1);
           }

           index = base + xiter % rulesize;
           for(vector<double> nextnn : nns[index]){
               x2 = x + nextnn[0];
               y2 = y + nextnn[1];
               if(x2>=bnds[0] && x2<=bnds[2] && y2>=bnds[1] && y2<=bnds[3]){
                   
                   p2 = Point(x2, y2);
                   piter = pmap.find(p2);
                   if(piter != pmap.end()){
                       p2.x = piter->first.x;
                       p2.y = piter->first.y;
                   }
                   else{
                       pmap.insert(make_pair(p2, pmap.size()));
		       nc.points.push_back(p2);
                   }

                   nc.edges.push_back(Edge(pmap[p1], pmap[p2]));
               }
           }
           x += rule[xiter++%rulesize + 2];
        }
        y += rule[1];
        base = (base + rulesize) % numrules;
        ruleiter++;
    }
}

double get_v_offset(vector<vector<double>> rules){
    double voffset = 0;

    for(vector<double> rule : rules){
        voffset += rule[1];
    }

    return voffset;
}

double get_h_offset(vector<vector<double>> rules){
    double hoffset = 0;
    int index;

    for(index = 2; index < rules[0].size(); index++){
        hoffset += rules[0][index];
    }

    return hoffset;
}

set<Edge> top_bottom(NetworkComplex nc){
    Point p1, p2, bleft, bright, tleft, tright, curr;
    set<Edge> topbottom;
    vector<Edge> perimeter;
    map<Point,Edge> neighbor_map;
    Edge next;
    double midx, midy, dx, dy, slope, xcent, ycent, miny, maxy;
    bool onborder;
    int iter;

    xcent = ycent = 0;
    miny = INF;
    maxy = -INF;

    bleft = Point(INF,INF);
    bright = Point(-INF,INF);
    tleft = Point(INF,-INF);
    tright = Point(-INF,-INF);

    for(Point p : nc.points){
        xcent += p.x;
        ycent += p.y;
    }
    xcent /= nc.points.size();
    ycent /= nc.points.size();

    for(Edge e : nc.edges){
        p1 = nc.points[e.idx1];
        p2 = nc.points[e.idx2];
        midx = (p1.x + p2.x)/2;
        midy = (p1.y + p2.y)/2;
        dx = midx - xcent;
        dy = midy - ycent;

        onborder = true;
        iter = 0;
        do{
            next = nc.edges[iter];
            if(!(next == e)){
                onborder = !intersection(Point(midx,midy), Point(midx+dx,midy+dy), p1, p2);
            }
            iter ++;
        }while(onborder && iter < nc.edges.size());
        if(onborder) perimeter.push_back(e);
    }

    for(Edge e : perimeter){
        p1 = nc.points[e.idx1];
        p2 = nc.points[e.idx2];

        if(p1.y <= miny){
            if(p1.y < miny) miny = p1.y;
            if(p1.x < bleft.x) bleft = p1;
            if(p1.x > bright.x) bright = p1;
        }
        if(p1.y >= maxy){
            if(p1.y > maxy){
                maxy = p1.y;
                tleft = p1;
                tright = p2;
            }
            if(p1.x < tleft.x) tleft = p1;
            if(p1.x > tright.x) tright = p1;
        }
        
        if(p2.y <= miny){
            if(p2.y < miny) miny = p2.y;
            if(p2.x < bleft.x) bleft = p2;
            if(p2.x > bright.x) bright = p2;
        }
        if(p2.y >= maxy){
            if(p2.y > maxy){
                maxy = p2.y;
                tleft = p2;
                tright = p2;
            }
            if(p2.x < tleft.x) tleft = p2;
            if(p2.x > tright.x) tright = p2;
        }
        
        if(neighbor_map.find(p1) == neighbor_map.end()){
            neighbor_map.insert(make_pair(p1,e));
        }
        else if(p2.x > nc.points[neighbor_map[p1].idx2].x) neighbor_map[p1] = e;
    }

    curr = bleft;
    while(!(curr == bright)){
        next = neighbor_map[curr];
        topbottom.insert(next);
        curr = nc.points[next.idx2];
    }
    curr = tleft;
    while(!(curr == tright)){
        next = neighbor_map[curr];
        topbottom.insert(next);
        curr = nc.points[next.idx2];
    }
    return topbottom;
}

//After edges have been removed from a network, cull unused points and reassign
//edge endpoint indices
void reassign_points_edges(NetworkComplex &nc, vector<Edge> kept_edges){

    set<int> kept_indices;
    map<int, int> reassignments;
    vector<Point> kept_points;

    //Figure out which points are retained
    for(Edge e : kept_edges){
        kept_indices.insert(e.idx1);
        kept_indices.insert(e.idx2);
    }

    //Make a map from indices of retained points to their new indices
    for(int point : kept_indices){
        reassignments.insert(make_pair(point, reassignments.size()));
	kept_points.push_back(nc.points[point]);
    }   

    //Reassign edges' endpoint indices based on updated point indices
    for(int idx = 0; idx < kept_edges.size(); idx++){
        kept_edges[idx].idx1 = reassignments[kept_edges[idx].idx1];
        kept_edges[idx].idx2 = reassignments[kept_edges[idx].idx2];
    }

    nc.points = kept_points;
    nc.edges = kept_edges;
}

void true_random(NetworkComplex &nc, double to_keep){
    int num_needed, index;
    vector<Edge> kept;
    unsigned seed = get_random_seed();
    mt19937 gen(seed);

    //Randomly retain a subset of all edges
    num_needed = (int) (to_keep * nc.edges.size());
    shuffle(nc.edges.begin(), nc.edges.end(), gen);

    index = 0;
    do{
        kept.push_back(nc.edges[index++]);
    }while(kept.size() < num_needed);

    reassign_points_edges(nc, kept);
}

/*vector<int> network_groups(vector<Edge> edges, set<Point> points){
    vector<int> canonical(points.size(), -1);
    vector<int> results(points.size(), 0);
    int root1, root2, max_zero;

    for(Edge e : edges){
        root1 = find_root(canonical, points.find(e.p1)->index);
        root2 = find_root(canonical, points.find(e.p2)->index);
        if(root1 != root2) makeunion(canonical, root1, root2);
    }
    
    for(int i = 0; i < points.size(); i ++){
        results[find_root(canonical,i)]++;
    }

    sort(results.begin(), results.end());
    for(int i = 0; i < results.size(); i++){
        if(!results[i]) max_zero = i;
        else break;
    }
    results.erase(results.begin(), results.begin() + max_zero + 1);
    return results;
}*/

void getangles(vector<double> angvec, double angle, double &low, double &high){
    int index = 0, size = angvec.size();

    while(index < size-1 && angle-angvec[index] > FLOAT_TOL) index++;
    
    low = angvec[(index+size-1)%size];
    high = angvec[(index+size+1)%size];
}

void changes(double low, double high, double hwidth, double &dx, double &dy){

    double diff = (high-low)/2;
    if(diff < 0) diff += M_PI;

    if(diff == 0){
        fprintf(stderr, "Illegal argument to function changes.\n");
        cerr << "Low: " << low << " High: " << high << "\n";
        return;
    }

    dx = hwidth * (cos(low)/tan(diff) - sin(low));
    dy = hwidth * (sin(low)/tan(diff) + cos(low));
}

void flush_with_edge(Point &pnt, double slope, double y_ext){
    double x_old, y_old, x_new;

    x_old = pnt.x;
    y_old = pnt.y;
    x_new = slope != 0 ? (y_ext - y_old) / slope + x_old : x_old;

    pnt.x = x_new;
    pnt.y = y_ext;
}

void add_if_missing(Point p, map<Point, int> &pmap, vector<Point> &plist){
    if(pmap.find(p) == pmap.end()){
        pmap.insert(make_pair(p, pmap.size()));
	plist.push_back(p);
    }
}

void add_thickness(NetworkComplex &current, double thickness, PolygonComplex &pc, bool makepoly, bool level){

    map<Point, vector<double>> angmap;
    map<Point, int> point_map, poly_pmap;
    Point p1, p2, key;
    Point p1f, p2f, p3f, p4f, p1fb, p2fb, p3fb, p4fb;
    double low, high, ang1, ang2, dx, dy, hwidth, slope;
    double ymin = FLT_MAX, ymax = FLT_MIN, ylow, yhigh, y_ext;
    bool p1fflag, p2fflag, p3fflag, p4fflag, p1_is_end, p2_is_end;
    vector<Point> replace_points;
    vector<Edge> replace_edges;
    hwidth = thickness/2;

    for(Edge e : current.edges){
        p1 = current.points[e.idx1];
	p2 = current.points[e.idx2];

	if(p1.y < ymin) ymin = p1.y;
	if(p2.y < ymin) ymin = p2.y;
	if(p1.y > ymax) ymax = p1.y;
	if(p2.y > ymax) ymax = p2.y;

	ang1 = atan2(p2.y - p1.y, p2.x - p1.x);
	if(ang1 < 0) ang1 += 2*M_PI;
	ang2 = ang1 < M_PI ? ang1 + M_PI : ang1 - M_PI;

	if(angmap.find(p1) == angmap.end()){
            angmap.insert(make_pair(p1, vector<double>()));
        }

	if(angmap.find(p2) == angmap.end()){
            angmap.insert(make_pair(p2, vector<double>()));
        }

	angmap[p1].push_back(ang1);
	angmap[p2].push_back(ang2);
    }

    ylow = ymin - hwidth;
    yhigh = ymax + hwidth;

    for(auto iter = angmap.begin(); iter != angmap.end(); iter ++){
        key = iter->first;
	sort(angmap[key].begin(), angmap[key].end());
    }

    for(Edge e : current.edges){
        p1 = current.points[e.idx1];
	p2 = current.points[e.idx2];

        ang1 = atan2(p2.y - p1.y, p2.x - p1.x);
	if(ang1 < 0) ang1 += 2*M_PI;
	ang2 = ang1 < M_PI ? ang1 + M_PI : ang1 - M_PI;

	slope = abs(p2.x - p1.x) > FLOAT_TOL ? (p2.y - p1.y)/(p2.x - p1.x) : BIG_SLOPE;
	p1fflag = false;
        p2fflag = false;
        p3fflag = false;
        p4fflag = false;
        p1_is_end = false;
        p2_is_end = false;

	getangles(angmap[p1], ang1, low, high);
        if(abs(ang1 - low) < FLOAT_TOL){
            p1_is_end = true;
            p1f = Point(p1.x + hwidth*sin(ang1), p1.y - hwidth*cos(ang1));
            p3f = Point(p1.x - hwidth*sin(ang1), p1.y + hwidth*cos(ang1));

            if((abs(p1.y - ymin) < FLOAT_TOL || abs(p1.y - ymax) < FLOAT_TOL) && abs(slope) >= FLOAT_TOL){
                y_ext = abs(p1.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p1f, slope, y_ext);
                flush_with_edge(p3f, slope, y_ext);
            }

            add_if_missing(p1f, point_map, replace_points);
            add_if_missing(p3f, point_map, replace_points);

            replace_edges.push_back(Edge(point_map[p1f], point_map[p3f]));
        }
        else{
            changes(low, ang1, hwidth, dx, dy);
            p1f = Point(p1.x + dx, p1.y + dy);
            changes(ang1, high, hwidth, dx, dy);
            p3f = Point(p1.x + dx, p1.y + dy);

            if((p1f.y < ylow - FLOAT_TOL || p1f.y > yhigh + FLOAT_TOL) && level){
                y_ext = abs(p1.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p1f, slope, y_ext);
                add_if_missing(p1f, point_map, replace_points);
                p1fb = Point(p1.x, y_ext);
                add_if_missing(p1fb, point_map, replace_points);
                p1fflag = true;
            }
	    else add_if_missing(p1f, point_map, replace_points);
            if((p3f.y < ylow - FLOAT_TOL || p3f.y > yhigh + FLOAT_TOL) && level){
                y_ext = abs(p1.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p3f, slope, y_ext);
                add_if_missing(p3f, point_map, replace_points);
                p3fb = Point(p1.x, y_ext);
                add_if_missing(p3fb, point_map, replace_points);
                p3fflag = true;
            }
	    else add_if_missing(p3f, point_map, replace_points);
        }

	getangles(angmap[p2], ang2, low, high);
        if(abs(ang2 - low) < FLOAT_TOL){
            p2_is_end = true;
            p2f = Point(p2.x + hwidth*sin(ang1), p2.y - hwidth*cos(ang1));
            p4f = Point(p2.x - hwidth*sin(ang1), p2.y + hwidth*cos(ang1));

            if((abs(p2.y - ymin) < FLOAT_TOL || abs(p2.y - ymax) < FLOAT_TOL) && abs(slope) >= FLOAT_TOL){
                y_ext = abs(p2.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p2f, slope, y_ext);
                flush_with_edge(p4f, slope, y_ext);
            }

            add_if_missing(p2f, point_map, replace_points);
            add_if_missing(p4f, point_map, replace_points);

            replace_edges.push_back(Edge(point_map[p2f], point_map[p4f]));
        }
        else{
            changes(ang2, high, hwidth, dx, dy);
            p2f = Point(p2.x + dx, p2.y + dy);
            changes(low, ang2, hwidth, dx, dy);
            p4f = Point(p2.x + dx, p2.y + dy);

            if((p2f.y < ylow - FLOAT_TOL || p2f.y > yhigh + FLOAT_TOL) && level){
                y_ext = abs(p2.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p2f, slope, y_ext);
                add_if_missing(p2f, point_map, replace_points);
                p2fb = Point(p2.x, y_ext);
                add_if_missing(p2fb, point_map, replace_points);
                p2fflag = true;
            }
	    else add_if_missing(p2f, point_map, replace_points);
            if((p4f.y < ylow - FLOAT_TOL || p4f.y > yhigh + FLOAT_TOL) && level){
                //cerr << "Current p4f: " << p4f.x << "\t" << p4f.y << "\n";
                y_ext = abs(p2.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p4f, slope, y_ext);
                add_if_missing(p4f, point_map, replace_points);
                p4fb = Point(p2.x, y_ext);
                add_if_missing(p4fb, point_map, replace_points);
                p4fflag = true;
                //cerr << "New p4f: " << p4f.x << "\t" << p4f.y << "\n\n";
            }
	    else add_if_missing(p4f, point_map, replace_points);
        }

        if(p1fflag) replace_edges.push_back(Edge(point_map[p1fb], point_map[p1f]));
        replace_edges.push_back(Edge(point_map[p1f], point_map[p2f]));
        if(p2fflag) replace_edges.push_back(Edge(point_map[p2f], point_map[p2fb]));
        if(p3fflag) replace_edges.push_back(Edge(point_map[p3fb], point_map[p3f]));
        replace_edges.push_back(Edge(point_map[p3f],point_map[p4f]));
        if(p4fflag) replace_edges.push_back(Edge(point_map[p4f], point_map[p4fb]));

        if(makepoly){

            //Ensure necessary points are added to the polycomplex object
            add_if_missing(p1, poly_pmap, pc.points);
            add_if_missing(p2, poly_pmap, pc.points);
            add_if_missing(p1f, poly_pmap, pc.points);
            add_if_missing(p2f, poly_pmap, pc.points);
            add_if_missing(p3f, poly_pmap, pc.points);
            add_if_missing(p4f, poly_pmap, pc.points);
	    if(p1fflag) add_if_missing(p1fb, poly_pmap, pc.points);
	    if(p2fflag) add_if_missing(p2fb, poly_pmap, pc.points);
	    if(p3fflag) add_if_missing(p3fb, poly_pmap, pc.points);
	    if(p4fflag) add_if_missing(p4fb, poly_pmap, pc.points);

            vector<int> vertices({poly_pmap[p1f], poly_pmap[p2f]});

            if(!p2_is_end){
                if(p2fflag){
		    vertices.push_back(poly_pmap[p2fb]);
                }
                vertices.push_back(poly_pmap[p2]);
            
                if(p4fflag){
		    vertices.push_back(poly_pmap[p4fb]);
                }
            }
            vertices.push_back(poly_pmap[p4f]);
            vertices.push_back(poly_pmap[p3f]);

	    if(! p1_is_end){
                if(p3fflag){
                    vertices.push_back(poly_pmap[p3fb]);
                }
                vertices.push_back(poly_pmap[p1]);
                if(p1fflag){
                    vertices.push_back(poly_pmap[p1fb]);
                }
            }

            //Make sure the vertices are traversed in CCW order
	    p1f = pc.points[vertices[0]];
	    p2f = pc.points[vertices[1]];
	    p3f = pc.points[vertices[2]];
	    ang1 = atan2(p1f.y - p2f.y, p1f.x - p2f.x);
	    if(ang1 < 0) ang1 += 2*M_PI;
	    ang2 = atan2(p3f.y - p2f.y, p3f.x - p2f.x);
	    if(ang2 < 0) ang2 += 2*M_PI;
	    if(ang1 < ang2) reverse(vertices.begin(), vertices.end());

	    //Add the polygon and its edges to the data structure describing
	    //the polygonal tiling at the large length scale.
            pc.polygons.push_back(NetPolygon(pc.points, vertices, pc.polygons.size()));
	    for(int idx = 0; idx < vertices.size(); idx++){
                pc.edges.insert(Edge(vertices[idx], vertices[(idx+1)%vertices.size()]));
            }
        }

    }

    current.clear();
    current.points.insert(current.points.begin(), replace_points.begin(), replace_points.end());
    current.edges.insert(current.edges.begin(), replace_edges.begin(), replace_edges.end());
}

vector<vector<double>> getrules(){
    vector<double> rule;
    vector<vector<double>> rules;
    bool newRule;

    while(true){
        newRule = yesno("Enter a new rule? ");
        if(!newRule){
            if(rules.size()){
                break;
            }
            else{
                fprintf(stderr, "Enter at least one rule.\n");
            }
        }
        else{
            rule = getdoubles("Enter rule: ");
            if(rule.size() < 3){
                fprintf(stderr, "Enter at least three numbers.\n");
            }
            else{
                rules.push_back(rule);
            }
        }
    }

    return rules;
}

map<int,vector<vector<double>>> getnns(vector<vector<double>> rules){
    map<int,vector<vector<double>>> nns;
    vector<double> nn;
    int ruleiter = 0, count = 0, pointiter;
    string prompt, base;

    base = "Add a nearest neighbor vector for rule ";

    for(vector<double> nextrule : rules){
       ruleiter ++;
       for(pointiter = 1; pointiter <= nextrule.size() - 2; pointiter ++){
           vector<vector<double>> nextset;
           while(true){
               prompt = base + to_string(ruleiter) + " point " + to_string(pointiter) + "?";
               if(!yesno(prompt)) break;
               nn = getdoubles("Enter the vector: ");
               if(nn.size() != 2) fprintf(stderr, "Enter two numbers.\n");
               else{
                   nextset.push_back(nn);
               }
           }
           nns.insert(pair<int, vector<vector<double>>>(count++, nextset));
       }
    }

    return nns;
}

void scale_vector(vector<double>& in, double scale){
    int index;
    for(index = 0; index < in.size(); index++){
        in.at(index) = in.at(index)*scale;
    }
}

bool import_lattice(vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns, double scale){
    ifstream latfile;
    string name, nextline;
    bool again, blank, fileopen = false;
    vector<double> rule, nn;
    int lcount = 0, pointiter, nncount = 0, ruleiter = 0;

    //Prompt for file name
    do{
        printf("Enter the lattice file name: ");
        getline(cin, nextline);
        name = split(nextline, ' ')[0];

        if(!name.empty()) latfile.open(name);

        if(!latfile.is_open()){
            again = yesno("No file was read. Try again? ");
            if(!again) return false;
        }
    }while(!latfile.is_open());

    //Process rules until a blank line is reached
    blank = false;
    while(!latfile.eof()){
        lcount ++;
        getline(latfile, nextline);
        if(nextline.empty()) break;

        rule = parse_doubles(split(nextline, ' '));
        if(rule.size() < 3){
            fprintf(stderr, "Too few numbers in rule on line line %d\n",lcount);
            latfile.close();
            return false;
        }
        scale_vector(rule, scale);
        rules.push_back(rule);
    }
    
    //Read nearest neighbor rules
    for(vector<double> nextrule : rules){
       ruleiter ++;
       for(pointiter = 1; pointiter <= nextrule.size() - 2; pointiter ++){
           vector<vector<double>> nextset;
           while(!latfile.eof()){
               lcount ++;
               getline(latfile, nextline);
               if(nextline.empty()) break;
               nn = parse_doubles(split(nextline, ' '));
               if(nn.size() != 2) fprintf(stderr, "Insufficient information for nearest neighbor rule on line %d.\n", lcount);
               else{
                   scale_vector(nn, scale);
                   nextset.push_back(nn);
               }
           }
           nns.insert(pair<int, vector<vector<double>>>(nncount++, nextset));
       }
    }


    latfile.close();
    return true;
}

void get_lattice_info(vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns){
    bool success;
    double scale;

    if(yesno("Read lattice data from file?")){
        do{
            scale = getdoubles("Enter the scale factor: ").at(0);
            success = import_lattice(rules, nns, scale);
            if(!success) if(!yesno("Read failed. Try again?")) break;
        }while(!success);
    }

    if(!success){
        rules = getrules();
        nns = getnns(rules);
    }
}

void loadNetStack(stack<NetData>& net_stack){
    //vector<vector<double>> rules;
    //map<int,vector<vector<double>>> nns;
    double width, to_keep;
    bool valid, protect, displace;
    vector<double> response;
    double sdev;

    do{
        vector<vector<double>> rules;
        map<int,vector<vector<double>>> nns;
        get_lattice_info(rules, nns);
        valid = true;
        do{
            response = getdoubles("Enter the bond width: ");
            if(!response.size() || response.at(0) < 0){
                fprintf(stderr,"Enter a non-negative number.\n");
                valid = false;
            }
            else if(response.at(0) == 0 && net_stack.size()){
                fprintf(stderr, "The width must be greater than 0.\n");
                valid = false;
            }
            else{
                width = response.at(0);
                valid = true;
            }
        }while(!valid);
       
        valid = true;
        do{
            response = getdoubles("Enter portion of bonds to keep: ");
            if(!response.size() || response.at(0) < 0){
                fprintf(stderr, "Enter a non-negative number.\n");
                valid = false;
            }
            else{
                to_keep = response.at(0);
                valid = true;
            }
        }while(!valid); 
        protect = yesno("Protect top and bottom edges?");

        displace = yesno("Displace edges?");
        if(displace){
            do{
                response = getdoubles("Enter displacement standard deviation: ");
                if(response.size() == 0 || response[0] < 0){
                    cerr << "Enter a non-negative number.\n";
                    valid = false;
                }
                else{
                    sdev = response[0];
                    valid = true;
                }
            }while(! valid);
        }

        if(!displace) net_stack.push(NetData(rules, nns, width, to_keep, protect));
        else net_stack.push(NetData(rules, nns, width, to_keep, sdev, protect, displace));

    }while(yesno("Add another layer?"));
}

//Given a set of points and edges connecting those points, partition edges into
//two subsets: a random set forming a minimum spanning tree, and all others.
void randomMST(NetworkComplex nc, vector<Edge> &keep, vector<Edge>& rejects){
    int index = 0, numkept = 0, numneeded;
    int root1, root2;
    Edge next;
    unsigned seed = get_random_seed();
    mt19937 gen(seed);

    vector<int> union_find_table(nc.points.size(), -1);
    numneeded = nc.points.size() - 1;

    shuffle(nc.edges.begin(), nc.edges.end(), gen);
    while(!nc.edges.empty()){
        next = nc.edges[0];
        nc.edges.erase(nc.edges.begin());
        root1 = find_root(union_find_table, next.idx1);
        root2 = find_root(union_find_table, next.idx2);
        if(root1 != root2){
            makeunion(union_find_table, root1, root2);
            keep.push_back(next);
            if(++numkept == numneeded) break;
        }
        else rejects.push_back(next);
    }

    for(Edge next_edge : nc.edges){
        rejects.push_back(next_edge);
    }
}

//This function is similar to the one above, but considers the case in which
//the edges are a proper subset of all edges in a network.
void subsetRandomMST(vector<Edge> in, vector<Edge> &keep, set<Edge>& rejects){
    int numkept = 0, numneeded, root1, root2;
    Edge next;
    map<int, int> pmap;
    unsigned seed = get_random_seed();
    mt19937 gen(seed);

    for(Edge e : in){
        if(pmap.find(e.idx1) == pmap.end()){
            pmap.insert(make_pair(e.idx1, pmap.size()));
        }
        if(pmap.find(e.idx2) == pmap.end()){
            pmap.insert(make_pair(e.idx2, pmap.size()));
        }
    }

    vector<int> union_find_table(pmap.size(), -1);
    numneeded = pmap.size() - 1;

    shuffle(in.begin(), in.end(), gen);
    while(!in.empty()){
        next = in[0];
        in.erase(in.begin());
        root1 = find_root(union_find_table, pmap[next.idx1]);
        root2 = find_root(union_find_table, pmap[next.idx2]);
        if(root1 != root2){
            makeunion(union_find_table, root1, root2);
            keep.push_back(next);
            if(++numkept == numneeded) break;
        }
        else rejects.insert(next);
    }

    for(Edge next_edge : in){
        rejects.insert(next_edge);
    }
}

NetworkComplex random_connected(NetworkComplex in, double to_keep){
    int num_needed, reject_index;
    vector<Edge> kept_edges, rejects;
    unsigned seed = get_random_seed();
    mt19937 gen(seed);
    set<int> kept_indices;
    map<int, int> pmap;
    vector<Point> kept_points;

    num_needed = (int) (to_keep * in.edges.size());
    randomMST(in, kept_edges, rejects);
    shuffle(rejects.begin(), rejects.end(), gen);

    reject_index = 0;
    do{
        kept_edges.push_back(rejects[reject_index++]);
    }while(reject_index < rejects.size() && kept_edges.size() < num_needed);

    //Identify kept points, and reassign indices in each kept edge
    for(Edge e : kept_edges){
        kept_indices.insert(e.idx1);
        kept_indices.insert(e.idx2);
    }
    for(auto iter = kept_indices.begin(); iter != kept_indices.end(); iter++){
        pmap.insert(make_pair(*iter, pmap.size()));
	kept_points.push_back(in.points[*iter]);
    }
    for(Edge e : kept_edges){
        e.idx1 = pmap[e.idx1];
        e.idx2 = pmap[e.idx2];
    }

    return NetworkComplex(kept_points, kept_edges);
}

/*
This function adds Gaussian random noise to the location of each point.
The function takes as arguments the original set of points and edges describing
the network, and a standard deviation for Gaussian random noise. Points are 
shifted, and edges are updated accordingly.
*/
void displace_points_grn(NetworkComplex &nc, double sdev){

    unsigned seed = get_random_seed();
    mt19937 gen(seed);
    normal_distribution<double> dist(0, sdev);

    for(int iter = 0; iter < nc.points.size(); iter++){
        nc.points[iter].x += dist(gen);
	nc.points[iter].y += dist(gen);
    }
}

//Given a set of vectors to nearest neighbors, find the smallest length
//among these vectors
double get_min_dist(map<int, vector<vector<double>>> nns){

    double dist_sq, min_dist_sq;
    min_dist_sq = FLT_MAX;

    for(auto iter = nns.begin(); iter != nns.end(); iter ++){
        for(vector<double> disp : iter->second){
            dist_sq = disp[0]*disp[0] + disp[1]*disp[1];
            if(dist_sq < min_dist_sq) min_dist_sq = dist_sq;
        }
    }

    return sqrt(min_dist_sq);
}

/*After the edges at one length scale that fall within the edges of the next
larger length scale have been determined, this function stitches together the 
network at the smaller length scale, given a knowledge of edges at the small 
length scale contained within a polygonal tiling of edges at the large length 
scale.
*/
vector<Edge> stitch_network(vector<vector<Edge>> collection, map<int, vector<Edge>> adj_map, set<int> in_points, double to_keep){

    int canonical_count, root1, root2;
    int numNeeded, discardIter;
    unordered_map<int, int> point_map;
    vector<Edge> pool, next_al, retain, edgeMST, discard;
    set<Edge> rejects;
    Edge front;
    vector<Edge>::iterator pool_iter;
    unsigned seed = get_random_seed();
    mt19937 gen(seed);

    canonical_count = in_points.size();

    for(int next_point : in_points){
        point_map.insert(make_pair(next_point, point_map.size()));
    }

    vector<int> mst_table(point_map.size(), -1);

    for(vector<Edge> next_set : collection){
        edgeMST.clear();
        subsetRandomMST(next_set, edgeMST, rejects);
        for(Edge next_edge : edgeMST){
            retain.push_back(next_edge);
            root1 = find_root(mst_table, point_map[next_edge.idx1]);
            root2 = find_root(mst_table, point_map[next_edge.idx2]);
            if(root1 != root2){
                makeunion(mst_table, root1, root2);
                canonical_count --;
            }
        }
    }

    for(auto it = adj_map.begin(); it != adj_map.end(); it++){
        next_al = it->second;

        shuffle(next_al.begin(), next_al.end(), gen);
        retain.push_back(next_al[0]);

        if(next_al.size() > 1){
            for(auto al_iter = next_al.begin()+1; al_iter != next_al.end(); al_iter++){
                rejects.insert(*al_iter);
            }
	}
    }

    numNeeded = (int) ((retain.size() + rejects.size()) * to_keep);
    pool.insert(pool.begin(), rejects.begin(), rejects.end());
    shuffle(pool.begin(), pool.end(), gen);
    for(pool_iter = pool.begin(); pool_iter != pool.end(); pool_iter ++){
        if(canonical_count == 1) break;
        front = *pool_iter;
        root1 = find_root(mst_table, point_map[front.idx1]);
        root2 = find_root(mst_table, point_map[front.idx2]);
        if(root1 != root2){
            makeunion(mst_table, root1, root2);
            canonical_count --;
            retain.push_back(front);
        }
        else discard.push_back(front);
    }

    discard.insert(discard.end(), pool_iter, pool.end());

    shuffle(discard.begin(), discard.end(), gen);

    discardIter = 0;
    while(retain.size() < numNeeded && discardIter < discard.size()){
        retain.push_back(discard[discardIter ++]);
    }

    return retain;
}

//Given a subset of all edges in a network, find out which points are retained,
//and construct a new network complex with the reduced point set and edges
//with relabled endpoint indices.
NetworkComplex reduced_network(vector<Point> all_points, vector<Edge> kept_edges){

    NetworkComplex reduced;
    map<int, int> pmap;

    for(Edge e : kept_edges){
        if(pmap.find(e.idx1) == pmap.end()){
            pmap.insert(make_pair(e.idx1, pmap.size()));
	    reduced.points.push_back(all_points[e.idx1]);
        }
        if(pmap.find(e.idx2) == pmap.end()){
            pmap.insert(make_pair(e.idx2, pmap.size()));
	    reduced.points.push_back(all_points[e.idx2]);
        }
    }

    for(Edge e : kept_edges){
        reduced.edges.push_back(Edge(pmap[e.idx1], pmap[e.idx2]));
    }

    return reduced;
}

//Determine whether a query edge intersects an interior edge of a polygon, i.e.,
//an edge that polygon shares with another polygon.
bool interior_crossing(vector<Point> pts, Edge e, PolygonComplex pc, int idx){

    Point p1, p2, p3, p4;
    int idx3, idx4, num_verts;
   
    num_verts = pc.polygons[idx].vertices.size();
    p1 = pts[e.idx1];
    p2 = pts[e.idx2];

    for(int iter = 0; iter < num_verts; iter++){
        idx3 = pc.polygons[idx].vertices[iter];
	idx4 = pc.polygons[idx].vertices[(iter+1) % num_verts];
        p3 = pc.points[idx3];
        p4 = pc.points[idx4];
        if(intersection(p1, p2, p3, p4)){
            if(pc.edges.find(Edge(idx4, idx3)) != pc.edges.end()) return true;
        }
    }

    return false;
}

/*
Helper function for sieve_edges that calculates indices for use as keys
to map pairs of adjacent polygons from a tiling of a network to edges straddling
those polygons
*/
int adj_map_index(int index1, int index2, int total){
    int low = index1 <= index2 ? index1 : index2;
    int high = low == index1 ? index2 : index1;

    return low * (total + total + 1 - low) / 2 + high;
}

/*This function is called when the edges describing the network at a given
length scale have been given non-zero width and stored in the structure top,
and the resulting geometrical object has been broken into a tiling of polygons
represented by the members of the structure pc. The edges describing the 
network at the next lowest length scale, stored in the structure bottom, are
examined to determine whether they are contained entirely in one of these
polygons. Those edges meeting this criterion are kept and returned in the
structure retain.
*/
NetworkComplex sieve_edges(NetworkComplex top, double length, NetworkComplex bottom, PolygonComplex pc, double to_keep){

    vector<Edge> kept_edges;
    vector<vector<Edge>> edge_collection;
    map<int, vector<Edge>> adj_map;
    map<int, vector<int>> pmap;
    set<int> in_points;
    int adj_idx;
    bool included;

    //Index polygons to enable rapid tests of membership
    shared_ptr<HNetBVH> hnbvh = make_shared<HNetBVH>(pc.points, pc.polygons);

    for(int i = 0; i < pc.polygons.size(); i ++){
        edge_collection.push_back(vector<Edge>());
    }

    //Make a record of whether a point is contained within any large-scale
    //tiles, and if so, which ones.
    for(int pindex = 0; pindex < bottom.points.size(); pindex++){
        vector<int> hits = hnbvh->within(bottom.points[pindex]);
	if(hits.size() > 0){
            sort(hits.begin(), hits.end());
            pmap.insert(make_pair(pindex, hits));
        }
    }

    //Check whether each edge is either wholly in a polygonal tile, or
    //contained in to adjacent polygonal tiles.
    for(Edge e : bottom.edges){

        included = false;

        if(pmap.find(e.idx1) != pmap.end() && pmap.find(e.idx2) != pmap.end()){
            vector<int> hits1 = pmap[e.idx1];
            vector<int> hits2 = pmap[e.idx2];

	    for(int poly1 : hits1){
                for(int poly2 : hits2){

                    //If the indices are equal, add a record of this edge to
		    //the edges within a given polygonal tile.
		    if(poly1 == poly2){
                        if(to_keep < 1){
                            edge_collection[poly1].push_back(e);
                            in_points.insert(e.idx1);
                            in_points.insert(e.idx2);
                        }
			else kept_edges.push_back(e);
			included = true;
			break;
                    }

		    //If the indices are different, determine whether the edge
		    //straddles an interior edge shared between two tiles. Only
		    //in this case should the edge be kept.
		    else{
                        if(interior_crossing(bottom.points, e, pc, poly1)){
                            if(to_keep < 1){
                                adj_idx = adj_map_index(poly1, poly2, pc.polygons.size());
                                if(adj_map.find(adj_idx) == adj_map.end()){
                                    adj_map.insert(make_pair(adj_idx, vector<Edge>()));
                                }
                                adj_map[adj_idx].push_back(e);
			        in_points.insert(e.idx1);
			        in_points.insert(e.idx2);
                            }
			    else kept_edges.push_back(e);
			    included = true;
			    break;
                        }
                    }
                }
		if(included) break;
            }
        }
    }

    if(to_keep < 1){
        kept_edges = stitch_network(edge_collection, adj_map, in_points, to_keep);
    }

    return reduced_network(bottom.points, kept_edges);
}

//Given a set of polygonal tiles, find all edges common to two tiles and
//arrange them in a map for the purpose of stitching together networks
//contained in different tiles
map<IdxPair, GrainFusion> get_border_map(PolygonComplex pc){
    map<IdxPair, GrainFusion> border_map;
    map<IdxPair, vector<int>> neighbor_map;

    for(int piter = 0; piter < pc.polygons.size(); piter ++){
        vector<int> vertices = pc.polygons[piter].vertices;
        for(int viter = 0; viter < vertices.size(); viter ++){
            IdxPair pair = IdxPair(vertices[viter], vertices[(viter+1)%vertices.size()]);
	    if(neighbor_map.find(pair) == neighbor_map.end()){
                neighbor_map.insert(make_pair(pair, vector<int>()));
            }
	    neighbor_map[pair].push_back(piter);
        }
    }

    for(auto iter = neighbor_map.begin(); iter != neighbor_map.end(); iter++){
        if(iter->second.size() < 2) continue;
        GrainFusion gf;
	gf.p1 = pc.points[iter->first.idx1];
	gf.p2 = pc.points[iter->first.idx2];
        border_map.insert(make_pair(IdxPair(iter->second[0],iter->second[1]), gf));
    }

    return border_map;
}

/*
Given a set of tiles describing a network with disorder in node position,
fill tiles with small-scale networks, then stitch tiles together.
*/
void make_edges_deformed(vector<vector<double>> rules, map<int, vector<vector<double>>> nns, NetworkComplex top, PolygonComplex pc, double width, double to_keep, double cutoff_sq, NetworkComplex &final_network){

    double hwidth = width / 2;
    //Map edges defining grain boundaries to points on either side of boundaries
    map<IdxPair, GrainFusion> border_map;

    //Bounds for creating a set of points and edges within a grain
    double h_offset = get_h_offset(rules);
    double starting_bounds[] = {-h_offset,0,0,width};
    vector<double> bounds(starting_bounds, starting_bounds + 4);
    double length, angle, dx, dy, distance, min;
    Point midpoint, p1, p2, p3, p4, closest, pivot;
    int index, pt_index, total_edges, num_needed;
    vector<Edge> retain, mst;
    vector<vector<Edge>> grain_lists;
    set<Edge> pool;
    map<int, vector<Edge>> adj_map;
    double left, right, low, high;
    map<Point, int> pmap;
    vector<Point> points;
    unsigned seed = get_random_seed();
    mt19937 gen(seed);
    bool contains_1, contains_2, no_cross;

    //Find borders between neighboring tiles
    border_map = get_border_map(pc);

    //Create an index of polygons to support rapid point-in-polygon queries
    shared_ptr<HNetBVH> hnbvh = make_shared<HNetBVH>(pc.points, pc.polygons);

    if(to_keep < 1){
        for(int iter = 0; iter < pc.polygons.size(); iter++){
            grain_lists.push_back(vector<Edge>());
        }
    }

    //Create a set of nodes and edges for each tile, retaining interior points 
    //from edges that cross borders with neighboring tiles
    for(index = 0; index < top.edges.size(); index++){
        p1 = top.points[top.edges[index].idx1];
        p2 = top.points[top.edges[index].idx2];
        angle = atan2(p2.y - p1.y, p2.x - p1.x);
        NetPolygon grain = pc.polygons[index];
        vector<Point> copy;

	//Find the vertices of the grain and rotate them so that the large-scale
	//bond from which the grain was made is oriented horizontally
        for(int iter = 0; iter < grain.vertices.size(); iter++){
            copy.push_back(pc.points[grain.vertices[iter]]);
        }

        pivot = copy[0];
        rotate_points(copy, -angle, pivot);
        get_extremes(copy, left, low, right, high); 
        length = right - left;
        midpoint = Point((left + right)/2, (high+low)/2);
        rotate_point(midpoint, angle, pivot);
        bounds[2] = length + h_offset;
	NetworkComplex grain_nc;
        makeedges(rules, nns, bounds, grain_nc);
        get_extremes(grain_nc.points, left, low, right, high);
        displace(grain_nc.points, midpoint.x - (left+right)/2, midpoint.y - (high+low)/2);
        rotate_points(grain_nc.points, angle, midpoint);
        
        for(Edge e : grain_nc.edges){

            p1 = grain_nc.points[e.idx1];
            p2 = grain_nc.points[e.idx2];
            contains_1 = point_in_polygon(pc.points, pc.polygons[index], p1);
            contains_2 = point_in_polygon(pc.points, pc.polygons[index], p2);

            if(contains_1){
                if(pmap.find(p1) == pmap.end()){
                    pmap.insert(make_pair(p1, pmap.size()));
		    points.push_back(p1);
                }
            }

            if(contains_2){
                if(pmap.find(p2) == pmap.end()){
                    pmap.insert(make_pair(p2, pmap.size()));
		    points.push_back(p2);
                }
            }

            //Identify edges entirely within a grain
            if(contains_1 && contains_2){
                if(to_keep == 1) retain.push_back(Edge(pmap[p1], pmap[p2]));
		else grain_lists[index].push_back(Edge(pmap[p1], pmap[p2]));
            }

            //Identify edges with one point in the current grain and the other
            //point in a neighboring grain
            else if(contains_1 || contains_2){
                vector<int> hits; 
		if(contains_1){
                    hits = hnbvh->within(p2);
		    pt_index = pmap[p1];
                }
		else{
                    hits = hnbvh->within(p1);
		    pt_index = pmap[p2];
                }
		for(int next_hit : hits){
                    IdxPair pair = IdxPair(index, next_hit);
                    if(border_map.find(pair) != border_map.end()){
                        if(index < next_hit){
                            border_map[pair].list1.insert(points[pt_index]);
                        }
			else{
                            border_map[pair].list2.insert(points[pt_index]);
                        }
                    }
                }
            }
        }
    }

    //Create edges joining adjacent grains
    for(auto iter = border_map.begin(); iter != border_map.end(); iter++){
        if(iter->second.list1.size() == 0 || iter->second.list2.size() == 0){
            continue;
        }

	vector<Edge> stitch_edges;
	vector<Point> list1, list2;
	list1.assign(border_map[iter->first].list1.begin(), border_map[iter->first].list1.end());
	list2.assign(border_map[iter->first].list2.begin(), border_map[iter->first].list2.end());
	shuffle(list1.begin(), list1.end(), gen);
        shuffle(list2.begin(), list2.end(), gen);

	//Propose edges connecting points on opposite sides of the grain
	//boundary. Make sure they neither exceed the maximum length nor
	//intersect already created edges.
        for(Point p1 : list1){
            for(Point p2 : list2){
		if(distance_sq(p1, p2) > cutoff_sq) continue;
                no_cross = true;
		for(Edge next_edge : stitch_edges){
                    p3 = points[next_edge.idx1];
                    p4 = points[next_edge.idx2];
		    if(intersection(p1, p2, p3, p4)){
                        no_cross = false;
			break;
                    }
                }
		if(no_cross){
		    stitch_edges.push_back(Edge(pmap[p1], pmap[p2]));
                }
            }
        }

	if(to_keep < 1){
            int key = adj_map_index(iter->first.idx1, iter->first.idx2, pc.polygons.size());
            adj_map.insert(make_pair(key, stitch_edges));
        }
	else{
            retain.insert(retain.end(), stitch_edges.begin(), stitch_edges.end());
	}
    }

    //If the network is to be diluted, use the procedure for creating a random,
    //diluted network that preserves large-scale connectivity. Otherwise,
    //create a network complex with all edges.
    if(to_keep < 1){
        set<int> in_points;
	for(int iter = 0; iter < points.size(); iter++) in_points.insert(iter);
        retain = stitch_network(grain_lists, adj_map, in_points, to_keep);
        final_network = reduced_network(points, retain);
    }

    else final_network = NetworkComplex(points, retain);
}

//Group points according to the large-scale bonds in which they lie. Also find
//all edges within a given large-scale bond, and produce lists of edges
//on the small scale that span large-scale bonds.
void sort_edges(double length, NetworkComplex bottom, PolygonComplex pc, vector<vector<Edge>> &edge_collection){

    int prev;
    vector<int> curr_list;

    for(int i = 0; i < pc.polygons.size(); i ++){
        edge_collection.push_back(vector<Edge>());
    }

    //Data structure for efficient point-in-polygon tests
    shared_ptr<HNetBVH> hnb = make_shared<HNetBVH>(pc.points, pc.polygons);

    for(Edge next_edge : bottom.edges){
        curr_list.clear();
        vector<int> hits1 = hnb->within(bottom.points[next_edge.idx1]);
        vector<int> hits2 = hnb->within(bottom.points[next_edge.idx2]);
	curr_list.insert(curr_list.end(), hits1.begin(), hits1.end());
	curr_list.insert(curr_list.end(), hits2.begin(), hits2.end());
	prev = -1;
	for(int iter = 0; iter < curr_list.size(); iter ++){
            if(curr_list[iter] != prev){
                edge_collection[curr_list[iter]].push_back(next_edge);
            }
	    prev = curr_list[iter];
        }
    }
}

double calc_alignment(vector<Point> small_points, vector<Edge> small_edges, vector<Point> large_points, Edge skel_edge){
    double alignment_sum = 0, mag, dx, dy, skel_nx, skel_ny;

    dx = large_points[skel_edge.idx2].x - large_points[skel_edge.idx1].x;
    dy = large_points[skel_edge.idx2].y - large_points[skel_edge.idx1].y;
    mag = sqrt(dx*dx + dy*dy);
    skel_nx = dx / mag;
    skel_ny = dy / mag;

    for(Edge next_edge : small_edges){
        dx = small_points[next_edge.idx2].x - small_points[next_edge.idx1].x;
        dy = small_points[next_edge.idx2].y - small_points[next_edge.idx1].y;
        mag = sqrt(dx*dx + dy*dy);
        alignment_sum += abs((dx * skel_nx + dy * skel_ny) / mag);
    }

    return alignment_sum / small_edges.size();
}

Point find_match(Point p1, Point p2, double y){
    if(p1.y == y) return p1;
    else return p2;
}

void prepare_grips_displaced(NetworkComplex &nc, bool connect_top_bottom){
    double low, high, miny, maxy;
    vector<double> ycuts;
    bool valid, curr_defined;
    set<Point> bps, tps;
    Point p1, p2, new_point, curr;
    vector<Edge> replace;
    double intersect_x;
    Edge next;
    map<Point, int> pmap;

    do{
        ycuts = getdoubles("Enter the heights for the lower and upper grips: ");

        if(ycuts.size() < 2){
            cerr << "Enter two numbers.\n";
            continue;
        }

        low = ycuts[0];
        high = ycuts[1];

        if(high < low){
            cerr << "The upper bound must not be less than the lower bound.\n";
        }

        else valid = true;
    }while(! valid);

    for(Point p : nc.points) pmap.insert(make_pair(p, pmap.size()));

    while(! nc.edges.empty()){
        next = nc.edges.front();
        p1 = nc.points[next.idx1];
        p2 = nc.points[next.idx2];
        nc.edges.erase(nc.edges.begin());
	miny = p1.y < p2.y ? p1.y : p2.y;
	maxy = p1.y > p2.y ? p1.y : p2.y;
        if(miny >= low && maxy <= high){
            replace.push_back(next);
            if(miny == low) bps.insert(find_match(p1, p2, low));
            if(maxy == high) tps.insert(find_match(p1, p2, high));
        }

        else if(miny < low && maxy >= low){
            intersect_x = x_intersect(nc.points, next, low);
            new_point = Point(intersect_x, low);
	    pmap.insert(make_pair(new_point, pmap.size()));
	    nc.points.push_back(new_point);
            replace.push_back(Edge(pmap[new_point], pmap[find_match(p1, p2, maxy)]));
            bps.insert(new_point);
        }
        
        else if(miny <= high && maxy > high){
            intersect_x = x_intersect(nc.points, next, high);
            new_point = Point(intersect_x, high);
	    pmap.insert(make_pair(new_point, pmap.size()));
	    nc.points.push_back(new_point);
            replace.push_back(Edge(pmap[find_match(p1, p2, miny)], pmap[new_point]));
            tps.insert(new_point);
        }
    }

    if(connect_top_bottom){
        curr_defined = false;
        for(Point next : bps){
            if(! curr_defined){
                curr = next;
                curr_defined = true;
            }
            else{
                replace.push_back(Edge(pmap[curr], pmap[next]));
                curr = next;
            }
        }

        curr_defined = false;
        for(Point next : tps){
            if(! curr_defined){
                curr = next;
                curr_defined = true;
            }
            else{
                replace.push_back(Edge(pmap[curr], pmap[next]));
                curr = next;
            }
        }
    }

    nc = reduced_network(nc.points, replace);
}

//Struct to sort edges by the lowest point
struct bottom_sort {

    vector<Point> points;

    bottom_sort(vector<Point> my_points){
        points = my_points;
    }

    bool operator() (Edge e1, Edge e2){
        double low1 = min(points[e1.idx1].y, points[e1.idx2].y);
        double low2 = min(points[e2.idx1].y, points[e2.idx2].y);
	return low1 < low2;
    }
};

//Struct to sort edges by the highest point
struct top_sort {

    vector<Point> points;

    top_sort(vector<Point> my_points){
        points = my_points;
    }

    bool operator() (Edge e1, Edge e2){
        double top1 = max(points[e1.idx1].y, points[e1.idx2].y);
        double top2 = max(points[e2.idx1].y, points[e2.idx2].y);
	return top1 < top2;
    }
};

void prepare_grips_simple(NetworkComplex &nc){
    vector<Edge> edges;
    set<Edge> low_edges, high_edges;
    set<Point> low_points, high_points;
    double miny, maxy, y1, y2;
    Edge proposed;
    map<Point, int> pmap;

    if(nc.edges.size() == 0) return;
    edges = nc.edges;

    sort(edges.begin(), edges.end(), bottom_sort(nc.points));
    miny = min(nc.points[edges[0].idx1].y, nc.points[edges[0].idx2].y);
    for(Edge next_edge : edges){
        y1 = nc.points[next_edge.idx1].y;
        y2 = nc.points[next_edge.idx2].y;
        if(min(y1, y2) > miny) break;
        low_edges.insert(next_edge);
        if(y1 == miny) low_points.insert(nc.points[next_edge.idx1]); 
        if(y2 == miny) low_points.insert(nc.points[next_edge.idx2]); 
    }

    sort(edges.begin(), edges.end(), top_sort(nc.points));
    maxy = max(nc.points[(*edges.rbegin()).idx1].y, nc.points[(*edges.rbegin()).idx2].y);
    for(auto eiter = edges.rbegin(); eiter != edges.rend(); eiter ++){
        y1 = nc.points[(*eiter).idx1].y;
        y2 = nc.points[(*eiter).idx2].y;
        if(min(y1, y2) < maxy) break;
        high_edges.insert(*eiter);
        if(y1 == maxy) high_points.insert(nc.points[(*eiter).idx1]); 
        if(y2 == maxy) high_points.insert(nc.points[(*eiter).idx2]); 
    }

    for(Point p : nc.points) pmap.insert(make_pair(p, pmap.size()));

    for(auto piter = low_points.begin(); piter != prev(low_points.end()); piter++){
        proposed = Edge(pmap[*piter], pmap[*(next(piter))]);
        if(low_edges.find(proposed) == low_edges.end()){
            nc.edges.push_back(proposed);
        }
    }

    for(auto piter = high_points.begin();piter != prev(high_points.end()); piter++){
        proposed = Edge(pmap[*piter], pmap[*(next(piter))]);
        if(high_edges.find(proposed) == high_edges.end()){
            nc.edges.push_back(proposed);
        }
    }
}

void prepare_grips(NetworkComplex &nc, bool displaced, bool connect){
    if(! displaced) prepare_grips_simple(nc);
    else prepare_grips_displaced(nc, connect);
}

void adjust_bounds(vector<double> &bounds, vector<vector<double>> rules, vector<Point> points){
    double hoffset, voffset, minx, miny, maxx, maxy;

    get_extremes(points, minx, miny, maxx, maxy);
    hoffset = get_h_offset(rules);
    voffset = get_v_offset(rules);

    if(minx < bounds[0]){
        bounds[0] -= ceil((bounds[0] - minx) / hoffset) * hoffset;
    }
    if(miny < bounds[1]){
        bounds[1] -= ceil((bounds[1] - miny) / voffset) * voffset;
    }
    if(maxx > bounds[2]){
        bounds[2] += ceil((maxx - bounds[2]) / hoffset) * hoffset;
    }
    if(maxy > bounds[3]){
        bounds[3] += ceil((maxy - bounds[3]) / voffset) * voffset;
    }
}

NetworkComplex edge_hierarchy(vector<double> bounds, int polyflag, bool getAlign, bool connect, bool verbose){
    stack<NetData> net_stack;
    vector<vector<Edge>> edge_collection;
    NetworkComplex top, bottom, backup;
    PolygonComplex pc;
    bool displacement = false;
    int iter;
    double length, alignment, cutoff;
    FILE *align_report, *poly_report;
    string response;

    loadNetStack(net_stack);
    NetData tdat = net_stack.top();
    length = get_min_dist(tdat.nns);
    net_stack.pop();
    makeedges(tdat.rules, tdat.nns, bounds, top);

    if(tdat.to_keep < 1){
        if(connect) top = random_connected(top, tdat.to_keep);
        else true_random(top, tdat.to_keep);
    }

    if(tdat.displace){
        displace_points_grn(top, tdat.sdev);
        displacement = true;
    }

    if(yesno("Report skeleton?")){
        print_network(top, "Enter a file name or enter to decline: ");
    }

    if(tdat.width > 0){
        if(tdat.displace || getAlign){
            backup.points.insert(backup.points.begin(), top.points.begin(), top.points.end());
            backup.edges.insert(backup.edges.begin(), top.edges.begin(), top.edges.end());
        }
        add_thickness(top, tdat.width, pc, net_stack.size() > 0 || polyflag, !tdat.displace);
    }

    if(net_stack.empty() && yesno("Prepare for grips?")){
        prepare_grips(top, displacement, connect);
    }

    if(polyflag == 1){

	do{
            cout << "Enter the name of the polygon file: ";
            getline(cin, response);
	}while(response.compare("") == 0);

	poly_report = fopen(response.c_str(), "w");

	if(poly_report != NULL){
            if(verbose) print_polygons_verbose(poly_report, pc);
	    else print_polygons_concise(poly_report, pc);
	}
	else cerr << "The file could not be opened.\n";
    }

    if(yesno("Report top network?")){
        print_network(top, "Enter a file name or enter to decline: ");
    }
    
    while(!net_stack.empty()){
        NetData bdat = net_stack.top();
        net_stack.pop();

        adjust_bounds(bounds, bdat.rules, top.points);
        if(displacement && yesno("Use grain-based approach?")){

            while(true){
                cout << "Enter the cutoff length for stitch edges: ";
		getline(cin, response);
		if(sscanf(response.c_str(), "%lf", &cutoff) == 0){
                    cerr << "Enter a number.\n";
                }
		else if(cutoff < 0) cerr << "Enter a positive number.\n";

		else break;
            }

            make_edges_deformed(bdat.rules, bdat.nns, backup, pc, tdat.width, bdat.to_keep, cutoff*cutoff, bottom);
        }
        else{
            makeedges(bdat.rules, bdat.nns, bounds, bottom);
            //cout << "Check 1\n";
            bottom = sieve_edges(top, length, bottom, pc, bdat.to_keep);
        }

        //cout << "Check2\n";

        if(bdat.displace){
            displace_points_grn(bottom, bdat.sdev);
            displacement = true;
        }

        //cout << "Check3\n";
        if(net_stack.empty()){
            if(yesno("Prepare for grips?")){
		double xlow, ylow, xhigh, yhigh;
		get_extremes(bottom.points, xlow, ylow, xhigh, yhigh);
                prepare_grips(bottom, displacement, connect);
		get_extremes(bottom.points, xlow, ylow, xhigh, yhigh);
            }

            if(getAlign){

                cout << "Enter the name for the alignment report: ";
                getline(cin, response);
                align_report = fopen(response.c_str(), "w");

                if(align_report != NULL){
                    sort_edges(length, bottom, pc, edge_collection);

                    for(iter = 0; iter < edge_collection.size(); iter++){
                        alignment  = calc_alignment(bottom.points, edge_collection[iter], backup.points, backup.edges[iter]);
                        fprintf(align_report, "%lf\t%ld\n", alignment, edge_collection[iter].size());
                    }
                    fclose(align_report);
                }
            }
        }

        /*map<Edge, int> edge_tallies;
	int duplicate_count = 0;

	for(Edge e : bottom.edges){
            if(edge_tallies.find(e) == edge_tallies.end()){
                edge_tallies.insert(make_pair(e, 0));
            }
	    edge_tallies[e] ++;
        }

	for(auto et_iter = edge_tallies.begin(); et_iter != edge_tallies.end(); et_iter++){
            if(et_iter->second > 1){
                duplicate_count ++;
                fprintf(stderr, "%lf\t%lf\n", bottom.points[et_iter->first.idx1].x, bottom.points[et_iter->first.idx1].y);
                fprintf(stderr, "%lf\t%lf\n\n", bottom.points[et_iter->first.idx2].x, bottom.points[et_iter->first.idx2].y);
            }
        }

	cout << "Number of edges: " << bottom.edges.size() << endl;
	cout << "Number of duplicates: " << duplicate_count << endl;*/

        if(bdat.width > 0){
            pc.clear();
            if(!net_stack.empty() && (displace || getAlign)){
                backup.clear();
                backup.points.insert(backup.points.begin(), top.points.begin(), top.points.end());
                backup.edges.insert(backup.edges.begin(), top.edges.begin(), top.edges.end());
            }
            add_thickness(bottom, bdat.width, pc, net_stack.size() > 0, !bdat.displace);
        }

	top.assign(bottom);
        length = get_min_dist(bdat.nns);
        //cout << "Check5\n";
    }

    return top;
}

int main(int argc, char **argv){
    vector<vector<double>> rules;
    map<int,vector<vector<double>>> nns;
    vector<double> rule, bounds, nn;
    NetworkComplex nc;
    string nextline;
    int ruleiter, pointiter, count;
    string filename;
    FILE *out;
    double width;
    bool flag, success = false, align = false, connect = true, compact = true;
    int c, polyflag = 0;

    opterr = 0;

    while((c = getopt(argc, argv, "adpv")) != -1){
        switch(c) {
            case 'a':
                align = true;
                break;
            case 'd':
                connect = false;
                break;
            case 'p':
                polyflag = 1;
                break;
            case 'v':
		compact = false;
		break;
            case '?':
                if(isprint(optopt)){
                    fprintf(stderr, "Unknown option: -%c.\n", optopt);
                }
                else{
                    fprintf(stderr, "Unknown option character.\n");
                }
            default:
                break;
         }
    }

    //Read in network bounds
    while(true){
        bounds = getdoubles("Enter bottom left and top right bounds: ");
        if(bounds.size() == 4){
            if((bounds[0] > bounds[2]) || (bounds[1] > bounds[3])){
                fprintf(stderr,"Bounds are improperly ordered.\n");
            }
            else break;
        }
        else{
            fprintf(stderr, "Enter four numbers.\n");
        }
    }

    nc = edge_hierarchy(bounds, polyflag, align, connect, !compact);

    if(! compact){
        print_network(nc, "Enter a file name for output: ");
    }
    else{
	print_network_compact(nc, "Enter a file name for output: ");
    }

    return 0;
}
