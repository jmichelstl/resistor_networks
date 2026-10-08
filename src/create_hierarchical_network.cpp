/*
 * Author: Jonathan MichenyThis program prompts a user for boundaries and rules
 * for producing lattice sites and nearest neighbor connections. From here, a
 * set of line segments joining nearest neighbor locations is produced.
 * Connections can be removed from the network in two ways; one method is
 * completely random, while the other uses randomness but also guarantees that
 * the network remains fully connected.
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
#include <tuple>
#include <memory>
#include "hnet_job.hpp"

using namespace std;

#define MAX_ATTEMPTS 10

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

//Struct to represent the data for constructing a lattice
struct LatticeRecipe {

    LatticeRecipe(){}

    LatticeRecipe(vector<vector<double>> my_rules,  map<int,vector<vector<double>>> my_nns){
        rules = my_rules;;
        nns = my_nns;
    }

    vector<vector<double>> rules;
    map<int,vector<vector<double>>> nns;
};

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
            cerr << "Out of bounds edge: (" << e.idx1 << ", " << e.idx2 <<")\n";
        }
        Point p1 = nc.points[e.idx1];
        Point p2 = nc.points[e.idx2];
        fprintf(out, "%10.8lf %10.8lf \n%10.8lf %10.8lf \n\n", p1.x, p1.y, p2.x, p2.y);
    }

    fclose(out);
}

void print_network_automated(NetworkComplex nc, string name){    
    FILE *out = NULL;

    out = fopen(name.c_str(), "w");
    if(out == NULL){
        cerr << "The specified output file could not be opened.\n";
        return;
    }

    for(Edge e : nc.edges){
        if(e.idx1 >= nc.points.size() || e.idx2 >= nc.points.size()){
            cerr << "Out of bounds edge: (" << e.idx1 << ", " << e.idx2 <<")\n";
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

void print_network_compact_automated(NetworkComplex nc, string name){ 
    FILE *out = NULL;

    out = fopen(name.c_str(), "w");
    if(out == NULL){
        cerr << "The specified output file could not be opened.\n";
        return;
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

void get_lattice_info(vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns){

    bool success = false;
    double scale;
    string nextline, name;
    int file_name_attempts = 0, import_attempts = 0;

    do{
        scale = getdoubles("Enter the scale factor: ").at(0);
        //Prompt for file name
        do{
            printf("Enter the lattice file name: ");
            getline(cin, nextline);
            name = split(nextline, ' ')[0];

            if(!name.empty()) break;

            if(! yesno("No file name was read. Try again? ")){
                return;
            }
        }while(file_name_attempts++ < MAX_ATTEMPTS);

        success = import_lattice(name, rules, nns, scale);
        if(!success) if(!yesno("Read failed. Try again?")) break;
    }while(!success && import_attempts < MAX_ATTEMPTS);
}

void loadNetStack(stack<NetData>& net_stack){
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

void prepare_grips_displaced_interactive(NetworkComplex &nc, bool connect_top_bottom){
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

void prepare_grips_displaced_automated(NetworkComplex &nc, bool connect_top_bottom, double low, double high){
    double miny, maxy;
    bool curr_defined;
    set<Point> bps, tps;
    Point p1, p2, new_point, curr, in_pt;
    vector<Point> replace_points;
    vector<Edge> replace_edges;
    double intersect_x;
    Edge next;
    map<Point, unsigned int> pmap;

    while(! nc.edges.empty()){
        next = nc.edges.front();
        p1 = nc.points[next.idx1];
        p2 = nc.points[next.idx2];
        nc.edges.erase(nc.edges.begin());
        miny = p1.y < p2.y ? p1.y : p2.y;
        maxy = p1.y > p2.y ? p1.y : p2.y;
        if(miny >= low && maxy <= high){
            if(pmap.find(p1) == pmap.end()){
                replace_points.push_back(p1);
                pmap.insert(make_pair(p1, pmap.size()));
            }
            if(pmap.find(p2) == pmap.end()){
                replace_points.push_back(p2);
                pmap.insert(make_pair(p2, pmap.size()));
            }
            replace_edges.push_back(Edge(pmap[p1], pmap[p2]));
            if(connect_top_bottom){
                if(miny == low) bps.insert(find_match(p1, p2, low));
                if(maxy == high) tps.insert(find_match(p1, p2, high));
            }
        }

        else if(miny < low && maxy >= low){
            intersect_x = x_intersect(nc.points, next, low);
            new_point = Point(intersect_x, low);
            if(pmap.find(new_point) == pmap.end()){
                pmap.insert(make_pair(new_point, pmap.size()));
                replace_points.push_back(new_point);
            }
            in_pt = find_match(p1, p2, maxy);
            if(pmap.find(in_pt) == pmap.end()){
                pmap.insert(make_pair(in_pt, pmap.size()));
                replace_points.push_back(in_pt);
            }
            
            replace_edges.push_back(Edge(pmap[new_point], pmap[in_pt]));
            if(connect_top_bottom) bps.insert(new_point);
        }
        
        else if(miny <= high && maxy > high){
            intersect_x = x_intersect(nc.points, next, high);
            new_point = Point(intersect_x, high);
            if(pmap.find(new_point) == pmap.end()){
                pmap.insert(make_pair(new_point, pmap.size()));
                replace_points.push_back(new_point);
            }
            in_pt = find_match(p1, p2, miny);
            if(pmap.find(in_pt) == pmap.end()){
                pmap.insert(make_pair(in_pt, pmap.size()));
                replace_points.push_back(in_pt);
            }
            replace_edges.push_back(Edge(pmap[in_pt],pmap[new_point]));
            if(connect_top_bottom) tps.insert(new_point);
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
                replace_edges.push_back(Edge(pmap[curr], pmap[next]));
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
                replace_edges.push_back(Edge(pmap[curr], pmap[next]));
                curr = next;
            }
        }
    }

    nc.points.assign(replace_points.begin(), replace_points.end());
    nc.edges.assign(replace_edges.begin(), replace_edges.end());
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
        if(max(y1, y2) < maxy) break;
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

void prepare_grips(NetworkComplex &nc, bool displaced){
    if(! displaced) prepare_grips_simple(nc);
    else {
        prepare_grips_displaced_interactive(nc, yesno("Connect the top and bottom? "));
    }
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

NetworkComplex edge_hierarchy_interactive(vector<double> bounds, int polyflag, bool getAlign, bool connect, bool verbose){
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
            backup.points.assign(top.points.begin(), top.points.end());
            backup.edges.assign(top.edges.begin(), top.edges.end());
        }
        add_thickness(top, tdat.width, pc, net_stack.size() > 0 || polyflag, !tdat.displace);
    }

    if(net_stack.empty() && yesno("Prepare for grips?")){
        prepare_grips(top, displacement);
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
            bottom = sieve_edges(top, length, bottom, pc, bdat.to_keep);
        }

        if(bdat.displace){
            displace_points_grn(bottom, bdat.sdev);
            displacement = true;
        }

        if(net_stack.empty()){
            if(yesno("Prepare for grips?")){
                prepare_grips(bottom, displacement);
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

        if(bdat.width > 0){
            pc.clear();
            if(!net_stack.empty() && (displace || getAlign)){
                backup.clear();
                backup.points.assign(top.points.begin(), top.points.end());
                backup.edges.assign(top.edges.begin(), top.edges.end());
            }
            add_thickness(bottom, bdat.width, pc, net_stack.size() > 0, !bdat.displace);
        }

        top.assign(bottom);
        length = get_min_dist(bdat.nns);
    }

    return top;
}

void create_network_interactive(bool align, bool connect, bool compact, char polyflag){

    vector<double> bounds;
    NetworkComplex nc;

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

    nc = edge_hierarchy_interactive(bounds, polyflag, align, connect, !compact);

    if(! compact){
        print_network(nc, "Enter a file name for output: ");
    }
    else{
        print_network_compact(nc, "Enter a file name for output: ");
    }
}

NetworkComplex edge_hierarchy_automated(HNTask task, vector<LatticeRecipe> recipes, bool verbose){
    vector<vector<Edge>> edge_collection;
    NetworkComplex top, bottom, backup;
    PolygonComplex pc;
    int iter, idx;
    double length, alignment;
    FILE *align_report, *poly_report;
    bool displaced = false;
    vector<double> bounds({task.minx, task.miny, task.maxx, task.maxy});

    HNLayer tlayer = task.layers[0];
    length = get_min_dist(recipes[0].nns);
    makeedges(recipes[0].rules, recipes[0].nns, bounds, top);

    if(tlayer.bond_occupation < 1){
        if(task.connected) top = random_connected(top, tlayer.bond_occupation);
        else true_random(top, tlayer.bond_occupation);
    }

    if(tlayer.displacement > 0){
        displace_points_grn(top, tlayer.displacement);
        displaced = true;
    }

    if(task.skeleton_file.compare("") != 0){
        print_network_automated(top, task.skeleton_file);
    }

    if(tlayer.width > 0){
        if(tlayer.displacement > 0 || task.alignment_file.compare("") != 0){
            backup.points.insert(backup.points.begin(), top.points.begin(), top.points.end());
            backup.edges.insert(backup.edges.begin(), top.edges.begin(), top.edges.end());
        }
        add_thickness(top, tlayer.width, pc, task.layers.size() > 1 || task.poly_file.compare("") != 0, !(tlayer.displacement > 0));
    }

    if(task.layers.size() == 1 && task.add_grips){
        if(task.displaced) {
            prepare_grips_displaced_automated(top, task.connect_top_bottom, task.grip_lower, task.grip_upper);
        }
        else prepare_grips_simple(top);
    }

    if(task.poly_file.compare("") != 0){

        poly_report = fopen(task.poly_file.c_str(), "w");

        if(poly_report != NULL){
            if(verbose) print_polygons_verbose(poly_report, pc);
            else print_polygons_concise(poly_report, pc);
        }
        else cerr << "The polygon file could not be opened.\n";
    }

    if(task.top_file.compare("") != 0){
        print_network_automated(top, task.top_file);
    }

    for(iter = 1; iter < task.layers.size(); iter ++){
        HNLayer blayer = task.layers[iter];

        adjust_bounds(bounds, recipes[iter].rules, top.points);
        if(displaced && blayer.grain_based){
            make_edges_deformed(recipes[iter].rules, recipes[iter].nns, backup, pc, tlayer.width, blayer.bond_occupation, blayer.stitch_cutoff*blayer.stitch_cutoff, bottom);
        }
        else{
            makeedges(recipes[iter].rules, recipes[iter].nns, bounds, bottom);
            bottom = sieve_edges(top, length, bottom,pc,blayer.bond_occupation);
        }

        if(blayer.displacement > 0){
            displace_points_grn(bottom, blayer.displacement);
            displaced = true;
        }

        if(iter == task.layers.size() - 1){
            if(task.add_grips){
                if(task.displaced) {
                    prepare_grips_displaced_automated(bottom, task.connect_top_bottom, task.grip_lower, task.grip_upper);
                }
                else prepare_grips_simple(bottom);
            }

            if(task.alignment_file.compare("") != 0){

                align_report = fopen(task.alignment_file.c_str(), "w");

                if(align_report != NULL){
                    sort_edges(length, bottom, pc, edge_collection);

                    for(idx = 0; idx < edge_collection.size(); idx++){
                        alignment  = calc_alignment(bottom.points, edge_collection[idx], backup.points, backup.edges[idx]);
                        fprintf(align_report, "%lf\t%ld\n", alignment, edge_collection[idx].size());
                    }
                    fclose(align_report);
                }
            }
        }

        if(blayer.width > 0){
            pc.clear();
            if(iter < task.layers.size() - 1 && (task.displaced || task.alignment_file.compare("") != 0)){
                backup.clear();
                backup.points.insert(backup.points.begin(), top.points.begin(), top.points.end());
                backup.edges.insert(backup.edges.begin(), top.edges.begin(), top.edges.end());
            }
            add_thickness(bottom, blayer.width, pc, iter < task.layers.size()-1, blayer.displacement == 0);
        }

        top.assign(bottom);
        length = get_min_dist(recipes[iter].nns);
        tlayer = blayer;
    }

    return top;
}

void create_networks_automated(HNJob job, bool compact){

    vector<LatticeRecipe> recipes;

    //Find instructions for creating points and edges for each layer
    for(int iter = 0; iter < job.lattice_files.size(); iter ++){
        LatticeRecipe r;
        if(! import_lattice(job.lattice_files[iter], r.rules, r.nns, job.scales[iter])){
            cerr << "Importation of lattice instructions failed.\n";
            return;
        }
        recipes.emplace(recipes.begin(), r);
    }

    //Build each realization for each combination of bond portions
    for(HNTask next_task : job.tasks){
        NetworkComplex nc = edge_hierarchy_automated(next_task, recipes, !compact);
        if(! compact){
            print_network_automated(nc, next_task.output_file);
        }
        else{
            print_network_compact_automated(nc, next_task.output_file);
        }
    }
}

int main(int argc, char **argv){
    bool align = false, connect = true, compact = true, interactive = true;
    int c, polyflag = 0;
    string job_file_name;

    opterr = 0;

    while((c = getopt(argc, argv, "adj:pv")) != -1){
        switch(c) {
            case 'a':
                align = true;
                break;
            case 'd':
                connect = false;
                break;
            case 'j':
                interactive = false;
                job_file_name.assign(string(optarg));
                break;
            case 'p':
                polyflag = 1;
                break;
            case 'v':
                compact = false;
                break;
            case '?':
                if(optopt == 'j'){
                    cerr << "Option \"j\" requires a job file name.\n";
                }
                else if(isprint(optopt)){
                    fprintf(stderr, "Unknown option: -%c.\n", optopt);
                }
                else{
                    fprintf(stderr, "Unknown option character.\n");
                }
            default:
                break;
         }
    }

    if(interactive) {
        create_network_interactive(align, connect, compact, polyflag);
    }

    else {
        HNJob job;
        if(! parse_hnet_job(job_file_name, job)){
            cerr << "The specified job file could not be parsed.\n";
            return -1;
        }
        create_networks_automated(job, compact);
    }

    return 0;
}
