/*
This program reads data describing a planar line graph by specifying pairs
of end points for edges, and creates an stl file describing a 3D volume produced
by extruding the 2D network.
*/

#define REAL double
#define VOID int

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <list>
#include <sstream>
#include <algorithm>
#include "triangle.hpp"
#include <stdio.h>
#include <string.h>
#include <cfloat>
#include "network_utils.hpp"

using namespace std;

#define tuple_3d tuple<double, double, double>
#define tuple_2d tuple<double, double>
#define NEG1 -1
#define POS1 1

/*void open_dat_file(string prompt, ifstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}*/

void write_header(FILE *stream, unsigned int num_triangles){
    const char *message = "imaginaryguitarnotesandimaginaryvocalsexistonlyintheimaginationoftheimaginerjoe\n";

    fwrite((const void *)message, sizeof(char), 80, stream);
    fwrite((const void *) &num_triangles, sizeof(int), 1, stream);
}

bool in_poly(vector<Point> plist, vector<int> vertices, Point p){
    bool inpoly = false;
    Point p1, p2;
    double x = p.x, y = p.y;
    int iter, vcount = vertices.size();

    for(iter = 0; iter < vcount; iter++){
        p1 = plist[vertices[iter]];
        p2 = plist[vertices[(iter+1)%vcount]];
        if(p1.y < y && p2.y >= y || p2.y < y && p1.y >= y){
            if(p1.x + (y - p1.y)/(p2.y-p1.y)*(p2.x-p1.x) < x){
                inpoly = !inpoly;
            }
        }
    }

    return inpoly;
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
datfile - an open stream to a file describing the planar line graph
imap - a map from vertices in the planar line graph to integer indices
emap - a map in which each pair is a mapping from a point's integer index to
the indices of all other points to which that point is connected
*/
/*bool read_edges(ifstream &datfile, map<int, Point> &imap, map<int, list<int>> &emap){

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
}*/

/*
*Produce a more compact list of points for the planar line graph that eliminates
groups of three or more colinear points, and generate the corresponding edges.
*/
void get_polygons(NetworkComplex nc, vector<vector<int>> cycle_list, PolygonComplex &pc, int &border_idx){
    int start_index, curr_index, next_index, cycle_index = 0;
    Point prev_point, curr_point, next_point;
    double slope1, slope2, dx1, dy1, dx2, dy2, angle1, angle2;
    double min = FLT_MAX, cycle_min;
    size_t length;

    for(vector<int> next_cycle : cycle_list){
        //list of points describing the next polygon in a set
        vector<int> vertices;
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
                vertices.push_back(pc.points.size());
                pc.points.push_back(curr_point);
            }

            prev_point = curr_point;
            curr_index = next_index;
        }while(curr_index != start_index);

        //Create edges associated with the polygon that has been created.
        for(int iter = 0; iter < vertices.size(); iter ++){
            pc.edges.insert(Edge(vertices[iter],vertices[(iter+1)%vertices.size()]));
        }

        //Determine whether this cycle has the lowest y coordinate of any cycle
        //so far. The border will contain the lowest point in the planar line
        //graph.
        if(cycle_min < min){
            min = cycle_min;
            border_idx = cycle_index;
        }

        pc.polygons.push_back(NetPolygon(pc.points, vertices, cycle_index));
        cycle_index ++;
    }
}

//Don't bother coallescing points. Just translate indices into points, and
//find the border.
void get_polygons_simple(vector<vector<int>> cycle_list, PolygonComplex &pc, int &border_idx){
    double miny = FLT_MAX;
    int c_iter, v_iter;

    for(c_iter = 0; c_iter < cycle_list.size(); c_iter++){
        vector<int> cycle = cycle_list[c_iter];
        for(int v_iter = 0; v_iter < cycle.size(); v_iter ++){
            if(pc.points[cycle[v_iter]].y < miny){
                miny = pc.points[cycle[v_iter]].y;
                border_idx = c_iter;
            }
            pc.edges.insert(Edge(cycle[v_iter], cycle[(v_iter+1)%cycle.size()]));
        }
        pc.polygons.push_back(NetPolygon(pc.points, cycle, c_iter));
    }
}

//Helper function for get_normals
void process_polygon(PolygonComplex pc, int poly_idx, vector<vector<tuple_3d>> &normals, int case1, int case2){

    double sign, midx, midy, dx, dy, inv_mag;
    int iter, pcount;
    Point p1, p2, guess1, guess2;
    bool is1, is2;
    vector<tuple_3d> normlist;
    NetPolygon poly = pc.polygons[poly_idx];

    pcount = poly.vertices.size();
    if(pcount < 3){
        cerr << "That's no polygon!\n";
    }

    p1 = pc.points[poly.vertices[0]];
    p2 = pc.points[poly.vertices[1]];
    dx = p2.x - p1.x;
    dy = p2.y - p1.y;
    midx = (p1.x + p2.x) / 2;
    midy = (p1.y + p2.y) / 2;
    //Find the normal to the first edge pointing into the polygon
    guess1 = Point(midx - dy, midy + dx);
    guess2 = Point(midx + dy, midy - dx);
    do{
        guess1 = Point(midx - dy, midy + dx);
        guess2 = Point(midx + dy, midy - dx);
        is1 = in_poly(pc.points, poly.vertices, guess1);
        is2 = in_poly(pc.points, poly.vertices, guess2);
        if(! (is1 || is2)){
            dx /= 2;
            dy /= 2;
        }
    }while(! (is1 || is2));

    //Set the procedure for finding normals for successive edges
    sign = is1 ? case1 : case2;
    inv_mag = 1 / sqrt(dx*dx + dy*dy);
    normlist.push_back(make_tuple((float)(-dy*sign*inv_mag),(float)(dx*sign*inv_mag),0));

    for(iter = 1; iter != pcount; iter ++){
        p1 = pc.points[poly.vertices[iter]];
        p2 = pc.points[poly.vertices[(iter + 1)%pcount]];
        dx = p2.x - p1.x;
        dy = p2.y - p1.y;
        normlist.push_back(make_tuple((float)(-dy*sign*inv_mag),(float)(dx*sign*inv_mag),0));
    }

    normals.push_back(normlist);
}

void get_normals(PolygonComplex pc, int border_idx, vector<vector<tuple_3d>> &normals){

    double sign, midx, midy, dx, dy, inv_mag;

    for(int iter = 0; iter < pc.polygons.size(); iter ++){
        //Process the border separately. Normals should point out of the border,
        //and no interior point is desired.
        if(iter == border_idx) {
            process_polygon(pc, border_idx, normals, NEG1, POS1);
        }
        else process_polygon(pc, iter, normals, POS1, NEG1);
    }
}

void get_holes(PolygonComplex pc, int border_idx, vector<tuple_2d> &holes){
    double midx, midy, dx, dy;
    Point p1, p2, guess1, guess2;
    bool is1, is2;

    for(int iter = 0; iter < pc.polygons.size(); iter ++){
        if(iter == border_idx) continue;
        NetPolygon next_poly = pc.polygons[iter];
        if(next_poly.vertices.size() < 3){
            cerr << "That's no polygon!\n";
        }

        p1 = pc.points[next_poly.vertices[0]];
        p2 = pc.points[next_poly.vertices[1]];
        dx = p2.x - p1.x;
        dy = p2.y - p1.y;
        midx = (p1.x + p2.x) / 2;
        midy = (p1.y + p2.y) / 2;
        //Find the normal to the first edge pointing into the polygon
        do{
            guess1 = Point(midx - dy, midy + dx);
            guess2 = Point(midx + dy, midy - dx);
            is1 = in_poly(pc.points, next_poly.vertices, guess1);
            is2 = in_poly(pc.points, next_poly.vertices, guess2);
            if(is1){
                holes.push_back(make_tuple(guess1.x, guess1.y));
            }
            else if(is2){
                holes.push_back(make_tuple(guess2.x, guess2.y));
            }
            else{
                dx /= 2;
                dy /= 2;
            }
        }while(! (is1 || is2));
    }
}

void add_segments(int *segmentlist, int pcount, int &plabel, int &seg_index){
    int iter, start_plabel = plabel;

    for(iter = 0; iter < pcount - 1; iter++){
        segmentlist[seg_index] = plabel;
        seg_index++;
        segmentlist[seg_index] = plabel + 1;
        seg_index++;
        plabel ++;
    }
    segmentlist[seg_index] = plabel;
    seg_index++;
    segmentlist[seg_index] = start_plabel;
    seg_index++;
    plabel++;
}

void create_input(triangulateio *in, PolygonComplex pc, int border_idx, vector<tuple_2d> holes){
    int pindex = 0, plabel = 0, seg_index = 0, hole_index = 0, point_count;

    point_count = pc.points.size();
    in->numberofpoints = point_count;
    in->numberofsegments = point_count;
    in->pointlist = (double *) malloc(point_count * 2 * sizeof(double));
    in->segmentlist = (int *) malloc(point_count * 2 * sizeof(int));
    in->numberofregions = 0;

    for(int iter = 0; iter < pc.polygons.size(); iter ++){
        if(iter == border_idx) continue;
        NetPolygon next_poly = pc.polygons[iter];
        for(int next_vertex : next_poly.vertices){
            in->pointlist[pindex++] = pc.points[next_vertex].x;
            in->pointlist[pindex++] = pc.points[next_vertex].y;
        }
        add_segments(in->segmentlist, next_poly.vertices.size(), plabel, seg_index);
    }

    NetPolygon border = pc.polygons[border_idx];
    for(int next_vertex : border.vertices){
        in->pointlist[pindex++] = pc.points[next_vertex].x;
        in->pointlist[pindex++] = pc.points[next_vertex].y;
    }
    add_segments(in->segmentlist, border.vertices.size(), plabel, seg_index);

    in->segmentmarkerlist = (int *) malloc(sizeof(int) * in->numberofsegments);
    in->numberofpointattributes = 0;
    in->pointmarkerlist = (int *) NULL;
    in->numberofedges = 0;
    in->edgelist = (int *) NULL;
    in->numberofholes = holes.size();
    in->holelist = (double *) malloc(in->numberofholes * 2 * sizeof(double));
    for(tuple_2d nexthole : holes){
        in->holelist[hole_index++] = get<0>(nexthole);
        in->holelist[hole_index++] = get<1>(nexthole);
    }
}

void setup_output(triangulateio *out){
    out->pointlist = (double *) NULL;
    out->numberofpointattributes = 0;
    out->numberoftriangleattributes = 0;
    out->trianglelist = (int *) NULL;
    out->edgelist = (int *) NULL;
    out->edgemarkerlist = (int *) NULL;
    out->segmentlist = (int *) NULL;
    out->pointmarkerlist = (int *) NULL;
    out->segmentmarkerlist = (int *) NULL;
}

void write_point(FILE *stlfile, tuple_3d point){
    float px = (float) get<0>(point), py = (float) get<1>(point), pz = (float) get<2>(point);
    fwrite((const void *) &px, sizeof(float), 1, stlfile);
    fwrite((const void *) &py, sizeof(float), 1, stlfile);
    fwrite((const void *) &pz, sizeof(float), 1, stlfile);
}

void to_farr(tuple_3d normal, float *farr){
    farr[0] = (float) get<0>(normal);
    farr[1] = (float) get<1>(normal);
    farr[2] = (float) get<2>(normal);
}

/*Ensure the order of the vertices of a triangle is consistent with its normal
*vector. This is accomplished by first taking the cross product of two vectors:
*the difference between point 2 and point 1, and the difference betwen point
*3 and point 2. This cross product is dotted with the normal. If the result is
*positive, the order of the points is consistent with the normal vector
*/
double check_orientation(tuple_3d normal, tuple_3d p1, tuple_3d p2, tuple_3d p3){
    double cross_x, cross_y, cross_z;

    cross_x = (get<1>(p2)-get<1>(p1))*(get<2>(p3)-get<2>(p2))*(get<2>(p2)-get<2>(p1))*(get<1>(p3)-get<1>(p2));
    cross_y = (get<2>(p2)-get<2>(p1))*(get<0>(p3)-get<0>(p2))*(get<0>(p2)-get<0>(p1))*(get<2>(p3)-get<2>(p2));
    cross_z = (get<0>(p2)-get<0>(p1))*(get<1>(p3)-get<1>(p2))*(get<1>(p2)-get<1>(p1))*(get<0>(p3)-get<0>(p2));
    return cross_x*get<0>(normal)+cross_y*get<1>(normal)+cross_z*get<2>(normal);
}

void write_side(FILE *stlfile, PolygonComplex pc, int index, vector<tuple_3d> normal_set, float height){

    int vertex_iter, poly_size;
    short unsigned int nada = 0;
    tuple_3d v1, v2, v3;
    double ocheck;
    float farr[3];
    Point p1, p2;
    NetPolygon polygon = pc.polygons[index];
    poly_size = polygon.vertices.size();

    for(vertex_iter = 0; vertex_iter < poly_size; vertex_iter ++){
        p1 = pc.points[polygon.vertices[vertex_iter]];
        p2 = pc.points[polygon.vertices[(vertex_iter + 1) % poly_size]];

        //Bottom triangle
        to_farr(normal_set[vertex_iter], farr);
        fwrite((const void*) farr,sizeof(float),3,stlfile);
        v1 = make_tuple(p1.x, p1.y, 0);
        v2  = make_tuple(p2.x, p2.y, 0);
        v3  = make_tuple(p2.x, p2.y, height);

        //Check orientation
        ocheck = check_orientation(normal_set[vertex_iter], v1, v2, v3);

        if(ocheck >0){
            write_point(stlfile, v1);
            write_point(stlfile, v2);
        }
        else{
            write_point(stlfile, v2);
            write_point(stlfile, v1);
        }
        write_point(stlfile, v3);
        fwrite((const void *)&nada, sizeof(short unsigned int), 1, stlfile);

        //Top triangle
        fwrite((const void*) farr,sizeof(float),3,stlfile);
        v1 = make_tuple(p2.x, p2.y, height);
        v2 = make_tuple(p1.x, p1.y, height);
        v3 = make_tuple(p1.x, p1.y, 0);

        ocheck = check_orientation(normal_set[vertex_iter], v1, v2, v3);

        if(ocheck >0){
            write_point(stlfile, v1);
            write_point(stlfile, v2);
        }
        else{
            write_point(stlfile, v2);
            write_point(stlfile, v1);
        }
        write_point(stlfile, v3);
        fwrite((const void *)&nada, sizeof(short unsigned int), 1, stlfile);
    }
}

void write_as_float(FILE *stream, double *arr){
    float fvals[2];
    fvals[0] = arr[0];
    fvals[1] = arr[1];
    fwrite((const void *) fvals, sizeof(float), 2, stream);
}

//void write_stl(FILE *stlfile, triangulateio out, vector<vector<Point>> polygons, vector<Point> border, vector<vector<tuple_3d>> normals,  float height){
void write_stl(FILE *stlfile, triangulateio out, PolygonComplex pc, int border_idx, vector<vector<tuple_3d>> normals,  float height){
    int tri_iter, vertex_iter, v_index, poly_iter;
    unsigned short int nada = 0;
    tuple_3d bnorm = make_tuple(0,0,-1), tnorm = make_tuple(0,0,1);
    float barr[3] = {0,0,-1}, tarr[3] = {0,0,1};
    Point p1, p2;
    vector<tuple_3d> normal_set;
    double ocheck;
    tuple_3d verts[3];
    float farr[3];

    //report top and bottom triangulations
    for(tri_iter = 0; tri_iter < out.numberoftriangles; tri_iter++){
        //Report the bottom
        fwrite((const void *) barr, sizeof(float), 3, stlfile);
        for(vertex_iter = 0; vertex_iter < 3; vertex_iter++){
            v_index = out.trianglelist[tri_iter * 3 + vertex_iter];
            verts[vertex_iter] = make_tuple(out.pointlist[2*v_index],out.pointlist[2*v_index+1],0);
        }
        ocheck = check_orientation(bnorm, verts[0], verts[1], verts[2]);
        if(ocheck >0){
            write_point(stlfile, verts[0]);
            write_point(stlfile, verts[1]);
        }
        else{
            write_point(stlfile, verts[1]);
            write_point(stlfile, verts[0]);
        }
        write_point(stlfile, verts[2]);
        fwrite((const void *) &nada, sizeof(short unsigned int), 1, stlfile);

        //Now write the top
        fwrite((const void *) tarr, sizeof(float), 3, stlfile);
        for(vertex_iter = 0; vertex_iter < 3; vertex_iter++){
            v_index = out.trianglelist[tri_iter * 3 + vertex_iter];
            verts[vertex_iter] = make_tuple(out.pointlist[2*v_index],out.pointlist[2*v_index+1], (double) height);
        }
        ocheck = check_orientation(tnorm, verts[0], verts[1], verts[2]);
        if(ocheck >0){
            write_point(stlfile, verts[0]);
            write_point(stlfile, verts[1]);
        }
        else{
            write_point(stlfile, verts[1]);
            write_point(stlfile, verts[0]);
        }
        write_point(stlfile, verts[2]);
        fwrite((const void *) &nada, sizeof(short unsigned int), 1, stlfile);
    }

    //Report triangulation for the sides

    //Start with interior polygons
    for(poly_iter = 0; poly_iter < pc.polygons.size(); poly_iter ++){
        if(poly_iter == border_idx) continue;
        normal_set = normals[poly_iter];
        write_side(stlfile, pc, poly_iter, normal_set, height);
    }

    //Finally, make the border
    write_side(stlfile, pc, border_idx, *normals.rbegin(), height);

}

void report_triangles(FILE *stream, triangulateio out){
    int tri_iter, index1, index2, index3;
    double *p1, *p2, *p3;

    for(tri_iter = 0; tri_iter < out.numberoftriangles; tri_iter++){
        index1 = out.trianglelist[tri_iter * 3];
        index2 = out.trianglelist[tri_iter * 3 + 1];
        index3 = out.trianglelist[tri_iter * 3 + 2];
        p1 = out.pointlist + 2*index1;
        p2 = out.pointlist + 2*index2;
        p3 = out.pointlist + 2*index3;
        fprintf(stream, "%f\t%f\n%f\t%f\n\n", p1[0], p1[1], p2[0], p2[1]);
        fprintf(stream, "%f\t%f\n%f\t%f\n\n", p2[0], p2[1], p3[0], p3[1]);
        fprintf(stream, "%f\t%f\n%f\t%f\n\n", p3[0], p3[1], p1[0], p1[1]);
    }
}

void report_segments(FILE *stream, triangulateio in){
    int seg_iter, p1, p2;

    for(seg_iter = 0; seg_iter < in.numberofsegments; seg_iter++){
        p1 = in.segmentlist[seg_iter * 2];
        p2 = in.segmentlist[seg_iter * 2 + 1];
        fprintf(stream, "%lf\t%lf\n%lf\t%lf\n\n", in.pointlist[2*p1], in.pointlist[2*p1+1], in.pointlist[2*p2], in.pointlist[2*p2+1]);
    }
}

void report_holes(FILE *stream, triangulateio in){
    int hole_iter;

    for(hole_iter = 0; hole_iter < in.numberofholes; hole_iter ++){
        fprintf(stream, "%f", in.holelist[hole_iter * 2]);
        fprintf(stream, "\t%f\n", in.holelist[hole_iter * 2 + 1]);
    }
}

void read_boundary_segments(triangulateio *out, vector<Point> &points, map<int, list<int>> &emap){
    int iter, pindex = 0, p1, p2;
    double x,y;

    map<int, int> imap;

    for(iter = 0; iter < out->numberofsegments; iter++){
        p1 = out->segmentlist[2*iter];
        p2 = out->segmentlist[2*iter + 1];

        if(imap.find(p1) == imap.end()){
            emap.insert(make_pair(pindex, list<int>()));
            imap.insert(make_pair(p1, pindex));
            x = out->pointlist[2*p1];
            y = out->pointlist[2*p1 + 1];
            points.push_back(Point(x,y));
            pindex ++;
        }
        if(imap.find(p2) == imap.end()){
            emap.insert(make_pair(pindex, list<int>()));
            imap.insert(make_pair(p2, pindex));
            x = out->pointlist[2*p2];
            y = out->pointlist[2*p2 + 1];
            points.push_back(Point(x,y));
            pindex ++;
        }
        emap[imap[p1]].push_back(imap[p2]);
        emap[imap[p2]].push_back(imap[p1]);
    }
}

/*
* Given the name of a data file and an output stl file, create
* an open input stream, obtain a map of indices to lists of other indices.
* Use the map to build a table of cycles, then  convert the table of cycles to
* polygon data. Extrude the polygonal tiling, and report the resulting surface
* to a binary .stl file.
*/
void export_extruded_svg(NetworkComplex nc, float height, string file_name){
    ifstream datfile;
    FILE *stl_file;
    //FILE *trifile, *segment_file, *hole_file;
    map<int, list<int>> emap;
    vector<vector<int>> cycle_list;
    //Specify one point inside each interior polygon for the triangulation
    //process
    vector<tuple_2d> holes;
    //Normals for the triangles perpendicular to the horizontal plane
    vector<vector<tuple_3d>> normals;
    int numread, border_idx, iter, p1, p2;
    unsigned int trianglecount;
    string response;
    char options[7] = "pqQze\0";
    PolygonComplex pc;

    //Structures to hold input and output data for triangulation of the network
    struct triangulateio in, out;

    //Use maps to build cycles described by integer indices in cyclical order
    build_plc_map(nc, emap);
    cycle_list = build_cycles(emap);
    //cout << "Cycles built.\n";

    //Resolve the planar line graph describing the network into a set of 
    //polygons
    //get_polygons(nc, cycle_list, polygons, border, point_count);
    get_polygons(nc, cycle_list, pc, border_idx);
    //cout << "Polygons made.\n";

    //Find holes for triangulation and normal vectors for stl file creation
    get_holes(pc, border_idx, holes);
    //cout << "Holes found.\n";

    //Set up input to Triangle's triangulate function and perform triangulation
    create_input(&in, pc, border_idx, holes);
    setup_output(&out);
    /*if(yesno("Write segments?")){
        do{
            cout << "Enter the file name: ";
            getline(cin, response);
            segment_file = fopen(response.c_str(), "w");
        }while(segment_file == NULL);

        report_segments(segment_file, in);
        fclose(segment_file);
    }
    if(yesno("Write holes?")){
        do{
            cout << "Enter the file name: ";
            getline(cin, response);
            hole_file = fopen(response.c_str(), "w");
        }while(hole_file == NULL);

        report_holes(hole_file, in);
        fclose(hole_file);
    }*/

    triangulate(options, &in, &out, (struct triangulateio *) NULL);

    //Now that new points have been added by Triangle, find the new boundary
    //segments, create a new map from boundary points to other boundary points,
    //and create new ordered cycles of boundary points.
    emap.clear();
    cycle_list.clear();
    pc.clear();
    read_boundary_segments(&out, pc.points, emap);
    cycle_list = build_cycles(emap);
    get_polygons_simple(cycle_list, pc, border_idx);
    get_normals(pc, border_idx, normals);

    //Find the total number of triangles, and report triangles to an stl file
    trianglecount = 2*out.numberofsegments + 2*out.numberoftriangles;

    //Write STL file
    stl_file = fopen(file_name.c_str(), "wb");
    if(stl_file == NULL){
        cerr << "The desired file could not be opened.\n";
        return;
    }
    write_header(stl_file, trianglecount);
    write_stl(stl_file, out, pc, border_idx, normals,  height);
    fclose(stl_file);

    /*if(yesno("Report triangles?")){
        do{
            cout << "Enter the file name: ";
            getline(cin, response);
            trifile = fopen(response.c_str(), "w");
        }while(trifile == NULL);

        report_triangles(trifile, out);
        fclose(trifile);
    }*/

    //Clean up triangulation data structures
    free(in.pointlist);
    free(in.segmentlist);
    free(in.segmentmarkerlist);
    free(out.pointlist);
    free(out.edgelist);
    free(out.edgemarkerlist);
    free(out.segmentlist);
    free(out.segmentmarkerlist);
    free(out.trianglelist);
}
