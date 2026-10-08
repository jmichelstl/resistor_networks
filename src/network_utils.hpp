#ifndef NET_UTILS_H
#define NET_UTILS_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cfloat>
#include <unordered_map>
#include <cstdint>
#include <random>

using namespace std;

#define PRIME1 (int) 1021297
#define PRIME2 (int) 1021651 //Primes for computing hash functions
#define FLOAT_TOL 1E-10 //Tolerance for deeming floating point numbers different
#define BIG_SLOPE 1E10
#define INF 1E12

//Function prototypes for Robust Geometric Predicates functionality, implemented
//by Jonathan Richard Shewchuck.
extern "C" double  orient2d(double *pa, double *pb, double *pc);
extern "C" void exactinit();

/*
 * The following are data structures for fundamental geometric entities used
 * in creating hierarchical filamentous networks.
 */

/*
 * 2D Cartesian point 
 * Points are considered equal if their coordinates are within  a given
 * numerical tolerance.
 */
struct Point {

    Point(){
        this->tol = FLOAT_TOL;
    }

    Point(double xval, double yval) : x(xval), y(yval) {
        this->tol = FLOAT_TOL;
    }

    Point(double xval, double yval, int i) : x(xval), y(yval), index(i) {
        this->tol = FLOAT_TOL;
    }

    double x, y, tol;
    int index;

    //Utility functions establishing equality and ordering relations
    bool operator == (const Point& point) const {
        return(abs(point.x - x ) < tol && abs(point.y - y) < FLOAT_TOL);
    }

    friend bool operator < (const Point &p1, const Point &p2){
        if(p2.y > p1.y + p1.tol) return true;
        else if(p1.y > p2.y + p2.tol) return false;
        else if(p2.x > p1.x + p1.tol) return true;
        else return false;
    }
};

//Description of the bounding box of a geometric object or group of such objects
struct Bounds2D {

    double minx, miny, maxx, maxy, cx, cy, xrange, yrange;

    Bounds2D(){}

    Bounds2D(double left, double low, double right, double high){
        minx = left;
        miny = low;
        maxx = right;
        maxy = high;
        cx = (minx + maxx) / 2;
        cy = (miny + maxy) / 2;
        xrange = maxx - minx;
        yrange = maxy - miny;
    }

    bool in_bounds(Point query){
        return query.x >= minx && query.x <= maxx && query.y >= miny && query.y <= maxy;
    }
};

/*
 * An edge specified as ordered pairs of unsigned integer indices.
 */
struct Edge {

    Edge(){}

    Edge(unsigned int i1, unsigned int i2) : idx1(i1), idx2(i2) {}

    unsigned int idx1, idx2;

    bool operator == (const Edge& edge) const {
        return idx1 == edge.idx1 && idx2 == edge.idx2;
    }

    friend bool operator < (const Edge &e1, const Edge &e2){
        if(e2.idx1 > e1.idx1) return true;
        else if(e1.idx1 > e2.idx1) return false;
        else if(e2.idx2 > e1.idx2) return true;
        else return false;
    }
};

/*
 * A 2D network specified as a set of points and a set of edges, with edges
 * specified by ordered pairs of indices into the list of points
 */
struct NetworkComplex {

    vector<Point> points;
    vector<Edge> edges;

    NetworkComplex(){}

    NetworkComplex(vector<Point> my_points, vector<Edge> my_edges){
        points = my_points;
        edges = my_edges;
    }

    void clear(){
        points.clear();
        edges.clear();
    }

    void assign(NetworkComplex other){
        points.assign(other.points.begin(), other.points.end());
        edges.assign(other.edges.begin(), other.edges.end());
    }
};

/*
 * A polygon described as a set of edges arranged in CCW order
 * This class is only designed to handle convex polygons with no holes.
 */
struct NetPolygon {

    vector<int> vertices;
    double left, right, low, high;
    int index;

    NetPolygon(vector<Point> points, vector<int> my_verts, int my_idx){
        left = INF;
        right = -INF;
        low = INF;
        high = -INF;
        vertices = my_verts;
        index = my_idx;
        double x, y;

        for(int v : vertices){
            x = points[v].x;
            y = points[v].y;

            left = x < left ? x : left;
            right = x > right ? x : right;
            low = y < low ? y : low;
            high = y > high ? y : high;
        }
    }

    Bounds2D get_bounds(){
        return Bounds2D(left, low, right, high);
    }
};

/*
 * A collection of interlocking polygons, some of which share vertices and
 * edges
 */
struct PolygonComplex {

    vector<Point> points;
    set<Edge> edges;
    vector<NetPolygon> polygons;

    PolygonComplex(){}

    void clear(){
        points.clear();
        edges.clear();
        polygons.clear();
    }
};

/*
 * Given a polygon, specified as a set of ordered edges, determine whether
 * a query point is contained within that polygon.
 */
bool point_in_polygon(vector<Point> points, NetPolygon poly, Point query);

namespace std {

    template<> struct hash<Point>{
        typedef size_t result_type;
        typedef Point argument_type;
        size_t operator() (const Point& p) const;
    };
}

bool yesno(string message);

vector<string> split(string input, char delim);

vector<double> parse_doubles(vector<string> numstring);

vector<double> getdoubles(string prompt);

void makeunion(vector<int>& setvec, int root1, int root2);

int find_root(vector<int> setvec, int elem);

char get_choice(string message, map<char, string> c_map, vector<char> choices);

string enter_decline(string message);

void open_dat_file(string prompt, ifstream& file_stream);

void open_output_file(string prompt, ofstream& file_stream);

//Perform a uniform displacement of all points in a set
void displace(vector<Point> &points, double xdisp, double ydisp);

//Determine whether a line segment with endpoints p1 and p2 intersects a line
//segment with endpoints p3 and p4
bool intersection(Point p1, Point p2, Point p3, Point p4);

//Rotate each point in a set by a common angle about a common pivot
void rotate_points(vector<Point> &points, double angle, Point pivot);

void rotate_point(Point &point, double angle, Point pivot);

void get_extremes(vector<Point> points, double &minx, double &miny, double &maxx, double &maxy);

double x_intersect(vector<Point> points, Edge e, double yval);

double y_intersect(vector<Point> points, Edge e, double xval);

//Compute the squared distance between two points
double distance_sq(Point p1, Point p2);

//Utility function to generate a random seed
unsigned get_random_seed();

//Means of importing rules to create a lattice-based network
bool import_lattice(string name, vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns, double scale);

//Given a set of points and edges connecting those points, partition edges into
//two subsets: a random set forming a minimum spanning tree, and all others.
void randomMST(NetworkComplex nc, vector<Edge> &keep, vector<Edge>& rejects);

/*
This function adds Gaussian random noise to the location of each point.
The function takes as arguments the original set of points and edges describing
the network, and a standard deviation for Gaussian random noise. Points are 
shifted according to a normal distribution.
*/
void displace_points_grn(NetworkComplex &nc, double sdev);

//Given a set of rules for the geometry and topology of a lattice-based network,
//construct the network row-by-row.
void makeedges(vector<vector<double>> rules, map<int,vector<vector<double>>> nns, vector<double> bnds, NetworkComplex &nc);

template <typename T>
void add_if_missing(const T &t, map<T, int> &idx_map, vector<T> &list);

void add_thickness(NetworkComplex &current, double thickness, PolygonComplex &pc, bool makepoly, bool level);

NetworkComplex random_connected(NetworkComplex in, double to_keep);

//After edges have been removed from a network, cull unused points and reassign
//edge endpoint indices
void reassign_points_edges(NetworkComplex &nc, vector<Edge> kept_edges);
#endif
