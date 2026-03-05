/*
 * This code is meant to test the BVH utility for testing which polygons, if
 * any, contain a given query point.
 */

#include <iostream>
#include <stdio.h>
#include "network_utils.hpp"
#include "hnet_bvh.hpp"
#include <boost/random.hpp>
#include <memory>

int main(int argc, char **argv){

    unsigned seed;
    FILE *ranfile, *results;
    vector<Point> poly_points, all_queries, hits1, hits2, outside;
    set<Edge> square_edges;
    vector<NetPolygon> squares;
    PolygonComplex pc;

    exactinit();

    ranfile = fopen("/dev/urandom", "r");

    if(ranfile == NULL){
        cerr << "Opening /dev/urandom failed. The RNG could not be seeded.\n";
        return -1;
    }

    //Random number generating infrastructure
    size_t size = fread(&seed, sizeof(unsigned), 1, ranfile);
    fclose(ranfile);
    if(size < 1){
        cerr << "Reading from /dev/urandom failed. No seed was obtained.\n";
        return -2;
    }

    boost::random::lagged_fibonacci9689 dgen{static_cast<uint32_t>(seed)};

    //Create two squares to test membership queries

    //Vertices and edges for square 1
    poly_points.push_back(Point(1, 1));
    poly_points.push_back(Point(2, 1));
    poly_points.push_back(Point(1, 2));
    poly_points.push_back(Point(2, 2));
    square_edges.insert(Edge(0, 1));
    square_edges.insert(Edge(1, 3));
    square_edges.insert(Edge(3, 2));
    square_edges.insert(Edge(2, 0));

    //Vertices and edges for square 2
    poly_points.push_back(Point(7, 7));
    poly_points.push_back(Point(8, 7));
    poly_points.push_back(Point(7, 8));
    poly_points.push_back(Point(8, 8));
    square_edges.insert(Edge(4, 5));
    square_edges.insert(Edge(5, 7));
    square_edges.insert(Edge(7, 6));
    square_edges.insert(Edge(6, 4));

    squares.push_back(NetPolygon(poly_points, vector<int>({0, 1, 3, 2}), 0));
    squares.push_back(NetPolygon(poly_points, vector<int>({4, 5, 7, 6}), 1));

    pc.points = poly_points;
    pc.edges = square_edges;
    pc.polygons = squares;

    //Randomly generate a set of query points
    for(int i = 0; i < 100000; i++){
        all_queries.push_back(Point(10*dgen(), 10*dgen()));
    }

    //Create the spatial index and determine point membership
    shared_ptr<HNetBVH> hnb = make_shared<HNetBVH>(poly_points, squares);

    for(Point p : all_queries){
        vector<int> membership = hnb->within(p);
	if(membership.size() == 0){
            outside.push_back(p);
            continue;	    
	}
	if(membership[0] == 0){
            hits1.push_back(p);
        }
	else hits2.push_back(p);
    }

    //Report membership results
    results = fopen("in1.txt", "w");
    for(Point p : hits1){
        fprintf(results, "%1.8lf\t%1.8lf\n", p.x, p.y);
    }
    fclose(results);

    results = fopen("in2.txt", "w");
    for(Point p : hits2){
        fprintf(results, "%1.8lf\t%1.8lf\n", p.x, p.y);
    }
    fclose(results);

    results = fopen("outside.txt", "w");
    for(Point p : outside){
        fprintf(results, "%1.8lf\t%1.8lf\n", p.x, p.y);
    }
    fclose(results);

    return 0;
}
