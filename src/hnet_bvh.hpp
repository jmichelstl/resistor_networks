#ifndef HNET_BVH
#define HNET_BVH

#include "network_utils.hpp"
#include <cstdint>
#include <cmath>

//Maximum quantized index along the x or y direction for computing a Morton
//code
#define MAX_IDX 0xFFFF

using namespace std;

//Node for building the BVH
struct BVHBuildNode;

struct MortonInfo;

//Compact node for searching the BVH after construction
struct BVHSearchNode {

    BVHSearchNode(){}

    Bounds2D bounds;
    union{
        int first_poly;
        int secondChildOffset;
    };
    unsigned int n_polys;
};

class HNetBVH {

    public:
        HNetBVH(vector<Point> points, vector<NetPolygon> polys);
	    vector<int> within(Point p);

    private:

        void index_polygons();

        BVHBuildNode *bvh_build(vector<NetPolygon> &sorted, vector<MortonInfo> &minfo, int start, int end, int bit, int *node_count);

        int flatten_tree(BVHBuildNode *node, int *offset);

        void search_recursive(vector<int> &hits, Point query, int node_idx);

        double minx, miny, xrange, yrange;
        vector<Point> points;
	    vector<NetPolygon> polygons;
        vector<BVHSearchNode> s_nodes;
	    Bounds2D bounds;
};

#endif
