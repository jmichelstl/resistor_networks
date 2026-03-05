#include "hnet_bvh.hpp"
#include <cmath>
#include <algorithm>

using namespace std;

//Find a union between two bounding boxes
Bounds2D box_union(Bounds2D box1, Bounds2D box2){

    double minx, maxx, miny, maxy;

    minx = fmin(box1.minx, box2.minx);
    maxx = fmax(box1.maxx, box2.maxx);
    miny = fmin(box1.miny, box2.miny);
    maxy = fmax(box1.maxy, box2.maxy);

    return Bounds2D(minx, maxx, miny, maxy);
}

struct BVHBuildNode {

    void init_internal(BVHBuildNode *left_child, BVHBuildNode *right_child){
        n_polys = 0;
        children[0] = left_child;
        children[1] = right_child;
        bounds = box_union(left_child->bounds, right_child->bounds);
    }
    
    void init_leaf(int my_first_poly, int my_n_polys, const Bounds2D &b) {
        first_poly = my_first_poly;
        n_polys = my_n_polys;
        children[0] = nullptr;
        children[1] = nullptr;
        bounds = b;
    }

    Bounds2D bounds;
    BVHBuildNode *children[2];
    int first_poly, n_polys;
};

int32_t left_shift_2(int32_t val){

    val = (val | (val << 8)) & 0b00000000111111110000000011111111;
    val = (val | (val << 4)) & 0b00001111000011110000111100001111;
    val = (val | (val << 2)) & 0b00110011001100110011001100110011;
    val = (val | (val << 1)) & 0b01010101010101010101010101010101;
    return val;
}

//Compute Morton codes for sorting of polygons
int32_t morton_code(Bounds2D bounds, NetPolygon poly){

    double cx, cy;
    int32_t xcode, ycode;

    cx = (poly.left + poly.right) / 2;
    cy = (poly.low + poly.high) / 2;
    xcode = (int32_t) floor(cx - bounds.minx) / bounds.xrange;
    if(xcode > MAX_IDX) xcode --;
    ycode = (int32_t) floor(cy - bounds.miny) / bounds.yrange;
    if(ycode > MAX_IDX) ycode --;
    xcode = left_shift_2(xcode);
    ycode = left_shift_2(ycode);

    //Spread x and y codes' bits out so that they are interspersed by zeros,
    //shift each y code bit one spot to the left, then combine the result
    //to obtain the final code.
    return xcode | (ycode << 1);
}

//Struct to handle sorting of primitives according to their Morton codes
struct MortonInfo {

    MortonInfo(){}

    MortonInfo(int poly_idx, uint32_t morton_code) : poly_idx(poly_idx), morton_code(morton_code) {}

    int poly_idx;
    uint32_t morton_code;
};

//Struct to support partition polygons according to the value of a given bit in
//a Morton encoding
struct bitcomp { 

    bitcomp(uint32_t mask) : mask(mask) {}
    
    inline bool operator()(const MortonInfo &arg) {
        return arg.morton_code & mask;
    }

    uint32_t mask;
};  

//Partition a sorted list of primitives recursively to build a tree
BVHBuildNode * HNetBVH::bvh_build(vector<NetPolygon> &sorted, vector<MortonInfo> &m_info, int start, int end, int bit, int *node_count){

    BVHBuildNode *node = (BVHBuildNode *) malloc(sizeof(BVHBuildNode));
    *node_count = *node_count + 1;

    //If the range between the start and end indices is 1, make a leaf
    if(end - start == 1){
        node->init_leaf(sorted.size(), 1, polygons[m_info[start].poly_idx].get_bounds());
        sorted.push_back(polygons[m_info[start].poly_idx]);
    }

    //Try to split the current subset of primitives along some axis. If
    //splitting is not achieved, place all nodes in a leaf.
    else{
        auto iter1 = next(m_info.begin(), start);
        auto iter2 = next(m_info.begin(), end);
        do{
            auto split = partition(iter1, iter2, bitcomp(1 << (bit-1)));
            if( split != iter1 && split != iter2){
                node->init_internal(bvh_build(sorted, m_info, start, start+distance(iter1, split), bit-1, node_count), bvh_build(sorted, m_info, start+distance(iter1,split), end, bit-1, node_count));
                return node;
            }

            bit --;
        }while(bit > 0);

        //If the loop terminates, no split was achieved, and all primitives
        //should be placed in a new leaf
        Bounds2D bounds = polygons[m_info[start].poly_idx].get_bounds();
        sorted.push_back(polygons[m_info[start].poly_idx]);
        for(auto iter = next(iter1); iter != iter2; iter ++){
            bounds = box_union(bounds,polygons[(*iter).poly_idx].get_bounds());
            sorted.push_back(polygons[(*iter).poly_idx]);
        }
        node->init_leaf(sorted.size()+start-end, end - start, bounds);
    }

    return node;
}

int HNetBVH::flatten_tree(BVHBuildNode *node, int *offset){

    BVHSearchNode *s_node = &s_nodes[*offset];
    s_node->bounds = node->bounds;
    int this_offset = (*offset)++;

    if(node->n_polys > 0){
        s_node->n_polys = node->n_polys;
        s_node->first_poly = node->first_poly;
    }

    else{
        s_node->n_polys = 0;
        flatten_tree(node->children[0], offset);
        s_node->secondChildOffset = flatten_tree(node->children[1], offset);
    }

    free(node);
    return this_offset;
}

void HNetBVH::index_polygons() {

    int node_count = 0, offset = 0, iter;
    vector<NetPolygon> sorted;
    vector<MortonInfo> m_info(polygons.size(), MortonInfo());

    //Compute Morton encodings
    #pragma omp parallel for
    for(iter = 0; iter < polygons.size(); iter ++){
        m_info[iter].poly_idx = iter;
	m_info[iter].morton_code = morton_code(bounds, polygons[iter]);
    }

    //Create the initial search tree with build nodes
    BVHBuildNode *root = bvh_build(sorted, m_info, 0, m_info.size(), 32, &node_count);
    polygons = sorted;

    //Create the more efficient tree representation with search nodes
    for(int iter = 0; iter < node_count; iter++){
        s_nodes.push_back(BVHSearchNode());
    }
    flatten_tree(root, &offset);
}

void HNetBVH::search_recursive(vector<int> &hits, Point query, int node_idx){

    Bounds2D bounds = s_nodes[node_idx].bounds;
    if(s_nodes[node_idx].n_polys > 0){
        int poly_start = s_nodes[node_idx].first_poly;
        int poly_end = poly_start + s_nodes[node_idx].n_polys;
        for(int poly_idx = poly_start; poly_idx < poly_end; poly_idx++){
            if(point_in_polygon(points, polygons[poly_idx], query)){
                hits.push_back(polygons[poly_idx].index);
            }
        }
    }

    else{
        if(s_nodes[node_idx+1].bounds.in_bounds(query)){
            search_recursive(hits, query, node_idx+1);
        }
        int offset = s_nodes[node_idx].secondChildOffset;
        if(s_nodes[offset].bounds.in_bounds(query)){
            search_recursive(hits, query, offset);
        }
    }
}

HNetBVH::HNetBVH(vector<Point> my_points, vector<NetPolygon> my_polys){

    points = my_points;
    polygons = my_polys;
    double maxx, maxy;
    get_extremes(points, minx, miny, maxx, maxy);
    xrange = maxx - minx;
    yrange = maxy - miny;
    bounds = Bounds2D(minx, miny, maxx, maxy);

    index_polygons();
}

vector<int> HNetBVH::within(Point query){

    vector<int> hits;

    if(s_nodes.size() == 0) return hits;

    search_recursive(hits, query, 0);

    return hits;
}
