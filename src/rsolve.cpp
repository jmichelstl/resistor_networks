#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include "network_utils.hpp"
#include <cmath>
#include <suitesparse/cholmod.h>
#include <suitesparse/SuiteSparseQR.hpp>
#include <unistd.h>
#include <stdlib.h>
#include <iomanip>
#include <stdlib.h>
#include <string.h>

using namespace std;

#define RSUCCESS 1
#define RFAILURE 0
#define TRUE 1
#define FALSE 0

//This is a data structure to represent a single ohmic resistor. The structure
//specifies a resistance and the neighboring node to which the resistor
//connects.
struct RPair{

    RPair(int idx, double rval) : index(idx), resistance(rval){}

    RPair(int idx) : index(idx) {
        resistance  = 1;
    }

    int index;
    double resistance;
};

//This structure supports a means of storing and compairing pairs of
//integer indices, for the purpose of creating and updating entries in a
//coefficient matrix used to solve for the potential field of a resistor
//network.
struct IPair{

    IPair(int i1, int i2) : index1(i1), index2(i2) {}

    int index1, index2;

    friend bool operator < (const IPair &ip1, const IPair &ip2){
        if(ip1.index1 < ip2.index1) return true;
        if(ip1.index1 > ip2.index1) return false;
        if(ip1.index2 < ip2.index2) return true;
        return false;
    }
};

inline bool in_range(int idx, int max){
    return idx >= 0 && idx < max;
}

//Read a description of a network as a planar line graph, and create a
//representation of a network as a planar line graph.
bool read_resistors(ifstream &datfile, vector<Point>& point_list, map<int, vector<RPair>>& resistor_map, double &ymin, double &ymax){

    string nextline;
    Point p1, p2;
    int num_pts, num_edges, edge_count = 0, index1, index2, mindex, maxdex;
    double x, y;

    ymin = FLT_MAX;
    ymax = FLT_MIN;

    //First, try to read a header specifying the number of vertices and edges
    if(! datfile.eof()){
        getline(datfile, nextline);
	if(sscanf(nextline.c_str(), "%d %d", &num_pts, &num_edges) < 2){
            cerr << "Error: could not read a valid header.\n";
	    datfile.close();
	    return false;
        }
	if(num_pts <= 0 || num_edges <= 0){
            cerr << "Error: an invalid header was encountered.\n";
	    datfile.close();
	    return false;
        }
    }
    else{
        cerr << "The network description file was empty.\n";
	datfile.close();
	return false;
    }

    //Next, attempt to read the vertices of the graph
    while(! datfile.eof() && point_list.size() < num_pts){
        getline(datfile, nextline);

        if(sscanf(nextline.c_str(), "%lf %lf", &x, &y) < 2){
            cerr << "An valid point could not be read.\n";
	    datfile.close();
	    return false;
        }

	if(y < ymin) ymin = y;
        if(y > ymax) ymax = y;

        point_list.push_back(Point(x, y));
	resistor_map.insert(make_pair(resistor_map.size(), vector<RPair>()));
    }

    if(resistor_map.size() < num_pts){
        cerr << "Too few points were read.\n";
	datfile.close();
	return false;
    }

    //Finally, attempt to read the edges
    while(! datfile.eof() && edge_count < num_edges){

        getline(datfile, nextline);
        if(sscanf(nextline.c_str(), "%d %d", &index1, &index2) < 2){
            cerr << "A pair of indices could not be read.\n";
	    datfile.close();
	    return false;
        }
	else if(! in_range(index1, num_pts) || ! in_range(index2, num_pts)){
            cerr << "An out-of-bounds index was discovered.\n";
	    datfile.close();
	    return false;
        }

        //If two valid indices have been found, add a resistor
        resistor_map[min(index1, index2)].push_back(RPair(max(index1, index2)));
	edge_count ++;
    }

    datfile.close();

    if(edge_count < num_edges){
        cerr << "Too few edges were read.\n";
	return false;
    }

    return true;
}

//Manage the indices within the arrays backing a cholmod_triplet structure
//at which information for matrix elements is stored
int get_triplet_index(IPair index_pair, int itype, void *Gi, void *Gj, void *Gx, map<IPair, int> &index_map, int &count){

    if(index_map.find(index_pair) == index_map.end()){
        if(itype == CHOLMOD_LONG){
            ((SuiteSparse_long *) Gi)[count] = index_pair.index1;
            ((SuiteSparse_long *) Gj)[count] = index_pair.index2;
        }
        else{
            ((int *) Gi)[count] = index_pair.index1;
            ((int *) Gj)[count] = index_pair.index2;
        }
        ((double *) Gx)[count] = 0;
        index_map.insert(make_pair(index_pair, count));
        count ++;
        return count - 1;
    }

    return index_map[index_pair];
}

//Helper function to update arrays of row and column indices within triplet
//sparse matrix data structures
void assign_index(void *index_list, int itype, int index, int value){
    if(itype == CHOLMOD_LONG){
        ((SuiteSparse_long *) index_list)[index] = value;
    }
    else{
        ((int *) index_list)[index] = value;
    }
}

//Helper function to update an array of matrix elements within a triplet sparse
//matrix structure
void assign_mat_elem(void *elem_list, int index, double value){
    ((double *) elem_list)[index] = value;
}

//Helper function to update members of a list of matrix elements within a
//triplet sparse matrix strucutre
void update_mat_elem(void *elem_list, int index, double value){
    ((double *) elem_list)[index] += value;
}

//Create a conductivity matrix that relates a potential field to a current
//distribution. Also make a note of the boundary conditions for the final
//equation used to solve for the potential field of the resistor network, and
//the nodes on the top and bottom of the network
int load_cond_mat(cholmod_sparse **cond_mat, cholmod_sparse **int_proj, cholmod_sparse **b_proj, cholmod_sparse **t_proj, cholmod_dense **bc_vec, cholmod_common *common, vector<Point> point_list, map<int, vector<RPair>> resistor_map, double ymin, double ymax){

    size_t point_count = point_list.size(), nzmax;
    int count = 0, index1, index2, mat_index, diag_index1, diag_index2;
    int int_count = 0, b_count = 0, t_count = 0, g_itype, b_itype, t_itype;
    int int_itype;
    Point p1, p2;
    double mat_elem, resistance;
    map<IPair, int> index_map;
    //Conductivity matrix, in sparse triplet form
    cholmod_triplet *cond_mat_triplet, *b_proj_triplet, *int_proj_triplet;
    cholmod_triplet *t_proj_triplet;
    //References to the arrays of the conductivity matrix, in sparse triplet
    //format
    void *Gi, *Gj, *Gx, *Bi, *Bj, *Bx, *Ii, *Ij, *Ix, *Ti, *Tj, *Tx, *BCx;

    //Find the total number of non-zero elements in the conductivity matrix
    nzmax = point_count;
    for(auto iter = resistor_map.begin(); iter != resistor_map.end(); iter++){
        nzmax += iter->second.size();
    }

    //Create the stiffness matrix and projection operators onto the interior
    //and boundary points
    cond_mat_triplet = cholmod_l_allocate_triplet(point_count, point_count, nzmax, -1, CHOLMOD_REAL, common);
    if(common->status < CHOLMOD_OK){
        cerr << "Could not allocate conductivity matrix.\n";
        return RFAILURE;
    }
    Gi = cond_mat_triplet->i;
    Gj = cond_mat_triplet->j;
    Gx = cond_mat_triplet->x;
    cond_mat_triplet->nnz = nzmax;
    g_itype = cond_mat_triplet->itype;

    b_proj_triplet = cholmod_l_allocate_triplet(point_count, point_count, point_count, -1, CHOLMOD_REAL, common);
    if(common->status < CHOLMOD_OK){
        cerr << "Could not allocate bottom projection matrix.\n";
        cholmod_l_free_triplet(&cond_mat_triplet, common);
        return RFAILURE;
    }
    Bi = b_proj_triplet->i;
    Bj = b_proj_triplet->j;
    Bx = b_proj_triplet->x;
    b_itype = b_proj_triplet->itype;
 
    t_proj_triplet = cholmod_l_allocate_triplet(point_count, point_count, point_count, -1, CHOLMOD_REAL, common);
    if(common->status < CHOLMOD_OK){
        cerr << "Could not allocate top projection matrix.\n";
        cholmod_l_free_triplet(&cond_mat_triplet, common);
        cholmod_l_free_triplet(&b_proj_triplet, common);
        return RFAILURE;
    }
    Ti = t_proj_triplet->i;
    Tj = t_proj_triplet->j;
    Tx = t_proj_triplet->x;
    t_itype = t_proj_triplet->itype;

    int_proj_triplet = cholmod_l_allocate_triplet(point_count, point_count, point_count, -1, CHOLMOD_REAL, common);
    if(common->status < CHOLMOD_OK){
        cerr << "Could not allocate interior projection matrix.\n";
        cholmod_l_free_triplet(&cond_mat_triplet, common);
        cholmod_l_free_triplet(&b_proj_triplet, common);
        cholmod_l_free_triplet(&t_proj_triplet, common);
        return RFAILURE;
    }
    Ii = int_proj_triplet->i;
    Ij = int_proj_triplet->j;
    Ix = int_proj_triplet->x;
    int_itype = int_proj_triplet->itype;

    *bc_vec = cholmod_l_zeros(point_count, 1, CHOLMOD_REAL, common);
    BCx = (*bc_vec)->x;

    //For each node in the network, create an entry along the diagonal of the
    //conductivity matrix for that node, as well as elements in the upper
    //diagonal for any resistor it shares with a point with a greater index
    for(auto it = resistor_map.begin(); it != resistor_map.end(); it++){
        index1 = it->first;
        p1 = point_list[index1];
        diag_index1 = get_triplet_index(IPair(index1, index1), g_itype, Gi, Gj, Gx, index_map, count);

        //Determine whether a point is on the top, bottom, or interior
        if(p1.y == ymin){
            assign_index(Bi, b_itype, b_count, index1);
            assign_index(Bj, b_itype, b_count, index1);
            assign_mat_elem(Bx, b_count, 1);
            b_count ++;

            //Stipulate that the potential on the bottom should be 1
            assign_mat_elem(BCx, index1, 1);
        }
        else if(p1.y == ymax){
            assign_index(Ti, t_itype, t_count, index1);
            assign_index(Tj, t_itype, t_count, index1);
            assign_mat_elem(Tx, t_count, 1);
            t_count ++;
        }
        else{
            assign_index(Ii, int_itype, int_count, index1);
            assign_index(Ij, int_itype, int_count, index1);
            assign_mat_elem(Ix, int_count, 1);
            int_count ++;
        }

        for(RPair rp : it->second){
            index2 = rp.index;
            resistance = rp.resistance;
            diag_index2 = get_triplet_index(IPair(index2, index2), g_itype, Gi, Gj, Gx, index_map, count);

            //According to the general formula for the conductivity matrix,
            //assign the element Gij the value 1 / Rij
            mat_elem = 1 / resistance;
            assign_index(Gi, g_itype, count, index1);
            assign_index(Gj, g_itype, count, index2);
            assign_mat_elem(Gx, count, mat_elem);
            count ++;

            //Along the main diagonal, subtract 1 / Rij from both Gii and Gjj
            update_mat_elem(Gx, diag_index1, -mat_elem);
            update_mat_elem(Gx, diag_index2, -mat_elem);
        }
    }

    int_proj_triplet->nnz = int_count;
    b_proj_triplet->nnz = b_count;
    t_proj_triplet->nnz = t_count;

    //Convert matrices from triplet to cholmod_sparse form, then free triplet
    //matrices
    *cond_mat = cholmod_l_triplet_to_sparse(cond_mat_triplet, nzmax, common);
    *int_proj = cholmod_l_triplet_to_sparse(int_proj_triplet, int_count,common);
    *b_proj = cholmod_l_triplet_to_sparse(b_proj_triplet, b_count, common);
    *t_proj = cholmod_l_triplet_to_sparse(t_proj_triplet, t_count, common);
    cholmod_l_free_triplet(&cond_mat_triplet, common);
    cholmod_l_free_triplet(&int_proj_triplet, common);
    cholmod_l_free_triplet(&b_proj_triplet, common);
    cholmod_l_free_triplet(&t_proj_triplet, common);

    return RSUCCESS;
}

//Prompt the user for a file describing a network as a planar line graph,
//read the contents of the file, build a conductivity matrix, and solve for
//the potential field of the nodes in the network. Finally, compute the current
//at the grounded end of the network, and report the effective resistance.
void resistance_run(ifstream &datfile, const char *report_name){

    vector<Point> point_list;
    map<int, vector<RPair>> resistor_map;
    cholmod_sparse *cond_mat, *int_proj, *b_proj, *t_proj, *masked_gmat, *lhs;
    cholmod_sparse *sparse_pfield, *all_current, *total_boundary, *tcurr;
    cholmod_dense *bc_vec, *pfield, *dense_tcurr;
    cholmod_common common;
    int point_count, status, iter;
    double ymin, ymax, current = 0;
    double unity[2] = {1,1};
    FILE *report;
    size_t total_mem, available_mem;

    //Read information about the nodes and connections in the network from
    //the specified file
    if(! read_resistors(datfile, point_list, resistor_map, ymin, ymax)){
        cerr << "Reading of the resistor network failed.\n";
	return;
    }
    point_count = point_list.size();

    //Initialize cholmod_common structure
    cholmod_l_start(&common);
    /*common.useGPU = true;
    cholmod_l_gpu_memorysize(&total_mem, &available_mem, &common);
    common.gpuMemorySize = available_mem;*/

    if(point_count == 0){
        cerr << "No edges were read.\n";
        return;
    }

    //Create the conductivity matrix
    status = load_cond_mat(&cond_mat, &int_proj, &b_proj, &t_proj, &bc_vec, &common, point_list, resistor_map, ymin, ymax);
    //If not all needed structures were successfully constructed, abort
    if(status == RFAILURE){
        cholmod_l_finish(&common);
        return;
    }

    //Set up the matrix equation to solve for the potential field
    masked_gmat = cholmod_l_ssmult(int_proj, cond_mat, 0, TRUE, TRUE, &common);
    total_boundary = cholmod_l_add(b_proj, t_proj, unity, unity, TRUE, TRUE, &common);
    lhs = cholmod_l_add(masked_gmat, total_boundary, unity, unity, TRUE, TRUE, &common);

    //Free matrices no longer needed for the rest of the computation
    cholmod_l_free_sparse(&masked_gmat, &common);
    cholmod_l_free_sparse(&int_proj, &common);
    cholmod_l_free_sparse(&b_proj, &common);
    cholmod_l_free_sparse(&total_boundary, &common);

    //Solve for the potential field of the nodes in the resistor network
    pfield = SuiteSparseQR <double>(lhs, bc_vec, &common);

    //Find the current distribution, and isolate for the net current flowing
    //into the top of the network
    sparse_pfield = cholmod_l_dense_to_sparse(pfield, TRUE, &common);
    all_current = cholmod_l_ssmult(cond_mat, sparse_pfield, 0, TRUE, TRUE, &common);
    tcurr = cholmod_l_ssmult(t_proj, all_current, 0, TRUE, TRUE, &common);
    dense_tcurr = cholmod_l_sparse_to_dense(tcurr, &common);

    //Sum up currents running into top nodes to find the total current flow
    //and total resistance
    for(iter = 0; iter < dense_tcurr->nrow; iter++){
        current += ((double *) (dense_tcurr->x))[iter];
    }

    //Optionally report the potential field to a file
    if(report_name){
        report = fopen(report_name, "w");
        if(! report){
            cerr << "The specified report file could not be opened.\n";
        }
        else{
            for(int iter = 0; iter < pfield->nzmax; iter++){
                fprintf(report, "%1.10le\n", ((double *) pfield->x)[iter]);
            }
            fclose(report);
        }
    }

    //TODO: Optionally report the conductivity matrix to a file. I'm not sure
    //I see the use in this option at the moment.

    //Clean house
    cholmod_free_sparse(&sparse_pfield, &common);
    cholmod_free_sparse(&tcurr, &common);
    cholmod_free_sparse(&cond_mat, &common);
    cholmod_free_sparse(&all_current, &common);
    cholmod_free_sparse(&t_proj, &common);
    cholmod_free_sparse(&lhs, &common);
    cholmod_l_free_dense(&bc_vec, &common);
    cholmod_l_free_dense(&pfield, &common);
    cholmod_l_free_dense(&dense_tcurr, &common);
    cholmod_l_finish(&common);

    printf("%1.10le\n", current);
}

int main(int argc, char **argv){
    ifstream netfile;
    bool interact = true;
    char *file_name, *output_name = NULL;
    double resistance;
    char c;

    while((c = getopt(argc, argv, "f:p:")) != -1){
        switch(c) {
            case 'f':
                netfile.open(string(optarg));
                if(! netfile.is_open()){
                    cerr << "The specified file could not be read.\n";
                }
                else interact = false;
                break;
            case 'p':
                output_name = (char *) malloc(sizeof(char)*(1+strlen(optarg)));
                strcpy(output_name, optarg);
                break;
            case '?':
                if(optopt == 'f'){
                    cerr << "Option \"f\" requires a file name. Interactive mode will be used.\n";
                }
                else if(optopt == 'p'){
                    cerr << "Option \"p\" requires an output file name.\n";
                }
                else if(isprint(optopt)){
                    fprintf(stderr, "Unknown option -%c.\n", optopt);
                }
                else cerr << "Unknown option character.\n";
                break;
            default:
                break;
        }
    }

    if(interact){
        do{
            string response;
            open_dat_file("Enter the name of a network file: ", netfile);
            if(yesno("Report the potential field?")){
                cout << "Enter the report file name: ";
                getline(cin, response);
            }
            if(response.compare("")) resistance_run(netfile,response.c_str());
            else resistance_run(netfile, NULL);
        }while(yesno("Analyze another network?"));
    }

    else{
        resistance_run(netfile, (const char *) output_name);
    }
    if(output_name) free(output_name);

    return 0;
}
