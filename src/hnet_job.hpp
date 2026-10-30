#ifndef HNET_JOB
#define HNET_JOB

#include <vector>
#include <string>
#include "isle.hpp"
#include <memory>
#include <string.h>
#include <cstdint>
#include <cfloat>
#include <sstream>
#include <cmath>
#include <climits>
#include <typeinfo>

//Default number of digits of precision to use in reporting bond portions
#define DEFAULT_BP_PRECISION 3

using namespace std;

//Type to distinguish between three types of output: exportation as a 
//network, specified by point geometry and topology, exportation as an SVG
//file, in which the network has grips added for conductivity tests, or
//exportation as an STL file, in which a surface is created by extruding the
//network out of the plane.
enum class output_mode {network, svg, stl_extrude};

//Utility function for validating an integer argument
bool read_int(string name, string arg, int &value, int min, int max){

    int num_read;

    if(! (num_read = sscanf(arg.c_str(), "%d", &value))){
        cerr << "Parameter " << name << " must be an integer.\n";
        return false;
    }
    if( value < min || value > max){
        cerr << "Out of bounds argument for parameter \"" << name << "\".\n";
        return false;
    }

    return true;
}

//Utility function for validating a float argument
bool read_float(string name, string arg, float &value, float min, float max){

    int num_read;

    if(! (num_read = sscanf(arg.c_str(), "%f", &value))){
        cerr << "Parameter " << name << " must be a real number.\n";
        return false;
    }
    if( value < min || value > max){
        cerr << "Out of bounds argument for parameter \"" << name << "\".\n";
        return false;
    }

    return true;
}

//Utility function for validating a double argument
bool read_double(string name, string arg, double &value, double min, double max){

    int num_read;

    if(! (num_read = sscanf(arg.c_str(), "%lf", &value))){
        cerr << "Parameter " << name << " must be a real number.\n";
        return false;
    }
    if( value < min || value > max){
        cerr << "Out of bounds argument for parameter \"" << name << "\".\n";
        return false;
    }

    return true;
}

//Utility function for validating a boolean argument
bool read_bool(string name, string arg, bool &value){

    if(arg.compare("false") == 0){
        value = false;
        return true;
    }

    if(arg.compare("true") == 0){
        value = true;
        return true;
    }

    return false;
}

bool check_for_attribute(shared_ptr<isle::Element> elem, string name, string &value, bool needed){

    if(! elem->query_first_attribute(name, value)){
        if(needed){
            cerr << "Error in parsing a layer: missing " << name << ".\n";
        }
        return false;
    }

    return true;
}

//Make sure a parameter is specified, and that its value is within a
//specified range.
bool validate_parameter(shared_ptr<isle::Element> elem, string name, double &param, double min_val, double max_val){

    string result;

    if(! check_for_attribute(elem, name, result, true)){
        cerr << "The parameter " << name << " was not specified.\n";
        return false;
    }
    else {
        if(! read_double(name, result, param, min_val, max_val)){
            cerr << name << " must be between " << min_val << " and " << max_val << "." << endl;
            return false;
        }
    }

    return true;
}

bool validate_parameter(shared_ptr<isle::Element> elem, string name, float &param, float min_val, float max_val){

    string result;

    if(! check_for_attribute(elem, name, result, true)){
        cerr << "The parameter " << name << " was not specified.\n";
        return false;
    }
    else {
        if(! read_float(name, result, param, min_val, max_val)){
            cerr << name << " must be between " << min_val << " and " << max_val << "." << endl;
            return false;
        }
    }

    return true;
}

//Representation of a single stage in a hierarchical network creation process
struct HNLayer {

    HNLayer (){
        lattice_file = "";
        scale = 0;
        width = 0;
        stitch_cutoff = 0;
        bond_occupation = 0;
        displacement = 0;
        grain_based = false;
    }

    string lattice_file;
    double scale, width, bond_occupation, displacement, stitch_cutoff;
    bool protect_top_bottom, grain_based;
};

//Representation of a set of instructions for creating a hierarchical network
struct HNTask {

    HNTask (){
        output_file = "";
        poly_file = "";
        skeleton_file = "";
        top_file = "";
        alignment_file = "";
        bp_prod = 1;
        displaced = false;
        connected = true;
        connect_top_bottom = false;
        add_grips = false;
    }

    string poly_file, skeleton_file, top_file, alignment_file, output_file;
    bool displaced, connected, connect_top_bottom, add_grips;
    double bp_prod, minx, miny, maxx, maxy, grip_lower, grip_upper, taper_width;
    double taper_height, corner_radius, tab_height;
    float ext_height;
    vector<HNLayer> layers;
};

//Representation of a set of networks to be created, sharing common lattice
//file names and scale factors, but differing in bond portions
struct HNJob {

    HNJob(){}

    vector<string> lattice_files;
    vector<double> scales;
    vector<HNTask> tasks;
};

//Append descriptive tags to output file names to indicate bond portions
void augment_file_names(HNTask &task, double bp, int num_digits){

    ostringstream oss;
    int bp_rep, quotient, padding = num_digits;

    bp_rep = (int) (round(bp * pow(10, num_digits - 1)));
    quotient = bp_rep;
    while(quotient){
        padding --;
        quotient /= 10;
    }
    oss << "_";
    for(int iter = 0; iter < padding; iter ++) oss << "0";
    oss << bp_rep;

    if(task.poly_file.compare("") != 0) task.poly_file += oss.str();
    if(task.skeleton_file.compare("") != 0) task.skeleton_file += oss.str();
    if(task.top_file.compare("") != 0) task.top_file += oss.str();
    if(task.alignment_file.compare("") != 0) task.alignment_file += oss.str();
    if(task.output_file.compare("") != 0) task.output_file += oss.str();
}

//Attempt to obtain the instructions for a single layer of a hierarchical
//network. Upon success, add a representation of the layer to a list.
bool parse_hnet_layer(shared_ptr<isle::Element> layer, HNJob &job, int num_digits){

    int num_read, iter;
    double min_bp, max_bp, bp_inc, scale;
    HNLayer new_layer;
    string result;

    if(! check_for_attribute(layer, string("lattice_file"), result, true)){
        return false;
    }
    else {
        job.lattice_files.push_back(result);
    }

    if(! check_for_attribute(layer, string("scale"), result, true)){
        return false;
    }
    else {
        if(! read_double(string("scale"), result, scale, FLT_MIN, FLT_MAX)){
            cerr << "The scale must be a positive real number.\n";
            return false;
        }
        else job.scales.push_back(scale);
    }

    if(! check_for_attribute(layer, string("width"), result, true)){
        return false;
    }
    else {
        if(! read_double(string("width"), result, new_layer.width, 0, FLT_MAX)){
            cerr << "The width must be a non-negative real number.\n";
            return false;
        }
    }

    if(check_for_attribute(layer, string("displacement"), result, false)){
        if(! read_double(string("displacement"), result, new_layer.displacement, FLT_MIN, FLT_MAX)){
            cerr << "The displacement must be a positive real number.\n";
            return false;
        }
    }

    if(check_for_attribute(layer, string("protect_top_bottom"), result, false)){
        if(! read_bool(string("protect_top_bottom"), result, new_layer.protect_top_bottom)){
            cerr << "The field protect_top_bottom must be true or false.\n";
            return false;
        }
    }

    if(check_for_attribute(layer, string("grain_based"), result, false)){
        if(! read_bool(string("grain_based"), result, new_layer.grain_based)){
            cerr << "The field grain_based must be true or false.\n";
            return false;
        }
    }

    //If a grain-based construction should be used, a cutoff length for edges
    //stitching together neighboring grains should be given.
    if(new_layer.grain_based){
        if(! check_for_attribute(layer, string("stitch_cutoff"), result, true)){
            return false;
        }
        else {
            if(! read_double(string("stitch_cutoff"), result, new_layer.stitch_cutoff, FLT_MIN, FLT_MAX)){
                cerr << "stitch_cutoff must be a positive real number.\n";
                return false;
            }
        }
    }

    if(! check_for_attribute(layer, string("bond_occupation"), result, true)){
        return false;
    }
    else {
        num_read = sscanf(result.c_str(), "%lf %lf %lf", &min_bp, &max_bp, &bp_inc);
        if(num_read == 1){
            if(min_bp < 0 || min_bp > 1){
                cerr << "The bond portion must be between 0 and 1, inclusive.\n";
                return false;
            }
            new_layer.bond_occupation = min_bp;
            for(iter = 0; iter < job.tasks.size(); iter ++) { 
                job.tasks[iter].layers.emplace(job.tasks[iter].layers.begin(),new_layer);
                job.tasks[iter].bp_prod *= min_bp;
                if(new_layer.displacement) job.tasks[iter].displaced = true;
            }
        }
        else if(num_read == 3){
            if(min_bp < 0 || max_bp <= min_bp || max_bp > 1 || bp_inc <= 0){
                cerr << "The bond_occupation argument is invalid.\n";
                return false;
            }
            vector<HNTask> expanded_tasks;
            for(new_layer.bond_occupation = min_bp; new_layer.bond_occupation <= max_bp; new_layer.bond_occupation += bp_inc){
                for(HNTask next_task : job.tasks) {
                    next_task.layers.emplace(next_task.layers.begin(),new_layer);
                    next_task.bp_prod *= new_layer.bond_occupation;
                    if(new_layer.displacement) next_task.displaced = true;
                    augment_file_names(next_task, new_layer.bond_occupation, num_digits);
                    expanded_tasks.push_back(next_task);
                }
            }
            job.tasks.assign(expanded_tasks.begin(), expanded_tasks.end());
        }
        else {
            cerr << "The argument for the bond_occupation field should be one or three non-negative real numbers, with the maximum value no greater than 1.\n";
            return false;
        }
    }
    return true;
}

string sequence_tag(int seq_num, int total_digits){

    ostringstream tag;
    int dec_digits = 0;

    for(int quotient = seq_num; quotient > 0; quotient /= 10){
        dec_digits ++;
    }

    for(int iter = 0; iter < total_digits - dec_digits; iter ++) tag << "0";
    tag << seq_num;

    return tag.str();
}

//Attempt to parse instructions for creating one or more hierarchical networks
//from a file, using the Input Serialization LanguagE(ISLE). If only
//one bond portion is to be used and one realization is to be made, file names
//will be used as they are. Otherwise, it will be assumed that file names are
//stubs, to which information should be appended specifying bond portions and
//realization numbers.
bool parse_hnet_job(string job_file_name, HNJob &job, output_mode mode){

    int num_read, realizations = 1, report_digits = DEFAULT_BP_PRECISION, lnum;
    double min_h, max_h, taper_w, taper_h, c_rad, tab_h;
    float height;
    shared_ptr<isle::Element> root;
    HNTask task;
    string result;
    vector<shared_ptr<isle::Element>> layers;
    string extension;

    //Choose the correct extension for reporting networks
    if(mode == output_mode::network){
        extension = ".dat";
    }
    else if(mode == output_mode::svg){
        extension = ".svg";
    }
    else if(mode == output_mode::stl_extrude){
        extension = ".stl";
    }

    //Attempt to parse a valid element tree, and make sure the root element is
    //of type hnet_job
    if(! isle::parse_isle(job_file_name, root)){
        return false;
    }

    if(root->get_name().compare("hnet_job")){
        cerr << "The root-level element must be of type hnet_job.\n";
        return false;
    }

    //Look for overall attributes

    //The output file and bounds, at a minimum, must be specified
    if(! root->query_first_attribute("output_name", task.output_file)){
        cerr << "No output name was specified.\n";
        return false;
    }

    if(root->query_first_attribute("bounds", result)){
        if((num_read = sscanf(result.c_str(), "%lf %lf %lf %lf", &task.minx, &task.miny, &task.maxx, &task.maxy)) != 4){
            cerr << "The bounds attribute should specify four real numbers.\n";
            return false;
        }
        if(task.minx >= task.maxx || task.miny >= task.maxy){
            cerr << "The bounds were mis-ordered.\n";
            return false;
        }
    }
    else {
        cerr << "No bounds were specified.\n";
        return false;
    }

    //Query for optional attributes
    root->query_first_attribute("poly_file", task.poly_file);
    root->query_first_attribute("skeleton_file", task.skeleton_file);
    root->query_first_attribute("top_file", task.top_file);
    root->query_first_attribute("alignment_file", task.alignment_file);

    if(root->query_first_attribute("report_digits", result)){
        if(! read_int(string("report_digits"), result, report_digits, 1, INT_MAX)){
            cerr << "Parsing failed.\n";
            return false;
        }
    }

    if(root->query_first_attribute("realizations", result)){
        if(! read_int(string("realizations"), result, realizations, 1, INT_MAX)){
            cerr << "Parsing failed.\n";
            return false;
        }
    }

    if(root->query_first_attribute("connected", result)){
        if(! read_bool(string("connected"), result, task.connected)){
            cerr << "Parsing failed.\n";
            return false;
        }
    }

    if(root->query_first_attribute("connect_top_bottom", result)){
        if(! read_bool(string("connect_top_bottom"), result, task.connect_top_bottom)){
            cerr << "Parsing failed.\n";
            return false;
        }
    }

    //Determine whether the network should be prepared for the addition of top 
    //and bottom grips.
    if(check_for_attribute(root, string("add_grips"), result, false)){
        if(! read_bool(string("add_grips"), result, task.add_grips)){
            cerr << "The field add_grips must be true or false.\n";
            return false;
        }
    }

    job.tasks.push_back(task);

    //Attempt to parse a set of layers. For each instance in which a range
    //of bond portions is specified, create multiple copies of all jobs.
    if(! root->query_children("layer", layers)){
        cerr << "No layers were specified.\n";
        return false;
    }

    lnum = 1;
    for(shared_ptr<isle::Element> next_layer : layers){
        if(! parse_hnet_layer(next_layer, job, report_digits)){
            cerr << "Parsing of layer " << lnum << " failed.\n";
            return false;
        }
        lnum ++;
    }

    //Make sure that, if large-scale or small-scale bonds have been displaced
    //and grips are desired, grip heights are spcified
    if((mode != output_mode::network || job.tasks[0].add_grips) && job.tasks[0].displaced){
        if(check_for_attribute(root, string("grip_heights"), result, true)){
            if((num_read = sscanf(result.c_str(), "%lf %lf", &min_h, &max_h)) != 2){
                cerr << "The grip heights field takes two numbers as its argument.\n";
                return false;
            }
            if(max_h <= min_h){
                cerr << "The maximum grip height must be greater than the minimum.\n";
                return false;
            }
            for(int iter = 0; iter < job.tasks.size(); iter ++){
                job.tasks[iter].grip_lower = min_h;
                job.tasks[iter].grip_upper = max_h;
            }
        }
        else {
            cerr << "If grips are desired for a network with vertex displacement, grip heights must be specified.\n";
            return false;
        }
    }

    //If the SVG mode is desired, information about grip taper sections and
    //tabs should be specified.
    if(mode == output_mode::svg) {
        if(! validate_parameter(root, "taper_width", taper_w, DBL_MIN, DBL_MAX)){
            return false;
        }
        if(! validate_parameter(root, "taper_height", taper_h, DBL_MIN, DBL_MAX)){
            return false;
        }
        if(! validate_parameter(root, "corner_radius", c_rad, DBL_MIN, DBL_MAX)){
            return false;
        }
        if(! validate_parameter(root, "tab_height", tab_h, DBL_MIN, DBL_MAX)){
            return false;
        }
        for(int iter = 0; iter < job.tasks.size(); iter ++){
            job.tasks[iter].taper_width = taper_w;
            job.tasks[iter].taper_height = taper_h;
            job.tasks[iter].corner_radius = c_rad;
            job.tasks[iter].tab_height = tab_h;
        }
    }

    else if(mode == output_mode::stl_extrude) {
        if(! validate_parameter(root, "extrusion_height", height, FLT_MIN, FLT_MAX)){
            return false;
        }
        for(int iter = 0; iter < job.tasks.size(); iter ++) {
            job.tasks[iter].ext_height = height;
        }
    }

    //Prepare for multiple realizations, if desired
    if(realizations > 1) {

        vector<HNTask> expanded_tasks;

        int total_digits = 0;
        for(int quotient = realizations; quotient > 0; quotient /= 10){
            total_digits ++;
        }
        vector<string> seq_tags;
        for(int iter = 1; iter <= realizations; iter ++){
            seq_tags.push_back(sequence_tag(iter, total_digits));
        }

        for(HNTask next_task : job.tasks){
            if(next_task.bp_prod < 1 || next_task.displaced){
                for(int iter = 0; iter < realizations; iter ++){
                    expanded_tasks.push_back(next_task);
                    if(next_task.poly_file.compare("") != 0){
                        (*expanded_tasks.rbegin()).poly_file += seq_tags[iter];
                    }
                    if(next_task.skeleton_file.compare("") != 0) {
                        (*expanded_tasks.rbegin()).skeleton_file += seq_tags[iter];
                    }
                    if(next_task.top_file.compare("") != 0) {
                        (*expanded_tasks.rbegin()).top_file += seq_tags[iter];
                    }
                    if(next_task.alignment_file.compare("") != 0) {
                        (*expanded_tasks.rbegin()).alignment_file += seq_tags[iter];
                    }
                    if(next_task.output_file.compare("") != 0) {
                        (*expanded_tasks.rbegin()).output_file += seq_tags[iter];
                    }
                }
            }
            else expanded_tasks.push_back(next_task);
        }
        job.tasks.assign(expanded_tasks.begin(), expanded_tasks.end());
    }

    //Finalize output file names.
    for(int iter = 0; iter < job.tasks.size(); iter ++){
        if(job.tasks[iter].poly_file.compare("") != 0){
            job.tasks[iter].poly_file += ".dat";
        }
        if(job.tasks[iter].skeleton_file.compare("") != 0) {
            job.tasks[iter].skeleton_file += ".dat";
        }
        if(job.tasks[iter].top_file.compare("") != 0) {
            job.tasks[iter].top_file += ".dat";
        }
        if(job.tasks[iter].alignment_file.compare("") != 0) {
            job.tasks[iter].alignment_file += ".dat";
        }
        if(job.tasks[iter].output_file.compare("") != 0) {
            job.tasks[iter].output_file += extension;
        }
    }

    return true;
}

#endif
