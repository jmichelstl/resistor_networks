#include "network_utils.hpp"
#include <time.h>
#include <algorithm>

using namespace std;

/*
 * Given a polygon, specified as a set of ordered edges, determine whether
 * a query point is contained within that polygon.
 */
bool point_in_polygon(vector<Point> points, NetPolygon poly, Point query){

    double qx, qy, x1, y1, x2, y2;
    double p1[2], p2[2], q[2];

    q[0] = query.x;
    q[1] = query.y;

    //This code assumes edges are ordered in CCW order. Given this convention,
    //this test uses the Robust Geometric Predicates orient2d method, developed
    //by Jonathan Richard Shewchuck, to determine whether the query point is
    //to the "left", in a generalized sense, of each polygon edge.
    for(int i = 0; i < poly.vertices.size(); i++){
        p1[0] = points[poly.vertices[i]].x;
        p1[1] = points[poly.vertices[i]].y;
        p2[0] = points[poly.vertices[(i+1)%poly.vertices.size()]].x;
        p2[1] = points[poly.vertices[(i+1)%poly.vertices.size()]].y;
        if(orient2d(p1, p2, q) < -FLOAT_TOL) return false;
    }

    return true;
}

namespace std {

    size_t hash<Point>::operator() (const Point& p) const {
        return (size_t) ((p.x+p.y)*(p.x+p.y)*(p.x*p.y+p.x*p.x+p.y*p.y*PRIME1))^PRIME2;
    }
}

bool yesno(string message){
    string response;

    cout << message << "(y/n): ";

    while(true){
        getline(cin, response);
        if(!response.compare("y")) return true;
        else if(!response.compare("n")) return false;
        else printf("Respond with \"y\" or \"n\": ");
    }
}

vector<string> split(string input, char delim){
    vector<string> result;
    size_t start, iter, len;
    bool reading_token = false;

    for(iter = 0; iter < input.size(); iter++){
        if(delim != input[iter]){
            if(!reading_token){
                reading_token = true;
                start = iter;
            }
        }

        else{
            if(reading_token){
                reading_token = false;
                result.push_back(input.substr(start, iter - start));
            }
        }
    }

    if(reading_token) result.push_back(input.substr(start, iter - start));


    return result;
}

vector<double> parse_doubles(vector<string> numstring){
    double nextnum;
    vector<double> result;

    for(string s : numstring){
        if(sscanf(s.c_str(), "%lf", &nextnum)) result.push_back(nextnum);
    }

    return result;
}

vector<double> getdoubles(string prompt){
    double nextnum;
    vector<double> result;
    string response;

    printf("%s", prompt.c_str());
    getline(cin, response);

    return parse_doubles(split(response, ' '));
}

void makeunion(vector<int>& setvec, int root1, int root2){
    if(setvec[root2] < setvec[root1]) setvec[root1] = root2;
    else{
        if(setvec[root1] == setvec[root2]) setvec[root1]--;

        setvec[root2] = root1;
    }
}

int find_root(vector<int> setvec, int elem){
    int pos = elem;

    if(elem >= setvec.size()) return -1;

    while(setvec[pos] >= 0) pos = setvec[pos];

    return pos;
}

char get_choice(string message, map<char, string> c_map, vector<char> choices){

    bool valid = false;
    string response;
    char choice;

    do{
        cout << message;
        for(auto it = c_map.begin(); it != c_map.end(); it++){
            cout << it->first << ": " << it->second << "\n";
        }
        getline(cin, response);
        if(response.compare("") == 0){
            cerr << "Please make a choice.\n";
        }

        else{
            choice = response[0];
            for(char next_choice : choices){
                if(choice == next_choice) return choice;
            }
            cerr << "Choice \"" << choice << "\"" << "not recognized.\n";
        }
    }while(!valid);

    return choice;
}

string enter_decline(string message){
    string response;
    cout << message << ", or return to decline: ";
    getline(cin, response);
    return response;
}

void open_dat_file(string prompt, ifstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}

void open_output_file(string prompt, ofstream& file_stream){
    string response;

    do{
        cout << prompt;
        getline(cin, response);
        file_stream.open(response, ofstream::out);
        if(! file_stream.is_open()){
            if(! yesno("The file could not be read. Try again?")) break;
        }
    }while(! file_stream.is_open());
}

//Perform a uniform displacement of all points in a set
void displace(vector<Point> &points, double xdisp, double ydisp){
    for(int iter = 0; iter < points.size(); iter ++){
        points[iter].x += xdisp;
        points[iter].y += ydisp;
    }
}

//Helper function to load point coordinates into C-style arrays
inline void assign_pt(Point p, double arr[2]){
    arr[0] = p.x;
    arr[1] = p.y;
}

//Determine whether a line segment with endpoints p1 and p2 intersects a line
//segment with endpoints p3 and p4. This method is based on an approach
//in Real-Time Collision Detection, by Christer Ericson, Morgan Kaufmann, 2004.
bool intersection(Point p1, Point p2, Point p3, Point p4){

    //Signed areas for performing crossing tests
    double a1, a2, a3, a4;
    //C-style arrays for holding points
    double cpt_1[2], cpt_2[2], cpt_3[2], cpt_4[2];

    assign_pt(p1, cpt_1);
    assign_pt(p2, cpt_2);
    assign_pt(p3, cpt_3);
    assign_pt(p4, cpt_4);

    a1 = orient2d(cpt_1, cpt_2, cpt_4);
    a2 = orient2d(cpt_1, cpt_2, cpt_3);
    if(a1 * a2 < 0){
        a3 = orient2d(cpt_3, cpt_4, cpt_1);
        a4 = a3 + a2 - a1;
        return a3*a4 < 0;
    }
    return false;
}

//Rotate each point in a set by a common angle about a common pivot
void rotate_points(vector<Point> &points, double angle, Point pivot){

    double sine = sin(angle), cosine = cos(angle);
    double x, y, newx, newy;

    for(int iter = 0; iter < points.size(); iter++){
        x = points[iter].x;
        y = points[iter].y;
        newx = cosine*x - sine*y + pivot.x*(1 - cosine) + pivot.y*sine;
        newy = cosine*y + sine*x + pivot.y*(1 - cosine) - pivot.x*sine;
        points[iter].x = newx;
        points[iter].y = newy;
    }
}

void rotate_point(Point &point, double angle, Point pivot){

    double sine = sin(angle), cosine = cos(angle);
    double x, y, newx, newy;

    x = point.x;
    y = point.y;
    newx = cosine*x - sine*y + pivot.x*(1 - cosine) + pivot.y*sine;
    newy = cosine*y + sine*x + pivot.y*(1 - cosine) - pivot.x*sine;
    point.x = newx;
    point.y = newy;
}

void get_extremes(vector<Point> points, double &minx, double &miny, double &maxx, double &maxy){

    minx = FLT_MAX;
    maxx = -FLT_MAX;
    miny = FLT_MAX;
    maxy = -FLT_MAX;

    for(Point p : points){
        if(p.x < minx) minx = p.x;
        if(p.x > maxx) maxx = p.x;
        if(p.y < miny) miny = p.y;
        if(p.y > maxy) maxy = p.y;
    }
}

double x_intersect(vector<Point> points, Edge e, double yval){
    double slope;
    Point p1 = points[e.idx1], p2 = points[e.idx2];

    slope = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y) / (p2.x - p1.x);
    return p1.x + (yval - p1.y) / slope;
}

double y_intersect(vector<Point> points, Edge e, double xval){
    double slope;
    Point p1 = points[e.idx1], p2 = points[e.idx2];

    slope = p1.x == p2.x ? BIG_SLOPE : (p2.y - p1.y) / (p2.x - p1.x);
    return p1.y + slope * (xval - p1.x);
}

double distance_sq(Point p1, Point p2){
    double xdiff, ydiff;
    xdiff = p1.x - p2.x;
    ydiff = p1.y - p2.y;
    return xdiff*xdiff + ydiff*ydiff;
}


//Utility function to create a random seed
unsigned get_random_seed(){
    
    unsigned seed;
    FILE *ranfile;
    int num_read = 0;
    time_t the_time;

    //Attempt to create a random seed by reading from /dev/urandom. If this
    //fails, use the current time since the beginning of the Unix epoch, in
    //seconds.    
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

void scale_vector(vector<double> &vec, double scale){
    for(int iter = 0; iter < vec.size(); iter ++) vec[iter] *= scale;
}

bool import_lattice(string name, vector<vector<double>>& rules, map<int,vector<vector<double>>>& nns, double scale){
    ifstream latfile;
    string nextline;
    bool again, blank, fileopen = false;
    vector<double> rule, nn;
    int lcount = 0, pointiter, nncount = 0, ruleiter = 0;

    //Attempt to open lattice file
    latfile.open(name);
        
    if(!latfile.is_open()){
        cerr << "The specified file could not be read read\n";
        return false;
    }   

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

/*
This function adds Gaussian random noise to the location of each point.
The function takes as arguments the original set of points and edges describing
the network, and a standard deviation for Gaussian random noise. Points are 
shifted according to a normal distribution.
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

template <typename T> 
void add_if_missing(const T &t, map<T, int> &idx_map, vector<T> &list){

    if(idx_map.find(t) == idx_map.end()){
        idx_map.insert(make_pair(t, idx_map.size()));
    }
}

void flush_with_edge(Point &pnt, double slope, double y_ext){
    double x_old, y_old, x_new;

    x_old = pnt.x;
    y_old = pnt.y;
    x_new = slope != 0 ? (y_ext - y_old) / slope + x_old : x_old;

    pnt.x = x_new;
    pnt.y = y_ext;
}

void add_thickness(NetworkComplex &current, double thickness, PolygonComplex &pc, bool makepoly, bool level){

    map<Point, vector<double>> angmap;
    map<Point, int> point_map, poly_pmap;
    Point p1, p2, key;
    Point p1f, p2f, p3f, p4f, p1fb, p2fb, p3fb, p4fb;
    double low, high, ang1, ang2, dx, dy, hwidth, slope;
    double ymin = FLT_MAX, ymax = -FLT_MAX, ylow, yhigh, y_ext;
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

            if((p1f.y < ylow - FLOAT_TOL || p1f.y > yhigh + FLOAT_TOL) &&level){
                y_ext = abs(p1.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p1f, slope, y_ext);
                add_if_missing(p1f, point_map, replace_points);
                p1fb = Point(p1.x, y_ext);
                add_if_missing(p1fb, point_map, replace_points);
                p1fflag = true;
            }

            else add_if_missing(p1f, point_map, replace_points);
            if((p3f.y < ylow - FLOAT_TOL || p3f.y > yhigh + FLOAT_TOL) &&level){
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

            if((p2f.y < ylow - FLOAT_TOL || p2f.y > yhigh + FLOAT_TOL) &&level){
                y_ext = abs(p2.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p2f, slope, y_ext);
                add_if_missing(p2f, point_map, replace_points);
                p2fb = Point(p2.x, y_ext);
                add_if_missing(p2fb, point_map, replace_points);
                p2fflag = true;
            }

            else add_if_missing(p2f, point_map, replace_points);
            if((p4f.y < ylow - FLOAT_TOL || p4f.y > yhigh + FLOAT_TOL) &&level){
                y_ext = abs(p2.y - ymin) < FLOAT_TOL ? ylow : yhigh;
                flush_with_edge(p4f, slope, y_ext);
                add_if_missing(p4f, point_map, replace_points);
                p4fb = Point(p2.x, y_ext);
                add_if_missing(p4fb, point_map, replace_points);
                p4fflag = true;
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
