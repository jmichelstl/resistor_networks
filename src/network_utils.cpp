#include "network_utils.hpp"

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
    maxx = FLT_MIN;
    miny = FLT_MAX;
    maxy = FLT_MIN;

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
