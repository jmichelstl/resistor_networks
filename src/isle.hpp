#ifndef ISLE
#define ISLE

#include <memory>
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <vector>
#include <regex>
#include <tuple>
#include <map>

namespace isle {

    /*
     * This class is meant to parse files in which hierarhical information is
     * organized using nested sets of curly braces. A set of curly braces must
     * be preceded by an alphanumeric name string. The permissible contents of
     * a pair of curly braces are comma-delimited strings and child elements,
     * which themselves are indicated by an alphanumeric string and a pair of
     * curly braces. Comma-delimted attribute strings should have the form
     * of an alphanumeric attribute string, an equal sign padded by arbitrary
     * whitespace, and a value string, enclosed by quotation marks. This 
     * standard defines ISLE: The Input Serialization LanguagE
     */
    class Element {

        public:
            Element(){
                name = "";
                parent = nullptr;
            }

            Element(std::shared_ptr<Element> my_parent, std::string my_name) {
                parent = my_parent;
                name = my_name;
            }

            std::string get_name() {
                return name;
            }

            std::vector<std::string> get_attribute_names(){
                std::vector<std::string> names;
                for(auto it = attr_map.begin(); it != attr_map.end(); it++){
                    names.push_back(it->first);
                }
                return names;
            }

            void add_attribute(std::string name, std::string value){
                if(attr_map.find(name) == attr_map.end()){
                    attr_map.insert(make_pair(name,std::vector<std::string>()));
                }
                attr_map[name].push_back(value);
            }

            void add_child(std::shared_ptr<Element> child){
                if(cmap.find(child->get_name()) == cmap.end()){
                    cmap.insert(make_pair(child->get_name(),std::vector<std::shared_ptr<Element>>()));
                }
                cmap[child->get_name()].push_back(child);
            }

            std::shared_ptr<Element> get_parent(){
                return parent;
            }

            std::vector<std::string> get_child_names() {
                std::vector<std::string> names;
                for(auto it = cmap.begin(); it != cmap.end(); it++){
                    names.push_back(it->first);
                }
                return names;
            }

            bool query_attributes(std::string query, std::vector<std::string> &results){
                if(attr_map.find(query) == attr_map.end()) return false;
                results.assign(attr_map[query].begin(), attr_map[query].end());
                return true;
            }

            bool query_children(std::string query,std::vector<std::shared_ptr<Element>> &results){
                if(cmap.find(query) == cmap.end()) return false;
                results.assign(cmap[query].begin(), cmap[query].end());
                return true;
            }

            bool query_first_attribute(std::string query, std::string &result){
                if(attr_map.find(query) == attr_map.end()){
                    result.assign("");
                    return false;
                }
                result.assign(attr_map[query][0]);
                return true;
            }

            bool query_first_child(std::string query, std::shared_ptr<Element> &result){
                if(cmap.find(query) == cmap.end()) {
                    result = nullptr;
                    return false;
                }
                result = cmap[query][0];
                return true;
            }

        private:
            std::string name;
            std::map<std::string, std::vector<std::string>> attr_map;
            std::shared_ptr<Element> parent;
            std::map<std::string, std::vector<std::shared_ptr<Element>>> cmap;
    };

    //Read the contents of a file and combine all lines that are neither empty
    //nor start with a '#' sign into a single string.
    std::string file_to_string(std::ifstream &input){

        std::string file_string = "", nextline;

        //Regular expression to recognized lines beginning with optional
        //whitespace and a '#' sign
        std::regex comment("^\\s*#.*");

        while(! input.eof()){
            std::getline(input, nextline);
            //Igore blank lines and comments
            if(nextline.compare("")!= 0 && !std::regex_match(nextline.c_str(), comment)){
                file_string += nextline;
            }
        }

        return file_string;
    }

    //Evaluate a string meant to hold a name-value pair, in which the name is an
    //alphanumeric string padded by optional whitespace, and the value is a
    //string enclosed by quotation marks. Assignment is indicated by an equal
    //sign. Upon successful parsing of the string, a name and value are
    //initialized and true is returned. Upon failure, an error message is
    //printed and false is returned.
    bool parse_attribute(std::string argument, std::string &name, std::string &value ){

        std::cmatch match;

        std::regex_match(argument.c_str(), match, std::regex("^\\s*(\\w+)\\s*=\\s*\"(.*)\"\\s*$"));

        if(! match.size()){
            std::cerr << "The attribute string " << argument << " is misformatted.\n";
            return false;
        }

        name = match[1];
        value = match[2];

        return true;
    }

    //Attempt to read a file in which information about attributes and
    //sub-elements is organized hierarchically using nested sets of curly
    //braces.
    bool parse_isle(std::string file_name, std::shared_ptr<Element> &root){

        std::string file_contents, curr_string, attr_name, attr_value;
        std::stack<std::vector<std::string>> attr_stack;
        std::ifstream input;
        int pos;
        char c;
        std::shared_ptr<Element> curr_div = nullptr;

        //Recognize a string of alphanumeric characters padded on the left and
        //right by zero or more whitespace characters
        std::regex name_regex("\\s*(\\w+)\\s*");
        //Recognize a block of whitespace
        std::regex ws_regex("\\s+");

        input.open(file_name);
        if(! input.is_open()){
            std::cerr << "The specified file could not be opened.\n";
            return false;
        }

        file_contents = file_to_string(input);
        input.close();
        if(file_contents.compare("") == 0){
            std::cerr << "The specified file has no lines of input.\n";
            return false;
        }

        curr_string = "";
        attr_stack.push(std::vector<std::string>());
        for(pos = 0; pos < file_contents.size(); pos ++){

            c = file_contents[pos];

            //New elements are opened by left curly braces
            if(c == '{'){
                std::cmatch name_match;
                std::regex_match(curr_string.c_str(), name_match, name_regex);

                //If a legal name cannot be read, parsing has failed.
                if(! name_match.size()){
                    std::cerr << "An invalid element name was encountered.\n";
                    return false;
                }

                std::shared_ptr<Element> new_div = std::make_shared<Element>(curr_div, name_match[1]);
                if(curr_div != nullptr){
                    curr_div->add_child(new_div);
                }

                curr_div = new_div;

                curr_string = "";
                attr_stack.push(std::vector<std::string>());
            }

            //Elements are closed by right curly braces
            else if(c == '}'){

                //If the current Element is null, the input is misformatted
                if(curr_div == nullptr){
                    std::cerr << "'}' encountered without an opening '{'\n";
                    return false;
                }

                //Add any attributes to the current element
                if(curr_string.compare("") != 0 && ! std::regex_match(curr_string.c_str(), ws_regex)) {
                    attr_stack.top().push_back(curr_string);
                }
                for(std::string attribute : attr_stack.top()){
                    if(! parse_attribute(attribute, attr_name, attr_value)){
                        std::cerr << "Parsing failed due to a misformatted attribute.\n";
                        return false;
                    }
                    curr_div->add_attribute(attr_name, attr_value);
                }
                attr_stack.pop();

                //If the parent is null, parsing is finished. Any further input
                //is superfluous and should be disregarded.
                if(curr_div->get_parent() == nullptr){
                    root = curr_div;
                    return true;
                }

                curr_div = curr_div->get_parent();
                curr_string = "";
            }

            //Strings are delimited by commas
            else if(c == ','){

                //If the comma is escaped, do not terminate the string.
                if(curr_string.size() > 0 && *curr_string.rbegin() == '\\'){
                    curr_string.pop_back();
                    curr_string += ',';
                    continue;
                }

                //If comma-delimited strings exist outside of a block of curly
                //braces, the document is invalid.
                if(curr_div == nullptr){
                    std::cerr << "Invalid content was encountered outside of an element.\n";
                    return false;
                }

                attr_stack.top().push_back(curr_string);
                curr_string = "";
            }

            else {
                curr_string += c;
            }
        }

        //If execution reaches here, either no parent Element has been parsed,
        //or there is an unclosed element. Either way, the document is
        //misformatted and parsing fails.
        if(curr_div == nullptr){
            std::cerr << "No root-level element could be read.\n";
        }
        else {
            std::cerr << "There was an open element when the end-of-file was encountered.\n";
        }

        return false;
    }

    void print_element(std::shared_ptr<Element> e, std::string indent){

        std::vector<std::string> values, child_names;
        std::vector<std::shared_ptr<Element>> children;
        int iter;

        std::cout << indent << "name: " << e->get_name() << std::endl;
        std::cout << indent << "attributes:" << std::endl;
        for(std::string next_attr : e->get_attribute_names()){
            e->query_attributes(next_attr, values);
            std::cout << indent << "    " << next_attr << ": " << values[0];
            for(iter = 1; iter < values.size(); iter++){
                std::cout << ", " << values[iter];
            }
            std::cout << "\n";
        }

        child_names = e->get_child_names();
        if(child_names.size()){
            std::cout << indent << "children:" << std::endl;
            for(std::string child_name : child_names){
                e->query_children(child_name, children);
                std::cout << indent << "    " << child_name << std::endl;
                for(std::shared_ptr<Element> child : children){
                    print_element(child, indent + "        ");
                }
            }
        }
    }
}

#endif
