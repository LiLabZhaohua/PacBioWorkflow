#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>
#include <cstdlib>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;



int main (int argc, char **argv) {

    cerr << "SearchLine <STRING, key_file> <INT, key_position> <STRING, search_file> <INT, search_position>" << endl;

    ifstream kinf(argv[1]);
    int key_pos = atoi(argv[2]) - 1;
    ifstream inf(argv[3]);
    int search_pos = atoi(argv[4]) - 1;


    set<string> key_set;
    set<string>::iterator it;


    while(kinf){

        string strInput;
        getline(kinf, strInput);

        if(strInput.length() > 0){
            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            key_set.insert(vec[key_pos]);
        }
    }


    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            it = key_set.find(vec[search_pos]);
            if(it == key_set.end())  cout << strInput << endl;
        }
    }


    kinf.close(); inf.close();

    return 0;
}
