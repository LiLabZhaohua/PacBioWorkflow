#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


int main (int argc, char **argv) {

    cerr << "DeleteRedundancyInGPD <GPD file> <redundant_list.txt>" << endl;

    ifstream annoinf(argv[1]);
    ifstream rddinf(argv[2]);

    set<string> rdd_set;
    set<string>::iterator it;

    while(rddinf){

        string strInput;
        getline(rddinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            rdd_set.insert(vec[0]);
        }
    }


    while(annoinf){

        string strInput;
        getline(annoinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            it = rdd_set.find(vec[1]);
            if(it == rdd_set.end())  cout << strInput << endl;
        }
    }

    annoinf.close(); rddinf.close();

    return 0;
}
