#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


string CutTail (string str) {

    int L = str.length();

    for(int i = str.length() - 1; i > 0; --i){
        if(str[i] == '_'){
            L = i;
            break;
        }
    }

    return str.substr(0,L);
}





int main (int argc, char **argv) {

    cerr << "CorrectRedundant <read_subisoform_fltr.txt> <redundant_list.txt>" << endl;

    ifstream subinf(argv[1]);
    ifstream rddinf(argv[2]);

    map<string,string> rdd_map;
    map<string,string>::iterator it;

    while(rddinf){

        string strInput;
        getline(rddinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            rdd_map[CutTail(vec[0])] = CutTail(vec[1]);
        }
    }

    //for(it = rdd_map.begin(); it != rdd_map.end(); ++it)  cout << (it->first) << '\t' << (it->second) << endl;

    while(subinf){

        string strInput;
        getline(subinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            it = rdd_map.find(vec[2]);
            if(it != rdd_map.end()){
                vec[2] = (it->second);
                cout << vec[0];
                for(int i = 1; i < vec.size(); ++i)  cout << '\t' << vec[i];
                cout << endl;
            }
            else  cout << strInput << endl;
        }
    }

    subinf.close(); rddinf.close();

    return 0;
}
