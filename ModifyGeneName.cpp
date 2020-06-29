#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


int main (int argc, char **argv) {

    ifstream linf(argv[1]);
    ifstream inf(argv[2]);

    map<string,string> gname_map;

    while(linf){

        string strInput;
        getline(linf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            gname_map[vec[0]] = vec[1];
        }
    }


    map<string,string>::iterator it;


    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            it = gname_map.find(vec[1]);
            if(it == gname_map.end())  cout << "Alert" << endl;
            else{
                cout << vec[0] << '\t' << it->second;
                for(int i = 2; i < vec.size(); ++i){
                    cout << '\t' << vec[i];
                }
                cout << endl;
            }
        }
    }


    linf.close(); inf.close();

    return 0;
}
