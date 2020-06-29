#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


map<string,vector<pair<int,int> > > BuildEndsMapFromFile (ifstream &inf) {

    map<string,vector<pair<int,int> > > ends_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[2] != "-1"){
                ends_map[vec[0]].push_back(pair<int,int>(lexical_cast<int>(vec[2]), lexical_cast<int>(vec[3])));
            }
        }
    }

    return ends_map;
}



int main (int argc, char **argv) {

    cerr << "IntactEndsIsoform <ends_isoform.txt> <intact-rd_ends.txt>" << endl;

    ifstream endsinf(argv[1]);
    map<string,vector<pair<int,int> > > ends_map = BuildEndsMapFromFile(endsinf);
    endsinf.close();

    map<string,vector<pair<int,int> > >::iterator it;
    ifstream inf(argv[2]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            it = ends_map.find(vec[1]);
            if(it != ends_map.end()){
                int tss = lexical_cast<int>(vec[3]);
                int polya = lexical_cast<int>(vec[4]);
                for(vector<pair<int,int> >::iterator it_vec = (it->second).begin(); it_vec != (it->second).end(); ++it_vec){
                    if(abs((*it_vec).first - tss) < 30 && abs((*it_vec).second - polya) < 30){
                        cout << strInput << '\t' << (*it_vec).first << '\t' << (*it_vec).second << endl;
                        break;
                    }
                }
            }
        }
    }

    inf.close();

    return 0;
}
