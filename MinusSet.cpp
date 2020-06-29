#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


int main (int argc, char **argv) {

    ifstream inf(argv[1]);
    ifstream sinf(argv[2]);
    ofstream ouf(argv[3], ios::trunc);

    set<string> tname_set;
    set<string>::iterator it;

    while(sinf){

        string strInput;
        getline(sinf, strInput);
        if(strInput.length() > 0){
            it = tname_set.find(strInput);
            if(it == tname_set.end())  tname_set.insert(strInput);
        }
    }

    cout << tname_set.size() << endl;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            //it = tname_set.find(vec[1]);
            it = tname_set.find(vec[0]);
            if(it == tname_set.end())  ouf << strInput << endl;
        }
    }

    inf.close(); sinf.close(); ouf.close();

    return 0;
}
