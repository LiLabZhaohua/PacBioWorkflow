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


int main (int argc, char **argv) {

    ifstream inf(argv[1]);

    map<string,string> isoform_map;
    map<string,string>::iterator it;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            string isoform = "";
            for(int i = 2; i < 11; ++i)  isoform = isoform + vec[i] + '/';

            it = isoform_map.find(isoform);
            if(it != isoform_map.end()){
                if((it->second).substr(0,4) == "CAGE" && vec[1].substr(0,4) != "CAGE")  cout << (it->second) << '\t' << vec[1] << endl;
                else  cout << vec[1] << '\t' << (it->second) << endl;
            }
            else  isoform_map[isoform] = vec[1];
        }
    }

    inf.close();

    return 0;
}
