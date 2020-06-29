#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


int main (int argc, char **argv) {

    ifstream inf(argv[1]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("_"));

            for(int i = 0; i < vec.size() - 1; ++i){
                cout << vec[i] << '_';
            }
            cout << endl;
        }
    }

    inf.close();

    return 0;
}
