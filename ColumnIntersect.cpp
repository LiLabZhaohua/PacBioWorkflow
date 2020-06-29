#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


set<string> BuildColumnSet (ifstream &inf, int column) {

    set<string> column_set;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            column_set.insert(vec[column - 1]);
        }
    }

    return column_set;
}



int main (int argc, char **argv) {

    ifstream inf1(argv[1]);
    int column1 = atoi(argv[2]);
    ifstream inf2(argv[3]);
    int column2 = atoi(argv[4]);

    set<string> set1 = BuildColumnSet(inf1, column1);
    set<string> set2 = BuildColumnSet(inf2, column2);

    set<string>::iterator it;
    set<string>::iterator it_find;

    for(it = set2.begin(); it != set2.end(); ++it){
        it_find = set1.find(*it);
        if(it_find != set1.end())  cout << (*it) << endl;
    }

    set1.clear(); set2.clear();

    inf1.close(); inf2.close();

    return 0;
}
