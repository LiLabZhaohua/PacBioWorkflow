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
    map<string,int>  gene_abundance;
    map<string,int>::iterator it;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            string gname = vec[0];
            int abundance = lexical_cast<int>(vec[7]);

            it = gene_abundance.find(gname);
            if(it == gene_abundance.end())  gene_abundance[gname] = abundance;
            else  (it->second) += abundance;
        }
    }

    inf.close();
    inf.open(argv[1]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            string gname = vec[0];
            int abundance = lexical_cast<int>(vec[7]);

            it = gene_abundance.find(gname);
            if(it == gene_abundance.end())  cout << "Alert" << endl;
            else{
                double proportion = abundance / ((it->second) + 0.0);
                if(proportion > 0.05 || abundance >= 5)  cout << strInput << endl;
            }
        }
    }

    inf.close();

    return 0;
}
