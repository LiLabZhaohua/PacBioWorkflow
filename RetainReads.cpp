#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct ISOFORM {

    int tss;
    int polya;
    string tname;
    string strand;
};


map<string,vector<ISOFORM> > BuildIsoformMap (ifstream &inf) {

    map<string,vector<ISOFORM> > isoform_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            ISOFORM isoform;
            isoform.tname = vec[1];
            isoform.strand = vec[2];
            isoform.tss = lexical_cast<int>(vec[4]);
            isoform.polya = lexical_cast<int>(vec[6]);

            isoform_map[isoform.tname].push_back(isoform);
        }
    }

    return isoform_map;
}




int main (int argc, char **argv) {

    ifstream isoforminf(argv[1]);
    map<string,vector<ISOFORM> > isoform_map = BuildIsoformMap(isoforminf);
    map<string,vector<ISOFORM> >::iterator it;
    isoforminf.close();

    ifstream inf(argv[2]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            string tname = vec[2];
            string strand = vec[3];
            int tss = lexical_cast<int>(vec[5]);
            int polya = lexical_cast<int>(vec[7]);

            it = isoform_map.find(tname);

            if(it != isoform_map.end()){

                for(int i = 0; i < (it->second).size(); ++i){
                    if((it->second)[i].strand == strand && (it->second)[i].tss == tss && (it->second)[i].polya == polya){
                        cout << strInput << endl;
                        break;
                    }
                }
            }
        }
    }

    inf.close();

    return 0;
}
