#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


using namespace std;
using namespace boost;


struct ENDS {

    string strand;
    int CAGE_pos;
    int PAS_pos;
};


bool SimilarEnds (ENDS a, ENDS b) {

    if(a.strand == b.strand && abs(a.CAGE_pos - b.CAGE_pos) < 30 && abs(a.PAS_pos - b.PAS_pos) < 30)  return true;
    else  return false;
}


vector<ENDS> ReduceEndsVectorII (vector<ENDS> &vec) {

    vector<ENDS> output_vec;
    output_vec.push_back(vec[0]);

    for(int i = 1; i < vec.size(); ++i){
        bool SIMILAR = false;
        for(int j = 0; j < output_vec.size(); ++j){
            SIMILAR = SimilarEnds(vec[i], output_vec[j]);
            if(SIMILAR)  break;
        }
        if(!SIMILAR)  output_vec.push_back(vec[i]);
    }

    return output_vec;
}


vector<ENDS> DealWithMinusOne (vector<ENDS> &vec) {

    vector<ENDS> output_vec;

    for(int i = 0; i < vec.size(); ++i){
        if(vec[i].CAGE_pos == -1){
            bool compatible = false;
            for(int j = 0; j < vec.size(); ++j){
                if((vec[j].strand == vec[i].strand) && (vec[j].CAGE_pos != -1) && (abs(vec[j].PAS_pos - vec[i].PAS_pos) < 30)){
                    compatible = true;
                    break;
                }
            }
            if(!compatible)  output_vec.push_back(vec[i]);
        }
        else  output_vec.push_back(vec[i]);
    }

    return output_vec;
}




int main (int argc, char **argv) {

    cerr << "EndsIsoform <isoform_CAGE_PAS_ends.txt>" << endl;

    ifstream endsinf(argv[1]);
    map<string,vector<ENDS> > isoform_ends_map;
    set<string> undetermined_set;
    map<string,vector<ENDS> >::iterator it_map;
    set<string>::iterator it_set;

    while(endsinf){

        string strInput;
        getline(endsinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[2] != "*"){
                ENDS ends;
                ends.strand = vec[1];
                ends.CAGE_pos = lexical_cast<int>(vec[3]);
                ends.PAS_pos = lexical_cast<int>(vec[5]);
                isoform_ends_map[vec[0]].push_back(ends);
            }
            else{
                ENDS ends;
                ends.strand = vec[1];
                ends.CAGE_pos = -1;
                ends.PAS_pos = lexical_cast<int>(vec[5]);
                isoform_ends_map[vec[0]].push_back(ends);
            }
        }
    }



    for(it_map = isoform_ends_map.begin(); it_map != isoform_ends_map.end(); ++it_map){
        vector<ENDS> reduced_vec = ReduceEndsVectorII(it_map->second);
        vector<ENDS> final_vec = DealWithMinusOne(reduced_vec);
        for(int i = 0; i < final_vec.size(); ++i){
            if(final_vec.size() == 0)  cout << "Alert" << endl;
            cout << (it_map->first) << '\t' << final_vec[i].strand << '\t' << final_vec[i].CAGE_pos << '\t' << final_vec[i].PAS_pos << endl;
        }
    }

            
    endsinf.close();

    return 0;
}
