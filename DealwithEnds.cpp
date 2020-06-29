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


struct GPD {

    string gname;
    string tname;
    string chro;
    string strand;
    int start;
    int end;
    int cstart;
    int cend;
    int nblock;
    vector<int> vec_exon_start;
    vector<int> vec_exon_end;
};


void BuildGPDFromString (string &str, GPD &gpd) {

    vector<string> vec_temp;
    boost::split(vec_temp, str, boost::is_any_of("\t"));
    gpd.gname = vec_temp[0];
    gpd.tname = vec_temp[1];
    int L = vec_temp[2].length();
    gpd.chro = vec_temp[2];
    gpd.strand = vec_temp[3];
    gpd.start = lexical_cast<int>(vec_temp[4]);
    gpd.end = lexical_cast<int>(vec_temp[5]);
    gpd.cstart = lexical_cast<int>(vec_temp[6]);
    gpd.cend = lexical_cast<int>(vec_temp[7]);
    gpd.nblock = lexical_cast<int>(vec_temp[8]);
    vector<string> vec1;
    boost::split(vec1, vec_temp[9], boost::is_any_of(","));
    for(int i = 0; i < vec1.size() - 1; ++i)  (gpd.vec_exon_start).push_back(lexical_cast<int>(vec1[i]));
    vector<string> vec2;
    boost::split(vec2, vec_temp[10], boost::is_any_of(","));
    for(int i = 0; i < vec2.size() - 1; ++i)  (gpd.vec_exon_end).push_back(lexical_cast<int>(vec2[i]));
}


void DisplayGPD (GPD gpd) {

    cout << gpd.gname << endl;
    cout << gpd.tname << endl;
    cout << gpd.chro << endl;
    cout << gpd.strand << endl;
    cout << gpd.start << endl;
    cout << gpd.end << endl;
    cout << gpd.cstart << endl;
    cout << gpd.cend << endl;
    cout << gpd.nblock << endl;
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << '\t';
    cout << endl;
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << '\t';
    cout << endl;
}


map<string,GPD> BuildIsoformMapFromGPDFile (ifstream &inf) {

    map<string,GPD> isoform_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            isoform_map.insert(pair<string,GPD>(gpd.tname, gpd));
        }
    }

    return isoform_map;
}


struct ENDS {

    int CAGE_pos;
    int PAS_pos;
};


bool SimilarEnds (ENDS a, ENDS b) {

    if(abs(a.CAGE_pos - b.CAGE_pos) < 10 && abs(a.PAS_pos - b.PAS_pos) < 10)  return true;
    else  return false;
}


void ReduceEndsVector (vector<ENDS> &vec) {

    if(vec.size() > 1){

        int last = 0;
        vector<int> to_delete;

        for(int i = 1; i < vec.size(); ++i){
            bool check = SimilarEnds(vec[last], vec[i]);
            if(check)  to_delete.push_back(i);
            else  last = i;
        }

        if(to_delete.size() > 0){
            cout << "To delete" << endl;
            cout << "=============" << endl;
            for(int k = 0; k < vec.size(); ++k)  cout << vec[k].CAGE_pos << '\t' << vec[k].PAS_pos << endl;
        }

        for(int j = 0; j < to_delete.size(); ++j)  vec.erase(vec.begin() + to_delete[j] - j);

        if(to_delete.size() > 0){
            cout << "=============" << endl;
            for(int k = 0; k < vec.size(); ++k)  cout << vec[k].CAGE_pos << '\t' << vec[k].PAS_pos << endl;
        }
    }
}





int main (int argc, char **argv) {

    cerr << "DealwithEnds <assembly.gpd> <isoform_CAGE_PAS_ends.txt>" << endl;

    ifstream asseminf(argv[1]);
    map<string,GPD> isoform_map = BuildIsoformMapFromGPDFile(asseminf);
    map<string,GPD>::iterator it_isoform;
    asseminf.close();

    map<string,vector<ENDS> > isoform_ends_map;
    ifstream endsinf(argv[2]);

    while(endsinf){

        string strInput;
        getline(endsinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            ENDS ends;

            if(vec[1] == "*"){
                it_isoform = isoform_map.find(vec[0]);
                if((it_isoform->second).strand == "+"){
                    ends.CAGE_pos = (it_isoform->second).vec_exon_start[0];
                    ends.PAS_pos = (it_isoform->second).vec_exon_end[(it_isoform->second).nblock - 1];
                }
                else{
                    ends.CAGE_pos = (it_isoform->second).vec_exon_end[(it_isoform->second).nblock - 1];
                    ends.PAS_pos = (it_isoform->second).vec_exon_start[0];
                }
            }
            else{
                ends.CAGE_pos = lexical_cast<int>(vec[2]);
                ends.PAS_pos = lexical_cast<int>(vec[4]);
            }

            isoform_ends_map[vec[0]].push_back(ends);
        }
    }

    map<string,vector<ENDS> >::iterator it;

    for(it = isoform_ends_map.begin(); it != isoform_ends_map.end(); ++it){
        cout << "~~~~~~~~~~~~~~~" << endl;
        cout << (it->first) << endl;
        ReduceEndsVector(it->second);
    }

            
    endsinf.close();

    return 0;
}
