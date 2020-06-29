#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct EXON {

    string chro;
    int start;
    int end;
    string gname;
};


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


int GPD_Length (GPD gpd) {

    int gpd_len = 0;
    for(int i = 0; i < gpd.nblock; ++i)  gpd_len += (gpd.vec_exon_end[i] - gpd.vec_exon_start[i] + 1);
    return gpd_len;
}


map<string,vector<EXON> > BuildExonMapFromFile (ifstream &inf) {

    map<string,vector<EXON> > exon_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            for(int i = 0; i < gpd.nblock; ++i){
                EXON exon_rc;
                exon_rc.chro = gpd.chro;
                exon_rc.start = gpd.vec_exon_start[i];
                exon_rc.end = gpd.vec_exon_end[i];
                exon_rc.gname = gpd.gname;
                exon_map[gpd.gname].push_back(exon_rc);
            }
        }
    }

    return exon_map;
}


int Overlap (EXON a, EXON b) {

    if(a.chro != b.chro)  return -1;
    else{
        if(a.end <= b.start || a.start >= b.end)  return -1;
        else{
            int real_start, real_end;
            if(a.start > b.start)  real_start = a.start;
            else  real_start = b.start;
            if(a.end > b.end)  real_end = b.end;
            else  real_end = a.end;
            return (real_end - real_start);
        }
    }
}


bool RealFusion (GPD gpd, string gene1, string gene2, map<string,vector<EXON> > &exon_map) {

    map<string,vector<EXON> >::iterator it;
    bool single1 = false;
    bool single2 = false;

    for(int i = 0; i < gpd.nblock; ++i){

        EXON exon_rc;
        exon_rc.chro = gpd.chro;
        exon_rc.start = gpd.vec_exon_start[i];
        exon_rc.end = gpd.vec_exon_end[i];
        exon_rc.gname = gpd.gname;

        //cout << exon_rc.start << '\t' << exon_rc.end << endl;

        bool overlap1 = false;
        it = exon_map.find(gene1);
        for(int j = 0; j < (it->second).size(); ++j){
            int overlap = Overlap(exon_rc, (it->second)[j]);
            if(overlap > 0){
                overlap1 = true;
                break;
            }
        }

        //if(overlap1)  cout << "overlap1" << endl;

        bool overlap2 = false;
        it = exon_map.find(gene2);
        for(int j = 0; j < (it->second).size(); ++j){
            int overlap = Overlap(exon_rc, (it->second)[j]);
            if(overlap > 0){
                overlap2 = true;
                break;
            }
        }

        //if(overlap2)  cout << "overlap2" << endl;

        if(overlap1 && (!overlap2))  single1 = true;
        if((!overlap1) && overlap2)  single2 = true;
    }

    return  single1 && single2;
}






int main (int argc, char **argv) {

    ifstream anno_inf(argv[1]);
    map<string,vector<EXON> > exon_map = BuildExonMapFromFile(anno_inf);
    anno_inf.close();

    ifstream inf(argv[2]);
    string strInput1, strInput2;

    while(inf){

        getline(inf, strInput1);
        if(strInput1.length() == 0)  break;  
        getline(inf, strInput2);

        GPD gpd;
        BuildGPDFromString(strInput1, gpd);
        vector<string> vec;
        split(vec, strInput2, is_any_of("\t"));
        string gene1 = vec[1];
        string gene2 = vec[3];

        bool real_fusion = RealFusion(gpd, gene1, gene2, exon_map);

        if(real_fusion)  cout << strInput2 << endl;
    }

    inf.close();

    return 0;
}
