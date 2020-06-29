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


struct ENDS {

    string strand;
    int tss;
    int polya;
};


map<string,vector<ENDS> > BuildEndsMapFromFile (ifstream &inf) {

    map<string,vector<ENDS> > ends_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            ENDS ends;
            ends.strand = vec[2];
            ends.tss = lexical_cast<int>(vec[4]);
            ends.polya = lexical_cast<int>(vec[6]);
            ends_map[vec[1]].push_back(ends);
        }
    }

    return ends_map;
}


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


void DisplayGPDInLine (GPD gpd) {

    cout << gpd.gname << '\t';
    cout << gpd.tname << '\t';
    cout << gpd.chro << '\t';
    cout << gpd.strand << '\t';
    cout << gpd.start << '\t';
    cout << gpd.end << '\t';
    cout << gpd.cstart << '\t';
    cout << gpd.cend << '\t';
    cout << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << ',';
    cout << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << ',';
    cout << endl;
}


map<string,GPD> BuildGPDMapFromFile (ifstream &inf) {

    map<string,GPD> output_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            output_map[gpd.tname] = gpd;
        }
    }

    return output_map;
}






int main (int argc, char **argv) {

    ifstream annoinf(argv[1]);
    ifstream endsinf(argv[2]);

    map<string,GPD> gpd_map = BuildGPDMapFromFile(annoinf);
    map<string,GPD>::iterator gpd_it;
    map<string,vector<ENDS> > ends_map = BuildEndsMapFromFile(endsinf);
    map<string,vector<ENDS> >::iterator ends_it;

    annoinf.close();
    endsinf.close();

    ifstream rdinf(argv[3]);

    while(rdinf){

        string strInput;
        getline(rdinf, strInput);

        if(strInput.length() > 0){
            //cout << "~~~~~~~~~~~~~~~~~~~~" << endl;
            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[4] == "*"){

                string tname = vec[2];
                string strand = vec[3];
                int polya = lexical_cast<int>(vec[7]);
                int tss = -1;
                bool ToReferAnno = true;

                ends_it = ends_map.find(tname);

                if(ends_it != ends_map.end()){
                    for(int i = 0; i < (ends_it->second).size(); ++i){
                        if((ends_it->second)[i].strand == strand && (ends_it->second)[i].polya == polya){
                            tss = (ends_it->second)[i].tss;
                            ToReferAnno = false;
                            cout << vec[0] << '\t' << vec[1] << '\t' << tname << '\t' << strand << '\t' << vec[6] << '\t' << tss << '\t' << vec[6] << '\t' << polya << "\tB" << endl;
                            break;
                        }
                    }
                }

                if(ToReferAnno){

                    gpd_it = gpd_map.find(tname);

                    if(gpd_it == gpd_map.end()){
                        cout << "Alert1" << endl;
                        exit(1);
                    }
                    else{
                        int nblock = (gpd_it->second).nblock;
                        if(strand == "+")  tss = (gpd_it->second).vec_exon_start[0];
                        else if(strand == "-")  tss = (gpd_it->second).vec_exon_end[nblock - 1];
                        else  ;
                        cout << vec[0] << '\t' << vec[1] << '\t' << tname << '\t' << strand << '\t' << vec[6] << '\t' << tss << '\t' << vec[6] << '\t' << polya << "\tC" << endl;
                    }
                }
            }
            else  cout << strInput << "\tA" << endl;
        }
    }

    rdinf.close();

    return 0;
}
