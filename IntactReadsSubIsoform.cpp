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


struct GPD {

    string gname;
    string tname;
    string chro;
    string strand;
    int start;
    int end;
    string cstart;
    string cend;
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
    //gpd.cstart = lexical_cast<int>(vec_temp[6]);
    //gpd.cend = lexical_cast<int>(vec_temp[7]);
    gpd.nblock = lexical_cast<int>(vec_temp[8]);
    vector<string> vec1;
    boost::split(vec1, vec_temp[9], boost::is_any_of(","));
    for(int i = 0; i < vec1.size() - 1; ++i)  (gpd.vec_exon_start).push_back(lexical_cast<int>(vec1[i]));
    vector<string> vec2;
    boost::split(vec2, vec_temp[10], boost::is_any_of(","));
    for(int i = 0; i < vec2.size() - 1; ++i)  (gpd.vec_exon_end).push_back(lexical_cast<int>(vec2[i]));
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


string TrimTranscriptName (string str) {

    string str_trimmed;
    for(int i = str.length() - 1; i > 0; --i){
        if(str[i] == '_'){
            str_trimmed = str.substr(0,i);
            break;
        }
    }

    return str_trimmed;
}


map<string,vector<GPD> > BuildSubIsoformMap (ifstream &inf) {

    map<string,vector<GPD> > isoform_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            string trimmed_name = TrimTranscriptName(gpd.tname);
            isoform_map[trimmed_name].push_back(gpd);
        }
    }

    return isoform_map;
}




int main (int argc, char **argv) {

    ifstream annoinf(argv[1]);
    map<string,vector<GPD> > isoform_map = BuildSubIsoformMap(annoinf);
    annoinf.close();

    ifstream inf(argv[2]);
    map<string,vector<GPD> >::iterator it;

    //for(it = isoform_map.begin(); it != isoform_map.end(); ++it)  cout << (it->first) << endl;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[8] == "A"){

                string strand = vec[3];
                string chro = vec[4];
                int tss_pos = lexical_cast<int>(vec[5]);
                int polya_pos = lexical_cast<int>(vec[7]);

                it = isoform_map.find(vec[2]);

                if(it == isoform_map.end()){
                    cout << "Alert" << endl;
                    //cout << strInput << endl;
                }
                else{
                    int matched = 0;
                    for(int i = 0; i < (it->second).size(); ++i){
                        if(strand == (it->second)[i].strand && chro == (it->second)[i].chro){
                            if(strand == "+"){
                                if(tss_pos == (it->second)[i].vec_exon_start[0] && polya_pos == (it->second)[i].vec_exon_end[(it->second)[i].nblock - 1]){
                                    cout << vec[0] << '\t' << vec[1] << '\t' << (it->second)[i].tname << endl;
                                    ++matched;
                                }
                            }
                            else if(strand == "-"){
                                if(tss_pos == (it->second)[i].vec_exon_end[(it->second)[i].nblock - 1] && polya_pos == (it->second)[i].vec_exon_start[0]){
                                    cout << vec[0] << '\t' << vec[1] << '\t' << (it->second)[i].tname << endl;
                                    ++matched;
                                }
                            }
                            else  cout << "Alert" << endl;
                        }
                    }
                }
            }
        }
    }

    annoinf.close(); inf.close();

    return 0;
}
