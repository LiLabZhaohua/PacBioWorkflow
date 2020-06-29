#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>
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

    ifstream gpdinf(argv[1]);
    ifstream endsinf(argv[2]);

    map<string,GPD> gpd_map = BuildGPDMapFromFile(gpdinf);
    map<string,GPD>::iterator it;
    map<string,int> count_map;
    for(it = gpd_map.begin(); it != gpd_map.end(); ++it)  count_map[(it->second).tname] = 0;
    map<string,int>::iterator it_count;

    while(endsinf){

        string strInput;
        getline(endsinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            string tname = vec[0];
            int cage_pos = lexical_cast<int>(vec[2]);
            int pas_pos = lexical_cast<int>(vec[3]);

            it = gpd_map.find(tname);
            it_count = count_map.find(tname);
            if(it == gpd_map.end() || it_count == count_map.end()){
                cout << "Alert" << endl;
                exit(1);
            }

            GPD corr_gpd = (it->second);
            int num = (it_count->second) + 1;
            ++(it_count->second);

            corr_gpd.tname = corr_gpd.tname + "_" + lexical_cast<string>(num);
            corr_gpd.strand = vec[1];
            int block_num = corr_gpd.nblock;

            if(corr_gpd.strand == "+"){
                if(cage_pos == -1)  corr_gpd.tname = corr_gpd.tname + "_u";
                else{
                    if(cage_pos < corr_gpd.vec_exon_end[0])  corr_gpd.vec_exon_start[0] = cage_pos;
                }
                corr_gpd.vec_exon_end[block_num - 1] = pas_pos;
            }
            else{
                if(cage_pos == -1)  corr_gpd.tname = corr_gpd.tname + "_u";
                else{
                    if(cage_pos > corr_gpd.vec_exon_start[block_num - 1])  corr_gpd.vec_exon_end[block_num - 1] = cage_pos;
                }
                corr_gpd.vec_exon_start[0] = pas_pos;
            }

            //if(corr_gpd.vec_exon_start[0] > corr_gpd.vec_exon_end[block_num - 1])  cout << "Alert" << endl;

            DisplayGPDInLine(corr_gpd);
        }
    }

    gpdinf.close(); endsinf.close();

    return 0;
}
