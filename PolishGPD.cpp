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
    //string cstart;
    //string cend;
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
    cout << "*" << '\t';
    cout << "*" << '\t';
    cout << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << ',';
    cout << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << ',';
    cout << endl;
}


map<string,GPD> BuildGPDMap (ifstream &inf) {

    map<string,GPD> gpd_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            gpd_map[gpd.tname] = gpd;
        }
    }

    return gpd_map;
}




int main (int argc, char **argv) {

    ifstream gpdinf(argv[1]);
    ifstream inf(argv[2]);

    map<string,GPD> gpd_map = BuildGPDMap(gpdinf);
    map<string,GPD>::iterator it;


    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            string read_name = vec[0];
            int tss = lexical_cast<int>(vec[5]);
            int polya = lexical_cast<int>(vec[7]);

            it = gpd_map.find(read_name);

            if(it == gpd_map.end()){
                cout << "Alert" << endl;
                continue;
            }

            if((it->second).strand == "+"){
                (it->second).vec_exon_start[0] = tss;
                (it->second).vec_exon_end[(it->second).nblock-1] = polya;
            }
            else{
                (it->second).vec_exon_start[0] = polya;
                (it->second).vec_exon_end[(it->second).nblock-1] = tss;
            }

            DisplayGPDInLine(it->second);
        }
    }

    gpdinf.close(); inf.close();

    return 0;
}
