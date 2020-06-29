#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
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
    //int cstart;
    //int cend;
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


void ChangeToBEDExon (GPD gpd) {

    if(gpd.strand == "+"){
        for(int i = 0; i < gpd.nblock; ++i){
            int left = gpd.vec_exon_start[i] - 200 * (i == 0);
            if(left < 1)  left = 1;
            int right = gpd.vec_exon_end[i];
            cout << gpd.chro << '\t' << left << '\t' << right << "\t.\t0\t" << gpd.strand << endl;
        }
    }
    else{
        for(int i = 0; i < gpd.nblock; ++i){
            int left = gpd.vec_exon_start[i];
            int right = gpd.vec_exon_end[i] + 200 * (i == gpd.nblock - 1);
            cout << gpd.chro << '\t' << left << '\t' << right << "\t.\t0\t" << gpd.strand << endl;
        }
    }
}




int main (int argc, char **argv) {

    ifstream inf(argv[1]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            ChangeToBEDExon(gpd);
        }
    }

    inf.close();

    return 0;
}
