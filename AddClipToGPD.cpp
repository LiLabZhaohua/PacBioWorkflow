#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
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


void DisplayGPDInLine (GPD gpd) {

    cout << gpd.gname << '\t';
    cout << gpd.tname << '\t';
    cout << gpd.chro << '\t';
    cout << gpd.strand << '\t';
    cout << gpd.start << '\t';
    cout << gpd.end << '\t';
    cout << "*\t";
    cout << "*\t";
    cout << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << ',';
    cout << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << ',';
    cout << endl;
}


pair<int,int> ClipLen (string &cigar) {

    char ch[5000];
    int index = 0;
    int aln_len = 0;
    int S_num = 0;
    int left_S = 0;
    int right_S = 0;

    for(int i = 0; i < cigar.length(); ++i){

        if(cigar[i] > 'A' && cigar[i] < 'Z'){

            for(int j = index; j < 5000; ++j)  ch[j] = '\0';
            int number = atoi(ch);
            index = 0;

            if(cigar[i] == 'M' || cigar[i] == 'D' || cigar[i] == 'N')  aln_len += number;
            else if(cigar[i] == 'S' || cigar[i] == 'H'){
                ++S_num;
                if(S_num == 1)  left_S = number;
                else if(S_num == 2)  right_S = number;
                else  ;
            }
        }
        else{
            ch[index] = cigar[i];
            ++index;
        }
    }

    return pair<int,int>(-left_S, right_S);
}





int main (int argc, char **argv) {

    ifstream gpdinf(argv[1]);
    ifstream saminf(argv[2]);

    map<string,GPD> gpd_map;
    map<string,GPD>::iterator it;

    while(gpdinf){

        string strInput;
        getline(gpdinf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        string name = vec[0];
        GPD gpd;
        BuildGPDFromString(strInput,gpd);

        gpd_map[name] = gpd;
    }


    while(saminf){

        string strInput;
        getline(saminf, strInput);
        if(strInput.length() == 0)  break;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        string name = vec[0];

        it = gpd_map.find(name);
        if(it == gpd_map.end())  cout << "Alert" << endl;

        pair<int,int> clip_len = ClipLen(vec[5]);
        (it->second).vec_exon_start[0] = (it->second).vec_exon_start[0] + clip_len.first;
        if((it->second).vec_exon_start[0] <= 0)  (it->second).vec_exon_start[0] = 1;
        (it->second).vec_exon_end[(it->second).nblock-1] = (it->second).vec_exon_end[(it->second).nblock-1] + clip_len.second;

        (it->second).start = (it->second).vec_exon_start[0];
        (it->second).end = (it->second).vec_exon_end[(it->second).nblock-1];
    }


    for(it = gpd_map.begin(); it != gpd_map.end(); ++it){
        DisplayGPDInLine(it->second);
    }

    gpdinf.close(); saminf.close();

    return 0;
}
