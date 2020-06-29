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


class ENDS {

public:

    int tss;
    int polya;
    string gname;
    string tname;
    string rname;
    string strand;
    string chro_tss;
    string chro_polya;

    void display () {
        //if(tss == -1)  cout << rname << '\t' << gname << '\t' << tname << '\t' << strand << '\t' << chro_tss << "\t*\t" << chro_polya << '\t' << polya << endl;
        if(tss == -1)  cout << rname << '\t' << gname << '\t' << tname << '\t' << strand << "\t*\t*\t" << chro_polya << '\t' << polya << endl;
        else  cout << rname << '\t' << gname << '\t' << tname << '\t' << strand << '\t' << chro_tss << '\t' << tss << '\t' << chro_polya << '\t' << polya << endl;
    }
};


int compare_ends_sort (const void *a, const void *b) {

    if( (*(ENDS*) a).tss >= (*(ENDS*) b).tss )  return 1;
    else  return -1;
}


map<string,vector<ENDS> > BuildEndsMap (ifstream &inf) {

    map<string,vector<ENDS> > ends_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            ENDS ends;
            ends.rname = vec[0];
            ends.gname = vec[1];
            ends.tname = vec[2];
            ends.strand = vec[3];
            ends.chro_tss = vec[4];
            ends.chro_polya = vec[6];
            if(vec[5] == "*")  ends.tss = -1;
            else  ends.tss = lexical_cast<int>(vec[5]);
            ends.polya = lexical_cast<int>(vec[7]);

            ends_map[ends.gname].push_back(ends);
        }
    }

    map<string,vector<ENDS> >::iterator it;

    for(it = ends_map.begin(); it != ends_map.end(); ++it){
        qsort(&((it->second)[0]), (it->second).size(), sizeof(ENDS), compare_ends_sort);
    }

    return ends_map;
}


struct FREQ {

    int key;
    int freq;
};


vector<FREQ> BuildFreqVector (vector<ENDS> &ends_vec) {

    vector<FREQ> freq_vec;
    int pre_value = -100;

    for(int i = 0; i < ends_vec.size(); ++i){
        if(ends_vec[i].tss == pre_value)  ++freq_vec[freq_vec.size() - 1].freq;
        else{
            FREQ freq_rc;
            freq_rc.key = ends_vec[i].tss;
            freq_rc.freq = 1;
            freq_vec.push_back(freq_rc);
            pre_value = freq_rc.key;
        }
    }

    return freq_vec;
}


int FindMaxFreq (vector<FREQ> &freq_vec) {

    int max_freq = -100;
    int max_pos = -1;

    for(int i = 0; i < freq_vec.size(); ++i){
        if(freq_vec[i].freq > max_freq){
            max_freq = freq_vec[i].freq;
            max_pos = freq_vec[i].key;
        }
    }

    return max_pos;
}


void PolishEnds (vector<ENDS> &ends_vec, vector<FREQ> &freq_vec) {

    int traversed_num = 0;

    while(traversed_num < freq_vec.size()){

        int max_key = FindMaxFreq(freq_vec);
        //cout << "max_key = " << max_key << endl;
        for(int i = 0; i < freq_vec.size(); ++i){
            if(abs(freq_vec[i].key - max_key) < 100 && freq_vec[i].freq > 0){
                freq_vec[i].freq = -1;
                ++traversed_num;
            }
        }

        //cout << "traversed_num = " << traversed_num << endl;

        for(int i = 0; i < ends_vec.size(); ++i){
            if(abs(ends_vec[i].tss - max_key) < 100)  ends_vec[i].tss = max_key;
        }
    }
}





int main (int argc, char **argv) {

    ifstream inf(argv[1]);
    map<string,vector<ENDS> > ends_map = BuildEndsMap(inf);
    inf.close();

    map<string,vector<ENDS> >::iterator it;

    for(it = ends_map.begin(); it != ends_map.end(); ++it){

        //cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

        /*for(int i = 0; i < (it->second).size(); ++i){
            cout << (it->second)[i].rname << '\t' << (it->second)[i].tss << endl;
        }*/

        vector<FREQ> freq_vec = BuildFreqVector(it->second);

        /*for(int i = 0; i < freq_vec.size(); ++i){
            cout << freq_vec[i].key << '\t' << freq_vec[i].freq << endl;
        }*/

        PolishEnds(it->second, freq_vec);

        for(int i = 0; i < (it->second).size(); ++i){
            (it->second)[i].display();
        }
    }

    return 0;
}
