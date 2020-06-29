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


struct PEAK {

    string chro;
    string strand;
    int pos;
};


int compare_peak_sort (const void *a, const void *b) {

    if( (*(PEAK*) a).chro > (*(PEAK*) b).chro )  return 1;
    else if ( (*(PEAK*) a).chro == (*(PEAK*) b).chro ){
        if( (*(PEAK*) a).pos >= (*(PEAK*) b).pos )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_peak (PEAK a, PEAK b) {

    if(a.chro > b.chro)  return 1;
    else if(a.chro == b.chro){
        if(a.pos >= b.pos)  return 1;
        else  return -1;
    }
    else  return -1;
}


vector<PEAK> ReadInPeakAndSort (ifstream &inf) {

    vector<PEAK> pk_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            PEAK pk;
            pk.chro = vec[0];
            pk.strand = vec[1];
            pk.pos = lexical_cast<int>(vec[2]);
            pk_vec.push_back(pk);
        }
    }

    //qsort(&pk_vec[0], pk_vec.size(), sizeof(PEAK), compare_peak_sort);

    return pk_vec;
}


int BinarySearch (vector<PEAK> &pk_vec, PEAK pk_query) {

    int up = 0;
    int down = pk_vec.size() - 1;

    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_peak(pk_vec[middle],pk_query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}





int main (int argc, char **argv) {

    cerr << "FilterTranscriptByPAS <peaks.txt> <assembly.gpd>" << endl;

    ifstream pinf(argv[1]);
    ifstream annoinf(argv[2]);
    vector<PEAK> pk_vec = ReadInPeakAndSort(pinf);

    while(annoinf){

        string strInput;
        getline(annoinf, strInput);
        if(strInput.length() == 0)  break;

        GPD gpd;
        BuildGPDFromString(strInput, gpd);

        PEAK three_prime_end;
        three_prime_end.chro = gpd.chro;
        three_prime_end.strand = gpd.strand;
        if(gpd.strand == "+")  three_prime_end.pos = gpd.vec_exon_end[gpd.nblock - 1];
        else  three_prime_end.pos = gpd.vec_exon_start[0];

        int locate = BinarySearch(pk_vec, three_prime_end);

        int dist_up = 10000;
        int dist_down = 10000;
        if(three_prime_end.chro == pk_vec[locate].chro && three_prime_end.strand == pk_vec[locate].strand)  dist_up = abs(three_prime_end.pos - pk_vec[locate].pos);
        if(locate < pk_vec.size() - 1 && three_prime_end.chro == pk_vec[locate + 1].chro && three_prime_end.strand == pk_vec[locate + 1].strand)  dist_down = abs(three_prime_end.pos - pk_vec[locate + 1].pos);

        if(dist_up < 30 || dist_down < 30)  cout << strInput << endl;
    }


    pinf.close(); annoinf.close();

    return 0;
}
