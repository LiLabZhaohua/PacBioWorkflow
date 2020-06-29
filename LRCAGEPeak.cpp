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


struct PEAK {

    string chro;
    //string orien;
    //string name;
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
            //pk.orien = vec[1];
            //pk.name = vec[2];
            //pk.pos = lexical_cast<int>(vec[3]);
            pk.pos = lexical_cast<int>(vec[1]);
            pk_vec.push_back(pk);
        }
    }

    qsort(&pk_vec[0], pk_vec.size(), sizeof(PEAK), compare_peak_sort);

    return pk_vec;
}


int AddLen (string &cigar, char orien) {

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

    int add_len = 0;
    if(orien == '+')  add_len = -left_S;
    else if(orien == '-')  add_len = aln_len + right_S;
    else  ;

    return add_len;
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

    cerr << "LRCAGEPeak <peaks.txt> <alignment.sam>" << endl;
    cerr << "The <peaks.txt> file contains two columns: chro & pos" << endl;

    ifstream pinf(argv[1]);
    ifstream sinf(argv[2]);
    vector<PEAK> pk_vec = ReadInPeakAndSort(pinf);
    int count[pk_vec.size()];
    for(int i = 0; i < pk_vec.size(); ++i)  count[i] = 0;

    while(sinf){

        string strInput;
        getline(sinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            int end_point = 0;
            if(vec[1] == "0"){
                end_point = lexical_cast<int>(vec[3]) + AddLen(vec[5], '+');
            }
            else if(vec[1] == "16"){
                end_point = lexical_cast<int>(vec[3]) + AddLen(vec[5], '-');
            }
            else  ;

            //cout << end_point << endl;

            PEAK pk_query;
            pk_query.chro = vec[2];
            pk_query.pos = end_point;

            int locate = BinarySearch(pk_vec, pk_query);

            int dist_up = 10000;
            int dist_down = 10000;
            if(pk_query.chro == pk_vec[locate].chro)  dist_up = abs(pk_query.pos - pk_vec[locate].pos);
            if(locate < pk_vec.size() - 1 && pk_query.chro == pk_vec[locate + 1].chro)  dist_down = abs(pk_query.pos - pk_vec[locate + 1].pos);

            if(dist_up < 30 && dist_up < dist_down)  ++count[locate];
            else if(dist_down < 30 && dist_down < dist_up)  ++count[locate + 1];
            else  ;
        }
    }


    for(int i = 0; i < pk_vec.size(); ++i){
        if(count[i] > 0)  cout << pk_vec[i].chro << '\t' << pk_vec[i].pos << '\t' << count[i] << endl;
    }


    pk_vec.clear();

    pinf.close(); sinf.close();

    return 0;
}
