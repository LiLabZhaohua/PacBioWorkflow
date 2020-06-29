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


pair<int,int> AddLen (string &cigar) {

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

    return pair<int,int>(-left_S, aln_len + right_S);
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


int MatchPeak (vector<PEAK> &pk_vec, PEAK pk_query) {

    int result = -1;
    int locate = BinarySearch(pk_vec, pk_query);

    int dist_up = 10000;
    int dist_down = 10000;
    if(pk_query.chro == pk_vec[locate].chro)  dist_up = abs(pk_query.pos - pk_vec[locate].pos);
    if(locate < pk_vec.size() - 1 && pk_query.chro == pk_vec[locate + 1].chro)  dist_down = abs(pk_query.pos - pk_vec[locate + 1].pos);

    if(dist_up < dist_down && dist_up < 30 && pk_query.strand == pk_vec[locate].strand)  result = locate;
    else if(dist_down < 30 && pk_query.strand == pk_vec[locate + 1].strand)  result = locate + 1;
    else  ;

    return result;
}





int main (int argc, char **argv) {

    cerr << "SupportByCAGEOutputPeak <cage_peak.txt> <sam_file>" << endl;

    ifstream cageinf(argv[1]);
    ifstream saminf(argv[2]);

    vector<PEAK> cage_pkvec = ReadInPeakAndSort(cageinf);


    while(saminf){

        string strInput;
        getline(saminf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            int end_point;
            PEAK cage_query;
            cage_query.chro = vec[2];

            pair<int,int> add_len = AddLen(vec[5]);
            if(vec[1] == "0"){
                cage_query.strand = "+";
                cage_query.pos = lexical_cast<int>(vec[3]) + add_len.first;
                end_point = lexical_cast<int>(vec[3]) + add_len.second;
            }
            else if(vec[1] == "16"){
                cage_query.strand = "-";
                cage_query.pos = lexical_cast<int>(vec[3]) + add_len.second;
                end_point = lexical_cast<int>(vec[3]) + add_len.first;
            }
            else  ;

            int cage_match = MatchPeak(cage_pkvec, cage_query);

            if(cage_match >= 0){
                cout << vec[0] << '\t';
                cout << cage_pkvec[cage_match].chro << '\t' << cage_pkvec[cage_match].pos << '\t';
                cout << vec[2] << '\t' << end_point << endl;
            }
            else  cout << vec[0] << "\t*\t*\t" << vec[2] << '\t' << end_point << endl;
        }
    }


    cage_pkvec.clear();

    cageinf.close(); saminf.close();

    return 0;
}
