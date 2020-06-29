#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cmath>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct END {

    string chro;
    string strand;
    int pos;
};


int compare_END_sort (const void *a, const void *b) {

    if( (*(END*) a).chro > (*(END*) b).chro )  return 1;
    else if ( (*(END*) a).chro == (*(END*) b).chro ){
        if( (*(END*) a).pos > (*(END*) b).pos )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_two_END (END a, END b) {

    if(a.chro > b.chro)  return 1;
    else if (a.chro == b.chro){
        if(a.pos > b.pos)  return 1;
        else  return -1;
    }
    else  return -1;
}


vector<END> BuildENDVector (ifstream &inf) {

    vector<END> end_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            END end;
            end.chro = vec[1];
            end.strand = vec[2];
            end.pos = lexical_cast<int>(vec[3]);

            end_vec.push_back(end);
        }
    }

    qsort(&end_vec[0], end_vec.size(), sizeof(END), compare_END_sort);

    return end_vec;
}


int BinarySearch (END query, vector<END> &target_vec) {

    int up = 0;
    int down = target_vec.size() - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_two_END(target_vec[middle],query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}





int main (int argc, char **argv) {

    ifstream annoinf(argv[1]);
    vector<END> end_vec = BuildENDVector(annoinf);
    annoinf.close();

    ifstream inf(argv[2]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            END query;
            query.chro = vec[1];
            query.strand = vec[2];
            query.pos = lexical_cast<int>(vec[3]);

            int locate = BinarySearch (query, end_vec);
            bool match = false;

            if(query.chro == end_vec[locate].chro && query.strand == end_vec[locate].strand && abs(query.pos - end_vec[locate].pos) < 30)  match = true;
            else if(locate < end_vec.size() - 1 && query.chro == end_vec[locate + 1].chro && query.strand == end_vec[locate + 1].strand && abs(query.pos - end_vec[locate + 1].pos) < 30)  match = true;
            else  ;

            if(!match)  cout << strInput << endl;
        }
    }

    annoinf.close(); inf.close();

    return 0;
}
