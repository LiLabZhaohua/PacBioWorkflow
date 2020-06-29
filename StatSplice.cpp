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


struct SPLICE {

    string chro;
    int pos;
};


int compare_splice_sort (const void *a, const void *b) {

    if( (*(SPLICE*) a).chro > (*(SPLICE*) b).chro )  return 1;
    else if ( (*(SPLICE*) a).chro == (*(SPLICE*) b).chro ){
        if( (*(SPLICE*) a).pos >= (*(SPLICE*) b).pos )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_splice (SPLICE a, SPLICE b) {

    if(a.chro > b.chro)  return 1;
    else if(a.chro == b.chro){
        if(a.pos >= b.pos)  return 1;
        else  return -1;
    }
    else  return -1;
}


vector<SPLICE> ReadInPeakAndSort (ifstream &inf) {

    vector<SPLICE> sp_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            SPLICE sp;
            sp.chro = vec[0];
            sp.pos = lexical_cast<int>(vec[1]);
            sp_vec.push_back(sp);
        }
    }

    qsort(&sp_vec[0], sp_vec.size(), sizeof(SPLICE), compare_splice_sort);

    return sp_vec;
}


int BinarySearch (vector<SPLICE> &sp_vec, SPLICE sp_query) {

    int up = 0;
    int down = sp_vec.size() - 1;

    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_splice(sp_vec[middle],sp_query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}







int main (int argc, char **argv) {

    ifstream ainf(argv[1]);
    ifstream rinf(argv[2]);
    vector<SPLICE> asp_vec = ReadInPeakAndSort(ainf);
    vector<SPLICE> rsp_vec = ReadInPeakAndSort(rinf);

    for(int i = 0; i < rsp_vec.size(); ++i){

        int locate = BinarySearch(asp_vec, rsp_vec[i]);
        int dist_up = 10000;
        int dist_down = 10000;
        if(rsp_vec[i].chro == asp_vec[locate].chro)  dist_up = abs(rsp_vec[i].pos - asp_vec[locate].pos);
        if(locate < asp_vec.size() - 1 && rsp_vec[i].chro == asp_vec[locate + 1].chro)  dist_down = abs(rsp_vec[i].pos - asp_vec[locate + 1].pos);

        if(dist_up < dist_down)  cout << dist_up << endl;
        else  cout << dist_down << endl;
    }


    asp_vec.clear(); rsp_vec.clear();
    ainf.close(); rinf.close();

    return 0;
}
