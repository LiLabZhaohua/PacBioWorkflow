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

#define THRESHOLD 30

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


int compare_gpd_sort (const void *a, const void *b) {

    if( (*(GPD*) a).chro > (*(GPD*) b).chro )  return 1;
    else if ( (*(GPD*) a).chro == (*(GPD*) b).chro ){
        if( (*(GPD*) a).start > (*(GPD*) b).start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_two_gpd (GPD a, GPD b) {

    if(a.chro > b.chro)  return 1;
    else if (a.chro == b.chro){
        if(a.start > b.start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int BinarySearch (GPD gpd_query, GPD *gpd_array, int total) {

    int up = 0;
    int down = total - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_two_gpd(gpd_array[middle],gpd_query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}





int main (int argc, char **argv) {

    cerr << "FilterStringtieResult <annotation.gpd> <stringtie.gpd> <filtered_result.gtf>" << endl;

    ifstream inf(argv[1]);
    int line = 0;

    while(inf){
        string strInput;
        getline(inf, strInput);
        if(strInput.length() > 0)  ++line;
    }

    inf.close();
    inf.open(argv[1]);

    GPD *gpd_array = new GPD[line];

    int index = 0;

    while(inf){
        string strInput;
        getline(inf, strInput);
        if(strInput.length() > 0){
            BuildGPDFromString(strInput, gpd_array[index]);
            ++index;
        }
    }

    inf.close();


    qsort(&gpd_array[0],line,sizeof(GPD),compare_gpd_sort);


    ifstream sinf(argv[2]);
    ofstream ouf(argv[3], ios::trunc);
    vector<string> str_vec;


    while(sinf){

        string strInput;
        getline(sinf, strInput);

        if(strInput.length() > 0){

            GPD gpd_query;
            BuildGPDFromString(strInput, gpd_query);

            if(gpd_query.nblock == 1){

                int locate = BinarySearch(gpd_query, gpd_array, line);

                int down_locate = locate;
                bool overlap_down = false;
                while(down_locate < line - 1){
                    if(gpd_array[down_locate].chro != gpd_query.chro || gpd_array[down_locate].start > gpd_query.end)  break;
                    else{
                        for(int i = 0; i < gpd_array[down_locate].nblock; ++i){
                            bool overlap = false;
                            if(abs(gpd_query.start - gpd_array[down_locate].vec_exon_start[i]) < THRESHOLD && abs(gpd_query.end - gpd_array[down_locate].vec_exon_end[i]) < THRESHOLD)  overlap = true;
                            overlap_down = overlap_down | overlap;
                        }
                    }
                    ++down_locate;
                }

                int up_locate = locate;
                bool overlap_up = false;
                while(up_locate > 0){
                    --up_locate;
                    if(gpd_array[up_locate].chro != gpd_query.chro || locate - up_locate > 100 || gpd_query.start - gpd_array[up_locate].end > 50000)  break;
                    else{
                        for(int i = 0; i < gpd_array[up_locate].nblock; ++i){
                            bool overlap = false;
                            if(abs(gpd_query.start - gpd_array[up_locate].vec_exon_start[i]) < THRESHOLD && abs(gpd_query.end - gpd_array[up_locate].vec_exon_end[i]) < THRESHOLD)  overlap = true;
                            overlap_up = overlap_up | overlap;
                        }
                    }
                }

                if(overlap_up || overlap_down)  ouf << strInput << endl;
            }
            else  ouf << strInput << endl;
        }
    }

    sinf.close();
    ouf.close();

    return 0;
}
