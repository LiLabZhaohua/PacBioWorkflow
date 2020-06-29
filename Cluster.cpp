#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#define THRESHOLD 10

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


void DisplayGPDinLine (GPD gpd) {

    cout << gpd.gname << '\t' << gpd.tname << '\t' << gpd.chro << '\t' << gpd.strand << '\t' << gpd.start << '\t' << gpd.end << '\t' << gpd.cstart << '\t';
    cout << gpd.cend << '\t' << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << ',';
    cout << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << ',';
    //cout << endl;
}

void DisplayGPDtoFile (GPD gpd, ofstream &ouf) {

    ouf << gpd.gname << '\t' << gpd.tname << '\t' << gpd.chro << '\t' << gpd.strand << '\t' << gpd.start << '\t' << gpd.end << '\t' << gpd.cstart << '\t';
    ouf << gpd.cend << '\t' << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  ouf << gpd.vec_exon_start[i] << ',';
    ouf << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  ouf << gpd.vec_exon_end[i] << ',';
    ouf << endl;
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


pair<int,int> AlignGPD (GPD gpd_query, GPD gpd_base) {

    int score = 0;
    int divergence = 0;
    if(gpd_query.chro != gpd_base.chro || gpd_query.nblock != gpd_base.nblock || gpd_query.strand != gpd_base.strand)  score = 0; // Add the third condition on 2017-10-09.
    else{
        if(gpd_query.nblock == 1){
            bool left_aligned = false;
            bool right_aligned = false;
            if(abs(gpd_query.start - gpd_base.start) < 100 || abs(gpd_query.start - gpd_base.cstart) < 100)  left_aligned = true;
            if(abs(gpd_query.end - gpd_base.end) < 100 || abs(gpd_query.end - gpd_base.cend) < 100)  right_aligned = true;
            divergence = abs(gpd_query.start - gpd_base.start) + abs(gpd_query.start - gpd_base.cstart);
            if(left_aligned && right_aligned)  score = 1;
            else  score = 0;
        }
        else{
            bool aligned = true;
            if(abs(gpd_query.vec_exon_end[0] - gpd_base.vec_exon_end[0]) > THRESHOLD)  aligned = false;
            divergence += abs(gpd_query.vec_exon_end[0] - gpd_base.vec_exon_end[0]);
            for(int i = 1; i <= gpd_query.nblock - 2; ++i){
                if(abs(gpd_query.vec_exon_start[i] - gpd_base.vec_exon_start[i]) > THRESHOLD)  aligned = false;
                if(abs(gpd_query.vec_exon_end[i] - gpd_base.vec_exon_end[i]) > THRESHOLD)  aligned = false;
                divergence = divergence + abs(gpd_query.vec_exon_start[i] - gpd_base.vec_exon_start[i]) + abs(gpd_query.vec_exon_end[i] - gpd_base.vec_exon_end[i]);
            }
            int last = gpd_query.nblock - 1;
            if(abs(gpd_query.vec_exon_start[last] - gpd_base.vec_exon_start[last]) > THRESHOLD)  aligned = false;
            divergence += abs(gpd_query.vec_exon_start[last] - gpd_base.vec_exon_start[last]);
            if(aligned)  score = 2;
            else  score = 0;
        }
    }

    return pair<int,int>(score,divergence);
}





int main (int argc, char **argv) {

    ifstream inf(argv[1]);
    ofstream ouf1(argv[2]);
    ofstream ouf2(argv[3]);
    vector<GPD> gpd_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            gpd_vec.push_back(gpd);
        }
    }

    inf.close();


    qsort(&gpd_vec[0],gpd_vec.size(),sizeof(GPD),compare_gpd_sort);


    vector<GPD> cluster;
    vector<int> indicator;
    vector<int> abundance;

    cluster.push_back(gpd_vec[0]);
    indicator.push_back(0);
    abundance.push_back(1);
    int index;


    for(int i = 1; i < gpd_vec.size(); ++i){

        index = cluster.size() - 1;
        bool Success = false;
        int match_pos = -1;

        while(index >= 0){
            if(cluster[index].chro != gpd_vec[i].chro || gpd_vec[i].start - cluster[index].end > 1000000)  break;
            else{
                pair<int,int> result = AlignGPD(gpd_vec[i], cluster[index]);
                if(result.first > 0){
                    Success = true;
                    match_pos = index;
                    break;
                }
            }
            --index;
        }

        if(Success){
            indicator.push_back(index);
            ++abundance[index];
        }
        else{
            cluster.push_back(gpd_vec[i]);
            indicator.push_back(cluster.size() - 1);
            abundance.push_back(1);
        }
    }



    int acculmulate = 0;
    for(int i = 0; i < cluster.size(); ++i){
        if(abundance[i] >= 2){
            ++acculmulate;
            cluster[i].tname = "CAGE_novel_" + lexical_cast<string>(acculmulate);
            DisplayGPDtoFile(cluster[i], ouf1);
        }
    }

    for(int i = 0; i < gpd_vec.size(); ++i){
        int n = indicator[i];
        if(abundance[n] >= 2){
            ouf2 << gpd_vec[i].tname << '\t' << cluster[n].gname << '\t' << cluster[n].tname << '\t' << cluster[n].strand << '\t';
            ouf2 << gpd_vec[i].chro << ':' << gpd_vec[i].start << '-' << gpd_vec[i].end << '\t';
            ouf2 << cluster[n].chro << ':' << cluster[n].start << '-' << cluster[n].end << '\t' << cluster[n].nblock << endl;
        }
    }


    inf.close();

    return 0;
}
