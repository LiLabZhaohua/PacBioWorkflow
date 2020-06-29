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
    cout << endl;
}


void DisplayGPDToFile (GPD gpd, ofstream &ouf) {

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


map<string,GPD> BuildGPDMap (ifstream &inf) {

    map<string,GPD> gpd_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            gpd_map[gpd.tname] = gpd;
        }
    }

    return gpd_map;
}


pair<int,int> ModifyEnds (string &cigar, int mapped_pos) {

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
            else if(cigar[i] == 'S'){
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

    return pair<int,int>(mapped_pos - left_S, mapped_pos + aln_len + right_S);
}


void Polish (map<string,GPD> &gpd_map, ifstream &sinf) {

    map<string,GPD>::iterator it;

    while(sinf){

        string strInput;
        getline(sinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            pair<int,int> polished_ends = ModifyEnds(vec[5], lexical_cast<int>(vec[3]));
            string read_name = vec[0];

            it = gpd_map.find(read_name);
            if(it == gpd_map.end())  cout << "Alert" << endl;

            (it->second).start = polished_ends.first;
            (it->second).cstart = polished_ends.first;
            (it->second).end = polished_ends.second;
            (it->second).cend = polished_ends.second;
            int block_num = (it->second).nblock;
            (it->second).vec_exon_start[0] = polished_ends.first;
            (it->second).vec_exon_end[block_num - 1] = polished_ends.second;
        }
    }
}


vector<GPD> TransferMapToVec (map<string, GPD> &gpd_map) {

    vector<GPD> gpd_vec;
    map<string,GPD>::iterator it;
    for(it = gpd_map.begin(); it != gpd_map.end(); ++it)  gpd_vec.push_back(it->second);
    return gpd_vec;
}


bool StrictAlignGPD (GPD gpd_query, GPD gpd_base) {

    bool aligned = true;

    if(gpd_query.nblock != gpd_base.nblock)  aligned = false;
    else{
        for(int i = 0; i < gpd_query.nblock; ++i){
            if(i == 0 && abs(gpd_query.vec_exon_start[i] - gpd_base.vec_exon_start[i]) > 30)  aligned = false;
            else if(i > 0 && abs(gpd_query.vec_exon_start[i] - gpd_base.vec_exon_start[i]) > THRESHOLD)  aligned = false;
            else  ;
            if(i == gpd_query.nblock - 1 && abs(gpd_query.vec_exon_end[i] - gpd_base.vec_exon_end[i]) > 30)  aligned = false;
            else if(i < gpd_query.nblock - 1 && abs(gpd_query.vec_exon_end[i] - gpd_base.vec_exon_end[i]) > THRESHOLD)  aligned = false;
            else  ;
        }
     }

    return aligned;
}






int main (int argc, char **argv) {

    cerr << "ClusterUniqueReads <read_aln.gpd> <read_aln.sam>" << endl;

    ifstream inf(argv[1]);
    map<string,GPD> gpd_map = BuildGPDMap(inf);
    inf.close();

    ifstream sinf(argv[2]);
    Polish(gpd_map, sinf);
    sinf.close();

    vector<GPD> gpd_vec = TransferMapToVec(gpd_map);
    //qsort(&gpd_vec[0], gpd_vec.size(), sizeof(GPD), compare_gpd_sort);

    ofstream ouf("POLISHED");
    for(int i = 0; i < gpd_vec.size(); ++i)  DisplayGPDToFile(gpd_vec[i], ouf);
    ouf.close();

    system("sort -k 3,3 -k 5,5n POLISHED > POLISHED_sort");
    //exit(1);
    vector<GPD> gpd_vec_sort;
    ifstream pinf("POLISHED_sort");

    while(pinf){

        string strInput;
        getline(pinf, strInput);
        if(strInput.length() == 0)  break;

        GPD gpd;
        BuildGPDFromString(strInput, gpd);
        gpd_vec_sort.push_back(gpd);
    }


    vector<GPD> cluster;
    vector<int> indicator;
    vector<int> abundance;

    cluster.push_back(gpd_vec_sort[0]);
    indicator.push_back(0);
    abundance.push_back(1);

    for(int i = 1; i < gpd_vec_sort.size(); ++i){

        int index = cluster.size() - 1;
        bool Success = false;
        int matched_cluster = -1;

        while(index >= 0){
            if(cluster[index].chro != gpd_vec_sort[i].chro || gpd_vec_sort[i].start - cluster[index].end > 50000)  break;
            else{
                bool strict_aln_result = StrictAlignGPD(gpd_vec_sort[i], cluster[index]);
                if(strict_aln_result){
                    Success = true;
                    matched_cluster = index;
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
            cluster.push_back(gpd_vec_sort[i]);
            indicator.push_back(cluster.size() - 1);
            abundance.push_back(1);
        }
    }

    if(cluster.size() != abundance.size()){
        cout << "Alert" << endl;
        exit(1);
    }


    //for(int i = 0; i < gpd_vec_sort.size(); ++i){
    //    cout << gpd_vec_sort[i].tname << "\tCluster_" << indicator[i] << '\t' << abundance[indicator[i]] << endl;
    //}

    for(int i = 0; i < cluster.size(); ++i){
        cout << cluster[i].gname << "\tcluster_" << i << endl;
    }

    return 0;
}
