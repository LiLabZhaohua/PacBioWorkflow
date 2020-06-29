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
    //int cstart;
    //int cend;
    string cstart;
    string cend;
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
    gpd.cstart = vec_temp[6];
    gpd.cend = vec_temp[7];
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


vector<GPD> BuildGPDVec (ifstream &inf) {

    vector<GPD> gpd_array;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            gpd_array.push_back(gpd);
        }
    }

    //qsort(&gpd_array[0],gpd_array.size(),sizeof(GPD),compare_gpd_sort);

    return gpd_array;
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


int BinarySearch (GPD gpd_query, vector<GPD> &gpd_array) {

    int up = 0;
    int down = gpd_array.size() - 1;
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
    if(gpd_query.chro != gpd_base.chro || gpd_query.nblock != gpd_base.nblock || gpd_query.strand != gpd_base.strand)  score = 0;
    else{
        if(gpd_query.nblock == 1){
            bool left_aligned = false;
            bool right_aligned = false;
            if(abs(gpd_query.start - gpd_base.start) < 130)  left_aligned = true;
            if(abs(gpd_query.end - gpd_base.end) < 130)  right_aligned = true;
            divergence = abs(gpd_query.start - gpd_base.start) + abs(gpd_query.end - gpd_base.end);
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


bool PartialAlignGPD (GPD gpd_query, GPD gpd_base) {

    if(gpd_query.chro != gpd_base.chro || gpd_query.strand != gpd_base.strand)  return false;
    else{
        if(gpd_query.nblock == 1){
            if(gpd_base.nblock == 1){
                if(gpd_query.start > gpd_base.vec_exon_start[0] - 130 && gpd_query.end < gpd_base.vec_exon_end[0] + 130)  return true;
                else  return false;
            }
            else{
                bool aligned = false;
                for(int i = 0; i < gpd_base.nblock; ++i){
                    if(i == 0){
                        if(gpd_query.start > gpd_base.vec_exon_start[i] - 130 && gpd_query.end <= gpd_base.vec_exon_end[i] + THRESHOLD){
                            aligned = true;
                            break;
                        }
                    }
                    else if(i == gpd_base.nblock - 1){
                        if(gpd_query.end < gpd_base.vec_exon_end[i] + 130 && gpd_query.start >= gpd_base.vec_exon_start[i] - THRESHOLD){
                            aligned = true;
                            break;
                        }
                    }
                    else{
                        if(gpd_query.end <= gpd_base.vec_exon_end[i] + THRESHOLD && gpd_query.start >= gpd_base.vec_exon_start[i] - THRESHOLD){
                            aligned = true;
                            break;
                        }
                    }
                }
                return aligned;
            }
        }
        else{
            bool anchor = false;
            bool aligned = true;
            for(int i = 0; i < gpd_base.nblock; ++i){
                if(abs(gpd_query.vec_exon_end[0] - gpd_base.vec_exon_end[i]) <= THRESHOLD){
                    anchor = true;
                    if(i + gpd_query.nblock <= gpd_base.nblock){
                        if((i == 0 && gpd_query.vec_exon_start[0] < gpd_base.vec_exon_start[i] - 130) || (i > 0 && gpd_query.vec_exon_start[0] < gpd_base.vec_exon_start[i] - THRESHOLD))  aligned = false;
                        for(int j = 1; j <= gpd_query.nblock - 2; ++j){
                            if(abs(gpd_query.vec_exon_start[j] - gpd_base.vec_exon_start[j+i]) > THRESHOLD)  aligned = false;
                            if(abs(gpd_query.vec_exon_end[j] - gpd_base.vec_exon_end[j+i]) > THRESHOLD)  aligned = false;
                        }
                        int last = gpd_query.nblock - 1;
                        if(abs(gpd_query.vec_exon_start[last] - gpd_base.vec_exon_start[last+i]) > THRESHOLD)  aligned = false;
                        if((last+i == gpd_base.nblock-1 && gpd_query.vec_exon_end[last] > gpd_base.vec_exon_end[last+i] + 130) || (last+i < gpd_base.nblock-1 && gpd_query.vec_exon_end[last] > gpd_base.vec_exon_end[last+i] + THRESHOLD))  aligned = false;
                    }
                    else  aligned = false;
                }
            }

            if(anchor && aligned)  return true;
            else  return false;
        }
    }
}


string FullLengthMap (GPD gpd_query, vector<GPD> &gpd_array, int locate) {

    pair<int,int> align_gpd_result;
    vector<pair<GPD,int> > vec_record;

    align_gpd_result = AlignGPD(gpd_query,gpd_array[locate]);
    if(align_gpd_result.first > 0)  vec_record.push_back(pair<GPD,int>(gpd_array[locate],align_gpd_result.second));

    int down_locate = locate;
    while(down_locate < gpd_array.size() - 1){
        ++down_locate;
        if(gpd_array[down_locate].chro != gpd_query.chro || gpd_array[down_locate].start > gpd_query.end)  break;
        else{
            align_gpd_result = AlignGPD(gpd_query,gpd_array[down_locate]);
            if(align_gpd_result.first > 0)  vec_record.push_back(pair<GPD,int>(gpd_array[down_locate],align_gpd_result.second));
        }
    }

    int up_locate = locate;
    while(up_locate > 0){
        --up_locate;
        if(gpd_array[up_locate].chro != gpd_query.chro || locate - up_locate > 100 || gpd_query.start - gpd_array[up_locate].end > 1000000)  break;
        else{
            align_gpd_result = AlignGPD(gpd_query,gpd_array[up_locate]);
            if(align_gpd_result.first > 0)  vec_record.push_back(pair<GPD,int>(gpd_array[up_locate],align_gpd_result.second));
        }
    }

    int min_score = 100000;
    int min_index = 0;
    if(vec_record.size() > 0){
        for(int i = 0; i < vec_record.size(); ++i){
            if(vec_record[i].second < min_score){
                min_score = vec_record[i].second;
                min_index = i;
            }
        }
        return vec_record[min_index].first.tname;
    }
    else  return "None";
}


vector<GPD> PartialLengthMap (GPD gpd_query, vector<GPD> &gpd_array, int locate) {

    bool align_gpd_result;
    vector<GPD> vec_record;

    align_gpd_result = PartialAlignGPD(gpd_query,gpd_array[locate]);
    if(align_gpd_result)  vec_record.push_back(gpd_array[locate]);

    int down_locate = locate;
    while(down_locate < gpd_array.size() - 1){
        ++down_locate;
        if(gpd_array[down_locate].chro != gpd_query.chro || gpd_array[down_locate].start > gpd_query.end)  break;
        else{
            align_gpd_result = PartialAlignGPD(gpd_query,gpd_array[down_locate]);
            if(align_gpd_result)  vec_record.push_back(gpd_array[down_locate]);
        }
    }

    int up_locate = locate;
    while(up_locate > 0){
        --up_locate;
        if(gpd_array[up_locate].chro != gpd_query.chro || locate - up_locate > 100 || gpd_query.start - gpd_array[up_locate].end > 1000000)  break;
        else{
            align_gpd_result = PartialAlignGPD(gpd_query,gpd_array[up_locate]);
            if(align_gpd_result)  vec_record.push_back(gpd_array[up_locate]);
        }
    }


    return vec_record;
}








int main (int argc, char **argv) {

    cerr << "MapReadToIsoform <annotation.gpd> <read_annotate.gpd>" << endl;

    ifstream inf(argv[1]);
    vector<GPD> gpd_array = BuildGPDVec(inf);
    inf.close();
    inf.open(argv[1]);
    map<string,GPD> gpd_map = BuildGPDMap(inf);
    inf.close();

    map<string,GPD>::iterator it;
    ifstream rinf(argv[2]);


    while(rinf){

        string strInput;
        getline(rinf, strInput);

        if(strInput.length() > 0){

            GPD gpd_query;
            BuildGPDFromString(strInput, gpd_query);
            int locate = BinarySearch(gpd_query, gpd_array);

            /*string full_length_map = FullLengthMap(gpd_query, gpd_array, locate);

            if(full_length_map == "None"){

                vector<GPD> partial_length_map = PartialLengthMap(gpd_query, gpd_array, locate);

                if(partial_length_map.size() > 0){
                    for(int i = 0; i < partial_length_map.size(); ++i){
                    //for(int i = 0; i < 1; ++i){
                        cout << gpd_query.tname << '\t' << partial_length_map[i].tname << "\tPartial";
                        it = gpd_map.find(partial_length_map[i].tname);
                        cout << '\t' << (it->second).strand << '\t' << gpd_query.strand << endl;
                    }
                }
                else  cout << gpd_query.tname << "\tUnfound\t*\t*\t*" << endl;
            }
            else{
                cout << gpd_query.tname << '\t' << full_length_map << "\tFull-length";
                it = gpd_map.find(full_length_map);
                cout << '\t' << (it->second).strand << '\t' << gpd_query.strand << endl;
            }*/

            vector<GPD> partial_length_map = PartialLengthMap(gpd_query, gpd_array, locate);

            if(partial_length_map.size() > 0){
                for(int i = 0; i < partial_length_map.size(); ++i){
                    cout << gpd_query.tname << '\t' << partial_length_map[i].tname << "\tCompatible";
                    it = gpd_map.find(partial_length_map[i].tname);
                    cout << '\t' << (it->second).strand << '\t' << gpd_query.strand << endl;
                }
            }
            else  cout << gpd_query.tname << "\tUnfound\t*\t*\t*" << endl;
        }
    }

    rinf.close();

    return 0;
}
