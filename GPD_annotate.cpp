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
    if(gpd_query.chro != gpd_base.chro || gpd_query.nblock != gpd_base.nblock || gpd_query.strand != gpd_base.strand)  score = 0;
    else{
        if(gpd_query.nblock == 1){
            bool left_aligned = false;
            bool right_aligned = false;
            if(abs(gpd_query.start - gpd_base.start) < 100 || abs(gpd_query.start - gpd_base.cstart) < 100)  left_aligned = true;
            if(abs(gpd_query.end - gpd_base.end) < 100 || abs(gpd_query.end - gpd_base.cend) < 100)  right_aligned = true;
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






int main (int argc, char **argv) {

    cerr << "GPD_annotate <annotation.gpd> <read_annotate.gpd> <prefix -STR>" << endl;

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


    //qsort(&gpd_array[0],line,sizeof(GPD),compare_gpd_sort);

    /*for(int i = 0; i < line; ++i){
        DisplayGPD(gpd_array[i]);
    }*/


    ifstream rinf(argv[2]);
    char sfilename[100];
    sprintf(sfilename, "%s_intact.txt", argv[3]);
    char ffilename[100];
    sprintf(ffilename, "%s_nonintact.gpd", argv[3]);
    ofstream souf(sfilename, ios::trunc);
    ofstream fouf(ffilename, ios::trunc);
    int line_num = 0;

    while(rinf){

        string strInput;
        getline(rinf, strInput);

        if(strInput.length() > 0){

            ++line_num;
            GPD gpd_query;
            BuildGPDFromString(strInput, gpd_query);
            int locate = BinarySearch(gpd_query, gpd_array, line);
            pair<int,int> align_gpd_result;
            vector<pair<GPD,int> > vec_record;

            //cout << strInput << endl;
            align_gpd_result = AlignGPD(gpd_query,gpd_array[locate]);
            if(align_gpd_result.first > 0){
                //cout << "Intact" << endl;
                //DisplayGPD(gpd_array[locate]);
                vec_record.push_back(pair<GPD,int>(gpd_array[locate],align_gpd_result.second));
            }

            int down_locate = locate;
            while(down_locate < line - 1){
                ++down_locate;
                if(gpd_array[down_locate].chro != gpd_query.chro || gpd_array[down_locate].start > gpd_query.end)  break;
                else{
                    align_gpd_result = AlignGPD(gpd_query,gpd_array[down_locate]);
                    if(align_gpd_result.first > 0){
                        //cout << "Intact" << endl;
                        //DisplayGPD(gpd_array[down_locate]);
                        vec_record.push_back(pair<GPD,int>(gpd_array[down_locate],align_gpd_result.second));
                    }
                }
            }

            int up_locate = locate;
            while(up_locate > 0){
                --up_locate;
                if(gpd_array[up_locate].chro != gpd_query.chro || locate - up_locate > 100 || gpd_query.start - gpd_array[up_locate].end > 1000000)  break;
                else{
                    align_gpd_result = AlignGPD(gpd_query,gpd_array[up_locate]);
                    if(align_gpd_result.first > 0){
                        //cout << "Intact" << endl;
                        //DisplayGPD(gpd_array[up_locate]);
                        vec_record.push_back(pair<GPD,int>(gpd_array[up_locate],align_gpd_result.second));
                    }
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

                souf << gpd_query.gname << '\t' << vec_record[min_index].first.gname << '\t' << vec_record[min_index].first.tname << '\t';
                souf << vec_record[min_index].first.strand << '\t' << gpd_query.strand << '\t';
                souf << gpd_query.chro << ':' << gpd_query.start << '-' << gpd_query.end << '\t';
                souf << vec_record[min_index].first.chro << ':' << vec_record[min_index].first.start << '-' << vec_record[min_index].first.end << '\t' << gpd_query.nblock << endl;
            }
            else  fouf << strInput << endl;
        }
    }

    rinf.close();

    return 0;
}
