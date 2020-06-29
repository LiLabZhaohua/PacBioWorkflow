#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
//#include <unordered_set>
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


void OutputGPDinLine (ofstream &ouf, GPD gpd) {

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


/*void PrintGPDinLine (GPD gpd) {

    cout << gpd.gname << '\t' << gpd.tname << '\t' << gpd.chro << '\t' << gpd.strand << '\t' << gpd.start << '\t' << gpd.end << '\t' << gpd.cstart << '\t';
    cout << gpd.cend << '\t' << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << ',';
    cout << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << ',';
    cout << endl;
}*/


pair<int,int> AlignGPD (GPD gpd_query, GPD gpd_base) {

    int score = 0;
    int divergence = 0;
    if(gpd_query.chro != gpd_base.chro || gpd_query.nblock != gpd_base.nblock)  score = 0;
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


bool IsSubset (GPD gpd_query, GPD gpd_base) {

    if(gpd_query.chro != gpd_base.chro)  return false;
    else if(gpd_query.nblock >= gpd_base.nblock)  return false;
    else{

        bool Subset = false;

        for(int shift = 0; shift <= gpd_base.nblock - gpd_query.nblock; ++shift){

            bool aligned = true;
            if(abs(gpd_query.vec_exon_end[0] - gpd_base.vec_exon_end[shift]) > THRESHOLD)  aligned = false;
            for(int i = 1; i <= gpd_query.nblock - 2; ++i){
                if(abs(gpd_query.vec_exon_start[i] - gpd_base.vec_exon_start[i+shift]) > THRESHOLD)  aligned = false;
                if(abs(gpd_query.vec_exon_end[i] - gpd_base.vec_exon_end[i+shift]) > THRESHOLD)  aligned = false;
            }
            int last = gpd_query.nblock - 1;
            if(abs(gpd_query.vec_exon_start[last] - gpd_base.vec_exon_start[last+shift]) > THRESHOLD)  aligned = false;

            if(aligned){
                Subset = true;
                break;
            }
        }

        return Subset;
    }
}





int main (int argc, char **argv) {

    cerr << "AssignTranscript <ref_inf.gpd> <query_inf.gpd> <output.gpd>" << endl;

    ifstream rinf(argv[1]);
    ifstream qinf(argv[2]);
    ofstream ouf(argv[3], ios::trunc);
    map<string,vector<GPD> > ref_map;
    map<string,vector<GPD> >::iterator it_ref;
    map<string,vector<GPD> > q_map;
    map<string,vector<GPD> >::iterator it_q;

    while(rinf){

        string strInput;
        getline(rinf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            ref_map[gpd.gname].push_back(gpd);
        }
    }


    /*for(it_ref = ref_map.begin(); it_ref != ref_map.end(); ++it_ref){
        qsort(&((it_ref->second)[0]),(it_ref->second).size(),sizeof(GPD),compare_gpd_sort);
    }*/


    while(qinf){

        string strInput;
        getline(qinf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            q_map[gpd.gname].push_back(gpd);
        }
    }


    /*for(it_q = q_map.begin(); it_q != q_map.end(); ++it_q){
        qsort(&((it_q->second)[0]),(it_q->second).size(),sizeof(GPD),compare_gpd_sort);
    }*/


    for(it_q = q_map.begin(); it_q != q_map.end(); ++it_q){

        it_ref = ref_map.find(it_q->first);

        if(it_ref == ref_map.end()){

            vector<GPD>::iterator pnt;
            for(pnt = (it_q->second).begin(); pnt != (it_q->second).end(); ++pnt){
                OutputGPDinLine(ouf, *pnt);
            }
        }
        else{

            vector<GPD>::iterator pnt_q;

            for(pnt_q = (it_q->second).begin(); pnt_q != (it_q->second).end(); ++pnt_q){

                bool aligned = false;
                vector<int> divergence_vec;
                for(int i = 0; i < (it_ref->second).size(); ++i)  divergence_vec.push_back(1000000);

                for(int j = 0; j < (it_ref->second).size(); ++j){
                    pair<int,int> aln_result = AlignGPD(*pnt_q, (it_ref->second)[j]);
                    if(aln_result.first > 0){
                        aligned = true;
                        divergence_vec[j] = aln_result.second;
                        //break;
                    }
                }

                if(aligned){
                    int min_value = 1000000;
                    int min_index = -1;
                    for(int i = 0; i < (it_ref->second).size(); ++i){
                        if(divergence_vec[i] < min_value){
                            min_value = divergence_vec[i];
                            min_index = i;
                        }
                    }
                    //(*pnt_q).tname = (it_ref->second)[min_index].tname;
                    OutputGPDinLine(ouf, (it_ref->second)[min_index]);
                }
                else  OutputGPDinLine(ouf, *pnt_q);

                /*if(!aligned){
                    if((*pnt_q).tname.substr(0,8) == "LR_novel")  (*pnt_q).tname = "ad_" + (*pnt_q).tname;
                    OutputGPDinLine(ouf, *pnt_q);
                }*/
            }
        }
    }


    rinf.close();
    qinf.close();
    ouf.close();

    return 0;
}
