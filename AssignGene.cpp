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


struct LOCUS {

    string gname;
    string chro;
    int start;
    int end;
};


void DisplayLocus (LOCUS locus) {

    cout << locus.gname << '\t' << locus.chro << '\t' << locus.start << '\t' << locus.end << endl;
}


void BuildLocusVector (ifstream &inf, vector<LOCUS> *vec_locus) {

    LOCUS locus;
    string pre_gname;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[0] != pre_gname){
                if(locus.gname != "")  (*vec_locus).push_back(locus);
                locus.gname = vec[0];
                locus.chro = vec[2];
                locus.start = lexical_cast<int>(vec[4]);
                locus.end = lexical_cast<int>(vec[5]);
            }
            else{
                int temp_start = lexical_cast<int>(vec[4]);
                int temp_end = lexical_cast<int>(vec[5]);
                if(temp_start < locus.start)  locus.start = temp_start;
                if(temp_end > locus.end)  locus.end = temp_end;
            }

            pre_gname = vec[0];
        }
    }

    (*vec_locus).push_back(locus);
}


int compare_locus_sort (const void *a, const void *b) {

    if( (*(LOCUS*) a).chro > (*(LOCUS*) b).chro )  return 1;
    else if ( (*(LOCUS*) a).chro == (*(LOCUS*) b).chro ){
        if( (*(LOCUS*) a).start > (*(LOCUS*) b).start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_two_locus (LOCUS a, LOCUS b) {

    if(a.chro > b.chro)  return 1;
    else if (a.chro == b.chro){
        if(a.start > b.start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int BinarySearchLOCUS (LOCUS locus, vector<LOCUS> *vec_locus, int total) {

    int up = 0;
    int down = total - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_two_locus((*vec_locus)[middle],locus) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}


pair<int,double> AlignLOCUS (LOCUS query, LOCUS base) {

    int score = 0;
    double r = 0.0;
    if(query.chro != base.chro)  score = 0;
    else{
        if(query.end < base.start || query.start > base.end)  score = 0;
        else{
            int a, b;
            if(query.start > base.start)  a = query.start;
            else  a = base.start;
            if(query.end < base.end)  b = query.end;
            else  b = base.end;

            if(b < a)  cout << "Alert" << endl;
            else{
                int overlap = b - a;
                //cout << overlap << endl;
                double r1 = overlap / (query.end - query.start + 0.0);
                double r2 = overlap / (base.end - base.start + 0.0);
                score = 1;
                r = r1 + r2;
            }
        }
    }

    return pair<int,double>(score,r);
}







int main (int argc, char **argv) {

    ifstream inf1(argv[1]);
    ifstream inf2(argv[2]);

    vector<LOCUS> vec_locus1;
    vector<LOCUS> vec_locus2;

    BuildLocusVector(inf1, &vec_locus1);

    //qsort(&vec_locus1[0],vec_locus1.size(),sizeof(LOCUS),compare_locus_sort);

    /*for(int i = 0; i < vec_locus1.size(); ++i){
        DisplayLocus(vec_locus1[i]);
    }*/

    int success_assign = 0;
    int novel_gene_num = 0;


    while(inf2){

        string strInput;
        getline(inf2, strInput);
        if(strInput.length() == 0)  break;

        LOCUS query;
        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        query.gname = vec[0];
        query.chro = vec[2];
        query.start = lexical_cast<int>(vec[4]);
        query.end = lexical_cast<int>(vec[5]);

        int locate = BinarySearchLOCUS(query, &vec_locus1, vec_locus1.size());
        //DisplayLocus(vec_locus1[locate]);
        pair<int,double> align_result;
        vector<pair<LOCUS,double> > vec_record;

        align_result = AlignLOCUS(query, vec_locus1[locate]);
        //cout << align_result.first << '\t' << align_result.second << endl;
        if(align_result.first > 0){
            vec_record.push_back(pair<LOCUS,double>(vec_locus1[locate],align_result.second));
        }

        int down_locate = locate;
        while(down_locate < vec_locus1.size() - 1){
            ++down_locate;
            if(vec_locus1[down_locate].chro != query.chro || vec_locus1[down_locate].start > query.end)  break;
            else{
                //DisplayLocus(vec_locus1[down_locate]);
                align_result = AlignLOCUS(query, vec_locus1[down_locate]);
                //cout << align_result.first << '\t' << align_result.second << endl;
                if(align_result.first > 0){
                    vec_record.push_back(pair<LOCUS,double>(vec_locus1[down_locate],align_result.second));
                }
            }
        }

        int up_locate = locate;
        while(up_locate > 0){
            --up_locate;
            if(vec_locus1[up_locate].chro != query.chro || locate - up_locate > 100 || query.start - vec_locus1[up_locate].end > 1000000)  break;
            else{
                //DisplayLocus(vec_locus1[up_locate]);
                align_result = AlignLOCUS(query, vec_locus1[up_locate]);
                //cout << align_result.first << '\t' << align_result.second << endl;
                if(align_result.first > 0){
                    vec_record.push_back(pair<LOCUS,double>(vec_locus1[up_locate],align_result.second));
                }
            }
        }

        double max_r = 0.0;
        int max_index = 0;

        if(vec_record.size() > 0){
            for(int i = 0; i < vec_record.size(); ++i){
                if(vec_record[i].second > max_r){
                    max_r = vec_record[i].second;
                    max_index = i;
                }
            }
        }


        if(max_r > 1.0){
            //cout << vec[0] << '\t' << vec[1] << '\t' << vec_record[max_index].first.gname << endl;
            cout << vec_record[max_index].first.gname;
            for(int i = 1; i < vec.size(); ++i)  cout << '\t' << vec[i];
            cout << endl;
            ++success_assign;
        }
        else{
            //cout << strInput << endl;
            ++novel_gene_num;
            cout << "Novel_gene_" << novel_gene_num;
            for(int i = 1; i < vec.size(); ++i)  cout << '\t' << vec[i];
            cout << endl;
        }
    }


    inf1.close();
    inf2.close();

    return 0;
}
