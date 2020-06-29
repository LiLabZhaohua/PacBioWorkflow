#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct EXON {

    string chro;
    int start;
    int end;
    string gname;
};


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


int GPD_Length (GPD gpd) {

    int gpd_len = 0;
    for(int i = 0; i < gpd.nblock; ++i)  gpd_len += (gpd.vec_exon_end[i] - gpd.vec_exon_start[i] + 1);
    return gpd_len;
}


void ParseExon (ifstream &inf, ofstream &ouf) {

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            for(int i = 0; i < gpd.nblock; ++i)  ouf << gpd.chro << '\t' << gpd.vec_exon_start[i] << '\t' << gpd.vec_exon_end[i] << '\t' << gpd.gname << endl;
        }
    }
}


int compare_exon_sort (const void *a, const void *b) {

    if( (*(EXON*) a).chro > (*(EXON*) b).chro )  return 1;
    else if ( (*(EXON*) a).chro == (*(EXON*) b).chro ){
        if( (*(EXON*) a).start > (*(EXON*) b).start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_two_exon (EXON a, EXON b) {

    if(a.chro > b.chro)  return 1;
    else if (a.chro == b.chro){
        if(a.start > b.start )  return 1;
        else  return -1;
    }
    else  return -1;
}


vector<EXON> BuildExonVectorFromFile (ifstream &inf) {

    vector<EXON> exon_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            EXON exon_rc;
            exon_rc.chro = vec[0];
            exon_rc.start = lexical_cast<int>(vec[1]);
            exon_rc.end = lexical_cast<int>(vec[2]);
            exon_rc.gname = vec[3];
            exon_vec.push_back(exon_rc);
        }
    }

    qsort(&exon_vec[0], exon_vec.size(), sizeof(EXON), compare_exon_sort);

    return exon_vec;
}


int BinarySearchExon (EXON query, vector<EXON> &array) {

    int up = 0;
    int down = array.size() - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_two_exon(array[middle],query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}


int Overlap (EXON a, EXON b) {

    if(a.chro != b.chro)  return -1;
    else{
        if(a.end <= b.start || a.start >= b.end)  return -1;
        else{
            int real_start, real_end;
            if(a.start > b.start)  real_start = a.start;
            else  real_start = b.start;
            if(a.end > b.end)  real_end = b.end;
            else  real_end = a.end;
            return (real_end - real_start);
        }
    }
}





int main (int argc, char **argv) {

    cerr << "ExonCompare <annotation.gpd> <reads.gpd>" << endl;

    ifstream annoinf(argv[1]);
    ofstream psouf("PARSEEXON", ios::trunc);
    ParseExon(annoinf, psouf);
    annoinf.close(); psouf.close();

    char cmd_sort[100];
    sprintf(cmd_sort, "sort -u PARSEEXON > PARSEEXON_uniq");
    system(cmd_sort);

    ifstream uniqinf("PARSEEXON_uniq");
    vector<EXON> exon_vec = BuildExonVectorFromFile(uniqinf);
    uniqinf.close();
    ifstream gpdinf(argv[2]);

    while(gpdinf){

        string strInput;
        getline(gpdinf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            int total_overlap_len = 0;
            map<string,int> gene_overlap_map;
            map<string,int>::iterator it;

            for(int i = 0; i < gpd.nblock; ++i){

                EXON exon_rc;
                exon_rc.chro = gpd.chro;
                exon_rc.start = gpd.vec_exon_start[i];
                exon_rc.end = gpd.vec_exon_end[i];
                exon_rc.gname = gpd.gname;

                int locate = BinarySearchExon(exon_rc, exon_vec);
                int max_overlap = 0;
                int max_index = -1;

                for(int j = -10; j <= 10; ++j){
                    int index = locate + j;
                    if(index >= 0 && index < exon_vec.size()){
                        int overlap = Overlap(exon_rc, exon_vec[index]);
                        if(overlap > max_overlap){
                            max_overlap = overlap;
                            max_index = index;
                        }
                    }
                }

                if(max_index != -1){
                    total_overlap_len += max_overlap;
                    it = gene_overlap_map.find(exon_vec[max_index].gname);
                    if(it == gene_overlap_map.end())  gene_overlap_map[exon_vec[max_index].gname] = max_overlap;
                    else  (it->second) += max_overlap;
                }
            }

            int gpd_len = GPD_Length(gpd);
            double ratio = total_overlap_len / (gpd_len + 0.0);

            int gene_max_overlap = 0;
            string max_gene_name = "*";
            for(it = gene_overlap_map.begin(); it != gene_overlap_map.end(); ++it){
                if((it->second) > gene_max_overlap){
                    gene_max_overlap = (it->second);
                    max_gene_name = (it->first);
                }
            }

            cout << gpd.gname << '\t' << max_gene_name << '\t' << ratio << endl;
        }
    }

    gpdinf.close();

    return 0;
}
