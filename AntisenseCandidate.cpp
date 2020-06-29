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


struct GENE {

    string name;
    string strand;
    string chro;
    int start;
    int end;
};


map<string,GENE> BuildGeneMapFromFile (ifstream &inf) {

    map<string,GENE> gene_map;
    map<string,GENE>::iterator it;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            it = gene_map.find(gpd.gname);

            if(it == gene_map.end()){
                GENE gene;
                gene.name = gpd.gname;
                gene.strand = gpd.strand;
                gene.chro = gpd.chro;
                gene.start = gpd.start;
                gene.end = gpd.end;
                gene_map[gpd.gname] = gene;
            }
            else{
                if(gpd.strand != (it->second).strand){
                    (it->second).start = -1;
                    (it->second).end = -1;
                }
                else if((it->second).start > 0 && (it->second).end > 0){
                    if(gpd.start < (it->second).start)  (it->second).start = gpd.start;
                    if(gpd.end > (it->second).end)  (it->second).end = gpd.end;
                }
                else  ;
            }
        }
    }

    return gene_map;
}


int compare_gene_sort (const void *a, const void *b) {

    if( (*(GENE*) a).chro > (*(GENE*) b).chro )  return 1;
    else if ( (*(GENE*) a).chro == (*(GENE*) b).chro ){
        if( (*(GENE*) a).start > (*(GENE*) b).start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_gene (GENE a, GPD b) {

    if(a.chro > b.chro)  return 1;
    else if(a.chro == b.chro){
        if(a.start > b.start)  return 1;
        else  return -1;
    }
    else  return -1;
}


int BinarySearch (GPD gpd, vector<GENE> &gene_vec) {

    int up = 0;
    int down = gene_vec.size() - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_gene(gene_vec[middle], gpd) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}


int ReverseOverlapLength (GENE gene, GPD gpd) {

    int overlap_len = 0;

    if(gene.strand != gpd.strand && gene.chro == gpd.chro){
        int left, right;
        if(gene.start > gpd.start)  left = gene.start;
        else  left = gpd.start;
        if(gene.end > gpd.end)  right = gpd.end;
        else  right = gene.end;
        if(right > left)  overlap_len = right - left;
    }

    return overlap_len;
}





int main (int argc, char **argv) {

    ifstream annoinf(argv[1]);
    map<string,GENE> gene_map = BuildGeneMapFromFile(annoinf);
    annoinf.close();

    vector<GENE> gene_vec;
    map<string,GENE>::iterator it;

    for(it = gene_map.begin(); it != gene_map.end(); ++it){
        GENE rc;
        rc.name = (it->second).name;
        rc.strand = (it->second).strand;
        rc.chro = (it->second).chro;
        rc.start = (it->second).start;
        rc.end = (it->second).end;
        gene_vec.push_back(rc);
    }

    qsort(&gene_vec[0], gene_vec.size(), sizeof(GENE), compare_gene_sort);

    /*for(int i = 0; i < gene_vec.size(); ++i){
        cout << gene_vec[i].name << '\t' << gene_vec[i].strand << '\t' << gene_vec[i].chro << '\t' << gene_vec[i].start << '\t' << gene_vec[i].end << endl;
    }*/


    ifstream rdinf(argv[2]);

    while(rdinf){

        string strInput;
        getline(rdinf, strInput);

        if(strInput.length() > 0){
            
            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            int locate = BinarySearch(gpd, gene_vec);
            int overlap_len = ReverseOverlapLength(gene_vec[locate], gpd);
            int max_overlap_len = overlap_len;
            int max_pos = locate;

            int down_locate = locate;
            while(down_locate < gene_vec.size()){
                ++down_locate;
                if(gene_vec[down_locate].chro != gpd.chro || gene_vec[down_locate].start > gpd.end)  break;
                else{
                    overlap_len = ReverseOverlapLength(gene_vec[down_locate], gpd);
                    if(overlap_len > max_overlap_len){
                        max_overlap_len = overlap_len;
                        max_pos = down_locate;
                    }
                }
            }

            int up_locate = locate;
            while(up_locate > 0){
                --up_locate;
                if(gene_vec[up_locate].chro != gpd.chro || locate - up_locate > 200 || gpd.start - gene_vec[up_locate].end > 50000)  break;
                else{
                    overlap_len = ReverseOverlapLength(gene_vec[down_locate], gpd);
                    if(overlap_len > max_overlap_len){
                        max_overlap_len = overlap_len;
                        max_pos = down_locate;
                    }
                }
            }

            if(max_overlap_len > 0){
                cout << gpd.gname << '\t' << gpd.strand << '\t' << gene_vec[max_pos].name << '\t' << gene_vec[max_pos].strand << '\t';
                cout << gpd.chro << '\t' << gpd.start << '\t' << gpd.end << '\t';
                cout << gene_vec[max_pos].chro << '\t' << gene_vec[max_pos].start << '\t' << gene_vec[max_pos].end << '\t';
                cout << max_overlap_len << '\t';
                double ratio = max_overlap_len / (gpd.end - gpd.start + 0.0);
                cout << ratio << endl;
            }
            else  cout << gpd.gname << '\t' << gpd.strand << "\t*\t*\t" << gpd.chro << '\t' << gpd.start << '\t' << gpd.end << "\t*\t*\t*\t0\t0" << endl;
        }
    }


    return 0;
}
