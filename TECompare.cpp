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


struct TE {

    string chro;
    int start;
    int end;
    string strand;
    string name;
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


int compare_te_sort (const void *a, const void *b) {

    if( (*(TE*) a).chro > (*(TE*) b).chro )  return 1;
    else if ( (*(TE*) a).chro == (*(TE*) b).chro ){
        if( (*(TE*) a).start > (*(TE*) b).start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_two_te (TE a, TE b) {

    if(a.chro > b.chro)  return 1;
    else if (a.chro == b.chro){
        if(a.start > b.start )  return 1;
        else  return -1;
    }
    else  return -1;
}


vector<TE> BuildTEVectorFromFile (ifstream &inf) {

    vector<TE> te_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            TE te_rc;
            te_rc.chro = vec[0];
            te_rc.start = lexical_cast<int>(vec[1]);
            te_rc.end = lexical_cast<int>(vec[2]);
            te_rc.strand = vec[3];
            te_rc.name = vec[4];
            te_vec.push_back(te_rc);
        }
    }

    qsort(&te_vec[0], te_vec.size(), sizeof(TE), compare_te_sort);

    return te_vec;
}


int BinarySearchTE (TE query, vector<TE> &array) {

    int up = 0;
    int down = array.size() - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_two_te(array[middle],query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}


int Overlap (TE a, TE b) {

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

    cerr << "TECompare <RP_anno-TE.txt> <reads.gpd>" << endl;

    ifstream bedinf(argv[1]);
    vector<TE> te_vec = BuildTEVectorFromFile(bedinf);
    bedinf.close();

    ifstream gpdinf(argv[2]);

    while(gpdinf){

        string strInput;
        getline(gpdinf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            int max_overlap = 0;
            int max_index = -1;

            if(gpd.nblock == 1){

                TE psdo_te;
                psdo_te.chro = gpd.chro;
                psdo_te.start = gpd.vec_exon_start[0];
                psdo_te.end = gpd.vec_exon_end[0];
                psdo_te.strand = gpd.strand;
                psdo_te.name = gpd.gname;

                int locate = BinarySearchTE(psdo_te, te_vec);

                for(int j = -50; j <= 50; ++j){
                    int index = locate + j;
                    if(index >= 0 && index < te_vec.size()){
                        int overlap = Overlap(psdo_te, te_vec[index]);
                        if(overlap > max_overlap){
                            max_overlap = overlap;
                            max_index = index;
                        }
                    }
                }

                if(max_index >= 0){

                    int gpd_len = GPD_Length(gpd);
                    int te_len = te_vec[max_index].end - te_vec[max_index].start + 1;
                    double ratio1 = max_overlap / (gpd_len + 0.0);
                    double ratio2 = max_overlap / (te_len + 0.0);
                    cout << gpd.gname << '\t' << te_vec[max_index].name << '\t' << gpd.strand << '\t' << te_vec[max_index].strand << '\t' << ratio1 << '\t' << ratio2 << endl;
                }
                else  cout << gpd.gname << "\t*\t*\t*\t0\t0" << endl;
            }
            else  cout << gpd.gname << "\t*\t*\t*\t0\t0" << endl;
        }
    }

    gpdinf.close();

    return 0;
}
