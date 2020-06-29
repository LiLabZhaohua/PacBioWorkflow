#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <map>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct JUNCTION {

    string chro;
    int start;
    int end;
};


struct EXON {

    string chro;
    int start;
    int end;
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


void OutputJunction (ifstream &inf, ofstream &ouf) {

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);

            for(int i = 0; i <= gpd.nblock - 2; ++i)  ouf << gpd.gname << '\t' << gpd.chro << '\t' << (gpd.vec_exon_end[i] + 1) << '\t' << (gpd.vec_exon_start[i + 1] - 1) << endl;
        }
    }
}


int compare_junction_sort (const void *a, const void *b) {

    if( (*(JUNCTION*) a).chro > (*(JUNCTION*) b).chro )  return 1;
    else if ( (*(JUNCTION*) a).chro == (*(JUNCTION*) b).chro ){
        if( (*(JUNCTION*) a).start > (*(JUNCTION*) b).start )  return 1;
        else  return -1;
    }
    else  return -1;
}


int compare_junction_exon (JUNCTION a, EXON b) {

    if(a.chro > b.chro)  return 1;
    else if (a.chro == b.chro){
        if(a.start > b.start )  return 1;
        else  return -1;
    }
    else  return -1;
}


vector<JUNCTION> BuildJunctionVector (ifstream &inf) {

    vector<JUNCTION> junction_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            JUNCTION junction;
            junction.chro = vec[1];
            junction.start = lexical_cast<int>(vec[2]);
            junction.end = lexical_cast<int>(vec[3]);

            junction_vec.push_back(junction);
        }
    }

    qsort(&junction_vec[0], junction_vec.size(), sizeof(JUNCTION), compare_junction_sort);

    return junction_vec;
}


int BinarySearch (vector<JUNCTION> &junction_vec, EXON query) {

    int up = 0;
    int down = junction_vec.size() - 1;
    while(down - up > 1){
        int middle = (up + down) / 2;
        if(compare_junction_exon(junction_vec[middle], query) == 1)  down = middle;
        else  up = middle;
    }

    return up;
}


bool Covered (JUNCTION a, EXON b) {

    bool covered = false;

    if(a.chro == b.chro && a.start > b.start && a.end < b.end)  covered = true;

    return covered;
}






int main (int argc, char **argv) {

    cerr << "IntronRetention <RefSeq.gpd> <query.gpd>" << endl;

    ifstream refseqinf(argv[1]);
    ifstream inf(argv[2]);

    ofstream ouf("ALL_JUNCTIONS", ios::trunc);
    OutputJunction(refseqinf, ouf);
    char command_merge[1000];
    sprintf(command_merge, "sort -u ALL_JUNCTIONS > JUNCTIONS");
    system(command_merge);

    ifstream mergedinf("JUNCTIONS");
    vector<JUNCTION> junction_vec = BuildJunctionVector(mergedinf);


    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);

            for(int i = 0; i < gpd.nblock; ++i){

                EXON query;
                query.chro = gpd.chro;
                query.start = gpd.vec_exon_start[i];
                query.end = gpd.vec_exon_end[i];

                int locate = BinarySearch(junction_vec,query);
                bool covered = Covered(junction_vec[locate], query);

                int down_locate = locate;
                while(down_locate < junction_vec.size() - 1){
                    ++down_locate;
                    if(junction_vec[down_locate].chro != query.chro || junction_vec[down_locate].start > query.end + 10000)  break;
                    else{
                        covered = Covered(junction_vec[down_locate], query);
                    }
                }

                int up_locate = locate;
                while(up_locate > 0){
                    --up_locate;
                    if(junction_vec[up_locate].chro != query.chro || query.start - junction_vec[up_locate].end > 10000)  break;
                    else{
                        covered = Covered(junction_vec[down_locate], query);
                    }
                }

                if(covered)  cout << gpd.gname << '\t' << gpd.tname << '\t' << (i + 1) << endl;
            }
        }
    }

    refseqinf.close(); inf.close(); ouf.close();

    return 0;
}
