#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <omp.h>
#include <cstdlib>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct PEAK {

    string chro;
    string strand;
    int pos;
};


vector<PEAK> BuildPeakVecFromFile (ifstream &inf) {

    vector<PEAK> pk_vec;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            PEAK pk;
            pk.chro = vec[0];
            pk.strand = vec[1];
            pk.pos = lexical_cast<int>(vec[2]);

            pk_vec.push_back(pk);
        }
    }

    return pk_vec;
}


bool CheckPeak (PEAK pk, char *bamfile) {

    char samfile[1000];
    sprintf(samfile, "%s%s%d.sam", pk.chro.c_str(), pk.strand.c_str(), pk.pos);

    char cmd_generate_samfile[1000];
    sprintf(cmd_generate_samfile, "samtools view %s %s:%d-%d > %s", bamfile, pk.chro.c_str(), pk.pos - 50, pk.pos + 50, samfile);
    system(cmd_generate_samfile);

    int fwd_count = 0;
    int rev_count = 0;
    ifstream inf(samfile);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[1] == "0")  ++fwd_count;
            else if(vec[1] == "16")  ++rev_count;
            else  ;
        }
    }

    char cmd_rm_samfile[1000];
    sprintf(cmd_rm_samfile, "rm %s", samfile);
    system(cmd_rm_samfile);

    if(fwd_count > rev_count && pk.strand == "+")  return true;
    else if(fwd_count < rev_count && pk.strand == "-")  return true;
    else  return false;
}




int main (int argc, char **argv) {

    ifstream inf(argv[1]);
    char *bamfile = argv[2];
    vector<PEAK> pk_vec = BuildPeakVecFromFile(inf);

    bool *check_array = new bool[pk_vec.size()];

    #pragma omp parallel for num_threads(5)
    for(int i = 0; i < pk_vec.size(); ++i){
        check_array[i] = CheckPeak(pk_vec[i], bamfile);
    }

    for(int i = 0; i < pk_vec.size(); ++i){
        if(check_array[i])  cout << pk_vec[i].chro << '\t' << pk_vec[i].strand << '\t' << pk_vec[i].pos << endl;
    }

    inf.close();

    return 0;
}
