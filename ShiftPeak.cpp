#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#define RANGE 100

using namespace std;
using namespace boost;


int AlignLen (string &cigar) {

    char ch[5000];
    int index = 0;
    int aln_len = 0;
    int S_num = 0;

    for(int i = 0; i < cigar.length(); ++i){

        if(cigar[i] > 'A' && cigar[i] < 'Z'){

            for(int j = index; j < 5000; ++j)  ch[j] = '\0';
            int number = atoi(ch);
            index = 0;

            if(cigar[i] == 'M' || cigar[i] == 'D' || cigar[i] == 'N')  aln_len += number;
        }
        else{
            ch[index] = cigar[i];
            ++index;
        }
    }

    return aln_len;
}


int FindArgMax (string str, int start, int end) {

    ifstream inf(str.c_str());
    map<int,int> freq_map;
    map<int,int>::iterator it;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            int pos = lexical_cast<int>(vec[3]);
            if(vec[1] == "16")  pos = pos + AlignLen(vec[5]) - 1; // For cage
            //if(vec[1] == "0")  pos = pos + AlignLen(vec[5]) - 1; // For PAS

            it = freq_map.find(pos);
            if(it == freq_map.end())  freq_map.insert(pair<int,int>(pos,1));
            else  ++(it->second);
        }
    }

    inf.close();


    int max = 0;
    int arg_max = 0;

    for(it = freq_map.begin(); it != freq_map.end(); ++it){
        if(it->first > start - RANGE && it->first < end + RANGE && it->second > max){
            max = it->second;
            arg_max = it->first;
        }
    }

    return arg_max;
}


int Shift (char *bam_file, string chro, int start, int end) {

    char command[1000];
    sprintf(command, "samtools view %s %s:%d-%d > SAMFILE", bam_file, chro.c_str(), start - RANGE, end + RANGE);
    system(command);

    int shift_peak = FindArgMax("SAMFILE", start, end);

    string rm_file = "rm SAMFILE";
    system(rm_file.c_str());

    return shift_peak;
}





int main (int argc, char **argv) {

    cerr << "ShiftPeak <ASPeak_results.txt> <sorted_bam_file> <output_file>" << endl;

    ifstream inf(argv[1]);
    char *bam_file = argv[2];
    ofstream ouf(argv[3], ios::trunc);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            string chro = vec[0];
            int start = lexical_cast<int>(vec[7]);
            int end = lexical_cast<int>(vec[8]);

            int shift_peak = Shift(bam_file, chro, start, end);

            ouf << vec[0] << '\t' << vec[1] << '\t' << vec[2] << '\t' << shift_peak << endl;
        }
    }


    inf.close(); ouf.close();

    return 0;
}







/*int main (int argc, char **argv) {

    //cerr << "ShiftPeak <cage_peak_maxdepth.txt> <sorted_bam_file> <output_file>" << endl;
    cerr << "ShiftPeak <ASPeak_results.txt> <sorted_bam_file> <output_file>" << endl;

    ifstream inf(argv[1]);
    char *bam_file = argv[2];
    ofstream ouf(argv[3], ios::trunc);
    int line = 0;
    string chro;
    int gene_start, gene_end;

    while(inf){

        string strInput;
        getline(inf, strInput);
        ++line;

        if(strInput.length() > 0){

            if(strInput[0] == '@'){

                ouf << strInput << endl;
                line = 0;
                vector<string> vec;
                split(vec, strInput, is_any_of("\t"));
                chro = vec[1];
                gene_start = lexical_cast<int>(vec[3]);
                gene_end = lexical_cast<int>(vec[4]);
            }
            else if(line == 1)  ouf << strInput << endl;
            else{

                vector<string> vec;
                split(vec, strInput, is_any_of("\t"));
                int pos = lexical_cast<int>(vec[0]);
                //cout << chro << '\t' << pos << endl;

                int shift_peak = Shift(bam_file, chro, gene_start, gene_end, pos);

                ouf << shift_peak << '\t' << vec[1] << '\t' << vec[2] << endl;
            }
        }
    }


    inf.close(); ouf.close();

    return 0;
}*/
