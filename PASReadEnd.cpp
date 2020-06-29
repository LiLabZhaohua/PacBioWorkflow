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


bool CigarContainS (string &cigar) {

    bool output = false;
    for(int i = 0; i < cigar.length(); ++i){
        if(cigar[i] == 'S'){
            output = true;
            break;
        }
    }

    return output;
}


int AddLen (string &cigar) {

    char ch[5000];
    int index = 0;
    int aln_len = 0;

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



int main (int argc, char **argv) {

    ifstream inf(argv[1]);

    while(inf){

        string strInput;
        getline(inf, strInput);
        if(strInput.length() == 0)  break;

        if(strInput[0] == '@')  continue;

        vector<string> vec;
        split(vec, strInput, is_any_of("\t"));
        string cigar = vec[5];
        if(CigarContainS(cigar))  continue;

        string chro = vec[2];
        string flag = vec[1];
        int pos = lexical_cast<int>(vec[3]);

        if(flag == "0"){
            pos = pos + AddLen(cigar);
            cout << chro << "\t+\t" << pos << endl;
        }
        else if(flag == "16"){
            cout << chro << "\t-\t" << pos << endl;
        }
    }

    inf.close();

    return 0;
}
