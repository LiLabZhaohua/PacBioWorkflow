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


pair<int,int> AddLen (string &cigar) {

    char ch[5000];
    int index = 0;
    int aln_len = 0;
    int S_num = 0;
    int left_S = 0;
    int right_S = 0;

    for(int i = 0; i < cigar.length(); ++i){

        if(cigar[i] > 'A' && cigar[i] < 'Z'){

            for(int j = index; j < 5000; ++j)  ch[j] = '\0';
            int number = atoi(ch);
            index = 0;

            if(cigar[i] == 'M' || cigar[i] == 'D' || cigar[i] == 'N')  aln_len += number;
            else if(cigar[i] == 'S' || cigar[i] == 'H'){
                ++S_num;
                if(S_num == 1)  left_S = number;
                else if(S_num == 2)  right_S = number;
                else  ;
            }
        }
        else{
            ch[index] = cigar[i];
            ++index;
        }
    }

    return pair<int,int>(-left_S, aln_len + right_S);
}





int main (int argc, char **argv) {

    cerr << "GetReadEnds <alignment.sam>" << endl;

    ifstream inf(argv[1]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            int begin_point, end_point;
            pair<int,int> add_len = AddLen(vec[5]);
            if(vec[1] == "0"){
                begin_point = lexical_cast<int>(vec[3]) + add_len.first;
                end_point = lexical_cast<int>(vec[3]) + add_len.second;
            }
            else if(vec[1] == "16"){
                begin_point = lexical_cast<int>(vec[3]) + add_len.second;
                end_point = lexical_cast<int>(vec[3]) + add_len.first;
            }
            else  ;

            cout << vec[0] << '\t' << vec[2] << '\t' << begin_point << '\t' << end_point << endl;
        }
    }

    inf.close();

    return 0;
}
