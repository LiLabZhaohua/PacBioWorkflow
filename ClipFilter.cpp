#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


bool ShortClip (string &cigar) {

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

    if(left_S < 100 && right_S < 100)  return true;
    else  return false;
}


int main (int argc, char **argv) {

    ifstream inf(argv[1]);

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            string cigar = vec[5];
            bool short_clip = ShortClip(cigar);
            if(short_clip)  cout << strInput << endl;
        }
    }

    inf.close();

    return 0;
}
