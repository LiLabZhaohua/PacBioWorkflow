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

string chro_name[22] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM"};


int main (int argc, char **argv) {

    cerr << "DivideBEDFile <input.bed> <directory>" << endl;

    ifstream inf(argv[1]);
    int record_num = 0;

    for(int i = 0; i < 22; ++i){
        char cmd_mkdir[1000];
        sprintf(cmd_mkdir, "mkdir %s/annotation_%s", argv[2], chro_name[i].c_str());
        system(cmd_mkdir);
    }
        

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            ++record_num;

            char outputfile[1000];
            sprintf(outputfile, "%s/annotation_%s/%d.bed", argv[2], vec[0].c_str(), record_num);
            //sprintf(outputfile, "%s/%d.bed", argv[2], record_num);

            ofstream ouf(outputfile, ios::trunc);
            ouf << strInput << endl;
            ouf.close();
        }
    }

    inf.close();

    return 0;
}
