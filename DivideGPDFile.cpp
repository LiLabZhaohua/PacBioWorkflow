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

string chro_name[44] = {"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM",
"chr1_GL456210_random","chr1_GL456211_random","chr1_GL456212_random","chr1_GL456221_random","chr4_GL456216_random","chr4_GL456350_random","chr4_JH584292_random","chr4_JH584293_random",
"chr4_JH584294_random","chr4_JH584295_random","chr5_GL456354_random","chr5_JH584296_random","chr5_JH584297_random","chr5_JH584298_random","chr5_JH584299_random","chr7_GL456219_random",
"chrUn_GL456372","chrUn_GL456378","chrUn_GL456389","chrUn_GL456392","chrUn_JH584304","chrX_GL456233_random"};


int main (int argc, char **argv) {

    cerr << "DivideGPDFile <input.gpd> <directory>" << endl;
    cerr << "The input GPD file should be sorted according to gene name." << endl;

    ifstream inf(argv[1]);
    ofstream ouf;
    ofstream nouf("NAMELIST", ios::trunc);
    string gname;
    string chro;

    for(int i = 0; i < 44; ++i){
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

            if(vec[0] != gname){

                char gpdfile[100];
                chro = vec[2];
                sprintf(gpdfile, "%s/annotation_%s/%s.gpd", argv[2], chro.c_str(), vec[0].c_str());
                nouf << gpdfile << endl;
                gname = vec[0];
                ouf.close();
                ouf.open(gpdfile);
            }

            ouf << strInput << endl;
        }
    }

    inf.close(); ouf.close(); nouf.close();


    ifstream ninf("NAMELIST");

    while(ninf){

        string strInput;
        getline(ninf, strInput);

        if(strInput.length() > 0){

            int L = strInput.length();
            string bedfile = strInput.substr(0,L-4) + "_presort.bed";
            string bedfile_sort = strInput.substr(0,L-4) + "_sorted.bed";
            string bedfile_merge = strInput.substr(0,L-4) + "_merge.bed";
            char command[100];
            sprintf(command, "./code/GPDToBED %s %s", strInput.c_str(), bedfile.c_str());
            system(command);
            char cm_sort[100];
            sprintf(cm_sort, "sort -n -k 2,2 %s > %s", bedfile.c_str(), bedfile_sort.c_str());
            system(cm_sort);
            char cm_merge[100];
            //cout << bedfile_sort << endl;
            sprintf(cm_merge, "bedtools merge -i %s > %s", bedfile_sort.c_str(), bedfile_merge.c_str());
            system(cm_merge);
            char cm_rm1[100];
            sprintf(cm_rm1, "rm %s", bedfile.c_str());
            system(cm_rm1);

            string outputfile = strInput.substr(0,L-3) + "bed";

            ifstream inf1(strInput.c_str());
            string str;
            getline(inf1, str);
            vector<string> vec;
            split(vec, str, is_any_of("\t"));
            int num = 0;

            ifstream inf2(bedfile_merge.c_str());
            ofstream oufbed(outputfile.c_str());

            while(inf2){

                string strInput;
                getline(inf2, strInput);

                if(strInput.length() > 0){
                    ++num;
                    oufbed << strInput << '\t' << vec[0] << "_mergeexon_" << num << "\t0\t" << vec[3] << endl;
                }
            }

            inf1.close(); inf2.close(); oufbed.close();

            char cm_rm2[100];
            sprintf(cm_rm2, "rm %s", bedfile_sort.c_str());
            system(cm_rm2);
            char cm_rm3[100];
            sprintf(cm_rm3, "rm %s", bedfile_merge.c_str());
            system(cm_rm3);
            char cm_rm4[100];
            sprintf(cm_rm4, "rm %s", strInput.c_str());
            system(cm_rm4);
        }
    }

    ninf.close();

    return 0;
}
