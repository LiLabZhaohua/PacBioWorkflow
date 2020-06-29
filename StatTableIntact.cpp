#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


struct STAT {

    string gname;
    int intact_num;
    int compatible_num;
};


map<string,STAT> BuildIsoformStatMap (ifstream &inf) {

    map<string,STAT> stat_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            string tname = vec[1];
            stat_map[tname].gname = vec[0];
            stat_map[tname].intact_num = 0;
            stat_map[tname].compatible_num = 0;
        }
    }

    return stat_map;
}


int main (int argc, char **argv) {

    cerr << "StatTableIntact <assembly_final.gpd> <intact_reads_subisoforms_*.txt> <map_to_isoform_*.txt>" << endl;

    ifstream annoinf(argv[1]);
    ifstream intactinf(argv[2]);
    ifstream compatibleinf(argv[3]);

    map<string,STAT> stat_map = BuildIsoformStatMap(annoinf);
    map<string,STAT>::iterator it;


    while(intactinf){

        string strInput;
        getline(intactinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            it = stat_map.find(vec[2]);
            if(it == stat_map.end()){
                cout << "Alert" << endl;
            }
            else{
                ++(it->second).intact_num;
            }
        }
    }


    while(compatibleinf){

        string strInput;
        getline(compatibleinf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));

            if(vec[2] == "Compatible"){

                it = stat_map.find(vec[1]);
                if(it == stat_map.end())  cout << "Alert" << endl;
                else{
                    ++(it->second).compatible_num;
                }
            }
        }
    }


    for(it = stat_map.begin(); it != stat_map.end(); ++it){
        cout << (it->second).gname << '\t' << (it->first) << '\t' << (it->second).intact_num << '\t' << ((it->second).compatible_num - (it->second).intact_num) << endl;
    }

    annoinf.close(); intactinf.close(); compatibleinf.close();

    return 0;
}
