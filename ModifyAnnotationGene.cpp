#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;


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


struct GENE {

    string name;
    string strand;
    string chro;
    int start;
    int end;
};


bool Overlap (GPD a, GPD b) {

    if(a.chro == b.chro){
        if(a.start > b.end || a.end < b.start)  return false;
        else  return true;
    }
    else  return false;
}


map<string,vector<GPD> > ChangeGeneNameMap (ifstream &inf) {

    map<string,vector<GPD> > change_gene_name_map;
    map<string,vector<GPD> >::iterator it;
    map<string,int> acc_map;
    map<string,int>::iterator it_acc;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);

            it = change_gene_name_map.find(gpd.gname);

            if(it == change_gene_name_map.end()){
                change_gene_name_map[gpd.gname].push_back(gpd);
                acc_map[gpd.gname] = 1;
            }
            else{
                bool overlap = false;
                int index = -1;
                for(int i = 0; i < (it->second).size(); ++i){
                    if(Overlap(gpd, (it->second)[i])){
                        overlap = true;
                        index = i;
                        break;
                    }
                }
                if(!overlap){
                    it_acc = acc_map.find(gpd.gname);
                    if(it_acc == acc_map.end())  cout << "Alert" << endl;
                    ++(it_acc->second);
                    gpd.gname = gpd.gname + "#" + lexical_cast<string>(it_acc->second);
                }
                else{
                    gpd.gname = (it->second)[index].gname;
                }
                (it->second).push_back(gpd);
            }
        }
    }

    return change_gene_name_map;
}


void OutputGPD (GPD gpd) {

    cout << gpd.gname << '\t';
    cout << gpd.tname << '\t';
    cout << gpd.chro << '\t';
    cout << gpd.strand << '\t';
    cout << gpd.start << '\t';
    cout << gpd.end << '\t';
    cout << gpd.cstart << '\t';
    cout << gpd.cend << '\t';
    cout << gpd.nblock << '\t';
    for(int i = 0; i < gpd.vec_exon_start.size(); ++i)  cout << gpd.vec_exon_start[i] << ',';
    cout << '\t';
    for(int i = 0; i < gpd.vec_exon_end.size(); ++i)  cout << gpd.vec_exon_end[i] << ',';
    cout << endl;
}








int main (int argc, char **argv) {

    ifstream annoinf(argv[1]);
    map<string,vector<GPD> > change_gene_name_map = ChangeGeneNameMap(annoinf);
    annoinf.close();

    map<string,vector<GPD> >::iterator it;

    for(it = change_gene_name_map.begin(); it != change_gene_name_map.end(); ++it){
        for(int i = 0; i < (it->second).size(); ++i)  OutputGPD((it->second)[i]);
    }

    return 0;
}
