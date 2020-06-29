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


map<string,vector<GPD> > BuildGeneMap (ifstream &inf) {

    map<string,vector<GPD> > gene_map;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            GPD gpd;
            BuildGPDFromString(strInput, gpd);
            gene_map[gpd.gname].push_back(gpd);
        }
    }

    return gene_map;
}


void MakeBEDFile (vector<GPD> &vec) {

    ofstream ouf("pre-intervals.bed", ios::trunc);

    for(int i = 0; i < vec.size(); ++i){
        for(int j = 0; j < vec[i].nblock; ++j){
            ouf << vec[i].chro << '\t' << vec[i].vec_exon_start[j] << '\t' << vec[i].vec_exon_end[j] << "\t*\t*\t" << vec[i].strand << endl;
        }
    }

    string cmd_sort = "sort -k 1,1 -k 2,2n pre-intervals.bed > pre-intervals-sorted.bed";
    system(cmd_sort.c_str());

    string cmd_merge = "bedtools merge -i pre-intervals-sorted.bed > intervals.bed";
    system(cmd_merge.c_str());
}


pair<int,double> StatDepth (ifstream &inf) {

    int essential_len = 0;
    int max_depth = 0;
    int total_depth = 0;

    while(inf){

        string strInput;
        getline(inf, strInput);

        if(strInput.length() > 0){

            vector<string> vec;
            split(vec, strInput, is_any_of("\t"));
            int current_depth = lexical_cast<int>(vec[2]);

            ++essential_len;
            total_depth += current_depth;
            if(current_depth > max_depth)  max_depth = current_depth;
        }
    }

    double ratio;
    if(max_depth * essential_len == 0)  ratio = 0.0;
    else  ratio = total_depth / (max_depth * essential_len + 0.0);

    return pair<int,double>(max_depth, ratio);
}




int main (int argc, char **argv) {

    ifstream annoinf(argv[1]);
    map<string,vector<GPD> > gene_map = BuildGeneMap(annoinf);
    annoinf.close();

    map<string,vector<GPD> >::iterator it;

    for(it = gene_map.begin(); it != gene_map.end(); ++it){

        MakeBEDFile(it->second);

        ifstream bedinf("intervals.bed");

        while(bedinf){

            string strInput;
            getline(bedinf, strInput);

            if(strInput.length() > 0){

                vector<string> vec;
                split(vec, strInput, is_any_of("\t"));

                string chro = vec[0];
                int start = lexical_cast<int>(vec[1]);
                int end = lexical_cast<int>(vec[2]);

                char cmd_get_depth[1000];
                sprintf(cmd_get_depth, "samtools depth -r %s:%d-%d %s >> depth.txt", chro.c_str(), start, end, argv[2]);
                system(cmd_get_depth);
            }
        }


        ifstream depthinf("depth.txt");
        pair<int,double> stat_rslt = StatDepth(depthinf);
        cout << (it->first) << '\t' << stat_rslt.first << '\t' << stat_rslt.second << endl;
        depthinf.close();

        //break;

        string cmd_rm = "rm pre-intervals.bed pre-intervals-sorted.bed intervals.bed depth.txt";
        system(cmd_rm.c_str());
    }


    annoinf.close();

    return 0;
}
