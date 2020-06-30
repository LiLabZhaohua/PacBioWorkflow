#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <algorithm>
//Yu Sun, 2020-05-15
//This script filters GPD annotation file using PAS peaks
//Only the GPD record with PAS peak (up or down 30nt) supported will be output
//peaks.txt should have 3 columns: chr, strand, location
//Usage: PASFilterGPDByPeaks.cpp <peaks.txt> <assembly.gpd>

using namespace std;

//A split function
vector<string> SplitStringByDelim (string s, char delimitor){
  istringstream stringstream(s);   //read string as a stream
  string token;
  vector<string> tokens;
  while(getline(stringstream, token, delimitor))
  	tokens.push_back(token);
  return tokens;
}

struct GPD {
  string gene;
  string tx;
  string chr;
  string strand;
  int start;
  int end;
  int cstart;
  int cend;
  int nblock;
  vector<string> vec_exon_start;
  vector<string> vec_exon_end;
  GPD(vector<string> anno){
  	gene=anno[0];tx=anno[1];chr=anno[2];strand=anno[3];start=stoi(anno[4]);end=stoi(anno[5]);
  	cstart=stoi(anno[6]);cend=stoi(anno[7]);nblock=stoi(anno[8]);
  	vec_exon_start=SplitStringByDelim(anno[9],',');
  	vec_exon_end=SplitStringByDelim(anno[10],',');
  }
};

struct PEAK {
  string chr;
  string strand;
  int peakloc;
  PEAK(string c, string s, int p){
  	chr=c; strand=s; peakloc=p;
  }
};

bool comparePeaks (PEAK a, PEAK b){         //Sort two PEAK object in the default small->great order
  if (a.chr != b.chr){						//This is a internal criteria of string comparison
    return (a.chr < b.chr);
  }else{
    if (a.strand != b.strand){
      return (a.strand < b.strand);
    }else{
      return (a.peakloc < b.peakloc);
    }
  }
}

PEAK SearchCloestPeak (vector<PEAK> peak_database, PEAK curr_peak){
  int start=0;
  int end=peak_database.size()-1;
  int breakpoint;
  while (end - start >1){                //After this binary search loop, start and end will be consecutive numbers: [start, start+1]
  	breakpoint=int((start+end)/2);
  	if (comparePeaks(peak_database[breakpoint], curr_peak)){
  	  start=breakpoint;
  	}else{
  	  end=breakpoint; //+1?
  	}
  }
  //Measure distance of the two boundaries
  int final = end;  //initializing
  if ((peak_database[start].chr == curr_peak.chr) && (peak_database[end].chr == curr_peak.chr)){
    int dist_start=abs(peak_database[start].peakloc - curr_peak.peakloc);
    int dist_end=abs(peak_database[end].peakloc - curr_peak.peakloc);
    if (dist_start > dist_end){    //select closest
      final=end;
    }else{final = start;}}
  else if ((peak_database[start].chr != curr_peak.chr) && (peak_database[end].chr == curr_peak.chr)){
     final=end;}
  else if ((peak_database[start].chr == curr_peak.chr) && (peak_database[end].chr != curr_peak.chr)){
     final=start;}
  else if ((peak_database[start].chr != curr_peak.chr) && (peak_database[end].chr != curr_peak.chr)){
     final=end;}
  return peak_database[final];
}

int main (int argc, char *argv[]) {
  if (argc == 3){
    ifstream peaks(argv[1]);
    ifstream annotation(argv[2]);
  	
  	//Loop through the file to store peaks into a vector: peak_database
  	string peak_line;
  	vector<PEAK> peak_database;
	while(getline(peaks, peak_line)){            //Put getline into the while loop, and when reaching EOF, it will return false
		vector<string> peak_vector=SplitStringByDelim(peak_line,'\t');  //Split by tab, peack_vector[i]
  		PEAK current_peak(peak_vector[0],peak_vector[1],stoi(peak_vector[2]));         //stoi: converts string to int
		peak_database.push_back(current_peak);
	}
  	
  	//Sort the vector, for binary search
  	sort(peak_database.begin(), peak_database.end(), comparePeaks);    //Use the criteria defined by compareVector to sort the PEAK object vector
  	
  	//Loop through the transcript GPD file and compare each 3end with peak_database
  	string anno_line;
  	while(getline(annotation, anno_line)){
  		vector<string> anno_vector = SplitStringByDelim(anno_line,'\t');
  		GPD anno_gpd=GPD(anno_vector);
  		PEAK anno_3end=PEAK("chr1","+",0);
  		if (anno_gpd.strand == "+"){
  		    anno_3end=PEAK(anno_gpd.chr,anno_gpd.strand,anno_gpd.end);
  		}else{
  			anno_3end=PEAK(anno_gpd.chr,anno_gpd.strand,anno_gpd.start);
  		}
  		PEAK closest_peak = SearchCloestPeak(peak_database, anno_3end);
  		//cout << "closest_peak: " << closest_peak.chr << "\t"<<closest_peak.strand<<"\t"<<closest_peak.peakloc<< endl;
  		
  		//closest peak distance, up means the peak is before the current 3end
  		int dist_up=10000;
  		int dist_down=10000;
  		if ((anno_3end.chr == closest_peak.chr) && (anno_3end.strand == closest_peak.strand)){
  			if (anno_3end.strand == "+"){
  				if (anno_3end.peakloc > closest_peak.peakloc){
  					dist_up=anno_3end.peakloc - closest_peak.peakloc;
  				}else{
  					dist_down=closest_peak.peakloc - anno_3end.peakloc;
  				}
  			}else{
  				if (anno_3end.peakloc > closest_peak.peakloc){
  					dist_down=anno_3end.peakloc - closest_peak.peakloc;
  				}else{
  					dist_up=closest_peak.peakloc - anno_3end.peakloc;
  				}
  			}
  		}
  		//cout << "up and down: "<<dist_up<<", "<<dist_down<<endl;      //output this if you want to know the peak distance
  		if (dist_up < 30 || dist_down < 30) cout << anno_line << endl;
  	}
  	
  	peaks.close();
  	annotation.close();
  }else{
    cout << "Welcome to use PASFilterGPDByPeaks" << endl;
    cout << "Yu Sun, 2020-05-15" << endl;
    cout << "Usage: test <peaks.txt> <assembly.gpd>" << endl;
    cerr << "Please input 2 arguments !"<<endl;
  }
}
