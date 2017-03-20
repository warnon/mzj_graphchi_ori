#include<iostream>
#include<fstream>
//#include "api/chifilenames.hpp"

//using namespace graphchi{
struct Degree{
	int indegree;
	int outdegree;
};

std::string degree_file_name(std::string filename){
	return filename+"_degs.bin";
}

unsigned GetMaxDegreeVertex(std::string basefilename){
	std::string fname = degree_file_name(basefilename);
	FILE* fp = fopen(fname.c_str(), "r");	
	Degree degarray [1024];
	int len = 0;
	unsigned maxvid = 0;
	unsigned long product = 0;
	unsigned  count = 0;
	while((len = fread(&degarray, sizeof(Degree), 1024, fp)) != 0){
		for(int i=0; i<len; i++){
			if(product < (unsigned long)(degarray[i].indegree * degarray[i].outdegree) ){
				product = (unsigned long) degarray[i].indegree * degarray[i].outdegree;
				maxvid = count;	
			}
			count++;
		}
	}
	fclose(fp);
	return maxvid;
}

/*
int main(int argc, const char** argv){
	printf("maxvid = %u\n", GetMaxDegreeVertex(std::string(argv[1])));	
}
*/
//}
