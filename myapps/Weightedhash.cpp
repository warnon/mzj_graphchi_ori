/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Simple pagerank implementation. Uses the basic vertex-based API for
 * demonstration purposes. A faster implementation uses the functional API,
 * "pagerank_functional".
 */

#include <string>
#include <fstream>
#include <cmath>

//#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"
#include "MSBFS.cpp"
using namespace graphchi;
 
#define THRESHOLD 1e-2    
#define RANDOMRESETPROB 0.15
/*
struct VertexInfo{
	vid_t vid;
	int deg;
};
*/
struct EdgeInfo{
	vid_t largelabel;
	vid_t smalllabel;   
	float weight;
	vid_t& my_label(vid_t myid, vid_t nbid){
		if(myid < nbid) return smalllabel;	
		else return largelabel;
	}
	vid_t& nb_label(vid_t myid, vid_t nbid){
		if(myid < nbid) return largelabel;	
		else return smalllabel;
	}
	EdgeInfo(){
		largelabel = smalllabel = (vid_t)-1;	
	}
};

typedef vid_t VType;
typedef EdgeInfo EType;

bool flag_weight = false;
//static void VARIABLE_IS_NOT_USED parse<EdgeDataType>(EdgeDataType& edata, const char* s){
static void  parse(EdgeInfo& edata, const char* s){
	if(!flag_weight){
		flag_weight = true;
		std::cout<<"edata with weight"<<std::endl;
	}
	edata.weight = atof(s);	
}
/*
bool sortFunc(const VertexInfo& v1, const VertexInfo& v2){
	if(v1.deg < v2.deg){
		return false;	
	}else if(v1.deg > v2.deg){
		return true;	
	}else{
		if(v1.vid < v2.vid)	return true;	
		else return false;
	}	
}
*/

//mutex lock;
FILE* fp_list = NULL;
//FILE* fp_vt = NULL;
//float epsilon = 0.001;
//int* array = NULL;
size_t num_vertices = 0;

//vid_t start_id = 1;
//vid_t start_id_vld = 1;
size_t num_edges = 0;
//vid_t num_vertices = 0;
//bool flag = false;//denote whether the smallest vid is 0 or 1
//bool num_tasks_print = false;

//std::vector<vid_t> idmap;
//std::vector<vid_t> vnumbers;
vid_t* idmap = NULL;
std::string pfilename = ""; 

void fixline(char* str){
	if(str == NULL) return;
	int len = strlen(str);
	if(str[len-1] == '\n' )
		str[len-1] = 0;	
}

void initIdMap(){
	//idmap.resize(num_vertices);
	assert(num_vertices > 0);
	idmap = (vid_t*)malloc(sizeof(vid_t)*num_vertices);
	memset(idmap, 0, sizeof(vid_t)*num_vertices);

	FILE* fpp = fopen(pfilename.c_str(), "r");			
	assert(fpp != NULL);
	char buffer[1024];
	char delimiter[] = "\t ,";
	//int count = 0;
	//int sum = 0;
	while(fgets(buffer, 1024, fpp) != NULL){
		fixline(buffer);	
		if(buffer[0] == '#'|| buffer[0] == '%') continue;
		char* str = strtok(buffer, delimiter);		
		int old_vid = atoi(str);	
		
		str = strtok(NULL, delimiter); 
		int new_vid = atoi(str);

		assert(old_vid < num_vertices && new_vid < num_vertices);	
		assert(idmap[old_vid] == 0);
		idmap[old_vid] = (vid_t)new_vid;	
	}

	/*
	for(int i=0; i<(int)vnumbers.size(); i++){
		std::cout<<i<<"-th partition size="<<vnumbers[i]<<std::endl;
	}
	
	for(int i=0; i<(int)vnumbers.size(); i++){
		int tmp = sum;			
		sum += vnumbers[i];
		vnumbers[i] = tmp;
	}	
	for(int i=0; i<vnumbers.size(); i++){
		std::cout<<i<<"-th partition start vid="<<vnumbers[i]<<std::endl;
	}
	*/
	fclose(fpp);
	fpp = NULL;
}

vid_t getNewId(vid_t vid){
	assert(vid < num_vertices);
	return idmap[vid];
}

void freeMem(){
	free(idmap);	
}
/*
int getPId(vid_t vid){
	assert(vid < idmap.size());
	return idmap[vid];
}
vid_t getNewId(int partid){
	assert(partid < (int)vnumbers.size());
	vid_t newid = 0;
	for(int i=0; i<vnumbers.size(); i++){
		if(i == partid){
			lock.lock();
			newid = vnumbers[i]++;	
			lock.unlock();
		}		
	}
	return newid;
}
*/

struct ConvertProgram : public GraphChiProgram<VType, EType> {
   	bool converged;
	bool interval_converged; 
	mutex lock;
	/*
	std::vector<VType> degvector;
	std::vector<vid_t> idmap;
	*/
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &info) {
			/*
			if(iteration == 0)
				assert(NULL != array);
			*/
			converged = iteration > 0;
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
		if(iteration == 0){

		}else if(iteration == 1){
			fflush(fp_list);
			ginfo.set_last_iteration(iteration);	
		}
    }
    
    /**
      * Pagerank update function.
	  */
	void update(graphchi_vertex<VType, EType> &v, graphchi_context &ginfo) {
		//array[v.id()]++;		
		if(v.num_edges() == 0)	return;
		if (ginfo.iteration == 0) {
			//int partid = getPId(v.id());	
			vid_t newid = getNewId(v.id()); 	
			v.set_data(newid);
			for(int i=0; i<v.num_edges(); i++){
				graphchi_edge<EType> * edge = v.edge(i);
				EType edata = edge->get_data();
				edata.my_label(v.id(), edge->vertex_id()) = newid;
				edge->set_data(edata);
			}	
		} else if(ginfo.iteration == 1){
			/*
			if(v.id() == 0){
				fprintf(fp_list, "%u %u\n", num_vertices, num_edges);	
			}
			*/
			if(v.num_outedges() > 0){	
				vid_t mylabel = v.get_data();
				for(int i=0; i<v.num_outedges(); i++){
					graphchi_edge<EType> * edge = v.outedge(i);
					EType edata = edge->get_data();
					vid_t nblabel = edata.nb_label(v.id(), edge->vertex_id());
					//vid_t nb_id = edge->vertex_id();
					assert(mylabel != nblabel);
					if(!flag_weight){
						lock.lock();
						fprintf(fp_list, "%u\t%u\n", mylabel, nblabel);		
						lock.unlock();
					}else{
						lock.lock();
						fprintf(fp_list, "%u\t%u\t%.3f\n", mylabel, nblabel, edata.weight);		
						lock.unlock();
					}
					//edge->set_data(edata);	
				}
			}/*else{
				fprintf(fp_list, "\n");
			}*/
		}
	}
};


int main(int argc, const char ** argv) {

	//msbfsMain(argc, argv); 
	//std::cout<<"===============================MSBFS has finished================================="<<std::endl;
	
    graphchi_init(argc, argv);
    metrics m("Weightedhash-partitioning");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 100000);
    bool scheduler          = get_option_int("scheduler", false);;                    // Non-dynamic version of pagerank.
	//std::string uwfilename  = get_option_string("file");//unweighed file name 
    //int ntop                = get_option_int("top", 50);
    //epsilon	 				= get_option_float("epsilon", 0.001);
	//pfilename   			= get_option_string("pfile");
	//int execthreads		    = get_option_int("execthreads", 1);
	//num_tasks_print			= get_option_int("print", false);
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EType>(filename, get_option_string("nshards", "auto"));
	assert(0 != nshards);
	num_vertices = get_num_vertices(filename);
	assert(num_vertices > 0);

	std::string orig_filename = filename.substr(0, filename.find("_weight")); 	
	orig_filename += ".txt";
	pfilename = orig_filename + ".hash.vmap";	
	initIdMap();
	//array = (int*)malloc(sizeof(int)*num_vertices);
	//memset(array, 0, sizeof(int)*num_vertices);
	//return 0;
    /* Run */
    graphchi_engine<VType, EType> engine(filename, nshards, scheduler, m); 
	//engine.set_exec_threads(1);
    //engine.set_modifies_inedges(false); // Improves I/O performance.
   	fp_list = fopen((filename+".hash").c_str(), "w+");	 
   	//fp_vt = fopen((filename+".sd.vt").c_str(), "w+");	 
	//assert(fp_list != NULL);
	//assert(fp_vt != NULL);
	//fprintf(fp_vt, "%new_vid, old_vid, degree\n");
    bool inmemmode = false;//engine.num_vertices() * sizeof(EType) < (size_t)engine.get_membudget_mb() * 1024L * 1024L;
	std::cout<<"===============================ReMap is started================================="<<std::endl;
    if (inmemmode) {
		/*
        logstream(LOG_INFO) << "Running Pagerank by holding vertices in-memory mode!" << std::endl;
        engine.set_modifies_outedges(false);
        engine.set_disable_outedges(true);
        engine.set_only_adjacency(true);
        SortProgramInmem program(engine.num_vertices());
        engine.run(program, niters);
		*/
    } else {
        ConvertProgram program;
        engine.run(program, niters);
    }
   	fclose(fp_list); 
	//fclose(fp_vt);
	//fp_vt = NULL;
	fp_list = NULL;
    /* Output top ranked vertices */
	/*
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
   	*/ 
	freeMem();


    metrics_report(m);    
	//copy unweighted file interval to weighted file
	//std::cout<<"ouput interval file to dst"<<std::endl;
	FILE* srcfp_interval = fopen((orig_filename+".hash.interval").c_str(), "r");
	assert(srcfp_interval != NULL);
	FILE* dstfp_interval = fopen((filename+".hash.interval").c_str(), "w+");
	assert(dstfp_interval != NULL);
	char tmp[1024];
	std::cout<<"output file is "<<filename+".hash.interval"<<std::endl;
	
	//while(int bytes = read(srcfp_interval, tmp, 1024) > 0){	
	while(int bytes = fgets(tmp, 1024, srcfp_interval) > 0){	
		//std::cout<<tmp<<std::endl;		
		//write(tmp, sizeof(char), bytes, dstfp_interval);	
		//write(dstfp_interval, tmp, bytes);	
		fputs(tmp, dstfp_interval);
		//std::cout<<bytes<<" has read"<<std::endl;
	}
	fclose(srcfp_interval);
	fclose(dstfp_interval);
	
	//std::cout<<"end ouput interval file to dst"<<std::endl;
    //metrics_report(m);    
	/*
	size_t total_updates = 0;
	for(int i=0; i<num_vertices; i++){
		total_updates += array[i];	
	}
	std::cout<<"num_vertices: "<<num_vertices<<"\ntotal updates: "<<total_updates<<"\naverage updates per vertex: "
		<<total_updates/num_vertices<<std::endl;
	*/
    return 0;
}

