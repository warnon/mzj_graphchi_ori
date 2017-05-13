
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
 * Template for GraphChi applications. To create a new application, duplicate
 * this template.
 */



#include <string>
#include <time.h>
#include <stdlib.h>
#include <map>
#include <set>
//#include <ext/hash_set>

#include "graphchi_basic_includes.hpp"
#include "util/active_analysis.hpp"
//#include "DAG.cpp"
//#include "util/labelanalysis.hpp"

using namespace graphchi;

/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */
struct Bilabel{
	int larger;
	int smaller;
	float weight;
	bool propagate;
	Bilabel(){
		//larger = smaller = -1;
		//propagate = true;	
	}

	Bilabel(int small, int large){
		smaller = small;
		larger = large;
		//propagate = true;
	}
	
	int& neighbor_label(vid_t myid, vid_t nbid){
		if(myid < nbid){
			return larger;
		}else{
			return smaller;
		}	
	}	
	int& my_label(vid_t myid, vid_t nbid){
		if(myid < nbid){
			return smaller;
		}else{
			return larger;
		}
	}

	bool labels_positive(){
		return  smaller >=0 && larger >=0;
	}
	bool can_propagate(){
		return propagate;
	}
	/*
	bool& propagate(){
		return propagate;
	}	
	*/
	void enable_propagate(){
		propagate = true;
	}
	void disable_propagate(){
		propagate = false;
	}
};

struct Vertexinfo{
	int label;
	bool inbfs;
	Vertexinfo(){
		//inbfs = false;
	}	
	Vertexinfo(int value){
		label = value;
		inbfs = false;
	}
	bool is_active(){
		return label < 0;
	}
	friend std::ostream& operator <<(std::ostream& out, Vertexinfo& vinfo){
		out << vinfo.label;	
		return out;
	}
};

bool operator < (const Vertexinfo& a, const Vertexinfo& b){
	return a.label < b.label;	
}

bool operator > (const Vertexinfo& a, const Vertexinfo& b){
	return a.label > b.label;	
}

bool operator != (const Vertexinfo& a, const Vertexinfo& b){
	return a.label != b.label;	
}

bool operator == (const Vertexinfo& a, const Vertexinfo& b){
	return a.label == b.label;	
}
//typedef int VertexDataType;
typedef Vertexinfo VertexDataType;
typedef Bilabel  EdgeDataType;

//unsigned single_source = 0;
bool converged = false;
//unsigned maxlevel = 100000;
bool scheduler = false;
bool weight_flag = false;

mutex lock;
static void parse(EdgeDataType& edata, const char* s){
	weight_flag = true;	
	edata.weight = atof(s);	
}

struct label_count{
	int label;
	vid_t count;
	//vid_t start_vid;
	label_count(int lb, vid_t ct){
		label = lb;
		count = ct;
	}	
	label_count(){
		label = -1;
		count = 0;
	}
};

FILE* vfout = NULL;
FILE* efout = NULL;
FILE* fp_interval = NULL;
std::vector<vid_t> start_vid;
//storing the mapping between block id with array index 
std::map<int, int> block_to_idx;

int get_new_id(int lb){
	std::map<int, int>::iterator miter = block_to_idx.find(lb); 
	if(miter == block_to_idx.end()){
		std::cout<<"unfound label: "<<lb<<std::endl;	
		assert(miter != block_to_idx.end());
	}
	int idx = miter->second; 
	int new_id = 0;
	lock.lock();	
	new_id = start_vid[idx];
	start_vid[idx]++;
	lock.unlock();
	return new_id;
}
/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */

struct initBFS : public GraphChiProgram<VertexDataType, EdgeDataType> {
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		if (gcontext.iteration == 0) {
			int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			VertexDataType vdata = vertex.get_data();
			vdata.label = -(vertex_id+1);
			vdata.inbfs = false;
			vertex.set_data(vdata);
			for(int id = 0; id < vertex.num_edges(); id++)
			{
				if(scheduler)
				{
					gcontext.scheduler->add_task(vertex.edge(id)->vertex_id());
				}
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = -(vertex_id+1);
				edata.enable_propagate();
				e->set_data(edata);
						
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			//assert(vertex.num_edges() ==4 );	

		} /*else {

		}*/
    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			if(!scheduler)
				converged = iteration > 0;
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		gcontext.set_last_iteration(iteration);
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};


double sample_rate = 0.001;
double sample_max = 0.1;
double multiply = 2.0;
double stop_ratio = 0.9;
int max_block_size = 1000;
int max_iterations = 10;

/*
float get_random(){
	return std::rand();	
}
*/
//pair<source_id, count>
std::vector< std::pair<vid_t, int> > sources;
//hash_set<std::pair<vid_t, int>> sources;

int get_index(vid_t sid){
	for(int i=0; i<(int)sources.size(); i++){
		if(sources[i].first == sid)
			return i;			
	}
		return -1;
}

/*
int get_block_size(vid_t sid){
		
}
*/

void insert_source(vid_t sid){
	std::pair<vid_t, int> sr(sid, 1);
	lock.lock();
	sources.push_back(sr);	
	//sources.insert(sr);	
	lock.unlock();
}

bool add_to_source(vid_t sid){
	int idx = get_index(sid);
	//only this round BFS 
	if(idx < 0){
		return false;
	}
	if(sources[idx].second < max_block_size ){
		lock.lock();	
		if(sources[idx].second >= max_block_size ){
			lock.unlock();	
			return false;
		}
		//lock.lock();	
		sources[idx].second++;
		lock.unlock();
		return true;
	}else{
		return false;
	}
}

int sum_block(){
	int sum = 0;
	for(int i=0; i<(int)sources.size(); i++){
		sum += sources[i].second;
	}
	return sum;
}

void show_sources(){
	for(int i=0; i<(int)sources.size(); i++){
		std::cout<<"source "<<sources[i].first<<"\t size "<<sources[i].second<<std::endl;				
	}
}

struct msBFS : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		//vertices explored are no longer considered	
		if((vertex.get_data()).label >= 0){
			assert((vertex.get_data()).inbfs);
			return;
		}

		if (gcontext.iteration == 0) {
			//never consider ioslated vertices when sampling
			if(vertex.num_edges() > 0){
				double prob = (double)rand()/RAND_MAX;
				if(prob <= sample_rate){//sampled	
					//vertex.set_data((int)vertex.id());
					VertexDataType vdata = vertex.get_data();
					vdata.label = (int)vertex.id();
					vdata.inbfs = true;
					vertex.set_data(vdata);

					insert_source(vertex.id());	
					for(int id = 0; id < vertex.num_edges(); id++)
					{
						if(scheduler)
						{
							gcontext.scheduler->add_task(vertex.edge(id)->vertex_id());
						}
						graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
						EdgeDataType edata = e->get_data();
						edata.my_label(vertex.id(), e->vertex_id()) = (int)vertex.id();
						e->set_data(edata);
						//std::cout<<"vid"<<vertex.edge(id)->vertexid
					}
				}
			}
		} else {
			/* Do computation */ 
			//hash_set<int> nblabels;	
			std::set<int> nblabels;	
			for(int i=0; i < (int)vertex.num_edges(); i++){
				graphchi_edge<EdgeDataType> * e = vertex.edge(i);	
				EdgeDataType edata = e->get_data();
				int nblabel = edata.neighbor_label(vertex.id(), e->vertex_id());
				if(nblabel > 0){
					nblabels.insert(nblabel);	
				}
				//e->set_data(edata);
			}

			if(nblabels.size() > 0){
				int random_idx = rand() % nblabels.size();					
				std::set<int>::iterator riter = nblabels.begin();				
				//hash_set<int>::iterator riter = nblabels.begin();				
				//int source = *(riter + random_idx);
				std::advance(riter , random_idx);
				int source = *riter;//std::advance(riter , random_idx);
				/* add the vertex to the small CC with cid source
					if the block size reached the limit, nothing happens
				*/
				//VertexDataType vdata = vertex.get_data();
				if(add_to_source((vid_t)source)){
					converged = false;
					VertexDataType vdata = vertex.get_data();
					vdata.label = source;
					vdata.inbfs = true;
					vertex.set_data(vdata);

					//vertex.set_data(source);	
					for(int i=0; i < vertex.num_edges(); i++) {
						graphchi_edge<EdgeDataType> * e = vertex.edge(i);	
						EdgeDataType edata = e->get_data();
						edata.my_label(vertex.id(), e->vertex_id()) = source;
						e->set_data(edata);
					}	
				}/*else{
				//	std::cout<<"add to source failed"<<std::endl;
					//assert(vdata.inbfs == false);
						
				}*/
					
				VertexDataType vdata = vertex.get_data();
				if(vdata.label >= 0)
					assert(vdata.inbfs == true);
			}	
		}
    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			/*
			if(!scheduler)
				converged = iteration > 0;
			*/
		srand((unsigned)time(NULL));
		converged = iteration > 0;	
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		/*
		if(!scheduler){
			if(converged)
			{
				logstream(LOG_INFO)<<"bfs program has converged"<<std::endl;
				gcontext.set_last_iteration(iteration);
			}
		}
		*/
		if(gcontext.iteration >= max_iterations || converged){
			logstream(LOG_INFO)<<"max number of iterations is reached, terminate now!"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}	
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

struct checkBFS : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		VertexDataType vdata = vertex.get_data();
		if(vdata.label < 0){
			assert(vdata.inbfs == false);
			return;
		}else{
			assert(vdata.inbfs == true);
		}
		return;

	}
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			/*
			if(!scheduler)
				converged = iteration > 0;
			*/
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		gcontext.set_last_iteration(iteration);
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

struct initWCC : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		VertexDataType vdata = vertex.get_data();
		if(vdata.label < 0){
			assert(vdata.inbfs == false);
			return;
		}

		if (gcontext.iteration == 0) {
			assert((vertex.get_data()).inbfs == true);
			//int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			//vertex.set_data(-vertex_id);
			//VertexDataType vdata = vertex.get_data();
			//vdata.label = vertex_id;
			//vdata.inbfs = false;
			//vertex.set_data(vdata);

			for(int id = 0; id < vertex.num_edges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				edata.propagate = false;
				//edata.propagate() = false;
				e->set_data(edata);
			}
			//assert(vertex.num_edges() ==4 );	

		
		} else {
			for(int id = 0; id < vertex.num_edges(); id++)
			{
			
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				//assert(edata.can_propagate() == false);
				assert(edata.propagate == false);
				//assert(edata.propagate == true);
			}

		}
    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			/*
			if(!scheduler)
				converged = iteration > 0;
			*/
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		//gcontext.set_last_iteration(iteration);
		if(iteration == 1)
			gcontext.set_last_iteration(iteration);
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

struct checkWCC : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		VertexDataType vdata = vertex.get_data();
		if(vdata.label < 0){
			assert(vdata.inbfs == false);
			return;
		}

		if (gcontext.iteration == 0) {
			assert((vertex.get_data()).inbfs == true);
			//int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			//vertex.set_data(-vertex_id);
			//VertexDataType vdata = vertex.get_data();
			//vdata.label = vertex_id;
			//vdata.inbfs = false;
			//vertex.set_data(vdata);

			for(int id = 0; id < vertex.num_edges(); id++)
			{
			
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				//assert(edata.can_propagate() == false);
				assert(edata.propagate == false);
				//assert(edata.propagate == true);
			}

		} /*else {

		}*/
    }
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			/*
			if(!scheduler)
				converged = iteration > 0;
			*/
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		gcontext.set_last_iteration(iteration);
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

struct WCC : public GraphChiProgram<VertexDataType, EdgeDataType> {
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		//don't consider verties in BFS
		if((vertex.get_data()).inbfs)
			return;

		if (gcontext.iteration == 0) {
			int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			//vertex.set_data(-vertex_id);
			VertexDataType vdata = vertex.get_data();
			vdata.label = vertex_id;
			//vdata.inbfs = false;
			vertex.set_data(vdata);

			for(int id = 0; id < vertex.num_edges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = vertex_id;
				//if(vdata.label >= 0)
				//edata.disable_propagate();
				e->set_data(edata);
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			//assert(vertex.num_edges() ==4 );	

		}else{
			VertexDataType vdata = vertex.get_data();
			int min_label = vdata.label;
			assert(min_label >= 0);
			for(int id = 0; id < vertex.num_edges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				if(edata.can_propagate()){
					//edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = vertex_id;
					int tmp =  edata.neighbor_label(vertex.id(), vertex.edge(id)->vertex_id());
					min_label = std::min(tmp, min_label);
					assert(tmp >= 0);
					//e->set_data(edata);
				}
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			if(vertex.num_edges()>0)	
				assert(min_label >= 0);	

			if(vdata.label > min_label){
				for(int id = 0; id < vertex.num_edges(); id++)
				{
					graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
					EdgeDataType edata = e->get_data();
					edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = min_label;
					e->set_data(edata);
				}
				vdata.label = min_label;
				vertex.set_data(vdata);
				converged = false;
			}
		}
	}
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			//if(!scheduler)
				converged = iteration > 0;
			
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		if(converged){	
			std::cout<<"No vertices are left, WCC is finished"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};


struct ReMap : public GraphChiProgram<VertexDataType, EdgeDataType> {
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		if(vertex.num_edges() == 0)
			return;
		if (gcontext.iteration == 0) {
			///int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			//vertex.set_data(-vertex_id);
			VertexDataType vdata = vertex.get_data();
			int new_id = get_new_id(vdata.label);
			//vdata.label = vertex_id;
			//vdata.inbfs = false;
			//vertex.set_data(vdata);

			for(int id = 0; id < vertex.num_edges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = new_id;
				//if(vdata.label >= 0)
				//edata.disable_propagate();
				e->set_data(edata);
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			//assert(vertex.num_edges() ==4 );	
			lock.lock();
			fprintf(vfout, "%d\t%d\n", (int)vertex.id(), new_id);	
			lock.unlock();

		}else{
			//VertexDataType vdata = vertex.get_data();
			//int min_label = vdata.label;
			for(int id = 0; id < vertex.num_outedges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.outedge(id);	
				EdgeDataType edata = e->get_data();
				vid_t my_id = edata.my_label(vertex.id(), e->vertex_id());	
				vid_t nb_id = edata.neighbor_label(vertex.id(), e->vertex_id());	
				assert(my_id != nb_id);
				if(!weight_flag){
					lock.lock();
					fprintf(efout, "%u\t%u\n", my_id, nb_id);
					lock.unlock();
				}else{
					lock.lock();
					fprintf(efout, "%u\t%u\t%.3f\n", my_id, nb_id, edata.weight);
					lock.unlock();
				}	
				/*	
				if(edata.can_propagate()){
					//edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = vertex_id;
					int tmp =  edata.neighbor_label(vertex.id(), vertex.edge(id)->vertex_id());
					min_label = std::min(tmp, min_label);
					//e->set_data(edata);
				}
				*/
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			
		}
	}
    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &gcontext) {
			
			//if(!scheduler)
				//converged = iteration > 0;
			
    }
    
    /**
     * Called after an iteration has finished.
     */
	void after_iteration(int iteration, graphchi_context &gcontext) {
		/*
		if(converged){	
			std::cout<<"WCC on left vertices is finished"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}
		*/
		if(iteration == 1)
			gcontext.set_last_iteration(iteration);
	}
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {        
    }
    
};

int main(int argc, const char ** argv) {
	//std::vector<std::pair<vid_t, vid_t> > range = DAGmain(argc, argv);	
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);
    
    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("ms-BFS");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 100000); // Number of iterations
    scheduler       	 = get_option_int("scheduler", false); // Whether to use selective scheduling
	max_block_size 		 = get_option_int("max_block", 100000);
	sample_rate			 = (double)get_option_float("sampling_rate", 0.001);
	int num_par		     = get_option_int("num_par", 5);
	max_iterations		 = get_option_int("max_iterations", 10);

   	//single_source = get_option_int("root", 0); 
    /* Detect the number of shards or preprocess an input to create them */
    int nshards          = convert_if_notexists<EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));
    
    /* Run */
    initBFS init_program;
    //bfs program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
	//save edges attributes after execution in in-memory mode
	engine.set_save_edgesfiles_after_inmemmode(true);
    engine.run(init_program, niters);
	/*
	graphchi_engine<VertexDataType, EdgeDataType> enginexx(filename, nshards, scheduler, m); 
	checkBFS check_program1;	
	enginexx.run(check_program1, niters);
	std::cout<<"check initBFS is finished!"<<std::endl;	
	*/
	int active_curr = 0;	
	int active_prev = active_vertices_count<VertexDataType>(filename);		

	msBFS msbfs_program;	
	int round = 0;
	int block_size_sum = 0;	
	int total_blocks = 0;
	double ratio = 0;
	while(true){
		std::cout<<"running msBFS round "<<round++<<"------------"<<std::endl;
    	//engine2.run(msbfs_program, niters);
    	engine.run(msbfs_program, niters);
		active_curr = active_vertices_count<VertexDataType>(filename);		
		ratio = (double)active_curr/active_prev;
		active_prev = active_curr;
		block_size_sum += sum_block();
		total_blocks += sources.size();
		if(ratio > stop_ratio){
			break;	
		}			
		sample_rate *= multiply;	
		if(sample_rate > sample_max){
			break;
		}	
		sources.clear();
	}

	//graphchi_engine<VertexDataType, EdgeDataType> enginex(filename, nshards, scheduler, m); 
	//checkBFS check_program;	
	//enginex.run(check_program, niters);

	//graphchi_engine<VertexDataType, EdgeDataType> enginewcc(filename, nshards, scheduler, m); 
	//checkWCC check_wcc_program;	

	if(active_curr > 0){
		//graphchi_engine<VertexDataType, EdgeDataType> engine3(filename, nshards, scheduler, m); 
		initWCC initwcc_program;	
		engine.run(initwcc_program, niters);
		/*
		enginewcc.run(check_wcc_program, niters);
		std::cout<<"check wcc init finished"<<std::endl;	
		*/
		//graphchi_engine<VertexDataType, EdgeDataType> engine4(filename, nshards, scheduler, m); 
		WCC wcc_program;	
		engine.run(wcc_program, niters);
		std::cout<<"WCC is finished"<<std::endl;	
	}

   /*analyze result*/
	//show_sources();
	/*
	std::cout<<"total blocks: "<<total_blocks<<"\t visited vertices "<<block_size_sum<<"\t active: "
			<<active_curr<<"\t total vertices: "<< get_num_vertices(filename)<<std::endl;
	std::cout<<"sample_rate "<<sample_rate<<"\t sample_max "<<sample_max
			<<"\t round: "<<round<<"\t a_ratio: "<<ratio<<"\tstop_ratio: "<<stop_ratio<<std::endl;
	*/
	//m.start_time("label-analysis");
	//analyze_labels<int>(filename, 40); 
	//analyze_labels<unsigned>(filename, 40); 
	std::vector<label_count> collect_result;
	count_labels<VertexDataType, label_count>(filename, 40, collect_result); 
	//std::cout<<"size: "<<collect_result.size()<<std::endl;
		
	//std::vector<int> partitions(nshards, 0);
	if(num_par == 0)
		num_par = nshards;
	std::vector<int> partitions(num_par, 0);
	std::vector<int> block_to_par(collect_result.size(), -1);
	//std::vector<vid_t> start_vid(collect_result.size(), 0);
	start_vid.resize(collect_result.size());
	
	std::cout<<"assigning blocks to partitions"<<std::endl;
	
	for(int i=0; i<(int)collect_result.size(); i++){
		int min_idx = 0;		
		for(int j=1; j<(int)partitions.size(); j++){
			if(partitions[j] < partitions[min_idx])
				min_idx = j;
		}	
		partitions[min_idx] += collect_result[i].count;				
		block_to_par[i] = min_idx;	
		block_to_idx.insert(std::pair<int,int>(collect_result[i].label, i));
	}						

	fp_interval = fopen((filename+".blogel.interval").c_str(), "w+");
	assert(fp_interval != NULL);
	//std::cout<<"assigning blocks to partitions finished!"<<std::endl;
	//calculate the prefix sum of all partitions and blocks		
	int prefix_tmp = 0;
	for(int i=0; i<(int)partitions.size(); i++){
		int tmp = partitions[i];	
		partitions[i] = prefix_tmp;
		prefix_tmp += tmp;
		//write the interval range to the file (interval is closed [])
		fprintf(fp_interval, "%d\t%d\n", partitions[i], prefix_tmp-1);
	}	
	fclose(fp_interval);
	//std::cout<<"assigning blocks to partitions finished!1"<<std::endl;
	//start_vid store the starting vertex id for each block;
	for(int i=0; i<(int)collect_result.size(); i++){
		start_vid[i] = partitions[block_to_par[i]];	
		partitions[block_to_par[i]] += collect_result[i].count;	
	}
	//std::cout<<"assigning blocks to partitions finished!2"<<std::endl;
	//vfout = fopen((filename+".bv").c_str(), "w+");
	vfout = fopen((filename+".blogel.vmap").c_str(), "w+");
	//efout = fopen((filename+".be").c_str(), "w+");
	efout = fopen((filename+".blogel").c_str(), "w+");
	assert(vfout != NULL && efout != NULL);
	fprintf(vfout, "# old_vid  new_vid\n");
	fprintf(efout, "# new_src  new_dst\n");
	fflush(vfout);
	fflush(efout);
	std::cout<<"ReMap is started"<<std::endl;
	//graphchi_engine<Vertexinfo, EdgeDataType> engine5(filename, nshards, scheduler, m);	
	ReMap remap;
	engine.run(remap, 2);
	fclose(vfout);
	fclose(efout);
	//std::cout<<"block number is "<<block_num<<std::endl;
    /* Report execution metrics */
	std::cout<<"total blocks: "<<total_blocks<<"\t visited vertices "<<block_size_sum<<"\t active: "
			<<active_curr<<"\t total vertices: "<< get_num_vertices(filename)<<std::endl;
	std::cout<<"sample_rate "<<sample_rate<<"\t sample_max "<<sample_max
			<<"\t round: "<<round<<"\t a_ratio: "<<ratio<<"\tstop_ratio: "<<stop_ratio<<std::endl;
    metrics_report(m);
    return 0;
}
