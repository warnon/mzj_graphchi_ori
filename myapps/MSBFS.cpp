
#include <string>
#include <time.h>
#include <stdlib.h>
#include <map>
#include <set>
//#include <ext/hash_set>

#include "graphchi_basic_includes.hpp"
#include "util/active_analysis.hpp"
#include "util/toplist.hpp"
#include "DAG.cpp"
//#include "util/labelanalysis.hpp"

using namespace graphchi;
int MAX_LEVEL = 10000000;
/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */
struct Bilabel{
	int larger;
	int smaller;
	float weight;
	bool propagate;
	//for BFS level constraint
	int level;
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
	int indeg;
	Vertexinfo(){
		//inbfs = false;
	}	
	Vertexinfo(int value){
		label = value;
		inbfs = false;
		indeg = 0;
	}

	bool is_active(){
		//label > 0 && inbfs == false means isolated vertices 
		return label < 0 && !inbfs;
	}

	int get_degree() const {
		return indeg;
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
int NumRoots = 0;
//int NumLevels = 0;
std::set<vid_t> roots;
std::map<vid_t, int> sizecount;

//int numPartitions = 0;

static void parse(EdgeDataType& edata, const char* s){
	if(!weight_flag)
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
std::vector<vid_t> start_vid;
//storing the mapping between block id with array index 
std::map<int, int> block_to_idx;

/*
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
*/

/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */

struct initBFS : public GraphChiProgram<VertexDataType, EdgeDataType> {
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		/*
		if(vertex.num_edges() == 0)
			return;
		*/
		if (gcontext.iteration == 0) {
			int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			VertexDataType vdata = vertex.get_data();
			vdata.label = -(vertex_id+1);
			assert(vdata.label < 0);
			//ignore isolated vertices
			if(vertex.num_edges() == 0) vdata.inbfs = true;
			else vdata.inbfs = false;
			vdata.indeg = vertex.num_inedges();
			//vdata.level = MAX_LEVEL;
			vertex.set_data(vdata);
			for(int id = 0; id < vertex.num_edges(); id++)
			{
				/*
				if(scheduler)
				{
					gcontext.scheduler->add_task(vertex.edge(id)->vertex_id());
				}
				*/
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = -(vertex_id+1);
				edata.enable_propagate();
				//initialize level on outedges
				if( id >= vertex.num_inedges()){
					edata.level = MAX_LEVEL;
				}
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

bool isSources(vid_t vid){
	return roots.find(vid) != roots.end();	
}

int MaxDeg = 0;
void increaseSize(vid_t vid){
	std::map<vid_t, int>::iterator it = sizecount.find(vid);	
	assert(it != sizecount.end());		
	lock.lock();
	it->second++;
	lock.unlock();
}

struct msBFS : public GraphChiProgram<VertexDataType, EdgeDataType> {
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		if(vertex.num_edges() == 0)
			return;
		//assert(MaxDeg >= (int)vertex.num_inedges());
		//vertices explored are no longer considered	
		if((vertex.get_data()).label >= 0){
			assert((vertex.get_data()).inbfs);
			return;
		}
			
		if (gcontext.iteration == 0) {
			if(isSources(vertex.id())){
				VertexDataType vdata = vertex.get_data();
				vdata.label = (int)vertex.id();
				vdata.inbfs = true;
				vertex.set_data(vdata);
				
				for(int id = 0; id < vertex.num_edges(); id++)
				{
					graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
					EdgeDataType edata = e->get_data();
					edata.my_label(vertex.id(), e->vertex_id()) = (int)vertex.id();
					//write level info to outedges
					if( id >= (int)vertex.num_inedges()){
						edata.level = 0;
						if(scheduler)
						{
							gcontext.scheduler->add_task(vertex.edge(id)->vertex_id());
						}
					}
					e->set_data(edata);
				}
				increaseSize(vertex.id());	
			}
			//never consider ioslated vertices when sampling
			/*
			if(vertex.num_edges() > 0){
				double prob = (double)rand()/RAND_MAX;
				if(prob <= sample_rate){//sampled	
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
					}
				}
			}
			*/
		} else {
			/* Do computation */ 
			//hash_set<int> nblabels;	
			//std::set<int> nblabels;	
			int minlevel = MAX_LEVEL;
			int minlabel = 0;
			for(int i=0; i < (int)vertex.num_inedges(); i++){
				graphchi_edge<EdgeDataType> * e = vertex.inedge(i);	
				EdgeDataType edata = e->get_data();
				int nblabel = edata.neighbor_label(vertex.id(), e->vertex_id());
				if(isSources(nblabel) && edata.level+1 < minlevel){	
					minlevel = edata.level+1;			
					//minlabel = edata.neighbor_label(vertex.id(), e->vertex_id());
					minlabel = nblabel;//edata.neighbor_label(vertex.id(), e->vertex_id());
				}
				/*
				int nblabel = edata.neighbor_label(vertex.id(), e->vertex_id());
				if(nblabel > 0){
					nblabels.insert(nblabel);	
				}
				*/
				//e->set_data(edata);
			}

			if(minlevel < MAX_LEVEL){
				converged = false;
				VertexDataType vdata = vertex.get_data();
				vdata.label = minlabel;
				vdata.inbfs = true;
				vertex.set_data(vdata);
				assert(isSources((vid_t)minlabel));
				for(int i=0; i < vertex.num_edges(); i++) {
					graphchi_edge<EdgeDataType> * e = vertex.edge(i);	
					EdgeDataType edata = e->get_data();
					edata.my_label(vertex.id(), e->vertex_id()) = minlabel;
					//write the level to outedges
					if(i >= (int)vertex.num_inedges()){
						edata.level = minlevel;
						if(scheduler)
						{
							gcontext.scheduler->add_task(e->vertex_id());
						}
					}
					e->set_data(edata);
				}	

				increaseSize((vid_t)minlabel);	
				/*
				VertexDataType vdata = vertex.get_data();
				if(vdata.label >= 0)
					assert(vdata.inbfs == true);
				*/
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
		//srand((unsigned)time(NULL));
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
		//if(gcontext.iteration >= max_iterations || converged){
		if(converged){
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
		if(vertex.num_edges() == 0)
			return;
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
		if(vertex.num_edges() == 0)
			return;
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
					int tmp =  edata.neighbor_label(vertex.id(), e->vertex_id());
					min_label = std::min(tmp, min_label);
					assert(tmp >= 0);
					//e->set_data(edata);
				}
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			/*
			if(vertex.num_edges()>0)	
				assert(min_label >= 0);	
			*/
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


//number of total blocks
int Size = 0;
std::map<int, int> label2idx;
//std::vector< std::vector<int> > adjmatrix;
int** adjmatrix;

int totalvertices = 0;//total number of non-isolated vertices
std::vector< bool > merged;

std::vector<label_count> collect_result;
std::vector<std::vector<int> > partitions;	
std::vector<int> partition_size;
//std::vector< std::set<int> > nbblocks; 
//number of cut-edges of each bfs tree
std::vector<int> cutedges;
//use nbblocks to record the neighboring blocks of a partition 
std::vector< std::set<int> > nbblocks; 
std::vector<vid_t> newid;
vid_t startvid = 0;

//get the corresponding index of the label_id
int label_index(int label_id){
	
	std::map<int, int>::iterator mit = label2idx.find(label_id);	
	if(mit == label2idx.end()){
		std::cout<<"label id: "<<label_id<<" doesn't exists in map"<<std::endl;		
	}
	assert(mit != label2idx.end());
	return mit->second;	
}

void increaseMatrixCell(int i, int j){
	assert(i < Size && j < Size);
	lock.lock();
	adjmatrix[i][j]++;
	lock.unlock();
}

struct initAdjMatrix : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
 
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		if(vertex.num_edges() == 0) return;
		VertexDataType vdata = vertex.get_data();

		if (gcontext.iteration == 0) {
			//assert((vertex.get_data()).inbfs == true);
			//int vertex_id = (int)vertex.id();
			// use negative value to denote inactive state
			//vertex.set_data(-vertex_id);
			//VertexDataType vdata = vertex.get_data();
			//vdata.label = vertex_id;
			//vdata.inbfs = false;
			//vertex.set_data(vdata);
			int ith = label_index(vdata.label);
				
			for(int id = 0; id < vertex.num_outedges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.outedge(id);	
				EdgeDataType edata = e->get_data();
				int nblabel = edata.neighbor_label(vertex.id(), e->vertex_id());	
				if(vdata.label != nblabel){
					int jth = label_index(nblabel);
					increaseMatrixCell(ith, jth);	
				}
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

vid_t get_new_id(int label_id){
	int idx = label_index(label_id);							
	vid_t tmpid = 0;				
	lock.lock();
	tmpid = newid[idx]++;
	lock.unlock();
	return tmpid;	
}
//remap vertex id 
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
			vid_t new_id = get_new_id(vdata.label);
			//vdata.label = vertex_id;
			//vdata.inbfs = false;
			//vertex.set_data(vdata);
			for(int id = 0; id < vertex.num_edges(); id++)
			{
				graphchi_edge<EdgeDataType> * e = vertex.edge(id);	
				EdgeDataType edata = e->get_data();
				edata.my_label(vertex.id(), vertex.edge(id)->vertex_id()) = (int)new_id;
				//if(vdata.label >= 0)
				//edata.disable_propagate();
				e->set_data(edata);
				//std::cout<<"vid"<<vertex.edge(id)->vertexid
			}
			//assert(vertex.num_edges() ==4 );	
			lock.lock();
			fprintf(vfout, "%u\t%u\n", vertex.id(), new_id);	
			lock.unlock();
		}/*else{
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
			}
			
		}*/
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

//add a block to a partition, we need to insert the block's 
//unmerged neighbors to the corresponding set 
void add_neighbors(int partitionid, int blockid){
	assert(partitionid < (int)nbblocks.size());		
	if(blockid >= Size)
		std::cout<<"block id: "<<blockid<<"\t Size: "<<Size<<std::endl;
	assert(blockid < Size);
	for(int i=0; i<Size; i++){
		if(!merged[i] && adjmatrix[blockid][i] > 0){
			nbblocks[partitionid].insert(i);
		}else if(!merged[i] && adjmatrix[i][blockid] > 0){
			nbblocks[partitionid].insert(i);	
		}	
	}
}

float compute_block_score(int partitionid, int blockid){
	float score = (float)0;
	for(int i=0; i<(int)partitions[partitionid].size(); i++){
		int blkid = partitions[partitionid][i];		
		score += (adjmatrix[blkid][blockid]+adjmatrix[blockid][blkid])/(float)cutedges[blockid];	
	}						
	return score;
}
//compute the score between the partition and the block
int highest_block(int partitionid){
	std::set<int> nbsets = nbblocks[partitionid];			
	float high_score = (float)0;
	int bid = -1;//the corresponding block id of the block the highest score 
	std::set<int>::iterator sit;			
	for(sit = nbsets.begin(); sit != nbsets.end(); sit++){
		//int blockid = *sit;
		if(!merged[*sit]){
			float tmp_score = compute_block_score(partitionid, *sit); 	
			if(high_score < tmp_score){
				high_score = tmp_score;	
				bid = *sit;
			}
			//high_score = high_score < tmp_score ? tmp_score : high_score;
		}	
	}
	return bid;	
	//return high_score;
}

int largest_block(){
	int blksize = 0;	
	int blkid = -1;
	for(int i=0; i<(int)collect_result.size(); i++){
		if(!merged[i] && collect_result[i].count > blksize){
			blksize = collect_result[i].count;		
			//blkid = collect_result[i].label;
			blkid = i;
		}	
	}
	return blkid;
}
//choose the neighbor block with highest neighboring score
int best_nbblock(int partitionid){
	assert(partitionid < (int)nbblocks.size());	
	int blkid = highest_block(partitionid);	
	if(blkid != -1){
		return blkid;									
	}else{
		//this partition currently has no unmerged neighbors, choose the largest block from the rest
		blkid = largest_block();	 					
		assert(blkid != -1);
		return blkid;
	}	
	return blkid;
} 

int least_loaded_partition(){
	int blksize = totalvertices + 1;
	int blkid = 0;
	//choose the least loaded partition								
	for(int i=0; i<(int)partition_size.size(); i++){
		if(blksize > partition_size[i]){
			blksize = partition_size[i];
			blkid = i;
		}		
	}
	return blkid;	
}

void free_matrix(){
	if(adjmatrix == NULL)
		return;
	for(int i=0; i<Size; i++){
		free(adjmatrix[i]);
	}
	free(adjmatrix);
	adjmatrix = NULL;
}

//the vertex id range of each partition [start, end]
std::vector<std::pair<vid_t, vid_t> > interval;

int msbfsMain(int argc, const char ** argv) {
	//get the id range of the giant-SCC [start, end)
	std::vector<std::pair<vid_t, vid_t> > range = DAGmain(argc, argv);	
	bool left_empty = false;
	bool right_empty = false;
	assert((int)range.size() == 3);
	std::cout<<"DAG finished!============================"<<std::endl;
	startvid =	range[1].first; 
	vid_t endvid = range[1].second;
	//return 0;
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);
    
    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("giantSCC-msBFS");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 100000); // Number of iterations
    scheduler       	 = get_option_int("scheduler", false); // Whether to use selective scheduling
	max_block_size 		 = get_option_int("max_block", 100000);
	sample_rate			 = (double)get_option_float("sampling_rate", 0.001);
	int num_par		     = get_option_int("num_par", 5);
	max_iterations		 = get_option_int("max_iterations", 10);
	NumRoots			 = get_option_int("roots", 1000);	
	MAX_LEVEL		     = get_option_int("levels", 10);
	stop_ratio			 = (double)get_option_float("stopratio", 0.1);	
	//numPartitions	     = get_option_int("");
	
	//exclude left and right part of the Graph	
	//check left part is empty!
	if(range[0].first < range[0].second){
		num_par = num_par-1;
	}else{
		left_empty = true;	
	}
	
	//check if right part is empty
	if(range[2].first < range[2].second){
		num_par = num_par-1;
	}else{
		right_empty = true;		
	}
	assert(num_par > 0);
	max_iterations = MAX_LEVEL;
   	//single_source = get_option_int("root", 0); 
    /* Detect the number of shards or preprocess an input to create them */
	get_option_string("filetype", "auto");	
	std::string orig_filename = filename;
	filename += ".bigscc";
    int nshards          = convert_if_notexists<EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));
    
    /* Run */
    initBFS init_program;
    //bfs program;
	
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
	engine.set_save_edgesfiles_after_inmemmode(true);
	std::cout<<"---------------------------start init program----------------------------------"<<std::endl;
    engine.run(init_program, niters);
	/*
	graphchi_engine<VertexDataType, EdgeDataType> enginexx(filename, nshards, scheduler, m); 
	checkBFS check_program1;	
	enginexx.run(check_program1, niters);
	std::cout<<"check initBFS is finished!"<<std::endl;	
	*/
	int active_curr = 0;	
	int active_prev = active_vertices_count<VertexDataType>(filename);		

	scheduler = false;
    //graphchi_engine<VertexDataType, EdgeDataType> engine2(filename, nshards, scheduler, m); 
	msBFS msbfs_program;	
	int round = 0;
	int block_size_sum = 0;	
	int total_blocks = 0;
	double ratio = 0;
	int totalcount = 0;
	totalvertices = active_prev;//total number of non-isolated vertices
	int tmpcount = 0;	
	while(true){
		std::cout<<"------------------------------running msBFS round "<<round++<<"------------"<<std::endl;
		//roots.clear();
		roots = get_top_degree_vertices<VertexDataType>(filename, NumRoots); 	
		sizecount.clear();
		for(std::set<vid_t>::iterator it = roots.begin(); it != roots.end(); it++){
			sizecount.insert(std::make_pair(*it, 0));
		}
    	//engine2.run(msbfs_program, niters);
    	engine.run(msbfs_program, niters);
		active_curr = active_vertices_count<VertexDataType>(filename);		
		//ratio = (double)active_curr/active_prev;
		ratio = (double)active_curr/totalvertices;
		for(std::map<vid_t, int>::iterator it=sizecount.begin(); it != sizecount.end(); it++){
			if(tmpcount++ < 5)
				std::cout<<"vid="<<it->first<<"\t size="<<it->second<<std::endl;
			totalcount += it->second;
		}
		std::cout<<"active_prev: "<<active_prev<<"\t active_curr: "
			<<active_curr<<"\t totalcount: "<<totalcount<<"\t ratio: "<<ratio<<std::endl;
		active_prev = active_curr;
		if(ratio < stop_ratio){
			std::cout<<"============================stop ratio is reached, break!"<<std::endl;
			break;
		}
		//block_size_sum += sum_block();
		//total_blocks += sources.size();
		/*
		if(ratio > stop_ratio){
			std::cout<<"============================stop ratio is reached, break!"<<std::endl;
			break;	
		}			
		*/

		/*
		sample_rate *= multiply;	
		if(sample_rate > sample_max){
			break;
		}	
		*/
		//sources.clear();
	}
	
	scheduler = false;
	//graphchi_engine<VertexDataType, EdgeDataType> enginex(filename, nshards, scheduler, m); 
	//checkBFS check_program;	
	//enginex.run(check_program, niters);

	//graphchi_engine<VertexDataType, EdgeDataType> enginewcc(filename, nshards, scheduler, m); 
	//checkWCC check_wcc_program;	
	
	if(active_curr > 0){
		//graphchi_engine<VertexDataType, EdgeDataType> engine3(filename, nshards, scheduler, m); 
		initWCC initwcc_program;	
		engine.run(initwcc_program, niters);

		//enginewcc.run(check_wcc_program, niters);
		//std::cout<<"check wcc init finished"<<std::endl;	

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



	//std::vector<label_count> collect_result;
	//count_labels<VertexDataType, label_count>(filename, 40, collect_result); 
	msbfs_count_labels<VertexDataType, label_count>(filename, 40, collect_result); 
	std::cout<<"num blocks: "<<collect_result.size()<<"\t partition number: "<<num_par<<std::endl;
	assert((int)collect_result.size() > 0);
	merged.resize(collect_result.size(), false);
	Size = (int)collect_result.size();
	adjmatrix = (int **) malloc(sizeof(int*)*Size);
	for(int i=0; i<(int)collect_result.size(); i++){
		adjmatrix[i] = (int*)malloc(sizeof(int)*Size);	
		memset(adjmatrix[i], 0, sizeof(int)*Size);
	}
	
	for(int i=0; i<(int)collect_result.size(); i++){
		label2idx.insert(std::pair<int, int>(collect_result[i].label, i));	
	}	
	std::cout<<"------------------start construct adjacent Matrix--------------"<<std::endl;
	//graphchi_engine<VertexDataType, EdgeDataType> engine5(filename, nshards, scheduler, m); 
	engine.set_save_edgesfiles_after_inmemmode(false);
	initAdjMatrix matrix_program;	
	engine.run(matrix_program, niters);
	
	std::cout<<"-------------------start assigning blocks to partitions--------------------"<<std::endl;
	partitions.resize(num_par);	
	partition_size.resize(num_par, 0);
	nbblocks.resize(num_par);
	//std::vector< std::set<int> > nbblocks; 
	//number of cut-edges of each bfs tree
	cutedges.resize(Size, 0);
	//initialize cutedges array
	for(int i=0; i<Size; i++){
		assert(adjmatrix[i][i] == 0);
		//sum i-th row and i-th column
		for(int j=0; j<Size; j++){
			cutedges[i] += adjmatrix[i][j];	
			cutedges[i] += adjmatrix[j][i];
		}	
	}
	//number of unmerged blocks
	int active_blocks = (int)collect_result.size();
	for(int i=0; i<num_par; i++){
		//choose the top num_par blocks as the initial vertex
		//label_count lcount = collect_result[i];

		partitions[i].push_back(i);		
		//partitions[i].push_back(collect_result[i].label);		
		partition_size[i] = collect_result[i].count;	
		//merged[label_index(collect_result[i].label)] = true;	
		merged[i] = true;	
		add_neighbors(i, i);
	}

	active_blocks = active_blocks - num_par;				

	while(active_blocks-- > 0){
		//choose the least loaded partition to merge blocks	
		int partid = least_loaded_partition();		
		//choose the most connected neighbor block to merge 
		int blkid = best_nbblock(partid);			
		//partitions[partid].push_back(collect_result[blkid].label);
		partitions[partid].push_back(blkid);
		partition_size[partid] += collect_result[blkid].count;
		merged[blkid] = true;
		add_neighbors(partid, blkid);
	}
					
	//FILE* finterval = fopen((orig_filename+".dag.interval").c_str(), "w+");
	//interval.resize(num_par+2);
	//vertex id range of the left part of the DAG
	if(!left_empty)	
		interval.push_back(std::pair<vid_t, vid_t>(range[0].first, range[0].second-1));	

	//calculate the vertex id range of each block			
	newid.resize(collect_result.size());	
	//vid_t startid = ;				
	for(int i=0; i<(int)partitions.size(); i++){
		int start_tmp = startvid;
		for(int j=0; j<(int)partitions[i].size(); j++){
			int blkid = partitions[i][j];
			newid[blkid] = startvid;	
			startvid += collect_result[blkid].count; 	
		}
		//the vertex id range of each partition in  the middle part(the giant-SCC) 
		interval.push_back(std::pair<vid_t, vid_t>(start_tmp, startvid-1));			
	}
	assert(endvid == startvid);	
	//vertex id range of the right part of DAG
	if(!right_empty)
		interval.push_back(std::pair<vid_t, vid_t>(range[2].first, range[2].second-1));

	std::cout<<"---------------ReMap is started-------------------"<<std::endl;
	//open vertex map file for appending	
	vfout = fopen((orig_filename+".vmap").c_str(), "a");
	assert(vfout != NULL);
	//efout = fopen((filename+".be").c_str(), "w+");
	//assert(vfout != NULL && efout != NULL);
	fprintf(vfout, "# old_vid\tnew_vid\n");
	//fprintf(efout, "# new_src\tnew_dst\n");

	//graphchi_engine<Vertexinfo, EdgeDataType> engine(filename, nshards, scheduler, m);	
	ReMap remap;
	//engine.set_save_edgesfiles_after_inmemmode(true);
	engine.run(remap, 1);
	
	fclose(vfout);
	free_matrix();
	//fclose(efout);
	FILE* finterval = fopen((orig_filename+".dag.interval").c_str(), "w+");
	assert(finterval != NULL);	
	for(int i=0; i<(int)interval.size(); i++){
		fprintf(finterval, "%u\t%u\n", interval[i].first, interval[i].second);		
	}
	fclose(finterval);
	return 0;
	/*
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
	
	//std::cout<<"assigning blocks to partitions finished!"<<std::endl;
	//calculate the prefix sum of all partitions and blocks		
	int prefix_tmp = 0;
	for(int i=0; i<(int)partitions.size(); i++){
		int tmp = partitions[i];	
		partitions[i] = prefix_tmp;
		prefix_tmp += tmp;
	}	
	//std::cout<<"assigning blocks to partitions finished!1"<<std::endl;
	//start_vid store the starting vertex id for each block;
	for(int i=0; i<(int)collect_result.size(); i++){
		start_vid[i] = partitions[block_to_par[i]];	
		partitions[block_to_par[i]] += collect_result[i].count;	
	}
	//std::cout<<"assigning blocks to partitions finished!2"<<std::endl;
	vfout = fopen((filename+".bv").c_str(), "w+");
	efout = fopen((filename+".be").c_str(), "w+");
	assert(vfout != NULL && efout != NULL);
	fprintf(vfout, "# new_vid  old_vid\n");
	fprintf(efout, "# new_src  new_dst\n");
	fflush(vfout);
	fflush(efout);
	std::cout<<"ReMap is started"<<std::endl;
	graphchi_engine<Vertexinfo, EdgeDataType> engine5(filename, nshards, scheduler, m);	
	ReMap remap;
	engine5.run(remap, 2);
	fclose(vfout);
	fclose(efout);
	//std::cout<<"block number is "<<block_num<<std::endl;
	std::cout<<"total blocks: "<<total_blocks<<"\t visited vertices "<<block_size_sum<<"\t active: "
			<<active_curr<<"\t total vertices: "<< get_num_vertices(filename)<<std::endl;
	std::cout<<"sample_rate "<<sample_rate<<"\t sample_max "<<sample_max
			<<"\t round: "<<round<<"\t a_ratio: "<<ratio<<"\tstop_ratio: "<<stop_ratio<<std::endl;
	*/

    /* Report execution metrics */
    metrics_report(m);
    return 0;
}
