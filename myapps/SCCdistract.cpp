
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
 * Strongly Connected Components. Based on technical report (2012):
 @article{salihoglucomputing,
 title={Computing Strongly Connected Components in Pregel-like Systems},
 author={Salihoglu, Semih and Widom, Jennifer},
 publisher={Stanford InfoLab}
 }
 */


//#define SUPPORT_DELETIONS 1

#include <string>
#include <ostream>

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"
#include "util/maxdeg.cpp"

using namespace graphchi;


/**
 * Unlike in weakly connected components, we need
 * to ensure that neighbors do not overwrite each
 * others values. This is achieved by keeping two values
 * in an edge. In this struct, smaller_one is the id of the
 * vertex that has smaller id, and larger_one the others.
 * This complexity is due to us ignoring the direction of an edge.
 */
struct bidirectional_label {
    vid_t smaller_one;
    vid_t larger_one;
    
   	bool is_equal(){
		 assert(larger_one != 0xffffffffu);
        assert(smaller_one != 0xffffffffu);
		return larger_one == smaller_one;
	} 
    vid_t & neighbor_label(vid_t myid, vid_t nbid) {
        assert(larger_one != 0xffffffffu);
        assert(smaller_one != 0xffffffffu);
        
        if (myid < nbid) {
            return larger_one;
        } else {
            return smaller_one;
        }
    }
    
    vid_t & my_label(vid_t myid, vid_t nbid) {
        assert(larger_one != 0xffffffffu);
        assert(smaller_one != 0xffffffffu);
        
        if (myid < nbid) {
            return smaller_one;
        } else {
            return larger_one;
        }
    }
    
    // Annoying hack
    bool deleted() {
        return smaller_one == 0xffffffffu;
    }
};

int super_step = 0;


// Id for the output stream for contracted graph
int CONTRACTED_GRAPH_OUTPUT;

struct SCCinfo {
    vid_t color;
    bool confirmed;
    SCCinfo() : color(0), confirmed(false) {}
    SCCinfo(vid_t color) : color(color), confirmed(false) {}
    SCCinfo(vid_t color, bool confirmed) : color(color), confirmed(confirmed) {}
    
    friend std::ostream& operator<< (std::ostream &out, SCCinfo &scc) {
        out << scc.color;
        return out;
    }
    
};

/* Overloaded operators to help with labelanalysis.hpp */

bool operator<(const SCCinfo &a, const SCCinfo &b);
bool operator<(const SCCinfo &a, const SCCinfo &b) {
    return a.color < b.color;
}

bool operator==(const SCCinfo &a, const SCCinfo &b);
bool operator==(const SCCinfo &a, const SCCinfo &b) {
    return a.color == b.color;
}
bool operator!=(const SCCinfo &a, const SCCinfo &b);
bool operator!=(const SCCinfo &a, const SCCinfo &b) {
    return a.color != b.color;
}

typedef SCCinfo VertexDataType;
typedef bidirectional_label EdgeDataType;

static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(bidirectional_label val);
static inline bool VARIABLE_IS_NOT_USED is_deleted_edge_value(bidirectional_label val) {
    return 0xffffffffu == val.smaller_one;
}

static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<bidirectional_label> * e);
static void VARIABLE_IS_NOT_USED remove_edgev(graphchi_edge<bidirectional_label> * e) {
    bidirectional_label deletedlabel;
    deletedlabel.smaller_one = 0xffffffffu;
    deletedlabel.larger_one = 0xffffffffu;
    e->set_data(deletedlabel);
}

bool first_iteration = true;
bool remainingvertices = true;
bool scheduler = false;
vid_t root = 0;
bool converged = false;
mutex lock;
unsigned long product = 0;

/* use Fw-Bw SCC finding algorithm to extract the SCC that contains the root vertex
 * created by mzj 2016/3/13
 */

/**
 * FORWARD-PHASE
 */
struct SCCForward : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    
    /**
     *  Vertex update function.
     */
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		VertexDataType vertexdata; //= vertex.get_data();
		bool propagate = false;
		if (gcontext.iteration == 0) {

			vertex.set_data(SCCinfo(vertex.id()));
			vertexdata = vertex.get_data();
			vertexdata.color = vertex.id();


			if(vertex.id() == root){
				product = vertex.num_inedges() * vertex.num_outedges(); 

				vertexdata.confirmed = true;
				vertex.set_data(vertexdata);
				for(int i=0; i<vertex.num_outedges(); i++){
					bidirectional_label edgedata = vertex.outedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid) = vertex.id();
					if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
					vertex.outedge(i)->set_data(edgedata);
				}	
			}else{
				vertexdata.confirmed = false;
				vertex.set_data(vertexdata);
				for(int i=0; i<vertex.num_outedges(); i++){
					bidirectional_label edgedata = vertex.outedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid) = vertex.id();
					//if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
					vertex.outedge(i)->set_data(edgedata);
				}	
				// initialize labels on in and out edges
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = vertex.id();
					//if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
					vertex.inedge(i)->set_data(edgedata);
				}	
			}
		} else {
			if(true == vertexdata.confirmed)
				return ;
			vertexdata = vertex.get_data();
			vid_t min_color = vertexdata.color;
			for(int i=0; i<vertex.num_inedges(); i++){
				//min_color = std::min(min_color, vertexdata.inedge(i)->get_data().neighbor_label(vertex.id(), vertex.inedge(i)->vertexid));		
				if(root == (vertex.inedge(i)->get_data()).neighbor_label(vertex.id(), vertex.inedge(i)->vertexid)){
						min_color = root;
						break;
				}
			}
			if(min_color != vertexdata.color){
				converged = false;
				vertexdata.confirmed = true;
				vertexdata.color = min_color;
				for(int i=0; i<vertex.num_outedges(); i++){
					bidirectional_label edgedata = vertex.outedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid) = min_color;
					if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
					vertex.outedge(i)->set_data(edgedata);
				}	
				vertex.set_data(vertexdata);
			}
		}
	}
    
	void before_iteration(int iteration, graphchi_context &gcontext) {
        //first_iteration = false;
		converged = iteration > 0;
    }

    
    void after_iteration(int iteration, graphchi_context &gcontext) {
        //first_iteration = false;
		if(converged){
			logstream(LOG_INFO)<<"scc_forward has finished!"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}
    }
};

/**
 * BACKWARD phase
 */
struct SCCBackward : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
              
		VertexDataType vertexdata; //= vertex.get_data();
		//bool propagate = false;
		if (gcontext.iteration == 0) {
			//vertex.set_data(SCCinfo(vertex.id()));
			vertexdata = vertex.get_data();
			/* vertices that is not visited in Fw phase is not in the giant SCC!
			 * minor improve by mzj 2016/3/13
			 */
			if(!vertexdata.confirmed)
				return;
			//assert(vertexdata.color == root);
			if(vertex.id() == root){
				//vertexdata.confirmed = true;
				vertexdata.color = vertex.id(); 
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = vertex.id();
					vertex.inedge(i)->set_data(edgedata);
					if(scheduler) gcontext.scheduler->add_task(vertex.inedge(i)->vertexid);	
					vertex.inedge(i)->set_data(edgedata);
				}	
				vertex.set_data(vertexdata);	
			}else{
				//vertexdata.confirmed = false;
				vertexdata.color = vertex.id(); 
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = vertex.id();
					vertex.inedge(i)->set_data(edgedata);
					//if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
				}	
				vertex.set_data(vertexdata);	
			}
			//vertex.set_data(vertexdata);	
		} else {
			vertexdata = vertex.get_data();
			if(!vertexdata.confirmed)
				return ;
			vid_t min_color = vertexdata.color;
			for(int i=0; i<vertex.num_outedges(); i++){
				//min_color = std::min(min_color, vertexdata.inedge(i)->get_data().neighbor_label(vertex.id(), vertex.inedge(i)->vertexid));		
				if(root == (vertex.outedge(i)->get_data()).neighbor_label(vertex.id(), vertex.outedge(i)->vertexid)){
						min_color = root;
						break;
				}
			}
			if(min_color != vertexdata.color){
				converged = false;
				//vertexdata.confirmed = true;
				vertexdata.color = min_color;
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = min_color;
					if(scheduler) gcontext.scheduler->add_task(vertex.inedge(i)->vertexid);	
					vertex.inedge(i)->set_data(edgedata);
				}	
				vertex.set_data(vertexdata);
			}
			/*
			for(int i=0; i<vertex.num_inedges(); i++){
				bidirectional_label edgedata = vertex.inedge(i)->get_data();
				edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = root;
				vertex.inedge(i)->set_data(edgedata);	
			}
			*/
			/*
			if(false == vertexdata.confirmed){
				converged = false;
				vertexdata.confirmed = true;
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = root;
					if(scheduler) gcontext.scheduler->add_task(vertex.inedge(i)->vertexid);	
					vertex.inedge(i)->set_data(edgedata);	
				}	
			}
			vertex.set_data(vertexdata);	
			*/
		}
    }
    
    void before_iteration(int iteration, graphchi_context &gcontext) {
		converged = iteration > 0;
	}
    
    void after_iteration(int iteration, graphchi_context &gcontext) {
		if(converged){
			logstream(LOG_INFO)<<"scc_backward has finished!"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}
	}
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) { }
    
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {}
};

graphchi_engine<VertexDataType, EdgeDataType> * gengine = NULL;
FILE* fpout = NULL;
FILE* fpout1 = NULL;
/* Simple contraction step that just outputs the non-deleted edges. Would be better
 done automatically, but the dynamic engine is a bit flaky. */
struct ContractionStep : public GraphChiProgram<VertexDataType, EdgeDataType> {

	void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
		
		assert(vertex.num_inedges() * vertex.num_outedges() <= product);

		for(int i=0; i<vertex.num_outedges(); i++){
			bidirectional_label edgedata = vertex.outedge(i)->get_data();
			if(edgedata.is_equal()){		
				if(root == edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid)){
					lock.lock();
						fprintf(fpout, "%u\t%u\n", vertex.id(), vertex.outedge(i)->vertexid);
					lock.unlock();
					continue;
				}
			}
			lock.lock();
			fprintf(fpout1, "%u\t%u\n", vertex.id(), vertex.outedge(i)->vertexid);
			lock.unlock();
		}
	}
	void before_iteration(int iteration, graphchi_context &gcontext) {
		//converged = iteration > 0;
		assert(fpout != NULL);
		assert(fpout1 != NULL);
		fflush(fpout);
		fflush(fpout1);
	}
	 void after_iteration(int iteration, graphchi_context &gcontext) {
	//	if(converged){
	//		logstream(LOG_INFO)<<"scc_backward has finished!"<<std::endl;
			fflush(fpout);
			fflush(fpout1);
			gcontext.set_last_iteration(iteration);
	//	}
	}
};

int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line
     arguments and the configuration file. */
    graphchi_init(argc, argv);
    global_logger().set_log_level(LOG_DEBUG);
    
    /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
    metrics m("giant-SCC-distract");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    //bool scheduler       = true;
    
    /* Detect the number of shards or preprocess an input to create them */
    
    int nshards = find_shards<EdgeDataType>(filename);
	/*
    if (nshards > 0) {
        delete_shards<EdgeDataType>(filename, nshards);
    }
   	*/ 
    nshards          = convert_if_notexists<EdgeDataType>(filename,
                                                          get_option_string("nshards", "auto"));
    
    
	root = get_option_int("root", -1); 
	scheduler = get_option_int("scheduler", false); 
	int niters = get_option_int("niters", 1000);
    /* Run */
    fpout = fopen((filename+".bigscc").c_str(), "w+"); 
    fpout1 = fopen((filename+".smallscc").c_str(), "w+"); 
	assert(fpout != NULL);
	assert(fpout1 != NULL);
	//use max outdegree*indegree as the pivot to find the largest SCC
	if(root == -1) root = GetMaxDegreeVertex(filename);	
	graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m);
	//graphchi_engine<VertexDataType, EdgeDataType> engine1(filename, nshards, scheduler, m);
	//graphchi_engine<VertexDataType, EdgeDataType> engine2(filename, nshards, scheduler, m);
	//engine.set_save_edgesfiles_after_inmemmode(true);
	SCCForward forwardscc;	
	engine.run(forwardscc, niters);	
	SCCBackward backward;
	engine.run(backward, niters);
	ContractionStep cstep;
	engine.run(cstep, niters);

    analyze_labels<VertexDataType>(filename);
    
    //delete_shards<EdgeDataType>(filename, nshards);
   	fclose(fpout); 
   	fclose(fpout1); 
    
    /* Report execution metrics */
    metrics_report(m);
	std::cout<<"root vertex is: "<<root<<std::endl;
    return 0;
}
