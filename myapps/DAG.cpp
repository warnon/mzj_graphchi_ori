

#define SUPPORT_DELETIONS 1

#include <string>
#include <ostream>

#include "graphchi_basic_includes.hpp"
//#include "util/labelanalysis.hpp"
#include "util/maxdeg.cpp"

using namespace graphchi;


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

//int super_step = 0;


// Id for the output stream for contracted graph
//int CONTRACTED_GRAPH_OUTPUT;
/*
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
*/
struct SCCinfo {
    vid_t color;
    bool confirmed;
	bool reconfirmed;
    SCCinfo() : color(0), confirmed(false), reconfirmed(false) {}
    SCCinfo(vid_t color) : color(color), confirmed(false), reconfirmed(false) {}
    SCCinfo(vid_t color, bool confirmed) : color(color), confirmed(confirmed), reconfirmed(false) {}
    
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

typedef SCCinfo VDataType;
typedef bidirectional_label EDataType;


bool first_iteration = true;
bool remainingvertices = true;
//bool scheduler = false;
vid_t root = 0;
bool isconverged = false;
unsigned long product = 0;

/* use Fw-Bw SCC finding algorithm to extract the SCC that contains the root vertex
 * created by mzj 2016/3/13
 */

/**
 * FORWARD-PHASE
 */
struct SCCForward : public GraphChiProgram<VDataType, EDataType> {
    
    
	//mutex lock;
    /**
     *  Vertex update function.
     */
	void update(graphchi_vertex<VDataType, EDataType> &vertex, graphchi_context &gcontext) {

		VDataType vertexdata; //= vertex.get_data();
		//bool propagate = false;
		if (gcontext.iteration == 0) {

			vertex.set_data(SCCinfo(vertex.id()));
			vertexdata = vertex.get_data();
			vertexdata.color = vertex.id();


			if(vertex.id() == root){
				product = (unsigned long)vertex.num_inedges() * vertex.num_outedges(); 

				vertexdata.confirmed = true;
				vertex.set_data(vertexdata);
				for(int i=0; i<vertex.num_outedges(); i++){
					bidirectional_label edgedata = vertex.outedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid) = vertex.id();
					//if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
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
				isconverged = false;
				vertexdata.confirmed = true;
				vertexdata.color = min_color;
				for(int i=0; i<vertex.num_outedges(); i++){
					bidirectional_label edgedata = vertex.outedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid) = min_color;
					//if(scheduler) gcontext.scheduler->add_task(vertex.outedge(i)->vertexid);	
					vertex.outedge(i)->set_data(edgedata);
				}	
				vertex.set_data(vertexdata);
			}
		}
	}
    
	void before_iteration(int iteration, graphchi_context &gcontext) {
        //first_iteration = false;
		isconverged = iteration > 0;
    }

    
    void after_iteration(int iteration, graphchi_context &gcontext) {
        //first_iteration = false;
		if(isconverged){
			logstream(LOG_INFO)<<"scc_forward has finished!"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}
    }
};

/**
 * BACKWARD phase
 */
struct SCCBackward : public GraphChiProgram<VDataType, EDataType> {
    
    /**
     *  Vertex update function.
     */
    void update(graphchi_vertex<VDataType, EDataType> &vertex, graphchi_context &gcontext) {
              
		VDataType vertexdata; //= vertex.get_data();
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
				vertexdata.reconfirmed = true; 
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = vertex.id();
					vertex.inedge(i)->set_data(edgedata);
					//if(scheduler) gcontext.scheduler->add_task(vertex.inedge(i)->vertexid);	
					vertex.inedge(i)->set_data(edgedata);
				}	
				vertex.set_data(vertexdata);	
			}else{
				vertexdata.reconfirmed = false;
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
				isconverged = false;
				//vertexdata.confirmed = true;
				vertexdata.reconfirmed = true;
				vertexdata.color = min_color;
				for(int i=0; i<vertex.num_inedges(); i++){
					bidirectional_label edgedata = vertex.inedge(i)->get_data();
					edgedata.my_label(vertex.id(), vertex.inedge(i)->vertexid) = min_color;
					//if(scheduler) gcontext.scheduler->add_task(vertex.inedge(i)->vertexid);	
					vertex.inedge(i)->set_data(edgedata);
				}	
				vertex.set_data(vertexdata);
			}
		}
    }
    
    void before_iteration(int iteration, graphchi_context &gcontext) {
		isconverged = iteration > 0;
	}
    
    void after_iteration(int iteration, graphchi_context &gcontext) {
		if(isconverged){
			logstream(LOG_INFO)<<"scc_backward has finished!"<<std::endl;
			gcontext.set_last_iteration(iteration);
		}
	}
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) { }
    
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &gcontext) {}
};

//graphchi_engine<VDataType, EDataType> * gengine = NULL;
FILE* fpout = NULL;
//FILE* fpout1 = NULL;
/* Simple contraction step that just outputs the non-deleted edges. Would be better
 done automatically, but the dynamic engine is a bit flaky. */
struct ContractionStep : public GraphChiProgram<VDataType, EDataType> {

	mutex lock;

	void update(graphchi_vertex<VDataType, EDataType> &vertex, graphchi_context &gcontext) {
		if(gcontext.iteration == 0){
			VDataType vertexdata = vertex.get_data();
			if(!vertexdata.confirmed || !vertexdata.reconfirmed)
				return ;	
			//assert((unsigned long)(vertex.num_inedges() * vertex.num_outedges()) <= product);

			for(int i=0; i<vertex.num_outedges(); i++){
				bidirectional_label edgedata = vertex.outedge(i)->get_data();
				if(edgedata.is_equal()){		
					if(root == edgedata.my_label(vertex.id(), vertex.outedge(i)->vertexid)){
						//std::cout<<"root is "<<root<<std::endl;
						lock.lock();
						fprintf(fpout, "%u\t%u\n", vertex.id(), vertex.outedge(i)->vertexid);
						lock.unlock();
						//continue;
					}
				}
				/*
				   lock.lock();
				   fprintf(fpout1, "%u\t%u\n", vertex.id(), vertex.outedge(i)->vertexid);
				   lock.unlock();
				   */
			}
		}
	}

	void before_iteration(int iteration, graphchi_context &gcontext) {
		//isconverged = iteration > 0;
		assert(fpout != NULL);
		//assert(fpout1 != NULL);
		fflush(fpout);
		//fflush(fpout1);
	}
	 void after_iteration(int iteration, graphchi_context &gcontext) {
	//	if(isconverged){
	//		logstream(LOG_INFO)<<"scc_backward has finished!"<<std::endl;
			fflush(fpout);
			//fflush(fpout1);
			gcontext.set_last_iteration(iteration);
	//	}
	}
};

vid_t left = 0;
vid_t middle = 0; 
vid_t right = 0;
vid_t left_bound = 0;
vid_t middle_bound = 0;
vid_t right_bound = 0;
struct CountingStep : public GraphChiProgram<VDataType, EDataType> {

	mutex lock;

	void update(graphchi_vertex<VDataType, EDataType> &vertex, graphchi_context &gcontext) {
		
	//	assert(vertex.num_inedges() * vertex.num_outedges() <= product);
		if(gcontext.iteration == 0){
			if(vertex.num_edges() == 0)
				return;
			VDataType vertexdata = vertex.get_data();
			if(!vertexdata.confirmed){
				lock.lock();
				left++;
				lock.unlock();
				return;
			}

			if(vertexdata.confirmed && vertexdata.reconfirmed){
				lock.lock();
				middle++;	
				lock.unlock();
			}else{
				lock.lock();
				right++;
				lock.unlock();
			}
		}	
	}
	void before_iteration(int iteration, graphchi_context &gcontext) {
		//isconverged = iteration > 0;
	}
	 void after_iteration(int iteration, graphchi_context &gcontext) {
			gcontext.set_last_iteration(iteration);
	}
};

FILE* vmap = NULL;
mutex lk;
vid_t getNewIdLeft(){
	vid_t tmp = 0;
	lk.lock();
	tmp = left++;	
	lk.unlock();		
	assert(tmp < left_bound);	
	return tmp;
} 

vid_t getNewIdRight(){
	vid_t tmp = 0;
	lk.lock();
	tmp = right++;	
	lk.unlock();
	assert(tmp < right_bound);
	return tmp;
}

struct RemapStep : public GraphChiProgram<VDataType, EDataType> {

	mutex lock;

	void update(graphchi_vertex<VDataType, EDataType> &vertex, graphchi_context &gcontext) {
		if(vertex.num_edges() == 0)
			return ;
		VDataType vertexdata = vertex.get_data();
		if(vertexdata.confirmed && vertexdata.reconfirmed)
			return ;	
		//assert(vertex.num_inedges() * vertex.num_outedges() <= product);
		if (gcontext.iteration == 0){
			if(vertexdata.confirmed){
				vertexdata.color =	getNewIdRight(); 
			}else{
				vertexdata.color = getNewIdLeft();
			}	
			vertex.set_data(vertexdata);
		}else if(gcontext.iteration == 1){ 
			lock.lock();
			fprintf(vmap, "%u\t%u\n", vertex.id(), vertexdata.color);
			lock.unlock();
		}
	}

	void before_iteration(int iteration, graphchi_context &gcontext) {
		//isconverged = iteration > 0;
		/*
		assert(vmap != NULL);
		fflush(vmap);
		*/
	}

	void after_iteration(int iteration, graphchi_context &gcontext) {
		if(iteration == 1){
			//fflush(vmap);
			gcontext.set_last_iteration(iteration);
		}
	}
};

std::vector<std::pair<vid_t, vid_t> > DAGmain(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line
     arguments and the configuration file. */
    graphchi_init(argc, argv);
    global_logger().set_log_level(LOG_DEBUG);
    
    /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
    metrics m("DAG-msBFS");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    //bool scheduler       = true;
    
    /* Detect the number of shards or preprocess an input to create them */
    
    int nshards = find_shards<EDataType>(filename);
	/*
    if (nshards > 0) {
        delete_shards<EDataType>(filename, nshards);
    }
   	*/ 
    nshards          = convert_if_notexists<EDataType>(filename,
                                                          get_option_string("nshards", "auto"));
    
    
	root = (vid_t)get_option_int("root", -1); 
	bool scheduler = get_option_int("scheduler", false); 
	int niters = get_option_int("niters", 1000);
    /* Run */
    fpout = fopen((filename+".bigscc").c_str(), "w+"); 
    vmap = fopen((filename+".vmap").c_str(), "w+"); 

	assert(fpout != NULL);
	assert(vmap != NULL);

	fprintf(vmap, "#old_vid\tnew_vid\n");
	fprintf(fpout, "#src_vid\tdst_vid\n");
	//use max outdegree*indegree as the pivot to find the largest SCC
	if(root == (vid_t)-1) root = GetMaxDegreeVertex(filename);	
	graphchi_engine<VDataType, EDataType> engine(filename, nshards, scheduler, m);
	engine.set_save_edgesfiles_after_inmemmode(true);
	//graphchi_engine<VDataType, EDataType> engine1(filename, nshards, scheduler, m);
	//graphchi_engine<VDataType, EDataType> engine2(filename, nshards, scheduler, m);
	SCCForward forwardscc;	
	engine.run(forwardscc, niters);	
	SCCBackward backward;
	engine.run(backward, niters);
	CountingStep count;
	engine.run(count, niters);

	std::cout<<"root vertex is: "<<root<<std::endl;
	std::cout<<"left: "<<left<<"\tmiddle: "<<middle<<"\tright: "<<right<<std::endl;
	int sum = left + middle + right;
	assert(sum != 0);
	std::cout<<"ratio left: "<<left/(float)sum<<"\tmiddel: "<<middle/(float)sum<<"\tright: "<<right/(float)sum<<std::endl;

	std::cout<<"---------------------running contraction step-------------------"<<std::endl;
	engine.set_save_edgesfiles_after_inmemmode(false);
	ContractionStep cstep;
	engine.run(cstep, niters);
    //analyze_labels<VDataType>(filename);
    
    //delete_shards<EDataType>(filename, nshards);
   	fclose(fpout); 
	//remap the left and right part of the graph

	std::vector<std::pair<vid_t, vid_t> > range;

	left_bound = left;	
	middle_bound = left+middle;
	right_bound = left + middle + right; 

	right = left+middle;
	middle = left;
	left = 0;

	range.push_back(std::pair<vid_t, vid_t>(left, left_bound));	
	range.push_back(std::pair<vid_t, vid_t>(middle, middle_bound));	
	range.push_back(std::pair<vid_t, vid_t>(right, right_bound));	

	std::cout<<"---------------------running remap step-------------------"<<std::endl;
	//engine.set_save_edgesfiles_after_inmemmode(true);
	RemapStep remap;	
	engine.run(remap, niters);
	/*	
	std::cout<<"---------------------running contraction step-------------------"<<std::endl;
	ContractionStep cstep;
	engine.run(cstep, niters);
    //analyze_labels<VDataType>(filename);
    
    //delete_shards<EDataType>(filename, nshards);
   	fclose(fpout); 
	*/
   	fclose(vmap); 
    
    /* Report execution metrics */
    //metrics_report(m);
	
	std::cout<<"DAGmain finished==================="<<std::endl;
    return range;
}
