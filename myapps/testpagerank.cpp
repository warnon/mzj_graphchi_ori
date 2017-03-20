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

#define GRAPHCHI_DISABLE_COMPRESSION


#include "graphchi_basic_includes.hpp"
#include "util/toplist.hpp"

using namespace graphchi;
 
#define THRESHOLD 1e-2    
#define RANDOMRESETPROB 0.15


typedef float VertexDataType;
typedef float EdgeDataType;

float epsilon = 0.001;
//int* array = NULL;
bool scheduler;
size_t num_vertices = 0;
std::vector<float> pr_value1;
std::vector<float> pr_value8;
//bool num_tasks_print = false;
struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
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
		//std::cout<<"delta sum: "<<ginfo.get_delta()<<std::endl;
		if(converged){
			logstream(LOG_INFO)<<"converged!"<<std::endl;
			ginfo.set_last_iteration(iteration);
		}
    }
    
    /**
      * Called before an execution interval is started. Not implemented.
      */
	void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
		interval_converged = true;
		/*
		if(num_tasks_print){
			vid_t sum = ginfo.scheduler->total_tasks(window_st, window_en);
			logstream(LOG_INFO)<<"num of vertices scheduled="<<sum<<"/"<<window_en-window_st+1<<std::endl;			
		}
		*/
	}
    
  	bool repeat_updates(graphchi_context &gcontext){
		/*
		if(gcontext.iteration == 0)
			return false;
		else{
			logstream(LOG_INFO)<<"repeated_updates isconverged="<<interval_converged<<"\t iteration="<<gcontext.iteration<<std::endl;
			return !interval_converged;
		}
		*/
		return false;
	}  
    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {

		//array[v.id()]++;		
        if (ginfo.iteration == 0) {
            /* On first iteration, initialize vertex and out-edges. 
               The initialization is important,
               because on every run, GraphChi will modify the data in the edges on disk. 
             */
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                edge->set_data(1.0 / v.num_outedges());
            }
            v.set_data(RANDOMRESETPROB); 
			if(scheduler){
				ginfo.scheduler->add_task(v.id());	
			}
        } else {
			float old_value = v.get_data();
            /* Compute the sum of neighbors' weighted pageranks by
               reading from the in-edges. */
        	float sum=0;
            for(int i=0; i < v.num_inedges(); i++) {
                //float val = v.inedge(i)->get_data();
                sum  += v.inedge(i)->get_data();
            }
            
            /* Compute my pagerank */
            float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
			/*
			if(v.id() == 0){
				//std::cout<<"id="<<v.id()<<" newRank="<<pagerank<<" oldRank="<<old_value<<std::endl;		
				pr_value.push_back(old_value);
			}*/
           	/*	
			//check for fluctuate
			if(ginfo.iteration > 1 && pagerank > old_value){
				std::cout<<"id="<<v.id()<<" pagerank="<<pagerank<<" old_rank="<<old_value<<" indeg="<<v.num_inedges()<<" outdeg="<<v.num_outedges()<<std::endl;
				assert(false);
			}	 
			*/
            /* Write my pagerank divided by the number of out-edges to
               each of my out-edges. */
            if (v.num_outedges() > 0) {
                float pagerankcont = pagerank / v.num_outedges();
                for(int i=0; i < v.num_outedges(); i++) {
                    graphchi_edge<float> * edge = v.outedge(i);
                    edge->set_data(pagerankcont);
					if(scheduler){
						vid_t eid = edge->vertexid;
						if(eid <= ginfo.interval_en)
							ginfo.scheduler->add_task(eid, false);
						else
							ginfo.scheduler->add_task(eid, true);
					}
                }
            }
                
            /* Keep track of the progression of the computation.
               GraphChi engine writes a file filename.deltalog. */
            ginfo.log_change(std::abs(pagerank - v.get_data()));
            
            /* Set my new pagerank as the vertex value */
            v.set_data(pagerank); 
			if(std::abs(pagerank-old_value) > epsilon){
				if(converged)
					converged = false;
				if(interval_converged)
					interval_converged = false;
			}
        }
		
		if(v.id() == 1 || v.id() == 8){
			if(v.id() == 1){
				pr_value1.push_back(v.get_data());	
			}else{
				pr_value8.push_back(v.get_data());
			}
			
		}
    }
    
};

/**
  * Faster version of pagerank which holds vertices in memory. Used only if the number
  * of vertices is small enough.
  */
struct PagerankProgramInmem : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    std::vector<EdgeDataType> pr;
    PagerankProgramInmem(int nvertices) :   pr(nvertices, RANDOMRESETPROB) {}
    
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
        if (ginfo.iteration > 0) {
            float sum=0;
            for(int i=0; i < v.num_inedges(); i++) {
              sum += pr[v.inedge(i)->vertexid];
            }
            if (v.outc > 0) {
                pr[v.id()] = (RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum) / v.outc;
            } else {
                pr[v.id()] = (RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum);
            }
        } else if (ginfo.iteration == 0) {
            if (v.outc > 0) pr[v.id()] = 1.0f / v.outc;
        }
        if (ginfo.iteration == ginfo.num_iterations - 1) {
            /* On last iteration, multiply pr by degree and store the result */
            v.set_data(v.outc > 0 ? pr[v.id()] * v.outc : pr[v.id()]);
        }
    }
    
};

int main(int argc, const char ** argv) {
    graphchi_init(argc, argv);
    metrics m("pagerank");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 100000);
    scheduler          		= get_option_int("scheduler", false);;                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 50);
    epsilon	 				= get_option_float("epsilon", 0.001);
	//num_tasks_print			= get_option_int("print", false);
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
	assert(0 != nshards);
	num_vertices = get_num_vertices(filename);
	assert(num_vertices > 0);
	//array = (int*)malloc(sizeof(int)*num_vertices);
	//memset(array, 0, sizeof(int)*num_vertices);
    /* Run */
    graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
    
    bool inmemmode = false;//engine.num_vertices() * sizeof(EdgeDataType) < (size_t)engine.get_membudget_mb() * 1024L * 1024L;
    if (inmemmode) {
        logstream(LOG_INFO) << "Running Pagerank by holding vertices in-memory mode!" << std::endl;
        engine.set_modifies_outedges(false);
        engine.set_disable_outedges(true);
        engine.set_only_adjacency(true);
        PagerankProgramInmem program(engine.num_vertices());
        engine.run(program, niters);
    } else {
        PagerankProgram program;
        engine.run(program, niters);
    }
    
    /* Output top ranked vertices */
    std::vector< vertex_value<float> > top = get_top_vertices<float>(filename, ntop);
    std::cout << "Print top " << ntop << " vertices:" << std::endl;
    for(int i=0; i < (int)top.size(); i++) {
        std::cout << (i+1) << ". " << top[i].vertex << "\t" << top[i].value << std::endl;
    }
    
    metrics_report(m);    
	/*
	size_t total_updates = 0;
	for(int i=0; i<num_vertices; i++){
		total_updates += array[i];	
	}
	std::cout<<"num_vertices: "<<num_vertices<<"\ntotal updates: "<<total_updates<<"\naverage updates per vertex: "
		<<total_updates/num_vertices<<std::endl;
	*/
	FILE* fp = fopen((filename+".specific.csv").c_str(), "w+");
	assert(fp != NULL);	
	fprintf(fp, "pr_value1,pr_value8\n");
	for(int i=0; i<pr_value1.size(); i++){
		fprintf(fp, "%f,%f\n", pr_value1[i], pr_value8[i]);
	}
	fclose(fp);
	fp == NULL;
    return 0;
}

