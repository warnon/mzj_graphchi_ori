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

struct PR_Count{
int count;
int tasks;
PR_Count(){
	count = tasks = 0;
}
};

float epsilon = 0.001;
//int* array = NULL;
PR_Count* array = NULL;
int repeat_count = 0;
//int max_repeat = 999999999;
bool scheduler = false;
//bool num_tasks_print = false;
vid_t sum_tasks= 0;
vid_t curr_sum = 0;

struct PagerankProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
    /**
      * Called before an iteration starts. Not implemented.
      */
    void before_iteration(int iteration, graphchi_context &ginfo) {
			/*
			if(ginfo.iteration == 1){
				ginfo.scheduler->add_task_to_all();	
			}
			*/
			converged = iteration > 0;
    }
    
    /**
      * Called after an iteration has finished. Not implemented.
      */
    void after_iteration(int iteration, graphchi_context &ginfo) {
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
	}
    
  	bool repeat_updates(graphchi_context &gcontext){
		return false;
		//skip	
		if(gcontext.iteration == 0)
			return false;
		else{
			//if(repeat_count % 100 == 0 || interval_converged)
			if(repeat_count % 100 == 0 || interval_converged)
			logstream(LOG_INFO)<<"interval="<<gcontext.exec_interval<<" iter="<<gcontext.iteration
				<<"\trepeat_count= "<<repeat_count<<std::endl;
			repeat_count = (interval_converged) ? 0 : repeat_count+1; 
			//repeat_count = (interval_converged) ? 0 : repeat_count+1; 
			/*
			if(num_tasks_print){
				sum_tasks = interval_converged ? 0 : sum_tasks; 
			}
			*/	
			//return !interval_converged ;
			return !(interval_converged);
		}
	}  
    /**
      * Pagerank update function.
      */
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
		//if(scheduler) ginfo.scheduler->remove_tasks_now(v.id(), v.id());
        //float sum=0;
        if (ginfo.iteration == 0) {
            for(int i=0; i < v.num_outedges(); i++) {
                graphchi_edge<float> * edge = v.outedge(i);
                //edge->set_data( RANDOMRESETPROB/v.num_outedges());
                edge->set_data(1.0 / v.num_outedges());
            }
            //v.set_data(RANDOMRESETPROB); 
			//////
			if(v.num_edges() != 0 && v.num_inedges() == 0){
				//std::cout<<"vid= "<<v.id()<<"\t deg= "<<v.num_edges()<<std::endl;		
			}
			//////
			//schedule this vertex for next iteration
			//if(scheduler) ginfo.scheduler->add_task(v.id(), false);
        } else {
        	float sum=0;
			float old_value = v.get_data();
            for(int i=0; i < v.num_inedges(); i++) {
                sum += v.inedge(i)->get_data();
            }
            float pagerank = RANDOMRESETPROB + (1 - RANDOMRESETPROB) * sum;
            
            if (v.num_outedges() > 0) {
                float pagerankcont = pagerank / v.num_outedges();
                for(int i=0; i < v.num_outedges(); i++) {
                    graphchi_edge<float> * edge = v.outedge(i);
					edge->set_data(pagerankcont);
				}
            }
            ginfo.log_change(std::abs(pagerank - v.get_data()));
            
            v.set_data(pagerank); 
			if(std::abs(pagerank-old_value) > epsilon){
				if(converged)
					converged = false;
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
    metrics m("block-pagerank");
    global_logger().set_log_level(LOG_DEBUG);

    /* Parameters */
    std::string filename    = get_option_string("file"); // Base filename
    int niters              = get_option_int("niters", 100000);
    scheduler          		= get_option_int("scheduler", false);                    // Non-dynamic version of pagerank.
    int ntop                = get_option_int("top", 50);
    epsilon					= get_option_float("epsilon", 0.001);
	bool parallel			= get_option_int("parallel", 0);
	//num_tasks_print			= get_option_int("print", false);
	//max_repeat				= get_option_int("repeat", 999999999); // max number of repeat allowed for each subgraph
    /* Process input file - if not already preprocessed */
    int nshards             = convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
	assert(0 != nshards);
	//array = (int*)malloc(sizeof(int)*nshards);

    /* Run */
    graphchi_engine<float, float> engine(filename, nshards, scheduler, m); 
    engine.set_modifies_inedges(false); // Improves I/O performance.
	if(parallel){
		std::cout<<"full parallelism is enabled!!!"<<std::endl;
    	engine.set_enable_deterministic_parallelism(false);
	}
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
	if(	num_tasks_print ){
		array = (PR_Count*)malloc(sizeof(PR_Count)*nshards);
		assert(array != NULL);
		for(int i=0; i<nshards; i++){
			array[i].count = 0;
			array[i].tasks = 0;	
		}
		FILE* pr_out = NULL;
		if(scheduler)
			pr_out = fopen((filename+".schd_pr_ana").c_str(), "w+");
		else
			pr_out = fopen((filename+".pr_ana").c_str(), "w+");
		assert(pr_out != NULL);
		int avg_count = 0;
		int avg_tasks = 0;
		int nvertices = engine.num_vertices();
		for(int i=0; i<nshards; i++){
			fprintf(pr_out, "%d\t%d\t%d\n", i, array[i].count, array[i].tasks);		
			avg_count += array[i].count;
			avg_tasks += array[i].tasks;
		}		
		fprintf(pr_out, "avg\t%.4f\t%.4f\n", (float)avg_count/nshards, (float)avg_tasks/nvertices);	
		fclose(pr_out);
		free(array);
		printf("avg\t%.4f\t%.4f\n", (float)avg_count/nshards, (float)avg_tasks/nvertices);	
	}
	*/
    return 0;
}

