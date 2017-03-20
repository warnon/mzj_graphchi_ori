
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

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"

using namespace graphchi;

/**
  * Type definitions. Remember to create suitable graph shards using the
  * Sharder-program. 
  */
typedef int VertexDataType;
typedef int  EdgeDataType;


unsigned single_source = 0;
bool converged = false;
unsigned maxlevel = 100000;
bool scheduler = false;
mutex lock;
FILE* fpout = NULL;
int level = 3;//number of levels we want to dump
/**
  * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
  * class. The main logic is usually in the update function.
  */
struct bfs : public GraphChiProgram<VertexDataType, EdgeDataType> {
    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {

		if (gcontext.iteration == 0) {
			if(vertex.id() == single_source)
			{
				vertex.set_data(level);
				for(int id = 0; id < vertex.num_outedges(); id++)
				{
					if(scheduler)
					{
						gcontext.scheduler->add_task(vertex.outedge(id)->vertex_id());
					}
					//vertex.outedge(id)->set_data(vertex.get_data());
					vertex.outedge(id)->set_data(level-1);
					lock.lock();
					fprintf(fpout, "%d,%d\n", vertex.id(), vertex.outedge(id)->vertex_id());
					lock.unlock();
					//std::cout<<"vid"<<vertex.edge(id)->vertexid
				}
				//assert(vertex.num_edges() ==4 );	
			}
			else
			{
				vertex.set_data(0);
				for(int id=0; id < vertex.num_outedges(); id++){
					vertex.outedge(id)->set_data(0);	
				}
				
			}        

		} else {
			//VertexDataType tmpval = vertex.get_data();
			VertexDataType tmpval = vertex.get_data();
			for(int i=0; i < vertex.num_inedges(); i++) {
				// Do something
				//    value += vertex.inedge(i).get_data();
				//unsigned tmpval = vertex.inedge(i)->get_data() + 1;
				tmpval = std::max(tmpval, vertex.inedge(i)->get_data());
			}
			if(tmpval > vertex.get_data())
			{
				converged =false;
				vertex.set_data(tmpval);
				for(int i=0; i < vertex.num_outedges(); i++)
				{
					// Do something
					// vertex.outedge(i).set_data(x)
					vertex.outedge(i)->set_data(tmpval-1);
					vid_t eid = vertex.outedge(i)->vertex_id();
					lock.lock();
					fprintf(fpout, "%d,%d\n", vertex.id(), eid);	
					lock.unlock();
					if(scheduler)
					{
						if(eid <= gcontext.interval_en)
							gcontext.scheduler->add_task(eid, false);
						else
							gcontext.scheduler->add_task(eid, true);
					}
				}
			}

		}
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
		if(!scheduler){
			if(converged)
			{
				logstream(LOG_INFO)<<"bfs program has converged"<<std::endl;
				gcontext.set_last_iteration(iteration);
				fclose(fpout);
			}
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


int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
       arguments and the configuration file. */
    graphchi_init(argc, argv);
    
    /* Metrics object for keeping track of performance counters
       and other information. Currently required. */
    metrics m("BFS-level");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 100000); // Number of iterations
    scheduler       	 = get_option_int("scheduler", false); // Whether to use selective scheduling
   	single_source 		 = get_option_int("root", 0); 
   	level 				 = get_option_int("level", 3); 
    /* Detect the number of shards or preprocess an input to create them */
    int nshards          = convert_if_notexists<EdgeDataType>(filename, 
                                                            get_option_string("nshards", "auto"));
   
	fpout = fopen((filename+".level").c_str(), "w+"); 
	assert(fpout != NULL);
    /* Run */
    bfs program;
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    engine.run(program, niters);
   /*analyze result*/
	m.start_time("label-analysis");
	analyze_labels<unsigned>(filename, 40); 
    /* Report execution metrics */
    metrics_report(m);
    return 0;
}
