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
struct VertexInfo{
	vid_t vid;
	int deg;
};

struct EdgeInfo{
	vid_t largelabel;
	vid_t smalllabel;   
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
typedef VertexInfo VertexDataType;
typedef EdgeInfo EdgeDataType;

mutex lock;
FILE* fp_edgelist = NULL;
FILE* fp_vt = NULL;
float epsilon = 0.001;
//int* array = NULL;
bool scheduler;
size_t num_vertices = 0;
//bool num_tasks_print = false;
struct SortProgram : public GraphChiProgram<VertexDataType, EdgeDataType> {
   	bool converged;
	bool interval_converged; 
	std::vector<VertexDataType> degvector;
	std::vector<VertexDataType> nonzerovector;
	std::vector<vid_t> idmap;
	SortProgram(int nvertices):degvector(nvertices), idmap(nvertices), nonzerovector(nvertices){
		/*
		degvector(nvertices);	
		idmap(nvertices);	
		*/
	}
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
			vid_t count = 1;
			std::sort(degvector.begin(), degvector.end(), sortFunc);		
			int idx = 0;
			//fprintf(fp_vt, "new_vid\told_vid\n");
			for(int i=0; i<degvector.size(); i++){
				//ignore zero degree vertex
				if(degvector[i].deg == 0) continue;
				//if(degvector[i].deg == 0) break;
				assert(idx < degvector.size());
				while(nonzerovector[idx].deg == 0) idx++;
				//idmap[degvector[i].vid] = count++;	
				idmap[degvector[i].vid] = idx++;	
				//fprintf(fp_vt, "%d\t%d\t%d\n", i, degvector[i].vid, degvector[i].deg);
				fprintf(fp_vt, "%d\t%d\t%d\n", idmap[degvector[i].vid], degvector[i].vid, degvector[i].deg);
			}
			fflush(fp_vt);		
			//fp_edgelist = fopen();
		}else if(iteration == 2){
			ginfo.set_last_iteration(iteration);	
		}
    }
    
    /**
      * Pagerank update function.
	  */
	void update(graphchi_vertex<VertexDataType, EdgeDataType> &v, graphchi_context &ginfo) {
		//array[v.id()]++;		
		if (ginfo.iteration == 0) {
			VertexDataType vdata = v.get_data();
			vdata.vid = v.id();
			vdata.deg = (int)v.num_edges();
			v.set_data(vdata);
			//assert(v.id() < degvector.size());
			degvector[v.id()] = vdata;
			nonzerovector[v.id()] = vdata;
		} else if(ginfo.iteration == 1){
			//VertexDataType vdata = v.get_data();
			for(int i=0; i<v.num_edges(); i++){
				graphchi_edge<EdgeDataType> * edge = v.edge(i);
				EdgeDataType edata = edge->get_data();
				edata.my_label(v.id(), edge->vertex_id()) = idmap[v.id()];
				edge->set_data(edata);	
			}
		}else{	
			for(int i=0; i<v.num_outedges(); i++){
				graphchi_edge<EdgeDataType> * edge = v.outedge(i);
				EdgeDataType edata = edge->get_data();
				vid_t mylabel = edata.my_label(v.id(), edge->vertex_id());
				vid_t nblabel = edata.nb_label(v.id(), edge->vertex_id());	
				assert(mylabel != nblabel);
				lock.lock();	
				fprintf(fp_edgelist, "%d\t%d\n", mylabel, nblabel);
				lock.unlock();

			}	
			/*
			lock.lock();
			fprintf(fp_vt, "%d\t%d\t%d\n", v.id(), idmap[v.id()], v.num_edges());
			lock.unlock();	
			*/
		}
	}
};

/**
  * Faster version of pagerank which holds vertices in memory. Used only if the number
  * of vertices is small enough.
  */
/*
struct SortProgramInmem : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    std::vector<EdgeDataType> pr;
    SortProgramInmem(int nvertices) :   pr(nvertices, RANDOMRESETPROB) {}
    
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
            v.set_data(v.outc > 0 ? pr[v.id()] * v.outc : pr[v.id()]);
        }
    }
    
};
*/
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
    graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
    //engine.set_modifies_inedges(false); // Improves I/O performance.
   	fp_edgelist = fopen((filename+".sdegree").c_str(), "w+");	 
   	fp_vt = fopen((filename+".sd.vt").c_str(), "w+");	 
	assert(fp_edgelist != NULL);
	assert(fp_vt != NULL);
	fprintf(fp_vt, "new_vid\told_vid\tdegree\n");
    bool inmemmode = false;//engine.num_vertices() * sizeof(EdgeDataType) < (size_t)engine.get_membudget_mb() * 1024L * 1024L;
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
        SortProgram program(num_vertices);
        engine.run(program, niters);
    }
   	fclose(fp_edgelist); 
	fclose(fp_vt);
	fp_vt = NULL;
	fp_edgelist = NULL;
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
    return 0;
}

