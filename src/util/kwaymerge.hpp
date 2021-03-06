
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
 * Generic k-way merge. Could reuse existing solutions, but as a graduate student I reserve the
 * right to do my own implementations for the sake of it :).
 */

#ifndef DEF_KWAYMERGE_GRAPHCHI
#define DEF_KWAYMERGE_GRAPHCHI

#include <assert.h>
#include <stdlib.h>

#include <vector>
#include "binary_minheap.hpp"


template <typename T>
class merge_source {
public:
    virtual bool has_more() = 0;
    virtual T next() = 0;
};
template <typename T>
class merge_sink {
public:
    virtual void add(T val) = 0;
    virtual void done() = 0;
};

template <typename T>
struct value_source {
    int sourceidx;
    T value;
    value_source(int sourceidx, T value) : sourceidx(sourceidx), value(value) {}
    
    bool operator< (value_source &x2)
    {
        return value < x2.value;
    }
	//T operator ()
};


template <typename T, class F>
class kway_merge {
    std::vector<merge_source<T> *> sources;
    merge_sink<T> * sink;
    int K;
    binary_minheap<value_source<T>, F> tip;

public:
    kway_merge(std::vector<merge_source<T> *> sources, merge_sink<T> * sink, F f): sources(sources), sink(sink), tip((int)sources.size(), f) {
        K = (int) sources.size();
    }
    
    ~kway_merge() {
        sink = NULL;
    }
    
    void merge() {
        int active_sources =(int)sources.size();
        
        for(int i=0; i<active_sources; i++) {
            tip.insert(value_source<T>(i, sources[i]->next()));
        }
        
        while(active_sources > 0 || !tip.empty()) {
            value_source<T> vv = tip.min();
            tip.extractMin();
            if (sources[vv.sourceidx]->has_more()) {
                tip.insert(value_source<T>(vv.sourceidx, sources[vv.sourceidx]->next()));
			//	logstream(LOG_INFO)<<"in merge if sec sourceid:"<<vv.sourceidx<<std::endl;
            } else {
                active_sources--;
            }
            sink->add(vv.value);
        }
        sink->done();
    }
};
template <typename T, class F>
class myway_merge {
    std::vector<merge_source<T> *> sources;
   // std::vector<merge_source<T> *> sources2;
    //merge_sink<T> * sink;
    int K;
	//int K2;
	int active_sources;
	//int active_sources2;
    binary_minheap<value_source<T>, F> tip;
	

public:
    myway_merge(std::vector<merge_source<T> *> sources1_, F f)
	: sources(sources1_), tip((int)sources1_.size(), f){
        K = (int) sources.size();
       // K2 = (int) sources2.size();
		construct_heap();
    }
	
	T*  getedge( ){
		T* ptr = NULL; 
        if(active_sources > 0 || !tip.empty()) {
            value_source<T> vv = tip.min();
            tip.extractMin();
            if (sources[vv.sourceidx]->has_more()) {
                tip.insert(value_source<T>(vv.sourceidx, sources[vv.sourceidx]->next()));
            } else {
                active_sources--;
            }
          	ptr =  &(vv.value);
        }
			return ptr;
	} 

	
	bool getedge2(T& edge){
			bool ptr = false; 
        if(active_sources > 0 || !tip.empty()) {
            value_source<T> vv = tip.min();
            tip.extractMin();
            if (sources[vv.sourceidx]->has_more()) {
                tip.insert(value_source<T>(vv.sourceidx, sources[vv.sourceidx]->next()));
            } else {
                active_sources--;
            }
          	edge =  (vv.value);
			ptr = true;
        }
			return ptr;

	}


 	void construct_heap(){
		
        active_sources =(int)sources.size();
      //  active_sources2 =(int)sources2.size();
        for(int i=0; i<K; i++) {
            tip.insert(value_source<T>(i, sources[i]->next()));
		}  
	/*	for(int j=0; j<K2; j++){
			 tip2.insert(value_source<T>(j, sources[j]->next()));
		}
	*/
     }

};
#endif











