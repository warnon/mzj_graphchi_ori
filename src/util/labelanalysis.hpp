
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
 * Analyses output of label propagation algorithms such as connected components
 * and community detection. Memory efficient implementation.
 *
 * @author Aapo Kyrola
 */


#include <vector>
#include <algorithm>
#include <errno.h>
#include <assert.h>
#include <fstream>
#include <iostream>
//#include<fiostream.h>

#include "io/stripedio.hpp"
#include "logger/logger.hpp"
#include "util/merge.hpp"
#include "util/ioutil.hpp"
#include "util/qsort.hpp"
#include "api/chifilenames.hpp"
#include "engine/auxdata/vertex_data.hpp"

#ifndef DEF_GRAPHCHI_LABELANALYSIS
#define DEF_GRAPHCHI_LABELANALYSIS

using namespace graphchi;

template <typename LabelType>
struct labelcount_tt {
    LabelType label;
    unsigned int count;  // Count excludes the vertex which has its own id as the label. (Important optimization)
    labelcount_tt(LabelType l, int c) : label(l), count(c) {}
    labelcount_tt() {}
};

template <typename LabelType>
bool label_count_greater(const labelcount_tt<LabelType> &a, const labelcount_tt<LabelType> &b) {
    return a.count > b.count;
}

template <typename LabelType>
//void analyze_labels(std::string basefilename, int printtop = 20) {    
int analyze_labels(std::string basefilename, int printtop = 20) {    
    typedef labelcount_tt<LabelType> labelcount_t;
    /**
     * NOTE: this implementation is quite a mouthful. Cleaner implementation
     * could be done by using a map implementation. But STL map takes too much
     * memory, and I want to avoid Boost dependency - which would have boost::unordered_map.
     */
    metrics m("labelanalysis");
    stripedio * iomgr = new stripedio(m);
    
    /* Initialize the vertex-data reader */
    vid_t readwindow = 1024 * 1024;
    vid_t numvertices = (vid_t) get_num_vertices(basefilename);
    vertex_data_store<LabelType> * vertexdata =
    new vertex_data_store<LabelType>(basefilename, numvertices, iomgr);
    
    std::vector<labelcount_t> curlabels;
    bool first = true;
    vid_t curvid = 0;
    LabelType * buffer = (LabelType*) calloc(readwindow, sizeof(LabelType));
    
    /* Iterate the vertex values and maintain the top-list */
    vid_t st = 0;
    vid_t en = numvertices - 1;
    
    while(st <= numvertices - 1) {
        en = st + readwindow - 1;
        if (en >= numvertices - 1) en = numvertices - 1;
        
        /* Load the vertex values */
        vertexdata->load(st, en);
        
        int nt = en - st + 1;
        
        /* Mark vertices with its own label with 0xffffffff so they will be ignored */
        for(int i=0; i < nt; i++) { 
            LabelType l = *vertexdata->vertex_data_ptr(i + st);
            if (l == curvid) buffer[i] = 0xffffffff;
            else buffer[i] = l;
			
            curvid++;
        }
        
        /* First sort the buffer */
        quickSort(buffer, nt, std::less<LabelType>());
        
        /* Then collect */
        std::vector<labelcount_t> newlabels;
        newlabels.reserve(nt);
        LabelType lastlabel = LabelType(0xffffffff);
        for(int i=0; i < nt; i++) {
            if (buffer[i] != LabelType(0xffffffff)) {
                if (buffer[i] != lastlabel) {
                    newlabels.push_back(labelcount_t(buffer[i], 1));
                } else {
                    newlabels[newlabels.size() - 1].count ++;
                }
                lastlabel = buffer[i];
            }
        }
        
        if (first) {
            for(int i=0; i < (int)newlabels.size(); i++) {
                curlabels.push_back(newlabels[i]);
            }
            
        } else {
            /* Merge current and new label counts */
            int cl = 0;
            int nl = 0;
            std::vector< labelcount_t > merged;
            merged.reserve(curlabels.size() + newlabels.size());
            while(cl < (int)curlabels.size() && nl < (int)newlabels.size()) {
                if (newlabels[nl].label == curlabels[cl].label) {
                    merged.push_back(labelcount_t(newlabels[nl].label, newlabels[nl].count + curlabels[cl].count));
                    nl++; cl++;
                } else {
                    if (newlabels[nl].label < curlabels[cl].label) {
                        merged.push_back(newlabels[nl]);
                        nl++;
                    } else {
                        merged.push_back(curlabels[cl]);
                        cl++;
                    }
                }
            }
            while(cl < (int)curlabels.size()) merged.push_back(curlabels[cl++]);
            while(nl < (int)newlabels.size()) merged.push_back(newlabels[nl++]);
            
            curlabels = merged;
        }
        
        first = false;
        st += readwindow;
    }
    
    /* Sort */
    std::sort(curlabels.begin(), curlabels.end(), label_count_greater<LabelType>);
 
    /* Write output file */
    std::string outname = basefilename + ".components";
    std::ofstream resf;
    resf.open(outname.c_str());
    if (resf == NULL) {
        logstream(LOG_ERROR) << "Could not write label outputfile : " << outname << std::endl;
        return 0;
    }
	int valid_count = 0;
    for(int i=0; i < (int) curlabels.size(); i++) {
        resf << curlabels[i].label << "," << curlabels[i].count + 1 << std::endl;
		// only count labels whose size greater than 1
		if(curlabels[i].count > 1)
			valid_count++;
    }
    resf.close();
    
    std::cout << "Total number of different labels (components/communities): " << curlabels.size() << std::endl;
    std::cout << "List of labels was written to file: " << outname << std::endl;
    
    for(int i=0; i < (int)std::min((size_t)printtop, curlabels.size()); i++) {
        std::cout << (i+1) << ". label: " << curlabels[i].label << ", size: " << curlabels[i].count + 1 << std::endl;
    }
   	/* 
	//add by mzj
	std::string degreename = basefilename+".csv";
	//FILE* degreef = fopen(degreename.c_str(), "w");
	std::ofstream degreef(degreename.c_str());
	//assert(NULL != degreef);
	//if(degreef.is_open()){
	for(int i=0; i < (int) curlabels.size(); i++){
		//fprintf(degreef, "%u,%u\n", curlabels[i].label, curlabels[i].count);
		degreef<<curlabels[i].label<<","<<curlabels[i].count<<std::endl;
	}
		degreef.close();
	*/
	//}
	//fclose(degreef);
	//end
	
    free(buffer);
    
    delete vertexdata;
    delete iomgr;

	//return curlabels.size();
	return valid_count;
}


/*
int get_digits(int x){
	int num_digitable[]={9, 99, 999, 9999, 99999,
						 999999, 9999999, 99999999, 999999999, };
}
*/
std::string itoa(int n){
//	int len = 0;
 const int max_size = std::numeric_limits<int>::digits10 + 1 /*sign*/ + 1 /*0-terminator*/;
   char buffer[max_size] = {0};
   sprintf(buffer, "%d", n);
   return std::string(buffer);
	//fprintf();
}

template <typename LabelType>
struct labelcount_T {
    int  label;
    int count;  // Count excludes the vertex which has its own id as the label. (Important optimization)
    labelcount_T(int la, int c) : label(la), count(c) {}
    //labelcount_tt(int la) : label(la), count(0) {}
    //labelcount_tt() {}
};

template <typename LabelType>
bool label_ct_greater(const labelcount_T<LabelType> &a, const labelcount_T<LabelType> &b) {
    return a.count > b.count;
}

template <typename LabelType>
bool label_greater(const labelcount_T<LabelType> &a, const labelcount_T<LabelType> &b) {
    return a.label > b.label;
}

template <typename LabelType>
bool label_less(const labelcount_T<LabelType> &a, const labelcount_T<LabelType> &b) {
    return a.label < b.label;
}



template <typename LabelType>
bool operator !=(const labelcount_T<LabelType> &a, const labelcount_T<LabelType> &b) {
    return a.label != b.label;
}


// imitate analyze_lables to count num of labels and select Roots
template <typename LabelType>
int count_labels(std::string basefilename, int ntop =20, int NP=0) {    
    typedef labelcount_T<LabelType> labelcount_t;
    /**
     * NOTE: this implementation is quite a mouthful. Cleaner implementation
     * could be done by using a map implementation. But STL map takes too much
     * memory, and I want to avoid Boost dependency - which would have boost::unordered_map.
     */
    metrics m("label_count");
    stripedio * iomgr = new stripedio(m);
    
    /* Initialize the vertex-data reader */
    vid_t readwindow = 1024 * 1024;
    vid_t numvertices = (vid_t) get_num_vertices(basefilename);
    vertex_data_store<LabelType> * vertexdata =
    new vertex_data_store<LabelType>(basefilename, numvertices, iomgr);
    
    std::vector<labelcount_t> curlabels;
    bool first = true;
    //vid_t curvid = 0;
    labelcount_t * buffer = (labelcount_t*) calloc(readwindow, sizeof(LabelType));
    
    /* Iterate the vertex values and maintain the top-list */
    vid_t st = 0;
    vid_t en = numvertices - 1;
    int negcount = 0;//sum these vertices comp_id negative.
    while(st <= numvertices - 1) {
        en = st + readwindow - 1;
        if (en >= numvertices - 1) en = numvertices - 1;
        
        /* Load the vertex values */
        vertexdata->load(st, en);
        
        int nt = en - st + 1;
        int ki=0;
        /* Mark vertices with its own label with 0xffffffff so they will be ignored */
        for(int i=0; i < nt; i++) { 
            LabelType lab = *vertexdata->vertex_data_ptr(i + st);
			//assert();
			if(!NP){
			if(lab.get_value() >= 0)
				buffer[ki++] = labelcount_t(lab.get_value(), 1);
			else
			{
				//buffer[i] = labelcount_t(lab.get_value());
				negcount++;	
			}
			}else{
			if(lab.get_value() <= 0)
				buffer[ki++] = labelcount_t(lab.get_value(), 1);
			else
			{
				//buffer[i] = labelcount_t(lab.get_value());
				//negcount++;	
			}

			}
        }
        nt = ki;
        /* First sort the buffer */
        quickSort(buffer, nt, label_less<LabelType>);// important modify!!        
        //quickSort(buffer, nt, std::less<labelcount_t>());

        /* Then collect */
        std::vector<labelcount_t> newlabels;
        newlabels.reserve(nt);
		// just for initialize
        labelcount_t lastlabel = labelcount_t(-2, 0);
        for(int i=0; i < nt; i++) {
           // if (buffer[i] != lastlabel) {
                if (buffer[i] != lastlabel) {
                    newlabels.push_back(buffer[i]);
                } else {
                    newlabels[newlabels.size() - 1].count ++;
                }
                lastlabel = buffer[i];
           // }
        }
        if (first) {
            for(int i=0; i < (int)newlabels.size(); i++) {
                curlabels.push_back(newlabels[i]);
            }
            
        } else {
            /* Merge current and new label counts */
            int cl = 0;
            int nl = 0;
            std::vector< labelcount_t > merged;
            merged.reserve(curlabels.size() + newlabels.size());
            while(cl < (int)curlabels.size() && nl < (int)newlabels.size()) {
                if (newlabels[nl].label == curlabels[cl].label) {
                    merged.push_back(labelcount_t(newlabels[nl].label, newlabels[nl].count + curlabels[cl].count));
                    nl++; cl++;
                } else {
                    if (newlabels[nl].label < curlabels[cl].label) {
                        merged.push_back(newlabels[nl]);
                        nl++;
                    } else {
                        merged.push_back(curlabels[cl]);
                        cl++;
                    }
                }
            }
            while(cl < (int)curlabels.size()) merged.push_back(curlabels[cl++]);
            while(nl < (int)newlabels.size()) merged.push_back(newlabels[nl++]);
            
            curlabels = merged;
        }
        
        first = false;
        st += readwindow;
    }
	
 
	//std::cout<<"after sorting ================================"<<std::endl;
    /* Sort */
    std::sort(curlabels.begin(), curlabels.end(), label_ct_greater<LabelType>);

	//std::cout<<"output end=============="<<std::endl;	
	//assert(false);    
    /* Write output file */
	/*
    std::string outname = basefilename + ".components";
    std::ofstream resf;
    resf.open(outname.c_str());
    if (resf == NULL) {
        logstream(LOG_ERROR) << "Could not write label outputfile : " << outname << std::endl;
        return;
    }
    for(int i=0; i < (int) curlabels.size(); i++) {
        resf << curlabels[i].label << "," << curlabels[i].count + 1 << std::endl;
    }
    resf.close();
    
    std::cout << "Total number of different labels (components/communities): " << curlabels.size() << std::endl;
    std::cout << "List of labels was written to file: " << outname << std::endl;
    */
   // for(int i=0; i < (int)std::min((size_t)printtop, curlabels.size()); i++) {
	
	//int index = 0;
	static int icount = 0;
	int cur_size =(int) curlabels.size();
	//total size should include size of unvisited vertices!
	cur_size += negcount;
	std::cout<<"total number of different labels is :"<<cur_size<<"\t num vertices:"<<numvertices<<"\t first c_num"<<curlabels[0].count<<std::endl;
		
	for(int i=0; i < (int)std::min((size_t)ntop, curlabels.size()); i++) {
        std::cout << (i+1) << ". label: " << curlabels[i].label << ", size: " << curlabels[i].count + 1 << std::endl;
    }
	
	/*
	for(int i=0; i<cur_size && index < Rsize; i++){
		LabelType label = curlabels[i].label;
		if(! label.confirm )
		{
			Rt[index++]=(int)label.get_value();
		}
    }
	*/
    std::string outname = basefilename + ".count";
	FILE* fc = NULL;
	fc = fopen(outname.c_str(), "a+");
	fprintf(fc, "iteration-%d \t %d\n", icount++, cur_size);
	fclose(fc);
    free(buffer);
    
    delete vertexdata;
    delete iomgr;
	return cur_size;	
}


template <typename LabelType>
void partition_labels(std::string basefilename, int printtop = 20) {    
    metrics m("label partition");
    stripedio * iomgr = new stripedio(m);
    
    /* Initialize the vertex-data reader */
    vid_t readwindow = 1024 * 1024;
    vid_t numvertices = (vid_t) get_num_vertices(basefilename);
    vertex_data_store<LabelType> * vertexdata =
    new vertex_data_store<LabelType>(basefilename, numvertices, iomgr);
    
    //std::vector<labelcount_t> curlabels;
    //bool first = true;
    //vid_t curvid = 0;
    //LabelType * buffer = (LabelType*) calloc(readwindow, sizeof(LabelType));
    
    /* Iterate the vertex values and maintain the top-list */
    vid_t st = 0;
    vid_t en = numvertices - 1;
    
    while(st <= numvertices - 1) {
        en = st + readwindow - 1;
        if (en >= numvertices - 1) en = numvertices - 1;
        
        /* Load the vertex values */
        vertexdata->load(st, en);
        
        int nt = en - st + 1;
        
        /* Mark vertices with its own label with 0xffffffff so they will be ignored */
		for(int i=0; i < nt; i++) { 
		LabelType* L; //=LabelType(13);
		L = vertexdata->vertex_data_ptr(i + st);
				L->comp_id =300;
		}
		vertexdata->save();
		std::cout<<"vertex id"<<st<<"\t comp_id"<<vertexdata->vertex_data_ptr(st)->comp_id<<std::endl;
        //first = false;
        st += readwindow;
    }
    //free(buffer);
    
    delete vertexdata;
    delete iomgr;
}
#endif


