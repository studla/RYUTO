/* 
 * File:   node_count_raw.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 15, 2015, 9:59 AM
 */

#include <vector>

#include "count_raw_edge.h"
#include "../../Logger/logger.h"
#include "../../Datatype_Templates/move.h"

count_raw_edge::count_raw_edge() : size(0), initialized(false){
}

count_raw_edge::count_raw_edge(unsigned int size) : size(size), initialized(false) {
    
    starts.assign(size, {});
    ends.assign(size, {});
}

count_raw_edge::~count_raw_edge() {
}

void count_raw_edge::add_edge(count_raw_edge* edge) {
    
    splits[0] += edge->splits[0];
    for (unsigned int i = 0; i < size; ++i) {

        for (std::map< rpos,rcount >::iterator li = edge->starts[i]->begin(); li != edge->starts[i]->end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = starts[i]->find(li->first);
            if (hit == starts[i]->end()) {
                starts[i]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = edge->ends[i]->begin(); li !=  edge->ends[i]->end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = ends[i]->find(li->first);
            if (hit == ends[i]->end()) {
                ends[i]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }

        splits[i+1] += edge->splits[i+1];
    }        
}


void count_raw_edge::add_initial_count(exon_group* exons) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Create Initial. " + exons->bin_mask.to_string() + " " + std::to_string(exons->size) + " " + std::to_string(exons->set_exons) + "\n");
    #endif
    
    initialized = true;
    
    size = exons->set_exons - 2;
    
    starts.assign(size, {});
    ends.assign(size, {});
    
    splits.assign(size + 1, 0);
    
    splits[0] = exons->total_lefts + exons->hole_end_counts[0] - exons->hole_start_counts[0];
    if (size != 0) {
        splits[size] = exons->total_rights;
    }
    for (unsigned int i = 1; i < exons->set_exons - 1; ++i) {
        _MOVE_RANGE(exons->hole_starts[i].begin(), exons->hole_starts[i].end(), std::inserter(ends[i-1].ref(), ends[i-1]->begin()));
        _MOVE_RANGE(exons->hole_ends[i].begin(), exons->hole_ends[i].end(), std::inserter(starts[i-1].ref(), starts[i-1]->begin()));
        splits[i] = splits[i-1] + exons->hole_end_counts[i] - exons->hole_start_counts[i];
    }
     
    #ifdef ALLOW_DEBUG
    for (unsigned int i = 0; i < size; ++i) {
       
        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
        
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
            ee += it->second;
        }
        
        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
    }
    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
    #endif
}


rcount count_raw_edge::add_sub_counts_eg(exon_group* exons, int from, int to, rcount split_v) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Subcount RAW. " + exons->bin_mask.to_string() + " " + std::to_string(exons->set_exons) + " " + std::to_string(from) + " " + std::to_string(to) + "\n");
    #endif
    
    if (!initialized) { // unitialized arc
        
        initialized = true;
        
        size = to - from + 1;
        
        starts.assign(size, {});
        ends.assign(size, {});

        splits.assign(size + 1, 0);
    
        splits[0] += split_v;
        unsigned int g = 1;
        for (unsigned int i = from; i <= to; ++i, ++g) {

            std::copy(exons->hole_starts[i].begin(), exons->hole_starts[i].end(), std::inserter(ends[g-1].ref(), ends[g-1]->begin()));
            std::copy(exons->hole_ends[i].begin(), exons->hole_ends[i].end(), std::inserter(starts[g-1].ref(), starts[g-1]->begin()));
            split_v = split_v + exons->hole_end_counts[i] - exons->hole_start_counts[i];
            splits[g] += split_v;
        } 
    } else {
        
        splits[0] += split_v;
        unsigned int g = 1;
        for (unsigned int i = from; i <= to; ++i, ++g) {

            for (std::map< rpos,rcount >::iterator li = exons->hole_ends[i].begin(); li != exons->hole_ends[i].end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = starts[g-1]->find(li->first);
                if (hit == starts[g-1]->end()) {
                    starts[g-1]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }
            for (std::map< rpos,rcount >::iterator li = exons->hole_starts[i].begin(); li != exons->hole_starts[i].end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = ends[g-1]->find(li->first);
                if (hit == ends[g-1]->end()) {
                    ends[g-1]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }
            
            split_v = split_v + exons->hole_end_counts[i] - exons->hole_start_counts[i];
            splits[g] += split_v;
        }        
    }
    
    
    #ifdef ALLOW_DEBUG
    for (unsigned int i = 0; i < size; ++i) {
       
        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
        
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
            ee += it->second;
        }
        
        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
    }
    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
    #endif
    
    return split_v;
}


rcount count_raw_edge::add_sub_counts_start(exon_group* exons, int g) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcount start. " + exons->bin_mask.to_string() + " " + std::to_string(g) + "\n");
    #endif
    
    for (std::map< rpos,rcount >::iterator li = exons->lefts->begin(); li != exons->lefts->end(); ++li) { 
        std::map< rpos,rcount >::iterator hit = starts[g]->find(li->first);
        if (hit == starts[g]->end()) {
            starts[g]->insert(*li);
        } else {
            hit->second += li->second;
        }
    }    
    for (std::map< rpos,rcount >::iterator li = exons->hole_ends[0].begin(); li != exons->hole_ends[0].end(); ++li) { 
                
        std::map< rpos,rcount >::iterator hit = starts[g]->find(li->first);
        if (hit == starts[g]->end()) {
            starts[g]->insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_starts[0].begin(); li != exons->hole_starts[0].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = ends[g]->find(li->first);
        if (hit == ends[g]->end()) {
            ends[g]->insert((*li));
        } else {
            hit->second += li->second;
        }
    }

    splits[g+1] += exons->total_lefts + exons->hole_end_counts[0] - exons->hole_start_counts[0];

    if (exons->bin_mask.id.count() == 1) { // just one long, lefts are rights
        for (std::map< rpos,rcount >::iterator li = exons->rights->begin(); li != exons->rights->end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = ends[g]->find(li->first);
            if (hit == ends[g]->end()) {
                ends[g]->insert(*li);
            } else {
                hit->second += li->second;
            }
        }
         // there can't be an imbalance for single exon evidence
        splits[g+1] -= exons->total_rights;
    }
    
    #ifdef ALLOW_DEBUG
    for (unsigned int i = 0; i < size; ++i) {
       
        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
        
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
            ee += it->second;
        }
        
        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
    }
    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
    #endif

    return exons->total_lefts + exons->hole_end_counts[0] - exons->hole_start_counts[0];
}


rcount count_raw_edge::add_sub_counts_end(exon_group* exons, int g) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcount end. " + exons->bin_mask.to_string() + " " + std::to_string(g) + "\n");
    #endif
    
    unsigned int last_i = exons->set_exons -1;
    
    for (std::map< rpos,rcount >::iterator li = exons->rights->begin(); li != exons->rights->end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = ends[g]->find(li->first);
        if (hit == ends[g]->end()) {
            ends[g]->insert(*li);
        } else {
            hit->second += li->second;
        }
    }    
    for (std::map< rpos,rcount >::iterator li = exons->hole_ends[last_i].begin(); li != exons->hole_ends[last_i].end(); ++li) { 
        std::map< rpos,rcount >::iterator hit = starts[g]->find(li->first);
        if (hit == starts[g]->end()) {
            starts[g]->insert((*li));
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_starts[last_i].begin(); li != exons->hole_starts[last_i].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = ends[g]->find(li->first);
        if (hit == ends[g]->end()) {
            ends[g]->insert((*li));
        } else {
            hit->second += li->second;
        }
    }
    
 //   splits[g] += exons->total_rights + exons->hole_start_counts[last_i] - exons->hole_end_counts[last_i]; by range!

    #ifdef ALLOW_DEBUG
    for (unsigned int i = 0; i < size; ++i) {
       
        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
        
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
            ee += it->second;
        }
        
        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
    }
    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
    #endif

    return splits[g+1];
}

rcount count_raw_edge::add_sub_counts_eg_range(exon_group* exons, int from, int to, unsigned int g1, rcount split_v) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcount raw range. " + exons->bin_mask.to_string() + " " + std::to_string(split_v) + " " + std::to_string(from) + " " + std::to_string(to) + " " + std::to_string(g1) + "\n");
    #endif
    
    if (to < from) {
        return split_v;
    }
    
    //splits[g1-1] += split_v; we explicitly don't do that here! set externally!
    unsigned int g = g1;
    for (unsigned int i = from; i <= to; ++i, ++g) {
                
        for (std::map< rpos,rcount >::iterator li = exons->hole_ends[i].begin(); li != exons->hole_ends[i].end(); ++li) {  
            std::map< rpos,rcount >::iterator hit = starts[g-1]->find(li->first);
            if (hit == starts[g-1]->end()) {
                starts[g-1]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }

        for (std::map< rpos,rcount >::iterator li = exons->hole_starts[i].begin(); li != exons->hole_starts[i].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = ends[g-1]->find(li->first);
            if (hit == ends[g-1]->end()) {
                ends[g-1]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }
        
        split_v = split_v + exons->hole_end_counts[i] - exons->hole_start_counts[i];
        splits[g] += split_v;
    } 
    
    
    #ifdef ALLOW_DEBUG
    for (unsigned int i = 0; i < size; ++i) {
       
        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
        
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
            ee += it->second;
        }
        
        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
    }
    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
    #endif
    
    return split_v;
}

void count_raw_edge::add_sub_counts_re(count_raw_edge &a, int from, int to) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcounts existing. " + std::to_string(a.size) + " " + std::to_string(from) + " " + std::to_string(to) + "\n");
    #endif
    
    if (from > to) {
        
        if (!initialized) {
            initialized = true;
            //just one split site remains, after that un-subable
            splits.push_back(a.splits[from]);
        } else {
            splits[0] += a.splits[from];
        }
        return;
    }
    
    if (!initialized) { // copy in
        
        initialized = true;
        
        size = to - from + 1;
        starts.reserve(size);
        ends.reserve(size);
        splits.reserve(size + 1);

        std::copy(a.starts.begin() + from, a.starts.begin() + to + 1, std::back_inserter(starts));
        std::copy(a.ends.begin() + from, a.ends.begin() + to + 1, std::back_inserter(ends));
        std::copy(a.splits.begin() + from, a.splits.begin() + to + 2, std::back_inserter(splits));
                
    } else { // merge with existing!
                
        unsigned g = 0;
        unsigned int i = from;
        for (; i <= to; ++i, ++g) {
            splits[g] += a.splits[i]; 
                        
            for (std::map< rpos,rcount >::iterator li = a.starts[i]->begin(); li != a.starts[i]->end(); ++li) { 
                std::map< rpos,rcount >::iterator hit = starts[g]->find(li->first);
                if (hit == starts[g]->end()) {
                    starts[g]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }
            for (std::map< rpos,rcount >::iterator li = a.ends[i]->begin(); li != a.ends[i]->end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = ends[g]->find(li->first);
                if (hit == ends[g]->end()) {
                    ends[g]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }
        }
        splits[g] += a.splits[i];
    }
    
    #ifdef ALLOW_DEBUG
    for (unsigned int i = 0; i < size; ++i) {
       
        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
        
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
            ee += it->second;
        }
        
        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
    }
    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
    #endif
}

capacity_type count_raw_edge::get_max() {
    
    if (!initialized) {
        return 0;
    }
    
    capacity_type count = splits[0];
    capacity_type max = count; // the initial value
        
    for (unsigned int i = 0; i < size; ++i) {
        
        std::map< rpos,rcount >::iterator li = starts[i]->begin();
        std::map< rpos,rcount >::iterator ri = ends[i]->begin();
           
        while (li != starts[i]->end() && ri != ends[i]->end()) {
            
            if ( li == starts[i]->end() || (ri != ends[i]->end() && ri->first + 1 < li->first) ) {
                count -= ri->second;
                ++ri;
            } else if (ri == ends[i]->end() || ri->first +1 > li->first) {
                count += li->second;
                ++li;
            } else {
               count -= ri->second;
               count += li->second;
               ++ri;
               ++li;
            }
            
            if (count > max) {
                max = count;
            }
        }
    }
    return max;
}

capacity_type count_raw_edge::get_min() {
    
    if (!initialized) {
        return 0;
    }
    
    capacity_type count = splits[0];
    capacity_type min = count; // the initial value
        
    for (unsigned int i = 0; i < size; ++i) {
        
        std::map< rpos,rcount >::iterator li = starts[i]->begin();
        std::map< rpos,rcount >::iterator ri = ends[i]->begin();
           
        while (li != starts[i]->end() && ri != ends[i]->end()) {
            
            if ( li == starts[i]->end() || (ri != ends[i]->end() && ri->first + 1 < li->first) ) {
                count -= ri->second;
                ++ri;
            } else if (ri == ends[i]->end() || ri->first +1 > li->first) {
                count += li->second;
                ++li;
            } else {
               count -= ri->second;
               count += li->second;
               ++ri;
               ++li;
            }
            
            if (count < min ) {
                min = count;
            }
        }
    }
    return min;
}
