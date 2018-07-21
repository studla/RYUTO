/* 
 * File:   count_raw_node.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 15, 2017, 9:47 AM
 */

#include <map>
#include "count_raw_node.h"
#include "../../Datatype_Templates/move.h"

count_raw_node::count_raw_node() : total_lefts(0), total_rights(0) {
}

count_raw_node::~count_raw_node() {
}

void count_raw_node::add_node(count_raw_node* node) {
    
    for (std::map< rpos,rcount >::iterator li = node->lefts.begin(); li != node->lefts.end(); ++li) {  
        
        std::map< rpos,rcount >::iterator hit = lefts.find(li->first);
        if (hit == lefts.end()) {
            lefts.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator ri = node->rights.begin(); ri != node->rights.end(); ++ri) {      
        std::map< rpos,rcount >::iterator hit = rights.find(ri->first);
        if (hit == rights.end()) {
            rights.insert(*ri);
        } else {
            hit->second += ri->second;
        }
    }
    total_lefts += node->total_lefts;
    total_rights += node->total_rights;
}

void count_raw_node::add_node_start(exon_group* exons) {
    // copy evidence into lefts

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Node Start at. " + exons->bin_mask.to_string() + "\n");
    
    rcount ss = 0;
    for (std::map< rpos,rcount >::iterator it = lefts.begin(); it != lefts.end(); ++it) {
        ss += it->second;
    }
    rcount ee = 0;
    for (std::map< rpos,rcount >::iterator it = rights.begin(); it != rights.end(); ++it) {
        ee += it->second;
    }
   
   logger::Instance()->debug("Start AT  " + std::to_string(total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(total_lefts) + "\n");
   #endif 
   
    for (std::map< rpos,rcount >::iterator li = exons->lefts->begin(); li != exons->lefts->end(); ++li) {  
        
        std::map< rpos,rcount >::iterator hit = lefts.find(li->first);
        if (hit == lefts.end()) {
            lefts.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_ends[0].begin(); li != exons->hole_ends[0].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = lefts.find(li->first);
        if (hit == lefts.end()) {
            lefts.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_starts[0].begin(); li != exons->hole_starts[0].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = rights.find(li->first);
        if (hit == rights.end()) {
            rights.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    
    total_lefts += exons->total_lefts + exons->hole_end_counts[0] - exons->hole_start_counts[0];

    if (exons->bin_mask.id.count() == 1) { // just one long, lefts are rights
        for (std::map< rpos,rcount >::iterator li = exons->rights->begin(); li != exons->rights->end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = rights.find(li->first);
            if (hit == rights.end()) {
                rights.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
         // there can't be an imbalance for single exon evidence
        total_lefts -= exons->total_rights;
    }
    
    #ifdef ALLOW_DEBUG
     ss = 0;
    for (std::map< rpos,rcount >::iterator it = lefts.begin(); it != lefts.end(); ++it) {
        ss += it->second;
    }
     ee = 0;
    for (std::map< rpos,rcount >::iterator it = rights.begin(); it != rights.end(); ++it) {
        ee += it->second;
    }
   
   logger::Instance()->debug("Set TO  " + std::to_string(total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(total_lefts) + "\n");
   #endif
    
}


void count_raw_node::add_node_end(exon_group* exons) {
    
    unsigned int last_i = exons->set_exons -1;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Node END at. " + exons->bin_mask.to_string() + " " + std::to_string(last_i) + "\n");
    logger::Instance()->debug("Start AT  " + std::to_string(total_rights) + " " + std::to_string(total_lefts) + "\n");
    #endif
    
    // copy evidence into lefts
    for (std::map< rpos,rcount >::iterator li = exons->rights->begin(); li != exons->rights->end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = rights.find(li->first);
        if (hit == rights.end()) {
            rights.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_ends[last_i].begin(); li != exons->hole_ends[last_i].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = lefts.find(li->first);
        if (hit == lefts.end()) {
            lefts.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_starts[last_i].begin(); li != exons->hole_starts[last_i].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = rights.find(li->first);
        if (hit == rights.end()) {
            rights.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
        
    total_rights += exons->total_rights + exons->hole_start_counts[last_i] - exons->hole_end_counts[last_i];
    
    #ifdef ALLOW_DEBUG
    rcount ss = 0;
    for (std::map< rpos,rcount >::iterator it = lefts.begin(); it != lefts.end(); ++it) {
        ss += it->second;
    }
    rcount ee = 0;
    for (std::map< rpos,rcount >::iterator it = rights.begin(); it != rights.end(); ++it) {
        ee += it->second;
    }
   
   logger::Instance()->debug("Set TO  " + std::to_string(total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(total_lefts) + "\n");
   #endif
 
}


rcount count_raw_node::add_node_initial_index(exon_group* exons, unsigned int i, rcount split_v) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Node Middle at. " + exons->bin_mask.to_string() + " " + std::to_string(i) + " " + std::to_string(split_v) + "\n");
    logger::Instance()->debug("Start AT  " + std::to_string(total_rights) + " " + std::to_string(total_lefts) + "\n");
    #endif
    
    total_rights += split_v;
    split_v = split_v + exons->hole_end_counts[i] - exons->hole_start_counts[i];
    total_lefts += split_v;
    
    for (std::map< rpos,rcount >::iterator li = exons->hole_ends[i].begin(); li != exons->hole_ends[i].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = lefts.find(li->first);
        if (hit == lefts.end()) {
            lefts.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = exons->hole_starts[i].begin(); li != exons->hole_starts[i].end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = rights.find(li->first);
        if (hit == rights.end()) {
            rights.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    
    #ifdef ALLOW_DEBUG
    rcount ss = 0;
    for (std::map< rpos,rcount >::iterator it = lefts.begin(); it != lefts.end(); ++it) {
        ss += it->second;
    }
    rcount ee = 0;
    for (std::map< rpos,rcount >::iterator it = rights.begin(); it != rights.end(); ++it) {
        ee += it->second;
    }
   
   logger::Instance()->debug("Set TO  " + std::to_string(total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(total_lefts) + "\n");
   #endif
    
    return split_v;
}

void count_raw_node::add_node_index(count_raw_edge& edge, unsigned int i) {
    
  #ifdef ALLOW_DEBUG  
  logger::Instance()->debug("Split NODE from Existing edge " + std::to_string(i) + "\n");
  logger::Instance()->debug("Start AT  " + std::to_string(total_rights) + " " + std::to_string(total_lefts) + "\n");  
  #endif
  
   total_rights += edge.splits[i];
   total_lefts += edge.splits[i+1];

   rcount count_l = 0;
   rcount count_r = 0;
   
   // merge the maps!
   for (std::map< rpos,rcount >::iterator li = edge.starts[i]->begin(); li != edge.starts[i]->end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = lefts.find(li->first);
        count_l += li->second;
        if (hit == lefts.end()) {
            lefts.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = edge.ends[i]->begin(); li != edge.ends[i]->end(); ++li) {      
        std::map< rpos,rcount >::iterator hit = rights.find(li->first);
         count_r += li->second;
        if (hit == rights.end()) {
            rights.insert(*li);
        } else {
            hit->second += li->second;
        }
    }
    
   
   #ifdef ALLOW_DEBUG
    rcount ss = 0;
    for (std::map< rpos,rcount >::iterator it = lefts.begin(); it != lefts.end(); ++it) {
        ss += it->second;
    }
    rcount ee = 0;
    for (std::map< rpos,rcount >::iterator it = rights.begin(); it != rights.end(); ++it) {
        ee += it->second;
    }
   
   logger::Instance()->debug("Set TO  " + std::to_string(total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(total_lefts) + "\n");
   #endif
}


capacity_type count_raw_node::get_max() {
    

    capacity_type count = total_rights;
    capacity_type max = count; // the initial value
        
    std::map< rpos,rcount >::iterator li = lefts.begin();
    std::map< rpos,rcount >::iterator ri = rights.begin();

    while (li != lefts.end() && ri != rights.end()) {

        if ( li == lefts.end() || (ri != rights.end() && ri->first + 1 < li->first) ) {
            count -= ri->second;
            ++ri;
        } else if (ri == rights.end() || ri->first +1 > li->first) {
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
    
    return max;
}

capacity_type count_raw_node::get_min() {
    

    capacity_type count = total_rights;
    capacity_type min = count; // the initial value
        
    std::map< rpos,rcount >::iterator li = lefts.begin();
    std::map< rpos,rcount >::iterator ri = rights.begin();

    while (li != lefts.end() && ri != rights.end()) {

        if ( li == lefts.end() || (ri != rights.end() && ri->first + 1 < li->first) ) {
            count -= ri->second;
            ++ri;
        } else if (ri == rights.end() || ri->first +1 > li->first) {
            count += li->second;
            ++li;
        } else {
           count -= ri->second;
           count += li->second;
           ++ri;
           ++li;
        }

        if (count < min) {
            min = count;
        }
    }
    
    return min;
}
