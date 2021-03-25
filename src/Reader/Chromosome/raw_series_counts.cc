/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   raw_series_counts.cc
 * Author: thomas
 * 
 * Created on October 18, 2018, 3:24 PM
 */

#include "raw_series_counts.h"
#include "../../Logger/logger.h"

raw_series_counts::raw_series_counts() : total_lefts(0), total_rights(0), count(0), paired_count(0) {
}

raw_series_counts::~raw_series_counts() {
}


void raw_series_counts::add_other_max_min(raw_series_counts& rsc, rpos min, rpos max) {
    
//    logger::Instance()->info("Pre " + std::to_string(total_lefts) + " " + std::to_string(total_rights) +  "\n");
//    logger::Instance()->info("Add Min Max " + std::to_string(rsc.total_lefts) + " " + std::to_string(rsc.total_rights) +  "\n");
    
    for (std::map< rpos,rcount >::iterator li = rsc.lefts->begin(); li != rsc.lefts->end(); ++li) {  
        rpos left = std::min(std::max(li->first, min), max);
        
        std::map< rpos,rcount >::iterator hit = lefts->find(left);
        if (hit == lefts->end()) {
            lefts->insert(std::make_pair(left, li->second));
        } else {
            hit->second += li->second;
        }
    }
    total_lefts += rsc.total_lefts;
    
    for (std::map< rpos,rcount >::iterator ri = rsc.rights->begin(); ri != rsc.rights->end(); ++ri) {  
        rpos right = std::max(std::min(ri->first, max), min);
        
        std::map< rpos,rcount >::iterator hit = rights->find(right);
        if (hit == rights->end()) {
            rights->insert(std::make_pair(right, ri->second));
        } else {
            hit->second += ri->second;
        }
    }
    total_rights += rsc.total_rights;
    
    for (std::map< rpos,rcount >::iterator li = rsc.hole_starts->begin(); li != rsc.hole_starts->end(); ++li) {  
        rpos hs = std::max(std::min(li->first, max), min);
        
        std::map< rpos,rcount >::iterator hit = hole_starts->find(hs);
        if (hit == hole_starts->end()) {
            hole_starts->insert(std::make_pair(hs, li->second));
        } else {
            hit->second += li->second;
        }
    }
    for (std::map< rpos,rcount >::iterator li = rsc.hole_ends->begin(); li != rsc.hole_ends->end(); ++li) {  
        rpos hs = std::min(std::max(li->first, min), max);
        
        std::map< rpos,rcount >::iterator hit = hole_ends->find(hs);
        if (hit == hole_ends->end()) {
            hole_ends->insert(std::make_pair(hs, li->second));
        } else {
            hit->second += li->second;
        }
    } 

    count += rsc.count;
    paired_count += rsc.paired_count;
    
}
