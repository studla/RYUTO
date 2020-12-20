/* 
 * File:   count_raw_node.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on May 15, 2017, 9:47 AM
 */

#include <map>
#include "count_raw_node.h"
#include "../../../Datatype_Templates/move.h"

count_raw_node::count_raw_node() {
}

count_raw_node::~count_raw_node() {
}

void count_raw_node::add_node(count_raw_node* node) {
    
    logger::Instance()->debug("ADD NODE\n");

    for(gmap<int, series_struct>::iterator ssi = node->series.begin(); ssi != node->series.end(); ++ssi) {
        
        int id = ssi->first;
        
        for (std::map< rpos,rcount >::iterator li = ssi->second.lefts.begin(); li != ssi->second.lefts.end(); ++li) {  

            std::map< rpos,rcount >::iterator hit = series[id].lefts.find(li->first);
            if (hit == series[id].lefts.end()) {
                series[id].lefts.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator ri = ssi->second.rights.begin(); ri != ssi->second.rights.end(); ++ri) {      
            std::map< rpos,rcount >::iterator hit = series[id].rights.find(ri->first);
            if (hit == series[id].rights.end()) {
                series[id].rights.insert(*ri);
            } else {
                hit->second += ri->second;
            }
        }
        
        logger::Instance()->debug("NODE CHECK A " + std::to_string(ssi->second.total_lefts) + " - " + std::to_string(ssi->second.total_rights) + "\n");
        
        series[id].total_lefts += ssi->second.total_lefts;
        series[id].total_rights += ssi->second.total_rights;
        
        logger::Instance()->debug("NODE OUT " + std::to_string(series[id].total_lefts) + " - " + std::to_string(series[id].total_rights) + "\n");
    }
}

void count_raw_node::add_node_start(exon_group* exons) {
    // copy evidence into lefts

    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Node Start at. " + exons->bin_mask.to_string() + "\n");
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
            rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.lefts->begin(); it != ssi->second.lefts->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.rights->begin(); it != ssi->second.rights->end(); ++it) {
            ee += it->second;
        }
        logger::Instance()->debug("Start AT  " + std::to_string(ssi->second.total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(ssi->second.total_lefts) + "\n");
    }
   #endif 
   
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {

        int id = ssi->first;
        
        for (std::map< rpos,rcount >::iterator li = ssi->second.lefts->begin(); li != ssi->second.lefts->end(); ++li) {  

            std::map< rpos,rcount >::iterator hit = series[id].lefts.find(li->first);
            if (hit == series[id].lefts.end()) {
                series[id].lefts.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[0].begin(); li != ssi->second.hole_ends[0].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].lefts.find(li->first);
            if (hit == series[id].lefts.end()) {
                series[id].lefts.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[0].begin(); li != ssi->second.hole_starts[0].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].rights.find(li->first);
            if (hit == series[id].rights.end()) {
                series[id].rights.insert(*li);
            } else {
                hit->second += li->second;
            }
        }

        series[id].total_lefts += ssi->second.total_lefts + ssi->second.hole_end_counts[0] - ssi->second.hole_start_counts[0];

        logger::Instance()->debug("NODE CHECK B " + std::to_string(ssi->second.total_lefts) + " - " + std::to_string(ssi->second.hole_end_counts[0]) + " - "+ std::to_string(ssi->second.hole_start_counts[0]) + "\n");

        if (exons->bin_mask.id.count() == 1) { // just one long, lefts are rights
            for (std::map< rpos,rcount >::iterator li = ssi->second.rights->begin(); li != ssi->second.rights->end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = series[id].rights.find(li->first);
                if (hit == series[id].rights.end()) {
                    series[id].rights.insert(*li);
                } else {
                    hit->second += li->second;
                }
            }
             // there can't be an imbalance for single exon evidence
             logger::Instance()->debug("NODE CHECK B Single " + std::to_string(series[id].total_lefts) + " :- " + std::to_string(ssi->second.total_rights) + "\n");

            series[id].total_lefts -= ssi->second.total_rights;
        }
        
        logger::Instance()->debug("NODE OUT " + std::to_string(series[id].total_lefts) + " - " + std::to_string(series[id].total_rights) + "\n");
    }
    
    #ifdef ALLOW_DEBUG
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.lefts->begin(); it != ssi->second.lefts->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.rights->begin(); it != ssi->second.rights->end(); ++it) {
            ee += it->second;
        }
        logger::Instance()->debug("Set TO  " + std::to_string(ssi->second.total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(ssi->second.total_lefts) + "\n");

    }
   
   #endif
    
}


void count_raw_node::add_node_end(exon_group* exons) {
    
    unsigned int last_i = exons->set_exons -1;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Node END at. " + exons->bin_mask.to_string() + " " + std::to_string(last_i) + "\n");

    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        logger::Instance()->debug("Start AT  " + std::to_string(ssi->second.total_rights) + " " + std::to_string(ssi->second.total_lefts) + "\n");
    }
    
    #endif
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {

        int id = ssi->first;
    
        // copy evidence into lefts
        for (std::map< rpos,rcount >::iterator li = ssi->second.rights->begin(); li != ssi->second.rights->end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].rights.find(li->first);
            if (hit == series[id].rights.end()) {
                series[id].rights.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[last_i].begin(); li != ssi->second.hole_ends[last_i].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].lefts.find(li->first);
            if (hit == series[id].lefts.end()) {
                series[id].lefts.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[last_i].begin(); li != ssi->second.hole_starts[last_i].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].rights.find(li->first);
            if (hit == series[id].rights.end()) {
                series[id].rights.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        
        logger::Instance()->debug("NODE CHECK C " + std::to_string(series[id].total_rights) + " :" + std::to_string(ssi->second.total_rights) + " - " + std::to_string(ssi->second.hole_start_counts[last_i]) + " - "+ std::to_string(ssi->second.hole_end_counts[last_i]) + "\n");

        series[id].total_rights += ssi->second.total_rights + ssi->second.hole_start_counts[last_i] - ssi->second.hole_end_counts[last_i];
        
        logger::Instance()->debug("NODE OUT " + std::to_string(series[id].total_lefts) + " - " + std::to_string(series[id].total_rights) + "\n");
    }
    
    #ifdef ALLOW_DEBUG
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.lefts->begin(); it != ssi->second.lefts->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.rights->begin(); it != ssi->second.rights->end(); ++it) {
            ee += it->second;
        }
        logger::Instance()->debug("Set TO  " + std::to_string(ssi->second.total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(ssi->second.total_lefts) + "\n");
    }
    #endif
 
}


void count_raw_node::add_node_initial_index(exon_group* exons, unsigned int i, gmap<int, rcount> &next_value) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Node Middle at. " + exons->bin_mask.to_string() + " " + std::to_string(i) + "\n");
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        logger::Instance()->debug("Start AT  " + std::to_string(ssi->second.total_rights) + " " + std::to_string(ssi->second.total_lefts) + "\n");
    }
    #endif
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {

        int id = ssi->first;
    
        logger::Instance()->debug("NODE CHECK D " + std::to_string(next_value[id]) + " - " + std::to_string( ssi->second.hole_end_counts[i]) + " - "+ std::to_string(ssi->second.hole_start_counts[i]) + "\n");

        
        series[id].total_rights += next_value[id];
        next_value[id] = next_value[id] + ssi->second.hole_end_counts[i] - ssi->second.hole_start_counts[i];
        series[id].total_lefts += next_value[id];

        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[i].begin(); li != ssi->second.hole_ends[i].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].lefts.find(li->first);
            if (hit == series[id].lefts.end()) {
                series[id].lefts.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[i].begin(); li != ssi->second.hole_starts[i].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].rights.find(li->first);
            if (hit == series[id].rights.end()) {
                series[id].rights.insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        
        logger::Instance()->debug("NODE OUT " + std::to_string(series[id].total_lefts) + " - " + std::to_string(series[id].total_rights) + "\n");
    }
    
    #ifdef ALLOW_DEBUG
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        rcount ss = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.lefts->begin(); it != ssi->second.lefts->end(); ++it) {
            ss += it->second;
        }
        rcount ee = 0;
        for (std::map< rpos,rcount >::iterator it = ssi->second.rights->begin(); it != ssi->second.rights->end(); ++it) {
            ee += it->second;
        }
        logger::Instance()->debug("Set TO  " + std::to_string(ssi->second.total_rights) + " S " + std::to_string(ss) + " E " + std::to_string(ee) + " " + std::to_string(ssi->second.total_lefts) + "\n");
    }
    #endif
    
}

void count_raw_node::add_node_index(count_raw_edge& edge, unsigned int i) {
    
  #ifdef ALLOW_DEBUG  
  logger::Instance()->debug("Split NODE from Existing edge " + std::to_string(i) + "\n");
  #endif
  
   for(gmap<int, count_raw_edge::series_struct>::iterator ssi = edge.series.begin(); ssi != edge.series.end(); ++ssi) {

        int id = ssi->first;
  
       logger::Instance()->debug("NODE CHECK E " + std::to_string(ssi->second.splits[i]) + " - " + std::to_string(ssi->second.splits[i+1]) + "\n");

        
        series[id].total_rights += ssi->second.splits[i];
        series[id].total_lefts += ssi->second.splits[i+1];

        rcount count_l = 0;
        rcount count_r = 0;

        // merge the maps!
        for (std::map< rpos,rcount >::iterator li = ssi->second.starts[i]->begin(); li != ssi->second.starts[i]->end(); ++li) {      
             std::map< rpos,rcount >::iterator hit = series[id].lefts.find(li->first);
             count_l += li->second;
             if (hit == series[id].lefts.end()) {
                 series[id].lefts.insert(*li);
             } else {
                 hit->second += li->second;
             }
         }
         for (std::map< rpos,rcount >::iterator li = ssi->second.ends[i]->begin(); li != ssi->second.ends[i]->end(); ++li) {      
             std::map< rpos,rcount >::iterator hit = series[id].rights.find(li->first);
              count_r += li->second;
             if (hit == series[id].rights.end()) {
                 series[id].rights.insert(*li);
             } else {
                 hit->second += li->second;
             }
         }
        
        logger::Instance()->debug("NODE OUT " + std::to_string(series[id].total_lefts) + " - " + std::to_string(series[id].total_rights) + "\n");

   }
   
   #ifdef ALLOW_DEBUG
    for(gmap<int, count_raw_edge::series_struct>::iterator ssi = edge.series.begin(); ssi != edge.series.end(); ++ssi) {
        int id = ssi->first;
        logger::Instance()->debug("Set TO  " + std::to_string(series[id].total_rights) +  " " + std::to_string(series[id].total_lefts) + "\n");
    }
    #endif
}



std::ostream& operator<<(std::ostream& os, const count_raw_node& crn)
{
    for (gmap<int, count_raw_node::series_struct>::const_iterator iss = crn.series.begin(); iss != crn.series.end(); ++iss) {
        os << std::to_string(iss->first) << ":" << std::to_string(iss->second.total_lefts) << "-" << std::to_string(iss->second.total_rights) << ";";
    }
    return os;
}
