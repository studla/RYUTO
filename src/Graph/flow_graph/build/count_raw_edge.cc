/* 
 * File:   node_count_raw.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 15, 2015, 9:59 AM
 */

#include <vector>
#include <boost/concept_check.hpp>

#include "count_raw_edge.h"
#include "../../../Logger/logger.h"
#include "../../../Datatype_Templates/move.h"

count_raw_edge::count_raw_edge() : size(0), initialized(false){
}

count_raw_edge::count_raw_edge(unsigned int size) : size(size), initialized(false) {
    
}

count_raw_edge::~count_raw_edge() {
}

void count_raw_edge::add_edge(count_raw_edge* edge) {
    
    for(gmap<int, series_struct>::iterator ssi = edge->series.begin(); ssi != edge->series.end(); ++ssi) {
        int id = ssi->first;
        
        if (!series[id].initialized) {
            series[id].starts.assign(size, {});
            series[id].ends.assign(size, {});
            series[id].splits.assign(size + 1, 0);
            series[id].initialized = true;
        }
        
        series[id].splits[0] += ssi->second.splits[0];
        for (unsigned int i = 0; i < size; ++i) {

            for (std::map< rpos,rcount >::iterator li = ssi->second.starts[i]->begin(); li != ssi->second.starts[i]->end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = series[id].starts[i]->find(li->first);
                if (hit == series[id].starts[i]->end()) {
                    series[id].starts[i]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }
            for (std::map< rpos,rcount >::iterator li = ssi->second.ends[i]->begin(); li !=  ssi->second.ends[i]->end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = series[id].ends[i]->find(li->first);
                if (hit == series[id].ends[i]->end()) {
                    series[id].ends[i]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }

            series[id].splits[i+1] += ssi->second.splits[i+1];
        }   
    }
    
}


void count_raw_edge::add_initial_count(exon_group* exons) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Create Initial. " + exons->bin_mask.to_string() + " " + std::to_string(exons->size) + " " + std::to_string(exons->set_exons) + "\n");
    #endif
    
    initialized = true;
    size = exons->set_exons - 2;
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        int id = ssi->first;
        
        if (!series[id].initialized) {
            series[id].starts.assign(size, {});
            series[id].ends.assign(size, {});
            series[id].splits.assign(size + 1, 0);
            series[id].initialized = true;
        }
    
        series[id].splits[0] = ssi->second.total_lefts + ssi->second.hole_end_counts[0] - ssi->second.hole_start_counts[0];
        if (size != 0) {
            series[id].splits[size] = ssi->second.total_rights;
        }
        for (unsigned int i = 1; i < exons->set_exons - 1; ++i) {
            _MOVE_RANGE(ssi->second.hole_starts[i].begin(), ssi->second.hole_starts[i].end(), std::inserter(series[id].ends[i-1].ref(), series[id].ends[i-1]->begin()));
            _MOVE_RANGE(ssi->second.hole_ends[i].begin(), ssi->second.hole_ends[i].end(), std::inserter(series[id].starts[i-1].ref(), series[id].starts[i-1]->begin()));
            series[id].splits[i] = series[id].splits[i-1] + ssi->second.hole_end_counts[i] - ssi->second.hole_start_counts[i];
        }
    }
}


void count_raw_edge::add_sub_counts_eg(exon_group* exons, int from, int to, gmap<int, rcount> &next_value) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Subcount RAW. " + exons->bin_mask.to_string() + " " + std::to_string(exons->set_exons) + " " + std::to_string(from) + " " + std::to_string(to) + "\n");
    #endif
    
    if (!initialized) { // unitialized arc
        
        initialized = true;
        size = to - from + 1;
    }
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        int id = ssi->first;

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("ID " + std::to_string(id) + "\n");
            #endif
        
        if (!series[id].initialized) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Init Size " + std::to_string(size) + "\n");
            #endif
            
            series[id].starts.assign(size, {});
            series[id].ends.assign(size, {});
            series[id].splits.assign(size + 1, 0);
            series[id].initialized = true;


            series[id].splits[0] += next_value[id];
            unsigned int g = 1;
            for (unsigned int i = from; i <= to; ++i, ++g) {

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("I G " + std::to_string(i) + " " + std::to_string(g) + "\n");
            #endif
                
                std::copy(ssi->second.hole_starts[i].begin(), ssi->second.hole_starts[i].end(), std::inserter(series[id].ends[g-1].ref(), series[id].ends[g-1]->begin()));
                std::copy(ssi->second.hole_ends[i].begin(), ssi->second.hole_ends[i].end(), std::inserter(series[id].starts[g-1].ref(), series[id].starts[g-1]->begin()));
                next_value[id] = next_value[id] + ssi->second.hole_end_counts[i] - ssi->second.hole_start_counts[i];
                series[id].splits[g] += next_value[id];
            } 
        } else {

            series[id].splits[0] += next_value[id];
            unsigned int g = 1;
            for (unsigned int i = from; i <= to; ++i, ++g) {

                for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[i].begin(); li != ssi->second.hole_ends[i].end(); ++li) {      
                    std::map< rpos,rcount >::iterator hit = series[id].starts[g-1]->find(li->first);
                    if (hit == series[id].starts[g-1]->end()) {
                        series[id].starts[g-1]->insert((*li));
                    } else {
                        hit->second += li->second;
                    }
                }
                for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[i].begin(); li != ssi->second.hole_starts[i].end(); ++li) {      
                    std::map< rpos,rcount >::iterator hit = series[id].ends[g-1]->find(li->first);
                    if (hit == series[id].ends[g-1]->end()) {
                        series[id].ends[g-1]->insert((*li));
                    } else {
                        hit->second += li->second;
                    }
                }

                next_value[id] = next_value[id] + ssi->second.hole_end_counts[i] - ssi->second.hole_start_counts[i];
                series[id].splits[g] += next_value[id];
            }  
        }
    }
    
//    #ifdef ALLOW_DEBUG
//    for (unsigned int i = 0; i < size; ++i) {
//       
//        logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
//        
//        rcount ss = 0;
//        for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
//            ss += it->second;
//        }
//        rcount ee = 0;
//        for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
//            ee += it->second;
//        }
//        
//        logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
//    }
//    logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
//    #endif
    
}


void count_raw_edge::add_sub_counts_start(exon_group* exons, int g, gmap<int, rcount> &next_value) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcount start. " + exons->bin_mask.to_string() + " " + std::to_string(g) + "\n");
    #endif
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        int id = ssi->first;

        if (!series[id].initialized) {
            series[id].starts.assign(size, {});
            series[id].ends.assign(size, {});
            series[id].splits.assign(size + 1, 0);
            series[id].initialized = true;
        }
    
        for (std::map< rpos,rcount >::iterator li = ssi->second.lefts->begin(); li != ssi->second.lefts->end(); ++li) { 
            std::map< rpos,rcount >::iterator hit = series[id].starts[g]->find(li->first);
            if (hit == series[id].starts[g]->end()) {
                series[id].starts[g]->insert(*li);
            } else {
                hit->second += li->second;
            }
        }    
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[0].begin(); li != ssi->second.hole_ends[0].end(); ++li) { 

            std::map< rpos,rcount >::iterator hit = series[id].starts[g]->find(li->first);
            if (hit == series[id].starts[g]->end()) {
                series[id].starts[g]->insert(*li);
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[0].begin(); li != ssi->second.hole_starts[0].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].ends[g]->find(li->first);
            if (hit == series[id].ends[g]->end()) {
                series[id].ends[g]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }

        series[id].splits[g+1] += ssi->second.total_lefts + ssi->second.hole_end_counts[0] - ssi->second.hole_start_counts[0];

        if (exons->bin_mask.id.count() == 1) { // just one long, lefts are rights
            for (std::map< rpos,rcount >::iterator li = ssi->second.rights->begin(); li != ssi->second.rights->end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = series[id].ends[g]->find(li->first);
                if (hit == series[id].ends[g]->end()) {
                    series[id].ends[g]->insert(*li);
                } else {
                    hit->second += li->second;
                }
            }
             // there can't be an imbalance for single exon evidence
            series[id].splits[g+1] -= ssi->second.total_rights;
        }

//        #ifdef ALLOW_DEBUG
//        for (unsigned int i = 0; i < size; ++i) {
//
//            logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
//
//            rcount ss = 0;
//            for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
//                ss += it->second;
//            }
//            rcount ee = 0;
//            for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
//                ee += it->second;
//            }
//
//            logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
//        }
//        logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
//        #endif
        
        next_value[id] = ssi->second.total_lefts + ssi->second.hole_end_counts[0] - ssi->second.hole_start_counts[0];
        
    }
        
}


void count_raw_edge::add_sub_counts_end(exon_group* exons, int g, gmap<int, rcount> &next_value) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcount end. " + exons->bin_mask.to_string() + " " + std::to_string(g) + "\n");
    #endif
    
    unsigned int last_i = exons->set_exons -1;
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        int id = ssi->first;

        if (!series[id].initialized) {
            series[id].starts.assign(size, {});
            series[id].ends.assign(size, {});
            series[id].splits.assign(size + 1, 0);
            series[id].initialized = true;
        }
    
        for (std::map< rpos,rcount >::iterator li = ssi->second.rights->begin(); li != ssi->second.rights->end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].ends[g]->find(li->first);
            if (hit == series[id].ends[g]->end()) {
                series[id].ends[g]->insert(*li);
            } else {
                hit->second += li->second;
            }
        }    
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[last_i].begin(); li != ssi->second.hole_ends[last_i].end(); ++li) { 
            std::map< rpos,rcount >::iterator hit = series[id].starts[g]->find(li->first);
            if (hit == series[id].starts[g]->end()) {
                series[id].starts[g]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }
        for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[last_i].begin(); li != ssi->second.hole_starts[last_i].end(); ++li) {      
            std::map< rpos,rcount >::iterator hit = series[id].ends[g]->find(li->first);
            if (hit == series[id].ends[g]->end()) {
                series[id].ends[g]->insert((*li));
            } else {
                hit->second += li->second;
            }
        }
    
        //splits[g] += exons->total_rights + exons->hole_start_counts[last_i] - exons->hole_end_counts[last_i]; by range!

//        #ifdef ALLOW_DEBUG
//        for (unsigned int i = 0; i < size; ++i) {
//
//           logger::Instance()->debug("Split " + std::to_string(splits[i]) + "\n"); 
//
//           rcount ss = 0;
//           for (std::map< rpos,rcount >::iterator it = starts[i]->begin(); it != starts[i]->end(); ++it) {
//               ss += it->second;
//           }
//           rcount ee = 0;
//           for (std::map< rpos,rcount >::iterator it = ends[i]->begin(); it != ends[i]->end(); ++it) {
//               ee += it->second;
//           }
//
//           logger::Instance()->debug("S " + std::to_string(ss) + " EE " + std::to_string(ee) + "\n"); 
//        }
//        logger::Instance()->debug("Split " + std::to_string(splits[size]) + "\n"); 
//        #endif
        
         next_value[id] = series[id].splits[g+1];
    }

}

void count_raw_edge::add_sub_counts_eg_range(exon_group* exons, int from, int to, unsigned int g1, gmap<int, rcount> &next_value) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcount raw range. " + exons->bin_mask.to_string() + " " + std::to_string(from) + " " + std::to_string(to) + " " + std::to_string(g1) + "\n");
    #endif
    
    if (to < from) {
        return ;
    }
    
    for(gmap<int, exon_group_count>::iterator ssi = exons->count_series.begin(); ssi != exons->count_series.end(); ++ssi) {
        int id = ssi->first;

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Splits Before " + std::to_string(id) + ": ");
        for(std::vector<rcount>::const_iterator is = series[id].splits.begin(); is != series[id].splits.end(); ++is) {
            logger::Instance()->debug(std::to_string(*is)+ ", ");
        }
        logger::Instance()->debug("\n");
        #endif
        
        if (!series[id].initialized) {
            series[id].starts.assign(size, {});
            series[id].ends.assign(size, {});
            series[id].splits.assign(size + 1, 0);
            series[id].initialized = true;
        }
    
        //splits[g1-1] += split_v; we explicitly don't do that here! set externally!
        unsigned int g = g1;
        for (unsigned int i = from; i <= to; ++i, ++g) {

            for (std::map< rpos,rcount >::iterator li = ssi->second.hole_ends[i].begin(); li != ssi->second.hole_ends[i].end(); ++li) {  
                std::map< rpos,rcount >::iterator hit = series[id].starts[g-1]->find(li->first);
                if (hit == series[id].starts[g-1]->end()) {
                    series[id].starts[g-1]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }

            for (std::map< rpos,rcount >::iterator li = ssi->second.hole_starts[i].begin(); li != ssi->second.hole_starts[i].end(); ++li) {      
                std::map< rpos,rcount >::iterator hit = series[id].ends[g-1]->find(li->first);
                if (hit == series[id].ends[g-1]->end()) {
                    series[id].ends[g-1]->insert((*li));
                } else {
                    hit->second += li->second;
                }
            }

            next_value[id] = next_value[id] + ssi->second.hole_end_counts[i] - ssi->second.hole_start_counts[i];
            series[id].splits[g] += next_value[id];
        } 


        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Splits After " + std::to_string(id) + ": ");
        for(std::vector<rcount>::const_iterator is = series[id].splits.begin(); is != series[id].splits.end(); ++is) {
            logger::Instance()->debug(std::to_string(*is)+ ", ");
        }
        logger::Instance()->debug("\n");
        #endif
    }
}

void count_raw_edge::add_sub_counts_re(count_raw_edge &a, int from, int to) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add subcounts existing. " + std::to_string(a.size) + " " + std::to_string(from) + " " + std::to_string(to) + "\n");
    #endif
    
    if (from > to) {
        
        if (!initialized) {
            initialized = true;
            //just one split site remains, after that un-subable
        }   
        for(gmap<int, series_struct>::iterator ssi = a.series.begin(); ssi != a.series.end(); ++ssi) {
            int id = ssi->first;

            if (!series[id].initialized) {

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Splits Before " + std::to_string(id) + ": -\n");
                #endif

                series[id].splits.push_back(a.series[id].splits[from]);
                series[id].initialized = true;
            } else { 

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Splits Before " + std::to_string(id) + ": " + std::to_string(series[id].splits[0]) + "\n");
                #endif  

                series[id].splits[0] += a.series[id].splits[from];
            }

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Splits After " + std::to_string(id) + ": " + std::to_string(series[id].splits[0]) + "\n");
            #endif  
        }
        return;
    }
    
    if (!initialized) { // copy in
        
        initialized = true;
        size = to - from + 1;
    }
    for(gmap<int, series_struct>::iterator ssi = a.series.begin(); ssi != a.series.end(); ++ssi) {
        int id = ssi->first;

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Splits Before " + std::to_string(id) + ": ");
        for(std::vector<rcount>::const_iterator is = series[id].splits.begin(); is != series[id].splits.end(); ++is) {
            logger::Instance()->debug(std::to_string(*is)+ ", ");
        }
        logger::Instance()->debug("\n");
        #endif
        
        if (!series[id].initialized) {
            std::copy(ssi->second.starts.begin() + from, ssi->second.starts.begin() + to + 1, std::back_inserter(series[id].starts));
            std::copy(ssi->second.ends.begin() + from, ssi->second.ends.begin() + to + 1, std::back_inserter(series[id].ends));
            std::copy(ssi->second.splits.begin() + from, ssi->second.splits.begin() + to + 2, std::back_inserter(series[id].splits));
            series[id].initialized = true;
        } else {

            unsigned g = 0;
            unsigned int i = from;
            for (; i <= to; ++i, ++g) {
                series[id].splits[g] += ssi->second.splits[i]; 

                for (std::map< rpos,rcount >::iterator li = ssi->second.starts[i]->begin(); li != ssi->second.starts[i]->end(); ++li) { 
                    std::map< rpos,rcount >::iterator hit = series[id].starts[g]->find(li->first);
                    if (hit == series[id].starts[g]->end()) {
                        series[id].starts[g]->insert((*li));
                    } else {
                        hit->second += li->second;
                    }
                }
                for (std::map< rpos,rcount >::iterator li = ssi->second.ends[i]->begin(); li != ssi->second.ends[i]->end(); ++li) {      
                    std::map< rpos,rcount >::iterator hit = series[id].ends[g]->find(li->first);
                    if (hit == series[id].ends[g]->end()) {
                        series[id].ends[g]->insert((*li));
                    } else {
                        hit->second += li->second;
                    }
                }
            }
            series[id].splits[g] += ssi->second.splits[i];
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Splits After " + std::to_string(id) + ": ");
        for(std::vector<rcount>::const_iterator is = series[id].splits.begin(); is != series[id].splits.end(); ++is) {
            logger::Instance()->debug(std::to_string(*is)+ ", ");
        }
        logger::Instance()->debug("\n");
        #endif
        
    }
    
}


std::ostream& operator<<(std::ostream& os, const count_raw_edge& cre)
{
    for (gmap<int, count_raw_edge::series_struct>::const_iterator iss = cre.series.begin(); iss != cre.series.end(); ++iss) {
        os << std::to_string(iss->first) << ":"; 
        if (iss->second.initialized) {    
            for(std::vector<rcount>::const_iterator is = iss->second.splits.begin(); is != iss->second.splits.end(); ++is) {
                os << std::to_string(*is) << ",";
            }
        } else { 
            os << "-;";
        }
    }
    return os;
}