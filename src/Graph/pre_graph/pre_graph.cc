/* 
 * File:   pre_graph.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on November 17, 2015, 9:48 AM
 */

#include <algorithm>
#include <deque>

#include "pre_graph.h"
#include "../../Datatype_Templates/move.h"

pre_graph::pre_graph() : size(0), paired(false), frag_count_chrom(0), read_count_chrom(0), 
    frag_count_region(0), read_count_region(0) {
    
    pairs_for_exon = NULL;
    single_for_exon = NULL;
    
}

pre_graph::~pre_graph() {
    
    if (pairs_for_exon != NULL)
        delete [] pairs_for_exon;
    if (single_for_exon != NULL)
        delete [] single_for_exon;
         
}

void pre_graph::set_size(const unsigned int size)  {
    this->size = size;
    
    pairs_for_exon =  new graph_list<paired_exon_group *>[size]();
    single_for_exon = new graph_list<std::pair<exon_group *, rcount> >[size]();
}

void pre_graph::initialize_exon_gaps_paired(std::deque<overlap_node> &nodes, std::deque<contained_node> &contained_nodes) {

    std::deque<overlap_node>::iterator oli = nodes.begin(); //because of inconsistencies, we use only maximal nodes
    std::deque<contained_node>::iterator coi = contained_nodes.begin();
    
    for (graph_list<paired_exon_group>::iterator pi = paired_bin_list.begin(); pi != paired_bin_list.end(); ++pi) {
        
        while (oli!=nodes.end() && oli->exons->bin_mask.true_smaller(pi->left_read->bin_mask)) {
            ++oli;
        }
        if (oli == nodes.end()) {
            break;
        }
        while (coi!=contained_nodes.end() && coi->exons->bin_mask.true_smaller(pi->left_read->bin_mask)) {
            ++coi; 
        }

        if (oli->exons->bin_mask == pi->left_read->bin_mask) {
           
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Set paired " + pi->left_read->bin_mask.to_string() + " " + pi->right_read->bin_mask.to_string() +"\n");
            #endif 

            // gives evidence for all exons in between pairs
            for (unsigned int i = pi->left_read->range_end + 1; i < pi->right_read->range_start; ++i) {
                pairs_for_exon[i].push_back(&*pi);
            }

            if (pi->left_read->range_end - pi->left_read->range_start >= 1) {
                pairs_for_exon[pi->left_read->range_end].push_back(&*pi);
            }

            if (pi->right_read->range_end - pi->right_read->range_start >= 1) {
                pairs_for_exon[pi->right_read->range_start].push_back(&*pi);
            }
            
        } else if ( size < 2000 && coi!=contained_nodes.end() && coi->contained_in.size() == 1 && coi->exons->bin_mask == pi->left_read->bin_mask){

            pi->left_read = coi->contained_in.front()->exons;
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Insert paired " + pi->left_read->bin_mask.to_string() + " " + pi->right_read->bin_mask.to_string() +"\n");
            #endif 
            
            // unfortunately this needs log n insert, because overlap is restored and earlier in sorting!
            for (unsigned int i = pi->left_read->range_end + 1; i < pi->right_read->range_start; ++i) {    
                graph_list<paired_exon_group *>::iterator ins_point = std::lower_bound(pairs_for_exon[i].begin(), pairs_for_exon[i].end(), &*pi, pointer_comp_paired_exon_group());
                pairs_for_exon[i].insert(ins_point, &*pi);
            } 
            
            if (pi->left_read->range_end - pi->left_read->range_start >= 1) {
                graph_list<paired_exon_group *>::iterator ins_point = std::lower_bound(pairs_for_exon[pi->left_read->range_end].begin(), pairs_for_exon[pi->left_read->range_end].end(), &*pi, pointer_comp_paired_exon_group());
                pairs_for_exon[pi->left_read->range_end].insert(ins_point, &*pi);
            }

            if (pi->right_read->range_end - pi->right_read->range_start >= 1) {
                graph_list<paired_exon_group *>::iterator ins_point = std::lower_bound(pairs_for_exon[pi->right_read->range_start].begin(), pairs_for_exon[pi->right_read->range_start].end(), &*pi, pointer_comp_paired_exon_group());
                pairs_for_exon[pi->right_read->range_start].insert(ins_point, &*pi);
            }
        }   
    }
    
//    // quick test output to see sorting!
//    for (unsigned i=0; i < size; ++i){
//        logger::Instance()->debug("-------------------------- " + std::to_string(i) +"\n");
//        for ( graph_list<paired_exon_group *>::iterator it = pairs_for_exon[i].begin(); it != pairs_for_exon[i].end(); ++it) {
//            logger::Instance()->debug("Insert paired " + (*it)->left_read->bin_mask.to_string() + " " + (*it)->right_read->bin_mask.to_string() +"\n");
//        }        
//    } 
     
}

void pre_graph::initialize_exon_gaps_single( std::deque<overlap_node> &nodes, std::deque<contained_node> &contained_nodes) {
    
    for (std::deque<overlap_node>::iterator oli = nodes.begin() ;  oli != nodes.end(); ++oli) {
        if (oli->exons->range_end - oli->exons->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between

            boost::dynamic_bitset<>::size_type  i =  oli->exons->bin_mask.id.find_next(oli->exons->range_start);
            while (i < oli->exons->range_end) {
                
                rcount additional_count = 0;
                for (graph_list<contained_node*>::iterator coi = oli->contains.begin(); coi !=  oli->contains.end(); ++coi) {
                    if ((*coi)->exons->range_start < i && (*coi)->exons->range_end > i) { // there is a real overlap, index is guaranteed to exist in contains
                        additional_count += (*coi)->exons->frag_count;
                    }
                }
                               
                single_for_exon[i].push_back(std::make_pair(oli->exons, additional_count)); // we do the splitting for matches later when needed
                i = oli->exons->bin_mask.id.find_next(i);
            }
        }
    }
}

void pre_graph::initialize_exon_gaps_single_all( std::deque<overlap_node> &nodes, std::deque<contained_node> &contained_nodes) {
    
    for (std::deque<overlap_node>::iterator oli = nodes.begin() ;  oli != nodes.end(); ++oli) {
        if (oli->exons->range_end - oli->exons->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between

            boost::dynamic_bitset<>::size_type  i =  oli->exons->bin_mask.id.find_next(oli->exons->range_start);
            while (i < oli->exons->range_end) {
                                               
                single_for_exon[i].push_back(std::make_pair(oli->exons, oli->exons->frag_count)); // we do the splitting for matches later when needed
                i = oli->exons->bin_mask.id.find_next(i);
            }
        }
    }
    for (std::deque<contained_node>::iterator oli = contained_nodes.begin() ;  oli != contained_nodes.end(); ++oli) {
        if (oli->exons->range_end - oli->exons->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between

            boost::dynamic_bitset<>::size_type  i =  oli->exons->bin_mask.id.find_next(oli->exons->range_start);
            while (i < oli->exons->range_end) {
            
                single_for_exon[i].push_back(std::make_pair(oli->exons, oli->exons->frag_count)); // we do the splitting for matches later when needed
                i = oli->exons->bin_mask.id.find_next(i);
            }
        }
    }
}


void pre_graph::initialize_exon_gaps_single_raw() {
    
    for(graph_list<exon_group>::iterator it = singled_bin_list.begin(); it != singled_bin_list.end(); ++it) {
        if ( it->reference_atom && it->frag_count == 0 ||  it->extended) {
            continue;
        }
                
        if (it->range_end - it->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between

            boost::dynamic_bitset<>::size_type  i =  it->bin_mask.id.find_next(it->range_start);
            while (i < it->range_end) {
                                               
                single_for_exon[i].push_back(std::make_pair(&*it, it->frag_count)); // we do the splitting for matches later when needed
                i = it->bin_mask.id.find_next(i);
            }
        }
    }
}



void pre_graph::initialize_exon_gaps_single_raw2() {
    
    graph_list< std::pair<exon_group*, unsigned int> > active;
    for(graph_list<exon_group>::iterator it = singled_bin_list.begin(); it != singled_bin_list.end(); ++it) {
        if ( it->reference_atom && it->frag_count == 0) {
            continue;
        }
        
        bool matched = false;
        for (graph_list< std::pair<exon_group*, unsigned int> >::iterator at = active.begin();  at != active.end() && it->range_start <= (at->first)->range_end ; ++at) {
             if ( *at->first > *it ) {
                 matched = true;
                 at->second += it->frag_count;
             }
        }
        
        if (!matched) {
            active.push_back(std::make_pair(&*it, it->frag_count)); 
        }
    }
    
    for (graph_list< std::pair<exon_group*, unsigned int> >::iterator at = active.begin();  at != active.end(); ++at) {
    //for(graph_list<exon_group>::iterator it = singled_bin_list.begin(); it != singled_bin_list.end(); ++it) {
        
        if (at->first->range_end - at->first->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between

            boost::dynamic_bitset<>::size_type  i =  at->first->bin_mask.id.find_next(at->first->range_start);
            while (i < at->first->range_end) {
                                               
                single_for_exon[i].push_back(std::make_pair(at->first, at->second)); // we do the splitting for matches later when needed
                i = at->first->bin_mask.id.find_next(i);
            }
        }
    }
}


void pre_graph::initialize_exon_gaps_paired_raw() {
    
    for (graph_list<paired_exon_group>::iterator pi = paired_bin_list.begin(); pi != paired_bin_list.end(); ++pi) {
       
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Set paired " + pi->left_read->bin_mask.to_string() + " " + pi->right_read->bin_mask.to_string() +"\n");
        #endif 

        // gives evidence for all exons in between pairs
        for (unsigned int i = pi->left_read->range_end + 1; i < pi->right_read->range_start; ++i) {
            pairs_for_exon[i].push_back(&*pi);
        }

        if (pi->left_read->range_end - pi->left_read->range_start >= 1) {
            pairs_for_exon[pi->left_read->range_end].push_back(&*pi);
        }

        if (pi->right_read->range_end - pi->right_read->range_start >= 1) {
            pairs_for_exon[pi->right_read->range_start].push_back(&*pi);
        }
        
//       if (pi->left_read->range_end - pi->left_read->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between
//
//            boost::dynamic_bitset<>::size_type  i =  pi->left_read->bin_mask.id.find_next(pi->left_read->range_start);
//            while (i < pi->left_read->range_end) {
//                                               
//                pairs_for_exon[i].push_back(&*pi); // we do the splitting for matches later when needed
//                i = pi->left_read->bin_mask.id.find_next(i);
//            }
//        }
//        
//        if (pi->right_read->range_end - pi->right_read->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between
//
//            boost::dynamic_bitset<>::size_type  i =  pi->right_read->bin_mask.id.find_next(pi->right_read->range_start);
//            while (i < pi->right_read->range_end) {
//                                               
//                pairs_for_exon[i].push_back(&*pi); // we do the splitting for matches later when needed
//                i = pi->right_read->bin_mask.id.find_next(i);
//            }
//        }      
    }
    
//    // quick test output to see sorting!
//    for (unsigned i=0; i < size; ++i){
//        logger::Instance()->debug("-------------------------- " + std::to_string(i) +"\n");
//        for ( graph_list<paired_exon_group *>::iterator it = pairs_for_exon[i].begin(); it != pairs_for_exon[i].end(); ++it) {
//            logger::Instance()->debug("Insert paired " + (*it)->left_read->bin_mask.to_string() + " " + (*it)->right_read->bin_mask.to_string() +"\n");
//        }        
//    } 
     
}


void pre_graph::initialize_exon_gaps_paired_raw2() {

    graph_list<paired_exon_group* > active;
    for (graph_list<paired_exon_group>::iterator pi = paired_bin_list.begin(); pi != paired_bin_list.end(); ++pi) {
        
        bool matched = false;
        for (graph_list< paired_exon_group* >::iterator at = active.begin();  at != active.end(); ++at) {
            
             if ( ( *(*at)->left_read > *pi->left_read || *(*at)->left_read == *pi->left_read ) 
                     && ( *(*at)->right_read > *pi->right_read || *(*at)->right_read == *pi->right_read ) ) {
                 matched = true;
             }
        }
        
        if (!matched) {
            active.push_back(&*pi); 
        }
    }

    for (graph_list<paired_exon_group* >::iterator pii = active.begin(); pii != active.end(); ++pii) {
       
        paired_exon_group* pi = *pii;
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Set paired " + pi->left_read->bin_mask.to_string() + " " + pi->right_read->bin_mask.to_string() +"\n");
        #endif 

        // gives evidence for all exons in between pairs
        for (unsigned int i = pi->left_read->range_end + 1; i < pi->right_read->range_start; ++i) {
            pairs_for_exon[i].push_back(&*pi);
        }

        if (pi->left_read->range_end - pi->left_read->range_start >= 1) {
            pairs_for_exon[pi->left_read->range_end].push_back(&*pi);
        }

        if (pi->right_read->range_end - pi->right_read->range_start >= 1) {
            pairs_for_exon[pi->right_read->range_start].push_back(&*pi);
        }
            
    }
        
}

//    for (graph_list<exon_group>::iterator si = singled_bin_list.begin(); si != singled_bin_list.end(); ++si) {
//        if (si->range_end - si->range_start + 1 > 2) {/
//            
//        }
//        
//    }

// creates the lists for single read data and paired end info
//void pre_graph::reduce_to_bins() {
//    
//    // should not happen if used correctly, however
//    if (singled_bin_list.empty()) {
//        return;
//    }
//    
//    // sort in N log N to find clusters in afterwards N
//    std::sort(singled_bin_list.begin(), singled_bin_list.end(), exon_group_sorter);
//    
//    // list with doubles removed
//    graph_list<exon_group> new_bins;
//    
//    graph_list<exon_group>::iterator it = singled_bin_list.begin();
//    new_bins.push_back(_MOVE(*it)); 
//    exon_group* base = new_bins.back();
//    update_parent(base, &(*it));
//    ++it;
//    
//    for(; it != singled_bin_list.end(); ++it) {
//        if (*it == *base) {
//            base->read_count += it->read_count;
//            base->frag_count += it->frag_count;
//            update_parent(base, &(*it));
//            
//        } else {
//            new_bins.push_back(_MOVE(*it)); 
//            base = new_bins.back();
//            update_parent(base, &(*it));
//        }
//    }
//}

//// helper function to add reads to assorted lists
//void pre_graph::update_parent(const exon_group* base, const exon_group* merged) {
//    
//    if (merged->parent == NULL) {
//        return;
//    }
//    
//    if (merged->parent->left_read == merged) {
//        merged->parent->left_read = base;
//    } else {
//        merged->parent->right_read = base;
//    }
//    
//}