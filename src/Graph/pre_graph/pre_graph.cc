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

pre_graph::pre_graph() : size(0), paired(false), average_fragment_length(1) {
    
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
    single_for_exon = new graph_list<std::pair<exon_group *, gmap<int, rcount> > >[size]();
}

void pre_graph::initialize_exon_gaps_single_raw() {
    
    for(graph_list<exon_group>::iterator it = singled_bin_list.begin(); it != singled_bin_list.end(); ++it) {
        if ( it->reference_atom && !it->has_coverage || it->extended) {
            continue;
        }
        
        gmap<int, rcount> sfe_map;
        for (gmap<int, exon_group_count>::iterator egci = it->count_series.begin(); egci != it->count_series.end(); ++egci) {
            sfe_map.insert(std::make_pair( egci->first, egci->second.frag_count));
        }
        
        if (it->range_end - it->range_start + 1 > 2) { // we could have > 3 exons, so we have overarching evidence :) for all exons in between

            boost::dynamic_bitset<>::size_type  i =  it->bin_mask.id.find_next(it->range_start);
            while (i < it->range_end) {
                         
                single_for_exon[i].push_back(std::make_pair(&*it, sfe_map)); // we do the splitting for matches later when needed
                i = it->bin_mask.id.find_next(i);
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
