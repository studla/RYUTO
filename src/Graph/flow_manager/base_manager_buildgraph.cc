/* 
 * File:   base_manager.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on December 1, 2015, 3:17 PM
 */

#include <deque>
#include <set>
#include <boost/unordered/unordered_map.hpp>

#include "base_manager.h"
#include "../pre_graph/pre_graph.h"
#include "../../Datatype_Templates/graph_list.h"
#include "../overlap_graph/overlap_node.h"
#include "../overlap_graph/contained_node.h"
#include "../flow_graph/exon_edge.h"
#include "../../Datatype_Templates/move.h"
#include "../overlap_graph/range_helper.h"
#include "../../Logger/logger.h"
#include "../../Options/options.h"
#include "../flow_graph/coverage/flow_series.h"

#include <deque>

typedef lemon::ListGraph Graph;


base_manager::base_manager( pre_graph* raw,  exon_meta* meta,  const std::string &chromosome, std::set<int> &ids) : meta(meta), raw(raw), chromosome(chromosome), ai(g), node_index(g), fs(g), input_ids(ids), input_count(ids.size())  {
    
}

base_manager::~base_manager() {
}

bool base_manager::build_basic_splitgraph() {  //TODO: REWORK COVERAGE
    
    logger::Instance()->info("Build Splitgraph, Bin List: "+std::to_string(raw->singled_bin_list.size()) +" Region: " + std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[meta->size-1].right) + ".\n");
    
    if (raw->singled_bin_list.empty()) {
        logger::Instance()->warning("All Reads filtered on "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[meta->size-1].right) +"\n");
        return false;
    }
    
    if (meta->size == 1) {
        logger::Instance()->warning("Single exon on "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[0].right) +"\n");
        
        std::map<int, float> capacity;
        for ( std::set<int>::iterator iii = input_ids.begin(); iii !=  input_ids.end(); ++iii) {
            
            if (raw->singled_bin_list.front().count_series.find(*iii) == raw->singled_bin_list.front().count_series.end()) {
                continue;
            }
            
            for (std::map< rpos,rcount >::iterator i = raw->singled_bin_list.front().count_series[*iii].hole_ends[0].begin(); i != raw->singled_bin_list.front().count_series[*iii].hole_ends[0].end(); ++i) {
                raw->singled_bin_list.front().count_series[*iii].lefts.ref()[i->first] += i->second;
            }

            for (std::map< rpos,rcount >::iterator i = raw->singled_bin_list.front().count_series[*iii].hole_starts[0].begin(); i !=  raw->singled_bin_list.front().count_series[*iii].hole_starts[0].end(); ++i) {
                raw->singled_bin_list.front().count_series[*iii].rights.ref()[i->first]+= i->second;
            }

            region r;
            create_region(0, 0, raw->singled_bin_list.front().count_series[*iii].lefts.ref(), raw->singled_bin_list.front().count_series[*iii].rights.ref(), r);
            if (r.get_average() > 3.0) capacity[*iii] = r.get_average();
        }
        
        if (!raw->singled_bin_list.empty() 
                && (options::Instance()->vote(capacity.size(), input_ids.size(), options::delete_on::group_vote_high) ) ) {
            single_exons.insert(single_exon(0, capacity, raw->singled_bin_list.front().reference_atom));
        }
        return false;
    }
    
    const unsigned int size = raw->size;
    
    ListDigraph::ArcMap<count_raw_edge> edge_counts(g);
    ListDigraph::NodeMap<count_raw_node> node_counts(g);
    bool* sourcelinks =  new bool[size]();
    bool* drainlinks =  new bool[size]();
    bool* gnode_set = new bool[raw->size]();
    
    // ##### init graph nodes
    
    s = g.addNode();
    t = g.addNode();
    
    // room to temporarily storing all nodes for the graph
    ListDigraph::Node* gnodes = new ListDigraph::Node[size];
    for (unsigned int i = 0; i < size; ++i) {
        gnodes[i] = g.addNode();
        gnode_set[i] = true;
    }

    // we can make the basic graph all in just one loop
    for(graph_list<exon_group>::iterator it = raw->singled_bin_list.begin();it != raw->singled_bin_list.end(); ++it) {
        
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Add " + it->bin_mask.to_string() + "------------------------------------\n");
            #endif
        
        if (it->source_evidence) {
            sourcelinks[it->range_start] = true;
        }
        if (it->drain_evidence) {
            drainlinks[it->range_end] = true;
        }
            
        count_raw_node* rc = &node_counts[gnodes[it->range_start]];
        rc->add_node_start(&*it);
        
        // special case for single exon evidence!
        if (it->range_start == it->range_end) {
            // just fix up node
            continue;
        }

        unsigned int offset = 0;
        gmap<int, rcount> next_value;
        for (gmap<int, exon_group_count>::iterator egi = it->count_series.begin(); egi != it->count_series.end(); ++egi) {
            next_value[egi->first] = egi->second.total_lefts + egi->second.hole_end_counts[0] - egi->second.hole_start_counts[0];
        }
    
        unsigned int start = it->range_start;
        for (unsigned int i = it->range_start+1; i<= it->range_end; ++i) {
            
            if ( !(*it)[i] ) {
                continue;
            }
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("i " + std::to_string(i) + "\n");
            #endif
            
            // so this index is an exon
            ListDigraph::OutArcIt e(g, gnodes[start]);
            while(e!=INVALID && g.target(e)!=gnodes[i]) ++e;
            
            
            ListDigraph::Arc arc;
            if (e==INVALID) {
                // this arc does not exist so far, we should add it
                arc = g.addArc(gnodes[start], gnodes[i]);
                
                exon_edge* edge = &ai[arc].edge_specifier;
                
                edge->reserve(raw->size);
                ai[arc].edge_specifier.add_exon(start);
                ai[arc].edge_specifier.add_exon(i);
            } else {
                arc = e;
            }
            
            count_raw_edge* rec = &edge_counts[arc];
            rec->add_sub_counts_eg(&*it, offset+1, offset, next_value);
            offset += 1;
              
            if (i < it->range_end) {
                count_raw_node* rcl = &node_counts[gnodes[i]];
                rcl->add_node_initial_index(&*it, offset, next_value);
            } else {
                count_raw_node* rc = &node_counts[gnodes[i]];
                rc->add_node_end(&*it);
            }
            
            start = i;
        }

    }
    
    ListDigraph::ArcIt a(g);
    if(a == INVALID) {
        return false;
    }
    
    // ##### basic graph done already, add sources and sinks
    
    for(unsigned int i=0; i < raw->size;++i) {
        // connect to source
        ListDigraph::OutArcIt e(g, gnodes[i]);
        if (e == INVALID) {
            drainlinks[i] = true;
        }
        
        ListDigraph::InArcIt r(g, gnodes[i]);
        if (r == INVALID) {
            sourcelinks[i] = true;
        }
    }
    
    for(unsigned int i=0; i < raw->size;++i) {
        
        if (sourcelinks[i] && drainlinks[i]) {
            g.erase(gnodes[i]);
            gnode_set[i] = false;
            continue;
        }
        
        if (sourcelinks[i]) {
            // add source to artifical source
            ListDigraph::Arc source_arc = g.addArc(s, gnodes[i]);
            initialize_source_drain_arc(source_arc);
        }

        // test if end also qualifies as a drain
        if (drainlinks[i]) {
            // add source to artifical source
            ListDigraph::Arc drain_arc = g.addArc(gnodes[i], t);
            initialize_source_drain_arc(drain_arc);
        }
    }
        
    // ##### graph done already, now create capacities
    
    create_final_capacities(
    node_counts,
    edge_counts,
    gnodes,
    gnode_set);
    
    
//    graph_list<exon_group>::iterator it = raw->singled_bin_list.begin();
//    std::deque<overlap_node> nodes;
//    std::deque<overlap_node*> active_nodes;
//    std::deque<contained_node> contained_nodes;
//    std::deque<contained_node> filtered_contained_nodes;
//    
//    // init values
//    unsigned int exon_index = 0;
//    
//    bool found = false;
//    while (!found) {
//        if ( it->reference_atom) {
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Added guide " + it->bin_mask.to_string() + ".\n");
//            #endif
//            
//            raw->guide_transcripts.push_back( &*it); 
//            if(it->frag_count == 0) {
//                ++it;
//                if (it == raw->singled_bin_list.end()) {
//                    logger::Instance()->warning("No split reads mapping on/only guides for "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[0].right) +" - Skip Graph\n");
//                    return false;
//                }
//                continue; // don't use this for building graph if no other frags exist
//            }
//        }
//        if (it->range_start == it->range_end) { // we don't want singles!
////            
////            if (it->frag_count >= options::Instance()->get_min_single_coverage() ) {
////                single_exons.insert(std::make_pair(it->range_start, it->frag_count));
////            }
////            
//            ++it;
//            if (it == raw->singled_bin_list.end()) {
//                    logger::Instance()->warning("No split reads mapping on/only guides for "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[0].right) +" - Skip Graph\n");
//                    return false;
//            }
//            continue;
//        }
//        found = true;
//    }
//    
//    nodes.push_back( overlap_node(&(*it)) );
//    active_nodes.push_back( &nodes.back() ); // add first node as active
//    ++it;
//    
//    for(;it != raw->singled_bin_list.end(); ++it) { // 
//        
//        if ( it->reference_atom) {
//            
//            #ifdef ALLOW_DEBUG
//            logger::Instance()->debug("Added guide " + it->bin_mask.to_string() + ".\n");
//            #endif
//            
//            raw->guide_transcripts.push_back( &*it); 
//            if(it->frag_count == 0) {
//                continue; // don't use this for building graph if no other frags exist
//            }
//        }
//        
//        // increase exon_index if needed and update active list
//        if (it->range_start > exon_index) {
//            exon_index = it->range_start;
//            update_active(active_nodes, exon_index); // remove nodes no longer in range
//        }
//        
//        // search all active for containment
//        bool contained = false;
//        // moving the if out may look shitty but improves branch prediction and optimization
//        if (it->length_filterd) {
//            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
//                if ( *(*ac)->exons > *it ) { // active node contains current next node
//                    if (!contained) {
//                        filtered_contained_nodes.push_back(contained_node( &(*it) )); // add to list
//                        contained = true;
//                    }
//
//                    filtered_contained_nodes.back().add_contained( *ac );
//                }
//            }
//        } else {
//            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
//                if ( *(*ac)->exons > *it ) { // active node contains current next node
//                    if (!contained) {
//                        contained_nodes.push_back(contained_node( &(*it) )); // add to list
//                        contained = true;
//                    }
//                    contained_nodes.back().add_contained( *ac );
//                    (*ac)->contains.push_back(&contained_nodes.back());
//                }
//            }
//        }
//        if (!contained) {
// 
//            nodes.push_back( overlap_node(&(*it)) );
//            
//            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
//                if ( *(*ac)->exons >= *it ) {
//                    (*ac)->push_link(&nodes.back());
//                    nodes.back().push_backlink(*ac);
//                }
//            }
//
//            // add to active nodes
//            active_nodes.push_back(&nodes.back());
//        } 
//    }
//    
//    reduce_transitive(nodes);
//    reduce_transitive_contains(contained_nodes);
//    reduce_transitive_contains(filtered_contained_nodes);
//    

    raw->initialize_exon_gaps_single_raw();
    raw->initialize_exon_gaps_paired_raw();
    
    return true;
}

// This function takes the raw reads and builds up a graph for
// inheriting flow modules. Multi-Exon evidence is left intact whenever possible
// without ambigous flow!
bool base_manager::build_extended_splitgraph() {
    
    logger::Instance()->info("Build Splitgraph, Bin List: "+std::to_string(raw->singled_bin_list.size()) +" Region: " + std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[meta->size-1].right) + ".\n");
    
    if (raw->singled_bin_list.empty()) {
        logger::Instance()->warning("All Reads filtered on "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[meta->size-1].right) +"\n");
        return false;
    }
    
    if (meta->size == 1) { // TODO: REAL handling
        logger::Instance()->warning("Single exon on "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[0].right) +"\n");
        
        std::map<int, float> capacity;
        for ( std::set<int>::iterator iii = input_ids.begin(); iii !=  input_ids.end(); ++iii) {
            
            if (raw->singled_bin_list.front().count_series.find(*iii) == raw->singled_bin_list.front().count_series.end()) {
                continue;
            }
            
            for (std::map< rpos,rcount >::iterator i = raw->singled_bin_list.front().count_series[*iii].hole_ends[0].begin(); i != raw->singled_bin_list.front().count_series[*iii].hole_ends[0].end(); ++i) {
                raw->singled_bin_list.front().count_series[*iii].lefts.ref()[i->first] += i->second;
            }

            for (std::map< rpos,rcount >::iterator i = raw->singled_bin_list.front().count_series[*iii].hole_starts[0].begin(); i !=  raw->singled_bin_list.front().count_series[*iii].hole_starts[0].end(); ++i) {
                raw->singled_bin_list.front().count_series[*iii].rights.ref()[i->first]+= i->second;
            }

            region r;
            create_region(0, 0, raw->singled_bin_list.front().count_series[*iii].lefts.ref(), raw->singled_bin_list.front().count_series[*iii].rights.ref(), r);
            if ( (raw->singled_bin_list.front().reference_atom && r.get_max() > 0) || r.get_average() > 3.0) capacity[*iii] = r.get_average();
        }
        
        if (!raw->singled_bin_list.empty() 
                && options::Instance()->vote(capacity.size(), input_ids.size(), options::delete_on::group_vote_high) ) {
            single_exons.insert(single_exon(0, capacity, raw->singled_bin_list.front().reference_atom));
        }
        return false;
    }
     
    // basic datastructures
    graph_list<exon_group>::iterator it = raw->singled_bin_list.begin();
    std::deque<overlap_node> nodes;
    std::deque<overlap_node*> active_nodes;
    std::deque<contained_node> contained_nodes;
    std::deque<contained_node> filtered_contained_nodes;
    
    // circular references to simplify insertion 
    gmap<exon_group*, overlap_node*> circle_frag_to_overlap;
    gmap<exon_group*, contained_node*> circle_frag_to_contained;
    
    // init values
    unsigned int exon_index = 0;
    
    bool found = false;
    while (!found) {
        if ( it->reference_atom) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Added guide " + it->bin_mask.to_string() + ".\n");
            #endif
            
            raw->guide_transcripts.push_back( &*it); 
            if(!it->has_coverage) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("With no Coverage.\n");
                #endif
                ++it;
                if (it == raw->singled_bin_list.end()) {
                    logger::Instance()->warning("No split reads mapping on/only guides for "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[0].right) +" - Skip Graph\n");
                    return false;
                }
                continue; // don't use this for building graph if no other frags exist
            }
        }
        if (it->range_start == it->range_end) { // we don't want singles!
//            
//            if (it->frag_count >= options::Instance()->get_min_single_coverage() ) {
//                single_exons.insert(std::make_pair(it->range_start, it->frag_count));
//            }
//            
            ++it;
            if (it == raw->singled_bin_list.end()) {
                    logger::Instance()->warning("No split reads mapping on/only guides for "+ chromosome + " " +std::to_string(meta->exons[0].left)+" - " + std::to_string(meta->exons[0].right) +" - Skip Graph\n");
                    return false;
            }
            continue;
        }
        found = true;
    }
    
    nodes.push_back( overlap_node(&(*it)) );
    active_nodes.push_back( &nodes.back() ); // add first node as active
    ++it;
    
    for(;it != raw->singled_bin_list.end(); ++it) { // 
        
        if ( it->reference_atom) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Added guide " + it->bin_mask.to_string() + ".\n");
            #endif
            
            raw->guide_transcripts.push_back( &*it); 
            if(!it->has_coverage) {
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("With no Coverage.\n");
                #endif
                continue; // don't use this for building graph if no other frags exist
            }
        }
        
        // increase exon_index if needed and update active list
        if (it->range_start > exon_index) {
            exon_index = it->range_start;
            update_active(active_nodes, exon_index); // remove nodes no longer in range
        }
        
        // search all active for containment
        bool contained = false;
        // moving the if out may look shitty but improves branch prediction and optimization
        if (it->length_filterd) {
            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
                if ( *(*ac)->exons > *it ) { // active node contains current next node
                    if (!contained) {
                        filtered_contained_nodes.push_back(contained_node( &(*it) )); // add to list
                        contained = true;
                    }

                    filtered_contained_nodes.back().add_contained( *ac );
                }
            }
             if (!contained) {
                 // if this is not contained, this was falsely flagged
                 // only elements that are too short should be removed
                 // if a longer element existed at the same place, this HAS to 
                 // fall into it
                 it->length_filterd = false; // remove flag for further reference!
             }
        } else {
            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
                if ( *(*ac)->exons > *it ) { // active node contains current next node
                    if (!contained) {
                        contained_nodes.push_back(contained_node( &(*it) )); // add to list
                        contained = true;
                        
                        if ( it->backlink_fragment ) circle_frag_to_contained[&(*it)] = &contained_nodes.back();
                    }

                    contained_nodes.back().add_contained( *ac );
                    (*ac)->contains.push_back(&contained_nodes.back());
                }
            }
        }
        if (!contained) {
            
//            // we filter out to lowly evidenced fragments
//            // since count it low we ignore them fully!
//            if (it->frag_count <= 2) {
//                continue;
//            }
            
            // if this is a not contained element search for normal overlaps

            if (it->range_start == it->range_end) {
                // this is a standalone exon with no splits
                // single exons cannot form edges, as they have no split
                // and thus constitute a direct link from source to drain,
                // causing a variety of problems: don't add

                continue;
                
            }
            
            nodes.push_back( overlap_node(&(*it)) );
            if ( it->backlink_fragment )  circle_frag_to_overlap[&(*it)] = &nodes.back();
            
            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
                if ( *(*ac)->exons >= *it ) {
                    (*ac)->push_link(&nodes.back());
                    nodes.back().push_backlink(*ac);
                }
            }

            // add to active nodes
            active_nodes.push_back(&nodes.back());
        } else {
            
//            if (it->range_start == it->range_end) {
//                continue;
//            }
//            
//            for (std::deque<overlap_node*>::iterator ac = active_nodes.begin(); ac != active_nodes.end(); ++ac ) { 
//                if ( *(*ac)->exons >= *it ) {
//                    contained_nodes.back().add_contained( *ac );
//                    (*ac)->contains.push_back(&contained_nodes.back());
//                }
//            }
            
        } 
    }
    
    //filter_in_graph(nodes, contained_nodes);
    
    #ifdef ALLOW_DEBUG
    print_raw(nodes, contained_nodes, filtered_contained_nodes);
    #endif
    
    //reduce transitive edges
    reduce_transitive(nodes);
    reduce_transitive_contains(contained_nodes);
    reduce_transitive_contains(filtered_contained_nodes);
    
    #ifdef ALLOW_DEBUG
    print_raw(nodes, contained_nodes, filtered_contained_nodes);
    #endif
    
    // this is for later processing to disambiguate nodes
    raw->initialize_exon_gaps_single_raw();
    raw->initialize_exon_gaps_paired_raw();
    
    // this invalidates contains lists
    reduce_single_nodes(nodes, contained_nodes, raw->singled_bin_list, circle_frag_to_overlap, circle_frag_to_contained);
    
    #ifdef ALLOW_DEBUG
    print_raw(nodes, contained_nodes, filtered_contained_nodes);
    #endif
    
//  logger::Instance()->info("Transitive And Single Done "+ chromosome + " " + std::to_string(nodes.size()) + " " + std::to_string(contained_nodes.size()) + " " + std::to_string(filtered_contained_nodes.size()) +"\n");
    
    create_final_graph(nodes, contained_nodes, filtered_contained_nodes, circle_frag_to_overlap, circle_frag_to_contained);
    
    return true;
}

void base_manager::update_active( std::deque<overlap_node*> &active_nodes, const unsigned int exon_index) {
    
    std::deque<overlap_node*>::iterator it = active_nodes.begin();
    while(it != active_nodes.end()) {
        if((*it)->exons->range_end < exon_index) {
            it = active_nodes.erase(it);
        } else {
            ++it;
        }
    }
   
}


void base_manager::filter_in_graph( std::deque<overlap_node> &nodes, std::deque<contained_node> &contained_nodes) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("In Filter\n");
    #endif

    for(std::deque<overlap_node>::iterator node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
        
        if (node_it->back_links.empty() && node_it->exons->range_start != 0 
                || node_it->links.empty() && node_it->exons->range_end != node_it->exons->size-1 ) {
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Filter Out " + node_it->exons->bin_mask.to_string() + "\n");
            #endif
            
            follow_back_erase(&*node_it); // also does contain removes!

            for (graph_list<overlap_node*>::iterator it = node_it->links.begin(); it != node_it->links.end(); ++it) {
                (*it)->erase_backlink(&*node_it);
            }
            node_it->links.clear();
        }
        
    }
    
    for(std::deque<overlap_node>::iterator node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
        if (node_it->back_links.empty() && node_it->exons->range_start != 0 
                || node_it->links.empty() && node_it->exons->range_end != node_it->exons->size-1 ) {
            node_it->activated = false;
        }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Overlaps Removed\n"); 
    #endif
}

void base_manager::follow_back_erase(overlap_node* on) {
    
    for (graph_list<contained_node*>::iterator it = on->contains.begin(); it!= on->contains.end(); ++it) {
        (*it)->remove_contained(on);
    }
    on->contains.clear();
    
    for (graph_list<overlap_node*>::iterator it = on->back_links.begin(); it != on->back_links.end(); ++it) {
        (*it)->erase_link(on);
        if( (*it)->links.empty() ) {
            follow_back_erase(*it);
        }
    }
    on->back_links.clear();
}

void base_manager::reduce_transitive( std::deque<overlap_node> &nodes) {
    // best option seems to be: E.W. Myers. The fragment assembly string graph
    // if applied while building the graph fragments are lost
    // OK for normal assembly, but probably deadly here
    
    for(std::deque<overlap_node>::iterator node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
        
        // mark adjoining nodes
        for(typename graph_list<overlap_node*>::iterator link_it = node_it->links.begin(); link_it != node_it->links.end(); ++link_it) {
            (*link_it)->marker = overlap_node::INPLAY;
        }
        
        // mark nodes to which edges go that can be eliminated
        for(typename graph_list<overlap_node*>::iterator link_it = node_it->links.begin(); link_it != node_it->links.end(); ++link_it) {
            if ( (*link_it)->marker == overlap_node::INPLAY ) {
                for(typename graph_list<overlap_node*>::iterator link2_it = (*link_it)->links.begin(); link2_it != (*link_it)->links.end(); ++link2_it) {
                    if ( (*link2_it)->marker == overlap_node::INPLAY ) {
                        (*link2_it)->marker = overlap_node::ELIMINATED;
                    }
                }
            }
        }
        
        // now finally mark eliminated edges
        typename graph_list<bool>::iterator edge_mark_it = node_it->edge_marker.begin();
        typename graph_list<overlap_node*>::iterator link_it = node_it->links.begin();
        for(; link_it != node_it->links.end(); ++link_it, ++edge_mark_it) {
            if ( (*link_it)->marker == overlap_node::ELIMINATED ) {
                *edge_mark_it = true; // this edge is to be remove!
            }
            (*link_it)->marker = overlap_node::VACANT;
        }
    }
 
    
// now remove all marked edges!
    for(std::deque<overlap_node>::iterator node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
        
        typename graph_list<bool>::iterator edge_mark_it = node_it->edge_marker.begin();
        typename graph_list<overlap_node*>::iterator link_it = node_it->links.begin();
        
        for(; link_it != node_it->links.end(); ) {
            if ( *edge_mark_it == true ) { // remove both to keep consistent
                
                (*link_it)->erase_backlink(&*node_it);
                
                link_it = node_it->links.erase(link_it);
                edge_mark_it = node_it->edge_marker.erase(edge_mark_it);
            } else {
                ++link_it;
                ++edge_mark_it;
            }
        }
    }
}


void base_manager::reduce_transitive_contains( std::deque<contained_node> &contained) {
    
    // this is a modification of E.W. Myers. to only look at neighbouring nodes
    // of contained ones and mark down containment counters
    
    // loop over nodes marked as contained
    for( std::deque<contained_node>::iterator c_int = contained.begin(); c_int !=  contained.end(); ++c_int) {
        
        // if it is only contained in one, thing to noreduce!
        if ( c_int->contained_in.size() == 1 ) {
            continue;
        }
        
        // mark all nodes neighbouring the contained node
        for(typename graph_list<overlap_node*>::iterator node_it = c_int->contained_in.begin(); node_it != c_int->contained_in.end(); ++node_it) {
            (*node_it)->marker = overlap_node::INPLAY;
        }
        
        // loop through all neighbouring nodes again and test if they reference 
        // an already marked node
        for(typename graph_list<overlap_node*>::iterator node_it = c_int->contained_in.begin(); node_it != c_int->contained_in.end(); ++node_it) {
            for(typename graph_list<overlap_node*>::iterator edge_it = (*node_it)->links.begin(); edge_it != (*node_it)->links.end(); ++edge_it) {
                if ( (*edge_it)->marker == overlap_node::INPLAY ) {
                    (*edge_it)->marker = overlap_node::ELIMINATED;
                }
            }
        }
        
        // now decrease active counter for all edges connecting to eliminated nodes
        for(typename graph_list<overlap_node*>::iterator node_it = c_int->contained_in.begin(); node_it != c_int->contained_in.end(); ) {

            if ( (*node_it)->marker == overlap_node::ELIMINATED) {
                (*node_it)->marker = overlap_node::VACANT;
                node_it = c_int->contained_in.erase(node_it);
            } else {
                (*node_it)->marker = overlap_node::VACANT;
                ++node_it;
            }       
        }
    }  
}

void base_manager::reduce_single_nodes( std::deque<overlap_node> &nodes,  std::deque<contained_node> &contained, graph_list<exon_group> &raw_exons,
        gmap<exon_group*, overlap_node*> &circle_frag_to_overlap, gmap<exon_group*, contained_node*> &circle_frag_to_contained) {
    

    //unsigned int size = nodes.size();
    //unsigned int i = 0;
    for(std::deque<overlap_node>::iterator it = nodes.begin() ; it < nodes.end(); ++it) {

        if ( it->activated && it->links.size() == 1 && it->links.front()->back_links.size() == 1) {
            
            raw_exons.push_back(*it->exons);
            raw_exons.back().extended = true;
            
            // create a new node representing all now joined nodes
            // for this we create a raw exon
          //  nodes.push_back( overlap_node( &raw_exons.back() ) ) ; this distorts order, we need to us current node!
            
            //the old exons are kept as contained nodes;
            contained.push_back( contained_node( &raw_exons.back() ) );
                     
            if (it->exons->backlink_fragment) {
                circle_frag_to_overlap.erase(it->exons);
                circle_frag_to_contained[it->exons] = &contained.back();
            }
                
            // note that active counter is set automatically to 1
            contained.back().add_contained(&*it);
         //   nodes.back().contains.push_back(&contained.back()); // not really needed but for debug
        
            // update contain_in for contained nodes,  node contains unchanged
//            overlap_node* joined = &*it;
//            joined->contains.push_back(&contained.back());
//            for( typename graph_list<contained_node*>::iterator c_it = it->contains.begin(); c_it != it->contains.end(); ++c_it) {
//                (*c_it)->replace_contained(&*it, joined);
//                
//                #ifdef ALLOW_DEBUG
//                joined->contains.push_back(*c_it); // probably not needed, but for debugging
//                #endif
//            }
          
              // don't need this block anymore since we change the node inplace!
//            // update links of previous
//            for( typename graph_list<overlap_node*>::iterator c_it = it->back_links.begin(); c_it != it->back_links.end(); ++c_it) {
//                (*c_it)->replace_link(&*it, &nodes.back()); 
//            }
//            nodes.back().back_links = it->back_links;
            
            raw_exons.push_back(*it->exons); // we need to keep the references to the original intact with all values!
            
            it->exons = &raw_exons.back();
            it->exons->reset_maps();
            it->exons->extended = true;
                   
            // make the current overlap node a new joined node 
            follow_single_reduce(*it->links.back(), *it, contained, circle_frag_to_overlap, circle_frag_to_contained);
            
            
        }
    }
    
}

void base_manager::follow_single_reduce( overlap_node& current, overlap_node& join,  std::deque<contained_node> &contained,
        gmap<exon_group*, overlap_node*> &circle_frag_to_overlap, gmap<exon_group*, contained_node*> &circle_frag_to_contained) {
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Follow reduce\n");
    #endif
    
    // update the exon mask of the joined node
    join.exons->update_mask(*current.exons);
    // create contained node for next element
    contained.push_back( contained_node( current.exons ) );
    
    if (current.exons->backlink_fragment) {
                circle_frag_to_overlap.erase(current.exons);
                circle_frag_to_contained[current.exons] = &contained.back();
    }
    
    contained.back().add_contained(&join);
//    join.contains.push_back(&contained.back()); // probably not needed
    
    for( typename graph_list<contained_node*>::iterator c_it = current.contains.begin(); c_it != current.contains.end(); ++c_it) {
        (*c_it)->replace_contained(&current, &join);
        #ifdef ALLOW_DEBUG
        if (std::find(join.contains.begin(), join.contains.end(), *c_it) == join.contains.end())
        join.contains.push_back(*c_it); // probably not needed, but for debugging
        #endif
    }
      
    current.activated = false;
    
    if (current.links.size() != 1 || current.links.front()->back_links.size() != 1) {
        // end condition for reducing, but this node is still part
        // we can copy over links
        join.links = current.links;
        
        
        for( graph_list<overlap_node*>::iterator link_it = current.links.begin(); link_it != current.links.end(); ++link_it) {
            (*link_it)->replace_backlink(&current, &join);
        }
        return;  
    } 
    
    follow_single_reduce(*current.links.back(), join, contained, circle_frag_to_overlap, circle_frag_to_contained);
}

void base_manager::create_final_graph( std::deque<overlap_node> &nodes,  std::deque<contained_node> &contained,  std::deque<contained_node> &filtered_contained,
        gmap<exon_group*, overlap_node*> &circle_frag_to_overlap, gmap<exon_group*, contained_node*> &circle_frag_to_contained) {
    
   
    logger::Instance()->info("In Create; "+ chromosome + "\n");
    
    // we mark sources and sinks by missing backlinks
    
    s = g.addNode();
    t = g.addNode();
        
    // since lemon does not support node pointer we need additional structure
    bool* gnode_set = new bool[raw->size]();
    // using pointer array is faster than looking into a lemon graph
    bool* source_set = new bool[raw->size]();
        // using pointer array is faster than looking into a lemon graph
    bool* drain_set = new bool[raw->size]();
    
     // room to temporarily storing all nodes for the graph
    ListDigraph::Node* gnodes = new ListDigraph::Node[raw->size];
    
    gmap<exon_edge, ListDigraph::Arc> junction_to_arc;
    
    // lemon maps here work on dynamic allocated arrays, so this should be pretty good
    ListDigraph::ArcMap<count_raw_edge> edge_counts(g);
    ListDigraph::ArcMap<arc_range*> edge_range(g);
    ListDigraph::NodeMap<count_raw_node> node_counts(g);
    
    // for easier memory management, we just keep all ranges in a local list here
    // we need no mystical singleton memory object or slow pointer counters
    std::deque<arc_range> range_list;
    
    // ####### First pass, add overlap_nodes to graph
     
    for(std::deque<overlap_node>::iterator node_it = nodes.begin(); node_it != nodes.end(); ++node_it) {
        
        // we add all normal overlaps first
        if (!node_it->activated) {
            // ignore deactivated nodes
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Add Node " + node_it->exons->bin_mask.to_string() + " "  + std::to_string(node_it->back_links.empty()) + " " + std::to_string(node_it->exons->source_evidence) + " " + std::to_string(node_it->exons->drain_evidence)  +"\n");
        for(gmap<int, exon_group_count>::iterator ssi = node_it->exons->count_series.begin(); ssi != node_it->exons->count_series.end(); ++ssi) {
            logger::Instance()->debug("Id " +  std::to_string(ssi->first) + ":" + std::to_string(ssi->second.total_lefts) + "-" + std::to_string(ssi->second.total_rights)  +"\n");
        }
        #endif
        
        
        
//        for( typename graph_list<overlap_node*>::iterator b_it = node_it->back_links.begin(); b_it != node_it->back_links.end(); ++b_it) {
//             logger::Instance()->debug("Backlink ");
//
//             for(arc_range* r_it = (*b_it)->ranges.begin() ; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {
//                logger::Instance()->debug( std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + ";");
//             }
//             logger::Instance()->debug("\n");
//        }
          
        // if this is a source node, insert it completely
        if (node_it->back_links.empty()) {

            unsigned int start = node_it->exons->range_start;
            unsigned int end = node_it->exons->range_end;
            
            // first create node and source link to it
            
            // two disjoint nodes CAN start at the same graph node
            if (!gnode_set[start]) {
                gnodes[start] = g.addNode();
                gnode_set[start] = true;
            }
            
            if (!source_set[start]) {
                // add source to artifical source
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added Source1\n");
                #endif
                
                ListDigraph::Arc source_arc = g.addArc(s, gnodes[start]);
                initialize_source_drain_arc(source_arc);
                source_set[start] = true;
            }
            
            // test if the end node exists and connect to drain
            if (!gnode_set[end]) {
                gnodes[end] = g.addNode();
                gnode_set[end] = true;
            }
             
            // test if this also qualifies as a drain
            if ((node_it->links.empty() || node_it->exons->drain_evidence) && !drain_set[end]) {
                // add source to artifical source
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added Drain1 "+ std::to_string(node_it->exons->drain_evidence)+ "\n");
                #endif
                
                ListDigraph::Arc drain_arc = g.addArc(gnodes[end], t);
                initialize_source_drain_arc(drain_arc);
                drain_set[end] = true;
            }
            
            // now translate the contained node to a binary exon_edge
            ListDigraph::Arc added_arc = g.addArc(gnodes[start], gnodes[end]); // this also creates the map entries already!
            exon_edge* edge = &ai[added_arc].edge_specifier;
            *edge = node_it->exons->bin_mask;
            
            junction_to_arc.insert(std::make_pair(*edge, added_arc));
            
            range_list.push_back( arc_range(start, end, added_arc));
                      
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Add Range " + std::to_string(range_list.back().start) + " " + std::to_string(range_list.back().end) + "\n");
            #endif
            
            node_it->ranges.addRange(&range_list.back());
            edge_range[added_arc] = &range_list.back();
            
            // add up counts
            count_raw_edge* rec = &edge_counts[added_arc];
            rec->add_initial_count(node_it->exons);
            
            count_raw_node* rc = &node_counts[gnodes[start]];
            rc->add_node_start(node_it->exons);
            
            if (start!= end) {
                rc = &node_counts[gnodes[end]];
                rc->add_node_end(node_it->exons);
            }
            
        } else {
                        
            
            typename graph_list<overlap_node*>::iterator b_it = node_it->back_links.begin();
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Backlink Node " + (*b_it)->exons->bin_mask.to_string() + "\n");
            #endif
            
            // ########## add in the edge for the first overlap
            
            unsigned int i = node_it->exons->range_start; // the index of the first split
            unsigned int j = (*b_it)->exons->range_end; // the index of the second split
            
            arc_range* r_it = (*b_it)->ranges.begin();
            
            #ifdef ALLOW_DEBUG
             logger::Instance()->debug("Split at i " + std::to_string(i) + " " + std::to_string(j)  + "\n");
            #endif
            
            for(; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Range1 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
                #endif
                
                if (r_it->start <= i && r_it->end > i) {
                    //we found the first arc to possibly separate
                    break;
                }
            }
            
            #ifdef ALLOW_DEBUG
            if (r_it!=NULL) logger::Instance()->debug("Range found1 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n"); else logger::Instance()->debug("Range found1 NULL \n");
            #endif
             
            // if end == i we already hit a node, isn't that nice, no splitting
            if (r_it!=NULL && r_it->start != i) {
                
                // split up graph and edge objects
                split_edge(i, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
                
            }
            
            // ########## from i to range end counts need to be updated
            
            // since we are handling overlaps, the current node on insert MUST end after current range
            count_raw_node* rcl = &node_counts[gnodes[i]];
            rcl->add_node_start(node_it->exons); 
           
            unsigned int offset = 0;
            gmap<int, rcount> next_value;
            for (gmap<int, exon_group_count>::iterator egi = node_it->exons->count_series.begin(); egi != node_it->exons->count_series.end(); ++egi) {
                next_value[egi->first] = egi->second.total_lefts + egi->second.hole_end_counts[0] - egi->second.hole_start_counts[0];
            }
            
            for(; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {
                
                count_raw_edge* rec = &edge_counts[r_it->arc];
                rec->add_sub_counts_eg(node_it->exons, offset+1, offset + rec->size, next_value);
                offset += rec->size + 1;
                
                count_raw_node* rcl = &node_counts[gnodes[r_it->end]];
                rcl->add_node_initial_index(node_it->exons, offset, next_value);
                
                node_it->ranges.addRange(r_it);
            }
            
           // ########## from j to node end add new edge 
            
            // start node j HAS to exist
            // test if the end node exists and connect to drain
            unsigned int end = node_it->exons->range_end;
            if (!gnode_set[end]) {
                gnodes[end] = g.addNode();
                gnode_set[end] = true;
            }
            // if marked this is also a source end also qualifies as a drain
            if (node_it->exons->source_evidence && !source_set[i]) {
                // add source to artifical source
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added Source2\n");
                #endif
                
                ListDigraph::Arc source_arc = g.addArc(s, gnodes[i]);
                initialize_source_drain_arc(source_arc);
                source_set[i] = true;
            }
            
            // test if end also qualifies as a drain
            if ((node_it->links.empty() || node_it->exons->drain_evidence) && !drain_set[end]) {
                // add source to artifical source
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added Drain2 " +std::to_string(node_it->exons->drain_evidence) +"\n");
                #endif
                
                ListDigraph::Arc drain_arc = g.addArc(gnodes[end], t);
                initialize_source_drain_arc(drain_arc);
                drain_set[end] = true;
            }
            
            // now translate the area node to a binary exon_edge
            exon_edge edge;
            node_it->exons->bin_mask.right_split(j, edge); // creates completely empty edge evidence
            
            ListDigraph::Arc new_arc;
            arc_range* new_range;
            // test if last arc exits, if not add it!
            typename gmap<exon_edge, ListDigraph::Arc>::iterator e = junction_to_arc.find(edge);
            if (e == junction_to_arc.end()) {
                // so the arc really does not exist already, so make it
                new_arc = g.addArc(gnodes[j], gnodes[end]);
                range_list.push_back(arc_range(j, end, new_arc));
                new_range = &range_list.back();
                
                e = junction_to_arc.insert(std::make_pair(edge, new_arc)).first;
                ai[new_arc].edge_specifier = edge;
                edge_range[new_arc] = new_range;
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added new Arc\n");
                #endif
                
            } else {
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Found Arc\n");
                #endif
                
                new_arc = e->second;
                new_range = edge_range[new_arc];
            }
            
            count_raw_edge* rec = &edge_counts[new_arc];
            rec->add_sub_counts_eg(node_it->exons, offset+1, offset + edge.id.count() - 2, next_value);
            
            rcl = &node_counts[gnodes[end]];
            rcl->add_node_end(node_it->exons);
            
           // range_list.push_back( arc_range(j, end, new_arc)); double add, why was this here?
            node_it->ranges.addRange(new_range);
            
            // ########## new node has been inserted into first overlap, now snap arcs to other overlaps! 
            ++b_it;
            for (; b_it != node_it->back_links.end(); ++b_it) {
                
                // ########## find arc of new b that contains first index (i is still valid and the same)
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Split at i " + std::to_string(i) + " " + std::to_string(j) + " " + (*b_it)->exons->bin_mask.to_string() + "\n");
                #endif
                
                arc_range* r_it = (*b_it)->ranges.begin();            
                for(; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {
                
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Test Range2 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
                    #endif
                    
                    if (r_it->start <= i && r_it->end > i) {
                        //we found the first arc to possibly separate
                        break;
                    }
                }

                #ifdef ALLOW_DEBUG
                if (r_it!=NULL) logger::Instance()->debug("Range found2 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n"); else logger::Instance()->debug("Range found2 NULL \n");
                #endif
                
                // if start == i we already hit a node, isn't that nice, no splitting
                if (r_it!=NULL && r_it->start != i) {
                    // split up graph and edge objects
                    split_edge(i, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
                    // r_it now points at the new created edge starting at i!
                }
                
                // ########## now we move both ranges down in tandem to find stuff that needs joining!
                arc_range* i_it = node_it->ranges.begin();
                while ( r_it != NULL) { // overlap HAS to end first
                                        
                    // at this point the start of both iterators have to have the same start
                    if (i_it->end == r_it->end) {
                        
                        // ranges are the same (i_it == r_it), no need for any changes
                        // e.g. while splitting all counts have already be joined
                        i_it = node_it->ranges.next();
                        r_it = (*b_it)->ranges.next();
                    } else {
                        // so we need a split, but which edge
                        if (i_it->end < r_it->end) {
                            // split overlap on inserted node
                            if (i_it->end > r_it->start)
                                split_edge(i_it->end, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
                            i_it = node_it->ranges.next();
                            
                        } else {
                            // split the inserted node on overlap
                            if (r_it->end > i_it->start)
                                split_edge(r_it->end, node_it->ranges, i_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
                            r_it = (*b_it)->ranges.next();
                        }
                    }
                }
            }
        } // end if back_links
    } // end iterator over overlaps
    
    
   logger::Instance()->info("First Pass Done ++++++++++++++++++++++++++++++++ "+ chromosome + "\n");

    #ifdef ALLOW_DEBUG
    if (options::Instance()->is_debug())
        print_graph_debug(std::cout, node_counts, edge_counts);
    #endif
    
    // ####### second pass, add all contained with active count over 1 (graph changing contains)
    
    for(std::deque<contained_node>::iterator node_it = contained.begin(); node_it != contained.end(); ++node_it) {
        
        if (node_it->contained_in.size() < 2) {
            continue;
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Contained Node " + node_it->exons->bin_mask.to_string() + "\n");
        for (typename graph_list<overlap_node* >::iterator t = node_it->contained_in.begin(); t!= node_it->contained_in.end(); ++t ) {
            logger::Instance()->debug("Contained In: " + (*t)->exons->bin_mask.to_string() + "\n");
        }
        #endif
        
        // ####### add in contain to first overlap, afterwards snap others to created arcs
        
        typename graph_list<overlap_node*>::iterator b_it = node_it->contained_in.begin();
        
        // ########## add in the edge for the first overlap
            
        unsigned int i = node_it->exons->range_start; // the index of the first split
        unsigned int j = node_it->exons->range_end;   // the index of the second split
    
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Split at i " + std::to_string(i) + " " + std::to_string(j)  + "\n");
        #endif
        
        arc_range* r_it = (*b_it)->ranges.begin();

        for(; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Range3 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
            #endif
            
            if (r_it->start <= i && r_it->end > i) {
                //we found the first arc to possibly separate
                break;
            }
        }
        
        #ifdef ALLOW_DEBUG
        if (r_it!=NULL) logger::Instance()->debug("Range found3 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n"); else logger::Instance()->debug("Range found3 NULL \n");
        #endif
        
        // if start == i we already hit a node, isn't that nice, no splitting
        if (r_it!=NULL && r_it->start != i) {
            // split up graph and edge objects
            split_edge(i, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
        }
            
         if (node_it->exons->source_evidence && !source_set[i]) {
                // add source to artifical source
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added Source c3\n");
                #endif
                
                ListDigraph::Arc source_arc = g.addArc(s, gnodes[i]);
                initialize_source_drain_arc(source_arc);
                source_set[i] = true;
        }
        
        // ########## from i to j counts need to be updated, with possible split at j
        
        // since we are handling contains, the current node on end in b_it
        count_raw_node* rcl = &node_counts[gnodes[i]];
        rcl->add_node_start(node_it->exons); 

        // we can have single exon templates here, no further is needed for them    
        if (i!=j) { 
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Split at j " + std::to_string(i) + " " + std::to_string(j)  + "\n");
            #endif
            
            unsigned int offset = 0;
            gmap<int, rcount> next_value;
            for (gmap<int, exon_group_count>::iterator egi = node_it->exons->count_series.begin(); egi != node_it->exons->count_series.end(); ++egi) {
                next_value[egi->first] = egi->second.total_lefts + egi->second.hole_end_counts[0] - egi->second.hole_start_counts[0];
            }
            
            for(; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Range4 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
                #endif

                if (r_it->start <= j && r_it->end > j) {
                    //we found the first arc to possibly separate
                    break;
                }

               count_raw_edge* rec = &edge_counts[r_it->arc];
               rec->add_sub_counts_eg(node_it->exons, offset+1, offset + rec->size, next_value);
               offset += rec->size + 1;
               
                count_raw_node* rcl = &node_counts[gnodes[r_it->end]];
                if (r_it->end == j) {
                    rcl->add_node_end(node_it->exons);
                } else {
                    rcl->add_node_initial_index(node_it->exons, offset, next_value);
                }

//                logger::Instance()->debug(" Range ADD 1 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
                node_it->ranges.addRange(r_it);
            } 
            
            #ifdef ALLOW_DEBUG
            if (r_it!=NULL) logger::Instance()->debug("Range found4 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n"); else logger::Instance()->debug("Range found4 NULL \n");
            #endif
            
            // if start == j we already hit a node, isn't that nice, no splitting
            if (r_it!=NULL && r_it->start != j) { 
                    
                // split up graph and edge objects
                split_edge(j, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);

                //since r_it now points to second half of inserted edge, go one back
                r_it = (*b_it)->ranges.previous();  
            
                count_raw_edge* rec = &edge_counts[r_it->arc];
                rec->add_sub_counts_eg(node_it->exons, offset+1, offset + rec->size, next_value);
                offset += rec->size + 1;
                
                count_raw_node* rcl = &node_counts[gnodes[r_it->end]];
                rcl->add_node_end(node_it->exons);
                
//                logger::Instance()->debug(" Range ADD 2 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
                node_it->ranges.addRange(r_it);
            }
        }
        
        if (node_it->exons->drain_evidence && !drain_set[j]) {
                // add source to artifical source
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Added Drain c3 "+ std::to_string(node_it->exons->drain_evidence)+ "\n");
                #endif
                
                ListDigraph::Arc drain_arc = g.addArc(gnodes[j], t);
                initialize_source_drain_arc(drain_arc);
                drain_set[j] = true;
        }
        
        // now snap others to this
        ++b_it;
        for (; b_it != node_it->contained_in.end(); ++b_it) {
            
            arc_range* r_it = (*b_it)->ranges.begin();

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Split at i " + std::to_string(i) + " " + std::to_string(j)  + "\n");
            #endif
            
            for(; r_it!=NULL ; r_it = (*b_it)->ranges.next() ) {
           
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Range5 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n");
                #endif
                
                if (r_it->start <= i && r_it->end > i) {
                    //we found the first arc to possibly separate
                    break;
                }
            }

            #ifdef ALLOW_DEBUG
            if (r_it!=NULL) logger::Instance()->debug("Range found5 " + std::to_string(r_it->start) + " " + std::to_string(r_it->end)  + "\n"); else logger::Instance()->debug("Range found5 NULL \n");
            #endif
            
            // if start == i we already hit a node, isn't that nice, no splitting
            if (r_it!=NULL && r_it->start != i) {
                // split up graph and edge objects
                split_edge(i, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
            }
            
            // ########## now we move both ranges down in tandem to find stuff that needs joining!
            
            arc_range* i_it = node_it->ranges.begin();
            for(; i_it!=NULL ; i_it = node_it->ranges.next() ) {
           
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Before Op " + std::to_string(i_it->start) + " " + std::to_string(i_it->end)  + "\n");
                #endif
                
            }
            
            i_it = node_it->ranges.begin();
            while ( i_it != NULL) { // contain HAS to end first, if i=j this is NULL

//                logger::Instance()->debug("Loop" + std::to_string(r_it->start) + " " + std::to_string(r_it->end) + " - " + std::to_string(i_it->start) + " " + std::to_string(i_it->end)+"\n");

                // at this point the start of both iterators have to have the same start
                if (i_it->end == r_it->end) {
                    // ranges are the same (i_it == r_it), no need for any changes
                    // e.g. while splitting all counts have already be joined
                    i_it = node_it->ranges.next();
                    r_it = (*b_it)->ranges.next();
//                    logger::Instance()->debug(" op1.1 " + std::to_string(r_it->start) + "," + std::to_string(r_it->end) +"\n");
//                    if (i_it != NULL) logger::Instance()->debug(" op1.2 " + std::to_string(i_it->start) + " " + std::to_string(i_it->end) + "\n");
//                    
                } else {
                    // so we need a split, but which edge
                    if (i_it->end < r_it->end) {
                        // split overlap on inserted node
                        
//                        logger::Instance()->debug(" op2 \n");
                        
                        if (i_it->end > r_it->start)
                            split_edge(i_it->end, (*b_it)->ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
                        i_it = node_it->ranges.next();

                    } else {
                        // split the inserted node on overlap
                        
//                        logger::Instance()->debug(" op3 \n");
                        
                        if (r_it->end > r_it->start)
                            split_edge(r_it->end, node_it->ranges, i_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
                        r_it = (*b_it)->ranges.next();
                    }
                }
            }
        } 
    }
    
    
    logger::Instance()->info("Second Pass Done ++++++++++++++++++++++++++++++++ " + chromosome + "\n");
    
    #ifdef ALLOW_DEBUG
    if (options::Instance()->is_debug()) print_graph_debug(std::cout, node_counts, edge_counts);;
    #endif
    
    // ####### third pass, add all contained with active count 1 (graph unchanging contains)
    
    for(std::deque<contained_node>::iterator node_it = contained.begin(); node_it != contained.end(); ++node_it) {
        
        if (node_it->contained_in.size() != 1) {
            continue;
        }
               
        overlap_node* overlap = node_it->contained_in.back();
        add_unaltering_contained_counts(overlap, node_it, gnodes, edge_counts, node_counts);
        if (options::Instance()->is_debug()) print_graph_debug(std::cout, node_counts, edge_counts);
    }

    // ####### fourth pass, add all nodes that where excluded from the direct overlap test (graph unchanging smaller, e.g clipped or differently sequenced, reads)
    
    // by definition these fragments might not actually fit a unique edge
    // we will add it to all as a heuristic, bordering values should take care of it
    // TODO: is this really true? Seems like a nice guess at least <.< 
    
    for(std::deque<contained_node>::iterator node_it = filtered_contained.begin(); node_it != filtered_contained.end(); ++node_it) {
        
        for(typename graph_list<overlap_node*>::iterator overlap = node_it->contained_in.begin(); overlap != node_it->contained_in.end(); ++overlap) { 
             add_unaltering_contained_counts(*overlap, node_it, gnodes, edge_counts, node_counts);
        }
    }
  
    logger::Instance()->info("Third and Fourth Pass Done ++++++++++++++++++++++++++++++++\n");
    
    // ####### the graph is basically done here in it's basics
    // what is left to do is:   - correct the starting and or end bias wherever possible (later in flow model)
    //                          - split on a) known starts/end in sub b) cycles c) large jumps in coverage
     
    
    #ifdef ALLOW_DEBUG
    if (options::Instance()->is_debug()) print_graph_debug(std::cout, node_counts, edge_counts);
    #endif
   
    // we can kill no longer needed stuff
    //delete junction_to_arc;
    //delete nodes;
    //delete contained;
    
    
    // ####### 5th pass add in CIRCULAR INFORMATION -> breaks the DAG ########
    
//    for ( graph_list<paired_exon_group>::iterator b_it = raw->chim_circle_bin_list.begin(); b_it != raw->chim_circle_bin_list.end(); ++b_it) {
//        
//        // get the left read (actually mapped right of both)
//        
//        ListDigraph::Node jump_start, jump_end;
//        
//        // read can be contained or overlap
//        typename gmap<exon_group*, overlap_node*>::iterator e = circle_frag_to_overlap.find(b_it->left_read);
//        if (e == circle_frag_to_overlap.end()) {
//             // since this is the case, we have a guarantee that the node exists and can directly use it
//             jump_start = gnodes[e->second->exons->range_end];
//        } else {
//            // this has to be in contained
//            contained_node* cont = circle_frag_to_contained.find(b_it->left_read)->second;
//            if (cont->contained_in.size() > 2) {
//                 // compacted down other edges, so node has to exist!
//                 jump_start = gnodes[e->second->exons->range_end];
//            } else {
//                 // complicated testcase, could be within edge
//                 
//                 unsigned int i = cont->exons->range_end;
//                 
//                 arc_range* r_it = cont->contained_in.front()->ranges.begin();
//                 for(; r_it!=NULL ; r_it = cont->contained_in.front()->ranges.next() ) {
//
//                    if (r_it->start <= i && r_it->end >= i) {
//                        //we found the first arc to possibly separate
//                        break;
//                    }
//                }
//                
//                if (i == r_it->start) {
//                    // again, node exists, but we can't be sure if it's a normal node or not
//                    
//                    jump_start = g.source(r_it->arc);  
//                } else if (i ==  r_it->end) {
//                    
//                    jump_start = g.target(r_it->arc);
//                } else {
//                    // we need to insert special node into otherwise nice edge, what a shame
//                    split_edge_without_compacting(i, cont->contained_in.front()->ranges, r_it, edge_specifier, edge_range, edge_counts, node_counts, range_list);
//                    jump_start = g.source(r_it->arc);
//                }
//                 
//            }
//        }
//        
//        // read can be contained or overlap
//        e = circle_frag_to_overlap.find(b_it->right_read);
//        if (e == circle_frag_to_overlap.end()) {
//             // since this is the case, we have a guarantee that the node exists and can directly use it
//             jump_end = gnodes[e->second->exons->range_start];
//        } else {
//            // this has to be in contained
//            contained_node* cont = circle_frag_to_contained.find(b_it->right_read)->second;
//            if (cont->contained_in.size() > 2) {
//                 // compacted down other edges, so node has to exist!
//                 jump_end = gnodes[e->second->exons->range_start];
//            } else {
//                 // complicated testcase, could be within edge
//                 
//                 unsigned int i = cont->exons->range_start;
//                 
//                 arc_range* r_it = cont->contained_in.front()->ranges.begin();
//                 for(; r_it!=NULL ; r_it = cont->contained_in.front()->ranges.next() ) {
//
//                    if (r_it->start <= i && r_it->end >= i) {
//                        //we found the first arc to possibly separate
//                        break;
//                    }
//                }
//                
//                if (i == r_it->start) {
//                    // again, node exists, but we can't be sure if it's a normal node or not
//                    
//                    jump_start = g.source(r_it->arc);  
//                } else if (i ==  r_it->end) {
//                    
//                    jump_start = g.target(r_it->arc);
//                } else {
//                    // we need to insert special node into otherwise nice edge, what a shame
//                    split_edge_without_compacting(i, cont->contained_in.front()->ranges, r_it, edge_specifier, edge_range, edge_counts, node_counts, range_list);
//                    jump_start = g.source(r_it->arc);
//                }
//            }
//        }
//        
//        // we have the nodes, test if arc already exist between them
//        ListDigraph::OutArcIt a(g, jump_start);
//        for (; a != INVALID; ++a) {
//            if ( g.target(a) == jump_end ) {
//                break;
//            } 
//        }
//        
//        ListDigraph::Arc arc;
//        if (a == INVALID) {
//            // does not exist, just add and have fun
//            arc = g.addArc(jump_start, jump_end);
//        } else {
//            arc = a;
//        } 
//        
//        //MASSIV TODO
//       // edge_counts[arc].add_count(b_it->count, b_it->count, 0);
//        
//    }
//    
//    // we need to identify all cycles with a unique id
//    unsigned int id = 1; // we start at 0
//    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
//        ListDigraph::Arc arc(a);
//        if (edge_type[arc] == edge_types::BACKLINK ) {
//            cycle_id_in[arc] = id;      
//            ++id;
//        }
//    }
//            
    //delete edge_range;
       
    // now look for areas of high coverage jumps that are NOT at a junction
    std::set<int> force_split;
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        if(edge_counts[a].size < 1) { // so at least 2 splits!
            continue;
        }
        
        #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Splitting " + std::to_string(g.id(a)) + "\n");
        #endif
        
        range_helper rh;   // we have no overlap nodes here anymore, so we create a new one
        arc_range* ar = edge_range[a];
        rh.addRange(ar);

        count_raw_edge* re = &edge_counts[a];
        // good arc, test for sudden increases!
        
        exon_edge edge = ai[a].edge_specifier;
        
        for (gmap<int, count_raw_edge::series_struct>::iterator ssi = re->series.begin(); ssi != re->series.end(); ++ssi) {
        
            std::deque<unsigned int> jump_end;
            std::deque<unsigned int> split_index_end;
            std::deque<unsigned int> jump_start;
            std::deque<unsigned int> split_index_start;
            
            std::vector<rcount>::iterator ri = ssi->second.splits.begin();
            unsigned int c = 0;
            unsigned int index = edge.id.find_first();
            rcount lastval = *ri;
            unsigned int last_index = index;
            ++ri;
            ++c;
            index = edge.id.find_next(index);
            for (; ri != ssi->second.splits.end(); ++ri, ++c) {

                rcount n1 = std::max(lastval, *ri);
                rcount n2 = std::min(lastval, *ri);

                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Test Split " + std::to_string(n1) + " " + std::to_string(n2) + " " + std::to_string((n1 - n2) * 100 / float(n2))+ "\n");
                #endif

                if (n2 <= options::Instance()->get_low_edge_mark() && n1 >= options::Instance()->get_low_edge_mark() + 1) {
                    force_split.insert(index);
                }

                if ( (n1 - n2) * 100 / float(n2) >= options::Instance()->get_coverage_change_limit() ) {

                    if (lastval > *ri) {
                        // down trend
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("add end " + std::to_string(c) + "\n");
                        #endif

                        jump_end.push_back(c);
                        split_index_end.push_back(index);
                    } else {

                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("add start " + std::to_string(c-1) + "\n");
                        #endif

                        jump_start.push_back(c-1);
                        split_index_start.push_back(last_index);
                    }
                }
                lastval = *ri;
                last_index = index;
                index = edge.id.find_next(index);
            }

            // post processing found sites
            unsigned int last = ssi->second.splits.size() + 1; // can't possibly be in there
            std::deque<unsigned int>::iterator spi = split_index_start.begin();
            for (std::deque<unsigned int>::iterator  si = jump_start.begin(); si != jump_start.end(); ++si, ++spi) {
                if (*si == 0 || *si - 1 != last ) {
                    // serious start evidence!

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("CSTART" + std::to_string(*spi) + "\n");
                    #endif

                    force_split.insert(*spi);
                }
                last = *si;
            }
            // post processing found sites
            last = ssi->second.splits.size() + 5; // can't possibly be in there
            std::deque<unsigned int>::reverse_iterator spr = split_index_end.rbegin();
            for (std::deque<unsigned int>::reverse_iterator  si = jump_end.rbegin(); si != jump_end.rend(); ++si, ++spr) {
                if (*si + 1 != last ) {
                    // serious end evidence!

                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("CEND" + std::to_string(*spr) + "\n");
                    #endif

                    force_split.insert(*spr);
                }
                last = *si;
            }
        }
    }
    
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
         
        rcount left = 0;
        bool has_left = false;
        for (ListDigraph::InArcIt a(g,n); a != INVALID; ++a) {
            if (g.source(a) == s || g.target(a) == t) {
                continue;
            }
            
            count_raw_edge* re = &edge_counts[a];
            for (gmap<int, count_raw_edge::series_struct>::iterator ssi = re->series.begin(); ssi != re->series.end(); ++ssi) {
                rcount v = std::max(ssi->second.splits.front(), ssi->second.splits.back()); 
                if (!has_left || v > left) {
                    left = v;
                    has_left = true;
                }
            }
        }
        
        rcount right = 0;
        bool has_right = false;
        unsigned int index;
        for (ListDigraph::OutArcIt a(g,n); a != INVALID; ++a) {
            
            if (g.source(a) == s || g.target(a) == t) {
                continue;
            }
            
            count_raw_edge* re = &edge_counts[a];
            for (gmap<int, count_raw_edge::series_struct>::iterator ssi = re->series.begin(); ssi != re->series.end(); ++ssi) {
                rcount v = std::max(ssi->second.splits.front(), ssi->second.splits.back()); 
                if (!has_right || v > right) {
                    right = v;
                    if (!has_right) {
                        exon_edge edge = ai[a].edge_specifier;
                        index = edge.id.find_first();
                        has_right = true;
                    }
                }
            }
        }
        
        rcount n1 = std::max(left, right);
        rcount n2 = std::min(left, right);
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Node Force Test " + std::to_string(g.id(n)) + " " + std::to_string(n1) + " " + std::to_string(n2) + "\n");
        #endif
        
        if (has_left && has_right && n2 <= options::Instance()->get_low_edge_mark() && n1 >= options::Instance()->get_low_edge_mark()+1) {
           
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Force " + std::to_string(index) + "\n");
            #endif
            
            force_split.insert(index);
        }
    }
    
    // force snaps
    while (true) {   
	    bool cont = false;       
	    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
		
		ListDigraph::Arc arc(a);
		++a;
		
		if(edge_counts[arc].size < 1) { // so at least 2 splits!
		    continue;
		}

	 	#ifdef ALLOW_DEBUG
		logger::Instance()->debug("Snap Edge " + std::to_string(g.id(arc)) +  " "+ ai[arc].edge_specifier.to_string() + "\n");
		if (options::Instance()->is_debug()) print_graph_debug(std::cout, node_counts, edge_counts);;
		#endif        

		exon_edge edge = ai[arc].edge_specifier;
		unsigned int index = edge.id.find_first();
		index = edge.id.find_next(index);
		
		range_helper rh;   // we have no overlap nodes here anymore, so we create a new one
		arc_range* ar = edge_range[arc];
		rh.addRange(ar);
		arc_range* r_it = rh.begin();
		
		unsigned int size = edge_counts[arc].size;
		
		bool was_split = false;
		for (unsigned int i = 0; i < size; ++i) {

		    #ifdef ALLOW_DEBUG
		    logger::Instance()->debug("splitnew " + edge.to_string() + " " + std::to_string(index) + " " + std::to_string(gnode_set[index]) + " s " + std::to_string(source_set[index]) + " d " + std::to_string(drain_set[index]) + "\n");
		    #endif

		    bool is_force = force_split.find(index) != force_split.end();
		    was_split = was_split || is_force;		            
    
		    if (is_force ) {
		        split_edge(index, rh, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);
		    }
		    index = edge.id.find_next(index);
		}
		if (was_split) {
		    cont = true;
		    break;
		}
	    } 
	    if (!cont) {
		break;
	    }
    }  
    
     #ifdef ALLOW_DEBUG
    if (options::Instance()->is_debug()) print_graph_debug(std::cout, node_counts, edge_counts);
    #endif

    //  use known splits
    for(std::deque<contained_node>::iterator node_it = contained.begin(); node_it != contained.end(); ++node_it) {
    
        if (node_it->contained_in.size() != 1) {
            continue;
        }
         
        if (node_it->exons->source_evidence && node_it->exons->range_start != node_it->contained_in.front()->exons->range_start) {
            overlap_node* overlap = node_it->contained_in.back();
            add_contained_start(overlap->ranges, node_it->exons->range_start, node_counts, edge_range, edge_counts, range_list, true);
        }
        if (node_it->exons->drain_evidence && node_it->exons->range_end != node_it->contained_in.front()->exons->range_end) {
            overlap_node* overlap = node_it->contained_in.back();              
            add_contained_end(overlap->ranges, node_it->exons->range_end, node_counts, edge_range, edge_counts, range_list, true);
        }
    }
   
      
    // snap the rest
    for (ListDigraph::ArcIt a(g); a != INVALID; ) {
        
        ListDigraph::Arc arc(a);
        ++a;
        
        if(edge_counts[arc].size < 1) { // so at least 2 splits!
            continue;
        }
        
        exon_edge edge = ai[arc].edge_specifier;
        unsigned int index = edge.id.find_first();
        index = edge.id.find_next(index);
        
        range_helper rh;   // we have no overlap nodes here anymore, so we create a new one
        arc_range* ar = edge_range[arc];
        rh.addRange(ar);
        arc_range* r_it = rh.begin();
        
        unsigned int size = edge_counts[arc].size;
        
        for (unsigned int i = 0; i < size; ++i) {

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("splitnew " + edge.to_string() + " " + std::to_string(index) + " " + std::to_string(gnode_set[index]) + " s " + std::to_string(source_set[index]) + " d " + std::to_string(drain_set[index]) + "\n");
            #endif

            bool is_source = false;
            bool is_drain = false;
            
            if (gnode_set[index] && ( source_set[index] || drain_set[index]) ) {
                if (source_set[index] ) {
                    is_source = true;
                    ListDigraph::InArcIt a(g, gnodes[index]);
                    if (g.source(a) == s) ++a;
                    if (a != INVALID) is_source = false;
                }
                if (drain_set[index] ) {
                    is_drain = true;
                    ListDigraph::OutArcIt a(g, gnodes[index]);
                    if (g.target(a) == t) ++a;
                    if (a != INVALID) is_drain = false;
                }
            }
            
            if (gnode_set[index] && ( is_source || is_drain) ) {
                split_edge_without_compacting(index, rh, r_it, edge_range, edge_counts, node_counts, gnodes, range_list, true);
            } else { 
                split_edge_without_compacting(index, rh, r_it, edge_range, edge_counts, node_counts, gnodes, range_list, false);
            }
            index = edge.id.find_next(index);
        }
    }   
    

    
    // now join up edges to minimum!
    while (true) {
        bool changed = false;
        
        for (ListDigraph::NodeIt n(g); n != INVALID;) {
       
            #ifdef ALLOW_DEBUG
            if (options::Instance()->is_debug()) print_graph_debug(std::cout, node_counts, edge_counts);
            #endif
            
            ListDigraph::Node node(n);
            ++n;
            
            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("test compact at node " + std::to_string(g.id(node)) + "\n");
            #endif
            
            for (ListDigraph::OutArcIt o1(g, node); o1 != INVALID; ) {
                
                ListDigraph::Arc a1(o1);
                ++o1;
                
                for (ListDigraph::OutArcIt o2(g, node); o2 != INVALID; ) {

                    ListDigraph::Arc a2(o2);
                    ++o2;

                    if (a1 == a2 || g.target(a1) == t || g.target(a2) == t || g.source(a1) == s || g.source(a2) == s) {
                        continue;
                    }
                    
                    if (ai[a1].edge_specifier == ai[a2].edge_specifier) {
                        // join those!
                        
                        unsigned int ic1 = 0;
                        for (ListDigraph::InArcIt it(g, g.target(a1)); it != INVALID; ++it) ++ic1;
                        unsigned int ic2 = 0;
                        for (ListDigraph::InArcIt it(g, g.target(a2)); it != INVALID; ++it) ++ic2;

                        if (ic2 > 1 || ic1 > 1) {
                            continue;
                        }
                        
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("compact R " + std::to_string(g.id(a1)) + " " + std::to_string(g.id(a2)) + "\n");
                        #endif  
                        if (g.target(a1) != g.target(a2)) {
                            node_counts[g.target(a1)].add_node(&node_counts[g.target(a2)]);
                            g.contract(g.target(a1), g.target(a2), false); // dels a2!
                        }
                        edge_counts[a1].add_edge(&edge_counts[a2]);
                        
                        if (o1 == a2) {
                            ++o1;
                        }
                        
                        g.erase(a2);
                        changed = true;
                    } 
                }
            }
            for (ListDigraph::InArcIt o1(g, node); o1 != INVALID; ) {
                
                ListDigraph::Arc a1(o1);
                ++o1;
                
                for (ListDigraph::InArcIt o2(g, node); o2 != INVALID; ) {

                    ListDigraph::Arc a2(o2);
                    ++o2;
                    
                    if (a1 == a2 || g.source(a1) == s || g.source(a2) == s || g.target(a1) == t || g.target(a2) == t) {
                        continue;
                    }
                    
                    if (ai[a1].edge_specifier == ai[a2].edge_specifier) {
                        // join those!
                        
                        unsigned int ic1 = 0;
                        for (ListDigraph::OutArcIt it(g, g.source(a1)); it != INVALID; ++it) ++ic1;
                        unsigned int ic2 = 0;
                        for (ListDigraph::OutArcIt it(g, g.source(a2)); it != INVALID; ++it) ++ic2;

                        if (ic2 > 1 || ic1 > 1) {
                            continue;
                        }
                        
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("compact L " + std::to_string(g.id(a1)) + " " + std::to_string(g.id(a2)) + "\n");
                        #endif  
                        if (g.source(a1) != g.source(a2)) {
                            node_counts[g.source(a1)].add_node(&node_counts[g.source(a2)]);
                            g.contract(g.source(a1), g.source(a2), false); // dels a2!
                        }
                        edge_counts[a1].add_edge(&edge_counts[a2]);

                        if (o1 == a2) {
                            ++o1;
                        }
                        
                        g.erase(a2);
                        changed = true;
                    } 
                }
            }
        }
        
        if (!changed) {
            break;
        }
    }
     
    create_final_capacities(node_counts, edge_counts, gnodes, gnode_set);
    
    delete [] gnode_set;
    delete [] gnodes;
    delete [] source_set;
    delete [] drain_set;
  
    
    // ####### done :)
}

void base_manager::create_final_capacities(
    ListDigraph::NodeMap<count_raw_node> &node_counts,
    ListDigraph::ArcMap<count_raw_edge> &edge_counts,
    ListDigraph::Node* gnodes,
    bool* gnode_set) {
       
    // expand nodes if this is wanted by the implementing manager
    if(expand_exon_nodes()) {
        
        for(unsigned int i=0; i < raw->size;++i) {
            if (gnode_set[i]) {
                node_index[gnodes[i]] = i;
            }
        }
        
        for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {        
                
            if (n == s || n == t) {
                continue;
            } 
            
            unsigned int i = node_index[n];
            
            ListDigraph::Node nn = g.split(n, false);
            ListDigraph::Arc na = g.addArc(n, nn);
                
            ai[na].edge_type = edge_types::NODE;
            
            create_region_from_node(i, node_counts[n], fs[na]);
            create_node_capacities(na, fs[na]);
                   
            node_index[nn] = i;    
            ai[na].edge_specifier.node_index = i;
            ai[na].edge_lengths.middle = meta->exons[i].exon_length;
                               
            // we do not add a specifier to exon edge
        }
        
    } else {
        for(unsigned int i=0; i < raw->size;++i) {
            if (gnode_set[i]) {
                node_index[gnodes[i]] = i;
            }
        }
    }
      
    // now create arc capacities
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        ListDigraph::Arc arc(a);
        
        if (ai[arc].edge_type == edge_types::NODE ) {
            continue;
        }
        
        // test if the edge is a source edge
        if(g.source(arc) == s || g.target(arc) == t) {
            ai[arc].edge_type = edge_types::HELPER;
            continue;
        }
        
        if (ai[arc].edge_type == edge_types::BACKLINK ) {
            create_region_from_edge(ai[arc].edge_specifier, edge_counts[arc], fs[arc]);
            create_edge_capacities(arc, fs[arc]);
            continue;
        }
        
        ai[arc].edge_type = edge_types::EXON;
        compute_edge_length(arc);
        create_region_from_edge(ai[arc].edge_specifier, edge_counts[arc], fs[arc]);
        create_edge_capacities(arc, fs[arc]);
    }

}

void base_manager::compute_edge_length( ListDigraph::Arc arc) {
    
    exon_edge* edge = &ai[arc].edge_specifier;
    edge_length* edge_length = &ai[arc].edge_lengths;
    
    unsigned int index = edge->id.find_first();
    edge_length->first_exon = meta->exons[index].exon_length;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Comp Edgelength " +  edge->to_string() + "\n");
    #endif
    
    // we have at least 2 entries
    unsigned int pre_index = edge->id.find_next(index);
    
    if (meta->exons[index].right + 1 == meta->exons[pre_index].left) {
        ai[arc].edge_specifier.left_consecutive = true;
        ai[arc].edge_specifier.right_consecutive = true;
    }
    
    index = edge->id.find_next(pre_index);
    
    while(index < edge->id.size() ) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Length Index " +  std::to_string(pre_index) + " " + std::to_string(index) + "\n");
        #endif
        edge_length->middle += meta->exons[pre_index].exon_length;
        
        if (meta->exons[pre_index].right + 1 == meta->exons[index].left) {
            ai[arc].edge_specifier.right_consecutive = true;
        } else {
            ai[arc].edge_specifier.right_consecutive = false;
        }
        
        pre_index = index;
        index = edge->id.find_next(index);
    }
    
    edge_length->last_exon = meta->exons[pre_index].exon_length;
}

void base_manager::add_unaltering_contained_counts( overlap_node* overlap,
    std::deque<contained_node>::iterator &node_it,
    ListDigraph::Node* gnodes,
    ListDigraph::ArcMap<count_raw_edge> &edge_counts,
    ListDigraph::NodeMap<count_raw_node> &node_counts) {
    // just one object
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Add Unaltering " +  node_it->exons->bin_mask.to_string() + " in " + overlap->exons->bin_mask.to_string() + "\n");
    #endif    

    unsigned int i = node_it->exons->range_start;
    unsigned int j = node_it->exons->range_end;

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("indices " + std::to_string(i) + " " + std::to_string(j) + "\n");
    #endif 
    
    unsigned int offset_new_exons = 0;
    unsigned int offset_in_arc = 0;
    gmap<int, rcount> next_value;

    arc_range* r_it = overlap->ranges.begin();
    for(; r_it!=NULL ; r_it = overlap->ranges.next() ) {

        if (r_it->start <= i && r_it->end >= i) {

             #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Test Arc 1 " +  std::to_string(r_it->start) + " " + std::to_string(r_it->end) + " \n");
            #endif 
            // start area
            if (i == r_it->start ) {
                
                count_raw_node* rcl = &node_counts[gnodes[i]];
                rcl->add_node_start(node_it->exons); 
                
                if (i == j) {
                    // just this node left
                    break;
                }
                
                count_raw_edge* rec = &edge_counts[r_it->arc];
                for (gmap<int, exon_group_count>::iterator egi = node_it->exons->count_series.begin(); egi != node_it->exons->count_series.end(); ++egi) {
                    next_value[egi->first] = egi->second.total_lefts + egi->second.hole_end_counts[0] - egi->second.hole_start_counts[0];
                    if (!rec->series[egi->first].initialized) {
                        rec->series[egi->first].starts.assign(rec->size, {});
                        rec->series[egi->first].ends.assign(rec->size, {});
                        rec->series[egi->first].splits.assign(rec->size + 1, 0);
                        rec->series[egi->first].initialized = true;
                    }
                    rec->series[egi->first].splits[0] += next_value[egi->first];
                    
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("iS+ " +  std::to_string(egi->first) + " " + std::to_string(next_value[egi->first]) + " to " + std::to_string(rec->series[egi->first].splits[0]) + " \n");
                    #endif 
                }
            } else if (i == r_it->end) {
                
                count_raw_node* rcl = &node_counts[gnodes[i]];
                rcl->add_node_start(node_it->exons); // does both counts in case i == j
                
                if (i == j) {
                    // just this node left
                    break;
                }
                // otherwise there IS another range!
                r_it = overlap->ranges.next();
                // so we take this a the current actual first range, same node applies though
                
                logger::Instance()->debug("Switch Arc " +  std::to_string(r_it->start) + " " + std::to_string(r_it->end) + " \n");
                
                count_raw_edge* rec = &edge_counts[r_it->arc];
                for (gmap<int, exon_group_count>::iterator egi = node_it->exons->count_series.begin(); egi != node_it->exons->count_series.end(); ++egi) {
                    next_value[egi->first] = egi->second.total_lefts + egi->second.hole_end_counts[0] - egi->second.hole_start_counts[0];
                    if (!rec->series[egi->first].initialized) {
                        rec->series[egi->first].starts.assign(rec->size, {});
                        rec->series[egi->first].ends.assign(rec->size, {});
                        rec->series[egi->first].splits.assign(rec->size + 1, 0);
                        rec->series[egi->first].initialized = true;
                    }
                    rec->series[egi->first].splits[0] += next_value[egi->first];
                }
            } else  {
                
                unsigned int c = find_index_global_to_sub(i, &ai[r_it->arc].edge_specifier);
                // the start hits in the middle of existing edge
                // so we add in there
                
                logger::Instance()->debug("Mid Start " +  std::to_string(edge_counts[r_it->arc].size) + " " + std::to_string(c)  + " \n");
                
                unsigned int index = edge_counts[r_it->arc].size + 1 - c;
                
                count_raw_edge* rec = &edge_counts[r_it->arc];
                rec->add_sub_counts_start(node_it->exons, index, next_value);
                
                offset_in_arc = index+1;
            }
            
            // middle and end!
            
            if (j < r_it->end && i != j) {
                
  
                logger::Instance()->debug("Special CASE \n");
                
                ++offset_new_exons;
                
                unsigned int c1 = find_index_global_to_sub(i, &ai[r_it->arc].edge_specifier);
                unsigned int c2 = find_index_global_to_sub(j, &ai[r_it->arc].edge_specifier);
                
                unsigned int l = edge_counts[r_it->arc].size + 1 - c2;
                
                count_raw_edge* rec = &edge_counts[r_it->arc];
                rec->add_sub_counts_eg_range(node_it->exons, offset_new_exons, offset_new_exons + c1 - c2 - 2 , offset_in_arc + 1, next_value);
                rec->add_sub_counts_end(node_it->exons, l, next_value);
                
            } else if (i != j && i != r_it->end) {  // j >= r_it_end 
                
                unsigned int c = find_index_global_to_sub(i, &ai[r_it->arc].edge_specifier);
                 
                --c; // start gone
                ++offset_new_exons;

                if (j >= r_it->end) {
                    --c;
                }

                logger::Instance()->debug("Special CASE 2\n");
                
                // we used up the first one
                count_raw_edge* rec = &edge_counts[r_it->arc];
                rec->add_sub_counts_eg_range(node_it->exons, offset_new_exons, offset_new_exons + c - 1, edge_counts[r_it->arc].size + 1 - c , next_value);
                offset_new_exons += c; 

                if (j == r_it->end) {
                    count_raw_node* rcl = &node_counts[gnodes[j]];
                    rcl->add_node_end(node_it->exons); 
                } else { // if (j > r_it->end) {
                    count_raw_node* rcl = &node_counts[gnodes[r_it->end]];
                    rcl->add_node_initial_index(node_it->exons, offset_new_exons, next_value);
                }
                ++offset_new_exons;

            }
            
            break;
        }
    }

    if (r_it==NULL) {
        return;
    }

    if (r_it->end < j) {
        r_it = overlap->ranges.next();
        
         logger::Instance()->debug("Next Value c \n");
        
        for(; r_it!=NULL ; r_it = overlap->ranges.next() ) {

            logger::Instance()->debug("Test Arc 2 " +  std::to_string(r_it->start) + " " + std::to_string(r_it->end) + " \n");

            if (j > r_it->start && j < r_it->end) {
                
                logger::Instance()->debug("End condition  \n");
                
                unsigned int c = find_index_global_to_sub(j, &ai[r_it->arc].edge_specifier);
                unsigned int l = edge_counts[r_it->arc].size + 1 - c;
                
                count_raw_edge* rec = &edge_counts[r_it->arc];
                for (gmap<int, exon_group_count>::iterator egi = node_it->exons->count_series.begin(); egi != node_it->exons->count_series.end(); ++egi) {
                    if (!rec->series[egi->first].initialized) {
                        rec->series[egi->first].starts.assign(rec->size, {});
                        rec->series[egi->first].ends.assign(rec->size, {});
                        rec->series[egi->first].splits.assign(rec->size + 1, 0);
                        rec->series[egi->first].initialized = true;
                    }
                    rec->series[egi->first].splits[0] += next_value[egi->first];
                }
                rec->add_sub_counts_eg_range(node_it->exons, offset_new_exons, offset_new_exons + l - 1, 1 , next_value);
                rec->add_sub_counts_end(node_it->exons, l, next_value);
                
                offset_new_exons += l + 1;
                
                break;
            } else {
                
                 logger::Instance()->debug("Fully covered  \n");
                
                // so we have a full arc to add!
                count_raw_edge* rec = &edge_counts[r_it->arc];
                for (gmap<int, exon_group_count>::iterator egi = node_it->exons->count_series.begin(); egi != node_it->exons->count_series.end(); ++egi) {
                    if (!rec->series[egi->first].initialized) {
                        rec->series[egi->first].starts.assign(rec->size, {});
                        rec->series[egi->first].ends.assign(rec->size, {});
                        rec->series[egi->first].splits.assign(rec->size + 1, 0);
                        rec->series[egi->first].initialized = true;
                    }
                    rec->series[egi->first].splits[0] += next_value[egi->first];
                }
                rec->add_sub_counts_eg_range(node_it->exons, offset_new_exons, offset_new_exons + edge_counts[r_it->arc].size-1, 1 , next_value);
                offset_new_exons += edge_counts[r_it->arc].size + 1;
                
                if (j == r_it->end) {
                    count_raw_node* rcl = &node_counts[gnodes[j]];
                    rcl->add_node_end(node_it->exons); 
                    break;
                } else if (j > r_it->end) {
                    count_raw_node* rcl = &node_counts[gnodes[r_it->end]];
                    rcl->add_node_initial_index(node_it->exons, offset_new_exons-1, next_value);
                } 
            }
            
        }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("------------------- \n");
    #endif  
    
}

unsigned int base_manager::find_index_global_to_sub(unsigned int start, exon_edge* edge) {
    unsigned int c = 0;
    boost::dynamic_bitset<>::size_type index = start;
    while(index != boost::dynamic_bitset<>::npos) {
        ++c;
        index = edge->id.find_next(index);
    }
    return c;
}

void base_manager::add_contained_start(range_helper &ranges,
    unsigned int i,
    ListDigraph::NodeMap<count_raw_node> &node_counts,
    ListDigraph::ArcMap<arc_range*> &edge_range,
    ListDigraph::ArcMap<count_raw_edge> &edge_counts,
    std::deque<arc_range> &range_list, bool add_source) {
    // just one object
    
        ListDigraph::Node start;
        
        arc_range* r_it = ranges.begin();
        for(; r_it!=NULL ; r_it = ranges.next() ) {

            if (r_it->start <= i && r_it->end >= i) {
                
                if (r_it->start == i) {
                    // we have an existing node
                    start = g.source(r_it->arc);
                } else if (i ==  r_it->end) {
                    start = g.target(r_it->arc);
                } else {
                    // we need to add a new node
                    split_edge_without_compacting(i, ranges, r_it, edge_range, edge_counts, node_counts, NULL, range_list, false);
                    start = g.source(r_it->arc);
                }
                
                break;
            }
        }
        
        ListDigraph::InArcIt a(g, start);
        for (; a != INVALID; ++a) {
            if ( g.source(a) == s ) {
                break;
            } 
        }
        
        ListDigraph::Arc arc;
        if (a == INVALID && add_source) {
            // does not exist, just add and have fun
            arc = g.addArc(s, start);
            initialize_source_drain_arc(arc);
        }  
}


void base_manager::snap_contained(range_helper &ranges,
    unsigned int i,
    ListDigraph::NodeMap<count_raw_node> &node_counts,
    ListDigraph::ArcMap<arc_range*> &edge_range,
    ListDigraph::ArcMap<count_raw_edge> &edge_counts,
    std::deque<arc_range> &range_list,
    ListDigraph::Node* gnodes, bool* gnode_set,
    gmap<exon_edge, ListDigraph::Arc> &junction_to_arc) {
    // just one object
    
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("snap contained " + std::to_string(i) + " " + std::to_string(gnode_set[i]) + "\n");
    #endif
    
    arc_range* r_it = ranges.begin();
    for(; r_it!=NULL ; r_it = ranges.next() ) {

        if (r_it->start <= i && r_it->end >= i) {
            // we need to add a new node;
            if (r_it->start != i && i != r_it->end) {
                split_edge(i, ranges, r_it, edge_range, edge_counts, node_counts, gnode_set, gnodes, junction_to_arc, range_list);      
            }

            break;
        }
    }
}

void base_manager::add_contained_end(range_helper &ranges,
    unsigned int i,
    ListDigraph::NodeMap<count_raw_node> &node_counts,
    ListDigraph::ArcMap<arc_range*> &edge_range,
    ListDigraph::ArcMap<count_raw_edge> &edge_counts,
    std::deque<arc_range> &range_list, bool add_drain) {
    // just one object
    
        ListDigraph::Node start;   
                
        arc_range* r_it = ranges.begin();
        for(; r_it!=NULL ; r_it = ranges.next() ) {
            
            if (r_it->start <= i && r_it->end >= i) {
                
                if (r_it->start == i) {
                    // we have an existing node
                    start = g.source(r_it->arc);
                } else if (i ==  r_it->end) {
                    start = g.target(r_it->arc);
                } else {
                    // we need to add a new node
                    split_edge_without_compacting(i, ranges, r_it, edge_range, edge_counts, node_counts, NULL, range_list, false);
                    start = g.source(r_it->arc);
                }
                
                break;
            }
        }
                
        ListDigraph::OutArcIt a(g, start);
        for (; a != INVALID; ++a) {
            if ( g.target(a) == t ) {
                break;
            } 
        }
        
        ListDigraph::Arc arc;
        if (a == INVALID && add_drain) {
            // does not exist, just add and have fun
            arc = g.addArc(start, t);
            initialize_source_drain_arc(arc);
        }  
}

void base_manager::split_edge(unsigned int i,
        range_helper &ranges, // backlinked node that overlaps
        arc_range*& r_it, // range of overlap

         ListDigraph::ArcMap<arc_range*> &edge_range,
         ListDigraph::ArcMap<count_raw_edge> &edge_counts,
         ListDigraph::NodeMap<count_raw_node> &node_counts,
        
        bool* gnode_set, ListDigraph::Node* gnodes, 
         gmap<exon_edge, ListDigraph::Arc> &junction_to_arc,
         std::deque<arc_range> &range_list) {
    
    // this arc needs to be split up
    ListDigraph::Arc arc_save = r_it->arc;
    exon_edge edge = ai[r_it->arc].edge_specifier; 
    exon_edge left;
    edge.left_split(i, left);
    exon_edge right;
    edge.right_split(i, right);

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Split Norm " + edge.to_string() + " at pos " + std::to_string(i) + " to " + left.to_string() + " and " + right.to_string() + ".\n");
    #endif
    
    // does the current node exist?
    if (!gnode_set[i]) {
        gnodes[i] = g.addNode();
        gnode_set[i] = true;
    }

    ListDigraph::Arc left_a, right_a;
    arc_range* rl, *rr;
    
    // test if left split arc already exist!
    typename gmap<exon_edge, ListDigraph::Arc>::iterator lm = junction_to_arc.find(left);
    if (lm == junction_to_arc.end()) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("New Left.\n");
        #endif
        
        // so the arc really does not exist already, so make it
        left_a = g.addArc(gnodes[r_it->start], gnodes[i]);
        range_list.push_back(arc_range(r_it->start, i, left_a));
        rl = &range_list.back();
        
        lm = junction_to_arc.insert(std::make_pair(left, left_a)).first;
        ai[left_a].edge_specifier = left;
        edge_range[left_a] = rl;
        
    } else {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Found Left.\n");
        #endif
        
        left_a = lm->second;
        rl = edge_range[left_a];
    }
    // test if right split arc already exist!
    typename gmap<exon_edge, ListDigraph::Arc>::iterator lr = junction_to_arc.find(right);
    if (lr == junction_to_arc.end()) {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("New Right.\n");
        #endif
        
        // so the arc really does not exist already, so make it
        right_a = g.addArc(gnodes[i], gnodes[r_it->end]);
        range_list.push_back(arc_range(i, r_it->end, right_a));
        rr = &range_list.back();
               
        lr = junction_to_arc.insert(std::make_pair(right, right_a)).first;
        ai[right_a].edge_specifier = right;
        edge_range[right_a] = rr;
        
    } else {
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Found Right.\n");
        #endif
        
        right_a = lr->second;
        rr = edge_range[right_a];
    }

    // update counts to newly created stuff
    unsigned int lc = left.id.count() - 2;
    edge_counts[left_a].add_sub_counts_re(edge_counts[arc_save], 0, lc-1); 
    edge_counts[right_a].add_sub_counts_re(edge_counts[arc_save], lc+1, lc + right.id.count() - 2);
     
    node_counts[gnodes[i]].add_node_index(edge_counts[arc_save], lc);
    
    // now split up the range        
    r_it = ranges.split(i, rl, rr);
            
    // now we can safely remove all wrong nodes
    g.erase(arc_save); // this also clean out the map
    junction_to_arc.erase(edge);
     
}

void base_manager::split_edge_without_compacting(unsigned int i,
        range_helper &ranges, // backlinked node that overlaps
        arc_range*& r_it, // range of overlap
         ListDigraph::ArcMap<arc_range*> &edge_range,
         ListDigraph::ArcMap<count_raw_edge> &edge_counts,
         ListDigraph::NodeMap<count_raw_node> &node_counts,
         ListDigraph::Node* gnodes,
         std::deque<arc_range> &range_list, bool snap_node) {
    
    // this arc needs to be split up
    ListDigraph::Arc arc_save = r_it->arc;
    exon_edge edge = ai[r_it->arc].edge_specifier;
    exon_edge left;
    edge.left_split(i, left);
    exon_edge right;
    edge.right_split(i, right);

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Split WOC " + edge.to_string() + " at pos " + std::to_string(i) + " to " + left.to_string() + " and " + right.to_string() + ".\n");
    #endif
    
    ListDigraph::Arc left_a, right_a;
    arc_range* rl, *rr;
    
    ListDigraph::Node node;
    
    if (!snap_node) {
        node = g.addNode();
        node_index[node] = i;
    } else {
        node = gnodes[i];
    }
    
    // so the arc really does not exist already, so make it
    left_a = g.addArc(g.source(arc_save), node);
    range_list.push_back(arc_range(r_it->start, i, left_a));
    rl = &range_list.back();

    ai[left_a].edge_specifier = left;
    edge_range[left_a] = rl;

   
    // so the arc really does not exist already, so make it
    right_a = g.addArc(node, g.target(arc_save));
    range_list.push_back(arc_range(i, r_it->end, right_a));
    rr = &range_list.back();

    ai[right_a].edge_specifier = right;
    edge_range[right_a] = rr;
        
    
    // update counts to newly created stuff
    unsigned int lc = left.id.count() - 2;
    edge_counts[left_a].add_sub_counts_re(edge_counts[arc_save], 0, lc-1); 
    edge_counts[right_a].add_sub_counts_re(edge_counts[arc_save], lc+1, lc + right.id.count() - 2);
    
    node_counts[node].add_node_index(edge_counts[arc_save], lc);
        
    // now split up the range        
    r_it = ranges.split(i, rl, rr);
            
    // now we can safely remove all wrong nodes
    g.erase(arc_save); // this also clean out the map
     
}


void base_manager::create_region(unsigned int i, rcount basevalue, std::map< rpos,rcount > &lefts, std::map< rpos,rcount > &rights, region &r) {

    int rstate = 0; // 0 undefined, 1 up, 2 down
    rpos start = meta->exons[i].left;
    rpos end = meta->exons[i].right;
    
    rpos pos = start; 
    rcount value = basevalue;
    
    rcount pre_value;
    rcount start_value = value;
    rcount basecount = 0;
    rpos pre_pos;
    
    std::map< rpos,rcount >::const_iterator li = lefts.begin();
    std::map< rpos,rcount >::const_iterator ri = rights.begin();
    
    if (li != lefts.end() && li->first == pos) { // take up reads starting directly exon start region
        start_value += li->second; 
        value += li->second;
        ++li;
    }
    
    while (li != lefts.end() || ( ri != rights.end() && ri->first + 1 <= end )) {
        
//        logger::Instance()->debug("Value: " + std::to_string(value) + "\n");
//        logger::Instance()->debug("STEP: " + std::to_string(li->first) + "-" +std::to_string(li->second) + " " + std::to_string(ri->first) + "-" +std::to_string(ri->second) + "\n");
//        
        pre_value = value;
        pre_pos = pos;
                
        if (li == lefts.end() || (ri != rights.end() && ri->first + 1 < li->first ) ) {            
            value -= ri->second;
            pos = ri->first + 1;
            ri++;
        } else if (ri == rights.end() || ri->first + 1 > li->first ) {
            value += li->second;
            pos = li->first;
            li++;
        } else {
            value -= ri->second;
            value += li->second;
            pos = li->first;
            ri++;
            li++;
        }

        switch(rstate) {
            case 0:
//              #ifdef ALLOW_DEBUG
//                  logger::Instance()->debug("State 0 "+ std::to_string(pre_value) + " " + std::to_string(value) +"\n");
//              #endif
                if (pre_value < value) {
                    rstate = 1;
                } else {
                    rstate = 2;
                }
                break;
            case 1:
//               #ifdef ALLOW_DEBUG
//                  logger::Instance()->debug("State 1 "+ std::to_string(pre_value) + " " + std::to_string(value) +"\n");
//              #endif
                if (pre_value > value) {
                    rstate = 2;
                    basecount += pre_value;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Region 1 " + std::to_string(start) + " " + std::to_string(pre_pos) + " " + std::to_string(start_value) + " " + std::to_string(pre_value) + " " + std::to_string(basecount)  + "\n");
                    #endif
                    r.subregions.push_back( region::subregion(start, pre_pos, start_value, pre_value, basecount));
                    start_value = pre_value;
                    start = pre_pos+1;
                    ++pre_pos;
                    basecount = 0;
                } else {
                    if (pos - start > 50) {
                        // split to keep accurate borders
                        basecount += pre_value;
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Region 1S " + std::to_string(start) + " " + std::to_string(pre_pos) + " " + std::to_string(start_value) + " " + std::to_string(pre_value) + " " + std::to_string(basecount)  + "\n");
                        #endif
                        r.subregions.push_back( region::subregion(start, pre_pos, start_value, pre_value, basecount));
                        start_value = pre_value;
                        start = pre_pos+1;
                        ++pre_pos;
                        basecount = 0;
                    }
                }
                break;
            case 2:
//              #ifdef ALLOW_DEBUG
//                  logger::Instance()->debug("State 2 "+ std::to_string(pre_value) + " " + std::to_string(value) +"\n");
//              #endif
                if (pre_value < value) {
                    rstate = 1;
                    basecount += pre_value;
                    #ifdef ALLOW_DEBUG
                    logger::Instance()->debug("Region 2 " + std::to_string(start) + " " + std::to_string(pre_pos) + " " + std::to_string(start_value) + " " + std::to_string(pre_value) + " " + std::to_string(basecount) + "\n");
                    #endif
                    r.subregions.push_back( region::subregion(start, pre_pos, start_value, pre_value, basecount));
                    start_value = pre_value;
                    start = pre_pos +1;
                    ++pre_pos;
                    basecount = 0;
                    
                    if (pre_value == 0 && start + 1 < pos) {
                        // insert a zero region!
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Region 0 " + std::to_string(start) + " " + std::to_string(pos-2) + "\n");
                        #endif
                        r.subregions.push_back( region::subregion(start, pos-2, 0, 0, 0));
                        start = pos-1;
                    }
                    
                } else {
                    if (pos - start > 50 ) {
                        // split to keep accurate borders
                        basecount += pre_value;
                        #ifdef ALLOW_DEBUG
                        logger::Instance()->debug("Region 2S " + std::to_string(start) + " " + std::to_string(pre_pos) + " " + std::to_string(start_value) + " " + std::to_string(pre_value) + " " + std::to_string(basecount)  + "\n");
                        #endif
                        r.subregions.push_back( region::subregion(start, pre_pos, start_value, pre_value, basecount));
                        start_value = pre_value;
                        start = pre_pos+1;
                        ++pre_pos;
                        basecount = 0;
                    }
                }
                break;
            default:
                break;
        }
        basecount += pre_value * (pos - pre_pos);
    }
    basecount += value * (end - pos + 1);
    r.subregions.push_back( region::subregion(start, end, start_value, value, basecount));
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Region " + std::to_string(rstate) + " " + std::to_string(start) + " " + std::to_string(end) + " " + std::to_string(start_value) + " " + std::to_string(value) + " " + std::to_string(basecount) + "\n");
    #endif
}

void base_manager::create_region_from_node(unsigned int i, count_raw_node &node_counts, flow_series &fs) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Create Node Region " + std::to_string(i) + "\n");
    #endif
    for ( std::set<int>::iterator iii = input_ids.begin(); iii !=  input_ids.end(); ++iii) {

      gmap<int, count_raw_node::series_struct>::iterator ssi = node_counts.series.find(*iii);
      region r;
      r.total_length = meta->exons[i].right - meta->exons[i].left + 1;

      if(ssi != node_counts.series.end()) {
    
		create_region(i, ssi->second.total_rights, ssi->second.lefts, ssi->second.rights, r);
		
		fs.series[ssi->first].capacity = std::max(r.get_max(), 1ul);
		fs.series[ssi->first].mean = capacity_mean(std::max(r.get_average(), 0.1f), r.total_length);
		
		fs.series[ssi->first].average_to_first_zero_from_left = r.get_average_to_first_zero_from_left();
		fs.series[ssi->first].average_to_first_zero_from_right = r.get_average_to_first_zero_from_right();
		fs.series[ssi->first].deviation = r.get_deviation();
		fs.series[ssi->first].left = r.get_left();
		fs.series[ssi->first].right = r.get_right();
		fs.series[ssi->first].length = r.total_length;
		fs.series[ssi->first].length_to_first_zero_left = r.get_length_to_first_zero_from_left();
		fs.series[ssi->first].length_to_first_zero_right = r.get_length_to_first_zero_from_right();
		fs.series[ssi->first].region_average = r.get_average();
		fs.series[ssi->first].region_max = r.get_max();
		fs.series[ssi->first].region_min = r.get_min();
		fs.series[ssi->first].right = r.get_right();   

       } else {
                fs.series[*iii].mean = capacity_mean(0, r.total_length); // leave rest unitialized
       }
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("------- \n");
    #endif
}

                
void base_manager::create_region_from_edge(exon_edge& edge, count_raw_edge &edge_counts, flow_series &fs) {
    
    for(gmap<int, count_raw_edge::series_struct>::iterator ssi = edge_counts.series.begin(); ssi != edge_counts.series.end(); ++ssi) {
    
        region r;
        unsigned int c = 0; // index in count
        unsigned int index = edge.id.find_first();
        index = edge.id.find_next(index);
        r.total_length = 0;

        r.subregions.push_back( region::subregion(meta->exons[index].left, meta->exons[index].left, ssi->second.splits[0], ssi->second.splits[0], ssi->second.splits[0]));
        r.total_length += 1;
        
        for(; c < edge_counts.size; ++c) {

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("Create Edge Region " + std::to_string(index) + "\n");
            #endif

            create_region(index, ssi->second.splits[c], ssi->second.starts[c].ref(), ssi->second.ends[c].ref(), r);

            r.total_length += meta->exons[index].right - meta->exons[index].left + 2;
            r.subregions.push_back( region::subregion(meta->exons[index].right, meta->exons[index].right, ssi->second.splits[c+1], ssi->second.splits[c+1], ssi->second.splits[c+1]));

            index = edge.id.find_next(index);

            #ifdef ALLOW_DEBUG
            logger::Instance()->debug("------- \n");
            #endif
        }

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Region Split \n");
        #endif

        fs.series[ssi->first].capacity = r.get_max();
        fs.series[ssi->first].mean = capacity_mean(std::max(r.get_average(), 0.1f), r.total_length);
        
        fs.series[ssi->first].average_to_first_zero_from_left = r.get_average_to_first_zero_from_left();
        fs.series[ssi->first].average_to_first_zero_from_right = r.get_average_to_first_zero_from_right();
        fs.series[ssi->first].deviation = r.get_deviation();
        fs.series[ssi->first].left = r.get_left();
        fs.series[ssi->first].right = r.get_right();
        fs.series[ssi->first].length = r.total_length;
        fs.series[ssi->first].length_to_first_zero_left = r.get_length_to_first_zero_from_left();
        fs.series[ssi->first].length_to_first_zero_right = r.get_length_to_first_zero_from_right();
        fs.series[ssi->first].region_average = r.get_average();
        fs.series[ssi->first].region_max = r.get_max();
        fs.series[ssi->first].region_min = r.get_min();
        fs.series[ssi->first].right = r.get_right();
    }

}


// ##################### Print #####################
void base_manager::print_raw( std::deque<overlap_node> &nodes,  std::deque<contained_node> &contained,  std::deque<contained_node> &filtered_contained) {
    
    logger::Instance()->debug("====================================================\n");
    logger::Instance()->debug("Nodes\n");
    unsigned int i = 0;
    for (std::deque<overlap_node>::iterator it = nodes.begin(); it!= nodes.end(); ++it, ++i) {
         if (!it->activated) {
             continue;
         }
         logger::Instance()->debug("Node "+std::to_string(i)+ " " + it->exons->bin_mask.to_string() + "\n");
         logger::Instance()->debug("Contains ");
         for (graph_list<contained_node*>::iterator c = it->contains.begin(); c!= it->contains.end(); ++c) {
             
             // find to print only those that are no transitive
             if (std::find( (*c)->contained_in.begin(), (*c)->contained_in.end(), &*it ) != (*c)->contained_in.end() ) {
                 logger::Instance()->debug((*c)->exons->bin_mask.to_string()+",");
             } 
             
         }
         logger::Instance()->debug("\n");
         logger::Instance()->debug("Links ");
         for(graph_list<overlap_node*>::iterator l = it->links.begin(); l != it->links.end(); ++l) {
             logger::Instance()->debug( std::to_string( std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), **l ))) + ", " );
         }
         logger::Instance()->debug("\n");
         logger::Instance()->debug("BackLinks ");
         for(graph_list<overlap_node*>::iterator l = it->back_links.begin(); l != it->back_links.end(); ++l) {
             logger::Instance()->debug( std::to_string( std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), **l ))) + ", " );
         }
         logger::Instance()->debug("\n");
    }
    logger::Instance()->debug("====================================================\n");
}

void base_manager::print_graph_debug(std::ostream &os, ListDigraph::NodeMap<count_raw_node> &node_counts, ListDigraph::ArcMap<count_raw_edge> &edge_counts) {
    
    if (meta->size != 1) {
        
        ListDigraph::NodeMap<capacity_type> nc(g);
        ListDigraph::ArcMap<capacity_type> ec(g);
        
        // there is a graph!
        digraphWriter(g, os)
            .arcMap("ai", ai)
            .arcMap("edge_coverage", edge_counts)
            .nodeMap("node_coverage", node_counts)
            .node("source", s)
            .node("drain", t)
            .run();  
    
    }
    
}

bool base_manager::has_single_exons() {
    return !single_exons.empty();
}

void base_manager::add_single_exons() {
    
    for (std::set<single_exon >::iterator it = single_exons.begin(); it != single_exons.end(); ++it) {
        float combined = 0;
        for (std::set<int>::iterator iii = input_ids.begin(); iii != input_ids.end(); ++iii) {
            int id = *iii;
            
            if (input_ids.size() > 1 && !options::Instance()->is_compute_all_singles() && id != -1) {
                    continue;
            }
            
            if (id == -1) {
                continue;
            }
            
            if (it->guide || (it->capacity.find(id) != it->capacity.end() && it->capacity.at(id) >= 10)) {
            
                transcripts[id].transcripts.push_back(lazy<transcript>());

                transcripts[id].transcripts.back()->found_edge =  exon_edge(meta->size);
                transcripts[id].transcripts.back()->found_edge.set(it->meta, true);
                for (std::map<int, float>::const_iterator ci = it->capacity.begin(); ci != it->capacity.end(); ++ci) {
		    //capacity_type cap = ci->second;
                    transcripts[id].transcripts.back()->series[ci->first].flow = ci->second;
                    transcripts[id].transcripts.back()->series[ci->first].mean = ci->second;
                    transcripts[id].transcripts.back()->series[ci->first].score = ci->second;
                }
                transcripts[id].transcripts.back()->guided = it->guide;
                
                if (input_ids.size() > 1) { // we have a join!
                    combined += it->capacity.at(id);
                }  
            }
        }
        if (combined > 0) {
            // add to joined transcripts
            transcripts[-1].transcripts.push_back(lazy<transcript>());

            transcripts[-1].transcripts.back()->found_edge =  exon_edge(meta->size);
            transcripts[-1].transcripts.back()->found_edge.set(it->meta, true);
            for (std::map<int, float>::const_iterator ci = it->capacity.begin(); ci != it->capacity.end(); ++ci) {
                //capacity_type cap = ci->second;
                transcripts[-1].transcripts.back()->series[ci->first].flow = ci->second;
                transcripts[-1].transcripts.back()->series[ci->first].mean = ci->second;
                transcripts[-1].transcripts.back()->series[ci->first].score = ci->second;
            }
            transcripts[-1].transcripts.back()->series[-1].flow = combined;
            transcripts[-1].transcripts.back()->series[-1].mean = combined;
            transcripts[-1].transcripts.back()->series[-1].score = combined;
                
            transcripts[-1].transcripts.back()->guided = it->guide;
        }
    }
}
