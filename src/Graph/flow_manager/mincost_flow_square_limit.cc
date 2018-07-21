/* 
 * File:   mincost_flow_square_limit.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on January 9, 2017, 2:39 PM
 */

#include "mincost_flow_square_limit.h"

mincost_flow_square_limit::mincost_flow_square_limit(pre_graph* raw, exon_meta* meta, const std::string &chromosome) : mincost_flow_square(raw, meta, chromosome) { 
}


mincost_flow_square_limit::~mincost_flow_square_limit() {
}


void mincost_flow_square_limit::add_offset_edge(capacity_type capacity, capacity_type orig_cap, int exon_count,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference) {
    
    if (orig_cap == 0) return;
    
    // we create the overflow edge first
    capacity_type percent_mark = orig_cap*1/options::Instance()->get_coverage_bias_limit();
    capacity_type overflow_capacity = capacity - percent_mark;
    
    ListDigraph::Arc na = og.addArc(sn, tn);
    upper[na] = overflow_capacity;
    cost[na] = overflow_cost;
    reference.push_back(na);
    
//    logger::Instance()->debug("Overflow = " + std::to_string(overflow_capacity) + " " + std::to_string(overflow_cost) + "\n");
    
    capacity_type i = step_interval;
    for (; i <= percent_mark; i+= step_interval) {
        
        ListDigraph::Arc na = og.addArc(sn, tn);
        upper[na] = step_interval;
        
        capacity_type c1 = i;
        capacity_type c2 = i-step_interval;
        
        cost[na] = (c1*c1 - c2*c2 ) / step_interval * int_scaling_factor / orig_cap * (exon_count-1);
        reference.push_back(na);
        
    }
    
    i = i - step_interval;
    capacity_type leftover = percent_mark - i;
    if (leftover > 0) {
        // add one last edge
        ListDigraph::Arc na = og.addArc(sn, tn);
        upper[na] = leftover;
        
        capacity_type c1 = percent_mark;
        capacity_type c2 = i;
        
        cost[na] = (c1*c1 - c2*c2 ) / leftover * int_scaling_factor / orig_cap * (exon_count-1);
        reference.push_back(na);
    } 
}

void mincost_flow_square_limit::init_variables(capacity_type max_supply) {
    
    step_interval = max_supply / options::Instance()->get_max_arc_split();
    if (step_interval == 0) {
        step_interval = 1;
    } 
    
    capacity_type max_edge = max_supply*1/options::Instance()->get_coverage_bias_limit();
    
    if (max_edge == 0) {
        max_edge = 1;
    }
    
    overflow_cost = (max_edge*max_edge - (max_edge-step_interval) * (max_edge-step_interval) ) / step_interval * 100 * 100000 / max_edge; // maximal regular cost in any arc times 100
    
    int_scaling_factor = pow(ceil(log10(max_supply)),10);
}

bool mincost_flow_square_limit::modify_repeat(
            ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            capacity_type &max_supply,
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap) {
   
    
    bool new_sources = false;
    for (ListDigraph::NodeIt node(g); node != INVALID; ++node) {

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Modify Repeat " + std::to_string(g.id(node)) + "\n");
        #endif
        
            if (node == s || node == t) {
                continue;
            }

            // we always expand nodes, hence we need to look at all special cases
            ListDigraph::InArcIt i(g, node);
            ListDigraph::OutArcIt o(g, node);
            // note that they are 
            if (edge_type[o] == edge_types::NODE) {
                // the extension is to the RIGHT, we can have a SOURCE connection here
                
                bool has_left_open = false;
                bool is_source_connected = false;
                for( ;i != INVALID; ++i) {
                    // the 20 percentile overflow arcs are on first position, if any of them are filled we have overfilling!
                    if ( !arc_ref_backward[i].empty() && (flowmap[arc_ref_forward[i].front()] > 0 || flowmap[arc_ref_backward[i].front()] > 0) ) { 
                        has_left_open = true;
                    }
                    if ( g.source(i) == s) {
                        is_source_connected = true;
                    } 
                }
                if (is_source_connected || has_left_open) {
                    continue;
                }
                
                // if we are here left is NOT open, we need to find if there is a path going out from here
                if (arc_ref_backward[o].empty() || ( !( flowmap[arc_ref_forward[o].front()] > 0 || flowmap[arc_ref_backward[o].front()] > 0 ) )) {
                    continue;
                }
                // so the node has leftover capacity
                // we only add if we have a full node edge node
                ListDigraph::OutArcIt o2(g, g.target(o));
                bool add_source = false;
                for( ;o2 != INVALID; ++o2) {
                    if ( !arc_ref_backward[o2].empty() && (flowmap[arc_ref_forward[o2].front()] > 0 || flowmap[arc_ref_backward[o2].front()] > 0) ) {
                        ListDigraph::OutArcIt o3(g, g.target(o2));
                        
                        if (o3 == INVALID) {
                            continue;
                        }
                        
                        if (!arc_ref_backward[o3].empty() && ( flowmap[arc_ref_forward[o3].front()] > 0 || flowmap[arc_ref_backward[o3].front()] > 0) ) {
                            add_source = true;
                        }
                        break;
                    }
                }
                
                if (add_source) {
                    // we now have to create the source arc
                    ListDigraph::Arc na = g.addArc(s, node);
                    capacity[na] = 0;
                    edge_type[na] = edge_types::HELPER;
                    
                    ListDigraph::Arc in_copy = og.addArc(node_ref[s], node_ref[node]);
                    upper[in_copy] = max_supply;
                    cost[in_copy] = 0;
                    arc_ref_forward[na].push_back(in_copy);
                    
                    new_sources = true;
                }
            } else if (edge_type[i] == edge_types::NODE) {
                // the extension is to the LEFT, we can have a DRAIN connection here
                
                bool has_right_open = false;
                bool is_drain_connected = false;
                for( ;o != INVALID; ++o) {
                    if ( !arc_ref_backward[o].empty() && ( flowmap[arc_ref_forward[o].front()] > 0 || flowmap[arc_ref_backward[o].front()] > 0 )) {
                        has_right_open = true;
                    }
                    if ( g.target(o) == t) {
                        is_drain_connected = true;
                    } 
                }
                if (is_drain_connected || has_right_open) {
                    continue;
                }
                
                // if we are here right is NOT open, we need to find if there is a path going out from here
                if ( arc_ref_backward[i].empty() || !( flowmap[arc_ref_forward[i].front()] > 0 || flowmap[arc_ref_backward[i].front()] > 0) ) {
                    continue;
                }
                // so the node has leftover capacity
                // we only add if we have a full node edge node
                ListDigraph::InArcIt i2(g, g.source(i));
                bool add_drain = false;
                for( ;i2 != INVALID; ++i2) {
                    if (!arc_ref_backward[i2].empty() && (flowmap[arc_ref_forward[i2].front()] > 0 || flowmap[arc_ref_backward[i2].front()] > 0 ) ) {
                        ListDigraph::InArcIt i3(g, g.source(i2));
                        
                        if (i3 == INVALID) {
                            continue;
                        }
                        
                        if (!arc_ref_backward[i3].empty() && (flowmap[arc_ref_forward[i3].front()] > 0 || flowmap[arc_ref_backward[i3].front()] > 0) ) {
                            add_drain = true;
                        }
                        break;
                    }
                }
                
                if (add_drain) {
                    
                    // we now have to create the drain arc
                    ListDigraph::Arc na = g.addArc(node, t);
                    capacity[na] = 0;
                    edge_type[na] = edge_types::HELPER;
                    
                    ListDigraph::Arc in_copy = og.addArc(node_ref[node], node_ref[t]);
                    upper[in_copy] = max_supply;
                    cost[in_copy] = 0;
                    arc_ref_forward[na].push_back(in_copy);
                    
                    new_sources = true;
                }
            }
    } 
    
    return new_sources;
}