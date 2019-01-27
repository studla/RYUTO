/* 
 * File:   mincost_flow_base.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on January 4, 2017, 4:35 PM
 */

#include "mincost_flow_base.h"
#include <math.h>

mincost_flow_base::mincost_flow_base(pre_graph* raw, exon_meta* meta, const std::string &chromosome, std::set<int> &ids) : base_manager(raw, meta, chromosome, ids) { 
}


mincost_flow_base::~mincost_flow_base() {
}

// ################## interface functions for graph creation ##################


// we want to expand nodes
bool const mincost_flow_base::expand_exon_nodes() {
    return true;
}

void mincost_flow_base::initialize_source_drain_arc(const ListDigraph::Arc &arc) {
   // nothing needed to do, we still keep this stub
}

void mincost_flow_base::create_node_capacities(ListDigraph::Arc &arc, flow_series &r) {
       
}

void mincost_flow_base::create_edge_capacities(ListDigraph::Arc& arc, flow_series &r) {
    // nothing needed to do, we still keep this stub
}


void mincost_flow_base::finalize_flow_graph(int id) {
    // nothing needed to do, we still keep this stub     
}

// ################## functions for flow creation ##################


void mincost_flow_base::compute_flow(int id) {
    
    // we try to find the maximum flow that needs to be pushed through all
    ListDigraph::NodeMap<unsigned_capacity_type> ex_flow(g);
    capacity_type max_supply = find_exogenous_flow(ex_flow, id);
    init_variables(max_supply);
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Max Supply " + std::to_string(max_supply) + ".\n");
    #endif

    // we first create the offset graph
    ListDigraph::NodeMap<ListDigraph::Node > node_ref(g);
    ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > arc_ref_forward(g);
    ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > arc_ref_backward(g); 
    
    ListDigraph og;
    ListDigraph::ArcMap<capacity_type> upper(og);
    ListDigraph::ArcMap<unsigned_capacity_type> cost(og); 
    ListDigraph::NodeMap<unsigned_capacity_type> supply(og);
    ListDigraph::ArcMap<unsigned_capacity_type> flowmap(og);
    
    create_initial_offset_graph(ex_flow, node_ref, arc_ref_forward, arc_ref_backward, max_supply, og, upper, cost, supply, id);
    
    bool cont = true;
    while (cont) {
                
        #ifdef ALLOW_DEBUG
        if (options::Instance()->is_debug()) {
            ListDigraph::ArcMap<std::string> afwd(g);
            ListDigraph::ArcMap<std::string> abwd(g);

            for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
                for (std::deque< ListDigraph::Arc>::iterator it = arc_ref_forward[a].begin(); it!=arc_ref_forward[a].end(); ++it) {
                    afwd[a] += std::to_string(g.id(*it)) + ", ";
                }
                for (std::deque< ListDigraph::Arc>::iterator it = arc_ref_backward[a].begin(); it!=arc_ref_backward[a].end(); ++it) {
                    abwd[a] += std::to_string(g.id(*it)) + ", ";
                }
            }

            logger::Instance()->debug("Enter Flow Loop.\n");
            digraphWriter(g, std::cout)
                    .arcMap("ai", ai)
                    .arcMap("fwd", afwd)
                    .arcMap("bwd", abwd)
                    .node("source", s)
                    .node("drain", t)
                    .run();

//            digraphWriter(og, std::cout)
//                    .arcMap("Cap", upper)
//                    .arcMap("Cost", cost)
//                    .arcMap("Flow", flowmap)
//                    .nodeMap("Supply", supply)
//                    .run();
        }
        #endif
        
        // we start each cycle with computing the flow in the offset graph
        push_mincost_flow(og, upper, cost, supply, flowmap);
        
//         #ifdef ALLOW_DEBUG
//        if (options::Instance()->is_debug()) {
//            logger::Instance()->debug("Flow Step.\n");
//            digraphWriter(og, std::cout)
//                    .arcMap("Cap", upper)
//                    .arcMap("Cost", cost)
//                    .arcMap("Flow", flowmap)
//                    .nodeMap("Supply", supply)
//                    .run();
//         }
//        #endif
        
        cont = modify_repeat(node_ref, arc_ref_forward, arc_ref_backward, max_supply, og, upper, cost, supply, flowmap);
    }
    
    #ifdef ALLOW_DEBUG
        if (options::Instance()->is_debug()) {
            logger::Instance()->debug("Flow Done.\n");
            digraphWriter(og, std::cout)
                    .arcMap("Cap", upper)
                    .arcMap("Cost", cost)
                    .arcMap("Flow", flowmap)
                    .nodeMap("Supply", supply)
                    .run();
         }
        #endif
    
    transform_flow_to_orig(node_ref, arc_ref_forward, arc_ref_backward, og, flowmap, id);

}

capacity_type mincost_flow_base::find_exogenous_flow(ListDigraph::NodeMap<unsigned_capacity_type> &ex_flow_map, int id) {
    
    capacity_type max_supply;
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        
        capacity_type in = 0;
        for (ListDigraph::InArcIt a(g, n); a!=INVALID; ++a) {
            in += fs[a].series[id].capacity;
        }
        
        capacity_type out = 0;
        for (ListDigraph::OutArcIt a(g, n); a!=INVALID; ++a) {
            out += fs[a].series[id].capacity;
        }
        
        unsigned_capacity_type ex_flow = in - out;
        ex_flow_map[n] = ex_flow;
        
        if (ex_flow > 0) {
            max_supply += ex_flow;
        }
    }
    return max_supply;
}

void mincost_flow_base::create_initial_offset_graph(
            ListDigraph::NodeMap<unsigned_capacity_type> &ex_flow,
            ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            capacity_type &max_supply,
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            int id) {
    
    // we need star supply nodes, create them first
    ListDigraph::Node sstar = og.addNode();
    ListDigraph::Node tstar = og.addNode();
    supply[sstar] = max_supply;
    supply[tstar] = -max_supply;
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("star " + std::to_string(og.id(sstar))+ " "+ std::to_string(og.id(tstar)) + "\n");
    #endif
    
    for (ListDigraph::NodeIt n(g); n != INVALID; ++n) {
        ListDigraph::Node nn = og.addNode();
        node_ref[n] = nn;
        
        #ifdef ALLOW_DEBUG
        if (n == s)  logger::Instance()->debug("s " + std::to_string(og.id(nn))+ "\n");
        if (n == t)  logger::Instance()->debug("t " + std::to_string(og.id(nn))+ "\n");
        logger::Instance()->debug("Node " + std::to_string(g.id(n)) + " " + std::to_string(og.id(nn))+ " " + std::to_string(ex_flow[n]) + "\n");
        #endif

        if (n != s && n!= t) {
            if (ex_flow[n] < 0) {   // IN < OUT
                
                ListDigraph::Arc star_arc = og.addArc(nn, tstar);
                upper[star_arc] = -ex_flow[n];
                cost[star_arc] = 0;
                
            } else if(ex_flow[n] > 0) {   // IN > OUT
                
                ListDigraph::Arc star_arc = og.addArc(sstar, nn);
                upper[star_arc] = ex_flow[n];
                cost[star_arc] = 0;
            }
        }
    }
    
    #ifdef ALLOW_DEBUG
//        logger::Instance()->debug("Mid Create.\n");
//        digraphWriter(og, std::cout)
//                .arcMap("Cap", upper)
//                .arcMap("Cost", cost)
//                .nodeMap("Supply", supply)
//                .node("source", s)
//                .node("drain", t)
//                .run();
        #endif
    
    // now we iterate over the original graph and create the arc according to offset definition
    
    // backcyling for end and start supply
    ListDigraph::Arc gts_arc = og.addArc(node_ref[t], node_ref[s]);
    upper[gts_arc] = max_supply;
    cost[gts_arc] = 0;
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        
        ListDigraph::Node sn = node_ref[g.source(a)];
        ListDigraph::Node tn = node_ref[g.target(a)];
        
        if (ai[a].edge_type == edge_types::HELPER) {
    
//            logger::Instance()->debug("Arc Helper " + std::to_string(og.id(sn)) + " " + std::to_string(og.id(tn)) + "\n");
            
            add_offset_helper(max_supply, sn, tn, og, upper, cost, arc_ref_forward[a]);
            
        } else {
            
            logger::Instance()->debug("Arc Norm " + std::to_string(g.id(a)) + ": " + std::to_string(og.id(sn)) + " " + std::to_string(og.id(tn)) + "\n");
            // this is the interesting case of the actual arcs
            // note that star connections have already been made
//            logger::Instance()->debug("fwd");
            add_offset_edge(max_supply, fs[a].series[id].capacity, ai[a].edge_specifier.id.count(), sn, tn, og, upper, cost, arc_ref_forward[a]);
//            logger::Instance()->debug("bwd");
            add_offset_edge(fs[a].series[id].capacity, fs[a].series[id].capacity, ai[a].edge_specifier.id.count(), tn, sn, og, upper, cost, arc_ref_backward[a]);
        } 
    }
}

void mincost_flow_base::add_offset_edge(capacity_type capacity, capacity_type orig_cap, int exon_count,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference) {
    
    // this is the base implementation, so we just add on linear edgecosts
    
    ListDigraph::Arc na = og.addArc(sn, tn);
    upper[na] = capacity;
    cost[na] = 1;
    reference.push_back(na);
    
}

void mincost_flow_base::add_offset_helper(capacity_type capacity,
        ListDigraph::Node &sn, ListDigraph::Node &tn,
        ListDigraph &og, ListDigraph::ArcMap<capacity_type> &upper, ListDigraph::ArcMap<unsigned_capacity_type> &cost,
        std::deque< ListDigraph::Arc> &reference) {
    
        // helper get single arc with INF cap and no cost
        ListDigraph::Arc na = og.addArc(sn, tn);
        upper[na] = capacity;
        cost[na] = 0;
        reference.push_back(na);
}

bool mincost_flow_base::modify_repeat(
            ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            capacity_type &max_supply,
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap) {
    
    return false;
}

void mincost_flow_base::push_mincost_flow(
            ListDigraph &og,
            ListDigraph::ArcMap<capacity_type> &upper,
            ListDigraph::ArcMap<unsigned_capacity_type> &cost,
            ListDigraph::NodeMap<unsigned_capacity_type> &supply,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap) {
    
    NetworkSimplex<ListDigraph, unsigned_capacity_type, unsigned_capacity_type> ns(og);
    ns.upperMap(upper).costMap(cost);
    ns.supplyMap(supply);
    ns.run(NetworkSimplex<ListDigraph, unsigned_capacity_type, unsigned_capacity_type>:: BLOCK_SEARCH);
    
//    CostScaling<ListDigraph, unsigned_capacity_type> ns(og);
//    ns.upperMap(upper).costMap(cost);
//    ns.supplyMap(supply);
//    ns.run();
    ns.flowMap(flowmap);
    
}

void mincost_flow_base::transform_flow_to_orig(ListDigraph::NodeMap<ListDigraph::Node > &node_ref,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_forward,
            ListDigraph::ArcMap<std::deque< ListDigraph::Arc> > &arc_ref_backward,
            ListDigraph &og,
            ListDigraph::ArcMap<unsigned_capacity_type> &flowmap, int id) {
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        fs[a].series[id].flow = fs[a].series[id].capacity;
        for (std::deque< ListDigraph::Arc>::iterator fi = arc_ref_forward[a].begin(); fi != arc_ref_forward[a].end(); ++fi) {
            fs[a].series[id].flow += flowmap[*fi];
        }
        for (std::deque< ListDigraph::Arc>::iterator fi = arc_ref_backward[a].begin(); fi != arc_ref_backward[a].end(); ++fi) {
            fs[a].series[id].flow -= flowmap[*fi];
        }
        
        // we are done for normal nodes
        // but we need to DAG-ify backlinks
        if (ai[a].edge_type == edge_types::BACKLINK) {
            
            capacity_type max_cap = fs[a].series[id].flow;
            ListDigraph::Node in = g.source(a);
            ListDigraph::Node out = g.target(a);
            
            unsigned int cycle_id = ai[a].cycle_id_in;
            g.erase(a);
            ListDigraph::Arc ns = g.addArc(s, out);
            fs[ns].series[id].flow = max_cap;
            ListDigraph::Arc nt = g.addArc(in, t);
            fs[nt].series[id].flow = max_cap;
            ai[ns].edge_type = edge_types::HELPER;
            ai[ns].cycle_id_in = cycle_id;
            ai[nt].edge_type = edge_types::HELPER;
            ai[nt].cycle_id_out = cycle_id;
            
        }
        
    }
    
}

void mincost_flow_base::init_variables(capacity_type max_supply) {
    // nothing to do
}


// ################## simple printing ##################

void mincost_flow_base::print_graph(std::ostream &os) {
   
}