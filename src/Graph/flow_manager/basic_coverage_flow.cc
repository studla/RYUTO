/* 
 * File:   basic_coverage_flow.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on April 28, 2016, 2:56 PM
 */

#include "basic_coverage_flow.h"

#include <math.h>

basic_coverage_flow::basic_coverage_flow(pre_graph* raw, exon_meta* meta, const std::string &chromosome) : base_manager(raw, meta, chromosome) { 
}


basic_coverage_flow::~basic_coverage_flow() {
}


// ################## interface functions for graph creation ##################


// we want to expand nodes
bool const basic_coverage_flow::expand_exon_nodes() {
    return true;
}

void basic_coverage_flow::initialize_source_drain_arc(const ListDigraph::Arc &arc) {
    drain_source_arcs.push_back(arc); // just add for now, we init later
}

void basic_coverage_flow::create_node_capacities(ListDigraph::Arc &arc, region &r) {
    // logger::Instance()->debug("Set_node Cov " + std::to_string(count.coverage) + " " + std::to_string(count.read_count) + "\n");
    capacity[arc] = r.get_average();
}

void basic_coverage_flow::create_edge_capacities(ListDigraph::Arc& arc, region &r) {
    capacity[arc] = r.get_average();
    
}


void basic_coverage_flow::finalize_flow_graph() {
    // the source and drain nodes have no capacity yet, set it!
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("finalize graph  \n");
    #endif
    
    for ( std::deque<ListDigraph::Arc>::iterator it = drain_source_arcs.begin(); it != drain_source_arcs.end(); ++it ) {
        
        if (!g.valid(*it)) {
            continue;
        }
        
        rcount cap = 0;
        if ( g.source(*it) == s) {  // edge from source
            
            // we need to test the source if this drain is in between
            ListDigraph::InArcIt i(g, g.target(*it));
            ++i;// at least one, the source itself
            if (i != INVALID) {
                // this is an evil in between sink...
                ListDigraph::OutArcIt o(g, g.target(*it));
                slow_add_sources.insert(std::make_pair( edge_specifier[o].node_index , g.id(*it)));
                cap = 0;
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Slow Source.  \n");
                #endif
            } else {
                for (ListDigraph::OutArcIt a(g, g.target(*it)); a != INVALID; ++a) {
                    cap += capacity[a];
                } 
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Normal Source.  \n");
                #endif
            }
        } else { // edge to drain
            
            // we need to test the source if this darin is in between
            ListDigraph::OutArcIt o(g, g.source(*it));
            ++o;// at least one, the drain itself
            if (o != INVALID) {
                // this is an evil in between sink...
                ListDigraph::InArcIt i(g, g.source(*it));
                slow_add_drains.insert(std::make_pair( edge_specifier[i].node_index , g.id(*it)));
                cap = 0;
                
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Slow Drain.  \n");
                #endif
                
            } else {
                for (ListDigraph::InArcIt a(g, g.source(*it)); a != INVALID; ++a) {
                    cap += capacity[a];
                } 
                #ifdef ALLOW_DEBUG
                logger::Instance()->debug("Normal Drain.  \n");
                #endif
            }
        }
        capacity[*it] = cap;
    } 
    
}

// ################## functions for flow creation ##################

void basic_coverage_flow::compute_flow() {
    
    // we iterate the computations until we have found all potential sources/sinks or limit reached
    
    //first now we create the working copy
    ListDigraph wc; // working copy for modification
    ListDigraph::ArcMap<capacity_type> fc(wc);

    ListDigraph::ArcMap<exon_edge> ces(wc); 
    ListDigraph::ArcMap<edge_types::edge_type> cet(wc);
    ListDigraph::ArcMap<edge_length> cel(wc);
    ListDigraph::NodeMap<unsigned int> cni(wc);
    ListDigraph::Node cs,ct;

    ListDigraph::ArcMap<capacity_type> cap(wc);

    ListDigraph::NodeMap<ListDigraph::Node> node_ref(wc);
    
    ListDigraph::ArcMap<ListDigraph::Arc> arc_ref(wc);
    ListDigraph::ArcMap<ListDigraph::Arc> arc_ref_cop(g);

    digraphCopy(g,wc).arcMap(flow, fc).arcMap(edge_lengths, cel)
            .arcMap(edge_specifier, ces).arcMap(edge_type, cet).arcMap(capacity, cap)
            .nodeMap(node_index, cni).node(s,cs).node(t,ct)
            .nodeCrossRef(node_ref)
            .arcCrossRef(arc_ref).arcRef(arc_ref_cop).run();
    
    // first we need to create the actual flow
    // we have to put flow on the original graph, as this gives the base values for resolving everything
        
    #ifdef ALLOW_DEBUG
    print_graph_debug_copy(std::cout, wc, cap, ces, cet, cel, cs, ct);
    #endif

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Push Linear Flow.\n");
    #endif
    push_graph_flow(wc, cs, ct, cap, fc);   // inner sources and drains are set to 0 by finalize!
    transfer(wc, cet, cs, ct, cap, fc, arc_ref);
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Push Additional Sources.\n");  
    print_graph_debug_copy(std::cout, wc, fc, ces, cet, cel, cs, ct);    
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    for (std::set<std::pair<unsigned int, int> >::iterator it = slow_add_sources.begin(); it !=  slow_add_sources.end(); ++it) {
        
        ListDigraph::Arc arc = arc_ref_cop[wc.arcFromId(it->second)];
        for (ListDigraph::OutArcIt a(wc, wc.target(arc)); a != INVALID; ++a) {
            cap[arc] += cap[a];
        } 
        
        push_graph_flow(wc, cs, ct, cap, fc);
        transfer(wc, cet, cs, ct, cap, fc, arc_ref);
    }
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Push Additional Sinks.\n");
    print_graph_debug_copy(std::cout, wc, fc, ces, cet, cel, cs, ct);
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif

    for (std::set<std::pair<unsigned int, int> >::reverse_iterator it = slow_add_drains.rbegin(); it !=  slow_add_drains.rend(); ++it) {
        
        ListDigraph::Arc arc = arc_ref_cop[wc.arcFromId(it->second)];
        for (ListDigraph::InArcIt a(wc, wc.source(arc)); a != INVALID; ++a) {
            cap[arc] += cap[a];
        } 
        
        push_graph_flow(wc, cs, ct, cap, fc);
        transfer(wc, cet, cs, ct, cap, fc, arc_ref);
    }
    // TODO: erasing them might not have an effecr, but likely hinders unwanted flow
    for (std::set<std::pair<unsigned int, int> >::iterator it = slow_add_sources.begin(); it !=  slow_add_sources.end(); ++it) {
        // remove from copy what could not be used
        wc.erase(arc_ref_cop[wc.arcFromId(it->second)]);
    }
    for (std::set<std::pair<unsigned int, int> >::iterator it = slow_add_drains.begin(); it !=  slow_add_drains.end(); ++it) {
        // remove from copy what could not be used
        wc.erase(arc_ref_cop[wc.arcFromId(it->second)]);
    }
    
    ListDigraph::ArcMap<unsigned int> cii(wc);
    ListDigraph::ArcMap<unsigned int> cio(wc);

    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Push Circular Flow.\n");
    print_graph_debug_copy(std::cout, wc, fc, ces, cet, cel, cs, ct);
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    // the backlink edges are treated special
    expand_cycles(wc, fc, cap, ces, cet, cel, cni, cii, cio, cs, ct, node_ref, arc_ref);
    transfer(wc, cet, cs, ct, cap, fc, arc_ref);
    
    // at this point we have maximal flow and a DAG with expanded cycles
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Before LOOP.\n");
    print_graph_debug_copy(std::cout, wc, fc, ces, cet, cel, cs, ct);
    print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
    #endif
    
    bool new_sources = true;
    unsigned int iterations = 0;
    while (new_sources && iterations <= options::Instance()->get_max_flow_iteration()) { 
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Enter Flow Loop.\n");
        print_graph_debug_copy(std::cout, wc, cap, ces, cet, cel, cs, ct);
        #endif
        
        // remove exon edges who are missing either side node
        for (ListDigraph::ArcIt a(wc); a != INVALID; ) {
            ListDigraph::Arc arc(a);
            ++a; 
        
            if (cet[arc] == edge_types::EXON ) {
                ListDigraph::OutArcIt o(wc, wc.target(arc));
                ListDigraph::InArcIt i(wc, wc.source(arc));
                if (o == INVALID || i == INVALID) {
                     
                    ListDigraph::Node target = wc.target(arc);
                    ListDigraph::Node source = wc.source(arc);

                    wc.erase(arc);

                    ListDigraph::OutArcIt ot(wc,target);
                    ListDigraph::InArcIt it(wc,target);
                    if (it == INVALID && ot == INVALID) {
                        wc.erase(target);
                    }

                    ListDigraph::OutArcIt os(wc,source);
                    ListDigraph::InArcIt is(wc,source);
                    if (is == INVALID && os == INVALID) {
                        wc.erase(source);
                    }     
                }
            }
        }
        // we need to find and remove nodes without any exon edges left on either side
        for (ListDigraph::ArcIt a(wc); a != INVALID; ) {
            ListDigraph::Arc arc(a);
            ++a; 
        
            if (cet[arc] == edge_types::NODE ) {
                ListDigraph::OutArcIt o(wc, wc.target(arc));
                ListDigraph::InArcIt i(wc, wc.source(arc));
                while (o != INVALID && cet[o] == edge_types::HELPER) {
                    ++o;
                }
                while (i != INVALID && cet[i] == edge_types::HELPER) {
                    ++i;
                }
                
                if (o == INVALID && i == INVALID) {
                    
                    ListDigraph::Node target = wc.target(arc);
                    ListDigraph::Node source = wc.source(arc);

                    wc.erase(arc);

                    ListDigraph::OutArcIt ot(wc,target);
                    ListDigraph::InArcIt it(wc,target);
                    if (it == INVALID && ot == INVALID) {
                        wc.erase(target);
                    }

                    ListDigraph::OutArcIt os(wc,source);
                    ListDigraph::InArcIt is(wc,source);
                    if (is == INVALID && os == INVALID) {
                        wc.erase(source);
                    } 
                }
            }
        }
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Reduction step.\n");
        print_graph_debug_copy(std::cout, wc, cap, ces, cet, cel, cs, ct);
        #endif
        
        // now find stretched that are still intact and add sources and drains
        //  we need to add new sources and sinks at start and end stretches for the next iteration
        new_sources = false;
        for (ListDigraph::NodeIt node(wc); node != INVALID; ++node) {

            if (node == cs || node == ct) {
                continue;
            }

            ListDigraph::OutArcIt o(wc, node);
            ListDigraph::InArcIt i(wc, node);

            if (i == INVALID) {
                
                capacity_type count = 0;
                bool forbidden = false;
                for (; o!= INVALID; ++o) {
                    if (wc.target(o) == ct) {
                        forbidden = true;
                        continue;
                    }
                    count += cap[o];
                }
                if (forbidden) {
                    continue;
                }
                
                ListDigraph::Arc na = wc.addArc(cs, node);
                ListDigraph::Arc nag;
                
                bool exists = false;
                ListDigraph::InArcIt a(g, node_ref[node]); 
                for (; a != INVALID; ++a) {
                    if (g.source(a) == node_ref[cs]) {
                        // exists
                        exists = true;
                        break;
                    }
                }
                if (exists) {
                    nag = a;
                } else {
                    nag = g.addArc(node_ref[cs], node_ref[node]);
                }
                
                cap[na] += count;
                
                edge_type[nag] = edge_types::HELPER;
                cet[na] = edge_types::HELPER;
                arc_ref[na] = nag;
                new_sources = true;
                continue;
            }

            if (o == INVALID) {
                
                capacity_type count = 0;
                bool forbidden = false;
                for (; i!= INVALID; ++i) {
                    if (wc.source(i) == cs) {
                        forbidden = true;
                        continue;
                    }
                    count += cap[i];
                }
                if (forbidden) {
                    continue;
                }
                
                ListDigraph::Arc na = wc.addArc(node, ct);
                ListDigraph::Arc nag;
                
                bool exists = false;
                ListDigraph::OutArcIt a(g, node_ref[node]); 
                for (; a != INVALID; ++a) {
                    if (g.target(a) == node_ref[ct]) {
                        // exists
                        exists = true;
                        break;
                    }
                }
                if (exists) {
                    nag = a;
                } else {
                    nag = g.addArc(node_ref[node], node_ref[ct]);
                }

                cap[na] += count;
                
                edge_type[nag] = edge_types::HELPER;
                cet[na] = edge_types::HELPER;
                arc_ref[na] = nag;
                new_sources = true;
            }
        }
  
        // if we have new sources/sinks, push it
        push_graph_flow(wc, cs, ct, cap, fc);
        transfer(wc, cet, cs, ct, cap, fc, arc_ref);

        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Loop Done.\n");
        print_graph_debug_copy(std::cout, g, flow, edge_specifier, edge_type, edge_lengths, s, t);
        #endif
        
        ++iterations;
    }
}


capacity_type basic_coverage_flow::push_graph_flow(ListDigraph &g, ListDigraph::Node s, ListDigraph::Node t, ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<capacity_type> &f) {
    
    Preflow<ListDigraph, ListDigraph::ArcMap<rcount> > preflow(g, cap, s, t);
    preflow.run();
    
    for (ListDigraph::ArcIt a(g); a != INVALID; ++a) {
        f[a] = preflow.flowMap()[a];
    }
    return preflow.flowValue();
}

void basic_coverage_flow::transfer(ListDigraph &wc, ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::Node cs, ListDigraph::Node ct, ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<ListDigraph::Arc> &arc_ref) {
    
    for (ListDigraph::ArcIt a(wc); a != INVALID;) { // loop through current flowgraph
        
        ListDigraph::Arc arc(a);
        ++a;
        
        flow[arc_ref[arc]] += fc[arc];  // add the flow in base graph
        cap[arc] -= fc[arc]; // remove from capacity in copy graph we work on
        fc[arc] = 0;

    }
            
}

void basic_coverage_flow::expand_cycles(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &fc,
        ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int> &cni,
        ListDigraph::ArcMap<unsigned int> &cycle_id_in, ListDigraph::ArcMap<unsigned int> &cycle_id_out,
        ListDigraph::Node s, ListDigraph::Node t,
        ListDigraph::NodeMap<ListDigraph::Node> &node_ref,
        ListDigraph::ArcMap<ListDigraph::Arc> &arc_ref) {

    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        ListDigraph::Arc arc(a);
        if (cet[arc] == edge_types::BACKLINK ) {
            // we found a backlink
            
            ListDigraph::Node in = wc.source(arc);
            ListDigraph::Node out = wc.target(arc);
            
            // protect against node_edges
            ListDigraph::InArcIt ia(wc, in);
            if (cet[ia] == edge_types::NODE) {
                in = wc.source(ia);
            }
            ListDigraph::OutArcIt oa(wc, out);
            if (cet[oa] == edge_types::NODE) {
                out = wc.target(oa);
            }
            
            unsigned int in_index = cni[in];
            unsigned int out_index = cni[out];
            
            // initial filling of the 
            gmap< std::pair<int,int> , flow_paired_path > circle_links;
            std::set<std::pair<int,int> > inserts;
            for (graph_list<paired_exon_group>::iterator circle = raw->chim_circle_bin_list.begin(); circle != raw->chim_circle_bin_list.end(); ++circle) {
                
                // bounds
                unsigned int in_right = circle->left_read->range_end; // jump
                unsigned int in_left = circle->left_read->range_start;
                unsigned int out_left = circle->right_read->range_start; // jump
                unsigned int out_right = circle->right_read->range_end;
                
                if ( in_index != in_right || out_index != out_left) {
                    continue;
                } 
                
                ListDigraph::Arc in_arc, out_arc;
                bool found_arc = false;
                for (ListDigraph::InArcIt ia(wc, in); ia!=INVALID; ++ia) {
                    ListDigraph::Arc inarc(ia);
                    
                    if (cet[inarc] != edge_types::EXON ) {
                        continue;
                    }
                    
                    if ( (in_left > cni[wc.source(inarc)] &&  circle->left_read->bin_mask.is_contained_in(ces[inarc], in_left,  cni[wc.target(inarc)]) ) 
                            || (in_left <= cni[wc.source(inarc)] && ces[inarc].is_contained_in(circle->left_read->bin_mask, cni[wc.source(inarc)], cni[wc.target(inarc)]) ) ) {
                        
                            // we have a hit :)
                            found_arc = true;
                            in_arc = inarc;
                            break;
                    }
                }
                
                if(!found_arc) {
                    continue;
                }
                found_arc = false;
                for (ListDigraph::OutArcIt oa(wc, out); oa!=INVALID; ++oa) {
                    ListDigraph::Arc outarc(oa);
                    if ( (out_right > cni[wc.target(outarc)] &&  ces[outarc].is_contained_in(circle->right_read->bin_mask, cni[wc.source(outarc)], cni[wc.target(outarc)]) ) 
                            || (out_right <= cni[wc.target(outarc)] && circle->right_read->bin_mask.is_contained_in(ces[outarc], cni[wc.source(outarc)], out_right) ) ) {
                        
                            // we have a hit :)
                            found_arc = true;
                            out_arc = outarc;
                            break;
                    }
                }
                if(!found_arc) {
                    continue;
                }
                
                // at this point there is a guarantee we hit a valid pair
                std::pair<int, int> key = std::make_pair(wc.id(in_arc), wc.id(out_arc));
                flow_paired_path& flow = circle_links[key];
                inserts.insert(key);
                flow.frags.push_back(&*circle);
                if (flow.in.size() == 0 ) {
                    flow.in.push_back(wc.id(in_arc));
                    flow.out.push_back(wc.id(out_arc));
                }
            }
            
            for(std::set<std::pair<int,int> >::iterator it = inserts.begin(); it != inserts.end(); ++it) {
                    simplify_ambigous_cycles(wc, ces, cet, cni, *it, circle_links[*it].frags, circle_links);
            }
               
            std::vector<std::pair<int,int> > cap_flow;
            std::vector<std::pair<int,int> > flow_cap;
            std::vector<std::pair<int,int> > cap_cap;
            std::vector< std::pair<std::pair<int,int> , bool> > flow_flow;
            
            // now we have a map of ALL pairs to resolve
            for (gmap< std::pair<int,int> , flow_paired_path >::iterator circ = circle_links.begin(); circ != circle_links.end(); ++circ) {
                
                flow_paired_path& flow = circ->second;
               
                int c_in_id = circ->first.first;
                int c_out_id = circ->first.second;
                
                int c_in_left = cni[wc.source(wc.arcFromId(c_in_id))];
                int c_out_right = cni[wc.target(wc.arcFromId(c_out_id))];
                
                // test IN first
                bool is_in_flow_limiting, is_in_circle;
                if (c_in_left <= out_index) {
                    // we search for the edge adjourning to the jump end and test if jump is in found path
                    for(graph_list<int>::reverse_iterator a_id = flow.in.rbegin(); a_id != flow.in.rend(); ++a_id) {
                        ListDigraph::Arc arc = wc.arcFromId(*a_id);
                        if(cni[wc.target(arc)] > out_index) {
                            // we have the actual overlapping arc
                            if(wc.source(arc) == out) {
                                is_in_flow_limiting = true;
                                is_in_circle = true;
                            } else {
                                is_in_flow_limiting = true;
                                is_in_circle = false;
                            }
                            break;
                        }
                    }
                } else {
                    // this means we have to test if this could be in a cycle
                    ListDigraph::Node in_start = wc.source(wc.arcFromId(c_in_id));
                    
                    Bfs<ListDigraph> bfs(wc);
                    bfs.init();
                    bfs.addSource(out);
                    bfs.start(in_start);
                    
                    if (bfs.reached(in_start)) {
                        is_in_flow_limiting = false;
                        is_in_circle = true;
                    } else {
                        is_in_flow_limiting = true;
                        is_in_circle = false;
                    }
                }
                // test OUT second
                bool is_out_flow_limiting, is_out_circle;
                if (c_out_right >= in_index) {
                    // we search for the edge adjourning to the jump end and test if jump is in found path
                    for(graph_list<int>::reverse_iterator a_id = flow.out.rbegin(); a_id != flow.out.rend(); ++a_id) {
                        ListDigraph::Arc arc = wc.arcFromId(*a_id);
                        if(cni[wc.source(arc)] < in_index) {
                            // we have the actual overlapping arc
                            if(wc.target(arc) == in) {
                                is_out_flow_limiting = true;
                                is_out_circle = true;
                            } else {
                                is_out_flow_limiting = true;
                                is_out_circle = false;
                            }
                            break;
                        }
                    }
                } else {
                    // this means we have to test if this could be in a cycle
                    ListDigraph::Node out_end = wc.target(wc.arcFromId(c_in_id));
                    
                    Bfs<ListDigraph> bfs(wc);
                    bfs.init();
                    bfs.addSource(out_end); // source
                    bfs.start(in); // drain
                    
                    if (bfs.reached(in)) {
                        is_out_flow_limiting = false;
                        is_out_circle = true;
                    } else {
                        is_out_flow_limiting = true;
                        is_out_circle = false;
                    }
                }
                
                if (!is_in_circle && !is_out_circle) {
                    continue; // there is no possible circulation to add under this condition
                }

                if (is_in_flow_limiting && is_out_flow_limiting) {
                    flow_flow.push_back(std::make_pair(circ->first, is_in_circle));
                } else if (is_in_flow_limiting && !is_out_flow_limiting) {
                    flow_cap.push_back(circ->first);
                } else if ( !is_in_flow_limiting && is_out_flow_limiting) {
                    cap_flow.push_back(circ->first);
                } else {
                    
                    bool true_cycle;
                    // we need to test if this is actually an endless cycle
                    if ( c_in_left > c_out_right) {
                        // we have a proper separation, wee need BFS
                        Bfs<ListDigraph> bfs(wc);
                        bfs.init();
                        bfs.addSource(wc.target(wc.arcFromId(c_out_id))); // source
                        bfs.start(wc.source(wc.arcFromId(c_in_id))); // drain

                        true_cycle = bfs.reached(wc.source(wc.arcFromId(c_in_id)));

                    } else {
                        // in this case there is an overlap in arcs
                        
                        true_cycle = true;
                        //make get initial overlap
                        int id = flow.in.back();
                        graph_list<int>::reverse_iterator a_id = flow.out.rbegin();
                        for (; a_id != flow.out.rend() && *a_id!= id; ++a_id) { }
                        // this is the first overlap
                        
                        graph_list<int>::reverse_iterator in_it = flow.out.rbegin();
                        graph_list<int>::iterator out_it = --a_id.base();
                        for (; out_it != flow.out.end(); ++out_it, ++in_it) {
                            if (*in_it!= *out_it) {
                                true_cycle = false;
                                break;
                            }
                        }
                    }
                    
                    if (true_cycle) {
                        cap_cap.push_back(circ->first);
                    } else {
                        // this qualifies as both as no circle through both can be established directly
                        cap_flow.push_back(circ->first);
                        flow_cap.push_back(circ->first); 
                    } 
                }
            }
            
            // now we have everything populated as expected
            
            // we have the maximal flow supported by the known flow in and flow out particles
            // time to push the actual circular flow (limited further by rest capacities)

            // starting at exactly known circles
            for(std::vector<std::pair<std::pair<int,int>, bool> >::iterator ff = flow_flow.begin(); ff != flow_flow.end(); ++ff) {
                
                flow_paired_path& f = circle_links[ff->first];
                if (ff->second) { // the in path covers everything
                    
                    capacity_type lowest = cap[arc] - fc[arc];
                    for(graph_list<int>::iterator p = f.in.begin(); p!= f.in.end(); ++p) {
                        if (cni[wc.source(wc.arcFromId(*p))] >= out_index) {
                            capacity_type a = cap[wc.arcFromId(*p)] - fc[wc.arcFromId(*p)];
                            if (a < lowest) {
                                lowest = a;
                            }
                        } else {
                            capacity_type a = fc[wc.arcFromId(*p)];
                            if (a < lowest) {
                                lowest = a;
                            }
                        }
                    }
                    for(graph_list<int>::iterator p = f.out.begin(); p!= f.out.end(); ++p) {
                        capacity_type a = fc[wc.arcFromId(*p)];
                        if (a < lowest) {
                            lowest = a;
                        }
                    }
                    
                    fc[arc] += lowest;
                    for(graph_list<int>::iterator p = f.in.begin(); p!= f.in.end() && wc.target(wc.arcFromId(*p)) != out; ++p) {
                        fc[wc.arcFromId(*p)] += lowest;
                    }
                    
                } else {
                    
                    capacity_type lowest = cap[arc] - fc[arc];
                    for(graph_list<int>::iterator p = f.out.begin(); p!= f.out.end(); ++p) {
                        if (cni[wc.target(wc.arcFromId(*p))] <= in_index) {
                            capacity_type a = cap[wc.arcFromId(*p)] - fc[wc.arcFromId(*p)];
                            if (a < lowest) {
                                lowest = a;
                            }
                        } else {
                            capacity_type a = fc[wc.arcFromId(*p)];
                            if (a < lowest) {
                                lowest = a;
                            }
                        }
                    }
                    for(graph_list<int>::iterator p = f.in.begin(); p!= f.in.end(); ++p) {
                        capacity_type a = fc[wc.arcFromId(*p)];
                        if (a < lowest) {
                            lowest = a;
                        }
                    }
                    
                    fc[arc] += lowest;
                    for(graph_list<int>::iterator p = f.out.begin(); p!= f.out.end() && wc.source(wc.arcFromId(*p)) != in; ++p) {
                        fc[wc.arcFromId(*p)] += lowest;
                    }
                    
                }
            }
            
            // we establish the maximum flow first

            //we need to find the minimal supporting flow each side
            capacity_type min_in = 0;
            capacity_type min_out = 0;

            // we iterate over 
            for(std::vector<std::pair<int,int> >::iterator cf = cap_flow.begin(); cf != cap_flow.end(); ++cf) {
                // 
                flow_paired_path& f = circle_links[*cf];
                graph_list<int>::iterator a_id = f.out.begin();
                capacity_type min = fc[wc.arcFromId(*a_id)];
                ++a_id;
                for (; a_id != f.out.end(); ++a_id) {
                    capacity_type af = fc[wc.arcFromId(*a_id)];
                    if (af < min) {
                        min = af;
                    }
                }
                min_out += min;
            }
            for(std::vector<std::pair<int,int> >::iterator fci = flow_cap.begin(); fci != flow_cap.end(); ++fci) {
                // 
                flow_paired_path& f = circle_links[*fci];
                graph_list<int>::iterator a_id = f.in.begin();
                capacity_type min = fc[wc.arcFromId(*a_id)];
                ++a_id;
                for (; a_id != f.in.end(); ++a_id) {
                    capacity_type af = fc[wc.arcFromId(*a_id)];
                    if (af < min) {
                        min = af;
                    }
                }
                min_in += min;
            }
            
            // now we need to process the actual flow problems
            
            // Flow_Cap
            ListDigraph::ArcMap<capacity_type> alt_cap(wc);
            ListDigraph::ArcMap<capacity_type> new_flow(wc);
            create_leftover_capacities(wc, fc, cap, alt_cap);
            marked_capacities_in(wc, alt_cap, flow_cap, circle_links);
            min_in = std::min(min_in, cap[arc] - fc[arc]);
            capacity_type val = push_circle_flow(wc, alt_cap, new_flow, min_in, out, in);
            add_flow(wc,fc,new_flow);
            fc[arc] += val;
            
            // Cap_Flow
            //alt_cap = ListDigraph::ArcMap<capacity_type>(wc);  // Automatically reset
            //new_flow = ListDigraph::ArcMap<capacity_type>(wc);
            create_leftover_capacities(wc, fc, cap, alt_cap);
            marked_capacities_out(wc, alt_cap, cap_flow, circle_links);
            min_out = std::min(min_out, cap[arc] - fc[arc]);
            val = push_circle_flow(wc, alt_cap, new_flow, min_out, out, in);
            add_flow(wc,fc,new_flow);
            fc[arc] += val;
            
            // cap_cap
            //alt_cap = ListDigraph::ArcMap<capacity_type>(wc); // Automatically reset
            //new_flow = ListDigraph::ArcMap<capacity_type>(wc);
            create_leftover_capacities(wc, fc, cap, alt_cap);
            marked_capacities_both(wc, alt_cap, cap_cap, circle_links);
            capacity_type min = std::min(std::min(min_out, min_in), cap[arc] - fc[arc]);
            val = push_circle_flow(wc, alt_cap, new_flow, min, out, in);
            add_flow(wc,fc,new_flow);
            fc[arc] += val;
            
            // we have the max flow, create edges in such a way that we have again a DAG
            capacity_type max_cap = fc[arc];
            unsigned int cycle_id = cycle_id_in[arc];
            g.erase(arc_ref[arc]);
            wc.erase(arc);
            ListDigraph::Arc ns = wc.addArc(s, out);
            ListDigraph::Arc gns = g.addArc(node_ref[s],  node_ref[out]);
            arc_ref[ns] = gns;
            fc[ns] = max_cap;
            ListDigraph::Arc nt = wc.addArc(in, t);
            ListDigraph::Arc gnt = g.addArc(node_ref[in],  node_ref[t]);
            arc_ref[nt] = gnt;
            fc[nt] = max_cap;
            cet[ns] = edge_types::HELPER;
            cycle_id_in[ns] = cycle_id;
            cet[nt] = edge_types::HELPER;
            cycle_id_out[nt] = cycle_id;
        }
    }
    
}


void basic_coverage_flow::marked_capacities_in(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        std::vector<std::pair<int,int> > &flow_cap, gmap< std::pair<int,int> , flow_paired_path > &circle_links) {
    
    // mark arcs that can be deleted
    ListDigraph::ArcMap<circle_marker> markers(wc);
    for(std::vector<std::pair<int,int> >::iterator ai = flow_cap.begin(); ai != flow_cap.end(); ++ai){
        flow_paired_path& f = circle_links[*ai];
        for (graph_list<int>::iterator p = f.in.begin(); p!= f.in.end(); ++p) {
            ListDigraph::Arc arc = wc.arcFromId(*p);
            markers[arc] = SAVE;
            
            for (ListDigraph::InArcIt oa(wc, wc.target(arc)); oa!=INVALID; ++oa) {
                if (markers[arc] != SAVE) { 
                    markers[arc] = ELIMINATE;
                }
            }
        }
    }
    
    // remove capacity for eliminated edges
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        if(markers[a] == ELIMINATE) {
            alt_cap[a] = 0;
        }
    }
}


void basic_coverage_flow::marked_capacities_out(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        std::vector<std::pair<int,int> > &cap_flow, gmap< std::pair<int,int> , flow_paired_path > &circle_links) {
    
    // mark arcs that can be deleted
    ListDigraph::ArcMap<circle_marker> markers(wc);
    for(std::vector<std::pair<int,int> >::iterator ai = cap_flow.begin(); ai != cap_flow.end(); ++ai){
        flow_paired_path& f = circle_links[*ai];
        for (graph_list<int>::iterator p = f.out.begin(); p!= f.out.end(); ++p) {
            ListDigraph::Arc arc = wc.arcFromId(*p);
            markers[arc] = SAVE;
            
            for (ListDigraph::OutArcIt oa(wc, wc.source(arc)); oa!=INVALID; ++oa) {
                if (markers[arc] != SAVE) { 
                    markers[arc] = ELIMINATE;
                }
            }
        }
    }
    
    // remove capacity for eliminated edges
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        if(markers[a] == ELIMINATE) {
            alt_cap[a] = 0;
        }
    }
}

void basic_coverage_flow::marked_capacities_both(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        std::vector<std::pair<int,int> > &cap_cap, gmap< std::pair<int,int> , flow_paired_path > &circle_links) {
    
    // mark arcs that can be deleted
    ListDigraph::ArcMap<circle_marker> markers(wc);
    for(std::vector<std::pair<int,int> >::iterator ai = cap_cap.begin(); ai != cap_cap.end(); ++ai){
        flow_paired_path& f = circle_links[*ai];
        for (graph_list<int>::iterator p = f.out.begin(); p!= f.out.end(); ++p) {
            ListDigraph::Arc arc = wc.arcFromId(*p);
            markers[arc] = SAVE;
            
            for (ListDigraph::OutArcIt oa(wc, wc.source(arc)); oa!=INVALID; ++oa) {
                if (markers[arc] != SAVE) { 
                    markers[arc] = ELIMINATE;
                }
            }
        }
    }
    
    for(std::vector<std::pair<int,int> >::iterator ai = cap_cap.begin(); ai != cap_cap.end(); ++ai){
        flow_paired_path& f = circle_links[*ai];
        for (graph_list<int>::iterator p = f.in.begin(); p!= f.in.end(); ++p) {
            ListDigraph::Arc arc = wc.arcFromId(*p);
            markers[arc] = SAVE;
            
            for (ListDigraph::InArcIt oa(wc, wc.target(arc)); oa!=INVALID; ++oa) {
                if (markers[arc] != SAVE) { 
                    markers[arc] = ELIMINATE;
                }
            }
        }
    }
    
    // remove capacity for eliminated edges
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        if(markers[a] == ELIMINATE) {
            alt_cap[a] = 0;
        }
    }
    
}

capacity_type basic_coverage_flow::push_circle_flow(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &alt_cap,
        ListDigraph::ArcMap<capacity_type> &new_flow, capacity_type max_flow, ListDigraph::Node &s, ListDigraph::Node &t) {
    
    ListDigraph::Node new_start = wc.addNode();
    ListDigraph::Node new_drain = wc.addNode();
    
    ListDigraph::Arc new_start_arc = wc.addArc(new_start, s);
    ListDigraph::Arc new_drain_arc = wc.addArc(t, new_drain);
    
    alt_cap[new_start_arc] = max_flow;
    alt_cap[new_drain_arc] = max_flow;
    
    capacity_type ret = push_graph_flow(wc, new_start, new_drain, alt_cap, new_flow);
    
    wc.erase(new_drain_arc);
    wc.erase(new_drain);
    wc.erase(new_start_arc);
    wc.erase(new_start);
    
    return ret;
}

void basic_coverage_flow::simplify_ambigous_cycles(ListDigraph &wc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::NodeMap<unsigned int> &cni,
        std::pair<int,int> arc_pair, 
        graph_list<paired_exon_group* > &evidences,
        gmap< std::pair<int,int> , flow_paired_path > & circles) {
    
    ListDigraph::Arc in_arc = wc.arcFromId(arc_pair.first);
    ListDigraph::Arc out_arc = wc.arcFromId(arc_pair.second);
    
    unsigned int in_left = cni[ wc.source(in_arc)];
    unsigned int out_right = cni[ wc.target(out_arc)];
    
    // test if arcs can be extended
    bool extend_in = true;
    bool extend_out = true;
    for (graph_list<paired_exon_group* >::iterator cit = evidences.begin(); cit != evidences.end(); ++cit) {
        
        if ( (*cit)->left_read->range_start >= in_left) {
               extend_in = false;
        }
        if ( (*cit)->right_read->range_end <= out_right) {
               extend_out = false;
        }
    }
    
    // in first
    if (extend_in && extend_out) {
        
        // extend on both counts!
        
        // protect against node edges
        ListDigraph::InArcIt ia(wc, wc.source(in_arc));
        if (cet[ia] == edge_types::NODE) {
            in_arc = ia;
        }
        ListDigraph::OutArcIt oa(wc, wc.target(out_arc));
        if (cet[oa] == edge_types::NODE) {
            out_arc = oa;
        }
        
        std::set< std::pair<int, int> > new_pairs;
        for (graph_list<paired_exon_group* >::iterator mem_it = evidences.begin(); mem_it != evidences.end(); ++mem_it) {
            
            
            ListDigraph::Arc new_left, new_right;
            bool init_left = false;
            bool init_right = false;
            // test left advancement
            for (ListDigraph::InArcIt ia(wc, wc.source(in_arc)); ia!=INVALID; ++ia) {
                
                if (cet[ia] != edge_types::EXON) {
                    continue;
                }
                
                unsigned int left_border = (*mem_it)->left_read->range_start;
                if ( left_border < cni[wc.source(ia)]) {
                    left_border = cni[wc.source(ia)];
                } 
                
                if ( (*mem_it)->left_read->bin_mask.is_contained_in(ces[ia], left_border, cni[wc.target(ia)]) ) {
                    new_left = ia;
                    init_left = true;
                }
            }
            
            // test right advancement
            for (ListDigraph::OutArcIt oa(wc, wc.target(out_arc)); oa!=INVALID; ++oa) {
                
                if (cet[oa] != edge_types::EXON) {
                   continue;
                }
                
                unsigned int right_border = (*mem_it)->right_read->range_end;
                if ( right_border > cni[wc.target(oa)]) {
                    right_border = cni[wc.target(oa)];
                } 
                
                if ( (*mem_it)->right_read->bin_mask.is_contained_in(ces[oa], cni[wc.source(out_arc)], right_border) ) {
                    new_right = oa;
                    init_right = true;
                }
                
            }
            
            if(!init_left || !init_right) {
                continue;
            }
            
            // so here we have the new arc
            std::pair<int, int> new_pair = std::make_pair(wc.id(new_left), wc.id(new_right));
            new_pairs.insert(new_pair);
            flow_paired_path & flow = circles[new_pair];
            flow.frags.push_back(*mem_it);
            if (flow.in.empty()) { // first time use
                flow.in = circles[arc_pair].in;
                flow.out = circles[arc_pair].out;
                flow.in.push_back(wc.id(new_left));
                flow.out.push_back(wc.id(new_right));
            }
        }
        
        circles.erase(arc_pair); // old pair
        
        for (std::set<std::pair<int, int> >::iterator np = new_pairs.begin(); np != new_pairs.end(); ++np) {
            simplify_ambigous_cycles(wc, ces, cet, cni, *np, circles[*np].frags,  circles);
        }
        
    } else if (extend_in && !extend_out) {
        // we tested such that EVERY in_edge can be closer specified
        
        // protect against node edges
        ListDigraph::InArcIt ia(wc, wc.source(in_arc));
        if (cet[ia] == edge_types::NODE) {
            in_arc = ia;
        }
        
        std::set<int> new_lefts;
        bool init_left = false;
        for (ListDigraph::InArcIt ia(wc, wc.source(in_arc)); ia!=INVALID; ++ia) {
            
            if (cet[ia] != edge_types::EXON) {
                continue;
            }
            
            int new_left = wc.id(ia);
            for (graph_list<paired_exon_group* >::iterator mem_it = evidences.begin(); mem_it != evidences.end(); ++mem_it) {
                
                unsigned int left_border = (*mem_it)->left_read->range_start;
                if ( left_border < cni[wc.source(ia)]) {
                    left_border = cni[wc.source(ia)];
                } 
                
                if ( (*mem_it)->left_read->bin_mask.is_contained_in(ces[ia], left_border, cni[wc.target(ia)]) ) {
                    // we have an extended hit
                    flow_paired_path & flow = circles[std::make_pair(new_left, arc_pair.second)];
                    flow.frags.push_back(*mem_it);
                    if (flow.in.empty()) { // first time use
                        flow.in = circles[arc_pair].in;
                        flow.out = circles[arc_pair].out;
                        flow.in.push_back(new_left);
                    }
                    
                    new_lefts.insert(new_left);
                    init_left = true;
                }
            }
        }
        
        if(!init_left) {
                return;
        }
        
        circles.erase(arc_pair); // old pair
        
        // recursive further updates
        for (std::set<int>::iterator nl = new_lefts.begin(); nl != new_lefts.end(); ++nl) {
            std::pair<int, int> p = std::make_pair(*nl, arc_pair.second);
            simplify_ambigous_cycles(wc, ces, cet, cni, p, circles[p].frags,  circles);
        }
        
    } else if (!extend_in && extend_out) {
        // we tested such that EVERY out_edge can be closer specified
        
        // protect against node edges
        ListDigraph::OutArcIt oa(wc, wc.target(out_arc));
        if (cet[oa] == edge_types::NODE) {
            out_arc = oa;
        }
        
        std::set<int> new_rights;
        bool init_right = false;
        for (ListDigraph::OutArcIt oa(wc, wc.target(out_arc)); oa!=INVALID; ++oa) {
            
            if (cet[oa] != edge_types::EXON) {
               continue;
            }
            
            int new_right= wc.id(oa);
            for (graph_list<paired_exon_group* >::iterator mem_it = evidences.begin(); mem_it != evidences.end(); ++mem_it) {
                
                unsigned int right_border = (*mem_it)->right_read->range_end;
                if ( right_border > cni[wc.target(oa)]) {
                    right_border = cni[wc.target(oa)];
                } 
                
                if ( (*mem_it)->right_read->bin_mask.is_contained_in(ces[oa], cni[wc.source(out_arc)], right_border) ) {
                    // we have an extended hit
                    flow_paired_path & flow = circles[std::make_pair(arc_pair.first, new_right)];
                    flow.frags.push_back(*mem_it);
                    if (flow.in.empty()) { // first time use
                        flow.in = circles[arc_pair].in;
                        flow.out = circles[arc_pair].out;
                        flow.out.push_back(new_right);
                    }
                    
                    new_rights.insert(new_right);
                    init_right = true;
                }
            }
        }
        
        if(!init_right) {
                return;
        }
        
        circles.erase(arc_pair); // old pair
        
        // recursive further updates
        for (std::set<int>::iterator nl = new_rights.begin(); nl != new_rights.end(); ++nl) {
            std::pair<int, int> p = std::make_pair(arc_pair.first, *nl);
            simplify_ambigous_cycles(wc, ces, cet, cni, p, circles[p].frags,  circles);
        }
        
    }

} 

void basic_coverage_flow::create_leftover_capacities(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &flow, ListDigraph::ArcMap<capacity_type> &cap, ListDigraph::ArcMap<capacity_type> &new_cap) {
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        new_cap[a] = cap[a] - flow[a];
    }
}


void basic_coverage_flow::add_flow(ListDigraph &wc, ListDigraph::ArcMap<capacity_type> &orig_flow, ListDigraph::ArcMap<capacity_type> &add_flow) {
    
    for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
        orig_flow[a] += add_flow[a];
    }
}


// ################## simple printing ##################

void basic_coverage_flow::print_graph(std::ostream &os) {
   
    if (meta->size != 1) {   
        // there is a graph!
        digraphWriter(g, os)
            .arcMap("edge_specifier", edge_specifier)
            .arcMap("edge_type", edge_type)
            //.arcMap("edge_lengths", edge_lengths)
            .arcMap("coverage", capacity)
            .node("source", s)
            .node("drain", t)
            .run();  
    
    }
    
     for (int j = 0; j < meta->size; j++) {
        os << "Exon" << std::to_string(j) << " " << std::to_string(meta->exons[j].left) << "-" << std::to_string(meta->exons[j].right)  << "\n";
    }
    
    
    // also print singles
    for( std::set<single_exon >::iterator i = single_exons.begin(); i != single_exons.end(); ++i) {
        os << "Single Exon from " << std::to_string(meta->exons[i->meta].left) << "-" << std::to_string(meta->exons[i->meta].right)  << "\n";
    }
}