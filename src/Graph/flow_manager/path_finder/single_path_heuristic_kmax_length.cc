/* 
 * File:   multi_path_heuristic.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 8, 2016, 3:11 PM
 */

#include "single_path_heuristic_kmax_length.h"

#include <lemon/lgf_writer.h>

single_path_heuristic_kmax_length::single_path_heuristic_kmax_length(ListDigraph& wc,
        ListDigraph::Node& s,
        ListDigraph::Node& t,
        ListDigraph::ArcMap<capacity_type>& cfc,
        ListDigraph::ArcMap<capacity_mean> &mc,
        ListDigraph::ArcMap<exon_edge>& ces,
        ListDigraph::ArcMap<edge_types::edge_type>& cet,
        ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int>& cni,
        ListDigraph::ArcMap<arc_bridge>& kp,
        ListDigraph::ArcMap<arc_back_bridge>& kbp,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in,
        ListDigraph::ArcMap<unsigned int>& cycle_id_out,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        unsigned int size) 
: path_finder(wc, s, t, cfc, mc, ces, cet, cel, cni,kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size) {
}


single_path_heuristic_kmax_length::~single_path_heuristic_kmax_length() {
}

void single_path_heuristic_kmax_length::extract_transcripts(alternative_transcript_collection& results) {

   while (true) {
        
//       logger::Instance()->debug("Extract Path\n");
//       digraphWriter(wc, std::cout)
//            .arcMap("edge_specifier", ces)
//            .arcMap("edge_type", cet)
//            .arcMap("flow", cfc)
//            .arcMap("length", cel)
//            .node("source", s)
//            .node("drain", t)
//            .run();  
        
        ListDigraph::ArcIt has_arc(wc);
        if ( has_arc == INVALID ) {
            // no more arcs left in the graph, we are at the end!  
            return;
        }

        // we find the max min path of the possible paths
        std::deque<ListDigraph::Arc> path;
        find_min_max(path);
        
        #ifdef ALLOW_DEBUG
        logger::Instance()->debug("Path found " + std::to_string( path.size() )+ " ");
        for (std::deque<ListDigraph::Arc>::iterator it = path.begin(); it!= path.end(); ++it) {
              logger::Instance()->debug(std::to_string( wc.id(*it) )+ ", ");
        }
        logger::Instance()->debug("\n");
        #endif
        
        capacity_type cap = add_path_to_collection(path, results, wc, cfc, ces, cet, kp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId);
      
        // now we need to remove this path from the copy
        capacity_mean mean = remove_full_path(path, cap, wc, cfc, mc, ces, cet, cel, kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId); 
        results.transcripts.back()->mean = mean.mean; 
    }
    
}


void single_path_heuristic_kmax_length::find_min_max( std::deque<ListDigraph::Arc> &path) {
    
    // since we have a DAG we can operate in linear time using a topological order
    // we get such a order using DFS
    std::deque<ListDigraph::Node> top_order;
    pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > visitor(top_order);
    DfsVisit<ListDigraph, pff::TopologicalListVisitor<ListDigraph, std::deque<ListDigraph::Node> > > dfs(wc, visitor);
    dfs.init();
    dfs.addSource(s);
    dfs.start();
    
    ListDigraph::NodeMap<unsigned int> tsi(wc);
    pff::order_to_nodemap(top_order, tsi);
    std::deque<ListDigraph::Node>::iterator top_end = top_order.end();
    --top_end;
    
    // we need to create the dynamic programming tables (forward + back)
    // by design the drain is the last node in the list
    dyn_paths* min_max = new dyn_paths[top_order.size()]();
    
    // we compute nodes in order of the topological sorting (forward)
    int i = 0;
    std::deque<ListDigraph::Node>::iterator n = top_order.begin();
    
    // we handle source differently
    for (ListDigraph::OutArcIt a(wc,*n); a!= INVALID; ++a) {
        evidence_element se;
        se.arcs.push_back(a);
        se.cap = cfc[a];
        se.length =  cel[a].first_exon + cel[a].middle + cel[a].last_exon;
        
        unsigned int j = tsi[wc.target(a)];
        add_kmax(min_max[j][wc.id(a)], se);
    }
    
    ++i;
    ++n;
    for (; n!= top_end; ++n, ++i ) {
        
//        logger::Instance()->debug("Work on Node " + std::to_string(wc.id(*n)) + "\n");
        
        // for the rest we update klist lists for each arc, so triple loop: Arcs To Add, for this Arc potential previous edges, elements in Pareto front this edge
        for (ListDigraph::OutArcIt a(wc,*n); a!= INVALID; ++a) {
            
            for (dyn_paths::iterator p_it = min_max[i].begin(); p_it != min_max[i].end(); ++p_it) {
                
                ListDigraph::Arc arc = wc.arcFromId(p_it->first);
                if (kp[arc].is_evidenced_path() && !kp[arc].has_path_evidence(wc.id(a))) {
                    continue; // skip this one
                }
                
                rcount brigde_count = 0;
                if (kp[arc].is_evidenced_path()) {
                    brigde_count = kp[arc].get_path_evidence(wc.id(a));
                }
                
                rpos added_length = cel[a].middle + cel[a].last_exon;
                
                for (klist::iterator ee = p_it->second.begin(); ee != p_it->second.end(); ++ee) {
                    evidence_element new_element = *ee; // deep copy
                    new_element.arcs.push_back(a);
                    new_element.brigde += brigde_count;
                    new_element.cap = std::min(new_element.cap, cfc[a]);
                    
                    if (new_element.length == 0) {
                        new_element.length = cel[a].first_exon;
                    }
                    new_element.length += added_length;
                    
                    unsigned int j = tsi[wc.target(a)];
                    add_kmax(min_max[j][wc.id(a)], new_element);
                    
                }
            }
        }
    }
//      logger::Instance()->debug("End Node " + std::to_string(wc.id(*n)) + "\n");
    // in the end we merge all possible arcs to last node
    klist end_list;
    for (dyn_paths::iterator p_it = min_max[i].begin(); p_it != min_max[i].end(); ++p_it) {
        for (klist::iterator ee = p_it->second.begin(); ee != p_it->second.end(); ++ee) {
            add_kmax(end_list, *ee);
        }
    }
    
    
    
    klist::iterator ee = end_list.begin();
    klist::iterator max_cap = ee;
    ++ee;
    for (; ee != end_list.end(); ++ee) {
        if (ee->cap > max_cap->cap) {
            max_cap = ee;
        }
    }
    
    // prune too worse ones!
    klist::iterator filter = end_list.begin();
    ++filter;
    for (; filter != end_list.end(); ) {
        if (ee->cap * 100 / float(max_cap->cap) < 80) {
            filter = end_list.erase(filter);
        } else {
            ++filter;
        }
    }
    
    ee = end_list.begin();
    klist::iterator max = ee;
    ++ee;
    for (; ee != end_list.end(); ++ee) {
        if (ee->length > max->length) {
            max = ee;
        }
    }
    path = max->arcs;
    
    
//    ee = end_list.begin();
//    klist::iterator choice = ee;
//    rpos maxlen = 0;
//    for (; ee != end_list.end(); ++ee) {
//        if (ee->cap * 100 / float(max->cap) > 80 && ee->length > maxlen ) {
//            choice = ee;
//            maxlen = ee->length;
//        }
//    }
//    path = choice->arcs;
    
    // TODO: so far lists are copied a few times, but reference are not helpful as we extend intermediate list that need to stay intact... 
    
    delete[] min_max;
}

void single_path_heuristic_kmax_length::add_kmax(klist &k_list, evidence_element& insert) {
    
//    logger::Instance()->debug("Insert kmax: Cap " + std::to_string(insert.cap) + " Len " + std::to_string(insert.length)+"\n");
    
    if (k_list.size() < k) {
        k_list.push_back(insert);
    } else {
        
        klist::iterator min = k_list.begin();
        klist::iterator it = k_list.begin();
        ++it;
        for (; it != k_list.end(); ++it ) {
            if (it->cap < min->cap) {
                min = it;
            }
        }
        
        if (insert.cap > min->cap) {
            k_list.push_back(insert);
            k_list.erase(min);
        }
    }
}


//void single_path_heuristic_kmax_length::add_kmax(klist &k_list, evidence_element& insert) {
//    
//    logger::Instance()->debug("Insert kmax: Cap " + std::to_string(insert.cap) + " Len " + std::to_string(insert.length)+"\n");
//    
//    if (k_list.size() < k) {
//        k_list.push_back(insert);
//    } else {
//        
//        klist::iterator min = k_list.begin();
//        klist::iterator it = k_list.begin();
//        ++it;
//        for (; it != k_list.end(); ++it ) {
//            if (it->length < min->length) {
//                min = it;
//            }
//        }
//        
//        if (insert.length > min->length) {
//            k_list.push_back(insert);
//            k_list.erase(min);
//        }
//    }
//}