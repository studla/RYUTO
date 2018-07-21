/* 
 * File:   single_path_heuristic.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 8, 2016, 3:03 PM
 */

#include "single_path_heuristic_no_neighbour.h"

#include <lemon/lgf_writer.h>

single_path_heuristic_no_neighbour::single_path_heuristic_no_neighbour(ListDigraph& wc,
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
        unsigned int size, exon_meta* meta) 
: path_finder(wc, s, t, cfc, mc, ces, cet, cel, cni,kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size), meta(meta) {
}

single_path_heuristic_no_neighbour::~single_path_heuristic_no_neighbour() {
}



void single_path_heuristic_no_neighbour::extract_transcripts(alternative_transcript_collection& results) {
    
    while (true) {
        
//       logger::Instance()->debug("Extract Path\n");
//       digraphWriter(wc, std::cout)
//            .arcMap("edge_specifier", ces)
//            .arcMap("edge_type", cet)
//            .arcMap("flow", cfc)
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

void single_path_heuristic_no_neighbour::find_min_max(std::deque<ListDigraph::Arc> &path) {
    
    // logging stuff to remove
//        logger::Instance()->debug("Min Max Pre \n"); 
//        for (ListDigraph::ArcIt a(wc); a != INVALID; ++a) {
//            
//           logger::Instance()->debug("Arc " + std::to_string(wc.id(a))); 
//           for (arc_bridge::iterator ab = kp[a].begin(); ab != kp[a].end(); ++ab) {
//               logger::Instance()->debug(" f " + std::to_string(ab->first)); 
//           }
//           for (arc_back_bridge::iterator ab = kbp[a].begin(); ab != kbp[a].end(); ++ab) {
//               logger::Instance()->debug(" b " + std::to_string(*ab)); 
//           }
//           logger::Instance()->debug("\n" );  
//        }
//        
        // end logging stuff
    
    
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
    
//    for (ListDigraph::NodeIt a(wc); a!= INVALID; ++a) {
//        logger::Instance()->debug("Backmap " + std::to_string(wc.id(a))+ " " + std::to_string(tsi[a]) + " \n");
//    } 
    
    struct dyn_info {
        unsigned int back_node = 0;
        ListDigraph::Arc arc;
        capacity_type cap = 0;
        capacity_type last_cap = 0;
        unsigned int neighbour_count = 0;
    };
    typedef path_evidence_map<int, dyn_info> dyn_paths;
    
    // we need to create the dynamic programming tables (forward + back)
    // by design the drain is the last node in the list
    dyn_paths* min_max = new dyn_paths[top_order.size()]();
    
    // we compute nodes in order of the topological sorting (forward)
    int i = 0;
    for (std::deque<ListDigraph::Node>::iterator n = top_order.begin(); n!= top_end; ++n, ++i ) {
//        logger::Instance()->debug("====== " + std::to_string(i) + "\n");
        for (ListDigraph::OutArcIt a(wc,*n); a!= INVALID; ++a) {
//            logger::Instance()->debug("---- "+std::to_string(wc.id(a)) +"\n");
            
            // we get max min from all values of last nodes
            capacity_type max = 0;
            capacity_type last_cap_max = 0;
            ListDigraph::Arc max_arc;
            unsigned int neighbour_min = std::numeric_limits<unsigned int>::max();
            
            for (dyn_paths::iterator p_it = min_max[i].begin(); p_it != min_max[i].end(); ++p_it) {
                
//                logger::Instance()->debug("fmax " + std::to_string(p_it->first) + " \n");
                
                ListDigraph::Arc arc = wc.arcFromId(p_it->first);
                
                if (kp[arc].is_evidenced_path() && !kp[arc].has_path_evidence(wc.id(a))) {
//                    logger::Instance()->debug("Skip \n");
                    continue; // skip this one
                }
                
                if ( p_it->second.neighbour_count < neighbour_min 
                        || (p_it->second.neighbour_count == neighbour_min && p_it->second.cap > max)
                        || (p_it->second.neighbour_count == neighbour_min && p_it->second.cap == max && p_it->second.last_cap > last_cap_max) ) {
                    max = p_it->second.cap;
                    max_arc = arc;
                    last_cap_max = p_it->second.last_cap;
                    neighbour_min = p_it->second.neighbour_count;
//                    logger::Instance()->debug("Setmax \n");
                }
            }
            
            if (max == 0) {
                
                if (*n != s) {
//                     logger::Instance()->debug("Continue \n");
                    continue;
                }
                neighbour_min = find_neighbouring(a);
                max = cfc[a];
            } else {
                neighbour_min += find_neighbouring(a);
            }
            capacity_type min_cap = std::min(cfc[a], max);
            
            dyn_info new_path;
            new_path.arc = max_arc; // this is fine as empty for source node
            new_path.back_node = i;
            new_path.cap = min_cap;
	    new_path.last_cap = std::max(last_cap_max, cfc[a]);
            new_path.neighbour_count = neighbour_min;
            
            unsigned int j = tsi[wc.target(a)];
            min_max[j].insert(std::make_pair(wc.id(a), new_path));
        }
    }
     
    
//    for (int j = 0; j <= i; j++) {
//        
//        logger::Instance()->debug("Backnode " + std::to_string(j) + " \n");
//        for (dyn_paths::iterator p_it = min_max[j].begin(); p_it != min_max[j].end(); ++p_it) { 
//            logger::Instance()->debug("E " + std::to_string(p_it->first) + " " + std::to_string(p_it->second.back_node) + " " + std::to_string(wc.id(p_it->second.arc)) + " " + std::to_string(p_it->second.cap) + "\n");
//        }
//    }
    
    // now for backtracking
    unsigned int index = top_order.size() - 1;
    capacity_type max = 0;
    capacity_type last_cap_max = 0;
    unsigned int neighbour_min = std::numeric_limits<unsigned int>::max();
    int best;
    // which path do we take?
    
//    logger::Instance()->debug("Backtracing start \n");
    for (dyn_paths::iterator p_it = min_max[index].begin(); p_it != min_max[index].end(); ++p_it) { 
        
//        logger::Instance()->debug("Test " + std::to_string(p_it->first) + " " + std::to_string(p_it->second.back_node) + " " + std::to_string(wc.id(p_it->second.arc)) + " " + std::to_string(p_it->second.cap) + "\n");
        if ( p_it->second.neighbour_count < neighbour_min 
                        || (p_it->second.neighbour_count == neighbour_min && p_it->second.cap > max)
                        || (p_it->second.neighbour_count == neighbour_min && p_it->second.cap == max && p_it->second.last_cap > last_cap_max) ) {
//            logger::Instance()->debug("Set \n");
            max = p_it->second.cap;
            best = p_it->first;
	    last_cap_max = p_it->second.last_cap;
            neighbour_min = p_it->second.neighbour_count;
        }
    }
    
//    logger::Instance()->debug("Backtracing  \n");
    while (index != 0) {

        int best_temp = best;
        path.push_front(wc.arcFromId(best));
        best = wc.id(min_max[index][best].arc);
        index = min_max[index][best_temp].back_node;
        
//         logger::Instance()->debug("back " + std::to_string(best) + " " + std::to_string(index) + " pushed " + std::to_string(best_temp) + "\n");
    }
    
    delete[] min_max;
}

unsigned int single_path_heuristic_no_neighbour::find_neighbouring(ListDigraph::Arc a) {
    
    unsigned int neighbours = 0;
    
    if (ces[a].id.size() == 0) {
	return 0;
    }

    unsigned int left; 
    unsigned int right = ces[a].id.find_first();
    
    boost::dynamic_bitset<>::size_type index = right;
    index = ces[a].id.find_next(index);
    if (index == boost::dynamic_bitset<>::npos) {
        return 0;
    }
    
    while(index != boost::dynamic_bitset<>::npos) {
        left = right;
        right = index;
        
        if (meta->exons[left].right + 1 == meta->exons[right].left) {
            ++neighbours;
        }
        
        index = ces[a].id.find_next(index);
    } 
        
    return neighbours;    
    
}
