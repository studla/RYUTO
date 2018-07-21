/* 
 * File:   multi_path_heuristic.cc
 * Author: Thomas Gatter <thomas(at)bioinf.uni-leipzig.de>
 * 
 * Created on June 8, 2016, 3:11 PM
 */

#include "multi_path_heuristic.h"

#include <lemon/lgf_writer.h>

multi_path_heuristic::multi_path_heuristic(ListDigraph& wc,
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
: path_finder(wc, s, t, cfc,mc, ces, cet, cel, cni,kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, size) {
}


multi_path_heuristic::~multi_path_heuristic() {
}

void multi_path_heuristic::extract_transcripts(alternative_transcript_collection& results) {

    alternative_transcript_collection first;
    unsigned int min_count = std::numeric_limits<unsigned int>::max();
    
    multi_path_recursion(wc, s, cfc, ces, cet, cel, cni, kp, kbp, cycle_id_in, cycle_id_out, unsecurityArc, unsecurityId, first, results, &min_count); 
    
}

void multi_path_heuristic::multi_path_recursion(ListDigraph &wc,  ListDigraph::Node &s,
        ListDigraph::ArcMap<capacity_type> &fc, ListDigraph::ArcMap<exon_edge> &ces,
        ListDigraph::ArcMap<edge_types::edge_type> &cet, ListDigraph::ArcMap<edge_length> &cel,
        ListDigraph::NodeMap<unsigned int> &cni,
        ListDigraph::ArcMap<arc_bridge> &know_paths, ListDigraph::ArcMap<arc_back_bridge> &know_back_paths,
        ListDigraph::ArcMap<unsigned int>& cycle_id_in, ListDigraph::ArcMap<unsigned int>& cycle_id_out,
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > &unsecurityArc,
        ListDigraph::NodeMap< unsecurity_id > &unsecurityId,
        alternative_transcript_collection &collection,
        alternative_transcript_collection& best, unsigned int* min_count) {
    
    #ifdef ALLOW_DEBUG
    logger::Instance()->debug("Enter recursion " + std::to_string(collection.size()) + " " + std::to_string(*min_count) + "\n");
    #endif
    
    if (collection.size() >= *min_count) { // branch and bound condition
        return;
    }
   
    ListDigraph::ArcIt has_arc(wc);
    if ( has_arc == INVALID ) {
        // no more acrs left in the graph, we are at the end!
        if (*min_count > collection.size()) {
            *min_count = collection.size();
            best = collection; // this only works thanks to deep copy
        }
        
        return;
    }
    
    std::deque<std::deque<ListDigraph::Arc> > paths;
    find_min_max(paths, wc, s, fc, know_paths);
    
    for (std::deque<std::deque<ListDigraph::Arc> >::iterator p_it = paths.begin(); p_it != paths.end(); ++p_it) {
        
        // for every path we need to have a fresh copy
        ListDigraph wcc; // working copy for modification
        ListDigraph::ArcMap<capacity_type> fcc(wcc);
        ListDigraph::ArcMap<exon_edge> cces(wcc); 
        ListDigraph::ArcMap<edge_types::edge_type> ccet(wcc);
        ListDigraph::ArcMap<edge_length> ccel(wcc);
        ListDigraph::NodeMap<unsigned int> ccni(wcc);
        ListDigraph::ArcMap<arc_bridge> ckp(wcc);
        ListDigraph::ArcMap<arc_back_bridge> ckbp(wcc);
        ListDigraph::ArcMap<unsigned int> ccidin(wcc);
        ListDigraph::ArcMap<unsigned int> ccidout(wcc);
        ListDigraph::ArcMap< lazy<std::set<transcript_unsecurity> > > cua(wcc);
        ListDigraph::NodeMap< unsecurity_id > cui(wcc);
        
        ListDigraph::Node cs;
        
        ListDigraph::ArcMap<ListDigraph::Arc> arc_ref(wc);
        ListDigraph::ArcMap<ListDigraph::Arc> arc_ref_rev(wc);

        DigraphCopy<ListDigraph, ListDigraph> copy(wc,wcc);
        copy.arcMap(fc, fcc).arcMap(cel, ccel)
            .arcMap(ces, cces).arcMap(cet, ccet)
            .arcMap(cycle_id_in,ccidin).arcMap(cycle_id_out,ccidout)
            .arcMap(unsecurityArc,cua).nodeMap(unsecurityId, cui)
            .arcMap(know_paths, ckp).arcMap(know_back_paths, ckbp)
            .nodeMap(cni, ccni).node(s,cs)
            .arcCrossRef(arc_ref);
        
        std::deque<ListDigraph::Arc> path;
        path.resize(p_it->size());
        int i = 0;
        for ( std::deque<ListDigraph::Arc>::iterator a = p_it->begin(); a!= p_it->end(); ++a, ++i) {
            copy.arc(*a, path[i]);
        }
        copy.run();

        reverse_arc_cross_map(wc, wcc, arc_ref, arc_ref_rev);
        // unfortunately we now also need to update all known and back_known, since IDs are not consistent
        update_known(ckp, wc, wcc, arc_ref_rev);
        update_known_back(ckbp, wc, wcc, arc_ref_rev);
        
//       logger::Instance()->debug("Copy \n");
//       digraphWriter(wc, std::cout)
//            .arcMap("edge_specifier", ces)
//            .arcMap("edge_type", cet)
//            .arcMap("flow", fc)
//            .node("source", s)
//            .run();  
        
        // we need to add new path to a collection copy
        alternative_transcript_collection cc = collection;
      //MC_FIX  capacity_type cap = add_path_to_collection(path, cc, wcc, fcc, cces, ccet, ckp, ccidin, ccidout, cua, cui);
        
        // now we need to remove this path from the copy
        //MC_FIX remove_full_path(path, cap, wcc, fcc, cces, ccet, ccel, ckp, ckbp,  ccidin, ccidout, cua, cui); // problem copying in KNOWN path and backpaths
 
        // start the recursion
        multi_path_recursion(wcc, cs, fcc, cces, ccet, ccel, ccni, ckp, ckbp, ccidin, ccidout, cua, cui, cc, best, min_count);
        
    }
}



void multi_path_heuristic::find_min_max( std::deque<std::deque<ListDigraph::Arc> > &path,
        ListDigraph& wc, ListDigraph::Node& s, ListDigraph::ArcMap<capacity_type>& cfc,
        ListDigraph::ArcMap<arc_bridge> &kp) {
    
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
        
        unsigned int j = tsi[wc.target(a)];
        add_pareto(min_max[j][wc.id(a)], se);
    }
    
    ++i;
    ++n;
    for (; n!= top_end; ++n, ++i ) {
                
        // for the rest we update pareto lists for each arc, so triple loop: Arcs To Add, for this Arc potential previous edges, elements in Pareto front this edge
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
                
                for (pareto::iterator ee = p_it->second.begin(); ee != p_it->second.end(); ++ee) {
                    evidence_element new_element = *ee; // deep copy
                    new_element.arcs.push_back(a);
                    new_element.brigde += brigde_count;
                    new_element.cap = std::min(new_element.cap, cfc[a]);
                    
                    unsigned int j = tsi[wc.target(a)];
                    add_pareto(min_max[j][wc.id(a)], new_element);
                    
                }
            }
        }
    }
    // in the end we merge all possible arcs to last node
    pareto end_list;
    for (dyn_paths::iterator p_it = min_max[i].begin(); p_it != min_max[i].end(); ++p_it) {
        for (pareto::iterator ee = p_it->second.begin(); ee != p_it->second.end(); ++ee) {
            add_pareto(end_list, *ee);
        }
    }
    for (pareto::iterator ee = end_list.begin(); ee != end_list.end(); ++ee) {
        path.push_back(ee->arcs);
    }
    
    // TODO: so far lists are copied a few times, but reference are not helpful as we extend intermediate list that need to stay intact... 
    
    delete[] min_max;
}

void multi_path_heuristic::add_pareto(pareto &pareto_list, evidence_element& insert) {
        
    for (pareto::iterator el = pareto_list.begin(); el != pareto_list.end(); ) {
        if (insert.cap == el->cap && insert.brigde == el->brigde) {
            // if this happens we know we can insert as we have a copy to an actual value
            pareto_list.push_back(insert);
            return;
        } else if (insert.cap >= el->cap && insert.brigde >= el->brigde) {
            // so the new inserted element (co-optimal retaining) dominates this one
            el = pareto_list.erase(el);
        } else if (insert.cap <= el->cap && insert.brigde <= el->brigde) {
            // new element is dominated by existing
            return;
        } else {
            ++el;
        }
    }
    pareto_list.push_back(insert);
}